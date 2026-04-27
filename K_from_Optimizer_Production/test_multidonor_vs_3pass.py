"""
test_multidonor_vs_3pass.py

Side-by-side comparison of:
  (A) the new multi-donor reweighting:  mc_engine.find_beta_c_multidonor
  (B) the trusted 3-pass GC scan:       mc_engine.find_beta_c

on the same lattice geometry and bracket.  Prints β_c estimates and
saves a 2-panel figure:

  Left  — χ(β) curves: 3-pass scan points (with jackknife errors) +
          multi-donor reweighted curve (continuous, with σ band) +
          donor positions and the GC peak fits from each method.
  Right — wall-time bar comparison.

Run from K_from_Optimizer_Production/:
    python test_multidonor_vs_3pass.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from run import ensure_binary, CONFIG

# ----------------------------------------------------------------------
# Geometry + bracket — small enough that the 3-pass scan finishes in
# under a couple of minutes, large enough that the χ(β) peak is sharp.
# ----------------------------------------------------------------------
LX, LY = 12, 12
TX, TY = 0, 0
K1 = K2 = K3 = 1.0
BETA_LO, BETA_HI = 0.18, 0.40   # straddles 2D Ising criticality on a triangle

# Equal total-trajectory budget so the comparison is fair.
N_TRAJ_BUDGET = 30_000

OUT_DIR = os.path.join(_HERE, "results", "_test_multidonor_vs_3pass")
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    exe = os.path.join(_HERE, CONFIG["exe"])
    ensure_binary(CONFIG["exe"])
    print(f"[test] geometry  : Lx={LX} Ly={LY} Tx={TX} Ty={TY}")
    print(f"[test] bracket   : [{BETA_LO}, {BETA_HI}]")
    print(f"[test] traj budget per method: {N_TRAJ_BUDGET}")
    print(f"[test] output dir: {OUT_DIR}")
    print()

    # ------------------------------------------------------------------
    # (A) Multi-donor reweighting
    # ------------------------------------------------------------------
    print("[test] === Multi-donor reweighting ===")
    t0 = time.time()
    md_result = mc_engine.find_beta_c_multidonor(
        exe, LX, LY, TX, TY, K1, K2, K3,
        BETA_LO, BETA_HI,
        n_traj_parent=N_TRAJ_BUDGET,
        n_grid=201,
        n_eff_floor=0.10,
        donor_overlap_alpha=0.75,
        pilot_n_traj=2000,
        max_donors=12,
        data_dir=os.path.join(OUT_DIR, "md_scratch"),
        jackknife=True,
    )
    md_wall = time.time() - t0
    md_bc, md_bc_sig, md_chi_peak, md_betas, md_chis, md_errs = md_result
    md_betas = np.asarray(md_betas)
    md_chis  = np.asarray(md_chis)
    md_errs  = np.asarray(md_errs)
    print(f"[md ] β_c = {md_bc:.6f} ± {md_bc_sig:.2e}   "
          f"χ_peak = {md_chi_peak:.4g}   wall = {md_wall:.1f}s")

    # Recover donor β values from scratch dirs (donor_00, donor_01, ...).
    import reweight as rw_mod
    donor_betas = []
    md_root = os.path.join(OUT_DIR, "md_scratch")
    for name in sorted(os.listdir(md_root)):
        if not name.startswith("donor_"):
            continue
        sub = os.path.join(md_root, name)
        for child in sorted(os.listdir(sub)):
            child_path = os.path.join(sub, child)
            if not os.path.isdir(child_path):
                continue
            tp = rw_mod.find_traces_file(child_path)
            if tp is None:
                continue
            try:
                meta = rw_mod.parse_traces(tp)
                donor_betas.append(float(meta["beta_parent"]))
            except Exception:
                pass
            break
    donor_betas.sort()
    print(f"[md ] {len(donor_betas)} donors at β = "
          + ", ".join(f"{b:.4f}" for b in donor_betas))

    # ------------------------------------------------------------------
    # (B) 3-pass GC scan (the ground truth)
    # ------------------------------------------------------------------
    # Match the total trajectory budget: the 3-pass scan does
    # (n_coarse * n_traj_coarse + (n_refine + n_refine2 + n_refine3) *
    #  n_traj_fine).  Keep the default refine counts and just scale
    # n_traj_* so the totals are comparable.
    n_coarse, n_refine, n_refine2, n_refine3 = 11, 5, 5, 5
    n_fine_pts = n_refine + n_refine2 + n_refine3   # 15
    # Allocate ~⅓ of the budget to coarse, ~⅔ to fine.
    n_traj_coarse = max(2000, N_TRAJ_BUDGET // (3 * n_coarse))
    n_traj_fine   = max(2000, (2 * N_TRAJ_BUDGET) // (3 * n_fine_pts))
    total_3p = n_coarse * n_traj_coarse + n_fine_pts * n_traj_fine
    print(f"\n[test] === 3-pass GC scan ===")
    print(f"[3p ] n_traj_coarse={n_traj_coarse}, n_traj_fine={n_traj_fine}, "
          f"total≈{total_3p}")
    t0 = time.time()
    sp_result = mc_engine.find_beta_c(
        exe, LX, LY, TX, TY, K1, K2, K3,
        BETA_LO, BETA_HI,
        n_coarse=n_coarse, n_refine=n_refine,
        n_refine2=n_refine2, n_refine3=n_refine3,
        n_traj_coarse=n_traj_coarse, n_traj_fine=n_traj_fine,
        data_dir=os.path.join(OUT_DIR, "scan_scratch"),
        jackknife=True,
    )
    sp_wall = time.time() - t0
    sp_bc, sp_bc_sig, sp_chi_peak, sp_betas, sp_chis, sp_errs = sp_result
    sp_betas = np.asarray(sp_betas)
    sp_chis  = np.asarray(sp_chis)
    sp_errs  = np.asarray(sp_errs)
    print(f"[3p ] β_c = {sp_bc:.6f} ± {sp_bc_sig:.2e}   "
          f"χ_peak = {sp_chi_peak:.4g}   wall = {sp_wall:.1f}s")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, (ax_chi, ax_bar) = plt.subplots(
        1, 2, figsize=(14, 5.5),
        gridspec_kw={"width_ratios": [3, 1]},
    )

    # Left: χ(β)
    md_valid = np.isfinite(md_chis)
    ax_chi.fill_between(
        md_betas[md_valid],
        (md_chis - md_errs)[md_valid],
        (md_chis + md_errs)[md_valid],
        alpha=0.20, color="C0", label="multi-donor σ band",
    )
    ax_chi.plot(
        md_betas[md_valid], md_chis[md_valid],
        "-", color="C0", lw=2, label=f"multi-donor χ(β) ({sum(np.isfinite(md_chis))} pts)",
    )
    ax_chi.errorbar(
        sp_betas, sp_chis, yerr=sp_errs,
        fmt="s", ms=5, color="C3", lw=1, capsize=2,
        label=f"3-pass GC scan ({len(sp_betas)} pts)",
    )

    # β_c estimates as vertical lines.
    ax_chi.axvline(md_bc, color="C0", ls="--", lw=1.5,
                   label=f"β_c (MD) = {md_bc:.5f}±{md_bc_sig:.1e}")
    ax_chi.axvline(sp_bc, color="C3", ls=":",  lw=1.5,
                   label=f"β_c (3p) = {sp_bc:.5f}±{sp_bc_sig:.1e}")

    # Donor β positions (rug).
    if donor_betas:
        ymin, ymax = ax_chi.get_ylim()
        rug_y = ymin + 0.02 * (ymax - ymin)
        ax_chi.plot(donor_betas, [rug_y] * len(donor_betas),
                    "v", color="C0", ms=8, alpha=0.7,
                    label=f"{len(donor_betas)} donor positions")

    ax_chi.set_xlabel("β")
    ax_chi.set_ylabel("χ(β)")
    ax_chi.set_title(
        f"Multi-donor reweight vs. 3-pass GC scan\n"
        f"L_x={LX} L_y={LY} T_x={TX} T_y={TY}, "
        f"K=({K1},{K2},{K3}), budget={N_TRAJ_BUDGET} traj/method"
    )
    ax_chi.legend(loc="upper left", fontsize=9)
    ax_chi.grid(True, alpha=0.3)

    # Right: wall-time bar.
    methods = ["3-pass GC", f"multi-donor\n({len(donor_betas)} donors)"]
    walls   = [sp_wall, md_wall]
    colors  = ["C3", "C0"]
    bars = ax_bar.bar(methods, walls, color=colors, edgecolor="k")
    for b, w in zip(bars, walls):
        ax_bar.text(b.get_x() + b.get_width() / 2, w,
                    f"{w:.1f}s", ha="center", va="bottom", fontsize=10)
    ax_bar.set_ylabel("wall time (s)")
    ax_bar.set_title("Wall time")
    ax_bar.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    out_png = os.path.join(OUT_DIR, "comparison.png")
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"\n[test] wrote: {out_png}")

    # Summary table.
    delta_bc = md_bc - sp_bc
    print()
    print("=" * 64)
    print(f"  method        β_c              σ(β_c)         χ_peak     wall")
    print("-" * 64)
    print(f"  3-pass GC    {sp_bc:.6f}    {sp_bc_sig:.2e}    "
          f"{sp_chi_peak:.4g}    {sp_wall:5.1f}s")
    print(f"  multi-donor  {md_bc:.6f}    {md_bc_sig:.2e}    "
          f"{md_chi_peak:.4g}    {md_wall:5.1f}s")
    print("-" * 64)
    print(f"  Δβ_c (MD − 3p)  =  {delta_bc:+.6f}    "
          f"({abs(delta_bc) / max(sp_bc_sig, 1e-12):.2f} σ_3p)")
    print(f"  speedup  3p / MD =  {sp_wall / max(md_wall, 1e-9):.2f}×")
    print("=" * 64)


if __name__ == "__main__":
    main()
