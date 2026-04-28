"""
test_multidonor_2pass.py

Three-way comparison on the same lattice/bracket/budget:
  (A) 3-pass GC scan      — trusted ground truth
  (B) multi-donor (1-pass) — uniform tile across the bracket
  (C) multi-donor (2-pass) — coarse Pass-1 → narrow Pass-2 around peak

Saves a 2-panel figure (chi(beta) curves + wall-time bars) to
results/_test_multidonor_2pass/comparison.png.

Run from K_from_Optimizer_Production/:
    python test_multidonor_2pass.py
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
import reweight as rw_mod
from run import ensure_binary, CONFIG

LX, LY = 12, 12
TX, TY = 0, 0
K1 = K2 = K3 = 1.0
BETA_LO, BETA_HI = 0.18, 0.40

N_TRAJ_BUDGET = 30_000

OUT_DIR = os.path.join(_HERE, "results", "_test_multidonor_2pass")
os.makedirs(OUT_DIR, exist_ok=True)


def _collect_donor_betas(root, prefix="donor_"):
    """Walk a multidonor scratch dir and return parent betas of each donor."""
    out = []
    if not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root)):
        if not name.startswith(prefix):
            continue
        sub = os.path.join(root, name)
        if not os.path.isdir(sub):
            continue
        for child in sorted(os.listdir(sub)):
            cp = os.path.join(sub, child)
            if not os.path.isdir(cp):
                continue
            tp = rw_mod.find_traces_file(cp)
            if tp is None:
                continue
            try:
                meta = rw_mod.parse_traces(tp)
                out.append(float(meta["beta_parent"]))
            except Exception:
                pass
            break
    return sorted(out)


def main():
    exe = os.path.join(_HERE, CONFIG["exe"])
    ensure_binary(CONFIG["exe"])
    print(f"[test] geometry  : Lx={LX} Ly={LY} Tx={TX} Ty={TY}")
    print(f"[test] bracket   : [{BETA_LO}, {BETA_HI}]")
    print(f"[test] traj budget per method: {N_TRAJ_BUDGET}")
    print(f"[test] output dir: {OUT_DIR}")
    print()

    # ------------------------------------------------------------------
    # (A) 3-pass GC scan
    # ------------------------------------------------------------------
    n_coarse, n_refine, n_refine2, n_refine3 = 11, 5, 5, 5
    n_fine_pts = n_refine + n_refine2 + n_refine3
    n_traj_coarse = max(2000, N_TRAJ_BUDGET // (3 * n_coarse))
    n_traj_fine   = max(2000, (2 * N_TRAJ_BUDGET) // (3 * n_fine_pts))
    print("[test] === 3-pass GC scan ===")
    print(f"[3p ] n_traj_coarse={n_traj_coarse}, n_traj_fine={n_traj_fine}")
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
    sp_betas = np.asarray(sp_betas); sp_chis = np.asarray(sp_chis); sp_errs = np.asarray(sp_errs)
    print(f"[3p ] β_c = {sp_bc:.6f} ± {sp_bc_sig:.2e}   wall = {sp_wall:.1f}s")

    # ------------------------------------------------------------------
    # (B) multi-donor 1-pass
    # ------------------------------------------------------------------
    print("\n[test] === Multi-donor (1-pass) ===")
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
    md_betas = np.asarray(md_betas); md_chis = np.asarray(md_chis); md_errs = np.asarray(md_errs)
    md_donor_betas = _collect_donor_betas(os.path.join(OUT_DIR, "md_scratch"),
                                          prefix="donor_")
    print(f"[md ] β_c = {md_bc:.6f} ± {md_bc_sig:.2e}   wall = {md_wall:.1f}s   "
          f"donors={len(md_donor_betas)}")

    # ------------------------------------------------------------------
    # (C) multi-donor 2-pass
    # ------------------------------------------------------------------
    print("\n[test] === Multi-donor (2-pass) ===")
    t0 = time.time()
    md2_result = mc_engine.find_beta_c_multidonor_2pass(
        exe, LX, LY, TX, TY, K1, K2, K3,
        BETA_LO, BETA_HI,
        n_traj_parent=N_TRAJ_BUDGET,
        n_grid=201,
        n_eff_floor=0.10,
        pass1_n_donors=4,
        pass1_budget_frac=0.30,
        donor_overlap_alpha=0.75,
        max_donors=12,
        data_dir=os.path.join(OUT_DIR, "md2_scratch"),
        jackknife=True,
    )
    md2_wall = time.time() - t0
    md2_bc, md2_bc_sig, md2_chi_peak, md2_betas, md2_chis, md2_errs = md2_result
    md2_betas = np.asarray(md2_betas); md2_chis = np.asarray(md2_chis); md2_errs = np.asarray(md2_errs)
    md2_p1 = _collect_donor_betas(os.path.join(OUT_DIR, "md2_scratch"), prefix="p1_donor_")
    md2_p2 = _collect_donor_betas(os.path.join(OUT_DIR, "md2_scratch"), prefix="p2_donor_")
    print(f"[md2] β_c = {md2_bc:.6f} ± {md2_bc_sig:.2e}   wall = {md2_wall:.1f}s   "
          f"donors=p1:{len(md2_p1)} + p2:{len(md2_p2)}")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, (ax_chi, ax_bar) = plt.subplots(
        1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [3, 1]},
    )

    # 3-pass scan points.
    ax_chi.errorbar(sp_betas, sp_chis, yerr=sp_errs,
                    fmt="s", ms=5, color="C3", lw=1, capsize=2,
                    label=f"3-pass GC ({len(sp_betas)} pts)")

    # 1-pass MD curve.
    md_v = np.isfinite(md_chis)
    ax_chi.fill_between(md_betas[md_v], (md_chis - md_errs)[md_v], (md_chis + md_errs)[md_v],
                        alpha=0.15, color="C0")
    ax_chi.plot(md_betas[md_v], md_chis[md_v], "-", color="C0", lw=2,
                label=f"multi-donor 1-pass ({len(md_donor_betas)} donors)")

    # 2-pass MD curve.
    md2_v = np.isfinite(md2_chis)
    ax_chi.fill_between(md2_betas[md2_v], (md2_chis - md2_errs)[md2_v],
                        (md2_chis + md2_errs)[md2_v], alpha=0.15, color="C2")
    ax_chi.plot(md2_betas[md2_v], md2_chis[md2_v], "-", color="C2", lw=2,
                label=f"multi-donor 2-pass (p1:{len(md2_p1)}+p2:{len(md2_p2)})")

    # β_c lines.
    ax_chi.axvline(sp_bc, color="C3", ls=":",  lw=1.4,
                   label=f"β_c (3p)  = {sp_bc:.5f}±{sp_bc_sig:.1e}")
    ax_chi.axvline(md_bc, color="C0", ls="--", lw=1.4,
                   label=f"β_c (1p)  = {md_bc:.5f}±{md_bc_sig:.1e}")
    ax_chi.axvline(md2_bc, color="C2", ls="-.",  lw=1.4,
                   label=f"β_c (2p)  = {md2_bc:.5f}±{md2_bc_sig:.1e}")

    # Donor rugs.
    ymin, ymax = ax_chi.get_ylim()
    rug_y0 = ymin + 0.02 * (ymax - ymin)
    rug_y1 = ymin + 0.05 * (ymax - ymin)
    rug_y2 = ymin + 0.08 * (ymax - ymin)
    if md_donor_betas:
        ax_chi.plot(md_donor_betas, [rug_y0] * len(md_donor_betas),
                    "v", color="C0", ms=7, alpha=0.7)
    if md2_p1:
        ax_chi.plot(md2_p1, [rug_y1] * len(md2_p1),
                    "v", color="C2", ms=7, alpha=0.4, label="2-pass: pass-1 donors")
    if md2_p2:
        ax_chi.plot(md2_p2, [rug_y2] * len(md2_p2),
                    "^", color="C2", ms=7, alpha=0.9, label="2-pass: pass-2 donors")

    ax_chi.set_xlabel("β")
    ax_chi.set_ylabel("χ(β)")
    ax_chi.set_title(
        f"Multi-donor 1-pass vs 2-pass vs 3-pass GC\n"
        f"L={LX}×{LY}, K=({K1},{K2},{K3}), budget={N_TRAJ_BUDGET} traj/method"
    )
    ax_chi.legend(loc="upper left", fontsize=8)
    ax_chi.grid(True, alpha=0.3)

    methods = ["3-pass GC", "1-pass MD", "2-pass MD"]
    walls   = [sp_wall, md_wall, md2_wall]
    colors  = ["C3", "C0", "C2"]
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

    print()
    print("=" * 72)
    print(f"  method        β_c              σ(β_c)         χ_peak     wall")
    print("-" * 72)
    print(f"  3-pass GC    {sp_bc:.6f}    {sp_bc_sig:.2e}    {sp_chi_peak:.4g}    {sp_wall:5.1f}s")
    print(f"  1-pass MD    {md_bc:.6f}    {md_bc_sig:.2e}    {md_chi_peak:.4g}    {md_wall:5.1f}s")
    print(f"  2-pass MD    {md2_bc:.6f}    {md2_bc_sig:.2e}    {md2_chi_peak:.4g}    {md2_wall:5.1f}s")
    print("-" * 72)
    sig = max(sp_bc_sig, 1e-12)
    print(f"  Δβ_c (1p − 3p) = {md_bc - sp_bc:+.6f}    ({abs(md_bc - sp_bc) / sig:.2f} σ_3p)")
    print(f"  Δβ_c (2p − 3p) = {md2_bc - sp_bc:+.6f}    ({abs(md2_bc - sp_bc) / sig:.2f} σ_3p)")
    print(f"  σ(2p) / σ(1p)  = {md2_bc_sig / max(md_bc_sig, 1e-12):.3f}")
    print(f"  speedup 3p/1p  = {sp_wall / max(md_wall, 1e-9):.2f}×")
    print(f"  speedup 3p/2p  = {sp_wall / max(md2_wall, 1e-9):.2f}×")
    print("=" * 72)


if __name__ == "__main__":
    main()
