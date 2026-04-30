#!/usr/bin/env python3
"""
compare_betac_stability.py

Side-by-side stability test of:
  (A) 3-pass Gram-Charlier scan  — mc_engine.find_beta_c
  (B) Single-donor reweighting   — mc_engine.find_beta_c_reweight

Runs N_TRIALS independent trials of each method on each test geometry/coupling,
then reports bias, std(β_c), dome fraction, and wall time.

Usage:
    python compare_betac_stability.py [--n-trials N] [--no-vis]

All scratch data written to results/_betac_stability_test/ and never deleted.
"""
from __future__ import annotations

import argparse
import json
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

# ---------------------------------------------------------------------------
# Test cases: (label, Lx, Ly, Tx, Ty, k1, k2, k3, beta_lo, beta_hi, beta_exact)
# beta_exact = None when unknown; used only for bias computation.
# ---------------------------------------------------------------------------
# equilateral isotropic: exact β_c = ln(3)/4 ≈ 0.27465
_LN3_4 = 0.25 * np.log(3.0)

TEST_CASES = [
    # Near-equilateral — easy, reweighting should work fine here
    dict(label="equil_12",
         Lx=12, Ly=12, Tx=0, Ty=0, k1=1.0, k2=1.0, k3=1.0,
         beta_lo=0.20, beta_hi=0.35, beta_exact=_LN3_4),
    # Moderate anisotropy — typical early optimizer iterate
    dict(label="aniso_r2_12",
         Lx=12, Ly=12, Tx=0, Ty=0, k1=2.0, k2=1.0, k3=1.0,
         beta_lo=0.10, beta_hi=0.30, beta_exact=None),
    # High anisotropy — like k=(5,7) region the 4-5-6 optimizer explores
    dict(label="aniso_r5_12",
         Lx=12, Ly=12, Tx=0, Ty=0, k1=5.0, k2=1.0, k3=1.0,
         beta_lo=0.05, beta_hi=0.25, beta_exact=None),
]

# MC budget per trial (keep modest so the full test runs in ~10 min)
N_TRAJ_SCAN_COARSE = 4_000
N_TRAJ_SCAN_FINE   = 8_000
N_TRAJ_REWEIGHT    = 30_000   # single parent run for reweight
N_GRID             = 201
N_EFF_FLOOR        = 0.10
MAX_RECENTERS      = 4

SCAN_N_COARSE  = 11
SCAN_N_REFINE  = 5
SCAN_N_REFINE2 = 5
SCAN_N_REFINE3 = 5
SCAN_MAX_SHIFTS = 4


def _has_dome(scan_betas, scan_chis):
    """Return True if the chi peak is interior (flanked on both sides)."""
    chi = np.array(scan_chis, dtype=float)
    valid = np.isfinite(chi)
    if valid.sum() < 3:
        return False
    peak_idx = int(np.argmax(chi[valid]))
    return 0 < peak_idx < valid.sum() - 1


def run_trial_3pass(exe, case, scratch_dir, trial_idx):
    """Single trial using the 3-pass GC scan."""
    t0 = time.time()
    try:
        bc, bc_sig, chi_pk, sb, sc, sce = mc_engine.find_beta_c(
            exe, case["Lx"], case["Ly"], case["Tx"], case["Ty"],
            case["k1"], case["k2"], case["k3"],
            case["beta_lo"], case["beta_hi"],
            n_coarse=SCAN_N_COARSE,
            n_refine=SCAN_N_REFINE, n_refine2=SCAN_N_REFINE2,
            n_refine3=SCAN_N_REFINE3,
            n_traj_coarse=N_TRAJ_SCAN_COARSE,
            n_traj_fine=N_TRAJ_SCAN_FINE,
            max_shifts=SCAN_MAX_SHIFTS,
            data_dir=os.path.join(scratch_dir, f"3pass_t{trial_idx:02d}"),
            jackknife=True,
        )
        wall = time.time() - t0
        dome = _has_dome(sb, sc)
        return dict(ok=True, bc=bc, bc_sig=bc_sig, chi_pk=chi_pk,
                    wall=wall, dome=dome, n_valid=sum(np.isfinite(sc)))
    except Exception as exc:
        wall = time.time() - t0
        return dict(ok=False, bc=float("nan"), bc_sig=0.0, chi_pk=float("nan"),
                    wall=wall, dome=False, n_valid=0, error=str(exc))


def run_trial_reweight(exe, case, scratch_dir, trial_idx):
    """Single trial using single-donor Ferrenberg-Swendsen reweighting."""
    t0 = time.time()
    try:
        bc, bc_sig, chi_pk, sb, sc, sce = mc_engine.find_beta_c_reweight(
            exe, case["Lx"], case["Ly"], case["Tx"], case["Ty"],
            case["k1"], case["k2"], case["k3"],
            case["beta_lo"], case["beta_hi"],
            n_traj_parent=N_TRAJ_REWEIGHT,
            n_grid=N_GRID,
            n_eff_floor=N_EFF_FLOOR,
            max_recenters=MAX_RECENTERS,
            data_dir=os.path.join(scratch_dir, f"rw_t{trial_idx:02d}"),
            jackknife=True,
        )
        wall = time.time() - t0
        dome = _has_dome(sb, sc)
        return dict(ok=True, bc=bc, bc_sig=bc_sig, chi_pk=chi_pk,
                    wall=wall, dome=dome, n_valid=sum(np.isfinite(sc)))
    except Exception as exc:
        wall = time.time() - t0
        return dict(ok=False, bc=float("nan"), bc_sig=0.0, chi_pk=float("nan"),
                    wall=wall, dome=False, n_valid=0, error=str(exc))


def summarise(results, beta_exact=None):
    ok = [r for r in results if r["ok"]]
    n_ok = len(ok)
    n_tot = len(results)
    if n_ok == 0:
        return dict(n_ok=0, n_tot=n_tot, mean=float("nan"), std=float("nan"),
                    bias=float("nan"), dome_frac=0.0, wall_mean=float("nan"),
                    wall_std=float("nan"))
    bcs = np.array([r["bc"] for r in ok])
    walls = np.array([r["wall"] for r in ok])
    dome_frac = sum(r["dome"] for r in ok) / n_ok
    bias = float(np.mean(bcs) - beta_exact) if beta_exact is not None else float("nan")
    return dict(n_ok=n_ok, n_tot=n_tot,
                mean=float(np.mean(bcs)), std=float(np.std(bcs, ddof=1) if n_ok > 1 else 0.0),
                bias=bias, dome_frac=dome_frac,
                wall_mean=float(np.mean(walls)), wall_std=float(np.std(walls, ddof=1) if n_ok > 1 else 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-trials", type=int, default=5)
    ap.add_argument("--no-vis", action="store_true")
    args = ap.parse_args()

    N_TRIALS = args.n_trials

    ensure_binary(CONFIG["exe"])
    exe = os.path.join(_HERE, CONFIG["exe"])

    out_dir = os.path.join(_HERE, "results", "_betac_stability_test")
    os.makedirs(out_dir, exist_ok=True)

    all_results = {}

    for case in TEST_CASES:
        label = case["label"]
        scratch = os.path.join(out_dir, label)
        os.makedirs(scratch, exist_ok=True)
        print(f"\n{'='*60}")
        print(f"CASE: {label}  k=({case['k1']},{case['k2']},{case['k3']})  "
              f"β∈[{case['beta_lo']},{case['beta_hi']}]  "
              f"exact={'%.5f'%case['beta_exact'] if case['beta_exact'] else '?'}")
        print(f"{'='*60}")

        res3p = []
        resrw = []

        for trial in range(N_TRIALS):
            print(f"\n--- Trial {trial+1}/{N_TRIALS} ---")

            print(f"  [3-pass] ", end="", flush=True)
            r3 = run_trial_3pass(exe, case, scratch, trial)
            status = f"β_c={r3['bc']:.5f}±{r3['bc_sig']:.2e}  dome={r3['dome']}  {r3['wall']:.1f}s" if r3["ok"] else f"FAILED: {r3.get('error','')}"
            print(status)
            res3p.append(r3)

            print(f"  [reweight]", end="", flush=True)
            rr = run_trial_reweight(exe, case, scratch, trial)
            status = f"β_c={rr['bc']:.5f}±{rr['bc_sig']:.2e}  dome={rr['dome']}  n_valid={rr['n_valid']}/201  {rr['wall']:.1f}s" if rr["ok"] else f"FAILED: {rr.get('error','')}"
            print(status)
            resrw.append(rr)

        sum3p = summarise(res3p, case["beta_exact"])
        sumrw = summarise(resrw, case["beta_exact"])

        all_results[label] = {"3pass": sum3p, "reweight": sumrw,
                              "trials_3pass": res3p, "trials_rw": resrw,
                              "case": {k: v for k, v in case.items()}}

        print(f"\nSUMMARY [{label}]  ({N_TRIALS} trials each)")
        print(f"  {'Method':<12}  {'mean β_c':>10}  {'std':>8}  {'bias':>8}  "
              f"{'dome%':>6}  {'wall/trial':>10}  {'ok/n':>6}")
        for meth, s in [("3-pass GC", sum3p), ("reweight", sumrw)]:
            print(f"  {meth:<12}  {s['mean']:>10.5f}  {s['std']:>8.5f}  "
                  f"{s['bias']:>+8.5f}  {100*s['dome_frac']:>5.0f}%  "
                  f"{s['wall_mean']:>8.1f}s  {s['n_ok']}/{s['n_tot']}")

    # Save JSON
    json_path = os.path.join(out_dir, "results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved → {json_path}")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    if args.no_vis:
        return

    n_cases = len(TEST_CASES)
    fig, axes = plt.subplots(2, n_cases, figsize=(5 * n_cases, 8),
                             squeeze=False)

    for col, case in enumerate(TEST_CASES):
        label = case["label"]
        data = all_results[label]
        res3 = data["trials_3pass"]
        resrw = data["trials_rw"]

        # Top row: β_c per trial
        ax = axes[0, col]
        t_idx = list(range(1, len(res3) + 1))
        bcs3 = [r["bc"] for r in res3]
        bcsrw = [r["bc"] for r in resrw]
        ax.plot(t_idx, bcs3,  "s-", color="C3", label="3-pass GC", ms=7)
        ax.plot(t_idx, bcsrw, "o-", color="C0", label="reweight",  ms=7)
        if case["beta_exact"] is not None:
            ax.axhline(case["beta_exact"], ls="--", color="k", lw=1, label=f"exact={case['beta_exact']:.4f}")

        s3 = data["3pass"]
        sr = data["reweight"]
        ax.axhline(s3["mean"],  ls=":", color="C3", lw=1)
        ax.axhline(sr["mean"], ls=":", color="C0", lw=1)
        ax.set_title(f"{label}\nk=({case['k1']},{case['k2']},{case['k3']})", fontsize=10)
        ax.set_xlabel("trial")
        ax.set_ylabel("β_c")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Bottom row: dome fraction + wall time bar chart
        ax2 = axes[1, col]
        methods = ["3-pass GC", "reweight"]
        dome_fracs = [100 * s3["dome_frac"], 100 * sr["dome_frac"]]
        walls = [s3["wall_mean"], sr["wall_mean"]]
        x = np.array([0, 1])
        width = 0.35
        bars1 = ax2.bar(x - width/2, dome_fracs, width, label="dome %", color=["C3","C0"])
        ax2b = ax2.twinx()
        bars2 = ax2b.bar(x + width/2, walls, width, label="wall time (s)",
                         color=["C3","C0"], alpha=0.5, hatch="//")
        ax2.set_xticks(x); ax2.set_xticklabels(methods, fontsize=9)
        ax2.set_ylabel("dome %")
        ax2b.set_ylabel("wall time (s)")
        ax2.set_ylim(0, 110)
        ax2.set_title("dome fraction & wall time", fontsize=9)
        ax2.grid(True, alpha=0.2, axis="y")

    fig.suptitle(
        f"β_c stability: 3-pass GC scan vs single-donor reweighting\n"
        f"({N_TRIALS} trials each, L=12, Nt_3pass≈{N_TRAJ_SCAN_COARSE}c/{N_TRAJ_SCAN_FINE}f, "
        f"Nt_rw={N_TRAJ_REWEIGHT})",
        fontsize=11,
    )
    fig.tight_layout()
    png_path = os.path.join(out_dir, "betac_stability.png")
    fig.savefig(png_path, dpi=130)
    print(f"Plot saved → {png_path}")


if __name__ == "__main__":
    main()
