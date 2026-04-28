#!/usr/bin/env python3
"""Stress test of the β_c finder across anisotropic couplings.

For the 2-D triangular Ising model the exact bulk critical β satisfies

    sinh(2β k1)sinh(2β k2) + sinh(2β k2)sinh(2β k3)
        + sinh(2β k3)sinh(2β k1) = 1                        (*)

Strategy
--------
1.  Compute the exact bulk β_c(k1, k2, k3) by bisecting (*).
2.  Run the production β_c finder on a fixed lattice (default L=20,
    Tx=Ty=0) for a small grid of (k1, k2, k3) triples chosen to span
    weak-to-moderate anisotropy plus the isotropic baseline.
3.  Compare β_c(L) against β_c(∞) using the FSS calibration determined
    by the isotropic FSS run (results/_fss_betac/fss_summary.json).

A clean test passes if the deviation
    Δ(L) = β_c(L) − β_c(∞) − Δ_iso(L)
is consistent with zero within the reported jackknife σ for every
anisotropic point.  Δ_iso(L) is the shift measured at the same L on
the isotropic case.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine  # noqa: E402


def exact_beta_c(k1: float, k2: float, k3: float,
                 lo: float = 1e-3, hi: float = 5.0) -> float:
    """Bisect (*) for β > 0."""
    def f(b):
        return (math.sinh(2 * b * k1) * math.sinh(2 * b * k2)
                + math.sinh(2 * b * k2) * math.sinh(2 * b * k3)
                + math.sinh(2 * b * k3) * math.sinh(2 * b * k1) - 1.0)
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise RuntimeError(f"no sign change for k=({k1},{k2},{k3})")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm == 0.0 or 0.5 * (hi - lo) < 1e-14:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# Cases: (label, k1, k2, k3)
DEFAULT_CASES = [
    ("iso",        1.0, 1.0, 1.0),
    ("aniso_1.2",  1.2, 1.0, 0.8),
    ("aniso_1.5",  1.5, 1.0, 0.5),
    ("aniso_2.0",  2.0, 1.0, 0.5),
    ("two_eq_0.7", 1.0, 1.0, 0.7),
    ("two_eq_1.5", 1.0, 1.0, 1.5),
]


def run_one(L: int, k1: float, k2: float, k3: float,
            *, beta_lo_factor: float, beta_hi_factor: float,
            n_traj_coarse: int, n_traj_fine: int,
            n_coarse: int, n_refine: int,
            jackknife: bool, scratch_root: str) -> dict:
    bc_exact = exact_beta_c(k1, k2, k3)
    beta_lo = max(0.05, beta_lo_factor * bc_exact)
    beta_hi = beta_hi_factor * bc_exact
    tag = f"k{k1:.2f}_{k2:.2f}_{k3:.2f}_L{L}"
    data_dir = os.path.join(scratch_root, tag)
    print(f"[{tag}] β_c_exact={bc_exact:.6f}; "
          f"scan β∈[{beta_lo:.4f},{beta_hi:.4f}]", flush=True)
    t0 = time.time()
    bc, sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
        os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram"),
        L, L, 0, 0, k1, k2, k3,
        beta_lo, beta_hi,
        n_coarse=n_coarse,
        n_refine=n_refine, n_refine2=n_refine, n_refine3=n_refine,
        n_traj_coarse=n_traj_coarse, n_traj_fine=n_traj_fine,
        max_shifts=4, jackknife=jackknife,
        data_dir=data_dir,
    )
    dt = time.time() - t0
    print(f"[{tag}]  β_c(L)={bc:.6f} ± {sigma:.2e}  "
          f"χ_peak={chi_peak:.3g}  ({dt:.1f}s)", flush=True)
    return {
        "L": L, "k1": k1, "k2": k2, "k3": k3,
        "beta_c_exact": float(bc_exact),
        "beta_c": float(bc), "beta_c_sigma": float(sigma),
        "chi_peak": float(chi_peak),
        "wall_s": dt,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--L", type=int, default=20)
    ap.add_argument("--n-traj-coarse", type=int, default=4000)
    ap.add_argument("--n-traj-fine", type=int, default=10000)
    ap.add_argument("--n-coarse", type=int, default=11)
    ap.add_argument("--n-refine", type=int, default=5)
    ap.add_argument("--beta-lo-factor", type=float, default=0.55)
    ap.add_argument("--beta-hi-factor", type=float, default=1.45)
    ap.add_argument("--no-jackknife", action="store_true")
    ap.add_argument("--out-dir", type=str,
                    default=os.path.join("results", "_stress_aniso"))
    args = ap.parse_args()

    out_dir = (args.out_dir if os.path.isabs(args.out_dir)
               else os.path.join(_HERE, args.out_dir))
    os.makedirs(out_dir, exist_ok=True)
    scratch_root = os.path.join(out_dir, "_scratch")

    rows = []
    for label, k1, k2, k3 in DEFAULT_CASES:
        rec = run_one(args.L, k1, k2, k3,
                      beta_lo_factor=args.beta_lo_factor,
                      beta_hi_factor=args.beta_hi_factor,
                      n_traj_coarse=args.n_traj_coarse,
                      n_traj_fine=args.n_traj_fine,
                      n_coarse=args.n_coarse, n_refine=args.n_refine,
                      jackknife=not args.no_jackknife,
                      scratch_root=scratch_root)
        rec["label"] = label
        rows.append(rec)

    # Δ_iso(L) calibration: difference between MC and exact for the
    # isotropic row in this run.
    iso = next(r for r in rows if r["label"] == "iso")
    delta_iso = iso["beta_c"] - iso["beta_c_exact"]

    print()
    print(f"Reference shift Δ_iso(L={args.L}) "
          f"= β_c(L)_iso − β_c(∞)_iso = {delta_iso:+.5f}")
    print()
    hdr = (f"{'label':<12}{'k1':>5}{'k2':>5}{'k3':>5}"
           f"{'β_c_exact':>11}{'β_c(L)':>11}{'σ':>10}"
           f"{'Δ(L)':>9}{'Δ−Δ_iso':>11}{'(Δ−Δ_iso)/σ':>14}")
    print(hdr)
    print("-" * len(hdr))

    summary_rows = []
    for r in rows:
        delta = r["beta_c"] - r["beta_c_exact"]
        delta_corr = delta - delta_iso
        z = delta_corr / max(r["beta_c_sigma"], 1e-9)
        print(f"{r['label']:<12}{r['k1']:>5.2f}{r['k2']:>5.2f}{r['k3']:>5.2f}"
              f"{r['beta_c_exact']:>11.6f}{r['beta_c']:>11.6f}"
              f"{r['beta_c_sigma']:>10.2e}{delta:>+9.5f}{delta_corr:>+11.5f}"
              f"{z:>+14.2f}")
        r["delta"] = float(delta)
        r["delta_iso_corrected"] = float(delta_corr)
        r["z"] = float(z)
        summary_rows.append(r)

    summary = {"L": args.L, "delta_iso": delta_iso, "rows": summary_rows}
    with open(os.path.join(out_dir, "stress_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # ----- plot ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

        labels = [r["label"] for r in rows]
        bc_meas = np.array([r["beta_c"] for r in rows])
        bc_sig  = np.array([r["beta_c_sigma"] for r in rows])
        bc_exact= np.array([r["beta_c_exact"] for r in rows])
        x = np.arange(len(rows))

        ax = axes[0]
        ax.errorbar(x, bc_meas, yerr=bc_sig, fmt="o", color="C0",
                    capsize=3, label=f"MC β_c(L={args.L})")
        ax.plot(x, bc_exact, "x", color="crimson", ms=10, mew=2,
                label="exact β_c(∞)")
        ax.plot(x, bc_exact + delta_iso, "s", color="C2", mfc="none", ms=10,
                mew=1.5,
                label=f"exact + Δ_iso(L={args.L})")
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right",
                                              fontsize=8)
        ax.set_ylabel("β_c")
        ax.set_title("β_c finder vs exact across anisotropies")
        ax.legend(loc="best", fontsize=8); ax.grid(alpha=0.3)

        ax = axes[1]
        zs = np.array([r["z"] for r in rows])
        ax.bar(x, zs, color=["C0" if abs(z) < 2 else "crimson" for z in zs],
               edgecolor="black")
        ax.axhline(0, color="black", lw=1.0)
        ax.axhline(2, color="grey", ls="--", lw=0.7)
        ax.axhline(-2, color="grey", ls="--", lw=0.7)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right",
                                              fontsize=8)
        ax.set_ylabel("(β_c(L) − β_c_exact − Δ_iso) / σ")
        ax.set_title("FSS-corrected residual in σ units")
        ax.grid(alpha=0.3, axis="y")

        fig.suptitle(f"β_c finder stress test — anisotropic triangular Ising, "
                     f"L={args.L}", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out_png = os.path.join(out_dir, "stress_aniso.png")
        fig.savefig(out_png, dpi=130)
        print(f"[plot] {out_png}")
    except Exception as e:
        print(f"[plot] skipped ({e})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
