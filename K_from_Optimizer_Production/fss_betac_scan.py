#!/usr/bin/env python3
"""Finite-size-scaling test of the β_c finder.

Runs the 3-pass GC susceptibility scan (mc_engine.find_beta_c) on the
isotropic equilateral triangular Ising lattice (k1=k2=k3=1, no twist)
for several Lx=Ly=L sizes, then fits

    β_c(L) = β_c(∞) + a · L^{-1/ν}

For the 2-D Ising universality class on the triangular lattice the
exact answers are β_c(∞) = (1/4) ln 3 ≈ 0.27465307 and ν = 1.

Cached references at results/_reference_LxL_LyL_Tx0_Ty0/ are reused.
New sizes are scanned and a fresh cache file written under
results/_reference_LxL_LyL_Tx0_Ty0/reference_betac_scan.json.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

# ensure the package modules are importable when run from any cwd
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine  # noqa: E402

EXACT_BETAC = math.log(3.0) / 4.0  # ≈ 0.27465307
EXACT_NU = 1.0

DEFAULT_SIZES = [8, 10, 12, 16, 20, 24, 32]


def cache_dir(L: int) -> str:
    return os.path.join(_HERE, "results",
                        f"_reference_Lx{L}_Ly{L}_Tx0_Ty0")


def load_or_run(L: int, *, beta_lo: float, beta_hi: float,
                n_traj_coarse: int, n_traj_fine: int,
                n_coarse: int, n_refine: int,
                jackknife: bool, force: bool) -> dict:
    cdir = cache_dir(L)
    cache = os.path.join(cdir, "reference_betac_scan.json")
    if os.path.exists(cache) and not force:
        with open(cache) as f:
            d = json.load(f)
        d["L"] = L
        d["from_cache"] = True
        return d

    os.makedirs(cdir, exist_ok=True)
    data_dir = os.path.join(cdir, "_fss_scan")
    print(f"[L={L}] scanning β ∈ [{beta_lo}, {beta_hi}] "
          f"(coarse={n_coarse}×{n_traj_coarse}, refine=3×{n_refine}×{n_traj_fine})",
          flush=True)
    t0 = time.time()
    bc, bc_sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
        os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram"),
        L, L, 0, 0, 1.0, 1.0, 1.0,
        beta_lo, beta_hi,
        n_coarse=n_coarse,
        n_refine=n_refine, n_refine2=n_refine, n_refine3=n_refine,
        n_traj_coarse=n_traj_coarse, n_traj_fine=n_traj_fine,
        max_shifts=4, jackknife=jackknife,
        data_dir=data_dir,
    )
    dt = time.time() - t0
    print(f"[L={L}]  β_c={bc:.6f} ± {bc_sigma:.2e}  χ_peak={chi_peak:.3g}  "
          f"({dt:.1f}s)", flush=True)

    out = {
        "L": L,
        "beta_c": float(bc),
        "beta_c_sigma": float(bc_sigma),
        "chi_peak": float(chi_peak),
        "scan_betas": [float(b) for b in sb],
        "scan_chis":  [float(c) for c in sc],
        "scan_chi_errs": [float(e) for e in sce],
        "from_cache": False,
        "wall_s": dt,
    }
    with open(cache, "w") as f:
        json.dump(out, f, indent=2)
    return out


def fit_fss(Ls: np.ndarray, betas: np.ndarray, sigmas: np.ndarray,
            *, fix_nu: float | None = None):
    """Fit β_c(L) = β_c_inf + a · L^{-1/ν}.

    If fix_nu is not None, hold ν fixed and do a weighted linear fit.
    Otherwise nonlinear least squares with three free parameters.
    Returns dict with 'beta_c_inf', 'a', 'nu', their stderrs, chisq, dof.
    """
    Ls = np.asarray(Ls, dtype=float)
    betas = np.asarray(betas, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    # guard against zero sigma
    w = 1.0 / np.maximum(sigmas, 1e-6) ** 2

    if fix_nu is not None:
        x = Ls ** (-1.0 / fix_nu)
        # weighted linear: y = β_inf + a x
        S = w.sum()
        Sx = (w * x).sum(); Sy = (w * betas).sum()
        Sxx = (w * x * x).sum(); Sxy = (w * x * betas).sum()
        det = S * Sxx - Sx * Sx
        a = (S * Sxy - Sx * Sy) / det
        beta_inf = (Sxx * Sy - Sx * Sxy) / det
        var_a = S / det
        var_b = Sxx / det
        resid = betas - (beta_inf + a * x)
        chi2 = float((w * resid * resid).sum())
        dof = max(len(Ls) - 2, 1)
        return {"beta_c_inf": float(beta_inf),
                "beta_c_inf_err": float(math.sqrt(max(var_b, 0.0))),
                "a": float(a),
                "a_err": float(math.sqrt(max(var_a, 0.0))),
                "nu": float(fix_nu), "nu_err": 0.0,
                "chi2": chi2, "dof": dof,
                "chi2_per_dof": chi2 / dof}

    # nonlinear: SciPy if available, else simple grid+linear over ν.
    try:
        from scipy.optimize import curve_fit  # type: ignore

        def model(L, b_inf, a, inv_nu):
            return b_inf + a * L ** (-inv_nu)

        # initial guess: ν=1 linear fit
        lin = fit_fss(Ls, betas, sigmas, fix_nu=1.0)
        p0 = [lin["beta_c_inf"], lin["a"], 1.0]
        popt, pcov = curve_fit(model, Ls, betas, p0=p0,
                               sigma=np.maximum(sigmas, 1e-6),
                               absolute_sigma=True, maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
        resid = betas - model(Ls, *popt)
        chi2 = float((w * resid * resid).sum())
        dof = max(len(Ls) - 3, 1)
        return {"beta_c_inf": float(popt[0]),
                "beta_c_inf_err": float(perr[0]),
                "a": float(popt[1]),
                "a_err": float(perr[1]),
                "nu": float(1.0 / popt[2]),
                "nu_err": float(perr[2] / popt[2] ** 2),
                "chi2": chi2, "dof": dof,
                "chi2_per_dof": chi2 / dof}
    except Exception as e:
        print(f"[fit] scipy unavailable or failed ({e}); "
              "falling back to ν grid search", flush=True)
        best = None
        for inv_nu in np.linspace(0.3, 2.5, 221):
            r = fit_fss(Ls, betas, sigmas, fix_nu=1.0 / inv_nu)
            if (best is None) or (r["chi2"] < best["chi2"]):
                best = r
        return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=DEFAULT_SIZES)
    ap.add_argument("--beta-lo", type=float, default=0.20)
    ap.add_argument("--beta-hi", type=float, default=0.32)
    ap.add_argument("--n-traj-coarse", type=int, default=4000)
    ap.add_argument("--n-traj-fine", type=int, default=10000)
    ap.add_argument("--n-coarse", type=int, default=11)
    ap.add_argument("--n-refine", type=int, default=5)
    ap.add_argument("--no-jackknife", action="store_true",
                    help="Disable jackknife σ on β_c (faster).")
    ap.add_argument("--force", action="store_true",
                    help="Re-run even if a cache file exists.")
    ap.add_argument("--out-dir", type=str,
                    default=os.path.join("results", "_fss_betac"))
    args = ap.parse_args()

    out_dir = os.path.join(_HERE, args.out_dir) if not os.path.isabs(args.out_dir) else args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for L in args.sizes:
        rec = load_or_run(L,
                          beta_lo=args.beta_lo, beta_hi=args.beta_hi,
                          n_traj_coarse=args.n_traj_coarse,
                          n_traj_fine=args.n_traj_fine,
                          n_coarse=args.n_coarse, n_refine=args.n_refine,
                          jackknife=not args.no_jackknife,
                          force=args.force)
        rows.append(rec)

    # ----- summary table ---------------------------------------------------
    print()
    print(f"{'L':>4} {'β_c(L)':>11} {'σ(β_c)':>10} "
          f"{'shift to βc_exact':>20} {'cache':>6}")
    print("-" * 60)
    Ls, bcs, sbs = [], [], []
    for r in rows:
        L = int(r["L"])
        bc = float(r["beta_c"])
        sb = float(r.get("beta_c_sigma", 0.0))
        if sb <= 0.0:
            sb = 1e-3  # safety floor for fit weighting
        Ls.append(L); bcs.append(bc); sbs.append(sb)
        cache_tag = "yes" if r.get("from_cache") else "new"
        print(f"{L:>4} {bc:>11.6f} {sb:>10.2e} "
              f"{(bc - EXACT_BETAC):>+20.6f} {cache_tag:>6}")
    Ls = np.array(Ls); bcs = np.array(bcs); sbs = np.array(sbs)

    # ----- fits ------------------------------------------------------------
    fit_free = fit_fss(Ls, bcs, sbs)
    fit_nu1  = fit_fss(Ls, bcs, sbs, fix_nu=1.0)

    print()
    print("FSS fit:  β_c(L) = β_c(∞) + a · L^(-1/ν)")
    print("-" * 60)
    print(f"  free ν :  β_c(∞) = {fit_free['beta_c_inf']:.6f}"
          f" ± {fit_free['beta_c_inf_err']:.2e}"
          f"   ν = {fit_free['nu']:.3f} ± {fit_free['nu_err']:.3f}"
          f"   χ²/dof = {fit_free['chi2_per_dof']:.2f}"
          f"  (dof={fit_free['dof']})")
    print(f"  ν = 1  :  β_c(∞) = {fit_nu1['beta_c_inf']:.6f}"
          f" ± {fit_nu1['beta_c_inf_err']:.2e}"
          f"   a = {fit_nu1['a']:.4f}"
          f"   χ²/dof = {fit_nu1['chi2_per_dof']:.2f}"
          f"  (dof={fit_nu1['dof']})")
    print(f"  exact  :  β_c(∞) = {EXACT_BETAC:.6f}   ν = {EXACT_NU:.3f}")
    print()
    dev_free = (fit_free["beta_c_inf"] - EXACT_BETAC) / max(fit_free["beta_c_inf_err"], 1e-9)
    dev_nu1  = (fit_nu1["beta_c_inf"]  - EXACT_BETAC) / max(fit_nu1["beta_c_inf_err"],  1e-9)
    print(f"  β_c(∞) deviation from exact:  free-ν fit = {dev_free:+.2f}σ"
          f",  ν=1 fit = {dev_nu1:+.2f}σ")

    summary = {
        "exact_beta_c": EXACT_BETAC,
        "exact_nu": EXACT_NU,
        "sizes": Ls.tolist(),
        "beta_c": bcs.tolist(),
        "beta_c_sigma": sbs.tolist(),
        "fit_free_nu": fit_free,
        "fit_fixed_nu1": fit_nu1,
    }
    with open(os.path.join(out_dir, "fss_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # ----- plot ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

        # Panel 1: β_c vs L^{-1/ν}
        ax = axes[0]
        x_ref = Ls ** (-1.0 / EXACT_NU)
        x_grid = np.linspace(0.0, max(x_ref) * 1.05, 200)
        ax.errorbar(x_ref, bcs, yerr=sbs, fmt="o", color="C0",
                    capsize=3, label="MC β_c(L)")
        for L, x, b in zip(Ls, x_ref, bcs):
            ax.annotate(f"L={L}", (x, b), xytext=(4, 4),
                        textcoords="offset points", fontsize=8)
        ax.plot(x_grid,
                fit_nu1["beta_c_inf"] + fit_nu1["a"] * x_grid,
                "-", color="crimson",
                label=f"ν=1 fit: β∞={fit_nu1['beta_c_inf']:.5f}")
        x_grid_free = np.linspace(0.0, max(Ls ** (-1.0 / fit_free['nu'])) * 1.05, 200)
        ax.plot(x_grid_free,
                fit_free["beta_c_inf"] + fit_free["a"] * x_grid_free,
                "--", color="C2",
                label=(f"free ν={fit_free['nu']:.2f} fit:"
                       f" β∞={fit_free['beta_c_inf']:.5f}"))
        ax.axhline(EXACT_BETAC, color="black", ls=":", lw=1.0,
                   label=f"exact ln3/4 = {EXACT_BETAC:.5f}")
        ax.set_xlabel("L^(-1/ν)  (with ν used by each curve)")
        ax.set_ylabel("β_c(L)")
        ax.set_title("Pseudo-critical β_c vs L^(-1/ν)")
        ax.legend(loc="best", fontsize=8); ax.grid(alpha=0.3)

        # Panel 2: log-log of |β_c(L) - β_c(∞)|
        ax = axes[1]
        shifts = np.abs(bcs - fit_free["beta_c_inf"])
        ax.errorbar(Ls, shifts, yerr=sbs, fmt="o", color="C0",
                    capsize=3, label="|β_c(L) - β_c(∞)|")
        L_grid = np.linspace(min(Ls) * 0.9, max(Ls) * 1.1, 200)
        # slope -1/ν line through the largest L
        a_free = abs(fit_free["a"])
        ax.plot(L_grid, a_free * L_grid ** (-1.0 / fit_free["nu"]),
                "--", color="C2",
                label=f"L^(-1/{fit_free['nu']:.2f}) (free)")
        ax.plot(L_grid, abs(fit_nu1["a"]) * L_grid ** (-1.0),
                "-", color="crimson",
                label="L^(-1) (ν=1)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("L"); ax.set_ylabel("|β_c(L) − β_c(∞)|")
        ax.set_title("Shift magnitude (log-log)")
        ax.legend(loc="best", fontsize=8); ax.grid(alpha=0.3, which="both")

        fig.suptitle("β_c finder FSS check — 2-D triangular Ising", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out_png = os.path.join(out_dir, "fss_betac.png")
        fig.savefig(out_png, dpi=130)
        print(f"[plot] {out_png}")
    except Exception as e:
        print(f"[plot] skipped ({e})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
