#!/usr/bin/env python3
"""Per-anisotropy finite-size-scaling stress test for the β_c finder.

For each (k1, k2, k3) case, scan β_c(L) at several L and fit

    β_c(L) = β_c(∞) + a · L^{-1/ν}

with ν fixed to 1 (Ising universality).  Compare β_c(∞) to the exact
bulk value from sinh(2β k_i)sinh(2β k_j) sum = 1.

Caches are stored under
    results/_fss_aniso/<case>/L<L>/scan.json
so re-runs only do what is missing.
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


# (label, k1, k2, k3)
DEFAULT_CASES = [
    ("iso",        1.0, 1.0, 1.0),
    ("aniso_1.2",  1.2, 1.0, 0.8),
    ("aniso_1.5",  1.5, 1.0, 0.5),
    ("two_eq_0.7", 1.0, 1.0, 0.7),
    ("two_eq_1.5", 1.0, 1.0, 1.5),
]

DEFAULT_SIZES = [10, 12, 16, 20, 24]


def case_dir(out_root: str, label: str) -> str:
    return os.path.join(out_root, label)


def scan_one(label: str, L: int, k1: float, k2: float, k3: float,
             *, beta_lo_factor: float, beta_hi_factor: float,
             n_traj_coarse: int, n_traj_fine: int,
             n_coarse: int, n_refine: int,
             jackknife: bool, out_root: str, force: bool) -> dict:
    cdir = case_dir(out_root, label)
    os.makedirs(cdir, exist_ok=True)
    cache = os.path.join(cdir, f"L{L}.json")
    bc_exact = exact_beta_c(k1, k2, k3)
    if os.path.exists(cache) and not force:
        with open(cache) as f:
            d = json.load(f)
        d["from_cache"] = True
        return d

    beta_lo = max(0.05, beta_lo_factor * bc_exact)
    beta_hi = beta_hi_factor * bc_exact
    print(f"[{label} L={L}] β_c_exact={bc_exact:.6f}; "
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
        data_dir=os.path.join(cdir, f"_scratch_L{L}"),
    )
    dt = time.time() - t0
    print(f"[{label} L={L}]  β_c(L)={bc:.6f} ± {sigma:.2e}  "
          f"χ_peak={chi_peak:.3g}  ({dt:.1f}s)", flush=True)
    rec = {
        "label": label, "k1": k1, "k2": k2, "k3": k3,
        "L": L,
        "beta_c_exact": float(bc_exact),
        "beta_c": float(bc), "beta_c_sigma": float(sigma),
        "chi_peak": float(chi_peak),
        "wall_s": dt, "from_cache": False,
    }
    with open(cache, "w") as f:
        json.dump(rec, f, indent=2)
    return rec


def fit_nu1(Ls: np.ndarray, bcs: np.ndarray, sigmas: np.ndarray) -> dict:
    Ls = np.asarray(Ls, dtype=float)
    bcs = np.asarray(bcs, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    w = 1.0 / np.maximum(sigmas, 1e-6) ** 2
    x = 1.0 / Ls
    S   = w.sum()
    Sx  = (w * x).sum();  Sy  = (w * bcs).sum()
    Sxx = (w * x * x).sum(); Sxy = (w * x * bcs).sum()
    det = S * Sxx - Sx * Sx
    a       = (S * Sxy - Sx * Sy) / det
    beta_inf= (Sxx * Sy - Sx * Sxy) / det
    var_a   = S / det
    var_b   = Sxx / det
    resid   = bcs - (beta_inf + a * x)
    chi2    = float((w * resid * resid).sum())
    dof     = max(len(Ls) - 2, 1)
    return {"beta_c_inf": float(beta_inf),
            "beta_c_inf_err": float(math.sqrt(max(var_b, 0.0))),
            "a": float(a), "a_err": float(math.sqrt(max(var_a, 0.0))),
            "chi2": chi2, "dof": dof,
            "chi2_per_dof": chi2 / dof}


def fit_nu1_with_subleading(Ls, bcs, sigmas) -> dict:
    """β_c(L) = β_c(∞) + a/L + b/L^2 (3-parameter weighted lin lstsq)."""
    Ls = np.asarray(Ls, dtype=float)
    bcs = np.asarray(bcs, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    if len(Ls) < 3:
        return {}
    w = 1.0 / np.maximum(sigmas, 1e-6)
    A = np.stack([np.ones_like(Ls), 1.0 / Ls, 1.0 / Ls ** 2], axis=1)
    Aw = A * w[:, None]; yw = bcs * w
    sol, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    perr = np.sqrt(np.diag(cov))
    yhat = A @ sol
    chi2 = float(np.sum(((bcs - yhat) / np.maximum(sigmas, 1e-6)) ** 2))
    dof = max(len(Ls) - 3, 1)
    return {"beta_c_inf": float(sol[0]),
            "beta_c_inf_err": float(perr[0]),
            "a": float(sol[1]), "a_err": float(perr[1]),
            "b": float(sol[2]), "b_err": float(perr[2]),
            "chi2": chi2, "dof": dof, "chi2_per_dof": chi2 / dof}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=DEFAULT_SIZES)
    ap.add_argument("--cases", type=str, nargs="+", default=None,
                    help="Subset of case labels to run.")
    ap.add_argument("--n-traj-coarse", type=int, default=4000)
    ap.add_argument("--n-traj-fine", type=int, default=10000)
    ap.add_argument("--n-coarse", type=int, default=11)
    ap.add_argument("--n-refine", type=int, default=5)
    ap.add_argument("--beta-lo-factor", type=float, default=0.55)
    ap.add_argument("--beta-hi-factor", type=float, default=1.45)
    ap.add_argument("--no-jackknife", action="store_true")
    ap.add_argument("--out-dir", type=str,
                    default=os.path.join("results", "_fss_aniso"))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_root = (args.out_dir if os.path.isabs(args.out_dir)
                else os.path.join(_HERE, args.out_dir))
    os.makedirs(out_root, exist_ok=True)

    cases = DEFAULT_CASES
    if args.cases:
        sel = set(args.cases)
        cases = [c for c in cases if c[0] in sel]

    summary = []
    for label, k1, k2, k3 in cases:
        rows = []
        for L in args.sizes:
            r = scan_one(label, L, k1, k2, k3,
                         beta_lo_factor=args.beta_lo_factor,
                         beta_hi_factor=args.beta_hi_factor,
                         n_traj_coarse=args.n_traj_coarse,
                         n_traj_fine=args.n_traj_fine,
                         n_coarse=args.n_coarse, n_refine=args.n_refine,
                         jackknife=not args.no_jackknife,
                         out_root=out_root, force=args.force)
            rows.append(r)
        Ls = np.array([r["L"] for r in rows], dtype=float)
        bcs = np.array([r["beta_c"] for r in rows], dtype=float)
        sigmas = np.array([max(r["beta_c_sigma"], 1e-6) for r in rows], dtype=float)
        bc_exact = rows[0]["beta_c_exact"]

        fit1 = fit_nu1(Ls, bcs, sigmas)
        fit2 = fit_nu1_with_subleading(Ls, bcs, sigmas)
        z1 = (fit1["beta_c_inf"] - bc_exact) / max(fit1["beta_c_inf_err"], 1e-12)
        z2 = ((fit2["beta_c_inf"] - bc_exact) / max(fit2["beta_c_inf_err"], 1e-12)
              if fit2 else None)
        summary.append({
            "label": label, "k1": k1, "k2": k2, "k3": k3,
            "beta_c_exact": bc_exact,
            "rows": rows,
            "fit_nu1_a_only": {**fit1, "z_vs_exact": z1},
            "fit_nu1_a_and_b": ({**fit2, "z_vs_exact": z2} if fit2 else None),
        })

    # ----- print table -----------------------------------------------------
    print()
    hdr = (f"{'case':<11}{'k1':>5}{'k2':>5}{'k3':>5}"
           f"{'β_exact':>10}"
           f"{'β∞(a)':>10}{'σ':>8}{'χ²/dof':>8}{'z':>7}"
           f"{'β∞(a,b)':>11}{'σ':>8}{'χ²/dof':>8}{'z':>7}")
    print(hdr)
    print("-" * len(hdr))
    for s in summary:
        f1 = s["fit_nu1_a_only"]
        f2 = s["fit_nu1_a_and_b"]
        line = (f"{s['label']:<11}{s['k1']:>5.2f}{s['k2']:>5.2f}{s['k3']:>5.2f}"
                f"{s['beta_c_exact']:>10.5f}"
                f"{f1['beta_c_inf']:>10.5f}{f1['beta_c_inf_err']:>8.1e}"
                f"{f1['chi2_per_dof']:>8.2f}{f1['z_vs_exact']:>+7.2f}")
        if f2:
            line += (f"{f2['beta_c_inf']:>11.5f}"
                     f"{f2['beta_c_inf_err']:>8.1e}"
                     f"{f2['chi2_per_dof']:>8.2f}"
                     f"{f2['z_vs_exact']:>+7.2f}")
        print(line)

    with open(os.path.join(out_root, "fss_aniso_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # ----- plot ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ncase = len(summary)
        ncols = 2
        nrows = (ncase + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(5.5 * ncols, 3.4 * nrows),
                                 squeeze=False)
        for i, s in enumerate(summary):
            ax = axes[i // ncols][i % ncols]
            Ls = np.array([r["L"] for r in s["rows"]], dtype=float)
            bcs = np.array([r["beta_c"] for r in s["rows"]], dtype=float)
            sigmas = np.array([max(r["beta_c_sigma"], 1e-6)
                               for r in s["rows"]], dtype=float)
            x = 1.0 / Ls
            xg = np.linspace(0.0, max(x) * 1.05, 200)
            f1 = s["fit_nu1_a_only"]; f2 = s["fit_nu1_a_and_b"]
            ax.errorbar(x, bcs, yerr=sigmas, fmt="o", color="C0",
                        capsize=3, label="MC")
            ax.plot(xg, f1["beta_c_inf"] + f1["a"] * xg, "-",
                    color="crimson",
                    label=f"a/L: β∞={f1['beta_c_inf']:.5f}±"
                          f"{f1['beta_c_inf_err']:.1e}")
            if f2:
                ax.plot(xg, f2["beta_c_inf"] + f2["a"] * xg
                        + f2["b"] * xg ** 2, "--", color="C2",
                        label=f"a/L+b/L²: β∞={f2['beta_c_inf']:.5f}±"
                              f"{f2['beta_c_inf_err']:.1e}")
            ax.axhline(s["beta_c_exact"], color="black", ls=":", lw=1.0,
                       label=f"exact={s['beta_c_exact']:.5f}")
            for L, xi, b in zip(Ls, x, bcs):
                ax.annotate(f"L={int(L)}", (xi, b), xytext=(4, 4),
                            textcoords="offset points", fontsize=7)
            ax.set_xlabel("1/L"); ax.set_ylabel("β_c(L)")
            ax.set_title(f"{s['label']}  (k={s['k1']},{s['k2']},{s['k3']})",
                         fontsize=9)
            ax.legend(loc="best", fontsize=7); ax.grid(alpha=0.3)
        for j in range(ncase, nrows * ncols):
            axes[j // ncols][j % ncols].set_axis_off()
        fig.suptitle("Per-anisotropy FSS of β_c(L) — 2-D triangular Ising",
                     fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        out_png = os.path.join(out_root, "fss_aniso.png")
        fig.savefig(out_png, dpi=130)
        print(f"[plot] {out_png}")
    except Exception as e:
        print(f"[plot] skipped ({e})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
