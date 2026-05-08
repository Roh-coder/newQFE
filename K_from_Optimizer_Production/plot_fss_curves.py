"""
plot_fss_curves.py — per-point J vs 1/L FSS curves for the ladder.

Reads the precomputed ladder pkls and reference data, computes J_α for
each (r1, r2, L_test) cell, then plots J vs 1/L_test for every grid
point in a grid-of-axes figure with the WLS fit line and J_∞ intercept.

Usage:
    python plot_fss_curves.py [--tag _ladder_111] [--out results/...png]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import _tile_interp, _SQRT3_2, boundary_paths


# reuse helpers from landscape_l2_ladder
def _sample_xy_on_torus(Lx, Ly, Tx, Ty, N_samp):
    xs, ys = [], []
    for (dm, dn) in boundary_paths(Lx, Ly, Tx, Ty):
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        for k in range(1, N_samp):
            t = k / N_samp
            xs.append(t * ex); ys.append(t * ey)
    return np.asarray(xs), np.asarray(ys)


def _sample_data(data, Lx, Ly, Tx, Ty, N_samp, copies=2):
    iref = _tile_interp(data, Lx, Ly, Tx, Ty, "conn", copies)
    xs, ys = _sample_xy_on_torus(Lx, Ly, Tx, Ty, N_samp)
    return np.asarray(iref(np.column_stack([xs, ys])), dtype=float)


def _wls_intercept(L, J):
    """OLS fit J = a + b/L, returns (J_inf, J_inf_err)."""
    L = np.asarray(L, float); J = np.asarray(J, float)
    X = np.column_stack([np.ones_like(L), 1.0 / L])
    try:
        beta, res, _, _ = np.linalg.lstsq(X, J, rcond=None)
        cov = np.linalg.inv(X.T @ X) * (
            np.sum((J - X @ beta) ** 2) / max(len(J) - 2, 1))
        return float(beta[0]), float(np.sqrt(cov[0, 0]))
    except np.linalg.LinAlgError:
        return float(np.mean(J)), float(np.std(J))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="_ladder_111")
    ap.add_argument("--test-sizes", type=int, nargs="+", default=[8, 12, 16])
    ap.add_argument("--ref-sizes",  type=int, nargs="+", default=[16, 24, 32])
    ap.add_argument("--N-samp", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rungs = list(zip(args.test_sizes, args.ref_sizes))
    root      = os.path.join(_HERE, "results", args.tag)
    test_root = os.path.join(root, "test")
    ref_root  = os.path.join(root, "ref")

    # Load ref samples
    ref_vec = {}
    for L_ref in args.ref_sizes:
        path = os.path.join(ref_root, f"L{L_ref}", "two_point_all_to_all.dat")
        rd = mc_engine.load_all_to_all(path)
        ref_vec[L_ref] = _sample_data(rd, L_ref, L_ref, 0, 0, args.N_samp)

    # Discover test pkls
    grid_dir = os.path.join(test_root, "grid")
    pts_files: dict[tuple, dict[int, str]] = {}
    for fn in sorted(os.listdir(grid_dir)):
        if not fn.endswith(".pkl"): continue
        base = fn[:-4]
        try:
            parts = base.split("_")
            r1 = float(parts[1]); r2 = float(parts[3]); L = int(parts[4][1:])
        except Exception:
            continue
        pts_files.setdefault((r1, r2), {})[L] = os.path.join(grid_dir, fn)

    pt_list = sorted(pts_files.keys())
    r1s = sorted({p[0] for p in pt_list})
    r2s = sorted({p[1] for p in pt_list})
    nr1, nr2 = len(r1s), len(r2s)

    # Compute J per rung per point
    Js: dict[tuple, list] = {}   # (r1,r2) -> list of J per rung
    Ls_used: list[int] = args.test_sizes

    for pt in pt_list:
        js = []
        for (L_test, L_ref) in rungs:
            if L_test not in pts_files[pt]:
                js.append(np.nan); continue
            with open(pts_files[pt][L_test], "rb") as f:
                pkl = pickle.load(f)
            G_t = _sample_data(pkl["test_data"], L_test, L_test, 0, 0, args.N_samp)
            G_r = ref_vec[L_ref]
            mask = np.isfinite(G_t) & np.isfinite(G_r)
            res  = G_t[mask] - G_r[mask]
            js.append(float(np.sum(res * res)))
        Js[pt] = js

    # --- plot grid of J vs 1/L curves ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    invL = 1.0 / np.array(Ls_used, float)
    invL_fine = np.linspace(0, invL.max() * 1.15, 200)

    # Keep figure small enough to save quickly
    cell_w, cell_h = 1.6, 1.5
    fig, axes = plt.subplots(nr2, nr1,
                             figsize=(cell_w * nr1, cell_h * nr2),
                             sharex=True)
    # axes[row, col]: row = r2 index (top = large r2), col = r1 index
    fig.subplots_adjust(hspace=0.7, wspace=0.5)

    truth = (1.0, 1.0)

    for pt in pt_list:
        r1, r2 = pt
        col = r1s.index(r1)
        row = nr2 - 1 - r2s.index(r2)   # flip so r2 increases upward
        ax  = axes[row, col] if nr2 > 1 else axes[col]

        js = np.asarray(Js[pt], float)
        valid = np.isfinite(js)
        Lv = np.array(Ls_used, float)[valid]
        jv = js[valid]

        # scatter
        ax.scatter(1.0 / Lv, jv, color="steelblue", s=30, zorder=5)

        # fit + extrapolation
        if valid.sum() >= 2:
            ipt, sig = _wls_intercept(Lv, jv)
            # fit line
            b = np.polyfit(1.0 / Lv, jv, 1)  # linear in 1/L
            ax.plot(invL_fine, np.polyval(b, invL_fine),
                    "k--", lw=0.8, alpha=0.7)
            # J_inf marker at 1/L=0
            color = "green" if abs(ipt) < min(jv) else "tomato"
            ax.axhline(ipt, color=color, lw=1.0, ls=":")
            ax.plot(0, ipt, marker="D", ms=5, color=color, zorder=6)
            ipt_str = f"{ipt:.3f}"
        else:
            ipt_str = "n/a"

        ax.axvline(0, color="gray", lw=0.5, ls="--")
        ax.set_title(f"({r1:.2f},{r2:.2f})\nJ∞={ipt_str}",
                     fontsize=6.5, pad=2)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=5)

        # Highlight the truth cell
        if (r1, r2) == truth:
            for spine in ax.spines.values():
                spine.set_edgecolor("red"); spine.set_linewidth(2)

    # shared axis labels
    fig.text(0.5, 0.01, "1 / L_test", ha="center", fontsize=9)
    fig.text(0.01, 0.5, "J  (sum sq. residuals)", va="center",
             rotation="vertical", fontsize=9)

    # column headers = r1 values
    for col, r1 in enumerate(r1s):
        axes[0, col].set_xlabel(f"r₁={r1:.2f}", fontsize=6,
                                labelpad=18)
    # row labels = r2 values
    for row, r2 in enumerate(reversed(r2s)):
        axes[row, 0].set_ylabel(f"r₂={r2:.2f}", fontsize=6)

    fig.suptitle(
        f"FSS curves  J vs 1/L   {args.tag}\n"
        f"test L∈{args.test_sizes}  ref L∈{args.ref_sizes}  "
        f"truth=(1,1) [red border]",
        fontsize=10)

    out_png = args.out or os.path.join(root, "fss_curves.png")
    fig.savefig(out_png, dpi=90)
    print(f"→ {out_png}")


if __name__ == "__main__":
    main()
