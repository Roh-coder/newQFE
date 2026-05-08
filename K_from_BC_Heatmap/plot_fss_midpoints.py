#!/usr/bin/env python3
"""Plot only the t=1/2 midpoint of each boundary cycle vs 1/L."""
from __future__ import annotations

import argparse
import os
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from cost import boundary_paths, _tile_interp, _SQRT3_2  # noqa: E402
from precompute_grid import _point_path                  # noqa: E402


def _sample_mid(pkl_path):
    with open(pkl_path, "rb") as f:
        pkl = pickle.load(f)
    data = pkl["data"]
    Lx, Ly, Tx, Ty = pkl["Lx"], pkl["Ly"], pkl["Tx"], pkl["Ty"]
    iG = _tile_interp(data, Lx, Ly, Tx, Ty, "conn",     copies=2)
    iE = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    G  = np.zeros(3)
    sG = np.zeros(3)
    for c, (dm, dn) in enumerate(paths):
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        pt = np.array([[0.5 * ex, 0.5 * ey]])
        G[c]  = float(np.asarray(iG(pt)).ravel()[0])
        sG[c] = abs(float(np.asarray(iE(pt)).ravel()[0]))
    return G, sG, pkl


def _wls_poly(invL, y, sigma, deg):
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    x, yv, sv = invL[mask], y[mask], sigma[mask]
    if len(x) < deg + 1:
        return np.nan, np.nan, None
    X = np.vander(x, deg + 1, increasing=True)
    w = 1.0 / sv ** 2
    cov = np.linalg.inv((X.T * w) @ X)
    beta = cov @ ((X.T * w) @ yv)
    return float(beta[0]), float(np.sqrt(max(cov[0, 0], 0.0))), beta


def _fit_power(L, y, sigma):
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    Lm, ym, sm = L[mask], y[mask], sigma[mask]
    if len(Lm) < 3:
        return np.nan, np.nan, None
    a0 = float(ym[np.argmax(Lm)])
    b0 = float(ym[np.argmin(Lm)] - a0)

    def model(LL, a, b, A):
        return a + b * (1.0 / LL) ** A

    y_min, y_max = float(ym.min()), float(ym.max())
    pad = max(y_max - y_min, 1e-6)
    bounds = ([y_min - 5 * pad, -100 * pad, 0.1],
              [y_max + 5 * pad,  100 * pad, 6.0])
    try:
        popt, pcov = curve_fit(model, Lm, ym, p0=[a0, b0, 1.5],
                                sigma=sm, absolute_sigma=True,
                                bounds=bounds, maxfev=8000)
        return float(popt[0]), float(np.sqrt(max(pcov[0, 0], 0.0))), popt
    except Exception:  # noqa: BLE001
        return np.nan, np.nan, None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="fss_eighths")
    ap.add_argument("--r1",  type=float, default=1.0)
    ap.add_argument("--r2",  type=float, default=1.0)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[4, 6, 8, 12, 16, 24, 32])
    args = ap.parse_args()

    out_root = os.path.join(_HERE, "results", args.tag)
    sizes = sorted(args.sizes)
    Larr = np.asarray(sizes, dtype=float)
    invL = 1.0 / Larr

    G  = {fam: np.full((3, len(sizes)), np.nan) for fam in ("test", "ref")}
    sG = {fam: np.full((3, len(sizes)), np.nan) for fam in ("test", "ref")}
    geoms = {}
    for li, L in enumerate(sizes):
        for fam in ("test", "ref"):
            pkl_path = _point_path(
                os.path.join(out_root, "grid", f"L{L}", fam),
                args.r1, args.r2)
            if not os.path.exists(pkl_path):
                print(f"[skip] missing {pkl_path}")
                continue
            g, sg, pkl = _sample_mid(pkl_path)
            G[fam][:, li]  = g
            sG[fam][:, li] = sg
            geoms[(L, fam)] = (pkl["Lx"], pkl["Ly"], pkl["Tx"], pkl["Ty"])

    invL_dense = np.linspace(0.0, invL.max() * 1.05, 200)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    cyc_color = ["#1f77b4", "#d62728", "#2ca02c"]
    fam_title = {"test": "test (untwisted)", "ref": "ref (twisted T = round(L/4))"}

    for fi, fam in enumerate(("test", "ref")):
        ax = axes[fi]
        for c in range(3):
            y, sy = G[fam][c], sG[fam][c]
            ax.errorbar(invL, y, yerr=sy, fmt="o", ms=6, color=cyc_color[c],
                        ecolor=cyc_color[c], capsize=3, mec="k", mew=0.5,
                        label=f"cycle {c}", zorder=4)

        ax.axvline(0.0, color="k", lw=0.6, ls=":", alpha=0.5)
        ax.set_xlabel(r"$1/L$")
        ax.set_ylabel(r"$G_{\rm conn}(t = 1/2)$")
        ax.set_title(fam_title[fam])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc="upper left")

    fig.suptitle(
        f"Boundary midpoint G(t=1/2) vs 1/L  @ (r1, r2) = ({args.r1:g}, {args.r2:g})  "
        f"sizes {sizes}",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plot_dir = os.path.join(out_root, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    out_png = os.path.join(plot_dir,
        f"fss_midpoint_r1_{args.r1:.3f}_r2_{args.r2:.3f}.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"[plot] {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
