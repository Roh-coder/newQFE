#!/usr/bin/env python3
"""
plot_landscape.py — render heatmaps from compute_landscape.py output.

Produces:
  - one heatmap per finite L
  - one continuum (L→∞) heatmap
  - one "extrapolation uncertainty" heatmap (σ_inf)
  - one combined panel PNG

Examples
--------
python plot_landscape.py --tag default --cost-mode huber_log
python plot_landscape.py --tag default --cost-mode huber_log --log
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))


def _heatmap(ax, rs, Z, title, log=False, vmax=None, cmap="viridis"):
    Zp = Z.copy()
    if log:
        Zp = np.log10(np.where(Zp > 0, Zp, np.nan))
    finite = Zp[np.isfinite(Zp)]
    if vmax is None and finite.size > 0:
        vmax = float(np.nanpercentile(finite, 99))
    vmin = float(np.nanmin(finite)) if finite.size > 0 else 0.0
    extent = [rs[0], rs[-1], rs[0], rs[-1]]
    im = ax.imshow(Zp.T, origin="lower", extent=extent, aspect="auto",
                   cmap=cmap, vmin=vmin, vmax=vmax)
    if finite.size > 0:
        # Mark argmin.
        idx = np.unravel_index(np.nanargmin(Z), Z.shape)
        ax.plot(rs[idx[0]], rs[idx[1]], "w*", ms=10,
                markeredgecolor="k", label=f"min=({rs[idx[0]]:.2f},{rs[idx[1]]:.2f})")
        ax.legend(loc="upper right", fontsize=7)
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--cost-mode", default="huber_log")
    ap.add_argument("--log", action="store_true",
                    help="render log10(cost)")
    ap.add_argument("--vmax", type=float, default=None,
                    help="cap colour at this cost value (default: 99th pct)")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", args.tag)
    npz_path = os.path.join(root, f"cost_{args.cost_mode}.npz")
    if not os.path.exists(npz_path):
        print(f"ERROR: {npz_path} not found — run compute_landscape.py first.",
              file=sys.stderr)
        return 1
    z = np.load(npz_path)
    rs = z["rs"]
    sizes = list(z["sizes"])
    twist_frac = float(z["twist_frac"])
    cost_inf = z["cost_inf"]
    sigma_inf = z["sigma_inf"]

    plot_dir = os.path.join(root, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Per-L panels
    n_panels = len(sizes) + 2  # +continuum +sigma
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.6 * n_cols, 4.0 * n_rows),
                             squeeze=False)
    axes_flat = axes.flatten()
    for k, L in enumerate(sizes):
        Z = z[f"cost_L{L}"]
        _heatmap(axes_flat[k], rs, Z,
                 f"L={L}  cost={args.cost_mode}",
                 log=args.log, vmax=args.vmax)
    _heatmap(axes_flat[len(sizes)], rs, cost_inf,
             f"continuum (L→∞)  fit={args.cost_mode}",
             log=args.log, vmax=args.vmax)
    _heatmap(axes_flat[len(sizes) + 1], rs, sigma_inf,
             "σ(continuum extrapolation)",
             log=False, cmap="magma")
    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    fig.suptitle(
        f"Untwisted vs twisted cost landscape  "
        f"(tag={args.tag}, twist_frac={twist_frac}, sizes={sizes})",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, f"landscape_{args.cost_mode}"
                                 f"{'_log' if args.log else ''}.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[plot] {out}")

    # Standalone continuum heatmap (publication-ready)
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    _heatmap(ax, rs, cost_inf,
             f"Continuum cost (twist_frac={twist_frac})",
             log=args.log, vmax=args.vmax)
    fig.tight_layout()
    out2 = os.path.join(plot_dir, f"continuum_{args.cost_mode}"
                                  f"{'_log' if args.log else ''}.png")
    fig.savefig(out2, dpi=160)
    plt.close(fig)
    print(f"[plot] {out2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
