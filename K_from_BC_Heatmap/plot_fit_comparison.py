#!/usr/bin/env python3
"""
plot_fit_comparison.py — compare linear, quadratic, and exponential continuum
extrapolations side-by-side for a given (tag, cost-mode).

Requires a .npz produced by compute_landscape.py with ``--fit all``:

    python compute_landscape.py --tag smoke --cost-mode huber_log --fit all

Then:

    python plot_fit_comparison.py --tag smoke --cost-mode huber_log

Produces (under results/<tag>/plots/):
  fit_comparison_<mode>.png
    A grid of heatmaps:
      - Row for each fit method (linear / quadratic / exponential):
        cost_inf, sigma_inf
      - Final summary row:
        std(cost_inf) across the three fits, and
        range(cost_inf) = max − min across the three fits.
  fit_comparison_<mode>_diff.png
    Pairwise difference maps: quadratic−linear, exponential−linear,
    exponential−quadratic.
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

_FIT_LABELS = {
    "linear":      r"linear  ($a + b/L$)",
    "quadratic":   r"quadratic  ($a + b/L + c/L^2$)",
    "exponential": r"exponential  ($a + b\,e^{-L/\xi}$)",
}


def _heatmap(ax, rs, Z, title, vmin=None, vmax=None, cmap="viridis", sym=False):
    """Render one heatmap panel and return the AxesImage."""
    finite = Z[np.isfinite(Z)]
    if finite.size == 0:
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        return None
    if sym:
        lim = float(np.nanmax(np.abs(finite)))
        vmin, vmax = -lim, lim
        cmap = "RdBu_r"
    else:
        if vmin is None:
            vmin = float(np.nanmin(finite))
        if vmax is None:
            vmax = float(np.nanpercentile(finite, 99))
    extent = [rs[0], rs[-1], rs[0], rs[-1]]
    im = ax.imshow(Z.T, origin="lower", extent=extent, aspect="auto",
                   cmap=cmap, vmin=vmin, vmax=vmax)
    # argmin marker
    idx = np.unravel_index(np.nanargmin(Z), Z.shape)
    ax.plot(rs[idx[0]], rs[idx[1]], "w*", ms=9,
            markeredgecolor="k",
            label=f"min=({rs[idx[0]]:.2f},{rs[idx[1]]:.2f})")
    ax.legend(loc="upper right", fontsize=6)
    ax.set_xlabel("r1", fontsize=8)
    ax.set_ylabel("r2", fontsize=8)
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="smoke")
    ap.add_argument("--cost-mode", default="huber_log")
    ap.add_argument("--in-name", default=None,
                    help="override npz basename (default: cost_<mode>_all_fits)")
    args = ap.parse_args()

    root    = os.path.join(_HERE, "results", args.tag)
    in_name = args.in_name or f"cost_{args.cost_mode}_all_fits"
    npz_path = os.path.join(root, f"{in_name}.npz")
    if not os.path.exists(npz_path):
        print(f"ERROR: {npz_path} not found.\n"
              f"Run: python compute_landscape.py "
              f"--tag {args.tag} --cost-mode {args.cost_mode} --fit all",
              file=sys.stderr)
        return 1

    z = np.load(npz_path)
    rs         = z["rs"]
    sizes      = list(z["sizes"])
    twist_frac = float(z["twist_frac"])

    # Detect which fit modes are present
    fits_present = [m for m in ["linear", "quadratic", "exponential"]
                    if f"cost_inf_{m}" in z]
    if not fits_present:
        print("ERROR: npz has no fit-suffixed keys — was it built with --fit all?",
              file=sys.stderr)
        return 1

    plot_dir = os.path.join(root, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    cost_inf  = {m: z[f"cost_inf_{m}"]  for m in fits_present}
    sigma_inf = {m: z[f"sigma_inf_{m}"] for m in fits_present}

    # -----------------------------------------------------------------------
    # Figure 1 — cost_inf and sigma_inf for each fit, plus summary row
    # -----------------------------------------------------------------------
    n_fit_rows  = len(fits_present)
    n_summary   = 1  # spread row
    n_rows      = n_fit_rows + n_summary
    n_cols      = 2  # cost_inf | sigma_inf
    # extend last row to 2 columns (std, range)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 4.6 * n_rows),
                             squeeze=False)

    for row, m in enumerate(fits_present):
        label = _FIT_LABELS.get(m, m)
        _heatmap(axes[row, 0], rs, cost_inf[m],
                 f"{label}\ncost_inf")
        _heatmap(axes[row, 1], rs, sigma_inf[m],
                 f"{label}\nσ_inf", cmap="magma")

    # Summary row: std and range across fit methods
    stack = np.stack([cost_inf[m] for m in fits_present], axis=0)
    std_map   = np.nanstd(stack,  axis=0)
    range_map = np.nanmax(stack,  axis=0) - np.nanmin(stack, axis=0)

    _heatmap(axes[n_fit_rows, 0], rs, std_map,
             "std(cost_inf) across fit types", cmap="magma")
    _heatmap(axes[n_fit_rows, 1], rs, range_map,
             "range(cost_inf) = max−min across fit types", cmap="hot")

    fig.suptitle(
        f"Continuum extrapolation comparison\n"
        f"tag={args.tag}  cost={args.cost_mode}  "
        f"sizes={sizes}  twist_frac={twist_frac}",
        fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out1 = os.path.join(plot_dir, f"fit_comparison_{args.cost_mode}.png")
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"[plot] {out1}")

    # -----------------------------------------------------------------------
    # Figure 2 — pairwise difference maps
    # -----------------------------------------------------------------------
    pairs = []
    for ia, ma in enumerate(fits_present):
        for mb in fits_present[ia + 1:]:
            pairs.append((ma, mb))

    if pairs:
        n_pair_cols = min(3, len(pairs))
        n_pair_rows = (len(pairs) + n_pair_cols - 1) // n_pair_cols
        fig2, axes2 = plt.subplots(n_pair_rows, n_pair_cols,
                                   figsize=(5.5 * n_pair_cols,
                                            4.6 * n_pair_rows),
                                   squeeze=False)
        axes2_flat = axes2.flatten()
        for k, (ma, mb) in enumerate(pairs):
            diff = cost_inf[mb] - cost_inf[ma]
            la = _FIT_LABELS.get(ma, ma).split("(")[0].strip()
            lb = _FIT_LABELS.get(mb, mb).split("(")[0].strip()
            _heatmap(axes2_flat[k], rs, diff,
                     f"cost_inf:  {lb} − {la}", sym=True)
        for ax in axes2_flat[len(pairs):]:
            ax.axis("off")

        fig2.suptitle(
            f"Pairwise continuum-fit differences\n"
            f"tag={args.tag}  cost={args.cost_mode}  "
            f"sizes={sizes}  twist_frac={twist_frac}",
            fontsize=12)
        fig2.tight_layout(rect=[0, 0, 1, 0.96])
        out2 = os.path.join(plot_dir, f"fit_comparison_{args.cost_mode}_diff.png")
        fig2.savefig(out2, dpi=150)
        plt.close(fig2)
        print(f"[plot] {out2}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
