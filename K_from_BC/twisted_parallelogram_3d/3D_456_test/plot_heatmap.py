#!/usr/bin/env python3
"""Heatmaps of the all-to-all connected spatial correlator on the 4-5-6 torus.

Each site in the unit cell is mapped to its physical (x,y) position via the
triangular lattice basis:
  x = m + 0.5 * n
  y = (sqrt(3)/2) * n

The connected correlator value (corr_conn) is shown as a scatter plot coloured
by log10(corr_conn), giving a 2D heatmap of the correlator landscape.

One subplot per scale factor s found in --datadir.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def load_all_to_all(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            m, n = int(parts[1]), int(parts[2])
            corr_conn = float(parts[5])
            err_conn  = float(parts[6])
            rows.append((m, n, corr_conn, err_conn))
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def to_xy(m, n):
    x = m + 0.5 * n
    y = (np.sqrt(3.0) / 2.0) * n
    return x, y


def parse_geometry(dirname: str):
    mt = re.search(r"Lx(\d+)_Ly(\d+)_Tx(-?\d+)_Ty(-?\d+)_Nt(\d+)", dirname)
    if not mt:
        return None
    return {k: int(v) for k, v in zip(["Lx", "Ly", "Tx", "Ty", "Nt"], mt.groups())}


def draw_period_vectors(ax, geo, color="white", alpha=0.7):
    """Overlay the two spatial period vectors of the torus."""
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.5, np.sqrt(3.0) / 2.0])
    lx, ly, tx, ty = geo["Lx"], geo["Ly"], geo["Tx"], geo["Ty"]
    v1 = lx * e1 + ty * e2   # period along first boundary direction
    v2 = tx * e1 - ly * e2   # period along second boundary direction
    origin = np.array([0.0, 0.0])
    for v, label in [(v1, "v₁"), (v2, "v₂")]:
        ax.annotate("", xy=v, xytext=origin,
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.5, alpha=alpha))
        ax.text(v[0] * 0.55, v[1] * 0.55, label, color=color,
                fontsize=8, alpha=alpha, ha="center")


def main():
    ap = argparse.ArgumentParser(description="Heatmaps of all-to-all connected correlator")
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="heatmap_all_to_all.png", help="Output figure path")
    ap.add_argument("--vmin", type=float, default=None,
                    help="log10 colour scale minimum (default: auto)")
    ap.add_argument("--vmax", type=float, default=None,
                    help="log10 colour scale maximum (default: auto)")
    ap.add_argument("--marker-size", type=float, default=None,
                    help="Scatter marker size (default: auto from density)")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_all_to_all.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No two_point_all_to_all.dat files found: {pattern}")

    n_panels = len(files)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 6))
    if n_panels == 1:
        axes = [axes]

    # Shared colour limits across all panels for direct comparison
    all_log_vals = []
    datasets = []
    for fpath in files:
        fpath = Path(fpath)
        geo = parse_geometry(fpath.parent.name)
        m, n, conn, _ = load_all_to_all(fpath)
        x, y = to_xy(m, n)
        pos_mask = conn > 0
        log_conn = np.where(pos_mask, np.log10(np.where(pos_mask, conn, 1e-30)), np.nan)
        datasets.append((x, y, log_conn, geo, fpath))
        all_log_vals.append(log_conn[np.isfinite(log_conn)])

    all_vals = np.concatenate(all_log_vals)
    vmin = args.vmin if args.vmin is not None else float(np.percentile(all_vals, 2))
    vmax = args.vmax if args.vmax is not None else float(np.percentile(all_vals, 99))

    cmap = plt.cm.inferno

    for ax, (x, y, log_conn, geo, fpath) in zip(axes, datasets):
        if geo is None:
            continue
        lx, nt = geo["Lx"], geo["Nt"]
        s = lx // 39

        # Marker size: shrink for larger lattices so dots don't overlap
        n_pts = len(x)
        ms = args.marker_size if args.marker_size is not None else max(0.3, 12.0 / (s ** 1.2))

        sc = ax.scatter(x, y, c=log_conn, cmap=cmap, vmin=vmin, vmax=vmax,
                        s=ms, linewidths=0, rasterized=True)

        # Mark the origin (self-correlator)
        ax.scatter([0], [0], c="cyan", s=20, zorder=5, label="origin")

        # Draw period vectors
        draw_period_vectors(ax, geo)

        ax.set_aspect("equal")
        ax.set_title(f"s={s}  (Lx={lx}, Nt={nt})\n{n_pts} sites", fontsize=10)
        ax.set_xlabel("x (lattice units)")
        if ax is axes[0]:
            ax.set_ylabel("y (lattice units)")
        ax.tick_params(labelsize=8)

    # Shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.05)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(r"$\log_{10}\,G_{\rm conn}(r)$", fontsize=10)

    fig.suptitle("Equal-time connected spatial correlator: 4-5-6 torus, K=0.161",
                 fontsize=12, y=1.01)

    outpath = (datadir / args.output
               if not Path(args.output).is_absolute() else Path(args.output))
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
