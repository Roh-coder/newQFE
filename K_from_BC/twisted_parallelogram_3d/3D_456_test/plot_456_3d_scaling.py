#!/usr/bin/env python3
"""Phase 2a: Scaling collapse of equal-time spatial connected two-point function.

Reads two_point_all_to_all.dat from each scale factor s=1,2,3 directory and
overlays G_spatial(r) * r^(2*Delta_sigma) vs r/L for a scaling collapse plot.

3D Ising spin scaling dimension: Delta_sigma ≈ 0.5182
Triangular-lattice distance: r = sqrt(m^2 + n^2 + m*n)
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DELTA_SIGMA = 0.5182  # 3D Ising spin scaling dimension


def load_connected(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            m, n, cc, ec = int(parts[1]), int(parts[2]), float(parts[5]), float(parts[6])
            rows.append((m, n, cc, ec))
    if not rows:
        raise SystemExit(f"No data in {path}")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def parse_geometry(dirname: str):
    """Extract Lx, Nt from directory name like Lx39_Ly48_Tx9_Ty-9_Nt4_k0.161"""
    m = re.search(r"Lx(\d+).*Nt(\d+)", dirname)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def radial_bin(m, n, vals, errs, n_bins=40):
    r = np.sqrt(m * m + n * n + m * n)
    rmax = r.max()
    edges = np.linspace(0, rmax, n_bins + 1)
    centers, means, errors = [], [], []
    for i in range(n_bins):
        mask = (r >= edges[i]) & (r < edges[i + 1]) & (vals > 0)
        if mask.sum() == 0:
            continue
        w = 1.0 / (errs[mask] ** 2 + 1e-30)
        mean = np.sum(w * vals[mask]) / np.sum(w)
        err = 1.0 / np.sqrt(np.sum(w))
        centers.append(0.5 * (edges[i] + edges[i + 1]))
        means.append(mean)
        errors.append(err)
    return np.array(centers), np.array(means), np.array(errors)


def main():
    ap = argparse.ArgumentParser(description="Scaling collapse of 3D Ising spatial two-point function")
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="plot_456_3d_scaling.png", help="Output figure path")
    ap.add_argument("--delta", type=float, default=DELTA_SIGMA, help="Spin scaling dimension")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_all_to_all.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No data files found matching: {pattern}")

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(files)))

    for fpath, color in zip(files, colors):
        fpath = Path(fpath)
        dirname = fpath.parent.name
        Lx, Nt = parse_geometry(dirname)
        if Lx is None:
            continue
        # infer scale factor from Lx (base Lx=39 at s=1)
        s = Lx // 39

        m, n, conn, err = load_connected(fpath)
        # physical triangular-lattice distance
        r = np.sqrt(m * m + n * n + m * n)
        # linear spatial extent: L = Lx (number of unit cells along e1)
        L = float(Lx)

        # bin in r/L
        r_over_L = r / L
        rmax = r_over_L.max()
        n_bins = 40
        edges = np.linspace(0, rmax, n_bins + 1)
        x_centers, y_means, y_errs = [], [], []
        for i in range(n_bins):
            mask = (r_over_L >= edges[i]) & (r_over_L < edges[i + 1]) & (conn > 0)
            if mask.sum() == 0:
                continue
            w = 1.0 / (err[mask] ** 2 + 1e-30)
            rmid = 0.5 * (edges[i] + edges[i + 1])
            r_phys = rmid * L
            mean_g = np.sum(w * conn[mask]) / np.sum(w)
            err_g = 1.0 / np.sqrt(np.sum(w))
            # scaled: G(r) * r^(2*Delta) using the bin-center physical radius
            scaled = mean_g * r_phys ** (2 * args.delta)
            scaled_err = err_g * r_phys ** (2 * args.delta)
            x_centers.append(rmid)
            y_means.append(scaled)
            y_errs.append(scaled_err)

        x = np.array(x_centers)
        y = np.array(y_means)
        ye = np.array(y_errs)

        label = f"s={s}  (Lx={Lx}, Nt={Nt})"
        ax.errorbar(x, y, yerr=ye, fmt="o-", markersize=3, color=color, label=label, lw=1)

    ax.set_xlabel(r"$r / L_x$")
    ax.set_ylabel(r"$G_{\rm spatial}(r) \cdot r^{2\Delta_\sigma}$")
    ax.set_title(
        rf"Scaling collapse: 4-5-6 torus, 3D Ising, $\Delta_\sigma={args.delta}$, K=0.161"
    )
    ax.legend()
    ax.set_yscale("log")
    fig.tight_layout()
    outpath = datadir / args.output if not Path(args.output).is_absolute() else Path(args.output)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
