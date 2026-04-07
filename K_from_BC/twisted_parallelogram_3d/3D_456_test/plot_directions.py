#!/usr/bin/env python3
"""Plot spatial connected correlators along the three boundary directions for each scale s.

Reads two_point_typed.dat files. Directions:
  type 0 / 4: e1        (direction along first lattice vector, length Lx=39s)
  type 1 / 5: e2        (direction along second lattice vector, length Ly=48s)
  type 2 / 6: e2 - e1  (direction along third side of triangle, length ~sqrt(39^2+48^2-39*48)=s*...)

Raw correlator vs r (lattice steps) for types 0,1,2.
Connected correlator vs r for types 4,5,6.
Overlay s=1,2,3 per direction in separate subplots.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DIR_LABELS = {
    0: "e1 (raw)",
    1: "e2 (raw)",
    2: "e2−e1 (raw)",
    4: "e1 (connected)",
    5: "e2 (connected)",
    6: "e2−e1 (connected)",
}

# Triangular-lattice physical step length per unit in each direction
# e1 = (1,0), |e1|=1
# e2 = (1/2, sqrt(3)/2), |e2|=1
# e2-e1 = (-1/2, sqrt(3)/2), |e2-e1|=1
# All have unit length in the triangular lattice; r in each direction = step index * 1.


def load_typed(path: Path):
    """Return dict: type -> (r_array, corr_array, err_array)."""
    data: dict[int, list] = {}
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            t, r, c, e = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
            data.setdefault(t, []).append((r, c, e))
    result = {}
    for t, rows in data.items():
        arr = np.array(rows)
        result[t] = (arr[:, 0], arr[:, 1], arr[:, 2])
    return result


def parse_scale(dirname: str):
    m = re.search(r"Lx(\d+)", dirname)
    return int(m.group(1)) // 39 if m else None


def main():
    ap = argparse.ArgumentParser(description="Plot spatial correlators by boundary direction")
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="plot_directions.png", help="Output figure path")
    ap.add_argument("--connected-only", action="store_true", default=True,
                    help="Plot only connected correlators (types 4,5,6)")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_typed.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No two_point_typed.dat files found: {pattern}")

    types_to_plot = [4, 5, 6]  # connected
    type_labels = {4: "e1", 5: "e2", 6: "e2−e1"}

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(files)))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    for ax, t in zip(axes, types_to_plot):
        for fpath, color in zip(files, colors):
            fpath = Path(fpath)
            s = parse_scale(fpath.parent.name)
            if s is None:
                continue
            m = re.search(r"Lx(\d+).*Nt(\d+)", fpath.parent.name)
            lx, nt = (int(m.group(1)), int(m.group(2))) if m else (39 * s, 4 * s)

            data = load_typed(fpath)
            if t not in data:
                continue
            r, corr, err = data[t]

            # Physical distance: each step in these directions = 1 triangular lattice unit
            # Normalize x-axis by period length (= r.max() for that direction)
            period = r.max()
            x = r / float(period)

            # Only plot the first half (correlator is symmetric about r=period/2)
            half = r <= period // 2 + 1
            ax.errorbar(x[half], corr[half], yerr=err[half],
                        fmt="o-", markersize=3, lw=1, color=color,
                        label=f"s={s}  Lx={lx} Nt={nt}")

        ax.set_yscale("log")
        ax.set_xlabel(r"$r / L_{\rm period}$")
        ax.set_title(f"Direction: {type_labels[t]}")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", ls=":", alpha=0.4)

    axes[0].set_ylabel(r"$G_{\rm conn}(r)$")
    fig.suptitle("Connected spatial correlators along boundary directions\n"
                 "4-5-6 torus, K=0.161", fontsize=12)
    fig.tight_layout()
    outpath = datadir / args.output if not Path(args.output).is_absolute() else Path(args.output)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
