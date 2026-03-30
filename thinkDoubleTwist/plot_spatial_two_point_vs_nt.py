#!/usr/bin/env python3
"""Compare equal-time spatial connected two-point functions across N_t runs.

Input files are the `two_point_all_to_all.dat` outputs produced by
`ising_tri_twisted_parallelogram` with optional time extrusion.
Each row is expected to include:
  d m n corr err corr_conn err_conn

This script computes triangular-lattice spatial radius
  r = sqrt(m^2 + n^2 + m*n),
then bins connected correlators by rounded radius and overlays curves.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_connected(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows: list[tuple[int, int, float, float]] = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, _c_s, _e_s, cc_s, ec_s = s.split()
            rows.append((int(m_s), int(n_s), float(cc_s), float(ec_s)))

    if not rows:
        raise SystemExit(f"No data rows found in {path}")

    arr = np.array(rows, dtype=float)
    m = arr[:, 0]
    n = arr[:, 1]
    conn = arr[:, 2]
    err = arr[:, 3]
    return m, n, conn, err


def radial_bin_stats(m: np.ndarray, n: np.ndarray, conn: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = np.sqrt(m * m + n * n + m * n)
    r_int = np.rint(r).astype(int)

    uniq = np.unique(r_int)
    means: list[float] = []
    for rv in uniq:
        mask = r_int == rv
        means.append(float(np.mean(conn[mask])))

    return uniq.astype(float), np.array(means, dtype=float)


def label_from_path(path: Path) -> str:
    # Expected parent name contains `NtX`, e.g. .../Lx39_..._Nt4_.../
    name = path.parent.name
    for token in name.split("_"):
        if token.startswith("Nt"):
            return token
    return name


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Overlay spatial connected two-point evolution vs N_t")
    ap.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="List of two_point_all_to_all.dat files to compare",
    )
    ap.add_argument("--output", required=True, help="Output figure path")
    ap.add_argument("--r-max", type=float, default=None, help="Optional max radius to display")
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    fig, ax = plt.subplots(figsize=(8.2, 5.4), dpi=200, constrained_layout=True)

    for p in args.inputs:
        path = Path(p)
        m, n, conn, _err = load_connected(path)
        r, cmean = radial_bin_stats(m, n, conn)

        if args.r_max is not None:
            keep = r <= float(args.r_max)
            r = r[keep]
            cmean = cmean[keep]

        ax.plot(r, cmean, marker="o", ms=2.6, lw=1.2, label=label_from_path(path))

    ax.set_xlabel("spatial radius r (triangular metric)")
    ax.set_ylabel("connected two-point (radial mean)")
    ax.set_title("Equal-time Spatial Connected Two-Point vs Time Extent")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
