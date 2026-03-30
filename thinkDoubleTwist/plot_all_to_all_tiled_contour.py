#!/usr/bin/env python3
"""Contour plot of all-to-all correlator with BC tessellation overlays.

For twisted-parallelogram BC:
  v = (L_x, T_y)
  u = (T_x, -L_y)

This script loads all-to-all data rows (d, m, n, corr, err, corr_conn, err_conn),
converts lattice coordinates to Euclidean triangular-basis coordinates, tiles the
fundamental cell by integer shifts in (v, u), and draws a contour map with
parallelogram boundaries for each tiled copy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def lattice_to_xy(m: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = m + 0.5 * n
    y = (np.sqrt(3.0) / 2.0) * n
    return x, y


def rotate_xy(x: np.ndarray, y: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    c = np.cos(theta)
    s = np.sin(theta)
    xr = c * x - s * y
    yr = s * x + c * y
    return xr, yr


def load_full(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append((int(m_s), int(n_s), float(c_s), float(cc_s)))

    if not rows:
        raise SystemExit("No all-to-all rows found.")

    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def main() -> None:
    ap = argparse.ArgumentParser(description="Tiled contour plot for all-to-all correlator.")
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    ap.add_argument("--output", required=True, help="Output image path")
    ap.add_argument("--levels", type=int, default=24, help="Number of contour levels")
    ap.add_argument(
        "--connected-only",
        action="store_true",
        help="Use connected correlator column (default if set); otherwise raw correlator.",
    )
    ap.add_argument("--tile-a-min", type=int, default=0, help="Minimum tile index along v")
    ap.add_argument("--tile-a-max", type=int, default=1, help="Maximum tile index along v")
    ap.add_argument("--tile-b-min", type=int, default=0, help="Minimum tile index along u")
    ap.add_argument("--tile-b-max", type=int, default=1, help="Maximum tile index along u")
    ap.add_argument(
        "--align-v-x",
        action="store_true",
        help="Rotate coordinates so v=(L_x,T_y) points along +x.",
    )
    args = ap.parse_args()

    m, n, c_raw, c_conn = load_full(Path(args.full))
    z = c_conn if args.connected_only else c_raw

    # Fundamental vectors in lattice coordinates.
    v_m, v_n = args.L_x, args.T_y
    u_m, u_n = args.T_x, -args.L_y

    m_tiles = []
    n_tiles = []
    z_tiles = []
    for ia in range(args.tile_a_min, args.tile_a_max + 1):
        for ib in range(args.tile_b_min, args.tile_b_max + 1):
            m_tiles.append(m + ia * v_m + ib * u_m)
            n_tiles.append(n + ia * v_n + ib * u_n)
            z_tiles.append(z)

    m_all = np.concatenate(m_tiles)
    n_all = np.concatenate(n_tiles)
    z_all = np.concatenate(z_tiles)

    x_all, y_all = lattice_to_xy(m_all, n_all)

    # Optional frame rotation so the v direction lies on +x.
    if args.align_v_x:
        vx, vy = lattice_to_xy(np.array([v_m], dtype=float), np.array([v_n], dtype=float))
        v_angle = np.arctan2(vy[0], vx[0])
        x_all, y_all = rotate_xy(x_all, y_all, -v_angle)

    fig, ax = plt.subplots(figsize=(8.8, 7.4), dpi=220, constrained_layout=True)

    tri = mtri.Triangulation(x_all, y_all)
    cntr = ax.tricontourf(tri, z_all, levels=args.levels, cmap="viridis")
    cbar = fig.colorbar(cntr, ax=ax, shrink=0.92, pad=0.02)
    cbar.set_label("connected correlator" if args.connected_only else "raw correlator")

    # Draw boundaries for each tile copy so tessellation is explicit.
    base = np.array(
        [
            [0, 0],
            [v_m, v_n],
            [v_m + u_m, v_n + u_n],
            [u_m, u_n],
            [0, 0],
        ],
        dtype=float,
    )

    for ia in range(args.tile_a_min, args.tile_a_max + 1):
        for ib in range(args.tile_b_min, args.tile_b_max + 1):
            shift = np.array([ia * v_m + ib * u_m, ia * v_n + ib * u_n], dtype=float)
            poly = base + shift
            xb, yb = lattice_to_xy(poly[:, 0], poly[:, 1])
            if args.align_v_x:
                xb, yb = rotate_xy(xb, yb, -v_angle)
            ax.plot(xb, yb, color="white", lw=1.1, alpha=0.9)
            ax.plot(xb, yb, color="black", lw=0.35, alpha=0.75)

    geom = f"L_x={args.L_x}, L_y={args.L_y}, T_x={args.T_x}, T_y={args.T_y}"
    tile_txt = f"tiles: a={args.tile_a_min}..{args.tile_a_max}, b={args.tile_b_min}..{args.tile_b_max}"
    extra = "; frame rotated: v || +x" if args.align_v_x else ""
    ax.set_title(f"All-to-All {'Connected' if args.connected_only else 'Raw'} Contour with BC Tessellation\n{geom}; {tile_txt}{extra}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
