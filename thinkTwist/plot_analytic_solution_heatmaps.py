#!/usr/bin/env python3
"""Heatmaps for analytical centered-cosh compensation solutions.

For fixed Nx and ranges of Ny, shift, this script computes analytical
compensation lengths via solve_centered_cosh_lengths and visualizes:
- left_length
- diag_length
- left/diag ratio
- triangle slack = left + diag - x_length

triangle slack > 0  : non-degenerate Euclidean triangle possible
triangle slack = 0  : degenerate boundary
triangle slack < 0  : no Euclidean triangle with those side lengths
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analytic_compensation import solve_centered_cosh_lengths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Nx", type=int, default=32)
    parser.add_argument("--Ny-min", type=int, default=32)
    parser.add_argument("--Ny-max", type=int, default=48)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=32)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument("--out-prefix", default="thinkTwist/analytic_heatmaps_Nx32_Ny32to48_shift0to32")
    args = parser.parse_args()

    ny_vals = np.arange(args.Ny_min, args.Ny_max + 1)
    sh_vals = np.arange(args.shift_min, args.shift_max + 1)

    nY = len(ny_vals)
    nS = len(sh_vals)

    left = np.zeros((nY, nS), dtype=float)
    diag = np.zeros((nY, nS), dtype=float)
    ratio = np.zeros((nY, nS), dtype=float)
    slack = np.zeros((nY, nS), dtype=float)

    px_map = np.zeros((nY, nS), dtype=float)
    pl_map = np.zeros((nY, nS), dtype=float)
    pd_map = np.zeros((nY, nS), dtype=float)

    rows = []

    for iy, ny in enumerate(ny_vals):
        for isf, shift in enumerate(sh_vals):
            sol = solve_centered_cosh_lengths(
                nx=args.Nx,
                ny=int(ny),
                shift=int(shift),
                x_length=args.x_length,
            )

            left[iy, isf] = sol.left_length
            diag[iy, isf] = sol.diag_length
            ratio[iy, isf] = sol.left_length / sol.diag_length if sol.diag_length != 0 else np.nan
            slack[iy, isf] = sol.left_length + sol.diag_length - sol.x_length

            px_map[iy, isf] = sol.p_x
            pl_map[iy, isf] = sol.p_left
            pd_map[iy, isf] = sol.p_diag

            rows.append(
                {
                    "Nx": sol.nx,
                    "Ny": sol.ny,
                    "shift": sol.shift,
                    "p_x": sol.p_x,
                    "p_left": sol.p_left,
                    "p_diag": sol.p_diag,
                    "x_length": sol.x_length,
                    "left_length": sol.left_length,
                    "diag_length": sol.diag_length,
                    "left_over_diag": ratio[iy, isf],
                    "triangle_slack": slack[iy, isf],
                }
            )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    csv_path = out_prefix.with_suffix(".csv")
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "Nx",
                "Ny",
                "shift",
                "p_x",
                "p_left",
                "p_diag",
                "x_length",
                "left_length",
                "diag_length",
                "left_over_diag",
                "triangle_slack",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    extent = [sh_vals.min(), sh_vals.max(), ny_vals.min(), ny_vals.max()]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    im0 = axes[0, 0].imshow(left, origin="lower", aspect="auto", extent=extent, cmap="viridis")
    axes[0, 0].set_title("left_length")
    axes[0, 0].set_xlabel("shift")
    axes[0, 0].set_ylabel("Ny")
    fig.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(diag, origin="lower", aspect="auto", extent=extent, cmap="viridis")
    axes[0, 1].set_title("diag_length")
    axes[0, 1].set_xlabel("shift")
    axes[0, 1].set_ylabel("Ny")
    fig.colorbar(im1, ax=axes[0, 1])

    im2 = axes[1, 0].imshow(ratio, origin="lower", aspect="auto", extent=extent, cmap="magma")
    axes[1, 0].set_title("left_length / diag_length")
    axes[1, 0].set_xlabel("shift")
    axes[1, 0].set_ylabel("Ny")
    fig.colorbar(im2, ax=axes[1, 0])

    vmax = np.max(np.abs(slack))
    im3 = axes[1, 1].imshow(
        slack,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    axes[1, 1].set_title("triangle slack = left + diag - x_length")
    axes[1, 1].set_xlabel("shift")
    axes[1, 1].set_ylabel("Ny")
    fig.colorbar(im3, ax=axes[1, 1])

    fig.suptitle(
        f"Analytical compensation heatmaps (Nx={args.Nx}, x_length={args.x_length})",
        fontsize=14,
    )

    png_path = out_prefix.with_suffix(".png")
    svg_path = out_prefix.with_suffix(".svg")
    fig.savefig(png_path, dpi=180)
    fig.savefig(svg_path)
    plt.close(fig)

    # Period heatmaps in a second figure.
    fig2, ax2 = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    impx = ax2[0].imshow(px_map, origin="lower", aspect="auto", extent=extent, cmap="cividis")
    ax2[0].set_title("p_x")
    ax2[0].set_xlabel("shift")
    ax2[0].set_ylabel("Ny")
    fig2.colorbar(impx, ax=ax2[0])

    impl = ax2[1].imshow(pl_map, origin="lower", aspect="auto", extent=extent, cmap="cividis")
    ax2[1].set_title("p_left")
    ax2[1].set_xlabel("shift")
    ax2[1].set_ylabel("Ny")
    fig2.colorbar(impl, ax=ax2[1])

    impd = ax2[2].imshow(pd_map, origin="lower", aspect="auto", extent=extent, cmap="cividis")
    ax2[2].set_title("p_diag")
    ax2[2].set_xlabel("shift")
    ax2[2].set_ylabel("Ny")
    fig2.colorbar(impd, ax=ax2[2])

    fig2.suptitle(f"Orbit periods (Nx={args.Nx})", fontsize=14)
    p_png = out_prefix.parent / (out_prefix.stem + "_periods.png")
    p_svg = out_prefix.parent / (out_prefix.stem + "_periods.svg")
    fig2.savefig(p_png, dpi=180)
    fig2.savefig(p_svg)
    plt.close(fig2)

    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")
    print(f"Wrote {p_png}")
    print(f"Wrote {p_svg}")


if __name__ == "__main__":
    main()
