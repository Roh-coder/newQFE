#!/usr/bin/env python3
"""Draw an equilateral triangular lattice patch and twisted BC identifications."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def lattice_point(x: float, y: float):
    # Equilateral triangular lattice basis: e1=(1,0), e2=(1/2, sqrt(3)/2)
    return np.array([x + 0.5 * y, (np.sqrt(3.0) / 2.0) * y])


def draw_lattice(ax, nx: int, ny: int, step: int):
    xs = list(range(0, nx + 1, step))
    ys = list(range(0, ny + 1, step))

    # Draw basis-parallel segments to communicate triangular connectivity.
    for y in ys:
        for x in xs:
            p = lattice_point(x, y)
            if x + step <= nx:
                q = lattice_point(x + step, y)
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#9aa3ad", lw=0.8, alpha=0.8)
            if y + step <= ny:
                q = lattice_point(x, y + step)
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#9aa3ad", lw=0.8, alpha=0.8)
            if x + step <= nx and y - step >= 0:
                q = lattice_point(x + step, y - step)
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#9aa3ad", lw=0.8, alpha=0.8)

    pts = np.array([lattice_point(x, y) for y in ys for x in xs])
    ax.scatter(pts[:, 0], pts[:, 1], s=8, c="#28323c", alpha=0.8, zorder=3)


def draw_parallelogram_and_bc(ax, nx: int, ny: int, twist_shift: int):
    o = lattice_point(0, 0)
    ex = lattice_point(nx, 0)
    ey = lattice_point(0, ny)
    top = ex + ey

    # Fundamental domain edges, color-coded by periodic identification pair.
    x_pair_color = "#cc3311"  # right edge maps to left edge + shift
    y_pair_color = "#0077bb"  # top edge maps to bottom edge
    ax.plot([o[0], ex[0]], [o[1], ex[1]], color=y_pair_color, lw=2.4)
    ax.plot([ey[0], top[0]], [ey[1], top[1]], color=y_pair_color, lw=2.4)
    ax.plot([o[0], ey[0]], [o[1], ey[1]], color=x_pair_color, lw=2.4)
    ax.plot([ex[0], top[0]], [ex[1], top[1]], color=x_pair_color, lw=2.4)

    # Edge labels to make boundary identification explicit.
    mid_bottom = 0.5 * (o + ex)
    mid_top = 0.5 * (ey + top)
    mid_left = 0.5 * (o + ey)
    mid_right = 0.5 * (ex + top)
    ax.text(mid_bottom[0], mid_bottom[1] - 0.45, "bottom boundary (y=0)",
            color=y_pair_color, fontsize=9, ha="center")
    ax.text(mid_top[0], mid_top[1] + 0.45, "top boundary (y=Ny)",
            color=y_pair_color, fontsize=9, ha="center")
    ax.text(mid_left[0] - 0.55, mid_left[1], "left boundary (x=0)",
            color=x_pair_color, fontsize=9, rotation=60, va="center")
    ax.text(mid_right[0] + 0.45, mid_right[1], "right boundary (x=Nx)",
            color=x_pair_color, fontsize=9, rotation=60, va="center")

    # Period vectors.
    ax.annotate("", xy=ex, xytext=o, arrowprops=dict(arrowstyle="->", lw=2.2, color="#117733"))
    ax.text(*(0.55 * ex + 0.02 * (ey - o)), "Lx = Nx e1", color="#117733", fontsize=10)

    ax.annotate("", xy=ey, xytext=o, arrowprops=dict(arrowstyle="->", lw=2.2, color="#aa3377"))
    ax.text(*(0.45 * ey + np.array([0.25, 0.0])), "Ly = Ny e2", color="#aa3377", fontsize=10)

    # Twisted x-boundary map: right edge maps to left edge shifted by twist in y.
    y0 = int(round(0.25 * ny))
    p_right = lattice_point(nx, y0)
    p_left_target = lattice_point(0, y0 + twist_shift)
    ax.annotate("", xy=p_left_target, xytext=p_right,
                arrowprops=dict(arrowstyle="->", lw=2.4, color="#cc3311", linestyle="--"))

    # y-periodic map marker.
    x0 = int(round(0.22 * nx))
    p_top = lattice_point(x0, ny)
    p_bottom = lattice_point(x0, 0)
    ax.annotate("", xy=p_bottom, xytext=p_top,
                arrowprops=dict(arrowstyle="->", lw=2.2, color="#0077bb", linestyle=":"))

    ax.text(0.5 * (p_right[0] + p_left_target[0]),
            0.5 * (p_right[1] + p_left_target[1]) + 0.28,
            f"x-wrap map: (Nx,{y0}) -> (0,{(y0 + twist_shift) % ny})",
            color="#cc3311", fontsize=9, ha="center")

    ax.text(0.5 * (p_top[0] + p_bottom[0]) - 0.25,
            0.5 * (p_top[1] + p_bottom[1]) - 0.1,
            f"y-wrap map: ({x0},Ny) -> ({x0},0)",
            color="#0077bb", fontsize=9, rotation=90, va="center")

    label = (
        f"Twisted BC:\n"
        f"(x+Nx, y) ~ (x, y+{twist_shift})\n"
        f"(x, y+Ny) ~ (x, y)"
    )
    anchor = top + np.array([0.7, -0.2])
    ax.text(anchor[0], anchor[1], label, fontsize=10,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#666666", alpha=0.95))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, required=True)
    ap.add_argument("--ny", type=int, required=True)
    ap.add_argument("--twist_shift", type=int, required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--sample_step", type=int, default=2)
    ap.add_argument("--title", default="Equilateral triangular lattice with twisted BC")
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    draw_lattice(ax, args.nx, args.ny, max(1, args.sample_step))
    draw_parallelogram_and_bc(ax, args.nx, args.ny, args.twist_shift)

    ax.set_title(args.title)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.tight_layout()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=190)
    print(out)


if __name__ == "__main__":
    main()
