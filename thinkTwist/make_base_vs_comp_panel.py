#!/usr/bin/env python3
"""Build a 3x2 panel comparing base vs compensated directional plots.

Rows: base, left, diag
Cols: base (equal lengths), compensated
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def canonicalize(x: int, y: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    y0 = y % ny
    b = (y - y0) // ny
    x0 = (x - b * shift) % nx
    return x0, y0


def orbit_period(step: tuple[int, int], nx: int, ny: int, shift: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, nx, ny, shift)
    for n in range(1, nx * ny + 1):
        x, y = canonicalize(x + sx, y + sy, nx, ny, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period found for step={step}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=32)
    parser.add_argument("--Nx", type=int, default=None)
    parser.add_argument("--Ny", type=int, default=None)
    parser.add_argument("--shift", type=int, default=1)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument("--in-dir", default="thinkTwist")
    parser.add_argument("--out", default="thinkTwist/two_point_base_vs_comp_panel_L32_shift1.png")
    args = parser.parse_args()

    in_dir = Path(args.in_dir)

    nx = args.Nx if args.Nx is not None else args.L
    ny = args.Ny if args.Ny is not None else args.L

    rows = ["base", "left", "diag"]
    if nx == ny:
        prefix = f"twisted_base_L{nx}_shift{args.shift}_one_period_fractional_two_point"
    else:
        prefix = f"twisted_base_Nx{nx}_Ny{ny}_shift{args.shift}_one_period_fractional_two_point"

    base_paths = [
        in_dir / f"{prefix}_{r}_BASE_equal.png"
        for r in rows
    ]
    comp_paths = [
        in_dir / f"{prefix}_{r}_COMPENSATED.png"
        for r in rows
    ]

    for p in base_paths + comp_paths:
        if not p.exists():
            raise FileNotFoundError(f"Missing image: {p}")

    fig, axes = plt.subplots(4, 2, figsize=(12, 17))

    for i, r in enumerate(rows):
        img_base = mpimg.imread(base_paths[i])
        img_comp = mpimg.imread(comp_paths[i])

        ax_b = axes[i, 0]
        ax_c = axes[i, 1]

        ax_b.imshow(img_base)
        ax_c.imshow(img_comp)

        ax_b.set_title(f"{r} - BASE", fontsize=11)
        ax_c.set_title(f"{r} - COMPENSATED", fontsize=11)

        ax_b.axis("off")
        ax_c.axis("off")

    # Seventh panel: compensation triangle (or infeasibility diagnostic).
    tri_ax = axes[3, 0]
    px = orbit_period((1, 0), nx, ny, args.shift)
    pl = orbit_period((0, 1), nx, ny, args.shift)
    pd = orbit_period((-1, 1), nx, ny, args.shift)

    # Exact centered-cosh matching lengths with ell_x fixed.
    lx = args.x_length
    ly = args.x_length * px / pl
    lr = args.x_length * px / pd

    lower = abs(lx - ly)
    upper = lx + ly
    feasible = lower <= lr <= upper

    if feasible:
        cos_theta = (lx * lx + ly * ly - lr * lr) / (2.0 * lx * ly)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        theta = math.acos(cos_theta)
        ex = (lx, 0.0)
        ey = (ly * math.cos(theta), ly * math.sin(theta))
        tri_ax.plot([0.0, ex[0]], [0.0, ex[1]], color="#1f77b4", lw=3, label="|x|")
        tri_ax.plot([0.0, ey[0]], [0.0, ey[1]], color="#d62728", lw=3, label="|left|")
        tri_ax.plot([ex[0], ey[0]], [ex[1], ey[1]], color="#2ca02c", lw=3, label="|diag|")
        tri_ax.scatter([0.0, ex[0], ey[0]], [0.0, ex[1], ey[1]], c="black", s=20)
        tri_ax.set_aspect("equal")
        tri_ax.grid(alpha=0.25)
        tri_ax.legend(loc="upper right", fontsize=8)
        tri_ax.set_title("Compensation Triangle")
    else:
        # Draw the degenerate boundary triangle explicitly as collinear points.
        # Place O=(0,0), X=(lx,0), and Y on x-axis so |OY|=ly and |XY|=lr.
        # For this boundary case (lr=|lx-ly|) Y lies between O and X (or beyond).
        yx = lx - lr
        ex = (lx, 0.0)
        ey = (yx, 0.0)
        tri_ax.plot([0.0, ex[0]], [0.0, ex[1]], color="#1f77b4", lw=3, label="|x|")
        tri_ax.plot([0.0, ey[0]], [0.0, ey[1]], color="#d62728", lw=3, label="|left|")
        tri_ax.plot([ex[0], ey[0]], [ex[1], ey[1]], color="#2ca02c", lw=3, label="|diag|")
        tri_ax.scatter([0.0, ex[0], ey[0]], [0.0, ex[1], ey[1]], c="black", s=20)
        tri_ax.text(0.0, 0.02, "O", fontsize=9, ha="center")
        tri_ax.text(ex[0], 0.02, "X", fontsize=9, ha="center")
        tri_ax.text(ey[0], -0.035, "Y", fontsize=9, ha="center")
        tri_ax.set_aspect("equal")
        tri_ax.grid(alpha=0.25)
        tri_ax.legend(loc="upper right", fontsize=8)
        tri_ax.set_title("Compensation Triangle (Degenerate Boundary Case)")
        tri_ax.text(
            0.02,
            0.97,
            f"Degenerate: diag=|x-left|={lr:.5f}",
            transform=tri_ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="#f7f7f7", ec="#cccccc"),
        )

    # Eighth slot left intentionally as notes panel.
    note_ax = axes[3, 1]
    note_ax.axis("off")
    note_ax.text(
        0.03,
        0.95,
        "Seventh panel: compensation triangle/lengths\n"
        "(base: left column, compensated: right column)",
        ha="left",
        va="top",
        fontsize=10,
    )

    fig.suptitle(
        f"Two-point comparison panel (L={args.L}, shift={args.shift})\n"
        f"(Nx={nx}, Ny={ny})\n"
        "Rows 1-3: base/left/diag comparisons. Row 4: compensation shape panel",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)

    print(f"Wrote {out_path}")
    print(f"Wrote {out_path.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
