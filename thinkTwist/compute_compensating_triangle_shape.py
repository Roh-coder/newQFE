#!/usr/bin/env python3
"""Compute and plot a triangle deformation that compensates twisted BC periods.

We use twisted BC:
  (i + L, j) ~ (i, j)
  (i, j + L) ~ (i + shift, j)

Directions are interpreted as:
  x  : e1 = (1, 0)
  y  : e2 = (0, 1)
  r  : e3 = y - x = (-1, 1)

Goal: choose local metric lengths (|x|, |y|, |r|) so that half-orbit physical
scales mimic square-continuum expectations for x, y, y-x as closely as possible.
Here we match the half-orbit target exactly.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def is_equivalent_to_zero(x: int, y: int, L: int, shift: int) -> bool:
    if y % L != 0:
        return False
    b = y // L
    rem_x = x - b * shift
    return rem_x % L == 0


def period(step: tuple[int, int], L: int, shift: int, max_n: int = 10000) -> int:
    sx, sy = step
    for n in range(1, max_n + 1):
        if is_equivalent_to_zero(n * sx, n * sy, L, shift):
            return n
    raise RuntimeError(f"No period found for step={step}")


def folded_distance_steps(r: np.ndarray, p: int) -> np.ndarray:
    rr = np.mod(r, p)
    return np.minimum(rr, p - rr)


def main() -> None:
    # Problem setup
    L = 64
    shift = 32
    xi = 10.0
    r = np.arange(0, 321)

    steps = {
        "x": (1, 0),
        "y": (0, 1),
        "y-x": (-1, 1),
    }

    px = period(steps["x"], L, shift)
    py = period(steps["y"], L, shift)
    pr = period(steps["y-x"], L, shift)

    # Match condition to square-continuum half-orbit scales:
    # twisted: (p_dir/2)*|dir|, target ratio x:y:(y-x) = 1:1:sqrt(2)
    # Fix |x| = 1 (scale choice), then solve for |y| and |y-x|.
    lx = 1.0
    ly = (px / py) * lx
    lr = (math.sqrt(2.0) * px / pr) * lx

    # Infer triangle angle between x and y from |y-x|^2 = |x|^2 + |y|^2 - 2|x||y|cos(theta)
    cos_theta = (lx * lx + ly * ly - lr * lr) / (2.0 * lx * ly)
    if abs(cos_theta) > 1.0:
        raise RuntimeError("No geometric triangle realization: |cos(theta)| > 1")
    theta_deg = math.degrees(math.acos(cos_theta))

    # Compare against right-isosceles guess (known to fail here)
    lxi, lyi, lri = 1.0, 1.0, math.sqrt(2.0)

    # Square periodic reference for schematic comparison (no twist)
    p_sq = L

    out_dir = Path("thinkTwist")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write summary report
    report = out_dir / "twisted_L64_shape_fit_report.txt"
    with report.open("w") as f:
        f.write("Twisted BC compensating triangle fit (L=64, shift=32)\n")
        f.write("====================================================\n\n")
        f.write(f"Periods under twist: px={px}, py={py}, p(y-x)={pr}\n")
        f.write("Target half-orbit ratio (square continuum): x:y:(y-x) = 1:1:sqrt(2)\n\n")
        f.write("Chosen compensating local lengths (scale |x|=1):\n")
        f.write(f"  |x|      = {lx:.6f}\n")
        f.write(f"  |y|      = {ly:.6f}\n")
        f.write(f"  |y-x|    = {lr:.6f}\n")
        f.write(f"  cos(theta_xy) = {cos_theta:.6f}\n")
        f.write(f"  theta_xy       = {theta_deg:.6f} deg\n\n")
        f.write("Equivalent side ratio (x : y : y-x):\n")
        f.write(f"  1 : {ly/lx:.6f} : {lr/lx:.6f}\n\n")
        f.write("Right-isosceles local lengths for comparison:\n")
        f.write(f"  |x|=|y|={lxi:.6f}, |y-x|={lri:.6f}\n")

    # Build directional correlators
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # (name, period_twisted, length_twisted_comp, length_right_iso, color)
    cfg = [
        ("x", px, lx, lxi, "#1f77b4"),
        ("y", py, ly, lyi, "#d62728"),
        ("y-x", pr, lr, lri, "#2ca02c"),
    ]

    for ax, (name, p_tw, l_comp, l_iso, color) in zip(axes, cfg):
        d_comp = folded_distance_steps(r, p_tw) * l_comp
        c_comp = np.exp(-d_comp / xi)

        d_iso = folded_distance_steps(r, p_tw) * l_iso
        c_iso = np.exp(-d_iso / xi)

        if name == "x":
            d_ref = folded_distance_steps(r, p_sq) * 1.0
        elif name == "y":
            d_ref = folded_distance_steps(r, p_sq) * 1.0
        else:
            d_ref = folded_distance_steps(r, p_sq) * math.sqrt(2.0)
        c_ref = np.exp(-d_ref / xi)

        ax.plot(r, c_ref, "k--", lw=1.8, label="square periodic reference")
        ax.plot(r, c_iso, color="#999999", lw=2.0, label="twisted + right-isosceles")
        ax.plot(r, c_comp, color=color, lw=2.4, label="twisted + compensating shape")

        ax.set_title(f"Direction {name}: p_twisted={p_tw}")
        ax.set_ylabel("C(r)")
        ax.set_ylim(-0.02, 1.05)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")

    axes[-1].set_xlabel("trajectory separation r (steps)")
    fig.suptitle(
        "L=64, shift=32 twisted BC: compensating triangle-shape test\n"
        "C(r) = exp(-d_eff/xi), xi=10",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    png = out_dir / "twisted_L64_shape_fit_comparison.png"
    svg = out_dir / "twisted_L64_shape_fit_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)

    # Also draw the triangle itself in physical embedding coordinates
    tri_fig, tri_ax = plt.subplots(figsize=(6.5, 5.5))
    ex = np.array([lx, 0.0])
    ey = np.array([ly * math.cos(math.radians(theta_deg)), ly * math.sin(math.radians(theta_deg))])
    e3 = ey - ex

    O = np.array([0.0, 0.0])
    X = ex
    Y = ey

    tri_ax.plot([O[0], X[0]], [O[1], X[1]], "-", color="#1f77b4", lw=3, label="x edge")
    tri_ax.plot([O[0], Y[0]], [O[1], Y[1]], "-", color="#d62728", lw=3, label="y edge")
    tri_ax.plot([X[0], Y[0]], [X[1], Y[1]], "-", color="#2ca02c", lw=3, label="y-x edge")

    tri_ax.scatter([O[0], X[0], Y[0]], [O[1], X[1], Y[1]], c="black", s=35)
    tri_ax.text((O[0] + X[0]) * 0.5, (O[1] + X[1]) * 0.5 - 0.03, f"|x|={lx:.3f}", color="#1f77b4")
    tri_ax.text((O[0] + Y[0]) * 0.5 - 0.06, (O[1] + Y[1]) * 0.5, f"|y|={ly:.3f}", color="#d62728")
    tri_ax.text((X[0] + Y[0]) * 0.5 + 0.02, (X[1] + Y[1]) * 0.5, f"|y-x|={lr:.3f}", color="#2ca02c")

    tri_ax.set_title(f"Compensating triangle shape (theta_xy={theta_deg:.2f} deg)")
    tri_ax.set_aspect("equal")
    tri_ax.grid(alpha=0.2)
    tri_ax.legend(loc="upper right")

    tri_png = out_dir / "twisted_L64_compensating_triangle_shape.png"
    tri_svg = out_dir / "twisted_L64_compensating_triangle_shape.svg"
    tri_fig.savefig(tri_png, dpi=180)
    tri_fig.savefig(tri_svg)

    print(f"Wrote {report}")
    print(f"Wrote {png}")
    print(f"Wrote {svg}")
    print(f"Wrote {tri_png}")
    print(f"Wrote {tri_svg}")
    print(
        "Compensating shape summary: "
        f"|x|={lx:.6f}, |y|={ly:.6f}, |y-x|={lr:.6f}, theta={theta_deg:.6f} deg"
    )


if __name__ == "__main__":
    main()
