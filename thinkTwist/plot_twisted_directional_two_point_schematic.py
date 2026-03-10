#!/usr/bin/env python3
"""Minimal two-point decays on an LxL triangular lattice with twisted BC.

BC used here:
- (i + L, j) ~ (i, j)
- (i, j + L) ~ (i + shift, j)
Typical Mobius-like case uses shift=L/2.

For each triangular lattice direction, we compute the trajectory period on this
quotient lattice and then model a smooth cosh profile centered at half-period,
with direction-dependent length scale:

    C_d(k) = cosh((k-p_d/2)*ell_d/a) / cosh((p_d/2)*ell_d/a)

This gives C_d(0)=C_d(p_d)=1 and a smooth minimum at k=p_d/2.

where ell_d is the side length for that direction and a is the base length.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class Direction:
    name: str
    step: tuple[int, int]
    color: str


def is_equivalent_to_zero(vec: tuple[int, int], nx: int, ny: int, shift: int) -> bool:
    """Return True if vec is in the BC translation lattice.

    Translation lattice generators:
      g1 = (L, 0)
      g2 = (shift, L)

    Need integer a,b so vec = a*g1 + b*g2.
    """

    x, y = vec
    if y % ny != 0:
        return False
    b = y // ny
    rem_x = x - b * shift
    return rem_x % nx == 0


def orbit_period(step: tuple[int, int], nx: int, ny: int, shift: int, max_n: int | None = None) -> int:
    """Smallest n>0 such that n*step is equivalent to zero under BC."""

    sx, sy = step
    if max_n is None:
        # Quotient has at most nx*ny states, so period is bounded by nx*ny.
        max_n = nx * ny
    for n in range(1, max_n + 1):
        if is_equivalent_to_zero((n * sx, n * sy), nx, ny, shift):
            return n
    raise RuntimeError(f"No period found for step={step} up to n={max_n}")


def folded_distance(r: np.ndarray, period: int) -> np.ndarray:
    """Periodic folded 1D distance along a closed orbit of length period."""

    rem = np.mod(r, period)
    return np.minimum(rem, period - rem)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=4, help="Lattice linear size (legacy; sets Nx=Ny=L)")
    parser.add_argument("--Nx", type=int, default=None, help="Lattice size in x")
    parser.add_argument("--Ny", type=int, default=None, help="Lattice size in y")
    parser.add_argument("--shift", type=int, default=2, help="Twist shift in lattice units")
    parser.add_argument(
        "--rmax",
        type=int,
        default=None,
        help="Maximum trajectory separation shown (default: 5*L)",
    )
    parser.add_argument(
        "--fractional-one-period",
        action="store_true",
        help="Plot n/p_d over exactly one period (n=0..p_d-1)",
    )
    parser.add_argument(
        "--x-length",
        type=float,
        default=1.0,
        help="x-direction step length in units of a",
    )
    parser.add_argument(
        "--left-length",
        type=float,
        default=1.0,
        help="left-direction step length in units of a",
    )
    parser.add_argument(
        "--diag-length",
        type=float,
        default=np.sqrt(2.0),
        help="diag-direction step length in units of a",
    )
    parser.add_argument(
        "--decay-scale",
        type=float,
        default=1.0,
        help="Divide exponent by this factor (e.g. 10 makes decay 10x less aggressive)",
    )
    args = parser.parse_args()

    nx = args.Nx if args.Nx is not None else args.L
    ny = args.Ny if args.Ny is not None else args.L
    shift = args.shift

    if nx <= 0 or ny <= 0:
        raise ValueError("Nx and Ny must be positive")
    if shift < 0:
        raise ValueError("shift must be non-negative")
    if args.decay_scale <= 0:
        raise ValueError("decay-scale must be positive")

    length_over_a = {
        "base": args.x_length,
        "left": args.left_length,
        "diag": args.diag_length,
    }

    directions = [
        Direction("base e1 = (1,0)", (1, 0), "#1f77b4"),
        Direction("left e2 = (0,1)", (0, 1), "#d62728"),
        Direction("diag e3 = (-1,1)", (-1, 1), "#2ca02c"),
    ]

    # Global r grid for standard mode; in fractional one-period mode this is per direction.
    r_max = args.rmax if args.rmax is not None else 5 * max(nx, ny)
    r_global = np.arange(0, r_max + 1)

    out_dir = Path("thinkTwist")
    out_dir.mkdir(parents=True, exist_ok=True)

    for d in directions:
        p = orbit_period(d.step, nx, ny, shift)
        if args.fractional_one_period:
            # Include the endpoint n=p (same state as n=0) to show a full closed period.
            r = np.arange(0, p + 1)
            xvals = r / p
            xlabel = "fractional separation n/p_d (one period)"
        else:
            r = r_global
            xvals = r
            xlabel = "trajectory separation r (lattice steps)"

        # Use periodic separation over a full cycle k=0..p.
        # Decay rule remains C_d(n)=exp(-n*ell_d/a), with n taken as folded
        # periodic separation so the plotted one-period curve is complete.
        if args.fractional_one_period:
            n_eff = folded_distance(r, p)
        else:
            n_eff = r.astype(float)
        slug = d.name.split()[0].replace("=", "").replace("(", "").replace(")", "")
        ell = length_over_a[slug]
        phase = np.mod(r, p).astype(float)
        centered = phase - 0.5 * p
        alpha = ell / args.decay_scale
        c = np.cosh(centered * alpha) / np.cosh(0.5 * p * alpha)

        fig, axes = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)

        axes[0].plot(xvals, n_eff, "o-", color=d.color)
        axes[0].set_ylabel("separation n")
        axes[0].set_title(f"{d.name}  (period={p})")
        axes[0].grid(alpha=0.25)

        axes[1].plot(xvals, c, "o-", color=d.color)
        axes[1].set_xlabel(xlabel)
        axes[1].set_ylabel("C_d(k): cosh centered at k=p_d/2")
        axes[1].set_ylim(-0.02, 1.05)
        axes[1].grid(alpha=0.25)

        # Minimal single-scale form in folded step units.
        formula = (
            "C_d(k)=cosh(((k-p_d/2)*ell_d/a)/s) / cosh((p_d/2)*ell_d/a/s)"
            + "\n"
            + f"for this panel: ell_d/a={ell:.3g}, p_d={p}, s={args.decay_scale:.3g}"
        )
        axes[1].text(
            0.02,
            0.03,
            formula,
            transform=axes[1].transAxes,
            ha="left",
            va="bottom",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="#f7f7f7", ec="#cccccc", alpha=0.9),
        )

        fig.suptitle(
            f"Twisted BC schematic (Nx={nx}, Ny={ny}, shift={shift})\n"
            "rule: cosh profile centered at half-period with direction-dependent ell_d",
            fontsize=11,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        if args.fractional_one_period:
            if nx == ny:
                stem = f"twisted_base_L{nx}_shift{shift}_one_period_fractional_two_point_{slug}"
            else:
                stem = f"twisted_base_Nx{nx}_Ny{ny}_shift{shift}_one_period_fractional_two_point_{slug}"
        else:
            if nx == ny:
                stem = f"twisted_base_L{nx}_shift{shift}_two_point_{slug}"
            else:
                stem = f"twisted_base_Nx{nx}_Ny{ny}_shift{shift}_two_point_{slug}"
        png_out = out_dir / f"{stem}.png"
        svg_out = out_dir / f"{stem}.svg"

        fig.savefig(png_out, dpi=170)
        fig.savefig(svg_out)
        plt.close(fig)

        print(f"Wrote {png_out}")
        print(f"Wrote {svg_out}")


if __name__ == "__main__":
    main()
