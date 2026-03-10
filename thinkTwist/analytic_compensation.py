#!/usr/bin/env python3
"""Analytical compensation lengths for twisted boundary conditions.

Boundary conditions:
  (i + Nx, j) ~ (i, j)
  (i, j + Ny) ~ (i + shift, j)

For centered-cosh correlators in fractional separation t = k / p_d,
exact directional matching is obtained when

  p_x * ell_x = p_left * ell_left = p_diag * ell_diag,

where p_x, p_left, p_diag are the orbit periods of directions
(1,0), (0,1), and (-1,1) under the quotient BC.

Given ell_x, this yields
  ell_left = ell_x * p_x / p_left
  ell_diag = ell_x * p_x / p_diag
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class CompensationResult:
    nx: int
    ny: int
    shift: int
    p_x: int
    p_left: int
    p_diag: int
    x_length: float
    left_length: float
    diag_length: float
    invariant: float


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
    raise RuntimeError(
        f"No period found for step={step}, nx={nx}, ny={ny}, shift={shift}"
    )


def solve_centered_cosh_lengths(
    nx: int,
    ny: int,
    shift: int,
    x_length: float = 1.0,
) -> CompensationResult:
    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be positive")
    if x_length <= 0:
        raise ValueError("x_length must be positive")

    p_x = orbit_period((1, 0), nx, ny, shift)
    p_left = orbit_period((0, 1), nx, ny, shift)
    p_diag = orbit_period((-1, 1), nx, ny, shift)

    left_length = x_length * p_x / p_left
    diag_length = x_length * p_x / p_diag

    # Shared matching product for all directions.
    invariant = p_x * x_length

    return CompensationResult(
        nx=nx,
        ny=ny,
        shift=shift,
        p_x=p_x,
        p_left=p_left,
        p_diag=p_diag,
        x_length=x_length,
        left_length=left_length,
        diag_length=diag_length,
        invariant=invariant,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Nx", type=int, required=True)
    parser.add_argument("--Ny", type=int, required=True)
    parser.add_argument("--shift", type=int, required=True)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format",
    )
    args = parser.parse_args()

    result = solve_centered_cosh_lengths(
        nx=args.Nx,
        ny=args.Ny,
        shift=args.shift,
        x_length=args.x_length,
    )

    if args.format == "json":
        print(json.dumps(asdict(result), indent=2, sort_keys=True))
        return

    print("Centered-cosh analytical compensation")
    print("====================================")
    print(f"Nx={result.nx}, Ny={result.ny}, shift={result.shift}")
    print(f"periods: p_x={result.p_x}, p_left={result.p_left}, p_diag={result.p_diag}")
    print(f"x_length={result.x_length:.12g}")
    print(f"left_length={result.left_length:.12g}")
    print(f"diag_length={result.diag_length:.12g}")
    print(
        "check: p_x*x_length = p_left*left_length = p_diag*diag_length = "
        f"{result.invariant:.12g}"
    )


if __name__ == "__main__":
    main()
