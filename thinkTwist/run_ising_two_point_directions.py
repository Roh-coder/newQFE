#!/usr/bin/env python3
"""Measure Ising two-point functions in three triangular directions.

Lattice BC:
  (i + Nx, j) ~ (i, j)
  (i, j + Ny) ~ (i + shift, j)

Directions:
  x    = (1, 0)
  left = (0, 1)
  diag = (-1, 1)
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def canonicalize(i: int, j: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    j0 = j % ny
    b = (j - j0) // ny
    i0 = (i - b * shift) % nx
    return i0, j0


def add_bc(i: int, j: int, di: int, dj: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    return canonicalize(i + di, j + dj, nx, ny, shift)


def orbit_period(step: tuple[int, int], nx: int, ny: int, shift: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, nx, ny, shift)
    for n in range(1, nx * ny + 1):
        x, y = canonicalize(x + sx, y + sy, nx, ny, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period for step={step}")


def make_neighbors(nx: int, ny: int, shift: int) -> list[list[list[tuple[int, int]]]]:
    dirs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    neigh = [[None for _ in range(ny)] for _ in range(nx)]
    for i in range(nx):
        for j in range(ny):
            neigh[i][j] = [add_bc(i, j, di, dj, nx, ny, shift) for di, dj in dirs]
    return neigh


def metropolis_sweep(spins: list[list[int]], neigh: list[list[list[tuple[int, int]]]], K: float, nx: int, ny: int) -> None:
    nsite = nx * ny
    boltz = {4: math.exp(-2.0 * K * 4), 8: math.exp(-2.0 * K * 8), 12: math.exp(-2.0 * K * 12)}
    for _ in range(nsite):
        i = random.randrange(nx)
        j = random.randrange(ny)
        s = spins[i][j]
        h = 0
        for ni, nj in neigh[i][j]:
            h += spins[ni][nj]
        dE = 2 * s * h
        if dE <= 0 or random.random() < boltz.get(dE, math.exp(-K * dE)):
            spins[i][j] = -s


def measure_direction(
    spins: list[list[int]],
    anchors: list[tuple[int, int]],
    step: tuple[int, int],
    p: int,
    nx: int,
    ny: int,
    shift: int,
) -> np.ndarray:
    out = np.zeros(p + 1, dtype=float)
    dx, dy = step
    for n in range(p + 1):
        csum = 0.0
        for i, j in anchors:
            ni, nj = add_bc(i, j, n * dx, n * dy, nx, ny, shift)
            csum += spins[i][j] * spins[ni][nj]
        out[n] = csum / len(anchors)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Nx", type=int, default=32)
    ap.add_argument("--Ny", type=int, default=32)
    ap.add_argument("--shift", type=int, default=16)
    ap.add_argument("--K", type=float, default=None)
    ap.add_argument("--therm", type=int, default=2000)
    ap.add_argument("--meas", type=int, default=4000)
    ap.add_argument("--sample-every", type=int, default=10)
    ap.add_argument("--anchors", type=int, default=256)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out-prefix", default="thinkTwist/ising_two_point_Nx32_Ny32_shift16")
    args = ap.parse_args()

    random.seed(args.seed)
    nx, ny, shift = args.Nx, args.Ny, args.shift
    Kc_tri = 0.25 * math.log(3.0)
    K = Kc_tri if args.K is None else args.K

    steps = {
        "x": (1, 0),
        "left": (0, 1),
        "diag": (-1, 1),
    }
    periods = {k: orbit_period(v, nx, ny, shift) for k, v in steps.items()}

    spins = [[1 if random.random() < 0.5 else -1 for _ in range(ny)] for _ in range(nx)]
    neigh = make_neighbors(nx, ny, shift)

    all_sites = [(i, j) for i in range(nx) for j in range(ny)]
    if args.anchors >= len(all_sites):
        anchors = all_sites
    else:
        anchors = random.sample(all_sites, args.anchors)

    for _ in range(args.therm):
        metropolis_sweep(spins, neigh, K, nx, ny)

    acc = {k: np.zeros(periods[k] + 1, dtype=float) for k in steps}
    nsamp = 0

    for sweep in range(1, args.meas + 1):
        metropolis_sweep(spins, neigh, K, nx, ny)
        if sweep % args.sample_every != 0:
            continue
        nsamp += 1
        for k in steps:
            acc[k] += measure_direction(spins, anchors, steps[k], periods[k], nx, ny, shift)

    if nsamp == 0:
        raise RuntimeError("No samples collected; reduce --sample-every or increase --meas")

    avg = {k: acc[k] / nsamp for k in steps}

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # Write text summary
    txt = out_prefix.with_suffix(".txt")
    with txt.open("w") as f:
        f.write("MC two-point functions in triangular directions\n")
        f.write("============================================\n")
        f.write(f"Nx={nx} Ny={ny} shift={shift} K={K:.12f}\n")
        f.write(f"Kc_tri={(0.25 * math.log(3.0)):.12f}\n")
        f.write(f"therm={args.therm} meas={args.meas} sample_every={args.sample_every}\n")
        f.write(f"anchors={len(anchors)} nsamp={nsamp}\n")
        f.write(f"periods: x={periods['x']} left={periods['left']} diag={periods['diag']}\n")

    # Combined plot
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 10), sharex=False)
    colors = {"x": "#1f77b4", "left": "#d62728", "diag": "#2ca02c"}

    for ax, key in zip(axes, ["x", "left", "diag"]):
        p = periods[key]
        frac = np.arange(0, p + 1) / p
        ax.plot(frac, avg[key], "o-", color=colors[key], ms=3)
        ax.set_title(f"{key} direction, period={p}")
        ax.set_ylabel("G(n)=<s(0)s(n)>")
        ax.set_ylim(-1.05, 1.05)
        ax.grid(alpha=0.25)

    axes[-1].set_xlabel("fractional separation n/p")
    fig.suptitle(f"Ising two-point functions (Nx={nx}, Ny={ny}, shift={shift}, K={K:.6f})")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    png = out_prefix.with_suffix(".png")
    svg = out_prefix.with_suffix(".svg")
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    # Separate direction plots
    for key in ["x", "left", "diag"]:
        p = periods[key]
        frac = np.arange(0, p + 1) / p
        f2, a2 = plt.subplots(figsize=(8, 4.5))
        a2.plot(frac, avg[key], "o-", color=colors[key], ms=3)
        a2.set_title(f"Ising two-point: {key} direction (period={p})")
        a2.set_xlabel("fractional separation n/p")
        a2.set_ylabel("G(n)")
        a2.set_ylim(-1.05, 1.05)
        a2.grid(alpha=0.25)
        f2.tight_layout()
        f2.savefig(out_prefix.parent / f"{out_prefix.stem}_{key}.png", dpi=180)
        f2.savefig(out_prefix.parent / f"{out_prefix.stem}_{key}.svg")
        plt.close(f2)

    print(f"Wrote {txt}")
    print(f"Wrote {png}")
    print(f"Wrote {svg}")
    print(f"periods: x={periods['x']} left={periods['left']} diag={periods['diag']} nsamp={nsamp}")


if __name__ == "__main__":
    main()
