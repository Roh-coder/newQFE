#!/usr/bin/env python3
"""Metropolis Ising on twisted Nx x Ny triangular lattice.

Boundary conditions:
  (i + Nx, j) ~ (i, j)
  (i, j + Ny) ~ (i + shift, j)

Triangular nearest neighbors are taken along basis directions:
  +/-e1, +/-e2, +e1-e2, -e1+e2

Critical coupling for isotropic ferromagnetic triangular Ising:
  Kc = (1/4) * ln(3)
  equivalently sinh(2Kc) = 1/sqrt(3)
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path


def canonicalize(i: int, j: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    j0 = j % ny
    b = (j - j0) // ny
    i0 = (i - b * shift) % nx
    return i0, j0


def add_bc(i: int, j: int, di: int, dj: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    return canonicalize(i + di, j + dj, nx, ny, shift)


def neighbor_list(nx: int, ny: int, shift: int) -> list[list[tuple[int, int]]]:
    neigh = [[None for _ in range(ny)] for _ in range(nx)]
    dirs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    for i in range(nx):
        for j in range(ny):
            neigh[i][j] = [add_bc(i, j, di, dj, nx, ny, shift) for di, dj in dirs]
    return neigh


def total_energy(spins: list[list[int]], neigh: list[list[list[tuple[int, int]]]]) -> float:
    nx = len(spins)
    ny = len(spins[0])
    e = 0.0
    # Count each bond once via three forward directions.
    forward_idx = [0, 2, 4]  # (1,0), (0,1), (1,-1)
    for i in range(nx):
        for j in range(ny):
            s = spins[i][j]
            for idx in forward_idx:
                ni, nj = neigh[i][j][idx]
                e -= s * spins[ni][nj]
    return e


def run_mc(nx: int, ny: int, shift: int, K: float, therm: int, meas: int, seed: int) -> dict[str, float]:
    random.seed(seed)
    nsite = nx * ny

    spins = [[1 if random.random() < 0.5 else -1 for _ in range(ny)] for _ in range(nx)]
    neigh = neighbor_list(nx, ny, shift)

    boltz = {4: math.exp(-2.0 * K * 4), 8: math.exp(-2.0 * K * 8), 12: math.exp(-2.0 * K * 12)}

    def sweep() -> None:
        for _ in range(nsite):
            i = random.randrange(nx)
            j = random.randrange(ny)
            s = spins[i][j]
            h = 0
            for ni, nj in neigh[i][j]:
                h += spins[ni][nj]
            dE = 2 * s * h
            if dE <= 0 or random.random() < boltz.get(dE, math.exp(-2.0 * K * dE / 2.0)):
                spins[i][j] = -s

    for _ in range(therm):
        sweep()

    e_acc = 0.0
    e2_acc = 0.0
    m_acc = 0.0
    mabs_acc = 0.0
    m2_acc = 0.0

    for _ in range(meas):
        sweep()
        e = total_energy(spins, neigh)
        m = sum(sum(row) for row in spins)

        e_acc += e
        e2_acc += e * e
        m_acc += m
        mabs_acc += abs(m)
        m2_acc += m * m

    norm = 1.0 / meas
    e_mean = e_acc * norm
    e2_mean = e2_acc * norm
    m_mean = m_acc * norm
    mabs_mean = mabs_acc * norm
    m2_mean = m2_acc * norm

    e_site = e_mean / nsite
    m_site = m_mean / nsite
    mabs_site = mabs_mean / nsite

    C = (e2_mean - e_mean * e_mean) * (K * K) / nsite
    chi = (m2_mean - m_mean * m_mean) * K / nsite

    return {
        "Nx": float(nx),
        "Ny": float(ny),
        "shift": float(shift),
        "K": K,
        "therm_sweeps": float(therm),
        "meas_sweeps": float(meas),
        "E_per_site": e_site,
        "M_per_site": m_site,
        "|M|_per_site": mabs_site,
        "C_per_site": C,
        "chi_per_site": chi,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--Nx", type=int, default=32)
    p.add_argument("--Ny", type=int, default=32)
    p.add_argument("--shift", type=int, default=16)
    p.add_argument("--K", type=float, default=None)
    p.add_argument("--therm", type=int, default=2000)
    p.add_argument("--meas", type=int, default=8000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--out", default="thinkTwist/ising_triangular_twist_L32_shift16.txt")
    args = p.parse_args()

    Kc_tri = 0.25 * math.log(3.0)
    K = Kc_tri if args.K is None else args.K

    res = run_mc(args.Nx, args.Ny, args.shift, K, args.therm, args.meas, args.seed)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("Triangular Ising on twisted lattice (Metropolis)\n")
        f.write("==============================================\n")
        f.write(f"Nx={args.Nx} Ny={args.Ny} shift={args.shift}\n")
        f.write(f"K={K:.12f}\n")
        f.write(f"Kc_tri=(1/4)ln3={Kc_tri:.12f}\n")
        f.write(f"sinh(2K)={math.sinh(2.0*K):.12f}\n")
        f.write(f"target 1/sqrt(3)={1.0/math.sqrt(3.0):.12f}\n")
        f.write(f"therm={args.therm} meas={args.meas} seed={args.seed}\n\n")
        f.write(f"E_per_site={res['E_per_site']:.8f}\n")
        f.write(f"M_per_site={res['M_per_site']:.8f}\n")
        f.write(f"|M|_per_site={res['|M|_per_site']:.8f}\n")
        f.write(f"C_per_site={res['C_per_site']:.8f}\n")
        f.write(f"chi_per_site={res['chi_per_site']:.8f}\n")

    print(f"Wrote {out}")
    print(f"K={K:.12f}, Kc_tri={(0.25 * math.log(3.0)):.12f}")
    print(f"sinh(2K)={math.sinh(2.0*K):.12f}, target={1.0/math.sqrt(3.0):.12f}")
    print(
        "observables: "
        f"E/site={res['E_per_site']:.6f}, "
        f"|M|/site={res['|M|_per_site']:.6f}, "
        f"C/site={res['C_per_site']:.6f}, "
        f"chi/site={res['chi_per_site']:.6f}"
    )


if __name__ == "__main__":
    main()
