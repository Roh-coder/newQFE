#!/usr/bin/env python3
"""
preseed_betac.py — Bootstrap the persistent β_c cache with an initial grid.

Useful before launching an optimizer run so that even the very first
evaluation hits the cache instead of paying for a 3-pass scan.  Skips
points already present (within tol_r) so it is safe to re-run.

Usage
-----
    python tools/preseed_betac.py --geom 24 24 0 0 \\
        --r1-range 0.5 1.5 --r2-range 0.5 1.5 --grid 5 \\
        --beta-range 0.20 0.32

Run from the K_from_Optimization_pnp/ directory (or pass --root explicitly
so paths resolve relative to your results/ tree).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG  = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import mc_engine  # noqa: E402
from betac_cache import BetacCache  # noqa: E402


def _already_close(cache, r1, r2, tol):
    """Skip if any cached sample is within tol_r in either coord."""
    cache._refresh_samples()
    for s in cache._samples:
        if abs(s["r1"] - r1) < tol and abs(s["r2"] - r2) < tol:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geom", nargs=4, type=int, required=True,
                    metavar=("Lx", "Ly", "Tx", "Ty"))
    ap.add_argument("--r1-range", nargs=2, type=float, default=(0.5, 1.5))
    ap.add_argument("--r2-range", nargs=2, type=float, default=(0.5, 1.5))
    ap.add_argument("--grid", type=int, default=5,
                    help="N grid points per axis (default 5 → 25 evals)")
    ap.add_argument("--beta-range", nargs=2, type=float, default=(0.20, 0.32))
    ap.add_argument("--n-traj-coarse", type=int, default=4000)
    ap.add_argument("--n-traj-fine", type=int, default=10000)
    ap.add_argument("--root", default=os.path.join(PKG, "results"),
                    help="results/ directory (cache lives under here)")
    ap.add_argument("--exe", default=os.path.join(PKG,
                                                  "bin/ising_tri_twisted_parallelogram"))
    ap.add_argument("--tol-r", type=float, default=0.02,
                    help="skip points within this distance of an existing one")
    args = ap.parse_args()

    cache = BetacCache(tuple(args.geom), root=args.root, verbose=True)
    Lx, Ly, Tx, Ty = args.geom
    print(f"[preseed] geom=({Lx},{Ly},{Tx},{Ty})  cache has "
          f"{len(cache)} samples already")

    grid_r1 = np.linspace(*args.r1_range, args.grid)
    grid_r2 = np.linspace(*args.r2_range, args.grid)
    todo = []
    for r1 in grid_r1:
        for r2 in grid_r2:
            if not _already_close(cache, float(r1), float(r2), args.tol_r):
                todo.append((float(r1), float(r2)))
    print(f"[preseed] {len(todo)} new (r1, r2) points to evaluate")

    for k, (r1, r2) in enumerate(todo, start=1):
        print(f"[preseed {k}/{len(todo)}] r1={r1:.4f} r2={r2:.4f}")
        beta_c, beta_c_sigma, _, sb, _, _ = mc_engine.find_beta_c(
            args.exe, Lx, Ly, Tx, Ty, r1, r2, 1.0,
            args.beta_range[0], args.beta_range[1],
            n_traj_coarse=args.n_traj_coarse,
            n_traj_fine=args.n_traj_fine,
            data_dir=os.path.join(args.root, "_preseed_scratch"),
        )
        cache.add(r1, r2, beta_c, beta_c_sigma, source_run="preseed",
                  n_traj_total=len(sb) * (args.n_traj_coarse +
                                          args.n_traj_fine) // 2)
    print(f"[preseed] done; cache now has {len(cache)} samples")


if __name__ == "__main__":
    main()
