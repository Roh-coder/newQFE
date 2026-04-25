#!/usr/bin/env python3
"""
tools/validate_cache.py — Audit a populated β_c cache (Speedup 2 check 1.D).

Picks N random hits inside the convex hull of the cached samples,
re-runs `mc_engine.find_beta_c` at each, and asserts that the median
absolute deviation between the cached interpolant and a fresh MC β_c
is below `3 × cache.tol_beta`.

Exit codes:
    0   all checks pass (median |Δβ_c| < 3 · tol_beta)
    2   audit failed (median |Δβ_c| ≥ 3 · tol_beta)
    1   not enough samples to audit (cache too sparse)

Usage
-----
    python tools/validate_cache.py --geom 24 24 0 0 --n 10
    python tools/validate_cache.py --geom 24 24 0 0 --n 10 \
        --n-traj-coarse 4000 --n-traj-fine 10000

Notes
-----
* Auditing is intentionally cheap MC (defaults: 4 000 / 10 000 traj)
  because we only care about the *interpolant* error, not absolute
  precision of β_c.  If you want a tighter audit, raise the trajectory
  counts via `--n-traj-coarse` / `--n-traj-fine`.
* The audit set is sampled uniformly from the cache's interior (i.e.
  excluding hull-touching simplices), so the comparison is on the same
  region that real lookups would actually hit.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))   # so we can import the package modules

import mc_engine  # noqa: E402
from betac_cache import BetacCache  # noqa: E402


def _interior_query_points(cache: BetacCache, n: int, rng: np.random.Generator):
    """Pick `n` random points inside the cached interpolant.

    Each point is the centroid of a non-hull simplex jittered by a
    barycentric weight drawn from a Dirichlet — guaranteed to land
    strictly inside the simplex and so be both a hull-interior and a
    non-extrapolating query.  Returns an (n, 2) array of (r1, r2).
    """
    cache._refresh_samples()
    cache._maybe_rebuild()
    if cache._tri is None or cache._hull_simplex_mask is None:
        return None

    interior = np.flatnonzero(~cache._hull_simplex_mask)
    if interior.size == 0:
        return None

    out = []
    for _ in range(n):
        s_idx = int(rng.choice(interior))
        verts_idx = cache._tri.simplices[s_idx]
        verts = cache._u_pts[verts_idx]
        # Dirichlet weights so the point is strictly inside.
        w = rng.dirichlet(np.ones(len(verts)))
        out.append((w[:, None] * verts).sum(axis=0))
    return np.asarray(out, dtype=float)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--geom", nargs=4, type=int, required=True,
                   metavar=("Lx", "Ly", "Tx", "Ty"),
                   help="Test geometry whose β_c cache to audit.")
    p.add_argument("--n", type=int, default=10,
                   help="Number of audit points to draw (default: 10).")
    p.add_argument("--results-root", default="results",
                   help="Root containing the _betac_cache_* directories.")
    p.add_argument("--exe", default="bin/ising_tri_twisted_parallelogram",
                   help="Simulator binary path (default: bin/ising_tri_twisted_parallelogram).")
    p.add_argument("--n-traj-coarse", type=int, default=4000,
                   help="Trajectories per coarse scan point (default: 4000).")
    p.add_argument("--n-traj-fine", type=int, default=10000,
                   help="Trajectories per refinement point (default: 10000).")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for sample selection (default: random).")
    p.add_argument("--scratch", default=None,
                   help="Scratch dir for the audit MC (default: /tmp/validate_cache_<pid>).")
    p.add_argument("--tol-multiplier", type=float, default=3.0,
                   help="Pass threshold: median |Δβ_c| < tol_mult · cache.tol_beta "
                        "(default: 3).")
    args = p.parse_args()

    cache = BetacCache(tuple(args.geom), root=args.results_root, verbose=False)
    n_samples = len(cache)
    print(f"[validate] cache at {cache.dir}")
    print(f"[validate] {n_samples} samples on disk, tol_beta={cache.tol_beta:.2e}, "
          f"tol_r={cache.tol_r}, min_neighbours={cache.min_neighbours}")

    if n_samples < max(cache.min_neighbours + 1, 5):
        print(f"[validate] not enough samples to audit "
              f"(need ≥ {max(cache.min_neighbours + 1, 5)}, have {n_samples}); "
              "preseed the cache or run a few optimization iterations first.",
              file=sys.stderr)
        return 1

    rng = np.random.default_rng(args.seed)
    queries = _interior_query_points(cache, args.n, rng)
    if queries is None or len(queries) == 0:
        print("[validate] interpolant did not build (no interior simplices); "
              "the cache may be co-linear or too sparse.", file=sys.stderr)
        return 1

    scratch = args.scratch or f"/tmp/validate_cache_{os.getpid()}"
    os.makedirs(scratch, exist_ok=True)

    Lx, Ly, Tx, Ty = args.geom
    diffs = []
    started = time.time()
    for i, (r1, r2) in enumerate(queries, start=1):
        cached = cache.lookup(float(r1), float(r2))
        if cached is None:
            # Hit conditions failed; skip this point.
            print(f"[validate] {i:02d}/{len(queries)} (r1={r1:.4f}, r2={r2:.4f}) "
                  "→ cache miss (skipping)")
            continue
        beta_cache, sigma_loo = cached

        # Quick fresh scan around the cached estimate.
        eps = max(0.10 * beta_cache, 0.03)
        bracket = (max(0.01, beta_cache - eps), beta_cache + eps)
        beta_mc, *_ = mc_engine.find_beta_c(
            args.exe, Lx, Ly, Tx, Ty,
            float(r1), float(r2), 1.0,
            bracket[0], bracket[1],
            n_coarse=9, n_refine=4, n_refine2=4, n_refine3=3,
            n_traj_coarse=int(args.n_traj_coarse),
            n_traj_fine=int(args.n_traj_fine),
            data_dir=os.path.join(scratch, f"q{i:02d}"),
            max_shifts=2,
        )
        d = abs(float(beta_mc) - float(beta_cache))
        diffs.append(d)
        print(f"[validate] {i:02d}/{len(queries)} "
              f"r1={r1:.4f} r2={r2:.4f}  "
              f"β_cache={beta_cache:.6f}  β_MC={beta_mc:.6f}  "
              f"|Δ|={d:.2e}  σ_LOO={sigma_loo:.2e}")

    if not diffs:
        print("[validate] no audit points qualified for hits; "
              "loosen tol_beta / tol_r or add more samples.", file=sys.stderr)
        return 1

    diffs = np.asarray(diffs, dtype=float)
    median = float(np.median(diffs))
    p95 = float(np.quantile(diffs, 0.95))
    threshold = args.tol_multiplier * cache.tol_beta
    elapsed = time.time() - started

    print()
    print(f"[validate] {len(diffs)} hits audited in {elapsed:.1f}s")
    print(f"[validate]   median |Δβ_c| = {median:.2e}")
    print(f"[validate]   p95    |Δβ_c| = {p95:.2e}")
    print(f"[validate]   threshold     = {threshold:.2e}  "
          f"({args.tol_multiplier}× tol_beta)")

    if median < threshold:
        print("[validate] PASS")
        return 0
    print("[validate] FAIL — interpolant is too coarse for current tol_beta. "
          "Tighten tol_beta or add more cache samples.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
