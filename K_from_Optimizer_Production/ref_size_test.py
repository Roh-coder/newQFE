#!/usr/bin/env python3
"""
ref_size_test.py — Hypothesis test: does the small reference lattice cause
the cost-surface false minimum?

Three reference geometries with the same 4-5-6 side ratio at three sizes
are tested against the same equilateral test lattice (16, 16, 0, 0):

  small  (1×):  (Lx=13, Ly=16, Tx=-3, Ty=3)  → sides ≈ 11.79, 14.73, 17.69
  medium (2×):  (Lx=26, Ly=32, Tx=-6, Ty=6)  → sides ≈ 23.58, 29.46, 35.38
  large  (3×):  (Lx=39, Ly=48, Tx=-9, Ty=9)  → sides ≈ 35.37, 44.19, 53.07

For each reference we:
  1) build (or load cached) reference correlator at that ref's β_c
  2) evaluate the cost on a uniform grid in (r1, r2) ∈ [r1_min, r1_max]
                                                    × [r2_min, r2_max]
  3) optionally also run a CMA-ES with bounds and large σ₀

The TRUE optimum (in (r1, r2)) is determined by the test geometry alone
and is INDEPENDENT of the reference scale: r1≈5.07, r2≈7.74, β_c≈0.0628.

Hypothesis: the grid argmin should approach the true optimum monotonically
as the reference lattice grows.

Usage:
  python ref_size_test.py [--sizes 1 2 3] [--grid-n 9]
                          [--r1-range 1 9] [--r2-range 1 9]
                          [--n-traj 5000] [--ref-n-traj 20000]
                          [--n-workers 6] [--cma-evals 36]
                          [--no-cma] [--quick]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# ---------------------------------------------------------------------------
# True optimum (test=(16,16,0,0), 4-5-6 ref geometry)
# ---------------------------------------------------------------------------
R1_TRUE = 5.0652
R2_TRUE = 7.7429
BC_TRUE = 0.06283

# ---------------------------------------------------------------------------
# Reference geometries
# ---------------------------------------------------------------------------
REF_SIZES = {
    1: dict(geom=(13, 16, -3,  3), label="small (1×)"),
    2: dict(geom=(26, 32, -6,  6), label="medium (2×)"),
    3: dict(geom=(39, 48, -9,  9), label="large (3×)"),
}


def _side_lengths(Lx, Ly, Tx, Ty):
    s32 = math.sqrt(3) / 2
    def L(m, n):
        return math.hypot(m + 0.5 * n, s32 * n)
    return L(Lx, Ty), L(Tx, -Ly), L(-Lx - Tx, Ly - Ty)


# ---------------------------------------------------------------------------
# Reference build (delegates to run.build_reference)
# ---------------------------------------------------------------------------

def ensure_reference(ref_size: int, ref_n_traj: int, results_root: str,
                     beta_seed=(0.05, 0.50)):
    from run import build_reference, ensure_binary

    geom = REF_SIZES[ref_size]["geom"]
    Lx, Ly, Tx, Ty = geom

    cfg = {
        "ref_Lx": Lx, "ref_Ly": Ly, "ref_Tx": Tx, "ref_Ty": Ty,
        "ref_n_traj": ref_n_traj,
        "ref_scan_n_traj_coarse": 4000,
        "ref_scan_n_traj_fine":   10000,
        "ref_beta_seed": tuple(beta_seed),
        "ref_scan_n_coarse": 11,
        "ref_scan_n_refine":  5,
        "ref_scan_n_refine2": 5,
        "ref_scan_n_refine3": 5,
        "ref_scan_max_shifts": 4,
        "scan_jackknife": True,
        "exe": os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram"),
    }
    ensure_binary(cfg["exe"])

    ref_dir = os.path.join(results_root,
                           f"_reference_Lx{Lx}_Ly{Ly}_Tx{Tx}_Ty{Ty}")
    ref_data, meta = build_reference(cfg, ref_dir)
    return ref_data, meta, ref_dir


# ---------------------------------------------------------------------------
# Grid evaluation
# ---------------------------------------------------------------------------

def grid_points(grid_n, r1_range, r2_range):
    r1s = np.linspace(r1_range[0], r1_range[1], grid_n)
    r2s = np.linspace(r2_range[0], r2_range[1], grid_n)
    pts = [(float(r1), float(r2)) for r2 in r2s for r1 in r1s]
    return pts, r1s, r2s


def run_grid(pool, points, eval_id_base=1):
    """Dispatch all grid points as one 'generation' to the pool."""
    print(f"[grid] dispatching {len(points)} points to "
          f"{pool.n_workers} workers …", flush=True)
    t0 = time.perf_counter()
    results = pool.map_generation(points, eval_id_base=eval_id_base)
    wall = time.perf_counter() - t0
    print(f"[grid] done  {len(results)} evals  t={wall:.0f}s "
          f"({wall/len(results):.1f}s/eval)", flush=True)
    return results


# ---------------------------------------------------------------------------
# CMA-ES run
# ---------------------------------------------------------------------------

def run_cma_for_ref(pool, x0, sigma0, max_evals, lower_bounds, seed):
    from optimizer import run_cmaes

    class _Shim:
        n_traj_prod = 0   # set below from pool kwargs
        current_simplex = None
        current_gaussian = None
    shim = _Shim()
    shim.n_traj_prod = pool._eval_kwargs.get("n_traj_prod", 0)

    summary = run_cmaes(
        shim,
        x0=tuple(x0),
        max_evals=max_evals,
        sigma0=sigma0,
        popsize=6,
        tolx=0.005,
        tolfun=1e-6,
        snr_floor=0.0,
        indist_stop_snr=0.0,
        snr_target=2.0,
        snr_max_traj_factor=4,
        seed=seed,
        pool=pool,
        lower_bounds=lower_bounds,
    )
    return summary


# ---------------------------------------------------------------------------
# Per-size driver
# ---------------------------------------------------------------------------

def run_one_size(ref_size: int, args, out_root: str):
    label = REF_SIZES[ref_size]["label"]
    geom  = REF_SIZES[ref_size]["geom"]
    sides = _side_lengths(*geom)

    print(f"\n{'='*72}")
    print(f"REF SIZE {ref_size}× — {label}  geom={geom}")
    print(f"  sides v={sides[0]:.2f}  u={sides[1]:.2f}  w={sides[2]:.2f}  "
          f"(min={min(sides):.2f})")
    print('='*72, flush=True)

    out_dir = os.path.join(out_root, f"size{ref_size}")
    os.makedirs(out_dir, exist_ok=True)

    # 1) Reference
    ref_data, ref_meta, ref_dir = ensure_reference(
        ref_size, args.ref_n_traj, out_root)
    print(f"  ref β_c = {ref_meta['beta_c']:.6f}  "
          f"n_traj={ref_meta['n_traj']}", flush=True)

    # 2) Pool
    from prod_runtime import ProdGenerationPool

    eval_kwargs = dict(
        exe=os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram"),
        ref_data=ref_data,
        ref_geom=geom,
        test_geom=(16, 16, 0, 0),
        output_dir=os.path.join(out_dir, "mc_grid"),
        n_traj_prod=args.n_traj,
        n_traj_scan_coarse=2000,
        n_traj_scan_fine=4000,
        scan_n_coarse=11,
        scan_n_refine=5,
        scan_n_refine2=5,
        scan_n_refine3=5,
        scan_max_shifts=4,
        scan_jackknife=True,
        beta_seed=(0.05, 0.50),
        reweight=False,
    )

    points, r1s, r2s = grid_points(
        args.grid_n, args.r1_range, args.r2_range)

    pool = ProdGenerationPool(
        n_workers=args.n_workers,
        evaluator_kwargs=eval_kwargs,
        fb_kwargs={},
        master_seed=args.seed,
    )

    grid_results = []
    cma_summary  = None
    cma_results  = []
    try:
        # ----- GRID -----
        grid_results = run_grid(pool, points, eval_id_base=1)

        # Save grid
        cost_grid = np.full((len(r2s), len(r1s)), np.nan)
        sigma_grid = np.full_like(cost_grid, np.nan)
        bc_grid    = np.full_like(cost_grid, np.nan)
        per_dir_grid = np.full((len(r2s), len(r1s), 3), np.nan)
        for k, res in enumerate(grid_results):
            i = k // len(r1s)
            j = k % len(r1s)
            cost_grid[i, j]  = res["cost"]
            sigma_grid[i, j] = res["sigma_cost"]
            bc_grid[i, j]    = res["beta_c"]
            per_dir_grid[i, j, :] = res["per_dir"]
        np.savez(os.path.join(out_dir, "grid.npz"),
                 r1s=r1s, r2s=r2s,
                 cost=cost_grid, sigma=sigma_grid,
                 beta_c=bc_grid, per_dir=per_dir_grid)

        with open(os.path.join(out_dir, "grid_log.jsonl"), "w") as f:
            for res in grid_results:
                f.write(json.dumps(res) + "\n")

        # Grid argmin
        i_min, j_min = np.unravel_index(np.nanargmin(cost_grid), cost_grid.shape)
        r1_min, r2_min = float(r1s[j_min]), float(r2s[i_min])
        c_min = float(cost_grid[i_min, j_min])
        print(f"\n[grid] argmin: r1={r1_min:.3f}  r2={r2_min:.3f}  "
              f"cost={c_min:.5f}  β_c={bc_grid[i_min,j_min]:.4f}", flush=True)

        # Cost at the NEAREST grid point to the true optimum
        i_true = int(np.argmin(np.abs(r2s - R2_TRUE)))
        j_true = int(np.argmin(np.abs(r1s - R1_TRUE)))
        c_true = float(cost_grid[i_true, j_true])
        pd_true = per_dir_grid[i_true, j_true]
        print(f"[grid] near-true: r1={r1s[j_true]:.3f}  r2={r2s[i_true]:.3f}  "
              f"cost={c_true:.5f}  per_dir={pd_true.tolist()}", flush=True)

        # ----- CMA-ES -----
        if not args.no_cma:
            # New evaluator dir for CMA so logs are separate
            eval_kwargs["output_dir"] = os.path.join(out_dir, "mc_cma")
            os.makedirs(eval_kwargs["output_dir"], exist_ok=True)
            # Need to recreate pool to point at new output_dir
            pool.shutdown()
            pool = ProdGenerationPool(
                n_workers=args.n_workers,
                evaluator_kwargs=eval_kwargs,
                fb_kwargs={},
                master_seed=args.seed,
            )
            print(f"\n[cma] x0=(3,4) σ0=3 bounds=(0.5,0.5) "
                  f"max_evals={args.cma_evals}", flush=True)
            cma_summary = run_cma_for_ref(
                pool, x0=(3.0, 4.0), sigma0=3.0,
                max_evals=args.cma_evals,
                lower_bounds=(0.5, 0.5),
                seed=args.seed,
            )
            # Read back the eval log
            cma_log = os.path.join(eval_kwargs["output_dir"], "eval_log.jsonl")
            if os.path.exists(cma_log):
                with open(cma_log) as f:
                    cma_results = [json.loads(l) for l in f if l.strip()]
                cma_best = min(cma_results, key=lambda e: e["cost"])
                print(f"[cma] best:  r1={cma_best['r1']:.3f}  "
                      f"r2={cma_best['r2']:.3f}  cost={cma_best['cost']:.5f}  "
                      f"β_c={cma_best['beta_c']:.4f}", flush=True)

    finally:
        pool.shutdown()

    # Per-size summary
    summary = dict(
        ref_size=ref_size,
        ref_geom=list(geom),
        ref_sides=[round(s, 3) for s in sides],
        ref_beta_c=ref_meta["beta_c"],
        ref_n_traj=ref_meta["n_traj"],
        n_traj_eval=args.n_traj,
        grid_argmin=dict(r1=r1_min, r2=r2_min, cost=c_min),
        near_true=dict(r1=float(r1s[j_true]), r2=float(r2s[i_true]),
                       cost=c_true, per_dir=pd_true.tolist()),
        cma_best=(dict(r1=cma_best["r1"], r2=cma_best["r2"],
                       cost=cma_best["cost"],
                       beta_c=cma_best["beta_c"])
                  if cma_results else None),
        n_grid=int(args.grid_n) ** 2,
        n_cma=len(cma_results) if cma_results else 0,
    )
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sizes", nargs="+", type=int, default=[1, 2, 3],
                   choices=[1, 2, 3])
    p.add_argument("--grid-n", type=int, default=9,
                   help="Grid resolution per axis (default 9 → 81 pts)")
    p.add_argument("--r1-range", nargs=2, type=float, default=[1.0, 9.0])
    p.add_argument("--r2-range", nargs=2, type=float, default=[1.0, 9.0])
    p.add_argument("--n-traj", type=int, default=5000,
                   help="Production MC traj per grid point (default 5k)")
    p.add_argument("--ref-n-traj", type=int, default=20000,
                   help="Reference MC traj for each ref size (default 20k)")
    p.add_argument("--n-workers", type=int, default=6)
    p.add_argument("--cma-evals", type=int, default=36)
    p.add_argument("--no-cma", action="store_true",
                   help="Skip CMA-ES, grid only")
    p.add_argument("--quick", action="store_true",
                   help="Tiny test: sizes=[1,2], grid_n=5, n_traj=2000, "
                        "ref_n_traj=5000, no CMA")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--results-dir", default=os.path.join(
        _HERE, "results", "ref_size_test"))
    args = p.parse_args()

    if args.quick:
        args.sizes      = [1, 2]
        args.grid_n     = 5
        args.n_traj     = 2000
        args.ref_n_traj = 5000
        args.no_cma     = True

    args.r1_range = tuple(args.r1_range)
    args.r2_range = tuple(args.r2_range)
    os.makedirs(args.results_dir, exist_ok=True)

    print("=" * 72)
    print("Reference-size hypothesis test")
    print(f"  test geom: (16, 16, 0, 0)   (equilateral, side=16)")
    print(f"  TRUE optimum (test): r1={R1_TRUE} r2={R2_TRUE} β_c≈{BC_TRUE}")
    print(f"  sizes: {args.sizes}  grid: {args.grid_n}×{args.grid_n}")
    print(f"  r1∈[{args.r1_range[0]},{args.r1_range[1]}]  "
          f"r2∈[{args.r2_range[0]},{args.r2_range[1]}]")
    print(f"  n_traj_eval={args.n_traj}  ref_n_traj={args.ref_n_traj}  "
          f"workers={args.n_workers}  cma_evals="
          f"{0 if args.no_cma else args.cma_evals}")
    print("=" * 72, flush=True)

    all_summaries = []
    for ref_size in args.sizes:
        try:
            s = run_one_size(ref_size, args, args.results_dir)
            all_summaries.append(s)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            all_summaries.append(dict(ref_size=ref_size, error=str(exc)))

    # Final table
    print(f"\n{'='*72}\nVERDICT\n{'='*72}")
    print(f"{'size':<6} {'sides (min)':<10} "
          f"{'grid argmin':<22} {'near-true cost':<16} "
          f"{'cma best':<22}")
    print("-" * 72)
    for s in all_summaries:
        if "error" in s:
            print(f"{s['ref_size']:<6} ERROR: {s['error']}")
            continue
        ga = s["grid_argmin"]
        cb = s.get("cma_best")
        ga_s = f"({ga['r1']:.2f},{ga['r2']:.2f}) c={ga['cost']:.4f}"
        cb_s = (f"({cb['r1']:.2f},{cb['r2']:.2f}) c={cb['cost']:.4f}"
                if cb else "—")
        side_min = min(s["ref_sides"])
        nt_c = s["near_true"]["cost"]
        print(f"{s['ref_size']:<6} {side_min:<10.2f} "
              f"{ga_s:<22} {nt_c:<16.5f} {cb_s:<22}")
    print(f"\n  TRUE: ({R1_TRUE},{R2_TRUE})")
    print("=" * 72)
    with open(os.path.join(args.results_dir, "all_summary.json"), "w") as f:
        json.dump(all_summaries, f, indent=2)
    print(f"\nFull results in {args.results_dir}")


if __name__ == "__main__":
    main()
