#!/usr/bin/env python3
"""
strategy_test.py — Compare CMA-ES strategies for the 4-5-6 triangle problem.

True answer: r1 ≈ 5.0652, r2 ≈ 7.7429 (β_c ≈ 0.0628).
False minimum: r1 ≈ 1, r2 ≈ 0  (β_c ≈ 0.42, spurious scan peak).

Four strategies tested:
  S1 — failed run baseline:  x0=(3,3), σ0=0.5,   no bounds
  S2 — large sigma only:     x0=(3,4), σ0=3.0,   no bounds   [preset default]
  S3 — large sigma + bounds: x0=(3,4), σ0=3.0,   r2≥0.5
  S4 — closer start:         x0=(5,6), σ0=2.0,   r2≥0.5

Usage:
  python strategy_test.py [--strategies S1 S2 S3 S4] [--max-evals N]
                          [--n-traj N] [--n-workers N] [--no-vis]
                          [--seed N]

Each strategy writes its run to  results/strategy_test/<SN>/.
After all strategies finish a summary table is printed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# ---------------------------------------------------------------------------
# True solution for reference
# ---------------------------------------------------------------------------
R1_TRUE = 5.0652
R2_TRUE = 7.7429
BC_TRUE = 0.06283

# "False minimum" zone: r2 < 0.5 and β_c > 0.2
FALSE_MIN_R2_THRESH  = 0.5
FALSE_MIN_BC_THRESH  = 0.2

# A solution is "correct" if it lands within dist_thresh of (R1_TRUE, R2_TRUE)
# in the normalised sense, OR simply r1 > 3.5 and r2 > 5 (relaxed).
CORRECT_DIST_THRESH = 3.0   # L2 distance in (r1,r2) space


def _dist_to_true(r1, r2):
    return np.sqrt((r1 - R1_TRUE) ** 2 + (r2 - R2_TRUE) ** 2)


def _is_correct(r1, r2):
    return _dist_to_true(r1, r2) < CORRECT_DIST_THRESH


def _is_false_min(r2, beta_c):
    return r2 < FALSE_MIN_R2_THRESH and beta_c > FALSE_MIN_BC_THRESH


# ---------------------------------------------------------------------------
# Strategy definitions
# ---------------------------------------------------------------------------
STRATEGIES = {
    "S1": dict(
        label="Baseline (failed run):  x0=(3,3) σ0=0.5  no bounds",
        x0=(3.0, 3.0), sigma0=0.5,
        lower_bounds=None,
    ),
    "S2": dict(
        label="Large σ only:          x0=(3,4) σ0=3.0  no bounds",
        x0=(3.0, 4.0), sigma0=3.0,
        lower_bounds=None,
    ),
    "S3": dict(
        label="Large σ + bounds:      x0=(3,4) σ0=3.0  r≥(0.5,0.5)",
        x0=(3.0, 4.0), sigma0=3.0,
        lower_bounds=(0.5, 0.5),
    ),
    "S4": dict(
        label="Closer start + bounds: x0=(5,6) σ0=2.0  r≥(0.5,0.5)",
        x0=(5.0, 6.0), sigma0=2.0,
        lower_bounds=(0.5, 0.5),
    ),
}

# ---------------------------------------------------------------------------
# Run one strategy
# ---------------------------------------------------------------------------

def run_strategy(tag: str, cfg: dict, run_dir: str,
                 max_evals: int, n_traj: int,
                 n_traj_scan_coarse: int, n_traj_scan_fine: int,
                 n_workers: int,
                 seed: int | None, no_vis: bool) -> dict:
    """Run one strategy and return a summary dict."""
    import mc_engine
    import cost as cost_module
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # installs fallback patch

    os.makedirs(run_dir, exist_ok=True)

    # Geometry — same as the failed run
    ref_geom  = (13, 16, -3,  3)
    test_geom = (16, 16,  0,  0)
    popsize = 6

    # Reference data (use cached)
    ref_dir = os.path.join(_HERE, "results",
                           "_reference_Lx13_Ly16_Tx-3_Ty3")
    ref_a2a = os.path.join(ref_dir, "two_point_all_to_all.dat")
    if not os.path.exists(ref_a2a):
        raise FileNotFoundError(
            f"Reference cache not found at {ref_a2a}.\n"
            "Run `python run.py --preset 4_5_6` once to build it, or "
            "point ref_dir at the correct path.")
    ref_data = mc_engine.load_all_to_all(ref_a2a)

    # Evaluator
    mc_scratch = os.path.join(run_dir, "mc_scratch")
    exe = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        exe = os.path.join(_HERE, "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        raise FileNotFoundError(f"Simulator not found: {exe}")

    ev_kwargs = dict(
        exe=exe,
        ref_data=ref_data,
        ref_geom=ref_geom,
        test_geom=test_geom,
        output_dir=mc_scratch,
        n_traj_prod=n_traj,
        n_traj_scan_coarse=n_traj_scan_coarse,
        n_traj_scan_fine=n_traj_scan_fine,
        scan_n_coarse=11,
        scan_n_refine=5,
        scan_n_refine2=5,
        scan_n_refine3=5,
        scan_max_shifts=4,
        scan_jackknife=True,
        beta_seed=(0.05, 0.50),
        reweight=False,
    )

    if n_workers > 1:
        from prod_runtime import ProdGenerationPool
        fb_log = os.path.join(run_dir, "fallback_log.jsonl")
        pool = ProdGenerationPool(
            n_workers=n_workers,
            evaluator_kwargs=ev_kwargs,
            fb_kwargs={},
            master_seed=seed,
        )
        # Shim for optimizer (no visualiser in test mode)
        class _Shim:
            n_traj_prod = n_traj
            current_simplex = None
            current_gaussian = None
        ev = _Shim()
        pool_arg = pool
    else:
        ev = Evaluator(**ev_kwargs)
        pool_arg = None

    try:
        t0 = time.perf_counter()
        summary = run_cmaes(
            ev,
            x0=cfg["x0"],
            max_evals=max_evals,
            sigma0=cfg["sigma0"],
            popsize=popsize,
            tolx=0.005,
            tolfun=1e-6,
            snr_floor=0.0,
            indist_stop_snr=0.0,
            snr_target=2.0,
            snr_max_traj_factor=4,
            seed=seed,
            pool=pool_arg,
            lower_bounds=cfg["lower_bounds"],
        )
        wall = time.perf_counter() - t0
    finally:
        if n_workers > 1 and hasattr(pool, "shutdown"):
            pool.shutdown()

    # The evaluator writes to mc_scratch/eval_log.jsonl; copy to run_dir too.
    mc_log = os.path.join(mc_scratch, "eval_log.jsonl")
    log_path = os.path.join(run_dir, "eval_log.jsonl")
    if os.path.exists(mc_log) and mc_log != log_path:
        import shutil
        shutil.copy2(mc_log, log_path)

    # Load history from the jsonl file (works for serial and parallel)
    history = []
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                line = line.strip()
                if line:
                    history.append(json.loads(line))

    # Analyse
    if not history:
        return {"tag": tag, "n_evals": 0, "wall_s": wall,
                "best_r1": None, "best_r2": None,
                "best_cost": None, "best_beta_c": None,
                "verdict": "NO_DATA"}

    def _cost(r):
        return float(r.get("cost", float("inf")))

    best = min(history, key=_cost)
    n_evals = len(history)
    r1    = best.get("r1")
    r2    = best.get("r2")
    bc    = best.get("beta_c")
    cost_ = _cost(best)
    dist  = _dist_to_true(r1, r2) if r1 is not None else float("inf")

    if _is_correct(r1, r2):
        verdict = "CORRECT"
    elif _is_false_min(r2, bc):
        verdict = "FALSE_MIN"
    else:
        verdict = "INCONCLUSIVE"

    result = dict(
        tag=tag, n_evals=n_evals, wall_s=round(wall, 1),
        best_r1=round(r1, 4), best_r2=round(r2, 4),
        best_cost=round(cost_, 6), best_beta_c=round(bc, 5),
        dist_to_true=round(dist, 3),
        verdict=verdict,
    )
    with open(os.path.join(run_dir, "result.json"), "w") as f:
        json.dump(result, f, indent=2)
    print(f"[{tag}] done  verdict={verdict}  best=({r1:.3f},{r2:.3f})  "
          f"cost={cost_:.4f}  β_c={bc:.4f}  dist={dist:.2f}  "
          f"t={wall:.0f}s", flush=True)
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--strategies", nargs="+", default=list(STRATEGIES),
                   choices=list(STRATEGIES), metavar="SN",
                   help="Which strategies to run (default: all).")
    p.add_argument("--max-evals", type=int, default=72,
                   help="Evaluations per strategy (default 72 = 12 gens × 6).")
    p.add_argument("--n-traj", type=int, default=20_000,
                   help="Production MC trajectories per eval (default 20k).")
    p.add_argument("--n-traj-scan-coarse", type=int, default=4000,
                   help="Coarse scan trajectories (default 4k).")
    p.add_argument("--n-traj-scan-fine", type=int, default=10000,
                   help="Fine scan trajectories (default 10k).")
    p.add_argument("--n-workers", type=int, default=6,
                   help="Parallel workers per strategy (default 6).")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for reproducibility (default 42).")
    p.add_argument("--no-vis", action="store_true",
                   help="Suppress dashboard/PNG output.")
    p.add_argument("--results-dir", type=str,
                   default=os.path.join(_HERE, "results", "strategy_test"),
                   help="Parent directory for strategy sub-runs.")
    args = p.parse_args()

    print("=" * 70)
    print("4-5-6 triangle strategy test")
    print(f"  True answer: r1={R1_TRUE}  r2={R2_TRUE}  β_c={BC_TRUE}")
    print(f"  max_evals={args.max_evals}  n_traj={args.n_traj}  "
          f"n_workers={args.n_workers}  seed={args.seed}")
    print("=" * 70)

    results = []
    for tag in args.strategies:
        cfg = STRATEGIES[tag]
        print(f"\n{'─'*70}")
        print(f"[{tag}] {cfg['label']}")
        print(f"      x0={cfg['x0']}  σ0={cfg['sigma0']}  "
              f"bounds={cfg['lower_bounds']}", flush=True)
        run_dir = os.path.join(args.results_dir, tag)
        try:
            res = run_strategy(
                tag=tag, cfg=cfg, run_dir=run_dir,
                max_evals=args.max_evals,
                n_traj=args.n_traj,
                n_traj_scan_coarse=args.n_traj_scan_coarse,
                n_traj_scan_fine=args.n_traj_scan_fine,
                n_workers=args.n_workers,
                seed=args.seed,
                no_vis=args.no_vis,
            )
        except Exception as exc:
            import traceback
            traceback.print_exc()
            res = {"tag": tag, "verdict": f"ERROR: {exc}"}
        results.append(res)

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'tag':<4} {'verdict':<14} {'best r1':>8} {'best r2':>8} "
          f"{'cost':>10} {'β_c':>8} {'dist':>7}  label")
    print("-" * 70)
    for res in results:
        tag      = res.get("tag", "?")
        verdict  = res.get("verdict", "?")
        r1       = res.get("best_r1")
        r2       = res.get("best_r2")
        cost_    = res.get("best_cost")
        bc       = res.get("best_beta_c")
        dist_    = res.get("dist_to_true")
        label    = STRATEGIES.get(tag, {}).get("label", "")
        r1s  = f"{r1:8.3f}" if r1  is not None else "       ?"
        r2s  = f"{r2:8.3f}" if r2  is not None else "       ?"
        cs   = f"{cost_:10.5f}" if cost_ is not None else "         ?"
        bcs  = f"{bc:8.4f}"    if bc   is not None else "       ?"
        ds   = f"{dist_:7.2f}" if dist_ is not None else "      ?"
        print(f"{tag:<4} {verdict:<14} {r1s} {r2s} {cs} {bcs} {ds}  {label}")
    print("=" * 70)

    # Persist summary
    sum_path = os.path.join(args.results_dir, "summary.json")
    os.makedirs(args.results_dir, exist_ok=True)
    with open(sum_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSummary written to {sum_path}")


if __name__ == "__main__":
    main()
