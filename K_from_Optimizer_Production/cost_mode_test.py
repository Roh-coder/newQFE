#!/usr/bin/env python3
"""
cost_mode_test.py — Workflow-level comparison of cost functions.

Runs the SAME CMA-ES strategy on the 4-5-6 triangle problem twice:

  A — cost_mode="l4mean_both_interp"  (legacy default;
       both ref & test interpolated, per-direction L²,
       directions aggregated as p=4 power-mean)
  B — cost_mode="test_native"  (ref interpolated; test sampled at
       the integer lattice sites along each test boundary direction;
       per-direction Lᵖ residual mean, directions plain-mean)
       — power 2 by default; ``--cost-power 4`` to also probe L⁴.

True optimum (geometry-determined, independent of cost mode):
    r1 ≈ 5.0652   r2 ≈ 7.7429   β_c ≈ 0.0628

Both runs share the same seed, x0, σ0, budget, and reference cache, so
any difference in the converged best is attributable to the cost.

Usage:
  python cost_mode_test.py                          # default: 30 evals each
  python cost_mode_test.py --max-evals 18           # cheaper smoke test
  python cost_mode_test.py --cost-power 4           # B uses L⁴ residuals
  python cost_mode_test.py --modes A B              # subset
  python cost_mode_test.py --n-workers 1 --no-vis   # serial / quiet
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# True solution for the 4-5-6 triangle reference (test = 16×16 untwisted).
R1_TRUE = 5.0652
R2_TRUE = 7.7429
BC_TRUE = 0.06283


def _dist(r1, r2):
    return float(np.hypot(r1 - R1_TRUE, r2 - R2_TRUE))


# ---------------------------------------------------------------------------
def run_one(tag: str, cost_mode: str, cost_power: int,
            run_dir: str, max_evals: int, n_traj: int,
            n_workers: int, seed: int) -> dict:
    """One CMA-ES run with the requested cost; returns a summary dict."""
    import mc_engine
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # installs reweight fallback patch

    os.makedirs(run_dir, exist_ok=True)

    ref_geom  = (13, 16, -3,  3)
    test_geom = (16, 16,  0,  0)
    popsize   = 6

    # Reuse the strategy_test reference cache to avoid a long ref MC.
    ref_dir = os.path.join(_HERE, "results",
                           "_reference_Lx13_Ly16_Tx-3_Ty3")
    ref_a2a = os.path.join(ref_dir, "two_point_all_to_all.dat")
    if not os.path.exists(ref_a2a):
        raise FileNotFoundError(
            f"Reference cache missing: {ref_a2a}\n"
            "Build it first (e.g. run strategy_test.py once or "
            "`python run.py` with the 13×16 −3,3 reference).")
    ref_data = mc_engine.load_all_to_all(ref_a2a)

    exe = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        raise FileNotFoundError(f"Simulator not found: {exe}")

    mc_scratch = os.path.join(run_dir, "mc_scratch")
    ev_kwargs = dict(
        exe=exe,
        ref_data=ref_data,
        ref_geom=ref_geom,
        test_geom=test_geom,
        output_dir=mc_scratch,
        n_traj_prod=n_traj,
        n_traj_scan_coarse=4000,
        n_traj_scan_fine=10000,
        scan_n_coarse=11,
        scan_n_refine=5,
        scan_n_refine2=5,
        scan_n_refine3=5,
        scan_max_shifts=4,
        scan_jackknife=True,
        beta_seed=(0.05, 0.50),
        reweight=False,
        cost_mode=cost_mode,
        cost_power=cost_power,
    )

    if n_workers > 1:
        from prod_runtime import ProdGenerationPool
        pool = ProdGenerationPool(
            n_workers=n_workers,
            evaluator_kwargs=ev_kwargs,
            fb_kwargs={},
            master_seed=seed,
        )
        class _Shim:
            n_traj_prod = n_traj
            current_simplex = None
            current_gaussian = None
        ev = _Shim()
        pool_arg = pool
    else:
        ev = Evaluator(**ev_kwargs)
        pool_arg = None

    t0 = time.perf_counter()
    try:
        run_cmaes(
            ev,
            x0=(5.0, 6.0),
            max_evals=max_evals,
            sigma0=2.0,
            popsize=popsize,
            tolx=0.005,
            tolfun=1e-6,
            snr_floor=0.0,
            indist_stop_snr=0.0,
            snr_target=2.0,
            snr_max_traj_factor=4,
            seed=seed,
            pool=pool_arg,
            lower_bounds=(0.5, 0.5),
        )
    finally:
        if n_workers > 1 and hasattr(pool, "shutdown"):
            pool.shutdown()
    wall = time.perf_counter() - t0

    # Read evaluation history written by the (workers') Evaluators.
    log_path = os.path.join(mc_scratch, "eval_log.jsonl")
    history: list[dict] = []
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                line = line.strip()
                if line:
                    history.append(json.loads(line))

    if not history:
        return {"tag": tag, "cost_mode": cost_mode, "cost_power": cost_power,
                "n_evals": 0, "wall_s": round(wall, 1),
                "verdict": "NO_DATA"}

    def _c(r): return float(r.get("cost", float("inf")))
    best = min(history, key=_c)

    # First eval index that lands within 1.0 of true (r1,r2).
    first_close = None
    for i, r in enumerate(history, start=1):
        if r.get("r1") is None or r.get("r2") is None:
            continue
        if _dist(r["r1"], r["r2"]) < 1.0:
            first_close = i
            break

    res = dict(
        tag=tag,
        cost_mode=cost_mode,
        cost_power=cost_power,
        n_evals=len(history),
        wall_s=round(wall, 1),
        best_r1=round(best["r1"], 4),
        best_r2=round(best["r2"], 4),
        best_cost=float(_c(best)),
        best_beta_c=round(best["beta_c"], 5),
        dist_to_true=round(_dist(best["r1"], best["r2"]), 3),
        first_close_eval=first_close,
    )
    with open(os.path.join(run_dir, "result.json"), "w") as f:
        json.dump(res, f, indent=2)

    fc = res["first_close_eval"]
    fc_str = str(fc) if fc is not None else ">budget"
    print(f"[{tag}] {cost_mode}(p={cost_power})  "
          f"best=({res['best_r1']:.3f},{res['best_r2']:.3f})  "
          f"cost={res['best_cost']:.4e}  β_c={res['best_beta_c']:.4f}  "
          f"dist={res['dist_to_true']:.2f}  "
          f"first<1.0 at eval={fc_str}  t={wall:.0f}s",
          flush=True)
    return res


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--modes", nargs="+", default=["A", "B", "C", "D"],
                   choices=["A", "B", "C", "D"],
                   help="Which runs to perform (A=l4mean B=test_native C=affine_fit D=huber_log).")
    p.add_argument("--max-evals", type=int, default=30,
                   help="CMA-ES evaluations per run (default 30 = 5 gens × 6).")
    p.add_argument("--n-traj", type=int, default=10_000,
                   help="Production MC trajectories per eval (default 10k).")
    p.add_argument("--n-workers", type=int, default=6,
                   help="Parallel workers (default 6).")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--cost-power", type=int, default=2, choices=(2, 4),
                   help="Residual power for mode B / test_native (default 2).")
    p.add_argument("--results-dir", type=str,
                   default=os.path.join(_HERE, "results", "cost_mode_test"))
    args = p.parse_args()

    print("=" * 72)
    print("Cost-function workflow comparison (4-5-6 triangle problem)")
    print(f"  ref=(13,16,-3,3)  test=(16,16,0,0)  "
          f"true r=(5.0652, 7.7429)  β_c≈0.0628")
    print(f"  max_evals={args.max_evals}  n_traj={args.n_traj}  "
          f"n_workers={args.n_workers}  seed={args.seed}")
    print("=" * 72)

    runs = {
        "A": dict(cost_mode="l4mean_both_interp", cost_power=2),
        "B": dict(cost_mode="test_native",        cost_power=args.cost_power),
        "C": dict(cost_mode="affine_fit",         cost_power=2),
        "D": dict(cost_mode="huber_log",          cost_power=2),
    }

    results = []
    for tag in args.modes:
        cfg = runs[tag]
        run_dir = os.path.join(args.results_dir, tag)
        print(f"\n[--- {tag}: cost_mode={cfg['cost_mode']!r} "
              f"cost_power={cfg['cost_power']} ---]", flush=True)
        try:
            res = run_one(
                tag=tag,
                cost_mode=cfg["cost_mode"],
                cost_power=cfg["cost_power"],
                run_dir=run_dir,
                max_evals=args.max_evals,
                n_traj=args.n_traj,
                n_workers=args.n_workers,
                seed=args.seed,
            )
        except Exception as exc:
            import traceback
            traceback.print_exc()
            res = {"tag": tag, "cost_mode": cfg["cost_mode"],
                   "cost_power": cfg["cost_power"],
                   "verdict": f"ERROR: {exc}"}
        results.append(res)

    # Summary
    print(f"\n{'=' * 72}")
    print("SUMMARY")
    print(f"{'tag':<3} {'cost_mode':<22} {'p':>2} {'best r1':>8} "
          f"{'best r2':>8} {'best cost':>11} {'β_c':>8} "
          f"{'dist':>6} {'first<1.0':>10} {'wall_s':>7}")
    print("-" * 102)
    for r in results:
        if "best_r1" not in r:
            print(f"{r.get('tag','?'):<3} {r.get('cost_mode','?'):<22} "
                  f"{r.get('cost_power','?'):>2}  {r.get('verdict','?')}")
            continue
        fc = r["first_close_eval"]
        fc_str = str(fc) if fc is not None else ">budget"
        print(f"{r['tag']:<3} {r['cost_mode']:<22} {r['cost_power']:>2} "
              f"{r['best_r1']:>8.3f} {r['best_r2']:>8.3f} "
              f"{r['best_cost']:>11.4e} {r['best_beta_c']:>8.4f} "
              f"{r['dist_to_true']:>6.2f} {fc_str:>10} "
              f"{r['wall_s']:>7.0f}")
    print("=" * 72)
    print(f"True target: r1={R1_TRUE}  r2={R2_TRUE}  β_c={BC_TRUE}")

    # Persist
    os.makedirs(args.results_dir, exist_ok=True)
    sum_path = os.path.join(args.results_dir, "summary.json")
    with open(sum_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSummary written to {sum_path}")


if __name__ == "__main__":
    main()
