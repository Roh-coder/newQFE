#!/usr/bin/env python3
"""
realmc_cmaes_patch_test.py — REAL-MC CMA-ES test of the T3 patch.

No modification to existing code: we run one CMA-ES pass on a contrived
test geometry whose middle boundary direction collapses (gcd=1), saving
the per-direction cost contributions to eval_log.jsonl, then shadow-
recompute the proposed (T3-patched) cost from per_dir/per_dir_sigma
already in the log.

Setup
-----
  ref  = (12, 12, 0, 0)  — cached reference at exact β_c=ln(3)/4
  test = (12, 11, 0, 0)  — boundary path gcds = [12, 11, 1] → 1 collapsed
  CMA-ES: x0=(1.0, 1.0), σ0=0.4, popsize=6, max_evals=18 (3 gens)
  n_traj=4000  (cheap; fixed-geometry comparison only)

Once the run completes, for every entry in eval_log.jsonl:
  - current cost     = mean(per_dir)        (with collapsed dir = 0)
  - proposed cost    = mean(per_dir over valid)   (drop collapsed)
  - check selection invariance: rank within each gen of popsize entries
    must be IDENTICAL between current and proposed

Outputs
-------
  results/_native_triage/realmc_cmaes/  (run dir, mc_scratch + eval log)
  results/_native_triage/realmc_cmaes_summary.json
  results/_native_triage/realmc_cmaes.png
  results/_native_triage/realmc_cmaes.log
"""
from __future__ import annotations

import io
import json
import os
import sys
import time
from typing import List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)
RUN_DIR = os.path.join(OUT, "realmc_cmaes")
os.makedirs(RUN_DIR, exist_ok=True)


def _proposed_from_perdir(per_dir, per_sig):
    """Drop dirs where (C, σ) == (0, 0) ⇒ collapsed, then plain mean."""
    valid = [(c, s) for c, s in zip(per_dir, per_sig)
             if not (c == 0.0 and s == 0.0)]
    if not valid:
        return float("nan"), float("nan"), 0
    n = len(valid)
    cost = sum(c for c, _ in valid) / n
    sig = (sum(s * s for _, s in valid) ** 0.5) / n
    return cost, sig, n


def main():
    log_buf = io.StringIO()

    class Tee:
        def __init__(self, *s): self.s = s
        def write(self, x):
            for s in self.s: s.write(x)
        def flush(self):
            for s in self.s: s.flush()
    sys.stdout = Tee(sys.__stdout__, log_buf)

    print("=" * 76)
    print("realmc_cmaes_patch_test — REAL MC CMA-ES on collapsed-dir geometry")
    print("=" * 76)

    import mc_engine
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # installs reweight fallback patch
    import cost as cost_mod

    ref_geom = (12, 12, 0, 0)
    test_geom = (12, 11, 0, 0)
    paths = cost_mod.boundary_paths(*test_geom)
    import math
    gcds = [math.gcd(abs(a), abs(b)) for a, b in paths]
    print(f"  ref={ref_geom}  test={test_geom}")
    print(f"  test boundary paths={paths}  gcds={gcds}  "
          f"(N_valid expected = {sum(1 for g in gcds if g >= 2)})")

    ref_dir = os.path.join(HERE, "results", "_reference_Lx12_Ly12_Tx0_Ty0")
    ref_a2a = os.path.join(ref_dir, "two_point_all_to_all.dat")
    if not os.path.exists(ref_a2a):
        raise FileNotFoundError(ref_a2a)
    ref_data = mc_engine.load_all_to_all(ref_a2a)

    exe = os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        raise FileNotFoundError(exe)

    mc_scratch = os.path.join(RUN_DIR, "mc_scratch")
    ev_kwargs = dict(
        exe=exe, ref_data=ref_data,
        ref_geom=ref_geom, test_geom=test_geom,
        output_dir=mc_scratch,
        n_traj_prod=4000,
        n_traj_scan_coarse=2000,
        n_traj_scan_fine=4000,
        scan_n_coarse=9,
        scan_n_refine=5, scan_n_refine2=5, scan_n_refine3=5,
        scan_max_shifts=4, scan_jackknife=True,
        beta_seed=(0.20, 0.32),     # equilateral region
        reweight=False,
        cost_mode="test_native", cost_power=2,
    )
    ev = Evaluator(**ev_kwargs)

    t0 = time.perf_counter()
    run_cmaes(
        ev,
        x0=(1.0, 1.0),
        max_evals=18,
        sigma0=0.4,
        popsize=6,
        tolx=0.005, tolfun=1e-6,
        snr_floor=0.0, indist_stop_snr=0.0, snr_target=2.0,
        snr_max_traj_factor=4,
        seed=2026,
        pool=None,
        lower_bounds=(0.5, 0.5),
    )
    wall = time.perf_counter() - t0
    print(f"\nCMA-ES wall time: {wall:.0f}s")

    # ---- Replay eval_log.jsonl with proposed cost ----
    log_path = os.path.join(mc_scratch, "eval_log.jsonl")
    rows = []
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    print(f"  eval_log entries: {len(rows)}")

    POP = 6
    summary = {
        "n_evals": len(rows),
        "wall_s": round(wall, 1),
        "ref_geom": list(ref_geom), "test_geom": list(test_geom),
        "boundary_path_gcds": gcds,
        "per_eval": [],
        "selection_flips_per_gen": [],
    }

    print("\n  per-eval current vs proposed (T3) cost:")
    print(f"  {'eid':>3} {'r1':>6} {'r2':>6} {'cur':>10} {'prop':>10} "
          f"{'ratio':>6} {'N_v':>3}")
    for r in rows:
        cur = r["cost"]
        prop, sig_p, nv = _proposed_from_perdir(r["per_dir"], r["per_dir_sigma"])
        ratio = cur / prop if prop and prop != 0 else float("nan")
        summary["per_eval"].append(dict(
            eid=r["eval_id"], r1=r["r1"], r2=r["r2"],
            cost_current=cur, sigma_current=r["sigma_cost"],
            cost_proposed=prop, sigma_proposed=sig_p,
            n_valid=nv, ratio_cur_over_prop=ratio,
        ))
        print(f"  {r['eval_id']:>3} {r['r1']:>6.3f} {r['r2']:>6.3f} "
              f"{cur:>10.3e} {prop:>10.3e} {ratio:>6.3f} {nv:>3}")

    # Within each generation, rank-order entries by current cost vs proposed.
    n_gens = len(rows) // POP
    flips_total = 0
    for g in range(n_gens):
        chunk = rows[g * POP:(g + 1) * POP]
        cur = [(i, r["cost"]) for i, r in enumerate(chunk)]
        prop = []
        for i, r in enumerate(chunk):
            p, _, _ = _proposed_from_perdir(r["per_dir"], r["per_dir_sigma"])
            prop.append((i, p))
        cur_rank = [i for i, _ in sorted(cur, key=lambda x: x[1])]
        prop_rank = [i for i, _ in sorted(prop, key=lambda x: x[1])]
        flip = int(cur_rank != prop_rank)
        flips_total += flip
        summary["selection_flips_per_gen"].append(dict(
            gen=g, current_rank=cur_rank, proposed_rank=prop_rank,
            same=cur_rank == prop_rank))
        print(f"  gen{g}: rank current={cur_rank}  proposed={prop_rank}  "
              f"{'SAME' if cur_rank == prop_rank else 'FLIPPED'}")

    print(f"\n  CMA-ES selection flips across {n_gens} generations: {flips_total}")
    summary["total_selection_flips"] = flips_total

    # Verify ratio matches predicted N_valid/3
    obs_ratios = np.array([e["ratio_cur_over_prop"]
                           for e in summary["per_eval"]
                           if e["n_valid"] > 0])
    n_v = summary["per_eval"][0]["n_valid"]
    pred = n_v / 3
    print(f"\n  observed ratio current/proposed: "
          f"mean={obs_ratios.mean():.4f}  std={obs_ratios.std():.2e}")
    print(f"  predicted N_valid/3 = {n_v}/3 = {pred:.4f}")

    # ---- Plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))

        eids = [e["eid"] for e in summary["per_eval"]]
        cur = [e["cost_current"] for e in summary["per_eval"]]
        prop = [e["cost_proposed"] for e in summary["per_eval"]]
        cur_sig = [e["sigma_current"] for e in summary["per_eval"]]
        prop_sig = [e["sigma_proposed"] for e in summary["per_eval"]]
        ax = axs[0]
        ax.errorbar(eids, cur, yerr=cur_sig, fmt="o-", color="C0",
                    label="current", capsize=3)
        ax.errorbar(eids, prop, yerr=prop_sig, fmt="s--", color="C1",
                    label="proposed (T3)", capsize=3)
        ax.set_xlabel("eval id"); ax.set_ylabel("cost")
        ax.set_title(f"Real MC CMA-ES: per-eval cost\n"
                     f"ref={ref_geom}  test={test_geom} (gcds={gcds})")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axs[1]
        cur_best = np.minimum.accumulate(cur)
        prop_best = np.minimum.accumulate(prop)
        ax.plot(eids, cur_best, "o-", color="C0", label="current")
        ax.plot(eids, prop_best, "s--", color="C1", label="proposed (T3)")
        ax.set_xlabel("eval id"); ax.set_ylabel("best-so-far cost")
        ax.set_title(f"Best-so-far convergence (CMA selection identical:"
                     f" {flips_total}/{n_gens} flips)")
        ax.legend(); ax.grid(alpha=0.3)

        plt.tight_layout()
        png = os.path.join(OUT, "realmc_cmaes.png")
        plt.savefig(png, dpi=130)
        plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, "realmc_cmaes.log"), "w") as f:
        f.write(log_buf.getvalue())
    with open(os.path.join(OUT, "realmc_cmaes_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
