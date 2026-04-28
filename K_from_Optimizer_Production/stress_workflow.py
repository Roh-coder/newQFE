#!/usr/bin/env python3
"""End-to-end workflow stress test.

For several lattice sizes L ∈ {sizes}, run the full production workflow
with ref=test=(L,L,0,0).  Ground truth: (r1,r2) = (1,1).

For each L, report:
  - best eval (r1,r2,cost,β_c)
  - CMA mean at last generation
  - mean of trailing 5 evals (low-noise smoothed estimate)
  - distance from (1,1)
  - wall time and total trajectories

A clean pass: distance from (1,1) shrinks to within ~σ_r1, σ_r2 set by
trajectory budget, and shrinks (or stays flat) with L.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def run_one(L: int, *, max_evals: int, popsize: int, sigma0: float,
            n_traj_prod: int, reweight_n_traj: int,
            x0: tuple[float, float], n_workers: int,
            results_dir: str) -> dict:
    run_name = f"stress_L{L}_{datetime.now().strftime('%H%M%S')}"
    out_dir = os.path.join(results_dir, run_name)
    print(f"\n==== L={L}  →  {out_dir} ====", flush=True)
    cmd = [
        sys.executable, "-u", "run.py",
        "--test-Lx", str(L), "--test-Ly", str(L),
        "--ref-Lx",  str(L), "--ref-Ly",  str(L),
        "--max-evals", str(max_evals),
        "--cma-popsize", str(popsize),
        "--cma-sigma0", str(sigma0),
        "--x0", f"{x0[0]}", f"{x0[1]}",
        "--n-traj-prod", str(n_traj_prod),
        "--reweight-n-traj", str(reweight_n_traj),
        "--multidonor-2pass",
        "--indist-stop-snr", "0.0",
        "--n-workers", str(n_workers),
        "--save-every", "5",
        "--no-dashboard", "--no-vis",
        "--results-dir", results_dir,
        "--run-name", run_name,
    ]
    t0 = time.time()
    log_path = os.path.join(results_dir, f"{run_name}.log")
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, cwd=_HERE, stdout=logf,
                              stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if proc.returncode != 0:
        print(f"[L={L}] FAILED rc={proc.returncode}; see {log_path}",
              flush=True)
        return {"L": L, "status": "failed", "wall_s": dt, "log": log_path}

    summary_path = os.path.join(out_dir, "summary.json")
    with open(summary_path) as f:
        summ = json.load(f)
    eval_log = os.path.join(out_dir, "eval_log.jsonl")
    rows = [json.loads(line) for line in open(eval_log)]

    # Trailing-5 mean of (r1, r2)
    tail = rows[-5:]
    r1_tail = float(np.mean([r["r1"] for r in tail]))
    r2_tail = float(np.mean([r["r2"] for r in tail]))

    cma_mean = summ["final_gaussian"]["mean"]
    out = {
        "L": L,
        "status": summ["status"],
        "n_evals": summ["n_evals"],
        "best_eval_id": summ["best_eval_id"],
        "best_r1": summ["best_r1"], "best_r2": summ["best_r2"],
        "best_cost": summ["best_cost"],
        "best_sigma": summ["best_sigma"],
        "best_snr": summ["best_snr"],
        "best_beta_c": summ["best_beta_c"],
        "cma_mean": cma_mean,
        "cma_sigma": summ["final_gaussian"]["sigma"],
        "trail5_r1": r1_tail, "trail5_r2": r2_tail,
        "total_traj": summ["total_traj"],
        "wall_s": dt,
        "out_dir": out_dir,
    }
    out["dist_best_to_one"]  = math.hypot(out["best_r1"] - 1.0,
                                          out["best_r2"] - 1.0)
    out["dist_mean_to_one"]  = math.hypot(cma_mean[0]   - 1.0,
                                          cma_mean[1]   - 1.0)
    out["dist_trail5_to_one"]= math.hypot(out["trail5_r1"] - 1.0,
                                          out["trail5_r2"] - 1.0)
    print(f"[L={L}]  best (r1,r2)=({out['best_r1']:.4f},{out['best_r2']:.4f})"
          f"  CMA-mean=({cma_mean[0]:.4f},{cma_mean[1]:.4f})"
          f"  trail5=({r1_tail:.4f},{r2_tail:.4f})"
          f"  d_mean={out['dist_mean_to_one']:.4f}"
          f"  wall={dt:.1f}s", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[8, 12, 16, 20])
    ap.add_argument("--max-evals", type=int, default=60)
    ap.add_argument("--popsize", type=int, default=6)
    ap.add_argument("--sigma0", type=float, default=0.20)
    ap.add_argument("--n-traj-prod", type=int, default=10000)
    ap.add_argument("--reweight-n-traj", type=int, default=20000)
    ap.add_argument("--x0", type=float, nargs=2, default=[1.4, 0.7])
    ap.add_argument("--n-workers", type=int, default=2)
    ap.add_argument("--out-dir", type=str,
                    default=os.path.join("results", "_stress_workflow"))
    args = ap.parse_args()

    out_root = (args.out_dir if os.path.isabs(args.out_dir)
                else os.path.join(_HERE, args.out_dir))
    os.makedirs(out_root, exist_ok=True)

    rows = []
    t_start = time.time()
    for L in args.sizes:
        rec = run_one(L, max_evals=args.max_evals, popsize=args.popsize,
                      sigma0=args.sigma0,
                      n_traj_prod=args.n_traj_prod,
                      reweight_n_traj=args.reweight_n_traj,
                      x0=tuple(args.x0), n_workers=args.n_workers,
                      results_dir=out_root)
        rows.append(rec)

    total_wall = time.time() - t_start
    with open(os.path.join(out_root, "stress_workflow_summary.json"), "w") as f:
        json.dump({"x0": args.x0, "max_evals": args.max_evals,
                   "n_traj_prod": args.n_traj_prod,
                   "rows": rows, "total_wall_s": total_wall}, f, indent=2)

    print()
    print("=" * 90)
    print(f"{'L':>4}{'best (r1,r2)':>22}{'CMA mean':>22}{'trail5 (r1,r2)':>22}"
          f"{'d_best':>10}{'d_mean':>10}{'d_trail':>10}")
    print("-" * 100)
    for r in rows:
        if r.get("status") == "failed":
            print(f"{r['L']:>4}  FAILED")
            continue
        m = r["cma_mean"]
        print(f"{r['L']:>4}  ({r['best_r1']:.4f},{r['best_r2']:.4f})"
              f"   ({m[0]:.4f},{m[1]:.4f})"
              f"   ({r['trail5_r1']:.4f},{r['trail5_r2']:.4f})"
              f"   {r['dist_best_to_one']:>7.4f}"
              f"   {r['dist_mean_to_one']:>7.4f}"
              f"   {r['dist_trail5_to_one']:>7.4f}")
    print(f"\ntotal wall: {total_wall:.1f}s")

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        Ls = np.array([r["L"] for r in rows if r.get("status") != "failed"])
        d_b = np.array([r["dist_best_to_one"] for r in rows if r.get("status") != "failed"])
        d_m = np.array([r["dist_mean_to_one"] for r in rows if r.get("status") != "failed"])
        d_t = np.array([r["dist_trail5_to_one"] for r in rows if r.get("status") != "failed"])
        wall = np.array([r["wall_s"] for r in rows if r.get("status") != "failed"])
        ntraj = np.array([r["total_traj"] for r in rows if r.get("status") != "failed"])

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
        ax = axes[0]
        ax.plot(Ls, d_b, "o-", color="C0", label="distance(best, (1,1))")
        ax.plot(Ls, d_m, "s-", color="C1", label="distance(CMA mean, (1,1))")
        ax.plot(Ls, d_t, "^-", color="C2", label="distance(trail-5 mean, (1,1))")
        ax.set_xlabel("L"); ax.set_ylabel("Euclidean distance from (1,1)")
        ax.set_title("Workflow convergence to known answer")
        ax.set_yscale("log"); ax.grid(alpha=0.3, which="both")
        ax.legend(loc="best", fontsize=8)
        for L, d in zip(Ls, d_m):
            ax.annotate(f"L={L}", (L, d), xytext=(4, 4),
                        textcoords="offset points", fontsize=8)

        ax = axes[1]
        ax.plot(Ls, wall, "o-", color="C0", label="wall (s)")
        ax.set_xlabel("L"); ax.set_ylabel("wall time (s)", color="C0")
        ax.tick_params(axis="y", labelcolor="C0")
        ax.grid(alpha=0.3)
        ax2 = ax.twinx()
        ax2.plot(Ls, ntraj / 1e6, "s--", color="crimson", label="total trajectories (M)")
        ax2.set_ylabel("total trajectories (M)", color="crimson")
        ax2.tick_params(axis="y", labelcolor="crimson")
        ax.set_title("Cost per workflow run")

        fig.suptitle("Workflow stress test — ref=test, ground truth (1,1)",
                     fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out_png = os.path.join(out_root, "stress_workflow.png")
        fig.savefig(out_png, dpi=130)
        print(f"[plot] {out_png}")
    except Exception as e:
        print(f"[plot] skipped ({e})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
