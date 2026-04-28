#!/usr/bin/env python3
"""Workflow stress test on twisted-equilateral geometries.

Cycle lengths on the triangular torus (Lx, Ly, Tx, Ty):
    |v|^2 = Lx^2 + Lx*Ty + Ty^2
    |u|^2 = Tx^2 - Tx*Ly + Ly^2
    |w|^2 = (Lx+Tx)^2 - (Lx+Tx)(Ly-Ty) + (Ly-Ty)^2

All three are equal iff (Lx, Ly, Tx, Ty) = (L, L+T, T, T) (up to lattice
symmetry).  This script runs the full production workflow with
ref == test == one of these twisted-equilateral tori, so the ground truth
optimum is (r1, r2) = (1, 1).  Verifies that the workflow handles twist.
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


def cycle_len(Lx, Ly, Tx, Ty):
    v = math.sqrt(Lx * Lx + Lx * Ty + Ty * Ty)
    u = math.sqrt(Tx * Tx - Tx * Ly + Ly * Ly)
    w = math.sqrt((Lx + Tx) ** 2 - (Lx + Tx) * (Ly - Ty) + (Ly - Ty) ** 2)
    return v, u, w


# (label, Lx, Ly, Tx, Ty)
DEFAULT_GEOMS = [
    ("equi_T0_L12",   12, 12,  0,  0),   # untwisted baseline
    ("twist_T1_L4",    4,  5,  1,  1),   # smallest non-trivial twist (small)
    ("twist_T2_L5",    5,  7,  2,  2),
    ("twist_T3_L8",    8, 11,  3,  3),
    ("twist_Tm3_L12", 12,  9, -3, -3),   # negative twist
    ("twist_T3_L11",  11, 14,  3,  3),   # near L=12 scale
]


def run_one(label: str, Lx: int, Ly: int, Tx: int, Ty: int,
            *, max_evals: int, popsize: int, sigma0: float,
            n_traj_prod: int, reweight_n_traj: int,
            x0: tuple[float, float], n_workers: int,
            results_dir: str) -> dict:
    run_name = f"twst_{label}_{datetime.now().strftime('%H%M%S')}"
    out_dir = os.path.join(results_dir, run_name)
    sides = cycle_len(Lx, Ly, Tx, Ty)
    print(f"\n==== {label}  geom=({Lx},{Ly},{Tx},{Ty})  cycle≈{sides[0]:.3f}"
          f"  →  {out_dir} ====", flush=True)
    cfg_path = os.path.join(results_dir, f"{run_name}_cfg.json")
    with open(cfg_path, "w") as f:
        json.dump({
            "ref_Tx": int(Tx), "ref_Ty": int(Ty),
            "test_Tx": int(Tx), "test_Ty": int(Ty),
        }, f)
    cmd = [
        sys.executable, "-u", "run.py",
        "--config", cfg_path,
        "--test-Lx", str(Lx), "--test-Ly", str(Ly),
        "--ref-Lx",  str(Lx), "--ref-Ly",  str(Ly),
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
        print(f"[{label}] FAILED rc={proc.returncode}; see {log_path}",
              flush=True)
        return {"label": label, "geom": [Lx, Ly, Tx, Ty],
                "status": "failed", "wall_s": dt, "log": log_path}

    summary_path = os.path.join(out_dir, "summary.json")
    if not os.path.exists(summary_path):
        print(f"[{label}] no summary.json; see {log_path}", flush=True)
        return {"label": label, "geom": [Lx, Ly, Tx, Ty],
                "status": "no_summary", "wall_s": dt, "log": log_path}
    with open(summary_path) as f:
        summ = json.load(f)
    eval_log = os.path.join(out_dir, "eval_log.jsonl")
    rows = [json.loads(line) for line in open(eval_log)]

    tail = rows[-5:]
    r1_tail = float(np.mean([r["r1"] for r in tail]))
    r2_tail = float(np.mean([r["r2"] for r in tail]))
    cma_mean = summ["final_gaussian"]["mean"]

    # Verify that the simulator actually used the twisted geometry by
    # checking config from summary.
    sim_geom = (summ["config"]["test_Lx"], summ["config"]["test_Ly"],
                summ["config"].get("test_Tx", 0),
                summ["config"].get("test_Ty", 0))
    geom_ok = (sim_geom == (Lx, Ly, Tx, Ty))

    out = {
        "label": label,
        "geom": [Lx, Ly, Tx, Ty],
        "sim_geom": list(sim_geom),
        "geom_ok": geom_ok,
        "cycle_lens": list(sides),
        "status": summ["status"],
        "n_evals": summ["n_evals"],
        "best_eval_id": summ["best_eval_id"],
        "best_r1": summ["best_r1"], "best_r2": summ["best_r2"],
        "best_cost": summ["best_cost"], "best_sigma": summ["best_sigma"],
        "best_snr": summ["best_snr"], "best_beta_c": summ["best_beta_c"],
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
    out["dist_trail5_to_one"]= math.hypot(r1_tail - 1.0, r2_tail - 1.0)
    print(f"[{label}]  cycle≈({sides[0]:.2f},{sides[1]:.2f},{sides[2]:.2f})"
          f"  best=({out['best_r1']:.4f},{out['best_r2']:.4f})"
          f"  CMA-mean=({cma_mean[0]:.4f},{cma_mean[1]:.4f})"
          f"  d_mean={out['dist_mean_to_one']:.4f}"
          f"  geom_ok={geom_ok}"
          f"  wall={dt:.1f}s", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-evals", type=int, default=40)
    ap.add_argument("--popsize", type=int, default=6)
    ap.add_argument("--sigma0", type=float, default=0.20)
    ap.add_argument("--n-traj-prod", type=int, default=10000)
    ap.add_argument("--reweight-n-traj", type=int, default=20000)
    ap.add_argument("--x0", type=float, nargs=2, default=[1.3, 0.8])
    ap.add_argument("--n-workers", type=int, default=2)
    ap.add_argument("--geoms", type=str, nargs="+", default=None)
    ap.add_argument("--out-dir", type=str,
                    default=os.path.join("results", "_stress_twisted"))
    args = ap.parse_args()

    out_root = (args.out_dir if os.path.isabs(args.out_dir)
                else os.path.join(_HERE, args.out_dir))
    os.makedirs(out_root, exist_ok=True)

    geoms = DEFAULT_GEOMS
    if args.geoms:
        sel = set(args.geoms)
        geoms = [g for g in geoms if g[0] in sel]

    rows = []
    t_start = time.time()
    for label, Lx, Ly, Tx, Ty in geoms:
        rec = run_one(label, Lx, Ly, Tx, Ty,
                      max_evals=args.max_evals, popsize=args.popsize,
                      sigma0=args.sigma0,
                      n_traj_prod=args.n_traj_prod,
                      reweight_n_traj=args.reweight_n_traj,
                      x0=tuple(args.x0), n_workers=args.n_workers,
                      results_dir=out_root)
        rows.append(rec)

    total_wall = time.time() - t_start
    with open(os.path.join(out_root, "stress_twisted_summary.json"), "w") as f:
        json.dump({"x0": args.x0, "max_evals": args.max_evals,
                   "rows": rows, "total_wall_s": total_wall}, f, indent=2)

    print()
    print("=" * 110)
    hdr = (f"{'label':<15}{'geom':>17}{'cycle':>9}{'best (r1,r2)':>20}"
           f"{'CMA mean':>20}{'d_best':>9}{'d_mean':>9}{'geom_ok':>9}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        if r.get("status") in ("failed", "no_summary"):
            print(f"{r['label']:<15}  FAILED ({r.get('status')})")
            continue
        g = r["geom"]; m = r["cma_mean"]
        cyc = r["cycle_lens"][0]
        geom_str = f"({g[0]},{g[1]},{g[2]},{g[3]})"
        best_str = f"({r['best_r1']:.4f},{r['best_r2']:.4f})"
        mean_str = f"({m[0]:.4f},{m[1]:.4f})"
        print(f"{r['label']:<15}{geom_str:>17}{cyc:>9.2f}"
              f"{best_str:>20}{mean_str:>20}"
              f"{r['dist_best_to_one']:>9.4f}{r['dist_mean_to_one']:>9.4f}"
              f"{str(r['geom_ok']):>9}")
    print(f"\ntotal wall: {total_wall:.1f}s")

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keep = [r for r in rows if r.get("status") not in ("failed", "no_summary")]
        labels = [r["label"] for r in keep]
        d_b = np.array([r["dist_best_to_one"] for r in keep])
        d_m = np.array([r["dist_mean_to_one"] for r in keep])
        d_t = np.array([r["dist_trail5_to_one"] for r in keep])
        cyc = np.array([r["cycle_lens"][0] for r in keep])
        wall = np.array([r["wall_s"] for r in keep])

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
        x = np.arange(len(keep))
        ax = axes[0]
        ax.bar(x - 0.27, d_b, 0.27, label="best", color="C0")
        ax.bar(x        , d_m, 0.27, label="CMA mean", color="C1")
        ax.bar(x + 0.27, d_t, 0.27, label="trail-5", color="C2")
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right",
                                              fontsize=8)
        ax.set_ylabel("Euclidean distance from (1, 1)")
        ax.set_title("Workflow optimum vs ground truth across twisted geoms")
        ax.set_yscale("log"); ax.grid(alpha=0.3, which="both", axis="y")
        ax.legend(loc="best", fontsize=8)
        for i, c in enumerate(cyc):
            ax.annotate(f"|c|={c:.1f}", (i, max(d_b[i], d_m[i], d_t[i]) * 1.1),
                        ha="center", fontsize=7)

        ax = axes[1]
        ax.bar(x, wall, color="C0")
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right",
                                              fontsize=8)
        ax.set_ylabel("wall time (s)")
        ax.set_title("Run cost per geometry")
        ax.grid(alpha=0.3, axis="y")

        fig.suptitle("Workflow stress test — twisted-equilateral ref=test",
                     fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out_png = os.path.join(out_root, "stress_twisted.png")
        fig.savefig(out_png, dpi=130)
        print(f"[plot] {out_png}")
    except Exception as e:
        print(f"[plot] skipped ({e})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
