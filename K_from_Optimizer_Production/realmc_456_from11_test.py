#!/usr/bin/env python3
"""
realmc_456_from11_test.py — same setup as realmc_456_t3_test.py but
starts CMA-ES from x0=(1.0, 1.0) so the migration toward (5.07, 7.74)
is visible.

  x0       = (1.0, 1.0)
  σ0       = 2.0      (need wide spread; truth is ~7 away)
  popsize  = 8
  max_evals= 80       (10 generations)
  n_traj   = 8000
  bounds   = lower (0.5, 0.5)

T3 patch installed (identity on this geometry; gcds=[16,16,16]).
"""
from __future__ import annotations

import io
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)
RUN = os.path.join(OUT, "realmc_456_from11")
os.makedirs(RUN, exist_ok=True)

R1_TRUE, R2_TRUE, BC_TRUE = 5.0652, 7.7429, 0.06283
X0 = (1.0, 1.0)
SIGMA0 = 2.0
POPSIZE = 8
MAX_EVALS = 80


def _install_t3_patch():
    import cost as cost_mod
    _orig = cost_mod._l2_cost_test_native

    def _patched(ref_data, test_data,
                 test_Lx, test_Ly, test_Tx, test_Ty,
                 ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                 copies=2, power=2):
        c, s, pd, ps = _orig(ref_data, test_data,
                             test_Lx, test_Ly, test_Tx, test_Ty,
                             ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                             copies=copies, power=power)
        valid = [(C, S) for C, S in zip(pd, ps) if not (C == 0.0 and S == 0.0)]
        if not valid:
            return float("nan"), float("nan"), pd, ps
        N = len(valid)
        return (sum(C for C, _ in valid) / N,
                math.sqrt(sum(S * S for _, S in valid)) / N,
                pd, ps)

    cost_mod._l2_cost_test_native = _patched
    print("  T3 patch installed")


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
    print(f"realmc_456_from11_test — CMA-ES from x0={X0}  σ0={SIGMA0}")
    print("=" * 76)

    _install_t3_patch()

    import mc_engine
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # noqa: F401  (installs reweight fallback)

    ref_geom  = (13, 16, -3,  3)
    test_geom = (16, 16,  0,  0)
    print(f"  ref={ref_geom}  test={test_geom}  TRUE=({R1_TRUE}, {R2_TRUE})")

    ref_dir = os.path.join(HERE, "results", "_reference_Lx13_Ly16_Tx-3_Ty3")
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))

    exe = os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram")
    mc_scratch = os.path.join(RUN, "mc_scratch")

    ev = Evaluator(
        exe=exe, ref_data=ref_data,
        ref_geom=ref_geom, test_geom=test_geom,
        output_dir=mc_scratch,
        n_traj_prod=8000,
        n_traj_scan_coarse=3000, n_traj_scan_fine=8000,
        scan_n_coarse=11, scan_n_refine=5, scan_n_refine2=5, scan_n_refine3=5,
        scan_max_shifts=4, scan_jackknife=True,
        beta_seed=(0.05, 0.50),
        reweight=False,
        cost_mode="test_native", cost_power=2,
    )

    t0 = time.perf_counter()
    run_cmaes(
        ev,
        x0=X0,
        max_evals=MAX_EVALS,
        sigma0=SIGMA0,
        popsize=POPSIZE,
        tolx=0.01, tolfun=1e-6,
        snr_floor=0.0, indist_stop_snr=0.0, snr_target=2.0,
        snr_max_traj_factor=4,
        seed=2026,
        pool=None,
        lower_bounds=(0.5, 0.5),
    )
    wall = time.perf_counter() - t0
    print(f"\nwall: {wall:.0f}s")

    rows = []
    with open(os.path.join(mc_scratch, "eval_log.jsonl")) as f:
        for line in f:
            line = line.strip()
            if line: rows.append(json.loads(line))
    if not rows:
        print("NO EVALS"); return

    def _dist(r): return float(np.hypot(r["r1"] - R1_TRUE, r["r2"] - R2_TRUE))

    best = min(rows, key=lambda r: r["cost"])
    closest = min(rows, key=_dist)
    first_close = next((i for i, r in enumerate(rows, 1) if _dist(r) < 1.0),
                       None)

    print(f"\n  best by COST:     ({best['r1']:.3f}, {best['r2']:.3f})  "
          f"cost={best['cost']:.4e}  β={best['beta_c']:.4f}  dist={_dist(best):.3f}")
    print(f"  closest to TRUE: ({closest['r1']:.3f}, {closest['r2']:.3f})  "
          f"cost={closest['cost']:.4e}  β={closest['beta_c']:.4f}  dist={_dist(closest):.3f}")
    print(f"  first eval within dist<1.0: "
          f"{first_close if first_close else '> budget'}")

    n_gens = len(rows) // POPSIZE
    gen_best = []
    gen_means = []
    print()
    for g in range(n_gens):
        chunk = rows[g * POPSIZE:(g + 1) * POPSIZE]
        b = min(chunk, key=lambda r: r["cost"])
        mean_r1 = float(np.mean([r["r1"] for r in chunk]))
        mean_r2 = float(np.mean([r["r2"] for r in chunk]))
        gen_best.append(b)
        gen_means.append((mean_r1, mean_r2))
        print(f"  gen{g}: pop-mean ({mean_r1:.2f}, {mean_r2:.2f})  "
              f"best ({b['r1']:.2f}, {b['r2']:.2f}) cost={b['cost']:.3e} "
              f"dist={_dist(b):.2f}")

    summary = dict(
        x0=X0, sigma0=SIGMA0, popsize=POPSIZE, max_evals=MAX_EVALS,
        wall_s=round(wall, 1), n_evals=len(rows),
        best=dict(r1=best["r1"], r2=best["r2"],
                  cost=best["cost"], beta_c=best["beta_c"], dist=_dist(best)),
        closest=dict(r1=closest["r1"], r2=closest["r2"],
                     cost=closest["cost"], beta_c=closest["beta_c"],
                     dist=_dist(closest)),
        first_close_eval=first_close,
        true=dict(r1=R1_TRUE, r2=R2_TRUE, beta_c=BC_TRUE),
        per_gen=[dict(gen=i, mean_r1=m[0], mean_r2=m[1],
                      best_r1=b["r1"], best_r2=b["r2"],
                      best_cost=b["cost"], best_dist=_dist(b))
                 for i, (m, b) in enumerate(zip(gen_means, gen_best))],
    )

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        ax = axs[0]
        rs1 = np.array([r["r1"] for r in rows])
        rs2 = np.array([r["r2"] for r in rows])
        cs  = np.array([r["cost"] for r in rows])
        sc = ax.scatter(rs1, rs2, c=np.log10(cs + 1e-12),
                        cmap="viridis", s=35, alpha=0.7)
        # generation-mean trajectory
        gm = np.array(gen_means)
        ax.plot([X0[0], *gm[:, 0]], [X0[1], *gm[:, 1]],
                "-o", color="red", lw=2, ms=6, label="pop-mean trajectory")
        ax.scatter([X0[0]], [X0[1]], marker="s", s=200,
                   facecolor="yellow", edgecolor="black", zorder=11,
                   label=f"x0={X0}")
        ax.scatter([R1_TRUE], [R2_TRUE], marker="*", s=400,
                   color="red", edgecolor="black", zorder=12, label="TRUE")
        ax.scatter([best["r1"]], [best["r2"]], marker="X", s=160,
                   color="orange", zorder=10, label="best by cost")
        for g, (m1, m2) in enumerate(gen_means):
            ax.annotate(f"g{g}", (m1, m2), fontsize=9, color="black",
                        xytext=(4, 4), textcoords="offset points")
        plt.colorbar(sc, ax=ax, label="log₁₀ cost")
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_title(f"CMA-ES migration from (1, 1)  →  best ({best['r1']:.2f}, "
                     f"{best['r2']:.2f}), dist={_dist(best):.2f}")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)

        ax = axs[1]
        eids = np.arange(1, len(rows) + 1)
        dists = np.array([_dist(r) for r in rows])
        best_dist = np.minimum.accumulate(dists)
        ax.plot(eids, dists, "o", color="C0", alpha=0.5, label="per-eval")
        ax.plot(eids, best_dist, "-", color="C3", lw=2,
                label="best-so-far (by dist)")
        ax.axhline(1.0, color="gray", ls=":", label="dist=1.0")
        for g in range(1, n_gens):
            ax.axvline(g * POPSIZE + 0.5, color="lightgray", lw=0.5)
        ax.set_xlabel("eval id"); ax.set_ylabel("dist to TRUE")
        ax.set_yscale("log")
        ax.set_title(f"Convergence  (first dist<1: "
                     f"{first_close if first_close else '> budget'})")
        ax.legend(); ax.grid(alpha=0.3, which="both")

        plt.tight_layout()
        png = os.path.join(OUT, "realmc_456_from11.png")
        plt.savefig(png, dpi=130)
        plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, "realmc_456_from11_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    with open(os.path.join(OUT, "realmc_456_from11.log"), "w") as f:
        f.write(log_buf.getvalue())
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
