#!/usr/bin/env python3
"""
realmc_456_t3_test.py — real-MC CMA-ES on the 4-5-6 triangle problem with
the T3-patched ``test_native`` cost installed by monkey-patch (no source
edits).

NOTE on cost equivalence
------------------------
Test geometry is (16, 16, 0, 0) → boundary path gcds = [16, 16, 16] → all
3 directions are valid. The T3 patch (drop dirs where (C, σ) == (0, 0))
is mathematically a no-op here. So this run demonstrates that the patched
cost CAN find the 4-5-6 solution; it does not exercise the bias-correction
on this particular problem.

Setup
-----
  ref  = (13, 16, -3,  3)   ← 4:5:6 triangle reference (cached)
  test = (16, 16,  0,  0)
  TRUE = (r1≈5.0652, r2≈7.7429, β_c≈0.06283)
  CMA-ES: x0=(5.0, 7.5), σ0=0.8, popsize=8, max_evals=48 (6 gens)
  n_traj=8000  (cheap; serial)

Memory hints applied
--------------------
  - σ₀ = 0.8 (sharper than the failed previous run's σ₀=2.0)
  - x0 = (5.0, 7.5) (closer to truth than failed (5.0, 6.0))
  - popsize 8 (was 6); 48 evals (was 30)

Outputs
-------
  results/_native_triage/realmc_456/
  results/_native_triage/realmc_456_summary.json
  results/_native_triage/realmc_456.png
  results/_native_triage/realmc_456.log
"""
from __future__ import annotations

import contextlib
import io
import json
import math
import os
import sys
import time
from typing import List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)
RUN = os.path.join(OUT, "realmc_456")
os.makedirs(RUN, exist_ok=True)

R1_TRUE, R2_TRUE, BC_TRUE = 5.0652, 7.7429, 0.06283


# ---------------------------------------------------------------------------
# Monkey-patch test_native cost with T3 (drop-collapsed) logic.
# ---------------------------------------------------------------------------
def _install_t3_patch():
    import cost as cost_mod
    _orig = cost_mod._l2_cost_test_native

    def _patched_test_native(ref_data, test_data,
                             test_Lx, test_Ly, test_Tx, test_Ty,
                             ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                             copies=2, power=2):
        c, s, pd, ps = _orig(ref_data, test_data,
                             test_Lx, test_Ly, test_Tx, test_Ty,
                             ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                             copies=copies, power=power)
        # T3: drop dirs where (C, σ) == (0, 0); plain mean over survivors.
        valid = [(C, S) for C, S in zip(pd, ps) if not (C == 0.0 and S == 0.0)]
        if not valid:
            print("    [T3] all 3 dirs collapsed → returning NaN")
            return float("nan"), float("nan"), pd, ps
        N = len(valid)
        cost = sum(C for C, _ in valid) / N
        sig = math.sqrt(sum(S * S for _, S in valid)) / N
        if N < 3:
            print(f"    [T3] dropped {3-N} collapsed dir(s); cost {c:.3e} → "
                  f"{cost:.3e}  (factor {cost/c if c else 'inf':.3f})")
        return cost, sig, pd, ps

    cost_mod._l2_cost_test_native = _patched_test_native
    print(f"  T3 patch installed on cost._l2_cost_test_native")


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
    print("realmc_456_t3_test — CMA-ES on 4-5-6 triangle with T3-patched cost")
    print("=" * 76)

    _install_t3_patch()

    import mc_engine
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # installs reweight fallback patch

    ref_geom  = (13, 16, -3,  3)
    test_geom = (16, 16,  0,  0)

    paths = []
    for dm, dn in [(test_geom[0], test_geom[3]),
                   (test_geom[2], -test_geom[1]),
                   (-test_geom[0] - test_geom[2], test_geom[1] - test_geom[3])]:
        paths.append((dm, dn))
    gcds = [math.gcd(abs(a), abs(b)) for a, b in paths]
    print(f"  ref={ref_geom}  test={test_geom}")
    print(f"  test boundary path gcds = {gcds}  → "
          f"T3 effect: {'IDENTITY' if all(g >= 2 for g in gcds) else 'ACTIVE'}")
    print(f"  TRUE: r1={R1_TRUE}  r2={R2_TRUE}  β_c={BC_TRUE}")

    ref_dir = os.path.join(HERE, "results", "_reference_Lx13_Ly16_Tx-3_Ty3")
    ref_a2a = os.path.join(ref_dir, "two_point_all_to_all.dat")
    if not os.path.exists(ref_a2a):
        raise FileNotFoundError(ref_a2a)
    ref_data = mc_engine.load_all_to_all(ref_a2a)

    exe = os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram")
    mc_scratch = os.path.join(RUN, "mc_scratch")

    ev = Evaluator(
        exe=exe, ref_data=ref_data,
        ref_geom=ref_geom, test_geom=test_geom,
        output_dir=mc_scratch,
        n_traj_prod=8000,
        n_traj_scan_coarse=3000,
        n_traj_scan_fine=8000,
        scan_n_coarse=11,
        scan_n_refine=5, scan_n_refine2=5, scan_n_refine3=5,
        scan_max_shifts=4, scan_jackknife=True,
        beta_seed=(0.05, 0.50),     # 4-5-6 region (β_c ≈ 0.063)
        reweight=False,
        cost_mode="test_native", cost_power=2,
    )

    t0 = time.perf_counter()
    run_cmaes(
        ev,
        x0=(5.0, 7.5),
        max_evals=48,
        sigma0=0.8,
        popsize=8,
        tolx=0.01, tolfun=1e-6,
        snr_floor=0.0, indist_stop_snr=0.0, snr_target=2.0,
        snr_max_traj_factor=4,
        seed=2026,
        pool=None,
        lower_bounds=(0.5, 0.5),
    )
    wall = time.perf_counter() - t0
    print(f"\nCMA-ES wall time: {wall:.0f}s")

    # Replay log
    log_path = os.path.join(mc_scratch, "eval_log.jsonl")
    rows = []
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    if not rows:
        print("NO EVALS"); return

    def _dist(r): return float(np.hypot(r["r1"] - R1_TRUE, r["r2"] - R2_TRUE))

    best = min(rows, key=lambda r: r["cost"])
    closest = min(rows, key=_dist)

    print(f"\n  best by COST:      ({best['r1']:.3f}, {best['r2']:.3f})  "
          f"cost={best['cost']:.4e}  β_c={best['beta_c']:.4f}  "
          f"dist={_dist(best):.3f}")
    print(f"  closest to TRUE:  ({closest['r1']:.3f}, {closest['r2']:.3f})  "
          f"cost={closest['cost']:.4e}  β_c={closest['beta_c']:.4f}  "
          f"dist={_dist(closest):.3f}")

    first_close = next((i for i, r in enumerate(rows, 1) if _dist(r) < 1.0),
                       None)
    print(f"  first eval within dist<1.0:  "
          f"{first_close if first_close else '> budget'}")

    # Trajectory: per-generation best (popsize=8)
    POP = 8
    n_gens = len(rows) // POP
    gen_best = []
    for g in range(n_gens):
        chunk = rows[g * POP:(g + 1) * POP]
        b = min(chunk, key=lambda r: r["cost"])
        gen_best.append(b)
        print(f"  gen{g}: best ({b['r1']:.3f}, {b['r2']:.3f})  "
              f"cost={b['cost']:.4e}  dist={_dist(b):.3f}")

    summary = dict(
        ref_geom=ref_geom, test_geom=test_geom,
        boundary_gcds=gcds,
        wall_s=round(wall, 1),
        n_evals=len(rows),
        best=dict(r1=best["r1"], r2=best["r2"],
                  cost=best["cost"], beta_c=best["beta_c"],
                  dist=_dist(best)),
        closest=dict(r1=closest["r1"], r2=closest["r2"],
                     cost=closest["cost"], beta_c=closest["beta_c"],
                     dist=_dist(closest)),
        first_close_eval=first_close,
        true=dict(r1=R1_TRUE, r2=R2_TRUE, beta_c=BC_TRUE),
        per_gen_best=[dict(gen=i, r1=b["r1"], r2=b["r2"],
                           cost=b["cost"], dist=_dist(b))
                      for i, b in enumerate(gen_best)],
    )

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
        ax = axs[0]
        rs1 = np.array([r["r1"] for r in rows])
        rs2 = np.array([r["r2"] for r in rows])
        cs = np.array([r["cost"] for r in rows])
        sc = ax.scatter(rs1, rs2, c=np.log10(cs + 1e-12),
                        cmap="viridis", s=40)
        ax.scatter([R1_TRUE], [R2_TRUE], marker="*", s=300,
                   color="red", zorder=10, label="true")
        ax.scatter([best["r1"]], [best["r2"]], marker="X", s=140,
                   color="orange", zorder=9, label="best by cost")
        ax.scatter([closest["r1"]], [closest["r2"]], marker="o",
                   facecolor="none", edgecolor="cyan", s=180, lw=2,
                   zorder=9, label="closest to true")
        for i, r in enumerate(rows):
            ax.annotate(str(i + 1), (r["r1"], r["r2"]), fontsize=6,
                        color="white")
        plt.colorbar(sc, ax=ax, label="log₁₀ cost")
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_title(f"4-5-6 CMA-ES (T3-patched test_native)\n"
                     f"true=(5.07, 7.74)  best=({best['r1']:.2f}, "
                     f"{best['r2']:.2f})  dist={_dist(best):.2f}")
        ax.legend(loc="lower right")
        ax.grid(alpha=0.3)

        ax = axs[1]
        eids = np.arange(1, len(rows) + 1)
        dists = np.array([_dist(r) for r in rows])
        best_dist = np.minimum.accumulate(dists)
        ax.plot(eids, dists, "o", color="C0", alpha=0.5, label="per-eval")
        ax.plot(eids, best_dist, "-", color="C3", lw=2,
                label="best-so-far (by dist)")
        ax.axhline(1.0, color="gray", ls=":", label="threshold (1.0)")
        ax.set_xlabel("eval id"); ax.set_ylabel("dist to TRUE")
        ax.set_yscale("log")
        ax.set_title(f"Convergence  ({first_close if first_close else '>budget'} "
                     "evals to dist<1)")
        ax.legend(); ax.grid(alpha=0.3, which="both")

        plt.tight_layout()
        png = os.path.join(OUT, "realmc_456.png")
        plt.savefig(png, dpi=130)
        plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, "realmc_456_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    with open(os.path.join(OUT, "realmc_456.log"), "w") as f:
        f.write(log_buf.getvalue())
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
