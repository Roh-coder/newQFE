#!/usr/bin/env python3
"""
cost_zoo_driver.py — A/B-test any cost_zoo variant via real-MC CMA-ES
on the 4-5-6 problem from x0=(1, 1).

Usage:
    python cost_zoo_driver.py <variant> [--max-evals N] [--popsize P]
        [--n-traj T] [--sigma0 S] [--seed S]

Variants: relative | log | arclength | effmass | pmean4

Example:
    python cost_zoo_driver.py arclength --max-evals 64

Outputs:
    results/_native_triage/zoo_<variant>/         (MC scratch + eval_log.jsonl)
    results/_native_triage/zoo_<variant>_summary.json
    results/_native_triage/zoo_<variant>.png
    results/_native_triage/zoo_<variant>.log
"""
from __future__ import annotations

import argparse, io, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT  = os.path.join(HERE, "results", "_native_triage")
R1_TRUE, R2_TRUE, BC_TRUE = 5.0652, 7.7429, 0.06283


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant", choices=["relative", "log", "arclength",
                                        "effmass", "pmean4"])
    ap.add_argument("--max-evals", type=int, default=80)
    ap.add_argument("--popsize",   type=int, default=8)
    ap.add_argument("--sigma0",    type=float, default=2.0)
    ap.add_argument("--n-traj",    type=int, default=8000)
    ap.add_argument("--seed",      type=int, default=2026)
    ap.add_argument("--x0", type=float, nargs=2, default=[1.0, 1.0])
    args = ap.parse_args()

    tag = f"zoo_{args.variant}"
    os.makedirs(OUT, exist_ok=True)
    RUN = os.path.join(OUT, tag); os.makedirs(RUN, exist_ok=True)

    log_buf = io.StringIO()

    class Tee:
        def __init__(self, *s): self.s = s
        def write(self, x):
            for s in self.s: s.write(x)
        def flush(self):
            for s in self.s: s.flush()
    sys.stdout = Tee(sys.__stdout__, log_buf)

    print("=" * 76)
    print(f"cost_zoo_driver — variant={args.variant}  x0={tuple(args.x0)} "
          f"σ0={args.sigma0}  pop={args.popsize}  evals={args.max_evals}  "
          f"n_traj={args.n_traj}")
    print("=" * 76)

    import cost_zoo
    cost_zoo.install(args.variant)

    import mc_engine
    from optimizer import run_cmaes
    from evaluator import Evaluator
    import prod_runtime  # noqa: F401

    ref_geom  = (13, 16, -3, 3)
    test_geom = (16, 16,  0, 0)
    print(f"  ref={ref_geom}  test={test_geom}  TRUE=({R1_TRUE}, {R2_TRUE})")

    ref_data = mc_engine.load_all_to_all(os.path.join(
        HERE, "results", "_reference_Lx13_Ly16_Tx-3_Ty3",
        "two_point_all_to_all.dat"))

    ev = Evaluator(
        exe=os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram"),
        ref_data=ref_data,
        ref_geom=ref_geom, test_geom=test_geom,
        output_dir=os.path.join(RUN, "mc_scratch"),
        n_traj_prod=args.n_traj,
        n_traj_scan_coarse=max(2000, args.n_traj // 3),
        n_traj_scan_fine=args.n_traj,
        scan_n_coarse=11, scan_n_refine=5, scan_n_refine2=5, scan_n_refine3=5,
        scan_max_shifts=4, scan_jackknife=True,
        beta_seed=(0.05, 0.50),
        reweight=False,
        cost_mode="test_native",   # zoo replaces _l2_cost_test_native
        cost_power=2,
    )

    t0 = time.perf_counter()
    run_cmaes(ev,
              x0=tuple(args.x0),
              max_evals=args.max_evals,
              sigma0=args.sigma0,
              popsize=args.popsize,
              tolx=0.01, tolfun=1e-6,
              snr_floor=0.0, indist_stop_snr=0.0, snr_target=2.0,
              snr_max_traj_factor=4,
              seed=args.seed, pool=None,
              lower_bounds=(0.5, 0.5))
    wall = time.perf_counter() - t0
    print(f"\nwall: {wall:.0f}s")

    rows = [json.loads(l) for l in
            open(os.path.join(RUN, "mc_scratch", "eval_log.jsonl"))
            if l.strip()]
    if not rows:
        print("NO EVALS"); return

    def _dist(r): return float(np.hypot(r["r1"] - R1_TRUE, r["r2"] - R2_TRUE))
    rows_finite = [r for r in rows if math.isfinite(r["cost"])]
    if not rows_finite:
        print("All costs NaN!"); return
    best    = min(rows_finite, key=lambda r: r["cost"])
    closest = min(rows, key=_dist)
    first_close = next((i for i, r in enumerate(rows, 1) if _dist(r) < 1.0),
                       None)

    print(f"\n  best by COST:    ({best['r1']:.3f}, {best['r2']:.3f})  "
          f"cost={best['cost']:.4e}  β={best['beta_c']:.4f}  "
          f"dist={_dist(best):.3f}")
    print(f"  closest to TRUE: ({closest['r1']:.3f}, {closest['r2']:.3f})  "
          f"cost={closest['cost']:.4e}  β={closest['beta_c']:.4f}  "
          f"dist={_dist(closest):.3f}")
    print(f"  first eval within dist<1.0: "
          f"{first_close if first_close else '> budget'}")

    POP = args.popsize
    n_gens = len(rows) // POP
    gen_means, gen_best = [], []
    print()
    for g in range(n_gens):
        c = rows[g * POP:(g + 1) * POP]
        m1 = float(np.mean([r["r1"] for r in c]))
        m2 = float(np.mean([r["r2"] for r in c]))
        cf = [r for r in c if math.isfinite(r["cost"])]
        b = min(cf, key=lambda r: r["cost"]) if cf else c[0]
        gen_means.append((m1, m2)); gen_best.append(b)
        print(f"  gen{g}: mean ({m1:.2f}, {m2:.2f})  best "
              f"({b['r1']:.2f}, {b['r2']:.2f}) cost={b['cost']:.3e} "
              f"dist={_dist(b):.2f}")

    summary = dict(
        variant=args.variant, x0=tuple(args.x0), sigma0=args.sigma0,
        popsize=args.popsize, max_evals=args.max_evals,
        n_traj=args.n_traj, seed=args.seed,
        wall_s=round(wall, 1), n_evals=len(rows),
        best=dict(r1=best["r1"], r2=best["r2"], cost=best["cost"],
                  beta_c=best["beta_c"], dist=_dist(best)),
        closest=dict(r1=closest["r1"], r2=closest["r2"], cost=closest["cost"],
                     beta_c=closest["beta_c"], dist=_dist(closest)),
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
        cs  = np.array([abs(r["cost"]) if math.isfinite(r["cost"]) else 1.0
                        for r in rows])
        sc = ax.scatter(rs1, rs2, c=np.log10(cs + 1e-12),
                        cmap="viridis", s=35, alpha=0.7)
        gm = np.array(gen_means)
        ax.plot([args.x0[0], *gm[:, 0]], [args.x0[1], *gm[:, 1]],
                "-o", color="red", lw=2, ms=6, label="pop-mean trajectory")
        ax.scatter([args.x0[0]], [args.x0[1]], marker="s", s=200,
                   facecolor="yellow", edgecolor="black", zorder=11,
                   label=f"x0={tuple(args.x0)}")
        ax.scatter([R1_TRUE], [R2_TRUE], marker="*", s=400, color="red",
                   edgecolor="black", zorder=12, label="TRUE")
        ax.scatter([best["r1"]], [best["r2"]], marker="X", s=160,
                   color="orange", zorder=10, label="best by cost")
        for g, (m1, m2) in enumerate(gen_means):
            ax.annotate(f"g{g}", (m1, m2), fontsize=9,
                        xytext=(4, 4), textcoords="offset points")
        plt.colorbar(sc, ax=ax, label="log₁₀ |cost|")
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_title(f"{args.variant} — best ({best['r1']:.2f}, "
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
            ax.axvline(g * POP + 0.5, color="lightgray", lw=0.5)
        ax.set_xlabel("eval id"); ax.set_ylabel("dist to TRUE")
        ax.set_yscale("log")
        ax.set_title(f"Convergence (first dist<1: "
                     f"{first_close if first_close else '> budget'})")
        ax.legend(); ax.grid(alpha=0.3, which="both")
        plt.tight_layout()
        png = os.path.join(OUT, f"{tag}.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, f"{tag}_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    with open(os.path.join(OUT, f"{tag}.log"), "w") as f:
        f.write(log_buf.getvalue())
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
