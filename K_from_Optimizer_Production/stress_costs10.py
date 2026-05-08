#!/usr/bin/env python3
"""
stress_costs10.py — stress test 10 NEW cost functions (beyond the
original l2/relative/log/effmass/pmean4 quintet) under the same
modular-aware noise sweep used for cost_selector.

Cost functions tested (registered in modular_data._PER_DIR):
  1.  l1          mean |Δ|                         L1 robust
  2.  l1_log      mean |Δlog G|                    L1 on log
  3.  huber       Huber on raw residuals           soft outlier clip
  4.  huber_log   Huber on log residuals           soft outlier clip on log
  5.  chi2        Σ Δ² / |G_r|                     Poisson-like weighting
  6.  cosine      1 - <G_t,G_r>/(|G_t||G_r|)       scale-invariant shape
  7.  slope       (slope_t - slope_r)²             global eff-mass via fit
  8.  trim_rel    relative L2 with top-20% trimmed outlier-resistant
  9.  median_log  median Δlog² G                   very robust
  10. rel_l1      mean |Δ/G_r|                     relative L1

For each, we run:
  • CMA-ES with the modular testbed
  • 6 noise levels σ ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.10}
  • 8 scenarios × 5 seeds  =  40 runs per (cost, σ)  =  240 runs/cost

Composite score (same formula as cost_selector):
   S = 0.40·acc + 0.30·rob + 0.20·unf + 0.10·eff

Outputs:  results/_native_triage/stress_costs10_<tag>/
   noise_curves.png    success% vs σ for all 10 + the 5 baselines
   ranking_bars.png    composite score breakdown
   success_heatmap.png 15 × 6 grid
   recommendation.txt  text verdict
   scores.json         raw numbers
"""
from __future__ import annotations

import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import modular_data as md
from modular_testbed import build_or_load_cache, cmaes
from cost_selector import (SelectorEvaluator, DEFAULT_PANEL,
                            run_one as _selector_run_one,
                            aggregate as _selector_agg)

OUT_BASE = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT_BASE, exist_ok=True)

NEW_MODES = ["l1", "l1_log", "huber", "huber_log", "chi2",
             "cosine", "slope", "trim_rel", "median_log", "rel_l1"]
BASELINES = ["l2", "relative", "log", "effmass"]
ALL_MODES = NEW_MODES + BASELINES   # 14 total


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--K-max", type=int, default=8)
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.4)
    ap.add_argument("--center", type=float, nargs=2, default=[-0.45, 0.85])
    ap.add_argument("--w1", type=float, default=16.0)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--sigma0", type=float, default=0.35)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--noise-levels", type=float, nargs="+",
                    default=[0.0, 0.005, 0.01, 0.02, 0.05, 0.10])
    ap.add_argument("--scenarios", type=int, default=len(DEFAULT_PANEL))
    ap.add_argument("--cache", default=None)
    ap.add_argument("--tag", default="prod")
    ap.add_argument("--baselines", action="store_true",
                    help="also run l2/relative/log/effmass for context")
    args = ap.parse_args()

    OUT = os.path.join(OUT_BASE, f"stress_costs10_{args.tag}")
    os.makedirs(OUT, exist_ok=True)

    cre, cim = args.center
    re_lo = cre - args.span_re / 2; re_hi = cre + args.span_re / 2
    im_lo = max(0.05, cim - args.span_im / 2)
    im_hi = cim + args.span_im / 2
    bounds_lo = np.array([re_lo, im_lo]); bounds_hi = np.array([re_hi, im_hi])

    panel = DEFAULT_PANEL[:args.scenarios]
    modes = list(NEW_MODES)
    if args.baselines: modes += BASELINES
    print(f"  testing {len(modes)} cost modes: {modes}")
    print(f"  panel ({len(panel)}): " + ", ".join(p[0] for p in panel))
    print(f"  noise: {args.noise_levels}  seeds: {args.seeds}")

    cache_path = args.cache or os.path.join(
        OUT_BASE, f"modular_cache_n{args.n}_K{args.K_max}_w{int(args.w1)}.pkl")
    cache = build_or_load_cache(cache_path, n=args.n,
                                 re_lo=re_lo, re_hi=re_hi,
                                 im_lo=im_lo, im_hi=im_hi,
                                 K_max=args.K_max, w1_norm=args.w1)
    refs = {label: md.make_threedir_profile(complex(*tt),
                                              w1_norm=args.w1,
                                              K_max=args.K_max)
            for label, tt, _ in panel}

    # --------------------------------------------------------------
    results = {m: {sig: {} for sig in args.noise_levels} for m in modes}
    n_total = len(modes) * len(args.noise_levels) * len(panel)
    done = 0; t0 = time.perf_counter()
    for mode in modes:
        for sig in args.noise_levels:
            for label, tt, ts in panel:
                runs = _selector_run_one(cache, refs[label], mode,
                                          tau_truth=complex(*tt),
                                          tau_start=complex(*ts),
                                          noise=sig, frozen_noise=False,
                                          seeds=args.seeds,
                                          sigma0=args.sigma0,
                                          max_evals=args.max_evals,
                                          popsize=args.popsize,
                                          bounds_lo=bounds_lo,
                                          bounds_hi=bounds_hi)
                results[mode][sig][label] = dict(
                    agg=_selector_agg(runs, args.max_evals), runs=runs)
                done += 1
                if done % max(1, n_total // 20) == 0:
                    el = time.perf_counter() - t0
                    eta = el / done * (n_total - done)
                    print(f"    {done}/{n_total} ({100*done/n_total:.0f}%) "
                          f"el={el:.0f}s ETA={eta:.0f}s")
    wall = time.perf_counter() - t0
    print(f"\n  sweep done in {wall:.1f}s "
          f"({wall / (n_total * args.seeds):.2f}s/run)")

    # --------------------------------------------------------------
    summary = {m: {} for m in modes}
    for mode in modes:
        for sig in args.noise_levels:
            succ = [results[mode][sig][l]["agg"]["success"]
                    for l, _, _ in panel]
            dists = [results[mode][sig][l]["agg"]["mean_dist"]
                     for l, _, _ in panel]
            evals = [results[mode][sig][l]["agg"]["mean_evals"]
                     for l, _, _ in panel]
            summary[mode][sig] = dict(
                mean_success=float(np.mean(succ)),
                mean_dist=float(np.nanmean(dists)),
                std_dist=float(np.nanstd(dists)),
                worst_dist=float(np.nanmax(dists)),
                mean_evals=float(np.nanmean(evals)),
            )

    # composite score
    sig0 = args.noise_levels[0]
    composite = {}
    for mode in modes:
        accuracy = summary[mode][sig0]["mean_success"] / 100.0
        robustness = float(np.mean([summary[mode][s]["mean_success"]
                                     for s in args.noise_levels])) / 100.0
        md0 = summary[mode][sig0]["mean_dist"]
        sd0 = summary[mode][sig0]["std_dist"]
        uniformity = max(0.0, 1.0 - sd0 / max(md0, 1e-6)) if md0 else 1.0
        eff_vals = [1.0 - summary[mode][s]["mean_evals"] / args.max_evals
                    for s in args.noise_levels]
        efficiency = float(np.mean(eff_vals))
        score = (0.40 * accuracy + 0.30 * robustness
                 + 0.20 * uniformity + 0.10 * efficiency)
        composite[mode] = dict(accuracy=accuracy, robustness=robustness,
                                uniformity=uniformity, efficiency=efficiency,
                                score=score)
    ranked = sorted(composite.items(), key=lambda kv: -kv[1]["score"])

    # --------------------------------------------------------------
    print()
    print("  SUCCESS% vs NOISE  (modular-aware):")
    hdr = f"  {'mode':<11}  " + "  ".join(f"σ={s:>6.3g}"
                                            for s in args.noise_levels)
    print(hdr); print("  " + "-" * len(hdr))
    for mode, _ in ranked:
        cells = "  ".join(f"{summary[mode][s]['mean_success']:>7.1f}"
                          for s in args.noise_levels)
        flag = " ← BASELINE" if mode in BASELINES else ""
        print(f"  {mode:<11}  {cells}{flag}")

    print()
    print(f"  COMPOSITE RANKING (acc.40 rob.30 unf.20 eff.10):")
    hdr = (f"  {'mode':<11}  {'score':>6}  {'acc':>5}  {'rob':>5}  "
           f"{'unf':>5}  {'eff':>5}")
    print(hdr); print("  " + "-" * len(hdr))
    for mode, c in ranked:
        flag = " ← BASELINE" if mode in BASELINES else ""
        print(f"  {mode:<11}  {c['score']:>6.3f}  "
              f"{c['accuracy']:>5.3f}  {c['robustness']:>5.3f}  "
              f"{c['uniformity']:>5.3f}  {c['efficiency']:>5.3f}{flag}")

    # --------------------------------------------------------------
    # Recommendation
    rec = []
    rec.append("=" * 72)
    rec.append(f"10-COST STRESS REPORT  ({len(modes)} modes total)")
    rec.append("=" * 72)
    rec.append(f"  panel : {len(panel)} scenarios × "
               f"{len(args.noise_levels)} noise × {args.seeds} seeds = "
               f"{len(panel)*len(args.noise_levels)*args.seeds} runs/mode")
    rec.append(f"  config: n={args.n}  K_max={args.K_max}  |ω₁|={args.w1}")
    rec.append(f"  metric: modular-aware distance (T,S equivalences)")
    rec.append("")
    winner_mode, winner_c = ranked[0]
    rec.append(f"WINNER → {winner_mode}   "
               f"(score = {winner_c['score']:.3f})")
    rec.append("")
    rec.append(f"  acc      = {winner_c['accuracy']*100:5.1f}% (noise-free)")
    rec.append(f"  rob      = {winner_c['robustness']*100:5.1f}% (mean over σ)")
    rec.append(f"  uniform  = {winner_c['uniformity']:.3f}")
    rec.append(f"  efficiency = {winner_c['efficiency']:.3f}")
    rec.append("")
    rec.append("Per-noise-level winner (success%):")
    for sig in args.noise_levels:
        best = max(modes, key=lambda m: summary[m][sig]["mean_success"])
        rec.append(f"  σ={sig:<6.3g}  →  {best:<11}  "
                   f"({summary[best][sig]['mean_success']:.1f}%)")
    rec.append("")
    rec.append("Top 5 ranking:")
    for mode, c in ranked[:5]:
        flag = " (baseline)" if mode in BASELINES else ""
        rec.append(f"  {mode:<11}  score={c['score']:.3f}  "
                   f"acc={c['accuracy']*100:.0f}%  "
                   f"rob={c['robustness']*100:.0f}%{flag}")
    rec.append("")
    rec.append("Bottom 3 (avoid):")
    for mode, c in ranked[-3:]:
        flag = " (baseline)" if mode in BASELINES else ""
        rec.append(f"  {mode:<11}  score={c['score']:.3f}  "
                   f"acc={c['accuracy']*100:.0f}%  "
                   f"rob={c['robustness']*100:.0f}%{flag}")
    rec.append("=" * 72)
    text = "\n".join(rec)
    print(); print(text)
    with open(os.path.join(OUT, "recommendation.txt"), "w") as f:
        f.write(text)
    with open(os.path.join(OUT, "scores.json"), "w") as f:
        json.dump(dict(args=vars(args), summary=summary,
                        composite=composite,
                        ranking=[dict(mode=m, **c) for m, c in ranked]),
                  f, indent=2)

    # --------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # noise curves: top 10 by composite score for legibility
        fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.8))
        top10 = [m for m, _ in ranked[:10]]
        cmap = plt.cm.tab20(np.linspace(0, 1, len(top10)))
        for mode, col in zip(top10, cmap):
            ys = [summary[mode][s]["mean_success"]
                  for s in args.noise_levels]
            ds = [summary[mode][s]["mean_dist"]
                  for s in args.noise_levels]
            xs = [max(s, 1e-4) for s in args.noise_levels]
            ls = "--" if mode in BASELINES else "-"
            lbl = mode + (" *" if mode in BASELINES else "")
            axs[0].semilogx(xs, ys, ls + "o", color=col, label=lbl, lw=1.6)
            axs[1].loglog(xs, ds, ls + "o", color=col, label=lbl, lw=1.6)
        axs[0].set_xlabel("σ_rel"); axs[0].set_ylabel("success%")
        axs[0].set_title("noise robustness — top 10 by composite")
        axs[0].grid(alpha=0.3); axs[0].legend(fontsize=8, ncol=2)
        axs[0].axhline(50, color="r", ls=":", alpha=0.4)
        axs[1].set_xlabel("σ_rel"); axs[1].set_ylabel("mean modular dist")
        axs[1].set_title("distance vs noise (lower = better)")
        axs[1].grid(which="both", alpha=0.3); axs[1].legend(fontsize=8, ncol=2)
        fig.suptitle(f"stress_costs10 · n={args.n} K={args.K_max} "
                     f"seeds={args.seeds} ·  * = baseline", fontsize=11)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "noise_curves.png"), dpi=130)
        plt.close(fig)

        # ranking bars
        fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(modes)), 5))
        names = [m for m, _ in ranked]
        accs = [composite[m]["accuracy"] * 0.40 for m in names]
        robs = [composite[m]["robustness"] * 0.30 for m in names]
        unfs = [composite[m]["uniformity"] * 0.20 for m in names]
        effs = [composite[m]["efficiency"] * 0.10 for m in names]
        x = np.arange(len(names))
        ax.bar(x, accs, label="accuracy ×0.4", color="#2c7bb6")
        ax.bar(x, robs, bottom=accs, label="robustness ×0.3", color="#abd9e9")
        bot2 = [a + r for a, r in zip(accs, robs)]
        ax.bar(x, unfs, bottom=bot2, label="uniformity ×0.2", color="#fdae61")
        bot3 = [b + u for b, u in zip(bot2, unfs)]
        ax.bar(x, effs, bottom=bot3, label="efficiency ×0.1", color="#d7191c")
        for i, m in enumerate(names):
            tot = composite[m]["score"]
            color = "k" if m not in BASELINES else "blue"
            ax.text(i, tot + 0.01, f"{tot:.3f}", ha="center", fontsize=8,
                    fontweight="bold", color=color)
        ax.set_xticks(x); ax.set_xticklabels(
            [m + ("*" if m in BASELINES else "") for m in names],
            rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("composite score")
        ax.set_title(f"cost stress ranking · WINNER: {winner_mode}  "
                      "(blue = baseline)")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "ranking_bars.png"), dpi=130)
        plt.close(fig)

        # heatmap
        fig, ax = plt.subplots(figsize=(1.0 + 1.4 * len(args.noise_levels),
                                         0.7 + 0.4 * len(modes)))
        H = np.array([[summary[m][s]["mean_success"]
                       for s in args.noise_levels]
                      for m, _ in ranked])
        im = ax.imshow(H, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100)
        ax.set_xticks(range(len(args.noise_levels)))
        ax.set_xticklabels([f"{s:g}" for s in args.noise_levels], fontsize=9)
        ax.set_yticks(range(len(modes)))
        ax.set_yticklabels([m + ("*" if m in BASELINES else "")
                            for m, _ in ranked], fontsize=9)
        ax.set_xlabel("σ_rel"); ax.set_title("success% (modular-aware)")
        for i in range(len(modes)):
            for j in range(len(args.noise_levels)):
                ax.text(j, i, f"{H[i, j]:.0f}", ha="center", va="center",
                        fontsize=7)
        plt.colorbar(im, ax=ax, fraction=0.04)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "success_heatmap.png"), dpi=130)
        plt.close(fig)
        print(f"\n  artifacts → {os.path.relpath(OUT, HERE)}")
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"  plot failed: {e}")


if __name__ == "__main__":
    main()
