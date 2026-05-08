#!/usr/bin/env python3
"""
cost_selector.py — principled cost-function selection for the
coupling→geometry comparison workflow.

GOAL
----
Pick the BEST cost function for the production CMA-ES coupling sweep
based on data-driven, multi-axis criteria, on the noise-free *and*
noisy modular testbed.

WHAT IT DOES
------------
For each cost mode in {l2, relative, log, effmass, pmean4} and each
noise level σ ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.10}:

  • run multi-seed CMA-ES on N scenarios (different (τ_truth, τ_start))
  • measure success using **MODULAR-AWARE** distance — converging to a
    T-translate or S-image of truth counts as success (it IS the same
    physical theory).  This was the missing piece in §8.7.4 that
    artificially deflated success rates by ~30 pp.
  • aggregate per-mode statistics across scenarios

COMPOSITE SCORE — the key output
--------------------------------
For each mode we report a single scalar **selection score** in [0, 1]:

   S = w_acc · accuracy
     + w_rob · robustness
     + w_unf · uniformity
     + w_eff · efficiency

where
   accuracy   = noise-free success rate              (does it find truth?)
   robustness = mean success across all noise levels (does noise destroy it?)
   uniformity = 1 − std_dist / mean_dist             (is performance scenario-independent?)
   efficiency = 1 − mean_evals / max_evals           (how cheap is convergence?)

Default weights: 0.40, 0.30, 0.20, 0.10.

OUTPUTS (in results/_native_triage/cost_selector_<tag>/):
  scores.json                    raw + composite scores
  noise_curves.png               success% vs σ per mode
  ranking_bars.png               composite score breakdown (stacked bar)
  recommendation.txt             plain-text "use mode X because…"
"""
from __future__ import annotations

import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft  # noqa: F401
import modular_data as md
from modular_testbed import build_or_load_cache, cmaes

OUT_BASE = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT_BASE, exist_ok=True)


# A diverse panel of (truth, start) τ pairs spanning the upper half plane.
# Selected to exercise: (i) tall vs short Im, (ii) edge vs interior Re,
# (iii) various start→truth distances.  Modular-aware distance handles
# T-translate ambiguity, so "near edge" τ's are now legit test points.
DEFAULT_PANEL = [
    ("4-5-6 prod",    (-0.901, 0.794), (-0.500, 0.866)),
    ("equilateral",   (-0.500, 0.866), ( 0.000, 1.000)),
    ("oblique short", (-0.300, 0.500), (-0.700, 0.700)),
    ("tall Im≫Re",    ( 0.000, 1.300), (-0.300, 0.700)),
    ("flat Im~½",     (-0.500, 0.500), ( 0.500, 0.500)),
    ("near-edge -1",  (-0.950, 0.600), ( 0.000, 0.900)),
    ("near-edge +0",  ( 0.000, 0.900), (-0.950, 0.600)),
    ("interior",      (-0.250, 0.700), (-0.250, 1.100)),
]


# --------------------------------------------------------------------------
# Modular-aware GridEvaluator (overrides the base one to add frozen noise).
# --------------------------------------------------------------------------
class SelectorEvaluator:
    def __init__(self, cache, ref_profile, mode, *,
                 noise_sigma=0.0, noise_seed_base=0,
                 frozen_noise=False, pmean=1):
        self.cache = cache
        self.ref = ref_profile
        self.mode = mode
        self.pmean = pmean
        self.noise_sigma = float(noise_sigma)
        self.noise_seed_base = int(noise_seed_base)
        self.frozen = bool(frozen_noise)
        self.re = np.array(cache["re"])
        self.im = np.array(cache["im"])
        self.calls = 0
        self.history = []

    def __call__(self, x):
        re, im = float(x[0]), float(x[1])
        i = int(np.argmin(np.abs(self.re - re)))
        j = int(np.argmin(np.abs(self.im - im)))
        prof = self.cache["grid"].get((i, j))
        if prof is None:
            self.calls += 1
            self.history.append((re, im, float("inf")))
            return float("inf")
        if self.noise_sigma > 0:
            if self.frozen:
                cell_seed = (self.noise_seed_base
                             ^ (i * 100003 + j * 1009)) & 0xFFFFFFFF
                prof_eff = md.add_noise_frozen(prof, self.noise_sigma,
                                                cell_seed)
            else:
                h = hash((self.mode, self.noise_seed_base,
                          i, j, self.calls)) & 0xFFFFFFFF
                rng = np.random.default_rng(h)
                prof_eff = md.add_noise(prof, self.noise_sigma, rng)
        else:
            prof_eff = prof
        try:
            c, _ = md.threedir_cost(self.ref, prof_eff,
                                    mode=self.mode, pmean=self.pmean)
        except Exception:
            c = float("inf")
        if not np.isfinite(c): c = 1e30
        self.calls += 1
        self.history.append((re, im, c))
        return c


# --------------------------------------------------------------------------
def run_one(cache, ref_profile, mode, *,
            tau_truth, tau_start, noise, frozen_noise,
            seeds, sigma0, max_evals, popsize,
            bounds_lo, bounds_hi):
    pmean = 4 if mode == "pmean4" else 1
    base = "l2" if mode == "pmean4" else mode
    runs = []
    for seed in range(seeds):
        ev = SelectorEvaluator(cache, ref_profile, base,
                                noise_sigma=noise,
                                noise_seed_base=1000 * seed,
                                frozen_noise=frozen_noise, pmean=pmean)
        try:
            xf, hist, nev = cmaes(
                ev, [tau_start.real, tau_start.imag],
                sigma0=sigma0, max_evals=max_evals, popsize=popsize,
                seed=seed, bounds_lo=bounds_lo, bounds_hi=bounds_hi)
        except Exception as e:
            runs.append(dict(seed=seed, crash=str(e)))
            continue
        # Modular-aware distance to truth
        d_mod = md.modular_distance(complex(*xf),
                                     complex(tau_truth.real, tau_truth.imag),
                                     n_T=2, include_S=True)
        d_eucl = float(np.linalg.norm(
            np.array(xf) - np.array([tau_truth.real, tau_truth.imag])))
        runs.append(dict(seed=seed, x_final=list(xf),
                          dist_mod=float(d_mod),
                          dist_eucl=d_eucl,
                          n_evals=nev))
    return runs


def aggregate(runs, max_evals):
    ok = [r for r in runs if "dist_mod" in r]
    if not ok:
        return dict(success=0.0, near=0.0, mean_dist=float("nan"),
                    best_dist=float("nan"), mean_evals=float("nan"))
    dists = [r["dist_mod"] for r in ok]
    return dict(
        success=100.0 * sum(1 for d in dists if d < 0.1) / len(ok),
        near=100.0 * sum(1 for d in dists if d < 0.3) / len(ok),
        mean_dist=float(np.mean(dists)),
        best_dist=float(np.min(dists)),
        mean_evals=float(np.mean([r["n_evals"] for r in ok])),
    )


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
    ap.add_argument("--frozen-noise", action="store_true",
                    help="seed noise per cell only (static landscape)")
    ap.add_argument("--scenarios", type=int, default=len(DEFAULT_PANEL))
    ap.add_argument("--cache", default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--w-acc", type=float, default=0.40)
    ap.add_argument("--w-rob", type=float, default=0.30)
    ap.add_argument("--w-unf", type=float, default=0.20)
    ap.add_argument("--w-eff", type=float, default=0.10)
    args = ap.parse_args()

    suffix = (args.tag if args.tag
              else f"n{args.n}_K{args.K_max}_seeds{args.seeds}"
              + ("_frozen" if args.frozen_noise else "")
              ).replace(".", "p")
    OUT = os.path.join(OUT_BASE, f"cost_selector_{suffix}")
    os.makedirs(OUT, exist_ok=True)

    cre, cim = args.center
    re_lo = cre - args.span_re / 2; re_hi = cre + args.span_re / 2
    im_lo = max(0.05, cim - args.span_im / 2)
    im_hi = cim + args.span_im / 2

    print(f"  grid: Re τ ∈ [{re_lo:.3f}, {re_hi:.3f}], "
          f"Im τ ∈ [{im_lo:.3f}, {im_hi:.3f}], n={args.n}×{args.n}")
    print(f"  K_max={args.K_max}  |ω₁|={args.w1}  σ₀={args.sigma0}  "
          f"max_evals={args.max_evals}  seeds={args.seeds}")
    print(f"  noise levels: {args.noise_levels}  "
          f"frozen={args.frozen_noise}")
    print(f"  output → {os.path.relpath(OUT, HERE)}")

    cache_path = args.cache or os.path.join(
        OUT_BASE, f"modular_cache_n{args.n}_K{args.K_max}_w{int(args.w1)}.pkl")
    print("  building/loading τ-grid cache…")
    cache = build_or_load_cache(cache_path, n=args.n,
                                re_lo=re_lo, re_hi=re_hi,
                                im_lo=im_lo, im_hi=im_hi,
                                K_max=args.K_max, w1_norm=args.w1)
    bounds_lo = np.array([re_lo, im_lo])
    bounds_hi = np.array([re_hi, im_hi])

    panel = DEFAULT_PANEL[:args.scenarios]
    modes = ["l2", "relative", "log", "effmass", "pmean4"]
    print(f"  panel ({len(panel)} scenarios): "
          + ", ".join(p[0] for p in panel))
    print(f"  modes: {modes}")
    print()

    # Pre-build ref profiles per scenario (analytic, fast).
    refs = {label: md.make_threedir_profile(
                complex(*tt), w1_norm=args.w1, K_max=args.K_max)
            for label, tt, _ in panel}

    # ------------------------------------------------------------------
    # Sweep:  results[mode][noise][scenario] = aggregate dict
    # ------------------------------------------------------------------
    results = {m: {sig: {} for sig in args.noise_levels} for m in modes}
    t0 = time.perf_counter()
    n_total = len(modes) * len(args.noise_levels) * len(panel)
    done = 0
    for mode in modes:
        for sig in args.noise_levels:
            for label, tt, ts in panel:
                runs = run_one(cache, refs[label], mode,
                                tau_truth=complex(*tt),
                                tau_start=complex(*ts),
                                noise=sig, frozen_noise=args.frozen_noise,
                                seeds=args.seeds, sigma0=args.sigma0,
                                max_evals=args.max_evals,
                                popsize=args.popsize,
                                bounds_lo=bounds_lo, bounds_hi=bounds_hi)
                results[mode][sig][label] = dict(
                    agg=aggregate(runs, args.max_evals),
                    runs=runs)
                done += 1
                if done % max(1, n_total // 20) == 0:
                    elapsed = time.perf_counter() - t0
                    eta = elapsed / done * (n_total - done)
                    print(f"    {done}/{n_total}  "
                          f"({100 * done / n_total:.0f}%)  "
                          f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s")
    wall = time.perf_counter() - t0
    print(f"\n  sweep done in {wall:.1f}s  "
          f"({wall / (len(modes) * len(args.noise_levels) * len(panel) * args.seeds):.2f}s/run)")

    # ------------------------------------------------------------------
    # Aggregate per (mode, noise) over scenarios.
    # ------------------------------------------------------------------
    summary = {m: {} for m in modes}
    for mode in modes:
        for sig in args.noise_levels:
            succ = [results[mode][sig][l]["agg"]["success"]
                    for l, _, _ in panel]
            near = [results[mode][sig][l]["agg"]["near"]
                    for l, _, _ in panel]
            dists = [results[mode][sig][l]["agg"]["mean_dist"]
                     for l, _, _ in panel]
            evals = [results[mode][sig][l]["agg"]["mean_evals"]
                     for l, _, _ in panel]
            summary[mode][sig] = dict(
                mean_success=float(np.mean(succ)),
                mean_near=float(np.mean(near)),
                mean_dist=float(np.nanmean(dists)),
                std_dist=float(np.nanstd(dists)),
                worst_dist=float(np.nanmax(dists)),
                mean_evals=float(np.nanmean(evals)),
            )

    # Print noise-curve table.
    print()
    print("  SUCCESS% vs NOISE  (modular-aware distance):")
    hdr = f"  {'mode':<10}  " + "  ".join(f"σ={s:>6.3g}"
                                            for s in args.noise_levels)
    print(hdr); print("  " + "-" * len(hdr))
    for mode in modes:
        cells = "  ".join(f"{summary[mode][s]['mean_success']:>7.1f}"
                          for s in args.noise_levels)
        print(f"  {mode:<10}  {cells}")

    print()
    print("  MEAN DIST vs NOISE:")
    print(hdr); print("  " + "-" * len(hdr))
    for mode in modes:
        cells = "  ".join(f"{summary[mode][s]['mean_dist']:>7.3f}"
                          for s in args.noise_levels)
        print(f"  {mode:<10}  {cells}")

    # ------------------------------------------------------------------
    # COMPOSITE SCORE
    # ------------------------------------------------------------------
    sig0 = args.noise_levels[0]
    composite = {}
    for mode in modes:
        # accuracy: noise-free success / 100
        accuracy = summary[mode][sig0]["mean_success"] / 100.0
        # robustness: mean success across all noise levels / 100
        robustness = float(np.mean([summary[mode][s]["mean_success"]
                                     for s in args.noise_levels])) / 100.0
        # uniformity: 1 - normalized scenario-std at the noise-free level
        md0 = summary[mode][sig0]["mean_dist"]
        sd0 = summary[mode][sig0]["std_dist"]
        if md0 and md0 > 0:
            uniformity = max(0.0, 1.0 - sd0 / max(md0, 1e-6))
        else:
            uniformity = 1.0
        # efficiency: 1 - mean_evals / max_evals (averaged over noise)
        eff_vals = [1.0 - summary[mode][s]["mean_evals"] / args.max_evals
                    for s in args.noise_levels]
        efficiency = float(np.mean(eff_vals))
        score = (args.w_acc * accuracy + args.w_rob * robustness
                 + args.w_unf * uniformity + args.w_eff * efficiency)
        composite[mode] = dict(
            accuracy=accuracy, robustness=robustness,
            uniformity=uniformity, efficiency=efficiency,
            score=score,
        )

    ranked = sorted(composite.items(),
                     key=lambda kv: -kv[1]["score"])
    print()
    print("  COMPOSITE SELECTION SCORE  (weights: "
          f"acc={args.w_acc} rob={args.w_rob} "
          f"unf={args.w_unf} eff={args.w_eff}):")
    hdr = (f"  {'mode':<10}  {'score':>6}  {'acc':>5}  {'rob':>5}  "
           f"{'unf':>5}  {'eff':>5}")
    print(hdr); print("  " + "-" * len(hdr))
    for mode, c in ranked:
        print(f"  {mode:<10}  {c['score']:>6.3f}  "
              f"{c['accuracy']:>5.3f}  {c['robustness']:>5.3f}  "
              f"{c['uniformity']:>5.3f}  {c['efficiency']:>5.3f}")

    winner_mode, winner_c = ranked[0]
    runner_mode, runner_c = ranked[1]
    margin = winner_c["score"] - runner_c["score"]

    # ------------------------------------------------------------------
    # Recommendation text
    # ------------------------------------------------------------------
    rec_lines = []
    rec_lines.append("=" * 70)
    rec_lines.append("COST-FUNCTION SELECTION REPORT")
    rec_lines.append(f"  panel: {len(panel)} scenarios × {len(args.noise_levels)} noise levels"
                     f" × {args.seeds} seeds = {len(panel)*len(args.noise_levels)*args.seeds} runs/mode")
    rec_lines.append(f"  grid: {args.n}×{args.n}, K_max={args.K_max}, "
                     f"|ω₁|={args.w1}")
    rec_lines.append(f"  noise model: " +
                     ("FROZEN per-cell (static landscape)"
                      if args.frozen_noise else
                      "STOCHASTIC per-call (re-sampled)"))
    rec_lines.append(f"  metric: MODULAR-AWARE distance (T,S equivalences)")
    rec_lines.append("=" * 70)
    rec_lines.append("")
    rec_lines.append(f"WINNER → {winner_mode}   "
                     f"(score = {winner_c['score']:.3f}, "
                     f"margin over #2 ({runner_mode}) = {margin:+.3f})")
    rec_lines.append("")
    rec_lines.append("Per-axis breakdown of the winner:")
    rec_lines.append(f"  accuracy   (noise-free success%) = "
                     f"{winner_c['accuracy']*100:5.1f}%")
    rec_lines.append(f"  robustness (mean over all σ)     = "
                     f"{winner_c['robustness']*100:5.1f}%")
    rec_lines.append(f"  uniformity (cross-scenario)      = "
                     f"{winner_c['uniformity']:5.3f}")
    rec_lines.append(f"  efficiency (1 − evals/budget)    = "
                     f"{winner_c['efficiency']:5.3f}")
    rec_lines.append("")
    rec_lines.append("Noise-curve summary (success% per σ, modular-aware):")
    rec_lines.append("  " + " " * 12 +
                     "  ".join(f"σ={s:<6.3g}" for s in args.noise_levels))
    for mode, _ in ranked:
        rec_lines.append(f"  {mode:<10}  " +
                          "  ".join(f"{summary[mode][s]['mean_success']:>6.1f} "
                                     for s in args.noise_levels))
    rec_lines.append("")
    # Auto-generated rationale
    rec_lines.append("Why this winner:")
    if winner_c["accuracy"] > 0.85:
        rec_lines.append(f"  • dominant noise-free accuracy ({winner_c['accuracy']*100:.0f}%) "
                          "→ correct cost-surface basin")
    if winner_c["robustness"] > 0.5:
        rec_lines.append(f"  • holds up under noise (mean robustness {winner_c['robustness']*100:.0f}%)")
    if winner_c["uniformity"] > 0.5:
        rec_lines.append(f"  • consistent across τ-space (uniformity {winner_c['uniformity']:.2f})")
    # Risk flags
    rec_lines.append("")
    rec_lines.append("Caveats / runners-up:")
    for mode, c in ranked[1:]:
        gap = winner_c["score"] - c["score"]
        notes = []
        if c["accuracy"] >= winner_c["accuracy"] - 0.05:
            notes.append("equal noise-free accuracy")
        if c["robustness"] > winner_c["robustness"]:
            notes.append("BETTER noise robustness — consider for noisy regimes")
        if c["uniformity"] > winner_c["uniformity"]:
            notes.append("more uniform across scenarios")
        rec_lines.append(f"  {mode:<10}  Δscore={-gap:+.3f}  "
                          + ("; ".join(notes) if notes else "no special advantage"))
    rec_lines.append("")
    # Per-noise winner — sometimes different
    rec_lines.append("Per-noise-level winner (success%):")
    for sig in args.noise_levels:
        best = max(modes, key=lambda m: summary[m][sig]["mean_success"])
        rec_lines.append(f"  σ={sig:<6.3g}  →  {best:<10}  "
                          f"({summary[best][sig]['mean_success']:.1f}%)")
    rec_lines.append("")
    rec_lines.append("=" * 70)

    rec = "\n".join(rec_lines)
    rec_path = os.path.join(OUT, "recommendation.txt")
    with open(rec_path, "w") as f: f.write(rec)
    print()
    print(rec)
    print(f"\n  recommendation → {os.path.relpath(rec_path, HERE)}")

    # ------------------------------------------------------------------
    # Save scores
    # ------------------------------------------------------------------
    json_path = os.path.join(OUT, "scores.json")
    with open(json_path, "w") as f:
        json.dump(dict(
            args=vars(args),
            panel=[dict(label=l, truth=tt, start=ts) for l, tt, ts in panel],
            summary=summary,
            composite=composite,
            ranking=[dict(mode=m, **c) for m, c in ranked],
        ), f, indent=2)
    print(f"  scores → {os.path.relpath(json_path, HERE)}")

    # ------------------------------------------------------------------
    # Visualizers
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Plot 1: noise-robustness curves
        fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
        colors = plt.cm.tab10(np.linspace(0, 1, len(modes)))
        for mode, col in zip(modes, colors):
            ys_succ = [summary[mode][s]["mean_success"]
                       for s in args.noise_levels]
            ys_dist = [summary[mode][s]["mean_dist"]
                       for s in args.noise_levels]
            xs = [max(s, 1e-4) for s in args.noise_levels]   # avoid log(0)
            axs[0].plot(xs, ys_succ, "-o", color=col, label=mode, lw=1.6)
            axs[1].plot(xs, ys_dist, "-o", color=col, label=mode, lw=1.6)
        axs[0].set_xscale("log"); axs[0].set_xlabel("noise σ_rel")
        axs[0].set_ylabel("success% (modular-aware)")
        axs[0].set_title("noise robustness — success rate")
        axs[0].grid(True, alpha=0.3); axs[0].legend(fontsize=9)
        axs[1].set_xscale("log"); axs[1].set_xlabel("noise σ_rel")
        axs[1].set_ylabel("mean modular dist → truth")
        axs[1].set_yscale("log")
        axs[1].set_title("noise robustness — mean distance")
        axs[1].grid(True, which="both", alpha=0.3); axs[1].legend(fontsize=9)
        fig.suptitle(f"cost_selector — n={args.n} K_max={args.K_max} "
                     f"seeds={args.seeds}  ({'frozen' if args.frozen_noise else 'stochastic'} noise)",
                     fontsize=11)
        plt.tight_layout()
        png = os.path.join(OUT, "noise_curves.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  noise curves → {os.path.relpath(png, HERE)}")

        # Plot 2: composite-score breakdown stacked bars + ranking
        fig, ax = plt.subplots(figsize=(8, 5))
        names = [m for m, _ in ranked]
        accs = [composite[m]["accuracy"]   * args.w_acc for m in names]
        robs = [composite[m]["robustness"] * args.w_rob for m in names]
        unfs = [composite[m]["uniformity"] * args.w_unf for m in names]
        effs = [composite[m]["efficiency"] * args.w_eff for m in names]
        x = np.arange(len(names))
        ax.bar(x, accs, label=f"accuracy   ×{args.w_acc}", color="#2c7bb6")
        ax.bar(x, robs, bottom=accs, label=f"robustness ×{args.w_rob}",
               color="#abd9e9")
        bot2 = [a + r for a, r in zip(accs, robs)]
        ax.bar(x, unfs, bottom=bot2, label=f"uniformity ×{args.w_unf}",
               color="#fdae61")
        bot3 = [b + u for b, u in zip(bot2, unfs)]
        ax.bar(x, effs, bottom=bot3, label=f"efficiency ×{args.w_eff}",
               color="#d7191c")
        for i, m in enumerate(names):
            tot = composite[m]["score"]
            ax.text(i, tot + 0.01, f"{tot:.3f}", ha="center", fontsize=9,
                    fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels(names, fontsize=10)
        ax.set_ylabel("composite score")
        ax.set_title(f"cost-function ranking  —  WINNER: {winner_mode}")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        png = os.path.join(OUT, "ranking_bars.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  ranking bars → {os.path.relpath(png, HERE)}")

        # Plot 3: per-mode per-noise success heatmap
        fig, ax = plt.subplots(figsize=(1.0 + 1.4 * len(args.noise_levels),
                                         0.7 + 0.6 * len(modes)))
        H = np.array([[summary[m][s]["mean_success"]
                       for s in args.noise_levels] for m in modes])
        im = ax.imshow(H, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100)
        ax.set_xticks(range(len(args.noise_levels)))
        ax.set_xticklabels([f"{s:g}" for s in args.noise_levels], fontsize=9)
        ax.set_yticks(range(len(modes)))
        ax.set_yticklabels(modes, fontsize=10)
        ax.set_xlabel("noise σ_rel"); ax.set_title("success% (modular-aware)")
        for i in range(len(modes)):
            for j in range(len(args.noise_levels)):
                ax.text(j, i, f"{H[i, j]:.0f}", ha="center", va="center",
                        fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.04)
        plt.tight_layout()
        png = os.path.join(OUT, "success_heatmap.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  success heatmap → {os.path.relpath(png, HERE)}")
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"  plotting failed: {e}")


if __name__ == "__main__":
    main()
