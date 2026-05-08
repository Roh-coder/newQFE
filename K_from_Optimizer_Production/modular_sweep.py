#!/usr/bin/env python3
"""
modular_sweep.py — multi-truth × multi-start CMA-ES proxy sweep on a
fine (default 100×100) τ-grid of analytic Ising-CFT three-direction
profiles, with a cost-function comparison visualizer.

For each cost mode in {l2, relative, log, effmass, pmean4} we:

  1. Build the analytic three-direction profile cache once.
  2. For each (τ_truth, τ_start) scenario, run multi-seed CMA-ES.
  3. Aggregate success%, dist→truth, evals → produce per-cost heatmaps,
     a leaderboard, and a textual pros/cons summary.

Usage:
    python modular_sweep.py --n 100 --scenarios 6 --seeds 5
    python modular_sweep.py --n 100 --scenarios 6 --seeds 5 --noise 0.02

OUTPUTS (in results/_native_triage/):
    modular_sweep_<tag>.json            full data
    modular_sweep_<tag>_overview.png    leaderboard heatmap (cost × scenario)
    modular_sweep_<tag>_traj_<mode>.png trajectory grid for one cost
    modular_sweep_<tag>_proscons.txt    plain-text analysis
"""
from __future__ import annotations

import argparse, json, math, os, pickle, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft
import modular_data as md
from modular_testbed import (
    build_or_load_cache, GridEvaluator, cmaes,
)

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------
# Default scenarios — reasonable τ pairs across the upper half plane
# ----------------------------------------------------------------------
DEFAULT_SCENARIOS = [
    # (label,            tau_truth,             tau_start)
    ("4-5-6  prod",      (-0.901, 0.794),       (-0.500, 0.866)),
    ("equilateral",      (-0.500, 0.866),       (+0.000, 1.000)),
    ("oblique short",    (-0.300, 0.500),       (-0.700, 0.700)),
    ("tall  Im≫Re",      ( 0.000, 1.300),       (-0.300, 0.700)),
    ("flat  Im~½",       (-0.500, 0.500),       (+0.500, 0.500)),
    ("near-edge -1",     (-0.950, 0.600),       ( 0.000, 0.900)),
    ("near-edge +0",     ( 0.000, 0.900),       (-0.950, 0.600)),
    ("interior shift",   (-0.250, 0.700),       (-0.250, 1.100)),
]


# ----------------------------------------------------------------------
def run_scenario(cache, ref_profile, mode, *,
                 tau_truth, tau_start,
                 noise, seeds, sigma0, max_evals, popsize,
                 bounds_lo, bounds_hi):
    pmean = 4 if mode == "pmean4" else 1
    base_mode = "l2" if mode == "pmean4" else mode
    runs = []
    for seed in range(seeds):
        ev = GridEvaluator(cache, ref_profile, base_mode,
                           noise_sigma=noise, noise_seed=1000 * seed,
                           pmean=pmean)
        try:
            xf, hist, nev = cmaes(
                ev, [tau_start.real, tau_start.imag],
                sigma0=sigma0, max_evals=max_evals, popsize=popsize,
                seed=seed, bounds_lo=bounds_lo, bounds_hi=bounds_hi)
        except Exception as e:
            runs.append(dict(seed=seed, crash=str(e)))
            continue
        x_truth = np.array([tau_truth.real, tau_truth.imag])
        d = float(np.linalg.norm(xf - x_truth))
        status = ("found" if d < 0.1 else "near" if d < 0.3 else "lost")
        runs.append(dict(seed=seed, x_final=xf.tolist(),
                         dist_truth=d, n_evals=nev,
                         status=status, history=hist))
    return runs


def aggregate(runs):
    ok = [r for r in runs if "dist_truth" in r]
    if not ok:
        return dict(success=0, near=0, mean_dist=float("nan"),
                    best_dist=float("nan"), mean_evals=float("nan"))
    dists = [r["dist_truth"] for r in ok]
    return dict(
        success=100 * sum(1 for r in ok if r["status"] == "found") / len(ok),
        near=100 * sum(1 for r in ok if r["status"] in ("found", "near")) / len(ok),
        mean_dist=float(np.mean(dists)),
        best_dist=float(np.min(dists)),
        mean_evals=float(np.mean([r["n_evals"] for r in ok])),
    )


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--K-max", type=int, default=8)
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.4)
    ap.add_argument("--center", type=float, nargs=2, default=[-0.45, 0.85],
                    help="grid centre (Re τ, Im τ)")
    ap.add_argument("--scenarios", type=int, default=len(DEFAULT_SCENARIOS),
                    help="how many of DEFAULT_SCENARIOS to use")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--sigma0", type=float, default=0.25)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--cache", default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--w1", type=float, default=16.0,
                    help="|ω₁| in lattice units")
    args = ap.parse_args()

    cre, cim = args.center
    re_lo = cre - args.span_re / 2; re_hi = cre + args.span_re / 2
    im_lo = max(0.05, cim - args.span_im / 2)
    im_hi = cim + args.span_im / 2
    print(f"  grid: Re τ ∈ [{re_lo:.3f}, {re_hi:.3f}], "
          f"Im τ ∈ [{im_lo:.3f}, {im_hi:.3f}],  n={args.n}×{args.n}")
    print(f"  K_max = {args.K_max},   |ω₁| = {args.w1}")
    print(f"  noise σ_rel = {args.noise:.3g}")

    cache_path = args.cache or os.path.join(
        OUT, f"modular_cache_n{args.n}_K{args.K_max}_w{int(args.w1)}.pkl")
    print("  building/loading τ-grid cache…")
    cache = build_or_load_cache(cache_path, n=args.n,
                                re_lo=re_lo, re_hi=re_hi,
                                im_lo=im_lo, im_hi=im_hi,
                                K_max=args.K_max, w1_norm=args.w1)
    bounds_lo = np.array([re_lo, im_lo])
    bounds_hi = np.array([re_hi, im_hi])

    scenarios = DEFAULT_SCENARIOS[:args.scenarios]
    modes = ["l2", "relative", "log", "effmass", "pmean4"]

    # ------------------------------------------------------------------
    # Pre-build ref profiles per scenario  (analytic, fast)
    # ------------------------------------------------------------------
    ref_profiles = {}
    for label, t_truth, _ in scenarios:
        tau = complex(*t_truth)
        ref_profiles[label] = md.make_threedir_profile(
            tau, w1_norm=args.w1, K_max=args.K_max)

    # ------------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------------
    results = {m: {} for m in modes}
    print()
    hdr = (f"  {'scenario':<18}  {'mode':<10}  {'success%':>8}  "
           f"{'near%':>6}  {'mean_d':>7}  {'best_d':>7}  {'evals':>6}")
    print(hdr); print("  " + "-" * len(hdr))
    t0 = time.perf_counter()
    for label, t_truth, t_start in scenarios:
        tau_truth = complex(*t_truth); tau_start = complex(*t_start)
        ref_prof = ref_profiles[label]
        for mode in modes:
            runs = run_scenario(cache, ref_prof, mode,
                                tau_truth=tau_truth, tau_start=tau_start,
                                noise=args.noise, seeds=args.seeds,
                                sigma0=args.sigma0,
                                max_evals=args.max_evals,
                                popsize=args.popsize,
                                bounds_lo=bounds_lo, bounds_hi=bounds_hi)
            agg = aggregate(runs)
            results[mode][label] = dict(agg=agg, runs=runs,
                                        tau_truth=t_truth,
                                        tau_start=t_start)
            print(f"  {label:<18}  {mode:<10}  "
                  f"{agg['success']:>7.0f}%  {agg['near']:>5.0f}%  "
                  f"{agg['mean_dist']:>7.3f}  {agg['best_dist']:>7.3f}  "
                  f"{agg['mean_evals']:>6.1f}")
    wall = time.perf_counter() - t0
    print(f"\n  sweep done in {wall:.1f}s")

    # ------------------------------------------------------------------
    # Aggregate scoring per cost mode (across scenarios)
    # ------------------------------------------------------------------
    print()
    print("  COST-MODE LEADERBOARD (averaged over scenarios):")
    print(f"  {'mode':<10}  {'mean succ%':>10}  {'mean near%':>10}  "
          f"{'mean dist':>10}  {'worst dist':>10}  {'std dist':>9}")
    print("  " + "-" * 64)
    leader = []
    for mode in modes:
        succs = [results[mode][s[0]]["agg"]["success"] for s in scenarios]
        nears = [results[mode][s[0]]["agg"]["near"] for s in scenarios]
        dists = [results[mode][s[0]]["agg"]["mean_dist"] for s in scenarios]
        msucc = float(np.mean(succs))
        mnear = float(np.mean(nears))
        mdist = float(np.nanmean(dists))
        wdist = float(np.nanmax(dists))
        sdist = float(np.nanstd(dists))
        print(f"  {mode:<10}  {msucc:>9.0f}%  {mnear:>9.0f}%  "
              f"{mdist:>10.3f}  {wdist:>10.3f}  {sdist:>9.3f}")
        leader.append((mode, msucc, mnear, mdist, wdist, sdist))
    leader.sort(key=lambda r: (-r[1], -r[2], r[3]))
    print()
    print("  RANKED  (success% desc, near% desc, mean dist asc):")
    for i, (m, ms, mn, md_, wd, sd) in enumerate(leader, 1):
        print(f"    {i}. {m:<10}  succ={ms:>3.0f}%  near={mn:>3.0f}%  "
              f"mean_d={md_:.3f}  worst_d={wd:.3f}  σ_d={sd:.3f}")

    # ------------------------------------------------------------------
    # Pros/cons heuristic
    # ------------------------------------------------------------------
    proscons = []
    proscons.append("=" * 70)
    proscons.append("COST-FUNCTION PROS / CONS  (data-driven)")
    proscons.append(f"  grid {args.n}×{args.n}, K_max={args.K_max},"
                    f" |ω₁|={args.w1}, σ_rel={args.noise:.3g},"
                    f" {args.seeds} seeds, {len(scenarios)} scenarios")
    proscons.append("=" * 70)
    for mode, ms, mn, md_, wd, sd in leader:
        succs = [results[mode][s[0]]["agg"]["success"] for s in scenarios]
        nears = [results[mode][s[0]]["agg"]["near"] for s in scenarios]
        dists = [results[mode][s[0]]["agg"]["mean_dist"] for s in scenarios]
        worst = max(zip(dists, [s[0] for s in scenarios]))
        best  = min(zip(dists, [s[0] for s in scenarios]))
        line = []
        line.append(f"\n[{mode}]   success={ms:.0f}%  near={mn:.0f}%  "
                    f"mean_d={md_:.3f}  worst_d={wd:.3f}  σ_d={sd:.3f}")
        # PROS
        pros = []
        if ms >= 80:
            pros.append("very high success rate across scenarios")
        elif ms >= 50:
            pros.append("solid majority-of-scenarios convergence")
        if sd < 0.1:
            pros.append("uniform behaviour across τ-space (low variance)")
        if md_ < 0.1:
            pros.append("typically lands inside the truth basin")
        if mode == "log":
            pros.append("scale-invariant under any overall G-amplitude shift")
        if mode == "effmass":
            pros.append("amplitude-free; targets effective mass directly")
        if mode == "relative":
            pros.append("dimensionless; cancels per-direction magnitude")
        if mode == "pmean4":
            pros.append("anisotropy-penalising: a single bad direction "
                        "raises the cost")
        if mode == "l2":
            pros.append("simplest kernel; no log/division → robust to "
                        "near-zero G")
        # CONS
        cons = []
        if ms < 50:
            cons.append(f"convergence < 50% on average  ({ms:.0f}%)")
        if wd > 0.5:
            cons.append(f"WORST scenario: '{worst[1]}' "
                        f"with mean_d={worst[0]:.3f}")
        if sd > 0.2:
            cons.append("high variance across scenarios (cost is "
                        "scenario-sensitive)")
        if mode in ("log", "effmass") and args.noise > 0:
            cons.append("log/effmass blow up when noise drives G≤0; "
                        "fragile under noise")
        if mode == "pmean4":
            cons.append("p=4 aggregation amplifies noise on the "
                        "worst-fitting direction")
        if mode == "l2":
            cons.append("not scale-invariant; absolute G magnitude "
                        "differences leak into cost")
        if mode == "relative":
            cons.append("denominator G_r near 0 can spike the cost; "
                        "regularised here with +1e-30")
        proscons.append("  " + "\n  ".join(line + ["  + " + p for p in pros]
                                            + ["  − " + c for c in cons]))
    proscons.append("\n" + "=" * 70)

    # ------------------------------------------------------------------
    # Save artifacts
    # ------------------------------------------------------------------
    suffix = (args.tag if args.tag
              else f"n{args.n}_K{args.K_max}_noise{args.noise:.2g}"
              ).replace(".", "p").replace("=", "")
    json_path = os.path.join(OUT, f"modular_sweep_{suffix}.json")
    serialise = {m: {s: {"agg": v["agg"],
                         "tau_truth": v["tau_truth"],
                         "tau_start": v["tau_start"],
                         "runs": [{kk: vv for kk, vv in r.items()
                                   if kk != "history"} for r in v["runs"]]}
                     for s, v in d.items()}
                 for m, d in results.items()}
    with open(json_path, "w") as f:
        json.dump(dict(args=vars(args),
                       scenarios=[dict(label=l, truth=tt, start=ts)
                                  for l, tt, ts in scenarios],
                       leaderboard=[dict(mode=m, mean_success=ms,
                                         mean_near=mn, mean_dist=md_,
                                         worst_dist=wd, std_dist=sd)
                                    for m, ms, mn, md_, wd, sd in leader],
                       results=serialise), f, indent=2)
    print(f"\n  summary → {os.path.relpath(json_path, HERE)}")

    proscons_path = os.path.join(OUT, f"modular_sweep_{suffix}_proscons.txt")
    with open(proscons_path, "w") as f:
        f.write("\n".join(proscons))
    print(f"  pros/cons → {os.path.relpath(proscons_path, HERE)}")
    print()
    print("\n".join(proscons))

    # ------------------------------------------------------------------
    # Visualizer 1: leaderboard heatmap (mode × scenario, 2 panels)
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc_labels = [s[0] for s in scenarios]
        success_grid = np.array([[results[m][s]["agg"]["success"]
                                  for s in sc_labels] for m in modes])
        dist_grid = np.array([[results[m][s]["agg"]["mean_dist"]
                               for s in sc_labels] for m in modes])

        fig, axs = plt.subplots(1, 2, figsize=(7.0 + 1.5 * len(scenarios),
                                                3.0 + 0.5 * len(modes)))
        for ax, data, title, cmap, vmin, vmax, fmt in [
            (axs[0], success_grid, "success%  (CMA-ES converged < 0.1 to truth)",
             "RdYlGn", 0, 100, "{:.0f}"),
            (axs[1], np.log10(np.clip(dist_grid, 1e-3, None)),
             "log₁₀(mean dist → truth)", "RdYlGn_r", -2.5, 0.0, "{:.2f}"),
        ]:
            im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xticks(range(len(sc_labels)))
            ax.set_xticklabels(sc_labels, rotation=35, ha="right", fontsize=8)
            ax.set_yticks(range(len(modes)))
            ax.set_yticklabels(modes, fontsize=10)
            ax.set_title(title, fontsize=10)
            plt.colorbar(im, ax=ax, fraction=0.04)
            for i in range(len(modes)):
                for j in range(len(sc_labels)):
                    raw = (success_grid if data is success_grid
                           else dist_grid)[i, j]
                    txt = fmt.format(raw if data is success_grid else raw)
                    ax.text(j, i, txt, ha="center", va="center",
                            color="black", fontsize=7)
        fig.suptitle(f"modular_sweep — n={args.n}, K_max={args.K_max}, "
                     f"σ_rel={args.noise:g}, seeds={args.seeds}",
                     fontsize=11)
        plt.tight_layout()
        png = os.path.join(OUT, f"modular_sweep_{suffix}_overview.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  overview → {os.path.relpath(png, HERE)}")

        # --------------------------------------------------------------
        # Visualizer 2: per-mode trajectory grid
        # one column per scenario, one row per mode (full picture)
        # --------------------------------------------------------------
        re_arr = np.array(cache["re"])
        im_arr = np.array(cache["im"])

        # Pre-compute heatmap per (mode, scenario): cost surface depends
        # on ref_profile so it changes per scenario.
        ncols = min(len(scenarios), 4)
        nrows_per_mode = math.ceil(len(scenarios) / ncols)
        for mode in modes:
            base_mode = "l2" if mode == "pmean4" else mode
            pmean = 4 if mode == "pmean4" else 1
            fig, axs = plt.subplots(nrows_per_mode, ncols,
                                    figsize=(4.6 * ncols,
                                             4.0 * nrows_per_mode),
                                    squeeze=False)
            for k, (label, t_truth, t_start) in enumerate(scenarios):
                r, c = k // ncols, k % ncols
                ax = axs[r, c]
                ref_prof = ref_profiles[label]
                grid_cost = np.full((args.n, args.n), np.nan)
                for (i, j), prof in cache["grid"].items():
                    if prof is None: continue
                    try:
                        cv, _ = md.threedir_cost(ref_prof, prof,
                                                 mode=base_mode,
                                                 pmean=pmean)
                    except Exception:
                        cv = float("nan")
                    if np.isfinite(cv) and cv > 0:
                        grid_cost[i, j] = cv
                with np.errstate(invalid="ignore", divide="ignore"):
                    logg = np.log10(grid_cost)
                valid = logg[np.isfinite(logg)]
                vmin, vmax = (np.percentile(valid, [2, 98])
                              if len(valid) else (-1, 1))
                im = ax.imshow(logg.T, origin="lower",
                               extent=[re_arr[0], re_arr[-1],
                                       im_arr[0], im_arr[-1]],
                               aspect="auto", cmap="viridis",
                               vmin=vmin, vmax=vmax, alpha=0.85)
                tau_truth = complex(*t_truth); tau_start = complex(*t_start)
                ax.scatter([tau_truth.real], [tau_truth.imag],
                           marker="*", s=260, color="red",
                           edgecolor="white", lw=1.4, zorder=10)
                ax.scatter([tau_start.real], [tau_start.imag],
                           marker="o", s=110, color="cyan",
                           edgecolor="black", lw=1.0, zorder=9)
                colors = plt.cm.tab10(np.linspace(0, 1, args.seeds))
                for run, col in zip(results[mode][label]["runs"], colors):
                    if "history" not in run: continue
                    hist = np.array([h[:2] for h in run["history"]])
                    if len(hist):
                        ax.plot(hist[:, 0], hist[:, 1], "-", color=col,
                                lw=1.0, alpha=0.85)
                        ax.scatter([hist[-1, 0]], [hist[-1, 1]],
                                   marker="X", s=60, color=col,
                                   edgecolor="black", lw=0.5, zorder=8)
                agg = results[mode][label]["agg"]
                ax.set_xlim(re_lo, re_hi); ax.set_ylim(im_lo, im_hi)
                ax.set_title(f"{label}\n{mode}: succ={agg['success']:.0f}%  "
                             f"d={agg['mean_dist']:.2f}", fontsize=8)
                ax.set_xlabel("Re τ", fontsize=8)
                ax.set_ylabel("Im τ", fontsize=8)
            for k in range(len(scenarios), nrows_per_mode * ncols):
                axs[k // ncols, k % ncols].axis("off")
            fig.suptitle(f"modular_sweep — cost = {mode}  "
                         f"(σ_rel={args.noise:g})", fontsize=11)
            plt.tight_layout()
            png = os.path.join(OUT,
                               f"modular_sweep_{suffix}_traj_{mode}.png")
            plt.savefig(png, dpi=110); plt.close(fig)
            print(f"  trajectories[{mode}] → {os.path.relpath(png, HERE)}")
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"  plotting failed: {e}")


if __name__ == "__main__":
    main()
