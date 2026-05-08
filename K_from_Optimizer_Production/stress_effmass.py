#!/usr/bin/env python3
"""
stress_effmass.py — focused stress test of the `effmass` cost.

Why a separate harness?
-----------------------
The cost_selector run flagged `effmass` as the best **noise-free** cost
(77.5% success) but with a steep noise cliff (77.5 → 50 → 30 → 2.5%
across σ = 0, 0.005, 0.01, 0.02).  Effective mass is m_k = −log(G_{k+1}/G_k):
the log-ratio amplifies multiplicative noise dramatically as G shrinks
with k.  This script characterises that breakdown along five axes:

  AXIS A — fine-grain noise sweep  σ ∈ {0, 0.001, 0.002, 0.005, 0.008,
                                          0.01, 0.015, 0.02, 0.03, 0.05}
           → find the **σ_50** breakpoint where success drops below 50%.
  AXIS B — K_max sensitivity         K ∈ {3, 4, 6, 8, 10, 12}
           → does longer m_eff window help (more data) or hurt
             (deeper into noisy tails)?
  AXIS C — start-distance scaling   |τ_start − τ_truth| ∈ {0.1, 0.3,
                                                            0.5, 0.8, 1.2}
           → does a wider basin survive noise better?
  AXIS D — seed stability           20 seeds per (mode, cell)
           → variance of the converged τ; outlier failure modes
  AXIS E — frozen vs stochastic     compare per-cell-frozen vs per-call
                                     re-sampled noise on identical sweeps.

Outputs in results/_native_triage/stress_effmass_<tag>/:
  noise_fine.png           success%, mean dist, P95 dist vs σ (fine grid)
  Kmax_scan.png            success% vs K_max at three σ levels
  start_distance.png       success% vs |start - truth| at three σ levels
  seed_stability.png       per-cell scatter of converged τ (20 seeds)
  frozen_vs_stochastic.png overlay
  failure_atlas.png        where in τ-space do failed runs land?
  stress_report.txt        text summary with σ_50 and dominant failure modes
  stress.json              raw numbers
"""
from __future__ import annotations

import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft  # noqa: F401
import modular_data as md
from modular_testbed import build_or_load_cache, cmaes
from cost_selector import SelectorEvaluator   # reuse the noise-aware evaluator

OUT_BASE = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT_BASE, exist_ok=True)

MODE = "effmass"

# 6 truths spanning the upper half plane — picked from cost_selector panel
# minus the trivially solved "interior" case.
STRESS_TRUTHS = [
    ("4-5-6 prod",   (-0.901, 0.794)),
    ("equilateral",  (-0.500, 0.866)),
    ("oblique short",(-0.300, 0.500)),
    ("tall Im≫Re",   ( 0.000, 1.300)),
    ("flat Im~½",    (-0.500, 0.500)),
    ("near-edge -1", (-0.950, 0.600)),
]


# --------------------------------------------------------------------------
def random_start(truth, dist, rng):
    """Sample start at exact Euclidean distance `dist` from `truth`."""
    theta = rng.uniform(0, 2 * np.pi)
    re = truth[0] + dist * np.cos(theta)
    im = truth[1] + dist * np.sin(theta)
    im = max(0.15, min(im, 1.55))   # keep in upper half plane (clipped)
    return (re, im)


def run_one(cache, ref_profile, *, tau_truth, tau_start, noise,
            frozen_noise, sigma0, max_evals, popsize,
            bounds_lo, bounds_hi, seed):
    ev = SelectorEvaluator(cache, ref_profile, MODE,
                            noise_sigma=noise,
                            noise_seed_base=1000 * seed,
                            frozen_noise=frozen_noise, pmean=1)
    try:
        xf, _, nev = cmaes(ev, [tau_start[0], tau_start[1]],
                            sigma0=sigma0, max_evals=max_evals,
                            popsize=popsize, seed=seed,
                            bounds_lo=bounds_lo, bounds_hi=bounds_hi)
        d_mod = md.modular_distance(complex(*xf), complex(*tau_truth),
                                     n_T=2, include_S=True)
        return dict(seed=seed, x=list(xf), dist=float(d_mod),
                    n_evals=int(nev))
    except Exception as e:
        return dict(seed=seed, crash=str(e))


def aggregate(runs, succ_thresh=0.1):
    ok = [r for r in runs if "dist" in r]
    if not ok: return dict(success=0.0, mean=float("nan"),
                            p50=float("nan"), p95=float("nan"),
                            std=float("nan"), n_ok=0, n_total=len(runs))
    d = np.array([r["dist"] for r in ok])
    return dict(
        success=100.0 * float(np.mean(d < succ_thresh)),
        mean=float(np.mean(d)),
        p50=float(np.percentile(d, 50)),
        p95=float(np.percentile(d, 95)),
        std=float(np.std(d)),
        n_ok=len(ok), n_total=len(runs),
    )


# --------------------------------------------------------------------------
def axis_A_noise_fine(cache, refs, panel, args, rng_master,
                      noise_levels, label):
    """Per-noise sweep across all panel scenarios × seeds."""
    print(f"\n  === AXIS A · {label} noise · {len(noise_levels)} levels ===")
    out = {}
    for sig in noise_levels:
        runs = []
        for name, tt in panel:
            for s in range(args.seeds):
                # Use a default start — small offset from truth
                ts = random_start(tt, 0.4,
                                   np.random.default_rng(7919 * s + 13))
                r = run_one(cache, refs[name], tau_truth=tt, tau_start=ts,
                            noise=sig, frozen_noise=(label == "frozen"),
                            sigma0=args.sigma0, max_evals=args.max_evals,
                            popsize=args.popsize,
                            bounds_lo=args.bounds_lo, bounds_hi=args.bounds_hi,
                            seed=s)
                r["scenario"] = name; r["truth"] = tt; r["start"] = ts
                runs.append(r)
        out[sig] = dict(runs=runs, agg=aggregate(runs))
        print(f"    σ={sig:<7g}  success={out[sig]['agg']['success']:5.1f}%  "
              f"mean={out[sig]['agg']['mean']:.3f}  "
              f"P95={out[sig]['agg']['p95']:.3f}")
    return out


def axis_B_Kmax(cache_factory, panel, args, K_values, fixed_sigmas):
    """Re-build cache at different K_max."""
    print(f"\n  === AXIS B · K_max sweep · K ∈ {K_values} ===")
    out = {}
    for K in K_values:
        cache, refs, blo, bhi = cache_factory(K_max=K)
        out[K] = {}
        for sig in fixed_sigmas:
            runs = []
            for name, tt in panel:
                for s in range(args.seeds):
                    ts = random_start(tt, 0.4,
                                       np.random.default_rng(7919 * s + 13))
                    r = run_one(cache, refs[name], tau_truth=tt, tau_start=ts,
                                noise=sig, frozen_noise=False,
                                sigma0=args.sigma0,
                                max_evals=args.max_evals,
                                popsize=args.popsize,
                                bounds_lo=blo, bounds_hi=bhi, seed=s)
                    runs.append(r)
            out[K][sig] = aggregate(runs)
            print(f"    K={K:<3d}  σ={sig:<6g}  "
                  f"success={out[K][sig]['success']:5.1f}%  "
                  f"mean={out[K][sig]['mean']:.3f}")
    return out


def axis_C_start_distance(cache, refs, panel, args, distances, fixed_sigmas):
    print(f"\n  === AXIS C · start-distance sweep · d ∈ {distances} ===")
    out = {}
    for d in distances:
        out[d] = {}
        for sig in fixed_sigmas:
            runs = []
            for name, tt in panel:
                for s in range(args.seeds):
                    ts = random_start(tt, d,
                                       np.random.default_rng(7919 * s + 13))
                    r = run_one(cache, refs[name], tau_truth=tt, tau_start=ts,
                                noise=sig, frozen_noise=False,
                                sigma0=args.sigma0,
                                max_evals=args.max_evals,
                                popsize=args.popsize,
                                bounds_lo=args.bounds_lo,
                                bounds_hi=args.bounds_hi, seed=s)
                    runs.append(r)
            out[d][sig] = dict(runs=runs, agg=aggregate(runs))
            print(f"    d={d:<5g}  σ={sig:<6g}  "
                  f"success={out[d][sig]['agg']['success']:5.1f}%  "
                  f"mean={out[d][sig]['agg']['mean']:.3f}")
    return out


def axis_D_seed_stability(cache, refs, panel, args, sigmas, seeds_for_stab):
    print(f"\n  === AXIS D · seed stability · {seeds_for_stab} seeds ===")
    out = {}
    for sig in sigmas:
        out[sig] = {}
        for name, tt in panel:
            ts = random_start(tt, 0.4, np.random.default_rng(42))
            xs = []; ds = []
            for s in range(seeds_for_stab):
                r = run_one(cache, refs[name], tau_truth=tt, tau_start=ts,
                            noise=sig, frozen_noise=False,
                            sigma0=args.sigma0, max_evals=args.max_evals,
                            popsize=args.popsize,
                            bounds_lo=args.bounds_lo,
                            bounds_hi=args.bounds_hi, seed=s)
                if "dist" in r:
                    xs.append(r["x"]); ds.append(r["dist"])
            xs = np.array(xs); ds = np.array(ds)
            out[sig][name] = dict(
                truth=tt, start=ts,
                xs=xs.tolist(), dists=ds.tolist(),
                std_re=float(np.std(xs[:, 0])) if len(xs) else float("nan"),
                std_im=float(np.std(xs[:, 1])) if len(xs) else float("nan"),
                success=100 * float(np.mean(ds < 0.1)) if len(ds) else 0.0,
            )
            print(f"    σ={sig:<6g}  {name:<14}  "
                  f"std(τ)=({out[sig][name]['std_re']:.3f},"
                  f"{out[sig][name]['std_im']:.3f})  "
                  f"success={out[sig][name]['success']:.0f}%")
    return out


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=80)
    ap.add_argument("--K-max", type=int, default=8)
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.4)
    ap.add_argument("--center", type=float, nargs=2, default=[-0.45, 0.85])
    ap.add_argument("--w1", type=float, default=16.0)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seeds-stab", type=int, default=20)
    ap.add_argument("--sigma0", type=float, default=0.35)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--scenarios", type=int, default=len(STRESS_TRUTHS))
    ap.add_argument("--axes", default="ABCDE",
                    help="which axes to run (subset of ABCDE)")
    ap.add_argument("--tag", default="prod")
    args = ap.parse_args()

    OUT = os.path.join(OUT_BASE, f"stress_effmass_{args.tag}")
    os.makedirs(OUT, exist_ok=True)

    cre, cim = args.center
    re_lo = cre - args.span_re / 2; re_hi = cre + args.span_re / 2
    im_lo = max(0.05, cim - args.span_im / 2)
    im_hi = cim + args.span_im / 2
    args.bounds_lo = np.array([re_lo, im_lo])
    args.bounds_hi = np.array([re_hi, im_hi])

    panel = STRESS_TRUTHS[:args.scenarios]
    print(f"  cost mode: {MODE}")
    print(f"  grid: Re τ ∈ [{re_lo:.3f},{re_hi:.3f}]  "
          f"Im τ ∈ [{im_lo:.3f},{im_hi:.3f}]  n={args.n}")
    print(f"  K_max={args.K_max}  |ω₁|={args.w1}  σ₀={args.sigma0}  "
          f"max_evals={args.max_evals}  seeds={args.seeds}")
    print(f"  panel ({len(panel)}): " + ", ".join(p[0] for p in panel))
    print(f"  axes: {args.axes}")
    print(f"  output → {os.path.relpath(OUT, HERE)}")

    # ----- build the default cache & refs -----
    def cache_factory(K_max):
        path = os.path.join(OUT_BASE,
                             f"modular_cache_n{args.n}_K{K_max}_w{int(args.w1)}.pkl")
        cache = build_or_load_cache(path, n=args.n,
                                    re_lo=re_lo, re_hi=re_hi,
                                    im_lo=im_lo, im_hi=im_hi,
                                    K_max=K_max, w1_norm=args.w1)
        refs = {name: md.make_threedir_profile(complex(*tt),
                                                w1_norm=args.w1, K_max=K_max)
                for name, tt in panel}
        return cache, refs, args.bounds_lo, args.bounds_hi

    cache, refs, _, _ = cache_factory(args.K_max)

    fine_noise = [0.0, 0.001, 0.002, 0.005, 0.008,
                  0.01, 0.015, 0.02, 0.03, 0.05]
    K_values = [3, 4, 6, 8, 10, 12]
    distances = [0.1, 0.3, 0.5, 0.8, 1.2]
    fixed_sigmas_BC = [0.0, 0.005, 0.02]
    stab_sigmas = [0.0, 0.005, 0.02]

    out = dict(args={k: (v.tolist() if isinstance(v, np.ndarray) else v)
                      for k, v in vars(args).items()},
                panel=panel)

    t0 = time.perf_counter()

    if "A" in args.axes:
        out["axis_A_stochastic"] = {
            f"sig={s}": dict(agg=v["agg"]) for s, v in
            axis_A_noise_fine(cache, refs, panel, args,
                               np.random.default_rng(0),
                               fine_noise, "stochastic").items()}
    if "E" in args.axes:
        out["axis_E_frozen"] = {
            f"sig={s}": dict(agg=v["agg"]) for s, v in
            axis_A_noise_fine(cache, refs, panel, args,
                               np.random.default_rng(1),
                               fine_noise, "frozen").items()}
    if "B" in args.axes:
        b = axis_B_Kmax(cache_factory, panel, args, K_values, fixed_sigmas_BC)
        out["axis_B_Kmax"] = {str(K): {str(s): a for s, a in v.items()}
                              for K, v in b.items()}
    if "C" in args.axes:
        c = axis_C_start_distance(cache, refs, panel, args,
                                    distances, fixed_sigmas_BC)
        out["axis_C_start"] = {str(d): {str(s): v["agg"] for s, v in dd.items()}
                               for d, dd in c.items()}
    if "D" in args.axes:
        d = axis_D_seed_stability(cache, refs, panel, args,
                                    stab_sigmas, args.seeds_stab)
        out["axis_D_stability"] = {
            str(s): {n: {k: v for k, v in info.items()
                          if k not in ("xs",)}    # drop big arrays from json
                     for n, info in dd.items()}
            for s, dd in d.items()}
    else:
        d = None

    wall = time.perf_counter() - t0
    print(f"\n  total wall: {wall:.1f}s")

    # ----- σ_50 breakpoint per axis A -----
    sigma_50 = None
    if "A" in args.axes:
        prev_sig, prev_succ = 0.0, 100.0
        for s in fine_noise:
            sc = out["axis_A_stochastic"][f"sig={s}"]["agg"]["success"]
            if sc < 50.0 and prev_succ >= 50.0:
                # linear interp
                if prev_succ == sc:
                    sigma_50 = s
                else:
                    sigma_50 = prev_sig + (s - prev_sig) * (
                        prev_succ - 50.0) / (prev_succ - sc)
                break
            prev_sig, prev_succ = s, sc

    # ----- text report -----
    rep = []
    rep.append("=" * 70)
    rep.append(f"STRESS REPORT — cost = {MODE}")
    rep.append("=" * 70)
    rep.append(f"  config: n={args.n}, K_max={args.K_max}, |ω₁|={args.w1},"
               f" max_evals={args.max_evals}, seeds={args.seeds}")
    rep.append(f"  panel : {len(panel)} truths   "
               f"({', '.join(p[0] for p in panel)})")
    rep.append("")

    if "A" in args.axes:
        rep.append("AXIS A — fine noise sweep (stochastic per-call)")
        rep.append(f"  {'σ':>7}  {'success%':>9}  {'mean':>6}  "
                   f"{'P50':>6}  {'P95':>6}")
        for s in fine_noise:
            a = out["axis_A_stochastic"][f"sig={s}"]["agg"]
            rep.append(f"  {s:>7g}  {a['success']:>9.1f}  {a['mean']:>6.3f}  "
                       f"{a['p50']:>6.3f}  {a['p95']:>6.3f}")
        if sigma_50 is not None:
            rep.append(f"\n  ★ σ_50 (success drops below 50%) ≈ {sigma_50:.4f}")
        else:
            rep.append("\n  ★ σ_50 not reached within sweep range")
        rep.append("")

    if "E" in args.axes:
        rep.append("AXIS E — frozen vs stochastic")
        rep.append(f"  {'σ':>7}  {'stoch%':>7}  {'frozen%':>8}  Δ")
        for s in fine_noise:
            sa = out["axis_A_stochastic"][f"sig={s}"]["agg"]["success"]
            fa = out["axis_E_frozen"][f"sig={s}"]["agg"]["success"]
            rep.append(f"  {s:>7g}  {sa:>7.1f}  {fa:>8.1f}  "
                       f"{fa - sa:+.1f}")
        rep.append("")

    if "B" in args.axes:
        rep.append("AXIS B — K_max sweep")
        rep.append(f"  {'K_max':>5}  " + "  ".join(
            f"σ={s:>5g}" for s in fixed_sigmas_BC))
        for K in K_values:
            row = "  ".join(f"{out['axis_B_Kmax'][str(K)][str(s)]['success']:>6.1f}"
                            for s in fixed_sigmas_BC)
            rep.append(f"  {K:>5d}  {row}")
        rep.append("")

    if "C" in args.axes:
        rep.append("AXIS C — start distance |τ_start − τ_truth|")
        rep.append(f"  {'d':>5}  " + "  ".join(
            f"σ={s:>5g}" for s in fixed_sigmas_BC))
        for dd in distances:
            row = "  ".join(f"{out['axis_C_start'][str(dd)][str(s)]['success']:>6.1f}"
                            for s in fixed_sigmas_BC)
            rep.append(f"  {dd:>5g}  {row}")
        rep.append("")

    if "D" in args.axes:
        rep.append("AXIS D — seed stability (std of converged τ)")
        rep.append(f"  {'σ':>7}  {'scenario':<14}  "
                   f"{'std_Re':>7}  {'std_Im':>7}  {'success%':>8}")
        for s in stab_sigmas:
            for name, tt in panel:
                info = out["axis_D_stability"][str(s)][name]
                rep.append(f"  {s:>7g}  {name:<14}  "
                           f"{info['std_re']:>7.3f}  {info['std_im']:>7.3f}  "
                           f"{info['success']:>8.1f}")
        rep.append("")

    rep.append("=" * 70)
    rep.append("PASS / FAIL VERDICT")
    rep.append("=" * 70)
    verdict = []
    if sigma_50 is not None:
        if sigma_50 < 0.005:
            verdict.append(f"  ✗ FRAGILE: σ_50 ≈ {sigma_50:.4f} (< 0.5%)")
        elif sigma_50 < 0.02:
            verdict.append(f"  ⚠ MODERATE: σ_50 ≈ {sigma_50:.4f} (0.5–2%)")
        else:
            verdict.append(f"  ✓ ROBUST: σ_50 ≈ {sigma_50:.4f} (> 2%)")
    if "B" in args.axes:
        # Best K at σ=0
        s0 = "0.0"
        bestK = max(K_values,
                     key=lambda K: out["axis_B_Kmax"][str(K)][s0]["success"])
        verdict.append(f"  best K_max at σ=0: {bestK}")
        # K-trend with noise
        s_noisy = "0.02"
        succ_lowK = out["axis_B_Kmax"][str(min(K_values))][s_noisy]["success"]
        succ_hiK = out["axis_B_Kmax"][str(max(K_values))][s_noisy]["success"]
        if succ_lowK > succ_hiK + 5:
            verdict.append(f"  ★ short m_eff window (K={min(K_values)}) "
                           f"survives noise better ({succ_lowK:.0f}% vs "
                           f"{succ_hiK:.0f}% at K={max(K_values)})")
        elif succ_hiK > succ_lowK + 5:
            verdict.append(f"  ★ long m_eff window helps under noise "
                           f"(K={max(K_values)}: {succ_hiK:.0f}% > "
                           f"K={min(K_values)}: {succ_lowK:.0f}%)")
    if "C" in args.axes:
        # basin radius at noise-free
        s0 = "0.0"
        last_ok_d = None
        for dd in distances:
            if out["axis_C_start"][str(dd)][s0]["success"] >= 50.0:
                last_ok_d = dd
        if last_ok_d is not None:
            verdict.append(f"  basin radius (σ=0, success ≥ 50%): "
                           f"d ≈ {last_ok_d}")
    rep.extend(verdict)
    rep.append("=" * 70)

    text = "\n".join(rep)
    print()
    print(text)
    with open(os.path.join(OUT, "stress_report.txt"), "w") as f:
        f.write(text)

    # ----- save json -----
    with open(os.path.join(OUT, "stress.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)

    # ----- plots -----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # AXIS A
        if "A" in args.axes:
            fig, axs = plt.subplots(1, 3, figsize=(13.5, 4))
            xs = [max(s, 1e-4) for s in fine_noise]
            succ = [out["axis_A_stochastic"][f"sig={s}"]["agg"]["success"]
                    for s in fine_noise]
            mean = [out["axis_A_stochastic"][f"sig={s}"]["agg"]["mean"]
                    for s in fine_noise]
            p95 = [out["axis_A_stochastic"][f"sig={s}"]["agg"]["p95"]
                   for s in fine_noise]
            axs[0].semilogx(xs, succ, "-o", color="C2", lw=1.8)
            axs[0].axhline(50, color="r", ls="--", alpha=0.5, label="50%")
            if sigma_50:
                axs[0].axvline(sigma_50, color="r", ls=":", alpha=0.7,
                               label=f"σ₅₀={sigma_50:.4f}")
            axs[0].set_xlabel("σ_rel"); axs[0].set_ylabel("success%")
            axs[0].set_title("AXIS A · success vs noise")
            axs[0].grid(alpha=0.3); axs[0].legend()
            axs[1].loglog(xs, mean, "-o", color="C0", lw=1.8, label="mean")
            axs[1].loglog(xs, p95,  "-o", color="C3", lw=1.5, label="P95")
            axs[1].set_xlabel("σ_rel"); axs[1].set_ylabel("modular dist")
            axs[1].set_title("dist (mean & P95) vs noise")
            axs[1].grid(which="both", alpha=0.3); axs[1].legend()
            # AXIS E overlay
            if "E" in args.axes:
                succ_f = [out["axis_E_frozen"][f"sig={s}"]["agg"]["success"]
                          for s in fine_noise]
                axs[2].semilogx(xs, succ,   "-o", color="C2", label="stochastic")
                axs[2].semilogx(xs, succ_f, "-s", color="C4", label="frozen")
                axs[2].set_xlabel("σ_rel"); axs[2].set_ylabel("success%")
                axs[2].set_title("AXIS E · stoch vs frozen noise")
                axs[2].grid(alpha=0.3); axs[2].legend()
            else:
                axs[2].axis("off")
            fig.suptitle(f"stress_effmass · noise sweep · "
                         f"n={args.n} K={args.K_max} seeds={args.seeds}",
                         fontsize=11)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "noise_fine.png"), dpi=130)
            plt.close(fig)

        if "B" in args.axes:
            fig, ax = plt.subplots(figsize=(7, 4.5))
            for sig in fixed_sigmas_BC:
                ys = [out["axis_B_Kmax"][str(K)][str(sig)]["success"]
                      for K in K_values]
                ax.plot(K_values, ys, "-o", lw=1.6, label=f"σ={sig}")
            ax.set_xlabel("K_max (m_eff window depth)")
            ax.set_ylabel("success%")
            ax.set_title("AXIS B · success vs K_max")
            ax.grid(alpha=0.3); ax.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "Kmax_scan.png"), dpi=130)
            plt.close(fig)

        if "C" in args.axes:
            fig, ax = plt.subplots(figsize=(7, 4.5))
            for sig in fixed_sigmas_BC:
                ys = [out["axis_C_start"][str(d)][str(sig)]["success"]
                      for d in distances]
                ax.plot(distances, ys, "-o", lw=1.6, label=f"σ={sig}")
            ax.set_xlabel("|τ_start − τ_truth|")
            ax.set_ylabel("success%")
            ax.set_title("AXIS C · basin radius (success vs start distance)")
            ax.grid(alpha=0.3); ax.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "start_distance.png"), dpi=130)
            plt.close(fig)

        if "D" in args.axes and d is not None:
            ncols = len(panel)
            nrows = len(stab_sigmas)
            fig, axs = plt.subplots(nrows, ncols,
                                      figsize=(2.4 * ncols, 2.4 * nrows),
                                      squeeze=False)
            for i, sig in enumerate(stab_sigmas):
                for j, (name, tt) in enumerate(panel):
                    a = axs[i, j]
                    info = d[sig][name]
                    xs = np.array(info["xs"])
                    if len(xs):
                        ds_arr = np.array(info["dists"])
                        good = ds_arr < 0.1
                        a.scatter(xs[good, 0], xs[good, 1], s=20,
                                  c="C2", alpha=0.7, label="✓")
                        a.scatter(xs[~good, 0], xs[~good, 1], s=20,
                                  c="C3", alpha=0.7, label="✗")
                    a.scatter(*tt, marker="*", s=140, color="k", zorder=5)
                    a.scatter(*info["start"], marker="x", s=60, color="C0")
                    a.set_title(f"σ={sig} {name}", fontsize=8)
                    a.set_xlim(re_lo, re_hi); a.set_ylim(im_lo, im_hi)
                    a.tick_params(labelsize=6)
            fig.suptitle("AXIS D · seed stability "
                         "(★ truth, × start, ✓/✗ converged)",
                         fontsize=11)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "seed_stability.png"), dpi=130)
            plt.close(fig)

            # failure atlas — combine all converged points colored by σ
            fig, ax = plt.subplots(figsize=(7, 5))
            for sig, color in zip(stab_sigmas, ["#1a9850", "#fdae61", "#d73027"]):
                allx, ally = [], []
                for name, tt in panel:
                    info = d[sig][name]
                    xs = np.array(info["xs"])
                    if len(xs):
                        allx.extend(xs[:, 0]); ally.extend(xs[:, 1])
                ax.scatter(allx, ally, s=12, alpha=0.6, color=color,
                           label=f"σ={sig}")
            for name, tt in panel:
                ax.scatter(*tt, marker="*", s=160, color="k", zorder=5)
            ax.set_xlim(re_lo, re_hi); ax.set_ylim(im_lo, im_hi)
            ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
            ax.set_title("Failure atlas — converged τ across all seeds")
            ax.legend(); ax.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "failure_atlas.png"), dpi=130)
            plt.close(fig)
        print()
        print(f"  artifacts → {os.path.relpath(OUT, HERE)}")
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"  plotting failed: {e}")


if __name__ == "__main__":
    main()
