#!/usr/bin/env python3
"""
modular_testbed.py — closed-loop noise-free (or noise-controlled) test
of CMA-ES + cost-function variants on the analytic Ising-CFT modular plane.

WORKFLOW
--------
1. Pre-compute analytic three-direction correlator profiles
   (Re τ, Im τ, Im τ−Re τ — the three torus geodesics) on an
   n × n grid in the (Re τ, Im τ) plane.  Cache to disk.

2. Pick a target τ_truth and a starting τ_start (defaults: τ_truth =
   the τ of the production reference geometry; τ_start = the τ of the
   production test geometry).

3. For each cost mode in {l2, relative, log, effmass, pmean4}, run a
   minimal 2-D CMA-ES in the continuous (Re τ, Im τ) plane.  Each
   evaluation snaps to the nearest grid cell and computes
   threedir_cost(ref_profile, cell_profile).

4. Optional: --noise σ injects multiplicative Gaussian noise on every
   profile entry, with a deterministic seed per (i, j, eval_idx) so the
   cost is repeatable but reflects realistic MC scatter.

5. Outputs per cost mode: trajectory plot overlaid on cost heatmap,
   summary JSON, ranking table.

USAGE
-----
    # noise-free
    python modular_testbed.py --n 21

    # with 5% relative noise
    python modular_testbed.py --n 21 --noise 0.05 --seeds 5
"""
from __future__ import annotations

import argparse, json, math, os, pickle, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft
import modular_data as md

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------------
# Cache: precompute three-direction profiles on a τ grid
# ----------------------------------------------------------------------
def build_or_load_cache(path: str, *, n: int,
                        re_lo, re_hi, im_lo, im_hi,
                        K_max: int, w1_norm: float):
    if os.path.exists(path):
        with open(path, "rb") as f:
            cache = pickle.load(f)
        if (cache["n"] == n
                and cache["K_max"] == K_max
                and abs(cache["w1_norm"] - w1_norm) < 1e-12
                and cache["re_lo"] == re_lo and cache["re_hi"] == re_hi
                and cache["im_lo"] == im_lo and cache["im_hi"] == im_hi):
            print(f"  loaded cache from {os.path.relpath(path, HERE)}  "
                  f"({n}×{n} cells)")
            return cache
        print("  cache parameters mismatch → rebuilding")

    re_arr = np.linspace(re_lo, re_hi, n)
    im_arr = np.linspace(im_lo, im_hi, n)
    grid = {}
    t0 = time.perf_counter()
    n_total = n * n
    cell = 0
    for i, re in enumerate(re_arr):
        for j, im in enumerate(im_arr):
            cell += 1
            if cell % max(1, n_total // 10) == 0:
                print(f"    cache cell {cell}/{n_total}  "
                      f"({100 * cell / n_total:.0f}%)")
            try:
                grid[(i, j)] = md.make_threedir_profile(
                    complex(re, im), w1_norm=w1_norm, K_max=K_max)
            except Exception:
                grid[(i, j)] = None
    wall = time.perf_counter() - t0
    print(f"  cache built in {wall:.1f}s  "
          f"({wall / n_total * 1000:.1f} ms/cell)")
    cache = dict(n=n, re=re_arr.tolist(), im=im_arr.tolist(),
                 K_max=K_max, w1_norm=w1_norm,
                 re_lo=re_lo, re_hi=re_hi, im_lo=im_lo, im_hi=im_hi,
                 grid=grid)
    with open(path, "wb") as f:
        pickle.dump(cache, f)
    print(f"  cache saved → {os.path.relpath(path, HERE)}")
    return cache


# ----------------------------------------------------------------------
# Evaluator
# ----------------------------------------------------------------------
class GridEvaluator:
    def __init__(self, cache, ref_profile, mode, *,
                 noise_sigma=0.0, noise_seed=0, pmean=1):
        self.cache = cache
        self.ref = ref_profile
        self.mode = mode
        self.pmean = pmean
        self.noise_sigma = float(noise_sigma)
        self.noise_seed = int(noise_seed)
        self.re = np.array(cache["re"])
        self.im = np.array(cache["im"])
        self.calls = 0
        self.history = []      # (re, im, cost) per call

    def _maybe_noisy(self, profile, ij):
        if self.noise_sigma <= 0:
            return profile
        # Deterministic seed per (variant, seed-base, cell, call).
        h = hash((self.mode, self.noise_seed, ij, self.calls)) & 0xFFFFFFFF
        rng = np.random.default_rng(h)
        return md.add_noise(profile, self.noise_sigma, rng)

    def __call__(self, x):
        re, im = float(x[0]), float(x[1])
        i = int(np.argmin(np.abs(self.re - re)))
        j = int(np.argmin(np.abs(self.im - im)))
        prof = self.cache["grid"].get((i, j))
        if prof is None:
            self.calls += 1
            self.history.append((re, im, float("inf")))
            return float("inf")
        prof_eff = self._maybe_noisy(prof, (i, j))
        try:
            c, _ = md.threedir_cost(self.ref, prof_eff,
                                    mode=self.mode, pmean=self.pmean)
        except Exception:
            c = float("inf")
        if not np.isfinite(c):
            c = 1e30
        self.calls += 1
        self.history.append((re, im, c))
        return c


# ----------------------------------------------------------------------
# Minimal (μ/μ_w, λ)-CMA-ES (2-D, bounded)
# ----------------------------------------------------------------------
def cmaes(eval_fn, x0, *, sigma0, max_evals, popsize=8, seed=0,
          bounds_lo=None, bounds_hi=None):
    rng = np.random.default_rng(seed)
    n = len(x0)
    mean = np.array(x0, dtype=float)
    sigma = float(sigma0)
    C = np.eye(n)
    pc = np.zeros(n)
    ps = np.zeros(n)
    lam = popsize
    mu = lam // 2
    weights = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
    weights /= weights.sum()
    mu_eff = 1.0 / (weights ** 2).sum()
    cs = (mu_eff + 2) / (n + mu_eff + 5)
    cc = (4 + mu_eff / n) / (n + 4 + 2 * mu_eff / n)
    c1 = 2.0 / ((n + 1.3) ** 2 + mu_eff)
    cmu = min(1 - c1,
              2 * (mu_eff - 2 + 1 / mu_eff) / ((n + 2) ** 2 + mu_eff))
    damps = 1 + 2 * max(0, math.sqrt((mu_eff - 1) / (n + 1)) - 1) + cs
    chiN = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))

    history_mean = []
    n_evals = 0
    gen = 0
    while n_evals < max_evals:
        try:
            B, D2, _ = np.linalg.svd(C); D = np.sqrt(np.maximum(D2, 1e-30))
        except np.linalg.LinAlgError:
            B = np.eye(n); D = np.ones(n)
        offspring = []
        for _ in range(lam):
            z = rng.standard_normal(n)
            x = mean + sigma * (B @ (D * z))
            if bounds_lo is not None: x = np.maximum(x, bounds_lo)
            if bounds_hi is not None: x = np.minimum(x, bounds_hi)
            f = eval_fn(x)
            n_evals += 1
            offspring.append((f, x, z))
            if n_evals >= max_evals: break
        offspring.sort(key=lambda t: t[0])
        sel = offspring[:mu]
        x_sel = np.stack([t[1] for t in sel])
        z_sel = np.stack([t[2] for t in sel])
        mean = weights @ x_sel
        z_mean = weights @ z_sel
        ps = (1 - cs) * ps + math.sqrt(cs * (2 - cs) * mu_eff) * (B @ z_mean)
        norm_ps = np.linalg.norm(ps)
        hs = norm_ps / math.sqrt(1 - (1 - cs) ** (2 * (gen + 1))) / chiN \
             < 1.4 + 2 / (n + 1)
        pc = (1 - cc) * pc + (1 if hs else 0) * \
             math.sqrt(cc * (2 - cc) * mu_eff) * (B @ (D * z_mean))
        artmp = (x_sel - mean) / sigma
        C = (1 - c1 - cmu) * C \
            + c1 * (np.outer(pc, pc) + (0 if hs else cc * (2 - cc)) * C) \
            + cmu * (artmp.T @ np.diag(weights) @ artmp)
        sigma *= math.exp((cs / damps) * (norm_ps / chiN - 1))
        sigma = float(np.clip(sigma, 1e-6, 5.0))
        gen += 1
        history_mean.append((float(mean[0]), float(mean[1]),
                             float(offspring[0][0])))
        if sigma * float(D.max()) < 1e-3: break
        if len(history_mean) >= 5:
            recent = [h[2] for h in history_mean[-5:]]
            if max(recent) - min(recent) < 1e-9: break
    return mean, history_mean, n_evals


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=21)
    ap.add_argument("--K-max", type=int, default=8,
                    help="max # samples per direction")
    ap.add_argument("--ref",  type=int, nargs=4, default=[13, 16, -3, 3],
                    help="(Lx, Ly, Tx, Ty) used only to derive τ_truth")
    ap.add_argument("--test", type=int, nargs=4, default=[16, 16,  0,  0],
                    help="(Lx, Ly, Tx, Ty) used only to derive τ_start")
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.7)
    ap.add_argument("--start", type=float, nargs=2, default=None)
    ap.add_argument("--sigma0", type=float, default=0.25)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--noise", type=float, default=0.0,
                    help="multiplicative relative noise σ on every G entry")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    Rx, Ry, RTx, RTy = args.ref
    Tx, Ty, TTx, TTy = args.test
    w1_r, w2_r = cft.torus_periods(Rx, Ry, RTx, RTy)
    tau_truth = w2_r / w1_r
    if tau_truth.imag < 0:
        tau_truth = complex(tau_truth.real, -tau_truth.imag)
    w1_t, w2_t = cft.torus_periods(Tx, Ty, TTx, TTy)
    tau_t_nat = w2_t / w1_t
    if tau_t_nat.imag < 0:
        tau_t_nat = complex(tau_t_nat.real, -tau_t_nat.imag)
    tau_start = (complex(args.start[0], args.start[1])
                 if args.start is not None else tau_t_nat)
    print(f"  τ_truth  = {tau_truth:.4f}  (from ref geom)")
    print(f"  τ_start  = {tau_start:.4f}  "
          f"(natural test τ = {tau_t_nat:.4f})")

    # τ-grid centred on τ_truth
    re_lo = tau_truth.real - args.span_re / 2
    re_hi = tau_truth.real + args.span_re / 2
    im_lo = max(0.05, tau_truth.imag - args.span_im / 2)
    im_hi = tau_truth.imag + args.span_im / 2
    print(f"  grid: Re τ ∈ [{re_lo:.3f}, {re_hi:.3f}], "
          f"Im τ ∈ [{im_lo:.3f}, {im_hi:.3f}], n = {args.n}")
    print(f"  noise σ_rel = {args.noise:.3g}")

    cache_path = args.cache or os.path.join(
        OUT, f"modular_cache_n{args.n}_K{args.K_max}.pkl")

    print("  building reference profile (analytic) at τ_truth …")
    ref_profile = md.make_threedir_profile(
        tau_truth, w1_norm=abs(w1_r), K_max=args.K_max)
    print(f"  ref K = {ref_profile['K']}; "
          f"|G|_re sample: {ref_profile['profiles']['re'][:3]}")

    print("  building/loading τ-grid cache of test profiles …")
    cache = build_or_load_cache(cache_path, n=args.n,
                                re_lo=re_lo, re_hi=re_hi,
                                im_lo=im_lo, im_hi=im_hi,
                                K_max=args.K_max, w1_norm=abs(w1_r))

    bounds_lo = np.array([re_lo, im_lo])
    bounds_hi = np.array([re_hi, im_hi])
    x_truth = np.array([tau_truth.real, tau_truth.imag])

    modes = ["l2", "relative", "log", "effmass", "pmean4"]
    print(f"  cost modes: {modes}")

    results = {}
    print()
    hdr = (f"  {'mode':<10}  {'seed':>4}  {'final τ':>16}  "
           f"{'dist→truth':>10}  {'cost@truth':>10}  {'#evals':>6}  status")
    print(hdr); print("  " + "-" * len(hdr))
    for mode in modes:
        runs = []
        for seed in range(args.seeds):
            pmean = 4 if mode == "pmean4" else 1
            base_mode = "l2" if mode == "pmean4" else mode
            ev = GridEvaluator(cache, ref_profile, base_mode,
                               noise_sigma=args.noise,
                               noise_seed=1000 * seed,
                               pmean=pmean)
            try:
                xf, hist, nev = cmaes(
                    ev, [tau_start.real, tau_start.imag],
                    sigma0=args.sigma0, max_evals=args.max_evals,
                    popsize=args.popsize, seed=seed,
                    bounds_lo=bounds_lo, bounds_hi=bounds_hi)
            except Exception as e:
                print(f"  {mode:<10}  seed {seed}  CRASH: {e}")
                continue
            d_truth = float(np.linalg.norm(xf - x_truth))
            c_truth = ev(x_truth)
            ev.calls -= 1; ev.history.pop()
            status = ("✓ found" if d_truth < 0.1
                      else "near"  if d_truth < 0.3
                      else "✗ lost")
            print(f"  {mode:<10}  {seed:>4}  "
                  f"({xf[0]:+.3f},{xf[1]:+.3f})  "
                  f"{d_truth:>10.3f}  {c_truth:>10.2e}  "
                  f"{nev:>6d}  {status}")
            runs.append(dict(seed=seed, x_final=xf.tolist(),
                             dist_truth=d_truth, n_evals=nev,
                             status=status, history=hist))
        results[mode] = runs

    print()
    print("  PERFORMANCE SUMMARY:")
    print(f"  {'mode':<10}  {'mean dist':>10}  {'best dist':>10}  "
          f"{'success%':>9}  {'mean evals':>10}")
    print("  " + "-" * 56)
    rank_rows = []
    for mode in modes:
        runs = results.get(mode, [])
        if not runs: continue
        dists = [r["dist_truth"] for r in runs]
        evals = [r["n_evals"] for r in runs]
        succ = 100.0 * sum(1 for r in runs if "✓" in r["status"]) / len(runs)
        print(f"  {mode:<10}  {np.mean(dists):>10.3f}  "
              f"{np.min(dists):>10.3f}  {succ:>8.0f}%  "
              f"{np.mean(evals):>10.1f}")
        rank_rows.append((mode, np.mean(dists), np.min(dists), succ))
    rank_rows.sort(key=lambda r: (-r[3], r[1]))
    print()
    print("  RANKING (success% desc, mean dist asc):")
    for i, (m, md_, bd, s) in enumerate(rank_rows, 1):
        print(f"    {i}. {m:<10}  success={s:>3.0f}%  "
              f"mean_dist={md_:.3f}  best_dist={bd:.3f}")

    suffix = (args.tag if args.tag
              else f"noise{args.noise:.2g}").replace(".", "p")
    json_path = os.path.join(OUT, f"modular_testbed_{suffix}.json")
    with open(json_path, "w") as f:
        json.dump(dict(
            tau_truth=[tau_truth.real, tau_truth.imag],
            tau_start=[tau_start.real, tau_start.imag],
            args=vars(args),
            ranking=[dict(mode=m, mean_dist=md_, best_dist=bd,
                          success_pct=s) for m, md_, bd, s in rank_rows],
            results={m: [{k: v for k, v in r.items() if k != "history"}
                         for r in rs] for m, rs in results.items()},
        ), f, indent=2)
    print(f"\n  summary → {os.path.relpath(json_path, HERE)}")

    # ------------------------------------------------------------------
    # Plot: cost heatmap per mode + trajectories overlay
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        re_arr = np.array(cache["re"])
        im_arr = np.array(cache["im"])
        ncols = 3
        nrows = (len(modes) + ncols - 1) // ncols
        fig, axs = plt.subplots(nrows, ncols,
                                figsize=(5.2 * ncols, 4.6 * nrows),
                                squeeze=False)
        for k, mode in enumerate(modes):
            r, c = k // ncols, k % ncols
            ax = axs[r, c]
            grid_cost = np.full((args.n, args.n), np.nan)
            base_mode = "l2" if mode == "pmean4" else mode
            pmean = 4 if mode == "pmean4" else 1
            for (i, j), prof in cache["grid"].items():
                if prof is None: continue
                try:
                    cv, _ = md.threedir_cost(ref_profile, prof,
                                             mode=base_mode, pmean=pmean)
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
            plt.colorbar(im, ax=ax, label="log₁₀ cost")
            ax.scatter([tau_truth.real], [tau_truth.imag], marker="*",
                       s=320, color="red", edgecolor="white", lw=1.4,
                       zorder=10, label="truth")
            ax.scatter([tau_start.real], [tau_start.imag], marker="o",
                       s=140, color="cyan", edgecolor="black", lw=1.2,
                       zorder=9, label="start")
            colors = plt.cm.tab10(np.linspace(0, 1, args.seeds))
            for run, col in zip(results.get(mode, []), colors):
                hist = np.array([h[:2] for h in run["history"]])
                if len(hist):
                    ax.plot(hist[:, 0], hist[:, 1], "-",
                            color=col, lw=1.2, alpha=0.9)
                    ax.scatter([hist[-1, 0]], [hist[-1, 1]],
                               marker="X", s=80, color=col,
                               edgecolor="black", lw=0.6, zorder=8)
            # mark cost min on grid
            if np.isfinite(grid_cost).any():
                i_min, j_min = np.unravel_index(
                    np.nanargmin(grid_cost), grid_cost.shape)
                ax.scatter([re_arr[i_min]], [im_arr[j_min]], marker="P",
                           s=120, color="orange", edgecolor="black",
                           lw=0.6, zorder=11, label="grid min")
            ax.set_xlim(re_lo, re_hi); ax.set_ylim(im_lo, im_hi)
            ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
            md_dist = (np.mean([r["dist_truth"] for r in results.get(mode, [])])
                       if results.get(mode) else float("nan"))
            succ = (100.0 * sum(1 for r in results.get(mode, [])
                                if "✓" in r["status"])
                    / max(1, len(results.get(mode, []))))
            ax.set_title(f"{mode}  (σ_rel={args.noise:g})\n"
                         f"success={succ:.0f}%  mean dist={md_dist:.2f}",
                         fontsize=10)
            ax.legend(loc="upper right", fontsize=7)
        for k in range(len(modes), nrows * ncols):
            axs[k // ncols, k % ncols].axis("off")
        plt.tight_layout()
        png = os.path.join(OUT, f"modular_testbed_{suffix}.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  plot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"  plot failed: {e}")


if __name__ == "__main__":
    main()
