#!/usr/bin/env python3
"""
zoo_grid_cmaes_navigate.py — closed-loop noise-free test:

Pre-compute analytic Ising correlator data dicts on a (Re τ, Im τ) grid,
pick a target (truth) τ₁ and a starting τ₂, then for each cost variant
run CMA-ES in the continuous (Re τ, Im τ) plane.  The evaluator does a
nearest-grid lookup of the test data dict and computes cost vs the
fixed reference data dict at τ₁.

This separates the cost-surface bug from the optimizer bug:
    - only one source of noise:  the grid discretisation
    - the optimizer is the same hand-rolled CMA-ES family used in
      production (here: minimal 2-D variant)
    - convergence to the true τ depends ONLY on cost-surface topology

Outputs:
    * pickle cache of the grid (reusable)
    * one trajectory plot per variant overlaid on the cost heatmap
    * summary table: dist(final, truth), # evals, success
"""
from __future__ import annotations
import argparse, io, math, os, pickle, sys, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft
import cost as cost_mod
import cost_zoo

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------------
# 1. Grid cache
# ----------------------------------------------------------------------

def build_or_load_cache(path: str, *, n: int,
                        Tx: int, Ty: int, TTx: int, TTy: int,
                        w1_ref: complex,
                        re_lo: float, re_hi: float,
                        im_lo: float, im_hi: float):
    """Build (or load) {(i,j): ising_data_dict} on an n×n grid in τ.

    ω₂(τ) = τ · ω₁_ref so that τ_test_eff = τ exactly.
    """
    if os.path.exists(path):
        with open(path, "rb") as f:
            cache = pickle.load(f)
        if (cache["n"] == n
                and cache["test_geom"] == (Tx, Ty, TTx, TTy)
                and abs(cache["w1_ref"] - w1_ref) < 1e-12
                and cache["re_lo"] == re_lo and cache["re_hi"] == re_hi
                and cache["im_lo"] == im_lo and cache["im_hi"] == im_hi):
            print(f"  loaded cache from {os.path.relpath(path, HERE)}  "
                  f"({n}×{n} = {n * n} cells)")
            return cache
        print("  cache parameters mismatch → rebuilding")

    re_arr = np.linspace(re_lo, re_hi, n)
    im_arr = np.linspace(im_lo, im_hi, n)
    grid: dict[tuple[int, int], dict] = {}
    real_stdout = sys.stdout
    null = io.StringIO()
    t0 = time.perf_counter()
    n_total = n * n
    cell = 0
    for i, re in enumerate(re_arr):
        for j, im in enumerate(im_arr):
            cell += 1
            if cell % max(1, n_total // 10) == 0:
                real_stdout.write(f"    cache cell {cell}/{n_total}  "
                                  f"({100 * cell / n_total:.0f}%)\n")
                real_stdout.flush()
            tau = complex(re, im)
            sys.stdout = null
            try:
                d = cft.make_ising_data_dict(
                    Tx, Ty, TTx, TTy,
                    omega1_override=w1_ref,
                    omega2_override=tau * w1_ref)
                grid[(i, j)] = d
            finally:
                sys.stdout = real_stdout
    wall = time.perf_counter() - t0
    print(f"  cache built in {wall:.1f}s  ({wall / n_total * 1000:.1f} ms/cell)")

    cache = dict(n=n, re=re_arr.tolist(), im=im_arr.tolist(),
                 test_geom=(Tx, Ty, TTx, TTy), w1_ref=w1_ref,
                 re_lo=re_lo, re_hi=re_hi, im_lo=im_lo, im_hi=im_hi,
                 grid=grid)
    with open(path, "wb") as f:
        pickle.dump(cache, f)
    print(f"  cache saved → {os.path.relpath(path, HERE)}")
    return cache


# ----------------------------------------------------------------------
# 2. Evaluator: nearest-grid lookup + cost
# ----------------------------------------------------------------------

class GridEvaluator:
    def __init__(self, cache, ref_data, cost_fn, *,
                 Tx, Ty, TTx, TTy, Rx, Ry, RTx, RTy):
        self.cache = cache
        self.ref_data = ref_data
        self.fn = cost_fn
        self.Tx, self.Ty, self.TTx, self.TTy = Tx, Ty, TTx, TTy
        self.Rx, self.Ry, self.RTx, self.RTy = Rx, Ry, RTx, RTy
        self.re = np.array(cache["re"])
        self.im = np.array(cache["im"])
        self.calls = 0
        self.history: list[tuple[float, float, float]] = []

    def __call__(self, x):
        re, im = float(x[0]), float(x[1])
        # nearest grid cell
        i = int(np.argmin(np.abs(self.re - re)))
        j = int(np.argmin(np.abs(self.im - im)))
        d = self.cache["grid"].get((i, j))
        if d is None:
            return float("inf")
        null = io.StringIO()
        old = sys.stdout
        sys.stdout = null
        try:
            c, *_ = self.fn(self.ref_data, d,
                            self.Tx, self.Ty, self.TTx, self.TTy,
                            self.Rx, self.Ry, self.RTx, self.RTy,
                            copies=2, power=2)
        except Exception:
            c = float("inf")
        finally:
            sys.stdout = old
        if not np.isfinite(c):
            c = 1e30
        self.calls += 1
        self.history.append((re, im, c))
        return c


# ----------------------------------------------------------------------
# 3. Minimal (μ/μ_w, λ)-CMA-ES, 2-D, deterministic, bounded
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

    history_mean: list[tuple[float, float, float]] = []
    n_evals = 0
    gen = 0
    while n_evals < max_evals:
        # sample
        try:
            B, D2, _ = np.linalg.svd(C); D = np.sqrt(np.maximum(D2, 1e-30))
        except np.linalg.LinAlgError:
            B = np.eye(n); D = np.ones(n)
        offspring = []
        for _ in range(lam):
            z = rng.standard_normal(n)
            x = mean + sigma * (B @ (D * z))
            if bounds_lo is not None:
                x = np.maximum(x, bounds_lo)
            if bounds_hi is not None:
                x = np.minimum(x, bounds_hi)
            f = eval_fn(x)
            n_evals += 1
            offspring.append((f, x, z))
            if n_evals >= max_evals: break
        offspring.sort(key=lambda t: t[0])
        # Best mu
        sel = offspring[:mu]
        # New mean
        x_sel = np.stack([t[1] for t in sel])
        z_sel = np.stack([t[2] for t in sel])
        mean = weights @ x_sel
        # Update step path
        z_mean = weights @ z_sel
        ps = (1 - cs) * ps + math.sqrt(cs * (2 - cs) * mu_eff) * (B @ z_mean)
        norm_ps = np.linalg.norm(ps)
        hs = norm_ps / math.sqrt(1 - (1 - cs) ** (2 * (gen + 1))) / chiN \
             < 1.4 + 2 / (n + 1)
        pc = (1 - cc) * pc + (1 if hs else 0) * \
             math.sqrt(cc * (2 - cc) * mu_eff) * (B @ (D * z_mean))
        # Covariance
        artmp = (x_sel - mean) / sigma
        C = (1 - c1 - cmu) * C \
            + c1 * (np.outer(pc, pc) + (0 if hs else cc * (2 - cc)) * C) \
            + cmu * (artmp.T @ np.diag(weights) @ artmp)
        # Step size
        sigma *= math.exp((cs / damps) * (norm_ps / chiN - 1))
        sigma = float(np.clip(sigma, 1e-6, 5.0))
        gen += 1
        history_mean.append((float(mean[0]), float(mean[1]),
                             float(offspring[0][0])))
        # convergence
        if sigma * float(D.max()) < 1e-3:
            break
        if len(history_mean) >= 5:
            recent = [h[2] for h in history_mean[-5:]]
            if max(recent) - min(recent) < 1e-9:
                break
    return mean, history_mean, n_evals


# ----------------------------------------------------------------------
# 4. Driver
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=21,
                    help="grid cells per side")
    ap.add_argument("--ref",  type=int, nargs=4, default=[13, 16, -3, 3])
    ap.add_argument("--test", type=int, nargs=4, default=[16, 16,  0,  0])
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.7)
    ap.add_argument("--start", type=float, nargs=2, default=None,
                    help="starting τ (Re, Im).  Default = τ_test natural")
    ap.add_argument("--sigma0", type=float, default=0.25)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--cache", default=None,
                    help="path to grid cache pickle")
    args = ap.parse_args()

    Rx, Ry, RTx, RTy = args.ref
    Tx, Ty, TTx, TTy = args.test
    print(f"  ref geom : ({Rx}, {Ry}, {RTx}, {RTy})")
    print(f"  test geom: ({Tx}, {Ty}, {TTx}, {TTy})")
    w1_r, w2_r = cft.torus_periods(Rx, Ry, RTx, RTy)
    tau_r = w2_r / w1_r
    if tau_r.imag < 0: tau_r = complex(tau_r.real, -tau_r.imag)
    print(f"  τ_truth = {tau_r:.4f}")
    w1_t, w2_t = cft.torus_periods(Tx, Ty, TTx, TTy)
    tau_t_natural = w2_t / w1_t
    if tau_t_natural.imag < 0:
        tau_t_natural = complex(tau_t_natural.real, -tau_t_natural.imag)
    print(f"  τ_natural test = {tau_t_natural:.4f}")
    if args.start is None:
        tau_start = tau_t_natural
    else:
        tau_start = complex(args.start[0], args.start[1])
    print(f"  τ_start = {tau_start:.4f}")

    re_lo = tau_r.real - args.span_re / 2
    re_hi = tau_r.real + args.span_re / 2
    im_lo = max(0.05, tau_r.imag - args.span_im / 2)
    im_hi = tau_r.imag + args.span_im / 2
    print(f"  grid: Re ∈ [{re_lo:.3f}, {re_hi:.3f}], "
          f"Im ∈ [{im_lo:.3f}, {im_hi:.3f}], n={args.n}")

    cache_path = args.cache or os.path.join(
        OUT, f"grid_cache_n{args.n}_T{Tx}_{Ty}_{TTx}_{TTy}.pkl")

    print("  building ref Ising correlator…")
    ref_data = cft.make_ising_data_dict(Rx, Ry, RTx, RTy)
    print("  building/loading test grid cache…")
    cache = build_or_load_cache(cache_path, n=args.n,
                                Tx=Tx, Ty=Ty, TTx=TTx, TTy=TTy,
                                w1_ref=w1_r,
                                re_lo=re_lo, re_hi=re_hi,
                                im_lo=im_lo, im_hi=im_hi)

    variants = {
        "production":  cost_mod._l2_cost_test_native,
        **{f"zoo_{k}": v for k, v in cost_zoo.ZOO.items()},
    }
    print(f"  variants: {list(variants.keys())}")

    bounds_lo = np.array([re_lo, im_lo])
    bounds_hi = np.array([re_hi, im_hi])
    x_truth = np.array([tau_r.real, tau_r.imag])

    results = {}
    print()
    hdr = (f"  {'variant':<14}  {'seed':>4}  {'final τ':>16}  "
           f"{'dist→truth':>10}  {'cost→truth':>10}  "
           f"{'dist→start':>10}  {'#evals':>6}  status")
    print(hdr); print("  " + "-" * len(hdr))
    for name, fn in variants.items():
        runs = []
        for seed in range(args.seeds):
            ev = GridEvaluator(cache, ref_data, fn,
                               Tx=Tx, Ty=Ty, TTx=TTx, TTy=TTy,
                               Rx=Rx, Ry=Ry, RTx=RTx, RTy=RTy)
            try:
                xf, hist, nev = cmaes(
                    ev, [tau_start.real, tau_start.imag],
                    sigma0=args.sigma0, max_evals=args.max_evals,
                    popsize=args.popsize, seed=seed,
                    bounds_lo=bounds_lo, bounds_hi=bounds_hi)
            except Exception as e:
                print(f"  {name:<14}  seed {seed}  CRASH: {e}")
                continue
            d_truth = float(np.linalg.norm(xf - x_truth))
            d_start = float(np.hypot(tau_start.real - xf[0],
                                     tau_start.imag - xf[1]))
            # cost at truth (for reference)
            c_truth = ev(x_truth)
            ev.calls -= 1; ev.history.pop()
            status = "OK"
            if d_truth < 0.1:        status = "✓ found"
            elif d_truth < 0.3:      status = "near"
            else:                    status = "✗ lost"
            print(f"  {name:<14}  {seed:>4}  "
                  f"({xf[0]:+.3f},{xf[1]:+.3f})  "
                  f"{d_truth:>10.3f}  {c_truth:>10.2e}  "
                  f"{d_start:>10.3f}  {nev:>6d}  {status}")
            runs.append(dict(seed=seed, x_final=xf.tolist(),
                             dist_truth=d_truth,
                             dist_start=d_start, n_evals=nev,
                             status=status,
                             history=hist,
                             eval_history=ev.history))
        results[name] = runs

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    print("  PERFORMANCE SUMMARY (averaged over seeds):")
    print(f"  {'variant':<14}  {'mean dist':>10}  {'best dist':>10}  "
          f"{'success%':>9}  {'mean evals':>10}")
    print("  " + "-" * 60)
    rank_rows = []
    for name in variants:
        runs = results.get(name, [])
        if not runs:
            print(f"  {name:<14}  no runs"); continue
        dists = [r["dist_truth"] for r in runs]
        evals = [r["n_evals"] for r in runs]
        succ = 100.0 * sum(1 for r in runs if "✓" in r["status"]) / len(runs)
        print(f"  {name:<14}  {np.mean(dists):>10.3f}  "
              f"{np.min(dists):>10.3f}  {succ:>8.0f}%  "
              f"{np.mean(evals):>10.1f}")
        rank_rows.append((name, np.mean(dists), np.min(dists), succ))
    rank_rows.sort(key=lambda r: (-r[3], r[1]))   # success desc, dist asc
    print()
    print("  RANKING (success% desc, mean dist asc):")
    for i, (n, md, bd, s) in enumerate(rank_rows, 1):
        print(f"    {i}. {n:<14}  success={s:>3.0f}%  "
              f"mean_dist={md:.3f}  best_dist={bd:.3f}")

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    summary_path = os.path.join(OUT, "zoo_grid_cmaes_summary.json")
    with open(summary_path, "w") as f:
        json.dump(dict(
            tau_truth=[tau_r.real, tau_r.imag],
            tau_start=[tau_start.real, tau_start.imag],
            cache=cache_path,
            args=vars(args),
            results={k: [{kk: vv for kk, vv in r.items()
                          if kk not in ("history", "eval_history")}
                         for r in v]
                     for k, v in results.items()},
        ), f, indent=2)
    print(f"\n  summary → {os.path.relpath(summary_path, HERE)}")

    # ------------------------------------------------------------------
    # Plot trajectories overlaid on cost heatmap
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import colors as mcolors

        re_arr = np.array(cache["re"])
        im_arr = np.array(cache["im"])
        # Build heatmap once per variant from cached evaluations
        ncols = 3
        nrows = (len(variants) + ncols - 1) // ncols
        fig, axs = plt.subplots(nrows, ncols,
                                figsize=(5.2 * ncols, 4.6 * nrows),
                                squeeze=False)

        # Compute full heatmap per variant (re-evaluating costs from cache).
        for k, (name, fn) in enumerate(variants.items()):
            r, c = k // ncols, k % ncols
            ax = axs[r, c]
            grid_cost = np.full((args.n, args.n), np.nan)
            null = io.StringIO()
            for (i, j), d in cache["grid"].items():
                old = sys.stdout; sys.stdout = null
                try:
                    cval, *_ = fn(ref_data, d,
                                  Tx, Ty, TTx, TTy,
                                  Rx, Ry, RTx, RTy,
                                  copies=2, power=2)
                except Exception:
                    cval = float("nan")
                finally:
                    sys.stdout = old
                if np.isfinite(cval) and cval > 0:
                    grid_cost[i, j] = cval
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
            # truth & start
            ax.scatter([tau_r.real], [tau_r.imag], marker="*",
                       s=320, color="red", edgecolor="white",
                       lw=1.4, zorder=10, label="truth")
            ax.scatter([tau_start.real], [tau_start.imag], marker="o",
                       s=140, color="cyan", edgecolor="black",
                       lw=1.2, zorder=9, label="start")
            # trajectories
            colors = plt.cm.tab10(np.linspace(0, 1, args.seeds))
            for run, col in zip(results.get(name, []), colors):
                hist = np.array([h[:2] for h in run["history"]])
                if len(hist):
                    ax.plot(hist[:, 0], hist[:, 1], "-",
                            color=col, lw=1.2, alpha=0.9)
                    ax.scatter([hist[-1, 0]], [hist[-1, 1]],
                               marker="X", s=80, color=col,
                               edgecolor="black", lw=0.6, zorder=8)
            ax.set_xlim(re_lo, re_hi); ax.set_ylim(im_lo, im_hi)
            ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
            md = (np.mean([r["dist_truth"] for r in results.get(name, [])])
                  if results.get(name) else float("nan"))
            succ = (100.0 * sum(1 for r in results.get(name, [])
                                if "✓" in r["status"])
                    / max(1, len(results.get(name, []))))
            ax.set_title(f"{name}\nsuccess={succ:.0f}%  "
                         f"mean dist={md:.2f}", fontsize=10)
            ax.legend(loc="upper right", fontsize=7)
        for k in range(len(variants), nrows * ncols):
            axs[k // ncols, k % ncols].axis("off")
        plt.tight_layout()
        png = os.path.join(OUT, "zoo_grid_cmaes_trajectories.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"  plot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"  plot failed: {e}")


if __name__ == "__main__":
    main()
