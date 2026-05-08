#!/usr/bin/env python3
"""
viz_cmaes_huber.py — visualize CMA-ES converging to the truth on the
`huber_log` cost surface.

Produces:
  cmaes_huber_<scenario>_<noise>.png  — 4-panel figure
    (a) cost surface (log10) with CMA-ES trajectory & population overlay
    (b) zoom near the basin
    (c) cost vs evaluation #  (best-so-far)
    (d) σ (CMA-ES step size) vs generation, plus distance-to-truth
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import modular_data as md
from modular_testbed import build_or_load_cache, cmaes
from cost_selector import SelectorEvaluator

OUT_BASE = os.path.join(HERE, "results", "_native_triage")

PANEL = {
    "456":         (-0.901, 0.794),
    "equilateral": (-0.500, 0.866),
    "tall":        ( 0.000, 1.300),
    "edge_minus":  (-0.950, 0.600),
    "interior":    (-0.250, 0.700),
}
DEFAULT_STARTS = {
    "456":         ( 0.300, 0.500),
    "equilateral": ( 0.500, 0.500),
    "tall":        (-0.700, 0.700),
    "edge_minus":  ( 0.000, 1.100),
    "interior":    (-0.250, 1.150),
}


# Wrapper that records every evaluation (not just generation means).
class TracingEvaluator:
    def __init__(self, base):
        self.base = base
        self.points = []   # list of (re, im, cost)
    def __call__(self, x):
        c = self.base(x)
        self.points.append((float(x[0]), float(x[1]), float(c)))
        return c


def compute_cost_grid(refs_profile, cache, mode, noise, seed_base):
    """Evaluate cost on the cached τ-grid (no extra CFT calls)."""
    re = np.array(cache["re"]); im = np.array(cache["im"])
    H = np.full((len(im), len(re)), np.nan)
    for (i, j), prof in cache["grid"].items():
        if prof is None: continue
        if noise > 0:
            h = hash(("grid_viz", seed_base, i, j)) & 0xFFFFFFFF
            rng = np.random.default_rng(h)
            prof_eff = md.add_noise(prof, noise, rng)
        else:
            prof_eff = prof
        try:
            c, _ = md.threedir_cost(refs_profile, prof_eff, mode=mode)
        except Exception:
            c = np.nan
        if not np.isfinite(c) or c <= 0: c = np.nan
        H[j, i] = c
    return re, im, H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=list(PANEL), default="456")
    ap.add_argument("--mode", default="huber_log")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--n", type=int, default=80)
    ap.add_argument("--K-max", type=int, default=8)
    ap.add_argument("--w1", type=float, default=16.0)
    ap.add_argument("--span-re", type=float, default=2.0)
    ap.add_argument("--span-im", type=float, default=1.4)
    ap.add_argument("--center", type=float, nargs=2, default=[-0.45, 0.85])
    ap.add_argument("--sigma0", type=float, default=0.35)
    ap.add_argument("--max-evals", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--start", type=float, nargs=2, default=None)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    cre, cim = args.center
    re_lo = cre - args.span_re/2; re_hi = cre + args.span_re/2
    im_lo = max(0.05, cim - args.span_im/2)
    im_hi = cim + args.span_im/2

    truth = PANEL[args.scenario]
    start = tuple(args.start) if args.start else DEFAULT_STARTS[args.scenario]

    cache_path = os.path.join(OUT_BASE,
        f"modular_cache_n{args.n}_K{args.K_max}_w{int(args.w1)}.pkl")
    cache = build_or_load_cache(cache_path, n=args.n,
                                 re_lo=re_lo, re_hi=re_hi,
                                 im_lo=im_lo, im_hi=im_hi,
                                 K_max=args.K_max, w1_norm=args.w1)
    ref = md.make_threedir_profile(complex(*truth),
                                    w1_norm=args.w1, K_max=args.K_max)

    print(f"  scenario : {args.scenario}  truth = {truth}")
    print(f"  start    : {start}")
    print(f"  mode     : {args.mode}   noise σ_rel = {args.noise}")

    print("  computing cost grid …")
    t0 = time.perf_counter()
    re_g, im_g, H = compute_cost_grid(ref, cache, args.mode,
                                       args.noise, args.seed)
    print(f"    grid {H.shape} in {time.perf_counter()-t0:.1f}s")

    base_eval = SelectorEvaluator(cache, ref, args.mode,
                                   noise_sigma=args.noise,
                                   noise_seed_base=args.seed,
                                   frozen_noise=False, pmean=1)
    tracer = TracingEvaluator(base_eval)

    print("  running CMA-ES …")
    bounds_lo = np.array([re_lo, im_lo]); bounds_hi = np.array([re_hi, im_hi])
    xf, hist, nev = cmaes(tracer, [start[0], start[1]],
                           sigma0=args.sigma0, max_evals=args.max_evals,
                           popsize=args.popsize, seed=args.seed,
                           bounds_lo=bounds_lo, bounds_hi=bounds_hi)
    d_mod = md.modular_distance(complex(*xf), complex(*truth))
    print(f"  done. final τ = ({xf[0]:.4f}, {xf[1]:.4f})  "
          f"modular d = {d_mod:.4f}  evals = {nev}")

    pts = np.array(tracer.points)         # (N, 3)
    means = np.array(hist)                # (G, 3)
    pop_size = args.popsize

    # build OUT
    tag = (args.tag if args.tag else
           f"{args.scenario}_{args.mode}_n{int(args.noise*1000):03d}")
    OUT = os.path.join(OUT_BASE, f"viz_cmaes_{tag}")
    os.makedirs(OUT, exist_ok=True)

    # ----------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.4, 1])

    # ---- (a) full surface with trajectory ----
    ax = fig.add_subplot(gs[0, :2])
    Hp = np.where(np.isnan(H) | (H <= 0), np.nan, H)
    vmin = max(np.nanmin(Hp), 1e-12)
    vmax = np.nanmax(Hp)
    pcm = ax.pcolormesh(re_g, im_g, Hp,
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading="auto", cmap="viridis")
    plt.colorbar(pcm, ax=ax, label=f"{args.mode} cost (log)")
    # contours
    try:
        levels = np.geomspace(vmin, vmax, 10)
        ax.contour(re_g, im_g, Hp, levels=levels,
                    colors="white", linewidths=0.4, alpha=0.5)
    except Exception:
        pass

    # population evaluations colored by generation
    n_ev = len(pts)
    gens_per_pt = np.arange(n_ev) // pop_size
    sc = ax.scatter(pts[:, 0], pts[:, 1], c=gens_per_pt,
                    cmap="autumn", s=15, alpha=0.7,
                    edgecolors="k", linewidths=0.3)
    plt.colorbar(sc, ax=ax, label="generation",
                 fraction=0.04, pad=0.10)

    # generation-mean trajectory
    ax.plot(means[:, 0], means[:, 1], "-", color="cyan", lw=2,
            alpha=0.9, label="CMA-ES mean")
    ax.plot(*start, "X", color="white", markersize=14,
             markeredgecolor="k", markeredgewidth=1.5, label="start")
    ax.plot(*truth, "*", color="lime", markersize=22,
             markeredgecolor="k", markeredgewidth=1.5, label="truth")
    ax.plot(xf[0], xf[1], "o", color="red", markersize=10,
             markeredgecolor="k", label="final τ")
    ax.set_xlim(re_lo, re_hi); ax.set_ylim(im_lo, im_hi)
    ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
    ax.set_title(f"(a) cost surface + CMA-ES trajectory  "
                  f"[{args.mode}, σ={args.noise}, scenario={args.scenario}]")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.2, color="white")

    # ---- (b) zoom near truth ----
    ax = fig.add_subplot(gs[0, 2])
    z = 0.4
    zlo_re, zhi_re = truth[0]-z, truth[0]+z
    zlo_im, zhi_im = max(0.1, truth[1]-z), truth[1]+z
    mask_re = (re_g >= zlo_re) & (re_g <= zhi_re)
    mask_im = (im_g >= zlo_im) & (im_g <= zhi_im)
    Hz = Hp[np.ix_(mask_im, mask_re)]
    pcm = ax.pcolormesh(re_g[mask_re], im_g[mask_im], Hz,
                         norm=LogNorm(vmin=max(1e-12, np.nanmin(Hz)),
                                       vmax=np.nanmax(Hz)),
                         shading="auto", cmap="viridis")
    plt.colorbar(pcm, ax=ax, fraction=0.05)
    in_zoom = ((pts[:, 0] >= zlo_re) & (pts[:, 0] <= zhi_re) &
               (pts[:, 1] >= zlo_im) & (pts[:, 1] <= zhi_im))
    ax.scatter(pts[in_zoom, 0], pts[in_zoom, 1],
                c=gens_per_pt[in_zoom], cmap="autumn",
                s=18, alpha=0.7, edgecolors="k", linewidths=0.3)
    ax.plot(means[:, 0], means[:, 1], "-", color="cyan", lw=2)
    ax.plot(*truth, "*", color="lime", markersize=20,
             markeredgecolor="k")
    ax.plot(xf[0], xf[1], "o", color="red", markersize=10,
             markeredgecolor="k")
    ax.set_xlim(zlo_re, zhi_re); ax.set_ylim(zlo_im, zhi_im)
    ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
    ax.set_title("(b) zoom near truth")
    ax.grid(alpha=0.2, color="white")

    # ---- (c) best-so-far cost vs evals ----
    ax = fig.add_subplot(gs[1, 0])
    best = np.minimum.accumulate(pts[:, 2])
    ax.semilogy(np.arange(1, n_ev+1), pts[:, 2], ".", color="0.6",
                ms=4, label="per-eval")
    ax.semilogy(np.arange(1, n_ev+1), best, "-", color="C3", lw=1.8,
                label="best-so-far")
    ax.set_xlabel("evaluation #")
    ax.set_ylabel(f"{args.mode} cost")
    ax.set_title("(c) cost convergence")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    # ---- (d) distance-to-truth + step size σ ----
    ax = fig.add_subplot(gs[1, 1])
    d_per_gen = []
    sigmas = []
    # reconstruct dist of mean per generation
    for k in range(len(means)):
        d_per_gen.append(np.linalg.norm(
            np.array([means[k][0], means[k][1]]) - np.array(truth)))
    gens = np.arange(1, len(d_per_gen)+1)
    ax.semilogy(gens, d_per_gen, "-o", color="C0", lw=1.6, ms=4,
                label="‖mean − truth‖")
    ax.axhline(0.1, color="g", ls="--", alpha=0.5, label="success thresh 0.1")
    ax.set_xlabel("generation")
    ax.set_ylabel("distance"); ax.set_title("(d) basin convergence")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    # ---- (e) summary text panel ----
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    txt = [
        f"scenario     : {args.scenario}",
        f"truth τ      : ({truth[0]:.4f}, {truth[1]:.4f})",
        f"start  τ     : ({start[0]:.4f}, {start[1]:.4f})",
        f"final  τ     : ({xf[0]:.4f}, {xf[1]:.4f})",
        f"modular dist : {d_mod:.4f}",
        f"euclid dist  : "
        f"{np.linalg.norm(np.array(xf)-np.array(truth)):.4f}",
        f"evals        : {nev}",
        f"generations  : {len(means)}",
        f"final cost   : {pts[-1,2]:.3e}",
        f"min  cost    : {pts[:,2].min():.3e}",
        "",
        f"VERDICT: " + ("✓ FOUND TRUTH" if d_mod < 0.1
                       else "✗ MISSED" if d_mod > 0.3 else "≈ NEAR"),
    ]
    ax.text(0.0, 0.95, "\n".join(txt), transform=ax.transAxes,
             family="monospace", fontsize=10, va="top")
    ax.set_title("(e) run summary", loc="left")

    plt.tight_layout()
    png = os.path.join(OUT, f"cmaes_{tag}.png")
    plt.savefig(png, dpi=130); plt.close(fig)
    print(f"\n  → {os.path.relpath(png, HERE)}")


if __name__ == "__main__":
    main()
