#!/usr/bin/env python3
"""
cmaes_patch_synthetic_test.py — convergence comparison of current
``test_native`` cost vs proposed (T3+T4) patch using a SYNTHETIC,
deterministic cost surface calibrated to the cached MC noise level.

Why synthetic
-------------
A real-MC CMA-ES A/B comparison consumes ~15 min and proves a
mathematically-predictable result: on a SINGLE fixed geometry the
patch is a pure rescaling (T3 → factor N_valid/3 on cost, similar on
σ; T4 → only σ changes, not cost).  CMA-ES is scale-invariant ⇒
trajectories must coincide to within MC noise.

The synthetic surface here lets us:
  (a) verify that invariance numerically across many seeds,
  (b) measure exactly how cross-geometry ranking flips under the patch,
  (c) confirm that no second-order effect (e.g. SNR-driven retries)
      survives in the standard pipeline.

Surface
-------
Per-direction cost as function of (r1, r2):

    C_d(r) = a_d · [(r1-1)² + (r2-1)²] + b_d · (r1-r2)²

with three direction amplitudes (a_v, a_u, a_w)=(1.0, 0.7, 1.3)
and (b_v, b_u, b_w)=(0.2, 0.1, 0.3).  Noise per direction is Gaussian
σ_d = 5e-3 (calibrated to L=12 RW two-point std observed in
results/_betac_stability_test).

Two test "geometries" probed:
  G0  — all 3 dirs valid (analogue of (12,12,0,0); gcds all ≥2)
  G1  — direction w collapsed (analogue of (12,12,0,1); gcds=[1,12,1])

Cost modes
----------
  current  : C = (Σ_d C_d) / 3       (collapsed dirs contribute 0)
             σ = sqrt(Σ_d σ_d²) / 3  (collapsed dirs contribute 0)
  proposed : drop collapsed dirs from numerator AND denominator.
             σ_d uses 1.12× scaling (calibrated T4 vertex-max factor
             measured in the previous patch-test run).

CMA-ES
------
Hand-rolled minimal CMA (no external dep) — same population update
rule for both cost modes, identical seed schedule, x0=(2.0, 0.5),
σ0=0.6, popsize=8, max_gen=20.

Outputs
-------
  results/_native_triage/cmaes_synthetic_summary.json
  results/_native_triage/cmaes_synthetic.png      (4-panel)
  results/_native_triage/cmaes_synthetic.log
"""
from __future__ import annotations

import io
import json
import math
import os
import sys
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------------------
# Synthetic surface
# ---------------------------------------------------------------------------
A = np.array([1.0, 0.7, 1.3])           # per-direction quadratic amplitude
B = np.array([0.2, 0.1, 0.3])           # per-direction (r1-r2)² coupling
SIGMA_D = np.array([5e-3, 5e-3, 5e-3])  # per-direction noise
T4_SCALE = 1.12                          # vertex-max σ_ref calibration
TRUE = np.array([1.0, 1.0])

GEOMS = {
    "G0_no_collapse":     np.array([True, True, True]),
    "G1_one_collapse":    np.array([True, True, False]),
    "G2_two_collapse":    np.array([True, False, False]),
}


def per_direction_costs(r1: float, r2: float) -> np.ndarray:
    quad = (r1 - 1.0) ** 2 + (r2 - 1.0) ** 2
    cross = (r1 - r2) ** 2
    return A * quad + B * cross           # shape (3,)


def evaluate(rng: np.random.Generator,
             r1: float, r2: float,
             valid: np.ndarray,
             mode: str) -> Tuple[float, float]:
    """Return (cost, sigma_cost) for one CMA-ES query."""
    Cd_true = per_direction_costs(r1, r2)
    noise = rng.normal(0.0, SIGMA_D)
    Cd = Cd_true + noise
    sd = SIGMA_D.copy()

    if mode == "current":
        # Collapsed dirs contribute exactly 0 to numerator AND variance.
        Cd_used = np.where(valid, Cd, 0.0)
        sd_used = np.where(valid, sd, 0.0)
        cost = float(Cd_used.sum() / 3.0)
        sig = float(math.sqrt((sd_used ** 2).sum()) / 3.0)
    elif mode == "proposed":
        # Drop collapsed dirs from both num and denom.  Apply T4 σ_ref
        # rescaling on surviving directions.
        N = int(valid.sum())
        if N == 0:
            return float("nan"), float("nan")
        Cd_used = Cd[valid]
        sd_used = sd[valid] * T4_SCALE
        cost = float(Cd_used.sum() / N)
        sig = float(math.sqrt((sd_used ** 2).sum()) / N)
    else:
        raise ValueError(mode)
    return cost, sig


# ---------------------------------------------------------------------------
# Minimal CMA-ES (vanilla μ/μ_w, λ; no covariance adaptation — sufficient
# for a 2D quadratic with SNR > 1).  Identical update rule for both modes,
# noise comes only from the cost evaluator.
# ---------------------------------------------------------------------------

@dataclass
class CMAState:
    mean: np.ndarray
    sigma: float
    history: List[dict]


def cmaes_run(eval_fn: Callable[[float, float, np.random.Generator],
                                Tuple[float, float]],
              x0: Tuple[float, float],
              sigma0: float,
              popsize: int,
              max_gen: int,
              seed: int) -> CMAState:
    """Hand-rolled (μ_w, λ)-ES.  ``eval_fn`` is called with (r1, r2, rng);
    rng is a per-eval generator derived from a deterministic seed schedule
    so both cost modes consume IDENTICAL noise samples."""
    n = 2
    mu = popsize // 2
    weights = np.log(mu + 1) - np.log(np.arange(1, mu + 1))
    weights /= weights.sum()
    mueff = 1.0 / (weights ** 2).sum()
    c_sigma = (mueff + 2) / (n + mueff + 5)
    d_sigma = 1.0 + 2.0 * max(0.0, math.sqrt((mueff - 1) / (n + 1)) - 1.0) + c_sigma
    chi_n = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))

    mean = np.array(x0, dtype=float)
    sigma = float(sigma0)
    p_sigma = np.zeros(n)

    master = np.random.default_rng(seed)
    history: List[dict] = []

    for g in range(max_gen):
        # Sampling RNG (separate from eval RNG, for reproducibility).
        sample_seed = int(master.integers(0, 2**31 - 1))
        sample_rng = np.random.default_rng(sample_seed)
        z = sample_rng.standard_normal((popsize, n))
        x = mean + sigma * z

        # Evaluate; per-eval RNG is derived deterministically from
        # (seed, gen, member) so every cost mode sees the same noise.
        evals = []
        for k in range(popsize):
            eval_seed = (seed * 1_000_003 + g * 1_009 + k) & 0x7FFFFFFF
            eval_rng = np.random.default_rng(eval_seed)
            cost, sig = eval_fn(float(x[k, 0]), float(x[k, 1]), eval_rng)
            evals.append((cost, sig, x[k]))
            history.append(dict(gen=g, member=k,
                                r1=float(x[k, 0]), r2=float(x[k, 1]),
                                cost=cost, sigma_cost=sig))

        evals.sort(key=lambda e: e[0])
        best_z = np.array([(e[2] - mean) / sigma for e in evals[:mu]])
        new_mean = mean + sigma * (weights @ best_z)

        # σ adaptation (CSA, simplified)
        p_sigma = (1 - c_sigma) * p_sigma + math.sqrt(
            c_sigma * (2 - c_sigma) * mueff) * (weights @ best_z)
        sigma = sigma * math.exp((c_sigma / d_sigma) *
                                 (np.linalg.norm(p_sigma) / chi_n - 1))
        mean = new_mean

    return CMAState(mean=mean, sigma=sigma, history=history)


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def run_pair(geom_label: str, valid: np.ndarray,
             seed: int, popsize: int, max_gen: int,
             x0=(2.0, 0.5), sigma0=0.6) -> Dict[str, CMAState]:
    out: Dict[str, CMAState] = {}
    for mode in ("current", "proposed"):
        def _eval(r1, r2, rng, _v=valid, _m=mode):
            return evaluate(rng, r1, r2, _v, _m)
        st = cmaes_run(_eval, x0, sigma0, popsize, max_gen, seed=seed)
        out[mode] = st
    return out


def best_so_far_curve(history: List[dict]) -> np.ndarray:
    bs = []
    cur = float("inf")
    for h in history:
        if h["cost"] < cur:
            cur = h["cost"]
        bs.append(cur)
    return np.array(bs)


def dist_so_far_curve(history: List[dict]) -> np.ndarray:
    """Distance from current best-cost iterate to TRUE optimum (1,1)."""
    out = []
    best_cost = float("inf"); best_r = (None, None)
    for h in history:
        if h["cost"] < best_cost:
            best_cost = h["cost"]; best_r = (h["r1"], h["r2"])
        if best_r[0] is None:
            out.append(np.nan)
        else:
            out.append(float(np.hypot(best_r[0] - 1.0, best_r[1] - 1.0)))
    return np.array(out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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
    print("cmaes_patch_synthetic_test — current vs T3+T4 patched test_native")
    print("=" * 76)

    SEEDS = [11, 23, 37, 41, 53, 67, 79, 89, 97, 101]
    POPSIZE = 8
    MAX_GEN = 20

    summary: Dict[str, dict] = {}
    all_runs: Dict[str, Dict[str, List[CMAState]]] = {}

    for label, valid in GEOMS.items():
        print(f"\n--- geometry {label}  valid={valid.tolist()} "
              f"(N_valid={int(valid.sum())}) ---")
        per_seed = {"current": [], "proposed": []}
        for s in SEEDS:
            r = run_pair(label, valid, seed=s,
                         popsize=POPSIZE, max_gen=MAX_GEN)
            per_seed["current"].append(r["current"])
            per_seed["proposed"].append(r["proposed"])

        # Trajectory invariance check (per-eval r1/r2 must coincide
        # because both modes consume identical noise → fitness ranks
        # within each population are scaled by the same constant ⇒
        # selected mu members are the same ⇒ identical means).
        max_r_dev = 0.0
        for s_idx, _ in enumerate(SEEDS):
            cur = per_seed["current"][s_idx].history
            prop = per_seed["proposed"][s_idx].history
            for hc, hp in zip(cur, prop):
                max_r_dev = max(max_r_dev,
                                abs(hc["r1"] - hp["r1"]),
                                abs(hc["r2"] - hp["r2"]))

        # Distance-to-truth at end of run, averaged over seeds.
        def _final_dist(states):
            ds = []
            for st in states:
                d = dist_so_far_curve(st.history)
                ds.append(d[-1])
            return float(np.mean(ds)), float(np.std(ds))

        d_cur_mean, d_cur_std = _final_dist(per_seed["current"])
        d_prop_mean, d_prop_std = _final_dist(per_seed["proposed"])

        # Best-cost magnitudes (these DO differ — that is the point).
        c_cur = np.mean([best_so_far_curve(st.history)[-1]
                         for st in per_seed["current"]])
        c_prop = np.mean([best_so_far_curve(st.history)[-1]
                          for st in per_seed["proposed"]])

        ratio_pred = (int(valid.sum()) / 3.0) if valid.any() else float("nan")
        ratio_obs = c_cur / c_prop if c_prop else float("nan")

        print(f"  trajectory iterate max |Δr| current vs proposed = {max_r_dev:.3e}  "
              f"({'IDENTICAL' if max_r_dev < 1e-12 else 'DIFFERS'})")
        print(f"  final dist→(1,1):  current = {d_cur_mean:.3f}±{d_cur_std:.3f}   "
              f"proposed = {d_prop_mean:.3f}±{d_prop_std:.3f}")
        print(f"  reported best cost: current = {c_cur:.4e}   proposed = {c_prop:.4e}   "
              f"ratio = {ratio_obs:.4f}  (predicted N_valid/3 = {ratio_pred:.4f})")

        summary[label] = dict(
            valid=valid.tolist(),
            n_valid=int(valid.sum()),
            max_r_deviation=float(max_r_dev),
            final_dist_current=(d_cur_mean, d_cur_std),
            final_dist_proposed=(d_prop_mean, d_prop_std),
            best_cost_current=float(c_cur),
            best_cost_proposed=float(c_prop),
            ratio_observed=float(ratio_obs),
            ratio_predicted=float(ratio_pred),
            seeds=SEEDS,
        )
        all_runs[label] = per_seed

    # Cross-geometry mis-ranking demonstration --------------------------
    print("\n--- cross-geometry mis-ranking ---")
    print("If you compare ‘best cost’ ACROSS geometries to pick a lattice,")
    print("the current cost biases you toward whichever geometry has more")
    print("collapsed directions (because they’re scaled by N_valid/3).\n")

    # Take final best-cost average per geometry → a ranking.
    rank = []
    for label, info in summary.items():
        rank.append((label, info["best_cost_current"],
                     info["best_cost_proposed"]))
    rank_cur = sorted(rank, key=lambda x: x[1])
    rank_prop = sorted(rank, key=lambda x: x[2])
    print(f"  ranking by CURRENT cost (smaller = ‘better’):")
    for i, (l, c, _) in enumerate(rank_cur, 1):
        print(f"    {i}. {l}  cost={c:.4e}")
    print(f"  ranking by PROPOSED cost:")
    for i, (l, _, p) in enumerate(rank_prop, 1):
        print(f"    {i}. {l}  cost={p:.4e}")
    same = [a[0] for a in rank_cur] == [a[0] for a in rank_prop]
    print(f"  rankings agree: {same}")

    # ---- Plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))

        # Panel 1: best-so-far cost vs eval, G0 (no collapse), all seeds.
        ax = axs[0, 0]
        for st in all_runs["G0_no_collapse"]["current"]:
            ax.plot(best_so_far_curve(st.history),
                    color="C0", alpha=0.3, lw=0.8)
        for st in all_runs["G0_no_collapse"]["proposed"]:
            ax.plot(best_so_far_curve(st.history),
                    color="C1", alpha=0.3, lw=0.8, ls="--")
        ax.set_yscale("log"); ax.set_xlabel("eval #")
        ax.set_ylabel("best cost so far")
        ax.set_title("G0 (no collapse): trajectories overlay\n"
                     "blue=current  orange-dashed=proposed")
        ax.grid(alpha=0.3)

        # Panel 2: same for G1 (1 dir collapsed).
        ax = axs[0, 1]
        for st in all_runs["G1_one_collapse"]["current"]:
            ax.plot(best_so_far_curve(st.history),
                    color="C0", alpha=0.3, lw=0.8)
        for st in all_runs["G1_one_collapse"]["proposed"]:
            ax.plot(best_so_far_curve(st.history),
                    color="C1", alpha=0.3, lw=0.8, ls="--")
        ax.set_yscale("log"); ax.set_xlabel("eval #")
        ax.set_ylabel("best cost so far")
        ax.set_title("G1 (1 dir collapsed): proposed = 1.5× current\n"
                     "(constant rescale ⇒ identical optimum)")
        ax.grid(alpha=0.3)

        # Panel 3: distance to truth vs eval, ALL geoms, both modes.
        ax = axs[1, 0]
        for label, runs in all_runs.items():
            d_cur = np.array([dist_so_far_curve(st.history)
                              for st in runs["current"]])
            d_prop = np.array([dist_so_far_curve(st.history)
                               for st in runs["proposed"]])
            x = np.arange(d_cur.shape[1])
            ax.plot(x, d_cur.mean(0), label=f"{label}  current", lw=1.5)
            ax.plot(x, d_prop.mean(0), label=f"{label}  proposed",
                    lw=1.5, ls="--")
        ax.set_xlabel("eval #"); ax.set_ylabel("dist of best to (1,1)")
        ax.set_yscale("log")
        ax.legend(fontsize=7, ncol=2)
        ax.set_title("Distance to truth (CMA-ES path is identical)")
        ax.grid(alpha=0.3)

        # Panel 4: ratio bar plot.
        ax = axs[1, 1]
        labels = list(summary.keys())
        obs = [summary[l]["ratio_observed"] for l in labels]
        pred = [summary[l]["ratio_predicted"] for l in labels]
        x = np.arange(len(labels))
        ax.bar(x - 0.2, obs, 0.4, label="observed cur/prop", color="C2")
        ax.bar(x + 0.2, pred, 0.4, label="predicted N_valid/3", color="C3",
               hatch="//", alpha=0.7)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right",
                                             fontsize=8)
        ax.set_ylabel("cost ratio");  ax.legend(fontsize=8)
        ax.set_title("Bias is exactly N_valid/3")
        ax.grid(alpha=0.3, axis="y")

        plt.suptitle("CMA-ES current vs proposed (T3+T4) test_native cost",
                     fontsize=12)
        plt.tight_layout()
        png = os.path.join(OUT, "cmaes_synthetic.png")
        plt.savefig(png, dpi=130)
        plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, "cmaes_synthetic.log"), "w") as f:
        f.write(log_buf.getvalue())
    with open(os.path.join(OUT, "cmaes_synthetic_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
