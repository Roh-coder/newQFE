#!/usr/bin/env python3
"""
spectral_vs_native_test.py
==========================

Controlled comparison of two cost functions on cached MC data:

  * test_native  — current production cost.  Test correlator read at
                   integer lattice sites along the three boundary
                   period vectors; reference interpolated at matching
                   physical positions.  Aggregate = mean over 3 dirs of
                   mean of squared residuals.

  * spectral     — new cost.  Both correlators NDFTed onto a polar
                   momentum grid (n_r × n_θ); cost = mean |Δ G̃(k)|².
                   No interpolation anywhere; uses ALL all-to-all data.

Datasets (under results/_betac_stability_test/):
  * equil_12     : 10 replicates at (k1,k2,k3) = (1,1,1)·k_c  (matched)
  * aniso_r2_12  : 10 replicates at (2,1,1)·k_c               (mismatched, mild)
  * aniso_r5_12  : 10 replicates at (5,1,1)·k_c               (mismatched, strong)

All three groups share geometry Lx=Ly=12, Tx=Ty=0, so the cost function
sees only physics differences.

Tests
-----
1. **MC noise floor (matched)** — pair each equil replicate with every
   other equil replicate (45 pairs).  Cost should be small; std gives
   the irreducible MC noise of the cost.

2. **Discrimination (mild)** — equil vs aniso_r2 (10×10 = 100 pairs).

3. **Discrimination (strong)** — equil vs aniso_r5 (100 pairs).

4. **Z-score** = (mean_mismatched − mean_matched) / σ_matched.  Higher
   is better for an optimizer (the gradient of cost w.r.t. true r1,r2
   is more visible above the MC noise).

5. **Hyperparameter sweep** for spectral — vary (n_r, n_θ, k_max) and
   plot discrimination Z vs each.

Outputs
-------
  results/_cost_compare/
    summary.json               — per-pair raw numbers + discrimination
    cost_compare.png           — 4-panel summary
    spectral_hyperparam.png    — Z vs (n_r, n_θ, k_max)

Run
---
  python spectral_vs_native_test.py
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cost as cost_mod
import mc_engine

ROOT     = Path(__file__).resolve().parent
STAB_DIR = ROOT / "results" / "_betac_stability_test"
OUT_DIR  = ROOT / "results" / "_cost_compare"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GEOM = dict(Lx=12, Ly=12, Tx=0, Ty=0)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_group(group_name: str):
    """Load all 10 replicates from one (k1,k2,k3) group."""
    base = STAB_DIR / group_name
    if not base.is_dir():
        raise SystemExit(f"missing dataset: {base}")
    out = []
    for sub in sorted(base.iterdir()):
        # subdir layout:  3pass_t00/12x12_t0x0_k1000_1000_1000_500/two_point_all_to_all.dat
        cand = list(sub.glob("*/two_point_all_to_all.dat"))
        if not cand:
            continue
        out.append(mc_engine.load_all_to_all(str(cand[0])))
    print(f"[load] {group_name}: {len(out)} replicates")
    return out


# ---------------------------------------------------------------------------
# Cost wrappers
# ---------------------------------------------------------------------------

def _native(ref, test, power=2):
    c, s, _, _ = cost_mod.l2_cost(
        ref, test,
        GEOM["Lx"], GEOM["Ly"], GEOM["Tx"], GEOM["Ty"],
        GEOM["Lx"], GEOM["Ly"], GEOM["Tx"], GEOM["Ty"],
        cost_mode="test_native", cost_power=power)
    return c, s


def _spectral(ref, test, **spec_kw):
    c, s, _, _ = cost_mod.l2_cost(
        ref, test,
        GEOM["Lx"], GEOM["Ly"], GEOM["Tx"], GEOM["Ty"],
        GEOM["Lx"], GEOM["Ly"], GEOM["Tx"], GEOM["Ty"],
        cost_mode="spectral", spectral_kwargs=spec_kw)
    return c, s


def _pair_costs(group_a, group_b, *, same_group: bool, spec_kw=None,
                want_native: bool = True):
    """Compute (native, spectral) cost for every cross-pair."""
    if spec_kw is None:
        spec_kw = {}
    nat, spec = [], []
    if same_group:
        idx = list(itertools.combinations(range(len(group_a)), 2))
        pairs = [(group_a[i], group_a[j]) for i, j in idx]
    else:
        pairs = [(a, b) for a in group_a for b in group_b]
    for ref, tst in pairs:
        with _silence():
            if want_native:
                n, _ = _native(ref, tst, power=2)
                nat.append(n)
            s, _ = _spectral(ref, tst, **spec_kw)
            spec.append(s)
    return np.array(nat), np.array(spec)


class _silence:
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        return self
    def __exit__(self, *a):
        sys.stdout.close(); sys.stdout = self._stdout


# ---------------------------------------------------------------------------
# Headline test: matched vs mild vs strong mismatch
# ---------------------------------------------------------------------------

def _summary(name, arr):
    return dict(name=name, n=len(arr),
                mean=float(np.mean(arr)),
                std=float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
                median=float(np.median(arr)),
                min=float(np.min(arr)), max=float(np.max(arr)))


def _z_score(matched, mismatched):
    """(mean_mis − mean_matched) / σ_matched."""
    sd = np.std(matched, ddof=1) if len(matched) > 1 else 0.0
    if sd <= 0:
        return float("inf") if np.mean(mismatched) > np.mean(matched) else 0.0
    return float((np.mean(mismatched) - np.mean(matched)) / sd)


def main():
    t0 = time.time()
    eq  = _load_group("equil_12")
    a2  = _load_group("aniso_r2_12")
    a5  = _load_group("aniso_r5_12")

    print("\n=== Test 1: MC noise floor (matched, 45 self-pairs) ===")
    nat_m, spec_m = _pair_costs(eq, eq, same_group=True)

    print("\n=== Test 2: mild mismatch (equil vs aniso_r2, 100 pairs) ===")
    nat_2, spec_2 = _pair_costs(eq, a2, same_group=False)

    print("\n=== Test 3: strong mismatch (equil vs aniso_r5, 100 pairs) ===")
    nat_5, spec_5 = _pair_costs(eq, a5, same_group=False)

    # Headline summary table.
    rows = []
    for label, nat_arr, spec_arr in [
            ("matched (eq–eq)",       nat_m, spec_m),
            ("mismatch r2 (eq–a2)",   nat_2, spec_2),
            ("mismatch r5 (eq–a5)",   nat_5, spec_5)]:
        rows.append(dict(label=label,
                         native=_summary("native",   nat_arr),
                         spectral=_summary("spectral", spec_arr)))

    z_native_2 = _z_score(nat_m,  nat_2)
    z_native_5 = _z_score(nat_m,  nat_5)
    z_spec_2   = _z_score(spec_m, spec_2)
    z_spec_5   = _z_score(spec_m, spec_5)

    print("\n" + "="*78)
    print(f"{'method':12s} {'matched μ±σ':>22s} {'r2 μ':>14s} {'r5 μ':>14s}"
          f" {'Z(r2)':>8s} {'Z(r5)':>8s}")
    print("-"*78)
    print(f"{'native':12s} "
          f"{nat_m.mean():>10.3e} ± {nat_m.std(ddof=1):>8.2e}"
          f" {nat_2.mean():>14.3e} {nat_5.mean():>14.3e}"
          f" {z_native_2:>8.2f} {z_native_5:>8.2f}")
    print(f"{'spectral':12s} "
          f"{spec_m.mean():>10.3e} ± {spec_m.std(ddof=1):>8.2e}"
          f" {spec_2.mean():>14.3e} {spec_5.mean():>14.3e}"
          f" {z_spec_2:>8.2f} {z_spec_5:>8.2f}")
    print("="*78)

    # ------------------------------------------------------------------
    # Hyperparameter sweep (Test 4)
    # ------------------------------------------------------------------
    print("\n=== Test 4: spectral hyperparam sweep ===")
    sweep = []
    base = dict(n_r=6, n_theta=12, k_max=1.6)
    grid = {
        "n_r":     [3, 4, 6, 8, 12, 16],
        "n_theta": [6, 8, 12, 16, 24, 32],
        "k_max":   [0.6, 1.0, 1.6, 2.4, 3.2, math.pi],
    }
    for axis, vals in grid.items():
        print(f"  axis={axis}")
        for v in vals:
            kw = dict(base); kw[axis] = v
            _, sm  = _pair_costs(eq, eq, same_group=True,  spec_kw=kw, want_native=False)
            _, sd2 = _pair_costs(eq, a2, same_group=False, spec_kw=kw, want_native=False)
            _, sd5 = _pair_costs(eq, a5, same_group=False, spec_kw=kw, want_native=False)
            row = dict(axis=axis, value=v,
                       matched_mean=float(sm.mean()),
                       matched_std=float(sm.std(ddof=1)),
                       r2_mean=float(sd2.mean()),
                       r5_mean=float(sd5.mean()),
                       Z_r2=_z_score(sm, sd2),
                       Z_r5=_z_score(sm, sd5))
            sweep.append(row)
            print(f"    {axis}={v!s:<8s}  Z(r2)={row['Z_r2']:7.2f}"
                  f"  Z(r5)={row['Z_r5']:7.2f}")

    # ------------------------------------------------------------------
    # Save raw + summary
    # ------------------------------------------------------------------
    summary = dict(geom=GEOM, base_spec=base, headline=rows,
                   z_scores=dict(native_r2=z_native_2,
                                 native_r5=z_native_5,
                                 spectral_r2=z_spec_2,
                                 spectral_r5=z_spec_5),
                   hyperparam_sweep=sweep,
                   raw=dict(matched_native=nat_m.tolist(),
                            matched_spectral=spec_m.tolist(),
                            r2_native=nat_2.tolist(),
                            r2_spectral=spec_2.tolist(),
                            r5_native=nat_5.tolist(),
                            r5_spectral=spec_5.tolist()),
                   wall_time_s=time.time() - t0)
    with open(OUT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[save] {OUT_DIR/'summary.json'}")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    _plot_headline(nat_m, nat_2, nat_5, spec_m, spec_2, spec_5,
                   z_native_2, z_native_5, z_spec_2, z_spec_5)
    _plot_sweep(sweep)
    print(f"[done] wall {time.time()-t0:.1f}s")


def _plot_headline(nat_m, nat_2, nat_5, spec_m, spec_2, spec_5,
                   zn2, zn5, zs2, zs5):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # Panel A: cost distributions, native (log-y)
    ax = axes[0, 0]
    ax.boxplot([nat_m, nat_2, nat_5],
               tick_labels=["matched", "r2 mismatch", "r5 mismatch"],
               showmeans=True)
    ax.set_yscale("log")
    ax.set_title("Native (test-native L²) — cost distributions")
    ax.set_ylabel("cost"); ax.grid(alpha=0.3, which="both")

    # Panel B: same for spectral
    ax = axes[0, 1]
    ax.boxplot([spec_m, spec_2, spec_5],
               tick_labels=["matched", "r2 mismatch", "r5 mismatch"],
               showmeans=True)
    ax.set_yscale("log")
    ax.set_title("Spectral (NDFT |ΔG̃|²) — cost distributions")
    ax.set_ylabel("cost"); ax.grid(alpha=0.3, which="both")

    # Panel C: discrimination Z (higher is better)
    ax = axes[1, 0]
    x = np.arange(2)
    ax.bar(x - 0.18, [zn2, zn5], width=0.36, label="native",   color="C0")
    ax.bar(x + 0.18, [zs2, zs5], width=0.36, label="spectral", color="C1")
    ax.set_xticks(x); ax.set_xticklabels(["r2 mismatch", "r5 mismatch"])
    ax.set_ylabel("Z = (μ_mis − μ_match) / σ_match")
    ax.set_title("Discrimination Z-score  (higher = optimizer sees signal\n"
                 " above MC noise more easily)")
    ax.axhline(0, color="black", lw=0.6)
    for i, (zn, zs) in enumerate(zip([zn2, zn5], [zs2, zs5])):
        ax.text(i - 0.18, zn, f"{zn:.0f}", ha="center", va="bottom", fontsize=9)
        ax.text(i + 0.18, zs, f"{zs:.0f}", ha="center", va="bottom", fontsize=9)
    ax.legend(); ax.grid(alpha=0.3, axis="y")

    # Panel D: scatter of native vs spectral cost on the same pairs
    ax = axes[1, 1]
    ax.scatter(nat_m, spec_m, s=14, alpha=0.7, label="matched",      color="C2")
    ax.scatter(nat_2, spec_2, s=14, alpha=0.7, label="r2 mismatch",  color="C0")
    ax.scatter(nat_5, spec_5, s=14, alpha=0.7, label="r5 mismatch",  color="C3")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("native cost"); ax.set_ylabel("spectral cost")
    ax.set_title("Per-pair correlation (one dot = one ref/test pair)")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)

    fig.suptitle("Spectral vs test-native cost — controlled MC stability test\n"
                 f"L=12 equilateral, 10 replicates per group", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = OUT_DIR / "cost_compare.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[save] {out}")


def _plot_sweep(sweep):
    axes_names = ["n_r", "n_theta", "k_max"]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4))
    for ax, axis in zip(axs, axes_names):
        rows = [r for r in sweep if r["axis"] == axis]
        xs = [r["value"] for r in rows]
        ax.plot(xs, [r["Z_r2"] for r in rows], "o-", label="Z(r2)", color="C0")
        ax.plot(xs, [r["Z_r5"] for r in rows], "s-", label="Z(r5)", color="C3")
        ax.set_xlabel(axis); ax.set_ylabel("discrimination Z")
        ax.set_title(f"sweep over {axis}"); ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Spectral cost — discrimination Z vs hyperparameters",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = OUT_DIR / "spectral_hyperparam.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
