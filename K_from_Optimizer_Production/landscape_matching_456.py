#!/usr/bin/env python3
"""
landscape_matching_456.py

8×8 cost-landscape heatmap for the 4-5-6 matching observable, on TWO test
sizes (L=8, L=16) followed by a 1/L continuum extrapolation.

Cost (correct matching condition, see algo.md §4):

  For each test point (r1, r2, L):
    1. Build tiled LinearND interpolant on test all-to-all G_conn(m,n).
    2. Sample G(t_k) at t_k = k/8 for k=1..7 along each of the 3 cycles.
    3. Form Z_σ-invariant log observable:
           L(t) = log|G(t)|  −  <log|G|>_{t ∈ 1/8..7/8}
       (subtracting the geometric mean cancels the global multiplicative
       normalisation Z_σ, and is far less noisy than dividing by a single
       sample G(1/2).)
    4. Pair test cycle ↔ reference cycle by PHYSICAL SIDE (truth pairing
       from algo.md §2.2):
            ref c0 ↔ test c0   (side 5)
            ref c1 ↔ test c2   (side 6)
            ref c2 ↔ test c1   (side 4)
    5. Cost = Σ_paired Σ_k (L_test(t_k) − L_ref(t_k))^2
       (UN-weighted residual, summed over k=1..7 and over the 3 sides.)

Reference: high-stats _fss_456/ref/a4 = (52, 64, -12, 12), iso K=1.

Continuum extrapolation (per (r1, r2)):
    cost(L) ≈ cost_∞ + a / L
    cost_∞ = 2·cost(L=16) - cost(L=8)
    Negative values are clipped to 0 (FSS noise on cells already at
    cost_∞ ≈ 0).

Output:
    results/_landscape/_match_456/
        match_cost_L8.png
        match_cost_L16.png
        match_cost_continuum.png
        match_cost_grid.png        ← 1×3 panel (L=8, L=16, ∞) side-by-side
        match_cost.npz             ← raw arrays
"""
from __future__ import annotations

import os
import sys
import math
import pickle
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine
from cost import boundary_paths, _SQRT3_2, _tile_interp

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
TEST_SIZES = [8, 16]  # FSS extrapolation pair

# 8x8 subgrid (step 1.0) over r1,r2 ∈ {2,3,...,9}
R_GRID = np.arange(2.0, 9.5, 1.0)  # 2,3,4,5,6,7,8,9 → 8×8 = 64 points
TRUTH = (5.06523, 7.74293)

T_SAMPLES = np.array([k / 8 for k in range(1, 8)], float)   # t = 1/8 .. 7/8

# Cycle ↔ physical side  (algo.md §2.2)
REF_CYCLE_SIDE  = {0: 5, 1: 6, 2: 4}    # twisted 4-5-6 ref
TEST_CYCLE_SIDE = {0: 5, 1: 4, 2: 6}    # untwisted, anisotropic-k test
SIDE_PAIR = {s: (next(c for c, v in REF_CYCLE_SIDE.items()  if v == s),
                 next(c for c, v in TEST_CYCLE_SIDE.items() if v == s))
             for s in (4, 5, 6)}

# Reference: high-stats α=4
REF_PATH = os.path.join(_HERE, "results", "_fss_456",
                        "ref", "a4", "two_point_all_to_all.dat")
REF_GEOM = (52, 64, -12, 12)

OUT_DIR  = os.path.join(_HERE, "results", "_landscape", "_match_456")
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def _sample_L(test_data, geom, cycle_idx, t_arr):
    """Sample log|G(t)| - <log|G|> along one cycle (Zσ-invariant)."""
    iG    = _tile_interp(test_data, *geom, "conn", copies=2)
    paths = boundary_paths(*geom)
    dm, dn = paths[cycle_idx]
    ex, ey = dm + 0.5*dn, _SQRT3_2*dn
    xy_t = np.column_stack([t_arr * ex, t_arr * ey])
    G    = np.asarray(iG(xy_t), float)
    LG   = np.log(np.maximum(np.abs(G), 1e-30))
    return LG - LG.mean()


def _ref_R_table():
    """Pre-compute L_ref(t_k) = log|G_ref| - <log|G_ref|> for cycles 0,1,2."""
    ref_data = mc_engine.load_all_to_all(REF_PATH)
    return {c: _sample_L(ref_data, REF_GEOM, c, T_SAMPLES) for c in range(3)}


def _test_cost(test_data, Lx, Ly, ref_table):
    """Un-weighted matching cost in log-mean-subtracted observable.

        cost = Σ_{paired sides} Σ_k (L_test(t_k) − L_ref(t_k))²
    """
    geom = (Lx, Ly, 0, 0)
    cost = 0.0
    for side, (rc, tc) in SIDE_PAIR.items():
        L_ref  = ref_table[rc]
        L_test = _sample_L(test_data, geom, tc, T_SAMPLES)
        m = np.isfinite(L_test) & np.isfinite(L_ref)
        if not np.any(m):
            return np.nan
        cost += float(np.sum((L_test[m] - L_ref[m])**2))
    return cost


def _load_landscape_pkl(L, r1, r2):
    p = os.path.join(_HERE, "results", "_landscape",
                     f"Lx{L}_Ly{L}_Tx0_Ty0", "grid",
                     f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
    if not os.path.exists(p):
        return None
    with open(p, "rb") as f:
        payload = pickle.load(f)
    return payload["test_data"]


# ─────────────────────────────────────────────────────────────────────────────
# Main computation
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("Loading reference (high-stats 4-5-6 α=4)...")
    ref_table = _ref_R_table()
    print(f"  L_ref samples: {[ref_table[c].round(3).tolist() for c in range(3)]}")

    nR = len(R_GRID)
    cost = {L: np.full((nR, nR), np.nan, dtype=float) for L in TEST_SIZES}

    for L in TEST_SIZES:
        print(f"\nComputing cost on L={L} grid ({nR}×{nR})...")
        n_done = n_miss = 0
        for i, r1 in enumerate(R_GRID):
            for j, r2 in enumerate(R_GRID):
                td = _load_landscape_pkl(L, r1, r2)
                if td is None:
                    n_miss += 1
                    continue
                c = _test_cost(td, L, L, ref_table)
                cost[L][i, j] = c
                n_done += 1
        print(f"  L={L}: {n_done} cells, {n_miss} missing")

    # ── 1/L extrapolation:  cost(L) = c_∞ + a/L  →  c_∞ = 2 c16 − c8 ────────
    cost_inf_raw = 2.0 * cost[16] - cost[8]
    a_coef       = 8.0 * (cost[8] - cost[16])     # cost(L) = c∞ + a/L
    # Negative extrapolations indicate cells already at zero cost up to FSS
    # noise (truth basin); clip to 0 so the heatmap is monotone in mismatch.
    cost_inf = np.where(cost_inf_raw > 0, cost_inf_raw, 0.0)

    np.savez(os.path.join(OUT_DIR, "match_cost.npz"),
             R_GRID=R_GRID,
             cost_L8=cost[8], cost_L16=cost[16],
             cost_inf=cost_inf, cost_inf_raw=cost_inf_raw, a_coef=a_coef,
             truth=np.array(TRUTH))

    # ─────────────────────────────────────────────────────────────────────────
    # Plot
    # ─────────────────────────────────────────────────────────────────────────
    def _plot(arr, title, fname, cmap="viridis", vmin=None, vmax=None):
        fig, ax = plt.subplots(figsize=(7.5, 6.0))
        # arr indexed (i=r1, j=r2); pcolormesh wants X=r2, Y=r1 if origin lower
        # Convention: x-axis = r1, y-axis = r2 → transpose
        Z = arr.T   # now shape (r2, r1) → imshow with origin lower works
        if vmin is None: vmin = float(np.nanmin(Z))
        if vmax is None: vmax = float(np.nanmax(Z))
        im = ax.imshow(Z, origin="lower", aspect="equal",
                        extent=[R_GRID[0]-0.5, R_GRID[-1]+0.5,
                                R_GRID[0]-0.5, R_GRID[-1]+0.5],
                        cmap=cmap, vmin=vmin, vmax=vmax)
        # mark truth
        ax.plot(TRUTH[0], TRUTH[1], marker="*", ms=22, mfc="red",
                mec="white", mew=1.5, label="truth (5.07, 7.74)")
        # mark argmin (finite cells only)
        if np.any(np.isfinite(arr)):
            ij = np.unravel_index(np.nanargmin(arr), arr.shape)
            ax.plot(R_GRID[ij[0]], R_GRID[ij[1]], marker="o", ms=14,
                    mfc="none", mec="cyan", mew=2.5,
                    label=f"argmin ({R_GRID[ij[0]]:.0f},{R_GRID[ij[1]]:.0f})")
        # cell-value annotations
        for i, r1 in enumerate(R_GRID):
            for j, r2 in enumerate(R_GRID):
                v = arr[i, j]
                if np.isfinite(v):
                    norm = (v - vmin) / max(vmax - vmin, 1e-12)
                    fc = "white" if norm < 0.55 else "black"
                    ax.text(r1, r2, f"{v:.3f}", ha="center", va="center",
                            color=fc, fontsize=8)
        ax.set_xlabel(r"$r_1$ (≡ $K_1$)")
        ax.set_ylabel(r"$r_2$ (≡ $K_2$)")
        ax.set_title(title)
        ax.set_xticks(R_GRID); ax.set_yticks(R_GRID)
        ax.legend(loc="upper left", framealpha=0.85, fontsize=9)
        plt.colorbar(im, ax=ax, label=r"match cost  $\Sigma(\Delta\mathcal{L})^2$")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, fname), dpi=140)
        plt.close(fig)

    _plot(cost[8],   "Match cost  L=8   (4-5-6 vs untwisted)",   "match_cost_L8.png")
    _plot(cost[16],  "Match cost  L=16  (4-5-6 vs untwisted)",   "match_cost_L16.png")
    _plot(cost_inf,  "Match cost  L→∞  (1/L extrapolated)",      "match_cost_continuum.png")

    # 1×3 side-by-side
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.2))
    vmin = float(min(np.nanmin(cost[8]), np.nanmin(cost[16]),
                     np.nanmin(cost_inf)))
    # cap vmax at 99th pct so the basin is visible (extreme corners dominate)
    all_finite = np.concatenate([
        cost[8].ravel()[np.isfinite(cost[8].ravel())],
        cost[16].ravel()[np.isfinite(cost[16].ravel())],
        cost_inf.ravel()[np.isfinite(cost_inf.ravel())],
    ])
    vmax = float(np.percentile(all_finite, 95))
    for ax, arr, ttl in zip(axes, [cost[8], cost[16], cost_inf],
                            ["L=8", "L=16", "L→∞ (1/L extrap.)"]):
        Z = arr.T
        im = ax.imshow(Z, origin="lower", aspect="equal",
                       extent=[R_GRID[0]-0.5, R_GRID[-1]+0.5,
                               R_GRID[0]-0.5, R_GRID[-1]+0.5],
                       cmap="viridis", vmin=vmin, vmax=vmax)
        ax.plot(TRUTH[0], TRUTH[1], marker="*", ms=22, mfc="red",
                mec="white", mew=1.5)
        if np.any(np.isfinite(arr)):
            ij = np.unravel_index(np.nanargmin(arr), arr.shape)
            ax.plot(R_GRID[ij[0]], R_GRID[ij[1]], marker="o", ms=14,
                    mfc="none", mec="cyan", mew=2.5)
        for i, r1 in enumerate(R_GRID):
            for j, r2 in enumerate(R_GRID):
                v = arr[i, j]
                if np.isfinite(v):
                    norm = (v - vmin) / max(vmax - vmin, 1e-12)
                    fc = "white" if norm < 0.55 else "black"
                    ax.text(r1, r2, f"{v:.3f}", ha="center", va="center",
                            color=fc, fontsize=7)
        ax.set_xlabel(r"$r_1$"); ax.set_ylabel(r"$r_2$"); ax.set_title(ttl)
        ax.set_xticks(R_GRID); ax.set_yticks(R_GRID)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(
        "Continuum-limit matching cost landscape\n"
        f"red ★ = analytic truth (5.07, 7.74);  cyan ○ = grid argmin   "
        f"(8×8 grid, step=1.0, ref = 4-5-6 α=4)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(OUT_DIR, "match_cost_grid.png"), dpi=140)
    plt.close(fig)

    # ─────────────────────────────────────────────────────────────────────────
    # Summary print
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("Summary (top-5 lowest-cost grid cells per L):")
    for label, arr in [("L=8", cost[8]), ("L=16", cost[16]),
                       ("L→∞", cost_inf)]:
        flat = []
        for i, r1 in enumerate(R_GRID):
            for j, r2 in enumerate(R_GRID):
                if np.isfinite(arr[i, j]):
                    d_truth = math.hypot(r1-TRUTH[0], r2-TRUTH[1])
                    flat.append((arr[i, j], r1, r2, d_truth))
        flat.sort()
        print(f"\n{label}:")
        for c, r1, r2, d in flat[:5]:
            print(f"  ({r1:.1f}, {r2:.1f})   cost={c:9.4f}   d_truth={d:.2f}")

    # honest diagnostic: rank-correlation cost vs distance-to-truth
    print("\nPearson(cost, d_truth)  — production target ≪ -0.5 :")
    for label, arr in [("L=8", cost[8]), ("L=16", cost[16]),
                       ("L→∞", cost_inf)]:
        d_arr = np.array([[math.hypot(r1-TRUTH[0], r2-TRUTH[1])
                           for r2 in R_GRID] for r1 in R_GRID])
        m = np.isfinite(arr)
        rho = float(np.corrcoef(arr[m], d_arr[m])[0, 1])
        print(f"  {label}:  rho = {rho:+.3f}")

    print("\nWrote:")
    for fn in ["match_cost_L8.png", "match_cost_L16.png",
               "match_cost_continuum.png", "match_cost_grid.png",
               "match_cost.npz"]:
        print(" ", os.path.join(OUT_DIR, fn))


if __name__ == "__main__":
    main()
