#!/usr/bin/env python3
"""
run_optimizer_test.py — Run only Phase 2 (optimization) using the existing
critical surface and reference data from the quick test.

Includes a live progress monitor that prints a formatted table after each
grid point and a per-level heatmap saved to disk.

Usage:
    python run_optimizer_test.py
"""

import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import mc_engine as mc
from betac_surface import BetacSurface

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =========================================================================
# Configuration
# =========================================================================
QUICK_TEST_DIR = "results/quick_test"
SURFACE_PATH   = os.path.join(QUICK_TEST_DIR, "betac_surfaces",
                               "surface_8x8_T0x0.json")
REF_META_PATH  = os.path.join(QUICK_TEST_DIR, "ref_data", "ref_metadata.json")
OUTPUT_DIR     = "results/optimizer_test"

EXE = "bin/ising_tri_twisted_parallelogram"

# Test lattice (must match the surface)
TEST_Lx = 8
TEST_Ly = 8

# Surface bounds (from quick_test config)
SURFACE_R1_LO = 0.6;  SURFACE_R1_HI = 1.6
SURFACE_R2_LO = 0.6;  SURFACE_R2_HI = 1.6

# Optimizer settings
OPT_R1_INIT    = 1.0
OPT_R2_INIT    = 1.0
OPT_HALF_SPAN  = 0.30
OPT_N_GRID     = 5      # 5×5 = 25 per level (denser than quick_test's 3×3)
OPT_MAX_LEVELS = 4
OPT_TOP_N      = 3
OPT_N_TRAJ     = 10000

# Reference lattice info (read from metadata but needed for cost function)
REF_Lx = 8;  REF_Ly = 8;  REF_Tx = 0;  REF_Ty = 0


# =========================================================================
# Monitor — formatted table + level heatmap
# =========================================================================

class ProgressMonitor:
    """Live progress monitor for the optimization grid search."""

    def __init__(self, output_dir):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.global_best_z2 = float("inf")
        self.global_best_pt = None
        self.start_time = time.time()
        self.all_results = []   # across all levels
        self._level_results = []  # current level only

    def _elapsed(self):
        s = time.time() - self.start_time
        m, s = divmod(int(s), 60)
        h, m = divmod(m, 60)
        return f"{h:d}h{m:02d}m{s:02d}s" if h else f"{m:d}m{s:02d}s"

    # ----- level lifecycle -----

    def begin_level(self, level, r1_vals, r2_vals, r1_center, r2_center,
                    hs_r1, hs_r2):
        self._level = level
        self._r1_vals = r1_vals
        self._r2_vals = r2_vals
        self._level_results = []
        n = len(r1_vals) * len(r2_vals)
        print()
        print("=" * 78)
        print(f"  LEVEL {level}   centre=({r1_center:.4f}, {r2_center:.4f})"
              f"   half-span=({hs_r1:.4f}, {hs_r2:.4f})")
        print(f"  r1 ∈ [{r1_vals[0]:.4f}, {r1_vals[-1]:.4f}]   "
              f"r2 ∈ [{r2_vals[0]:.4f}, {r2_vals[-1]:.4f}]   "
              f"({n} pts)")
        print("=" * 78)
        print(f"  {'#':>4s}  {'r1':>8s}  {'r2':>8s}  {'β_c':>10s}  "
              f"{'Z²':>10s}  {'time':>6s}  {'best_Z²':>10s}  status")
        print(f"  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*10}  "
              f"{'-'*10}  {'-'*6}  {'-'*10}  {'-'*8}")
        sys.stdout.flush()

    def log_point(self, idx, total, r1, r2, beta_c, z2, dt, inside_hull):
        """Log one completed grid point."""
        is_best = z2 < self.global_best_z2
        if is_best:
            self.global_best_z2 = z2
            self.global_best_pt = (r1, r2, beta_c, z2)

        status = ""
        if is_best:
            status = "★ NEW BEST"
        elif not inside_hull:
            status = "(extrap)"

        rec = {"r1": r1, "r2": r2, "beta_c": beta_c, "z2": z2}
        self._level_results.append(rec)
        self.all_results.append(rec | {"level": self._level})

        print(f"  {idx:>4d}  {r1:8.5f}  {r2:8.5f}  {beta_c:10.7f}  "
              f"{z2:10.5f}  {dt:5.0f}s  {self.global_best_z2:10.5f}  "
              f"{status}")
        sys.stdout.flush()

    def end_level(self, level, best_r1, best_r2, best_z2, best_z4,
                  top_results):
        """Print level summary and save heatmap."""
        print()
        print(f"  Level {level} best:  r1={best_r1:.6f}  r2={best_r2:.6f}  "
              f"Z²={best_z2:.5f}  Z⁴={best_z4:.5f}")

        if len(top_results) > 1:
            print(f"  Top {len(top_results)} Z⁴ re-scores:")
            for rank, tr in enumerate(top_results):
                print(f"    #{rank+1}  r1={tr['r1']:.5f}  r2={tr['r2']:.5f}"
                      f"  Z²={tr['z2']:.5f}  Z⁴={tr['z4']:.5f}")

        print(f"  Elapsed: {self._elapsed()}")
        print(f"  Global best so far:  r1={self.global_best_pt[0]:.6f}  "
              f"r2={self.global_best_pt[1]:.6f}  "
              f"Z²={self.global_best_pt[3]:.5f}")
        sys.stdout.flush()

        # Heatmap
        self._save_heatmap(level)

    def _save_heatmap(self, level):
        r1v = np.array(sorted(set(r["r1"] for r in self._level_results)))
        r2v = np.array(sorted(set(r["r2"] for r in self._level_results)))
        Z = np.full((len(r2v), len(r1v)), np.nan)
        for r in self._level_results:
            i = np.searchsorted(r1v, r["r1"])
            j = np.searchsorted(r2v, r["r2"])
            if i < len(r1v) and j < len(r2v):
                Z[j, i] = r["z2"]

        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(Z, origin="lower", aspect="auto",
                       extent=[r1v[0], r1v[-1], r2v[0], r2v[-1]],
                       cmap="viridis_r")
        cb = fig.colorbar(im, ax=ax, label="Z² cost")
        # Mark the best point
        best = min(self._level_results, key=lambda r: r["z2"])
        ax.plot(best["r1"], best["r2"], "r*", ms=15, label="best")
        ax.set_xlabel("r₁")
        ax.set_ylabel("r₂")
        ax.set_title(f"Level {level}: Z² cost heatmap  "
                     f"(best Z²={best['z2']:.5f})")
        ax.legend()
        fig.tight_layout()
        path = os.path.join(self.output_dir, f"heatmap_level{level:02d}.png")
        fig.savefig(path, dpi=120)
        plt.close(fig)
        print(f"  Heatmap saved: {path}")

    def final_summary(self):
        print()
        print("=" * 78)
        print("  OPTIMISATION COMPLETE")
        print("=" * 78)
        if self.global_best_pt:
            r1, r2, bc, z2 = self.global_best_pt
            print(f"  Best:  r1 = {r1:.6f}")
            print(f"         r2 = {r2:.6f}")
            print(f"         β_c = {bc:.8f}")
            print(f"         Z²  = {z2:.6f}")
        print(f"  Total time: {self._elapsed()}")
        print(f"  Points evaluated: {len(self.all_results)}")
        print("=" * 78)

        # Save convergence plot
        self._save_convergence()
        sys.stdout.flush()

    def _save_convergence(self):
        """Plot running-best Z² vs evaluation number."""
        z2s = [r["z2"] for r in self.all_results]
        running_best = np.minimum.accumulate(z2s)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(1, len(z2s) + 1), z2s, "o", ms=3, alpha=0.5,
                label="Z² per point")
        ax.plot(range(1, len(z2s) + 1), running_best, "r-", lw=2,
                label="running best")
        # Mark level boundaries
        levels = sorted(set(r["level"] for r in self.all_results))
        idx = 0
        for lv in levels:
            n = sum(1 for r in self.all_results if r["level"] == lv)
            if idx > 0:
                ax.axvline(idx + 0.5, color="gray", ls="--", lw=0.8)
                ax.text(idx + 1, ax.get_ylim()[1] * 0.95, f"L{lv}",
                        fontsize=8, va="top")
            idx += n
        ax.set_xlabel("Evaluation #")
        ax.set_ylabel("Z² cost")
        ax.set_title("Optimizer convergence")
        ax.legend()
        fig.tight_layout()
        path = os.path.join(self.output_dir, "convergence.png")
        fig.savefig(path, dpi=120)
        plt.close(fig)
        print(f"  Convergence plot: {path}")


# =========================================================================
# Optimizer (Phase 2 only, with monitor)
# =========================================================================

def run_optimization():
    t0_wall = time.time()

    # Load existing surface and reference
    surface = BetacSurface(SURFACE_PATH)
    surface.summary()

    with open(REF_META_PATH) as f:
        ref_meta = json.load(f)
    ref_data = mc.load_all_to_all(ref_meta["a2a_path"])
    print(f"Reference loaded: β_c={ref_meta['beta_c']:.10f}")

    opt_dir = OUTPUT_DIR
    run_dir = os.path.join(opt_dir, "runs")
    os.makedirs(run_dir, exist_ok=True)

    monitor = ProgressMonitor(opt_dir)

    r1_center = OPT_R1_INIT
    r2_center = OPT_R2_INIT
    hs_r1 = OPT_HALF_SPAN
    hs_r2 = OPT_HALF_SPAN

    for level in range(OPT_MAX_LEVELS):
        r1_lo = max(r1_center - hs_r1, SURFACE_R1_LO)
        r1_hi = min(r1_center + hs_r1, SURFACE_R1_HI)
        r2_lo = max(r2_center - hs_r2, SURFACE_R2_LO)
        r2_hi = min(r2_center + hs_r2, SURFACE_R2_HI)

        r1_vals = np.linspace(r1_lo, r1_hi, OPT_N_GRID)
        r2_vals = np.linspace(r2_lo, r2_hi, OPT_N_GRID)

        monitor.begin_level(level, r1_vals, r2_vals,
                            r1_center, r2_center, hs_r1, hs_r2)

        grid = {}
        idx = 0
        total = OPT_N_GRID * OPT_N_GRID

        for i, r1 in enumerate(r1_vals):
            for j, r2 in enumerate(r2_vals):
                idx += 1
                k3 = 1.0
                k1, k2 = r1 * k3, r2 * k3

                beta_c = surface(r1, r2)
                inside = surface.domain_contains(r1, r2)

                t_pt = time.time()
                label = f"L{level}_i{i}_j{j}"
                prod_dir = os.path.join(run_dir, f"level{level:02d}",
                                        f"prod_{label}")
                try:
                    stdout, subdir = mc.run_simulator(
                        EXE, TEST_Lx, TEST_Ly, 0, 0, k1, k2, k3, beta_c,
                        n_traj=OPT_N_TRAJ, n_therm=3000,
                        data_dir=prod_dir)
                except Exception as e:
                    print(f"  {idx:>4d}  r1={r1:.5f} r2={r2:.5f}  FAILED: {e}")
                    continue

                a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
                test_data = mc.load_all_to_all(a2a_path)

                chi2, ndof, chi2_err = mc.boundary_slices_normed(
                    ref_data, test_data,
                    TEST_Lx, TEST_Ly, 0, 0,
                    REF_Lx, REF_Ly, REF_Tx, REF_Ty)
                z2_cost = chi2 / max(ndof, 1)
                dt = time.time() - t_pt

                monitor.log_point(idx, total, r1, r2, beta_c, z2_cost,
                                  dt, inside)

                grid[(i, j)] = {
                    "r1": r1, "r2": r2, "beta_c": beta_c,
                    "z2": z2_cost, "a2a_path": a2a_path,
                }

        if not grid:
            print("  WARNING: no points evaluated!")
            continue

        # Re-score top N with Z⁴
        ranked = sorted(grid.keys(), key=lambda k: grid[k]["z2"])
        top_keys = ranked[:OPT_TOP_N]

        top_results = []
        for key in top_keys:
            pt = grid[key]
            test_data = mc.load_all_to_all(pt["a2a_path"])
            chi2_q, ndof_q, _ = mc.boundary_slices_quartic(
                ref_data, test_data,
                TEST_Lx, TEST_Ly, 0, 0,
                REF_Lx, REF_Ly, REF_Tx, REF_Ty)
            z4 = chi2_q / max(ndof_q, 1)
            pt["z4"] = z4
            top_results.append({"r1": pt["r1"], "r2": pt["r2"],
                                "z2": pt["z2"], "z4": z4})

        best_key = min(top_keys,
                       key=lambda k: grid[k].get("z4", grid[k]["z2"]))
        best = grid[best_key]

        monitor.end_level(level, best["r1"], best["r2"],
                          best["z2"], best.get("z4", best["z2"]),
                          top_results)

        # Save level JSON
        level_summary = {
            "level": level,
            "r1_center": float(r1_center), "r2_center": float(r2_center),
            "hs_r1": float(hs_r1), "hs_r2": float(hs_r2),
            "best_r1": float(best["r1"]), "best_r2": float(best["r2"]),
            "best_z2": float(best["z2"]),
            "best_z4": float(best.get("z4", best["z2"])),
            "best_beta_c": float(best["beta_c"]),
            "grid": {f"{k[0]}_{k[1]}": {
                "r1": v["r1"], "r2": v["r2"],
                "beta_c": v["beta_c"], "z2": v["z2"],
            } for k, v in grid.items()},
        }
        with open(os.path.join(opt_dir, f"grid_level_{level:02d}.json"),
                  "w") as f:
            json.dump(level_summary, f, indent=2)

        # --- Translate or refine ---
        bi, bj = best_key
        r1_center = float(best["r1"])
        r2_center = float(best["r2"])

        r1_border = (bi == 0) or (bi == OPT_N_GRID - 1)
        r2_border = (bj == 0) or (bj == OPT_N_GRID - 1)

        if r1_border:
            print(f"  r1: border → translate to {r1_center:.4f}")
        else:
            hs_r1 /= 2.0
            print(f"  r1: interior → refine, hs_r1={hs_r1:.4f}")

        if r2_border:
            print(f"  r2: border → translate to {r2_center:.4f}")
        else:
            hs_r2 /= 2.0
            print(f"  r2: interior → refine, hs_r2={hs_r2:.4f}")

    # Done
    monitor.final_summary()

    # Save final result
    fit = {
        "best_r1": monitor.global_best_pt[0],
        "best_r2": monitor.global_best_pt[1],
        "best_beta_c": monitor.global_best_pt[2],
        "best_z2": monitor.global_best_pt[3],
        "levels": OPT_MAX_LEVELS,
        "grid_size": OPT_N_GRID,
        "n_traj": OPT_N_TRAJ,
        "all_results": monitor.all_results,
    }
    with open(os.path.join(opt_dir, "fit_result.json"), "w") as f:
        json.dump(fit, f, indent=2, default=float)
    print(f"  Results saved: {opt_dir}/fit_result.json")


if __name__ == "__main__":
    run_optimization()
