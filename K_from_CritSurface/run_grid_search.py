#!/usr/bin/env python3
"""
run_grid_search.py — Grid search with pre-computed β_c surface
==============================================================

Uses a pre-built β_c(r₁, r₂) surface (from build_betac_surface.py) to
skip the expensive MC susceptibility scan at every grid point.  Each
evaluation only runs a production MC at the interpolated β_c, then scores
with Z² / Z⁴ boundary-slice cost functions.

Workflow:
  1. Run build_betac_surface.py to create the β_c surface (once per geometry).
  2. Edit the CONFIG section below.
  3. Run this script.

The β_c surface eliminates:
  - ~60% of compute time per grid point (no scan overhead)
  - All β_c finder jitter (the dominant noise source in the old approach)
"""

# ===========================================================================
# USER CONFIGURATION
# ===========================================================================

# --- Test lattice geometry ---
Lx = 32
Ly = 32
Tx = 0
Ty = 0

# --- Reference lattice geometry ---
REF_Lx = 32
REF_Ly = 32
REF_Tx = 0
REF_Ty = 0

# --- β_c surface file (from build_betac_surface.py) ---
BETAC_SURFACE = "betac_surfaces/surface_32x32_T0x0.json"

# --- Initial grid search parameters ---
R1_INIT    = 1.5
R2_INIT    = 1.5
HALF_SPAN  = 1.0

# --- MC statistics for production runs ---
N_TRAJ_PROD = 50000

# --- Grid search parameters ---
MAX_ITER = 6
N_GRID   = 5
TOP_N    = 3     # top Z² candidates re-scored with Z⁴ per level

# --- Output directory ---
OUTPUT_DIR = "results/equilateral_lut_from_1.5"

# --- Reference data ---
REFERENCE  = None     # None → auto-detect from K_from_TwoPoint_standalone/ref_data/
REF_N_TRAJ = 500000
REF_BETA_C = None     # e.g. 0.27465307216702745 for exact equilateral

# --- Simulator path ---
import sys as _sys
EXE = ("../K_from_TwoPoint_standalone/bin/ising_tri_twisted_parallelogram.exe"
       if _sys.platform == "win32"
       else "../K_from_TwoPoint_standalone/bin/ising_tri_twisted_parallelogram")

# ===========================================================================
# END OF USER CONFIGURATION
# ===========================================================================

import os, sys, json, time, warnings
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
_STANDALONE = os.path.join(_PARENT, "K_from_TwoPoint_standalone")
if _STANDALONE not in sys.path:
    sys.path.insert(0, _STANDALONE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
os.chdir(_HERE)

import optimise_couplings as oc
from betac_surface import BetacSurface

# ---------------------------------------------------------------------------
# Load β_c surface
# ---------------------------------------------------------------------------
if not os.path.isfile(BETAC_SURFACE):
    print(f"ERROR: β_c surface not found: {BETAC_SURFACE}")
    print("Run build_betac_surface.py first to generate the surface.")
    sys.exit(1)

surface = BetacSurface(BETAC_SURFACE)
surface.summary()
print()

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------
if REFERENCE is None:
    _ref_tag = f"ref_{REF_Lx}x{REF_Ly}" + (f"_t{REF_Tx}x{REF_Ty}" if REF_Tx or REF_Ty else "")
    REFERENCE = os.path.join(_STANDALONE, "ref_data", _ref_tag, "ref_metadata.json")
    _bundled = os.path.join(_STANDALONE, "ref_data", "ref_metadata.json")
    if (not os.path.isfile(REFERENCE) and os.path.isfile(_bundled)
            and REF_Lx == 32 and REF_Ly == 32 and REF_Tx == 0 and REF_Ty == 0):
        REFERENCE = _bundled

if not os.path.isfile(REFERENCE):
    print(f"ERROR: Reference not found: {REFERENCE}")
    print("Generate it first via K_from_TwoPoint_standalone (run_grid_search.py or gen_ref.py).")
    sys.exit(1)

with open(REFERENCE) as f:
    ref_meta = json.load(f)
ref_data = oc.load_all_to_all(ref_meta["a2a_path"])
print(f"Reference loaded: {REFERENCE}")
print(f"  β_c = {ref_meta['beta_c']:.10f}, Lx={ref_meta['Lx']}, Ly={ref_meta['Ly']}")

L_eff = min(Lx, Ly)
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Evaluate one grid point: LUT β_c → production → Z² cost
# ---------------------------------------------------------------------------
def evaluate_point_lut(r1, r2, level, i, j, data_dir):
    """Return (beta_c, z2_cost, z2_ndof, test_data_path)."""
    k3 = 1.0
    k1, k2 = r1 * k3, r2 * k3
    label = f"L{level}_i{i}_j{j}"

    level_dir = os.path.join(data_dir, f"level{level:02d}")
    os.makedirs(level_dir, exist_ok=True)

    plotter.set_current_point(r1, r2)

    # β_c from the pre-computed surface — instant, no MC scan, no jitter
    beta_c = surface(r1, r2)
    print(f"  [{label}] r1={r1:.6f} r2={r2:.6f}: β_c={beta_c:.8f} (LUT)",
          end="", flush=True)
    t0 = time.time()

    # Production run
    prod_dir = os.path.join(level_dir, f"prod_{label}")
    stdout, subdir = oc.run_simulator(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3, beta_c,
        n_traj=N_TRAJ_PROD, n_therm=3000, data_dir=prod_dir)
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    test_data = oc.load_all_to_all(a2a_path)

    # Z² cost (normed quadratic — most stable for grid search)
    chi2, ndof, chi2_err = oc.compute_chi2(
        ref_data, test_data,
        r_min=0.0, r_max_frac=0.33, L_eff=L_eff, Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
        ref_Lx=REF_Lx, ref_Ly=REF_Ly, ref_Tx=REF_Tx, ref_Ty=REF_Ty,
        cost="boundary_slices_normed",
    )
    z2_cost = chi2 / max(ndof, 1)
    t_total = time.time() - t0
    print(f"  Z²={z2_cost:.4f} ({t_total:.0f}s)")
    sys.stdout.flush()

    slices = None
    try:
        slices = oc.extract_boundary_slices(
            ref_data, test_data, Lx, Ly, Tx, Ty,
            ref_Lx=REF_Lx, ref_Ly=REF_Ly,
            ref_Tx=REF_Tx, ref_Ty=REF_Ty)
    except Exception:
        pass

    # Record for live plot (no scan data to show — pass empty lists)
    plotter.record_point(r1, r2, beta_c, chi2, ndof,
                         [], [], [],
                         parabola=None,
                         chi2_err=chi2_err, slices=slices)

    return beta_c, z2_cost, ndof, a2a_path


def rescore_z4(a2a_path):
    """Re-score an already-run production with Z⁴ (no extra MC needed)."""
    test_data = oc.load_all_to_all(a2a_path)
    chi2, ndof, chi2_err = oc.compute_chi2(
        ref_data, test_data,
        r_min=0.0, r_max_frac=0.33, L_eff=L_eff, Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
        ref_Lx=REF_Lx, ref_Ly=REF_Ly, ref_Tx=REF_Tx, ref_Ty=REF_Ty,
        cost="boundary_slices_normed_quartic",
    )
    return chi2 / max(ndof, 1)


# ---------------------------------------------------------------------------
# Grid search: LUT β_c + Z² grid / Z⁴ re-score
# ---------------------------------------------------------------------------
print(f"\n{'='*70}")
print(f"Grid Search with Pre-Computed β_c Surface")
print(f"  Lattice: {Lx}×{Ly}  Tx={Tx} Ty={Ty}")
print(f"  Reference: {REF_Lx}×{REF_Ly}  Tx={REF_Tx} Ty={REF_Ty}")
print(f"  Grid: {N_GRID}×{N_GRID}  Levels: {MAX_ITER}")
print(f"  Start: r1={R1_INIT}, r2={R2_INIT}, half_span={HALF_SPAN}")
print(f"  MC: prod={N_TRAJ_PROD} (no scan — β_c from surface)")
print(f"  Top-{TOP_N} Z⁴ re-scoring per level")
print(f"{'='*70}\n")

r1_center = R1_INIT
r2_center = R2_INIT
hs_r1 = HALF_SPAN
hs_r2 = HALF_SPAN

import matplotlib
matplotlib.use("Agg")
plotter = oc.LivePlotter(OUTPUT_DIR, cost="boundary_slices_normed")

all_results = []
run_dir = os.path.join(OUTPUT_DIR, "runs")

for level in range(MAX_ITER):
    r1_lo = max(r1_center - hs_r1, 0.1)
    r1_hi = r1_center + hs_r1
    r2_lo = max(r2_center - hs_r2, 0.1)
    r2_hi = r2_center + hs_r2

    r1_vals = np.linspace(r1_lo, r1_hi, N_GRID)
    r2_vals = np.linspace(r2_lo, r2_hi, N_GRID)

    print(f"\n{'='*70}")
    print(f"Level {level}: centre=({r1_center:.4f}, {r2_center:.4f})"
          f"  hs=({hs_r1:.4f}, {hs_r2:.4f})")
    print(f"  r1 ∈ [{r1_lo:.4f}, {r1_hi:.4f}]  r2 ∈ [{r2_lo:.4f}, {r2_hi:.4f}]")
    print(f"{'='*70}")
    sys.stdout.flush()

    plotter.start_level(level, r1_vals, r2_vals)

    # --- Phase 1: evaluate all grid points with Z² ---
    grid = {}
    for i, r1 in enumerate(r1_vals):
        for j, r2 in enumerate(r2_vals):
            beta_c, z2, ndof, a2a_path = evaluate_point_lut(
                r1, r2, level, i, j, run_dir)
            grid[(i, j)] = {
                "r1": r1, "r2": r2, "beta_c": beta_c,
                "z2": z2, "a2a_path": a2a_path,
            }
            all_results.append({
                "level": level, "i": i, "j": j,
                "r1": r1, "r2": r2, "beta_c": beta_c,
                "z2_cost": z2,
            })

    # --- Phase 2: re-score top N with Z⁴ ---
    ranked = sorted(grid.keys(), key=lambda k: grid[k]["z2"])
    top_keys = ranked[:TOP_N]

    print(f"\n  Re-scoring top {TOP_N} with Z⁴:")
    for rank, key in enumerate(top_keys):
        pt = grid[key]
        z4 = rescore_z4(pt["a2a_path"])
        pt["z4"] = z4
        print(f"    #{rank+1}: ({key[0]},{key[1]}) r1={pt['r1']:.4f} r2={pt['r2']:.4f}"
              f"  Z²={pt['z2']:.4f}  Z⁴={z4:.4f}")

    # Best by Z⁴ among top candidates
    best_key = min(top_keys, key=lambda k: grid[k].get("z4", grid[k]["z2"]))
    bi, bj = best_key
    best = grid[best_key]
    r1_best, r2_best = best["r1"], best["r2"]
    best_beta_c = best["beta_c"]
    best_z2 = best["z2"]
    best_z4 = best.get("z4", best_z2)

    print(f"\n  Level {level} best: ({bi},{bj})  r1={r1_best:.6f} r2={r2_best:.6f}"
          f"  Z²={best_z2:.4f}  Z⁴={best_z4:.4f}  β_c={best_beta_c:.8f}")

    # --- Border detection ---
    r1_border = (bi == 0) or (bi == N_GRID - 1)
    r2_border = (bj == 0) or (bj == N_GRID - 1)

    # --- Save level summary ---
    level_summary = {
        "level": level,
        "r1_center": r1_center, "r2_center": r2_center,
        "hs_r1": hs_r1, "hs_r2": hs_r2,
        "r1_vals": r1_vals.tolist(), "r2_vals": r2_vals.tolist(),
        "best_ij": [bi, bj],
        "best_r1": r1_best, "best_r2": r2_best,
        "best_z2": best_z2, "best_z4": best_z4,
        "best_beta_c": best_beta_c,
        "r1_border": r1_border, "r2_border": r2_border,
        "top_n_z4": [
            {"ij": list(k), "r1": grid[k]["r1"], "r2": grid[k]["r2"],
             "z2": grid[k]["z2"], "z4": grid[k].get("z4")}
            for k in top_keys
        ],
    }
    with open(os.path.join(OUTPUT_DIR, f"grid_level_{level:02d}.json"), "w") as f:
        json.dump(level_summary, f, indent=2)

    # --- Refine / translate ---
    # Centre on best; halve span only for interior axes.
    r1_center = r1_best
    if r1_border:
        print(f"  r1: border → translate to {r1_center:.4f} (span unchanged)")
    else:
        hs_r1 /= 2.0
        print(f"  r1: interior → refine, hs_r1={hs_r1:.4f}")

    r2_center = r2_best
    if r2_border:
        print(f"  r2: border → translate to {r2_center:.4f} (span unchanged)")
    else:
        hs_r2 /= 2.0
        print(f"  r2: interior → refine, hs_r2={hs_r2:.4f}")

    # Check if grid search will move outside the surface domain
    if not surface.domain_contains(r1_center - hs_r1, r2_center - hs_r2):
        print(f"  WARNING: next grid extends outside β_c surface domain")
    if not surface.domain_contains(r1_center + hs_r1, r2_center + hs_r2):
        print(f"  WARNING: next grid extends outside β_c surface domain")

    sys.stdout.flush()


# ---------------------------------------------------------------------------
# Final results
# ---------------------------------------------------------------------------
print(f"\n{'='*70}")
print("OPTIMISATION COMPLETE")
print(f"{'='*70}")

best_overall = min(all_results, key=lambda r: r["z2_cost"])
print(f"Best point: r1={best_overall['r1']:.6f} r2={best_overall['r2']:.6f}"
      f"  Z²={best_overall['z2_cost']:.6f}  β_c={best_overall['beta_c']:.8f}")

fit_result = {
    "best_r1": best_overall["r1"],
    "best_r2": best_overall["r2"],
    "best_k3": 1.0,
    "best_beta_c": best_overall["beta_c"],
    "best_z2_cost": best_overall["z2_cost"],
    "method": "LUT surface + Z²/Z⁴",
    "betac_surface": BETAC_SURFACE,
    "levels_completed": min(level + 1, MAX_ITER),
    "all_results": all_results,
}

with open(os.path.join(OUTPUT_DIR, "fit_result.json"), "w") as f:
    json.dump(fit_result, f, indent=2)

plotter.save_final()
print(f"\nResults saved: {OUTPUT_DIR}/fit_result.json")
print("Done.")
