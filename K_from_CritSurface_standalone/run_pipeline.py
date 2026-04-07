#!/usr/bin/env python3
"""
run_pipeline.py — Standalone 3-phase pipeline for finding anisotropic
couplings that reproduce a twisted equilateral lattice's two-point function.

  Phase 0:  Generate reference data on the twisted equilateral lattice.
  Phase 1:  Build β_c(r₁, r₂) surface on the untwisted test lattice.
  Phase 2:  Adaptive grid search using the surface as a LUT.

Edit the CONFIG section below, then run:

    python run_pipeline.py

The script is fully self-contained: it builds the C++ simulator if needed,
generates the reference if needed, builds the surface if needed, and then
runs the optimization.  Intermediate results are saved incrementally so
you can interrupt and resume at any point.

Designed to run overnight on a local PC.
"""

# ===========================================================================
# USER CONFIGURATION
# ===========================================================================

# ---------------------------------------------------------------------------
# REFERENCE LATTICE (equilateral triangles, possibly twisted)
# This defines the "target" physics you want to match.
# k1=k2=k3 (equilateral); twist via (Tx, Ty).
# ---------------------------------------------------------------------------
REF_Lx = 32          # reference lattice x-size
REF_Ly = 32          # reference lattice y-size
REF_Tx = 0           # reference twist x (e.g. 8 for 45-deg on 32×32)
REF_Ty = 0           # reference twist y

# Override the reference beta_c.  When set, the susceptibility scan is
# skipped for the reference lattice.  For the equilateral triangular
# Ising model on an infinite lattice: ln(3)/4 ≈ 0.27465307216702745.
# For finite lattices the pseudo-critical beta is slightly shifted;
# set to None to let the scan find it automatically.
REF_BETA_C = None     # e.g. 0.27465307216702745

# MC trajectories for the reference production run.
REF_N_TRAJ = 500000

# MC statistics for the reference beta_c scan (if REF_BETA_C is None).
REF_N_TRAJ_SCAN_COARSE = 30000
REF_N_TRAJ_SCAN_FINE   = 60000

# ---------------------------------------------------------------------------
# TEST LATTICE (untwisted, anisotropic couplings)
# This is the lattice you're optimizing over (r1, r2).
# k1 = r1*k3, k2 = r2*k3, k3 = 1.0;  Tx=Ty=0 always.
# ---------------------------------------------------------------------------
TEST_Lx = 32         # test lattice x-size
TEST_Ly = 32         # test lattice y-size

# ---------------------------------------------------------------------------
# SCANNING REGION for the β_c surface (Phase 1)
# The surface is built on a regular grid over this region.
# Make sure it covers the optimization search region.
# ---------------------------------------------------------------------------
SURFACE_R1_LO = 0.5
SURFACE_R1_HI = 2.5
SURFACE_R2_LO = 0.5
SURFACE_R2_HI = 2.5
SURFACE_N_R1  = 11   # grid density (11×11 = 121 scan points)
SURFACE_N_R2  = 11

# MC statistics for each surface scan point.
SURFACE_N_TRAJ_COARSE = 30000
SURFACE_N_TRAJ_FINE   = 60000

# ---------------------------------------------------------------------------
# OPTIMIZATION (Phase 2) — adaptive grid search
# ---------------------------------------------------------------------------
OPT_R1_INIT   = 1.0   # initial grid centre r1
OPT_R2_INIT   = 1.0   # initial grid centre r2
OPT_HALF_SPAN = 0.50   # half-width of grid in r-space
OPT_N_GRID    = 5      # grid points per side per level (5×5 = 25)
OPT_MAX_LEVELS = 6     # max refinement levels
OPT_TOP_N      = 3     # re-score top N with Z⁴

# Production MC statistics per grid point.
OPT_N_TRAJ_PROD = 50000

# Initial guess for beta_c (used by the surface builder as starting point).
BETA_INIT = 0.27

# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------
OUTPUT_DIR = "results/run_001"

# ---------------------------------------------------------------------------
# Path to the compiled simulator (auto-built if missing).
# ---------------------------------------------------------------------------
import sys as _sys
EXE = ("bin/ising_tri_twisted_parallelogram.exe"
       if _sys.platform == "win32"
       else "bin/ising_tri_twisted_parallelogram")

# ===========================================================================
# END OF USER CONFIGURATION
# ===========================================================================

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


# ===================================================================
# Auto-build simulator
# ===================================================================

def _build_exe():
    if os.path.isfile(EXE):
        return
    print(f"Simulator not found at '{EXE}' — building...")
    import subprocess
    for make_cmd in (["mingw32-make"] if sys.platform == "win32" else []) + ["make"]:
        try:
            r = subprocess.run([make_cmd], cwd=_HERE, capture_output=True, text=True)
            if r.returncode == 0 and os.path.isfile(EXE):
                print(f"Build successful via {make_cmd}")
                return
        except FileNotFoundError:
            pass
    # Fallback: direct g++ invocation
    src = os.path.join(_HERE, "src", "ising_tri_twisted_parallelogram.cc")
    inc = os.path.join(_HERE, "include")
    os.makedirs(os.path.join(_HERE, "bin"), exist_ok=True)
    cmd = ["g++", "-std=c++14", "-O3", "-Wall", "-Wno-sign-compare",
           "-Wno-unused-variable", f"-I{inc}", src, "-o",
           os.path.join(_HERE, EXE)]
    r = subprocess.run(cmd, cwd=_HERE, capture_output=True, text=True)
    if r.returncode == 0 and os.path.isfile(EXE):
        print(f"Build successful via g++")
        return
    print(f"Build failed:\n{r.stderr}")
    sys.exit(1)


# ===================================================================
# Phase 0: Generate reference data
# ===================================================================

def phase0_reference(output_dir):
    """Generate reference data on the equilateral twisted lattice."""
    ref_dir = os.path.join(output_dir, "ref_data")
    meta_path = os.path.join(ref_dir, "ref_metadata.json")

    if os.path.isfile(meta_path):
        print(f"Phase 0: Reference already exists at {meta_path}")
        with open(meta_path) as f:
            return json.load(f)

    os.makedirs(ref_dir, exist_ok=True)
    print(f"\n{'='*70}")
    print("Phase 0: Generating reference data")
    print(f"  Lattice: {REF_Lx}×{REF_Ly}  Tx={REF_Tx} Ty={REF_Ty}")
    print(f"  Equilateral couplings: k1=k2=k3=1.0")
    print(f"{'='*70}\n")

    # Find or use beta_c
    if REF_BETA_C is not None:
        beta_c = REF_BETA_C
        print(f"  Using provided β_c = {beta_c:.10f}")
    else:
        print(f"  Scanning for β_c ...")
        margin = max(0.20 * BETA_INIT, 0.04)
        beta_lo = max(0.01, BETA_INIT - margin)
        beta_hi = BETA_INIT + margin
        scan_dir = os.path.join(ref_dir, "beta_scan")
        beta_c, chi_peak, _, _, _ = mc.find_beta_c(
            EXE, REF_Lx, REF_Ly, REF_Tx, REF_Ty,
            k1=1.0, k2=1.0, k3=1.0,
            beta_lo=beta_lo, beta_hi=beta_hi,
            n_traj_coarse=REF_N_TRAJ_SCAN_COARSE,
            n_traj_fine=REF_N_TRAJ_SCAN_FINE,
            data_dir=scan_dir,
        )
        print(f"  Found β_c = {beta_c:.10f}  (χ_peak = {chi_peak:.6f})")

    # Production run at beta_c
    print(f"\n  Running production ({REF_N_TRAJ} trajectories) ...")
    t0 = time.time()
    prod_dir = os.path.join(ref_dir, "ref_production")
    stdout, subdir = mc.run_simulator(
        EXE, REF_Lx, REF_Ly, REF_Tx, REF_Ty,
        k1=1.0, k2=1.0, k3=1.0, beta=beta_c,
        n_traj=REF_N_TRAJ, n_therm=5000, data_dir=prod_dir)
    t_prod = time.time() - t0
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    print(f"  Production complete ({t_prod:.0f}s)")
    print(f"  Data: {a2a_path}")

    meta = {
        "beta_c": beta_c,
        "Lx": REF_Lx, "Ly": REF_Ly, "Tx": REF_Tx, "Ty": REF_Ty,
        "k1": 1.0, "k2": 1.0, "k3": 1.0,
        "n_traj": REF_N_TRAJ,
        "a2a_path": os.path.relpath(a2a_path, _HERE),
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"  Metadata saved: {meta_path}")
    return meta


# ===================================================================
# Phase 1: Build β_c surface
# ===================================================================

def phase1_surface(output_dir):
    """Build the β_c(r1, r2) surface on the untwisted test lattice."""
    surface_dir = os.path.join(output_dir, "betac_surfaces")
    os.makedirs(surface_dir, exist_ok=True)
    tag = f"surface_{TEST_Lx}x{TEST_Ly}_T0x0"
    surface_path = os.path.join(surface_dir, f"{tag}.json")

    # Check for existing complete surface
    if os.path.isfile(surface_path):
        with open(surface_path) as f:
            raw = json.load(f)
        expected = SURFACE_N_R1 * SURFACE_N_R2
        if len(raw.get("points", [])) >= expected:
            print(f"Phase 1: Surface already complete at {surface_path}")
            return surface_path

    print(f"\n{'='*70}")
    print("Phase 1: Building β_c surface")
    print(f"  Test lattice: {TEST_Lx}×{TEST_Ly}  Tx=0 Ty=0")
    print(f"  Scanning: r1 ∈ [{SURFACE_R1_LO}, {SURFACE_R1_HI}]  "
          f"r2 ∈ [{SURFACE_R2_LO}, {SURFACE_R2_HI}]")
    print(f"  Grid: {SURFACE_N_R1}×{SURFACE_N_R2} = "
          f"{SURFACE_N_R1 * SURFACE_N_R2} points")
    print(f"{'='*70}\n")

    r1_vals = np.linspace(SURFACE_R1_LO, SURFACE_R1_HI, SURFACE_N_R1)
    r2_vals = np.linspace(SURFACE_R2_LO, SURFACE_R2_HI, SURFACE_N_R2)

    # Load any previously completed points (incremental resume)
    completed = {}
    if os.path.isfile(surface_path):
        with open(surface_path) as f:
            raw = json.load(f)
        for p in raw.get("points", []):
            key = (round(p["r1"], 8), round(p["r2"], 8))
            completed[key] = p["beta_c"]

    total = SURFACE_N_R1 * SURFACE_N_R2
    done = len(completed)
    print(f"  Resuming: {done}/{total} points already completed.\n")

    scan_base = os.path.join(output_dir, "surface_scans")
    # Use the isotropic beta as initial guess; update adaptively
    last_beta_c = BETA_INIT

    for i, r1 in enumerate(r1_vals):
        for j, r2 in enumerate(r2_vals):
            key = (round(r1, 8), round(r2, 8))
            if key in completed:
                last_beta_c = completed[key]
                continue

            done += 1
            k3 = 1.0
            k1, k2 = r1 * k3, r2 * k3
            label = f"r1_{r1:.4f}_r2_{r2:.4f}"
            scan_dir = os.path.join(scan_base, label)

            margin = max(0.20 * last_beta_c, 0.04)
            beta_lo = max(0.01, last_beta_c - margin)
            beta_hi = last_beta_c + margin

            print(f"  [{done}/{total}] r1={r1:.4f} r2={r2:.4f}: scanning ...",
                  end="", flush=True)
            t0 = time.time()
            try:
                beta_c, chi_peak, _, _, _ = mc.find_beta_c(
                    EXE, TEST_Lx, TEST_Ly, 0, 0, k1, k2, k3,
                    beta_lo, beta_hi,
                    n_traj_coarse=SURFACE_N_TRAJ_COARSE,
                    n_traj_fine=SURFACE_N_TRAJ_FINE,
                    data_dir=scan_dir,
                )
            except Exception as e:
                print(f"  FAILED: {e}")
                continue
            elapsed = time.time() - t0
            print(f" β_c={beta_c:.8f} ({elapsed:.0f}s)")

            completed[key] = beta_c
            last_beta_c = beta_c

            # Save incrementally
            _save_surface(surface_path, r1_vals, r2_vals, completed)

    _save_surface(surface_path, r1_vals, r2_vals, completed)
    print(f"\n  Surface saved: {surface_path}  ({len(completed)} points)")
    return surface_path


def _save_surface(path, r1_vals, r2_vals, completed):
    """Write surface JSON with current data."""
    points = []
    for (r1, r2), bc in sorted(completed.items()):
        points.append({"r1": r1, "r2": r2, "beta_c": bc})
    data = {
        "Lx": TEST_Lx, "Ly": TEST_Ly, "Tx": 0, "Ty": 0,
        "r1_vals": [round(float(v), 8) for v in r1_vals],
        "r2_vals": [round(float(v), 8) for v in r2_vals],
        "n_r1": len(r1_vals), "n_r2": len(r2_vals),
        "n_traj_scan_coarse": SURFACE_N_TRAJ_COARSE,
        "n_traj_scan_fine": SURFACE_N_TRAJ_FINE,
        "points": points,
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


# ===================================================================
# Phase 2: Grid search optimisation using the surface
# ===================================================================

def phase2_optimize(surface_path, ref_meta, output_dir):
    """Adaptive grid search using the β_c surface as LUT."""
    print(f"\n{'='*70}")
    print("Phase 2: Grid search optimisation (surface LUT)")
    print(f"  Test lattice: {TEST_Lx}×{TEST_Ly}  Tx=0 Ty=0")
    print(f"  Reference: {REF_Lx}×{REF_Ly}  Tx={REF_Tx} Ty={REF_Ty}")
    print(f"  Grid: {OPT_N_GRID}×{OPT_N_GRID}  Levels: {OPT_MAX_LEVELS}")
    print(f"  Start: r1={OPT_R1_INIT}, r2={OPT_R2_INIT}, "
          f"half_span={OPT_HALF_SPAN}")
    print(f"  MC production: {OPT_N_TRAJ_PROD} trajectories per point")
    print(f"  Top-{OPT_TOP_N} Z⁴ re-scoring per level")
    print(f"{'='*70}\n")

    # Load surface
    surface = BetacSurface(surface_path)
    surface.summary()

    # Load reference data
    ref_data = mc.load_all_to_all(ref_meta["a2a_path"])
    print(f"\n  Reference loaded: β_c={ref_meta['beta_c']:.10f}")

    opt_dir = os.path.join(output_dir, "optimization")
    run_dir = os.path.join(opt_dir, "runs")
    os.makedirs(opt_dir, exist_ok=True)

    r1_center = OPT_R1_INIT
    r2_center = OPT_R2_INIT
    hs_r1 = OPT_HALF_SPAN
    hs_r2 = OPT_HALF_SPAN

    all_results = []

    for level in range(OPT_MAX_LEVELS):
        r1_lo = max(r1_center - hs_r1, SURFACE_R1_LO)
        r1_hi = min(r1_center + hs_r1, SURFACE_R1_HI)
        r2_lo = max(r2_center - hs_r2, SURFACE_R2_LO)
        r2_hi = min(r2_center + hs_r2, SURFACE_R2_HI)

        r1_vals = np.linspace(r1_lo, r1_hi, OPT_N_GRID)
        r2_vals = np.linspace(r2_lo, r2_hi, OPT_N_GRID)

        print(f"\n{'='*70}")
        print(f"Level {level}: centre=({r1_center:.4f}, {r2_center:.4f})"
              f"  hs=({hs_r1:.4f}, {hs_r2:.4f})")
        print(f"  r1 ∈ [{r1_lo:.4f}, {r1_hi:.4f}]  "
              f"r2 ∈ [{r2_lo:.4f}, {r2_hi:.4f}]")
        print(f"{'='*70}")
        sys.stdout.flush()

        # --- Evaluate all grid points ---
        grid = {}
        for i, r1 in enumerate(r1_vals):
            for j, r2 in enumerate(r2_vals):
                label = f"L{level}_i{i}_j{j}"
                k3 = 1.0
                k1, k2 = r1 * k3, r2 * k3

                # LUT lookup
                beta_c = surface(r1, r2)
                inside = surface.domain_contains(r1, r2)

                print(f"  [{label}] r1={r1:.6f} r2={r2:.6f}: "
                      f"β_c={beta_c:.8f} (LUT{'!' if not inside else ''})",
                      end="", flush=True)

                # Production run
                t0 = time.time()
                level_dir = os.path.join(run_dir, f"level{level:02d}")
                prod_dir = os.path.join(level_dir, f"prod_{label}")
                try:
                    stdout, subdir = mc.run_simulator(
                        EXE, TEST_Lx, TEST_Ly, 0, 0, k1, k2, k3, beta_c,
                        n_traj=OPT_N_TRAJ_PROD, n_therm=3000,
                        data_dir=prod_dir)
                except Exception as e:
                    print(f"  FAILED: {e}")
                    continue
                a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
                test_data = mc.load_all_to_all(a2a_path)

                # Z² cost
                chi2, ndof, chi2_err = mc.boundary_slices_normed(
                    ref_data, test_data,
                    TEST_Lx, TEST_Ly, 0, 0,
                    REF_Lx, REF_Ly, REF_Tx, REF_Ty)
                z2_cost = chi2 / max(ndof, 1)
                t_total = time.time() - t0
                print(f"  Z²={z2_cost:.4f} ({t_total:.0f}s)")
                sys.stdout.flush()

                grid[(i, j)] = {
                    "r1": r1, "r2": r2, "beta_c": beta_c,
                    "z2": z2_cost, "a2a_path": a2a_path,
                }
                all_results.append({
                    "level": level, "i": i, "j": j,
                    "r1": float(r1), "r2": float(r2),
                    "beta_c": beta_c, "z2_cost": z2_cost,
                })

        if not grid:
            print("  WARNING: no points evaluated at this level!")
            continue

        # --- Re-score top N with Z⁴ ---
        ranked = sorted(grid.keys(), key=lambda k: grid[k]["z2"])
        top_keys = ranked[:OPT_TOP_N]

        print(f"\n  Re-scoring top {OPT_TOP_N} with Z⁴:")
        for rank, key in enumerate(top_keys):
            pt = grid[key]
            test_data = mc.load_all_to_all(pt["a2a_path"])
            chi2_q, ndof_q, _ = mc.boundary_slices_quartic(
                ref_data, test_data,
                TEST_Lx, TEST_Ly, 0, 0,
                REF_Lx, REF_Ly, REF_Tx, REF_Ty)
            z4 = chi2_q / max(ndof_q, 1)
            pt["z4"] = z4
            print(f"    #{rank+1}: ({key[0]},{key[1]}) "
                  f"r1={pt['r1']:.4f} r2={pt['r2']:.4f}  "
                  f"Z²={pt['z2']:.4f}  Z⁴={z4:.4f}")

        # Best by Z⁴ among top candidates
        best_key = min(top_keys,
                       key=lambda k: grid[k].get("z4", grid[k]["z2"]))
        bi, bj = best_key
        best = grid[best_key]
        r1_best, r2_best = best["r1"], best["r2"]
        best_z2 = best["z2"]
        best_z4 = best.get("z4", best_z2)
        best_beta_c = best["beta_c"]

        print(f"\n  Level {level} best: ({bi},{bj})  "
              f"r1={r1_best:.6f} r2={r2_best:.6f}  "
              f"Z²={best_z2:.4f}  Z⁴={best_z4:.4f}  "
              f"β_c={best_beta_c:.8f}")

        # --- Border detect (translate-or-refine) ---
        r1_border = (bi == 0) or (bi == OPT_N_GRID - 1)
        r2_border = (bj == 0) or (bj == OPT_N_GRID - 1)

        # --- Save level summary ---
        level_summary = {
            "level": level,
            "r1_center": float(r1_center), "r2_center": float(r2_center),
            "hs_r1": float(hs_r1), "hs_r2": float(hs_r2),
            "r1_vals": [float(v) for v in r1_vals],
            "r2_vals": [float(v) for v in r2_vals],
            "best_ij": [bi, bj],
            "best_r1": float(r1_best), "best_r2": float(r2_best),
            "best_z2": float(best_z2), "best_z4": float(best_z4),
            "best_beta_c": float(best_beta_c),
            "r1_border": r1_border, "r2_border": r2_border,
            "grid": {f"{k[0]}_{k[1]}": {
                "r1": v["r1"], "r2": v["r2"],
                "beta_c": v["beta_c"], "z2": v["z2"],
            } for k, v in grid.items()},
        }
        with open(os.path.join(opt_dir, f"grid_level_{level:02d}.json"), "w") as f:
            json.dump(level_summary, f, indent=2)

        # --- Refine / translate ---
        r1_center = float(r1_best)
        if r1_border:
            print(f"  r1: border → translate to {r1_center:.4f} (span unchanged)")
        else:
            hs_r1 /= 2.0
            print(f"  r1: interior → refine, hs_r1={hs_r1:.4f}")

        r2_center = float(r2_best)
        if r2_border:
            print(f"  r2: border → translate to {r2_center:.4f} (span unchanged)")
        else:
            hs_r2 /= 2.0
            print(f"  r2: interior → refine, hs_r2={hs_r2:.4f}")

        sys.stdout.flush()

    # --- Final results ---
    print(f"\n{'='*70}")
    print("OPTIMISATION COMPLETE")
    print(f"{'='*70}")

    if all_results:
        best_overall = min(all_results, key=lambda r: r["z2_cost"])
        print(f"Best point: r1={best_overall['r1']:.6f} "
              f"r2={best_overall['r2']:.6f}  "
              f"Z²={best_overall['z2_cost']:.6f}  "
              f"β_c={best_overall['beta_c']:.8f}")

        fit_result = {
            "best_r1": best_overall["r1"],
            "best_r2": best_overall["r2"],
            "best_k3": 1.0,
            "best_beta_c": best_overall["beta_c"],
            "best_z2_cost": best_overall["z2_cost"],
            "method": "critical_surface_lut + Z²/Z⁴ hybrid",
            "ref_Lx": REF_Lx, "ref_Ly": REF_Ly,
            "ref_Tx": REF_Tx, "ref_Ty": REF_Ty,
            "test_Lx": TEST_Lx, "test_Ly": TEST_Ly,
            "levels_completed": min(level + 1, OPT_MAX_LEVELS),
            "all_results": all_results,
        }
        with open(os.path.join(opt_dir, "fit_result.json"), "w") as f:
            json.dump(fit_result, f, indent=2)
        print(f"\nResults saved: {opt_dir}/fit_result.json")
    else:
        print("No results collected.")

    return all_results


# ===================================================================
# Main
# ===================================================================

def main():
    t_start = time.time()

    _build_exe()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save configuration snapshot
    config = {
        "ref": {"Lx": REF_Lx, "Ly": REF_Ly, "Tx": REF_Tx, "Ty": REF_Ty,
                "beta_c_override": REF_BETA_C, "n_traj": REF_N_TRAJ},
        "test": {"Lx": TEST_Lx, "Ly": TEST_Ly, "Tx": 0, "Ty": 0},
        "surface": {
            "r1_range": [SURFACE_R1_LO, SURFACE_R1_HI],
            "r2_range": [SURFACE_R2_LO, SURFACE_R2_HI],
            "n_r1": SURFACE_N_R1, "n_r2": SURFACE_N_R2,
        },
        "opt": {
            "r1_init": OPT_R1_INIT, "r2_init": OPT_R2_INIT,
            "half_span": OPT_HALF_SPAN,
            "n_grid": OPT_N_GRID, "max_levels": OPT_MAX_LEVELS,
            "n_traj_prod": OPT_N_TRAJ_PROD,
        },
    }
    with open(os.path.join(OUTPUT_DIR, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # Phase 0
    ref_meta = phase0_reference(OUTPUT_DIR)

    # Phase 1
    surface_path = phase1_surface(OUTPUT_DIR)

    # Phase 2
    results = phase2_optimize(surface_path, ref_meta, OUTPUT_DIR)

    t_total = time.time() - t_start
    print(f"\nTotal wall time: {t_total/3600:.1f} hours ({t_total:.0f}s)")
    print("Done.")


if __name__ == "__main__":
    main()
