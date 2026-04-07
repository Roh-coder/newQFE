#!/usr/bin/env python3
"""
build_betac_surface.py — Pre-compute β_c(r₁, r₂) via MC susceptibility scans
=============================================================================

For a given test lattice (Lx, Ly, Tx, Ty), this script runs MC susceptibility
scans over a grid of (r₁, r₂) coupling-ratio values to find β_c at each
point.  The results are saved to a JSON file which can then be loaded by
the grid-search optimizer to look up β_c via FEM interpolation instead of
running a fresh susceptibility scan at every optimization grid point.

Usage:
    python build_betac_surface.py

Edit the CONFIG section below, then run.  The output is a JSON file containing:
  - The (r₁, r₂) grid coordinates
  - The measured β_c at each grid point
  - Lattice geometry metadata

The scan uses GC3 (3-pass Gram-Charlier) for each point, the same method
used in K_from_TwoPoint_standalone.  This is the expensive step done once;
the optimizer then benefits from instant, jitter-free β_c lookups.
"""

# ===========================================================================
# USER CONFIGURATION
# ===========================================================================

# --- Lattice geometry (must match the test lattice in the optimizer) ---
Lx = 32
Ly = 32
Tx = 0
Ty = 0

# --- Coupling-ratio grid for the β_c surface ---
# The grid should cover the region you expect to search.
# For an equilateral search near r₁ = r₂ = 1, a range of (0.5, 2.5)
# with ~11–21 points per axis gives a good surface.
R1_LO = 0.5
R1_HI = 2.5
R2_LO = 0.5
R2_HI = 2.5
N_R1  = 11       # grid points along r₁ axis
N_R2  = 11       # grid points along r₂ axis

# --- Initial beta guess (for the first scan point; subsequent points
#     use the β_c from the nearest completed neighbor) ---
BETA_INIT = 0.27

# --- MC scan statistics ---
N_TRAJ_SCAN_COARSE = 30000
N_TRAJ_SCAN_FINE   = 60000

# --- Output ---
OUTPUT_FILE = "betac_surfaces/surface_32x32_T0x0.json"

# --- Path to simulator ---
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
os.chdir(_HERE)

import optimise_couplings as oc
from scipy.optimize import curve_fit, minimize_scalar


# ---------------------------------------------------------------------------
# GC3 β_c finder (3-pass Gram-Charlier) — same as K_from_TwoPoint_standalone
# ---------------------------------------------------------------------------
def _gc_refit(sb, sc, beta_init):
    """Fit Gram-Charlier to scan data and return mode."""
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    sb_arr, sc_arr = np.array(sb), np.array(sc)
    b_rng = sb_arr.max() - sb_arr.min()
    p0 = [beta_init, b_rng / 2, 0.0, 0.0, float(sc_arr.max())]
    bounds = ([sb_arr.min(), 1e-8, -5.0, -3.0, 0.0],
              [sb_arr.max(), b_rng * 2, 5.0, 10.0, np.inf])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        popt, _ = curve_fit(_gram_charlier, sb_arr, sc_arr,
                            p0=p0, bounds=bounds, maxfev=8000)
    mu_fit, sigma_fit = popt[0], popt[1]
    mode_lo = max(sb_arr.min(), mu_fit - 3.0 * sigma_fit)
    mode_hi = min(sb_arr.max(), mu_fit + 3.0 * sigma_fit)
    res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                          bounds=(mode_lo, mode_hi), method='bounded')
    return float(np.clip(res.x, sb_arr.min(), sb_arr.max())), tuple(popt)


def find_beta_gc3(k1, k2, k3, beta_guess, scan_dir):
    """3-pass Gram-Charlier β_c finder."""
    margin = max(0.20 * beta_guess, 0.04)
    beta_lo = max(0.01, beta_guess - margin)
    beta_hi = beta_guess + margin

    beta_c_2, chi_peak, sb, sc, se, fit_params_12 = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3,
        beta_lo, beta_hi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE, n_traj_fine=N_TRAJ_SCAN_FINE,
        data_dir=scan_dir,
    )
    # Pass 3: very tight bracket around pass-2 result
    b_range = max(sb) - min(sb)
    half3 = b_range / (2 * 11) / 4
    fine3_betas = np.linspace(beta_c_2 - half3, beta_c_2 + half3, 7)
    for b in fine3_betas:
        scan_dir3 = scan_dir + "_p3"
        stdout_f, _ = oc.run_simulator(
            EXE, Lx, Ly, Tx, Ty, k1, k2, k3, float(b),
            n_traj=N_TRAJ_SCAN_FINE, data_dir=scan_dir3)
        chi_f, chi_f_err = oc.parse_stdout_value_with_err(stdout_f, "m_susc:")
        sb.append(float(b))
        sc.append(chi_f)
        se.append(chi_f_err)

    try:
        beta_c_3, _ = _gc_refit(sb, sc, beta_c_2)
    except Exception:
        beta_c_3 = beta_c_2
    return beta_c_3, chi_peak


# ---------------------------------------------------------------------------
# Build the surface
# ---------------------------------------------------------------------------
def main():
    os.makedirs(os.path.dirname(OUTPUT_FILE) or ".", exist_ok=True)

    r1_vals = np.linspace(R1_LO, R1_HI, N_R1)
    r2_vals = np.linspace(R2_LO, R2_HI, N_R2)

    print(f"Building β_c surface for {Lx}×{Ly} Tx={Tx} Ty={Ty}")
    print(f"  r₁ ∈ [{R1_LO}, {R1_HI}], {N_R1} points")
    print(f"  r₂ ∈ [{R2_LO}, {R2_HI}], {N_R2} points")
    print(f"  Total: {N_R1 * N_R2} scan points")
    print(f"  MC: coarse={N_TRAJ_SCAN_COARSE}, fine={N_TRAJ_SCAN_FINE}")
    print(f"  Output: {OUTPUT_FILE}")
    print()

    # Try to load existing partial results for resume support
    results = {}
    if os.path.isfile(OUTPUT_FILE):
        with open(OUTPUT_FILE) as f:
            existing = json.load(f)
        for entry in existing.get("points", []):
            key = (round(entry["r1"], 10), round(entry["r2"], 10))
            results[key] = entry["beta_c"]
        print(f"  Resuming: {len(results)} points already computed\n")

    scan_base = os.path.join("scan_data", f"{Lx}x{Ly}_T{Tx}x{Ty}")
    total = N_R1 * N_R2
    done = 0

    for i, r1 in enumerate(r1_vals):
        for j, r2 in enumerate(r2_vals):
            done += 1
            key = (round(float(r1), 10), round(float(r2), 10))
            if key in results:
                print(f"  [{done}/{total}] r1={r1:.4f} r2={r2:.4f}: "
                      f"cached β_c={results[key]:.8f}")
                continue

            k3 = 1.0
            k1, k2 = r1 * k3, r2 * k3

            # Use nearest completed neighbor as beta guess
            beta_guess = _best_guess(results, r1, r2, BETA_INIT)

            scan_dir = os.path.join(scan_base, f"r1_{r1:.4f}_r2_{r2:.4f}")
            print(f"  [{done}/{total}] r1={r1:.4f} r2={r2:.4f}: "
                  f"scanning (guess={beta_guess:.6f}) ...", end="", flush=True)
            t0 = time.time()

            try:
                beta_c, _ = find_beta_gc3(k1, k2, k3, beta_guess, scan_dir)
            except Exception as e:
                print(f" FAILED: {e}")
                continue

            dt = time.time() - t0
            results[key] = float(beta_c)
            print(f" β_c={beta_c:.8f} ({dt:.0f}s)")

            # Save incrementally
            _save(r1_vals, r2_vals, results, OUTPUT_FILE)

    _save(r1_vals, r2_vals, results, OUTPUT_FILE)
    print(f"\nDone. {len(results)} points saved to {OUTPUT_FILE}")


def _best_guess(results, r1, r2, default):
    """Find β_c from the nearest completed point as initial guess."""
    if not results:
        return default
    best_dist = float("inf")
    best_beta = default
    for (r1k, r2k), bc in results.items():
        d = (r1 - r1k) ** 2 + (r2 - r2k) ** 2
        if d < best_dist:
            best_dist = d
            best_beta = bc
    return best_beta


def _save(r1_vals, r2_vals, results, path):
    """Save surface data as JSON."""
    points = []
    for (r1, r2), bc in sorted(results.items()):
        points.append({"r1": r1, "r2": r2, "beta_c": bc})

    data = {
        "Lx": Lx, "Ly": Ly, "Tx": Tx, "Ty": Ty,
        "r1_vals": [round(float(v), 10) for v in r1_vals],
        "r2_vals": [round(float(v), 10) for v in r2_vals],
        "n_r1": len(r1_vals), "n_r2": len(r2_vals),
        "n_traj_scan_coarse": N_TRAJ_SCAN_COARSE,
        "n_traj_scan_fine": N_TRAJ_SCAN_FINE,
        "points": points,
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    main()
