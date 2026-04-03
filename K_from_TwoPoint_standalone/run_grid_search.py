#!/usr/bin/env python3
"""
run_grid_search.py
==================
User-facing entry point for the K_from_TwoPoint adaptive grid search.

Edit the CONFIG section below, then run:

    python run_grid_search.py

All results are written to the directory specified by OUTPUT_DIR.
A live plot (optimisation_live.png) is updated after each grid point.
"""

# ===========================================================================
# USER CONFIGURATION — edit the values in this section
# ===========================================================================

# ---------------------------------------------------------------------------
# Lattice geometry for the runs you want to optimise
# (the reference is always the equilateral 32x32 lattice included in ref_data/)
# ---------------------------------------------------------------------------
Lx = 32          # lattice x size
Ly = 32          # lattice y size
Tx = 0           # twist in x-direction (0 for untwisted)
Ty = 0           # twist in y-direction (e.g. Ty = Lx//4 for 45-deg tilt)

# ---------------------------------------------------------------------------
# Initial grid search parameters
# The optimizer searches over (r1, r2) = (k1/k3, k2/k3).
# Start centred on (1, 1) = isotropic couplings; widen HALF_SPAN if the
# answer is expected to be far from isotropic.
# ---------------------------------------------------------------------------
R1_INIT    = 1.0   # starting grid centre, r1 = k1/k3
R2_INIT    = 1.0   # starting grid centre, r2 = k2/k3
HALF_SPAN  = 0.15  # half-width of grid in r-space (e.g. 0.15 → grid ±15%)

# ---------------------------------------------------------------------------
# Initial guess for the inverse critical temperature beta_c.
# For a 32×32 isotropic equilateral lattice beta_c ≈ 0.267.
# Adjust if your geometry is expected to be significantly different.
# ---------------------------------------------------------------------------
BETA_INIT = 0.27

# ---------------------------------------------------------------------------
# Monte Carlo statistics
# More trajectories → more reliable fits, but slower runs.
# For a quick test, lower all three to ~5000 / 10000 / 10000.
# For production quality (like test_32x32_v7), use 50000 / 30000 / 60000.
# ---------------------------------------------------------------------------
N_TRAJ_PROD         = 50000   # trajectories for production (two-point function)
N_TRAJ_SCAN_COARSE  = 30000   # trajectories per beta point during coarse beta scan
N_TRAJ_SCAN_FINE    = 60000   # trajectories per beta point during fine beta refine

# ---------------------------------------------------------------------------
# Cost function — how the two-point functions are compared
#   "fem_integral"  recommended: exact integral of (G_test - G_ref)^2 on a
#                   bilinear mesh; preserves physical amplitudes.
#   "log_ratio"     amplitude-free GLS log-ratio chi2 (good general default).
#   "pair_ratio"    all pairs of log-ratios; amplitude cancels exactly.
#   "residuals"     plain L2, no amplitude correction.
#   "beta_deriv"    first-order correction for beta mismatch (advanced).
# ---------------------------------------------------------------------------
COST = "fem_integral"

# ---------------------------------------------------------------------------
# Range of displacements included in the chi2 sum (for non-fem_integral costs)
# R_MIN:     ignore displacements closer than this (lattice units)
# R_MAX_FRAC: ignore displacements farther than this fraction of L_eff = min(Lx,Ly)
# ---------------------------------------------------------------------------
R_MIN      = 0.0
R_MAX_FRAC = 0.33

# ---------------------------------------------------------------------------
# Maximum number of grid refinement / translation levels.
# Each level evaluates a full 5×5 grid, so MAX_ITER=4 → up to 100 simulator runs.
# ---------------------------------------------------------------------------
MAX_ITER = 4

# ---------------------------------------------------------------------------
# Output directory (relative to this script's location)
# A subdirectory will be created; all plots, JSON summaries, and MC data go here.
# ---------------------------------------------------------------------------
OUTPUT_DIR = "results/my_run"

# ---------------------------------------------------------------------------
# Reference data
# Leave as None to use the bundled 32×32 equilateral reference.
# If the reference does not exist yet it will be generated automatically
# before the grid search starts (takes ~10–30 min for REF_N_TRAJ=500000).
# Set to a path to an existing ref_metadata.json to skip generation.
# ---------------------------------------------------------------------------
REFERENCE = None   # None → auto-use/generate ref_data/ref_metadata.json

# Number of MC trajectories for the automatic reference generation run.
# 500000 gives good statistics; use fewer (e.g. 50000) for a quick test.
REF_N_TRAJ = 500000

# ---------------------------------------------------------------------------
# Path to the compiled simulator (relative to this script's location).
# On Windows, compile with MSYS2/MinGW (see README.md) then set this to
# "bin/ising_tri_twisted_parallelogram.exe"
# ---------------------------------------------------------------------------
import sys as _sys
EXE = ("bin/ising_tri_twisted_parallelogram.exe"
       if _sys.platform == "win32"
       else "bin/ising_tri_twisted_parallelogram")

# ===========================================================================
# END OF USER CONFIGURATION — no need to edit below this line
# ===========================================================================

import os
import sys

# Make sure we can find optimise_couplings.py even if the script is called from
# a different working directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Change to the package root so all relative paths work correctly.
os.chdir(_HERE)

# Validate the binary exists
if not os.path.isfile(EXE):
    print(f"ERROR: Simulator not found at '{EXE}'")
    print("Please build it first with:  make")
    sys.exit(1)

# Choose reference path
if REFERENCE is None:
    REFERENCE = os.path.join("ref_data", "ref_metadata.json")

# Auto-generate reference if it doesn't exist yet
if not os.path.isfile(REFERENCE):
    print(f"Reference not found at '{REFERENCE}' — generating it now.")
    print(f"  Lattice: {Lx}x{Ly}  Tx={Tx} Ty={Ty}")
    print(f"  Trajectories: {REF_N_TRAJ}  (edit REF_N_TRAJ in CONFIG to change)")
    print()
    import optimise_couplings as _oc
    ref_output_dir = os.path.dirname(REFERENCE)
    if not ref_output_dir:
        ref_output_dir = "."
    sys.argv = [
        "optimise_couplings.py",
        "--exe",        EXE,
        "--Lx",         str(Lx),
        "--Ly",         str(Ly),
        "--Tx",         str(Tx),
        "--Ty",         str(Ty),
        "--beta_init",  str(BETA_INIT),
        "--ref_n_traj", str(REF_N_TRAJ),
        "--output_dir", ref_output_dir,
        "--gen_ref",
    ]
    import importlib
    importlib.reload(_oc)
    _oc.main()
    print()
    if not os.path.isfile(REFERENCE):
        print(f"ERROR: Reference generation did not produce '{REFERENCE}'")
        sys.exit(1)
    print(f"Reference generated. Proceeding with grid search.\n")

# Build sys.argv so that optimise_couplings.main() picks up our settings
sys.argv = [
    "optimise_couplings.py",
    "--exe",               EXE,
    "--Lx",               str(Lx),
    "--Ly",               str(Ly),
    "--Tx",               str(Tx),
    "--Ty",               str(Ty),
    "--ref",              REFERENCE,
    "--r1_init",          str(R1_INIT),
    "--r2_init",          str(R2_INIT),
    "--half_span",        str(HALF_SPAN),
    "--beta_init",        str(BETA_INIT),
    "--n_traj_prod",      str(N_TRAJ_PROD),
    "--n_traj_scan_coarse", str(N_TRAJ_SCAN_COARSE),
    "--n_traj_scan_fine", str(N_TRAJ_SCAN_FINE),
    "--cost",             COST,
    "--r_min",            str(R_MIN),
    "--r_max_frac",       str(R_MAX_FRAC),
    "--max_iter",         str(MAX_ITER),
    "--output_dir",       OUTPUT_DIR,
]

import optimise_couplings
optimise_couplings.main()
