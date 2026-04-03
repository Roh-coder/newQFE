#!/usr/bin/env python3
"""
gen_ref.py
==========
Generate a new equilateral reference dataset for a given lattice size.

This only needs to be done ONCE per lattice size.  The bundled reference
in ref_data/ was generated for 32×32.  Use this script if you want to
work with a different lattice size (e.g. 16×16 or 48×48).

Edit the CONFIG block below, then run:

    python gen_ref.py

The script runs a susceptibility-peak scan followed by a long production
run at beta_c and writes:
    <OUTPUT_DIR>/ref_production/...../two_point_all_to_all.dat
    <OUTPUT_DIR>/ref_metadata.json

Point subsequent calls to run_grid_search.py at the new ref_metadata.json:
    REFERENCE = "<OUTPUT_DIR>/ref_metadata.json"
"""

# ===========================================================================
# USER CONFIGURATION
# ===========================================================================

# Lattice size for the reference (should match Lx, Ly in run_grid_search.py)
REF_Lx = 32
REF_Ly = 32
REF_Tx = 0
REF_Ty = 0

# Initial guess for beta_c (equilateral isotropic triangular Ising):
#   Kc = ln(3)/4 ≈ 0.2747  →  beta_c ≈ 0.267 for 32×32
BETA_INIT = 0.27

# Number of trajectories for the reference production run.
# 500000 gives good statistics; use fewer for a quick test.
REF_N_TRAJ = 500000

# Output location for the reference data
OUTPUT_DIR = "ref_data/custom_ref"

EXE = "bin/ising_tri_twisted_parallelogram"

# ===========================================================================
# END OF USER CONFIGURATION
# ===========================================================================

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
os.chdir(_HERE)

if not os.path.isfile(EXE):
    print(f"ERROR: Simulator not found at '{EXE}'")
    print("Please build it first with:  make")
    sys.exit(1)

sys.argv = [
    "optimise_couplings.py",
    "--exe",        EXE,
    "--Lx",         str(REF_Lx),
    "--Ly",         str(REF_Ly),
    "--Tx",         str(REF_Tx),
    "--Ty",         str(REF_Ty),
    "--beta_init",  str(BETA_INIT),
    "--ref_n_traj", str(REF_N_TRAJ),
    "--output_dir", OUTPUT_DIR,
    "--gen_ref",
]

import optimise_couplings
optimise_couplings.main()
