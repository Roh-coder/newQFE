#!/usr/bin/env python3
"""
run_quick_test.py — Fast overnight test with reduced statistics.

Uses a small lattice (8×8) and coarse grids to verify the full pipeline
works end-to-end in ~30 minutes on a modern laptop.

Edit only if needed, then run:

    python run_quick_test.py
"""

# ===========================================================================
# Quick test configuration (8×8, low stats)
# ===========================================================================

# Reference: 8×8 equilateral, untwisted
REF_Lx = 8
REF_Ly = 8
REF_Tx = 0
REF_Ty = 0
REF_BETA_C = None     # let the scan find it
REF_N_TRAJ = 50000
REF_N_TRAJ_SCAN_COARSE = 5000
REF_N_TRAJ_SCAN_FINE   = 10000

# Test: 8×8 untwisted, anisotropic
TEST_Lx = 8
TEST_Ly = 8

# Surface: small region, coarse grid
SURFACE_R1_LO = 0.6
SURFACE_R1_HI = 1.6
SURFACE_R2_LO = 0.6
SURFACE_R2_HI = 1.6
SURFACE_N_R1  = 5     # 5×5 = 25 scan points
SURFACE_N_R2  = 5
SURFACE_N_TRAJ_COARSE = 5000
SURFACE_N_TRAJ_FINE   = 10000

# Optimization: fast grid search
OPT_R1_INIT   = 1.0
OPT_R2_INIT   = 1.0
OPT_HALF_SPAN = 0.30
OPT_N_GRID    = 3     # 3×3 = 9 points per level
OPT_MAX_LEVELS = 3
OPT_TOP_N     = 2
OPT_N_TRAJ_PROD = 10000

BETA_INIT = 0.27
OUTPUT_DIR = "results/quick_test"

# ===========================================================================
# Import and override run_pipeline's config, then run
# ===========================================================================

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import run_pipeline as rp

# Patch the config variables
for attr in [
    "REF_Lx", "REF_Ly", "REF_Tx", "REF_Ty", "REF_BETA_C", "REF_N_TRAJ",
    "REF_N_TRAJ_SCAN_COARSE", "REF_N_TRAJ_SCAN_FINE",
    "TEST_Lx", "TEST_Ly",
    "SURFACE_R1_LO", "SURFACE_R1_HI", "SURFACE_R2_LO", "SURFACE_R2_HI",
    "SURFACE_N_R1", "SURFACE_N_R2",
    "SURFACE_N_TRAJ_COARSE", "SURFACE_N_TRAJ_FINE",
    "OPT_R1_INIT", "OPT_R2_INIT", "OPT_HALF_SPAN",
    "OPT_N_GRID", "OPT_MAX_LEVELS", "OPT_TOP_N", "OPT_N_TRAJ_PROD",
    "BETA_INIT", "OUTPUT_DIR",
]:
    setattr(rp, attr, globals()[attr])

rp.main()
