"""Shared pytest fixtures.

The tiny_geom fixture exposes a CONFIG-shaped dict tuned for fast smoke
tests (~5–15 s per evaluator call on a single core).  Tests that need
the real C++ simulator can use it; tests for pure-Python modules
(betac_cache, parallel scheduling) typically don't need to run any MC
at all and can monkey-patch the evaluator.
"""
from __future__ import annotations

import os
import sys

# Make the package modules importable when pytest is run from anywhere.
HERE = os.path.dirname(os.path.abspath(__file__))
PKG  = os.path.dirname(HERE)
if PKG not in sys.path:
    sys.path.insert(0, PKG)

import pytest


@pytest.fixture
def tiny_geom():
    """A small CONFIG dict; ~seconds per eval if the simulator is invoked."""
    return {
        "ref_Lx":        8, "ref_Ly":        8, "ref_Tx": 0, "ref_Ty": 0,
        "test_Lx":       8, "test_Ly":       8, "test_Tx": 0, "test_Ty": 0,
        "ref_n_traj":               2_000,
        "ref_scan_n_traj_coarse":     500,
        "ref_scan_n_traj_fine":     1_000,
        "ref_beta_seed":            (0.20, 0.40),
        "ref_scan_n_coarse":        7,
        "ref_scan_n_refine":        3,
        "ref_scan_n_refine2":       3,
        "ref_scan_n_refine3":       3,
        "ref_scan_max_shifts":      2,
        "n_traj_prod":              2_000,
        "n_traj_scan_coarse":         500,
        "n_traj_scan_fine":         1_000,
        "beta_seed":                (0.20, 0.40),
        "scan_n_coarse":            7,
        "scan_n_refine":            3,
        "scan_n_refine2":           3,
        "scan_n_refine3":           3,
        "scan_max_shifts":          2,
        "scan_jackknife":           False,
        "optimizer":        "nelder_mead",
        "x0":               (1.0, 1.0),
        "max_evals":        5,
        "nm_sigma0":        0.10,
        "nm_xatol":         0.005,
        "nm_fatol":         1e-6,
        "nm_shrink":        0.75,
        "cma_sigma0":       0.10,
        "cma_popsize":      4,
        "cma_max_gens":     0,
        "cma_tolx":         0.005,
        "cma_tolfun":       1e-6,
        "cma_seed":         12345,
        "snr_target":       0.0,
        "snr_max_traj_factor": 4,
        "snr_floor":        0.0,
        "indist_stop_snr":  0.0,
        # Speedup flags — defaults match production CONFIG (all OFF for
        # legacy parity).  Tests can override per-case.
        "betac_cache":      False,
        "betac_tol_r":      0.05,
        "betac_tol_beta":   2e-3,
        "betac_min_neighbours": 4,
        "n_workers":        1,
        "master_seed":      None,
        "backend":          "subprocess",
        "results_dir":      "results",
        "run_name":         None,
        "save_every":       1,
        "no_vis":           True,
        "dashboard":        False,
        "exe":              "bin/ising_tri_twisted_parallelogram",
    }
