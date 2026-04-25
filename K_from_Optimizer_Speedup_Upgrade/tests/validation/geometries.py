"""Geometry catalogue for the post-speedup validation panel.

Every test runs on the same fixed test (target) lattice:
    Lx = Ly = 12, Tx = Ty = 0   (square, untwisted, equilateral)

What varies is the reference (target-correlator) geometry.  These five
tags exercise: same lattice, scaled lattice, twisted, rectangular, and
pure-shear — covering the corners of the cost-surface stress space.

Each entry is a dict consumable by run.py CONFIG:
    ref_Lx, ref_Ly, ref_Tx, ref_Ty,
    test_Lx (=12), test_Ly (=12), test_Tx (=0), test_Ty (=0),
    expected_optimum: rough (r1*, r2*) ground-truth ball centre
    optimum_radius:   acceptance ball radius (passed to the panel cost check)
"""
from __future__ import annotations

REFERENCE_GEOMETRIES = {
    "R-equi": {
        "ref_Lx": 12, "ref_Ly": 12, "ref_Tx": 0, "ref_Ty": 0,
        "expected_optimum": (1.0, 1.0),
        "optimum_radius": 0.05,
        "description": "Self-consistency: same geometry as the test lattice; "
                        "optimum must be (1, 1) to within MC noise.",
    },
    "R-equiL": {
        "ref_Lx": 24, "ref_Ly": 24, "ref_Tx": 0, "ref_Ty": 0,
        "expected_optimum": (1.0, 1.0),
        "optimum_radius": 0.08,
        "description": "Cross-size, same shape: optimum still near (1, 1) "
                        "but MC noise scaled by the size mismatch.",
    },
    "R-twist": {
        "ref_Lx": 13, "ref_Ly": 16, "ref_Tx": 3, "ref_Ty": -3,
        "expected_optimum": None,   # unknown; recorded by Phase 0.1 baseline
        "optimum_radius": 0.10,
        "description": "Twisted parallelogram reference; pushes optimum off "
                        "(1, 1).  Acceptance ball centred on the legacy "
                        "baseline best.",
    },
    "R-rect": {
        "ref_Lx": 16, "ref_Ly": 12, "ref_Tx": 0, "ref_Ty": 0,
        "expected_optimum": None,
        "optimum_radius": 0.10,
        "description": "Rectangular zero-twist reference; r1 ≠ r2 expected.",
    },
    "R-shear": {
        "ref_Lx": 15, "ref_Ly": 15, "ref_Tx": 4, "ref_Ty": 0,
        "expected_optimum": None,
        "optimum_radius": 0.10,
        "description": "Pure shear, no rectangular stretch; off-axis ratio.",
    },
}

# Fixed test (target) lattice — same for all panel runs.
TEST_GEOMETRY = {
    "test_Lx": 12, "test_Ly": 12, "test_Tx": 0, "test_Ty": 0,
}

# CMA-ES starting distributions for the convergence sub-panel.
# Each entry is (start_id, x0, sigma0).  See README "Validation test
# panel" §C.
CMAES_START_DISTRIBUTIONS = [
    ("near",   (1.00, 1.00), 0.10),
    ("offset", (0.85, 1.15), 0.15),
    ("far_lo", (0.50, 0.50), 0.30),
    ("far_hi", (1.50, 1.50), 0.30),
    ("skew",   (1.50, 0.60), 0.30),
]

# Per-geometry × per-start panel cell: number of independent seeds.
PANEL_SEEDS_PER_CELL = 3

# CMA-ES eval budget for the panel runs (≈ 2× legacy budget).
PANEL_CMAES_MAX_EVALS = 50

# NM control test starts (cheap regression smoke test).
NM_CONTROL_STARTS = ["near", "far_lo"]
NM_CONTROL_MAX_EVALS = 30


def list_geometries() -> list:
    """Return list of geometry tags."""
    return sorted(REFERENCE_GEOMETRIES.keys())


def get_geometry(tag: str) -> dict:
    """Return the full geometry dict for ``tag`` (raises if unknown)."""
    if tag not in REFERENCE_GEOMETRIES:
        raise KeyError(f"unknown reference-geometry tag {tag!r}; "
                       f"choose from {list_geometries()}")
    g = dict(REFERENCE_GEOMETRIES[tag])
    g.update(TEST_GEOMETRY)
    g["tag"] = tag
    return g
