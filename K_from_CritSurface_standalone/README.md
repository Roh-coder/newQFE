# K_from_CritSurface_standalone — Standalone Critical-Surface Pipeline

Find the anisotropic couplings (k₁, k₂, k₃) on an **untwisted** triangular
lattice whose two-point function matches a **twisted equilateral** reference
lattice, using a pre-computed β_c surface to eliminate per-point
susceptibility scans.

Self-contained package: download this directory, build, and run overnight.

---

## The Problem

Given a triangular lattice with some twist (Lx, Ly, Tx, Ty) and equilateral
couplings (k₁ = k₂ = k₃), the Ising model at criticality produces a
specific two-point correlation function G_ref(m, n).

**Question**: Can we find anisotropic couplings (k₁, k₂, k₃) on an
*untwisted* rectangular lattice (Tx = Ty = 0) that reproduce the same
correlator along all three lattice boundary directions?

This is a coupling-matching or "universality-mapping" problem: different
microscopic couplings may yield the same macroscopic (CFT) correlator if
they flow to the same universality class under RG.

---

## Algorithm

The pipeline has three phases:

### Phase 0: Generate Reference Data

1. **Input**: Reference lattice parameters (Lx_ref, Ly_ref, Tx_ref, Ty_ref)
   with equilateral couplings k₁ = k₂ = k₃ = 1.
2. **Find β_c** for the reference lattice:
   - If the user provides `REF_BETA_C`, use it directly (skip scan).
   - Otherwise, run a 3-pass Gram-Charlier (GC3) susceptibility scan to
     locate the finite-size pseudo-critical temperature.
3. **Production run** at β_c: long MC simulation producing the all-to-all
   connected two-point function G_ref(m, n).
4. **Output**: `ref_metadata.json` + `two_point_all_to_all.dat`.

### Phase 1: Build the β_c Surface

For the **test lattice** (Lx_test × Ly_test, Tx = Ty = 0) with anisotropic
couplings parameterised by r₁ = k₁/k₃ and r₂ = k₂/k₃ (k₃ = 1):

1. **Create a regular grid** of (r₁, r₂) points over the user-specified
   scanning region.
2. **For each grid point**: run a GC3 susceptibility scan on the test
   lattice to find β_c(r₁, r₂).
3. **Build a Delaunay triangulation** of the (r₁, r₂) grid and construct a
   piecewise-linear (P1 FEM) interpolant over it.
4. **Output**: `surface_{Lx}x{Ly}_T0x0.json` — a smooth, deterministic
   β_c(r₁, r₂) lookup table.

This is the expensive phase (~2–4 hours for an 11×11 grid) but is done
**once per test geometry** and amortised across all optimisation runs.
Progress is saved incrementally — you can interrupt and resume.

```
r₂ ↑
   │  ●───●───●───●───●
   │  │ ╲ │ ╲ │ ╲ │ ╲ │
   │  ●───●───●───●───●      ● = MC β_c scan point
   │  │ ╲ │ ╲ │ ╲ │ ╲ │      ─ = Delaunay triangulation
   │  ●───●───●───●───●
   │  │ ╲ │ ╲ │ ╲ │ ╲ │
   │  ●───●───●───●───●
   └──────────────────────→ r₁
```

### Phase 2: Adaptive Grid Search

With the surface built, each grid point evaluation is fast (production MC
only, no scan):

1. **Initialise** a search grid centred on (r₁_init, r₂_init) with
   half-span hs.
2. **For each refinement level** (up to MAX_LEVELS):
   a. Create an N_GRID × N_GRID test grid.
   b. For each (r₁, r₂) point:
      - **Look up** β_c from the surface (instant, deterministic).
      - **Run production MC** at that β_c on the test lattice.
      - **Compute Z² cost**: variance-normalised boundary-slice integral
        comparing G_test(t) vs G_ref(t) along the three torus directions
        (v, u, w).
   c. **Re-score** the top N candidates with Z⁴ (quartic variant for
      steeper discrimination).
   d. **Select best** by Z⁴.
   e. **Refine or translate**:
      - Interior best → halve the half-span (zoom 2×), re-centre.
      - Border best → keep span, translate centre.
3. **Output**: `fit_result.json` with optimal (r₁*, r₂*, β_c*) and per-level
   grid summaries.

```
┌─────────────────────────────────────────────────────────────┐
│ Phase 0          Phase 1              Phase 2               │
│                                                             │
│ Reference    →   β_c surface    →   Adaptive grid search   │
│ MC @ β_c         MC scans on        LUT lookup + MC prod   │
│ (twisted eq.)    (r₁,r₂) grid       Z²/Z⁴ cost matching   │
│                  Delaunay FEM        → best (r₁*, r₂*)     │
└─────────────────────────────────────────────────────────────┘
```

---

## Cost Function: Boundary-Slice Matching

The cost measures how well the test lattice reproduces the reference
correlator along the three fundamental directions of the torus:

- **v** = (Lx, Ty) — the first periodicity vector
- **u** = (Tx, −Ly) — the second periodicity vector
- **w** = −(v + u) — the diagonal

For each direction, the connected correlator is interpolated along the
boundary path as a function of t ∈ [0, 1]:

$$\text{Z}^2_d = \int_0^1 \frac{(G_\text{test}(t) - G_\text{ref}(t))^2}{\sigma_\text{ref}^2(t) + \sigma_\text{test}^2(t)} \, dt$$

The total cost aggregates per-direction costs as a sum of squares:

$$\text{Cost} = \sum_{d \in \{v,u,w\}} (\text{Z}^2_d)^2$$

This penalises large deviations in any single direction more than moderate
misses spread across all three.  The Z⁴ variant uses $(G/\sigma)^4$ for even
steeper sensitivity.

---

## Directory Structure

```
K_from_CritSurface_standalone/
├── README.md                 ← this file
├── Makefile                  ← builds the C++ simulator
├── requirements.txt          ← Python dependencies
├── .gitignore
│
├── run_pipeline.py           ← MAIN ENTRY POINT (all 3 phases)
├── run_quick_test.py         ← fast 8×8 end-to-end test (~30 min)
├── mc_engine.py              ← MC engine (simulator, β_c finder, costs)
├── betac_surface.py          ← FEM interpolant class + CLI inspector
│
├── src/
│   └── ising_tri_twisted_parallelogram.cc   ← C++ Ising simulator
├── include/
│   ├── ising.h               ← Ising model (Wolff + Metropolis)
│   ├── lattice.h             ← Triangular lattice on twisted parallelogram
│   ├── rng.h                 ← MT19937 RNG
│   ├── statistics.h          ← Measurement accumulators
│   └── timer.h               ← Wall-time tracking
├── bin/                      ← compiled binary (auto-built)
├── ref_data/                 ← reference MC data (auto-generated)
├── betac_surfaces/           ← pre-computed surface files
└── results/                  ← optimisation output
```

---

## Quick Start

### 1. Prerequisites

- **C++ compiler**: g++ with C++14 support (Linux, macOS, WSL, or MSYS2)
- **Python 3.8+** with numpy, scipy, matplotlib:
  ```bash
  pip install numpy scipy matplotlib
  ```

### 2. Build

```bash
cd K_from_CritSurface_standalone
make
```

This compiles `bin/ising_tri_twisted_parallelogram`.  If `make` is not
available, the Python scripts will attempt to build automatically via g++.

### 3. Quick Test (~30 minutes)

Run the fast 8×8 end-to-end test to verify everything works:

```bash
python run_quick_test.py
```

This uses an 8×8 lattice with reduced statistics (5×5 surface grid, 3×3
optimisation grid, 3 levels).  The output goes to `results/quick_test/`.

### 4. Production Run (overnight)

Edit the CONFIG section in `run_pipeline.py` to set your geometry:

```python
# Reference: equilateral triangles with twist
REF_Lx = 32;  REF_Ly = 32;  REF_Tx = 8;  REF_Ty = 0

# Test: untwisted, anisotropic
TEST_Lx = 32;  TEST_Ly = 32

# Scanning region
SURFACE_R1_LO = 0.5;  SURFACE_R1_HI = 2.5
SURFACE_R2_LO = 0.5;  SURFACE_R2_HI = 2.5
SURFACE_N_R1  = 11;   SURFACE_N_R2  = 11

# Override reference beta_c if known exactly
REF_BETA_C = 0.27465307216702745  # equilateral exact (infinite volume)
```

Then run:

```bash
python run_pipeline.py
```

The script will:
1. Generate the reference (or load if already done)
2. Build the surface (or resume from partial progress)
3. Run the adaptive grid search
4. Save results to `results/run_001/`

---

## Configuration Guide

### Reference Lattice

| Parameter | Description |
|-----------|-------------|
| `REF_Lx, REF_Ly` | Lattice dimensions |
| `REF_Tx, REF_Ty` | Twist vector (0,0 for untwisted) |
| `REF_BETA_C` | Override β_c (None = auto-scan) |
| `REF_N_TRAJ` | Production trajectories (500k default) |

**Tip**: For exact equilateral on an infinite lattice, β_c K = ln(3)/4.
For finite lattices, set `REF_BETA_C = None` to let the scan find the
finite-size pseudo-critical temperature.

### Test Lattice & Surface

| Parameter | Description |
|-----------|-------------|
| `TEST_Lx, TEST_Ly` | Test lattice dimensions (Tx=Ty=0 always) |
| `SURFACE_R1_LO/HI` | r₁ range for the surface |
| `SURFACE_R2_LO/HI` | r₂ range for the surface |
| `SURFACE_N_R1/N_R2` | Grid density (11×11 recommended) |

### Optimisation

| Parameter | Description |
|-----------|-------------|
| `OPT_R1_INIT, OPT_R2_INIT` | Initial grid centre |
| `OPT_HALF_SPAN` | Half-width of search grid |
| `OPT_N_GRID` | Points per side per level (5 → 25 evals) |
| `OPT_MAX_LEVELS` | Max refinement levels |
| `OPT_N_TRAJ_PROD` | Production trajectories per grid point |

---

## Inspecting the Surface

After Phase 1 completes, inspect the surface:

```bash
# Summary + query a specific point
python betac_surface.py results/run_001/betac_surfaces/surface_32x32_T0x0.json \
    --r1 1.0 --r2 1.0

# Generate a heatmap plot
python betac_surface.py results/run_001/betac_surfaces/surface_32x32_T0x0.json \
    --plot
```

---

## Output Files

After a complete run, the output directory contains:

```
results/run_001/
├── config.json               ← snapshot of all configuration parameters
├── ref_data/
│   ├── ref_metadata.json     ← reference β_c, lattice params, data path
│   └── ref_production/...    ← reference MC data files
├── betac_surfaces/
│   └── surface_32x32_T0x0.json  ← pre-computed β_c surface
├── surface_scans/            ← raw MC data from surface construction
│   └── r1_X.XXXX_r2_Y.YYYY/ ← per-point scan directories
└── optimization/
    ├── fit_result.json       ← FINAL RESULT: best (r₁, r₂, β_c, Z² cost)
    ├── grid_level_00.json    ← per-level grid summaries with all point data
    ├── grid_level_01.json
    └── runs/                 ← raw MC data from production runs
```

---

## When to Rebuild the Surface

You need a new surface when:
- **Test lattice geometry changes** (different Lx, Ly)
- **You want a different (r₁, r₂) scanning region**
- **You want higher accuracy** (denser grid or more statistics)

You do **not** need a new surface when:
- Changing the reference lattice (changing twist, changing REF_BETA_C)
- Changing the cost function or optimisation parameters
- Changing production statistics (OPT_N_TRAJ_PROD)

---

## Comparison with Sister Packages

| Aspect | K_from_TwoPoint_standalone | K_from_CritSurface_standalone |
|--------|---------------------------|-------------------------------|
| β_c method | GC3 scan per grid point (~100s) | LUT interpolation (instant) |
| β_c noise | σ ≈ 0.0007 (dominant noise) | Zero (deterministic) |
| Cost variance | ~50% from β_c, ~50% from MC | 100% from MC (averages cleanly) |
| Wall time per point | ~120s (scan + production) | ~15s (production only) |
| Setup cost | None | One-time surface build (~2–4h) |
| Reference generation | Built-in | Built-in |
| Self-contained | Yes | Yes |
| Cross-geometry support | Yes | Yes (ref twisted, test untwisted) |

---

## Physics Background

The 2D triangular Ising model with couplings (k₁, k₂, k₃) has an exact
critical condition in the thermodynamic limit:

$$\sinh(2\beta k_1)\sinh(2\beta k_2) + \sinh(2\beta k_2)\sinh(2\beta k_3) + \sinh(2\beta k_3)\sinh(2\beta k_1) = 1$$

On a **finite lattice**, the susceptibility peak (pseudo-critical point) is
shifted by O(1/L) finite-size corrections.  The MC susceptibility scan
captures this finite-size effect, which is essential for fair correlator
comparisons between the reference and test lattices at the *same effective
criticality*.

The boundary-slice cost function tests whether the two lattices produce
the same correlator along the three fundamental directions of the torus.
By the universality hypothesis, if the anisotropic couplings on the
untwisted lattice correctly reproduce the twisted equilateral correlator
in all three directions, the two models are in the same universality class
with the same effective conformal structure.

---

## Troubleshooting

**"simulator not found"** — Run `make` or install g++ and let the script
auto-compile.

**Surface scan is slow** — For a quick test, reduce `SURFACE_N_R1/N_R2` to
5×5 and lower `SURFACE_N_TRAJ_COARSE/FINE` to 5000/10000.

**"outside the surface domain" warnings** — The optimisation grid moved
outside the surface region.  Either widen `SURFACE_R1_LO/HI` and rebuild,
or narrow `OPT_HALF_SPAN`.

**Interrupted run** — Just re-run `python run_pipeline.py`.  All three
phases check for existing results and resume from where they left off.

**NaN in cost function** — Usually means the interpolator couldn't evaluate
a point (too few tiling copies or the reference and test geometries are
incompatible).  Check that REF and TEST lattice sizes are compatible.
