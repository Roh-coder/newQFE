# K_from_TwoPoint — Standalone Package

Find the critical anisotropic couplings (k₁, k₂, k₃) of the 2D Ising model
on a triangular lattice by matching the two-point function to a
high-statistics reference via an adaptive grid search.

Supports both **same-geometry** matching (e.g. 32×32 equilateral ref vs
32×32 anisotropic test) and **cross-geometry** matching (e.g. 4-5-6 twisted
ref vs 8×8 untwisted test), where the optimizer finds couplings on the test
lattice that reproduce the boundary-direction correlators of the reference.

---

## What this package contains

```
K_from_TwoPoint_standalone/
├── README.md                      ← this file
├── Makefile                       ← builds the C++ simulator
├── requirements.txt               ← Python dependencies
│
├── src/
│   └── ising_tri_twisted_parallelogram.cc   ← C++ simulator (self-contained)
├── include/
│   ├── ising.h, lattice.h, rng.h, statistics.h, timer.h  ← C++ headers
│
├── ref_data/
│   ├── two_point_all_to_all.dat   ← 32×32 equilateral reference correlator
│   ├── ref_metadata.json          ← beta_c and path for the 32×32 reference
│   ├── ref_8x8/                   ← auto-generated 8×8 equilateral reference
│   └── ref_13x16_t3x-3/          ← bundled 4-5-6 triangle reference
│
├── optimise_couplings.py          ← optimisation engine (do not edit)
├── run_grid_search.py             ← USER ENTRY POINT — hybrid GC3 + Z²/Z⁴
├── run_grid_search_8x8.py         ← fast 8×8 test (same hybrid engine)
├── run_grid_search_hybrid.py      ← original hybrid prototype (superseded)
├── run_stability_test.py          ← stability test: N evals of the same point
├── run_stability_normed_compare.py ← side-by-side stability comparison of 3 costs
├── run_betac_diagnostic.py        ← beta_c finder method comparison (5 methods)
├── run_fss_betac.py               ← finite-size scaling β_c test (4 methods × 5 sizes)
├── gen_ref.py                     ← optional: generate your own reference data
├── plot_boundary_slices.py        ← visualise boundary-direction correlator slices
├── plot_betac_diagnostic.py       ← 6-panel β_c diagnostic plots
├── plot_score_comparison.py       ← 6-panel cost function comparison plots
└── results/                       ← output directory (created on first run)
```

---

## Quick start

### 1. Requirements

- **C++ compiler**: g++ supporting C++14 (any modern GCC or Clang)
- **Python**: 3.8+
- **Python packages**: numpy, scipy, matplotlib

Install Python dependencies:

```bash
pip install -r requirements.txt
```

> **Linux**: install g++ via `sudo apt install build-essential` if not present.  
> **macOS**: `make` works out of the box with Apple clang.  
> **Windows**: the bundled binary in `bin/` is Linux-only. See [Windows build instructions](#windows) below.

---

### 2. Build the simulator

```bash
make
```

This compiles `src/ising_tri_twisted_parallelogram.cc` and places the binary
at `bin/ising_tri_twisted_parallelogram`.

---

### 3. Configure and run

Open **`run_grid_search.py`** in any text editor.  There is a clearly marked
`USER CONFIGURATION` section at the top.  The script uses the **hybrid
strategy** — GC3 (3-pass Gram-Charlier) for the β_c finder, normalised Z²
for grid evaluation, and Z⁴ re-scoring of top candidates per level.  The
main parameters are:

| Parameter | What it does | Default |
|-----------|--------------|---------|
| `Lx`, `Ly` | Lattice size | 32, 32 |
| `Tx`, `Ty` | Twist vector components | 0, 0 |
| `R1_INIT`, `R2_INIT` | Starting grid centre (r₁ = k₁/k₃, r₂ = k₂/k₃) | 1.5, 1.5 |
| `HALF_SPAN` | Half-width of the search grid | 0.30 |
| `BETA_INIT` | Initial guess for critical β | 0.27 |
| `N_TRAJ_PROD` | MC trajectories for each production run | 50 000 |
| `N_TRAJ_SCAN_COARSE` | MC trajectories per β point (coarse scan) | 30 000 |
| `N_TRAJ_SCAN_FINE` | MC trajectories per β point (fine scan) | 60 000 |
| `TOP_N` | Top N candidates per level re-scored with Z⁴ | 3 |
| `MAX_ITER` | Maximum grid refinement levels | 6 |
| `REF_BETA_C` | Exact β_c for reference generation (skip scan if set) | `None` |
| `OUTPUT_DIR` | Where results are written | `"results/equilateral_from_1.5"` |

After editing, run:

```bash
python run_grid_search.py
```

A live plot (`results/my_run/optimisation_live.png`) is written after every
grid-point evaluation showing the search progress.  The final results are
saved in `results/my_run/fit_result.json`.

---

### 4. Quick test (fast, low statistics)

For a quick sanity check before committing to a full run, reduce the
trajectory counts in `run_grid_search.py`:

```python
N_TRAJ_PROD        = 10000
N_TRAJ_SCAN_COARSE = 10000
N_TRAJ_SCAN_FINE   = 20000
MAX_ITER           = 2
OUTPUT_DIR         = "results/quick_test"
```

Or use the 8×8 fast test which is ~16× faster per MC sweep:

```bash
python run_grid_search_8x8.py
```

Both use the same hybrid GC3 + Z²/Z⁴ engine.

---

### 5. Running in the background (Linux / macOS)

For long runs it is convenient to use `nohup` and redirect output:

```bash
nohup python run_grid_search.py > results/my_run.log 2>&1 &
```

Monitor progress:

```bash
tail -f results/my_run.log
```

---

## Generating a new reference dataset

The bundled reference (`ref_data/`) was generated for a **32×32 equilateral**
lattice.  If you want to work with a different lattice size, generate a new
reference first:

1. Edit the CONFIG block in **`gen_ref.py`** to match your target lattice size.
2. Run:

   ```bash
   python gen_ref.py
   ```

   This runs a susceptibility scan and a long production run (500 k trajectories
   by default) and writes `ref_data/custom_ref/ref_metadata.json`.

3. In `run_grid_search.py`, point `REFERENCE` at the new metadata file:

   ```python
   REFERENCE = "ref_data/custom_ref/ref_metadata.json"
   ```

### Using an exact β_c

If you know the exact critical β for your reference geometry, you can skip
the susceptibility scan entirely by setting `REF_BETA_C` in the CONFIG
section of `run_grid_search.py` (or `run_grid_search_8x8.py`):

```python
REF_BETA_C = 0.27465307216702745   # ln(3)/4, exact equilateral triangular Ising
```

This is passed as `--ref_beta_c` to the optimizer's `--gen_ref` mode.
The production run is executed at exactly that β, avoiding any finite-size
shift from the susceptibility peak.  Set to `None` (default) to scan as
before.

---

## Understanding the output

After a completed run you will find in `OUTPUT_DIR/`:

| File | Description |
|------|-------------|
| `fit_result.json` | Best-fit r₁, r₂, β_c and chi²/ndof |
| `grid_level_NN.json` | Per-level summary of grid search |
| `optimisation_live.png` | Live plot updated after each point |
| `optimisation_final.png` | Final snapshot of the live plot |
| `runs/levelNN/` | Raw MC data for each grid point |

The key result is in `fit_result.json`:

```json
{
  "best_r1": 1.075,       ← k1/k3 at the best-fit point
  "best_r2": 0.924,       ← k2/k3 at the best-fit point
  "best_k3": 1.0,
  "best_beta_c": 0.26711, ← inverse critical temperature
  "best_chi2_ndof": 0.06  ← goodness of fit (< 1 is good)
}
```

The critical couplings are then:  **k₁ = r₁ · k₃**, **k₂ = r₂ · k₃**, **k₃** chosen to satisfy
the criticality condition exp(−2β·k₁) + exp(−2β·k₂) + exp(−2β·k₃) = 1.

---

## How it works (brief overview)

1. **Reference**: a long, high-statistics MC run on the equilateral 32×32 lattice
   at its critical β produces the reference two-point function G_ref(m, n).

2. **Outer loop**: the optimizer searches over the coupling ratio plane
   (r₁ = k₁/k₃, r₂ = k₂/k₃) on a 5×5 grid.  **Each axis is handled
   independently**: after each level the best-point index is checked
   per-dimension:
   - If the minimum is interior in that axis → halve that axis's half-span
     (zoom in).
   - If the minimum is on the border in that axis → translate the centre
     in that direction by an extra half-span so the next grid extends well
     beyond the old boundary, keeping the span unchanged.
   This means a side-border minimum refines the interior axis while
   translating only the border axis, and a corner minimum translates both.

3. **Inner loop (GC3 β_c finder)**: for each grid point (r₁, r₂), the
   critical β is found via **3-pass Gram-Charlier**:
   - Pass 1–2: coarse scan (11 points) + two refinement passes (5 + 5
     fine-statistics points) using all accumulated scan data.
   - Pass 3: 7 ultra-tight points around the pass-2 result, then a
     Gram-Charlier refit over all data.  This extra pass tightens β_c
     by ~35% vs the 2-pass finder (see β_c diagnostic results below).
   A production run at the resulting β_c measures G_test(m, n).

4. **Cost (hybrid Z²/Z⁴)**: each grid point is scored with **normalised Z²**
   (∫Z² dt, lowest coefficient of variation for stable ranking).  At the
   end of each level, the **top 3 candidates** are re-scored with
   **normalised Z⁴** (∫Z⁴ dt, steeper landscape) using the same MC data
   (no additional runs).  The Z⁴ winner is used 3for grid refinement.
   See [Cost functions](#cost-functions) below.

---

## Cost functions

| Name | Description |
|------|-------------|
| `"boundary_slices"` | Sum of ∫₀¹ (G_test(t) − G_ref(t))⁴ dt along the three torus boundary paths v, u, w. The **quartic** (diff⁴) integrand massively punishes large pointwise deviations while preserving physical amplitudes. Per-direction costs are printed for diagnostics. Supports cross-geometry mode. |
| `"boundary_slices_normed"` | **Variance-normalised quadratic**: ∫₀¹ Z(t)² dt, where Z = ΔG/σ and σ² = σ_ref² + σ_test². Under the null (test = ref) the integrand has expectation 1 per sample point, making the total cost independent of MC statistics. Very stable (low CoV) but less sensitive to large deviations than the quartic. |
| `"boundary_slices_normed_quartic"` | **Variance-normalised quartic**: ∫₀¹ Z(t)⁴ dt with the same Z = ΔG/σ normalisation. Combines the stability of normalisation (no CoV blow-up) with the steep landscape of the quartic: for signal s = δ/σ, the excess over the null scales as 6s² + s⁴, giving ~44% better discrimination than the quadratic at s = 2 and ~116% at s = 3. **Recommended** for production runs. |
| `"fem_integral"` | Exact integral of (G_test − G_ref)² over the full 2D periodic domain via bilinear FEM. Requires same-geometry ref and test. |
| `"log_ratio"` | Amplitude-free GLS log-ratio χ². Good general default for same-geometry. |
| `"pair_ratio"` | All pairs of log-ratios; amplitude cancels exactly. |
| `"residuals"` | Plain L² of log-diffs, no amplitude correction. |
| `"beta_deriv"` | First-order β-mismatch correction (advanced). |

The `boundary_slices` cost works by:
1. Building a tiled interpolant of G_conn(x, y) for each dataset (ref and test)
   from their respective lattice geometries.
2. Sampling along each of the three boundary-direction paths as a function of
   the fractional parameter t ∈ [0, 1].
3. Computing ∫(G_test − G_ref)⁴ dt for each path (quartic penalty) and summing.
   The per-direction breakdown (v, u, w) is printed to the log so that
   pathological directions can be identified at a glance.

When the reference and test have **different geometries** (cross-geometry mode),
each interpolator is built from its own torus periodicity, and the three
directions are matched by index (v↔v, u↔u, w↔w).

---

## Cross-geometry mode

You can match correlators across **different lattice geometries**.  For example,
match a 4-5-6 triangular (13×16, Tx=3, Ty=−3) isotropic reference against an
8×8 untwisted test lattice with anisotropic couplings:

```bash
python optimise_couplings.py \
  --Lx 8 --Ly 8 --Tx 0 --Ty 0 \
  --ref ref_data/ref_13x16_t3x-3/ref_metadata.json \
  --ref_Lx 13 --ref_Ly 16 --ref_Tx 3 --ref_Ty -3 \
  --cost boundary_slices \
  --beta_init 0.275 \
  --n_traj_prod 20000 --n_grid 3 --max_iter 2 \
  --output_dir results/cross_geom_456_vs_8x8
```

The `--ref_Lx`, `--ref_Ly`, `--ref_Tx`, `--ref_Ty` flags tell the optimizer
that the reference lives on a different torus.  The boundary-slice interpolants
are built independently for each geometry.

---

## Known lattice geometries

| Triangle | Lx | Ly | Tx | Ty | Volume | Notes |
|----------|----|----|----|----|--------|-------|
| Equilateral (small) | 8 | 8 | 0 | 0 | 64 | Fast tests; 8×8 ref auto-generated |
| Equilateral (standard) | 32 | 32 | 0 | 0 | 1024 | Bundled reference in `ref_data/` |
| 4-5-6 (base) | 13 | 16 | 3 | −3 | 208 | Bundled reference in `ref_data/ref_13x16_t3x-3/` |
| 4-5-6 (3× scaled) | 39 | 48 | 9 | −9 | 1872 | Used in 3D Nt scans |

---

## Visualising boundary slices

Use `plot_boundary_slices.py` to inspect the three boundary-direction
correlator slices for a single dataset or compare reference vs test:

```bash
# Single dataset: view the three boundary slices
python plot_boundary_slices.py \
  --dat results/test_456_iso/.../two_point_all_to_all.dat \
  --Lx 13 --Ly 16 --Tx 3 --Ty -3 \
  --output results/slices_456_iso.png

# Comparison: reference (dashed) vs test (solid) with residuals
python plot_boundary_slices.py \
  --ref  results/test_456_iso/.../two_point_all_to_all.dat \
  --test results/test_456_aniso/.../two_point_all_to_all.dat \
  --Lx 13 --Ly 16 --Tx 3 --Ty -3 \
  --output results/slices_456_compare.png
```

---

## Troubleshooting

### Windows

The binary in `bin/` is a Linux ELF executable and **will not run on Windows**.
You must compile it first using one of these options:

**Option A — MSYS2/MinGW (recommended, native Windows binary)**

1. Install MSYS2 from https://www.msys2.org/
2. Open the **MSYS2 MINGW64** terminal and run:
   ```bash
   pacman -S mingw-w64-x86_64-gcc make
   cd /c/path/to/K_from_TwoPoint_standalone
   make
   ```
3. The binary will be at `bin/ising_tri_twisted_parallelogram.exe`
4. In `run_grid_search.py`, update the EXE line to:
   ```python
   EXE = "bin/ising_tri_twisted_parallelogram.exe"
   ```

**Option B — WSL (Windows Subsystem for Linux)**

1. Enable WSL and install Ubuntu from the Microsoft Store
2. Inside WSL, `cd` to the package directory and run `make` as normal
3. Run `python run_grid_search.py` from the WSL terminal

---

**`make` fails with "command not found"**  
Install a C++ compiler: `sudo apt install build-essential` (Linux) or install
Xcode Command Line Tools on macOS (`xcode-select --install`).

**`Error: simulator failed`**  
The binary panicked or exited non-zero.  Check that the geometry parameters
(Lx, Ly, Tx, Ty) are valid.  The constraint is: the parallelogram spanned by
(Lx, 0) and (Tx, Ly) must tile the plane (Lx and Ly must be positive integers;
Tx and Ty can be any integers).

**`No artists with labels found`** (matplotlib warning)  
Harmless; appears before the first grid point is evaluated.

**Susceptibility scan shifts window repeatedly**  
Your BETA_INIT is too far from the true beta_c.  Look at the log for a
`[scan] peak at lower/upper edge` message and adjust BETA_INIT.

**Run is very slow**  
Reduce `N_TRAJ_*` for a first pass, then increase for production quality.
Each MC trajectory takes O(V) spin flips where V = Lx × Ly × Nt.

---

## Additional scripts

### `run_grid_search_8x8.py`

Fast 8×8 equilateral convergence test.  Uses the same **hybrid GC3 + Z²/Z⁴
engine** as `run_grid_search.py` but configured for an 8×8 lattice (starts at
r₁ = r₂ = 1.5, half-span 0.30, 6 levels).  The 8×8 reference is auto-generated
on first run and cached in `ref_data/ref_8x8/`.  Output goes to
`results/equilateral_8x8_from_1.5/`.

Current configuration:

| Parameter | Value |
|-----------|-------|
| `Lx, Ly` | 8, 8 |
| `R1_INIT, R2_INIT` | 1.5, 1.5 |
| `HALF_SPAN` | 0.30 |
| `N_TRAJ_PROD` | 50 000 |
| `N_TRAJ_SCAN_COARSE` | 50 000 |
| `N_TRAJ_SCAN_FINE` | 100 000 |
| `MAX_ITER` | 6 |
| `TOP_N` | 3 |

### `run_stability_test.py`

Measures the **statistical stability** of the optimizer score by evaluating the
same grid point 100 times under independent MC noise.  Each trial runs the
full pipeline (β_c scan → production → boundary_slices cost).  Outputs:

- `results/stability_test_8x8/stability_results.json` — per-trial β_c,
  chi²/ndof, and timing (saved incrementally after each trial).
- `results/stability_test_8x8/stability_histograms.png` — histograms of β_c,
  cost, and log₁₀(cost) distributions with summary statistics (mean, std,
  min, max, coefficient of variation).

Current configuration: 100 trials at (r₁, r₂) = (1, 1) on the 8×8 lattice
with 50 k / 50 k / 100 k trajectories (prod / coarse / fine).

### `run_stability_normed_compare.py`

Side-by-side paired comparison of all three boundary-slice cost functions:

1. `boundary_slices` — raw quartic ∫(diff⁴) dt
2. `boundary_slices_normed` — normalised quadratic ∫(Z²) dt
3. `boundary_slices_normed_quartic` — normalised quartic ∫(Z⁴) dt

Each trial generates ONE MC dataset (same β_c scan + production run) and
evaluates all three costs on it, so the comparison is perfectly paired.
Produces histograms of each cost's distribution and paired scatter plots
showing inter-cost correlation.  Output goes to
`results/stability_normed_compare/`.

Current configuration: 100 trials at (r₁, r₂) = (1, 1) on the 8×8 lattice
with 50 k / 50 k / 100 k trajectories.

### `run_betac_diagnostic.py`

Diagnoses how much of the cost-function noise comes from the β_c **finder**
vs pure MC **production** noise.  Tests 5 β_c methods side-by-side, all on
the same (r₁=r₂=1, 8×8) point:

| Method | Description | Extra cost vs gc2 |
|--------|-------------|-------------------|
| `exact` | Use reference's own β_c (bypass finder entirely) — **noise floor control** | Cheapest (no scan) |
| `gc2` | Current 2-pass Gram-Charlier (status quo) | 1× (baseline) |
| `gc3` | 3-pass Gram-Charlier (extra tighter refinement layer) | ~1.3× |
| `gc2_histat` | 2-pass GC with 2× trajectories per scan point | ~2× |
| `gc2_boot` | 2-pass GC + bootstrap median (10 resamples of scan data) | ~1× (reuses same scan) |

Each trial evaluates all three cost functions (raw quartic, normed Z², normed Z⁴)
on the SAME production data, providing a perfectly controlled comparison.
Output goes to `results/betac_diagnostic/`.

Current configuration: 30 trials × 5 methods, 8×8 lattice,
50k prod / 50k coarse scan / 100k fine scan.

---

## Current tests and results

### Tests completed / in progress

| Test | Script | Status |
|------|--------|--------|
| β_c finder diagnostic (5 methods, 30 trials) | `run_betac_diagnostic.py` | **Complete** (150/150) |
| 8×8 equilateral convergence (6 levels) | `run_grid_search_8x8.py` | **Complete** — r₁=0.96, r₂=1.00 |
| FSS β_c test (4 methods × 5 sizes × 10 trials) | `run_fss_betac.py` | **152/200** (L=4,6,8 done; L=12 8/10; L=16 pending) |
| Stability (raw quartic, 100 evals) | `run_stability_test.py` | 99/100 |
| Normed cost comparison (3 costs, 100 evals) | `run_stability_normed_compare.py` | 75/100 |

### Stability test results (raw quartic, 83/100 trials)

The original `boundary_slices` (raw quartic diff⁴) stability test at the
equilateral point (r₁ = r₂ = 1, 8×8 lattice, 50k/50k/100k trajectories):

| Metric | Value |
|--------|-------|
| β_c mean ± std | 0.24082 ± 0.00070 |
| chi²/ndof mean | 5.2 × 10⁻⁸ |
| chi²/ndof std | 8.2 × 10⁻⁸ |
| **CoV (cost)** | **1.574** |
| log₁₀ range | −10.7 to −6.3 (4.4 decades) |

A CoV of 1.6 means two evaluations of the *same physical point* can differ
by ~1000×.  This motivated the normalised cost functions.

### Normed cost comparison (20/100 trials, paired)

All three costs evaluated on the *same* MC data per trial:

| Cost function | Mean | Std | **CoV** | log₁₀ span |
|---------------|------|-----|--------|------------|
| Raw quartic ∫(diff⁴)dt | 7.3 × 10⁻⁸ | 1.4 × 10⁻⁷ | **1.92** | 3.5 dec |
| Normed Z² ∫(Z²)dt | 4.5 | 4.5 | **1.01** | 2.0 dec |
| Normed Z⁴ ∫(Z⁴)dt | 45 | 85 | **1.87** | 3.8 dec |

**Key findings** (20 paired trials):
- Normed Z² has the lowest CoV ≈ 1.0, roughly 2× better than both raw
  quartic and normed Z⁴.
- **log₁₀(Z⁴) vs log₁₀(Z²) correlation = 0.999** — they rank-order grid
  points identically.  The quartic power adds no discriminative information
  at the null; it only amplifies β_c-driven noise.
- The Z⁴/Z² ratio has CoV = 0.91, confirming that the quartic exponent
  does not provide stable independent information.
- CoV > 1 persists even for normed Z², confirming that the dominant noise
  source is **β_c finder jitter** (std ≈ 0.0007), not the cost function.

### β_c finder diagnostic (round 1 complete, running 30 rounds)

First trial of each method (preliminary — need 30 for reliable statistics):

| Method | β_c found | raw quartic | Z⁴ | Time |
|--------|-----------|-------------|-----|------|
| **exact** (ref β_c = 0.2389) | 0.23888464 | 1.9 × 10⁻⁹ | 1.18 | 4s |
| **gc2** (current) | 0.24077417 | 9.0 × 10⁻⁹ | 7.67 | 102s |
| **gc3** (3-pass) | 0.23957956 | 1.9 × 10⁻¹⁰ | 0.08 | 138s |
| **gc2_histat** | 0.24127309 | — | — | ~200s |
| **gc2_boot** | 0.24154046 | — | — | ~102s |

**Key insight**: The "exact" control (bypassing the finder entirely) gives
cost ≈ 2 × 10⁻⁹, while the gc2 finder gives ≈ 9 × 10⁻⁹ — a 5× penalty
from β_c mismatch (finder found 0.2408 vs reference's 0.2389).  The gc3
3-pass finder found 0.2396 (closer) and achieved even lower cost.

---

### ✅ Adopted: Hybrid GC3 + Z²/Z⁴ (2025-04-05)

Based on the completed diagnostic tests, the **hybrid strategy is now the
default** in both `run_grid_search.py` and `run_grid_search_8x8.py`:

1. **GC3 (3-pass Gram-Charlier) for the β_c finder** — the extra refinement
   pass tightens β_c.  The β_c diagnostic (30 trials × 5 methods) confirmed
   that all methods converge to the same physical peak; GC3 is ~35% slower
   than GC2 but provides a tighter estimate.  FSS analysis (L=4–12) shows
   all four methods (gc2, gc3, gc2_histat, gc2_boot) give consistent
   extrapolations to β_c(∞) ≈ 0.2773–0.2778 (exact = 0.2747; the +0.003
   residual is a finite-size artefact from missing L=16 data).

2. **Normed Z² (`boundary_slices_normed`) for grid evaluation** — CoV ≈ 1.01,
   the most reliable ranking of grid points.  Correlation with Z⁴ = 0.999
   (they rank-order identically).

3. **Normed Z⁴ (`boundary_slices_normed_quartic`) re-scoring** of top 3
   candidates per level.  No extra MC runs (reuses existing two-point data).
   The Z⁴ winner determines grid refinement direction.

4. **Translate-or-refine grid logic** — border axis: re-centre on best,
   span unchanged; interior axis: re-centre on best, span halved.  No
   overshoot.

5. **Sum-of-squares per-direction aggregation** — penalises concentrated
   deviations in one direction more than spread deviations.

The **8×8 equilateral convergence test** completed 6 levels and converged to
r₁ = 0.956, r₂ = 1.003, β_c = 0.244 (expected: r₁ = r₂ ≈ 1.0).

### FSS β_c test (partial: L=4,6,8 complete, L=12 at 8/10, L=16 pending)

Finite-size scaling fit β_c(L) = β_∞ + a/L (ν=1, 2D Ising):

| Method | β_∞ | ± err | bias from ln(3)/4 |
|--------|-----|-------|--------------------|
| GC2 (standard) | 0.27732 | 0.00030 | +0.00267 |
| GC3 (3-pass) | 0.27779 | 0.00033 | +0.00314 |
| GC2 hi-stat (2×) | 0.27777 | 0.00022 | +0.00311 |
| GC2 + bootstrap | 0.27781 | 0.00041 | +0.00315 |

All methods overshoot by ~0.003 — expected from the 1/L fit without L=16.
The L=16 data (in progress) will constrain the extrapolation further.
Key finding: **all four methods agree** — the FSS shift dominates over
method differences.

---

## Changelog

### 2025-04-05

**25. Hybrid GC3 + Z²/Z⁴ becomes default in plug-and-play scripts**

`run_grid_search.py` and `run_grid_search_8x8.py` now use the hybrid
strategy directly (self-contained, no longer delegating to
`optimise_couplings.main()`):
- **GC3 β_c finder** (3-pass Gram-Charlier) for every grid point.
- **Z² (boundary_slices_normed)** for grid evaluation — lowest CoV.
- **Z⁴ (boundary_slices_normed_quartic)** re-scoring of top 3 per level.
- Translate-or-refine grid logic, sum-of-squares aggregation.
- Auto-build of the C++ simulator if binary is missing.
- Auto-generation of reference data if not found.

The `COST` and `R_MIN`/`R_MAX_FRAC` parameters have been removed from the
CONFIG section (the hybrid strategy handles cost selection internally).
The old `run_grid_search_hybrid.py` is retained as an archive.

**16. Grid refinement bug fix (translate-or-refine)**

The per-axis grid refinement logic had a bug: when the best point fell on a
border, the centre was pushed by an extra half-span *past* the best point
(overshoot), and the span was never reduced on border axes.  This caused
oscillation/overshoot visible in the 32×32 PC run (best at r₁≈0.75 instead
of 1.0, β_c drifting to 0.16).

**Fix (`optimise_couplings.py`)**: Simplified to a clean rule:
- **Border axis** → re-centre on best point, span **unchanged** (translate).
- **Interior axis** → re-centre on best point, span **halved** (refine).

No extra push, no asymmetric offsets.  The grid always centres on the
current best and only zooms in when the best is away from the edge.

**17. FSS β_c test (`run_fss_betac.py`)**

New finite-size scaling test for ground-truthing β_c finder methods.  Runs
4 methods (gc2, gc3, gc2_histat, gc2_boot) across 5 lattice sizes
(L = 4, 6, 8, 12, 16) with 10 trials each, then fits
β_c(L) = β_c(∞) + a/L^p to extrapolate to L→∞ vs the known exact
ln(3)/4 ≈ 0.2747.  Resume support included.

**18. Visualization scripts**

- `plot_betac_diagnostic.py`: 6-panel diagnostic of β_c finder comparison
  (distributions, bias, Z⁴ by method, scatter, win rates, efficiency).
- `plot_score_comparison.py`: 6-panel cost function comparison (histograms,
  CoV bars, Z² vs Z⁴ scatter, bias vs cost, β_c histogram, Z⁴/Z² ratio).

**19. Resume support for β_c diagnostic**

`run_betac_diagnostic.py` now loads existing `diagnostic_results.json` on
start and resumes from where it left off, rather than restarting from
trial 0.

**20. Hybrid GC3 + Z²/Z⁴ grid search (`run_grid_search_hybrid.py`)**

New grid search script implementing the full hybrid strategy:
- **GC3 β_c finder** (3-pass Gram-Charlier) for every grid point.
- **Z² (boundary_slices_normed)** as the grid evaluation cost — lowest CoV,
  most stable for ranking.
- **Z⁴ (boundary_slices_normed_quartic)** re-scoring of top 3 candidates
  per level — steepest landscape for final discrimination.
- Uses the fixed translate-or-refine grid logic (#16 above).
- No extra production runs needed for Z⁴ re-scoring (reuses existing
  two-point data from the Z² evaluation).

Config: 5×5 grid, 6 levels, 8×8 equilateral, starting at r₁=r₂=1.5.
Output: `results/hybrid_gc3_z2z4_8x8/`.

**21. Error bar scaling fix in LivePlotter**

The `record_point()` method in `LivePlotter` stored the raw `chi2_err`
(= √(2·total)) but plotted it as error bars on the `chi2/ndof` axis.
With `ndof = 3` (three boundary directions), the error bars were **3×
too large**.  Fixed by dividing `chi2_err` by `ndof` when storing.

**22. Variance decomposition analysis**

Regression analysis of the β_c diagnostic data (30 trials × 5 methods)
to decompose Z² cost variance without trusting any "exact" β_c value.
Using within-method linear regression Z² ~ β_c (trust-nothing approach):
- β_c jitter explains **47–66%** of Z² variance (gc2–gc3).
- MC statistics account for the remaining **34–53%**.
This is a roughly **50/50 split**, making β_c the single largest (but
not overwhelming) noise source.

**23. High-statistics hybrid test relaunch**

Restarted `run_grid_search_hybrid.py` with increased MC statistics:
- **Reference**: 100k → **1M trajectories** (10×), stored in `ref_data/ref_8x8_1M/`
- **Test production**: 50k → **100k trajectories** (2×)
- Scan statistics unchanged (50k coarse, 100k fine — used for β_c only)

Estimated impact: the normalised Z² denominator (σ²_ref + σ²_test) shrinks
by **2.7×**, making the cost landscape 2.7× sharper.  The null floor (Z² ≈ 1
at the correct point) is unchanged because the normalisation cancels the
extra statistics.  However, β_c jitter is also amplified 2.7× into Z²,
shifting the variance split to ~88% β_c / 12% MC.  GC3 is considered
adequate for the β_c finder at this stage.

Config: 5×5 grid, 6 levels, 8×8, r₁=r₂=1.5, ref=1M, test=100k.
Output: `results/hybrid_gc3_z2z4_8x8_histat/`.

**24. Sum-of-squares per-direction cost aggregation**

The Z² and Z⁴ cost functions previously summed per-direction costs linearly:
`total = v + u + w`.  This meant a large deviation in one direction (e.g.
v=30, u=0, w=0 → total=30) scored the same as moderate deviations in all
directions (v=10, u=10, w=10 → total=30).  In practice the optimizer needs
to penalise anisotropy: a large miss in one direction is a worse coupling
choice than a moderate miss spread evenly.

**Fix**: both `_boundary_slices_normed_cost` and `_boundary_slices_normed_quartic_cost`
now aggregate as **sum of squares**: `total = v² + u² + w²`.  With this:
- (10, 10, 10) → 300
- (30, 0, 0) → 900 (3× worse — concentrated deviation penalised)
- (5, 5, 20) → 450 vs (10, 10, 10) → 300

This makes the cost landscape steeper in the anisotropic directions,
improving the grid search's ability to distinguish genuinely isotropic
couplings from those that happen to average out across directions.

### 2025-04-04

**4. Gram-Charlier mode extraction bug fix**

The `_gc_fit` function used `minimize_scalar` over the full scan range to
find the susceptibility peak (mode) of the GC polynomial.  The unconstrained
search was sometimes attracted to **spurious secondary maxima** created by
the skew/kurtosis terms of the polynomial far from the true peak (e.g.
returning β_c ≈ 0.326 at the scan boundary instead of the correct ≈ 0.235).
**Fix**: the mode search is now constrained to **±3σ of the fitted μ**,
eliminating false peaks while still accommodating the genuine GC mode shift.

**5. Error bars on boundary_slices cost**

`_boundary_slices_cost` now computes and returns a proper uncertainty
estimate.  The MC errors on the reference and test two-point functions are
propagated through the quartic integrand via:

$$\sigma_{\mathrm{int}} = 4 \, |G_{\mathrm{test}} - G_{\mathrm{ref}}|^3 \; \sqrt{\sigma_{\mathrm{ref}}^2 + \sigma_{\mathrm{test}}^2}$$

The variance is accumulated through the trapezoidal integration and the
total error is returned as √(∑ var).  Previously the error was hard-coded
to 0.0.  The error bars appear automatically in the fit-quality (Panel 2)
plot when non-zero.

**6. Log-scaled fit-quality plot**

The Panel 2 y-axis ("Fit quality") now uses `ax.set_yscale("log")` in all
cases, including when error bars are present.  Previously `semilogy` was
only used in the no-error-bar branch, and `errorbar` plotted on a linear
axis which made the tiny boundary_slices costs (∼10⁻⁸) invisible.

**7. Scientific notation for tiny chi²/ndof values**

When chi²/ndof < 10⁻³ the log output now prints in scientific notation
(`:.4e`) instead of fixed-point (`:.4f`), so boundary_slices costs are
visible instead of rounding to `0.0000`.

**8. 8×8 fast test script (`run_grid_search_8x8.py`)**

New self-contained entry point for fast 8×8 tests.  Auto-generates and
caches its own 8×8 equilateral reference (100 k trajectories).  The 8×8
lattice is ~16× faster per MC sweep than 32×32, making it useful for rapid
convergence checks.

**9. Stability test script (`run_stability_test.py`)**

New script that evaluates the same (r₁, r₂) point N times (default 100)
under independent MC noise, running the full pipeline each time.  Produces
JSON results and histogram plots showing the distribution of β_c and cost
values, allowing quantification of the optimizer's signal-to-noise ratio.

**10. Exact β_c for reference generation (`--ref_beta_c`)**

New `--ref_beta_c` CLI argument (exposed as `REF_BETA_C` in the run
scripts).  When provided, reference generation skips the susceptibility
scan entirely and uses the supplied value as β_c.  This is useful when the
exact critical point is known analytically — e.g. ln(3)/4 ≈ 0.2747 for
the equilateral triangular Ising model — avoiding any finite-size shift
from the susceptibility peak.

**11. Variance-normalised quadratic cost (`boundary_slices_normed`)**

New cost function: ∫₀¹ Z(t)² dt, where Z = (G_test − G_ref) / σ and
σ² = σ_ref² + σ_test².  Under the null (test = ref) the integrand has
expectation 1 per sample point, making the total cost independent of MC
statistics.  Selected via `--cost boundary_slices_normed`.

**12. Variance-normalised quartic cost (`boundary_slices_normed_quartic`)**

New cost function: ∫₀¹ Z(t)⁴ dt with the same Z = ΔG/σ normalisation.
Combines variance normalisation with the steep quartic landscape.  For
signal s = δ/σ, the excess over null scales as 6s² + s⁴ (vs s² for the
quadratic), giving ~44% better discrimination at s = 2 and ~116% at s = 3.
Null expectation is E[Z⁴] = 3 per point with std = √96.  Selected via
`--cost boundary_slices_normed_quartic`.

**Motivation**: The raw quartic cost (diff⁴) was found to have
CoV ≈ 1.6 in stability tests — two evaluations of the same point can
differ by 1000×.  The normalisation was expected to fix this, but early
results (15 trials) show CoV ≈ 1.1 (normed Z²) and CoV ≈ 2.1 (normed Z⁴),
indicating that the dominant noise source is β_c finder jitter, not the
cost function per se.

**13. Normed cost comparison script (`run_stability_normed_compare.py`)**

Side-by-side paired comparison of all three boundary-slice cost functions
on the same MC data.  100 trials at (r₁=r₂=1, 8×8), producing histograms
and paired scatter plots.  Currently running.

**14. β_c finder diagnostic script (`run_betac_diagnostic.py`)**

Comprehensive diagnostic that tests 5 β_c finder methods side-by-side:
exact (control), gc2 (current), gc3 (3-pass), gc2_histat (2× scan stats),
and gc2_boot (bootstrap median).  30 trials per method interleaved for
fairness, evaluating all three cost functions on each trial.  Produces
β_c distribution plots, CoV bar charts, and cost distribution grids.
Currently running.

**15. Recommendation section (flagged for discussion)**

Added a "⚠️ RECOMMENDATION FOR DISCUSSION" section to the test-results
area of this README, proposing gc3 + normed Z² for grid search and normed
Z⁴ for final verification only.  Based on 20 paired trials showing
log₁₀(Z⁴) vs log₁₀(Z²) correlation = 0.999 and Z² CoV (1.01) roughly
half of Z⁴ (1.87).  The gc3 3-pass finder is not yet integrated into
`find_beta_c` — currently only implemented in `run_betac_diagnostic.py`.

### 2025-04-03

**1. Quartic boundary-slice cost (diff⁴)**

The `boundary_slices` cost function now uses a **quartic integrand**:
∫₀¹ (G_test − G_ref)⁴ dt, replacing the previous quadratic (diff²).
This strongly penalises large pointwise deviations: a deviation twice as
large now contributes 16× more to the cost instead of 4×.  Per-direction
costs (v, u, w) are printed after each evaluation for diagnostics.

**2. Per-dimension adaptive grid search**

The outer grid search now maintains **independent half-spans** `hs_r1` and
`hs_r2` for the two coupling-ratio axes.  At each level, the border/interior
decision is made independently per dimension:
- Interior axis → halve its half-span (zoom in).
- Border axis → translate the centre by an extra half-span in the border
  direction, keeping the span unchanged.

Previously a single `half_span` was shared and any border detection blocked
refinement in both axes simultaneously.  The new logic means a side-border
minimum correctly refines the interior axis while translating the border
axis; a corner minimum translates both.

**3. More robust β_c finder**

- **Wider initial scan window**: margin increased from 15% to 20% of β_guess
  (minimum 0.04 instead of 0.02), reducing the chance of the peak falling
  outside the initial range.
- **Denser coarse scan**: 11 points per sweep (was 7), increasing the chance
  of bracketing a narrow susceptibility peak.
- **More refinement points**: 5 + 5 fine-statistics Gram-Charlier refinement
  passes (was 3 + 3).
- **All-data GC fits**: both refinement passes now fit the Gram-Charlier
  envelope to ALL accumulated scan data (coarse + fine), not just the latest
  narrow window.  This dramatically reduces the probability of a "missed
  peak" when the fine window lands off-centre.
- **Weighted-top-3 parabolic fallback**: if the Gram-Charlier fit fails
  entirely, the susceptibility-weighted centroid of the top 3 scan points is
  used instead of the raw argmax of the fine window.

---

## Citation / acknowledgement

If you use this code, please cite the newQFE project (github.com/brower/newQFE)
and the relevant paper(s) on quantum finite elements / triangular Ising FEM.
