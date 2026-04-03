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
│   └── ref_metadata.json          ← beta_c and path for the reference
│
├── optimise_couplings.py          ← optimisation engine (do not edit)
├── run_grid_search.py             ← USER ENTRY POINT — edit CONFIG here
├── gen_ref.py                     ← optional: generate your own reference data
├── plot_boundary_slices.py        ← visualise boundary-direction correlator slices
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
`USER CONFIGURATION` section at the top.  The main parameters are:

| Parameter | What it does | Default |
|-----------|--------------|---------|
| `Lx`, `Ly` | Lattice size | 32, 32 |
| `Tx`, `Ty` | Twist vector components | 0, 0 |
| `R1_INIT`, `R2_INIT` | Starting grid centre (r₁ = k₁/k₃, r₂ = k₂/k₃) | 1.0, 1.0 |
| `HALF_SPAN` | Half-width of the search grid | 0.15 |
| `BETA_INIT` | Initial guess for critical β | 0.27 |
| `N_TRAJ_PROD` | MC trajectories for each production run | 50 000 |
| `N_TRAJ_SCAN_COARSE` | MC trajectories per β point (coarse scan) | 30 000 |
| `N_TRAJ_SCAN_FINE` | MC trajectories per β point (fine scan) | 60 000 |
| `COST` | Cost function (see [Cost functions](#cost-functions) below) | `"boundary_slices"` |
| `MAX_ITER` | Maximum grid refinement levels | 4 |
| `OUTPUT_DIR` | Where results are written | `"results/my_run"` |

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
N_TRAJ_PROD        = 5000
N_TRAJ_SCAN_COARSE = 5000
N_TRAJ_SCAN_FINE   = 10000
MAX_ITER           = 2
OUTPUT_DIR         = "results/quick_test"
```

This should complete in a few minutes on a modern laptop.

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
   (r₁ = k₁/k₃, r₂ = k₂/k₃) on a 5×5 grid.  After each level:
   - If the minimum is in the interior → grid is halved and re-centred.
   - If the minimum is on the border → grid is translated toward it.

3. **Inner loop**: for each grid point (r₁, r₂), the critical β is found by
   locating the susceptibility peak (coarse grid scan + two Gram-Charlier
   refinement passes).  A production run at that β_c then measures G_test(m, n).

4. **Cost**: the cost function quantifies how well the test correlator matches
   the reference.  See [Cost functions](#cost-functions) below.

---

## Cost functions

| Name | Description |
|------|-------------|
| `"boundary_slices"` | **Recommended.** Sum of ∫₀¹ (G_test(t) − G_ref(t))² dt along the three torus boundary paths v, u, w. Directly targets the physically meaningful signal. Supports cross-geometry mode. |
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
3. Computing ∫(G_test − G_ref)² dt for each path and summing.

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

## Citation / acknowledgement

If you use this code, please cite the newQFE project (github.com/brower/newQFE)
and the relevant paper(s) on quantum finite elements / triangular Ising FEM.
