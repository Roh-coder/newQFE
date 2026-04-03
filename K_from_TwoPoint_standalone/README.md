# K_from_TwoPoint — Standalone Package

Find the critical anisotropic couplings (k₁, k₂, k₃) of the 2D Ising model
on a triangular lattice by matching the all-to-all two-point function to a
high-statistics equilateral reference via an adaptive grid search.

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
| `COST` | Cost function (`"fem_integral"` recommended) | `"fem_integral"` |
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

4. **Cost**: the `fem_integral` cost computes the exact analytical integral
   of (G_test − G_ref)² over the periodic lattice domain using a bilinear
   FEM representation.  When this equals zero the two correlators agree
   everywhere and the couplings are critical and isotropic.

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
