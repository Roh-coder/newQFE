# boundary_opt_standalone

Self-contained zip-and-go boundary-continuum optimization workflow.

This folder is designed to be zipped, moved to another machine, unzipped,
opened in Spyder, and run offline. The Monte Carlo engine is built from the
local C++ source in `src/ising_tri_twisted_parallelogram.cc`.

## Contents

```text
Makefile
src/ising_tri_twisted_parallelogram.cc
include/
lib/mc_engine.py
lib/cost.py
run_boundary_opt.py
configs/
results/
```

## External requirements

- Python packages: `numpy`, `scipy`, `matplotlib`
- Toolchain: `make` and a C++14 compiler (`g++` or similar)

## Quick start in Spyder

1. Open `run_boundary_opt.py` in Spyder.
2. Edit the hard-coded variables near the top of the file (for example
  `RUN_TAG`, `REFERENCE_MODE`, `TEST_SIZES`, `REFERENCE_LARGE_GEOM`), or set
  `CONFIG_FILE` to one of the JSON files in `configs/`.
3. Press Run.
4. The script will build the local simulator with `make` if needed.
5. Watch the live progress monitor and inspect the final score and RMS
   z-score heatmaps.

## Quick start in a terminal

```bash
make
/workspaces/newQFE/.venv/bin/python run_boundary_opt.py --config configs/smoke_continuum.json
```

Useful flags:

```bash
/workspaces/newQFE/.venv/bin/python run_boundary_opt.py \
  --config configs/production_validation.json \
  --tag prod_validation \
  --no-monitor \
  --no-show-final
```

Reference-strategy override flags:

```bash
/workspaces/newQFE/.venv/bin/python run_boundary_opt.py \
  --config configs/production_run.json \
  --reference-mode continuum

/workspaces/newQFE/.venv/bin/python run_boundary_opt.py \
  --config configs/production_run.json \
  --reference-mode large \
  --reference-large-geom 80 80 0 0
```

## Configuration surface

The runner exposes all relevant parameters through the configuration surface.

Top-level sections:

- `run`
  - bookkeeping such as `tag`
- `paths`
  - `results_dir`, `resume`, progress log name
- `execution`
  - local executable override, auto-build policy, build command, worker count
- `grid`
  - `r_min`, `r_max`, `r_step`, optional explicit `r1_values` and `r2_values`, `k3`
- `test_family`
  - FSS sizes and either geometry defaults or explicit geometry map
- `reference_family`
  - reference mode (`continuum` or `large`), reference sizes, geometry map,
    large reference geometry, fixed isotropic `(r1, r2)` for reference payloads
- `mc`
  - production MC controls: `n_traj`, `n_skip`, `n_therm`
- `beta_c_finder`
  - scan bracket, point counts, coarse/fine trajectories, bracket shifts,
    jackknife toggle
- `analysis`
  - `k_values`, `k_denominator`, `weighted`, `fit_mode`
- `monitor`
  - live plotting, refresh cadence, snapshot saving, figure sizes, dpi

## Output layout

Outputs land under `results/<tag>/`.

Important files:

- `manifest_boundary_opt.json`
- `progress_log.jsonl`
- `dat/raw_test_channels.dat`
- `dat/raw_ref_channels.dat`
- `dat/continuum_test.dat`
- `dat/continuum_ref.dat`
- `dat/score_map.dat`
- `dat/zscore_map.dat`
- `dat/score_minimum.dat`
- `score_heatmap.png`
- `score_and_zscore_heatmaps.png`
- `monitor_snapshot.png`

## Reproducible configs shipped here

- `configs/smoke_continuum.json`
  - small continuum-reference smoke run
- `configs/smoke_large_ref.json`
  - small large-reference smoke run
- `configs/production_validation.json`
  - moderate production-style validation run used during development

## Choosing The Reference Workflow

The runner supports two reference workflows behind the same `reference_family`
section.

- `mode: "continuum"`
  - build one fixed-isotropic reference lattice for every `L` in
    `reference_family.sizes`
  - sample channels on each reference size
  - fit those channels in `1/L` and compare against the continuum intercept

- `mode: "large"`
  - build one fixed-isotropic reference lattice at
    `reference_family.large_geom = [Lx, Ly, Tx, Ty]`
  - sample that single lattice directly
  - compare test data against that single large reference without a reference-side
    continuum fit

You can choose the workflow in either of two ways:

1. Edit `reference_family.mode` in your JSON config.
2. Keep one config file and switch at launch time with `--reference-mode`.

Recommended pattern:

- keep your test grid, test sizes, MC controls, and channel set fixed in the
  JSON file
- set both `reference_family.sizes` and `reference_family.large_geom` in that
  same file
- choose `--reference-mode continuum` when you want a reference extrapolation
- choose `--reference-mode large` when you want a single large reference run

You can also define explicit per-size lattice tables (full parameters exposed)
directly in the script defaults or in JSON:

- test table format: `[[Lx, Ly, Tx, Ty], ...]`
- reference table format: `[[Lx, Ly, Tx, Ty], ...]`

For these explicit lattice tables, the internal fit-size `L` is inferred as
`max(|Lx|, |Ly|)`.

When exactly one reference lattice is defined, the runner automatically uses
large-reference mode and uses that single geometry as `large_geom`.

If you use `--reference-mode large`, you can also override the large lattice on
the command line with `--reference-large-geom LX LY TX TY`.

## Notes

- The runner does not import project code from outside this folder.
- The local source tree is authoritative for the MC engine.
- If you modify `src/` or `include/`, the runner will rebuild before running.