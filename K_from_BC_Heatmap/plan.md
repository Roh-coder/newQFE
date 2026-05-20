# Continuum Boundary-Correlator Optimization Plan

## 2026-05-12 — Zip-Folder Standalone Workflow Plan

### Goal
Build a self-contained folder that can be zipped, copied to another machine, unzipped, opened in Spyder, and run offline with no dependency on the surrounding repo tree.

The workflow should:
- use the local C++ Monte Carlo source file `ising_tri_twisted_parallelogram.cc` as the MC engine,
- compile that source inside the standalone folder,
- expose a single plug-and-play Python runner for Spyder,
- allow the user to reconfigure all relevant run parameters from that runner,
- show a live visual monitor during the run,
- finish by showing the final score heatmap and RMS z-score heatmap,
- save all outputs inside the standalone folder's own `results/` tree.

### Grounded starting point
The correct pattern is the existing zip-and-go layout already present in `fss_standalone/`.

Important evidence:
- `fss_standalone/README.md` already describes a standalone zipped workflow.
- `fss_standalone/Makefile` builds the local simulator from source.
- `fss_standalone/src/ising_tri_twisted_parallelogram.cc` is the MC source file to use.
- `run_boundary_continuum_opt.py` already computes the final boundary score and RMS z-score products we want.

So the plan is not to invent a new MC path. The plan is to package the existing boundary-continuum analysis into a self-contained standalone folder whose MC binary is built from the local `.cc` source.

### Deliverable
Add one new standalone folder, tentatively:
- `boundary_opt_standalone/`

This folder should be zip-safe and include everything needed except standard Python packages and a C++ compiler.

Proposed contents:
- `Makefile`
- `src/ising_tri_twisted_parallelogram.cc`
- `include/*.h`
- `lib/mc_engine.py`
- `lib/cost.py`
- `run_boundary_opt.py`
- `config.example.py` or top-level editable `CONFIG` block inside `run_boundary_opt.py`
- `README.md`
- optional prebuilt `bin/ising_tri_twisted_parallelogram` for convenience, but source remains authoritative.

### User-facing workflow
The intended offline workflow should be:

1. zip the standalone folder,
2. move it to another machine,
3. unzip it anywhere,
4. open `run_boundary_opt.py` in Spyder,
5. edit the top-level config,
6. press Run,
7. let the script build the MC executable from `src/ising_tri_twisted_parallelogram.cc` if needed,
8. watch the live monitor,
9. inspect the final score and z-score heatmaps on screen,
10. keep the saved outputs under `results/<tag>/`.

### What “uses the .cc file for the MC” means here
For this workflow, the source of truth for MC must be the local C++ source inside the standalone folder.

That means:
- the standalone script must not depend on the parent repo binary in `K_from_BC_Heatmap/bin/`,
- the standalone script must prefer the local executable built from local source,
- if the executable is missing, stale, or on a new platform, the workflow should rebuild locally.

Recommended behavior:
- on startup, `run_boundary_opt.py` checks whether `bin/ising_tri_twisted_parallelogram` exists,
- if not, or if the source is newer, it runs `make`,
- if build fails, it prints a clear actionable error in Spyder.

### Scope for v1
Included in v1:
- one self-contained standalone folder,
- one Spyder-friendly runner script,
- local Makefile-based build of the C++ simulator,
- no imports from the parent repo,
- full precompute plus continuum analysis,
- live matplotlib progress monitor,
- final on-screen score and z-score figures,
- saved `.dat`, `.json`, and `.png` outputs.

Explicitly out of scope for v1:
- GUI installers,
- Conda environment export,
- platform-specific packaged executables,
- deep C++ refactoring,
- in-scan per-sweep animation from inside the simulator.

### Standalone folder design
The standalone directory should mirror the successful `fss_standalone` structure, but target the boundary-optimization use case.

Suggested layout:

- `Makefile`
  builds `bin/ising_tri_twisted_parallelogram` from local source.

- `src/ising_tri_twisted_parallelogram.cc`
  the exact MC source file used for production sampling.

- `include/`
  local headers required by the simulator build.

- `lib/mc_engine.py`
  local beta-c finder and simulator wrapper.

- `lib/cost.py`
  local geometry utilities and boundary-path interpolation helpers.

- `run_boundary_opt.py`
  main Spyder entrypoint with editable config and live monitor.

- `README.md`
  short zip-folder instructions: build, Python packages, run in Spyder.

### Python dependency policy
The standalone folder should be self-contained with respect to project code, but it is acceptable to require standard Python packages installed on the destination machine.

Expected external Python packages:
- NumPy
- SciPy
- Matplotlib

Expected external toolchain:
- a C++14 compiler,
- `make` or an equivalent POSIX-style build tool.

The README should state these explicitly.

### Main Python runner design
The user still wants a plug-and-play Python file for Spyder, so the standalone workflow should keep a single obvious entrypoint:
- `run_boundary_opt.py`

This file should:
- contain a compact top-level `CONFIG` block,
- be executable as a normal script,
- also support `main(config=None)` for reuse,
- print clear progress lines in the Spyder console,
- show live plots using matplotlib interactive mode.

### Runner reconfiguration requirements
The standalone runner must let the user reconfigure all parameters that materially affect the run, without editing internal helper code.

At minimum the top-level configuration should expose:
- bookkeeping:
  - `tag`
  - `results_dir`
  - `resume`
  - `auto_build`
- coupling grid:
  - `r_min`, `r_max`, `r_step`
  - optional explicit `r1_values`, `r2_values`
  - `k3`
- test family:
  - `test_sizes`
  - test geometry policy or explicit map
- reference family:
  - `reference_mode`
  - `ref_sizes`
  - `ref_geom_map`
  - `ref_large`
  - `ref_fixed_r1`, `ref_fixed_r2`
- Monte Carlo controls:
  - `n_traj`
  - `n_traj_coarse`
  - `n_traj_fine`
  - `n_skip`
  - `n_therm`
  - `beta_seed`
  - beta scan point counts for each pass
  - `max_shifts`
  - `jackknife`
- execution:
  - `n_workers`
  - optional `exe`
  - `build_command`
- scoring and analysis:
  - `weighted`
  - sampled boundary fractions or `k_values`
  - continuum fit mode if multiple are supported
- monitor and plotting:
  - `show_live_monitor`
  - `monitor_refresh_seconds`
  - `save_monitor_snapshot`
  - figure sizing or dpi if needed.

Preferred design:
- keep a readable default `CONFIG` dictionary inside `run_boundary_opt.py`,
- optionally support loading an external config file for larger runs,
- validate config fields up front and fail early with readable messages.

### Live monitor design
The live monitor should stay matplotlib-only so it behaves naturally inside Spyder.

Recommended dashboard: 4 panels.

Panel 1: run summary
- stage name,
- completed jobs / total jobs,
- completed score cells / total cells,
- ok / cached / error counts,
- elapsed time and ETA.

Panel 2: per-size completion map
- test sizes by reference/test kind,
- a simple matrix or grouped bars showing what has completed.

Panel 3: partial score heatmap
- filled progressively as `(r1, r2)` analysis completes,
- current best point marked with a star.

Panel 4: partial RMS z-score heatmap
- same progressive fill,
- current best z-score point marked with a star.

Refresh cadence for v1:
- after each completed MC job,
- after each completed `(r1, r2)` analysis cell.

### Final figures to display
At the end of the run, the script should show:

1. final score heatmap,
2. final RMS z-score heatmap,
3. optionally the last monitor dashboard snapshot.

These should also be saved to disk.

### Output contract
The standalone workflow should preserve the useful existing boundary-analysis outputs.

Keep:
- `manifest_boundary_opt.json`
- `dat/raw_test_channels.dat`
- `dat/raw_ref_channels.dat`
- `dat/continuum_test.dat`
- `dat/continuum_ref.dat`
- `dat/score_map.dat`
- `dat/zscore_map.dat`
- `dat/score_minimum.dat`
- `score_heatmap.png`
- `score_and_zscore_heatmaps.png`

Add if useful:
- `monitor_snapshot.png`
- `progress_log.jsonl`

### Implementation approach
The build should use two existing anchors:
- the standalone packaging pattern from `fss_standalone/`,
- the boundary-continuum numerical logic from `run_boundary_continuum_opt.py`.

Step 1: create standalone folder shell
- copy or vendor the required `src/`, `include/`, and `lib/` files into the new standalone directory,
- add a local Makefile,
- ensure nothing imports from outside the new folder.

Step 2: port the boundary runner
- create `run_boundary_opt.py` in the standalone folder,
- reuse the validated analysis path for score and z-score maps,
- replace CLI-only assumptions with a Spyder-friendly config block.

Step 3: add build bootstrap
- add a helper that detects whether the simulator needs to be built,
- run `make` automatically when needed,
- surface build errors clearly.

Step 4: add live monitor
- add a small `LiveMonitor` helper using matplotlib,
- update after completed jobs and completed score cells,
- keep the monitor optional but on by default.

Step 5: zip-folder validation
- test from the standalone directory only,
- verify local build from source,
- verify no accidental dependency on the parent repo tree.

### Key engineering choices
#### Choice A: use local source build, not parent binary
This is required by the workflow request and makes the zip folder portable.

#### Choice B: keep one obvious Spyder entry script
The user should not need to orchestrate multiple scripts by hand.

#### Choice C: keep the monitor in matplotlib
This keeps the workflow simple and editor-friendly.

### Risks and mitigations
Risk: local build may fail on machines without compiler tooling.
Mitigation:
- detect this early,
- print a short build troubleshooting section in the README,
- optionally allow shipping a convenience binary while still keeping local source authoritative.

Risk: multiprocessing plus interactive plotting can be fragile in some Spyder setups.
Mitigation:
- keep `n_workers = 1` as a debug-safe mode,
- make refresh cadence modest,
- allow turning the monitor off.

Risk: accidental hidden imports from the main repo would break portability.
Mitigation:
- vendor all needed Python modules into the standalone folder,
- validate by running only from the standalone directory.

### Acceptance criteria
This task is complete when:
1. there is one standalone folder that can be zipped and moved independently,
2. that folder contains the local `ising_tri_twisted_parallelogram.cc` source and local build files,
3. the Python runner uses the executable built from that local source,
4. the runner can be opened and run directly in Spyder,
5. the runner exposes all relevant run parameters in one user-editable configuration surface,
6. a live progress monitor updates during the run,
7. the final displayed plots include the score heatmap and RMS z-score heatmap,
8. the workflow runs without importing code from the parent repo.

### Review gate before build
Before implementation, confirm these scope decisions:
- build a new zip-safe standalone folder, not just a new file inside the existing tree,
- compile and use the local `ising_tri_twisted_parallelogram.cc` source for MC,
- keep one Spyder-friendly Python entrypoint in that folder,
- preserve the current boundary score/z-score outputs and add live monitoring on top.

## 1. Objective
Build a production pipeline that scans coupling-ratio space `(r1, r2)` and computes a continuum-level mismatch score between:
- a test lattice family (untwisted, multiple sizes), and
- a user-defined twisted reference lattice family evaluated once at fixed isotropic couplings and reused across all `(r1, r2)` points.

The output is a smooth score manifold over `(r1, r2)` with a clear minimum for geometry-to-coupling calibration.

## 2. Parameter-Space and Geometry Definition

### 2.1 Coupling scan
- Scan a rectangular grid in `(r1, r2)`.
- Use fixed `k3` (default `1.0`) and set:
	- `k1 = r1`
	- `k2 = r2`

### 2.2 Test lattices (fixed production set)
- Test sizes:
	- `L_test = [8, 16, 24, 32, 40]`
- Test geometry default:
	- untwisted square: `(Lx, Ly, Tx, Ty) = (L, L, 0, 0)`

### 2.3 Reference lattices (user-configurable)
Allow user input for a list of reference geometry levels:
- `(Lx_ref_i, Ly_ref_i, Tx_ref_i, Ty_ref_i)` for each level `i`
- Optional independent reference size list `L_ref`

Reference coupling policy (important):
- Reference jobs are run at fixed isotropic couplings (`k1=k2=k3`, default ratio `1:1:1`).
- Reference channels depend on geometry and size, not on `(r1, r2)`.
- Therefore, reference payloads are generated once per reference level and reused for the full `(r1, r2)` scan.

Two operation modes:
1. **Reference continuum mode**: run multiple reference sizes and extrapolate to continuum.
2. **Reference large-L mode**: run a single very large twisted reference lattice and use it as a near-continuum proxy.

## 3. Observable Construction

### 3.1 Boundary-direction correlators
For each lattice realization, sample connected correlators along all 3 boundary cycles at 8 fractional positions:
- `t_k = k/8`, `k = 0..7`

Channel index set:
- cycle index `c in {0,1,2}`
- fraction index `k in {0,...,7}`

Stored channels:
- `G_test(c, k, L)`
- `G_ref(c, k, L)` (or `G_ref_large(c, k)` in large-L mode), reused across all `(r1, r2)`

### 3.2 Critical-point policy
For test lattices, at each size and `(r1, r2)` point:
- locate pseudo-critical `beta_c(L)` via the existing Gram-Charlier scan, then
- run production MC at `beta_c(L)` to measure correlators.

For reference lattices, once per reference size/geometry level:
- locate pseudo-critical `beta_c_ref(L_ref)` at isotropic couplings,
- run production MC once and cache channels for reuse across all scan points.

## 4. Continuum Extrapolation

For each channel `(c, k)`, extrapolate in `x = 1/L` using quadratic form:

`G(L) = a + b x + c x^2`

Continuum value is the intercept:
- `G_inf = a`

Apply this to:
- test channels: `G_test_inf(c,k)`
- reference channels (continuum mode): `G_ref_inf(c,k)`

In reference large-L mode:
- set `G_ref_inf(c,k) := G_ref_large(c,k)`

## 5. Optimization Score Definition

Residual per channel:

`R(c,k) = G_test_inf(c,k) - G_ref_inf(c,k)`

Point score at each `(r1, r2)`:

`S(r1,r2) = sum_{c=0..2} sum_{k=0..7} w(c,k) * R(c,k)^2`

Default:
- uniform weights `w(c,k)=1`

Optional production weighting:
- inverse-variance weights `w(c,k)=1/sigma(c,k)^2` when propagated continuum uncertainties are available.

## 6. Production Workflow

### Phase A: Precompute data grid
For every `(r1, r2)` point:
1. Run all test sizes in `L_test`.
2. Reuse the precomputed reference set (either multi-size continuum mode or one-shot large-L mode), generated once at isotropic couplings.
3. Cache all raw outputs in a resumable directory layout.

### Phase B: Build continuum datasets
1. For each channel `(c,k)`, fit test vs `1/L` quadratically and store `G_test_inf(c,k)`.
2. Build reference continuum channels similarly, or ingest large-L proxy channels.
3. Save continuum tensors and uncertainty fields.

### Phase C: Compute score manifold
1. Evaluate `S(r1,r2)` on the full grid.
2. Save `score_map` and minimum location.
3. Generate heatmaps and residual-diagnostic plots around low-score basin.

### Phase D: Refinement near best basin
1. Locally refine `(r1, r2)` grid near minimum.
2. Increase MC stats only in refined region.
3. Recompute and report final optimum with confidence intervals.

## 7. Data Products
For each production tag `results/<tag>/`:
- `manifest.json`: full run configuration
- `grid/...`: cached per-size per-geometry payloads
- `continuum_test.npz`: `G_test_inf(c,k,r1,r2)` and uncertainty
- `continuum_ref.npz`: `G_ref_inf(c,k)` and uncertainty (or large-L snapshots), independent of `(r1,r2)`
- `score_map.npz`: `S(r1,r2)` and metadata
- `plots/score_heatmap.png`: score manifold with marked minimum
- `plots/residual_channels_*.png`: channel-level residual structure

## 8. Operational Notes
- Keep caching and resumability as first-class behavior to support long runs.
- Maintain fixed MC policy across the scan to avoid introducing procedural bias into the manifold.
- Prefer reference continuum mode for highest fidelity; use large-L reference mode when budget is constrained.
- Use local refinement around the coarse minimum to improve final parameter precision at modest additional cost.

## 9. Success Criteria
The production run is complete when:
1. a stable minimum in `S(r1,r2)` is identified,
2. local refinement preserves the same basin,
3. final best-fit `(r1, r2)` and uncertainty are reported with full run metadata.
