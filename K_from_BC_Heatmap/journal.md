# Journal

## 2026-05-12 — Standalone zip-folder boundary optimizer implemented and validated

### Objective
Build a self-contained boundary-optimization package that can be zipped,
moved off-repo, opened in Spyder, and run against a local build of the
twisted-parallelogram Ising Monte Carlo source.

### What was implemented
Added a new standalone folder:
- `K_from_BC_Heatmap/boundary_opt_standalone/`

Key design points:
- Monte Carlo is built locally from:
   - `boundary_opt_standalone/src/ising_tri_twisted_parallelogram.cc`
- Local build is handled by:
   - `boundary_opt_standalone/Makefile`
- The main user entrypoint is:
   - `boundary_opt_standalone/run_boundary_opt.py`
- The runner is self-contained with local helper imports only:
   - `boundary_opt_standalone/lib/mc_engine.py`
   - `boundary_opt_standalone/lib/cost.py`
- The runner exposes a comprehensive configuration surface via:
   - top-level `CONFIG` in `run_boundary_opt.py`
   - optional JSON override files under:
      - `boundary_opt_standalone/configs/`

### Runner capabilities
The standalone runner now supports user reconfiguration of all run-critical knobs from one place, including:
- output tag and results location,
- resume behavior,
- auto-build and executable override,
- `(r1, r2)` grid definition,
- test and reference geometry families,
- reference mode (`continuum` or `large`),
- fixed isotropic reference couplings,
- production MC controls (`n_traj`, `n_skip`, `n_therm`),
- beta-scan controls (bracket, point counts, coarse/fine trajectories, bracket shifts, jackknife),
- channel sampling (`k_values`, denominator),
- fit mode (`linear`, `quadratic`, `auto`),
- weighted vs unweighted score,
- live monitor behavior and plot sizing.

### Monitoring and outputs
The runner now provides:
- a live matplotlib dashboard during the run when enabled,
- final score and RMS z-score heatmaps,
- a progress event log:
   - `results/<tag>/progress_log.jsonl`
- a monitor snapshot when the monitor path is enabled:
   - `results/<tag>/monitor_snapshot.png`

The standalone output contract matches the main boundary workflow:
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

### Configs added for validation and reuse
Reusable JSON configs were added under:
- `boundary_opt_standalone/configs/smoke_continuum.json`
- `boundary_opt_standalone/configs/smoke_large_ref.json`
- `boundary_opt_standalone/configs/production_validation.json`

### Validation run 1 — continuum-reference smoke
Command class:
- standalone runner with `configs/smoke_continuum.json`

Observed result:
- local build path used the standalone executable under:
   - `boundary_opt_standalone/bin/ising_tri_twisted_parallelogram`
- jobs total: `15`
- jobs ok: `15`
- jobs err: `0`
- analysis cells: `4/4`
- wall time: `14.191 s`

Best coarse point:
- `(r1, r2) = (0.9, 1.1)`
- `score_min = 24.7781902696`

Artifacts written under:
- `boundary_opt_standalone/results/smoke_continuum/`

### Validation run 2 — large-reference smoke
Command class:
- standalone runner with `configs/smoke_large_ref.json`

Observed result:
- jobs total: `9`
- jobs ok: `9`
- jobs err: `0`
- analysis cells: `4/4`
- wall time: `8.507 s`

Best coarse point:
- `(r1, r2) = (0.9, 1.1)`
- `score_min = 1.16646557147`

Artifacts written under:
- `boundary_opt_standalone/results/smoke_large_ref/`

### Validation run 3 — resume/caching smoke
Command class:
- same continuum smoke config run twice under tag:
   - `smoke_continuum_resume_check`

Observed result for the second invocation:
- jobs total: `15`
- jobs ok: `0`
- jobs skip: `15`
- jobs err: `0`
- analysis cells recomputed: `4/4`
- wall time: `1.966 s`

Interpretation:
- resumable caching works as intended for the precompute stage.

Artifacts written under:
- `boundary_opt_standalone/results/smoke_continuum_resume_check/`

### Validation run 4 — live-monitor path
Command class:
- large-reference smoke config with live monitor enabled and final display disabled

Observed result:
- jobs total: `9`
- jobs ok: `9`
- jobs err: `0`
- analysis cells: `4/4`
- wall time: `11.718 s`
- monitor snapshot successfully written:
   - `boundary_opt_standalone/results/smoke_large_ref_monitor/monitor_snapshot.png`

Best coarse point:
- `(r1, r2) = (1.1, 1.1)`
- `score_min = 5.17464849590`

Interpretation:
- the monitor code path executes successfully and produces a saved dashboard snapshot.

### Production-style validation run
Command class:
- standalone runner with `configs/production_validation.json`

Scope:
- grid: `3 x 3` points over `(r1, r2) in {0.8, 1.0, 1.2}^2`
- test sizes: `[4, 6, 8]`
- reference mode: `continuum`
- reference sizes: `[4, 6, 8]`
- production trajectories: `250`
- beta-scan trajectories: coarse `30`, fine `60`

Observed result:
- jobs total: `30`
- jobs ok: `30`
- jobs err: `0`
- analysis cells: `9/9`
- wall time: `25.329 s`

Best coarse point in this validation run:
- `(r1, r2) = (0.8, 1.2)`
- `score_min = 6.17549008262`

Artifacts written under:
- `boundary_opt_standalone/results/production_validation/`

### Conclusions
1. The standalone zip-folder workflow is functional.
2. The runner uses the local standalone C++ source build, not the parent repo binary.
3. Both reference modes (`continuum` and `large`) run successfully.
4. Resume/caching behavior is working.
5. The live monitor path is implemented and produces a dashboard snapshot.
6. A moderate production-style validation run completed cleanly end to end.

### Remaining caveat
The validation runs above are intentionally small enough to execute quickly in development. They validate packaging, configurability, build-from-source behavior, monitoring, and end-to-end analysis flow. They are not substitutes for a high-statistics physics production campaign.

---

## 2026-05-12 — Standalone production run executed

### Run context
Executed a real standalone production run using the new zip-folder workflow:
- runner:
   - `K_from_BC_Heatmap/boundary_opt_standalone/run_boundary_opt.py`
- config:
   - `K_from_BC_Heatmap/boundary_opt_standalone/configs/production_run.json`
- local MC executable built from:
   - `K_from_BC_Heatmap/boundary_opt_standalone/src/ising_tri_twisted_parallelogram.cc`

### Production configuration used
- reference mode: `continuum`
- grid:
   - `r1, r2 in {0.6, 0.8, 1.0, 1.2, 1.4}`
- test sizes:
   - `L = [8, 12, 16, 24]`
- reference sizes:
   - `L = [8, 12, 16, 24]`
- reference geometry:
   - quarter-twist family `(L, L, round(L/4), round(L/4))`
- production MC:
   - `n_traj = 500`
   - `n_skip = 10`
   - `n_therm = 500`
- beta_c finder:
   - coarse trajectories: `60`
   - fine trajectories: `120`
   - point counts: `11 / 5 / 5 / 5`
   - `max_shifts = 4`
- workers:
   - `4`

### Completion status
- jobs total: `104`
- jobs ok: `104`
- jobs err: `0`
- analysis cells: `25/25`
- wall time: `393.765 s` (`~6.56 min`)

Manifest:
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_20260512/manifest_boundary_opt.json`

### Key results
Global minimum from the standalone production run:
- `(r1, r2) = (1.4, 1.0)`
- `score_min = 0.06982707883`
- `z_rms = 0.2932348462`

This same point is also the best cell in the RMS z-score map.

### Top low-score cells
From `score_map.dat`, the leading cells are:
- `(1.4, 1.0) -> score = 0.0698271`
- `(0.8, 1.2) -> score = 0.0837523`
- `(0.8, 0.6) -> score = 0.1411193`
- `(1.2, 1.0) -> score = 0.1511940`
- `(1.2, 0.6) -> score = 0.1573089`

From `zscore_map.dat`, the leading low-z cells are:
- `(1.4, 1.0) -> z_rms = 0.293235`
- `(0.8, 0.6) -> z_rms = 0.438592`
- `(0.8, 1.2) -> z_rms = 0.446157`
- `(1.2, 1.0) -> z_rms = 0.470040`
- `(1.2, 0.8) -> z_rms = 0.499231`

### Basin shape
The score landscape is not flat.

Qualitative features:
- strong penalty bands are visible near `(0.6, 1.0)`, `(0.6, 1.2)`, `(1.0, 1.2)`, and `(1.4, 1.2-1.4)`;
- the best basin is anisotropic and favors larger `r1` with `r2` near `1.0`;
- a secondary low-score competitor exists near `(0.8, 1.2)` but remains above the winner in both score and RMS z-score.

---

## 2026-05-12 — Higher-stat standalone production rerun completed

### Run context
Reran the same standalone production grid at substantially higher statistics to test whether the apparent minimum near `(1.4, 1.0)` was stable or a moderate-stat fluctuation.

Runner and output roots:
- runner:
   - `K_from_BC_Heatmap/boundary_opt_standalone/run_boundary_opt.py`
- config:
   - `K_from_BC_Heatmap/boundary_opt_standalone/configs/production_run_high_stats.json`
- manifest:
   - `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/manifest_boundary_opt.json`

### What changed relative to the first production run
The grid, sizes, reference mode, fit mode, and channel set were held fixed.

The statistics were increased by a factor of about four:
- production MC:
   - `n_traj: 500 -> 2000`
   - `n_therm: 500 -> 1000`
- beta finder:
   - coarse trajectories: `60 -> 240`
   - fine trajectories: `120 -> 480`

### Completion status
- jobs total: `104`
- jobs ok: `104`
- jobs err: `0`
- analysis cells: `25/25`
- wall time: `460.052 s` (`~7.67 min`)

### Key result
The global minimum did not move.

Best cell in the higher-stat rerun:
- `(r1, r2) = (1.4, 1.0)`
- `score_min = 0.009216393701`
- `z_rms = 0.2657631678`

Comparison to the first production run at the same grid point:
- score improved from `0.06982707883` to `0.009216393701`
- RMS z-score improved from `0.2932348462` to `0.2657631678`

Interpretation:
- the earlier winner at `(1.4, 1.0)` is reinforced rather than overturned;
- the minimum is not just surviving the higher-stat rerun, it is becoming more sharply preferred in score;
- the best-fit basin still favors anisotropy with `r1 > r2` and `r2` close to `1.0`.

### Leading cells in the higher-stat rerun
From `score_map.dat`, the lowest-score cells are:
- `(1.4, 1.0) -> score = 0.00921639`
- `(1.4, 0.8) -> score = 0.0411699`
- `(1.2, 0.8) -> score = 0.0422568`
- `(0.8, 0.8) -> score = 0.0563991`
- `(0.6, 1.0) -> score = 0.0650945`

From `zscore_map.dat`, the lowest RMS-z cells are:
- `(1.4, 1.0) -> z_rms = 0.265763`
- `(1.4, 0.8) -> z_rms = 0.456523`
- `(1.2, 0.8) -> z_rms = 0.480068`
- `(0.8, 0.8) -> z_rms = 0.518612`
- `(0.6, 1.0) -> z_rms = 0.594942`

### Takeaway
At the current coarse `5 x 5` grid and with the higher-stat rerun completed, the standalone workflow now gives a consistent answer:
- the preferred basin remains centered on the cell `(r1, r2) = (1.4, 1.0)`;
- the main moderate-stat competitor near `(0.8, 1.2)` drops out of the top tier in the higher-stat rerun;
- the next rational production step is no longer “rerun the same grid,” but a local refinement around the surviving basin, especially with `r1` above `1.2` and `r2` near `0.9-1.1`.

### Output artifacts
Key files from this run:
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/dat/score_map.dat`
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/dat/zscore_map.dat`
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/dat/score_minimum.dat`
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/score_heatmap.png`
- `K_from_BC_Heatmap/boundary_opt_standalone/results/production_run_high_stats_20260512/score_and_zscore_heatmaps.png`

### Interpretation
Within this higher-stat standalone production configuration, the coarse-grid winner `(1.4, 1.0)` is clearly stabilized. The earlier moderate-stat competitor `(0.8, 1.2)` no longer sits near the minimum, so the remaining uncertainty is not whether the basin exists, but how to refine its location inside the surviving anisotropic neighborhood.

## 2026-05-08 — Continuum Boundary Optimization Smoke Run and Production Plan

### Context
Implemented and executed an end-to-end pipeline for continuum-calibrated boundary-correlator optimization over `(r1, r2)`:
- Script: `run_boundary_continuum_opt.py`
- Output tag: `boundary_opt_smoke`
- Output root: `results/boundary_opt_smoke/`

The pipeline performs:
1. MC precompute for test and reference lattices on size ladders.
2. Boundary-direction sampling at `t_k = k/8`, `k=0..7`, for all 3 cycles.
3. Continuum extrapolation per channel using quadratic fit in `1/L`.
4. Score construction:
   `S(r1,r2) = sum_{c=0..2} sum_{k=0..7} [G_test_inf(c,k)-G_ref_inf(c,k)]^2`
   (uniform weighting in this run).
5. Export of reusable `.dat` files.

---

### What was run
Smoke command (reduced scope for turnaround):
- `r-grid`: `r1,r2 in {0.9, 1.1}` (2x2 points)
- test sizes: `L = [8,16,24]`
- reference mode: `continuum`
- reference sizes: `L = [8,16,24]`
- reference geometries:
  - `L=8  -> (Lx,Ly,Tx,Ty)=(8,8,2,2)`
  - `L=16 -> (16,16,4,4)`
  - `L=24 -> (24,24,6,6)`
- `n_traj=300`, `n_traj_coarse=80`, `n_traj_fine=120`
- workers: `2`

Run status:
- jobs: `24/24 ok`
- errors: `0`
- precompute wall time: `~103.9 s`

---

### Key results
From `results/boundary_opt_smoke/dat/score_map.dat`:
- `(0.9, 0.9) -> 8.726826633`
- `(0.9, 1.1) -> 4.487569639`
- `(1.1, 0.9) -> 0.6815159369`  **minimum**
- `(1.1, 1.1) -> 3.39490116`

Minimum from `results/boundary_opt_smoke/dat/score_minimum.dat`:
- `r1_min = 1.1`, `r2_min = 0.9`, `score_min = 0.6815159369`

Data completeness:
- `n_channels_used = 24` at each point (3 cycles x 8 fractions), confirming full channel coverage.

Continuum fit usage:
- `n_used = 3` points per channel in this smoke run (expected, given 3 sizes).

---

### Interpretation
1. **Pipeline correctness is established** at smoke scale:
   - full precompute-compute-export loop ran successfully,
   - score map and minimum were produced,
   - all requested `.dat` products were emitted.

2. **Physics conclusion is preliminary**:
   - the run is intentionally low-stat and coarse-grid,
   - only 3 sizes were used for quadratic-in-`1/L` fits,
   - therefore uncertainty on channel intercepts is broad.

3. **Directionality signal exists**:
   - the minimum is not flat across the 4 points,
   - suggests structure is already emerging and worth refining with production settings.

---

### Reusable data generated
All reusable text datasets are under:
- `results/boundary_opt_smoke/dat/`

Files:
- `raw_test_channels.dat`
- `raw_ref_channels.dat`
- `continuum_test.dat`
- `continuum_ref.dat`
- `score_map.dat`
- `score_minimum.dat`

Metadata:
- `results/boundary_opt_smoke/manifest_boundary_opt.json`

Visualization:
- `results/boundary_opt_smoke/score_heatmap.png`

---

### Production plan (synthesized)

#### Phase 1 — Baseline production grid
- Test size ladder: `L_test = [8,16,24,32,40]`.
- Keep reference in `continuum` mode.
- Use matching reference size ladder initially (`[8,16,24,32,40]`) with user-defined twists.
- Increase statistics substantially (`n_traj >= 20000` minimum; higher if surface remains noisy).
- Start with moderate `(r1,r2)` mesh to map basin topology.

#### Phase 2 — Basin refinement
- Identify low-score basin from Phase 1.
- Re-run local subgrid around minimum with finer `r-step`.
- Increase `n_traj` only in local region for compute efficiency.

#### Phase 3 — Reference strategy cross-check
- Re-run selected points with `reference_mode=large` using a very large twisted reference lattice.
- Compare minimum location and basin shape against continuum-reference result.
- Use agreement as confidence check for final reported optimum.

#### Phase 4 — Final selection and handoff
- Report final `(r1,r2)` with neighboring score contrast.
- Archive manifest + all `.dat` files as authoritative reusable products.
- Preserve exact CLI command lines used for reproducibility.

---

### Operational notes for next run
- Keep output tags explicit (for example `boundary_opt_prod_v1`, `boundary_opt_prod_refined`).
- Do not overwrite smoke results; treat smoke as baseline benchmark.
- Continue using `.dat` as canonical exported data format.

---

## 2026-05-11 — Phase 1 Baseline Production Grid (in progress)

### Objective
Execute Phase 1 from the production plan:
- moderate `(r1,r2)` grid,
- higher statistics,
- continuum-reference strategy,
- expanded size ladder.

### Command launched
```bash
/workspaces/newQFE/.venv/bin/python K_from_BC_Heatmap/run_boundary_continuum_opt.py \
   --tag boundary_opt_prod_phase1 \
   --r-min 0.6 --r-max 1.4 --r-step 0.2 \
   --test-sizes 8 16 24 32 40 \
   --reference-mode continuum \
   --ref-sizes 8 16 24 32 40 \
   --ref-geom 8:8:8:2:2 \
   --ref-geom 16:16:16:4:4 \
   --ref-geom 24:24:24:6:6 \
   --ref-geom 32:32:32:8:8 \
   --ref-geom 40:40:40:10:10 \
   --n-traj 20000 --n-traj-coarse 3000 --n-traj-fine 6000 \
   --beta-seed 0.05 0.40 \
   --n-workers 4
```

### Planned workload
- r-grid points: `5 x 5 = 25`
- physics-distinct test jobs: `25 x 5 = 125`
- physics-distinct reference jobs (fixed geometry, isotropic): `5` (one per reference size)
- physics-distinct total: `130`

Implementation note:
- The current runner still schedules reference payload tasks with `(r1,r2)`-keyed filenames, so the precompute counter may report a larger implementation-level total (for example `250`).
- This does **not** mean reference physics is scanned over anisotropic couplings; the intended reference dataset is fixed-geometry isotropic and reused.

### Live status snapshot (latest poll)
- stage: `precompute`
- reported by runner: `total=250 cached=0 todo=250 workers=4` (implementation scheduler count)
- observed completed jobs in log: `33/250`, all `ok`, `0` errors observed so far
- current completed tranche: full test `L=8` grid completed (`25/25`), then test `L=16` in progress (`8/25` points complete)
- representative beta values (L=8 test):
   - `beta_c(r1=0.6,r2=0.6) = 0.32998`
   - `beta_c(r1=1.0,r2=1.0) = 0.23910`
   - `beta_c(r1=1.4,r2=1.4) = 0.18977`

### Output paths
- run root: `results/boundary_opt_prod_phase1/`
- intermediate caches: `results/boundary_opt_prod_phase1/grid/`
- analysis exports (`.dat`) will appear after precompute completes under:
   - `results/boundary_opt_prod_phase1/dat/`

### Interpretation (interim)
1. Run is healthy and advancing normally (no failures in the first 10% of jobs).
2. The expected monotone decrease of `beta_c` with increasing anisotropy ratio is visible in early `L=8` points.
3. No score-map physics conclusion is available yet; continuum channel fits and `score_map.dat` are produced only after full precompute completion.

---

## 2026-05-11 — Phase 1 Baseline Production Grid (completed analysis)

### Completion status
- Precompute and continuum analysis completed successfully for the full `5 x 5` grid.
- Final exports were written under `results/boundary_opt_prod_phase1/`, including:
   - `dat/score_map.dat`
   - `dat/score_minimum.dat`
   - `dat/zscore_map.dat`
   - `dat/fit_quality_map.dat`
   - `score_heatmap.png`
   - `score_and_zscore_heatmaps.png`
   - `plots/fss_with_fit_all_points/`
- The production score uses the nontrivial boundary channels only: `k=1..7` for each of the 3 cycles (`21` channels total). The `k=0` self-point was dropped from the scoring and FSS visualization because it is identically `1` and carries no discriminating continuum information.

### Key results
- Global minimum from `score_minimum.dat`:
   - `(r1, r2) = (1.0, 1.4)`
   - `score_min = 0.001914935309`
- The same point is also the minimum of the RMS test-vs-reference z-score map:
   - `z_rms(1.0, 1.4) = 0.5019656755`
- Best nearby competitors on the coarse grid are:
   - `(1.4, 1.2) -> score = 0.00233381682, z_rms = 0.5992526710`
   - `(0.6, 1.0) -> score = 0.002452757784, z_rms = 0.5617903094`
   - `(1.0, 1.0) -> score = 0.002581866268, z_rms = 0.5766450641`
   - `(0.8, 1.4) -> score = 0.003847987843, z_rms = 0.6888622808`
- The basin is clearly anisotropic rather than symmetric:
   - `(1.0, 1.4)` is the best point on the grid,
   - while the transposed point `(1.4, 1.0)` is poor: `score = 0.03228862015`, `z_rms = 2.023420327`.

### Fit diagnostics
- FSS-with-fit plots were generated for all `25` test points and the reference point, with:
   - weighted quadratic fits in `1/L`,
   - jackknife error bars on the measured channels,
   - continuum intercept markers with propagated fit uncertainty.
- The fit-quality map shows substantial variation across the grid.
   - clean extrapolation regions include `(1.0, 1.2)` with median reduced `chi2 ~ 0.048` and `(0.8, 1.0)` with median reduced `chi2 ~ 0.141`,
   - poor regions include `(1.2, 1.2)` with median reduced `chi2 ~ 13.15` and `(0.6, 1.4)` with median reduced `chi2 ~ 12.46`.
- The optimum `(1.0, 1.4)` is statistically favored by both score and z-score, but it is not the cleanest fit region in reduced-`chi2` terms:
   - median reduced `chi2(1.0, 1.4) ~ 4.03`.
- The correct reading is therefore: `(1.0, 1.4)` is the best coarse-grid candidate from the present production scan, not yet the final precision optimum.
- A further caution is that the top of the basin is only weakly resolved. The leading score values are very close:
   - `(1.0, 1.4) -> 0.001915`
   - `(1.4, 1.2) -> 0.002334`
   - `(0.6, 1.0) -> 0.002453`
   - `(1.0, 1.0) -> 0.002582`
- In particular, `(1.0, 1.0)` differs from `(1.0, 1.4)` by only `6.67e-4` in score, and the channel-by-channel continuum shifts between those two points are typically sub-`1 sigma`, with inter-point RMS `z ~ 0.68`. That makes it plausible that the ranking within the best few cells is noise-dominated.
- By contrast, the broader map is not purely noise: several cells have RMS test-vs-reference `z > 2`, so the data do still discriminate clearly against parts of parameter space far from the low-score basin.

### Synthesis blurb
The phase-1 production scan identifies a low-score basin whose best coarse-grid point is `(r1, r2) = (1.0, 1.4)`. This point minimizes both the continuum score and the channel-wise RMS z-score, so it is the correct provisional target for follow-up. At the same time, the top of the basin is only weakly resolved: `(1.0, 1.0)`, `(0.6, 1.0)`, and `(1.4, 1.2)` are all very close in score, and the channel-level differences between `(1.0, 1.0)` and `(1.0, 1.4)` are mostly below `1 sigma`. That suggests the fine ordering inside the best few cells is likely noise-dominated. The stronger conclusion is therefore not that `(1.0, 1.4)` has been decisively separated from `(1.0, 1.0)`, but that the current data disfavors large regions of the grid and points to an anisotropic candidate basin that should now be revisited with a higher-statistics, finer local scan.

### Next steps
1. Move `K_from_BC_Heatmap` into a plug-and-play production workflow so a high-statistics run can be launched overnight with one command rather than by hand stitching together analysis steps.
2. The plug-and-play version should include a single top-level config/CLI, resumable per-point caches, automatic export of `.dat` tables and diagnostic plots, and an unambiguous run manifest for reproducibility.
3. Use that overnight workflow to run a refined local scan centered on `(1.0, 1.4)`, for example with `r1` near `1.0` and `r2` in the `1.2-1.6` band, at substantially higher statistics than the current `n_traj = 20000` baseline.
4. Re-rank the refined basin using score, RMS z-score, and fit-quality together before promoting a final recommended boundary condition.
