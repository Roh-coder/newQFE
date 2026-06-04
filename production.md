# Production Plan

## Scope Assumed Here

- Dense coupling grid: $r_1, r_2 \in [0.3, 3.0]$ with step $0.01$.
- Inclusive grid count: $((3.0 - 0.3) / 0.01) + 1 = 271$ points per axis, so
  $271^2 = 73{,}441$ total coupling points.
- High-stat benchmark definition used for estimates: `n_traj=100000`, `n_therm=10000`, `n_skip=10`.
- Runtime estimate below is first anchored to a live `32x32` pilot, because the request specified a dense $(r_1,r_2)$ grid but did not specify a size ladder.

## What Was Measured

Two live pilot runs were executed on `2026-06-04` with the current simulator:

1. `32x32`, isotropic `(r1,r2)=(1,1)`, `n_traj=20000`, `n_therm=2000`, `n_skip=10`, `--dump_traces 1`
   - duration: `7.355173 s`

2. `32x32`, isotropic `(r1,r2)=(1,1)`, `n_traj=100000`, `n_therm=10000`, `n_skip=10`, `--dump_traces 1`
   - duration: `35.944629 s`
   - `two_point_all_to_all.dat`: `104812 B`
   - `traces_*.dat`: `580193 B`

Benchmark outputs are under:

- `/workspaces/newQFE/_benchmarks/reweight_32x32_20k/`
- `/workspaces/newQFE/_benchmarks/reweight_32x32_hi100k/`

The measured scaling from the live pilots is close to linear:

- $35.944629 / 7.355173 = 4.887$

That is the scaling factor used below when lifting existing `20k` ladder timings to the `100k` case.

## Runtime Estimate

### Case A: single `32x32` lattice only

Using the live `32x32` high-stat pilot:

- per coupling point: `35.944629 s`
- total serial wall time:
  $$73{,}441 \times 35.944629\,\text{s} = 30.55\,\text{days}$$
- total wall time on this 2-core dev container with 2 independent workers:
  $$15.28\,\text{days}$$

For comparison, the same dense `32x32` grid at the older `20k/2k/10` setting is still large:

- serial: `6.25 days`
- 2 workers: `3.13 days`

### Case B: the usual 7-size ladder `(4, 6, 8, 12, 16, 24, 32)`

Using the existing acute456 `20k` metadata at the center point for the 7-size ladder,
the summed per-coupling runtime is about `32.35 s` at `20k`. Applying the measured live
scaling factor `4.887` gives:

- estimated per coupling point at `100k`: `158.08 s`
- total serial wall time:
  $$73{,}441 \times 158.08\,\text{s} = 134.37\,\text{days}$$
- total wall time on 2 workers:
  $$67.19\,\text{days}$$

This excludes rescoring, plotting, retries, and any failed-job handling.

## Current Histogram-Reweighting Status

The current simulator already has a Ferrenberg-Swendsen trace dump path:

- `K_from_continuum/src/ising_tri_twisted_parallelogram.cc`
- opt-in CLI: `--dump_traces 1`

What it writes today is exactly one row per measured trajectory with three values:

- `E_per_site`
- `abs_m`
- `m^2`

That is enough for susceptibility-style reweighting,
$$\chi(\beta) = \langle m^2 \rangle_\beta - \langle |m| \rangle_\beta^2,$$
and it matches the expectations of the existing reweight implementation used elsewhere in the repo.

It is **not** enough to reweight correlation functions.

For correlation-function reweighting, the production run needs per-measured-trajectory data for the observables that will be reweighted. At minimum that means:

- `E_per_site` for the Boltzmann weight shift
- `abs_m` and `m^2` if the connected subtraction continues to use the current convention
- `corr(d)` for each tracked displacement or tracked observable channel
- stable block structure or one-row-per-measurement data so jackknife/block errors can be reconstructed offline

There is another practical issue: `K_from_continuum/lib/mc_engine.py` imports a local `reweight` module, but `K_from_continuum/lib/` in the current workspace only contains `cost.py` and `mc_engine.py`. A known-good `reweight.py` exists elsewhere in the repo, for example:

- `/workspaces/newQFE/K_from_Optimizer_Production/reweight.py`

So the current `K_from_continuum` tree is not yet self-contained for histogram reweighting work.

## Storage Estimate

### With the current trace format only

Using the live `32x32` `100k` pilot sizes:

- `two_point_all_to_all.dat` per point: `104812 B`
- `traces_*.dat` per point: `580193 B`

For the full dense grid of `73,441` points:

- all-to-all outputs: `7.70 GB`
- FS traces: `42.61 GB`
- subtotal: about `50.3 GB` before metadata, logs, and restarts

### If correlation-function traces are added

At `100k` trajectories with `n_skip=10`, there are about `10000` measured rows per point.

If correlation reweighting stores binary doubles for

- `E_per_site`
- `abs_m`
- `m^2`
- plus `k` tracked correlator observables

then the dense-grid storage is approximately:

- `k=5`: `0.047 TB`
- `k=25`: `0.165 TB`
- `k=100`: `0.605 TB`
- `k=271`: `1.61 TB`
- `k=1024` (full `32x32` all-to-all): `6.03 TB`

Those numbers are for binary payloads. Text output would be much larger and is not a sane format for the correlation-reweighting production path.

## Recommendation

Do **not** launch the full `0.01` dense grid yet.

The two reasons are structural, not cosmetic:

1. The wall time is already multiple weeks even for a single `32x32` lattice on this machine.
2. The current trace format only supports susceptibility reweighting, not correlation-function reweighting.

## Concrete Production Plan

### Phase 0: Make the reweighting path self-contained

1. Vendor a working `reweight.py` into `K_from_continuum/lib/`.
2. Add a small test that confirms the current 3-column trace format parses and reproduces the existing susceptibility reweighting path.

### Phase 1: Add correlation-trace support

1. Extend `K_from_continuum/src/ising_tri_twisted_parallelogram.cc` with a new optional correlation trace mode.
2. Default should remain off.
3. The new mode should support a **selected observable set**, not full all-to-all by default.
4. Preferred output format: binary or containerized binary (`npz`, `npy`, HDF5, or a compact custom format), not text.

Suggested trace payload per measurement row:

- `E_per_site`
- `abs_m`
- `m^2`
- `corr(obs_1)`
- `corr(obs_2)`
- ...

Where `obs_i` should be limited to the actual score observables used downstream. For the current boundary-anchor workflow, that likely means a very small set such as:

- `v` quarter-point
- `u` quarter-point
- `w` quarter-point
- anchor point

Ratios can then be reconstructed offline from reweighted means rather than stored directly.

### Phase 2: Validate the new trace path on a small pilot

1. Run an `11 x 11` pilot grid, for example `r1,r2 \in [0.5,1.5]` with step `0.1`.
2. Use `32x32`, `100k/10k/10`, and the new correlation trace mode.
3. Measure:
   - runtime per point
   - trace size per point
   - effective reweighting window in $\beta$
   - whether reweighted correlators are stable enough to replace reruns nearby in coupling space

### Phase 3: Decide between brute force and adaptive refinement

After the pilot, choose one of two paths:

1. If the score landscape is broad, do an adaptive scan or coarse-to-fine refinement rather than a full brute-force `0.01` grid.
2. Only if the score develops sharp localized structure that justifies it should the full `0.01` production grid be launched.

### Phase 4: Only then schedule production

If the final choice is still the full brute-force dense grid:

1. Move it off this 2-core dev container onto a larger machine or batch scheduler.
2. Keep the production observable set minimal.
3. Use resumable manifests and chunk the grid into tiles.
4. Separate raw trace storage from derived score products.

## Bottom Line

For a dense $0.01$ grid over $[0.3,3.0]^2$:

- single `32x32` high-stat run: about `15.3 days` on this machine with 2 workers
- 7-size ladder to `32x32`: about `67.2 days` on this machine with 2 workers

And the current codebase does **not** yet record enough per-trajectory data to do histogram reweighting of the correlation functions themselves.

So the correct immediate move is:

1. make the reweighting stack self-contained in `K_from_continuum`
2. add compact selected-observable correlation traces
3. run a small high-stat pilot grid
4. only then decide whether a full `0.01` brute-force production is justified

## Prototype Continuation Plan (2026-06-04 Update)

The sections above were written before the first local correlation-reweighting prototypes were run.
The current working-tree status is now more specific:

- selected-correlator nearby-coupling reweighting works at `32x32`, `100k/10k/10`
- the same idea extended over a midpoint FSS ladder remains good through `L <= 24`, but starts to miss at `L = 32`
- nearest-donor interpolation of the **full geometry-match proxy score surface** is not yet reliable, even when `N_eff / N` stays extremely high

So the next work should be framed as a **prototype-debugging program**, not as production scanning.

### Working Baseline To Reuse

Use the following scripts and artifact trees as the starting point for all follow-up work:

- nearby single-correlator control:
  `responsible_method_tests/scripts/test_correlator_reweight_nearby.py`
- nearby correlator plots:
  `responsible_method_tests/scripts/plot_correlator_reweight_nearby.py`
- geometry-match interpolation pilot:
  `responsible_method_tests/scripts/test_geometry_match_grid_interpolation.py`
- nearby single-correlator result tree:
  `responsible_method_tests/results/correlator_reweight_nearby_iso111_L32_hi100k_20260604/`
- nearby FSS ladder result tree:
  `responsible_method_tests/results/correlator_reweight_fss_iso111_sizes4_6_8_12_16_24_32_mid_hi100k_20260604/`
- geometry-match interpolation result tree:
  `responsible_method_tests/results/geometry_match_reweight_interp_iso111_grid5x5_donor3x3_sizes8_12_16_24_target32_20260604/`

Current headline numbers to keep in mind:

- nearby 3-point correlator control: all off-diagonal reweight checks were within `|z| < 1`
- nearby midpoint FSS ladder: worst `L=32` midpoint miss reached about `|z| = 2.46`
- geometry-match interpolation pilot: direct best `(0.99, 1.00)` vs reweighted best `(0.995, 1.00)`, `delta_z2_sum_rmse = 4.12`, worst `|delta_z2_sum| = 13.68`, while `min N_eff / N > 0.99985`

### Main Prototype Question

The immediate question is no longer "can observables be reweighted at all?"

That answer is already yes, locally.

The real question now is:

> At which stage does the geometry-match interpolation fail?

The likely failure points are:

1. raw reweighted connected observables
2. connected-ratio construction
3. finite-size fitting / extrapolated intercepts
4. final multi-panel score aggregation
5. the nearest-donor approximation itself, even if the raw observable reweight is fine

### Phase A: Diagnose The Existing Failure With No New MC

Do not spend more Monte Carlo first. Reuse the cached `5x5` / `3x3` iso111 interpolation data.

Add a diagnostic layer that dumps, for each fine-grid point and each fitted size:

- direct connected values for all stored observables
- reweighted connected values from the donor
- direct ratios for every score panel
- reweighted ratios for every score panel
- fitted continuum intercepts and fit errors per panel
- final panel `z` values and total `z2_sum`

Focus first on four sentinel cells:

- exact donor/self-check: `(0.99, 1.00)`
- argmin shift case: `(0.995, 1.00)`
- catastrophic mismatch case: `(0.99, 1.005)` from donor `(0.99, 1.00)`
- center baseline: `(1.00, 1.00)`

Deliverables for this phase:

- one machine-readable diagnostics file per sentinel case
- one compact comparison plot showing where the direct and reweighted pipelines first separate

Exit condition for Phase A:

- we can name the first stage where the large mismatch appears, rather than only observing the final `z2_sum` failure

### Phase B: Minimal High-Stat Subset To Separate Noise From Structural Failure

Once Phase A identifies the first unstable stage, rerun only the smallest subset needed to answer whether the current failure is mainly statistical or structural.

Recommended subset:

- donors and targets needed for the four sentinel cells above
- same size ladder `L = 8, 12, 16, 24`, plus `L = 32` holdout
- high-stat setting `100k/10k/10`

Interpretation rule:

- if the bad sentinel cells collapse toward the direct result at `100k`, the present failure is mostly noise amplification through the nonlinear pipeline
- if the same cells remain badly wrong at `100k`, the nearest-donor model itself is the problem

Do not scale this to the full `5x5` grid until that distinction is made.

### Phase C: Prototype Better Interpolants Only After A/B

If Phase B shows that nearest-donor is structurally too crude, test better interpolants on the same cached data before running any new broad grid.

Recommended order:

1. multidonor observable blend using a local `2x2` or `3x3` donor stencil in `(r1, r2)`
2. local linear interpolation of each **raw connected observable** before ratio formation
3. compare `interpolate observables -> build ratios -> fit` against `fit donor families first -> interpolate fitted intercepts`

The default should be to interpolate the most linear object available as early in the pipeline as possible. That means raw connected observables are a better starting point than final `z2_sum`.

### Provisional Acceptance Criteria

Do not promote any interpolant just because `N_eff / N` is large. Use explicit surface-level checks.

A candidate interpolation method should pass the iso111 control only if it does all of the following:

1. recovers the direct-grid argmin on the `5x5` iso111 control
2. achieves `delta_z2_sum_rmse < 1.0` over that grid
3. keeps `max |delta_z2_sum| < 2.0` on all cells with direct `z2_sum < 5`
4. preserves `N_eff / N > 0.99` as a sanity check, but not as the primary success metric

These thresholds are prototype gates, not final physics criteria. They are there to prevent another false sense of safety from overlap numbers alone.

### Phase D: Transfer Before Any Dense Scan

If and only if one interpolation method passes the iso111 local control, transfer it to one anisotropic local box before discussing any dense production grid.

Use exactly one of the following as the next transfer case:

- acute456 local box around the current standard-search candidate region
- quarter-twist anchored local box if that workflow is easier to instrument with the same selected observables

The transfer goal is not to find the production optimum. It is only to check whether the method survives a non-iso111 local landscape.

### Hard Stop Rule

Do not launch the full dense `0.01` geometry grid until all of the following are true:

1. the failure stage in the current iso111 pilot is identified
2. a replacement interpolant passes the iso111 control gate above
3. that interpolant also survives one anisotropic transfer test
4. the reweighting runtime / storage path is self-contained inside `K_from_continuum`

Until then, the correct scale of work is targeted prototype debugging, not production sampling.