# Standard Raw Two-Point Campaign

This subdirectory is the raw-data staging area for a simple, resumable source-built simulator campaign. The goal is to produce only large twisted references plus untwisted parameter-space sweeps and save the resulting all-to-all two-point data for later organization and analysis.

Nothing here runs the continuum-manifold comparison pipeline. The runner goes straight to the source-built simulator through `K_from_continuum/workflow_common.py` and copies only `two_point_all_to_all.dat` plus a small metadata sidecar into this directory tree.

## Scope

The campaign requested here is:

- `iso111` twisted reference: one large isotropic lattice with volume above `20000`
- `acute456` twisted reference: one large story-style acute `4-5-6` lattice with volume above `20000`
- `iso111` untwisted sweep: a `5 x 5` grid in `(r1, r2)` around `(1, 1)`
- `acute456` untwisted sweep: a `5 x 5` grid in `(r1, r2)` around the story acute `4-5-6` center `(4.702782819756, 7.353910143333)`
- for each untwisted point: square lattices `8x8`, `16x16`, `24x24`, `32x32`, `48x48`, `64x64`
- Monte Carlo production settings: `n_traj = 20000`, `n_skip = 10`, `n_therm = 2000`
- critical tuning: exact triangular Ising sinh rule at each coupling point

The explicit multiplier stencil is the one from the request, not a symmetric `+/-10%` shorthand:

- `r1` multipliers: `0.8, 0.9, 1.0, 1.1, 1.2`
- `r2` multipliers: `0.8, 0.9, 1.0, 1.1, 1.2`

That gives:

- `2` twisted reference jobs
- `25 x 6 = 150` untwisted jobs per geometry
- `302` total jobs in the full campaign

## Geometry Choices

These are the fixed choices encoded in `run_standard_campaign.py`.

### `iso111`

- twisted reference lattice: `(144, 144, 0, 0)`
- volume: `20736`
- twisted couplings: `(k1, k2, k3) = (1, 1, 1)`
- untwisted center: `(r1*, r2*) = (1, 1)`

### `acute456`

- twisted reference lattice: `(144, 144, 72, 24)`
- volume: `22464`
- twisted couplings: `(k1, k2, k3) = (1, 1, 1)`
- untwisted center: `(r1*, r2*) = (4.702782819756, 7.353910143333)`

This `acute456` center is the story-compatible realized untwisted match used around the existing acute `4-5-6` figures, rather than the exact `(13,16,-3,3)` continuum-side solution.

## Critical Beta Rule

Every run here uses the exact anisotropic triangular Ising critical point from the sinh rule,

`sinh(2 beta k1) sinh(2 beta k2) + sinh(2 beta k2) sinh(2 beta k3) + sinh(2 beta k3) sinh(2 beta k1) = 1`.

The runner computes `beta_c` from this equation for each untwisted candidate and for the isotropic twisted references.

## Output Layout

The final copied data live under:

- `responsible_method_tests/standard/data/iso111/...`
- `responsible_method_tests/standard/data/acute456/...`

Each finished job writes:

- `two_point_all_to_all.dat`
- `two_point_all_to_all.meta.json`

The `.meta.json` records the lattice, couplings, center point, multiplier offsets, exact beta, MC settings, and the original scratch-file path produced by the simulator.

The campaign-level bookkeeping files are:

- `campaign_plan.json`: full static job list and campaign settings
- `job_log.jsonl`: append-only per-job completion or failure log
- `status_summary.json`: rolling counts of completed, skipped, failed, and pending jobs
- `_mc_scratch/`: raw simulator output directories used during production

## Current Run State

As of the May 25, 2026 session that created this directory:

- all `300` untwisted jobs completed successfully
- the large `iso111` twisted reference completed successfully
- the large `acute456` twisted reference was the only stubborn case

Two separate issues occurred for the large twisted references during the initial batch:

- the original workflow_common/mc_engine path imposed a hard `600 s` subprocess timeout, which was too short for `144x144` twisted references at `n_traj=20000`
- after removing that wrapper timeout in `run_standard_campaign.py`, the `acute456` twisted reference still returned a non-zero exit code in one scripted rerun even though a tiny-stat direct probe on the same geometry succeeded

Because of that, the final acute `4-5-6` twisted reference was relaunched outside the Python batch runner, in its own detached session, to avoid any remaining wrapper or shell-lifetime effects:

```bash
setsid -f K_from_continuum/bin/ising_tri_twisted_parallelogram \
	--L_x 144 --L_y 144 --T_x 72 --T_y 24 \
		--k1 1.000000000000 --k2 1.000000000000 --k3 1.000000000000 \
			--beta 0.274653072167 \
				--n_traj 20000 --n_skip 10 --n_therm 2000 \
					--seed 246813579 \
						--data_dir responsible_method_tests/standard/_manual/acute456_twisted_full \
							> responsible_method_tests/standard/_manual/acute456_twisted_full/stdout.log \
								2> responsible_method_tests/standard/_manual/acute456_twisted_full/stderr.log \
									< /dev/null
									```

									The detached-session version matters here: a plain background job owned by the interactive shell was still vulnerable to shell cleanup, while the `setsid -f` launch leaves the simulator re-parented to init and still running after the shell exits.

									The tiny-stat direct geometry probe that confirmed the lattice itself is valid used:

									```bash
									K_from_continuum/bin/ising_tri_twisted_parallelogram \
										--L_x 144 --L_y 144 --T_x 72 --T_y 24 \
											--k1 1.000000000000 --k2 1.000000000000 --k3 1.000000000000 \
												--beta 0.274653072167 \
													--n_traj 1 --n_skip 1 --n_therm 10 \
														--seed 123 \
															--data_dir responsible_method_tests/standard/_probes/acute456_twisted_probe
															```

															If the manual direct run finishes successfully, the remaining housekeeping step is simple:

															- import it with:

															```bash
															python3 responsible_method_tests/standard/run_standard_campaign.py \
																--import-manual-job acute456__twisted__reference \
																	--manual-source-dir responsible_method_tests/standard/_manual/acute456_twisted_full
																	```

																	That command copies the finished `two_point_all_to_all.dat` into the standard data tree, writes the matching metadata, appends a `completed` row for `acute456__twisted__reference` to `job_log.jsonl`, and refreshes `status_summary.json`.

																	In the current live setup, that import step is also armed automatically by a detached close-write watcher rooted at the same manual directory:

																	- `import_watch.pid` records the watcher shell pid
																	- `import_watch.log` will capture the eventual import action once `two_point_all_to_all.dat` is closed after writing

																	## How To Launch Or Resume

																	From the repo root:

																	```bash
																	source /workspaces/newQFE/.venv/bin/activate
																	python responsible_method_tests/standard/run_standard_campaign.py --workers 2
																	```

																	The runner is resumable. If both the copied `.dat` file and its `.meta.json` already exist for a job, that job is skipped automatically.

																	Useful restricted restarts:

																	```bash
																	python responsible_method_tests/standard/run_standard_campaign.py --geometry iso111 --workers 2
																	python responsible_method_tests/standard/run_standard_campaign.py --geometry acute456 --method twisted --workers 1
																	python responsible_method_tests/standard/run_standard_campaign.py --method untwisted --max-jobs 10 --workers 2
																	python responsible_method_tests/standard/run_standard_campaign.py --dry-run
																	```

																	## Notes For Later Analysis

																	- This campaign intentionally does not generate manifold interpolation, continuum fits, or comparison plots.
																	- The acute `4-5-6` center here is the story-compatible realized match. If we later want the exact `4-5-6` continuum center `(5.065231..., 7.742930...)`, that should be a second standard campaign rather than silently mutating this one.
																	- The twisted references are both single large lattices above the `V > 20000` threshold and are meant as raw anchor data for later holdout or direct-comparison work.
																	- The untwisted ladders use the exact sizes requested by hand: `8, 16, 24, 32, 48, 64`.
																	- If a run fails, the traceback is captured in `job_log.jsonl`; the failed job can then be rerun with the same command because successful jobs are skipped.