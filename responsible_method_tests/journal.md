# Journal

## 2026-05-27: Acute456 Full 5x5 Anchor-Ratio Landscape

### Goal

Test the stronger version of the anchor-ratio hypothesis for the standard acute `4-5-6` campaign:

- run the full `5 x 5` acute456 candidate grid with `normalization_mode = anchor_ratio`
- use the existing orbit-reduced pointwise target score as the ranking metric
- check whether the target candidate `(r1, r2) = (4.702783, 7.353910)` is actually the global minimum of that landscape
- decide whether the origin-inclusive full-point aggregate should be part of the main score or only treated as a diagnostic

### Commands Used

Full acute456 anchor-ratio landscape:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py \
	--family acute456 \
	--normalization-mode anchor_ratio
```

Full-point anchor-ratio checks for the target and two challengers:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_all_base_points_individual_fits_anchor_ratio.tsv

python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--untwisted-dir responsible_method_tests/standard/data/acute456/untwisted/r1_5p643339__r2_8p089301 \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_r1_5p643339_r2_8p089301_all_base_points_anchor_ratio.tsv

python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--untwisted-dir responsible_method_tests/standard/data/acute456/untwisted/r1_5p643339__r2_8p824692 \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_r1_5p643339_r2_8p824692_all_base_points_anchor_ratio.tsv
```

### Artifacts Written

- landscape figure:
	- `responsible_method_tests/standard/standard_pointwise_score_landscapes_anchor_ratio.png`
- acute456 landscape table:
	- `responsible_method_tests/standard/acute456_pointwise_score_landscape_anchor_ratio.tsv`
- target full-point anchor-ratio table:
	- `responsible_method_tests/standard/acute456_all_base_points_individual_fits_anchor_ratio.tsv`
- challenger full-point anchor-ratio tables:
	- `responsible_method_tests/standard/acute456_r1_5p643339_r2_8p089301_all_base_points_anchor_ratio.tsv`
	- `responsible_method_tests/standard/acute456_r1_5p643339_r2_8p824692_all_base_points_anchor_ratio.tsv`

### Main Quantitative Result

The full acute456 `5 x 5` anchor-ratio landscape does **not** place the target at the global minimum.

Top five orbit-reduced anchor-ratio candidates from `acute456_pointwise_score_landscape_anchor_ratio.tsv`:

| rank | `(r1, r2)` | orbit RMS z | jackknife sigma |
| --- | --- | ---: | ---: |
| 1 | `(3.762226, 6.618519)` | `0.1603` | `0.0232` |
| 2 | `(4.702783, 7.353910)` | `0.1918` | `0.0253` |
| 3 | `(5.173061, 8.089301)` | `0.3374` | `0.0437` |
| 4 | `(5.643339, 8.089301)` | `0.3376` | `0.0309` |
| 5 | `(5.643339, 8.824692)` | `0.3548` | `0.0554` |

So the target is **rank 2 / 25**, not rank 1.

The result is therefore mixed but informative:

- anchor-ratio scoring does beat the two explicit challengers checked earlier
- but it still prefers a lower-coupling candidate `(3.762226, 6.618519)` over the target

### Full-Point Versus Orbit-Reduced Behavior

For the target candidate `(4.702783, 7.353910)`:

- full-point aggregate RMS z = `1.3129`
- orbit-reduced RMS z = `0.1918`
- worst point is the origin with `z = -10.3191`

For the two explicit off-target challengers that were checked with the full-point anchor-ratio fit script:

- `(5.643339, 8.089301)` gave full-point aggregate RMS z = `1.3536` and orbit-reduced RMS z = `0.3376`
- `(5.643339, 8.824692)` gave full-point aggregate RMS z = `1.4094` and orbit-reduced RMS z = `0.3548`

So the same origin-dominance pattern persists away from the target: once the origin is included, the full-point aggregate is inflated far above the orbit-reduced score and becomes a poor primary ranking statistic.

The operational lesson is the same as in the earlier pointwise tables:

- the origin carries a qualitatively different normalization signal
- once it is mixed into the full-point aggregate, that aggregate no longer cleanly tracks the geometry-realization question
- the orbit-reduced score is the more stable quantity for ranking candidates

### Score-Definition Decision

For the standard anchor-ratio workflow, keep the **orbit-reduced score as the primary geometry score** and treat the full-point aggregate as a **diagnostic only**.

Concretely:

- do **not** use the origin-inclusive full-point aggregate RMS z as the main ranking metric
- do use the orbit-reduced RMS z table/landscape as the primary selection surface
- keep the full-point aggregate, worst-point list, and origin residuals as debugging outputs that explain why a candidate looks bad

I did **not** change the landscape code for this, because `plot_standard_pointwise_score_landscapes.py` already ranks candidates by the orbit-reduced score. The new acute456 run supports keeping that design rather than switching the primary score to the origin-inclusive aggregate.

### Interpretation

- The anchor-ratio construction is useful: it separates the target from obvious off-target challengers much better than the raw all-point aggregate does.
- But it is not yet a complete realization score: even on the full acute456 `5 x 5` grid, the target is only the **second-best** candidate.
- The current evidence therefore supports the weaker statement:
	- anchor-ratio FSS contains real geometry information
	- orbit-reduced anchor-ratio is a reasonable primary score
	- but anchor-ratio alone does not yet recover the exact target as the unique minimum of the acute456 grid

### Current Takeaway

For standard acute456, the best-supported workflow is:

- use the orbit-reduced anchor-ratio landscape to rank candidates
- treat the full-point aggregate as a diagnostic for origin pathologies and other local failures
- do not claim yet that anchor-ratio FSS alone identifies the target geometry uniquely

### Follow-Up: L8-Ratio Cross-Check

I followed up with the second normalization choice `normalization_mode = l8_ratio`.

That cross-check shows that the anchor-ratio winner is **not** stable under a normalization change:

- `l8_ratio` landscape winner is `(5.643339, 8.089301)` with orbit RMS z = `0.2813 +/- 0.0275`
- target `(4.702783, 7.353910)` is rank `7 / 25` with orbit RMS z = `0.4382 +/- 0.0197`
- former anchor-ratio winner `(3.762226, 6.618519)` falls to rank `10 / 25` with orbit RMS z = `0.5430 +/- 0.0376`

So the specific anchor-ratio minimum at `(3.762226, 6.618519)` is normalization-dependent rather than a robust optimum.

### Follow-Up: Origin-Dropped Full-Point Summary

To make the full-point reporting less ambiguous, `fit_standard_acute456_all_base_points.py` now prints both:

- `aggregate_score` (origin included if present in the summary rows)
- `origin_dropped_aggregate_score` (same point set with the origin removed)

Under `l8_ratio`, dropping the origin barely changes the full-point target score:

- target `(4.702783, 7.353910)`: `0.4285 -> 0.4318`, with orbit-reduced RMS z = `0.4382`
- former anchor-ratio winner `(3.762226, 6.618519)`: `0.5374 -> 0.5417`, with orbit-reduced RMS z = `0.5430`

So for `l8_ratio`, the remaining mismatch is **not** mainly an origin pathology. The stronger interpretation is:

- the anchor-ratio off-target minimum is not normalization-stable
- once `l8_ratio` is used, the target miss is a broader shape/scoring bias, not just origin contamination

### Follow-Up: Direct L8-Ratio Winner Comparison

I then put the actual `l8_ratio` winner `(5.643339, 8.089301)` on the same direct full-point footing as the target.

For the `l8_ratio` winner `(5.643339, 8.089301)`:

- full-point aggregate RMS z = `0.2607`
- origin-dropped full-point aggregate RMS z = `0.2609`
- orbit-reduced RMS z = `0.2813`
- worst point is `(a,b) = (0.125, 0.000)` and its orbit partner, with `z = +0.6391`

For the target `(4.702783, 7.353910)` under the same normalization:

- full-point aggregate RMS z = `0.4285`
- origin-dropped full-point aggregate RMS z = `0.4318`
- orbit-reduced RMS z = `0.4382`
- worst point is `(a,b) = (0.250, 0.250)` and its orbit partner, with `z = -0.6666`

So the `l8_ratio` winner beats the target on every direct aggregate that was checked, not just on the orbit-reduced landscape score.

### Follow-Up: Landscape TSV Diagnostics

`plot_standard_pointwise_score_landscapes.py` now writes the pointwise aggregate diagnostics explicitly into the landscape TSV:

- `point_origin_dropped_rms_z`
- `point_origin_dropped_rms_z_jackknife_sigma`
- `point_origin_inclusive_rms_z`
- `point_origin_inclusive_rms_z_jackknife_sigma`

The older `point_rms_z` column was already the origin-dropped aggregate because the landscape builder excluded the origin before computing that pointwise RMS. The new columns make that fact explicit and add the origin-inclusive counterpart, so future sweeps can show immediately whether the origin is responsible for a bad ranking.

In the regenerated acute456 `l8_ratio` table, the difference is tiny for the current winner:

- `(5.643339, 8.089301)`: `0.260860` dropped-origin versus `0.260663` origin-inclusive

and likewise small for the target:

- `(4.702783, 7.353910)`: `0.431841` dropped-origin versus `0.428519` origin-inclusive

So the `l8_ratio` mis-ranking is not being driven by the origin.

### Status Of The Score Hypothesis

The current status of the hypothesis

- "we can build a score that picks the untwisted lattice which correctly realizes the target geometry from the twisted reference"

is: **not yet established**.

What the current evidence does support is the weaker statement:

- the score family contains real geometry information
- normalization matters a lot
- orbit-reduced anchor-ratio is currently the most stable primary acute456 metric among the tested choices

But the stronger target-selection claim is not yet supported, because:

- for acute456, `anchor_ratio` ranks the target only `2 / 25`
- for acute456, `l8_ratio` ranks the target only `7 / 25`
- for iso111, anchor-ratio fixes the raw `(0.8, 0.8)` inversion but still minimizes near `(1.1, 1.1)`, not exactly at `(1.0, 1.0)`
- the `l8_ratio` follow-up shows that the remaining acute456 failure is not mainly an origin artifact, so the scoring bias is deeper than one pathological point

So the honest summary is:

- we do have scores that are informative and nontrivial
- we do not yet have a normalization-stable score that reliably selects the correct untwisted target as the unique optimum
- the next scoring improvement has to address a broader shape-selection bias, not just the origin treatment

## 2026-05-21: Twisted-Reference 4x Trial

### Goal

Test the user hypothesis that the remaining twisted-versus-untwisted mismatch is mainly caused by the two continuum manifolds being sampled on different lattice point sets. The rerun changes the comparison workflow to mimic the intended coupling-to-geometry pipeline:

- build a denser twisted reference family
- fit the twisted continuum manifold first
- interpolate that twisted continuum manifold onto the untwisted continuum points
- use only the twisted-on-untwisted residuals to decide distinguishability

### Configuration

- Run config: `configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json`
- Comparison mode: `python scripts/align_manifolds.py --comparison-mode twisted_reference`
- Twisted family rule: base cell chosen to be the nearest integer multiple of the original pickGeo cell whose area is about `4 x` the untwisted `12x12` base, then scaled by `{1,2,3,4,5}`
- Untwisted family rule: keep the original size ladder `(12,24,36,48,60)`

### Run Status

- Completed benchmark: `geometry_acute_111`
- Completed comparison artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- Completed summary artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/SMOKE_SUMMARY.md`
- Interrupted benchmark: `geometry_acute_122`
- Interruption point: twisted scale 5, lattice `(160,100,40,20)`
- Interruption mode: manual `KeyboardInterrupt`, so no trustworthy campaign-wide conclusion exists yet for the 10-geometry set

### Measured Result For `geometry_acute_111`

This is the cleanest control point because both branches are isotropic and target the same equilateral torus.

| quantity | previous symmetric FSS run | 4x twisted-reference run |
| --- | ---: | ---: |
| run tag | `raw_manifold_fss_pickgeo10_fss_20260521` | `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` |
| comparison mode | `symmetric` | `twisted_reference` |
| twisted comparable points | `143` | `575` |
| untwisted comparable points | `143` | `143` |
| twisted alpha | `0.3047` | `0.4128` |
| twisted chi2/dof | `0.1843` | `0.2111` |
| untwisted alpha | `0.2863` | `0.4573` |
| untwisted chi2/dof | `0.0505` | `0.2859` |
| twisted-on-untwisted RMS z | `2.2363` | `3.2160` |
| untwisted-on-twisted RMS z | `2.1648` | `3.9717` |
| common-grid relative RMS | `0.0696` | `0.0907` |
| verdict | `distinguishable` | `distinguishable` |

### Interpretation

- The individual twisted and untwisted continuum fits still look healthy. Both chi2/dof values remain comfortably below `1` in the 4x trial.
- The asymmetric reference-style comparison does **not** reduce the disagreement for the equilateral control. The primary metric, twisted-on-untwisted RMS z, rises from `2.2363` to `3.2160`.
- The denser twisted family increases the twisted continuum point count from `143` to `575`, so this is not a case where the result is being driven by a lack of reference points.
- For `geometry_acute_111`, the hypothesis "the manifolds only look different because the point sets do not sit on top of each other" is therefore disfavored. Even when the twisted manifold is treated as the dense reference and sampled only at the untwisted continuum points, the mismatch stays well above the `Z = 2` marginal threshold.

### Operational Note

The 4x twisted-reference workflow is substantially more expensive than the previous symmetric run.

- `geometry_acute_111` twisted scale 5 already reached lattice `(120,120,0,0)` and took about `51 s` for one production job.
- `geometry_acute_122` was interrupted during twisted scale 5 at `(160,100,40,20)`.
- At `n_workers = 1`, the full 10-geometry campaign is practical but materially slower than the earlier FSS run.

### Current Scientific Takeaway

From the one completed control benchmark, denser twisted sampling by itself does not rescue indistinguishability. If the broader 10-geometry campaign is resumed, it should be treated as a stronger test of the continuum-equivalence claim, not as a formality expected to automatically collapse the manifolds.

### Next Useful Step

- Resume `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` from the current partial tree, ideally benchmark-by-benchmark, so runtime is easier to control.
- After at least a few additional geometries complete, regenerate the comparison plots with `--comparison-mode twisted_reference` and build a dedicated campaign summary for that run.

### Resume Later

Current saved state before stopping:

- run root: `responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/`
- completed benchmark manifest: `manifest_geometry_acute_111.json`
- completed comparison artifact: `geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- next blocked slice when interrupted: `geometry_acute_122`, twisted scale 5, lattice `(160,100,40,20)`
- config already has `"resume": true`, so rerunning the same generation command should reuse the completed payloads under the current run root

Resume command:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py \
	--config responsible_method_tests/configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json
```

After enough benchmarks finish, regenerate the asymmetric comparison artifacts with:

```bash
source /workspaces/newQFE/.venv/bin/activate
for manifest in responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_*.json; do
	python responsible_method_tests/scripts/align_manifolds.py \
		--benchmark-manifest "$manifest" \
		--comparison-mode twisted_reference
done
```

Then build the campaign summary with:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/summarize_benchmark.py \
	--comparison-mode twisted_reference \
	--output responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/PICKGEO10_TWISTED_REFERENCE4X_SUMMARY.md \
	--benchmark-manifest responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_111.json
```

and add the remaining `--benchmark-manifest` arguments once those benchmark manifests exist.

## 2026-05-21: Large-Twisted Fixed-Alpha Residual Collapse Test

### Goal

Test the simpler question behind the "large reference as continuum" idea without involving the untwisted branch: for one fixed twisted geometry, does the average per-site residual to a fixed-scale modular two-point function collapse as the twisted lattice is enlarged?

The residual was defined exactly as requested:

- choose one fixed geometry
- keep a single fixed modular normalization `alpha`
- compare the twisted MC connected correlator to `alpha * g_modular` at the corresponding lattice-point locations
- average the sitewise residuals over the matched non-origin lattice points

Because the CFT two-point function diverges near the origin while the MC correlator saturates at `1`, the comparison excluded the origin and used a single overall amplitude factor `alpha`.

### Benchmark Choice

I used `geometry_acute_567`, because it was already the hardest case for the earlier large-reference idea:

- continuum twisted vs continuum untwisted was **not distinguished**
- but the largest finite twisted lattice still sat visibly away from the untwisted continuum when sampled directly

That makes it a good stress test for the question "does a larger twisted lattice itself move toward the modular continuum?"

### Fixed-Alpha Definition

Two normalizations were checked.

1. The original twisted continuum-to-modular alignment from the corrected 5-scale run:
	 - `alpha = 0.21849098395920252`
2. The regenerated alignment after adding the larger twisted lattice:
	 - `alpha = 0.2189182337`

The qualitative conclusion was the same for both, so the result is not coming from retuning the amplitude.

### Practical Size Ceiling

The current all-to-all simulator cannot reach the user-suggested million-site twisted regime for this observable.

- `geometry_acute_567` scale `6` with lattice `(120,150,60,90)` runs successfully
- scale `7` already fails, even in an `n_traj=1` allocation probe
- the failure is therefore a simulator/storage ceiling, not a long-runtime issue

So the test could only be pushed to the **largest feasible twisted scale**, not to `O(10^6)` volume.

### Data Products

- Probe config: `configs/raw_manifold_large_ref_probe_567_20260521.json`
- High-stat twisted probe run root: `results/raw_manifold_large_ref_probe_567_20260521/`
- Fixed-alpha collapse script: `scripts/plot_fixed_alpha_modular_collapse.py`
- Final plot using regenerated run alpha:
	- `results/raw_manifold_large_ref_probe_567_20260521/geometry_acute_567/twisted/geometry_acute_567_twisted_fixed_alpha_modular_collapse.png`
- Strict fixed-alpha plot using the original pre-extension normalization:
	- `results/raw_manifold_large_ref_probe_567_20260521/geometry_acute_567/twisted/geometry_acute_567_twisted_fixed_alpha_modular_collapse_alpha021849.png`

### Main Result

Using the **strict original fixed normalization** `alpha = 0.21849098395920252`, the mean absolute residual per non-origin lattice site was:

| scale | family size | mean abs residual/site | RMS residual/site |
| --- | ---: | ---: | ---: |
| 1 | 1 | 0.19588758 | 0.19639174 |
| 2 | 2 | 0.11623277 | 0.11672905 |
| 3 | 3 | 0.08863659 | 0.08893889 |
| 4 | 4 | 0.06010132 | 0.06040861 |
| 5 | 5 | 0.05721346 | 0.05724360 |
| 6 | 6 | 0.04783519 | 0.04791435 |

This is a clear collapse with increasing twisted size. The scale-6 point lowers the mean absolute residual by about `16.4%` relative to scale 5:

- `0.05721346 -> 0.04783519`

Using the regenerated run alpha gives the same conclusion:

- scale 5 mean abs residual/site = `0.05668565`
- scale 6 mean abs residual/site = `0.04730739`

### Important Operational Detail

The first scale-6 attempt at only `n_traj=100` produced an apparent reversal (`~0.0858` at fixed old alpha), but a direct higher-stat check at `n_traj=500` showed that this was just a noisy outlier. After regenerating the scale-6 probe point at `n_traj=500`, the collapse resumed cleanly.

### Interpretation

- For this fixed geometry, the twisted lattice correlator does move toward the modular target as the twisted size is increased.
- The collapse is visible directly in the finite-lattice data using a single fixed amplitude normalization, without leaning on the untwisted branch.
- However, the current all-to-all implementation hits a hard ceiling already between scales `6` and `7` for this geometry, so this workflow cannot yet test the stronger hypothesis that an arbitrarily large twisted reference can simply replace the continuum limit.

### Current Takeaway

The finite-size trend supports the idea that a larger twisted reference is moving in the right direction, but with the present simulator the test stops well before the truly asymptotic regime. For non-integrable theories, this fixed-alpha twisted-vs-modular collapse is still a useful diagnostic, but it is not yet a substitute for a true continuum extrapolation.

## 2026-05-22: Matched-Volume Square-Untwisted Cross-Checks for 4-5-6 and 5-6-7

### Goal

Test the stricter matched-volume variant requested after the shared-coordinate 4-5-6 plots:

- keep the untwisted branch square, `L x L`
- choose twisted and untwisted base lattices with approximately matched areas above `100`
- handle physical side orientation explicitly through the correct untwisted `embedding_cycles`
- push both branches through the same `1..10` size ladder

This was also the first pass where I exposed the post-fit manifold difference directly instead of only as one panel inside the standard comparison plot.

### Configurations

Two acute scalene benchmarks were run.

1. `geometry_acute_456_matched_volume`
	- config: `responsible_method_tests/configs/raw_manifold_fss_acute_456_matched_volume_20260522.json`
	- twisted base: `[12,12,6,2]`, volume `156`
	- untwisted base: square `[12,12,0,0]`, volume `144`
	- untwisted `embedding_cycles=[0,1]`

2. `geometry_acute_567_matched_volume`
	- config: `responsible_method_tests/configs/raw_manifold_fss_acute_567_matched_volume_20260522.json`
	- twisted base: `[12,15,6,9]`, volume `234`
	- untwisted base: square `[15,15,0,0]`, volume `225`
	- untwisted `embedding_cycles=[2,1]`

Both runs used `n_traj=2000` and scales `1..10` on both branches.

### New Diagnostic

I added `responsible_method_tests/scripts/plot_continuum_difference.py` to render the post-fit difference explicitly from the same shared-cell comparison data already built by `align_manifolds.py`.

For each benchmark it writes three panels:

- common-grid heatmap of `twisted - untwisted`
- differences sampled on the untwisted fitted points
- differences sampled on the twisted fitted points

This makes it easier to see whether a result is just noisy pointwise scatter or a coherent smooth offset across the target cell.

### Data Products

For `geometry_acute_456_matched_volume`:

- run root: `responsible_method_tests/results/raw_manifold_fss_acute_456_matched_volume_20260522/`
- shared-coordinate FSS: `geometry_acute_456_matched_volume_shared_coordinate_twisted_untwisted_fss.{png,tsv}`
- shared-cell comparison: `geometry_acute_456_matched_volume_twisted_vs_untwisted_common_grid.{png,json}`
- dedicated post-fit difference: `geometry_acute_456_matched_volume_continuum_difference.png`

For `geometry_acute_567_matched_volume`:

- run root: `responsible_method_tests/results/raw_manifold_fss_acute_567_matched_volume_20260522/`
- shared-coordinate FSS: `geometry_acute_567_matched_volume_shared_coordinate_twisted_untwisted_fss.{png,tsv}`
- shared-cell comparison: `geometry_acute_567_matched_volume_twisted_vs_untwisted_common_grid.{png,json}`
- dedicated post-fit difference: `geometry_acute_567_matched_volume_continuum_difference.png`

### Measured Results

| quantity | 4-5-6 matched-volume | 5-6-7 matched-volume |
| --- | ---: | ---: |
| twisted chi2/dof | `0.0120` | `0.2011` |
| untwisted chi2/dof | `0.1547` | `0.1234` |
| common-grid RMS abs | `0.00754` | `0.01975` |
| common-grid relative RMS | `0.02877` | `0.07700` |
| twisted-on-untwisted RMS z | `0.7497` | `1.8487` |
| untwisted-on-twisted RMS z | `0.6789` | `1.9597` |
| verdict | `not distinguished within continuum point errors` | `marginally separated` |

### Shared-Coordinate FSS Cross-Check

The pointwise FSS plots tell the same story.

For `4-5-6`, the auto-selected points were `point_id = {1,2,4,6}`. The fitted continuum differences were small:

- `+0.00241` (`z=+0.31`)
- `-0.00437` (`z=-0.41`)
- `-0.00486` (`z=-0.35`)
- `-0.00558` (`z=-0.37`)

For `5-6-7`, the auto-selected points were `point_id = {1,3,5,7}`. All four fitted differences had the same sign and somewhat larger magnitude:

- `-0.00510` (`z=-0.75`)
- `-0.01505` (`z=-1.32`)
- `-0.01432` (`z=-1.04`)
- `-0.01523` (`z=-1.02`)

So the matched-volume `5-6-7` result is not just one bad point; it looks more like a coherent offset at the tested shared coordinates.

### Interpretation

- The earlier mild `4-5-6` separation was indeed a side-orientation artifact. Once the untwisted basis was corrected to the proper physical-side map, the benchmark moved to the clearly safe side of the continuum-error threshold.
- Repeating the same matched-volume protocol on `5-6-7` does **not** produce the same collapse. Both individual fits remain healthy, but the two continua retain a smooth, mostly positive `twisted - untwisted` offset across most of the shared target cell, with bidirectional RMS z values near `2`.
- The scientific lesson is therefore narrower than the optimistic `4-5-6` control alone suggested: correct orientation handling is necessary, but it is not by itself sufficient to make every matched-volume square-untwisted comparison collapse.

### Current Takeaway

At current statistics and with the present Taylor-2 continuum fits:

- `4-5-6` supports the claim "same continuum within current errors"
- `5-6-7` remains a borderline counterexample, landing in the `marginally separated` regime even after the orientation-aware matched-volume construction

The dedicated difference plots make that split especially clear: `4-5-6` shows a small residual field at the level of the fitted point errors, while `5-6-7` shows a larger smooth offset of order `2e-2` across much of the cell.

## 2026-05-22: Exact-Critical 1-1-1 Control Sweep

### Workflow Changes

Two continuum-workflow fixes landed before the control reruns:

- `K_from_continuum/workflow_common.py` now supports `beta_finder.mode = "exact_triangular_sinh_rule"` via `exact_triangular_ising_beta(...)`, so the triangular Ising critical point can be supplied exactly from
	`sinh(2 beta k1) sinh(2 beta k2) + sinh(2 beta k2) sinh(2 beta k3) + sinh(2 beta k3) sinh(2 beta k1) = 1`
- `K_from_continuum/02_generate_tests.py` stopped incorrectly forwarding `n_workers=` into `find_beta_for_lattice(...)`, which was breaking point jobs that used continuum-based beta strategies

The exact-beta path was validated immediately on the isotropic `1-1-1` case, giving `beta = ln(3)/4` as expected.

### Runs

I used that exact-critical path for three increasingly sharp `1-1-1` continuum-vs-reference controls:

- Pilot axial stencil: `K_from_continuum/results/iso111_distinguishabilityK_20260522/`
- Dense high-stat axial stencil: `K_from_continuum/results/iso111_distK_dense_hi5000_20260522/`
- Fresh narrow repeat control: `K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/`

All three used a twisted equilateral reference family and square untwisted test families at the exact critical couplings.

### Raw-Score Results

At the original raw continuum score, the equilateral control did **not** stabilize to the isotropic point as statistics increased.

- Pilot run: best point moved to `r1 = 1.01`, `r2 = 1.00`, with `z_rms = 0.6682`; isotropic `(1.00, 1.00)` still sat at `z_rms = 1.1588`
- Dense high-stat run: best point moved again to `r1 = 0.995`, `r2 = 1.00`, with `z_rms = 1.8069`; isotropic `(1.00, 1.00)` worsened to `z_rms = 3.2136`
- Fresh repeat run: the same off-isotropic preference repeated, with best point `r1 = 0.995`, `r2 = 1.00`, `z_rms = 2.2419`; isotropic `(1.00, 1.00)` still sat at `z_rms = 2.7362`

### Interpretation

- The failure mode sharpened with more statistics, which strongly argues against blaming the issue on the exact beta finder itself
- The repeat run showed that the bias was stable enough to reproduce, so the problem had to be in the reference/scoring construction rather than in one noisy outlier dataset
- That made the next question sharper: what part of the raw score is the workflow over-penalizing?

## 2026-05-23: Additive-Offset and Shape-Only Rescoring of the 1-1-1 Control

### Workflow Change

To answer that question, `K_from_continuum/03_score_continuum_vs_reference.py` now supports

- `analysis.score_mode = raw`
- `analysis.score_mode = additive_offset`
- `analysis.score_mode = shape_only`

and writes `score_diagnostics.dat` so the fitted nuisance offsets are visible instead of being implicit.

The idea is simple:

- `raw`: score the channels exactly as before
- `additive_offset`: subtract one best-fit common offset before scoring
- `shape_only`: subtract one best-fit offset per cycle before scoring

### Rescoring Results

Offline rescoring of the existing high-stat `1-1-1` artifacts changed the interpretation materially.

For the dense `hi5000` run:

- raw best: `(0.995, 1.000)` with `z_rms = 1.8069`
- additive-offset best: `(1.000, 1.000)` with `z_rms = 0.1999`
- shape-only best: `(1.000, 1.000)` with `z_rms = 0.0950`

For the repeat `hi5000` run:

- raw best: `(0.995, 1.000)` with `z_rms = 2.2419`
- additive-offset best: `(1.000, 1.000)` with `z_rms = 0.2839`
- shape-only best: `(0.995, 1.000)` with `z_rms = 0.1029`, but isotropic `(1.000, 1.000)` is still only `z_rms = 0.2145`

### Interpretation

- The dominant raw-score failure is a common-mode additive mismatch, not a large shape disagreement in the continuum channels
- Once that nuisance direction is removed, the isotropic point is either recovered outright or becomes essentially degenerate with the tiny `0.5%` offset candidate
- The `1-1-1` exact-critical problem is therefore much more about how the score is defined than about how `beta_c` is being chosen

## 2026-05-23: Responsible-Method Local 7x7 Iso111 Untwisted Scan

I also ran a narrow responsible-method local scan around the isotropic square target:

- run root: `responsible_method_tests/results/raw_manifold_fss_iso111_local_scan_grid7x7_hi5000_20260523/`
- target twisted manifest: `results/raw_manifold_fss_integer_multiples_20260520/geometry_111/twisted/manifest_geometry_111_twisted.json`

This did **not** recover isotropic couplings.

- best candidate: `(r1, r2) = (1.005, 1.020)`, score `1.3966`, interpretation `marginally separated`
- isotropic `(1.000, 1.000)`: score `3.2217`, interpretation `distinguishable at continuum-point level`

So on the responsible-method side, the near-isotropic local landscape was still favoring a noticeably shifted untwisted square candidate even after the continuum rescoring result above.

## 2026-05-24 to 2026-05-25: Wide Iso111 Responsible-Method Grid and Path-Length Fix

### Root Cause Fix

The wide untwisted square scan initially crashed at `(r1, r2) = (0.5, 1.0)`, but this turned out to be a simulator output-path bug rather than a Monte Carlo instability.

`K_from_continuum/src/ising_tri_twisted_parallelogram.cc` had been building runtime output subdirectories with fixed-size `char` buffers, and the generated scratch path hit that limit exactly for some of the widened-grid labels. That code now uses `std::string` for the run subdirectory and output paths.

### Production Run

After the path fix, the full widened responsible-method scan completed:

- twisted target run: `responsible_method_tests/results/raw_manifold_fss_iso111_twisted_target_hi10000_20260524/`
- untwisted square grid run: `responsible_method_tests/results/raw_manifold_fss_iso111_grid_r0p5to1p5_r0p5to1p5_step0p1_hi10000_20260524/`
- fit summary: `responsible_method_tests/results/raw_manifold_fss_iso111_grid_r0p5to1p5_r0p5to1p5_step0p1_hi10000_20260524/fit_vs_twisted_hi10000/geometry_111_coupling_fit_scan.json`

### Result

The best responsible-method square candidate on the full `11 x 11` grid was

- `iso111_scan_r11p200_r20p900`
- `(r1, r2) = (1.2, 0.9)`
- score `0.32237`
- twisted-on-untwisted RMS z `0.31967`
- untwisted-on-twisted RMS z `0.32237`
- common-grid relative RMS `0.01311`

The widened scan therefore still prefers a significantly anisotropic square candidate over isotropic `1-1-1`, even though the continuum rescoring analysis above showed that the raw continuum score was largely being driven by a nuisance offset.

## 2026-05-25: Fresh Responsible-Method Control Reruns and Larger-Lattice Follow-Up

### Fresh Small-Ladder Rerun

To check whether the wide-grid preference was just coming from a stale or unlucky reference, I rebuilt a fresh twisted `geometry_111` reference and compared it directly against two untwisted square candidates:

- run root: `responsible_method_tests/stupid_method/`
- config: `responsible_method_tests/configs/raw_manifold_fss_stupid_method_20260525.json`
- candidates: `(1.0, 1.0)` and `(1.2, 0.9)`

On the four selected shared-coordinate panels, the earlier dramatic isotropic self-mismatch did **not** reproduce.

- isotropic `(1.0, 1.0)` panel z values: `{-0.965, -0.832, -0.604, -0.441}`
- `(1.2, 0.9)` panel z values: `{-0.740, -0.809, -0.499, -0.199}`

So the fresh rerun puts both candidates inside `|z| < 1` on those selected shared points, with `(1.2, 0.9)` still somewhat closer.

### Larger-Lattice High-Stat Follow-Up

I then pushed the same control to larger ladders and higher statistics:

- run root: `responsible_method_tests/stupid_method/large_hi20000/`
- config: `responsible_method_tests/configs/raw_manifold_fss_stupid_method_large_hi20000_20260525.json`
- twisted scales: `1..7`
- untwisted sizes: `12..84`
- `n_traj = 20000`

On the same four selected shared-coordinate panels, the larger-lattice run still favors `(1.2, 0.9)` over isotropic, although isotropic improved materially.

- isotropic selected-panel RMS z: `0.3852`
- `(1.2, 0.9)` selected-panel RMS z: `0.1035`

### Small Workflow Extensions

Two small script changes supported these reruns:

- `responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py` now accepts explicit `entries=[{scale, family_size, lattice}, ...]`, so multiplicative scale and displayed family size can be decoupled when a family is not naturally expressed as one simple size ladder
- `responsible_method_tests/scripts/fit_target_geometry.py` now has `--annotate-points` to overlay candidate labels directly on the scan plot when that is useful

I also ran the size-`10..100` comparison bundle under `responsible_method_tests/stupid_method/iso111_and_far456_size10to100/`, which exported about `100` standalone shared-point panels per benchmark for slide triage.

## 2026-05-25: Far-Point Responsible-Method Control

To see whether the responsible-method ranking was just measuring distance from isotropic couplings, I compared the fresh twisted `1-1-1` reference against a square untwisted family at the **exact** `4-5-6` couplings:

- run root: `responsible_method_tests/stupid_method/far_456_vs_iso111/`
- config: `responsible_method_tests/configs/raw_manifold_fss_far_456_vs_iso111_20260525.json`
- square untwisted couplings: `(r1, r2) = (5.065231170850374, 7.742930488974055)`

The result was much less separated than naive “distance in coupling space” intuition would suggest.

- primary RMS z: `0.60139`
- twisted-on-untwisted RMS z: `0.59148`
- common-grid relative RMS: `0.02216`
- interpretation: `not distinguished within continuum point errors`

The selected shared-point z values were only `{1.64, 0.73, 0.50, 0.41}`.

So at current errors, the responsible-method comparison is clearly not just “how far are the couplings from `(1, 1)`?” A square untwisted family at the exact `4-5-6` point can still sit surprisingly close to the fresh twisted `1-1-1` reference.

## 2026-05-25: Acute 4-5-6 Literal-Size Blind Holdout

### Goal

Build a real blind holdout for the acute `4-5-6` shared-boundary midpoint using

- a literal-size untwisted square family
- a twisted family built from explicit multiples of the small `4-5-6` base `[6,6,3,1]`
- one held-out twisted lattice that is large enough to lie well outside the training window

### Setup

The final high-stat run used

- run root: `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_hi40000_notaylor/`
- twisted training family: `[6,6,3,1] * {2,4,6,8,10,12,16}`
- untwisted training sizes: `{8,16,24,32,48,64,96}`
- training MC: `n_traj = 40000`, `n_therm = 1000`
- blind holdout MC: `n_traj = 10000`, `n_therm = 200`

The blind datum was the smallest twisted multiple with the chosen boundary midpoint on-lattice and volume above `20000`:

- holdout lattice: `[144,144,72,24]`
- holdout volume: `22464`
- measured connected correlator: `0.28278461 +/- 0.00503771`
- holdout production runtime: about `54.54 s`

### Holdout Result

The blind-holdout comparison did **not** pick one universal best extrapolant across both branches.

For the twisted branch:

- best model: `power_free`
- holdout z: `-0.4325`
- training chi2/dof: `1.0090`

For the untwisted branch:

- best model: `pade21`
- holdout z: `-0.3269`
- training chi2/dof: `0.7234`

The untwisted `power_free` fit was still acceptable at the blind point (`z = +0.5831`), but `pade11` over-predicted badly (`z = +3.2404`).

### Interpretation

- Both branches can extrapolate to the held-out twisted datum within `|z| < 1`
- The preferred asymptotic ansatz is branch-dependent for this observable and this training window, not universal
- A lower-stat predecessor under `acute456_literal_sizes_holdout20k_hi20000/` landed on the same holdout geometry but gave noisier model rankings, so the `hi40000_notaylor` rerun is the one worth trusting

## 2026-05-25: Single-Displacement Holdout Mode for Large Blind Points

### Code Change

`K_from_continuum/src/ising_tri_twisted_parallelogram.cc` now accepts

- `--single_disp_m`
- `--single_disp_n`

and, when those are present, builds only the requested displacement map instead of the full `N_cell x N_cell` add table.

That matters because the full add table is the memory bottleneck for large blind-holdout points.

### Workflow

I used that new mode in the `responsible_method_tests/stupid_method/iso111_base4_size4to100_step4/` workflow:

- config: `responsible_method_tests/configs/raw_manifold_fss_iso111_base4_size4to100_step4_20260525.json`
- run root: `responsible_method_tests/stupid_method/iso111_base4_size4to100_step4/`

### Results

For the quick `L = 200` single-point holdout at point `pt002` (`m = 0`, `n = -100`):

- measured connected correlator: `0.27914208 +/- 0.00908915`
- the selected-model quick summary is written in `holdout_L200_pt002_quick_summary.json`

For the real `L = 1000` blind holdout at the corresponding scaled point (`m = 0`, `n = -500`):

- measured connected correlator: `0.17363803 +/- 0.01114078`
- output summary: `holdout_L1000_pt002_selected_models_logx.json`

At `L = 1000`, the asymptotic `power_free` model is clearly the best model on **both** branches at the blind point:

- untwisted `power_free`: `z = -0.4318`
- twisted `power_free`: `z = +0.8964`

while the Taylor and Padé models miss by much larger margins.

### Interpretation

- Single-point blind holdouts at volumes that would be awkward or impossible for the full all-to-all add-table workflow are now practical
- That gives a much cleaner direct test of continuum extrapolants against genuinely out-of-sample large lattices

## 2026-05-25: Acute 4-5-6 Beta-Finder Holdout Rerun Stopped Midway

### Context

I started a production rerun of the acute `4-5-6` literal-size blind-holdout workflow using the new per-lattice `beta_c` finder path:

- script: `responsible_method_tests/scripts/run_acute456_pow2_blind_holdout.py`
- run root: `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/`
- training MC: `n_traj = 40000`, `n_therm = 1000`
- holdout MC: `n_traj = 10000`, `n_therm = 200`
- beta finder: `11 + 7 + 7 + 7` scan points with `n_traj = 120` on the coarse pass and `240` on refinement passes, with jackknife enabled

The run was then intentionally killed before completion so I could switch focus to collaborator-facing project storytelling.

### What Completed Before The Stop

The production rerun did complete the held-out twisted lattice and most of the twisted training ladder:

- completed holdout beta scan and holdout production point
- completed twisted training beta scans and production points for
	- `[12,12,6,2]`
	- `[24,24,12,4]`
	- `[36,36,18,6]`
	- `[48,48,24,8]`
	- `[60,60,30,10]`
	- `[72,72,36,12]`
- started but did not finish the twisted `96x96_t48x16` beta scan
- did **not** reach any untwisted square beta scans
- therefore did **not** write the final blind-holdout summary JSON or any of the production PNG figures for this rerun

The last log line before termination was during the `96x96_t48x16` twisted scan:

- pass `2/4`, point `4/7`, `beta = 0.27024485`

### Partial Numerical Takeaway

The holdout beta scan itself was successful and peaked close to the exact isotropic triangular critical value `ln(3)/4 = 0.27465307`:

- holdout lattice: `[144,144,72,24]`
- holdout `beta_c = 0.27350859 +/- 0.00099267`
- holdout scan peak susceptibility: `0.0405321`
- holdout beta-scan wall time: about `908 s`

The completed twisted training pseudo-critical estimates showed a sensible finite-size trend toward the exact value:

- `12x12_t6x2`: `beta_c = 0.25699835 +/- 0.00175`
- `24x24_t12x4`: `beta_c = 0.26441420 +/- 0.0206`
- `36x36_t18x6`: `beta_c = 0.26985540 +/- 0.00216`
- `48x48_t24x8`: `beta_c = 0.26919364 +/- 0.00149`
- `60x60_t30x10`: `beta_c = 0.26905701 +/- 0.00108`
- `72x72_t36x12`: `beta_c = 0.27003786 +/- 0.00274`

So the beta-finder path looks operational on the twisted branch and the partial drift is physically reasonable, but the study is scientifically incomplete until the final twisted size, the full untwisted ladder, and the blind extrapolation fits are all finished.

### Artifacts Worth Reusing Later

- completed holdout beta summary JSON:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/beta_scans/holdout_twisted/acute456_holdout_smallbase024_m-108_n60_beta_scan.json`
- completed twisted beta summary JSONs currently exist for the six finished training sizes under:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/beta_scans/twisted/`
- completed MC production payloads currently exist under:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/mc_runs/holdout_twisted/`
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/mc_runs/twisted/`

Because the beta-scan summaries are cached per lattice, resuming later should not need to recompute the already-finished holdout and twisted sizes if the runner is pointed back at the same output tree without `--force`.

## 2026-05-27: Standard Iso111 Anchor-Ratio Landscape Note

I added a separate standard-workflow journal note for the iso111 pointwise target-score rerun with normalization turned on:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That note records the raw-score bias diagnosis, the new `anchor_ratio` scoring path, and the resulting ranking shift.

The short version is:

- raw scoring had favored `(0.8, 0.8)` over `(1.0, 1.0)`
- under anchor-ratio normalization, `(1.0, 1.0)` becomes clearly better than `(0.8, 0.8)`
- but the full anchor-ratio landscape still minimizes at `(1.1, 1.1)` rather than at isotropic `(1.0, 1.0)`

So the raw inversion was indeed a normalization-sensitive artifact, but normalization alone did not fully restore isotropic `1-1-1` as the unique optimum in the standard pointwise target score.

## 2026-05-28: Standard Boundary fit_method Follow-up Note

I added a substantial follow-up entry to the standard-workflow journal note:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That update records three concrete outcomes from the boundary-score session:

- the completed `iso111` quarter-point jackknife boundary summary and its surviving direction-ratio contribution
- the new `acute456` dual-target midpoint / quarter boundary sweeps, including separate large-target and small-target summary tables
- the new `fit_method` prototype, which fits the twisted and untwisted ladders separately before comparing continuum intercepts

The short version is:

- the new `fit_method` is not yet convincing; both the continuum-fit and blind-power variants fail to restore isotropic `1-1-1`, and the ratio-panel extrapolations look unstable
- the older direct boundary `chi2` summaries, including the direction-ratio channels, still look more promising and may contain real signal rather than pure noise

So the practical session conclusion is currently skeptical of the new fit-vs-fit score, but still interested in the old boundary `chi2` + direction-ratio line of attack.

## 2026-05-28: Acute456 Large-L Follow-up Hypothesis Note

I added another update to the standard-workflow note:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

This one records the current interpretation of the `acute456` small-target versus large-target behavior.

The short version is:

- the smaller twisted target may be looking better mainly because the present untwisted ladder only has to propagate from `64 -> 66`, which is almost a local holdout
- the large target currently asks the untwisted fit to bridge `64 -> 144`, which is a much longer extrapolation and may dominate the score more than the twisted/untwisted cutoff mismatch itself

The detailed note now includes a minimal next test:

- a `12`-job untwisted-only large-`L` campaign for the quarter-point `acute456` score
- shortlist = exact target, current large-target winner, current small-target winner
- new untwisted sizes = `96, 112, 128, 144`

The point of that campaign is to answer one narrow question before widening scope again:

- does the exact `acute456` target recover when the large-target comparison is made local near `L = 144`, instead of being extrapolated all the way from `L = 64`?

## 2026-05-28: Acute456 Ranking Diagnostics And Next-Step Freeze

I extended the same detailed standard-workflow note again:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That latest update records the quarter-point `acute456` ranking experiments and the current decision about what **not** to do next.

The short version is:

- a mixed squared-rank aggregation across the quarter-point boundary metrics recovers the exact target only for the smaller twisted target
- ratio-only and correlator-only rank aggregations do **not** reproduce that exact-target win, so the effect is hybrid and heuristic rather than cleanly tied to one observable block
- same-size untwisted-to-twisted reweighting may be interesting later, but it does not solve the present `64 -> 144` size-extrapolation problem

So the ranking pass was useful as diagnosis, but the current decision is to stop inventing new scores and move to the cleaner experiment:

- a 12-job `acute456` quarter-point untwisted large-`L` shortlist
- sizes `96,112,128,144`
- couplings = exact target, current large-target winner, current small-target winner
- start at `n_traj = 40000`
