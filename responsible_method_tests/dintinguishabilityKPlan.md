# Distinguishability In K-Space Plan

## Purpose

This document records a first pilot of the twisted-reference workflow on the
`1-1-1` geometry and lays out the next steps needed to turn it into a reliable
`\Delta K` distinguishability test.

The immediate goal is narrow:

> Build a reference from the twisted-family side, evaluate a small untwisted
> coupling stencil against it, and determine whether nearby points in coupling
> space are already distinguishable.

For this pilot we intentionally skip the Monte Carlo `beta_c` finder and use
the exact triangular Ising criticality condition

`sinh(2 beta k1) sinh(2 beta k2) + sinh(2 beta k2) sinh(2 beta k3) + sinh(2 beta k3) sinh(2 beta k1) = 1`

so that the test isolates the reference-comparison workflow rather than mixing
it with finite-size critical-tuning noise.

## Why The `1-1-1` Case Is Useful

The `1-1-1` geometry is the cleanest possible control for this workflow.

- the twisted reference family is the isotropic zero-twist triangular torus
- the untwisted family can be moved through nearby anisotropic coupling points
  in a controlled one-parameter stencil
- any failure here is easier to interpret than in strongly deformed targets
  such as `4-5-6` or `5-6-7`

This makes `1-1-1` the right place to answer a first question:

Can the twisted-reference score tell apart nearby points in untwisted coupling
space before we reintroduce the `beta_c` finder?

## Pilot That Was Run

### Reference Family

- workflow: `K_from_continuum` Step 1/2/3
- run tag: `iso111_distinguishabilityK_20260522`
- geometry: `(L,L,0,0)`
- reference sizes: `L = {24, 48, 72, 96}`
- reference couplings: `k1 = k2 = k3 = 1`
- production stats: `n_traj = 2000`, `n_skip = 10`, `n_therm = 1500`
- continuum channel fit: second-order Taylor form `A + B/L + C/L^2`

### Untwisted Test Stencil

We used a five-point axial slice through the isotropic point:

- `k3 = 1`
- `k2/k3 = 1`
- `k1/k3 in {0.98, 0.99, 1.00, 1.01, 1.02}`
- test sizes: `L = {12, 24, 36, 48, 60}`

### Three-Directional Score Object

The score uses three directional channels with seven fractions on each:

- cycles `0, 1, 2`
- fractions `t = {1/8, 2/8, ..., 7/8}`
- total scored continuum channels per point: `21`

This is the closest currently implemented proxy to the three-directional
twisted-reference workflow in `workflowPlan.md`.

## Pilot Result

### Reference-Vs-Test Ranking

The scored output is:

| `k1/k3` | exact `beta_c` | `Delta K1 = beta*k1 - beta_iso` | score | `z_rms` | interpretation vs twisted reference |
| --- | ---: | ---: | ---: | ---: | --- |
| `0.98` | `0.276500` | `-3.6833e-03` | `0.03919` | `1.1460` | marginal |
| `0.99` | `0.275572` | `-1.8363e-03` | `0.16255` | `2.0286` | distinguishable |
| `1.00` | `0.274653` | `0` | `0.04139` | `1.1588` | marginal |
| `1.01` | `0.273741` | `+1.8258e-03` | `0.01148` | `0.6682` | not distinguished |
| `1.02` | `0.272837` | `+3.6411e-03` | `0.04942` | `1.1678` | marginal |

The raw minimum did **not** land at the isotropic point. The best match to the
twisted reference on this stencil was `k1/k3 = 1.01`.

That matters because it means this pilot does **not** yet give a clean,
monotonic local distance-to-reference curve.

### Per-Cycle RMS Z Against The Reference

| `k1/k3` | cycle 0 | cycle 1 | cycle 2 |
| --- | ---: | ---: | ---: |
| `0.98` | `1.0857` | `1.0326` | `1.3020` |
| `0.99` | `1.8659` | `1.9512` | `2.2489` |
| `1.00` | `1.0298` | `1.3270` | `1.0987` |
| `1.01` | `0.7159` | `0.5598` | `0.7167` |
| `1.02` | `1.3171` | `1.1238` | `1.0459` |

The `0.99` point is not being flagged by a single pathological direction; the
excess is distributed across all three channels.

### Pairwise Test-Test Separation

The most relevant pairwise RMS z values on the five-point stencil are:

| pair | pairwise `z_rms` | interpretation |
| --- | ---: | --- |
| `0.98` vs `0.99` | `1.1329` | not distinguished |
| `0.99` vs `1.00` | `1.1338` | not distinguished |
| `1.00` vs `1.01` | `0.7325` | not distinguished |
| `1.01` vs `1.02` | `0.7472` | not distinguished |
| `0.99` vs `1.01` | `1.7385` | still below clear separation |

No pairwise test-test comparison on this stencil exceeded `z_rms = 2`.

## What This Means Right Now

The pilot supports three conclusions.

### 1. The workflow mechanics are working

The exact-critical shortcut, twisted-reference generation, untwisted test
family generation, and continuum channel scoring all completed successfully.

That means the basic workflow structure is viable.

### 2. A nearby point can already register against the reference

The `k1/k3 = 0.99` point reached `z_rms = 2.03` against the twisted reference,
corresponding to

- `Delta (k1/k3) = -0.01`
- `Delta K1 = -1.8363e-03`

So the current workflow can already see a percent-level offset on one side of
the isotropic point.

### 3. The pilot does not yet define a clean threshold `Delta K`

Because the minimum landed at `1.01` instead of `1.00`, the local response is
not monotonic.

In practice that means:

- `Delta(k1/k3) = 0.01` is **sometimes** visible against the reference
- the same nominal offset on the opposite side is **not** visible in the same
  way
- even `0.99` vs `1.01` only gave pairwise `z_rms = 1.74`

So the conservative reading is:

> This five-point pilot does not yet establish a stable one-number
> distinguishability threshold in `K` space.

The strongest safe statement is that the current implementation is sensitive to
offsets of order `|Delta K1| ~ 2e-3`, but not yet stably enough to quote that as
the real threshold.

## Dense High-Statistic Follow-Up

After the first five-point pilot, we ran the next exact-critical control that
was recommended before touching the `beta_c` finder.

### What Changed

- denser stencil in `k1/k3`
- higher production statistics for both the reference family and the untwisted
  test family
- same three-directional channel score
- same exact sinh-rule criticality

### Follow-Up Setup

- run tag: `iso111_distK_dense_hi5000_20260522`
- reference sizes: `L = {24, 48, 72, 96}`
- test sizes: `L = {12, 24, 36, 48, 60}`
- production stats: `n_traj = 5000`, `n_skip = 10`, `n_therm = 1500`
- stencil:
  `k1/k3 = {0.990, 0.995, 1.000, 1.005, 1.010, 1.015, 1.020}` with `k2/k3 = 1`

### Reference-Vs-Test Ranking At Higher Statistics

| `k1/k3` | `Delta K1 = beta*k1 - beta_iso` | score | `z_rms` | interpretation vs twisted reference |
| --- | ---: | ---: | ---: | --- |
| `0.990` | `-1.8363e-03` | `0.22201` | `3.8567` | distinguishable |
| `0.995` | `-9.1683e-04` | `0.03979` | `1.8069` | marginal |
| `1.000` | `0` | `0.13609` | `3.2136` | distinguishable |
| `1.005` | `+9.1420e-04` | `0.13036` | `3.1077` | distinguishable |
| `1.010` | `+1.8258e-03` | `0.17046` | `3.4712` | distinguishable |
| `1.015` | `+2.7348e-03` | `0.09321` | `2.6537` | distinguishable |
| `1.020` | `+3.6411e-03` | `0.22637` | `4.0139` | distinguishable |

The best point on the denser stencil moved again, this time to `k1/k3 = 0.995`,
but even that point was only marginal, not truly passing.

The most important fact is not the shifted minimum. It is that the exact
isotropic point `k1/k3 = 1.000` became **more** separated from the twisted
reference when the statistics were increased.

### Pairwise Test-Test Separation At Higher Statistics

The denser run finally produced some pairwise test-test separations above the
`z_rms = 2` level.

The clearest examples are:

| pair | pairwise `z_rms` | interpretation |
| --- | ---: | --- |
| `0.990` vs `0.995` | `2.5220` | distinguishable |
| `0.995` vs `1.010` | `2.0156` | distinguishable |
| `0.995` vs `1.020` | `2.6499` | distinguishable |

But this does **not** rescue the control, because the exact isotropic point is
still not the preferred point on the stencil.

### What The Follow-Up Means

This second run changes the scientific interpretation in an important way.

The first pilot could still be read as: the workflow sees a signal, but the
surface is noisy.

The denser high-stat pilot instead says:

> Increasing statistics and local resolution did not restore the isotropic
> control point. It sharpened a mismatch.

That is stronger evidence that the remaining problem is **not** just sparse
sampling in `k` space.

At this stage, the main unresolved issue is not the `Delta K` threshold itself.
It is the fact that the exact-critical `1-1-1` control still does not score its
own isotropic point as the best, or even as a passing, match.

## Explicit Repeat Control Around The Isotropic Point

To test whether the shifted minimum was just a one-off fluctuation, we ran a
fresh reference and rescored the narrow central window directly.

### Repeat Setup

- run tag: `iso111_distK_repeat3_hi5000_20260522`
- fresh twisted reference family with the same sizes `L = {24, 48, 72, 96}`
- fresh untwisted tests with the same sizes `L = {12, 24, 36, 48, 60}`
- production stats: `n_traj = 5000`, `n_skip = 10`, `n_therm = 1500`
- exact sinh-rule criticality again
- repeat stencil: `k1/k3 = {0.995, 1.000, 1.005}` with `k2/k3 = 1`

### Reference-Vs-Test Ranking In The Repeat

| `k1/k3` | `Delta K1 = beta*k1 - beta_iso` | score | `z_rms` | interpretation vs twisted reference |
| --- | ---: | ---: | ---: | --- |
| `0.995` | `-9.1683e-04` | `0.06819` | `2.2419` | distinguishable |
| `1.000` | `0` | `0.11378` | `2.7362` | distinguishable |
| `1.005` | `+9.1420e-04` | `0.18473` | `3.4685` | distinguishable |

The minimum stayed at `k1/k3 = 0.995` even with a completely fresh isotropic
reference family.

That is the cleanest result in this whole control campaign so far: the isotropic
point did **not** recover when the reference itself was rerun.

### Per-Cycle RMS Z Against The Fresh Reference

| `k1/k3` | cycle 0 | cycle 1 | cycle 2 |
| --- | ---: | ---: | ---: |
| `0.995` | `2.3510` | `1.8087` | `2.5058` |
| `1.000` | `2.6612` | `2.5239` | `3.0012` |
| `1.005` | `3.5058` | `3.3877` | `3.5108` |

Again the excess is broad across all three directional channels. This is not a
single-cycle pathology.

### Pairwise Test-Test Separation In The Repeat

| pair | pairwise `z_rms` | interpretation |
| --- | ---: | --- |
| `0.995` vs `1.000` | `0.7102` | not distinguished |
| `0.995` vs `1.005` | `1.5131` | not distinguished |
| `1.000` vs `1.005` | `0.7993` | not distinguished |

So the repeat campaign says two things simultaneously:

- the fresh reference still prefers the off-isotropic point `0.995`
- the narrow three-point window is not internally resolved at `z_rms > 2`

That combination matters. It means the present failure mode is not simply
"local points are indistinguishable." The more specific problem is that the
twisted-reference scoring surface is biased away from the exact isotropic point
in a repeat-stable way.

## Offline Rescoring With Nuisance Removal

After the repeat control made the raw-score bias unambiguous, we re-ran Step 3
only on the existing dense and repeat artifacts with two nuisance-removed score
definitions.

No new Monte Carlo data were generated in this step. The reference and test
payloads were unchanged; only the continuum-to-reference score was changed.

### The Two Alternative Scores

- `additive_offset`: fit one weighted additive offset across all `21`
  continuum channels at a point, then score only the residuals around that
  offset
- `shape_only`: fit one weighted additive offset separately on each of the
  three directional cycles, then score only the residuals around those
  per-cycle offsets

The point of these rescoring passes is simple: test whether the raw control
failure is really a geometry mismatch, or whether the raw score is mostly
penalizing a common-mode amplitude shift.

### Dense Run Under The New Scores

| `k1/k3` | raw `z_rms` | additive-offset `z_rms` | shape-only `z_rms` |
| --- | ---: | ---: | ---: |
| `0.990` | `3.8567` | `0.4481` | `0.4282` |
| `0.995` | `1.8069` | `0.3174` | `0.2238` |
| `1.000` | `3.2136` | `0.1999` | `0.0950` |
| `1.005` | `3.1077` | `0.2781` | `0.1063` |
| `1.010` | `3.4712` | `0.3267` | `0.1859` |
| `1.015` | `2.6537` | `0.2319` | `0.1345` |
| `1.020` | `4.0139` | `0.4132` | `0.3619` |

This is a qualitative change.

With the raw score, the dense run preferred `k1/k3 = 0.995` and strongly
rejected the exact isotropic point.

With either nuisance-removed score, the dense run instead picks
`k1/k3 = 1.000` as the best point on the stencil.

### Repeat Run Under The New Scores

| `k1/k3` | raw `z_rms` | additive-offset `z_rms` | shape-only `z_rms` |
| --- | ---: | ---: | ---: |
| `0.995` | `2.2419` | `0.3031` | `0.1029` |
| `1.000` | `2.7362` | `0.2839` | `0.2145` |
| `1.005` | `3.4685` | `0.4957` | `0.4942` |

The repeat run tells a slightly subtler version of the same story.

- under `additive_offset`, the best point moves to the exact isotropic point
  `k1/k3 = 1.000`
- under `shape_only`, there is still a mild preference for `0.995`, but all
  three points are now deep in the sub-unit `z_rms` regime

So even the strictest nuisance-removed score no longer sees the central
three-point window as strongly separated from the reference.

### What The Rescoring Says About The Raw Bias

The diagnostics from the new score runs show that the raw mismatch was indeed
dominated by a common-mode offset.

For example:

- on the dense run at `k1/k3 = 1.000`, the mean raw channel difference was
  `+0.0803`, and the fitted global additive nuisance was `+0.0808`
- on the repeat run at `k1/k3 = 1.000`, the mean raw channel difference was
  `+0.0732`, and the fitted global additive nuisance was `+0.0711`

In other words, the raw score was mostly measuring a coherent upward shift of
all continuum channels together, not a directional shape failure.

## Revised Interpretation

The current state of the workflow is now:

1. The raw continuum-channel score is strongly contaminated by a common-mode
  additive offset.
2. Once that nuisance direction is removed, the exact isotropic point is
  recovered as the best point in the dense run and in the repeat run under the
  global additive-offset score.
3. The stricter per-cycle `shape_only` score still leaves a small residual
  preference for `0.995` in the repeat run, but the whole central window drops
  to `z_rms < 0.5`.

This changes the diagnosis in an important way.

The dominant problem is no longer best described as "the twisted-reference
workflow fails the `1-1-1` control." It is better described as:

> the original raw score was not quotienting out a common-mode continuum
> offset, so it was over-calling mismatches that are not really directional
> shape failures.

That is a much better place to be, but it still does **not** mean the workflow
is ready for a `beta_c`-finder run. The new score definition should be frozen
and repeat-validated first.

## Recommended Next Run Sequence

To turn this into a publishable or collaborator-facing `Delta K` threshold, the
next sequence should now be:

1. Promote nuisance-removed scoring to the primary control diagnostic.
  At present, `additive_offset` is the cleanest candidate because it restores
  the isotropic point in both dense and repeat rescoring.
2. Keep the raw score as a secondary diagnostic only, because it is still
  useful for tracking the common-mode offset itself.
3. Add explicit reference-vs-reference repeats at `1-1-1` and measure the
  repeat floor under the nuisance-removed score, not just under the raw score.
4. Test one more score-cleanup step on the same existing data: prune mirrored
  fractions or otherwise avoid double-counting the symmetric `t` pairs.
5. Only after the exact-critical isotropic control is stable under that final
  score definition should the same stencil be rerun with the production
  `beta_c` finder.

This sequence separates three questions cleanly:

- workflow mechanics
- raw-score contamination from common-mode offsets
- repeat stability of the nuisance-removed local score surface
- extra noise from numerical critical tuning

## Ultimate Test After The Exact-Critical Pilot

Once the exact-critical `1-1-1` pilot is stable, the production version should
be:

1. Same three-directional twisted-reference workflow.
2. Same local untwisted coupling stencil.
3. Replace the exact sinh-rule criticality with the actual `beta_c` finder.
4. Re-measure the repeat floor and the `Delta K` threshold.

If the threshold survives that replacement, then the geometry-grading part of
the workflow is not specific to the exactly solved Ising case.

That is the key prerequisite for porting the method to `phi^4`, where the exact
criticality shortcut is unavailable but the reference-vs-test geometry workflow
should still be the same.

## Artifacts From This Pilot

- reference config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_reference_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_reference_20260522.json)
- tests config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_tests_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_tests_20260522.json)
- score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_score_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_score_20260522.json)
- scored results root:
  [../K_from_continuum/results/iso111_distinguishabilityK_20260522](../K_from_continuum/results/iso111_distinguishabilityK_20260522)
- score table:
  [../K_from_continuum/results/iso111_distinguishabilityK_20260522/comparison_analysis_data/score_map.dat](../K_from_continuum/results/iso111_distinguishabilityK_20260522/comparison_analysis_data/score_map.dat)
- score heatmap:
  [../K_from_continuum/results/iso111_distinguishabilityK_20260522/comparison_analysis_data/score_zscore_heatmaps.png](../K_from_continuum/results/iso111_distinguishabilityK_20260522/comparison_analysis_data/score_zscore_heatmaps.png)

## Artifacts From The Dense High-Statistic Follow-Up

- reference config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_reference_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_reference_20260522.json)
- tests config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_tests_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_tests_20260522.json)
- score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_20260522.json)
- scored results root:
  [../K_from_continuum/results/iso111_distK_dense_hi5000_20260522](../K_from_continuum/results/iso111_distK_dense_hi5000_20260522)
- dense score table:
  [../K_from_continuum/results/iso111_distK_dense_hi5000_20260522/comparison_analysis_data/score_map.dat](../K_from_continuum/results/iso111_distK_dense_hi5000_20260522/comparison_analysis_data/score_map.dat)
- dense score heatmap:
  [../K_from_continuum/results/iso111_distK_dense_hi5000_20260522/comparison_analysis_data/score_zscore_heatmaps.png](../K_from_continuum/results/iso111_distK_dense_hi5000_20260522/comparison_analysis_data/score_zscore_heatmaps.png)

## Artifacts From The Explicit Repeat Control

- reference config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_reference_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_reference_20260522.json)
- tests config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_tests_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_tests_20260522.json)
- score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_20260522.json](../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_20260522.json)
- scored results root:
  [../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522](../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522)
- repeat score table:
  [../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/comparison_analysis_data/score_map.dat](../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/comparison_analysis_data/score_map.dat)
- repeat score heatmap:
  [../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/comparison_analysis_data/score_zscore_heatmaps.png](../K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/comparison_analysis_data/score_zscore_heatmaps.png)

## Artifacts From The Offline Rescoring Pass

- dense additive-offset score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_additive_offset_20260523.json](../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_additive_offset_20260523.json)
- dense shape-only score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_shape_only_20260523.json](../K_from_continuum/configs/iso111_distinguishabilityK_dense_hi5000_score_shape_only_20260523.json)
- repeat additive-offset score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_additive_offset_20260523.json](../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_additive_offset_20260523.json)
- repeat shape-only score config:
  [../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_shape_only_20260523.json](../K_from_continuum/configs/iso111_distinguishabilityK_repeat3_hi5000_score_shape_only_20260523.json)
- dense additive-offset results:
  [../K_from_continuum/results/iso111_distK_dense_hi5000_additive_offset_20260523](../K_from_continuum/results/iso111_distK_dense_hi5000_additive_offset_20260523)
- dense shape-only results:
  [../K_from_continuum/results/iso111_distK_dense_hi5000_shape_only_20260523](../K_from_continuum/results/iso111_distK_dense_hi5000_shape_only_20260523)
- repeat additive-offset results:
  [../K_from_continuum/results/iso111_distK_repeat3_hi5000_additive_offset_20260523](../K_from_continuum/results/iso111_distK_repeat3_hi5000_additive_offset_20260523)
- repeat shape-only results:
  [../K_from_continuum/results/iso111_distK_repeat3_hi5000_shape_only_20260523](../K_from_continuum/results/iso111_distK_repeat3_hi5000_shape_only_20260523)

## Bottom Line

The denser, higher-stat exact-critical follow-up and the explicit repeat
campaign did expose a real problem, but the offline rescoring pass shows that
the dominant problem sits in the original raw score definition.

Once a common-mode additive nuisance is removed, the isotropic point is
recovered as the best match in the dense run and in the repeat run under the
`additive_offset` score, while the stricter `shape_only` score leaves only a
small residual preference for `0.995` in the repeat case. That means the raw
control failure was overstating the true geometry mismatch.

The workflow is therefore closer to a valid exact-critical control than the raw
tables suggested, but the improved score still needs to be frozen and checked
against repeat floors before moving on to the `beta_c` finder.