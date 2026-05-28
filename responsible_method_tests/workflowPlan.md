# Twisted-Reference Workflow Plan

## Purpose

This document proposes a practical workflow for using twisted lattices to build
a reference two-point function that anisotropic untwisted lattice data can be
graded against.

The workflow is designed for theories such as `phi^4`, where no analytical
modular-space solution is known.

The central idea is simple:

> Use the twisted lattices themselves to define a surrogate continuum reference,
> then grade untwisted candidates against that reference with a controlled
> uncertainty budget.

This plan uses the lessons from the responsible Ising runs as design
constraints, not just as motivation.

## Executive Summary

- Do **not** use a single largest twisted lattice as the reference.
- Do **not** rely on an analytical modular solution.
- Build a **twisted-only continuum surrogate** from the largest feasible
  twisted size family.
- Make the **primary comparison three-directional**, along the three physical
  boundary-path channels of the torus.
- Use a **shape-first** score so the workflow can transfer to theories like
  `phi^4`, where the overall continuum normalization is not known a priori.
- Keep the full-cell manifold comparison as a **secondary diagnostic**, not the
  primary grade.
- Carry an explicit uncertainty model that includes fit error, window
  dependence, and residual finite-size bias.

The responsible runs strongly suggest that this is the most defensible route.
They also show what must be avoided:

- orientation mistakes can create fake disagreement
- denser twisted point sets alone do not rescue a bad comparison
- more statistics can sharpen a disagreement instead of washing it out
- the simulator ceiling prevents us from replacing extrapolation with one giant
  twisted lattice

## Recommended Reference Object

For non-analytic theories, the primary reference should be a **set of three
normalized directional shape profiles**, not the raw correlator field by
itself.

These three channels are the physically corresponding torus directions:

- primitive cycle 0
- primitive cycle 1
- the third boundary path associated with the short diagonal cycle

In repository terms, this is the same three-channel structure already returned
by `sample_directional_channels(...)` in
[K_from_continuum/workflow_common.py](K_from_continuum/workflow_common.py).

Two robust options are:

1. Anchor-normalized field:

  `R^(c)(t) = G^(c)(t) / G^(c)(t_ref)`

2. Centered log-profile:

  `L^(c)(t) = log G^(c)(t) - average_t(log G^(c)(t))`

The first is simpler and usually easier to explain. The second is more robust
when overall amplitude drift is a concern.

The twisted reference is then the twisted-only continuum extrapolation of those
three shape profiles:

`R_ref_tw^(c)(t) = continuum limit of the twisted family, for c = 0, 1, 2`

The untwisted candidate is graded against them pointwise in the directional
parameter `t`:

`z_shape^(c)(t) = [R_un^(c)(t) - R_ref_tw^(c)(t)] / sqrt(sigma_un^(c)(t)^2 + sigma_ref^(c)(t)^2)`

An amplitude comparison can still be reported as a secondary diagnostic, but it
should not be the primary grading variable.

The full periodic cell field should still be built when possible, but only as a
secondary visualization and debugging object.

## Design Constraints From Responsible Runs

The responsible Ising campaigns already tell us what a credible non-analytic
workflow must respect.

| Evidence from responsible runs | Quantitative result | Consequence for workflow |
| --- | --- | --- |
| Correct 4-5-6 side orientation matters | `4-5-6` moved to `RMS z = 0.750 / 0.679`, `rel RMS = 0.0288`, verdict `not distinguished` after fixing `embedding_cycles` | Physical-side matching is mandatory. Never compare by raw cycle index alone. |
| Correct orientation is not enough in general | `5-6-7` matched-volume run at `n_traj=2000` gave `RMS z = 1.849 / 1.960`, `rel RMS = 0.0770`, verdict `marginally separated` | A valid workflow must allow genuine disagreement after the side map is fixed. |
| Higher stats can strengthen disagreement | High-stat `5-6-7` large-window rerun at `n_traj=5000` gave `RMS z = 2.092 / 2.268`, `rel RMS = 0.1282`, verdict `distinguishable` | Do not assume residual mismatch is just noise. The reference must carry a real uncertainty budget. |
| A larger twisted point cloud is not sufficient | In the twisted-reference 4x trial, `geometry_acute_111` worsened from `tw->un RMS z = 2.236` to `3.216` | Denser twisted sampling alone is not a solution. The reference must be an extrapolated field, not just a denser finite lattice. |
| Twisted large-size trend is smooth but finite | Fixed-alpha twisted probe on `geometry_acute_567` reduced mean abs residual/site from `0.05721` to `0.04784` from scale 5 to 6 | Twisted data do carry a meaningful continuum trend, so they can support a surrogate reference. |
| The simulator has a hard ceiling | `5-6-7` scale-11 twisted probe `(132,165,66,99)` failed even at `n_traj=1` | The reference must be extrapolated from the feasible family. One giant twisted lattice is not an option. |

## Proposed Workflow

### Stage 0. Choose The Geometry And Side Map

For each benchmark geometry:

- define the target torus shape through the twisted cell
- identify the physical side ordering from the geometry record
- assign the untwisted `embedding_cycles` by physical side, not by cycle index

This is the single most important bookkeeping rule. The 4-5-6 control already
showed that the wrong side map can fake a continuum mismatch.

### Stage 1. Build A Twisted Family Only

Generate a twisted family at the largest feasible multiplicative window.

Guidelines:

- prefer the largest feasible window over the broadest window
- use explicit family entries when the feasible sizes are irregular
- for non-integrable theories, tune each size consistently to the same
  critical surface using the usual numerical criteria

The new explicit-entry support in
[scripts/generate_pointwise_manifold_dataset.py](scripts/generate_pointwise_manifold_dataset.py)
already allows this style of family specification.

### Stage 2. Define The Twisted Master Point Set And The Three Boundary Paths

Use the smallest lattice in the twisted family as the master point set.

For each point on that smallest lattice:

- map it to wrapped cell coordinates `(a, b)`
- map it to the corresponding cell coordinate `nu`
- identify the matching point on each larger twisted lattice by the same
  multiplicative scaling rule already used in the current pointwise FSS code

This preserves point identity across the twisted family.

In addition, define a fixed fraction grid `t in [0,1]` along the three
physically corresponding boundary paths.

This yields three directional channels per payload:

- `G^(0)(t)`
- `G^(1)(t)`
- `G^(2)(t)`

These three channels are the primary objects for grading.

### Stage 3. Fit The Twisted Continuum Surrogate

For the primary workflow, fit the three directional observables across twisted
sizes using the existing continuum machinery.

Recommended observable order:

1. primary: normalized directional shape field `R^(c)(t)`
2. optional secondary: raw directional connected correlator `G^(c)(t)`
3. optional tertiary diagnostic: the full pointwise field `G(nu)`

The continuum intercept at each sampled `t` value defines the twisted reference
profile for each of the three channels.

Output conceptually:

- `twisted_reference_directional.dat`
- continuum value for each sampled `t` on each channel
- uncertainty for each sampled `t` on each channel
- optional `twisted_reference_points.dat` as a secondary whole-cell object

### Stage 4. Build The Three Directional Twisted Reference Profiles

Interpolate or fit the twisted directional continuum points so that each of the
three channels is available on a common fraction grid.

This produces the primary reference objects:

`R_ref_tw^(0)(t)`, `R_ref_tw^(1)(t)`, and `R_ref_tw^(2)(t)`

If desired, also build the full periodic cell surface as a secondary object for
debugging and visualization.

The existing repository already has the right ingredients:

- `sample_directional_channels(...)` in
  [K_from_continuum/workflow_common.py](K_from_continuum/workflow_common.py)
- periodic pointwise interpolation in
  [scripts/align_manifolds.py](scripts/align_manifolds.py)
- shared-coordinate FSS logic in
  [scripts/plot_shared_coordinate_twisted_untwisted_fss.py](scripts/plot_shared_coordinate_twisted_untwisted_fss.py)

### Stage 5. Quantify Reference Uncertainty

The twisted reference uncertainty should not be just the statistical fit error.

Use:

`sigma_ref^2 = sigma_fit^2 + sigma_window^2 + sigma_bias^2 + sigma_tuning^2`

Where:

- `sigma_fit`: pointwise continuum-fit uncertainty from the twisted family
- `sigma_window`: sensitivity to the chosen size window, for example full
  feasible family versus large-window-only family
- `sigma_bias`: residual finite-size bias estimated from the gap between the
  largest twisted lattice and the twisted continuum surrogate
- `sigma_tuning`: uncertainty from critical tuning, if the theory is not at an
  exact critical point

This is the place where the responsible `5-6-7` runs matter most. They show
that the largest twisted finite lattice can still sit far enough from the
continuum surrogate to matter scientifically.

### Stage 6. Build Untwisted Candidate Families

For each anisotropic untwisted candidate:

- generate an untwisted `L x L` family at the chosen coupling point
- use the physically correct `embedding_cycles`
- fit the same observable used for the twisted reference

This produces the untwisted continuum directional profiles on the same three
physical channels.

### Stage 7. Grade Untwisted Against The Twisted Reference

This comparison should be **three-directional** in physical-cycle space.

Operationally:

- evaluate the twisted directional surrogate on the untwisted directional
  continuum points for all three channels
- compare untwisted values to the twisted reference on each channel
- summarize the three channel residuals jointly as the primary grade

Recommended primary metrics:

- directional `z_shape^(c)(t)`
- per-channel RMS z for `c = 0, 1, 2`
- aggregate directional score, for example `max_c RMS_t z_shape^(c)(t)`
- mean signed directional offset on each channel
- max absolute directional z on each channel

Recommended secondary metrics:

- best-fit amplitude factor between untwisted and twisted reference
- amplitude mismatch significance
- full-cell difference heatmap for visualization and debugging

The full-cell manifold comparison should remain available, but it should be
treated as a secondary diagnostic layer rather than the main grading rule.

### Stage 8. Rank And Decide

Use the three-directional shape score as the primary ranking variable.

Initial interpretation rule, inherited from the current responsible tests:

- `max directional RMS z <= 1`: not distinguished within current point errors
- `1 < max directional RMS z <= 2`: marginal
- `max directional RMS z > 2`: distinguishable

These thresholds are good starting points for collaborator discussions. They can
be retuned later once non-Ising benchmark data exist.

## Supporting Data Snapshot

### Benchmark-Level Summary

The table below condenses the most relevant responsible-run evidence for this
workflow choice.

| Benchmark | Size window and stats | `tw->un RMS z` | `un->tw RMS z` | common-grid rel RMS | Verdict | Main lesson |
| --- | --- | ---: | ---: | ---: | --- | --- |
| `geometry_acute_456_matched_volume` | full ladder `1..10`, `n_traj=2000` | `0.7497` | `0.6789` | `0.0288` | not distinguished | Correct orientation can make the constructions agree. |
| `geometry_acute_567_matched_volume` | full ladder `1..10`, `n_traj=2000` | `1.8487` | `1.9597` | `0.0770` | marginal | Correct orientation is necessary but not sufficient. |
| `geom_acute_567_mv_hs` | feasible large window `{2,4,6,8,10}`, `n_traj=5000` | `2.0918` | `2.2683` | `0.1282` | distinguishable | Higher stats and larger feasible sizes sharpen the disagreement. |

### Selected Shared-Coordinate Pointwise Offsets

These pointwise summaries are useful because they show whether a result is being
driven by one bad point or by a coherent field-level offset.

| Benchmark | Representative `delta_A = A_un - A_ref` | Pointwise interpretation |
| --- | --- | --- |
| `4-5-6` matched-volume | `+0.00241`, `-0.00437`, `-0.00486`, `-0.00558` | All within about `0.4 sigma`; consistent with the control passing. |
| `5-6-7` matched-volume, `n_traj=2000` | `-0.00510`, `-0.01505`, `-0.01432`, `-0.01523` | Coherent sign within the run; not just one outlier. |
| `5-6-7` high-stat large window, `n_traj=5000` | `+0.01500`, `+0.02563`, `+0.03483`, `+0.02861` | Coherent sign within the run and larger magnitude; disagreement sharpened. |

### Key Artifacts

- `4-5-6` summary:
  [results/raw_manifold_fss_acute_456_matched_volume_20260522/CONTINUUM_BENCHMARK_SUMMARY.md](results/raw_manifold_fss_acute_456_matched_volume_20260522/CONTINUUM_BENCHMARK_SUMMARY.md)
- `5-6-7` matched-volume summary:
  [results/raw_manifold_fss_acute_567_matched_volume_20260522/CONTINUUM_BENCHMARK_SUMMARY.md](results/raw_manifold_fss_acute_567_matched_volume_20260522/CONTINUUM_BENCHMARK_SUMMARY.md)
- `5-6-7` high-stat large-window summary:
  [results/raw_mfss_567_mv_hs_20260522/CONTINUUM_BENCHMARK_SUMMARY.md](results/raw_mfss_567_mv_hs_20260522/CONTINUUM_BENCHMARK_SUMMARY.md)
- `5-6-7` high-stat difference plot:
  [results/raw_mfss_567_mv_hs_20260522/geom_acute_567_mv_hs_continuum_difference.png](results/raw_mfss_567_mv_hs_20260522/geom_acute_567_mv_hs_continuum_difference.png)

## Why This Transfers To `phi^4`

This workflow does not require any exact modular correlator.

For `phi^4` or another non-integrable theory:

- the reference is built entirely from twisted Monte Carlo data
- the continuum limit is still extracted by pointwise FSS
- the critical point can be found numerically, size by size or by a shared
  tuning workflow
- the grading variable can be the three directional normalized shape profiles
  rather than the raw correlator amplitude

In other words, the only thing the Ising responsible runs contribute here is
methodological discipline and calibration, not an analytical input.

## Implementation Plan In This Repository

The most natural implementation is to add two scripts.

### 1. `scripts/build_twisted_reference.py`

Inputs:

- twisted method manifest
- choice of observable (`raw`, `anchor-normalized`, or `centered-log`)
- fitting window selection

Outputs:

- twisted directional continuum profiles on the three channels
- directional uncertainty payload
- optional periodic interpolator payload for the secondary whole-cell field
- reference uncertainty payload
- optional diagnostic plots for window dependence and largest-lattice bias

### 2. `scripts/grade_against_twisted_reference.py`

Inputs:

- twisted reference payload
- untwisted method manifest
- comparison mode settings

Outputs:

- per-channel `z(t)` tables on the three physical directions
- per-channel and aggregate directional RMS z summary
- common-grid difference heatmap as a secondary diagnostic
- amplitude-mismatch diagnostic
- one-row summary JSON or Markdown for ranking many untwisted candidates

These can reuse most of the existing interpolation and plotting code already in
[scripts/align_manifolds.py](scripts/align_manifolds.py) and
[scripts/plot_continuum_difference.py](scripts/plot_continuum_difference.py).

## Recommended First Validation Sequence

1. Rebuild the current `4-5-6` control without using any modular quantities in
   the final grade.
2. Rebuild the current `5-6-7` matched-volume cases with the same three-directional twisted-reference grade.
3. Verify that the control still passes and that the hard case still fails or
   stays borderline.
4. Only then port the workflow to `phi^4`.

This keeps the first non-analytic workflow benchmark anchored to cases we
already understand.

## Bottom Line

The proposed collaborator-facing conclusion is:

> For non-analytic theories, the most defensible reference object is a
> twisted-only continuum surrogate built from the largest feasible twisted size
> family, with explicit uncertainty from fit error, window dependence, and
> residual finite-size bias. Untwisted anisotropic candidates should then be
> graded against that twisted reference through the three physically
> corresponding directional channels, using a shape-first directional score and
> amplitude as a secondary diagnostic.

That is the workflow most consistent with everything the responsible runs have
shown so far.