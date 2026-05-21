# Scripts

This folder is reserved for the lightweight comparison utilities that sit on
top of existing data-generation workflows.

## Intended First Scripts

- `build_modular_manifold.py`: sample the analytic modular correlator on the
  same normalized cell used by the lattice comparisons.
- `align_manifolds.py`: load modular, twisted, and untwisted all-to-all data;
  map them into common cell coordinates; and build interpolated comparison
  grids.
- `plot_pointwise_fss.py`: plot representative pointwise finite-size-scaling
  panels from a generated method manifest, including the continuum intercept
  and aligned modular reference.
- `plot_fss_interpolated_manifolds.py`: render non-extrapolated, periodically
  interpolated correlator manifolds for each finite-size scale as transparent
  overlays on the common target torus cell.
- `plot_continuum_vs_modular_manifold.py`: overlay the continuum-extrapolated
  twisted/untwisted manifolds with the aligned modular target manifold for a
  completed benchmark.
- `plot_modular_convergence.py`: summarize, per benchmark, how twisted and
  untwisted finite-size manifolds approach their shared modular target using
  per-scale amplitude-fitted residual metrics and largest-scale comparisons.
- `plot_residual_scaling.py`: show only the per-scale residual trends versus
  lattice-size scale factor for twisted and untwisted families.
- `align_manifolds.py`: reinterpolate the twisted and untwisted continuum
  manifolds onto the same target torus cell, emit direct cross-method
  agreement metrics, render a common-grid comparison figure, and optionally
  switch to a `twisted_reference` mode where only the twisted interpolant
  sampled on untwisted continuum points controls the distinguishability
  verdict.
- `fit_target_geometry.py`: rank untwisted coupling candidates against a
  chosen twisted target torus using the same common-cell continuum score,
  and emit a localization table plus a scan plot.
- `score_campaign.py`: compute whole-manifold and cycle-cut scores for one
  benchmark entry from `configs/campaign_template.json`.
- `summarize_benchmark.py`: emit a compact table of matched and deformed scores
  plus pass/fail status.

The currently implemented summary workflow is:

- `align_manifolds.py` for one benchmark manifest.
- `align_manifolds.py --comparison-mode twisted_reference` when the intended
  workflow is a dense twisted reference manifold compared against a sparser
  untwisted candidate manifold.
- `fit_target_geometry.py` for a local untwisted coupling scan around one
  target torus.
- `summarize_benchmark.py` for a small curated list of benchmark manifests,
  using the same `--comparison-mode` as the align step.

## Scope Rule

These scripts should consume completed manifests and data products. They should
not duplicate Monte Carlo production logic that already lives elsewhere in the
repository.