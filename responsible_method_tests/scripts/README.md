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
- `score_campaign.py`: compute whole-manifold and cycle-cut scores for one
  benchmark entry from `configs/campaign_template.json`.
- `summarize_benchmark.py`: emit a compact table of matched and deformed scores
  plus pass/fail status.

## Scope Rule

These scripts should consume completed manifests and data products. They should
not duplicate Monte Carlo production logic that already lives elsewhere in the
repository.