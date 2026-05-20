# K_from_continuum

A self-contained continuum workflow for coupling-to-geometry scoring.

This directory is the symmetric version of the large-reference workflow:
- Step 1 builds a reference size family at fixed reference couplings.
- Step 1 fits a shared continuum beta from the finite-size pseudo-critical sequence.
- Step 1 reruns production at that shared beta and extrapolates the reference channels to the continuum.
- Step 2 does the same shared-beta production workflow for every test coupling point.
- Step 3 fits each test channel to the continuum and compares test continuum channels directly to the reference continuum channels.

The goal is to reduce score-landscape roughness from per-size beta jitter and to avoid the asymmetry of comparing continuum-extrapolated test data against a single finite reference lattice.

## Layout

```text
K_from_continuum/
  01_generate_reference.py
  02_generate_tests.py
  03_score_continuum_vs_reference.py
  workflow_common.py
  Makefile
  bin/
  include/
  lib/
  src/
  configs/
  results/
```

## Quick start

Run from this directory:

```bash
cd /path/to/K_from_continuum
python -m pip install -r requirements.txt

python 01_generate_reference.py --config configs/reference_example.json
python 02_generate_tests.py --config configs/tests_example.json
python 03_score_continuum_vs_reference.py --config configs/score_example.json
```

For a fast mechanics check:

```bash
python 01_generate_reference.py --config configs/smoke_reference.json
python 02_generate_tests.py --config configs/smoke_tests.json
python 03_score_continuum_vs_reference.py --config configs/smoke_score.json
```

## Config summary

Step 1 uses:
- `reference_family.sizes`
- `reference_family.geometry_defaults` or `reference_family.geometry_map`
- `reference_couplings`
- `beta_strategy` (must be `free_power_continuum`)
- `beta_finder`
- `mc`
- `analysis.k_values`
- `analysis.k_denominator`
- `analysis.power_fit`

Main Step 1 outputs:
- `reference_data/grid/L*/reference_*.dat`
- `reference_data/continuum_beta_extrapolation.json`
- `reference_data/reference_raw_channels.dat`
- `reference_data/reference_continuum_channels.dat`
- `reference_data/reference_fss_data.dat`
- `reference_data/reference_fss.png`
- `reference_data/manifest_reference.json`

Step 2 supports three beta modes:
- `per_size_beta`
- `free_power_continuum`
- `taylor2_continuum`

For small test-size families, `taylor2_continuum` is often more stable than the free-power fit.

Main Step 2 outputs:
- `test_data/grid/L*/test_*.dat`
- `test_data/continuum_beta_extrapolations.json`
- `test_data/manifest_tests.json`

Step 3 uses:
- `input.reference_manifest`
- `input.tests_manifest`
- `analysis.k_values`
- `analysis.k_denominator`
- `analysis.power_fit`
- `analysis.fss_plots`
- `analysis.weighted`

Main Step 3 outputs:
- `reference_data/reference_channels_scored.dat`
- `test_data/raw_test_channels_scored.dat`
- `comparison_analysis_data/continuum_test_channels.dat`
- `comparison_analysis_data/channel_comparison.dat`
- `comparison_analysis_data/score_map.dat`
- `comparison_analysis_data/score_minimum.dat`
- `comparison_analysis_data/fss_plot_data.dat`
- `comparison_analysis_data/fss_plot_index.dat`
- `comparison_analysis_data/fss_plots/*.png`
- `comparison_analysis_data/score_heatmap.png`
- `comparison_analysis_data/zscore_heatmap.png`
- `comparison_analysis_data/score_zscore_heatmaps.png`
- `comparison_analysis_data/manifest_score.json`

## Notes

- Channel fits use the free-power form `A + B * (1/L)^C` with guardrails.
- `analysis.power_fit.min_sizes_for_free_C` controls when the scorer is allowed to fit the exponent `C` freely.
- With too few sizes, or for obviously unstable channels, the code falls back to fixed `C = 1` rather than returning pathological continuum intercepts.
- Relative config paths are resolved from this directory.

## Session notes

- 2026-05-18 Windows bundle debug: `dist/newQFE_K_from_continuum_L128_local.zip` is a reduced standalone package, not a full mirror of this directory. During debugging, the bundle diverged from the live repo and its packaged `02_generate_tests.py` had gone stale relative to `run_from_spyder.py` and `configs/local_L128_tests.json`.
- The Windows/Spyder reproduction used the bundle-local `configs/local_L128_reference.json` with isotropic couplings `k1=k2=k3=1`, sizes `[16, 32, 48, 64, 80, 96, 112, 128]`, `beta_finder` trajectories `(200 coarse, 400 fine)`, and production `mc.n_traj=2000`.
- The immediate Step 2 failure after a successful Step 1 run was a config-schema mismatch, not an MC/runtime problem: the stale bundled `02_generate_tests.py` was still expecting `reference_family`/`reference_couplings`/`analysis`, while the bundle test config correctly provided `test_family`/`couplings` and the newer Step 2 workflow. The packaged Step 2 driver under `dist/newQFE_K_from_continuum_L128_local/K_from_continuum/02_generate_tests.py` was resynced to the current repo version.
- Step 1 can look stuck even after the C++ simulator has finished. In the Spyder wrapper, `run_from_spyder.py` keeps printing the last parsed status snapshot every few seconds when the child process is quiet, so repeated lines such as `phase=production MC` may be stale status rather than evidence that the simulator is still running.
- The quiet tail after the final `run_simulator done` line is likely Python post-processing in `01_generate_reference.py`, not Monte Carlo. Step 1 reloads each `two_point_all_to_all.dat`, then repeatedly calls `sample_directional_channels(...)`, which rebuilds tiled `LinearNDInterpolator` objects for every payload during the channel/FSS fit pass.
- Practical resume rules: if the last visible simulator line is `run_simulator done` for the smallest reference lattice, wait for the analysis tail before killing the job; if Step 2 fails immediately on missing top-level sections, check that the bundled `02_generate_tests.py` still matches the test-family/couplings schema before rerunning. If this needs more debugging later, instrument the channel-analysis loop in `01_generate_reference.py` and/or add an explicit post-MC phase update in `run_from_spyder.py`.