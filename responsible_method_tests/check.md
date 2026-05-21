# Quarter-Twist Follow-Up Check

## Goal

Test whether the quarter-twist inverse-matching signal sharpens when the unstable slice is rerun at higher Monte Carlo statistics.

## Working Hypothesis

The current quarter-twist ambiguity is dominated by run-to-run Monte Carlo and continuum-fit noise, not by a genuine degeneracy between nearby untwisted couplings.

## Tests To Run

1. Generate one fresh higher-stat twisted target repeat for the quarter-twist geometry at the exact isotropic critical point.
2. Generate two independent higher-stat untwisted repeats at the current best coupling candidate `(r1, r2) = (1.4, 1.0)`.
3. Generate one fresh higher-stat untwisted neighbor at `(r1, r2) = (1.4, 0.8)`.
4. Rank the three untwisted candidates against the fresh twisted target with `fit_target_geometry.py`.
5. Compare the two same-coupling untwisted repeats directly on the shared modular cell.

## Pass / Fail Criteria

- Pass: the two `(1.4, 1.0)` repeats are no longer distinguished within current continuum errors and both score better against the twisted target than `(1.4, 0.8)`.
- Mixed: `(1.4, 1.0)` still ranks best, but same-coupling repeats remain distinguishable.
- Fail: the ranking of nearby couplings is unstable and same-coupling repeats remain clearly separated.

## Run Status

- Status: completed
- Config: `responsible_method_tests/configs/raw_manifold_fss_quarter_twist_followup_high_stats.json`
- Run root: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/`

## Observed Results

### 1. Target-vs-candidate ranking

Fresh higher-stat twisted target versus fresh higher-stat untwisted candidates:

| rank | candidate | r1 | r2 | score | tw->un RMS z | un->tw RMS z | rel RMS | interpretation |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | qtwist_repeatB_r140_r100_hs | 1.4 | 1.0 | 1.3014 | 1.2595 | 1.3014 | 0.0184 | marginally separated |
| 2 | qtwist_repeatA_r140_r100_hs | 1.4 | 1.0 | 1.4141 | 1.4141 | 1.3165 | 0.0194 | marginally separated |
| 3 | qtwist_repeat_r140_r080_hs | 1.4 | 0.8 | 3.5254 | 3.5254 | 3.3703 | 0.0585 | distinguishable |

### 2. Same-coupling repeat stability

Direct comparison between the two fresh `(1.4, 1.0)` untwisted repeats:

- common-grid relative RMS = `0.00345`
- repeatA->repeatB RMS z = `0.2678`
- repeatB->repeatA RMS z = `0.2735`
- interpretation = `not distinguished within continuum point errors`

## Verdict

Pass on the immediate follow-up test: with `n_traj = 10000`, the same-coupling untwisted repeats are stable against each other and both beat the nearby `(1.4, 0.8)` candidate by a wide margin. The inverse signal has sharpened enough to separate the nearest tested neighbor.

The remaining limitation is that the fresh twisted target versus fresh best untwisted candidates is still only marginal rather than fully collapsed, with best scores around `1.30` to `1.41`. That means the local map is now usable at coarse resolution near `(1.4, 1.0)`, but it is not yet at the point where sub-neighbor localization should be trusted.

## Artifacts

- Fit scan markdown: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/followup_fit_scan/geometry_quarter_twist_target_hs_coupling_fit_scan.md`
- Fit scan JSON: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/followup_fit_scan/geometry_quarter_twist_target_hs_coupling_fit_scan.json`
- Fit scan plot: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/followup_fit_scan/geometry_quarter_twist_target_hs_coupling_fit_scan.png`
- Twisted target method manifest: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/geometry_quarter_twist_target_hs/twisted/manifest_geometry_quarter_twist_target_hs_twisted.json`
- Untwisted repeat A manifest: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/qtwist_repeatA_r140_r100_hs/untwisted/manifest_qtwist_repeatA_r140_r100_hs_untwisted.json`
- Untwisted repeat B manifest: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/qtwist_repeatB_r140_r100_hs/untwisted/manifest_qtwist_repeatB_r140_r100_hs_untwisted.json`
- Untwisted neighbor manifest: `responsible_method_tests/results/raw_manifold_fss_quarter_twist_followup_high_stats_20260521/qtwist_repeat_r140_r080_hs/untwisted/manifest_qtwist_repeat_r140_r080_hs_untwisted.json`