# Continuum Benchmark Summary

Direct twisted-vs-untwisted comparison uses shared-cell interpolation and directional continuum-point residuals.
Interpretation rule: max directional RMS z <= 1.0 means not distinguished within current continuum point errors, 1.0-2.0 is marginal, and > 2.0 is distinguishable.

| benchmark | twisted alpha | twisted chi2/dof | untwisted alpha | untwisted chi2/dof | tw->un RMS z | un->tw RMS z | common-grid rel RMS | interpretation | comparison | metrics |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| geometry_111 (baseline) | 0.2126 | 0.0634 | 0.2097 | 0.2778 | 0.5524 | 0.5211 | 0.0231 | not distinguished within continuum point errors | n/a | n/a |
| geometry_456 (high stat) | 0.2700 | 0.1194 | 0.2689 | 0.0549 | 1.7524 | 1.0475 | 0.0084 | marginally separated | [geometry_456_twisted_vs_untwisted_common_grid.png](raw_manifold_fss_integer_multiples_high_stats_20260520/geometry_456_twisted_vs_untwisted_common_grid.png) | [geometry_456_twisted_vs_untwisted_common_grid.json](raw_manifold_fss_integer_multiples_high_stats_20260520/geometry_456_twisted_vs_untwisted_common_grid.json) |
| geometry_456 (taylor2 refit) | 0.2305 | 0.0609 | 0.2323 | 0.0625 | 0.7492 | 0.5014 | 0.0134 | not distinguished within continuum point errors | [geometry_456_twisted_vs_untwisted_common_grid.png](raw_manifold_fss_integer_multiples_high_stats_taylor2_refit_20260520/geometry_456_twisted_vs_untwisted_common_grid.png) | [geometry_456_twisted_vs_untwisted_common_grid.json](raw_manifold_fss_integer_multiples_high_stats_taylor2_refit_20260520/geometry_456_twisted_vs_untwisted_common_grid.json) |

