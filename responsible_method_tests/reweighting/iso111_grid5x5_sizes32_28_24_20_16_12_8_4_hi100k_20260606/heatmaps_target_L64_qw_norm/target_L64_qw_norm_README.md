# Reweighting q_w-Normalized Heatmaps

Main raw family root: /workspaces/newQFE/responsible_method_tests/reweighting/iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606

Target raw root: /workspaces/newQFE/responsible_method_tests/reweighting/iso111_target_L64_hi100k_20260606

Target holdout: L=64, (r1,r2)=(1.000,1.000)

Observable set:
- midpoint ratios normalized by q_w: mid_v/q_w, mid_u/q_w, mid_w/q_w
- quarter ratios normalized by q_w: q_v/q_w, q_u/q_w

Each grid cell is scored by fitting the raw finite-size ladder with A + B x^omega at x = 1/L, then comparing the predicted value at the target x = 1/64 to the direct L64 target using signed log residuals and chi2 = residual^2.

Important caveat: the truly unused point (L/4, 3L/4) is not present in the current saved bundles. These heatmaps therefore use q_w = (L/4, L/4) as the currently available interior normalization point.

Balance mode: none
