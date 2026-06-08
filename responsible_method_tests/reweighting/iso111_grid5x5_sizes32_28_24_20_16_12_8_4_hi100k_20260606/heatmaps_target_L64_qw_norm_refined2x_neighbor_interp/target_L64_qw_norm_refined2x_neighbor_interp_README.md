# Reweighting q_w-Normalized Refined Surface

Main raw family root: /workspaces/newQFE/responsible_method_tests/reweighting/iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606

Target raw root: /workspaces/newQFE/responsible_method_tests/reweighting/iso111_target_L64_hi100k_20260606

Target holdout: L=64, (r1,r2)=(1.000,1.000)

Base grid values: 0.900, 0.950, 1.000, 1.050, 1.100

Refined grid values: 0.900, 0.925, 0.950, 0.975, 1.000, 1.025, 1.050, 1.075, 1.100

Inserted midpoint nodes are evaluated by reweighting the enclosing base-grid corner donors to the target coupling, then blending those donor predictions with bilinear weights inside each base-grid cell. Original 5x5 nodes keep their direct scored values.

Gradient fields are rendered separately from a direct local reweight stencil at each original 5x5 node: for each arrow, the score is re-evaluated at r +/- delta e_i using that node's own donor payload with delta=0.012500. The gradient figure therefore does not finite-difference the 2x-refined neighbor-interpolated surface.

Observable set: midpoint and quarter ratios normalized by q_w = (L/4, L/4): mid_v/q_w, mid_u/q_w, mid_w/q_w, q_v/q_w, q_u/q_w.

Balance mode: none
