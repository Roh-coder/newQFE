# Mixed score candidates

Input TSV: /workspaces/newQFE/responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep_anchor_ratio/acute456_boundary_quarter_fss_anchor_ratio_small_target_zscores.tsv

All candidates are dimensionless RMS-like costs built from z-based chi^2 summaries, so they avoid mixing raw observables with incompatible units.

Definitions:

- corr_panel_mean_chi2 = corr_chi2 / 3
- ratio_panel_mean_chi2 = ratio_chi2 / 2

1. Balanced sectors: `sqrt(0.50*corr_chi2/3 + 0.50*ratio_chi2/2)`
   Equal weight to anchored-correlator and pure-ratio sectors after per-panel normalization.
2. Ratio-tilted: `sqrt(0.35*corr_chi2/3 + 0.65*ratio_chi2/2)`
   Favors direction-shape agreement while still keeping some anchored-correlator information.
3. Correlator-tilted: `sqrt(0.65*corr_chi2/3 + 0.35*ratio_chi2/2)`
   Favors the anchored-correlator sector while still penalizing directional-ratio mismatches.
4. Anchor + shape: `sqrt(0.50*u/4_chi2 + 0.25*(v/u)_chi2 + 0.25*(w/u)_chi2)`
   A nonredundant-ish amplitude-plus-shape blend: one anchored amplitude channel plus two shape ratios.
5. Sector guard: `max(sqrt(corr_chi2/3), sqrt(ratio_chi2/2))`
   A conservative guard metric that only looks good when both sectors look good.
