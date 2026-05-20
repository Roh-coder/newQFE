# Benchmark Matrix

This file turns the high-level README into an actionable benchmark list for the
first collaborator-facing validation campaign.

## Benchmarks

| ID | Name | Matched constructions | Nominal size family | Primary question |
| --- | --- | --- | --- | --- |
| B0 | Equilateral control | modular vs twisted `(L,L,0,0)` vs untwisted `(r1,r2)=(1,1)` | `L = 16, 32, 48, 64, 80` | Does the full pipeline agree at the symmetric control point? |
| B1 | Quarter-twist benchmark | modular vs twisted `(L,L,L/4,L/4)` vs tuned untwisted anisotropy | `L = 16, 32, 48, 64, 80, 96` | Does a nontrivial but moderate torus shape match across all three routes? |
| B2 | 4-5-6 stress benchmark | modular vs twisted `alpha*(13,16,-3,3)` vs untwisted `(r1,r2) ~= (5.0652, 7.7429)` | `alpha = 1, 2, 3, 4` for twisted; untwisted `L = 16, 24, 32, 48, 64` | Does the equivalence survive a strongly anisotropic shape? |

## Observables To Compare

- Whole-manifold interpolated RMS mismatch on a common normalized cell.
- Whole-manifold uncertainty-weighted mismatch where bootstrap or jackknife
  errors are available.
- Line-cut mismatch along the three physically corresponding short cycles.
- Shape-only mismatch after removing one fitted global amplitude.
- Continuum-fit residuals for each route separately.

## Local Deformation Stencils

Each benchmark should include a local wrong-answer stencil so that the score is
tested for separation power rather than only for zero-at-truth behavior.

### B0 deformations

- Coupling offsets: `(0.95, 1.00)`, `(1.05, 1.00)`, `(1.00, 0.95)`, `(1.00, 1.05)`.
- Twist offsets: `(L,L,L/8,0)` and `(L,L,L/8,L/8)`.

### B1 deformations

- Coupling offsets around the tuned untwisted point: start with `+-2%` and
  `+-5%` in each ratio separately.
- Twist offsets around the quarter-twist point: `(Tx,Ty) = (L/4 +- 1, L/4)` and
  `(L/4, L/4 +- 1)` on each lattice size.

### B2 deformations

- Coupling offsets around `(5.0652, 7.7429)`: start with multiplicative
  factors `0.97`, `1.03`, `0.95`, `1.05` applied one ratio at a time.
- Twist offsets around `alpha*(13,16,-3,3)`: test `Ty -> Ty +- 1` at fixed
  `Tx`, and `Tx -> Tx +- 1` at fixed `Ty`, then rescale by `alpha`.

## Exit Criteria

Each benchmark is considered passed only if all of the following are true.

1. The matched point is statistically consistent with zero mismatch on the
   primary whole-manifold score, with `|Z_match| <= 2`.
2. The matched point is the lowest-score point or statistically tied for the
   lowest within its local deformation stencil.
3. At least one wrong point in the coupling stencil and one wrong point in the
   twist stencil satisfy `Z_sep >= 3` relative to the matched point.
4. The pass/fail classification is unchanged when the smallest size is dropped.
5. For B1 and B2, the physically paired-cycle comparison beats naive
   index-paired matching.

## Deliverables Per Benchmark

- One summary note describing the matched geometry and the deformation stencil.
- One plot showing the aligned manifolds and their difference surface.
- One plot showing the three matched physical line cuts.
- One table of continuum scores with uncertainties.
- One table of deformation scores and `Z_sep` values.