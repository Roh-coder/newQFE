# Modular Analogue Plan

## Goal

Build an analytical analogue of the current acute456 anchor-ratio landscape by replacing Monte Carlo correlators with the modular Ising torus solution evaluated at different modular parameters.

The immediate question is:

- does the modular solution by itself produce the same kind of clean trough we see in the anchored-ratio score?
- if it does, is that trough already present in pure modular space, or only after mapping anisotropic couplings `(r1, r2)` into an effective modular parameter?

The cleanest way to answer that is a two-stage plan.

## Existing Reusable Pieces

The repo already contains almost everything needed for the analytic side.

1. Modular parameter and modular two-point function:
   - `responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py`
   - reusable functions:
     - `_modular_tau(lx, ly, tx, ty)`
     - `_ising_torus_shape(nu, tau, theta1p0, q)`
   - this script already writes `modular_raw.dat`, `modular_aligned.dat`, and `modular_alignment.json` for a given target tau.

2. Standalone modular comparison against data:
   - `thinkDoubleTwist/compare_modular_ising_tau.py`
   - useful as a reference implementation for direct `G(nu | tau)` evaluation and amplitude fitting.

3. Existing modular-manifold diagnostics:
   - `responsible_method_tests/scripts/plot_continuum_vs_modular_manifold.py`
   - `responsible_method_tests/scripts/plot_modular_convergence.py`
   - these already consume manifests containing `target_tau`, `modular_raw`, `modular_aligned`, and alignment summaries.

4. Current boundary-anchor score shape to mirror:
   - `responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py`
   - this provides the current quarter-point observable family and the anchor convention.

5. Current heatmap style to reuse:
   - `responsible_method_tests/scripts/plot_standard_iso111_boundary_zscore_heatmaps.py`
   - `responsible_method_tests/scripts/plot_mixed_score_candidates.py`

## Recommended Scope

Mirror the current small-target acute456 quarter-point anchored-ratio story first.

That means the analytic target should be the small twisted reference geometry

- `(Lx, Ly, Tx, Ty) = (66, 66, 33, 11)`

with modular parameter obtained from `_modular_tau(66, 66, 33, 11)`.

After that works, repeat for the large target

- `(144, 144, 72, 24)`

as a robustness check.

## Stage 1: Pure Tau-Space Analogue

### Objective

Construct a deterministic score landscape in modular space

- `tau = tau_real + i tau_imag`

using only the modular solution, with no MC and no coupling-to-geometry mapping yet.

This is the fastest way to see whether the trough is intrinsic to the modular observable family.

### Target Observables

Use the same quarter-point family as the current boundary FSS plots:

- `v/4`
- `u/4`
- `w/4`
- `(v/4)/(u/4)`
- `(w/4)/(u/4)`

In modular coordinates this is naturally the fixed point set

- `nu_v = 1/4`
- `nu_u = tau/4`
- `nu_w = (1 + tau)/4`

because those correspond to `(a, b)` values

- `(1/4, 0)`
- `(0, 1/4)`
- `(1/4, 1/4)`.

### Anchor Choice

There are two viable versions.

1. Fast manifold version:
   - choose a small fixed modular offset `nu_anchor = eps * tau` or `eps`
   - for example `eps = 1/64`
   - this gives a clean analytic anchored-ratio score without reproducing the lattice-size ladder exactly.

2. Faithful standard-workflow version:
   - reproduce the exact standard anchor `(m, n) = (0, -1)` on each untwisted lattice size `L in {8,16,24,32,48,64}`
   - convert that discrete displacement to wrapped `(a, b)` using the same lattice-coordinate logic already used in `plot_standard_iso111_boundary_midpoint_fss.py`
   - evaluate the modular solution at that size-dependent anchor point
   - build the same anchor-normalized finite-size sequence the MC workflow uses.

Recommendation:

- do the fast manifold version first to expose the structural geometry of the trough
- then do the faithful finite-size-anchor version if the fast version already shows something interesting.

### Score Definition

Do not introduce artificial chi-square language at this stage, because the modular values are exact and have no statistical error bars.

Use a deterministic shape score such as

$$
S_{\mathrm{AR}}(\tau; \tau_*)
=
\sqrt{\frac{1}{N} \sum_i \left[\log R_i(\tau) - \log R_i(\tau_*)\right]^2}
$$

where

- `R_i(tau)` are the anchored modular observables
- `tau_*` is the target modular parameter.

This is the analytic analogue of the anchored-ratio score because:

- it is dimensionless
- it quotients out the overall amplitude
- it compares only relative shape.

Also record two sector summaries:

- anchored-correlator sector: `v/4`, `u/4`, `w/4`
- pure-ratio sector: `(v/4)/(u/4)`, `(w/4)/(u/4)`

so the analytic trough direction can be compared directly to the MC `corr_chi2` and `ratio_chi2` panels.

### Heatmaps To Produce

For a rectangular tau grid centered on the small-target tau, write:

1. `corr_score(tau)`
2. `ratio_score(tau)`
3. `total_score(tau)`
4. optional signed residual panels for the five observables

Suggested first grid:

- `tau_real`: target real part plus/minus about `0.08`
- `tau_imag`: target imaginary part plus/minus about `0.20`
- grid size: `61 x 61`

Use the same visual conventions as the current heatmaps:

- star = best point for that panel
- x = target tau

### Minimal Implementation Path

Add a new script, for example:

- `responsible_method_tests/scripts/plot_modular_analogue_tau_landscape.py`

Core steps:

1. import `_modular_tau` and `_ising_torus_shape` from `generate_pointwise_manifold_dataset.py`
2. compute the target tau from the small target lattice `(66,66,33,11)`
3. evaluate the target modular observables once
4. loop over candidate tau values on a grid
5. compute anchored observables and deterministic residual scores
6. write a TSV plus heatmaps

Suggested outputs:

- `results/small_target_tau_landscape.tsv`
- `results/small_target_tau_landscape.png`

inside this folder or a sibling results directory.

## Stage 2: Coupling-Mapped Analogue In `(r1, r2)`

### Objective

Map the candidate couplings `(r1, r2)` to an effective modular parameter `tau_eff(r1, r2)` and then evaluate the Stage 1 analytic score on the same `5 x 5` or `9 x 9` coupling grid used by the MC campaign.

This is what will let the analytic trough be compared directly to the observed anchored-ratio trough in coupling space.

### Required Geometry Map

This stage needs the critical anisotropy-to-geometry relation for the triangular Ising model.

The intended mapping is:

1. for fixed critical couplings `(k1, k2, k3)` with `k3 = 1`, solve for the equivalent continuum triangle angles `theta_i`
2. construct the continuum metric / basis vectors associated with that triangle
3. from that basis, compute the effective torus modular parameter `tau_eff`

The right relation to use is the standard critical anisotropy-angle correspondence,

$$
\sinh(2 \beta_c k_i) = \cot(\theta_i),
\qquad
\theta_1 + \theta_2 + \theta_3 = \pi.
$$

This stage is harder because it is not the same as simply reading `target_tau` from an existing manifest. It needs a new explicit `couplings -> tau_eff` helper.

### Why This Stage Matters

If the current MC trough is really geometry-driven, then a deterministic landscape built from

- `tau_eff(r1, r2)`
- the modular solution at that tau

should show a similar valley orientation in `(r1, r2)`.

If it does not, then the observed trough is likely not explained by modular-shape geometry alone.

### Minimal Implementation Path

Add a second script, for example:

- `responsible_method_tests/scripts/plot_modular_analogue_rgrid.py`

Core steps:

1. implement `tau_eff_from_couplings(r1, r2, beta_c)`
2. evaluate `S_AR(tau_eff(r1, r2); tau_target)` on the same coupling grid as the MC scan
3. write a TSV with columns
   - `r1`
   - `r2`
   - `tau_eff_real`
   - `tau_eff_imag`
   - `corr_score`
   - `ratio_score`
   - `total_score`
4. render heatmaps with the same conventions as the current boundary plots

Recommended first target grid:

- the existing small-target acute456 `5 x 5`
- then rerun on the `9 x 9` backfill grid when available.

## Stage 3: Independent-Crossing Test

Once Stage 1 and Stage 2 exist, test the dream scenario directly.

Desired outcome:

- the MC anchor-ratio surface provides one trough
- a second, more independent analytic score provides a second trough
- the crossing lies near the true couplings.

The most promising analytic candidates are:

1. pure anchored modular ratios at a different anchor scale
2. absolute modular amplitude after single-parameter alignment
3. modular local-slope or directional-derivative observables at the same sampled points

The key rule is that the second score must not just be a trivial reweighting of the same anchored quarter-point observables, or it will inherit the same trough instead of cutting across it.

## Success Criteria

### Stage 1 succeeds if

- the code writes a stable tau-space heatmap
- the target tau is a local minimum of the deterministic anchored-ratio score
- the trough direction can be identified visually and numerically.

### Stage 2 succeeds if

- the analytic `(r1, r2)` trough direction can be compared directly to the MC trough direction
- the truth lies on or near the analytic trough
- the result is robust between small and large target tau.

### Stage 3 succeeds if

- a second analytic observable family produces a trough that is not parallel to the first one
- the intersection region is close to the true acute456 couplings.

## Recommended Execution Order

1. Implement Stage 1 fast manifold version around the small target tau.
2. Confirm whether the pure modular anchored-ratio score already shows a diagonal trough.
3. Add the faithful finite-size-anchor version only if the fast version is promising.
4. Implement `tau_eff(r1, r2)` and map the analytic score back onto the coupling grid.
5. Only then design a second analytic observable family intended to cross the first trough.

## Expected Interpretation

What this plan should tell us is not just whether the modular solution fits the data, but whether the very existence and direction of the current anchored-ratio trough are already encoded by modular geometry itself.

That is the main analytical analogue worth testing.