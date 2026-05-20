# Continuum-Beta Test Plan

## Goal

Test a common-beta workflow for the finite test lattices:

1. For each coupling point `(k1/k3, k2/k3)`, measure the finite-size pseudo-critical sequence `beta_c(L)` with the existing beta finder.
2. Fit that sequence with a free-power extrapolation
   `beta_c(L) = beta_c(infinity) + a * L^(-p)`.
3. Re-run production MC for every requested test size using the shared `beta_c(infinity)` instead of the per-size `beta_c(L)`.
4. Keep Step 3 scoring unchanged: sample directed correlators at every eighth, fit the continuum intercept of the two-point channels across test sizes, then compare those continuum channels to the single large-reference lattice.

## Implementation plan

### Step 2 changes

- Split the current `find beta -> production` path into reusable pieces:
  - finite-size beta finding for one lattice,
  - production MC at an explicitly supplied beta.
- Add a new Step 2 mode:
  - legacy: `per_size_beta`
  - new: `free_power_continuum`
- In `free_power_continuum` mode, process jobs by coupling point instead of by `(L, r1, r2)` payload:
  - gather `beta_c(L)` over all requested sizes,
  - fit one shared `beta_c(infinity)`,
  - run production at that shared beta for each size.
- Preserve the existing manifest contract used by Step 3.

### Data products

- Keep the existing per-payload `.dat` and `.meta.json` outputs.
- Record both betas in metadata:
  - `finite_size_beta_c`
  - `production_beta`
- Write a test-side extrapolation summary file so each coupling point has an auditable `beta_c(infinity)` fit.

### Validation

- Run a smoke-sized Step 2 job in `free_power_continuum` mode.
- Run Step 3 scoring against the existing smoke reference.
- Confirm:
  - Step 2 manifest remains consumable by Step 3.
  - Test metadata now distinguishes finite-size and production betas.
  - Scoring completes without code changes in Step 3.

## Implementation status

- Implemented in `workflow_common.py` and `02_generate_tests.py`.
- Added config surface: `beta_strategy.mode = "free_power_continuum"`.
- Added smoke and example configs for the new mode.
- Smoke validation completed:
  - Step 2: `results/smoke_tests_free_power_continuum_20260515/`
  - Step 3: `results/smoke_compare_free_power_continuum_20260515/`

## Expected strengths

- Removes per-size beta-finder jitter from the production two-point data.
- Should smooth FSS trends in the directed correlators if the roughness is mainly caused by noisy `beta_c(L)`.
- Keeps the large-reference side unchanged, so the experiment isolates the test-side beta prescription.
- Cheap compared with new observables or a full crossing-based rewrite.

## Expected weaknesses

- It is not true pseudo-critical sampling anymore; each finite lattice is run at a shared extrapolated beta, not at its own `beta_c(L)`.
- The quality of the method depends on the stability of the free-power extrapolation of `beta_c(L)`.
- Small-size points can still bias `beta_c(infinity)` if the finite-size beta sequence is noisy or non-asymptotic.
- This is still a diagnostic / workflow experiment, not a principled replacement for RG-invariant crossing analyses.