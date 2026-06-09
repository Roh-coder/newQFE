# Journal

## 2026-06-08: Proposal For Reweight-Augmented CMA-ES With Local Curvature

### Question

Could histogram reweighting be used to accelerate CMA-ES by attaching a local gradient or Hessian estimate to each evaluated point?

### Short Answer

Yes, but only as a local augmentation.

CMA-ES already learns an inverse-Hessian-like covariance from ranked population history, so feeding every candidate's raw gradient or Hessian directly into the core update is probably the wrong first design. The more promising use of reweighting is to provide cheap local curvature information around a small number of directly simulated points, especially the generation mean or the top `1-3` elites.

### Why This Could Help

Because a direct MC point already stores the per-trajectory decomposition needed for nearby reweighting, one expensive direct evaluation can support several cheap local score probes around the same donor point.

In the current `2d` control problem, a local quadratic model is affordable:

- center `f(r1, r2)`
- axial probes `f(x +/- delta e1)`, `f(x +/- delta e2)` for the gradient and diagonal curvature
- corner probes `f(x +/- delta e1 +/- delta e2)` for the mixed curvature

That would yield:

- local gradient `g`
- local Hessian `H`
- jackknife error bars on both, if the current block-replica machinery is extended from gradients to Hessian entries

If `g` and `H` are reliable, they can tell CMA-ES which directions are stiff, flat, or badly correlated before the rank-based covariance update has enough generations to learn that structure itself.

### Why Not Give Every Candidate Its Own Gradient And Hessian

That would likely help less than it sounds.

- CMA-ES is robust partly because it updates from ranks rather than trusting noisy derivatives.
- A per-candidate Hessian is only local to that candidate; across a spread-out population these local models can disagree strongly.
- If every individual contributes both a function value and a noisy curvature tensor, the algorithm stops behaving like plain CMA-ES but does not yet become a clean trust-region Newton method either.
- The expensive part is still the direct MC point. Reweighting is cheap relative to that, but not free enough to justify a full Hessian build for the entire population unless the population is tiny.
- Bad Hessians are dangerous: negative curvature, near-singularity, and mixed-observable compromise directions can produce overconfident steps.

So the right comparison is not

- plain CMA-ES
- versus CMA-ES where every point has a Hessian

but rather

- plain staged-fidelity CMA-ES
- versus CMA-ES plus a small number of curvature-guided local proposals.

### Proposed Hybrid

1. Keep direct `total_score` as the only quantity that decides selection and the CMA-ES rank update.
2. At each generation, choose a small curvature set: the current mean and the top `1-2` elites after cheap direct evaluation.
3. For each curvature point, use local histogram reweighting to estimate `g`, `H`, and jackknife uncertainties.
4. Regularize `H` aggressively:
	- symmetrize it
	- floor eigenvalues
	- reject or downweight it when curvature SNR is poor or the condition number is extreme
5. Use the regularized inverse Hessian only to augment proposal generation, not to replace CMA-ES:
	- initialize or warm-start the local sampling covariance near the mean
	- add one or two trust-region Newton proposals `x_newton = x - alpha (H + lambda I)^-1 g`
	- whiten mutation directions using `(H + lambda I)^-1` inside a bounded local radius
6. Evaluate those extra proposals with direct MC and let normal CMA-ES ranking decide whether they survive.
7. Fall back to plain staged-fidelity CMA-ES whenever curvature quality is poor.

### Expected Benefit

This should reduce the number of generations needed to learn the dominant anisotropy inside a basin.

It is most likely to help:

- after CMA-ES has already found the right basin
- when the local score surface is smooth enough for nearby reweighting to be trustworthy
- when the covariance learned from the population is still lagging behind the local ridge geometry

It is less likely to help:

- very far from the basin
- when the objective is still multimodal or strongly non-quadratic over the current population radius
- when the cheap-budget MC noise is already large enough to destabilize the current gradient estimates

### Minimal Experiment

The first version should be narrow.

1. Run staged-fidelity CMA-ES exactly as already proposed.
2. At the generation mean only, build a local reweight quadratic model from one direct payload.
3. Spawn one extra trust-region Newton candidate and one extra Hessian-whitened random candidate.
4. Score those two candidates with direct MC and include them in the same generation ranking.
5. Compare against baseline CMA-ES at matched direct-MC cost:
	- best score versus direct evaluations
	- time to reach the current basin floor
	- stability of repeated runs with different seeds
	- agreement between predicted and direct improvement for the curvature-guided candidates

### Recommendation

So the answer is:

- gradients would probably help CMA-ES find the minimum faster if they are used to bias a few local proposals
- Hessians could help even more by exposing the local ridge geometry
- but a full "every point gets a gradient and Hessian" design is probably too noisy and too complicated for the first pass
- the right first experiment is a curvature-augmented CMA-ES where reweighting informs elite local proposals while CMA-ES keeps the global, rank-based search logic

The good news is that this does not require inventing a new data path. The current local reweight score probes already exist, and the current jackknife gradient machinery already shows how to attach uncertainty estimates to local derivative objects.

## 2026-06-08: Launched 50-Step Far-Start Live-Covariance Runs And Sketched Two Next Optimizer Variants

### 50-Step Runs Launched

After validating the live-updating covariance-aware trajectory plots on the `5`-step farther-start runs from

- `(r1, r2) = (0.5, 5.0)`
- `n_traj = 10000`
- `n_therm = 1000`
- fit sizes `{32, 28, 24, 20, 16, 12, 8, 4}`

I launched matching `50`-step continuations under fresh roots:

- `responsible_method_tests/reweighting/iso111_gradient_flow_10k_start_r0p500_r5p000_scalar_snr_steps50_livecov_20260608/`
- `responsible_method_tests/reweighting/iso111_gradient_flow_10k_start_r0p500_r5p000_covariance_steps50_livecov_20260608/`

Both roots were preseeded from the completed `5`-step live-cov runs by copying

- `_target/`
- `raw/`

from the corresponding `steps5_livecov` roots. So the first five steps should be cache hits and the new work should begin at the first not-yet-sampled proposal beyond those earlier trajectories.

The driver now rewrites

- `trajectory.tsv`
- `trajectory.png`

after every attempted step, and the live figure includes

- the parameter-space path
- the score history
- covariance components by step
- the latest `2 x 2` gradient-covariance heatmap

which makes the long run much easier to inspect while it is still active.

### 50-Step Outcomes

Both `50`-step runs are now complete.

#### Scalar-SNR (`steps50_livecov`)

Root:

- `responsible_method_tests/reweighting/iso111_gradient_flow_10k_start_r0p500_r5p000_scalar_snr_steps50_livecov_20260608/`

Final summary:

- final point `(0.5934887434906397, 4.968642287164648)`
- final score `0.03850576492914759`
- accepted / attempted = `3 / 4`
- stop reason = `gradient_snr_below_threshold`

Important observation: giving scalar-SNR `50` allowed steps did **not** change the outcome relative to the earlier `5`-step run. It hit the same SNR wall at the same point. So in this regime the controlling limitation is the cheap-budget gradient noise floor, not `num_steps`.

#### Covariance (`steps50_livecov`)

Root:

- `responsible_method_tests/reweighting/iso111_gradient_flow_10k_start_r0p500_r5p000_covariance_steps50_livecov_20260608/`

Final summary:

- final point `(0.5745067081108214, 4.605038046221583)`
- final score `0.034821679542951706`
- accepted / attempted = `23 / 24`
- stop reason = `gradient_snr_below_threshold`

This is materially different from scalar-SNR. At the same `10k / 1k` budget and the same farther start `(0.5, 5.0)`, the covariance-preconditioned descent direction remained usable for a long time, continuing through `23` accepted steps before finally hitting the SNR gate.

#### Comparison / Continuation Context

The key result from this `50`-step test is:

- scalar-SNR stops early because the local gradient estimate becomes noise-dominated almost immediately
- covariance continues to make progress for much longer at the same budget
- covariance also reaches the lower final score (`0.03482` vs `0.03851`)

So the live conclusion to carry forward is not simply “more steps help.” It is narrower and more important:

- increasing `num_steps` alone is useless once the gradient SNR gate is the limiting factor
- covariance preconditioning substantially delays that failure mode
- the next useful improvement is therefore **adaptive budget escalation**, not just larger step budgets

If we resume later, the first place to look is the covariance `steps50_livecov` root, especially

- `trajectory.tsv`
- `trajectory.png`
- `trajectory_path.png` after rerendering if needed

because that run contains the longest successful cheap-budget descent we have so far from the farther start.

### Potential Plan A: CMA-ES On The Direct-MC Objective

The clean version of a CMA-ES follow-up would be:

1. Optimize directly in `(r1, r2)` using the same `total_score` objective now used by gradient flow.
2. Treat each candidate evaluation as a full direct MC point, reusing the existing cached `raw/L*/.../selected_observables_bundle.dat` outputs.
3. Start with a deliberately cheap evaluation budget for the whole population, then selectively reevaluate elites at higher stats before updating the CMA distribution.

The practical design I would propose is:

- state: `x = (r1, r2)` in the existing bounded box
- initial mean: either the user start point or the best point from a short gradient pre-run
- initial sigma: chosen from the desired exploration radius, not from the gradient covariance
- objective: direct `total_score`
- cache policy: exact reuse by `(r1, r2, sizes, n_traj, n_therm, n_skip)`

To keep CMA-ES sane under Monte Carlo noise, it should not update from single noisy scores alone. The minimum viable noise handling would be:

- evaluate the whole population at a cheap baseline budget
- reevaluate the top `k` candidates at the same budget with independent seeds to estimate rank stability
- optionally promote the best `1-3` candidates to a higher budget before the generation update

That creates a staged-fidelity CMA-ES rather than a naive one-shot noisy CMA-ES.

The main reason this is attractive is that CMA-ES can still move when the local gradient estimate becomes too noisy or too anisotropic to trust. The main reason to be careful is cost: a population-based method can burn through many direct points quickly if we do not aggressively cache and stage the fidelity.

### Potential Plan B: Gradient Flow With Adaptive Budget Escalation

The second follow-up should keep the current local-reweight gradient logic but adapt the Monte Carlo budget as the walker approaches the floor and the signal-to-noise ratio degrades.

The key idea is:

- far from the minimum, use cheap MC because the gradient is large and the direction is easy
- near the minimum, increase trajectories and/or the lattice ladder because the gradient norm shrinks and the acceptance decision becomes noise-limited

The clean control variables are:

- `grad_norm_snr`
- `gradient_vector_snr`
- gradient covariance anisotropy / correlation
- disagreement between predicted and direct proposal scores
- actual accepted step norm

A practical tiered scheme would be:

1. Start with a cheap budget, for example `n_traj = 10000`, `n_therm = 1000`, fit sizes `{32, 28, 24, 20, 16, 12, 8, 4}`.
2. If `grad_norm_snr` falls below a target band, rerun the current point at a higher budget before giving up on the step.
3. If the gradient remains noisy, promote the fit ladder first toward larger lattices and then promote `n_traj` if needed.
4. Near the floor, require the local gradient to be resolved at the promoted budget before accepting another move.
5. Stop only after the largest allowed budget still cannot recover a robust descent signal.

There are two natural versions of this:

- trajectory-first escalation: keep the same sizes, but increase `n_traj`
- ladder-first escalation: drop the noisiest small sizes or add emphasis to the largest sizes before increasing trajectories

My current expectation is that a mixed schedule will work best:

- use fewer / cheaper sizes far away
- broaden and stabilize the ladder near the minimum
- only pay for high `n_traj` once the walker is already in a narrow basin

This should preserve the basin-local efficiency of the gradient method while avoiding the current failure mode where the walker stops simply because the cheap-budget local gradient becomes noise-dominated.

## 2026-06-08: Local Reweight Gradient Stepping Looks Plausible As A Basin-Local Optimizer, But Not As Uncontrolled Fixed-Step Descent

### Question

Evaluate whether the geometry-fit workflow could shift from a grid search to a loop of the form:

- choose one coupling point `(r1, r2)`
- run direct Monte Carlo there
- use reweighting from that same run to estimate the local score gradient
- take a downhill step
- run Monte Carlo again at the new point
- repeat

The comparison point is the existing local-grid workflow, where many points are run directly and then ranked.

### Data Used

I used the new direct local-reweight gradient outputs for the iso111 control:

- gradient PNG:
	- `responsible_method_tests/reweighting/iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606/gradient_fields_target_L64_qw_norm_local_reweight/target_L64_qw_norm_local_reweight_aggregate_gradient_fields.png`
- gradient TSV:
	- `responsible_method_tests/reweighting/iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606/gradient_fields_target_L64_qw_norm_local_reweight/target_L64_qw_norm_local_reweight_aggregate_gradient_fields.tsv`
- refined `2x` score landscape:
	- `responsible_method_tests/reweighting/iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606/heatmaps_target_L64_qw_norm_refined2x_neighbor_interp/target_L64_qw_norm_refined2x_neighbor_interp_landscape.tsv`

Important distinction: this gradient field is **not** the old `np.gradient` of an interpolated surface. Each arrow now comes from local reweighting around the current base donor point at

$$
(r_1, r_2) \to (r_1 \pm \delta, r_2), (r_1, r_2 \pm \delta)
$$

with `delta = 0.0125`, using the exact per-trajectory energy decomposition already saved in the selected-observable bundles.

So this test is really about a **local reweight-derived optimizer step**, not about trusting a globally interpolated score surface.

### Visual Read Of The Field

The total-score field is not random-looking. Over most of the `5 x 5` base grid, the arrows point coherently toward the lower-score basin around

- `r1 ~ 0.95 - 1.00`
- `r2 ~ 0.95 - 1.00`

The large-score corners `(1.10, 0.90)` and `(0.90, 1.10)` both point back toward that basin. That is already qualitatively different from a useless noisy vector field.

### Quantitative One-Step Check

To turn that into something harder than eyeballing, I took every base-grid point, moved one normalized steepest-descent step of length `0.025`, snapped that proposal to the nearest refined `2x` grid node, and compared the refined total score before and after.

Result over all `25` base points:

- improved: `23 / 25`
- unchanged: `1 / 25`
- worsened: `1 / 25`

The unchanged case was the lower-left boundary point `(0.90, 0.90)`, where the downhill direction points outside the box and the snapped proposal stays on the boundary.

Representative cases:

| start | start score | downhill proposal | proposal score | change |
| --- | ---: | --- | ---: | ---: |
| `(1.00, 1.00)` | `1.74612e-05` | `(0.975, 0.975)` | `5.2265e-06` | `-1.22347e-05` |
| `(1.10, 0.90)` | `9.36150e-05` | `(1.075, 0.925)` | `4.46532e-05` | `-4.89618e-05` |
| `(0.90, 1.10)` | `6.72256e-05` | `(0.925, 1.075)` | `3.48513e-05` | `-3.23743e-05` |
| `(0.95, 0.95)` | `9.12240e-06` | `(0.950, 0.975)` | `7.1705e-06` | `-1.95190e-06` |
| `(0.95, 1.00)` | `3.80040e-06` | `(0.975, 1.000)` | `5.9690e-06` | `+2.16860e-06` |

This is the key evidence.

Far from the basin, the local reweight gradient gives a very usable downhill direction. Near the apparent minimum, a fixed step can already overshoot and get worse.

### Interpretation

This makes the gradient-step idea look **materially better** than the earlier nearest-donor surface interpolation idea.

The earlier failure mode was:

- use a donor grid to infer missing points elsewhere
- trust that interpolated/reweighted surrogate enough to rank the whole surface

That failed badly even when histogram overlap stayed excellent.

What is different here is that the reweighting is used only **locally around the current Monte Carlo point**. That is a much narrower and more defensible use case.

The present evidence supports the following statement:

- local reweighting is probably good enough to tell you a downhill direction inside a basin

but it does **not** yet support the stronger statement:

- a naive fixed-step gradient descent can safely replace grid search everywhere

There are at least four reasons for that caution.

First, the near-minimum counterexample is real. At the current base-grid best point `(0.95, 1.00)`, the local gradient is not small enough to stop a fixed `0.025` step from moving uphill on the refined landscape.

Second, the objective is a compromise of multiple sectors, and those sectors do not have the same preferred location. On the refined landscape,

- `midpoint_score` best is at `(0.95, 1.00)`
- `quarter_score` best is at `(1.075, 1.050)`
- `total_score` best is again at `(0.95, 1.00)`

So the total score is not a simple isotropic bowl; it is a negotiated multi-observable landscape. That makes overconfident raw gradient descent risky.

Third, the current gradient field is sampled only on the coarse `5 x 5` base nodes. The proposed workflow itself avoids that limitation, because after each accepted move you would run a new direct MC and recompute the gradient there. But it does mean that the present test validates **one-step local usefulness**, not full multi-step convergence.

Fourth, the gradient magnitude changes a lot across the box. For the total score it runs from about

- `5.31e-05` near `(0.95, 1.00)`
- up to `8.94e-04` at `(1.10, 0.90)`

So a single fixed spatial step is not likely to be optimal everywhere.

### What This Says About Workflow Design

The practical conclusion is:

- do **not** replace the grid search with uncontrolled fixed-step gradient descent
- do consider replacing it with a **local optimizer loop** driven by direct MC plus local reweight probes

In other words, the promising object here is not “gradient descent” in the abstract. It is a **trust-region / line-search workflow** whose local model comes from exact reweighting around the current direct run.

The sensible version would look like:

1. run direct MC at the current point `x_k`
2. use reweighting from that same payload to estimate the gradient and evaluate a few candidate points along `-grad`, for example `alpha in {delta/2, delta, 2 delta}`
3. accept a step only if the predicted score decrease is clear and the local overlap / effective sample size remain healthy
4. rerun direct MC at the accepted point `x_{k+1}`
5. repeat until no tested step gives a robust decrease

That is a much more defensible workflow than either of these extremes:

- full brute-force dense grid search everywhere
- blind raw gradient descent with one hard-coded step size

### Recommendation

My current recommendation is a **hybrid**:

- use a coarse grid search only to locate the basin and guard against missing disconnected valleys
- once inside a basin, switch to direct-MC plus local-reweight stepping
- use line search or a trust radius, not a fixed step
- near the floor, validate each accepted move with another direct run before trusting the direction further

So the answer is:

- **yes**, this looks promising as a replacement for the expensive inner part of a grid search
- **no**, it is not yet convincing as a full drop-in replacement for global search without safeguards

The present iso111 field is strong evidence that local reweight-derived gradients contain real optimization signal. The same field also shows why the method should be wrapped in a conservative step-acceptance rule rather than used as pure unconstrained gradient descent.

## 2026-06-08: Implemented A Direct-MC Plus Local-Reweight Gradient-Flow Driver

I added

- `responsible_method_tests/scripts/run_reweight_gradient_flow.py`

to turn the previous gradient-field reasoning into an executable workflow.

### What It Does

The script runs a discrete optimizer loop:

1. start from a user-chosen `(r1, r2)`
2. run direct MC for all requested fit sizes at that point
3. estimate the local score gradient by reweighting that same direct payload to
	- `(r1 - delta, r2)`
	- `(r1 + delta, r2)`
	- `(r1, r2 - delta)`
	- `(r1, r2 + delta)`
4. propose a downhill step
5. run direct MC at the proposal point
6. either accept the step or stop, depending on the acceptance rule

So this is the actual “direct MC + local reweight gradient + direct MC again” workflow, not just a post-hoc plotter.

### User Controls

The main user-facing controls are now exposed directly on the CLI:

- `--start-r1`, `--start-r2`
- `--step-size`
- `--num-steps`

plus supporting knobs:

- `--step-mode {normalized,raw}`
- `--gradient-step`
- `--accept-rule {always,direct_decrease}`
- `--metric-key`
- target size / target point / fit sizes / MC statistics / bounds / output root

Default objective is the total score for the chosen balance mode.

Default acceptance is conservative:

- `direct_decrease`

meaning the script only accepts a step if the **next direct MC score** is no worse than the current direct score.

### Outputs

For each run it writes:

- `trajectory.tsv`
- `trajectory.png`
- `summary.json`

under the chosen output root, while reusing cached direct raw bundles if they already exist.

### Smoke Validation

A small low-stat smoke run completed successfully with

- start `(1.0, 1.0)`
- `step_size = 0.02`
- `num_steps = 2`
- fit sizes `{12, 8, 4}`
- target size `16`
- `n_traj = 200`, `n_therm = 50`, `n_skip = 10`

Output root:

- `responsible_method_tests/results/_smoke_reweight_gradient_flow/`

In that smoke run the first step was accepted, the second was rejected by the conservative rule, and the script stopped with

- `final_reason = direct_score_increase`

which is exactly the intended safeguard behavior near an overshooting proposal.

## 2026-06-04: Nearest-Donor Histogram Reweighting Does Not Yet Preserve The Iso111 Geometry-Match Grid

### Goal

Check whether the nearby-coupling histogram-reweighting machinery is already good enough
to interpolate a small geometry-matching search grid, so that a coarse set of directly-run
donor points could stand in for a denser direct scan.

The concrete question here was narrower than the full production problem:

- target geometry fixed at iso111 `(r1, r2) = (1.0, 1.0)`
- compare a direct `5 x 5` local score surface against a `3 x 3` donor grid
- use nearest-donor reweighting to predict the missing fine-grid points
- ask whether the final geometry-match proxy score is preserved well enough to keep the same valley/minimum structure

### Test Script

I added

- `responsible_method_tests/scripts/test_geometry_match_grid_interpolation.py`

which runs the selected-observable bundles on a local coupling grid, builds both

- the **direct** finite-size family for each fine-grid point, and
- the **reweighted** family obtained from the nearest donor point,

then compares the two after the same blind finite-size extrapolation to a `32 x 32`
target holdout.

### Setup

- fine grid: `r1, r2 in {0.99, 0.995, 1.0, 1.005, 1.01}`
- donor grid: `r1, r2 in {0.99, 1.0, 1.01}`
- fit sizes: `L in {8, 12, 16, 24}`
- target holdout size: `L = 32`
- statistics at every run point: `n_traj = 20000`, `n_therm = 2000`, `n_skip = 10`
- observables stored per trajectory:
	- midpoint `v, u, w`
	- quarter-point `v, u, w`
	- anchor `(0, -1)`
- score panels built from those observables:
	- midpoint-to-anchor: `v/a`, `u/a`, `w/a`
	- midpoint ratios: `v/u`, `w/u`
	- quarter-to-anchor: `v/a`, `u/a`, `w/a`
	- quarter ratios: `v/u`, `w/u`
- total proxy score: `z2_sum = sum(panel_z^2)` against the `L=32`, `(1,1)` target

Output tree:

- `responsible_method_tests/results/geometry_match_reweight_interp_iso111_grid5x5_donor3x3_sizes8_12_16_24_target32_20260604/`

Key files:

- `summary.json`
- `grid_scores.tsv`
- `geometry_match_grid_interpolation.png`

### Main Result

The raw histogram overlap stayed excellent everywhere, but the **interpolated score surface did not**.

Across all reweighted panel fits, the minimum reported effective sample-size fraction stayed above

- `N_eff / N > 0.99985`

so this is **not** a simple overlap-collapse failure.

At donor points themselves, direct and reweighted scores agree essentially exactly, as they should.

But away from donor points, the final geometry-match proxy can shift a lot after the ratio-building
and finite-size fit stages:

- RMS difference over the `5 x 5` grid: `delta_z2_sum_rmse = 4.116`
- worst absolute discrepancy: `max |delta_z2_sum| = 13.680`

The direct and reweighted best points on the same `5 x 5` surface were:

| surface | best point `(r1, r2)` | score |
| --- | --- | ---: |
| direct | `(0.99, 1.00)` | `0.5691` |
| reweighted | `(0.995, 1.00)` | `0.5355` |

So even in this tiny local control window, the minimum moved by one half-step in `r1`.

More importantly, some individual cells changed by an amount large enough to completely alter the local ranking. The clearest failure in this pilot was

- point `(0.99, 1.005)` using donor `(0.99, 1.00)`
- direct `z2_sum = 14.2847`
- reweighted `z2_sum = 0.6050`
- `delta_z2_sum = -13.6797`

despite the same point still having

- `N_eff / N = 0.999861...`

on every reweighted panel.

### Interpretation

This sharply limits what can be concluded from the earlier nearby-correlator control.

That earlier test showed that a **single selected correlator** can be reweighted reliably between
very nearby coupling points. This new test shows that this is **not yet enough** to guarantee a stable
interpolation of the downstream geometry-match objective.

The failure is happening after several nonlinear steps are composed:

- connected-correlator subtraction
- ratio formation
- multi-panel comparison to the target
- blind finite-size extrapolation on each panel
- final aggregation into `z2_sum`

So excellent histogram overlap at the observable level does not automatically imply that the
derived score landscape is preserved.

### Takeaway

Nearest-donor histogram reweighting is **not yet safe as a drop-in interpolant for the geometry-matching grid search**, even in this small iso111 control box.

The current result supports a narrower statement instead:

- nearby observable reweighting works as a local measurement tool
- nearest-donor patch interpolation does not yet preserve the fine-grid score ranking

If this idea is pursued further, the next thing to test should be something stricter than nearest-donor patching, for example a multidonor blend or a more local direct-vs-reweighted validation rule before trusting any interpolated minimum.

## 2026-06-04: Three-Point Nearby-Coupling Correlator Reweighting Control Works At 32x32

### Goal

Test whether a single selected correlator can be reweighted between nearby points in
`(r1, r2)` space, rather than only in `beta`, and check that the reweighted predictions
agree with direct Monte Carlo runs.

### What Had To Change

The existing Ferrenberg-Swendsen trace dump in

- `K_from_continuum/src/ising_tri_twisted_parallelogram.cc`

already wrote per-measurement

- `E_per_site`
- `abs_m`
- `m^2`

which is enough for susceptibility reweighting in `beta`, but not enough for
reweighting a correlator across nearby coupling points.

For this test, the existing single-displacement sample dump was extended so that when

- `--single_disp_m`
- `--single_disp_n`
- `--single_disp_samples_name`

are used, the per-measurement sample file now also records the uncoupled directional
bond energies per site:

- `e1_per_site`
- `e2_per_site`
- `e2me1_per_site`

with

$$
E_{\mathrm{site}} = K_1 e_1 + K_2 e_2 + K_3 e_{2-e1}.
$$

This is enough to reweight a selected correlator between nearby `(r1, r2)` points at
their exact triangular critical `beta` values.

### Test Script

I added

- `responsible_method_tests/scripts/test_correlator_reweight_nearby.py`

which runs three nearby points directly, loads the selected-displacement trace file,
and compares direct connected-correlator estimates against all pairwise reweighted
predictions.

### Setup

- lattice: `(Lx, Ly, Tx, Ty) = (32, 32, 0, 0)`
- selected displacement: `(m, n) = (8, 0)`
- statistics: `n_traj = 100000`, `n_therm = 10000`, `n_skip = 10`
- three nearby coupling points:
	- `center = (1.000, 1.000)`
	- `r1_plus_0p01 = (1.010, 1.000)`
	- `r2_plus_0p01 = (1.000, 1.010)`
- exact critical `beta` values from the triangular sinh-rule:
	- center: `0.274653072167`
	- both anisotropic neighbors: `0.273741433831`

Output tree:

- `responsible_method_tests/results/correlator_reweight_nearby_iso111_L32_hi100k_20260604/`

Key files:

- `summary.json`
- `reweight_vs_direct.tsv`

### Direct Results

The direct connected correlator at `(8,0)` was:

| point | connected | sigma | wall time |
| --- | ---: | ---: | ---: |
| center | `0.4521014` | `0.0017147` | `36.34 s` |
| `r1 + 0.01` | `0.4508243` | `0.0017368` | `33.76 s` |
| `r2 + 0.01` | `0.4500686` | `0.0017325` | `35.93 s` |

So the nearby anisotropic points really do move the correlator by a measurable but small
amount, of order `1e-3` to `2e-3`.

### Reweighting Result

All off-diagonal reweight predictions agreed with the direct target run within `1 sigma`.

Off-diagonal `z = (prediction - direct) / sigma_combined` values were:

| source -> target | `z` | `N_eff / N` |
| --- | ---: | ---: |
| center -> `r1 + 0.01` | `+0.584` | `0.9990` |
| center -> `r2 + 0.01` | `+0.833` | `0.9990` |
| `r1 + 0.01` -> center | `-0.480` | `0.9990` |
| `r1 + 0.01` -> `r2 + 0.01` | `+0.183` | `0.9970` |
| `r2 + 0.01` -> center | `-0.938` | `0.9990` |
| `r2 + 0.01` -> `r1 + 0.01` | `-0.258` | `0.9970` |

The worst discrepancy was still only about `0.94 sigma`, and all effective sample-size
fractions remained above `0.9969`.

### Interpretation

This is a clean positive control for **nearby-coupling correlator reweighting**.

At least for this `32x32`, `100k` iso111 control test and this selected displacement,
the direct runs and the reweighted predictions are statistically consistent. So the
basic parameter-space reweighting idea is not blocked at the first sanity-check level.

The test does **not** yet establish that this will remain valid for:

- larger coupling moves
- sharper anisotropies
- larger lattices
- multi-observable score panels
- full all-to-all or boundary-quarter selected sets

but it does show that the minimal nearby-point version is working.

### Takeaway

The next sensible step is not a full dense production grid yet. It is a slightly larger
pilot using the same reweighting idea on a small selected observable set, for example the
quarter-boundary `v`, `u`, `w`, and anchor observables, to see how rapidly the reweight
agreement degrades as the move in `(r1, r2)` grows.

## 2026-05-29: Iso111 Modular-Analogue Ratio Trough Directions Split Away From Quarter Point

### Goal

Check the precise meaning of the apparent ratio-sector trough in the iso111 modular-analogue control scans.

The narrow question was not whether the `v/u` and `w/u` panels both have low-score regions, but whether those low-score regions point in the **same direction in `(r1, r2)` space**.

### Setup

I used the exact modular anchored-ratio r-grid script

- `modular analogue/plot_modular_anchored_ratio_rgrid_heatmaps.py`

for the iso111 control case

- target geometry `(12, 12, 0, 0)`
- truth `(r1, r2) = (1.0, 1.0)`
- scan box `r1, r2 in [0.6, 1.4]`

at four sample fractions along the correlator direction:

- `1/8`
- `1/4`
- `3/8`
- `1/2`

To define the trough direction, I took the low-score region of each ratio surface and fit its principal axis in `(r1, r2)` space. Using the lowest `2%`, `5%`, and `10%` of each score surface gave essentially the same answer, so the direction estimate is stable.

### Main Result

The `v/u` trough direction is stable and diagonal, but the `w/u` trough direction is **not**.

At the `5%` low-score level, the fitted trough slopes are:

| fraction | `v/u` trough slope | `w/u` trough slope |
| --- | ---: | ---: |
| `1/8` | `1.000` | `1.539` |
| `1/4` | `1.000` | `1.005` |
| `3/8` | `1.000` | `0.288` |
| `1/2` | `1.000` | `0.000` |

So the correct statement is:

- `v/u` always prefers the near-diagonal family `r2 ~ r1`
- `w/u` aligns with that diagonal only near the quarter point
- away from quarter point, `w/u` rotates strongly

By `1/2`, the `w/u` trough is essentially horizontal, i.e. much closer to `r2 ~ const` than to `r2 ~ r1`.

### Interpretation

This sharpens the earlier aggregate-trough observation.

The combined ratio sector and total score can still show a common soft valley through `(1,1)`, but that does **not** mean the two constituent ratio observables are probing the same geometric direction.

What is really happening is:

- `v/u` supplies a robust diagonal trough across all four fractions
- `w/u` agrees with it at `1/4`
- `w/u` becomes a different directional probe at `3/8` and especially at `1/2`

So quarter point is special here: it is the point where the two ratio sectors are closest to sharing the same parameter-space trough direction.

### Takeaway

For the iso111 modular-analogue control, it is wrong to say that `v/u` and `w/u` generically have the same trough direction.

The more accurate statement is:

- they are directionally aligned at `1/4`
- they are already different at `1/8`
- they clearly split by `3/8`
- by `1/2`, `w/u` is nearly horizontal while `v/u` remains diagonal

## 2026-05-29: Sigma-Operator Normalization Control Via `G_conn / <sigma>^2`

### Idea

There is another clean way to control for the arbitrary lattice normalization of the spin operator.

If the lattice spin field is related to the continuum operator by

$$
\sigma_{\mathrm{lat}} = Z_\sigma \, \sigma_{\mathrm{cont}},
$$

then:

- the one-point function scales like

$$
\langle \sigma_{\mathrm{lat}} \rangle = Z_\sigma \, \langle \sigma_{\mathrm{cont}} \rangle,
$$

- while the connected two-point function scales like

$$
G_{\mathrm{conn,lat}}(x,y)
=
\langle \sigma_{\mathrm{lat}}(x) \sigma_{\mathrm{lat}}(y) \rangle_c
=
Z_\sigma^2 \, G_{\mathrm{conn,cont}}(x,y).
$$

Therefore the ratio

$$
\frac{G_{\mathrm{conn,lat}}(x,y)}{\langle \sigma_{\mathrm{lat}} \rangle^2}
$$

is independent of the unknown multiplicative normalization `Z_sigma`.

### Why This Is Useful

This gives a second normalization-control strategy, distinct from anchor-ratio normalization.

- anchor-ratio normalization divides one two-point observable by another two-point observable
- this alternative divides the connected two-point function by the square of a one-point function

So it removes any implicit overall scale in the lattice sigma operator without referring to a separate anchor displacement.

Operationally, this could provide a more amplitude-aware normalization control than the current anchored-ratio construction, because it ties the two-point amplitude directly to the one-point normalization of the same operator.

### Main Caveat

This only makes sense when the one-point function is actually nonzero and under control.

On a finite torus in a symmetry-preserving ensemble, the plain one-point function usually vanishes by symmetry, so this ratio would be undefined or dominated by noise.

So this method is only viable in a setting where at least one of the following is true:

- the ensemble has an explicitly nonzero one-point function
- the symmetry is broken or fixed by boundary conditions
- a consistent substitute for the one-point amplitude is defined and measured with good precision

If those conditions are met, then plotting

$$
G_{\mathrm{conn}} / \langle \sigma \rangle^2
$$

is a principled way to remove the normalization ambiguity of the lattice sigma operator.

### Status

This is now recorded as an available normalization-control idea for future score construction.

It is potentially useful exactly because it is not just another reweighting of the existing anchored-ratio panels, so it may offer a genuinely more independent comparison surface if a reliable nonzero one-point function is available.

## 2026-05-29: Why The Acute456 Anchor-Ratio Landscape Shows A Trough

### Goal

Understand why the acute456 anchor-ratio scores form a clean diagonal trough that passes near the story target

- `(r1, r2) = (4.702783, 7.353910)`

and decide whether a mixed correlator-plus-ratio score can provide a genuinely independent crossing surface.

### Main Observation

The acute456 anchor-ratio landscape is not just a noisy set of isolated low cells. It has a real soft direction.

From the full `5 x 5` orbit-reduced pointwise anchor-ratio table

- `responsible_method_tests/standard/acute456_pointwise_score_landscape_anchor_ratio.tsv`

the top six candidates are:

| rank | `(r1, r2)` | orbit RMS z |
| --- | --- | ---: |
| 1 | `(3.762226, 6.618519)` | `0.1603` |
| 2 | `(4.702783, 7.353910)` | `0.1918` |
| 3 | `(5.173061, 8.089301)` | `0.3374` |
| 4 | `(5.643339, 8.089301)` | `0.3376` |
| 5 | `(5.643339, 8.824692)` | `0.3548` |
| 6 | `(5.173061, 7.353910)` | `0.3608` |

This is already the signature of a trough rather than a single isolated minimum: the low-score points march roughly along the direction where both `r1` and `r2` increase together.

### Local Evidence Around The True Point

The target is a genuine local minimum on the coarse `5 x 5` grid, but it sits inside a soft diagonal valley rather than at the bottom of a sharply isolated bowl.

Target score:

- `(4.702783, 7.353910)` -> orbit RMS z = `0.1918`

Immediate neighbors:

- left `(4.232505, 7.353910)` -> `0.7024`
- right `(5.173061, 7.353910)` -> `0.3608`
- down `(4.702783, 6.618519)` -> `0.4348`
- up `(4.702783, 8.089301)` -> `1.3304`
- up-right `(5.173061, 8.089301)` -> `0.3374`
- down-left `(4.232505, 6.618519)` -> `0.7154`

So the truth is not being bypassed by the trough. The trough really does intersect the true point on this coarse grid. The soft direction is just visibly tilted: moving up-right or down-left hurts the score much less than moving purely left/right/up/down, especially much less than moving upward in `r2` alone.

### Quarter-Point Boundary Ratios Show The Same Direction

The same pattern shows up in the quarter-point small-target boundary ratio sector

- `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep_anchor_ratio/acute456_boundary_quarter_fss_anchor_ratio_small_target_zscores.tsv`

If the ratio-sector minimum in each `r1` row is listed, the best `r2` values march upward with `r1`:

| fixed `r1` | best `r2` in `ratio_chi2` | ratio chi^2 |
| --- | --- | ---: |
| `3.762226` | `6.618519` | `0.2924` |
| `4.232505` | `7.353910` | `0.2302` |
| `4.702783` | `8.089301` | `0.0584` |
| `5.173061` | `8.824692` | `0.0624` |
| `5.643339` | `8.824692` | `0.9195` |

Column minima show the same diagonal drift:

| fixed `r2` | best `r1` in `ratio_chi2` | ratio chi^2 |
| --- | --- | ---: |
| `6.618519` | `3.762226` | `0.2924` |
| `7.353910` | `4.232505` | `0.2302` |
| `8.089301` | `4.702783` | `0.0584` |
| `8.824692` | `5.173061` | `0.0624` |

So the trough direction is not an artifact of the full pointwise scorer. The boundary ratio sector by itself already prefers the same rising diagonal family.

### Interpretation

The trough makes physical sense once the anchor-ratio construction is viewed as a projection.

What anchor-ratio normalization does is:

- divide out a common amplitude factor by comparing each sampled point to the same anchor point
- suppress finite-target common-mode normalization drift
- emphasize relative directional shape rather than absolute correlator magnitude

At criticality, with `beta_c` re-tuned at every candidate, that means the score is mostly probing the shape of the emergent anisotropic metric seen by the sampled points, not the absolute overall scale of the correlator.

That naturally creates a soft direction. A whole one-parameter family of nearby couplings can preserve almost the same relative distances from the anchor to the sampled quarter/boundary points, even while the absolute correlator amplitudes move noticeably. In other words:

- the anchor-ratio observables strongly constrain one combination of anisotropy deformations
- they are softer along the combination that looks more like a common rescaling of the sampled geometry after the exact `beta_c(r1, r2)` retuning has absorbed part of the overall scale change

On the current grid, that soft combination is the approximate down-left / up-right diagonal.

This is also why varying only one coupling at fixed partner is expensive. Pure horizontal or vertical moves spoil the relative directional balance quickly, while coordinated moves in both couplings can keep the anchor-normalized pattern more nearly intact.

### Why The Trough Passes Near The Truth

The truth should lie on the trough if the trough is really encoding the target geometry, because the true couplings are one member of the family that reproduces the right relative shape.

The fact that the truth is a local minimum but not always the global minimum on the coarse grid is then naturally explained by the remaining nuisance effects:

- finite target size
- anchor choice `(m, n) = (0, -1)`
- limited observable set (all built from the same anchor-normalized family)
- current MC noise / coarse `5 x 5` resolution

So the right picture is not

- “the trough is wrong and should disappear”

but rather

- “the trough is the structural geometry signal, and the remaining task is to find another observable family that cuts across it.”

### Mixed-Score Follow-Up: Not Independent Enough

I tested five mixed correlator-plus-ratio candidates built from the same quarter-point small-target TSV and wrote the artifacts under

- `responsible_method_tests/standard/mixed_scores/`

with summary table

- `responsible_method_tests/standard/mixed_scores/acute456_boundary_quarter_fss_anchor_ratio_small_target_zscores_mixed_score_best_points.tsv`

All five mixed candidates chose the same minimum:

- `(r1, r2) = (5.173061, 8.824692)`

That is strong evidence that these mixtures are **not** an independent second metric. They only reweight the same underlying anchored-correlator and ratio residuals, so they inherit the same diagonal trough instead of producing a new crossing trough.

### What A Genuinely Independent Second Trough Should Use

If the goal is a second trough that intersects the anchor-ratio trough near the truth, the second score should not be another recombination of the same anchor-normalized panels.

It should retain information that the anchor-ratio construction intentionally throws away, for example:

- an absolute-amplitude channel tied to the anchor itself
- a different spatial fraction or a different point set, not the same quarter-point family
- an observable family with a different continuum projection, such as a fit-based amplitude, local slope, or another nonredundant shape descriptor

The present mixed scores are useful as diagnostics, but not as independent geometry measurements.

### Current Takeaway

- the acute456 anchor-ratio trough is structural, not a plotting accident
- it reflects the fact that anchor-ratio scoring largely projects onto relative-shape information and quotients out common amplitude information
- the trough intersects the true point, which is exactly what we would want from a real geometry signal
- the remaining miss of the global minimum is likely a finite-target / anchor / observable-set issue, not evidence that the trough itself is meaningless
- mixed scores built from the same anchored panels do not produce a second independent crossing surface

So the next real scoring improvement should aim for a deliberately different observable class, not just a different weighting of the current anchored-ratio ingredients.

## 2026-05-27: Acute456 Full 5x5 Anchor-Ratio Landscape

### Goal

Test the stronger version of the anchor-ratio hypothesis for the standard acute `4-5-6` campaign:

- run the full `5 x 5` acute456 candidate grid with `normalization_mode = anchor_ratio`
- use the existing orbit-reduced pointwise target score as the ranking metric
- check whether the target candidate `(r1, r2) = (4.702783, 7.353910)` is actually the global minimum of that landscape
- decide whether the origin-inclusive full-point aggregate should be part of the main score or only treated as a diagnostic

### Commands Used

Full acute456 anchor-ratio landscape:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py \
	--family acute456 \
	--normalization-mode anchor_ratio
```

Full-point anchor-ratio checks for the target and two challengers:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_all_base_points_individual_fits_anchor_ratio.tsv

python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--untwisted-dir responsible_method_tests/standard/data/acute456/untwisted/r1_5p643339__r2_8p089301 \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_r1_5p643339_r2_8p089301_all_base_points_anchor_ratio.tsv

python responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py \
	--untwisted-dir responsible_method_tests/standard/data/acute456/untwisted/r1_5p643339__r2_8p824692 \
	--normalization-mode anchor_ratio \
	--output responsible_method_tests/standard/acute456_r1_5p643339_r2_8p824692_all_base_points_anchor_ratio.tsv
```

### Artifacts Written

- landscape figure:
	- `responsible_method_tests/standard/standard_pointwise_score_landscapes_anchor_ratio.png`
- acute456 landscape table:
	- `responsible_method_tests/standard/acute456_pointwise_score_landscape_anchor_ratio.tsv`
- target full-point anchor-ratio table:
	- `responsible_method_tests/standard/acute456_all_base_points_individual_fits_anchor_ratio.tsv`
- challenger full-point anchor-ratio tables:
	- `responsible_method_tests/standard/acute456_r1_5p643339_r2_8p089301_all_base_points_anchor_ratio.tsv`
	- `responsible_method_tests/standard/acute456_r1_5p643339_r2_8p824692_all_base_points_anchor_ratio.tsv`

### Main Quantitative Result

The full acute456 `5 x 5` anchor-ratio landscape does **not** place the target at the global minimum.

Top five orbit-reduced anchor-ratio candidates from `acute456_pointwise_score_landscape_anchor_ratio.tsv`:

| rank | `(r1, r2)` | orbit RMS z | jackknife sigma |
| --- | --- | ---: | ---: |
| 1 | `(3.762226, 6.618519)` | `0.1603` | `0.0232` |
| 2 | `(4.702783, 7.353910)` | `0.1918` | `0.0253` |
| 3 | `(5.173061, 8.089301)` | `0.3374` | `0.0437` |
| 4 | `(5.643339, 8.089301)` | `0.3376` | `0.0309` |
| 5 | `(5.643339, 8.824692)` | `0.3548` | `0.0554` |

So the target is **rank 2 / 25**, not rank 1.

The result is therefore mixed but informative:

- anchor-ratio scoring does beat the two explicit challengers checked earlier
- but it still prefers a lower-coupling candidate `(3.762226, 6.618519)` over the target

### Full-Point Versus Orbit-Reduced Behavior

For the target candidate `(4.702783, 7.353910)`:

- full-point aggregate RMS z = `1.3129`
- orbit-reduced RMS z = `0.1918`
- worst point is the origin with `z = -10.3191`

For the two explicit off-target challengers that were checked with the full-point anchor-ratio fit script:

- `(5.643339, 8.089301)` gave full-point aggregate RMS z = `1.3536` and orbit-reduced RMS z = `0.3376`
- `(5.643339, 8.824692)` gave full-point aggregate RMS z = `1.4094` and orbit-reduced RMS z = `0.3548`

So the same origin-dominance pattern persists away from the target: once the origin is included, the full-point aggregate is inflated far above the orbit-reduced score and becomes a poor primary ranking statistic.

The operational lesson is the same as in the earlier pointwise tables:

- the origin carries a qualitatively different normalization signal
- once it is mixed into the full-point aggregate, that aggregate no longer cleanly tracks the geometry-realization question
- the orbit-reduced score is the more stable quantity for ranking candidates

### Score-Definition Decision

For the standard anchor-ratio workflow, keep the **orbit-reduced score as the primary geometry score** and treat the full-point aggregate as a **diagnostic only**.

Concretely:

- do **not** use the origin-inclusive full-point aggregate RMS z as the main ranking metric
- do use the orbit-reduced RMS z table/landscape as the primary selection surface
- keep the full-point aggregate, worst-point list, and origin residuals as debugging outputs that explain why a candidate looks bad

I did **not** change the landscape code for this, because `plot_standard_pointwise_score_landscapes.py` already ranks candidates by the orbit-reduced score. The new acute456 run supports keeping that design rather than switching the primary score to the origin-inclusive aggregate.

### Interpretation

- The anchor-ratio construction is useful: it separates the target from obvious off-target challengers much better than the raw all-point aggregate does.
- But it is not yet a complete realization score: even on the full acute456 `5 x 5` grid, the target is only the **second-best** candidate.
- The current evidence therefore supports the weaker statement:
	- anchor-ratio FSS contains real geometry information
	- orbit-reduced anchor-ratio is a reasonable primary score
	- but anchor-ratio alone does not yet recover the exact target as the unique minimum of the acute456 grid

### Current Takeaway

For standard acute456, the best-supported workflow is:

- use the orbit-reduced anchor-ratio landscape to rank candidates
- treat the full-point aggregate as a diagnostic for origin pathologies and other local failures
- do not claim yet that anchor-ratio FSS alone identifies the target geometry uniquely

### Follow-Up: L8-Ratio Cross-Check

I followed up with the second normalization choice `normalization_mode = l8_ratio`.

That cross-check shows that the anchor-ratio winner is **not** stable under a normalization change:

- `l8_ratio` landscape winner is `(5.643339, 8.089301)` with orbit RMS z = `0.2813 +/- 0.0275`
- target `(4.702783, 7.353910)` is rank `7 / 25` with orbit RMS z = `0.4382 +/- 0.0197`
- former anchor-ratio winner `(3.762226, 6.618519)` falls to rank `10 / 25` with orbit RMS z = `0.5430 +/- 0.0376`

So the specific anchor-ratio minimum at `(3.762226, 6.618519)` is normalization-dependent rather than a robust optimum.

### Follow-Up: Origin-Dropped Full-Point Summary

To make the full-point reporting less ambiguous, `fit_standard_acute456_all_base_points.py` now prints both:

- `aggregate_score` (origin included if present in the summary rows)
- `origin_dropped_aggregate_score` (same point set with the origin removed)

Under `l8_ratio`, dropping the origin barely changes the full-point target score:

- target `(4.702783, 7.353910)`: `0.4285 -> 0.4318`, with orbit-reduced RMS z = `0.4382`
- former anchor-ratio winner `(3.762226, 6.618519)`: `0.5374 -> 0.5417`, with orbit-reduced RMS z = `0.5430`

So for `l8_ratio`, the remaining mismatch is **not** mainly an origin pathology. The stronger interpretation is:

- the anchor-ratio off-target minimum is not normalization-stable
- once `l8_ratio` is used, the target miss is a broader shape/scoring bias, not just origin contamination

### Follow-Up: Direct L8-Ratio Winner Comparison

I then put the actual `l8_ratio` winner `(5.643339, 8.089301)` on the same direct full-point footing as the target.

For the `l8_ratio` winner `(5.643339, 8.089301)`:

- full-point aggregate RMS z = `0.2607`
- origin-dropped full-point aggregate RMS z = `0.2609`
- orbit-reduced RMS z = `0.2813`
- worst point is `(a,b) = (0.125, 0.000)` and its orbit partner, with `z = +0.6391`

For the target `(4.702783, 7.353910)` under the same normalization:

- full-point aggregate RMS z = `0.4285`
- origin-dropped full-point aggregate RMS z = `0.4318`
- orbit-reduced RMS z = `0.4382`
- worst point is `(a,b) = (0.250, 0.250)` and its orbit partner, with `z = -0.6666`

So the `l8_ratio` winner beats the target on every direct aggregate that was checked, not just on the orbit-reduced landscape score.

### Follow-Up: Landscape TSV Diagnostics

`plot_standard_pointwise_score_landscapes.py` now writes the pointwise aggregate diagnostics explicitly into the landscape TSV:

- `point_origin_dropped_rms_z`
- `point_origin_dropped_rms_z_jackknife_sigma`
- `point_origin_inclusive_rms_z`
- `point_origin_inclusive_rms_z_jackknife_sigma`

The older `point_rms_z` column was already the origin-dropped aggregate because the landscape builder excluded the origin before computing that pointwise RMS. The new columns make that fact explicit and add the origin-inclusive counterpart, so future sweeps can show immediately whether the origin is responsible for a bad ranking.

In the regenerated acute456 `l8_ratio` table, the difference is tiny for the current winner:

- `(5.643339, 8.089301)`: `0.260860` dropped-origin versus `0.260663` origin-inclusive

and likewise small for the target:

- `(4.702783, 7.353910)`: `0.431841` dropped-origin versus `0.428519` origin-inclusive

So the `l8_ratio` mis-ranking is not being driven by the origin.

### Status Of The Score Hypothesis

The current status of the hypothesis

- "we can build a score that picks the untwisted lattice which correctly realizes the target geometry from the twisted reference"

is: **not yet established**.

What the current evidence does support is the weaker statement:

- the score family contains real geometry information
- normalization matters a lot
- orbit-reduced anchor-ratio is currently the most stable primary acute456 metric among the tested choices

But the stronger target-selection claim is not yet supported, because:

- for acute456, `anchor_ratio` ranks the target only `2 / 25`
- for acute456, `l8_ratio` ranks the target only `7 / 25`
- for iso111, anchor-ratio fixes the raw `(0.8, 0.8)` inversion but still minimizes near `(1.1, 1.1)`, not exactly at `(1.0, 1.0)`
- the `l8_ratio` follow-up shows that the remaining acute456 failure is not mainly an origin artifact, so the scoring bias is deeper than one pathological point

So the honest summary is:

- we do have scores that are informative and nontrivial
- we do not yet have a normalization-stable score that reliably selects the correct untwisted target as the unique optimum
- the next scoring improvement has to address a broader shape-selection bias, not just the origin treatment

## 2026-05-21: Twisted-Reference 4x Trial

### Goal

Test the user hypothesis that the remaining twisted-versus-untwisted mismatch is mainly caused by the two continuum manifolds being sampled on different lattice point sets. The rerun changes the comparison workflow to mimic the intended coupling-to-geometry pipeline:

- build a denser twisted reference family
- fit the twisted continuum manifold first
- interpolate that twisted continuum manifold onto the untwisted continuum points
- use only the twisted-on-untwisted residuals to decide distinguishability

### Configuration

- Run config: `configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json`
- Comparison mode: `python scripts/align_manifolds.py --comparison-mode twisted_reference`
- Twisted family rule: base cell chosen to be the nearest integer multiple of the original pickGeo cell whose area is about `4 x` the untwisted `12x12` base, then scaled by `{1,2,3,4,5}`
- Untwisted family rule: keep the original size ladder `(12,24,36,48,60)`

### Run Status

- Completed benchmark: `geometry_acute_111`
- Completed comparison artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- Completed summary artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/SMOKE_SUMMARY.md`
- Interrupted benchmark: `geometry_acute_122`
- Interruption point: twisted scale 5, lattice `(160,100,40,20)`
- Interruption mode: manual `KeyboardInterrupt`, so no trustworthy campaign-wide conclusion exists yet for the 10-geometry set

### Measured Result For `geometry_acute_111`

This is the cleanest control point because both branches are isotropic and target the same equilateral torus.

| quantity | previous symmetric FSS run | 4x twisted-reference run |
| --- | ---: | ---: |
| run tag | `raw_manifold_fss_pickgeo10_fss_20260521` | `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` |
| comparison mode | `symmetric` | `twisted_reference` |
| twisted comparable points | `143` | `575` |
| untwisted comparable points | `143` | `143` |
| twisted alpha | `0.3047` | `0.4128` |
| twisted chi2/dof | `0.1843` | `0.2111` |
| untwisted alpha | `0.2863` | `0.4573` |
| untwisted chi2/dof | `0.0505` | `0.2859` |
| twisted-on-untwisted RMS z | `2.2363` | `3.2160` |
| untwisted-on-twisted RMS z | `2.1648` | `3.9717` |
| common-grid relative RMS | `0.0696` | `0.0907` |
| verdict | `distinguishable` | `distinguishable` |

### Interpretation

- The individual twisted and untwisted continuum fits still look healthy. Both chi2/dof values remain comfortably below `1` in the 4x trial.
- The asymmetric reference-style comparison does **not** reduce the disagreement for the equilateral control. The primary metric, twisted-on-untwisted RMS z, rises from `2.2363` to `3.2160`.
- The denser twisted family increases the twisted continuum point count from `143` to `575`, so this is not a case where the result is being driven by a lack of reference points.
- For `geometry_acute_111`, the hypothesis "the manifolds only look different because the point sets do not sit on top of each other" is therefore disfavored. Even when the twisted manifold is treated as the dense reference and sampled only at the untwisted continuum points, the mismatch stays well above the `Z = 2` marginal threshold.

### Operational Note

The 4x twisted-reference workflow is substantially more expensive than the previous symmetric run.

- `geometry_acute_111` twisted scale 5 already reached lattice `(120,120,0,0)` and took about `51 s` for one production job.
- `geometry_acute_122` was interrupted during twisted scale 5 at `(160,100,40,20)`.
- At `n_workers = 1`, the full 10-geometry campaign is practical but materially slower than the earlier FSS run.

### Current Scientific Takeaway

From the one completed control benchmark, denser twisted sampling by itself does not rescue indistinguishability. If the broader 10-geometry campaign is resumed, it should be treated as a stronger test of the continuum-equivalence claim, not as a formality expected to automatically collapse the manifolds.

### Next Useful Step

- Resume `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` from the current partial tree, ideally benchmark-by-benchmark, so runtime is easier to control.
- After at least a few additional geometries complete, regenerate the comparison plots with `--comparison-mode twisted_reference` and build a dedicated campaign summary for that run.

### Resume Later

Current saved state before stopping:

- run root: `responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/`
- completed benchmark manifest: `manifest_geometry_acute_111.json`
- completed comparison artifact: `geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- next blocked slice when interrupted: `geometry_acute_122`, twisted scale 5, lattice `(160,100,40,20)`
- config already has `"resume": true`, so rerunning the same generation command should reuse the completed payloads under the current run root

Resume command:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py \
	--config responsible_method_tests/configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json
```

After enough benchmarks finish, regenerate the asymmetric comparison artifacts with:

```bash
source /workspaces/newQFE/.venv/bin/activate
for manifest in responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_*.json; do
	python responsible_method_tests/scripts/align_manifolds.py \
		--benchmark-manifest "$manifest" \
		--comparison-mode twisted_reference
done
```

Then build the campaign summary with:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/summarize_benchmark.py \
	--comparison-mode twisted_reference \
	--output responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/PICKGEO10_TWISTED_REFERENCE4X_SUMMARY.md \
	--benchmark-manifest responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_111.json
```

and add the remaining `--benchmark-manifest` arguments once those benchmark manifests exist.

## 2026-05-21: Large-Twisted Fixed-Alpha Residual Collapse Test

### Goal

Test the simpler question behind the "large reference as continuum" idea without involving the untwisted branch: for one fixed twisted geometry, does the average per-site residual to a fixed-scale modular two-point function collapse as the twisted lattice is enlarged?

The residual was defined exactly as requested:

- choose one fixed geometry
- keep a single fixed modular normalization `alpha`
- compare the twisted MC connected correlator to `alpha * g_modular` at the corresponding lattice-point locations
- average the sitewise residuals over the matched non-origin lattice points

Because the CFT two-point function diverges near the origin while the MC correlator saturates at `1`, the comparison excluded the origin and used a single overall amplitude factor `alpha`.

### Benchmark Choice

I used `geometry_acute_567`, because it was already the hardest case for the earlier large-reference idea:

- continuum twisted vs continuum untwisted was **not distinguished**
- but the largest finite twisted lattice still sat visibly away from the untwisted continuum when sampled directly

That makes it a good stress test for the question "does a larger twisted lattice itself move toward the modular continuum?"

### Fixed-Alpha Definition

Two normalizations were checked.

1. The original twisted continuum-to-modular alignment from the corrected 5-scale run:
	 - `alpha = 0.21849098395920252`
2. The regenerated alignment after adding the larger twisted lattice:
	 - `alpha = 0.2189182337`

The qualitative conclusion was the same for both, so the result is not coming from retuning the amplitude.

### Practical Size Ceiling

The current all-to-all simulator cannot reach the user-suggested million-site twisted regime for this observable.

- `geometry_acute_567` scale `6` with lattice `(120,150,60,90)` runs successfully
- scale `7` already fails, even in an `n_traj=1` allocation probe
- the failure is therefore a simulator/storage ceiling, not a long-runtime issue

So the test could only be pushed to the **largest feasible twisted scale**, not to `O(10^6)` volume.

### Data Products

- Probe config: `configs/raw_manifold_large_ref_probe_567_20260521.json`
- High-stat twisted probe run root: `results/raw_manifold_large_ref_probe_567_20260521/`
- Fixed-alpha collapse script: `scripts/plot_fixed_alpha_modular_collapse.py`
- Final plot using regenerated run alpha:
	- `results/raw_manifold_large_ref_probe_567_20260521/geometry_acute_567/twisted/geometry_acute_567_twisted_fixed_alpha_modular_collapse.png`
- Strict fixed-alpha plot using the original pre-extension normalization:
	- `results/raw_manifold_large_ref_probe_567_20260521/geometry_acute_567/twisted/geometry_acute_567_twisted_fixed_alpha_modular_collapse_alpha021849.png`

### Main Result

Using the **strict original fixed normalization** `alpha = 0.21849098395920252`, the mean absolute residual per non-origin lattice site was:

| scale | family size | mean abs residual/site | RMS residual/site |
| --- | ---: | ---: | ---: |
| 1 | 1 | 0.19588758 | 0.19639174 |
| 2 | 2 | 0.11623277 | 0.11672905 |
| 3 | 3 | 0.08863659 | 0.08893889 |
| 4 | 4 | 0.06010132 | 0.06040861 |
| 5 | 5 | 0.05721346 | 0.05724360 |
| 6 | 6 | 0.04783519 | 0.04791435 |

This is a clear collapse with increasing twisted size. The scale-6 point lowers the mean absolute residual by about `16.4%` relative to scale 5:

- `0.05721346 -> 0.04783519`

Using the regenerated run alpha gives the same conclusion:

- scale 5 mean abs residual/site = `0.05668565`
- scale 6 mean abs residual/site = `0.04730739`

### Important Operational Detail

The first scale-6 attempt at only `n_traj=100` produced an apparent reversal (`~0.0858` at fixed old alpha), but a direct higher-stat check at `n_traj=500` showed that this was just a noisy outlier. After regenerating the scale-6 probe point at `n_traj=500`, the collapse resumed cleanly.

### Interpretation

- For this fixed geometry, the twisted lattice correlator does move toward the modular target as the twisted size is increased.
- The collapse is visible directly in the finite-lattice data using a single fixed amplitude normalization, without leaning on the untwisted branch.
- However, the current all-to-all implementation hits a hard ceiling already between scales `6` and `7` for this geometry, so this workflow cannot yet test the stronger hypothesis that an arbitrarily large twisted reference can simply replace the continuum limit.

### Current Takeaway

The finite-size trend supports the idea that a larger twisted reference is moving in the right direction, but with the present simulator the test stops well before the truly asymptotic regime. For non-integrable theories, this fixed-alpha twisted-vs-modular collapse is still a useful diagnostic, but it is not yet a substitute for a true continuum extrapolation.

## 2026-05-22: Matched-Volume Square-Untwisted Cross-Checks for 4-5-6 and 5-6-7

### Goal

Test the stricter matched-volume variant requested after the shared-coordinate 4-5-6 plots:

- keep the untwisted branch square, `L x L`
- choose twisted and untwisted base lattices with approximately matched areas above `100`
- handle physical side orientation explicitly through the correct untwisted `embedding_cycles`
- push both branches through the same `1..10` size ladder

This was also the first pass where I exposed the post-fit manifold difference directly instead of only as one panel inside the standard comparison plot.

### Configurations

Two acute scalene benchmarks were run.

1. `geometry_acute_456_matched_volume`
	- config: `responsible_method_tests/configs/raw_manifold_fss_acute_456_matched_volume_20260522.json`
	- twisted base: `[12,12,6,2]`, volume `156`
	- untwisted base: square `[12,12,0,0]`, volume `144`
	- untwisted `embedding_cycles=[0,1]`

2. `geometry_acute_567_matched_volume`
	- config: `responsible_method_tests/configs/raw_manifold_fss_acute_567_matched_volume_20260522.json`
	- twisted base: `[12,15,6,9]`, volume `234`
	- untwisted base: square `[15,15,0,0]`, volume `225`
	- untwisted `embedding_cycles=[2,1]`

Both runs used `n_traj=2000` and scales `1..10` on both branches.

### New Diagnostic

I added `responsible_method_tests/scripts/plot_continuum_difference.py` to render the post-fit difference explicitly from the same shared-cell comparison data already built by `align_manifolds.py`.

For each benchmark it writes three panels:

- common-grid heatmap of `twisted - untwisted`
- differences sampled on the untwisted fitted points
- differences sampled on the twisted fitted points

This makes it easier to see whether a result is just noisy pointwise scatter or a coherent smooth offset across the target cell.

### Data Products

For `geometry_acute_456_matched_volume`:

- run root: `responsible_method_tests/results/raw_manifold_fss_acute_456_matched_volume_20260522/`
- shared-coordinate FSS: `geometry_acute_456_matched_volume_shared_coordinate_twisted_untwisted_fss.{png,tsv}`
- shared-cell comparison: `geometry_acute_456_matched_volume_twisted_vs_untwisted_common_grid.{png,json}`
- dedicated post-fit difference: `geometry_acute_456_matched_volume_continuum_difference.png`

For `geometry_acute_567_matched_volume`:

- run root: `responsible_method_tests/results/raw_manifold_fss_acute_567_matched_volume_20260522/`
- shared-coordinate FSS: `geometry_acute_567_matched_volume_shared_coordinate_twisted_untwisted_fss.{png,tsv}`
- shared-cell comparison: `geometry_acute_567_matched_volume_twisted_vs_untwisted_common_grid.{png,json}`
- dedicated post-fit difference: `geometry_acute_567_matched_volume_continuum_difference.png`

### Measured Results

| quantity | 4-5-6 matched-volume | 5-6-7 matched-volume |
| --- | ---: | ---: |
| twisted chi2/dof | `0.0120` | `0.2011` |
| untwisted chi2/dof | `0.1547` | `0.1234` |
| common-grid RMS abs | `0.00754` | `0.01975` |
| common-grid relative RMS | `0.02877` | `0.07700` |
| twisted-on-untwisted RMS z | `0.7497` | `1.8487` |
| untwisted-on-twisted RMS z | `0.6789` | `1.9597` |
| verdict | `not distinguished within continuum point errors` | `marginally separated` |

### Shared-Coordinate FSS Cross-Check

The pointwise FSS plots tell the same story.

For `4-5-6`, the auto-selected points were `point_id = {1,2,4,6}`. The fitted continuum differences were small:

- `+0.00241` (`z=+0.31`)
- `-0.00437` (`z=-0.41`)
- `-0.00486` (`z=-0.35`)
- `-0.00558` (`z=-0.37`)

For `5-6-7`, the auto-selected points were `point_id = {1,3,5,7}`. All four fitted differences had the same sign and somewhat larger magnitude:

- `-0.00510` (`z=-0.75`)
- `-0.01505` (`z=-1.32`)
- `-0.01432` (`z=-1.04`)
- `-0.01523` (`z=-1.02`)

So the matched-volume `5-6-7` result is not just one bad point; it looks more like a coherent offset at the tested shared coordinates.

### Interpretation

- The earlier mild `4-5-6` separation was indeed a side-orientation artifact. Once the untwisted basis was corrected to the proper physical-side map, the benchmark moved to the clearly safe side of the continuum-error threshold.
- Repeating the same matched-volume protocol on `5-6-7` does **not** produce the same collapse. Both individual fits remain healthy, but the two continua retain a smooth, mostly positive `twisted - untwisted` offset across most of the shared target cell, with bidirectional RMS z values near `2`.
- The scientific lesson is therefore narrower than the optimistic `4-5-6` control alone suggested: correct orientation handling is necessary, but it is not by itself sufficient to make every matched-volume square-untwisted comparison collapse.

### Current Takeaway

At current statistics and with the present Taylor-2 continuum fits:

- `4-5-6` supports the claim "same continuum within current errors"
- `5-6-7` remains a borderline counterexample, landing in the `marginally separated` regime even after the orientation-aware matched-volume construction

The dedicated difference plots make that split especially clear: `4-5-6` shows a small residual field at the level of the fitted point errors, while `5-6-7` shows a larger smooth offset of order `2e-2` across much of the cell.

## 2026-05-22: Exact-Critical 1-1-1 Control Sweep

### Workflow Changes

Two continuum-workflow fixes landed before the control reruns:

- `K_from_continuum/workflow_common.py` now supports `beta_finder.mode = "exact_triangular_sinh_rule"` via `exact_triangular_ising_beta(...)`, so the triangular Ising critical point can be supplied exactly from
	`sinh(2 beta k1) sinh(2 beta k2) + sinh(2 beta k2) sinh(2 beta k3) + sinh(2 beta k3) sinh(2 beta k1) = 1`
- `K_from_continuum/02_generate_tests.py` stopped incorrectly forwarding `n_workers=` into `find_beta_for_lattice(...)`, which was breaking point jobs that used continuum-based beta strategies

The exact-beta path was validated immediately on the isotropic `1-1-1` case, giving `beta = ln(3)/4` as expected.

### Runs

I used that exact-critical path for three increasingly sharp `1-1-1` continuum-vs-reference controls:

- Pilot axial stencil: `K_from_continuum/results/iso111_distinguishabilityK_20260522/`
- Dense high-stat axial stencil: `K_from_continuum/results/iso111_distK_dense_hi5000_20260522/`
- Fresh narrow repeat control: `K_from_continuum/results/iso111_distK_repeat3_hi5000_20260522/`

All three used a twisted equilateral reference family and square untwisted test families at the exact critical couplings.

### Raw-Score Results

At the original raw continuum score, the equilateral control did **not** stabilize to the isotropic point as statistics increased.

- Pilot run: best point moved to `r1 = 1.01`, `r2 = 1.00`, with `z_rms = 0.6682`; isotropic `(1.00, 1.00)` still sat at `z_rms = 1.1588`
- Dense high-stat run: best point moved again to `r1 = 0.995`, `r2 = 1.00`, with `z_rms = 1.8069`; isotropic `(1.00, 1.00)` worsened to `z_rms = 3.2136`
- Fresh repeat run: the same off-isotropic preference repeated, with best point `r1 = 0.995`, `r2 = 1.00`, `z_rms = 2.2419`; isotropic `(1.00, 1.00)` still sat at `z_rms = 2.7362`

### Interpretation

- The failure mode sharpened with more statistics, which strongly argues against blaming the issue on the exact beta finder itself
- The repeat run showed that the bias was stable enough to reproduce, so the problem had to be in the reference/scoring construction rather than in one noisy outlier dataset
- That made the next question sharper: what part of the raw score is the workflow over-penalizing?

## 2026-05-23: Additive-Offset and Shape-Only Rescoring of the 1-1-1 Control

### Workflow Change

To answer that question, `K_from_continuum/03_score_continuum_vs_reference.py` now supports

- `analysis.score_mode = raw`
- `analysis.score_mode = additive_offset`
- `analysis.score_mode = shape_only`

and writes `score_diagnostics.dat` so the fitted nuisance offsets are visible instead of being implicit.

The idea is simple:

- `raw`: score the channels exactly as before
- `additive_offset`: subtract one best-fit common offset before scoring
- `shape_only`: subtract one best-fit offset per cycle before scoring

### Rescoring Results

Offline rescoring of the existing high-stat `1-1-1` artifacts changed the interpretation materially.

For the dense `hi5000` run:

- raw best: `(0.995, 1.000)` with `z_rms = 1.8069`
- additive-offset best: `(1.000, 1.000)` with `z_rms = 0.1999`
- shape-only best: `(1.000, 1.000)` with `z_rms = 0.0950`

For the repeat `hi5000` run:

- raw best: `(0.995, 1.000)` with `z_rms = 2.2419`
- additive-offset best: `(1.000, 1.000)` with `z_rms = 0.2839`
- shape-only best: `(0.995, 1.000)` with `z_rms = 0.1029`, but isotropic `(1.000, 1.000)` is still only `z_rms = 0.2145`

### Interpretation

- The dominant raw-score failure is a common-mode additive mismatch, not a large shape disagreement in the continuum channels
- Once that nuisance direction is removed, the isotropic point is either recovered outright or becomes essentially degenerate with the tiny `0.5%` offset candidate
- The `1-1-1` exact-critical problem is therefore much more about how the score is defined than about how `beta_c` is being chosen

## 2026-05-23: Responsible-Method Local 7x7 Iso111 Untwisted Scan

I also ran a narrow responsible-method local scan around the isotropic square target:

- run root: `responsible_method_tests/results/raw_manifold_fss_iso111_local_scan_grid7x7_hi5000_20260523/`
- target twisted manifest: `results/raw_manifold_fss_integer_multiples_20260520/geometry_111/twisted/manifest_geometry_111_twisted.json`

This did **not** recover isotropic couplings.

- best candidate: `(r1, r2) = (1.005, 1.020)`, score `1.3966`, interpretation `marginally separated`
- isotropic `(1.000, 1.000)`: score `3.2217`, interpretation `distinguishable at continuum-point level`

So on the responsible-method side, the near-isotropic local landscape was still favoring a noticeably shifted untwisted square candidate even after the continuum rescoring result above.

## 2026-05-24 to 2026-05-25: Wide Iso111 Responsible-Method Grid and Path-Length Fix

### Root Cause Fix

The wide untwisted square scan initially crashed at `(r1, r2) = (0.5, 1.0)`, but this turned out to be a simulator output-path bug rather than a Monte Carlo instability.

`K_from_continuum/src/ising_tri_twisted_parallelogram.cc` had been building runtime output subdirectories with fixed-size `char` buffers, and the generated scratch path hit that limit exactly for some of the widened-grid labels. That code now uses `std::string` for the run subdirectory and output paths.

### Production Run

After the path fix, the full widened responsible-method scan completed:

- twisted target run: `responsible_method_tests/results/raw_manifold_fss_iso111_twisted_target_hi10000_20260524/`
- untwisted square grid run: `responsible_method_tests/results/raw_manifold_fss_iso111_grid_r0p5to1p5_r0p5to1p5_step0p1_hi10000_20260524/`
- fit summary: `responsible_method_tests/results/raw_manifold_fss_iso111_grid_r0p5to1p5_r0p5to1p5_step0p1_hi10000_20260524/fit_vs_twisted_hi10000/geometry_111_coupling_fit_scan.json`

### Result

The best responsible-method square candidate on the full `11 x 11` grid was

- `iso111_scan_r11p200_r20p900`
- `(r1, r2) = (1.2, 0.9)`
- score `0.32237`
- twisted-on-untwisted RMS z `0.31967`
- untwisted-on-twisted RMS z `0.32237`
- common-grid relative RMS `0.01311`

The widened scan therefore still prefers a significantly anisotropic square candidate over isotropic `1-1-1`, even though the continuum rescoring analysis above showed that the raw continuum score was largely being driven by a nuisance offset.

## 2026-05-25: Fresh Responsible-Method Control Reruns and Larger-Lattice Follow-Up

### Fresh Small-Ladder Rerun

To check whether the wide-grid preference was just coming from a stale or unlucky reference, I rebuilt a fresh twisted `geometry_111` reference and compared it directly against two untwisted square candidates:

- run root: `responsible_method_tests/stupid_method/`
- config: `responsible_method_tests/configs/raw_manifold_fss_stupid_method_20260525.json`
- candidates: `(1.0, 1.0)` and `(1.2, 0.9)`

On the four selected shared-coordinate panels, the earlier dramatic isotropic self-mismatch did **not** reproduce.

- isotropic `(1.0, 1.0)` panel z values: `{-0.965, -0.832, -0.604, -0.441}`
- `(1.2, 0.9)` panel z values: `{-0.740, -0.809, -0.499, -0.199}`

So the fresh rerun puts both candidates inside `|z| < 1` on those selected shared points, with `(1.2, 0.9)` still somewhat closer.

### Larger-Lattice High-Stat Follow-Up

I then pushed the same control to larger ladders and higher statistics:

- run root: `responsible_method_tests/stupid_method/large_hi20000/`
- config: `responsible_method_tests/configs/raw_manifold_fss_stupid_method_large_hi20000_20260525.json`
- twisted scales: `1..7`
- untwisted sizes: `12..84`
- `n_traj = 20000`

On the same four selected shared-coordinate panels, the larger-lattice run still favors `(1.2, 0.9)` over isotropic, although isotropic improved materially.

- isotropic selected-panel RMS z: `0.3852`
- `(1.2, 0.9)` selected-panel RMS z: `0.1035`

### Small Workflow Extensions

Two small script changes supported these reruns:

- `responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py` now accepts explicit `entries=[{scale, family_size, lattice}, ...]`, so multiplicative scale and displayed family size can be decoupled when a family is not naturally expressed as one simple size ladder
- `responsible_method_tests/scripts/fit_target_geometry.py` now has `--annotate-points` to overlay candidate labels directly on the scan plot when that is useful

I also ran the size-`10..100` comparison bundle under `responsible_method_tests/stupid_method/iso111_and_far456_size10to100/`, which exported about `100` standalone shared-point panels per benchmark for slide triage.

## 2026-05-25: Far-Point Responsible-Method Control

To see whether the responsible-method ranking was just measuring distance from isotropic couplings, I compared the fresh twisted `1-1-1` reference against a square untwisted family at the **exact** `4-5-6` couplings:

- run root: `responsible_method_tests/stupid_method/far_456_vs_iso111/`
- config: `responsible_method_tests/configs/raw_manifold_fss_far_456_vs_iso111_20260525.json`
- square untwisted couplings: `(r1, r2) = (5.065231170850374, 7.742930488974055)`

The result was much less separated than naive “distance in coupling space” intuition would suggest.

- primary RMS z: `0.60139`
- twisted-on-untwisted RMS z: `0.59148`
- common-grid relative RMS: `0.02216`
- interpretation: `not distinguished within continuum point errors`

The selected shared-point z values were only `{1.64, 0.73, 0.50, 0.41}`.

So at current errors, the responsible-method comparison is clearly not just “how far are the couplings from `(1, 1)`?” A square untwisted family at the exact `4-5-6` point can still sit surprisingly close to the fresh twisted `1-1-1` reference.

## 2026-05-25: Acute 4-5-6 Literal-Size Blind Holdout

### Goal

Build a real blind holdout for the acute `4-5-6` shared-boundary midpoint using

- a literal-size untwisted square family
- a twisted family built from explicit multiples of the small `4-5-6` base `[6,6,3,1]`
- one held-out twisted lattice that is large enough to lie well outside the training window

### Setup

The final high-stat run used

- run root: `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_hi40000_notaylor/`
- twisted training family: `[6,6,3,1] * {2,4,6,8,10,12,16}`
- untwisted training sizes: `{8,16,24,32,48,64,96}`
- training MC: `n_traj = 40000`, `n_therm = 1000`
- blind holdout MC: `n_traj = 10000`, `n_therm = 200`

The blind datum was the smallest twisted multiple with the chosen boundary midpoint on-lattice and volume above `20000`:

- holdout lattice: `[144,144,72,24]`
- holdout volume: `22464`
- measured connected correlator: `0.28278461 +/- 0.00503771`
- holdout production runtime: about `54.54 s`

### Holdout Result

The blind-holdout comparison did **not** pick one universal best extrapolant across both branches.

For the twisted branch:

- best model: `power_free`
- holdout z: `-0.4325`
- training chi2/dof: `1.0090`

For the untwisted branch:

- best model: `pade21`
- holdout z: `-0.3269`
- training chi2/dof: `0.7234`

The untwisted `power_free` fit was still acceptable at the blind point (`z = +0.5831`), but `pade11` over-predicted badly (`z = +3.2404`).

### Interpretation

- Both branches can extrapolate to the held-out twisted datum within `|z| < 1`
- The preferred asymptotic ansatz is branch-dependent for this observable and this training window, not universal
- A lower-stat predecessor under `acute456_literal_sizes_holdout20k_hi20000/` landed on the same holdout geometry but gave noisier model rankings, so the `hi40000_notaylor` rerun is the one worth trusting

## 2026-05-25: Single-Displacement Holdout Mode for Large Blind Points

### Code Change

`K_from_continuum/src/ising_tri_twisted_parallelogram.cc` now accepts

- `--single_disp_m`
- `--single_disp_n`

and, when those are present, builds only the requested displacement map instead of the full `N_cell x N_cell` add table.

That matters because the full add table is the memory bottleneck for large blind-holdout points.

### Workflow

I used that new mode in the `responsible_method_tests/stupid_method/iso111_base4_size4to100_step4/` workflow:

- config: `responsible_method_tests/configs/raw_manifold_fss_iso111_base4_size4to100_step4_20260525.json`
- run root: `responsible_method_tests/stupid_method/iso111_base4_size4to100_step4/`

### Results

For the quick `L = 200` single-point holdout at point `pt002` (`m = 0`, `n = -100`):

- measured connected correlator: `0.27914208 +/- 0.00908915`
- the selected-model quick summary is written in `holdout_L200_pt002_quick_summary.json`

For the real `L = 1000` blind holdout at the corresponding scaled point (`m = 0`, `n = -500`):

- measured connected correlator: `0.17363803 +/- 0.01114078`
- output summary: `holdout_L1000_pt002_selected_models_logx.json`

At `L = 1000`, the asymptotic `power_free` model is clearly the best model on **both** branches at the blind point:

- untwisted `power_free`: `z = -0.4318`
- twisted `power_free`: `z = +0.8964`

while the Taylor and Padé models miss by much larger margins.

### Interpretation

- Single-point blind holdouts at volumes that would be awkward or impossible for the full all-to-all add-table workflow are now practical
- That gives a much cleaner direct test of continuum extrapolants against genuinely out-of-sample large lattices

## 2026-05-25: Acute 4-5-6 Beta-Finder Holdout Rerun Stopped Midway

### Context

I started a production rerun of the acute `4-5-6` literal-size blind-holdout workflow using the new per-lattice `beta_c` finder path:

- script: `responsible_method_tests/scripts/run_acute456_pow2_blind_holdout.py`
- run root: `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/`
- training MC: `n_traj = 40000`, `n_therm = 1000`
- holdout MC: `n_traj = 10000`, `n_therm = 200`
- beta finder: `11 + 7 + 7 + 7` scan points with `n_traj = 120` on the coarse pass and `240` on refinement passes, with jackknife enabled

The run was then intentionally killed before completion so I could switch focus to collaborator-facing project storytelling.

### What Completed Before The Stop

The production rerun did complete the held-out twisted lattice and most of the twisted training ladder:

- completed holdout beta scan and holdout production point
- completed twisted training beta scans and production points for
	- `[12,12,6,2]`
	- `[24,24,12,4]`
	- `[36,36,18,6]`
	- `[48,48,24,8]`
	- `[60,60,30,10]`
	- `[72,72,36,12]`
- started but did not finish the twisted `96x96_t48x16` beta scan
- did **not** reach any untwisted square beta scans
- therefore did **not** write the final blind-holdout summary JSON or any of the production PNG figures for this rerun

The last log line before termination was during the `96x96_t48x16` twisted scan:

- pass `2/4`, point `4/7`, `beta = 0.27024485`

### Partial Numerical Takeaway

The holdout beta scan itself was successful and peaked close to the exact isotropic triangular critical value `ln(3)/4 = 0.27465307`:

- holdout lattice: `[144,144,72,24]`
- holdout `beta_c = 0.27350859 +/- 0.00099267`
- holdout scan peak susceptibility: `0.0405321`
- holdout beta-scan wall time: about `908 s`

The completed twisted training pseudo-critical estimates showed a sensible finite-size trend toward the exact value:

- `12x12_t6x2`: `beta_c = 0.25699835 +/- 0.00175`
- `24x24_t12x4`: `beta_c = 0.26441420 +/- 0.0206`
- `36x36_t18x6`: `beta_c = 0.26985540 +/- 0.00216`
- `48x48_t24x8`: `beta_c = 0.26919364 +/- 0.00149`
- `60x60_t30x10`: `beta_c = 0.26905701 +/- 0.00108`
- `72x72_t36x12`: `beta_c = 0.27003786 +/- 0.00274`

So the beta-finder path looks operational on the twisted branch and the partial drift is physically reasonable, but the study is scientifically incomplete until the final twisted size, the full untwisted ladder, and the blind extrapolation fits are all finished.

### Artifacts Worth Reusing Later

- completed holdout beta summary JSON:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/beta_scans/holdout_twisted/acute456_holdout_smallbase024_m-108_n60_beta_scan.json`
- completed twisted beta summary JSONs currently exist for the six finished training sizes under:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/beta_scans/twisted/`
- completed MC production payloads currently exist under:
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/mc_runs/holdout_twisted/`
	- `responsible_method_tests/stupid_method/acute456_literal_sizes_holdout20k_betac_hi40000/mc_runs/twisted/`

Because the beta-scan summaries are cached per lattice, resuming later should not need to recompute the already-finished holdout and twisted sizes if the runner is pointed back at the same output tree without `--force`.

## 2026-05-27: Standard Iso111 Anchor-Ratio Landscape Note

I added a separate standard-workflow journal note for the iso111 pointwise target-score rerun with normalization turned on:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That note records the raw-score bias diagnosis, the new `anchor_ratio` scoring path, and the resulting ranking shift.

The short version is:

- raw scoring had favored `(0.8, 0.8)` over `(1.0, 1.0)`
- under anchor-ratio normalization, `(1.0, 1.0)` becomes clearly better than `(0.8, 0.8)`
- but the full anchor-ratio landscape still minimizes at `(1.1, 1.1)` rather than at isotropic `(1.0, 1.0)`

So the raw inversion was indeed a normalization-sensitive artifact, but normalization alone did not fully restore isotropic `1-1-1` as the unique optimum in the standard pointwise target score.

## 2026-05-28: Standard Boundary fit_method Follow-up Note

I added a substantial follow-up entry to the standard-workflow journal note:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That update records three concrete outcomes from the boundary-score session:

- the completed `iso111` quarter-point jackknife boundary summary and its surviving direction-ratio contribution
- the new `acute456` dual-target midpoint / quarter boundary sweeps, including separate large-target and small-target summary tables
- the new `fit_method` prototype, which fits the twisted and untwisted ladders separately before comparing continuum intercepts

The short version is:

- the new `fit_method` is not yet convincing; both the continuum-fit and blind-power variants fail to restore isotropic `1-1-1`, and the ratio-panel extrapolations look unstable
- the older direct boundary `chi2` summaries, including the direction-ratio channels, still look more promising and may contain real signal rather than pure noise

So the practical session conclusion is currently skeptical of the new fit-vs-fit score, but still interested in the old boundary `chi2` + direction-ratio line of attack.

## 2026-05-28: Acute456 Large-L Follow-up Hypothesis Note

I added another update to the standard-workflow note:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

This one records the current interpretation of the `acute456` small-target versus large-target behavior.

The short version is:

- the smaller twisted target may be looking better mainly because the present untwisted ladder only has to propagate from `64 -> 66`, which is almost a local holdout
- the large target currently asks the untwisted fit to bridge `64 -> 144`, which is a much longer extrapolation and may dominate the score more than the twisted/untwisted cutoff mismatch itself

The detailed note now includes a minimal next test:

- a `12`-job untwisted-only large-`L` campaign for the quarter-point `acute456` score
- shortlist = exact target, current large-target winner, current small-target winner
- new untwisted sizes = `96, 112, 128, 144`

The point of that campaign is to answer one narrow question before widening scope again:

- does the exact `acute456` target recover when the large-target comparison is made local near `L = 144`, instead of being extrapolated all the way from `L = 64`?

## 2026-05-28: Acute456 Ranking Diagnostics And Next-Step Freeze

I extended the same detailed standard-workflow note again:

- `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`

That latest update records the quarter-point `acute456` ranking experiments and the current decision about what **not** to do next.

The short version is:

- a mixed squared-rank aggregation across the quarter-point boundary metrics recovers the exact target only for the smaller twisted target
- ratio-only and correlator-only rank aggregations do **not** reproduce that exact-target win, so the effect is hybrid and heuristic rather than cleanly tied to one observable block
- same-size untwisted-to-twisted reweighting may be interesting later, but it does not solve the present `64 -> 144` size-extrapolation problem

So the ranking pass was useful as diagnosis, but the current decision is to stop inventing new scores and move to the cleaner experiment:

- a 12-job `acute456` quarter-point untwisted large-`L` shortlist
- sizes `96,112,128,144`
- couplings = exact target, current large-target winner, current small-target winner
- start at `n_traj = 40000`

## 2026-06-08: Background Fixed-Step 100k Run Launched From `(0.5, 5.0)`

### Why This Was Started

I wanted a live `100k` run advancing in the cloud while away from the editor, using the farther-start iso111 gradient-flow setup at

- start point `(r1, r2) = (0.5, 5.0)`
- `n_traj = 100000`
- `n_therm = 10000`
- `n_skip = 10`
- fit sizes `{32, 28, 24, 20, 16, 12, 8, 4}`
- `step_mode = normalized`
- `step_adaptation = fixed`
- `step_size = 0.025`
- `num_steps = 50`
- `min_gradient_snr = 1.0`

### Important Behavior Reminder

`responsible_method_tests/scripts/run_reweight_gradient_flow.py` now forces

- `accept_rule = "gradient_only"`

inside `main()`, regardless of the CLI value.

So this live fixed-step run is **not** using the old `direct_decrease` stopping logic from the earlier summaries. The active stopping conditions are now:

- gradient SNR falls below the threshold
- zero or clipped step
- requested step budget exhausted

This matters because the direct score at the proposal can worsen while the move is still accepted.

### Live Output Root And Process

Fresh detached run root:

- `responsible_method_tests/reweighting/iso111_gradient_flow_100k_start_r0p500_r5p000_fixed_steps50_bg_20260608_173456/`

Runner PID recorded in that root:

- `8116`

The run was launched with `nohup` from the workspace `.venv` Python and is intended to keep running without an attached terminal.

### Cache / Seeding Context

To avoid starting cold, I seeded the new background root by copying

- `_target/`
- `raw/`

from the earlier partial fixed-step root

- `responsible_method_tests/reweighting/iso111_gradient_flow_100k_start_r0p500_r5p000_fixed_steps50_livecov_20260608/`

I did **not** copy the old trajectory tables or PNGs, so the new root writes a clean live trajectory while still reusing any already-generated raw MC bundles.

### State At Handoff Time

The live files already existed when I stopped checking:

- `trajectory.tsv`
- `trajectory.png`
- `runner.pid`
- `nohup.log`

`nohup.log` was still empty at that moment, which is expected here because the script only prints the final JSON summary to stdout at the end.

The first written trajectory row was:

- current point `(0.5000000000, 5.0000000000)`
- proposal `(0.5237861689, 4.9923046656)`
- predicted score `5.2658145826e-02`
- direct proposal score `5.5600831763e-02`
- accepted = `True`

So the new background root had definitely started and was no longer just sitting at launch.

### Fast Resume Checklist

When picking this up later, the first things to inspect are:

- `ps -p 8116 -o pid,ppid,etimes,%cpu,%mem,cmd=`
- `responsible_method_tests/reweighting/iso111_gradient_flow_100k_start_r0p500_r5p000_fixed_steps50_bg_20260608_173456/trajectory.tsv`
- `responsible_method_tests/reweighting/iso111_gradient_flow_100k_start_r0p500_r5p000_fixed_steps50_bg_20260608_173456/trajectory.png`

Interpretation guide:

- if `summary.json` exists, the run finished cleanly
- if the PID is gone and `summary.json` is missing, inspect `nohup.log`
- if the PID is alive and `trajectory.tsv` is growing, just let it continue

This new `_bg_...` root is the active one to follow, not the older `...fixed_steps50_livecov_20260608/` root.

## 2026-06-08: Current-Stack CMA-ES Runner, Resume Support, And Acute456 Handoff

I added a new current-stack CMA-ES driver in

- `responsible_method_tests/scripts/run_reweight_cmaes.py`

The intent is to keep the optimizer on the same direct-MC / reweight-aware stack as the newer responsible-method scripts, rather than falling back to the old production optimizer code.

### What The New Runner Already Does

- optimizes the direct objective using the same current-stack evaluator path as `run_reweight_gradient_flow.py`
- writes `evals.tsv`, `generations.tsv`, `summary.json`, and `trajectory.png`
- renders the learned Gaussian covariance as `1 sigma` and `2 sigma` ellipses so the search distribution can be inspected generation by generation
- supports `--save-frames` for a per-generation history if needed
- now defaults its output root into `responsible_method_tests/reweighting/`

### Iso111 Smoke Validation

The first cheap smoke run succeeded at

- start `(r1, r2) = (0.5, 5.0)`
- target = iso111 square-untwisted current-stack target
- cheap budget: `target_size = 16`, fit sizes `{12, 8, 4}`, `n_traj = 200`, `n_therm = 50`, `n_skip = 1`

Smoke root:

- `responsible_method_tests/reweighting/iso111_cmaes_smoke_farstart_r0p500_r5p000_ntraj200_20260608/`

Smoke outcome:

- `18` direct evaluations
- `3` generations
- best score `0.004411952874`
- best point `(2.2636599282860415, 4.979669551335301)`

### Iso111 10k Far-Start Run

The first full current-stack iso111 run used

- start `(0.5, 5.0)`
- target size `64`
- fit sizes `{32, 28, 24, 20, 16, 12, 8, 4}`
- `n_traj = 10000`, `n_therm = 1000`, `n_skip = 10`
- `popsize = 8`

Run root:

- `responsible_method_tests/reweighting/iso111_cmaes_10k_farstart_r0p500_r5p000_total_score_pop8_eval120_20260608/`

The first `120`-evaluation run finished at `15` generations with

- best score `0.00044441459766`
- best point `(4.339072827435914, 4.172406504258914)`
- final mean `(3.30829212930538, 3.1153852978291643)`

That was still in the wrong basin, so I added true in-place resume support.

### Resume Support

`run_reweight_cmaes.py` now supports resuming a prior run root by reconstructing the CMA-ES state from

- `summary.json`
- `evals.tsv`
- `generations.tsv`

The new CLI path is

- `--resume-root <existing_root>`
- `--additional-gens <n>`

One replay bug mattered here: the selected offspring must be replayed in score-ranked order, otherwise the weighted CMA update is wrong even if the same set of points is present. That bug is fixed.

Resume was validated first on a copied smoke root and then on the real far-start iso111 run.

### Iso111 Continuation Results

After two successful resumptions, the same far-start iso111 run reached `35` generations / `280` evaluations.

Best-so-far progression:

- `15` generations: best `(4.339072827435914, 4.172406504258914)`, score `4.4441459766e-04`
- `25` generations: best `(1.738291649523858, 1.4716631387912569)`, score `1.2049723474e-04`
- `35` generations: best `(1.0531621969384206, 0.9342590147317011)`, score `4.8576189672e-05`

Final `35`-generation state:

- final mean `(1.0351051708550236, 0.9630145677144921)`
- final sigma `0.43798858407914937`

So the important scientific result is that the current-stack CMA-ES did eventually leave the far-start wrong basin and enter the neighborhood of the true iso111 target, but it needed more generations than the initial `120`-evaluation run allowed.

### Acute456 / 4-5-6 Handoff State

The next requested step was a detached `35`-generation CMA-ES run using the acute `4-5-6` target, but the original iso111 path could not truthfully do that because it depends on the square-untwisted `_run_point(...)` machinery in `test_geometry_match_grid_interpolation.py`.

To keep the work on the current stack, I started extending `run_reweight_cmaes.py` with a separate

- `--target-mode acute456`

path that reuses the existing acute456 helper code already present in the repo, especially

- `responsible_method_tests/scripts/run_acute456_pow2_blind_holdout.py`
- `responsible_method_tests/scripts/plot_standard_acute456_center_fss.py`

The new acute456 mode is meant to

- simulate untwisted families at candidate `(r1, r2)`
- compare them against the current acute456 twisted reference target
- score candidates with the existing aggregate target-score machinery instead of the iso111 square objective

Important handoff status:

- the acute456 mode has been patched into `run_reweight_cmaes.py`
- imports and parser wiring are in place
- file-level syntax checks passed
- but the acute456 path was **not yet smoke-tested** at handoff time
- the detached `35`-generation acute456 run was therefore **not launched yet**

So the next concrete step, before trusting any `4-5-6` detached production run, is:

1. run a focused cheap smoke test with `--target-mode acute456`
2. confirm that the family directories, per-size outputs, and scalar aggregate score are all written correctly
3. only then launch the detached `35`-generation run into `responsible_method_tests/reweighting/`

## 2026-06-09: Acute456 Detached CMA-ES Launch, Live-Refresh Patch, And Overnight Handoff

The `acute456` current-stack path has now been smoke-tested and the detached production-style run has been launched.

### Smoke Validation

Fresh smoke root:

- `responsible_method_tests/reweighting/_acute456_smoke_cmaes_20260609/`

That focused check confirmed that `run_reweight_cmaes.py --target-mode acute456`

- reads the current large twisted acute456 target
- writes per-candidate untwisted family directories
- produces `evals.tsv`, `generations.tsv`, `summary.json`, and `trajectory.png`
- returns a scalar aggregate `total_score` without falling back to the old production optimizer stack

### Detached Run Launched

Detached run root:

- `responsible_method_tests/reweighting/acute456_cmaes_10k_start_r3p000_r4p000_total_score_pop8_eval280_20260609/`

Launch configuration:

- `target_mode = acute456`
- start `(r1, r2) = (3, 4)`
- bounds `r1, r2 >= 0.5`
- `n_traj = 10000`, `n_therm = 1000`, `n_skip = 10`
- `popsize = 8`
- `max_evals = 280`
- `save_frames = true`

The run was left active in a detached terminal and was still writing output at the time of this note.

### Current Checkpoint At Handoff

At the checkpoint used for this handoff, the on-disk state had reached `2` completed generations / `16` completed evaluations.

Best-so-far progression visible in the live TSVs:

- generation `1`: best point `(2.8808346557, 4.5967903410)`, best score `3.7200530604e-01`
- generation `2`: best point `(3.6773289378, 5.1382313644)`, best score `1.5025178748e-01`

Current generation distribution row at that checkpoint:

- mean `(3.6096304676, 4.7899528540)`
- sigma `0.8716234865`

Primary live artifacts to inspect tomorrow:

- `responsible_method_tests/reweighting/acute456_cmaes_10k_start_r3p000_r4p000_total_score_pop8_eval280_20260609/evals.tsv`
- `responsible_method_tests/reweighting/acute456_cmaes_10k_start_r3p000_r4p000_total_score_pop8_eval280_20260609/generations.tsv`
- `responsible_method_tests/reweighting/acute456_cmaes_10k_start_r3p000_r4p000_total_score_pop8_eval280_20260609/trajectory.png`
- `responsible_method_tests/reweighting/acute456_cmaes_10k_start_r3p000_r4p000_total_score_pop8_eval280_20260609/launch.log`

### Live-Trajectory Update Note

During the run review I found that the original live plotting path only rewrote `trajectory.png` after a full CMA-ES generation completed. For acute456 that makes the figure appear frozen for long stretches because each generation evaluates `8` candidates and each candidate runs the full six-size lattice ladder.

`responsible_method_tests/scripts/run_reweight_cmaes.py` has now been patched so that new runs rewrite

- `evals.tsv`
- `generations.tsv`
- `trajectory.png`

after each completed evaluation, not only after each completed generation.

That per-evaluation refresh behavior was validated separately under:

- `responsible_method_tests/reweighting/_acute456_live_refresh_validate_20260609/`

Important caveat: the already-running detached acute456 production run was launched **before** this patch. So that specific run will continue to redraw only at generation boundaries. Future runs from the patched script should show genuine mid-generation live updates.

### Practical Overnight Expectation

Nothing in the current detached run is waiting for interactive input. So the run should continue on its own provided the codespace / container stays alive and nobody explicitly kills the terminal. That is an environment caveat rather than an algorithmic blocker.
