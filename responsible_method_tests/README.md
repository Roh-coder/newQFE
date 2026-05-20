# Responsible Method Tests

This directory is the staging area for a careful continuum-level validation of
three nominally equivalent ways of obtaining the all-to-all two-point function
of the 2D affine Ising model on a torus. The immediate goal is not to optimize
anything, but to establish a clean scientific baseline: when the geometry is
matched, do the three constructions converge to the same correlator manifold,
and when the geometry is perturbed, are the resulting differences large enough
to be resolved by a practical score?

## Scientific Motivation

At criticality, the long-distance spin-spin correlator on a torus is fixed by
the continuum conformal field theory up to the torus shape and an overall
normalization. In the present project, the torus shape is represented either by
its modular parameter `tau` in the analytic description or by a lattice
fundamental cell built from primitive cycles in the numerical descriptions. If
two constructions realize the same continuum torus, then after finite-size and
normalization effects are removed, they should define the same smooth
all-to-all correlator surface over the unit cell.

That statement is the foundation of the broader program. If it is true in a
controlled validation study, then twisted lattices and untwisted anisotropic
lattices can both be used as faithful lattice realizations of the same
continuum geometry. If it fails, then any later optimization or inference based
on matching correlators across geometries is on uncertain ground.

The validation problem therefore has two parts:

1. Equivalence: verify that all three constructions agree in the continuum
	 limit when they are tuned to represent the same torus.
2. Discriminability: verify that small but intentional deformations of either
	 the couplings or the twist move the correlator far enough in function space
	 that a reasonable score can tell the difference.

## The Three Constructions

### 1. Analytic modular-space solution

This is the continuum benchmark. For a torus with modular parameter `tau`, the
critical Ising two-point function can be written in terms of Jacobi theta
functions. In this route there is no finite-size extrapolation in principle:
once `tau` and the operator normalization convention are fixed, the correlator
is a continuum object.

Operationally, this construction provides the target manifold against which the
lattice realizations should converge. It also provides a clean way to generate
directional or full-cell correlator profiles at arbitrary points in modular
space without Monte Carlo noise.

### 2. Twisted lattice construction

In the twisted-lattice route, the torus is encoded directly in the lattice
fundamental domain through integers `(Lx, Ly, Tx, Ty)`. The twist changes which
sites are identified across the periodic boundaries and therefore changes the
shape of the effective torus even before any coupling anisotropy is introduced.

At criticality, the full all-to-all lattice correlator should approach the same
continuum manifold as the analytic modular solution for the corresponding
torus. This route is geometrically transparent, but care is required when
comparing directions or cycles: the physically corresponding sides of the torus
need not coincide with the raw lattice-axis labels.

### 3. Untwisted `L x L` lattice with anisotropic couplings

In the untwisted route, the boundary conditions remain simple but the couplings
`(k1, k2, k3)` are made anisotropic. The idea is that the anisotropy deforms
the effective continuum metric so that an apparently square lattice can realize
the same continuum torus as a twisted isotropic lattice.

This construction is the most attractive operationally because it avoids the
bookkeeping of twisted boundaries, but it is also the one that most needs to be
validated. The central claim is that there exists a map from anisotropic
couplings to continuum shape such that the continuum-limit all-to-all manifold
matches the modular and twisted descriptions.

## What Must Be Shown

The object of comparison is not a single line cut or a handful of correlation
lengths, but the full all-to-all correlator manifold over the torus unit cell.
Concretely, each numerical realization produces a table of correlator values on
all lattice separations, which can be mapped into normalized cell coordinates
and interpolated onto a common domain. The modular solution provides the same
kind of object directly in the continuum.

The main scientific claim to test is:

For matched geometry and matched critical tuning, the three constructions
converge to the same correlator manifold in the continuum limit, modulo the
expected overall normalization convention and controlled finite-size errors.

That claim splits into several concrete checks:

- Full-manifold agreement after alignment on a common normalized cell.
- Agreement of line cuts along the physically corresponding primitive cycles.
- Stable continuum extrapolation with increasing lattice size.
- No hidden dependence on whether the geometry is encoded by boundary twist or
	by coupling anisotropy.

## Why This Is Nontrivial

Several effects can fake agreement or fake disagreement if they are not treated
carefully.

First, critical tuning matters. The modular solution is already at continuum
criticality, while the lattice constructions require a finite-size estimate of
the pseudo-critical point and then a continuum extrapolation. If the beta value
is allowed to jitter independently at each lattice size, the resulting score
surface can become artificially rough. This is why the continuum workflows in
this repository use shared-beta extrapolation strategies across size families.

Second, cycle labels matter. Two tori can represent the same continuum shape
while assigning different lattice directions to the same physical side. Any
comparison that matches directions by index instead of by physical cycle can
generate an apparent disagreement that is purely kinematic.

Third, short-distance agreement is not enough. Many costs are dominated by the
near-origin region where all correlators are large and smooth. A successful test
must confirm agreement across the whole manifold, not just at the smallest
separations.

Fourth, a score is only useful if it separates the truth from nearby wrong
answers. It is not enough that the matched geometry gives a small score; the
wrong geometry must give a measurably larger one.

## Proposed Validation Program

### Stage A. Establish continuum equivalence

For a small set of representative torus shapes:

- Generate the analytic modular correlator manifold.
- Generate twisted-lattice all-to-all correlators across a size family and
	extrapolate them to the continuum.
- Generate untwisted anisotropic-lattice all-to-all correlators across the same
	size family and extrapolate them to the continuum.
- Align all three on the same normalized cell coordinates.
- Compare surfaces, line cuts, and continuum-fit residuals.

The minimal deliverable from this stage is a figure set showing that the
matched manifolds collapse onto one another within uncertainty bands.

### Stage B. Probe discriminability under controlled deformations

Starting from a validated matched point, introduce deliberate offsets:

- Perturb the anisotropic couplings away from the tuned point on the untwisted
	lattice.
- Perturb the twist parameters away from the matched torus on the twisted
	lattice.
- Optionally perturb both in directions that preserve one gross observable
	while changing the full manifold, to test whether the score is genuinely
	manifold-sensitive.

The aim is to measure how rapidly the score rises away from the matched point
and whether that rise dominates the statistical and extrapolation noise.

## Score Criteria

The score should be judged by two standards: physical correctness and practical
separation power. A good score is small for the matched continuum manifolds,
stable under continuum extrapolation, and clearly larger for nearby but wrong
geometries.

Candidate score families include:

- Pointwise manifold mismatch on a common interpolated grid, with uncertainty
	weighting where available.
- Directional mismatch along the primitive cycles and the third short cycle.
- Amplitude-robust shape scores, for example comparing normalized correlators
	or centered log-profiles so that the test is not driven only by an overall
	scale factor.
- Z-scores or related significance measures that compare matched and deformed
	cases relative to their combined uncertainty.

For collaborator discussions, the key operational criterion is simple:

The matched construction should produce a score statistically consistent with
zero, while intentional coupling or twist deformations should produce scores
that are separated from zero and from the matched case by a comfortable margin.

## Provisional Benchmark Matrix

The first campaign should start from a small number of named benchmarks rather
than a broad scan. The purpose of the benchmark set is to separate three very
different questions: whether the infrastructure works at a trivial control
point, whether it works for a nontrivial but moderate torus deformation, and
whether it still works in a strongly anisotropic regime.

| ID | Benchmark | Lattice realization to match | Role in campaign |
| --- | --- | --- | --- |
| B0 | Equilateral control | Twisted family `(L,L,0,0)` vs untwisted anisotropic target `(r1,r2)=(1,1)` | Null test for alignment, normalization, and continuum extrapolation machinery |
| B1 | Quarter-twist benchmark | Twisted family `(L,L,L/4,L/4)` with `L` chosen as multiples of 4; untwisted anisotropic target to be tuned from continuum matching | First nontrivial modular shape close to existing quarter-twist workflow conventions |
| B2 | 4-5-6 stress benchmark | Twisted family `alpha*(13,16,-3,3)` and untwisted anisotropic target near `(r1,r2)=(5.0652,7.7429)` | Strong-shape stress test with an exact continuum target already identified in repo notes |

The benchmark order matters. B0 should be used to debug the mechanics of the
comparison. B1 should be the first real physics test of the equivalence claim.
B2 should only be treated as decisive after B0 and B1 are stable, because it
is the most likely case to confuse true physics failure with analysis or tuning
failure.

## Provisional Pass/Fail Criteria

These criteria are intentionally phrased in significance language rather than
raw-score language so that they remain meaningful as the score definition is
refined.

- Match consistency: for a matched benchmark, the final continuum-extrapolated
	pairwise comparison scores among modular, twisted, and anisotropic routes
	should be statistically consistent with zero, with `|Z_match| <= 2` on the
	primary whole-manifold score.
- Ranking: within a local deformation stencil around the matched point, the
	matched configuration should be the lowest-score point or statistically tied
	for lowest.
- Separation: at least one deliberate twist deformation and one deliberate
	coupling deformation should satisfy `Z_sep >= 3` against the matched case on
	the primary score. A mature benchmark should aim for `Z_sep >= 5`.
- Robustness: removing the smallest lattice size or switching between nearby
	continuum-fit variants should not flip the pass/fail conclusion.
- Cycle matching sanity check: when a benchmark involves nontrivial cycle
	relabeling, the physically paired cycle comparison should outperform the
	naive index-paired comparison in a clear and reproducible way.

## Starter Scaffold

This directory now contains a minimal first-pass campaign scaffold:

- `plans/benchmark_matrix.md` for the benchmark table and exit criteria.
- `plans/first_validation_campaign.md` for the first campaign sequence.
- `configs/campaign_template.json` for a machine-readable starter spec.
- `configs/README.md` for how to map the starter spec onto existing workflows.
- `scripts/README.md` for the intended analysis entry points.
- `results/README.md` for artifact layout conventions.

## Current Status: 2026-05-20

The `geometry_456` benchmark is now the most useful recovery point for this
directory. The main findings from the latest pass are:

- The large untwisted mismatch was caused by a side-ordering bug in the
	normalized-cell embedding, not by a failure of continuum equivalence. For the
	4-5-6 benchmark, twisted/reference cycles `0,1,2` map to physical sides
	`5,6,4`, while untwisted/test cycles `0,1,2` map to `5,4,6`. The untwisted
	embedding therefore has to use physical sides `5,6`, which in this code path
	means `embedding_cycles: [0, 2]` instead of the default `[0, 1]`.
- The embedding fix is implemented in `scripts/generate_pointwise_manifold_dataset.py`.
	The working `geometry_456` configs now carry `embedding_cycles: [0, 2]` in the
	untwisted branch.
- After that fix, the high-stat benchmark under
	`results/raw_manifold_fss_integer_multiples_high_stats_20260520/` shows both
	twisted and untwisted approaching the aligned modular manifold well. The
	post-fix per-scale residuals are small for both methods: twisted has
	`chi2/dof ~ 0.04-0.20` and untwisted has `chi2/dof ~ 0.056-0.363`. The earlier
	untwisted `chi2/dof ~ 18-22` was kinematic fallout from the wrong embedding.
- The main high-stat artifact set to reuse is:
	`geometry_456_modular_convergence.png`,
	`geometry_456_residual_scaling.png`,
	`geometry_456_continuum_vs_modular_manifolds.png`,
	`geometry_456_twisted_fss_sampled_points.png`, and
	`geometry_456_untwisted_fss_sampled_points.png` in
	`results/raw_manifold_fss_integer_multiples_high_stats_20260520/`.
- A dedicated quadratic refit config now exists at
	`configs/raw_manifold_fss_campaign_more_large_sizes_high_stats_taylor2_refit.json`.
	It keeps the same high-stat data but lowers `min_sizes_for_free_C` from `8`
	to `3`, so the existing `taylor2` branch actually runs on the six-size
	`geometry_456` family.
- The quadratic-refit outputs live in
	`results/raw_manifold_fss_integer_multiples_high_stats_taylor2_refit_20260520/`.
	In those `pointwise_continuum.dat` tables the `fit_mode` field is `taylor2`
	for both twisted and untwisted. In this mode the columns `A, B, C` mean the
	continuum intercept, the coefficient of `1/scale`, and the coefficient of
	`1/scale^2`; `C` is not a power exponent.
- The sampled quadratic-refit plots are
	`geometry_456_twisted_fss_sampled_points_taylor2.png` and
	`geometry_456_untwisted_fss_sampled_points_taylor2.png` in the same refit
	directory. For the representative points checked so far, the quadratic-fit
	continuum intercepts remain close to the modular target within current
	uncertainties.
- Operational note: `generate_pointwise_manifold_dataset.py` resumes only within
	the results tree selected by `run.tag`. If the goal is to recompute only the
	derived continuum products under a new tag, first clone the finished results
	tree into the new tag directory (for example with `rsync -a --delete`) before
	re-running the generator; otherwise it will start a fresh MC campaign.
- Best next checks from here are:
	(1) regenerate the benchmark-level continuum-vs-modular manifold overlay from
	the quadratic-refit manifest,
	(2) drop the smallest lattice size and compare intercept stability between the
	linear-fallback and quadratic fits, and
	(3) script the ansatz-comparison step if this will be repeated on more than one
	benchmark.

## Success Conditions

This effort will be considered successful if we can demonstrate all of the
following:

- The modular, twisted, and anisotropic constructions agree on the continuum
	all-to-all manifold for the same target torus.
- The agreement survives more than one torus shape and is not confined to a
	special symmetric point.
- The preferred score ranks the matched point ahead of nearby deformations with
	clear statistical separation.
- The comparison protocol is explicit enough that collaborators can reproduce
	the result without relying on undocumented alignment choices.

## Immediate Deliverables For This Directory

This directory is intended to hold:

- Brief test plans for each target torus.
- Small reproducible configs for twisted and anisotropic size families.
- Analysis scripts that align manifolds and compute candidate scores.
- Summary notes and plots showing both equivalence and discriminability.

## Working Principles

- Prefer a few well-understood benchmark geometries over broad scans at the
	start.
- Match physical cycles, not raw coordinate labels.
- Treat continuum extrapolation and score design as part of the physics test,
	not as post-processing details.
- Keep all validation outputs separate from production optimization runs.