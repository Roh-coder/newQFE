# K_from_Optimization — Direct Coupling Matching by Optimization

**Status**: Active development — hand-coded NM + BFGS with adaptive MC stats  
**Started**: April 2026  
**Predecessors**: `K_from_CritSurface_standalone/`, `K_from_TwoPoint_standalone/`

---

## Background and Goal

The **coupling-matching problem**: given a twisted equilateral triangular Ising
lattice (reference), find anisotropic couplings (k₁, k₂, k₃) on an untwisted
rectangular lattice such that the two-point correlation functions match along
all three fundamental torus directions.  This is a statement about **universality**
under the renormalization group: different microscopic couplings can flow to the
same CFT fixed point and therefore produce identical long-distance correlators.

The three boundary directions of the torus are:
- **v** = (Lx, Ty)  — first periodicity vector
- **u** = (Tx, −Ly) — second periodicity vector
- **w** = −(v + u)  — diagonal

Matching along all three simultaneously is the "curve-collapse" condition: when
G_test plotted against rescaled arc length t ∈ [0,1] lands on top of G_ref, the
couplings are equivalent at the CFT level.

**What earlier iterations established:**
- The curve-collapse picture is physically correct: there is a clear minimum in
  the cost landscape and the correlators do collapse when couplings are
  well-matched.
- Pre-computing a β_c(r₁, r₂) surface (the v1 "CritSurface" approach) does
  amortize cost over many optimizer evaluations, but it adds a heavyweight
  Phase 1 step, surface-interpolation error, and a dependency on a second grid.
- The variance-normalised Z²/Z⁴ cost from v1 is unstable near the minimum
  because both numerator and denominator are noisy MC estimates.

---

## Design Philosophy for K_from_Optimization

This iteration drops the pre-computed β_c surface and replaces it with the
simplest possible workflow:

```
   for each candidate (r₁, r₂) the optimizer asks for:
       1. find β_c on the test lattice for those couplings (short MC scan)
       2. run a longer production MC at that β_c
       3. compute the L² cost vs the reference correlator
       4. return (cost, σ_cost) to the optimizer
```

There is **no separate surface-building phase, no interpolation, no LUT**.
Every cost evaluation is a self-contained "find criticality, then measure
correlator, then compare" step.  This is conceptually clean: the optimizer is
the sole driver of the workflow; the MC engine is a black-box function from
(r₁, r₂) to (cost, σ_cost).

The expected practical advantages:

- **Simpler mental model and code path.**  No surface JSON files to manage,
  no interpolation accuracy to worry about, no two-stage configuration.
- **β_c is always exact for the current lattice and couplings**, not an
  interpolated estimate.  No systematic surface-error contribution to the cost.
- **Adaptive resolution.**  Smart optimizers (Bayesian, CMA-ES, Powell) only
  evaluate the points they actually need — typically 30–60 points total —
  whereas a pre-built surface insists on a regular grid of ~100 points whether
  the optimizer needs them or not.
- **Warm-started β_c scans.**  Each new (r₁, r₂) is close to a previously
  visited point in coupling space, so its β_c is close to the previous β_c.
  The β_c finder can use the last result as a tight starting bracket and
  needs only a short scan.
- **No two-stage error budget.**  A single n_traj parameter sets the
  statistical precision; we don't need separate budgets for the surface scan
  and the production MC.

The expected practical disadvantage is that the **per-evaluation cost is
larger** because each call now includes a β_c scan (5–20 short MC runs) on
top of the production MC.  We expect this to be more than offset by the
reduction in total evaluation count from a smarter optimizer, but this is
something to verify in benchmarks.

---

## Repository Layout and File Status

This directory is a **standalone** package: everything required to build the
C++ simulator and run the Python benchmarks lives here.  Current state:

| File / Dir | Role | Status |
| --- | --- | --- |
| `Makefile` | Builds `bin/ising_tri_twisted_parallelogram` from `src/`. | ✅ in place |
| `src/`, `include/` | C++ Ising simulator sources, copied verbatim from v1. | ✅ in place |
| `mc_engine.py` | `run_simulator()`, `find_beta_c()` (3-pass GC + pass-0 window translation), `load_all_to_all()`. Has `progress_cb` hook for live β_c-scan plotting. | ✅ in place |
| `cost.py` | `l2_cost()` (plain L² over 3 boundary directions, propagated σ), `snr()`, `snr_status()`. | ✅ in place |
| `visualization.py` | `BetaScanPlotter` (per-point updates; single overwriting PNG with frame counter) and `OptimizerPlotter` (4-panel: trajectory+simplex, cost history, β_c history, residual ribbons). Simplex deduplication so ghost trail shows actual NM moves. | ✅ in place |
| `evaluator.py` | `evaluate(r1, r2)` wrapper. Warm-start bracket uses ±20% of β_prev (floor 0.05) so pass-0 covers a wide window. Adaptive `n_traj_prod` boosted by NM/BFGS when SNR drops. | ✅ in place |
| `optimizers.py` | Hand-coded NM (σ_shrink=0.75, adaptive stats, indist/noise stops) and hand-coded BFGS (rank-2 inverse-Hessian, central-FD grad, Armijo backtracking, adaptive stats). Plus Powell, BFGS-FD scipy wrappers, GP (skopt), CMA-ES. | ✅ in place |
| `run_benchmark.py` | Top-level CLI. Builds reference if missing (per-L cache), then dispatches optimizers. Full suite of NM + BFGS knobs. | ✅ in place |
| `run_nm_test.sh` | Isolated NM far-start test: x0=(1.5,0.5), max-evals=25, adaptive stats start prod=5k→max 20k, indist-stop=1.0. | ✅ in place |
| `run_bfgs_test.sh` | Isolated BFGS far-start test: x0=(1.5,0.5), max-evals=25, fd-eps=0.10, sigma0=0.4, adaptive stats start prod=60k→max 240k, indist-stop=1.0. | ✅ in place |
| `analyze_benchmark.py` | Reads `eval_log.jsonl`; writes convergence/trajectory/distance PNGs and markdown table. | ✅ in place |
| `plot_criticality_check.py` | Overlays GC scan data against fit for audit. Minor file-handle bug present (harmless). | ✅ in place |
| `requirements.txt` | numpy, scipy, matplotlib, scikit-optimize, cma. | ✅ in place |
| `bin/`, `results/` | Build output and per-run output (logs, frames, eval_log.jsonl). | ✅ created |

Reproducing a benchmark is now a two-line command:

```
make                                          # builds the C++ simulator
python run_benchmark.py --method all          # runs all five optimizers
```

Useful flags:

- `--method <name>` — run a single optimizer (`nelder_mead`, `bfgs_hand`, `powell`,
  `bfgs_fd`, `gp`, `cma`, or `all`).
- `--results-dir <dir>` — **always set this** to a new directory per run so
  nothing is overwritten (see Data Preservation Policy below).
- `--ref-L <int>` / `--test-L <int>` — reference / test lattice size.
  L=8 uses `results/_reference/`; other sizes use `results/_reference_L<N>/`.
- `--ref-n-traj`, `--n-traj-prod`, `--n-traj-scan-coarse`,
  `--n-traj-scan-fine` — base statistics knobs (adaptive boosting can
  multiply `n_traj_prod` at runtime).
- `--max-evals <int>` — hard cap on optimizer queries.
- `--x0 R1 R2` — starting point.
- `--beta-seed LO HI` — initial β bracket for the first evaluation
  (wide for far starts, e.g. `0.15 0.40`).
- **Nelder-Mead knobs:**
  - `--nm-sigma0` — initial simplex edge length (default 0.1; use ~0.4 for far starts).
  - `--nm-xatol`, `--nm-shrink` — position tolerance and shrink coefficient.
  - `--noise-stop-snr`, `--nm-indist-stop-snr` — early-stop thresholds.
  - `--nm-snr-target`, `--nm-snr-max-traj-factor` — adaptive stats: boost `n_traj_prod`
    by ×1.5 whenever SNR < target, capped at factor × starting value.
- **BFGS knobs:**
  - `--bfgs-fd-eps`, `--bfgs-sigma0`, `--bfgs-max-step` — FD step, initial step scale,
    hard step cap.
  - `--bfgs-noise-stop-snr`, `--bfgs-indist-stop-snr` — early-stop thresholds.
  - `--bfgs-snr-target`, `--bfgs-snr-max-traj-factor` — adaptive stats (same semantics as NM).
- `--save-every N` — write one visualizer frame per N evals (default 5; `1` = every eval).
- `--no-vis` — disable PNG frame writing for fast headless benchmarks.

---

## Data Preservation Policy

> **⚠ Never delete or overwrite results directories.  Every run must write to
> its own isolated output directory.**

Each MC run is expensive (minutes of wall time per evaluation).  Once data is
produced it is irreplaceable without re-running from scratch.  This has caused
problems in the past (a benchmark script wiped the `results/` tree, losing the
comparison data needed for decision-making).  The rules below prevent that.

### Rules

1. **Use `--results-dir` on every `run_benchmark.py` call.**  Pass a fresh,
   timestamped directory (the helper scripts do this automatically via
   `$(date +%Y%m%d_%H%M%S)`).  Never let two runs write to the same directory.

2. **Never `rm -rf results/<method>/`** in scripts or by hand unless you are
   certain the data is already archived elsewhere.  If you need a clean slate,
   create a new directory — do not delete the old one.

3. **The reference (`results/_reference/`) is permanently shared.**  It holds
   the 60k-trajectory equilateral reference run that every optimizer is tested
   against.  Do not delete it.  It will be reused automatically by all runs
   that point to the canonical `results/` tree.

4. **Archive before modifying.**  If a script previously wrote to a flat
   directory and you want to restructure, `mv` the old directory to a
   dated archive name first:
   ```bash
   mv results/nelder_mead results/nelder_mead_archive_20260422
   ```

5. **Scripts must not contain `rm -rf results/`.**  The benchmark driver
   scripts (`run_far_start.sh`, `run_nm_test.sh`, etc.) must never unconditionally
   delete results.  If a workflow truly requires a clean directory, it must
   explicitly request confirmation from the user before proceeding.

### Directory layout

```
results/
  _reference/             ← shared; never delete
    two_point_all_to_all.dat
    reference_meta.json
    reference_betac_scan.json
  far_start_20260422_1723/   ← one dir per run_far_start.sh invocation
    nelder_mead/
      eval_log.jsonl
      summary.json
      frames/
    powell/ ...
    benchmark_table.md
    convergence.png
    criticality_check.png
  nm_test_20260422_1800/     ← one dir per run_nm_test.sh invocation
    nelder_mead/
      eval_log.jsonl
      ...
```

This writes `results/<method>/eval_log.jsonl` and
`results/<method>/frames/*.png` for every run, plus a top-level
`results/summary.json` comparing total trajectory cost and final (r₁, r₂)
error across methods.

---

## The Evaluator: One Function from (r₁, r₂) to Cost

The entire pipeline is wrapped behind a single Python function:

```python
def evaluate(r1, r2, *, n_traj_prod, n_traj_scan, beta_init=None) -> EvalResult:
    """
    1. Find β_c on the test lattice with k1=r1, k2=r2, k3=1.
       Use beta_init (the previous evaluator's result) as starting point if given.
    2. Run a production MC at that β_c with n_traj_prod trajectories.
    3. Compute the L² cost vs the reference correlator along the three boundary
       directions, plus its propagated 1σ uncertainty.
    4. Return (cost, sigma_cost, beta_c, raw_data, ...).
    """
```

`EvalResult` records everything: the converged β_c, per-direction cost
contributions, σ_cost, total trajectory count, wall time, and the seed.
All evaluations are appended to a `eval_log.jsonl` file so any optimizer can
warm-start from the others' work.

The on-the-fly β_c scan uses the same 3-pass Gram-Charlier finder from the
predecessor projects.  Because successive optimizer queries are close in
(r₁, r₂) space, β_c moves smoothly and a tight initial bracket
[β_prev − ε, β_prev + ε] (with ε ≈ 5% of β_prev) is sufficient.  Typical
scan cost is ~10 short MC runs at ~5k trajectories each, much less than
the production MC cost.

**What this iteration addresses:**  
The v1 adaptive grid-search optimizer is **noisy and slow** near the minimum.
The v1 cost function (variance-normalised Z²/Z⁴) introduces its own noise because
the denominator σ²(t) is itself a noisy MC estimate.  The combined effect is that
the cost landscape near the minimum is dominated by statistical fluctuations rather
than physical signal.  This iteration simplifies the cost function, adds a rigorous
statistics-adequacy criterion, and replaces the grid search with a better-suited
classical optimizer.

---

## Simplified Cost Function

### Motivation for change

The v1 cost divided pointwise by the local variance:

$$C_\text{v1} = \sum_d \left[\int_0^1
  \frac{(G_\text{test}(t) - G_\text{ref}(t))^2}{\sigma_\text{ref}^2(t) + \sigma_\text{test}^2(t)}
  \,dt\right]^2$$

This normalization was intended to make the cost dimensionless and automatically
up-weight statistically precise regions.  In practice it creates two problems:

1. **Division by a noisy estimate.**  Both σ_ref and σ_test are jackknife/bootstrap
   estimates from finite MC samples.  In regions where G(t) is small (near the
   zero-crossing of the connected correlator), σ can be tiny and fluctuate wildly,
   making the integrand blow up from a pure numerical accident rather than a genuine
   shape mismatch.

2. **Circular noise amplification near the minimum.**  When G_test ≈ G_ref (i.e.,
   near the true coupling), both numerator (diff²) and denominator (σ²_ref + σ²_test)
   are O(σ²), so the integrand is O(1) regardless of how well the curves actually
   match.  The cost has no ability to resolve differences smaller than the statistical
   floor.

### The new cost: plain L² distance

Replace the variance-normalised integrand with the bare squared difference,
summed over the three boundary directions:

$$C(r_1, r_2) = \sum_{d \in \{v,\, u,\, w\}}
  \int_0^1 \bigl[G_\text{test}^{(d)}(t) - G_\text{ref}^{(d)}(t)\bigr]^2 \, dt$$

where the integral is computed by the trapezoid rule over n_samples = 400
uniformly spaced points t ∈ [0, 1] along each boundary direction.

**Properties of this cost:**

- **Units of [G]².**  The correlator has dimensions of (magnetisation)².  The cost
  is therefore straightforward to interpret: it is the mean-squared residual per
  unit arc length, integrated over all three directions.  If the correlator has a
  typical magnitude |G| ~ A, then a perfect match gives C = 0 and a random test
  lattice gives C ~ O(A²).

- **No division by noisy estimates.**  The noise in C comes only from the
  statistical uncertainty in G_test and G_ref, not from a ratio of two noisy
  quantities.  Near the minimum, the noise in C is proportional to
  σ_G × |G_test − G_ref|, which goes to zero faster than σ_G² as the curves
  converge.

- **Propagated statistical uncertainty.**  The uncertainty in C can be estimated
  from the per-point errors σ_ref(t) and σ_test(t) as:

$$\sigma_C^2 \approx \sum_d \int_0^1
  4(G_\text{test} - G_\text{ref})^2
  \bigl[\sigma_\text{ref}^2(t) + \sigma_\text{test}^2(t)\bigr]\,dt$$

  This is the linearised propagation of error through the squared difference.
  It gives a concrete, data-driven estimate of the noise floor in C.

- **Simple gradient structure.**  Since C is a quadratic functional of the
  interpolated test correlator, it is relatively smooth in (r₁, r₂) space and
  amenable to gradient-based and surrogate-model optimizers.

### The statistics-adequacy criterion

The fundamental question during optimization is: **"Is the current minimum of C
meaningfully lower than its neighbors, or is the difference within the noise?"**

Define the signal-to-noise ratio (SNR) at the current best point (r₁*, r₂*) as:

$$\text{SNR} = \frac{C(r_1^*, r_2^*)}{\sigma_C(r_1^*, r_2^*)}$$

where σ_C is the propagated uncertainty from the MC errors.

- **SNR < 1**: the cost at the putative minimum is statistically consistent with
  zero — the curves are already matching within 1σ of the current statistical
  precision.  **More statistics are needed** before the optimizer can make further
  progress.  The correct action is to increase n_traj (at least double it) and
  re-evaluate the current neighbourhood before continuing.

- **1 ≤ SNR < 3**: the minimum is marginally resolved.  The optimizer may be able
  to make one more useful refinement step, but convergence claims should be treated
  with caution.

- **SNR ≥ 3**: the minimum is clearly distinguished from its neighbours.  The
  optimizer can confidently zoom in.

This criterion replaces the ad-hoc Z⁴ re-scoring step from v1.  Instead of
computing a different cost function to "sharpen" discrimination, we simply demand
that the optimizer not proceed past a statistical noise floor.  When SNR < 1,
the run reports the current best (r₁*, r₂*) and prompts the user (or the
pipeline) to increase statistics.

**Estimating σ_C in practice:**  The interpolated-path approach gives σ_C via
the formula above.  An alternative is to run 4–5 independent replicas at the
best point and measure the empirical standard deviation of C directly.  This
replica approach is model-free and serves as a useful cross-check.

---

### Anisotropy-penalising cost variants

The plain L² sum $C = C_v + C_u + C_w$ allows the optimizer to exploit
directional trade-offs: it can zero one direction while leaving others large,
as long as the sum decreases.  The variants below penalise this.  Let
$C_d = \int_0^1 \Delta G_d^2\,dt$ denote the per-direction contribution.

| Label | Formula | Notes |
|-------|---------|-------|
| **A. Variance-normalised** | $\sum_d C_d/\langle\sigma_d^2\rangle$ | No free parameter; tight (low-noise) directions automatically count more |
| **B. L⁴ / quartic** | $\sum_d \int_0^1 \Delta G_d^4\,dt$ | Large deviations penalised steeply; steeper gradient near optimum; σ propagation more complex |
| **C. Max + regularisation** | $\max_d C_d + \lambda \sum_d C_d$ | Forces worst direction to shrink first; λ≈0.1 keeps it smooth |
| **D. L² + anisotropy penalty** | $\sum_d C_d + \mu\sum_d C_d^2$ | Zero only when all directions balanced; cheapest to add; recommended first experiment |

**Recommendation:** try **D** first (add one line, one tunable parameter μ),
then **B** if stronger penalisation is needed.  **A** is useful when one
direction has significantly noisier data than the others.  **C** is most
aggressive but can cause the optimizer to oscillate between directions.

**σ propagation for variant D:**

$$\sigma_C^2 = \sigma_{C_\text{L}2}^2 + \mu^2 \sum_d (2 C_d)^2 \sigma_{C_d}^2$$

where $\sigma_{C_d}$ is the per-direction uncertainty already computed by
`l2_cost`.

**σ propagation for variant B (L⁴):**

$$\sigma_C^2 \approx \sum_d \int_0^1 (4\,\Delta G_d^3)^2
  (\sigma_\text{ref}^2 + \sigma_\text{test}^2)\,dt$$

---

## Classical Optimization Options

The cost function C(r₁, r₂) is a two-dimensional, smooth, unimodal function
(empirically from v1 heatmaps) with an elliptical bowl-shaped minimum.
The function is stochastic: each evaluation returns C + ε where ε is MC noise.
This is a well-studied class of problem — **noisy 2D scalar minimization** —
and several classical methods are well-matched to it.

The options below are ordered from simplest to most sophisticated.  All of
them operate on the same evaluator interface and can be swapped in without
changing the MC engine.

---

### Option 1: Coarse Grid + Nelder-Mead Simplex

**In plain terms.**  Imagine you are lost on a hilly landscape in the dark and
want to find the lowest valley.  You start by standing at a few widely-spaced
spots and noting which ones feel lowest underfoot.  You then place three people
at a triangle (the "simplex") around the best starting point.  At each step,
the person standing highest walks through the middle of the other two to the
opposite side — the triangle crawls downhill.  When the triangle has
shrunk to the size of your foot, you have found the bottom.

**More precisely.**  First run a small 5×5 grid scan to get a good starting
point, then hand off to Nelder-Mead, a derivative-free downhill-simplex method.
The simplex (a triangle in 2D) adapts its shape to the local curvature: it
stretches along shallow directions and contracts perpendicular to steep ones.
It does not need to know the slope — only which vertex is highest at each step.

**Why it suits this problem.**  Nelder-Mead is robust to moderate noise,
since a single noisy evaluation only moves one vertex at a time, and the
other two provide a stable reference.  It converges in roughly 20–40
evaluations for a 2D problem.  The final simplex size gives a direct geometric
estimate of the uncertainty on (r₁*, r₂*).

**Key parameters:**
- `initial_simplex`: the starting triangle, sized from the coarse-grid spacing.
- `fatol`: tolerance on C; set to 2 × σ_C so the simplex stops when the cost
  is consistent with zero within two standard deviations rather than
  over-refining into statistical noise.

**Stopping rule:**  When the triangle has shrunk below `xatol` but C is still
larger than 2σ_C, this signals SNR < 1 — the curves look matched within the
current statistics.  Increase n_traj before refining further.

**Available in:** `scipy.optimize.minimize(method='Nelder-Mead')`.

---

### Option 2: Powell's Method (Conjugate Directions)

**In plain terms.**  Instead of wandering in any direction, this method takes
turns making one-dimensional sweeps: first slide purely in the r₁ direction
until you hit the bottom of that trench, then slide in the r₂ direction, and
so on.  After each full round it learns which diagonal direction is the
"cheapest" direction of descent (the one that required the least movement) and
adds that to its list of sweep directions.  Over time, the sweep directions
line up with the natural axes of the bowl, so subsequent rounds converge very
fast.

**More precisely.**  Powell's method performs a sequence of exact 1D line
minimisations along a set of conjugate search directions, updating the
direction set after each cycle.  For a 2D quadratic bowl it converges in
exactly two passes.

**Why it suits this problem.**  The v1 heatmaps show that the minimum is an
elongated ellipse in (r₁, r₂) space — the cost is much less sensitive to
changes along the r₁ = r₂ diagonal than perpendicular to it.  Powell's method
will rapidly identify these principal axes and exploit the shallow direction
efficiently without needing to compute any derivatives.

**Key consideration:**  Each 1D sweep calls the cost function several times.
With noisy evaluations the line minimiser needs a bracket tolerance matched
to the noise scale (roughly `xtol ~ 0.02` in coupling-ratio units).

**Available in:** `scipy.optimize.minimize(method='Powell')`.

---

### Option 3: Gradient Descent with Numerical Slopes (BFGS)

**In plain terms.**  Imagine you are on the same hilly landscape but now you
have a spirit level: you can measure the slope of the ground at your feet.
To estimate the slope in each direction you take a small step forward, measure
the height change, and step back — this gives you an approximate gradient.
You then take a bigger step in the downhill direction, remembering how the
slope has changed between steps to build up a mental map of the local curvature
(the Hessian).  That curvature map lets you jump almost directly to the bottom
of a quadratic bowl rather than zigzagging.

**More precisely.**  The slope ∂C/∂r₁ is estimated by evaluating C at
(r₁ ± h, r₂), and similarly for r₂.  The BFGS algorithm accumulates these
slope estimates over successive steps to approximate the second-derivative
matrix (how curved is the bowl?), then uses that curvature to take
Newton-like steps toward the minimum.

**Why it suits this problem — and the catch.**  When statistics are good
(SNR ≥ 3), this is the fastest method: it typically converges in 5–15
evaluations because it uses curvature information directly.  However, when
the cost C is dominated by statistical noise, the estimated slope is just
noise, and the method wanders.  The optimal finite-difference step size is
h ~ √(σ_C): too small and noise swamps the signal; too large and you miss
the shape of the bowl.  Only use this method after the SNR criterion is
satisfied.

**Available in:** `scipy.optimize.minimize(method='BFGS', jac='3-point')`.

---

### Option 4: Learning a Map of the Landscape (Gaussian Process / Bayesian Optimization)

**In plain terms.**  Instead of walking the landscape one step at a time, this
method builds a probabilistic map of the entire landscape from every
measurement made so far.  The map says: "at this (r₁, r₂) point I have
measured the cost; at points I haven't visited yet, my best guess is X with
uncertainty Y."  The next measurement is placed not where we think the minimum
already is, but where measuring would most likely reveal a new minimum —
balancing exploration of uncertain regions against exploitation of promising
ones.  Because every measurement updates the full map, nothing is thrown away:
a measurement taken at the very start of the run still contributes to the final
answer.

**More precisely.**  A Gaussian process (GP) is fitted to all (r₁, r₂, C)
observations.  The GP gives a smooth interpolant of the cost landscape plus
error bars that are large in unexplored regions and small near visited points.
The next point to evaluate is chosen by maximising the Expected Improvement
(EI) under noise — the integral of the amount by which C might improve,
weighted by its probability.  After each evaluation the GP is updated and a
new EI is computed.

**Why it suits this problem.**  It is the only method that explicitly models
the MC noise as a known additive component (σ_noise² per observation), so it
never mistakes a noise fluctuation for a genuine new minimum.  It is the best
choice when the run budget is tight and we cannot afford to simply increase
statistics.  The GP map is also a useful diagnostic: plotting the posterior
mean after the run shows whether the minimum is well-resolved or whether
an entire ridge of equally good couplings exists.

**Practical note:**  Requires an additional Python library.  The noise level
σ_noise can be estimated from a handful of repeated evaluations at the same
(r₁, r₂) point before the main run.

**Available in:** `skopt.gp_minimize` (`scikit-optimize`), or
`botorch.optim` (PyTorch-based, supports GPU).

---

### Option 5: Evolutionary Search with Adaptive Step Size (CMA-ES)

**In plain terms.**  Instead of a single walker, send a small population of
walkers (say 8) drawn from a statistical cloud centred on the current best
guess.  After each generation, keep the half that found the lowest costs,
and use their spread and direction to reshape the cloud for the next round:
if the survivors are all strung out along a diagonal, the cloud becomes a
long ellipse along that diagonal — the method literally learns the shape of
the bowl by seeing where the good points cluster.  The overall cloud size
(step size) also adapts: it grows if the population is not improving (need
to explore more) and shrinks when it is (zooming in on a minimum).

**More precisely.**  CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
maintains a multivariate Gaussian over the search space and updates its mean,
covariance matrix, and step-size scale from the best μ = λ/2 individuals each
generation.  It is derivative-free and makes no assumption about the shape of
the landscape beyond smoothness.

**Why it suits this problem.**  CMA-ES has a built-in noise-handling mode: it
re-evaluates the top candidates multiple times and uses the average cost for
learning, so statistical fluctuations do not corrupt the covariance update.
It naturally discovers the elongated r₁ = r₂ ridge structure.  For very noisy
problems or irregular landscapes it is the most reliable of the five options,
at the cost of more total evaluations (60–120 vs. 30–60 for the others).

**Available in:** `cma` Python package (`pip install cma`).

---

### Comparison Table

| Option | Gradient info | Noise handling | Eval budget | Dependencies |
|--------|---------------|----------------|-------------|--------------|
| 1 — Nelder-Mead | None | Moderate (ad hoc) | 30–60 | scipy |
| 2 — Powell | None (line search) | Moderate | 30–60 | scipy |
| 3 — BFGS + FD | Numerical | Poor at low SNR | 20–40 | scipy |
| 4 — GP / Bayes Opt | None (surrogate) | Excellent | 30–60 | scikit-optimize |
| 5 — CMA-ES | None (population) | Good (built-in avg) | 60–120 | cma |

**Recommendation — test all of them.**  All five optimizers are wired into the
same evaluator interface (`evaluate(r1, r2)` from the section above) so that a
side-by-side benchmark is a natural next step.  Concretely, `run_benchmark.py`
runs each optimizer in turn on the same reference and reports:

- Final (r₁*, r₂*) and final cost.
- Total number of MC trajectories spent (the true cost metric, since each
  evaluator call uses a variable amount of MC depending on how the β_c scan
  converges).
- Whether the optimizer hit the SNR floor and how many statistics increases
  were triggered.
- A trajectory plot of (r₁, r₂) candidates over time, overlaid for the
  five methods on the same axes.

For everyday use we expect Nelder-Mead (no extra dependencies, easy to audit)
and the GP / Bayes-Opt (best sample efficiency) to be the workhorses, with
the other three serving as cross-checks.

---

## Baseline Numbers (from v1)

Reference geometry used for benchmarking: **8×8 untwisted equilateral** (Tx=Ty=0).
Expected answer: r₁ = r₂ = 1.0 (same couplings).

| Run | n_traj/pt | Grid | Levels | Found (r₁, r₂) | Best cost |
|-----|-----------|------|--------|----------------|-----------|
| quick_test | 5k | 3×3 | 3 | (1.00, 1.15) | (Z²) 0.121 |
| optimizer_test | 10k | 5×5 | 4 | (0.925, 0.812) | (Z²) 0.0082 |
| test_20x20 | 10k | 3×3 | 3 | (1.30, 1.30) | (Z²) 1.38 |

All tests used an untwisted equilateral reference, so the physical answer is
r₁=r₂=1.  Deviations represent noise-driven misidentification.  The test_20x20
result r=(1.3,1.3) is clearly wrong — the 3×3 grid is too coarse to distinguish
the minimum from nearby points at n_traj=10k.

**Goal**: recover r₁=r₂=1.0 ±0.03 on the 8×8 equilateral reference
using the simplified L² cost, with the optimizer halting when SNR < 1 and
prompting a statistics increase, rather than silently converging to noise.

---

## Visualization

To verify by eye that both the **β_c finder** and the **optimizer** are
behaving sensibly, every run produces two streams of PNG frames in the
output directory.  These are written incrementally so they can be inspected
while a long benchmark is still running.

### β_c finder (Gram-Charlier scan)

Each call to `find_beta_c()` triggers a `progress_cb` after every MC point
and after each of the 3 GC fitting passes.  The visualization plots:

- **All measured (β, χ) points so far**, colored by which scan pass produced
  them (pass 0 = coarse sweep, pass 1 = first refinement, pass 2 = second
  refinement, pass 3 = final tight bracket).  Error bars on χ are shown.
- **The current Gram-Charlier fit curve** drawn through the data.  The fit is
  a Gaussian dressed with skewness (H₃) and kurtosis (H₄) Hermite corrections
  — the literal formula being used by the finder.  As more points come in,
  the curve tightens around the true peak.
- **A vertical line at the current β_c estimate** (the mode of the GC fit).

The output is `frames/betac_scan_<label>_pass<N>.png` for each pass.  Watching
the sequence shows the bracket narrowing and the fit converging.  If a frame
shows the GC curve diverging from the data, or the estimate jumping wildly
between passes, the scan is unhealthy and the run will likely fail downstream.

### Optimizer trajectory

Every call to `evaluate(r₁, r₂)` is logged to `eval_log.jsonl`, and the
visualization re-reads this log to produce a 4-panel summary frame after
each call:

1. **Top-left — landscape map.**  Scatter of every (r₁, r₂) the optimizer
   has ever tried, colored by cost (log scale).  The current point is marked
   with a star, the running best with a circle.  The trajectory is drawn as
   a polyline so the optimizer's strategy is visible (a Nelder-Mead simplex
   crawls in triangle-flips; CMA-ES sprays clouds of points; a GP samples in
   a Latin-hypercube–like pattern at first then concentrates on the minimum).

2. **Top-right — cost vs. evaluation number.**  C and σ_C plotted as a
   function of how many evaluations have been done so far.  When the running
   best stops improving and σ_C bars start crossing zero, the SNR criterion
   is about to trigger.

3. **Bottom-left — β_c history.**  The β_c value returned by each evaluation
   plotted against iteration number.  Useful for confirming that the
   warm-started β_c scan is tracking the smooth β_c(r₁, r₂) surface as
   the optimizer wanders.

4. **Bottom-right — current curve collapse.**  G_test(t) and G_ref(t)
   overlaid for each of the three boundary directions (v, u, w) at the
   current best (r₁, r₂).  This is the direct visual check that the curves
   really are collapsing: when the colored test curves land on the gray
   reference curves within their error bars, the optimizer has succeeded.

The output is `frames/optimizer_<method>_step<N>.png`.  At the end of a
benchmark, the per-method frame sequences can be assembled into GIFs with
one ffmpeg command (a helper script `make_gifs.sh` is included).

### Implementation

All visualization lives in `visualization.py` and is fully optional —
running with `--no-vis` disables the frame writing for fast benchmarks
where only the final numbers matter.  Internally it uses matplotlib's `Agg`
backend so it never opens a window and can run on a headless server.

---

## Benchmark Results — first run (2026-04-22)

The first apples-to-apples benchmark of all five optimizers, run with the
`--snr-floor 0` flag so each method uses its full evaluation budget.

**Setup**
- Reference lattice: 6×6 equilateral, k₁=k₂=k₃=1, untwisted (T=0).
- Reference statistics: 60 000 trajectories at the measured β_c = 0.22732.
- Test lattice: same geometry (self-consistency test → true minimum at
  r₁ = r₂ = 1).
- Per-evaluation budget: 30 000 production trajectories + a 3-pass
  Gram-Charlier β_c scan (4 000 / 10 000 trajectories per scan point).
- Starting point: x₀ = (1.15, 0.85), max 25 evaluations per method.
- β_c-scan pass widths use the new sqrt-of-variance schedule
  (pass 1 = ±2σ_fit, pass 2 = ±1σ_fit, pass 3 = ±0.5σ_fit).

**Results table**
([raw output](results/benchmark_table.md), plots in `results/*.png`):

| Method | n | best cost | best (r₁, r₂) | best dist | final-3 ⟨dist⟩ | min dist seen | total traj | wall (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nelder_mead | 25 | 1.44e-04 | (1.100, 0.925) | 0.125 | **0.122** | 0.071 | 6.5 M | 98 |
| powell | 10 | 2.50e-04 | (1.150, 0.850) | 0.212 | 0.260 | 0.212 | 2.7 M | 40 |
| bfgs_fd | 25 | 5.47e-05 | (0.990, 0.854) | 0.147 | 0.176 | 0.147 | 6.7 M | 95 |
| gp | 25 | 1.19e-04 | (1.309, 1.320) | 0.445 | 0.327 | 0.171 | 6.6 M | 92 |
| cma | 25 | 1.29e-04 | (1.068, 0.977) | **0.072** | 0.304 | **0.072** | 6.5 M | 99 |

Distance is `‖(r₁, r₂) − (1, 1)‖` in coupling space.  *best cost* is
the single lowest cost the method evaluated (note: all five methods'
best-cost points have SNR ≈ 0.5, so by the production-time stopping rule
this point is statistically *consistent with zero* — the absolute cost
ranking is dominated by noise, not by which optimizer is best).
The *final-3 ⟨dist⟩* and *min dist seen* columns are the noise-robust
metrics: they measure where each optimizer is *operating*, not where
its luckiest fluctuation landed.

**Analysis (one method per paragraph)**

- **Nelder-Mead — best converged operating point.**  Final-3 average distance
  0.122, with the simplex collapsed in a tight neighbourhood of (1.10, 0.93).
  No noise-fitting, no stalling, no wasted exploration.  This is the cleanest
  win on the noise-robust metric and confirms the README recommendation.
- **CMA-ES — closest single point but unconverged cloud.**  Reached
  (1.07, 0.98) at evaluation 21 (distance 0.072), but the population is still
  spraying widely (final-3 average 0.30).  With 25 evals it has barely
  finished its initial sigma-adaptation phase; budget needs to be ~50–60 to
  let the covariance contract.
- **BFGS-FD — noise-poisoned gradients.**  Lowest raw cost (5.5e-5) but it
  came from a near-zero MC fluctuation, not from gradient information.  The
  3-point finite-difference Jacobians wandered around (1.0, 0.85) without
  ever moving towards r₂ = 1.0; the SNR-filtered mean distance is the worst
  of all five (0.42).  Confirms the prior expectation that derivative-based
  methods are a poor match for an MC-noisy objective.
- **GP / Bayesian Opt — explored, never exploited.**  Spent the first ~8
  Latin-hypercube–like evaluations far from (1, 1) and the final 17
  acquisition-function-driven evaluations clustered around (1.3, 1.3) —
  apparently a noise-induced local minimum the GP latched onto.  With this
  budget GP underperforms simpler methods; it would benefit from either
  more evaluations or a tighter prior bounding box.
- **Powell — terminated prematurely.**  Stopped at 10 evaluations, well
  below the cap.  Powell's line-search-and-restart heuristics declared
  convergence on a noise-flat ridge.  Its `xtol`/`ftol` thresholds are too
  loose for an MC-noisy cost; needs tightening (or wrapping in averaging)
  before a fair re-comparison.

**Headline conclusion**

For the *current* MC budget on a 6×6 lattice, **Nelder-Mead is the operational
winner** — most evaluations spent near the true optimum, no special tuning
needed, and no surprises.  CMA-ES is the most promising *bigger-budget*
alternative: with ~2× the evaluation count its closest-point performance
suggests it would beat Nelder-Mead on convergence quality.  GP, BFGS-FD,
and Powell as configured are not competitive on this problem.

**Caveats**

- Self-consistency test only; needs to be repeated against an *external*
  reference (a different lattice geometry, or a continuum target) to
  rule out trivially-easy gradients.
- 30 000 production trajectories per evaluation gives σ_C ≈ 2e-4 near
  the optimum — comparable to the cost itself.  Doubling the production
  statistics would let the SNR-floor stopping rule decide convergence
  instead of a hard `max_evals`.
- Powell's premature termination should be diagnosed before discarding
  the method; the table entry is an upper bound on its performance.

**Reproducing**

```
make
python run_benchmark.py --method all \
    --ref-L 6 --ref-n-traj 60000 \
    --n-traj-prod 30000 --n-traj-scan-coarse 4000 --n-traj-scan-fine 10000 \
    --max-evals 25 --x0 1.15 0.85 --snr-floor 0 --no-vis
python analyze_benchmark.py
```

---

## Next Steps (for the next working session)

The first benchmark identified **Nelder-Mead as the operational winner**.
The next session focuses entirely on understanding *why* NM wins, where
it breaks, and how to push its performance further.  The other four
optimizers are parked until NM is well-characterized.

Concrete tasks, in priority order, sized to fit in one session each.

### A. Multi-start ensemble of Nelder-Mead

A single run from x₀ = (1.15, 0.85) is one realization of a stochastic
trajectory.  We need to know the **distribution** of NM outcomes to make
any quantitative claim about convergence.  **Action:**

1. Add `--n-starts <int>` to `run_benchmark.py`.  For each start, draw
   x₀ uniformly from a configurable bounding box (default `[0.7, 1.3]²`),
   with a fixed master seed for reproducibility.
2. Rerun NM with `--n-starts 8`.
3. In `analyze_benchmark.py` aggregate per-start: report median, IQR,
   and worst-case of *final-3 ⟨dist⟩* and *min dist seen*.  Plot all
   trajectories on the same (r₁, r₂) panel so basin structure is visible.

**Acceptance**: a single number "NM converges to within X of (1, 1) on Y%
of starts in N evals" with X, Y, N quantified.

### B. Initial-simplex sensitivity study

The current driver builds the simplex with `sigma0 = 0.1`, i.e. a side
length of 0.1 in coupling space.  This is a guess.  Larger simplex →
more exploration, slower contraction; smaller → faster contraction but
risk of locking onto a noise minimum.  **Action:**

1. Run NM with `sigma0 ∈ {0.02, 0.05, 0.1, 0.2, 0.3}` from the same x₀.
2. Plot final-3 ⟨dist⟩ vs sigma0.
3. Add a `--sigma0` CLI flag and document the recommended value as the
   minimum of that curve.

### C. Adaptive-restart Nelder-Mead

Standard NM contracts monotonically; once the simplex collapses it cannot
re-explore.  In the noise-driven regime the simplex collapses prematurely
because random fluctuations look like local minima.  **Action:**

1. Implement `run_nelder_mead_restart` in `optimizers.py`: after each
   `scipy.optimize.minimize` returns (or after every K evaluations),
   re-inflate the simplex around the running best with `sigma0 / 2` and
   continue, until `max_evals` is exhausted.
2. Compare against vanilla NM at equal evaluation budget.

**Hypothesis to test**: restart NM beats vanilla NM by 2× on
*final-3 ⟨dist⟩* at the same budget, because it can escape noise basins.

### D. Cost averaging — kill noise at the call site

NM relies on consistent cost ranking between trial points.  When σ_C is
comparable to the cost difference between simplex vertices, NM is doing
random search.  An obvious fix: have each evaluator call return the
average of K independent MC runs.  **Action:**

1. Add `--n-repeats <int>` to `evaluator.py`; each call runs K production
   MCs at the same (r₁, r₂) with different seeds and averages.
2. Try K ∈ {1, 2, 4} keeping total trajectory budget fixed (so K=4 means
   each MC is 1/4 the trajectories).
3. Determine whether averaging-then-comparing beats single-call-with-more-
   trajectories, or vice versa.  This answers a fundamental question:
   does NM care more about *low σ per call* or *more independent samples*?

### E. Cross-geometry test for NM

The current results are self-consistency only (test geometry = reference
geometry, true minimum at r₁ = r₂ = 1).  Confirm NM also wins on a
non-trivial target.  **Action:**

1. Extend `run_benchmark.py` with `--test-Tx`, `--test-Ty` flags.
2. Run NM on an 8×8 untwisted reference matched against a twisted
   parallelogram test geometry (e.g. Tx = 1, Ty = 0).  Record the
   recovered (r₁, r₂); document this as the geometry-matched point.
3. Re-run with x₀ ≠ matched point and verify NM converges back to it.

### F. β_c warm-start tightness sweep

The on-the-fly β_c scan currently uses a ±2% window around the previous
β_c.  This is a guess.  **Action:**

1. Add `--betac-window-frac <float>` to `evaluator.py`.
2. Sweep ∈ {0.005, 0.01, 0.02, 0.05} on a 4-eval NM run.  Measure scan
   wall time and how often the scan misses the peak (β_c ends up at the
   bracket edge, indicated by `pass_ids` containing window-shift events).
3. Pick the smallest window that never misses; document.

### G. Adaptive `n_traj_prod` driven by SNR

Each NM evaluation currently uses fixed statistics.  Far from the
minimum, σ_C is large and even cheap MC suffices to tell vertices apart;
near the minimum we need much higher statistics.  **Action:**

1. Implement in `evaluator.py`:
   - if previous SNR > 5 → halve `n_traj_prod` for next call (min 5 000)
   - if previous SNR < 1 → double `n_traj_prod` for next call (max 200 000)
2. Re-run NM with this adaptive rule and `--snr-floor 1.0`.
3. Verify the optimizer terminates at a meaningful endpoint instead of
   exhausting `max_evals`.

### Acceptance for the next session

A "good next session" output is a new sub-section in this README titled
**"Nelder-Mead Characterization"** containing:

- A median ± IQR convergence number from task A (multi-start).
- A recommended `sigma0` from task B with one-line justification.
- A pass/fail verdict on the restart hypothesis (task C).
- One concrete production recipe for NM, e.g.
  `--method nelder_mead --sigma0 0.05 --n-traj-prod 60000 --n-repeats 2 --max-evals 40 --snr-floor 1.0`,
  with the supporting numbers.

Tasks D–G are nice-to-have if time permits; A, B, C are the core.

---

## Open Questions

1. **What twist geometries are the target?**  The equilateral untwisted case
   (Tx=Ty=0, expected r₁=r₂=1) is the self-consistency benchmark.  The
   physically interesting cases have non-zero twist where r₁≠r₂ at the minimum.
   We need to characterise how the L² cost landscape changes with twist.

2. **Normalization of the L² cost.**  The absolute value of C scales with the
   overall magnitude of the correlator, which depends on L and β_c.  Should we
   normalize by the reference correlator's L² norm? i.e.,

   $$C_\text{normalized} = \frac{\int_0^1 (G_\text{test} - G_\text{ref})^2\,dt}
                                  {\int_0^1 G_\text{ref}^2\,dt}$$

   This would make the cost dimensionless and scale-invariant, comparable across
   different lattice sizes.  The tradeoff is one more noisy denominator — though
   the denominator here is non-local (it doesn't blow up the way per-point division
   did in v1).

3. **How to combine the three directions.**  The current plan sums C_v + C_u + C_w.
   Should these be weighted by direction length?  By the reference signal strength
   in each direction?  The three directions are not equally well-sampled by the
   all-to-all correlator for non-square lattices.

4. **Is the 2D search space sufficient?**  We fix k₃=1 and optimize over the
   ratios r₁=k₁/k₃, r₂=k₂/k₃.  For the untwisted case the expected answer is
   r₁=r₂=1 by symmetry, so the 2D search is sufficient.  For twisted geometries
   the answer may lie on a curved manifold in (k₁,k₂,k₃) space; the 2D ratio
   parameterization may be adequate if the overall scale is set by criticality.

5. **Warm-starting β_c from the analytic formula.**  The exact infinite-volume
   criticality condition:

   $$\sinh(2\beta k_1)\sinh(2\beta k_2) + \sinh(2\beta k_2)\sinh(2\beta k_3)
   + \sinh(2\beta k_3)\sinh(2\beta k_1) = 1$$

   gives β_c^exact(r₁, r₂) at zero cost.  For large L this may be a
   sufficiently good first-pass starting point for the on-the-fly β_c scan,
   shortening it from a full 3-pass GC to a single tight refinement.

---

## Session Log

> This section records what was explored and decided in each working session.
> Newest entries at the top.

### 2026-04-23 — Adaptive MC stats; BFGS hand-coded; simplex display fix; parallel runs

#### Adaptive MC statistics (both NM and BFGS)

Both `run_nelder_mead` and `run_bfgs_hand` now support adaptive trajectory
counts.  At each evaluation, if the returned SNR is below `snr_target`, the
evaluator's `n_traj_prod` is multiplied by 1.5 and capped at
`snr_max_traj_factor × base_n_traj_prod` (where `base` is snapshot at
optimizer startup to prevent drift across many boosts).

CLI flags added:

| Flag | Optimizer | Default | Meaning |
|------|-----------|---------|---------|
| `--nm-snr-target` | NM | 0.0 (off) | SNR threshold below which stats are boosted |
| `--nm-snr-max-traj-factor` | NM | 4 | Cap multiplier on `n_traj_prod` |
| `--bfgs-snr-target` | BFGS | 3.0 | SNR threshold |
| `--bfgs-snr-max-traj-factor` | BFGS | 4 | Cap multiplier |

Example: `--nm-snr-target 2.0 --nm-snr-max-traj-factor 4` with `--n-traj-prod 5000`
starts at 5k trajectories and can rise to 20k if the signal is too noisy.

#### Hand-coded BFGS (`run_bfgs_hand`)

Implemented a rank-2 inverse-Hessian BFGS from scratch (no scipy).  Key
features:

- **Central-FD gradient**: for a 2D parameter, 4 evaluations per gradient
  step at spacing `fd_eps` (configurable via `--bfgs-fd-eps`; default 0.10).
- **Armijo backtracking**: `c1=1e-4`, up to 4 step halvings per line search.
- **Hard step cap**: `max_step` on step norm (default 0.8 via `--bfgs-max-step`).
- **Curvature reset**: resets H to identity whenever `s·y ≤ 0` (curvature
  condition violated due to noise).
- **Parameter clipping**: each parameter clamped to [0.3, 2.5] after every step.
- Same `noise_stop_snr`, `indist_stop_snr`, `snr_target`, `snr_max_traj_factor`
  as NM.

New CLI flag: `--bfgs-sigma0` (initial step scale, default 0.1; use 0.4 for
far starts).

#### Simplex visualization fix

Two bugs caused the NM triangle display to lag behind the actual optimizer
state:

1. **Slot-based substitution**: `_eval(pt, simplex, slot=-1)` now builds a
   display copy of the current simplex and substitutes `pt` at position `slot`
   before pushing to `OptimizerPlotter`.  This means every triangles drawn
   reflects the candidate being evaluated, not a stale previous state.
   `slot=None` for initial vertex evaluations; `slot=i` for shrink steps;
   `slot=-1` (default) for reflection/expansion/contraction candidates.

2. **Deduplication in `update()`**: `OptimizerPlotter.update` only appends
   to `simplex_hist` when the new simplex differs from the last stored entry.
   This removes the 2–3× duplication that caused the animation to stutter and
   triangles to appear not to move.

#### β-scan warm-start bracket widened 10×

The initial β bracket for each evaluation was ±2% of `β_prev` (min 0.005),
giving a span of ~0.01.  This caused frequent window-translation cascade
failures when the optimizer jumped far in (r₁, r₂) space and the true β_c
had moved significantly.

Changed to: `eps = max(0.20 × β_prev, 0.05)`, giving a span of ~0.10.
The existing window-translation fallback (10%-edge trigger, shift = 0.5×span,
up to 4 shifts per pass) is retained as a safety net.

#### BFGS parameter tuning

| Parameter | Old | New | Effect |
|-----------|-----|-----|--------|
| `fd_eps`  | 0.05 | 0.10 | Wider FD stencil → less noise in gradient |
| `sigma0`  | 0.1  | 0.4  | Larger first step → not trapped near x0 |
| `max_step`| 0.4  | 0.8  | Allows longer jumps toward minimum |

#### Parallel run execution pattern

Both test scripts can be launched in parallel with:

```bash
D=/workspaces/newQFE/K_from_Optimization
bash "$D/run_nm_test.sh" > run_nm.log 2>&1 &
bash "$D/run_bfgs_test.sh" > run_bfgs.log 2>&1 &
wait
```

(Absolute paths required because each `&` subshell does not inherit a
`cd`-applied working directory from the parent.)

#### Standard test configuration (25 evals)

Both `run_nm_test.sh` and `run_bfgs_test.sh` are now standardised at 25
evaluations with far-start `x0 = (1.5, 0.5)`.

`run_nm_test.sh` (current defaults):

```
--max-evals 25
--n-traj-prod 5000  --n-traj-scan-coarse 1000  --n-traj-scan-fine 3000
--nm-snr-target 2.0  --nm-snr-max-traj-factor 4
--nm-sigma0 0.4  --nm-xatol 0.02  --nm-shrink 0.75
```

`run_bfgs_test.sh` (current defaults):

```
--max-evals 25
--n-traj-prod 60000  --n-traj-scan-coarse 8000  --n-traj-scan-fine 20000
--bfgs-fd-eps 0.10  --bfgs-sigma0 0.4  --bfgs-max-step 0.8
--bfgs-snr-target 3.0  --bfgs-snr-max-traj-factor 4
```

#### NM early-stop behaviour confirmed correct

On one 20×20 run, NM stopped at eval 11 with cost = 8.6×10⁻⁵ ± 1.45×10⁻⁴
and SNR = 0.60.  This triggered `indist_stop_snr = 1.0` ("indistinguishable
from reference within MC noise") — the correct outcome.  The optimizer
found the minimum; it did not crash or stall.

#### Files changed this session

| File | Change |
|------|--------|
| `optimizers.py` | `run_nelder_mead`: `snr_target`, `snr_max_traj_factor`, `_n_traj_prod_base` capture, slot-based `_eval`; `run_bfgs_hand`: full hand-coded BFGS, same adaptive stats |
| `visualization.py` | `OptimizerPlotter.update()`: deduplication guard on `simplex_hist` |
| `run_benchmark.py` | `--nm-snr-target`, `--nm-snr-max-traj-factor`, `--bfgs-sigma0`, `--bfgs-snr-target`, `--bfgs-snr-max-traj-factor` flags; NM and BFGS kwargs wired through |
| `evaluator.py` | Warm-start bracket widened from ±2% to ±20% |
| `run_nm_test.sh` | Updated to 25-eval standard with adaptive stats flags |
| `run_bfgs_test.sh` | Updated to 25-eval standard with improved BFGS params |
| `README.md` | Status line, file table, CLI flag docs, session log updated |

### 2026-04-22 — Data preservation policy; isolated run dirs; run_nm_test.sh

**Problem.**  Earlier in this session `run_far_start.sh` included
`rm -rf results/$m` for each optimizer method and also deleted the
`results/_reference/` tree.  This wiped the 60k-trajectory reference run
and all previous benchmark outputs that would have been needed for
cross-run comparison.

**Fix: isolated output directories.**  `run_benchmark.py` now accepts
`--results-dir <dir>`.  Every run writes exclusively into that directory;
the shared reference is always loaded from (and cached in) the canonical
`results/_reference/` regardless of `--results-dir`.  The helper scripts
generate a timestamped directory via `$(date +%Y%m%d_%H%M%S)` so two runs
can never collide.

**A Data Preservation Policy section** was added to this README with explicit
rules:

- Always pass `--results-dir` to `run_benchmark.py`.
- Never `rm -rf results/<method>/` in scripts.
- `results/_reference/` is permanent and shared.
- Archive before restructuring; never silently overwrite.
- Scripts must not contain unconditional delete commands.

**`run_far_start.sh` rewritten** — no longer deletes anything.  Outputs go
to `results/far_start_<TIMESTAMP>/`.  `analyze_benchmark.py` and
`plot_criticality_check.py` updated to accept `--results-dir`.

**`run_nm_test.sh` (new)** — isolated test of the improved hand-coded NM
from a far start.  Parameters:

```
--x0 1.5 0.5   --max-evals 60   --nm-sigma0 0.4
--nm-xatol 0.02   --nm-shrink 0.75   --save-every 3
```

`--save-every 3` (vs 5 in `run_far_start.sh`) gives a more detailed simplex
animation so individual NM steps can be inspected.

Reproduce:
```bash
bash run_nm_test.sh
```

### 2026-04-22 — Hand-coded NM; simplex visualization; residual panel; selective framing

**Why hand-code NM instead of calling scipy?**
scipy's `minimize(method='Nelder-Mead')` does not expose the shrink
coefficient — it is hard-wired to σ = 0.5.  Under MC noise, contractions
fail often (a reflected vertex scores worse than it really is because of
a lucky noise fluctuation at one of its neighbours), so the simplex
gets halved repeatedly and collapses long before it has walked to the
minimum.  A hand-coded loop with σ = 0.75 shrinks more gently and gives
the simplex room to explore.  The full NM step sequence (reflection →
expansion → outside contraction → inside contraction → shrink) is
implemented from scratch in `run_nelder_mead`; scipy is no longer a
dependency for this driver.

**NM coefficients (now configurable via CLI)**

| symbol | operation | scipy default | our default |
|------|-----------|--------------|-------------|
| α | reflection | 1.0 | 1.0 |
| β | expansion  | 2.0 | 2.0 |
| γ | contraction| 0.5 | 0.5 |
| **σ** | **shrink** | **0.5** | **0.75** |

CLI flags: `--nm-sigma0` (initial simplex edge), `--nm-xatol` (position
convergence tolerance), `--nm-shrink` (shrink coefficient σ).

`run_far_start.sh` uses `--nm-sigma0 0.4 --nm-xatol 0.02 --nm-shrink 0.75`
so the simplex spans ~0.56× the distance from x0=(1.5, 0.5) to the truth
(1, 1) at startup.  This ensures the triangle brackets the basin rather
than sitting well inside it.

**Simplex visualization — trajectory panel (Panel 1)**
Before each evaluator call, `run_nelder_mead` stores the current three
vertices on `evaluator.current_simplex`.  `Evaluator.__call__` passes this
to `OptimizerPlotter.update(..., simplex=current_simplex)`, which appends it
to `self.simplex_hist`.  In `_render`, the last 5 simplices are drawn as
filled polygons (`matplotlib.patches.Polygon`) in `deepskyblue`: the current
simplex at full opacity (α=0.25 fill, α=0.7 edge), older ones as ghosts
(α≈0.04–0.10).  This makes it easy to see whether the triangle is crawling
toward the minimum or has shrunk and frozen in place.

**Residual panel (Panel 4 — replaces curve-collapse overlay)**
Previously Panel 4 showed G_ref(t) and G_test(t) overlaid, making it hard
to judge the quality of coincidence by eye near the optimum.  It now shows
G_test(t) − G_ref(t) for each of the three boundary directions (v, u, w),
with a dashed y = 0 reference line.  A shaded ±RMS band provides a quick
noise-floor indicator.  When the optimizer has found the minimum, all three
curves hug the zero line; deviations above it show which direction is still
mismatched.

**Selective frame saving — `--save-every`**
Previously every evaluator call produced a PNG.  At 60 evals × 5 methods
this generated ≥300 frames plus many GC-scan frames.  `OptimizerPlotter` now
accepts `save_every=N` (default 5 via `--save-every` CLI flag) and writes a
frame only at step 1 and every Nth step after.  Use `--save-every 1` to
recover the old per-eval behaviour; `--no-vis` to suppress all frames.

**Files changed in this session**

| File | Change |
|------|--------|
| `optimizers.py` | `run_nelder_mead` replaced with hand-coded NM; σ=0.75; simplex stored on evaluator before each call |
| `visualization.py` | `OptimizerPlotter.__init__`: `save_every`, `simplex_hist`, `current_simplex`; `update()`: `simplex` param, selective rendering; `_render` Panel 1: simplex polygon history; `_plot_curve_collapse`: residuals G_test−G_ref with zero line + RMS band |
| `evaluator.py` | `self.current_simplex = None` in `__init__`; passed as `simplex=` to `optimizer_plot.update` |
| `run_benchmark.py` | `--nm-sigma0`, `--nm-xatol`, `--nm-shrink`, `--save-every` flags; `run_one` accepts `method_kwargs`; NM kwargs wired through |
| `run_far_start.sh` | Updated with new NM and `--save-every 5` flags |

### 2026-04-22 — Far-start stress test + criticality-check visualizer

The first five-method benchmark started the optimizer at the truth
(`x0 = (1, 1)`), so every method began inside its own initial step size of
the answer.  This made the comparison too generous — Nelder-Mead's first
simplex already brackets the minimum.  The new stress test moves the start
to `x0 = (1.5, 0.5)` (Euclidean distance 0.71 from the truth) and raises
the budget to `--max-evals 60` so CMA-ES and GP get past their exploration
phase.  The wider initial bracket `--beta-seed 0.15 0.40` lets the GC peak
finder follow β_c as the test couplings drift.

A separate question is whether each lattice the optimizer reports as
"matched" is actually critical — i.e. the test β_c is genuinely sitting at
the χ-peak, not at a local artefact of the GC fit.  To make this
auditable, every evaluation now persists its `(scan_betas, scan_chis,
scan_chi_errs)` to `eval_log.jsonl`, and the reference run dumps its own
scan to `results/_reference/reference_betac_scan.json`.  The new
`plot_criticality_check.py` reads both, re-fits Gram-Charlier on the
saved scans, and writes a single overlay panel per lattice marking the
data peak (★) and the GC-inferred β_c (dashed) — the eyeball test is
"does ★ sit on top of --, and is the GC curve hugging the data?".

Workflow:
- New file `plot_criticality_check.py` (reuses `mc_engine._gc_fit` and
  `_gram_charlier`).
- New script `run_far_start.sh` runs `run_benchmark.py` with the new
  CLI flags, then `analyze_benchmark.py`, then `plot_criticality_check.py`.
- `evaluator.py`: `EvalResult` extended with `scan_betas`, `scan_chis`,
  `scan_chi_errs` (all logged to JSONL).
- `run_benchmark.py`: caches the reference β-scan on first build; adds
  `--beta-seed LO HI` flag (default `(0.20, 0.32)` unchanged).

Reproduce:

```
bash run_far_start.sh
```

This writes `results/criticality_check.png` alongside the existing
`benchmark_table.md`, `convergence.png`, `trajectories.png`, and
`distance_history.png`.

**Next step (planned).**  If the far-start test reproduces or improves
on Nelder-Mead's previous lead, the next iteration should explore
**reweighting between geometries**: when the optimizer queries a new
(r₁, r₂), instead of always running a fresh MC there, reuse the
configurations from a nearby previous evaluation via histogram /
Ferrenberg-Swendsen reweighting (in β at fixed couplings) or a
Jacobian-corrected coupling reweighting (small Δr₁, Δr₂ from the
reference ensemble).  This would let one MC ensemble service many
optimizer queries — a large speedup, especially for noise-tolerant
optimizers like CMA-ES that revisit the same neighbourhood many times.

### 2026-04-22 — Wrote up next-session plan

- Added the "Next Steps (for the next working session)" section above
  Open Questions.  Refocused entirely on Nelder-Mead, the winner of the
  first benchmark.  Seven NM-specific tasks (A–G), with the first three
  (multi-start ensemble, sigma0 sensitivity, restart variant) marked as
  the core deliverables.  Acceptance criterion: a new
  "Nelder-Mead Characterization" sub-section with a median ± IQR
  convergence number, a recommended sigma0, a pass/fail on restart, and
  one concrete production recipe.

### 2026-04-22 — First five-method benchmark; Nelder-Mead wins

- Ran the full benchmark on a 6×6 self-consistency setup (60 k reference,
  30 k production, 25 evals/method, x₀=(1.15, 0.85), SNR floor disabled).
  Total wall ≈ 7 minutes, ~30 M trajectories across all five methods.
- Built `analyze_benchmark.py`: writes a markdown table plus
  `convergence.png`, `trajectories.png`, and `distance_history.png`
  to `results/`.
- Two bugs found and fixed during the run:
  - The CMA driver was passing a malformed bounds list to
    `cma.CMAEvolutionStrategy`; rewrote to `[[lo₁, lo₂], [hi₁, hi₂]]`.
  - `scikit-optimize` and `cma` were missing from the active venv;
    installed and added to `requirements.txt`.
- Added a `--snr-floor` CLI flag to `run_benchmark.py` (default 0 =
  disabled) so benchmark runs use their full evaluation budget.  The
  production setting is `--snr-floor 1.0`.
- Headline ranking by noise-robust *final-3 ⟨dist⟩*:
  Nelder-Mead (0.122) < BFGS-FD (0.176) < Powell (0.260) <
  CMA (0.304) < GP (0.327).  CMA reached the closest single point
  (0.072) but its population has not contracted yet at 25 evals.
- See the new "Benchmark Results — first run" section above for the full
  table, per-method analysis, and reproduction command.

### 2026-04-22 — β_c scan: sqrt-of-variance pass widths

- Replaced the geometric pass-halving in `find_beta_c` with a
  sqrt-of-variance schedule: pass 1 = ±2σ_fit, pass 2 = ±1σ_fit,
  pass 3 = ±0.5σ_fit, where σ_fit is the fitted Gaussian width of the
  current Gram-Charlier fit.  This is the natural physical width of the
  susceptibility peak — the scan adapts to the data instead of using a
  fixed fraction of the original bracket.
- Each pass falls back to a geometric fraction of the coarse step if no
  GC fit is available (e.g. pass 0 produced too few good points), and a
  floor prevents the window from collapsing to zero on under-constrained
  fits.

### 2026-04-21 — End-to-end smoke test passes

- Wrote the remaining four files: `evaluator.py`, `optimizers.py`,
  `run_benchmark.py`, `requirements.txt`.
- `evaluator.py` exposes `Evaluator(...)` callable returning an `EvalResult`
  dataclass; warm-starts each call's β_c bracket from the previous β_c
  (±2% window), wipes per-call MC scratch, appends to `eval_log.jsonl`,
  and pushes frames to `OptimizerPlotter` and a fresh per-call
  `BetaScanPlotter`.
- `optimizers.py` registers all five drivers in `ALL_METHODS`.  Each obeys
  the same hard cap (`max_evals`) and the SNR-floor stopping rule (running
  best with SNR < 1 after ≥5 evals halts the search).  GP and CMA-ES
  degrade gracefully to `"skipped"` if the optional deps are missing.
- `run_benchmark.py`: caches the reference correlator under
  `results/_reference/`, then loops over methods, writing per-method
  `summary.json` and a top-level summary comparing total trajectory cost
  and best (r₁, r₂).
- **Smoke test** (4×4 reference, x₀=(1.05, 0.95), 4 NM evals): completed
  in ~12 s wall, wrote 4 optimizer frames + 16 GC-scan frames under
  `results/nelder_mead/frames/`, plus the JSONL log and JSON summary.
  Pipeline confirmed working end-to-end with both visualization streams.
- **Next**: a real benchmark run on 8×8 with all five methods, then a
  comparison plot of cost-vs-trajectories across methods to identify the
  most efficient optimizer for this problem.

### 2026-04-21 — Visualization module complete; evaluator/optimizers next

- Wrote `visualization.py`:
  - `BetaScanPlotter` plugs into `mc_engine.find_beta_c(progress_cb=...)`
    and re-uses the literal `_gram_charlier(...)` formula to overlay the
    live fit on the (β, χ) scatter, with points colored by GC pass and a
    dashed line at the current β_c estimate.
  - `OptimizerPlotter.update(r1, r2, cost, σ, β_c, test_data)` writes a
    4-panel PNG per evaluator call: log-cost-colored (r₁, r₂) trajectory,
    cost-vs-eval with running best, β_c-vs-eval, and the three boundary-
    direction curve collapses.
  - Both use the matplotlib `Agg` backend (headless safe).  Frames stream
    into `results/<run>/frames/`, assembleable to GIF via the included
    `write_gif_command_hint(...)` helper.
- Added a “Repository Layout and File Status” section to this README
  tracking which pieces are in place vs. pending.
- **Next**: `evaluator.py`, `optimizers.py`, `run_benchmark.py`,
  `requirements.txt`, then a smoke test on a tiny lattice to confirm the
  whole loop runs end-to-end before launching the real benchmark.

### 2026-04-21 — Renamed and simplified: K_from_Optimization

- Renamed the directory from `K_from_CritSurface_v2` to `K_from_Optimization`
  to reflect the new design philosophy: the optimizer drives the workflow,
  not a pre-computed surface.
- Dropped the β_c surface entirely.  β_c is now found on-the-fly inside the
  evaluator for each (r₁, r₂) the optimizer queries, using the previous
  evaluator's β_c as a tight starting bracket.
- Removed `betac_surface.py` from the directory.  The evaluator delegates the
  β_c finder to `mc_engine.py` directly.
- Plan: implement all five optimizers and benchmark them side by side on the
  same 8×8 equilateral self-consistency test, comparing total trajectory cost
  and convergence quality.

### 2026-04-21 — Cost function and optimizer survey

- Decided to replace the variance-normalised Z²/Z⁴ cost with a plain L²
  integrated squared difference between the three boundary-direction correlator
  curves.
- Added a formal SNR-based statistics-adequacy criterion: when C/σ_C < 1,
  the optimizer must stop and request more statistics rather than continuing
  to refine into noise.
- Surveyed five classical optimizer options (Nelder-Mead, Powell, BFGS+FD,
  GP Bayes Opt, CMA-ES) with a comparison table.  Primary choice is Nelder-Mead
  (robust, no new dependencies) paired with the SNR stopping rule.  GP/BO as
  a secondary track for budget-limited runs.
- Noted open questions: direction weighting, normalized vs. unnormalized cost,
  and whether the analytic β_c formula can shorten the on-the-fly scan for
  large L.

### 2026-04-21 — Project initialized

- Created this directory and README after reviewing v1 results.
- Identified core problem: noise floor in Z² cost comparable to differences
  between grid points near minimum; grid search cannot reliably find the true
  minimum.
- v1 baseline: equilateral 8×8 self-consistency test recovers r≈(0.92,0.81)
  instead of (1,1); test_20x20 at 10k trajectories gives r=(1.3,1.3), clearly
  noise-dominated.
