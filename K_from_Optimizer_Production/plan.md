# QFE Cost-Function Development Plan

> **Project goal.** Build a *reliable, MC-only, theory-agnostic* cost
> function that maps lattice couplings (e.g. anisotropic plaquette
> weights `r₁, r₂, …` on a square test torus) to physical geometry by
> matching against MC data on a twisted reference torus. Once it works
> for 2D Ising it must port directly to φ⁴ (and other scalar QFTs)
> without bespoke modifications.

---

## 1. Problem Statement

We have:

* a **test lattice** — small, anisotropic, parametrized by couplings
  `θ = (r₁, r₂, …)`. (Currently a `(Lx,Ly,Tx,Ty) = (16,16,0,0)`
  square torus; about to upgrade to `(8,8,0,0)` at 50 k trajectories.)
* a **reference lattice** — large, twisted, with cycle-length triplet
  in known ratios. Currently `(39,48,−9,9)` whose three boundary
  cycles have lengths in ratio **4 : 5 : 6** (the “4-5-6 triangle”).
  Run at high statistics (~100 k trajectories).

Both lattices are simulated *at* their critical β. We measure the
connected two-point function `G_conn(m, n)` (and, if helpful,
susceptibility / higher-point functions). For a fixed candidate `θ`
on the test lattice we want a cost `J(θ)` whose minimum is the
coupling that makes the test lattice **physically equivalent** to the
reference.

**“Truth”** for the current 4-5-6 problem: `(r₁, r₂) ≈ (5.07, 7.74)`,
`β_c ≈ 0.0628`.

---

## 2. FSS Continuum-Limit Analysis (4-5-6 Geometry)

### Goal
Verify that the **twisted isotropic reference** (α=1,2,3,4; lattice scaled by
α from the `(13,16,−3,3)` base) and the **untwisted anisotropic test** (L=8,16,
24,32,48,64 at truth couplings r₁=5.065, r₂=7.743) both extrapolate to the
**same continuum-limit two-point function** G∞(t) along each of the three
boundary cycles of the 4-5-6 triangle.

### Method
Empirical power-law FSS:

$$
G(t;\, L_{\rm phys}) = G_\infty(t) + a(t) \cdot L_{\rm phys}^{-\omega}
$$

with G∞, a, ω all **free** (3-parameter fit per family per sample position t_k).

* **Test** physical scale: `L_phys = L` (cycles have |p| = L for square torus).
* **Ref** physical scale: `L_phys = α × |p_c^{α=1}|`, where the base cycle
  lengths are 14.731, 17.692, 11.790 for cycles 0, 1, 2.
* Fit uses the **N_FIT=4 largest available sizes** per family to stay in the
  asymptotic regime and have 1 dof.

### Datasets

| Family | Sizes available | Runs | Status |
|--------|----------------|------|--------|
| Test (untwisted) | L=8,16,24,32,48,64 | 50k–50k traj | α=1–3 done; α=4 running |
| Ref (twisted iso) | α=1,2,3,4 | 50k–100k traj | L=8–32 done; L=48,64 running |

Precompute scripts:
* **Original**: `precompute_456_fss.py` (α=1–3, L=8–32) → `results/_fss_456/`
* **Extended**: `precompute_456_fss_extend.py` (α=4, L=48, L=64) → same dir

### Output scripts
* `plot_fss_correlator_456.py` — 3-panel G(t) vs t overlay (FSS convergence check)
* `plot_fss_continuum_limit_456.py` — 3-row × 3-col continuum extrapolation figure:
  - Row 0: G vs 1/L_phys with power-law fit curves and G∞ intercept markers
  - Row 1: fitted ω(t_k) per family — shows whether corrections are ω≈2 or ω≈1
  - Row 2: G∞^test vs G∞^ref vs CFT envelope; % deviation (test−ref)/ref

### Current status
* Fixed ω=2 fits (earlier run): test G∞ sits ~15–22% **above** CFT, ref G∞
  ~13–30% **below** CFT — families disagree. Root cause uncertain (wrong ω or
  not yet asymptotic).
* Free-ω fits pending completion of extended MC (ref α=4, test L=48, L=64).
  All three jobs running in parallel via `ProcessPoolExecutor(max_workers=3)`.

---

## 3. What We Have Confirmed (Cost-Function Development)

| Finding | Evidence |
|---|---|
| Naïve L2 in log-G has a trivial bowl at the small-geometry corner. | All 14 huber/log/etc. kernels argmin at small r. |
| The bowl is caused by a **uniform multiplicative offset** `G_test ≈ c · G_ref` with `c ≈ 1.4`. | `viz_residuals_normalized.py` |
| Per-direction amplitude offsets `(c_u, c_v, c_w)` carry the anisotropy fingerprint: spread(c) is large near truth, small at iso/small-r. | `landscape_offset.py` |
| Dropping the contact-term sample (t≈0) and the wraparound sample (t≈1) is essential. | `offset_cost_drop_t0_t1.png` |
| Best cost so far: `shape/spread`, **d_truth ≈ 0.94 grid steps** at 5 k trajectories on a 16×16 test. | Same. |
| Per-site relative noise ≈ 5 % in the informative midrange t∈[0.25,0.75]; floor on a properly-averaged residual ≈ 0.7 %. | `landscape_snr.py` |
| `boundary_paths` and `_direction_lattice_steps` correctly construct u,v,w cycles and pick up the periodic-torus singularities at t=0 and t=1 on **both** test and reference. | `verify_cycles.py` |

**Constraint that drives the design.** Anything that bakes in
2D-Ising specifics (exact modular `τ`, theta functions, `Δ_σ = 1/8`,
analytic FSS amplitudes) won't generalize to φ⁴. The cost must use
**only** what MC produces.

---

## 4. Design Principles

1. **MC-only.** Inputs are `G_conn(m, n)` with errors; possibly
   `χ`, higher-point functions. No analytic 2D-Ising special structure.
2. **Reparameterization-invariant.** The matching task is
   `G_test(x) = c · G_ref(M·x)` for an unknown affine map `M` and
   unknown overall scale `c`. The cost must therefore be **scale-free
   and amplitude-free** so the trivial offset `c` does not create a
   spurious basin.
3. **Use cycle structure.** The torus has three fundamental boundary
   cycles `u, v, w` whose physical lengths form a triplet that *is*
   the geometric content of the lattice. Costs built from per-cycle
   observables capture anisotropy directly.
4. **Use error bars.** MC gives heteroscedastic uncertainties. A
   reliable cost must be a proper χ² (or use weighted fits) so that
   the resulting `θ̂` has a defensible covariance.
5. **Compose, don't stack.** Build small, orthogonal observables
   (decay-rate triplet, susceptibility triplet, amplitude triplet);
   combine them in a single weighted χ². Each piece must be
   theory-agnostic on its own.

---

## 5. Cost-Function Menu (MC-only)

### Tier 1 — Primary observables, scale-free by construction

#### A. **Cycle-length triplet from exponential decay**
For each cycle `i ∈ {u,v,w}` fit a single-state ansatz to **both** test
and ref correlators along that cycle:

$$
G(t \, \hat p_i) \;\approx\; A_i \left[\, e^{-E_i \cdot t \cdot |p_i|_\text{phys}} + e^{-E_i (1-t) |p_i|_\text{phys}} \,\right].
$$

Per cycle we recover the dimensionless decay rate
`λ_i = E_i · |p_i|_phys`. Take two scale-free ratios

$$
R_1^\lambda = \lambda_u / \lambda_w, \qquad R_2^\lambda = \lambda_v / \lambda_w
$$

and demand they match the reference's. Cost contribution:

$$
J_\text{decay} = \sum_{k=1,2} \left(\frac{R_k^{\lambda,\text{test}} - R_k^{\lambda,\text{ref}}}{\sigma_k}\right)^{2}.
$$

**Universality.** Only assumes lowest-state dominance along the
cycle. Holds for any massive scalar QFT near criticality.

#### B. **Susceptibility (zero-momentum) triplet**
Direction-projected zero-momentum susceptibilities

$$
\chi_i = \sum_{t \in \text{cycle } i} G_\text{conn}(t \, \hat p_i).
$$

Two scale-free ratios `χ_u/χ_w`, `χ_v/χ_w`. No fits, errors via
jackknife. Universal.

#### C. **Per-direction amplitude `c` triplet** (already partially built)
Fit a per-cycle multiplicative offset `c_i` such that `G_test(t·p̂_i) ≈
c_i · G_ref(t·p̂_i)`. This is what we currently call the offset
removal in `landscape_offset.py`. Promote `(c_u, c_v, c_w)` to two
scale-free ratios `c_u/c_w, c_v/c_w` and match.

**Three orthogonal probes.** A uses *decay rates*, B uses *integrated
amplitudes*, C uses *short-distance amplitudes*. All three must
agree at the truth, only there. Combine into a single 6-component
χ² for the final cost.

### Tier 2 — Robust statistical wrappers (theory-agnostic; drop-in)

#### D. **Bootstrap / jackknife the cost**
Resample MC trajectories → cost samples → cost mean and variance.
Either (a) minimize the lower confidence bound or (b) treat the cost
as a proper χ² so the recovered `θ̂` has a covariance matrix.

#### E. **GP smoothing of `G(m,n)`**
Fit a 2D Gaussian process to the test correlator with `conn_err` as
heteroscedastic observation noise. Evaluate the cost on the GP mean
rather than raw MC. Suppresses per-site speckle; sharpens any of
A–C.

#### F. **Heteroscedastic χ² in log-G**
Replace uniform-weight `(log G_test − log G_ref)²` with

$$
\frac{(\log G_\text{test} - \log G_\text{ref})^{2}}{\sigma^{2}_\text{test}/G_\text{test}^{2} + \sigma^{2}_\text{ref}/G_\text{ref}^{2}}.
$$

A no-brainer upgrade that doesn't change anything else.

### Tier 3 — Geometric, fully MC-driven

#### G. **Metric-tensor fit**
Posit `G_test(m,n) = F( √((m,n)ᵀ A (m,n)) )` for an unknown 2×2
positive-definite metric `A` and an unknown 1D function `F`. Fit `A`
(2 d.o.f.) jointly with `F` nonparametrically (e.g. cubic spline on
log-distance). Compare `A_test` vs `A_ref` after a common scale
convention. Most rigorous and the most general — works for any
scalar QFT whose two-point function is approximately a function of
one geodesic distance.

#### H. **Spectral axis ratio from 2D FFT**
2D FFT of `G_conn(m,n)` on the test torus → `C̃(k)`. Lowest-mass
momentum shell traces an ellipse whose axes encode the metric. Fit
the ellipse, compare ratio to ref's same ellipse.

#### I. **Higher-point function ratios** (for φ⁴ later)
Any scale- and field-rescaling-invariant n-point function ratio adds
constraints. Defer until Tier 1 stack works for Ising.

---

## 6. Recommended Build Order

| Step | Module | Output | Why now |
|---|---|---|---|
| 1 | `anisotropy_triplet.py` | computes `(λ_ratios, χ_ratios, c_ratios)` + jackknife errors on a single MC dataset | the three primary observables A, B, C in one place |
| 2 | `landscape_anisotropy.py` | replays each observable's heatmap on the existing precompute grid | first apples-to-apples comparison of A vs B vs C against the truth |
| 3 | `combined_chi2.py` | weighted χ² combining A+B+C with their jackknife covariance | the actual production cost |
| 4 | `bootstrap_wrap.py` | adds bootstrap error bars to *any* cost | tells us whether minima are real |
| 5 | wire `combined_chi2` into `cost.py` as a new `cost_mode='aniso_chi2'` | drop-in replacement for `huber_log` | production CMA-ES uses it |
| 6 | CMA-ES validation on real MC at 50 k from `x₀ = (1,1)` and `(5,6)` | recovered `(r̂₁, r̂₂)` with covariance | end-to-end test |
| 7 | metric-tensor fit (Tier 3 G) as a separate validator | `Â_test`, `Â_ref` independent recovery | catches errors in the cycle-based stack |
| 8 | port to φ⁴: same `combined_chi2`, same observables | recovers physical geometry from φ⁴ MC on the same lattices | proves portability |

---

## 7. Validation Protocol

For each new cost in this plan we report **all five** numbers:

1. **Argmin on the precompute grid** (current truth distance metric `d_truth`).
2. **Heatmap quality** (visual: smooth basin? speckle? trivial corner bowl?).
3. **CMA-ES from two seeds** (small `x₀=(1,1)` and near-truth `x₀=(5,6)`):
   does it converge to the same point?
4. **Bootstrap-implied uncertainty** on `θ̂`.
5. **Generalization check**: does the same cost run unchanged against
   reference `(13,16,−3,3)` and recover an answer consistent with
   the more-stats `(39,48,−9,9)` reference?

A cost only graduates to "production" after passing all five.

---

## 7. Current Status (snapshot)

* Test grid: `Lx16_Ly16_Tx0_Ty0`, 400 pkls, n_traj_prod = 5 000.
* Fine-grid zoom: `Lx16_Ly16_Tx0_Ty0_zoom25`, 625 pkls, step=0.25,
  r∈[2.5,8.5]², n_traj = 20 000. **Complete.**
* Reference: `_reference_Lx39_Ly48_Tx-9_Ty9`.
* Per-site SNR (midrange t): ≈ 20.
* Plan: launch 8×8 (or 16×16) test grid at 50 k trajectories
  while we build steps 1–4 above; rerun every kernel on the
  high-stats grid to see if the truth basin sharpens to the
  predicted ~3× level.

### 7.1 Cost-quality scoreboard

Truth at `(5.07, 7.74)`; `d_truth` is Euclidean distance from grid
argmin to truth in `(r₁, r₂)` space.

| Cost | Grid / stats | argmin | `d_truth` | Notes |
|---|---|---|---|---|
| huber_log dozen vs ref(13,16,−3,3) | coarse 5k | small-r corner | 8.56 | trivial bowl |
| huber_log dozen vs ref(39,48,−9,9) | coarse 5k | small-r corner | 5.84 | trivial bowl, less severe |
| half_ratio | coarse 5k | (1.5, 2.0) | 6.76 | diagonal ridge but wrong side |
| ellipse_residual | coarse 5k | – | 5.53 | half-decay contour ellipse fit |
| xcollapse_powerlaw | coarse 5k | – | 3.29 | cross-geometry curve collapse |
| **shape_only (drop t≈0,1)** | coarse 5k | (6.5, 8.5) | **1.62** | offset-removed L2 shape residual |
| **shape/spread (drop t≈0,1)** | coarse 5k L=16 | (4.5, 8.5) | **0.94** | best at coarse grid |
| shape/spread | coarse 5k L=8 | (varies) | 4.82 | FSS — much worse at small L |
| **shape/spread** | **fine 20k L=16** | **(4.50, 8.25)** | **0.76** | fine-grid zoom25; still speckled |
| A: half-decay ratio test↔ref | coarse 5k | (1.0, 10.0) | 4.65 | speckle dominated |
| B: susceptibility ratio test↔ref | coarse 5k | small-r corner | 8.56 | small-r limit gives trivial ratios |
| C: amplitude ratio test↔ref | coarse 5k | small-r corner | 8.56 | same — trivial ratios at small r |
| F: heteroscedastic χ²(log G) | coarse 5k | (2.5, 2.5) | 5.84 | mild improvement on plain log-L2 |
| Combined A + B + C | coarse 5k | small-r corner | 8.56 | inherits B/C corner |

### 7.2 Findings from `landscape_anisotropy.py`

* The "scale-free ratio" Tier-1 costs A, B, C **do not work as
  written**: at small `(r₁, r₂)` the test correlator decays so fast
  that all three per-cycle observables (decay rate, susceptibility,
  amplitude) become small and approximately equal across u,v,w.
  Their inter-direction ratios collapse toward unity — which happens
  to match the reference's ratios trivially in some normalisation,
  giving a false minimum at the small-r corner.
* The B and C heatmaps **do** show a clean diagonal valley running
  through the truth; the information is present but the cost
  normalisation lets the corner win. The required fix is to
  **regularise out the small-amplitude regime** (e.g. weight by
  signal magnitude, or take ratios *of differences from unity*
  rather than ratios themselves).
* The shape-residual cost (already known) at `d_truth ≈ 1.6` is
  still the best single observable; the `shape/spread` combination
  at `d_truth ≈ 0.94` is the best overall.

### 7.3 Implications for the plan

* The "compose three orthogonal scale-free ratios" recipe (Step 3
  Combined χ²) is *not* automatically reliable — the trivial-ratio
  failure mode at small r must be designed out before the combination
  is run.
* The empirical winner so far (`shape/spread`) is structurally an
  observable on the *anisotropy of c_i across directions*, i.e. a
  direct measurement of how non-uniform the per-cycle amplitude
  offset is — which is exactly the right physical signal but a
  different functional form from the "scale-free ratio" recipe.
* Next development priority shifts: take `shape/spread`-style
  observables seriously as a class, add bootstrap error bars to them,
  and re-derive the Tier-1 stack with the small-r failure mode
  excluded.

### 7.4 FSS behavior of shape/spread (L=8 vs L=16)

* `d_truth` improves dramatically with L: 4.82 (L=8, 50k traj) →
  0.94 (L=16, 5k traj) → 0.76 (L=16, 20k fine grid). FSS is the
  dominant effect — 5× larger than the statistical effect of going
  from 5k to 50k at fixed L.
* Visual: the basin in the shape/spread heatmap moves toward truth
  as L increases (`fss_compare.png`). The tilt is FSS-driven, not
  noise-driven.
* Even at 20k traj and fine grid (step=0.25), the heatmap is
  **highly speckled** — no smooth bowl. The shape/spread cost has
  high MC variance because its denominator (std of c_i across
  directions) is itself a noisy quantity. Finer grid reveals more
  noise rather than resolving the basin.
* Conclusion: shape/spread needs either (a) much larger L (≥ 32),
  or (b) a lower-variance replacement for the denominator.

### 7.5 Revised next steps

1. Add bootstrap / jackknife error bars to `shape/spread` so we can
   minimize a proper χ². (Step 4 of original plan, accelerated.)
2. Re-cast B and C as **deviation-from-isotropy** observables:
   instead of the raw triplet ratio `c_u/c_w`, use
   `Δ_c = (c_u − c̄)² + (c_v − c̄)² + (c_w − c̄)²` divided by
   `c̄²`. This vanishes whenever the test is locally isotropic
   (small-r AND truth-r both have characteristic anisotropy
   patterns; truth's matches ref's, small-r's doesn't).
3. Heatmap each new variant; only after the small-r corner is gone
   do we combine.
4. Run the 8 × 8 test grid at 50 k trajectories and re-evaluate
   every cost on it.
5. Investigate per-correlator FSS (§9.3) as the correct alternative
   to scalar-cost FSS (§9.10) — see §9.11.

---

## 8. Out-of-Scope (deliberately)

The following are tempting but excluded because they'd commit the
machinery to 2D-Ising specifics:

* exact modular-`τ` matching from the parallelogram parametrisation;
* theta-function fits to the all-to-all correlator;
* using `Δ_σ = 1/8` as a hardcoded conformal weight;
* SL(2,ℤ)-invariant cost on the full Z(τ) partition function.

These remain useful as **independent sanity checks** outside the
production cost path; we may use them once during validation but they
will not appear in `cost.py`.

---

## 9. Continuum-Extrapolated L2 Heatmap (FSS on both sides)

**Goal.** Replace the single-L `Σ (G_test − G_ref)²` with the same
cost evaluated on **continuum-extrapolated correlators on both sides**,
and produce one heatmap on the existing `(r₁, r₂)` grid. The cost
form is unchanged — the only new ingredient is FSS — so the heatmap
isolates "is the small-r corner bowl a finite-volume artifact?".

### 9.1 Geometry families

Scale linear size by integer factor `α` while keeping aspect and twist
fixed.

* **Test family** (square, no twist): `(8α, 8α, 0, 0)` for
  `α ∈ {1, 2, 3}` → L ∈ {8, 16, 24}. The existing 16×16 grid is α=2.
* **Reference family** (4-5-6 twist): `(13α, 16α, −3α, 3α)` for
  `α ∈ {1, 2, 3}` → identical aspect, identical twist fraction,
  identical cycle-length ratios. α=1 = `(13,16,−3,3)` and
  α=3 = `(39,48,−9,9)` are already cached; only α=2 = `(26,32,−6,6)`
  is new.

### 9.2 Sample lattice (fixed across α)

Pick once, before any MC, a finite set of physical sample positions
common to all `α`:

* For each cycle `i ∈ {u, v, w}`: `t_k = k/N_samp` for
  `k = 1, …, N_samp − 1` (drop t=0 and t=1; known-required from §7).
* `N_samp = 8` → `3 · (N_samp − 1) = 21` comparison samples regardless
  of L.

Each `(t_k, p̂_i)` maps to a physical Euclidean point on each torus;
both sides interpolate to that point with `_tile_interp`.

### 9.3 Per-side continuum extrapolation

For each candidate `θ = (r₁, r₂)` on the test side (and once for the
reference):

1. Run MC at every `α`, at that side's β_c(α; θ), with matched stats.
2. Build interpolators `G_α(x)` from the all-to-all data via the
   existing `_tile_interp`.
3. Evaluate at the fixed sample points `x_s` → `G_α(x_s)` with
   jackknife `σ_α(x_s)`.
4. Per-sample weighted least squares fit
   `G_α(x_s) = G_∞(x_s) + a(x_s)/L_α + b(x_s)/L_α²` (3 points; drop
   the `b/L²` term if it adds noise) → `Ĝ_∞(x_s)` with propagated
   `σ_∞(x_s)`.

Outputs: two length-21 vectors `Ĝ_∞^test(θ)` and `Ĝ_∞^ref` (cached).

### 9.4 The cost

$$
J(\theta) \;=\; \sum_{s=1}^{21} \bigl(\hat G_\infty^{\text{test}}(\theta; x_s) - \hat G_\infty^{\text{ref}}(x_s)\bigr)^{2}.
$$

Optional `chi2` variant divides each summand by
`σ_∞^test² + σ_∞^ref²`. Both reported on the same heatmap. No
reweighting, no offset removal, no per-cycle ratios — same form as
`l2_diff` from the dozen; FSS is the only new ingredient.

### 9.5 Implementation modules

| File | Role | Approx LOC |
|---|---|---|
| `precompute_landscape_multiL.py` | extends `precompute_landscape.py`: takes `--alphas 1 2 3`, runs the test scaling family per `(r₁,r₂)`, writes `grid/r1_X_r2_Y_a{α}.pkl` | 210 |
| `build_ref_multiL.py` | runs the three reference geometries once, writes `_reference_456_multiL/{a1,a2,a3}/two_point_all_to_all.dat` plus a manifest (α=1, α=3 already exist; only α=2 is new) | 120 |
| `continuum_extrap.py` | pure function: list of `(L, G_array, σ_array)` → `(G_∞, σ_∞)` per sample. Closed-form weighted least squares, no scipy | 80 |
| `landscape_l2_continuum.py` | replay-only: loads multi-L pkls + ref multi-L cache, builds the 21-vector on each, applies §9.3, computes `J(θ)`, writes `l2_continuum.png` and `.npz`. Also writes the single-L `l2_diff` heatmap from the same data for direct A/B comparison | 180 |

All four mirror the structure of `precompute_landscape.py` /
`landscape_dozen.py`.

### 9.6 Compute budget (rough)

Single 16×16 test point at 5k traj ≈ 1 unit.

* **Test family per `(r₁,r₂)`**: α=1 (8×8) ≈ 0.25, α=2 (16×16) ≈ 1,
  α=3 (24×24) ≈ 2.25 → ≈ 3.5 units per grid point.
* **Test grid** (20×20 = 400 pts): ≈ 1400 units, ~4× current grid
  cost. Overnight on 6 workers.
* **Reference**: only α=2 is new (≈ 4 units total). Negligible.

### 9.7 Validation panel (one figure, four subpanels)

1. `J_singleL(θ)` from α=2 only — reproduces today's "corner bowl".
2. `J_continuum(θ)` — the new cost.
3. `J_continuum_chi2(θ)` — σ-weighted version.
4. Per-θ mean of `|G_α=2 − Ĝ_∞|` over the 21 samples — diagnostic of
   how much FSS actually moved things.

Argmins, `d_truth`, `log10(min)` reported in the same scoreboard
format as `landscape_dozen.py`.

### 9.8 Order of operations

1. `build_ref_multiL.py` (only α=2 ref run is new).
2. `continuum_extrap.py` + 5-line unit test on synthetic
   `G_α = G_∞ + a/L`.
3. Pilot `precompute_landscape_multiL.py` on a 5×5 sub-grid —
   confirms wall time and disk.
4. Full 20×20 multi-L precompute overnight.
5. `landscape_l2_continuum.py` — produces the deciding figure.

Stop after step 5; only proceed to wire FSS into `cost.py` as a
production mode if the corner bowl is gone.

### 9.9 Decision rule

* `J_continuum` argmins **at or near (5.07, 7.74)** while `J_singleL`
  argmins at the corner on the *same* MC data → FSS is the entire fix
  and no cost-form changes are needed for production.
* `J_continuum` still argmins at the corner → the cost form is
  genuinely broken; §7.4 (Δ-from-isotropy etc.) is mandatory and FSS
  is at most a noise-reduction wrapper, not a fix.

### 9.10 First implementation: 1-1-1 (equilateral) ladder

Initial test problem before tackling 4-5-6. Truth: `(r₁, r₂) = (1, 1)`,
β_c = ln(3)/4 ≈ 0.27465 (exact). Picking the equilateral case lets
us check FSS behavior against an analytically known truth and an
isotropic correlator — anisotropy artifacts cannot mask FSS effects.

**Ladder.**

* Test sizes: `(Lx, Ly, Tx, Ty) ∈ {(8,8,0,0), (12,12,0,0), (16,16,0,0)}`
  (no twist, three points in `1/L`).
* Reference sizes (untwisted square at the same isotropic point, since
  for 1-1-1 the reference is a scaled copy of the test):
  `{(16,16,0,0), (24,24,0,0), (32,32,0,0)}`. "Effective size" = the
  linear extent that the test must FSS toward; here we pick reference
  L's roughly 2× the test L's.

The pairing is a "ladder" of `(L_test, L_ref)` cells:

| α | L_test | L_ref |
|---|---|---|
| 1 | 8     | 16    |
| 2 | 12    | 24    |
| 3 | 16    | 32    |

For each cell we compute the simple residual cost
`J_α(θ) = Σ_s (G_test_α(θ; x_s) − G_ref_α(x_s))²` on the fixed sample
lattice from §9.2. The **continuum cost** is the `1/L → 0`
extrapolation of `J_α(θ)` itself (linear or `1/L + 1/L²` fit on the
three rungs). This is the simplest possible FSS: extrapolate the
*scalar cost*, not each correlator sample, so there is no per-sample
fit complexity.

Heatmap: `J_∞(θ)` over the same `(r₁, r₂)` precompute grid, with
`J_α=1`, `J_α=2`, `J_α=3` shown as comparison panels.

**Why 1-1-1 first.**

* Truth at `(1,1)` is on a grid node and inside any reasonable bowl,
  so the small-r corner pathology cannot mimic the right answer.
* Ref and test families are both untwisted squares — minimal MC
  pipeline changes; reuses existing β_c finder seeds.
* If continuum extrapolation of `J(1,1) ≈ 0` while the off-truth
  cells stay finite, the ladder is working and we promote it to
  4-5-6 (twisted ref, anisotropic test).
* If `J_∞` argmins away from `(1,1)` on isotropic data, the residual
  cost is broken before FSS even enters — abort and revisit cost form.

**Modules to add (1-1-1 specific).**

| File | Role |
|---|---|
| `precompute_landscape_ladder.py` | computes the 9-cell ladder per `(r₁,r₂)` grid point; writes `grid/r1_X_r2_Y_α.pkl` × 3, plus reference cache for the three ref sizes |
| `landscape_l2_ladder.py` | replay-only: loads ladder, computes `J_α` per cell, fits `1/L → 0`, writes 4-panel heatmap (`J_α=1,2,3, J_∞`) and an `.npz` |

After the 1-1-1 heatmap looks clean, swap geometries to the 4-5-6
ladder (test families as in §9.1, ref families as in §9.1) and rerun
the same `landscape_l2_ladder.py` unchanged.

### 9.11 Scalar-cost FSS extrapolation: failure analysis

**Finding (from 1-1-1 ladder, 20k traj per rung).**
Extrapolating the *scalar cost* `J(L) = J_∞ + a/L` (§9.10) is
structurally broken for residual-style costs. Two independent
power-law errors:

1. **Wrong power at the truth.** At `θ = truth`, the test and ref
   correlators agree in the continuum limit, so `J_truth(L) → 0`
   as `L → ∞`. The leading FSS term is `J ~ c/L²` (volume
   suppression of the residual squared), not `c/L`. Fitting
   `J = J_∞ + a/L` forces `J_∞ < 0` at the truth:
     - 1-1-1 truth `(1,1)`:  `J_∞ ≈ −0.173`
     - (0.75, 0.75) near-isotropic:  `J_∞ ≈ −0.002`
   The extrapolated surface cannot distinguish truth from
   near-isotropic points; its minimum is dictated by extrapolation
   bias, not physics.

2. **Non-monotone J(L) off-truth.** Away from truth the residual
   `ΔG = G_test − G_ref` has a finite continuum limit `ΔG_∞ ≠ 0`
   and a subleading FSS correction `Δa/L`. The scalar cost is
   `J(L) = ΔG_∞² + 2 ΔG_∞ Δa/L + (Δa)²/L² + …`.
   The cross term `∝ ΔG_∞ · Δa / L` changes sign when `ΔG_∞` and
   `Δa` have opposite signs, making `J(L)` non-monotone. Measured
   example: `θ = (1.5, 1.5)` shows J rising L=8→12 before falling
   L=12→16. A one-parameter fit on a non-monotone curve returns
   garbage regardless of statistics.

**Conclusion.** Scalar-cost FSS cannot be salvaged by more statistics
or more rungs. Do **not** use §9.10 as a production wrapper.

**Correct alternative — per-correlator FSS (§9.3).** Extrapolate
each `G(x_s)` sample independently before squaring. This avoids
both pathologies: the linear fit `G_α(x_s) = G_∞(x_s) + a(x_s)/L`
is correctly specified at every grid point (with zero residual at
truth), and cross terms never appear because the fit precedes
squaring.

**Single-rung alternative.** For the 4-5-6 test at L=16, using
L=16 directly with a matched-scale reference is a viable production
path that bypasses FSS entirely. Current `d_truth` at L=16/20k
fine grid: 0.76. Estimated at L=32: ~0.3–0.5.

| Option | d_truth (est.) | Compute cost | When to use |
|---|---|---|---|
| L=16 single rung, shape/spread | 0.76 | 1× | optimizer seed |
| L=32 single rung, shape/spread | ~0.3–0.5 | 4× | tighter seed |
| Per-correlator FSS §9.3 | ~0.0 | 3× + complexity | production quality |
| Scalar-cost FSS §9.10 | broken | — | **do not use** |

### 9.12 FSS two-point correlator figures

Demonstrates that FSS is working for both geometries by showing
`G_conn(t)` along each boundary cycle converging to the CFT prediction
`G(t) = A/sin(πt)^{1/4}` (2D Ising Δ_σ = 1/8) as lattice size grows.

#### 9.12a  1-1-1 equilateral geometry
- **Script**: `plot_fss_correlator_111.py`
- **Output**: `results/_ladder_111_line20k/fss_correlator_111.png`
- **Data**: test L=8,16 from `_ladder_111_line20k/test/grid/` (20k traj PKL);
  ref L=24,32 from `_ladder_111/ref/` (20k traj all-to-all)
- **Method**: L ∈ {8,16,24,32} (all multiples of N_SAMP=8); sample positions
  t_k = k/8 hit exact integer lattice sites m = k·(L/8) — no interpolation.
  CFT amplitude fitted to L=32. m-index annotations on each marker.
- **Finding**: curves converge monotonically to CFT from above; shape clearly
  1/sin^{1/4}; cost sample positions well-placed along the curve.

#### 9.12b  4-5-6 anisotropic geometry
- **Script**: `plot_fss_correlator_456.py`
- **Output**: `results/_fss_456/fss_correlator_456.png`
- **Data** (generated by `precompute_456_fss.py`):
  - REF twisted iso (k1=k2=k3=1): α=1 (13×16, copied), α=2 (26×32, new 50k,
    β_c=0.26589), α=3 (39×48, copied) → `results/_fss_456/ref/a{α}/`
  - TEST untwisted at exact truth (r1≈5.0652, r2≈7.7429, k3=1):
    L=8 (100k, β_c=0.05507), L=16 (100k, β_c=0.05896),
    L=24 (50k, β_c=0.06030), L=32 (50k, β_c=0.06082)
    → `results/_fss_456/test/L{L}/`
- **Method**: 3-panel figure, one per boundary cycle. Both ref (blue, α=1,2,3)
  and test (orange, L=8,16,24,32) sampled at t_k=k/8 via `_tile_interp`.
  CFT amplitude fitted per cycle to α=3 ref. Cycle 0,1,2 amplitudes:
  A≈0.208, 0.196, 0.238 at α=3.
- **Finding**: both families converge to the same CFT curve in all three cycles,
  confirming that the 4-5-6 geometry matching and truth couplings are correct.
  The test (untwisted anisotropic) and ref (twisted isotropic) approach the
  identical continuum limit, validating the core methodology.

---

## 10. Production Matching-Cost Landscape (next milestone)

### Context

Landscape replay on the existing diagnostic precompute (5k traj/cell, L=8 & 16
only) using the §4 matching cost shows
**Pearson(cost, d_truth) ≈ 0** — the cost field is MC noise at this statistics
level. The matching-cost design is correct (validated to 0.28% in
`match_456_twisted_vs_untwisted.py`); what is missing is a dedicated
higher-stats precompute. See `algo.md §11` and `LANDSCAPE.md §8` for the
full diagnostic.

### Goal

Build a production landscape precompute on which the matching-cost basin at
the analytic truth `(r₁, r₂) = (5.07, 7.74)` is **clearly resolved**
(Pearson ≲ −0.5), then confirm it with the replay script before seeding CMA-ES.

### Step 1 — Production precompute

Run `precompute_landscape.py` (or a new wrapper) with the following knobs:

| knob              | diagnostic (current)  | production (target)               |
|-------------------|-----------------------|-----------------------------------|
| `n_traj_prod`     | 5 000                 | **≥ 50 000**                      |
| FSS sizes `L`     | 8, 16                 | **8, 16, 24, 32** (4 sizes)       |
| grid (r₁, r₂)    | step 0.5, r ∈ [0.5,10] | step 0.5 near truth; step 1.0 periphery |
| observables       | `G_conn(m,n)`         | same (all-to-all)                 |

Estimated cost: 400 cells × 4 sizes × 50k traj ≈ **80M sweeps** (~6–8 h on 8
workers). Alternatively run a coarse pass first (step 1.0, 20k traj) to
locate the basin, then refine at 0.5/50k near truth.

Script entry point: `precompute_landscape.py --n-traj 50000 --sizes 8 16 24 32`
(add `--sizes` CLI flag if not present).

### Step 2 — Matching-cost replay gate

After the precompute finishes, run:

```bash
python landscape_matching_456.py
```

Update `TEST_SIZES = [8, 16, 24, 32]` in the script header.

**Acceptance criterion**: Pearson(cost, d_truth) ≲ −0.5 for at least two
individual sizes. If the criterion is met:

- argmin of the continuum-extrapolated landscape should lie within
  1.5 r-units of truth (d < 1.5).
- The FSS extrapolation must use a **multi-L jackknife fit** (upgrade from
  the current 2-point `c_∞ = 2c16 − c8`).

If not met, double `n_traj_prod` and rerun.

### Step 3 — Jackknife multi-L FSS extrapolation

Replace the 2-point estimator with a proper linear fit:

$$
c(L) = c_\infty + a / L, \quad L \in \{8, 16, 24, 32\}
$$

Fit by ordinary least squares; use jackknife-over-L for error bars on
$c_\infty$. Clipping negative values is **not** done on the fitted intercept
(negative intercept = consistent with truth basin at that cell).

### Step 4 — CMA-ES seeding

Once the landscape basin is confirmed:

1. Use the continuum-extrapolated argmin `(r₁*, r₂*)` as CMA-ES initial
   point `--x0 r₁* r₂*`.
2. Set `--cma-sigma0` to half the basin half-width (estimated from the
   landscape heatmap).
3. Run CMA-ES with the §4 matching cost (not the old `test_native`):
   ```bash
   python run.py --config cfg_456.json \
       --x0 <r1*> <r2*> --cma-sigma0 0.5 \
       --max-evals 150 --multidonor-2pass \
       --n-workers 4 --save-frames \
       --run-name 4_5_6_prod_seed
   ```
4. Validate convergence with `match_456_twisted_vs_untwisted.py` at the
   recovered `(r₁, r₂)`.

### Step 5 — Acceptance criteria for the full pipeline

| check | pass condition |
|-------|----------------|
| Landscape Pearson | ≲ −0.5 (at least 2 individual L) |
| Landscape argmin  | d_truth < 1.5 r-units |
| CMA-ES best eval  | d_truth < 0.5 r-units after 150 evals |
| Matching validation | A_test / A_ref = 1 ± 2% (per side) |

### Open items

- [ ] Add `--sizes` CLI flag to `precompute_landscape.py`
- [ ] Implement jackknife multi-L FSS fit in `landscape_matching_456.py`
- [ ] Run production precompute (Step 1)
- [ ] Re-run replay gate with 4 sizes (Step 2–3)
- [ ] Seed CMA-ES run (Step 4)
- [ ] Document outcome in `algo.md §12` and `LANDSCAPE.md §9`
