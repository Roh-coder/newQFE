# K_from_Optimizer — Production

**Plug-and-play CMA-ES coupling optimizer with Ferrenberg–Swendsen
reweighting, automatic 3-pass fallback, optional process-pool
parallelism, and a live PNG viewer.**

This directory is **fully self-contained** — copy it to any machine,
build the C++ simulator, and run `python run.py`.  No sister-directory
dependencies.

---

## Contents

```
run.py                      ← single CLI entry point + the editable CONFIG dict
prod_runtime.py             ← reweight→3-pass fallback patch + parallel pool
live_viewer.py              ← Tk window that auto-reloads the optimizer PNG
run_overnight.sh            ← 4-run validation harness (legacy / RW / RW+FB / parallel)
check_overnight.py          ← post-run report generator (writes OVERNIGHT_REPORT.md)
Makefile                    ← builds bin/ising_tri_twisted_parallelogram
src/, include/              ← C++ simulator source
evaluator.py                ← per-eval (r1, r2) → (cost, σ_cost) [from upstream]
mc_engine.py                ← run_simulator, find_beta_c, find_beta_c_reweight,
                               find_beta_c_multidonor, find_beta_c_multidonor_2pass
optimizer.py                ← run_cmaes (hand-coded), run_nelder_mead, run_bo
parallel.py                 ← legacy GenerationPool (kept for reference)
reweight.py                 ← Reweighter, CombinedReweighter
visualization.py            ← OptimizerPlotter (4-panel PNG)
dashboard.py                ← rich-based terminal dashboard
cost.py                     ← L⁴ power-mean curve-collapse cost
betac_cache.py              ← optional persistent β_c interpolation cache
test_multidonor_vs_3pass.py ← side-by-side: multi-donor 1-pass vs 3-pass GC scan
test_multidonor_2pass.py    ← three-way: 1-pass MD vs 2-pass MD vs 3-pass GC scan
fss_betac_scan.py           ← isotropic FSS validation of β_c finder (ν fit)
fss_aniso_betac.py          ← per-anisotropy FSS validation of β_c finder
stress_aniso_betac.py       ← β_c finder stress test across anisotropic couplings
stress_twisted.py           ← end-to-end workflow stress test on twisted tori
stress_workflow.py          ← end-to-end workflow stress test (untwisted, vs L)
strategy_test.py            ← compares 4 CMA-ES strategies on the 4-5-6 problem
ref_size_test.py            ← grid + CMA scan over multiple reference lattice sizes
paraboloid_interp_test.py   ← analytic conceptual test of LinearND interp residuals
cft_torus.py                ← analytic Ising CFT torus correlator (Jacobi theta)
cost_zoo.py                 ← five candidate cost variants (relative, log,
                               arclength, effmass, pmean4) — see §8.6
modular_data.py             ← three-direction τ-aligned profile builder + costs
                               + noise injection (see §8.7)
modular_testbed.py          ← noise-free / noise-controlled τ-grid CMA-ES
                               testbed driver (see §8.7)
precompute_456_fss.py       ← runs MC for 4-5-6 FSS families: ref α=1,2,3 and
                               test L=8,16,24,32 → results/_fss_456/
precompute_456_fss_extend.py← extends the FSS dataset: ref α=4, test L=48,64
                               (adds dof for free-ω continuum extrapolation)
plot_fss_correlator_456.py  ← 3-panel FSS correlator convergence plot for 4-5-6
plot_fss_correlator_111.py  ← same for 1-1-1 (equilateral) ladder geometry
plot_fss_continuum_limit_456.py
                            ← empirical free-ω FSS fit G=G∞+a·L^(-ω), separate
                               fits for twisted-ref and untwisted-test families;
                               3-row figure (FSS bands, ω values, G∞ comparison)
                               → results/_fss_456/fss_continuum_limit_456.png
```

---

## 1. Quick start

```bash
# 1. Build the C++ simulator (one-time; auto-runs on first python run.py)
make

# 2. Install Python deps in a venv
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 3. Default 12×12 run, ~5–10 min on 4 cores
python run.py --n-workers 4

# 4. Same thing, but with a live PNG window that updates each eval
python run.py --n-workers 4 --live --save-frames

# 5. Print the resolved CONFIG (after CLI overrides) and exit
python run.py --print-config
```

Outputs land in `results/<run_name>/`:

```
results/<run>/
  summary.json            final best (r1, r2), β_c, cost, fallback stats
  eval_log.jsonl          one JSON object per evaluation
  fallback_log.jsonl      one record per β_c finder call (status + geom)
  run_meta.json           run name, geometries, full CONFIG snapshot
  frames/
    optimizer_cmaes.png   rolling 4-panel frame (overwritten each eval)
  frames_history/         only when --save-frames
    frame_0001.png … frame_NNNN.png
```

---

## 2. CONFIG vs CLI

`run.py` ships a single `CONFIG` dict at the top of the file — edit it
in-place for permanent defaults.  Every key is also exposed as a CLI
flag (dashed names: `cma_popsize` → `--cma-popsize`):

| Flag | CONFIG key | Meaning |
|---|---|---|
| `--config FILE.json` | — | Merge JSON keys on top of CONFIG |
| `--results-dir DIR` | `results_dir` | Where outputs live |
| `--run-name NAME` | `run_name` | Sub-directory name (auto if omitted) |
| `--max-evals N` | `max_evals` | Total evaluation budget |
| `--n-workers N` | `n_workers` | 1=serial, ≥2=ProcessPool |
| `--cma-popsize N` | `cma_popsize` | λ (evals/generation) |
| `--cma-sigma0 X` | `cma_sigma0` | Initial CMA-ES step size |
| `--n-traj-prod N` | `n_traj_prod` | MC trajectories per eval |
| `--reweight-n-traj N` | `reweight_n_traj` | Parent ensemble for reweight |
| `--save-every N` | `save_every` | Render PNG every N evals |
| `--ref-Lx`, `--ref-Ly`, `--test-Lx`, `--test-Ly` | … | Lattice sizes |
| `--x0 R1 R2` | `x0` | Optimizer starting point |
| `--r-lower R1 R2` | `r_lower` | Lower bound on `(r1, r2)`; out-of-bounds CMA-ES offspring receive a death penalty (see §8) |
| `--r-upper R1 R2` | `r_upper` | Upper bound on `(r1, r2)`; same semantics |
| `--no-reweight` | `reweight=False` | Force the legacy 3-pass scan |
| `--no-vis` | `no_vis=True` | Skip PNG rendering (fastest) |
| `--no-dashboard` | `dashboard=False` | Skip the rich terminal dashboard |
| `--live` | `live=True` | Pop a Tk window watching the rolling PNG |
| `--save-frames` | `save_frames=True` | Archive every PNG to `frames_history/` |
| `--live-poll-ms MS` | `live_poll_ms` | Viewer refresh interval |

Tuple keys not on the CLI (e.g. `beta_seed`, `ref_beta_seed`) are best
edited in `CONFIG` directly or supplied via `--config patch.json`.

---

## 3. The three speed-ups, and how they compose

The per-evaluation cost has two terms: the **β_c finder** and the
**production MC**.  We attack both:

| Speed-up | What it does | Effect | How to enable |
|---|---|---|---|
| **S4 — FS reweighting** | Replaces the 3-pass GC scan with one long parent MC + reweighted χ(β) on a dense grid | β_c finder time drops from ~3× to ~1× the production MC | `reweight: True` (default) |
| **S4a — multi-donor 1-pass** | Tiles `n_donors` parent ensembles across the bracket; `CombinedReweighter` merges them so every β has sufficient N_eff | Ends cold-start / narrow-bracket failures without triggering the 3-pass fallback | `--multidonor` |
| **S4b — multi-donor 2-pass** | Pass 1: coarse tile → GC fit locates peak interval, with automatic bracket translation if the peak lands on an edge. Pass 2: dense donors over ±1 gap around that interval. All donors combined for final fit. | Concentrates MC budget near the peak; ~1.5–2× faster than 3-pass at comparable accuracy; robust to far-from-equilateral starts where the seed bracket misses β_c | `--multidonor-2pass` |
| **3-pass fallback**     | Catches reweight failures transparently and re-runs the legacy GC scan | Makes S4/S4a/S4b safe to leave on for far-from-equilateral starts | always on, audited via `fallback_log.jsonl` |
| **S3 — parallel pop**   | Fans the λ candidates of each CMA-ES generation across N worker processes | ~min(λ, N)× wall-time reduction per generation | `--n-workers N` (≥2) |

### Audit the fallback rate

```bash
# Quick breakdown after a run
python -c "
import collections, json, sys
c = collections.Counter()
for line in open(sys.argv[1]):
    c[json.loads(line)['status']] += 1
print(c)
" results/<run>/fallback_log.jsonl
```

Healthy near-equilateral runs show ~100% `reweight`, ~0% `fallback_3pass`.
Far-start runs (e.g. `--x0 1.5 0.5`) typically show 5–20% fallback —
that's the safety net doing its job.

> **Parallel caveat:** `summary.json["fallback_stats"]` only counts
> events in the **main** process.  Worker processes maintain their own
> counters; for an aggregate count across all workers, parse
> `fallback_log.jsonl` (it is `flock`-coordinated and contains one row
> per call from every worker).

---

## 4. β_c finder hierarchy

`prod_runtime.py` patches `mc_engine.find_beta_c_reweight` with a
three-tier cascade.  At each call the fastest method that succeeds is
used; failures fall through to the next:

```
┌─────────────────────────────────────────────────────────────────┐
│  Tier 1 — multi-donor reweighting  (find_beta_c_multidonor)     │
│    • N donors auto-placed by action-variance heuristic          │
│    • CombinedReweighter guarantees N_eff coverage everywhere    │
│    • Validated: ≤0.21σ bias vs 3-pass on 12×12, 1.3× faster    │
│    enabled via --multidonor                                      │
└──────────────────────┬──────────────────────────────────────────┘
                       │ N_eff floor not met anywhere
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│  Tier 2 — single-donor FS reweighting  (find_beta_c_reweight)  │
│    • One parent MC at bracket midpoint                          │
│    • Fast when β_c is close to the starting point              │
└──────────────────────┬──────────────────────────────────────────┘
                       │ peak outside reweight validity window
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│  Tier 3 — 3-pass GC scan  (find_beta_c)                        │
│    • Trusted ground truth; most expensive                       │
│    • Always succeeds within the bracket                         │
└─────────────────────────────────────────────────────────────────┘
```

### Two-pass multi-donor (S4b)

`find_beta_c_multidonor_2pass` implements an interval-based refinement
strategy on top of the standard multi-donor tile:

1. **Pass 1** — `pass1_n_donors` donors evenly across `[beta_lo, beta_hi]`
   with `pass1_budget_frac` of the total trajectory budget.  A GC fit
   on the combined χ(β) identifies which gap `[p1_betas[i], p1_betas[i+1]]`
   contains the approximate peak β*.

2. **Pass 2** — dense donors over the interval-based window
   `[p1_betas[max(0, i−1)], p1_betas[min(N−1, i+2)]]` — the peak gap
   **plus one neighbouring gap on each side** — with the remaining budget.
   Donor spacing is controlled by the action-variance heuristic
   `|Δβ|_safe = sqrt(−ln(n_eff_floor) / Var(E))`.

3. **Combined fit** — all Pass-1 and Pass-2 donors are merged in a single
   `CombinedReweighter`.  The final GC fit is restricted to the Pass-2
   window for sharper peak resolution; Pass-1 wings prevent N_eff dropout
   in the tails.

Benchmark (12×12 equilateral lattice, 30 000 traj/method, single run):

| method       | β_c     | Δ vs 3-pass | wall |
|--------------|---------|-------------|------|
| 3-pass GC    | 0.2537  | —           | 5.5s |
| 1-pass MD    | 0.2492  | 3.4 σ       | 4.5s |
| **2-pass MD**| **0.2514** | **1.7 σ** | **3.3s** |

The 2-pass approach halves the bias relative to 1-pass while remaining
1.7× faster than the 3-pass scan.  The larger jackknife σ on β_c is a
budget-split artefact (see next-steps item 2 below).

### Bracket translation (off-isotropic robustness)

Both `find_beta_c_reweight` and `find_beta_c_multidonor_2pass` now
detect when the χ(β) maximum (or the GC-fit peak) lands within the
first/last 10 % of the current scan grid.  In that case the bracket is
translated by **half its span** in the direction of the peak (left if
the peak is at the low edge, right if at the high edge) and the scan
resumes:

* In `find_beta_c_reweight` the parent ensemble is re-spawned at the
  new bracket midpoint; the recenter loop continues up to
  `max_recenters` translations.
* In `find_beta_c_multidonor_2pass` Pass-1 donors are *accumulated*
  across translations (up to `MAX_TRANSLATIONS = 6` half-span hops),
  so wing coverage strictly grows.  Only when the combined χ(β)
  argmax lands in the interior does the routine proceed to Pass 2.

This makes the optimizer robust on strongly anisotropic geometries
(e.g. the 4-5-6 triangle) where the true β_c can sit far outside the
initial seed bracket.

### CLI flags for the 2-pass path

| Flag | CONFIG key | Default | Meaning |
|---|---|---|---|
| `--multidonor-2pass` | `multidonor_2pass` | False | Enable the 2-pass path (slots in *ahead* of 1-pass MD) |
| `--multidonor-2pass-pass1-n-donors` | `multidonor_2pass_pass1_n_donors` | 4 | Donors used in Pass 1 |
| `--multidonor-2pass-pass1-budget-frac` | `multidonor_2pass_pass1_budget_frac` | 0.30 | Trajectory share for Pass 1 |
| `--multidonor-2pass-alpha` | `multidonor_2pass_alpha` | 0.75 | Pass-2 donor spacing in units of `|Δβ|_safe` |
| `--multidonor-2pass-pilot-n-traj` | `multidonor_2pass_pilot_n_traj` | 2000 | Pilot traj. count when `Var(E)` unknown |

When `--multidonor-2pass` is set the runtime cascade in
`prod_runtime.py` becomes

    2-pass MD  →  1-pass MD (if --multidonor)  →  single-donor reweight  →  3-pass GC scan

failing through transparently and recording each step in
`fallback_log.jsonl`.

### Planned augmentations (next steps)

1. **Adaptive `pass1_n_donors`** — set gap width = δ_safe from a cheap
   pilot run, giving guaranteed single-gap resolution without user tuning.
2. **Budget recycling** — after the Pass-2 window is chosen, reallocate
   the Pass-1 budget from outside-window donors to Pass-2 donors, keeping
   the total trajectory count fixed and reducing the jackknife uncertainty.
3. **Three-pass extension** — optional Pass 3 placing a handful of donors
   within a single δ_safe of the Pass-2 β_c for sub-percent uncertainty
   at large lattice sizes.
4. **Unified test script** — merge `test_multidonor_2pass.py` into
   `test_multidonor_vs_3pass.py` as a third curve on the same shared-budget
   comparison plot.

---

## 5. The visualizer

The 4-panel PNG (`frames/optimizer_cmaes.png`) shows:

1. **(r1, r2) trajectory** with CMA-ES generation ellipses
2. **cost vs eval**
3. **β_c vs eval**
4. **per-direction residuals** of the test correlator vs the reference

### Three modes

| Mode | Flag(s) | Behaviour |
|---|---|---|
| Headless rolling | (default) | Single PNG `optimizer_<method>.png` overwritten every `save_every` evals |
| Live window | `--live` | Stdlib Tk window auto-reloads the rolling PNG every `--live-poll-ms` (default 500 ms).  Window stays open after the run finishes — close it to exit |
| Presentation archive | `--save-frames` | Each new rolling-PNG mtime is also copied to `frames_history/frame_NNNN.png`.  Combine with `ffmpeg` for an animation: `ffmpeg -framerate 4 -pattern_type glob -i 'frames_history/frame_*.png' -c:v libx264 -pix_fmt yuv420p run.mp4` |

`--live` requires Tk (preinstalled with Python on Windows/macOS; on Linux
install `python3-tk`).  Pillow is optional but improves PNG rendering
quality.

---

## 6. Overnight validation suite

```bash
bash run_overnight.sh                    # → results/_overnight/
bash run_overnight.sh /scratch/myroot    # custom root
```

Runs four 30-eval CMA-ES runs back-to-back on a 12×12 lattice:

| Run | Description | Expectation |
|---|---|---|
| **A_legacy** | 3-pass GC scan on every eval (`--no-reweight`) | baseline timing & β_c |
| **B_reweight** | Reweight on, fallback never fires | matches A on β_c, ~2× faster |
| **C_reweight_fallback** | Squeezed reweight bracket forces fallbacks | mix of `reweight` + `fallback_3pass`, β_c still matches A |
| **D_parallel** | Reweight + 4 workers | additional ~2–3× wall reduction |

`check_overnight.py` writes `OVERNIGHT_REPORT.md` with a comparison
table (n_eval, wall, s/eval, best (r1,r2), best β_c, speedup, fallback
counts).  Total wall ~30–90 min on a typical 4-core laptop.

---

## 7. Switching geometries / starting points

Edit the geometry block of `CONFIG`, or pass on the CLI:

```bash
# 24×24 production run from a far-from-equilateral start
python run.py \
    --ref-Lx 24 --ref-Ly 24 --test-Lx 24 --test-Ly 24 \
    --max-evals 60 --cma-popsize 8 --n-workers 4 \
    --n-traj-prod 50000 --reweight-n-traj 80000 \
    --x0 1.5 0.5 --save-frames
```

The reference correlator is cached per-geometry under
`results/_reference_Lx<..>_Ly<..>_Tx<..>_Ty<..>/`.  Switching back to a
previously-built geometry skips the ~1-minute reference build entirely.

### Cross-geometry runs (twisted ref vs untwisted test)

The ref and test lattices may have *different twists*.  Pass twist
parameters via a small JSON config:

```bash
cat > /tmp/cfg_456.json <<'EOF'
{"ref_Tx": -3, "ref_Ty": 3, "test_Tx": 0, "test_Ty": 0}
EOF

python run.py \
    --config /tmp/cfg_456.json \
    --ref-Lx 13 --ref-Ly 16 --test-Lx 13 --test-Ly 16 \
    --x0 3.0 4.0 --cma-sigma0 3.0 \
    --max-evals 150 --multidonor-2pass \
    --n-workers 2 --save-frames \
    --run-name 4_5_6_triangle --no-dashboard
```

This is the canonical setup for matching an *untwisted* test torus to
an asymmetrically twisted reference.  The (Lx,Ly,Tx,Ty) = (13,16,−3,3)
geometry above gives reference cycle lengths in ratio 4 : 5 : 6 to
within 0.04 %; the analytic optimum from the sinh rule is

* r1 = (2 ln 5 − ln 7)/(2 ln 3 − ln 7) ≈ 5.0652
* r2 = ln 7/(2 ln 3 − ln 7) ≈ 7.7429
* β_c = ln(3/√7)/2 ≈ 0.0628

Note that β_c here is far below the default seed bracket [0.20, 0.32].
The bracket-translation logic in the β_c finders (see §4) sweeps until
the χ(β) peak sits in the interior, so the seed bracket need not
contain β_c — but a very wide CMA-ES initial Gaussian (σ₀ ≳ 3 in
the r-space units of this example) is required to actually reach the
far-from-isotropic optimum.

---

## 8. Cross-geometry pitfalls — the 4-5-6 triangle case study

When the *test* and *reference* lattices have different shapes — the
canonical example being a 16×16 untwisted equilateral test compared
against a (13, 16, −3, 3) twisted reference whose three boundary
cycles are in 4 : 5 : 6 ratio — naive CMA-ES runs from a generic
starting point (e.g. `--x0 3 3 --cma-sigma0 0.5`) drift toward a
spurious low-cost region near `(r1 ≈ 1, r2 ≈ 0, β_c ≈ 0.42)` and
never reach the analytic optimum at `(r1, r2, β_c) ≈ (5.07, 7.74, 0.063)`.

This section documents the diagnostics performed on run `j002317_p40196`
and the resulting fixes.  The same mechanisms apply to *any*
cross-geometry run where the three reference boundary cycles differ
significantly in length.

### 8.1 Three independent failure modes

1. **Initial CMA-ES Gaussian too narrow.**  With `cma_sigma0 = 0.5` and
   `x0 = (3, 3)`, the 3σ exploration ellipse covers only
   `r ∈ [1.5, 4.5]²` — the true optimum at `(5.07, 7.74)` is ~4.7 σ
   away and is never sampled.

2. **Spurious cost minimum at `r2 → 0`.**  When `r2` collapses, the
   lattice becomes effectively one-dimensional; the χ(β) scanner
   latches onto a finite-size pseudo-peak ~6 % below the true T_c
   (`crit_lhs(β) = 0.853` instead of 1.0).  The reported cost at
   `(r1 ≈ 0.88, r2 ≈ 0.06)` is ~0.0009 — beating every legitimate
   eval, despite SNR of only 5.6 (vs ≥ 10 for genuine evals).
   This is a noise-driven artifact, not a true minimum.

3. **Cost surface is dominated by the shortest reference cycle.**
   For the (13, 16, −3, 3) reference, the three boundary cycles
   have lengths `v ≈ 14.73`, `u ≈ 17.69`, `w ≈ 11.79`.  The
   per-direction cost contributions at the *true* optimum are
   `Cv ≈ 0.0006`, `Cu ≈ 0.0022`, `Cw ≈ 0.0058` — `Cw` is 9.6× larger
   than `Cv` and dominates the total.  Combined with the L⁴-mean
   aggregation in `cost.py`, the optimizer is effectively minimising
   `Cw` alone, and any `(r1, r2)` configuration that makes `Cw` look
   small wins regardless of how it does on `Cv`/`Cu`.

### 8.2 Diagnosis of the cost-surface failure

Three independent tests confirmed that the small reference is the
primary issue:

* **`strategy_test.py`** (S1–S4).  Four CMA-ES configurations were
  run with 36 evals each at 5 000 production trajectories.  Even the
  strategy that *did* sample near the true basin (S2: `x0=(3,4)`,
  `σ0=3`, no bounds) recorded a *higher* cost there (0.00402 at
  `(4.82, 6.01)`) than in the middle region (0.00283 at
  `(2.42, 2.61)`).  The cost surface itself, as evaluated, has its
  minimum in the wrong place.

* **`ref_size_test.py`** (5×5 grid in `r ∈ [1, 9]²`, sizes 1× and 2×).
  Doubling the reference lattice from (13, 16, −3, 3) to
  (26, 32, −6, 6) — same 4 : 5 : 6 ratio, sides doubled — drops the
  per-direction `Cw` contamination at the near-true grid point
  `(5, 7)` from 0.00169 to **0.00054** (3× cleaner), inverts the
  `Cw / Cv` mean ratio from 10.5× to 0.25×, but does **not** by
  itself move the grid argmin away from the corner — the L⁴ mean
  and the 5 000-traj MC noise floor still dominate the surface.

* **`paraboloid_interp_test.py`** (analytic).  On the smooth target
  `f(x, y) = x² + y²` over the unit disk, comparing a 10 × 10 test
  grid against an `N_ref × N_ref` rotated reference grid:

  | N_ref | h = 2R/N | RMS_lin | RMS_lin / h² | RMS_cubic |
  |-------|----------|---------|--------------|-----------|
  |   5   | 0.40     | 8.9 e-2 | 0.55         | 1.6 e-2   |
  |  10   | 0.20     | 1.8 e-2 | 0.44         | 3.9 e-3   |
  |  20   | 0.10     | 4.1 e-3 | 0.41         | 6.0 e-4   |
  |  50   | 0.040    | 5.8 e-4 | 0.36         | 6.1 e-5   |
  | 100   | 0.020    | 1.2 e-4 | 0.30         | 9.9 e-6   |
  | 200   | 0.010    | 3.3 e-5 | 0.33         | 2.4 e-6   |

  `LinearNDInterpolator` residuals scale cleanly as `O(h²)` and
  `CloughTocher2DInterpolator` is ~10× better at every density.
  For the actual 4-5-6 reference, the shortest cycle has only
  ~12 lattice spacings, equivalent to `h ≈ 0.17` in the unit-disk
  normalisation, giving a predicted interpolation noise floor of
  `≈ 0.013` — comparable to or larger than the 1 e-3 signal at the
  true optimum.

  The same script also showed that **interpolating both reference
  and test** (the current `cost.py` behaviour) produces a *higher*
  L² value at the true minimum than at slightly perturbed points —
  a textbook spurious-minimum signature caused by correlated test
  interpolation error being worst at the actual data positions.
  Restricting interpolation to the reference and using the test
  values *exactly at the test lattice sites* eliminated this artifact
  entirely on the toy model.

### 8.3 Implemented fixes

* **Death-penalty hard bounds in `run_cmaes`**
  (`optimizer.py:run_cmaes`).  New keyword arguments
  `lower_bounds` / `upper_bounds` skip the MC for any candidate
  outside the box and assign a fixed penalty of 10.0.  Out-of-bounds
  offspring still participate in the covariance update so the mean
  is pushed away from the wall on subsequent generations.  Wired
  through `run.py` as `r_lower` / `r_upper` in CONFIG and
  `--r-lower R1 R2` / `--r-upper R1 R2` on the CLI.

* **Updated `4_5_6` preset.**  The shipped preset now uses
  `x0 = (3, 4)`, `cma_sigma0 = 3.0`, and
  `r_lower = (0.5, 0.5)` to combine all three cures by default.

### 8.4 Recommended workflow for cross-geometry runs

For any 4-5-6-style cross-geometry problem we now recommend a
**two-stage** run:

```bash
# Stage A — coarse exploration (~5 min): find the right basin
python run.py --preset 4_5_6 \
    --x0 3 4 --cma-sigma0 3.0 \
    --r-lower 0.5 0.5 --r-upper 12 12 \
    --n-traj-prod 5000 --max-evals 50 --run-name 456_stageA

# Stage B — refine in the winning basin (~30 min): tight convergence
# (replace 5.1 7.8 with the (r1, r2) Stage A reports as best)
python run.py --preset 4_5_6 \
    --x0 5.1 7.8 --cma-sigma0 0.7 \
    --r-lower 0.5 0.5 --r-upper 12 12 \
    --n-traj-prod 20000 --max-evals 40 --run-name 456_stageB
```

Stage A finds the basin reliably; Stage B converges with statistics
adequate to distinguish neighbouring `(r1, r2)` evaluations.

### 8.5 Open issues and next steps

The bounds + large σ₀ + two-stage workflow reliably keeps the
optimizer **out** of the spurious `r2 → 0` basin, but a dense grid
scan with the small reference still places the cost minimum in the
wrong corner.  The cost function itself needs improvement.  Three
fixes are slated based on the diagnostics above:

1. **Larger reference lattice.**  Doubling the reference linear
   size cleanly reduces the dominant `Cw` per-direction signal at
   the true optimum by ~3× (verified via `ref_size_test.py`).  The
   one-time MC cost is amortised across all subsequent runs against
   that reference.

2. **Test-native sampling in `cost.py`.**  Use the test correlator
   *exactly at the test lattice sites* along each boundary
   direction and only interpolate the reference correlator.  The
   paraboloid test demonstrated that this single change removes the
   spurious-minimum-at-truth artifact on smooth targets.  Side
   benefit: σ_cost becomes statistically honest (currently it is
   under-estimated by ~5× because 400 oversampled points are treated
   as independent).

3. **Restricted integration window** `t ∈ [0, 0.4]` in the
   per-direction L² integrals.  Avoids periodic-image contamination
   that arises whenever the sampled point is more than half a cycle
   length from the origin in the *short* reference direction.

4. **Cubic (Clough–Tocher) interpolation for the reference.**
   `CloughTocher2DInterpolator` gives ~10× lower interpolation
   residuals than `LinearNDInterpolator` at the same reference
   density — almost free upgrade once the reference is dense enough
   for cubic interpolation to be well-conditioned.

5. **Plain mean (`P = 1`) instead of L⁴ for direction aggregation.**
   The L⁴ mean amplifies any single direction's noise; on the
   small 4-5-6 reference it makes the optimizer essentially blind
   to `Cv` and `Cu`.  Plain-mean aggregation lets the smaller-cycle
   directions actually contribute to the cost.

These should be implemented as flag-controlled options on
`l2_cost(...)` so existing runs are reproducible while the new
behaviour can be A/B tested on cached MC data without re-running.

---

## 8.6 Cost-function zoo (2026-05-01)

`cost_zoo_test.py` is a standalone driver that runs nine candidate
cost functions plus a replicate-baseline post-processing step on the
cached `_betac_stability_test/` triple (10 replicates each at
`(k1,k2,k3) = (1,1,1)·k_c`, `(2,1,1)·k_c`, `(5,1,1)·k_c`, all
`L = 12` equilateral so the cost only sees physics differences).

The discrimination metric is the **Z-score**

$$Z \;=\; \frac{\mu_{\text{mismatched}} - \mu_{\text{matched}}}{\sigma_{\text{matched}}}$$

a higher `|Z|` means the optimizer can see the true `(r1, r2)`
mismatch above MC noise.

Crucially, the test loops over **three reference choices** (eq, a2,
a5).  An honest distance-like cost gives roughly symmetric `|Z|`
across all three; a cost that only fires when the reference is
isotropic (`ref = eq`) is exposed as an anisotropy probe rather than
a true distance.

### 8.6.1 Headline results

`|Z|` ranges across all 6 ref/mismatch combinations:

| candidate            | best `|Z|` | worst `|Z|` | verdict                                    |
|----------------------|-----------:|------------:|--------------------------------------------|
| **affine_fit**       |        158 |        12.6 | ★ winner — robust across all 3 ref choices |
| chi2_native          |         43 |        0.06 | ~2–3× boost over native, same MC cost      |
| native (baseline)    |         31 |        0.04 | reference                                  |
| cycle_loop           |         31 |        0.04 | tracks native (rescaled same signal)       |
| spectral             |         25 |        0.02 | rescaled native (perfect log-log corr.)    |
| angular_iso          |         25 |        0.10 | only fires at `ref=eq` → isotropy probe    |
| moment_r2            |        6.8 |       −0.40 | sign-flips, unreliable                     |
| eff_xi               |        5.2 |        0.20 | weak signal everywhere                     |
| **cycle_var**        |    220 000 |        1.6  | pure isotropy probe — fake winner          |
| baseline_subtract    |         29 |       −0.07 | does not fix r2 mismatch                   |

Two ideas are **N/A** from cached data: jackknife (needs per-config
block dumps) and off-diagonal covariance (needs per-config full
all-to-all matrix).  Both would require re-running the simulator
with extra output hooks.

Output artefacts: `results/_cost_compare/{zoo_summary.json,
zoo_compare.png, zoo_distributions.png}`.

### 8.6.2 `affine_fit` — detailed explanation

For each ref/test pair, take the intersection of `(m, n)` keys present
in both all-to-all dicts and stack into vectors `G_r`, `G_t`.  Solve
the 2-parameter linear least-squares problem

$$(a^\*, b^\*) \;=\; \arg\min_{a,b}\;\bigl\| a\,G_t + b\,\mathbf{1} - G_r \bigr\|_2^2$$

(one `np.linalg.lstsq` call on `A = [G_t, 1]` against `G_r`) and
report the residual mean-square as the cost:

$$C_{\text{affine}} \;=\; \frac{1}{N}\,\bigl\| a^\*\,G_t + b^\*\,\mathbf{1} - G_r \bigr\|_2^2.$$

The cost is **invariant under the two transformations that don't
carry information about `(r1, r2)`**:

1. **Multiplicative rescaling** (slope `a`).  Two MC replicates at
   identical couplings can disagree on the overall amplitude of
   `G(m, n)` because of: `β_c` jitter (a 10⁻³ shift in the GC-fit
   critical β rescales `G` by an O(1) factor at finite `L`),
   finite-volume normalisation differences, and reweighting-overlap
   wobble.  All of these live entirely in `a` and leave the *shape*
   of `G(m, n)` untouched.
2. **Additive offset** (intercept `b`).  The connected correlator
   should asymptote to zero, but at finite `L` each replicate has a
   residual `⟨m⟩²` subtraction error that varies replicate-to-
   replicate.  This is a constant offset in `G` and is absorbed by
   `b`.

By projecting both nuisance directions out **before** measuring
residual, the matched-pair distribution shrinks dramatically (matched
std drops from `1.4e-2` for native to `9.4e-6` for affine_fit — a
~1500× MC-noise suppression), while genuine geometric mismatch
(which is a *shape* difference, not a scale difference) survives.
That is why `|Z|` jumps by 1–2 orders of magnitude across **all six**
mismatch tests.

Caveats: it discards absolute amplitude (use `chi2_native` if you
want amplitude in the cost); it needs ≥ 4 shared `(m, n)` keys
(trivial for same-geometry runs); cross-geometry runs need overlapping
integer lattice points (rare failure case).

### 8.6.3 `cycle_var` — what it is and why it's a fake winner

For each replicate, sum `G_conn(m, n)` along the integer lattice
points on each of the three torus boundary period vectors
`(Lx, Ty)`, `(Tx, -Ly)`, `(-Lx-Tx, Ly-Ty)`:

$$\Phi_d \;=\; \sum_{(m,n)\,\in\,\text{cycle }d} G_{\text{conn}}(m,n),
   \qquad d \in \{v, u, w\}.$$

Then take the **sample variance over the three scalars**
`{Φ_v, Φ_u, Φ_w}`: large when the cycle integrals disagree
(anisotropic correlator), zero when they agree (perfectly isotropic).
The pair cost is the squared difference of this single number
between ref and test:

$$C_{\text{cycle\_var}} \;=\; \bigl( \operatorname{Var}\{\Phi_d^{\text{test}}\} - \operatorname{Var}\{\Phi_d^{\text{ref}}\} \bigr)^2.$$

This is a **single-sample, reference-free statistic** — each
replicate has its own anisotropy number; the cost only compares the
two numbers.  When `ref = eq` (true isotropy), the ref variance is
≈ 0 and any anisotropic test sticks out enormously (`|Z| ~ 10⁴–10⁵`).
When `ref = a2` (already anisotropic), both variances are macroscopic
and their difference is small — `|Z|` collapses to 1.6.  In production
the reference is almost never perfectly isotropic, so `cycle_var`
mostly measures *distance-to-isotropic-point*, not *distance to
the actual reference*.  Useless for navigating `(r1, r2)` space
toward an anisotropic target.

This is the canonical failure mode the 3-reference test was designed
to expose: any future cost-function candidate must be tested against
an anisotropic reference, not just the equilateral one.

### 8.6.4 Promotions to `cost.py`

1. **`cost_mode='affine_fit'`** — promoted (2026-05-02).  Wired into
   the `l2_cost(...)` dispatcher with an optional `affine_kwargs`
   argument:
   * `affine_kwargs={'mode': 'all_to_all'}` — intersects the
     `(m, n)` keys present in both ref and test (same-geometry case).
   * `affine_kwargs={'mode': 'test_native'}` — samples `G_t` on the
     integer test-boundary lattice sites and looks up `G_r` via
     `tile_interp` (cross-geometry; this is the **default** so the
     dispatcher "just works" for problems like the 4-5-6 triangle).

   Returns the same 4-tuple `(cost, σ_cost, per_dir, per_dir_σ)` as
   the other modes; per-direction wedges are angular bins (v, u, w)
   for diagnostic parity with `test_native`.

2. **`cost_mode='chi2_native'`** — *still pending*.  Drop-in upgrade
   of `test_native` at zero extra MC cost (divides each residual by
   `sqrt(σ_t² + σ_r²)`); ~2–3× boost across the board.

### 8.6.5 CMA-ES validation: `affine_fit` vs `test_native` (2026-05-02)

`cost_mode_test.py` extended with a third run (mode `C`,
`cost_mode='affine_fit'`) so the workflow-level CMA-ES comparison can
include the new dispatcher branch.  A single-seed run on the 4-5-6
problem (`--modes B C --max-evals 30 --n-workers 6 --seed 42`,
ref=(13,16,−3,3), test=(16,16,0,0), x0=(5.0, 6.0), σ0=2.0,
popsize=6, n_traj=10 000) gave:

| mode             | best (r1, r2)   | best cost   | β_c    | dist→true | first<1.0 |
|------------------|-----------------|-------------|--------|-----------|-----------|
| B `test_native`  | (3.136, 2.387)  | 2.020 × 10⁻³ | 0.121  | 5.69      | >budget   |
| C `affine_fit`   | (3.441, 3.135)  | 2.018 × 10⁻³ | 0.104  | 4.89      | >budget   |

**Both runs failed to converge** within budget — both stalled in a
`(r1, r2) ≈ (3, 3)` corner, never approaching the true
`(5.07, 7.74)`, with essentially identical `best_cost`.

The cost-landscape diagnostic, however, *does* show the expected
`affine_fit` advantage:

| metric                                | B (test_native) | C (affine_fit) |
|---------------------------------------|-----------------|----------------|
| dynamic range (max / min over 30 evals) | 2.31×          | 2.40×          |
| mean SNR (cost / σ)                   | 20.8            | 21.0           |
| dist of closest-to-true eval          | 1.05            | 1.34           |
| (cost at closest eval) / (global min) | **1.59×**       | **1.43×**      |

At the eval that landed nearest to the true solution, `affine_fit`
rates it as only 1.43× the run's cost minimum, vs `test_native`'s
1.59× — i.e. `affine_fit` is more honest about the true region,
matching the zoo-test prediction.  But on this geometry the cost
surface itself is nearly degenerate (a flat valley of ~2.0 × 10⁻³
spanning very different `(r1, r2)`), and neither cost can break the
symmetry at this budget.

**Take-aways:**

1. The zoo-test-level discrimination win (|Z| = 12–158 for
   `affine_fit` vs 0.04–30 for `test_native` on noise tests) **does
   not automatically translate** to a CMA-ES convergence advantage
   at production budgets when the cost surface is intrinsically flat.
2. Single-seed 30-eval comparisons are insufficient — the cost
   surface dominates over the cost choice at this scale.
3. To actually validate the affine_fit advantage end-to-end, need
   one of: (a) multi-seed sweep, (b) larger budget
   (`--max-evals 80` or `popsize 12`), or (c) a sharper start
   `x0=(5.0, 7.5)` with `σ0=0.8`.

Artifacts:
* [`results/cost_mode_test/B/result.json`](results/cost_mode_test/B/result.json)
* [`results/cost_mode_test/C/result.json`](results/cost_mode_test/C/result.json)
* [`results/cost_mode_test/analysis_BC_seed42.md`](results/cost_mode_test/analysis_BC_seed42.md)

---

## 8.7 Modular-plane testbed (2026-05-01)

A noise-free, fully analytic CMA-ES + cost-discrimination harness that
isolates the **cost surface** from the optimizer and from MC noise.
Builds on `cft_torus.ising_torus_F` (Jacobi-theta Ising spin-spin
two-point function) so the "truth" is a closed-form expression and
correlator data can be generated at arbitrary modular parameter τ
without running any C++ simulation.

### 8.7.1 Why this testbed

The 4-5-6 real-MC convergence test (§8.5) showed CMA-ES drifting
*away* from truth.  The 1-D s-scan and 2-D τ-scan
([`results/_native_triage/zoo_ising_tau_grid.png`](results/_native_triage/zoo_ising_tau_grid.png))
revealed the production cost (`_l2_cost_test_native`) and four of the
five [`cost_zoo`](cost_zoo.py) variants place their minimum at the
boundary of any τ-window, never at truth.  Diagnosis: the lattice
boundary-path / `gcd`-stride sampling makes the per-direction cost
depend on the **test geometry's** stride pattern rather than the
underlying CFT, so the cost minimum is decoupled from the true τ.

### 8.7.2 Three-direction τ-aligned profiles

[`modular_data.py`](modular_data.py) replaces lattice boundary paths
with **three torus geodesics aligned with the modular parameter**:

| label   | direction in z-plane | physical meaning      |
|---------|----------------------|-----------------------|
| `re`    | along `ω₁` axis      | the **Re τ** cycle    |
| `im`    | perpendicular axis (`i`) | the **Im τ** cycle    |
| `diag`  | along `ω₂ − ω₁`      | the **Im τ − Re τ** cycle (after one twist) |

For each τ on the modular grid we evaluate `ising_torus_F` along each
direction at integer physical distances `s = 1, 2, …, K_min` (with
`K_min` bounded by the shortest half-period).  The result is three
short vectors `G[s]` per cell, all referred to the **same** physical
units (lattice spacing 1).  No `gcd` strides, no boundary-path tuples,
no `t = k/N` reparameterization.

Five cost modes operate on these profiles
([`modular_data.threedir_cost`](modular_data.py)):

| mode        | per-direction kernel                                   |
|-------------|--------------------------------------------------------|
| `l2`        | `mean( (G_t − G_r)² )`                                  |
| `relative`  | `mean( ((G_t − G_r) / G_r)² )`                          |
| `log`       | `mean( (log G_t − log G_r)² )`                          |
| `effmass`   | `mean( (m_eff_t − m_eff_r)² )` , `m_eff_k = −log(G_{k+1}/G_k)` |
| `pmean4`    | `( mean_d C_d⁴ )^(1/4)` over directions of `l2` per-dir |

### 8.7.3 Closed-loop CMA-ES navigation

[`modular_testbed.py`](modular_testbed.py) implements the closed loop:

1. Pre-compute three-direction profiles on an `n × n` grid in
   `(Re τ, Im τ)` and pickle to
   [`results/_native_triage/modular_cache_n*_K*.pkl`](results/_native_triage/).
2. Pick `τ_truth` (default = `τ` of the production reference geometry,
   `(13, 16, −3, 3)` → `τ_truth ≈ −0.901 + 0.794 i`) and `τ_start`
   (default = `τ` of the test geometry, `(16, 16, 0, 0)` → `τ_start ≈
   −0.500 + 0.866 i`).
3. For each cost mode, run a minimal 2-D `(μ/μ_w, λ)`-CMA-ES
   (inline in [`modular_testbed.cmaes`](modular_testbed.py)) in the
   continuous `(Re τ, Im τ)` plane; each evaluation snaps to the
   nearest grid cell and computes `threedir_cost` against
   `ref_profile` at `τ_truth`.
4. Optional `--noise σ` injects multiplicative Gaussian on every
   profile entry with a deterministic seed per
   `(mode, seed-base, cell, call_idx)` so the cost is reproducible
   while reflecting realistic MC scatter (re-sampled per call).

```bash
# noise-free
python modular_testbed.py --n 21 --K-max 8 --seeds 3 --max-evals 120

# 5% relative noise per profile entry
python modular_testbed.py --n 21 --K-max 8 --seeds 5 --noise 0.05
```

Outputs:
[`results/_native_triage/modular_testbed_<tag>.{json,png}`](results/_native_triage/) —
2 × 3 grid of cost heatmaps with `truth`/`start`/`grid min` markers
and per-seed CMA-ES trajectories overlaid.

### 8.7.4 Headline results (smoke n=11, K_max=6, seeds=2-3)

**Noise-free:**

| mode       | success | mean dist → truth | best dist | rank |
|------------|--------:|------------------:|----------:|-----:|
| `l2`       |    100% |             0.037 |     0.023 |  ★ 1 |
| `relative` |     50% |             0.084 |     0.023 |    2 |
| `log`      |     50% |             0.084 |     0.023 |    3 |
| `effmass`  |     50% |             0.468 |     0.092 |    4 |
| `pmean4`   |      0% |             0.155 |     0.120 |    5 |

`l2` on three-direction profiles **finds the truth in every seed**.
This is the opposite of the lattice-boundary-path version (§8.6 and
§8.5), where every cost candidate placed its global minimum at the
edge of the τ-window.  The fix is the *sampling reformulation*, not
the cost kernel.

**With 5% multiplicative noise per entry** (3 seeds, max 80 evals):
all modes lose convergence (`success = 0%`, `mean dist ≈ 0.4`).  The
basin width shrinks below the noise floor at this `K_max=6` /
`n=11` setting; raising `K_max` and / or averaging multiple noise
realizations per call (a "fewer-but-cleaner-evals" budget split) is
the obvious next step.

### 8.7.5 Roadmap

1. **Production sweep**: re-run §8.7.4 at `n=21, K_max=8, seeds=10`
   to get statistically meaningful success rates per cost mode at
   noise levels `σ ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.10}`.
2. **Frozen-noise mode**: add a `--noise-frozen` flag that perturbs
   each cell *once* (seed = cell index) so the cost surface has a
   fixed roughness rather than re-sampled noise per call; this
   isolates "static landscape roughness" from "stochastic eval
   noise."
3. **Block-averaged eval**: add `--noise-blocks K` that averages
   `K` noise samples per call (effective σ → σ/√K), letting us
   trade evals for noise reduction inside the CMA-ES loop.
4. **Port back to lattice cost**: once a cost mode is validated on
   the modular testbed, wire its three-direction kernel into
   [`cost.py`](cost.py) by replacing the `boundary_paths()` /
   `_direction_lattice_steps()` machinery with τ-aligned sampling
   on the actual MC data dict (using bilinear interpolation of
   `_tile_interp` at the three direction unit vectors).
5. **Multi-truth survey**: instead of fixing `τ_truth` to the
   `(13, 16, −3, 3)` reference, pick 6–10 representative truths
   spanning the modular fundamental domain and rank costs by
   *aggregate* success.  This guards against cost modes that only
   work near one corner of `(Re τ, Im τ)`.

Artifacts:
* [`modular_data.py`](modular_data.py) — three-direction profile builder + costs + noise
* [`modular_testbed.py`](modular_testbed.py) — cache builder + CMA-ES navigation driver
* [`results/_native_triage/modular_testbed_noise0.png`](results/_native_triage/modular_testbed_noise0.png) — noise-free 5-panel landscape + trajectories
* [`results/_native_triage/modular_testbed_noise0p05.png`](results/_native_triage/modular_testbed_noise0p05.png) — 5% noise version

---

## 9. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `g++: command not found` during build | No C++ compiler | `apt install build-essential` / Xcode CLT / MSVC |
| `--live` raises `tkinter.TclError: no display name` | SSH without X forwarding | drop `--live`; use `--save-frames` and copy PNGs locally |
| Reweight fails on every eval | Bracket too narrow for current β_c | widen `beta_seed` in CONFIG, or accept the fallback path |
| Long path errors on Windows | `mc_scratch/` filename overflow | shorten `--run-name` (≤ ~15 chars); avoid OneDrive paths |
| Workers crash with `pickling error` | Custom callback object passed to pool | the production pool keeps workers headless — don't attach plotters/dashboards to the worker `Evaluator` |

---

## 10. Programmatic use

`run_optimizer(cfg) → summary_dict` in `run.py` is the library entry
point.  Build a CONFIG dict in your own script, optionally run a
parameter sweep, and call it directly:

```python
import os, sys
sys.path.insert(0, "K_from_Optimizer_Production")
from run import CONFIG, run_optimizer

for sigma0 in (0.05, 0.10, 0.20):
    cfg = dict(CONFIG, cma_sigma0=sigma0,
               run_name=f"sigma_{sigma0:.2f}", live=False)
    summary = run_optimizer(cfg)
    print(sigma0, summary["best_cost"])
```

---

## 11. Provenance

The C++ simulator and the bulk of the Python pipeline (`mc_engine.py`,
`evaluator.py`, `optimizer.py`, `reweight.py`, `visualization.py`,
`dashboard.py`, `cost.py`, `betac_cache.py`) are copied verbatim from
`K_from_Optimizer_Speedup_Upgrade/`.  This directory adds:

- `prod_runtime.py` — the reweight→3-pass-fallback patch and a worker
  pool that installs the patch in every child process,
- `live_viewer.py` — the Tk live PNG viewer,
- `run.py` — a single-file, fully-CLI-overridable driver,
- `run_overnight.sh` + `check_overnight.py` — the validation harness.
- `mc_engine.find_beta_c_multidonor` — multi-donor 1-pass reweighting
  (added here; not in upstream).
- `mc_engine.find_beta_c_multidonor_2pass` — interval-based 2-pass
  multi-donor reweighting (added here; not in upstream).
- `test_multidonor_vs_3pass.py` — validation script comparing 1-pass
  multi-donor against the 3-pass GC scan.
- `test_multidonor_2pass.py` — three-way comparison: 1-pass MD vs
  2-pass MD vs 3-pass GC scan on a shared trajectory budget.
- `fss_betac_scan.py` — isotropic finite-size-scaling validation of
  the β_c finder (fits β_c(L) = β_c(∞) + a·L^{−1/ν} on the equilateral
  triangular Ising model; verifies ν ≈ 1 and β_c(∞) ≈ ln(3)/4).
- `fss_aniso_betac.py` — per-anisotropy FSS scan that fits β_c(L)
  with a/L (and a/L + b/L²) for each (k1, k2, k3) case, comparing
  the extrapolated β_c(∞) to the exact triangular-Ising sinh-rule
  value.
- `stress_aniso_betac.py` — β_c finder stress test on a fixed L
  across anisotropic couplings, with FSS-corrected residuals in σ
  units.
- `stress_workflow.py` — end-to-end workflow stress test on
  isotropic (untwisted) ref = test = (L,L,0,0) lattices for several
  L; ground truth (r1,r2) = (1,1).  Reports best, CMA-mean and
  trailing-5 distances from the answer.
- `stress_twisted.py` — same harness as `stress_workflow.py` on
  twisted-equilateral tori (Lx, Ly, Tx, Ty) = (L, L+T, T, T) where
  the three cycle lengths are equal; ground truth still (1,1).
- `strategy_test.py` — comparative harness running four different
  CMA-ES strategies (baseline, large σ, large σ + bounds,
  closer-start + bounds) on the 4-5-6 cross-geometry problem.
  Used to demonstrate the failure modes documented in §8.
- `ref_size_test.py` — sweeps reference lattice size (1×, 2×, 3×
  the 4-5-6 base geometry) and runs both a uniform `(r1, r2)` grid
  evaluation *and* a CMA-ES optimisation against each reference.
  Used to quantify how much of the false minimum is caused by
  reference-side discretisation.
- `paraboloid_interp_test.py` — analytic conceptual test isolating
  the geometric contribution of `LinearNDInterpolator` /
  `CloughTocher2DInterpolator` residuals on the smooth target
  `f(x, y) = x² + y²` over the unit disk, with the reference
  rotated by 30° relative to the test grid.  Confirms the expected
  `O(h²)` scaling and motivates the cost-function changes proposed
  in §8.5.
- Hard bounds in `optimizer.run_cmaes` — new `lower_bounds` /
  `upper_bounds` keyword arguments implementing a death penalty
  for out-of-bounds CMA-ES offspring.  Wired through `run.py`
  CONFIG (`r_lower`, `r_upper`) and CLI (`--r-lower`, `--r-upper`).

Bug fixes that land in the upstream copies should be `cp`-mirrored back
into this directory (the upstream `parallel.py`, `evaluator.py`, etc.
remain reusable as-is — `prod_runtime.py` only monkey-patches
`mc_engine.find_beta_c_reweight`).
