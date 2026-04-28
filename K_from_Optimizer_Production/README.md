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

## 8. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `g++: command not found` during build | No C++ compiler | `apt install build-essential` / Xcode CLT / MSVC |
| `--live` raises `tkinter.TclError: no display name` | SSH without X forwarding | drop `--live`; use `--save-frames` and copy PNGs locally |
| Reweight fails on every eval | Bracket too narrow for current β_c | widen `beta_seed` in CONFIG, or accept the fallback path |
| Long path errors on Windows | `mc_scratch/` filename overflow | shorten `--run-name` (≤ ~15 chars); avoid OneDrive paths |
| Workers crash with `pickling error` | Custom callback object passed to pool | the production pool keeps workers headless — don't attach plotters/dashboards to the worker `Evaluator` |

---

## 9. Programmatic use

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

## 10. Provenance

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

Bug fixes that land in the upstream copies should be `cp`-mirrored back
into this directory (the upstream `parallel.py`, `evaluator.py`, etc.
remain reusable as-is — `prod_runtime.py` only monkey-patches
`mc_engine.find_beta_c_reweight`).
