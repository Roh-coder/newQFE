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
run.py             ← single CLI entry point + the editable CONFIG dict
prod_runtime.py    ← reweight→3-pass fallback patch + parallel pool
live_viewer.py     ← Tk window that auto-reloads the optimizer PNG
run_overnight.sh   ← 4-run validation harness (legacy / RW / RW+FB / parallel)
check_overnight.py ← post-run report generator (writes OVERNIGHT_REPORT.md)
Makefile           ← builds bin/ising_tri_twisted_parallelogram
src/, include/     ← C++ simulator source
evaluator.py       ← per-eval (r1, r2) → (cost, σ_cost) [from upstream]
mc_engine.py       ← run_simulator, find_beta_c, find_beta_c_reweight
optimizer.py       ← run_cmaes (hand-coded), run_nelder_mead, run_bo
parallel.py        ← legacy GenerationPool (kept for reference)
reweight.py        ← Reweighter, CombinedReweighter
visualization.py   ← OptimizerPlotter (4-panel PNG)
dashboard.py       ← rich-based terminal dashboard
cost.py            ← L⁴ power-mean curve-collapse cost
betac_cache.py     ← optional persistent β_c interpolation cache
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
| **3-pass fallback**     | Catches reweight failures (bracket too narrow, peak outside validity window) and re-runs the legacy GC scan transparently | Makes S4 safe to leave on for far-from-equilateral starts | always on, audited via `fallback_log.jsonl` |
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

## 4. The visualizer

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

## 5. Overnight validation suite

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

## 6. Switching geometries / starting points

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

---

## 7. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `g++: command not found` during build | No C++ compiler | `apt install build-essential` / Xcode CLT / MSVC |
| `--live` raises `tkinter.TclError: no display name` | SSH without X forwarding | drop `--live`; use `--save-frames` and copy PNGs locally |
| Reweight fails on every eval | Bracket too narrow for current β_c | widen `beta_seed` in CONFIG, or accept the fallback path |
| Long path errors on Windows | `mc_scratch/` filename overflow | shorten `--run-name` (≤ ~15 chars); avoid OneDrive paths |
| Workers crash with `pickling error` | Custom callback object passed to pool | the production pool keeps workers headless — don't attach plotters/dashboards to the worker `Evaluator` |

---

## 8. Programmatic use

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

## 9. Provenance

The C++ simulator and the bulk of the Python pipeline (`mc_engine.py`,
`evaluator.py`, `optimizer.py`, `reweight.py`, `visualization.py`,
`dashboard.py`, `cost.py`, `betac_cache.py`) are copied verbatim from
`K_from_Optimizer_Speedup_Upgrade/`.  This directory adds:

- `prod_runtime.py` — the reweight→3-pass-fallback patch and a worker
  pool that installs the patch in every child process,
- `live_viewer.py` — the Tk live PNG viewer,
- `run.py` — a single-file, fully-CLI-overridable driver,
- `run_overnight.sh` + `check_overnight.py` — the validation harness.

Bug fixes that land in the upstream copies should be `cp`-mirrored back
into this directory (the upstream `parallel.py`, `evaluator.py`, etc.
remain reusable as-is — `prod_runtime.py` only monkey-patches
`mc_engine.find_beta_c_reweight`).
