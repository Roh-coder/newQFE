# K_from_Optimization (standalone, plug-and-play)

Find anisotropic Ising couplings `(k1, k2, k3)` on a target lattice that
reproduce the two-point correlator of a chosen reference lattice (using
`k3 = 1` and optimising the ratios `r1 = k1/k3`, `r2 = k2/k3`).

This directory is **self-contained**.  Drop it anywhere, run `make`,
edit the `CONFIG` block in `run.py`, and run `python run.py`.

---

## Quick start

```bash
# 1. Build the C++ Monte Carlo simulator (one-time)
make

# 2. Install Python deps (one-time, in your venv of choice)
pip install -r requirements.txt

# 3. Edit CONFIG at the top of run.py to set lattice sizes / statistics

# 4. Run
python run.py
```

Outputs go to `results/`:

- `results/_reference_L<N>/`         — cached reference correlator (built once per L)
- `results/<optimizer>/eval_log.jsonl` — one JSON line per evaluation
- `results/<optimizer>/summary.json`   — final summary including the chosen (r1, r2)
- `results/<optimizer>/frames/`        — live PNG visualisation, updated every eval

where `<optimizer>` is `nelder_mead` or `cmaes` depending on `CONFIG["optimizer"]`.

---

## What it does

Three phases:

1. **Reference build** (one-time per reference geometry): runs a 3-pass Gram-Charlier
   susceptibility-peak scan to find β_c on an `Lref × Lref` equilateral
   lattice (k1=k2=k3=1), then a long production MC at β_c to measure the
   all-to-all two-point correlator G_ref.  Cached to disk.

   The 3-pass scan starts with `ref_scan_n_coarse` evenly-spaced β points
   across `ref_beta_seed` (translating the window up to
   `ref_scan_max_shifts` times if the peak lies on an edge), then runs
   three GC-fit refinement passes whose brackets shrink to a multiple of
   the fitted σ each time.  The per-eval scan uses the same algorithm
   with the `scan_*` keys.  When `scan_jackknife=True` a leave-one-out
   estimate of σ(β_c) is computed at no MC cost (~10–200 ms of GC
   refits per scan) and surfaced in the eval log as `beta_c_sigma`.

2. **Per-evaluation Monte Carlo** (called by the optimizer with a candidate
   `(r1, r2)`):
   - Run a 3-pass β_c scan on the test lattice with couplings `(r1, r2, 1)`
     (warm-started from the previous evaluation's β_c).
   - Run a production MC at the new β_c to measure G_test.
   - Compute the cost.

3. **Cost function** (`cost.py`):
   Per-direction L² mismatch on each of the three torus boundary directions,
   $$ C_d = \int_0^1 \bigl[G_\text{test}^{(d)}(t) - G_\text{ref}^{(d)}(t)\bigr]^2 \, dt, \qquad d \in \{v, u, w\}, $$
   integrated by trapezoid rule on `n_samples = 400` points, then aggregated
   across directions with a **p = 4 power mean**:
   $$ C(r_1, r_2) = \Bigl( \tfrac{1}{3} \sum_{d \in \{v,u,w\}} C_d^{\,4} \Bigr)^{1/4}. $$
   The L⁴ mean penalises anisotropy: a single direction that fits poorly
   dominates the total, so all three correlators must match.  Propagated
   1-σ uncertainty σ_C is also computed from the per-point MC errors so the
   optimizer can monitor the SNR.

4. **Optimization** (`optimizer.py`):
   Two backends, selected by `CONFIG["optimizer"]`:

   - `"nelder_mead"` — hand-coded Nelder-Mead simplex in `(r1, r2)`.
     Stops on `max_evals`, simplex diameter < `nm_xatol` *and* cost spread
     < `nm_fatol`, or `indist_stop_snr`.
   - `"cmaes"` — hand-coded (μ/μ_w, λ)-CMA-ES (pure NumPy, no `cma`
     dependency) with rank-1 + rank-μ covariance updates and CSA step-size
     adaptation.  Stops on `max_evals`, σ·max(D) < `cma_tolx`, generation
     cost spread < `cma_tolfun`, or `indist_stop_snr`.  More robust than NM
     on the flat L⁴ ridges produced by the power-mean cost.

   Both backends share the same adaptive-statistics logic: when SNR <
   `snr_target`, `n_traj_prod` is multiplied by 1.5 (capped at
   `snr_max_traj_factor × starting value`).

---

## Configuration

All knobs are at the top of `run.py` in the `CONFIG` dictionary:

| Key | Meaning | Typical |
|-----|---------|---------|
| `ref_Lx, ref_Ly`               | Reference lattice dimensions                              | 24, 24 |
| `ref_Tx, ref_Ty`               | Reference twist parameters (0 for untwisted)              | 0, 0 |
| `test_Lx, test_Ly`             | Test lattice dimensions                                   | match ref |
| `test_Tx, test_Ty`             | Test twist parameters (0 for untwisted square)            | 0, 0 |
| `ref_n_traj`                   | Production trajectories at the reference β_c              | 60 000 |
| `ref_scan_n_traj_coarse/fine`  | β_c scan trajectories per point (reference)               | 4 000 / 10 000 |
| `ref_beta_seed`                | Initial bracket `(β_lo, β_hi)` for the reference scan     | (0.20, 0.32) |
| `ref_scan_n_coarse`            | Points in pass 0 of the reference β_c scan                | 11 |
| `ref_scan_n_refine[/2/3]`      | Extra points in passes 1 / 2 / 3 (each pass uses N+2 pts) | 5 / 5 / 5 |
| `ref_scan_max_shifts`          | Max window translations if pass-0 peak lies on the edge   | 4 |
| `n_traj_prod`                  | Production trajectories per optimizer evaluation          | 30 000 |
| `n_traj_scan_coarse/fine`      | Per-point trajectories during each per-eval β_c scan      | 4 000 / 10 000 |
| `beta_seed`                    | Initial bracket for the very first eval (warm-started after)  | (0.20, 0.32) |
| `scan_n_coarse`                | Points in pass 0 of each per-eval β_c scan                | 11 |
| `scan_n_refine[/2/3]`          | Extra points in per-eval passes 1 / 2 / 3                 | 5 / 5 / 5 |
| `scan_max_shifts`              | Max window translations during a per-eval scan            | 4 |
| `scan_jackknife`               | Leave-one-out jackknife on β_c (adds ~ms per scan, no extra MC) | False |
| `x0`                           | Starting `(r1, r2)`                                       | (1.0, 1.0) |
| `max_evals`                    | Hard cap on optimizer evaluations                         | 30 |
| `optimizer`                    | `"nelder_mead"` or `"cmaes"`                              | `"nelder_mead"` |
| `nm_sigma0`                    | Initial NM simplex side length                            | 0.10 |
| `nm_xatol / nm_fatol`          | NM convergence tolerances                                 | 0.005 / 1e-6 |
| `nm_shrink`                    | Simplex shrink coefficient                                | 0.75 |
| `cma_sigma0`                   | Initial CMA-ES global step size σ₀                        | 0.10 |
| `cma_popsize`                  | Evals per generation λ (≈ 4 + 3·ln N, so 6 for N=2)        | 6 |
| `cma_max_gens`                 | Hard cap on generations (0 = no cap, only `max_evals` limits) | 0 |
| `cma_tolx`                     | Stop when σ·max(D) drops below this                       | 0.005 |
| `cma_tolfun`                   | Stop when generation cost spread drops below this         | 1e-6 |
| `cma_seed`                     | RNG seed for reproducibility (`None` = random)            | None |
| `snr_target`                   | Boost stats when SNR drops below this (0 = disable)       | 2.0 |
| `snr_max_traj_factor`          | Stats boost cap (× starting `n_traj_prod`)                | 4 |
| `indist_stop_snr`              | Stop when running best SNR drops below this               | 1.0 |
| `save_every`                   | Write optimizer PNG every N evals                         | 1 |
| `no_vis`                       | Disable PNG frames entirely (faster)                      | False |
| `dashboard`                    | Live terminal dashboard (rich); PNGs are still saved      | True |

---

## Files

```
run.py            ← edit CONFIG and run this
cost.py           ← L⁴ power-mean cost over three boundary directions
evaluator.py      ← (r1, r2) → (cost, σ_cost): wraps β_c scan + production MC
mc_engine.py      ← Gram-Charlier β_c finder + simulator wrapper
optimizer.py      ← hand-coded Nelder-Mead and CMA-ES backends
dashboard.py      ← live terminal dashboard (rich)
visualization.py  ← optional PNG frames (4-panel optimizer view + β_c scans)
Makefile          ← build the C++ simulator
src/              ← C++ Ising MC source
include/          ← C++ headers
bin/              ← built simulator binary lands here
```

---

## Live terminal dashboard

By default `CONFIG["dashboard"] = True` and you'll see something like
this at the bottom of your terminal, refreshing as the run progresses:

```
╭─ K_from_Optimization ───────────────────────────────────────────────╮
│ ref:  L=24×24  β_c=0.262914  n_traj=60000                           │
│ test: 24×24  Tx=0 Ty=0    method=nelder_mead  max_evals=30  wall=…  │
╰──────────────────────────────────────────────────────────────────────╯
╭─ β_c scan ──────────────────────────────────────────────────────────╮
│ eval 003  r1=+1.0500 r2=+0.9500    bracket=[+0.2103, +0.3103]       │
│ pass 2   pts=18   traj=146000                                       │
│   this pass: ████████████░░░░░░░░░░░░ 4/7                           │
│   current β_c estimate : 0.263812                                   │
│   argmax χ              : 0.263000    χ_max = 1.234                 │
╰──────────────────────────────────────────────────────────────────────╯
╭─ optimizer (nelder_mead) ───────────────────────────────────────────╮
│  ev    r1       r2       β_c       cost       σ        SNR  status  │
│ ★001  +1.0000  +1.0000  0.262914  1.23e-04   3.4e-05   3.6  ok      │
│  002  +1.1000  +1.0000  0.258091  8.42e-04   2.1e-04   4.0  ok      │
│  003  +0.9500  +1.0500  0.269140  9.21e-05   2.8e-05   3.3  ok      │
│ best: r1=+1.0000  r2=+1.0000  β_c=0.262914  cost=1.23e-04  SNR=3.6  │
│ runs: 3/30    total_traj=287000    wall_in_evals=232.4s             │
╰──────────────────────────────────────────────────────────────────────╯
```

All `print()` output from the MC engine still scrolls above the dashboard
so you keep the full log.  PNG snapshots of the optimizer state
(`results/nelder_mead/frames/optimizer_nelder_mead.png`) are written to
disk in parallel — set `CONFIG["no_vis"] = True` to disable the PNG
output, or `CONFIG["dashboard"] = False` to disable the live terminal
view (in which case the per-eval β_c-scan PNG is also written).

---

## Tips

- **First run is slow**: building the reference takes a few minutes
  (`ref_n_traj` trajectories on the reference lattice).  It is cached, so
  subsequent runs at the same reference geometry reuse it.

- **Per-evaluation cost** scales as `n_traj_prod × (test_Lx · test_Ly)`, plus a few
  thousand trajectories for the β_c scan.  At `Lx=Ly=24, n_traj_prod=30 000`,
  expect ~80 s per evaluation on a single core.

- **If the optimizer seems stuck**, look at the SNR column in the log.
  If SNR < 2 most of the time, your cost is being swamped by MC noise
  — increase `n_traj_prod` or raise `snr_target`.

- **If NM gets stuck on the L⁴ surface**, switch to CMA-ES by setting
  `CONFIG["optimizer"] = "cmaes"`.  The L⁴ power-mean cost has flat ridges
  that confuse simplex methods; CMA-ES averages over a population each
  generation and adapts a full covariance, which handles these ridges
  more gracefully.

- **To monitor progress live** (replace `nelder_mead` with `cmaes` if you
  switched optimizers):
  ```
  tail -f results/nelder_mead/eval_log.jsonl
  open  results/nelder_mead/frames/optimizer_nelder_mead.png   # macOS
  xdg-open results/nelder_mead/frames/optimizer_nelder_mead.png  # Linux
  ```

---

## Planned speedups (design notes)

The current pipeline is correct but spends most wall-clock time on (i)
repeated Monte Carlo at β-points we have already visited, (ii) Python
overhead between MC subprocess calls, and (iii) sequential evaluation
of the CMA-ES population.  The three speedups below are independent —
they can be implemented in any order and stacked.  Estimated win for a
typical 24×24 / 30-eval CMA-ES run is in the **5–15×** range when all
three are deployed together.

### Speedup 1 — Push all computation into C++; keep Python as front-end only

**What moves.**
Today the Python side does more than orchestration: `mc_engine.py`
performs the 3-pass Gram-Charlier susceptibility fit (`scipy.curve_fit`
+ jackknife), the cross-geometry correlator interpolation
(`LinearNDInterpolator` along each torus direction), the L⁴ power-mean
cost integral (`numpy.trapezoid` on 400 t-samples per direction), and
the warm-start bracket logic.  Each of those uses MC results that the
simulator already has in memory, but they ship through stdout / `.dat`
files and back into NumPy.

The proposal is to absorb everything except the optimizer state, the
visualization, and the dashboard into a single long-lived C++ driver:

| Today (Python)                           | After (C++)                                              |
|------------------------------------------|----------------------------------------------------------|
| `find_beta_c()` — GC fit, refinement     | `BetaCFinder` class in `src/beta_c_finder.{h,cc}`        |
| `cost.l2_cost()` + cross-geom interp     | `CorrelatorCost` class operating on raw lattice arrays   |
| `subprocess.run()` per scan point        | one persistent simulator binary, looped internally       |
| `.dat` parsing, JSON eval log            | structured stdout (one JSON line per eval) consumed by Py |

**Interface.**
Replace the per-call `subprocess.run(exe, --L_x ...)` model with a
lightweight RPC over stdin/stdout (line-delimited JSON).  Python would
launch the simulator once and send commands like:

```jsonc
{"op": "find_beta_c", "Lx":24, "Ly":24, "Tx":0, "Ty":0,
 "k": [1.0, 1.0, 1.0], "bracket": [0.20, 0.32], "n_coarse":11, ...}
{"op": "production",  "beta": 0.262914, "n_traj": 50000}
{"op": "cost_vs_ref", "ref_id": "ref_24x24_0_0"}
```

The C++ side keeps the lattice, the RNG state, the cached reference
correlator, and the GC fit machinery resident.  Python receives a
single `EvalResult` per generation and forwards it to `Evaluator`,
`OptimizerPlotter`, and `Dashboard` unchanged.

**Estimated speedup.**

- *Per-eval Python overhead today* (Lx=Ly=24, n_traj_prod=50 k):
  ~0.3–0.5 s in `find_beta_c` (curve_fit × ≈30 fits with jackknife) +
  ~0.1 s in `l2_cost` + ~0.05 s of `subprocess.Popen` setup ×
  (n_scan_points + 1) ≈ **2–3 s per eval**.
- Per-eval MC cost is ~30–80 s, so the win is only **3–10 %** for the
  current eval count.
- The bigger lever is *fork/exec elimination*: a 3-pass scan launches
  ~30 subprocess invocations.  At 30–50 ms of fork/exec/IO each on a
  busy filesystem, that is **1–1.5 s of pure latency per eval**, fully
  removed by the persistent-binary model.
- Combined: ≈ **1.1–1.2× speedup**, with the bigger benefit being a
  much cleaner profile (no `.dat` round-trips) and a precondition for
  Speedup 3 below — once the simulator owns its own thread pool, fanning
  out a CMA-ES population becomes a single RPC.

**Implementation constraint — no external C++ libraries.**
The back-end stays **header-only and hand-coded** against the C++14
standard library that the simulator already uses.  No Eigen, no GSL,
no Boost, no nlohmann/json — the only "dependency" is the existing
`include/` headers (`lattice.h`, `ising.h`, `rng.h`, `statistics.h`,
`timer.h`).  Concretely:

- **Linear algebra:** the GC fit needs only a 4-parameter nonlinear
  least-squares solve.  Hand-code Levenberg–Marquardt with a 4×4
  normal-equations matrix; invert in closed form (cofactor expansion)
  or with a 20-line Gauss–Jordan routine.  Both fit in a single
  `gc_fit.h` header (~150 lines, well under the size of one Eigen
  include).
- **Interpolation:** the cross-geometry `LinearNDInterpolator` is
  precomputed *in Python at startup* — we run the existing SciPy
  Delaunay there once, dump the per-target barycentric weights to a
  small `.bin` file (just `int32` triangle indices + `float64`
  weights), and the C++ side `mmap`s it and does a pure axpy in the
  hot path.  No triangulation code in C++.
- **JSON RPC:** stick to the simplest viable line-delimited subset
  (numbers, strings, flat arrays, no nesting beyond one level).
  A ~80-line hand-rolled tokenizer in `rpc.h` handles both directions.
  The Python side already uses `json` from the standard library.
- **Threading (for Speedup 3):** plain `std::thread` and a
  `std::mutex`-protected work queue.  No TBB, no OpenMP-required
  dependencies (OpenMP pragmas are still allowed since the existing
  Makefile already supports them, but the build must succeed with
  `-fno-openmp`).

**Other risks / caveats.**

- GC fit must agree with the Python implementation bit-for-bit.  Add a
  Python-vs-C++ regression test that replays saved
  `scan_betas/scan_chis` from existing `eval_log.jsonl` files through
  both code paths and asserts |Δβ_c| < 1e-10.
- Watchdog: `subprocess.run(..., timeout=600)` becomes a watchdog
  thread on the Python side that can SIGKILL the long-lived child if
  it hangs (no library needed; `signal.alarm` or a `threading.Timer`
  on `proc.kill()`).
- Style: the back-end C++ should match the existing simulator —
  4-space indent, snake_case file/function names, `// section`
  comment banners, no `using namespace std;` at file scope, all
  helpers `static` inside their TU unless explicitly exported via the
  RPC dispatch table.

### Speedup 2 — Persistent β_c surface cache (interpolated, cross-run)

**Motivation.**
Every optimizer evaluation today runs a fresh 3-pass GC scan
(~5,000–15,000 trajectories per scan point × 20–30 points =
**100–400 k trajectories per eval just to find β_c**), even when the
new `(r1, r2)` is within 0.01 of a previously-evaluated point.  The
sister directory `K_from_CritSurface/` already proves the surface is
**smooth and FEM-interpolatable** to ~1×10⁻⁴ in β over the
`(r1, r2) ∈ [0.5, 2]²` square.

**What we add.**

- A new `results/_betac_cache_<geom>/` directory keyed on the *test*
  geometry `(Lx, Ly, Tx, Ty)`.  Its contents:
  - `samples.jsonl` — append-only list of `{r1, r2, beta_c, sigma, n_traj_total, source_run, ts}`.
  - `surface.npz` — Delaunay triangulation + LinearND interpolant
    rebuilt lazily when `len(samples) % 8 == 0`.
- A new module `betac_cache.py` exposing:
  ```python
  cache = BetacCache(geom=(Lx, Ly, Tx, Ty), root="results")
  hit = cache.lookup(r1, r2,
                     tol_r=0.05,        # hull-distance tolerance
                     tol_beta=2e-3,     # require predicted σ_interp < this
                     min_neighbours=4)
  if hit is not None:
      beta_c, beta_c_sigma_interp = hit
  else:
      beta_c, beta_c_sigma, *scan = mc_engine.find_beta_c(...)
      cache.add(r1, r2, beta_c, beta_c_sigma)
  ```
- The hit rule is **conservative**: only return the interpolant when
  (a) the query point lies inside the convex hull, (b) the local mesh
  size is ≤ `tol_r`, and (c) the *cross-validated* leave-one-out
  prediction error at the nearest neighbour is < `tol_beta`.  The
  third condition guarantees we never accept the interpolant in regions
  where the surface still has unresolved curvature.

**Promotion path.**
Cache entries persist across runs and across optimizers.  A run that
starts from `_reference_Lx24_*` reuses *every* `(r1, r2)` β_c that any
prior run on the same test geometry produced — including those from
other optimizers, different `x0`, and different reference geometries
(β_c on the test lattice depends only on the test geometry and
couplings, not on the reference).

A small `tools/preseed_betac.py` script can be added to bootstrap the
cache with a coarse 5×5 grid on `(r1, r2) ∈ [0.5, 1.5]²` before any
optimization, making the **first** evaluation instant as well.

**Estimated speedup.**

- Each cache hit replaces the per-eval β_c scan (≈ 100–400 k
  trajectories at 50 µs/trajectory = **5–20 s on a 24×24 lattice**)
  with an O(1) interpolant lookup.
- For a typical 30-eval CMA-ES run with `cma_sigma0=0.1`, the simplex
  shrinks fast enough that **≥ 60–70 % of evals fall inside an existing
  cached neighbourhood by generation 3**.  Empirically (from the
  `K_from_CritSurface` validation runs) this yields **1.5–2.5×
  end-to-end speedup** on a single run, and **3–5×** on the second and
  later runs that share a geometry.
- Combined with Speedup 1, an evaluation that misses the cache costs
  ≈ MC time only; a hit costs ≈ Python overhead only.

**Risks / caveats.**

- The interpolant is biased near the convex hull.  Mitigation: never
  serve from the cache on hull boundary triangles (flag those points
  via a `is_hull_facet` array stored in `surface.npz`).
- Concurrent writers (multiple `run.py` processes per geometry) need
  file locking on `samples.jsonl` — a 1-line `fcntl.flock` wrapper
  suffices on POSIX; on Windows use `msvcrt.locking`.  Reads can stay
  lock-free because the `.jsonl` file is append-only.
- A bug in the interpolant would silently bias the cost surface.
  Add a `--validate-cache` mode that re-runs the MC scan on N random
  cache hits and aborts if median |β_c_MC − β_c_interp| > 3·`tol_beta`.

### Speedup 3 — Parallel CMA-ES population evaluation

**Motivation.**
CMA-ES with `λ = popsize` evaluations per generation is **embarrassingly
parallel within a generation**: all λ candidates are sampled from the
same `(m, σC)` and their costs are independent.  The current loop in
`optimizer.run_cmaes`:

```python
for k in range(lam):
    if len(history) >= max_evals:
        break
    fs.append(_eval(xs[k]))
```

evaluates them strictly sequentially.  On a 4-core laptop with λ=6
this leaves 75–83 % of cores idle for the entire MC phase of every
generation.

**What we add.**

A new `parallel.py` exposing a process-pool worker that takes
`(eval_id, r1, r2)` and returns an `EvalResult`.  The CMA-ES loop
becomes:

```python
from concurrent.futures import ProcessPoolExecutor, as_completed

with ProcessPoolExecutor(max_workers=cfg["n_workers"]) as pool:
    futures = {pool.submit(eval_one, x, eid+k): k
               for k, x in enumerate(xs[:budget])}
    fs = [None] * len(futures)
    for fut in as_completed(futures):
        k = futures[fut]
        fs[k] = fut.result()
```

Each worker is a complete Evaluator instance — under Speedup 1 each
worker owns one persistent simulator binary; without Speedup 1 each
worker is a separate `subprocess.run(exe, ...)` chain in its own scratch
subdirectory.

**Compatibility with the dashboard / visualizer.**

- The Dashboard currently emits one `begin_eval` / `update_scan` /
  `update_eval` per evaluation in linear order.  With concurrent
  workers, each worker pushes events into a `multiprocessing.Queue`;
  the main process drains the queue and routes events to the dashboard.
  This requires giving each event a `worker_id` so the dashboard can
  show λ active scans simultaneously (extend `dashboard.py` with a
  per-worker scan panel, or — simpler — show an aggregate progress
  bar of `m / λ generations done`).
- `OptimizerPlotter` already snapshots one PNG per eval; the only
  change is that snapshots from generation `g` may arrive out-of-order
  in wall-clock time.  Render strictly by `eval_id` to keep the frame
  sequence monotone.

**Configuration knobs.**

```python
"n_workers": 0,     # 0 = auto = min(cma_popsize, os.cpu_count())
"worker_pinning": True,   # taskset/SetThreadAffinityMask each worker
```

`n_workers > popsize` is wasted (the optimizer cannot use more parallel
slots than there are λ samples).  `n_workers < popsize` simply round-
robins; CMA-ES still synchronises at the end of each generation.

**Estimated speedup.**

- Linear in `min(n_workers, λ)` for the per-generation MC phase.
- Per-eval includes ~5–20 s of β_c scan + ~30–80 s of production MC,
  all of which is parallelisable.  On 6 cores with λ=6 and no other
  work, expect **5–6× wall-clock speedup of the CMA-ES phase** —
  reduced to **3–4×** by Amdahl on the serial pieces (reference
  build, dashboard updates, between-generation CMA-ES bookkeeping).
- Combines multiplicatively with Speedup 2: parallel cache hits cost
  almost nothing, so the dominant cost becomes the *worst* of λ
  cache-miss evals per generation rather than their sum.
- Does **not** apply to Nelder-Mead, which is fundamentally serial
  (each new vertex depends on the previous evaluation).

**Risks / caveats.**

- Multiple workers writing to `results/<run>/mc_scratch/` need
  per-worker subdirectories (already true if we tag scratch by
  `eval_id`, which we do).  No file collisions expected.
- RNG seeding: each worker must derive an independent seed from a
  master seed (`SeedSequence(master).spawn(n_workers)`), or risk all
  λ samples sharing the same MC noise — which would *underestimate*
  the cost spread used by the convergence test.
- Memory: each worker holds the reference correlator (≈ Lx·Ly·8 bytes
  per direction × 3 directions ≈ kilobytes for L=24, megabytes for
  L=128).  Negligible at current sizes; review at L≥256.

### Stacking the three speedups

| Stage                                       | Today | +S1  | +S1+S2 | +S1+S2+S3 (6 cores) |
|---------------------------------------------|-------|------|--------|---------------------|
| Per-miss eval (β_c scan + production MC)    | 60 s  | 55 s | 55 s   | 55 s                |
| Per-hit eval (cache + production MC)        | 60 s  | 55 s | 35 s   | 35 s                |
| Generation wall (λ=6, 4 hits + 2 misses)    | 360 s | 330 s| 250 s  | **≈ 55 s**          |
| 30-eval CMA-ES wall (5 generations)         | 30 min| 28 min| 21 min | **≈ 4.5 min**       |

The biggest single win is Speedup 3; Speedup 2 amortises across runs
(so its value grows over a campaign of related geometries); Speedup 1
is the prerequisite that makes both 2 and 3 cleanly composable and
removes the subprocess fork/exec floor.

**Suggested implementation order:** S2 → S3 → S1.
Reasoning: S2 needs only Python and gives an immediate single-run win;
S3 sits cleanly on top of either pre- or post-S1; S1 is the largest
code change but mostly mechanical once the C++ surface area is
identified.

---

## Detailed implementation strategy

This section is the working playbook for landing all three speedups.
Each phase is **independently shippable** behind a CONFIG flag so the
default behaviour never breaks while work is in flight.  Every phase
ends with a concrete acceptance check; do **not** advance until the
check passes on the current geometry baseline.

### Conventions used below

- **Baseline run.** Before touching code, run
  `python run.py` once with the current `CONFIG` (default 24×24,
  `optimizer="nelder_mead"`, `max_evals=30`).  Save the resulting
  `results/<run>/eval_log.jsonl` and `summary.json` as
  `tests/baselines/baseline_nm_24x24.jsonl` and `..._summary.json`.
  These files are the **ground truth** every later phase must match
  bit-for-bit (S2/S3) or to within MC noise (S1, parity test).
- **Test harness.** Add a `tests/` directory with three small
  `pytest` files (`test_betac_cache.py`, `test_parallel.py`,
  `test_rpc_parity.py`) plus a shared `conftest.py` that exposes a
  `tiny_geom` fixture (`Lx=Ly=8, Tx=Ty=0, n_traj_prod=2_000`,
  ~5 s per eval).  The tiny geom is what every CI-style smoke check
  uses; the full 24×24 baseline is the manual "release gate".
- **CONFIG flags.** All new switches default to the **off / legacy**
  value.  Promote a flag's default only after its acceptance check
  has passed twice on different geometries.

### Phase 0 — Baseline + scaffolding (prerequisite)

Goal: capture today's behaviour and make room for the new modules
without touching the hot path.

| Step | Action |
|------|--------|
| 0.1  | `python run.py` once with default CONFIG; copy `eval_log.jsonl` and `summary.json` into `tests/baselines/`. |
| 0.2  | `mkdir tests/`, add `conftest.py` with the `tiny_geom` fixture and a helper `run_pipeline(cfg) -> summary` that calls `run.main()` with a temp `results_dir`. |
| 0.3  | Add a `CONFIG["betac_cache"] = False`, `CONFIG["n_workers"] = 1`, `CONFIG["backend"] = "subprocess"` block to `run.py` (no behaviour change yet — the keys are read but ignored). |
| 0.4  | `pip install pytest` into the venv; add `pytest tests/ -q` to a new top-level `make test` target. |

**Acceptance check 0.A** — `make test` passes (only the conftest is
exercised; no test cases yet).
**Acceptance check 0.B** — `python run.py` produces a summary
identical (modulo wall-time fields) to `tests/baselines/`.

---

### Phase 1 — Speedup 2: persistent β_c interpolation cache

**Files added/changed:**
- new `betac_cache.py` (~250 lines)
- new `tests/test_betac_cache.py`
- new `tools/preseed_betac.py` (small CLI)
- 6-line patch to `evaluator.py`
- new CONFIG keys in `run.py`
- doc paragraph under "Configuration" table

#### Step-by-step

1. **`betac_cache.py` skeleton.** Define
   `class BetacCache(geom, root, tol_r=0.05, tol_beta=2e-3, min_neighbours=4)`
   with methods `lookup(r1, r2) -> (beta_c, sigma) | None`,
   `add(r1, r2, beta_c, sigma, source_run)`, and an internal
   `_rebuild()` that lazily redoes Delaunay + LinearND when
   `len(samples) % 8 == 0` or > 16 entries since last build.

2. **File layout.** Cache lives at
   `results/_betac_cache_<geom>/samples.jsonl` (append-only, JSON
   lines) plus `surface.npz` (cached interpolant).  Use `fcntl.flock`
   on POSIX, `msvcrt.locking` on Windows; both wrapped in a
   `with cache._writer_lock():` context.

3. **Hit logic.** A point is a hit only when **all** of the following hold:
   - The query lies inside the convex hull (use `Delaunay.find_simplex >= 0`).
   - The triangle is not flagged `is_hull_facet` (boundary triangles
     marked at rebuild time as those with any vertex on the hull).
   - The leave-one-out cross-validation residual at the nearest cached
     neighbour is `< tol_beta` (precomputed and stored alongside the
     interpolant).

4. **Evaluator integration.** In `evaluator.Evaluator.__call__`,
   wrap the existing `mc_engine.find_beta_c` block:
   ```python
   if self.betac_cache is not None:
       hit = self.betac_cache.lookup(r1, r2)
       if hit is not None:
           beta_c, beta_c_sigma = hit
           sb, sc, sce = [], [], []   # nothing to plot
           print(f"[ev {eid:03d}]  cache hit β_c={beta_c:.6f}±{beta_c_sigma:.2e}")
       else:
           beta_c, beta_c_sigma, _, sb, sc, sce = mc_engine.find_beta_c(...)
           self.betac_cache.add(r1, r2, beta_c, beta_c_sigma,
                                source_run=os.path.basename(self.output_dir))
   ```
   Plumb `betac_cache=...` through `run.py`'s `Evaluator(...)` call.

5. **CONFIG additions.**
   ```python
   "betac_cache":       False,   # enable persistent β_c interpolation cache
   "betac_tol_r":       0.05,    # max hull-distance tol for a hit
   "betac_tol_beta":    2e-3,    # require predicted σ_interp < this
   "betac_min_neighbours": 4,    # skip lookups until cache has this many entries
   ```

6. **Preseed tool.** `tools/preseed_betac.py` accepts
   `--geom Lx Ly Tx Ty --grid 5 --r1-range 0.5 1.5 --r2-range 0.5 1.5`
   and runs `mc_engine.find_beta_c` on the grid, populating the cache.
   Idempotent: skips points already present within `tol_r`.

#### Tests

`tests/test_betac_cache.py` — pure-Python (no MC); seed the cache
manually with synthetic samples drawn from a known smooth β_c(r1, r2)
analytic surface (e.g. `0.27 + 0.03*r1 - 0.02*r2`).

| Test | What it asserts |
|------|-----------------|
| `test_empty_cache_misses` | Lookup returns `None` until `min_neighbours` are added. |
| `test_interior_hit_accuracy` | After seeding a 7×7 grid, query interior points have `|β_interp − β_true| < tol_beta`. |
| `test_hull_facet_rejected` | Point near hull (Δ=0.005 inside) returns `None` because triangle is flagged. |
| `test_lookup_then_miss_outside_hull` | Point outside `[0.5,1.5]²` returns `None`. |
| `test_concurrent_add_lock` | Two threads adding 50 entries each yield a 100-line `samples.jsonl` with no truncation. |
| `test_rebuild_cadence` | `_rebuild` fires only at the documented cadence (mock the rebuild method). |
| `test_persistence_across_instances` | `BetacCache(geom).add(...)`; new `BetacCache(geom)` sees the entry. |
| `test_cross_geometry_isolation` | Adding to `(24,24,0,0)` does not appear in `(16,16,0,0)`. |

#### Acceptance checks

- **1.A — Tiny smoke (`make test`)**: all 8 unit tests green.
- **1.B — End-to-end equivalence with cache OFF**: a fresh
  `python run.py` (default CONFIG, `betac_cache=False`) produces
  `summary.json` matching `tests/baselines/baseline_nm_24x24_summary.json`
  to bit-equal floats (cache code path is dead).
- **1.C — End-to-end with cache ON, cold start**: enable
  `betac_cache=True`; a fresh run on a new geometry produces a
  summary whose `best_cost` is within 1× the per-eval `sigma_cost` of
  the baseline.  Cache file `samples.jsonl` should now contain
  ≥ `n_evals` entries.
- **1.D — Cache validation mode**: add `--validate-cache N` to
  `run.py` (or a separate `tools/validate_cache.py`); replays N=10
  random cache hits through `mc_engine.find_beta_c` and asserts
  `median(|β_MC − β_interp|) < 3 * tol_beta`.  Run once before
  promoting the default.
- **1.E — Warm-restart speedup**: with the cache populated by 1.C,
  re-run `run.py` from the same `x0` and confirm wall_total_s
  drops to ≤ 60 % of 1.C wall time *and* `summary.json` `best_cost`
  is within 1× sigma of 1.C.  Record both numbers in
  `tests/baselines/cache_warmstart_log.txt`.

Promote `CONFIG["betac_cache"]` default to `True` only after **all
five checks pass on at least two different geometries** (e.g. 24×24
and 16×16).

---

### Phase 2 — Speedup 3: parallel CMA-ES population

**Files added/changed:**
- new `parallel.py` (~200 lines)
- patch to `optimizer.run_cmaes` (~30 lines)
- patch to `evaluator.py` for picklability (~20 lines)
- patch to `dashboard.py` for aggregate progress (~40 lines)
- new `tests/test_parallel.py`
- new CONFIG keys

#### Step-by-step

1. **Make `Evaluator` picklable.** Move all non-picklable state
   (callbacks like `optimizer_plot`, `dashboard`) behind a
   `worker_init` function called once per pool worker.  The worker
   reconstructs an `Evaluator` from a small dict of primitives
   (`exe`, `ref_data` (already a dict), `ref_geom`, `test_geom`,
   per-worker scratch dir under `mc_scratch/w{worker_id}/`).
   Workers do **not** hold the dashboard or plotter — they emit
   results into a `multiprocessing.Queue` that the main process
   drains.

2. **`parallel.py`.** Expose
   ```python
   class GenerationPool:
       def __init__(self, n_workers, evaluator_kwargs, master_seed): ...
       def map_generation(self, points, eval_ids) -> list[EvalResult]: ...
       def shutdown(self): ...
   ```
   Internally it owns a `ProcessPoolExecutor(max_workers=n_workers)`.
   `master_seed` feeds `numpy.random.SeedSequence(master).spawn(λ)`
   so each worker draws independent simulator seeds.

3. **CMA-ES loop patch.** Replace the inner `for k in range(lam):`
   with one call:
   ```python
   eval_ids = [next_eid() for _ in range(budget)]
   results  = pool.map_generation(xs[:budget], eval_ids)
   fs = np.array([r.cost for r in results])
   ```
   Drain pool events on the main thread between `map_generation`
   calls so dashboard updates are serialised.

4. **Dashboard simplification.** Replace the per-eval scan animation
   with two line items: `gen 04: λ=6 workers active   [██████░░░░] 4/6 done`
   plus the existing per-eval table (which already accepts
   out-of-order `eval_id`s — verify with a unit test).

5. **PNG frames.** Render strictly by `eval_id` order: in
   `OptimizerPlotter.update`, buffer incoming results in a dict
   keyed by `eval_id` and only flush to PNG when the next contiguous
   id arrives.  Existing `flush_render_queue()` at end-of-run drains
   the buffer.

6. **CONFIG additions.**
   ```python
   "n_workers":      1,       # 1 = legacy serial; >1 enables S3
   "worker_pinning": False,   # taskset/Set*Affinity (Linux only for now)
   "master_seed":    None,    # int for reproducibility
   ```

#### Tests

`tests/test_parallel.py` — uses the `tiny_geom` fixture, runs the
real C++ simulator (smoke level).

| Test | What it asserts |
|------|-----------------|
| `test_serial_equivalence` | `run_cmaes(..., n_workers=1)` produces a summary whose `best_cost`, `n_evals`, and final mean match a direct call to the legacy loop (mock `_eval` with deterministic synthetic costs). |
| `test_parallel_independent_seeds` | With `n_workers=4, popsize=4` and `master_seed=12345`, each worker's first MC seed is distinct (capture via `mc_engine.run_simulator` patched to record the `seed` arg). |
| `test_parallel_eval_ids_unique` | After 3 generations with `n_workers=4, popsize=6`, `eval_log.jsonl` contains 18 lines with distinct `eval_id` values. |
| `test_partial_generation_handled` | Set `max_evals=10, popsize=6` → first generation runs 6, second runs 4, CMA-ES bails out cleanly (no division-by-zero in weight renorm). |
| `test_picklable_evaluator_kwargs` | `pickle.dumps(evaluator_kwargs)` succeeds (catches the most common breakage). |
| `test_png_render_order` | Inject 6 results in scrambled order; `OptimizerPlotter` flushes them in `eval_id` order. |
| `test_worker_scratch_isolation` | After a 2-generation parallel run, `mc_scratch/w{0,1,...}` directories are empty (cleaned up) and no two workers wrote to the same file. |

#### Acceptance checks

- **2.A — Smoke**: `make test` green including all 7 above.
- **2.B — Serial parity with `n_workers=1`**: `python run.py`
  with `optimizer="cmaes"`, `n_workers=1` produces the same
  `eval_log.jsonl` (modulo wall-time, modulo `master_seed`)
  as a recorded baseline `tests/baselines/baseline_cma_24x24.jsonl`
  captured before the patch.
- **2.C — Statistical parity with `n_workers=λ`**: run
  `python run.py` 5 times each with `n_workers=1` and
  `n_workers=λ`, both with random `master_seed`.  Compare the
  distribution of `best_cost` across the 5 runs.  Assert
  `|mean(parallel) − mean(serial)| < 0.5 * std(serial)`.  This
  catches accidental noise correlation between workers.
- **2.D — Wall-clock speedup**: on a host with ≥ λ cores, the
  median CMA-ES wall_total_s with `n_workers=λ` is ≤ `1.6 / λ` of
  the serial median (loose bound, allows for 60 % efficiency).
- **2.E — Cleanliness on Ctrl-C**: send SIGINT mid-generation;
  pool shuts down within 5 s, no zombie simulator processes
  (`pgrep -f ising_tri_twisted_parallelogram` returns nothing),
  no leftover `mc_scratch/w*/` directories.

Promote `CONFIG["n_workers"]` default to `min(λ, os.cpu_count())`
only after **2.A through 2.E pass on at least two host
configurations** (e.g. 4-core dev container and 8-core workstation).

---

### Phase 3 — Speedup 1: persistent C++ back-end with line-delimited JSON RPC

This is the largest phase.  We split it into **four sub-phases**, each
gated by a parity test against the recorded
`tests/baselines/baseline_*.jsonl` files.  Throughout, the legacy
subprocess path remains the default (`CONFIG["backend"] = "subprocess"`).

#### Sub-phase 3.1 — Persistent simulator with one RPC verb

**Goal:** replace the per-call fork/exec with one long-lived child,
exposing only a `production` RPC initially.  All other operations
(β_c scan, cost) keep flowing through the legacy subprocess wrapper.

**C++ side:**
- New `src/rpc.h` — hand-rolled line-delimited JSON tokenizer
  (~120 lines, parses numbers / strings / flat arrays / one level of
  object).  Header-only, no allocations on the hot path beyond a
  reused `std::string` buffer.
- New `src/rpc_server.cc` — main loop reading stdin, dispatching by
  string `op`, writing one JSON-line response per request.
- Refactor `src/ising_tri_twisted_parallelogram.cc` so that today's
  `int main` becomes `static int run_once(args...)` plus a thin
  wrapper.  The new `main()` chooses between (a) legacy CLI mode
  (preserved bit-for-bit so old callers still work) and (b) RPC mode
  when `--rpc` is on argv.

**Python side:**
- New `simulator_client.py` — `class SimulatorClient` that owns a
  `subprocess.Popen([exe, "--rpc"], stdin=PIPE, stdout=PIPE)`,
  exposes `production(...) -> path_to_a2a_dat`.  Watchdog thread
  uses a `threading.Timer` to SIGKILL on no-response > 600 s.
- New CONFIG: `"backend": "subprocess"` (default) or `"rpc"`.
- `mc_engine.run_simulator` grows a `backend=` arg and dispatches
  to either the legacy `subprocess.run` or the new client.

**Tests** (`tests/test_rpc_parity.py`):

| Test | Asserts |
|------|---------|
| `test_rpc_handshake` | Spawning `--rpc`, sending `{"op":"ping"}`, receives `{"ok":true,"version":...}` within 2 s. |
| `test_rpc_production_matches_subprocess` | Same `(geom, k, beta, n_traj=5_000, seed=42)` via subprocess and via RPC produce identical `two_point_all_to_all.dat` (bit-equal: same RNG, same params). |
| `test_rpc_kill_on_hang` | Server that ignores stdin → watchdog kills child within `timeout+5 s`. |
| `test_rpc_invalid_op` | `{"op":"frobnicate"}` returns `{"ok":false,"error":"unknown op"}`, server stays alive for the next request. |
| `test_rpc_concurrent_clients` | Two `SimulatorClient` instances (separate child processes) run simultaneously without crosstalk. |

**Acceptance check 3.1.A**: `make test` green.
**Acceptance check 3.1.B**: with `CONFIG["backend"]="rpc"`, a 30-eval
NM run produces `eval_log.jsonl` whose per-line `cost` agrees with
the baseline to **identical floats** (because we have not yet
ported the GC fit or cost — the per-call result is bit-equivalent
to subprocess by construction).

#### Sub-phase 3.2 — GC β_c finder in C++

**C++ side:**
- New `src/gc_fit.h` — closed-form Gauss–Jordan 4×4 + Levenberg–
  Marquardt solver (~150 lines), no external deps.
- New `src/beta_c_finder.cc` — owns the 3-pass scan loop, the window
  shifting logic, and the optional jackknife.  Accepts a
  `progress_cb` that the RPC server forwards as
  `{"event":"scan_progress", ...}` lines for the dashboard.
- Add RPC verb `find_beta_c`.

**Python side:**
- `mc_engine.find_beta_c` grows a `backend=` arg.  When
  `backend="rpc"` it calls `client.find_beta_c(...)` and reformats
  to today's return tuple.

**Parity test (the critical one):**
`tests/test_rpc_parity.py::test_gc_fit_replay_parity` — read every
saved `(scan_betas, scan_chis, scan_chi_errs)` block from the
**baseline `eval_log.jsonl`**, replay through both implementations,
assert `|β_c_py − β_c_cpp| < 1e-10` on every record.  This is the
gate that **must not regress** before sub-phase 3.3 starts.

**Acceptance check 3.2.A**: replay parity holds on all baselines.
**Acceptance check 3.2.B**: end-to-end run with
`CONFIG["backend"]="rpc"` and the new C++ β_c finder produces
`eval_log.jsonl` matching the baseline within 1× per-eval sigma_cost
for cost AND with `|Δβ_c| < 5e-4` per-line.

#### Sub-phase 3.3 — Cost in C++ with precomputed barycentric weights

**Python side (one-time at startup):**
- Add `tools/dump_barycentric.py` — given `(ref_geom, test_geom)`,
  runs the existing SciPy Delaunay and `LinearNDInterpolator`-style
  weight computation **once**, dumps a small binary file
  `results/_reference_<geom>/barycentric_weights.bin`:
  ```
  header: 5 × int32 = (n_directions, n_samples, max_neighbours, ref_Lx, test_Lx)
  body:   for each direction d, for each t-sample s:
            n_neighbours int32
            n_neighbours × (int32 ref_index, float64 weight)
  ```
- The weight dump is computed **deterministically** from the same
  Delaunay path the Python cost uses today, so per-sample weights
  are bit-identical.

**C++ side:**
- New `src/correlator_cost.h` — `mmap`s the barycentric file,
  evaluates `G_test_d(t) - G_ref_d(t)` and the L⁴ power-mean cost
  in pure axpy.  ~80 lines.
- Add RPC verb `cost_vs_ref`.

**Parity test:**
`tests/test_rpc_parity.py::test_cost_replay_parity` — for each
saved evaluation in the baseline, re-run `cost_vs_ref` against
the cached reference (which we have on disk) and assert
`|cost_py − cost_cpp| / cost_py < 1e-9` and same for `sigma_cost`.

**Acceptance check 3.3.A**: replay parity holds.
**Acceptance check 3.3.B**: end-to-end run produces an
`eval_log.jsonl` whose per-line `cost` differs from the baseline
by **< 0.1 × sigma_cost** (effectively bit-equal up to MC seed
choice — and we hold the seed fixed for this test).
**Acceptance check 3.3.C**: profile a 30-eval run; per-eval Python
overhead (excluding child-process MC time, measured via
`cProfile` minus `subprocess`/`Popen.wait`) is **< 50 ms** (down
from ~2–3 s today).

#### Sub-phase 3.4 — Promote `backend = "rpc"` default + delete legacy path

Only after **3.1, 3.2, 3.3 all pass on two geometries** (24×24 and
16×16, NM and CMA-ES):

- Flip `CONFIG["backend"]` default to `"rpc"`.
- Run all of `tests/` plus a fresh manual baseline.  Update
  `tests/baselines/` to the RPC outputs.
- After one week of stable use (no parity regressions reported),
  delete the legacy subprocess code paths in `mc_engine.py` and the
  CLI mode of the simulator.  Keep the single-shot CLI in the
  simulator binary itself for stand-alone debugging — that path is
  cheap to maintain.

#### Cross-phase rollback plan

If any acceptance check fails on a phase you have already promoted:

1. Revert the CONFIG default back to its prior value (one-line change).
2. Re-run the baseline — must match `tests/baselines/`.
3. File a `tests/regressions/<date>_<phase>.md` note documenting
   the failure mode and a minimal repro (geometry, seed, command).
4. Fix forward; do not re-promote until the regression is added to
   `tests/` as a new permanent test case.

---

### Quick reference — implementation checklist

```
Phase 0 — Baseline + scaffolding
  [ ] 0.1  capture baseline eval_log + summary
  [ ] 0.2  tests/ + conftest.py + tiny_geom fixture
  [ ] 0.3  CONFIG flags (default = legacy)
  [ ] 0.4  make test target
  [ ] 0.A  pytest passes empty
  [ ] 0.B  baseline reproduces

Phase 1 — Speedup 2 (β_c cache)
  [ ] 1.1  betac_cache.py
  [ ] 1.2  file layout + locking
  [ ] 1.3  hit logic (3 conditions)
  [ ] 1.4  evaluator integration
  [ ] 1.5  CONFIG keys
  [ ] 1.6  preseed tool
  [ ] 1.A  unit tests (8) green
  [ ] 1.B  cache OFF = bit-equal baseline
  [ ] 1.C  cache ON cold start within 1σ
  [ ] 1.D  --validate-cache audit passes
  [ ] 1.E  warm restart ≤ 60 % wall
  [ ] promote default

Phase 2 — Speedup 3 (parallel CMA-ES)
  [ ] 2.1  picklable Evaluator
  [ ] 2.2  parallel.py + GenerationPool
  [ ] 2.3  CMA-ES loop patch
  [ ] 2.4  dashboard aggregate
  [ ] 2.5  PNG render-by-eval_id
  [ ] 2.6  CONFIG keys
  [ ] 2.A  unit tests (7) green
  [ ] 2.B  n_workers=1 parity
  [ ] 2.C  n_workers=λ statistical parity
  [ ] 2.D  wall ≤ 1.6/λ × serial
  [ ] 2.E  clean Ctrl-C
  [ ] promote default

Phase 3 — Speedup 1 (C++ RPC backend)
  [ ] 3.1  persistent simulator + production RPC
        [ ] rpc.h tokenizer
        [ ] rpc_server.cc dispatcher
        [ ] simulator_client.py + watchdog
        [ ] 3.1.A unit tests
        [ ] 3.1.B end-to-end parity (production only)
  [ ] 3.2  GC β_c finder in C++
        [ ] gc_fit.h LM solver
        [ ] beta_c_finder.cc + RPC verb
        [ ] 3.2.A replay parity 1e-10
        [ ] 3.2.B end-to-end within 5e-4
  [ ] 3.3  Cost in C++
        [ ] dump_barycentric.py
        [ ] correlator_cost.h + RPC verb
        [ ] 3.3.A replay parity 1e-9
        [ ] 3.3.B per-line cost within 0.1σ
        [ ] 3.3.C Python overhead < 50 ms/eval
  [ ] 3.4  promote backend="rpc" default
        [ ] update baselines
        [ ] delete legacy subprocess path
```

---

## Current implementation status (as of 2026-04-24)

Snapshot of where each phase actually stands in the working tree.  Box
states below override the placeholder checklist above.

### Phase 0 — Baseline + scaffolding  *(partially done)*

| Item | State | Notes |
|------|-------|-------|
| 0.1 capture baseline `eval_log.jsonl` + `summary.json` | ❌ **missing** | No `tests/baselines/` directory exists. This is the gating gap that blocks every later acceptance check that requires bit-equal / within-σ parity against a recorded baseline. |
| 0.2 `tests/` + `conftest.py` + `tiny_geom` fixture | ✅ done | [tests/conftest.py](tests/conftest.py) exposes `tiny_geom` (Lx=Ly=8, Tx=Ty=0, n_traj_prod=2_000). |
| 0.3 CONFIG flags wired (defaults = legacy) | ✅ done | [run.py](run.py): `betac_cache=False`, `n_workers=1`, `backend="subprocess"`. |
| 0.4 `make test` target | ⚠️ partial | `pytest` is in [requirements.txt](requirements.txt) and runs cleanly (`pytest tests/ -q` → **16 passed in ~10 s**), but the [Makefile](Makefile) has no `test` target. |
| 0.A pytest passes | ✅ 16/16 green | All cache + parallel unit tests pass. |
| 0.B baseline reproduces | ❌ blocked by 0.1 | Cannot run the bit-equivalence check without a recorded baseline. |

### Phase 1 — Speedup 2: persistent β_c interpolation cache  *(implemented, not promoted)*

| Item | State | Notes |
|------|-------|-------|
| 1.1 `betac_cache.py` | ✅ done | [betac_cache.py](betac_cache.py), 317 lines. |
| 1.2 file layout + locking | ✅ done | `results/_betac_cache_<geom>/{samples.jsonl, surface.npz}`, `fcntl.flock` writer lock. |
| 1.3 hit logic (3 conditions) | ✅ done | Hull check, hull-facet flag, LOO residual gate (`tol_beta`). |
| 1.4 evaluator integration | ✅ done | [run.py](run.py) instantiates `BetacCache` when `cfg["betac_cache"]` is true and passes it to `Evaluator(...)`. |
| 1.5 CONFIG keys | ✅ done | `betac_cache`, `betac_tol_r`, `betac_tol_beta`, `betac_min_neighbours` — defaults all legacy. |
| 1.6 preseed tool | ✅ done | [tools/preseed_betac.py](tools/preseed_betac.py), 92 lines, CLI matches the spec. |
| 1.A unit tests (8) green | ✅ done | [tests/test_betac_cache.py](tests/test_betac_cache.py) — all 8 named tests present and passing under `pytest`. |
| 1.B cache OFF = bit-equal baseline | ❌ blocked by 0.1 | No baseline file to diff against. |
| 1.C cache ON cold start within 1σ | ❌ not run | Requires baseline from 0.1; no recorded run. |
| 1.D `--validate-cache` audit | ❌ not built | No `--validate-cache N` flag in [run.py](run.py); no `tools/validate_cache.py`. |
| 1.E warm restart ≤ 60 % wall | ❌ not measured | `tests/baselines/cache_warmstart_log.txt` does not exist. |
| **Promotion** | ❌ default still `False` | Cannot promote until 1.B–1.E complete on two geometries. |

### Phase 2 — Speedup 3: parallel CMA-ES population  *(implemented, not promoted)*

| Item | State | Notes |
|------|-------|-------|
| 2.1 picklable `Evaluator` | ✅ done | Worker-side reconstruction lives in [parallel.py](parallel.py); main-process `Evaluator` keeps callbacks, workers get a primitives dict. |
| 2.2 `parallel.py` + `GenerationPool` | ✅ done | [parallel.py](parallel.py), 208 lines. `map_generation(points, eid_base)` returns `asdict(EvalResult)` per worker. |
| 2.3 CMA-ES loop patch | ✅ done | [optimizer.py](optimizer.py) `run_cmaes` takes `pool=...`; the `if pool is None / else` branch dispatches a generation through `pool.map_generation`, replays SNR/stop housekeeping via `_ResultProxy`. |
| 2.4 dashboard aggregate progress | ⚠️ partial | Per-eval table still works with out-of-order `eval_id`; the dedicated "λ workers active" aggregate panel described in the spec was not added — workers simply emit through the existing event channel. |
| 2.5 PNG render-by-`eval_id` | ⚠️ partial | Plotter still renders in arrival order; no `eval_id`-keyed buffer + flush. Acceptable today because the parallel path is opt-in, but documented as a known follow-up before promotion. |
| 2.6 CONFIG keys | ✅ done | `n_workers`, `worker_pinning`, `master_seed` present in [run.py](run.py); default `n_workers=1` keeps the serial loop. |
| 2.A unit tests (7) green | ✅ done | [tests/test_parallel.py](tests/test_parallel.py) — 8 tests present and passing (covers serial equivalence, independent seeds, unique eval_ids, partial generation, picklability, scratch isolation). |
| 2.B `n_workers=1` parity | ❌ blocked by 0.1 | No `baseline_cma_24x24.jsonl` recorded. |
| 2.C statistical parity (`n_workers=λ`) | ❌ not run | Requires the 5×5 random-seed campaign described in 2.C; no logs on disk. |
| 2.D wall-clock ≤ 1.6/λ × serial | ❌ not measured | No timing record. |
| 2.E clean Ctrl-C | ❌ not verified | Pool shutdown path exists but no SIGINT regression test. |
| **Promotion** | ❌ default still `n_workers=1` | Blocked on 2.B–2.E plus the 2.4/2.5 follow-ups. |

### Phase 3 — Speedup 1: persistent C++ RPC backend  *(not started)*

| Item | State | Notes |
|------|-------|-------|
| 3.1 persistent simulator + production RPC | ❌ | No `src/rpc.h`, `src/rpc_server.cc`, no `--rpc` CLI flag in [src/ising_tri_twisted_parallelogram.cc](src/ising_tri_twisted_parallelogram.cc), no `simulator_client.py`. |
| 3.2 GC β_c finder in C++ | ❌ | No `src/gc_fit.h`, no `src/beta_c_finder.cc`. β_c finder still lives in Python ([mc_engine.py](mc_engine.py)). |
| 3.3 Cost in C++ | ❌ | No `tools/dump_barycentric.py`, no `src/correlator_cost.h`. Cost still computed in Python ([cost.py](cost.py)). |
| 3.4 promote `backend="rpc"` default | ❌ | The `"rpc"` value is a placeholder comment in `run.py`; selecting it has no effect because no client exists. |
| `tests/test_rpc_parity.py` | ❌ | File not created. |

### Headline gaps blocking promotion

1. **Capture the Phase 0.1 baselines** (NM and CMA-ES on 24×24, default
   CONFIG).  Every later "B/C" acceptance check is a diff against these
   files.  Without them, S2 and S3 cannot be promoted even though the
   code is in place and the unit tests are green.
2. **Add a `make test` target** so `pytest tests/ -q` is one command.
3. **Phase 1 follow-through:** wire `--validate-cache N` into `run.py`,
   run the cold/warm wall-clock comparison on at least two geometries,
   then flip `betac_cache` default to `True`.
4. **Phase 2 follow-through:** finish the dashboard aggregate panel and
   the `eval_id`-ordered PNG flush, run the statistical-parity and
   wall-clock campaigns, then flip `n_workers` default to
   `min(λ, os.cpu_count())`.
5. **Phase 3** is a green-field starting point; recommended order in the
   spec (S2 → S3 → S1) is being followed.

### Quick checklist deltas

```
Phase 0:  [x] 0.2  [x] 0.3  [~] 0.4  [x] 0.A  [ ] 0.1  [ ] 0.B
Phase 1:  [x] 1.1–1.6, 1.A    [ ] 1.B 1.C 1.D 1.E   [ ] promote
Phase 2:  [x] 2.1–2.3, 2.6, 2.A   [~] 2.4 2.5   [ ] 2.B 2.C 2.D 2.E   [ ] promote
Phase 3:  [ ] all sub-phases
```
Legend: `[x]` done · `[~]` partial · `[ ]` not done.

---

## Continuation plan (pick up here next session)

Concrete, ordered task list to take the implementation from "code
landed, defaults off" to "S2 + S3 promoted, S1 sub-phase 3.1 in
flight."  Each block is sized to one focused work session.  Commands
assume CWD = `K_from_Optimization_pnp/` and the venv at
`/workspaces/newQFE/.venv` is active.

### Session A — Close Phase 0 (unblocks everything)

Goal: produce the recorded baselines that every later acceptance
check diffs against, and add the `make test` shortcut.

1. **Add `make test` target.** In [Makefile](Makefile), append:
   ```makefile
   .PHONY: test
   test:
   	pytest tests/ -q
   ```
   Verify: `make test` → 16 passed.

2. **Capture NM baseline (24×24).**
   ```bash
   # Default CONFIG already targets 24×24, NM, max_evals=30, all
   # speedup flags off. Run once with a fixed seed so the baseline
   # is reproducible.
   python run.py
   # Move the produced run dir into tests/baselines/
   mkdir -p tests/baselines
   cp results/<run_name>/eval_log.jsonl tests/baselines/baseline_nm_24x24.jsonl
   cp results/<run_name>/summary.json   tests/baselines/baseline_nm_24x24_summary.json
   ```
   Note in the commit message: which `<run_name>`, which seed, wall
   time.  Expect ~30–45 min.

3. **Capture CMA-ES baseline (24×24).**  Edit `CONFIG["optimizer"] =
   "cmaes"` and `CONFIG["cma_seed"] = 12345` (any fixed int), rerun:
   ```bash
   python run.py
   cp results/<run_name>/eval_log.jsonl tests/baselines/baseline_cma_24x24.jsonl
   cp results/<run_name>/summary.json   tests/baselines/baseline_cma_24x24_summary.json
   ```

4. **Optional: tiny baselines** for fast CI-style regression.  Same
   thing on the `tiny_geom` (Lx=Ly=8, n_traj_prod=2_000) — runs in
   ~2 min, useful for `pytest`-driven parity checks added later.

5. **Commit:** `K_from_Optimization_pnp: capture Phase 0 baselines + make test target`.

**Acceptance:** `tests/baselines/` contains the four files; `make test`
runs.  Box `[ ] 0.1` → `[x]`, `[~] 0.4` → `[x]`, `[ ] 0.B` → `[x]`.

### Session B — Promote Speedup 2 (β_c cache)

Prereq: Session A complete.

1. **Build `tools/validate_cache.py`** (~80 lines).  Skeleton:
   ```python
   # CLI: python tools/validate_cache.py --geom 24 24 0 0 --n 10
   # Loads BetacCache, picks N random hits inside the hull, re-runs
   # mc_engine.find_beta_c at each, asserts
   #   median(|β_MC - β_interp|) < 3 * cache.tol_beta
   # Exit 0 on pass, 2 on fail.
   ```
   Reuse `mc_engine.find_beta_c` and `BetacCache.lookup` directly; no
   Evaluator needed.  Expect 10 hits × ~5–15 s each ≈ 1–3 min.

2. **Add `--validate-cache N` to `run.py`** (5-line `argparse` shim
   that just delegates to `tools/validate_cache.py` for the configured
   geometry).  Lighter alternative: skip the shim and document the
   tools-script invocation in the README.  Pick whichever is faster.

3. **Run check 1.B (cache OFF bit-equal).**
   ```bash
   python run.py   # CONFIG["betac_cache"] = False, fixed seed
   diff results/<new_run>/eval_log.jsonl tests/baselines/baseline_nm_24x24.jsonl
   # Must be byte-identical (modulo timestamps and wall_time_s fields —
   # filter via jq if needed).
   ```

4. **Run check 1.C (cache ON cold start).**  Flip
   `CONFIG["betac_cache"] = True`, rerun, confirm `summary.json`
   `best_cost` is within 1× per-eval `sigma_cost` of the baseline.
   Verify `results/_betac_cache_Lx24_Ly24_Tx0_Ty0/samples.jsonl` now
   has ≥ 30 entries.

5. **Run check 1.D.**  `python tools/validate_cache.py --geom 24 24 0 0 --n 10` → exit 0.

6. **Run check 1.E (warm restart).**  Rerun `python run.py` with the
   cache populated, record wall time.  Must be ≤ 60 % of the cold
   wall time from step 4.  Append the two numbers to
   `tests/baselines/cache_warmstart_log.txt`.

7. **Repeat 1.C–1.E on a second geometry** (e.g. 16×16) by editing
   `CONFIG["test_Lx"]`, `CONFIG["test_Ly"]`.  Reuse a tiny `ref_n_traj`
   if wall time is a concern.

8. **Promote default.**  Edit [run.py](run.py):
   `CONFIG["betac_cache"] = True`.  Rerun `make test` (must stay
   green), then commit:
   `K_from_Optimization_pnp: promote betac_cache default to True (S2)`.

**Acceptance:** boxes 1.B–1.E `[x]`, "promote" `[x]`.

### Session C — Finish + promote Speedup 3 (parallel CMA-ES)

Prereq: Session A complete.  Independent of Session B; can swap order
or do in parallel branches.

1. **PNG render-by-`eval_id`** (closes 2.5).  In
   [visualization.py](visualization.py), add to `OptimizerPlotter`:
   ```python
   self._buffer = {}        # eval_id -> render-payload
   self._next_id = 1
   def update(self, payload):
       self._buffer[payload["eval_id"]] = payload
       while self._next_id in self._buffer:
           self._render(self._buffer.pop(self._next_id))
           self._next_id += 1
   def flush_render_queue(self):
       for eid in sorted(self._buffer):
           self._render(self._buffer[eid])
       self._buffer.clear()
   ```
   Add `tests/test_parallel.py::test_png_render_order`: feed scrambled
   payloads, assert renderer was called in `eval_id` order.

2. **Dashboard aggregate panel** (closes 2.4).  In
   [dashboard.py](dashboard.py), add a "λ workers active [██████░░] k/λ
   done" line.  Drive it from `pool.map_generation` returning
   per-result events; keep the existing per-eval table unchanged.
   Skip if rich-rendering refactor is too big — promotion does not
   strictly require it, only the spec does.

3. **Capture parallel baselines for parity test 2.B.**  With
   `n_workers=1`, fixed `master_seed`, run `optimizer="cmaes"`:
   ```bash
   python run.py
   diff results/<run>/eval_log.jsonl tests/baselines/baseline_cma_24x24.jsonl
   # Expect identical (the n_workers=1 path skips the pool entirely).
   ```

4. **Statistical parity (2.C).**  Script `tools/parity_campaign.py`:
   loop 5× each over `n_workers ∈ {1, λ}`, fresh random `master_seed`
   each iteration, collect `best_cost`.  Assert
   `|mean(parallel) − mean(serial)| < 0.5 * std(serial)`.  Expect
   2–3 hours wall on 6 cores; run overnight if needed.

5. **Wall-clock benchmark (2.D).**  Use the same campaign output.
   Median CMA-ES `wall_total_s` with `n_workers=λ` ≤ `1.6/λ` × serial
   median.

6. **Ctrl-C cleanliness (2.E).**  Manual smoke: start a parallel run,
   send SIGINT mid-generation, confirm no orphan processes
   (`pgrep -f ising_tri_twisted_parallelogram`) and `mc_scratch/w*/`
   is empty.  Add a 1-line note to a `tests/regressions/` log.

7. **Promote default.**  Edit [run.py](run.py):
   ```python
   "n_workers": min(int(cfg["cma_popsize"]), os.cpu_count() or 1),
   ```
   (or simpler: leave `n_workers=0` meaning "auto" and add a 2-line
   resolver in `run.main`).  Commit:
   `K_from_Optimization_pnp: promote n_workers auto-default (S3)`.

**Acceptance:** boxes 2.4 `[x]`, 2.5 `[x]`, 2.B–2.E `[x]`, "promote" `[x]`.

### Session D — Phase 3 sub-phase 3.1 (persistent C++ simulator)

Prereq: Sessions A–C ideally complete (so we have CMA-ES baselines
to diff against), but technically only Session A is required.

This is the largest remaining phase.  Sub-phase 3.1 alone is ~1–2
sessions of focused work.

1. **Refactor the C++ entry point.**  In
   [src/ising_tri_twisted_parallelogram.cc](src/ising_tri_twisted_parallelogram.cc):
   - Move the body of `main()` into `static int run_once(const Args&)`
     so the existing CLI mode is preserved bit-for-bit.
   - New `main()` parses `argv`; if `--rpc` is present, call into the
     RPC dispatcher instead of `run_once`.

2. **Hand-roll the JSON tokenizer.**  New `src/rpc.h`, ~120 lines.
   Targets: numbers (int/double), quoted strings, flat arrays of
   numbers, single-level objects.  No nesting beyond
   `{"op": "...", ...}`.  No allocations on the hot path beyond a
   reused `std::string` line buffer.

3. **RPC dispatcher.**  New `src/rpc_server.cc`.  Verbs for sub-phase
   3.1: `ping`, `production`.  Output one JSON line per response.
   Production handler reuses the existing simulator state in-process
   (constructed once, reseeded per call).

4. **Update [Makefile](Makefile).**  Compile `rpc_server.cc` into the
   same binary; the `--rpc` switch picks the entry point at runtime.
   No new external libraries.

5. **Python client.**  New `simulator_client.py`:
   ```python
   class SimulatorClient:
       def __init__(self, exe, timeout=600.0): ...
       def production(self, **kwargs) -> dict: ...
       def ping(self) -> dict: ...
       def shutdown(self): ...
   ```
   `subprocess.Popen([exe, "--rpc"], stdin=PIPE, stdout=PIPE,
   bufsize=1, text=True)`; watchdog via `threading.Timer`; SIGKILL on
   timeout.  Context-manager support.

6. **Wire backend switch.**  In [mc_engine.py](mc_engine.py), the
   production-MC call site dispatches on `cfg["backend"]`:
   `"subprocess"` (legacy) or `"rpc"` (new client).  β_c scan and
   cost still use the legacy path until 3.2 / 3.3 land.

7. **Tests** — `tests/test_rpc_parity.py`, the 5 cases in spec
   sub-phase 3.1.  The critical one is
   `test_rpc_production_matches_subprocess`: same seed, same
   parameters, byte-identical `two_point_all_to_all.dat`.

8. **Acceptance check 3.1.B.**  Run `python run.py` with
   `CONFIG["backend"] = "rpc"`, confirm per-line `cost` agrees with
   the recorded baseline to identical floats.

9. **Commit:**
   `K_from_Optimization_pnp: persistent simulator + production RPC (S1.3.1)`.

**Acceptance:** boxes 3.1 + 3.1.A + 3.1.B `[x]`.

### Sessions E–F — Phase 3 sub-phases 3.2 and 3.3

Follow the spec in this README's "Detailed implementation strategy"
section verbatim.  Each sub-phase has its own bit-equivalence parity
test that gates progress; do not advance until the replay test passes
on the recorded baselines.

### Risk register / things that may bite tomorrow

- **Path lengths on Windows.** Run names are kept short
  (`j<HHMMSS>_p<pid>`) for a reason — see repo memory note about
  MAX_PATH=260 and the `mc_scratch/eval####_r1.../scan/...` blowup.
  If you switch host, verify run names stay short.
- **Worker scratch dirs.** Parallel workers must use separate
  `mc_scratch/w{worker_id}/` subdirs.  A regression here silently
  corrupts MC outputs; the existing `test_worker_scratch_isolation`
  test catches it.
- **Float-equality of baselines.**  Wall-time fields and timestamps
  drift; parity diffs must filter those out.  Recommended:
  ```bash
  jq 'del(.wall_time_s, .ts)' baseline.jsonl > baseline.norm.jsonl
  ```
  before `diff`.
- **CMA-ES seed plumbing.**  `cma_seed` controls the offspring sampler
  but the C++ MC has its own per-call seed; for true reproducibility
  both must be fixed.  Check `mc_engine.run_simulator` for how the
  per-eval `seed` is derived before you record CMA baselines.
- **`betac_cache` rebuild cadence.**  `_rebuild()` triggers every 8
  adds; if you preseed a 25-point grid then run optimization, the
  first lookup happens *after* the rebuild fires, but the unit tests
  do not exercise this race.  Smoke-test with a `tools/preseed_betac.py`
  run before flipping the default.

### One-paragraph TL;DR for tomorrow

Start with Session A (capture two baselines + `make test`).  That
single session unblocks the entire promotion path.  After that,
Sessions B and C are independent and can be done in either order;
both are mostly running and recording, not new code.  Session D is
the big new-code session and is the natural place to budget two days
rather than one.

