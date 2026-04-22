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
- `results/nelder_mead/eval_log.jsonl` — one JSON line per evaluation
- `results/nelder_mead/summary.json`   — final summary including the chosen (r1, r2)
- `results/nelder_mead/frames/`        — live PNG visualisation, updated every eval

---

## What it does

Three phases:

1. **Reference build** (one-time per `ref_L`): runs a 3-pass Gram-Charlier
   susceptibility-peak scan to find β_c on an `Lref × Lref` equilateral
   lattice (k1=k2=k3=1), then a long production MC at β_c to measure the
   all-to-all two-point correlator G_ref.  Cached to disk.

2. **Per-evaluation Monte Carlo** (called by the optimizer with a candidate
   `(r1, r2)`):
   - Run a 3-pass β_c scan on the test lattice with couplings `(r1, r2, 1)`
     (warm-started from the previous evaluation's β_c).
   - Run a production MC at the new β_c to measure G_test.
   - Compute the cost.

3. **Cost function** (`cost.py`):
   $$ C(r_1, r_2) = \sum_{d \in \{v, u, w\}}
        \int_0^1 \bigl[G_\text{test}^{(d)}(t) - G_\text{ref}^{(d)}(t)\bigr]^2 \, dt $$
   summed over the three torus boundary directions, integrated by trapezoid
   rule on `n_samples = 400` points.  Plain L²; no anisotropy penalty, no
   variance normalisation.  Propagated 1-σ uncertainty σ_C is also computed
   from the per-point MC errors so the optimizer can monitor the SNR.

4. **Optimization** (`optimizer.py`):
   Hand-coded Nelder-Mead simplex in `(r1, r2)` plane.  Stops on:
   - max_evals reached
   - simplex diameter < `nm_xatol` AND cost spread < `nm_fatol`
   - running best SNR < `indist_stop_snr` (test indistinguishable from ref)

   Adaptive statistics: when SNR < `snr_target`, `n_traj_prod` is multiplied
   by 1.5 (capped at `snr_max_traj_factor × starting value`).

---

## Configuration

All knobs are at the top of `run.py` in the `CONFIG` dictionary:

| Key | Meaning | Typical |
|-----|---------|---------|
| `ref_L`                        | Reference lattice side `L`                                | 8 / 16 / 24 |
| `test_Lx, test_Ly`             | Test lattice dimensions                                   | same as `ref_L` |
| `test_Tx, test_Ty`             | Twist parameters (0 for untwisted square)                 | 0, 0 |
| `ref_n_traj`                   | Production trajectories at the reference β_c              | 60 000 |
| `ref_scan_n_traj_coarse/fine`  | β_c scan trajectories per point (reference)               | 4 000 / 10 000 |
| `ref_beta_seed`                | Initial bracket `(β_lo, β_hi)` for the reference scan     | (0.20, 0.32) |
| `n_traj_prod`                  | Production trajectories per optimizer evaluation          | 30 000 |
| `n_traj_scan_coarse/fine`      | Per-point trajectories during each per-eval β_c scan      | 4 000 / 10 000 |
| `beta_seed`                    | Initial bracket for the very first eval (warm-started after)  | (0.20, 0.32) |
| `x0`                           | Starting `(r1, r2)`                                       | (1.0, 1.0) |
| `max_evals`                    | Hard cap on optimizer evaluations                         | 30 |
| `nm_sigma0`                    | Initial NM simplex side length                            | 0.10 |
| `nm_xatol / nm_fatol`          | NM convergence tolerances                                 | 0.005 / 1e-6 |
| `nm_shrink`                    | Simplex shrink coefficient                                | 0.75 |
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
cost.py           ← plain L² cost on three boundary directions
evaluator.py      ← (r1, r2) → (cost, σ_cost): wraps β_c scan + production MC
mc_engine.py      ← Gram-Charlier β_c finder + simulator wrapper
optimizer.py      ← hand-coded Nelder-Mead simplex
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
  (`ref_n_traj` trajectories at `ref_L²` sites).  It is cached, so
  subsequent runs at the same `ref_L` reuse it.

- **Per-evaluation cost** scales as `n_traj_prod × ref_L²`, plus a few
  thousand trajectories for the β_c scan.  At `ref_L=24, n_traj_prod=30 000`,
  expect ~80 s per evaluation on a single core.

- **If the optimizer seems stuck**, look at the SNR column in the log.
  If SNR < 2 most of the time, your cost is being swamped by MC noise
  — increase `n_traj_prod` or raise `snr_target`.

- **To monitor progress live**:
  ```
  tail -f results/nelder_mead/eval_log.jsonl
  open  results/nelder_mead/frames/optimizer_nelder_mead.png   # macOS
  xdg-open results/nelder_mead/frames/optimizer_nelder_mead.png  # Linux
  ```
