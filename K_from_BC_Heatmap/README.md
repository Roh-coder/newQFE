# K_from_BC_Heatmap — Plug-and-play continuum cost landscapes

Generate cost landscapes over the (`r1`, `r2`) coupling-ratio plane that
**compare an untwisted square torus to a twisted parallelogram torus** at
the same `(r1, r2)`.  Each landscape point is extrapolated to the
**continuum (L → ∞)** from at least two finite lattice sizes, so the
heatmap shows a physical, lattice-spacing-independent mismatch.

This directory is self-contained: it carries its own copy of the C++
simulator binary and the relevant Python modules (`mc_engine.py`,
`cost.py`).  No path-fragile imports from the rest of the repo.

## TL;DR

```bash
# 1. End-to-end smoke test (~1–2 minutes)
bash run_smoke.sh

# 2. Production via unified orchestrator (recommended)
python3 run.py --tag prod --sizes 12 16 24 \
               --ref-Tx-frac 0.25 --ref-Ty-frac 0.25 \
               --r-min 0.5 --r-max 3.0 --r-step 0.25 \
               --n-traj 50000 --n-workers 8

# 2b. Or call the three scripts individually:
python3 precompute_grid.py  --tag prod --sizes 12 16 24 \
                            --ref-Tx-frac 0.25 --ref-Ty-frac 0.25 \
                            --r-min 0.5 --r-max 3.0 --r-step 0.25 \
                            --n-traj 50000 --n-workers 8
python3 compute_landscape.py --tag prod --cost-mode huber_log
python3 plot_landscape.py    --tag prod --cost-mode huber_log
```

Outputs land under `results/<tag>/` (data) and `results/<tag>/plots/`
(figures).  The pickle grid in `results/<tag>/grid/L*/{test,ref}/` is
**resumable** — re-running `precompute_grid.py` skips finished points.

## What it does

For every `(r1, r2)` on the requested grid, and for every lattice size
`L ∈ --sizes`, the precompute step runs the C++ simulator
(`bin/ising_tri_twisted_parallelogram`) twice:

| label  | geometry                 | (Lx, Ly, Tx, Ty)                              |
|--------|--------------------------|-----------------------------------------------|
| `test` | configurable (default: untwisted square) | `(round(lx_mult·L), round(ly_mult·L), round(tx_frac·L), round(ty_frac·L))` |
| `ref`  | configurable (default: quarter-twist)    | same formula, separate multipliers/fracs      |

with the default parameters `test=(1,1,0,0)` and `ref=(1,1,0.25,0.25)`.
At each point we

1. find the finite-size pseudo-critical β_c via the multi-pass
   Gram–Charlier susceptibility-peak finder (`mc_engine.find_beta_c`),
2. run a production MC at β_c (`--n-traj` trajectories),
3. dump the all-to-all two-point correlator into a pickle.

`compute_landscape.py` then loads the matched test/ref pickles per `L`,
evaluates a chosen cost-mode (currently the `cost.py` zoo: `huber_log`
(default), `test_native`, `affine_fit`, `spectral`,
`l4mean_both_interp`), and **fits cost(L) versus 1/L** (linear by
default, quadratic optional) to obtain `cost(r1, r2; L → ∞)` and an
extrapolation uncertainty.

`plot_landscape.py` renders one heatmap per `L`, one for the continuum
extrapolation, and one for its uncertainty σ.  The argmin of each
heatmap is annotated with a white star.

## Output layout

```
results/<tag>/
├── manifest.json                 # all sweep parameters + per-L geometries
├── grid/
│   ├── L<L1>/{test,ref}/r1_X.XXX_r2_Y.YYY.pkl
│   ├── L<L2>/{test,ref}/...
│   └── ...
├── cost_<mode>.npz               # per-L + continuum cost arrays
├── cost_<mode>_fits.json         # per-(r1,r2) continuum-fit details
└── plots/
    ├── landscape_<mode>.png      # all-L panel + continuum + sigma
    └── continuum_<mode>.png      # standalone continuum heatmap
```

The `.npz` is keyed:

* `rs`               — 1-D coupling grid
* `R1`, `R2`         — meshgrids
* `sizes`            — list of `L`
* `test_Lx_mult`, `test_Ly_mult`, `test_Tx_frac`, `test_Ty_frac` — test geometry params
* `ref_Lx_mult`, `ref_Ly_mult`, `ref_Tx_frac`, `ref_Ty_frac`     — ref geometry params
* `cost_L<L>`        — finite-L cost arrays (shape `(N, N)`)
* `sigma_L<L>`       — corresponding 1-σ statistical uncertainty
* `betac_test_L<L>`, `betac_ref_L<L>` — pseudo-critical β arrays
* `cost_inf`         — continuum-extrapolated cost
* `sigma_inf`        — extrapolation uncertainty
* `slope_inf`        — leading 1/L slope from the fit
* `ok_inf`           — boolean mask (False if too few `L` to fit)

## CLI reference

### `precompute_grid.py`

| flag               | default              | meaning                                        |
|--------------------|----------------------|------------------------------------------------|
| `--sizes`          | `8 12 16`            | lattice sizes for FSS extrapolation             |
| `--test-Lx-mult`   | `1.0`                | test Lx = round(M·L)                           |
| `--test-Ly-mult`   | `1.0`                | test Ly = round(M·L)                           |
| `--test-Tx-frac`   | `0.0`                | test Tx = round(F·L) — 0 = untwisted           |
| `--test-Ty-frac`   | `0.0`                | test Ty = round(F·L)                           |
| `--ref-Lx-mult`    | `1.0`                | ref  Lx = round(M·L)                           |
| `--ref-Ly-mult`    | `1.0`                | ref  Ly = round(M·L)                           |
| `--ref-Tx-frac`    | `0.25`               | ref  Tx = round(F·L)                           |
| `--ref-Ty-frac`    | `0.25`               | ref  Ty = round(F·L)                           |
| `--r-min/max/step` | `0.5` / `3.0` / `0.5`| (r1, r2) grid                                  |
| `--n-traj`        | `20000`              | production trajectories per point    |
| `--n-traj-coarse` | `2000`               | β_c finder coarse-pass trajectories  |
| `--n-traj-fine`   | `4000`               | β_c finder fine-pass trajectories    |
| `--beta-seed`     | `0.05 0.40`          | initial β bracket                    |
| `--n-workers`     | `4`                  | parallel MC processes                |
| `--tag`           | `default`            | output dir name under `results/`     |
| `--exe`           | `bin/...`            | override simulator path              |

### `run.py` (unified orchestrator)

Wraps all three stages.  All parameters from the three scripts are exposed
as CLI flags and also as a `CONFIG` dict at the top of the file.

| flag group         | flags                                                                 |
|--------------------|-----------------------------------------------------------------------|
| bookkeeping        | `--tag`, `--exe`, `--stages precompute compute plot`                  |
| coupling grid      | `--r-min`, `--r-max`, `--r-step`                                      |
| FSS sizes          | `--sizes L [L ...]`                                                   |
| test geometry      | `--test-Lx-mult`, `--test-Ly-mult`, `--test-Tx-frac`, `--test-Ty-frac` |
| ref geometry       | `--ref-Lx-mult`, `--ref-Ly-mult`, `--ref-Tx-frac`, `--ref-Ty-frac`   |
| MC                 | `--n-traj`, `--n-traj-coarse`, `--n-traj-fine`, `--n-skip`, `--n-therm`, `--beta-seed LO HI`, `--n-workers` |
| cost replay        | `--cost-modes MODE [MODE ...]`, `--cost-power`, `--fit linear\|quadratic` |
| plotting           | `--plot-log`, `--plot-vmax`                                           |
| monitoring         | `--status-interval` (seconds between heartbeat lines)                 |

Run a subset of stages:
```bash
# overnight MC only
python3 run.py --stages precompute --tag prod --n-traj 50000
# re-plot with a different cost mode, no MC
python3 run.py --stages compute plot --tag prod --cost-modes spectral
```

### `compute_landscape.py`

| flag            | default     | meaning                                           |
|-----------------|-------------|---------------------------------------------------|
| `--tag`         | `default`   | sweep tag (must match precompute)                 |
| `--cost-mode`   | `huber_log` | one of the modes in `cost.py`                     |
| `--cost-power`  | `2`         | residual exponent for `test_native`               |
| `--fit`         | `linear`    | `linear` (a + b/L) or `quadratic` (+ c/L²)        |
| `--out-name`    | derived     | basename for `.npz` / `.json`                     |

### `plot_landscape.py`

| flag          | default     | meaning                              |
|---------------|-------------|--------------------------------------|
| `--tag`       | `default`   | sweep tag                            |
| `--cost-mode` | `huber_log` | which `.npz` to plot                 |
| `--log`       | off         | render `log10(cost)`                 |
| `--vmax`      | 99th pct    | cap colour scale                     |

## Production sizing

**Continuum quality** is set by `--sizes` (need ≥3 for a quadratic fit;
2 minimum for a linear fit) and `--n-traj` (cost noise floor scales like
`1/√n_traj`).  Earlier landscape work in
`K_from_Optimizer_Production/LANDSCAPE.md` showed that `n_traj = 5000`
at a single `L` left no resolvable basin; for production use **at least
`n_traj ≥ 50000` and ≥ 3 sizes**.

Rough wall-time budget at `n_traj = 50000` on 8 workers, 100-point grid,
3 sizes (test+ref ⇒ 600 MC runs) is on the order of one CPU-day.
Resumability means you can stop / restart safely.

## Sanity checks

* Set `--ref-Tx-frac 0 --ref-Ty-frac 0` (matching test defaults) to give
  test = ref geometry — the continuum cost should vanish inside statistical
  noise.  Use this to calibrate the noise floor for your `--n-traj` choice.
* `betac_test_L<L>` and `betac_ref_L<L>` arrays in the `.npz` should be
  smooth in `(r1, r2)` and decreasing (roughly) along the diagonal —
  they are a free diagnostic that the β_c finder did not blow up.

## Smoke test

`run_smoke.sh` runs the entire pipeline at tiny sizes (`L ∈ {6, 8}`,
`r ∈ {0.5, 1.0, 1.5}`, `n_traj = 800`, ~35 s of MC on 2 workers) and
exercises both `huber_log` and `test_native` cost modes.  The five
output PNGs land under `results/smoke/plots/`.

## Related work in this repo

* `K_from_Optimizer_Production/precompute_landscape.py` —
  single-`L`, single-cost landscape; this directory generalises that
  workflow to multiple `L` + continuum extrapolation + matched
  twisted/untwisted geometry pair.
* `K_from_Optimizer_Production/LANDSCAPE.md` — discussion of why the
  argmin of a noisy cost landscape can drift away from physical truth at
  low statistics; informs the production-sizing guidance above.
* `K_from_BC/` — the original boundary-condition coupling-calibration
  workflow; this directory shares no source files with it (vendored
  copies only) so it can be edited independently.
