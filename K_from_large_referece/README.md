# large_reference_workflow

A collaborator-friendly, fully standalone large-reference workflow for boundary-correlator scoring.

This directory is intentionally simple and explicit:
- Step 1 builds exactly one large reference all-to-all two-point payload.
- Step 2 builds test payloads over user-chosen lattice sizes and coupling-ratio grid.
- Step 3 fits each test channel to the continuum with
  A + B * (1 / L)^C,
  then compares those continuum channels to the same fractional channels of the large reference.

No reference-side continuum extrapolation is used here.

## Why this exists

The original runner in this repository is powerful but monolithic:
- one script handles build, run scheduling, multiple reference modes, monitor UI, and scoring;
- many nested config branches are valid, which increases accidental complexity for handoff.

This replacement workflow favors transparency over flexibility:
- each script has one job;
- each step writes a manifest with paths and summary stats;
- outputs are plain `.dat` data artifacts for auditability (no pickle payloads).

This directory is now self-contained:
- local simulator source, headers, Makefile, and optional binary are included;
- local Python helper modules (`lib/mc_engine.py`, `lib/cost.py`) are included;
- no parent-directory runtime files are required.

All results are organized into three explicit buckets under each run tag:
- `reference_data`
- `test_data`
- `comparison_analysis_data`

## Directory layout

```text
large_reference_workflow/
  01_generate_reference.py
  02_generate_tests.py
  03_score_continuum_vs_reference.py
  workflow_common.py
  Makefile
  bin/
  include/
  lib/
  src/
  configs/
    reference_example.json
    tests_example.json
    score_example.json
  results/
    <run_tag>/
      reference_data/
      test_data/
      comparison_analysis_data/
```

## Prerequisites

From this folder only:
- Python packages: `numpy`, `scipy`, `matplotlib`
- Build toolchain for simulator (`make`, C++14 compiler)

## Quick start

Run from this directory:

```bash
cd /path/to/large_reference_workflow

python -m pip install -r requirements.txt

python 01_generate_reference.py --config configs/reference_example.json
python 02_generate_tests.py --config configs/tests_example.json
python 03_score_continuum_vs_reference.py --config configs/score_example.json
```

Optional explicit local build:

```bash
make
```

If you use a pre-existing virtual environment:

```bash
/path/to/venv/bin/python 01_generate_reference.py --config configs/reference_example.json
/path/to/venv/bin/python 02_generate_tests.py --config configs/tests_example.json
/path/to/venv/bin/python 03_score_continuum_vs_reference.py --config configs/score_example.json
```

## Step-by-step config guide

### Step 1 config: `reference_example.json`

Edit only these sections:
- `reference_lattice`: `Lx`, `Ly`, `Tx`, `Ty`
- `reference_couplings`: `k3`, `k1_over_k3`, `k2_over_k3` (usually isotropic 1,1)
- `beta_finder`
- `mc`
- `run.tag`

Main output:
- `results/<tag>/reference_data/reference_all_to_all.dat`
- `results/<tag>/reference_data/reference_all_to_all.meta.json`
- `results/<tag>/reference_data/reference_channels.dat`
- `results/<tag>/reference_data/manifest_reference.json`

### Step 2 config: `tests_example.json`

Edit only these sections:
- `test_family.sizes`
- `test_family.geometry_defaults` or `test_family.geometry_map`
- `couplings.k1_over_k3_values`
- `couplings.k2_over_k3_values`
- optional `beta_strategy`
- `beta_finder`
- `mc`
- `execution.n_workers`
- `run.tag`

Main output:
- many `.dat` test correlator files under `results/<tag>/test_data/grid/L*/`
- sidecar metadata files `*.meta.json` alongside each `.dat` file
- `results/<tag>/test_data/continuum_beta_extrapolations.json` when using
  `beta_strategy.mode = "free_power_continuum"`
- `results/<tag>/test_data/manifest_tests.json`

Optional Step 2 beta strategy:
- default legacy mode: per-size pseudo-critical production (`beta_strategy`
  omitted or `mode = "per_size_beta"`)
- continuum-beta mode: find `beta_c(L)` for each requested size, fit
  `beta_c(infinity) = beta_c(L) + a L^{-p}` with free power scaling, then run
  production at the shared `beta_c(infinity)` for every size by setting
  `beta_strategy.mode = "free_power_continuum"`

In continuum-beta mode, each payload metadata file records both:
- `finite_size_beta_c`
- `production_beta`

### Step 3 config: `score_example.json`

Edit only these sections:
- `input.reference_data_file`
- `input.reference_metadata_file` (optional; defaults to sidecar `*.meta.json`)
- `input.tests_manifest`
- `analysis.k_values`
- `analysis.k_denominator`
- `analysis.power_fit` (`c_min`, `c_max`, `c_initial`)
- `analysis.fss_plots` (optional: `enabled`, `dpi`, `max_points`)
- `analysis.weighted`
- `run.tag`

Main output:
- `results/<tag>/reference_data/reference_channels_scored.dat`
- `results/<tag>/test_data/raw_test_channels_scored.dat`
- `results/<tag>/comparison_analysis_data/continuum_test_channels.dat`
- `results/<tag>/comparison_analysis_data/channel_comparison.dat`
- `results/<tag>/comparison_analysis_data/score_map.dat`
- `results/<tag>/comparison_analysis_data/score_minimum.dat` (if finite scores exist)
- `results/<tag>/comparison_analysis_data/fss_plot_data.dat`
- `results/<tag>/comparison_analysis_data/fss_plot_index.dat`
- `results/<tag>/comparison_analysis_data/fss_plots/*.png`
- `results/<tag>/comparison_analysis_data/score_heatmap.png`
- `results/<tag>/comparison_analysis_data/zscore_heatmap.png`
- `results/<tag>/comparison_analysis_data/score_zscore_heatmaps.png`
- `results/<tag>/comparison_analysis_data/manifest_score.json`

FSS plot behavior in Step 03:
- one FSS figure is saved per `(k1/k3, k2/k3)` point;
- each panel shows one channel `(cycle, k)`;
- test-size points, fitted curve `A + B(1/L)^C`, continuum intercept `A`, and the large-reference point are all shown together.

## Logic and formulas

For each boundary cycle and each fractional sampling point `t = k / denominator`:

1. Collect test-channel values across sizes `L`.
2. Fit the continuum model:

   `G_test(L) = A + B * (1 / L)^C`

3. Interpret `A` as the continuum test estimate at that `(cycle, t)`.
4. Compare to large-reference channel `G_ref(cycle, t)`.

Per-point score is a channel sum of squared differences:
- unweighted: `sum((A - G_ref)^2)`
- weighted (optional): `sum((A - G_ref)^2 / (sigma_A^2 + sigma_ref^2))`

`z_rms` is also reported from channel z-values where variances are finite.

## Guardrails for collaborators

- Strict required-field checks with readable errors.
- Resume support (`paths.resume`) for long runs.
- Atomic payload writes (`.tmp` then replace) to reduce corruption risk.
- Manifests at every step with paths and summary counts.
- Missing-size warnings at scoring time, with explicit `n_missing_sizes` in `score_map.dat`.

## Notes

- Relative paths in configs are resolved relative to this directory.
- Step 3 can score partial test sets, but continuum fits are better with at least 3 sizes.
- If only two sizes are available for a channel, the fitter falls back to fixed `C=1`.

## Journal

### 2026-05-12: high-stat production run (5 test sizes)

Production configs used:
- `configs/production_reference_20260512.json`
- `configs/production_tests_20260512.json`
- `configs/production_score_20260512.json`

Execution summary:
- Step 1 reference: completed (`run_tag=production_ref_highstats_20260512`)
  - lattice: `80x80`, `(Tx,Ty)=(0,0)`
  - couplings: `(k1/k3, k2/k3)=(1.0, 1.0)`, `k3=1.0`
  - result: `beta_c = 0.270654881862 ± 0.001153760683`
  - payload wall time: `182.536 s`
- Step 2 tests: completed (`run_tag=production_tests_highstats_20260512`)
  - sizes: `[8, 12, 16, 24, 32]`
  - coupling grid: `k1/k3, k2/k3 in {0.9, 1.0, 1.1}` (`3x3` points)
  - jobs: `45/45` completed, `0` errors, `0` skipped
  - summed job wall time: `1084.096 s` (`0.301 h`)
- Step 3 scoring: completed (`run_tag=production_compare_highstats_20260512`)
  - channels scored per point: `21`
  - grid points scored: `9/9`
  - missing-size count: `0` for all points
  - fit mode usage: all `189` channel fits used full power-law mode `A + B*(1/L)^C`

Primary scoring result (unweighted score map):
- best point: `(k1/k3, k2/k3) = (1.0, 1.0)`
- minimum score: `0.04608512592`
- corresponding `z_rms`: `1.444351808`
- second-best score: `0.05834294903` at `(1.1, 1.0)`
- score separation: `+0.01225782311` absolute (`1.266x` relative to minimum)

Interpretation:
- The production scan cleanly prefers the isotropic point `(1.0, 1.0)` under the current unweighted objective.
- The score landscape has strong contrast (worst/min ratio about `3.33e3`), indicating robust rejection of much of the scanned anisotropic grid.
- `z_rms` and raw score are not strictly aligned across points because channel uncertainties vary strongly by point; use both diagnostics together.

Artifacts of interest:
- `results/production_compare_highstats_20260512/comparison_analysis_data/score_map.dat`
- `results/production_compare_highstats_20260512/comparison_analysis_data/score_minimum.dat`
- `results/production_compare_highstats_20260512/comparison_analysis_data/score_zscore_heatmaps.png`
- `results/production_compare_highstats_20260512/comparison_analysis_data/fss_plot_index.dat`
- `results/production_compare_highstats_20260512/comparison_analysis_data/fss_plots/`

Data format check:
- Verified no `.pkl` files in `production_ref_highstats_20260512`, `production_tests_highstats_20260512`, or `production_compare_highstats_20260512`.

### 2026-05-13: FSS / autocorrelation diagnostic — root cause of rough score landscape

**TL;DR.** The roughness in the (r1, r2) score landscape is **not** caused by MC
autocorrelation, undertherma­lization, or interpolation. It is caused by the
per-cell `find_beta_c` step running with too few trajectories
(`n_traj_coarse=120`, `n_traj_fine=240` in `configs/production_tests_20260512.json`).
Two effects compound:
1. Genuine FSS shift of `β_χmax(L)`: at small L the χ-peak sits well below
   `β_c(∞)` (e.g. ≈ −16% at L=8 in our test). This is real physics.
2. Stochastic noise on `β_χmax(L)`: the 240-trajectory peak fit gives σ_β of
   order 1e-3–5e-3 per cell, which is enough to scatter `G(t, L)` between
   neighbouring (r1, r2) cells.

#### Diagnostic harness

New script:
- [diag_fss_autocorr.py](K_from_large_referece/diag_fss_autocorr.py) — per L, runs
  the simulator at the symmetric isotropic point K=(1,1,1) with
  `--dump_traces 1` so per-trajectory `(E, |m|, m²)` are dumped, then
  computes integrated autocorrelation time τ_int (Madras–Sokal automatic
  windowing), 5-block thermalization drift, and FSS of the connected
  two-point function at fractional separations `t = k/8` along all three
  cycles. Has two modes:
  - default: uses the analytic `β_c = ln 3 / 4` (no β-finder),
  - `--find-beta`: runs `mc_engine.find_beta_c` per L with configurable
    `--bf-n-traj-coarse` / `--bf-n-traj-fine` / `--bf-n-coarse` /
    `--bf-n-refine{,2,3}` / `--beta-{lo,hi}`,
  - `--find-beta --beta-from-largest-size`: finds `β_c` once at
    `L=max(sizes)` and reuses that same value for every smaller lattice.

Helper:
- [diag_replot_thermalization.py](K_from_large_referece/diag_replot_thermalization.py)
  — re-renders running-mean and 100-block thermalization plots from the
  scratch traces of an existing diagnostic run.
- [diag_fit_beta_continuum.py](K_from_large_referece/diag_fit_beta_continuum.py)
  — reads a diagnostic `raw/summary.dat`, fits `beta_used(L)` with a
  Taylor-style polynomial in `1/L`, writes JSON / `.dat` fit artifacts, and
  optionally reports fit-window drift via `--window-scan`.

#### Runs completed today

Both runs used L ∈ {8, 16, 24, 32, 40, 48, 56}, K=(1,1,1), `n_therm=4000`,
`n_traj=8000`, `n_skip=1`, 4 workers.

1. **Analytic β_c (control)** — `results/diag_fss_autocorr/`:
   - τ_int(|m|): 0.51 (L=8) → 1.12 (L=56). τ_int(E) similar.
   - therm drift ≤ 0.06 σ for every L.
   - FSS of G(t, L) at every (cycle, k) is **smooth and monotone in 1/L**;
     the three cycles agree as required by C₃ symmetry.
   - Conclusion: at fixed β, MC quality is excellent; statistical errors do
     not produce the observed roughness.

2. **β found per L with the production β-finder settings**
   (`--find-beta`, `--bf-n-traj-coarse 120`, `--bf-n-traj-fine 240`) —
   `results/diag_fss_autocorr_findbeta/`:
   - β-finder results vs analytic 0.27465:
     L=8 → 0.2296 (−16.4%), L=16 → 0.2522 (−8.2%), L=24 → 0.2611 (−4.9%),
     L=32 → 0.2688 (−2.1%), L=40 → 0.2693 (−1.9%), L=48 → 0.2688 (−2.1%),
     L=56 → 0.2688 (−2.1%).
   - τ_int and therm drift are essentially unchanged from the control.
   - FSS of G(t, L) is **non-monotone**: every (cycle, k) panel shows a bump
     centered around L ≈ 24–32. This is the same character as the
     score-landscape roughness.

Side-by-side comparison plots:
- [analytic-β FSS](K_from_large_referece/results/diag_fss_autocorr/plots/fss_two_point.png)
- [β-finder FSS](K_from_large_referece/results/diag_fss_autocorr_findbeta/plots/fss_two_point.png)
- [τ_int(L)](K_from_large_referece/results/diag_fss_autocorr/plots/tau_int_vs_L.png)

Raw outputs:
- `results/diag_fss_autocorr{,_findbeta}/raw/summary.dat`
  (columns: `L  n_meas  tau_int_E  tau_int_|m|  mean_|m|  sd_|m|  therm_drift_sigmas  beta_used  beta_sigma`)
- `results/diag_fss_autocorr{,_findbeta}/raw/fss_G.dat`
  (columns: `L  cycle  k  t  G  sG`)

### 2026-05-14: FSS diagnostic — shared `β_c` from the largest lattice

**TL;DR.** Finding `β_c` once at the largest test lattice and reusing it for
all smaller sizes removes the rough FSS bump completely at the same low
β-finder statistics. It is therefore a much cleaner diagnostic than tuning
each size to its own pseudo-critical point individually, but it is not the
same physical prescription.

Run completed:
- shared-largest-L β mode — `results/diag_fss_autocorr_findbeta_sharedL56/`
  - command:

    ```bash
    cd /workspaces/newQFE/K_from_large_referece
    /workspaces/newQFE/.venv/bin/python diag_fss_autocorr.py \
        --sizes 8 16 24 32 40 48 56 \
        --n-therm 4000 --n-traj 8000 --n-skip 1 \
        --workers 4 --find-beta --beta-from-largest-size \
        --out-root /workspaces/newQFE/K_from_large_referece/results/diag_fss_autocorr_findbeta_sharedL56
    ```
  - shared β-finder result at `L=56`:
    - `β_c(L=56) = 0.26947454`
    - shift from exact `ln(3)/4 = 0.27465307`: `−1.89%`
    - β-finder scan work: `32` scan points, `6360` β-finder trajectories total
  - production autocorrelation / therm check:
    - `τ_int(|m|)` runs from `0.55` (`L=8`) to `1.32` (`L=56`)
    - thermalization drift stays below `0.10 σ` for every `L`

Direct comparison to per-L pseudo-critical tuning (`results/diag_fss_autocorr_findbeta/`):
- FSS smoothness:
  - shared-largest-L β: all `21/21` `(cycle, k)` channels are monotone in `1/L`
  - per-L β tuning: all `21/21` channels show a sign flip (the same mid-size bump)
- roughness metrics from `raw/fss_G.dat`:
  - mean total-variation ratio `TV / |G(L_max)-G(L_min)|`:
    - shared-largest-L β: `1.000`
    - per-L β tuning: `2.554`
  - mean RMS second difference across channels:
    - shared-largest-L β: `0.0284`
    - per-L β tuning: `0.0633`
- closeness to the analytic-β control (`results/diag_fss_autocorr/`):
  - global RMSE of `G` values:
    - shared-largest-L β: `0.1271`
    - per-L β tuning: `0.2323`

Interpretation:
- Reusing `β_c(L=56)` for every size almost entirely suppresses the β-finder
  jitter that roughened the original FSS curves.
- This method is therefore better if the immediate goal is a smooth
  common-β diagnostic or a score-landscape sanity check.
- It is **not** a drop-in replacement for true pseudo-critical sampling,
  because small lattices are no longer run at their own `β_χmax(L)`.
  Instead they are forced onto a common β tied to the largest lattice.

Artifacts of interest:
- `results/diag_fss_autocorr_findbeta_sharedL56/raw/summary.dat`
- `results/diag_fss_autocorr_findbeta_sharedL56/raw/fss_G.dat`
- `results/diag_fss_autocorr_findbeta_sharedL56/plots/fss_two_point.png`

#### Next steps (resume here)

The remaining unresolved question is whether higher-stat per-L β-finding can
recover the same smoothness **without** giving up pseudo-critical tuning.
The next targeted run is therefore still the 10× per-L β-finder test:

```bash
cd /workspaces/newQFE/K_from_large_referece
/workspaces/newQFE/.venv/bin/python diag_fss_autocorr.py \
    --sizes 8 16 24 32 40 48 56 \
    --n-therm 4000 --n-traj 8000 --n-skip 1 \
    --workers 4 --find-beta \
    --bf-n-traj-coarse 1200 --bf-n-traj-fine 2400 \
    --out-root /workspaces/newQFE/K_from_large_referece/results/diag_fss_autocorr_findbeta_10x \
  2>&1 | grep -E '^\[diag\]|find_beta_c done' | tail -40
```

Notes for resuming:
- 10× β-finder stats (coarse 120 → 1200, fine 240 → 2400) — keeps n_refine
  passes at 5/5/5; only the per-point trajectory count goes up.
- Estimated wall: at 4 workers the bottleneck is L=56. The 1× β-finder run
  took 143 s of β-finder + 99 s of production at L=56. Scaling β-finder by
  10× ⇒ ~24 min at L=56; total wall ≈ 25–30 min.
- If 10× still leaves residual roughness, the next experiment is to either
  (a) raise `n_refine{,2,3}` so the third pass narrows further around the
  true peak, or (b) bypass per-cell β-finding entirely by reweighting onto
  a common reduced-temperature target (the simulator already supports
  `--dump_traces 1`; `lib/mc_engine.py` has reweight infrastructure used by
  `find_beta_c_reweight` and `find_beta_c_multidonor*`).
- After the new run, rerun `03_score_continuum_vs_reference.py` on a fresh
  `02_generate_tests.py` output configured with the higher β-finder stats
  to confirm the score landscape smooths.

Open questions deferred for the next session:
- Does the score landscape track the β-finder stats monotonically, or is
  there a knee where extra trajectories stop helping? (Run 5×, 10×, 20× to
  bracket.)
- Is the genuine FSS shift `β_χmax(L) − β_c(∞)` better absorbed by an FSS
  extrapolation of β_c(L) before sampling G, instead of brute-force higher
  stats?

### 2026-05-14: notes on `β_c` extrapolation from pseudo-critical `χ` peaks

This discussion clarified what is and is not legitimate when the end goal is
the continuum-limit critical two-point function.

Shared-`β_c` diagnostic vs production use:
- Reusing one `β_c(L_max)` for all sizes is legitimate as a **diagnostic** for
  suppressing β-finder jitter and checking whether rough FSS curves are caused
  by noisy per-size tuning.
- It is **not** a principled production prescription for the critical
  continuum two-point function, because it does not keep the thermal scaling
  variable fixed across sizes. Small lattices are then sampled at a β tied to
  `L_max`, not at their own pseudo-critical point and not necessarily at
  `β_c(∞)`.

How to use pseudo-critical betas properly:
- If `β_pc(L)` is measured consistently for many lattice sizes, it is valid to
  extrapolate those pseudo-critical couplings to an infinite-volume critical
  `β_c`.
- The shared-`β_c` run
  (`results/diag_fss_autocorr_findbeta_sharedL56/raw/summary.dat`) is **not**
  useful for this extrapolation because `beta_used` is constant by
  construction.
- The relevant input is the per-L β-finder run in
  `results/diag_fss_autocorr_findbeta/raw/summary.dat`.

Generic FSS ansatz for `χ`-peak pseudo-critical couplings:
- For a standard second-order critical point,

  ```text
  β_χ,peak(L) = β_c + A L^(-1/ν) + B L^(-1/ν-ω) + C L^(-2/ν) + ...
  ```

- Here `ν` is the thermal exponent and `ω` is the leading irrelevant
  correction exponent.
- In CFT language,

  ```text
  1/ν = d - Δ_ε
  ω   = Δ_irr - d
  ```

  where `Δ_ε` is the scaling dimension of the thermal operator and `Δ_irr` is
  the dimension of the leading symmetry-allowed irrelevant operator.
- For a generic model there is an infinite tower of irrelevant corrections in
  principle, but practical fits usually keep one leading thermal term plus one
  or two subleading corrections.

Practical minimal ansatz:
- For an Ising-like case with `ν = 1`, the standard first working model is

  ```text
  β_χ,peak(L) = β_c + a/L + b/L^2
  ```

- Treating `β_c` as the constant term of a plain polynomial in `1/L` is often
  acceptable as a pragmatic estimator in this setting, but the physics
  justification is finite-size scaling, not analyticity in `1/L` by itself.

Fit performed on the existing diag-FSS per-L β sequence (`L = 8,16,24,32,40,48,56`):
- using

  ```text
  β_pc(L) = β_c + a/L + b/L^2
  ```

- best-fit parameters:
  - `β_c = 0.277710204887 ± 0.003048617073`
  - `a   = -0.388027394432 ± 0.121968973427`
  - `b   =  0.014817478950 ± 0.833838125527`
- fit quality:
  - RMS residual: `1.642672139218e-03`
  - exact triangular Ising value: `β_c = ln(3)/4 = 0.274653072167`
  - fitted intercept shift from exact: `+0.003057132720`

Interpretation of that fit:
- The `8–56` data can be forced into the `a/L + b/L^2` ansatz, but the `b`
  coefficient is very poorly constrained.
- Earlier fit-window tests also showed strong drift in the extrapolated
  intercept when changing the minimum fitted size.
- Conclusion: the current low-stat per-L β-finder data are sufficient for a
  diagnostic extrapolation, but not yet for a robust production determination
  of `β_c(∞)`.

Recommended next use of this note:
- After the planned `10×` per-L β-finder run, repeat the same `β_pc(L)`
  extrapolation and check whether the fitted `β_c` stabilizes across fit
  windows and across `a/L` vs `a/L + b/L^2` ansätze.
- Command for the current helper:

  ```bash
  cd /workspaces/newQFE/K_from_large_referece
  /workspaces/newQFE/.venv/bin/python diag_fit_beta_continuum.py \
      --summary results/diag_fss_autocorr_findbeta/raw/summary.dat \
      --degree 2 --window-scan
  ```

- After the `10×` run completes, point the same helper at
  `results/diag_fss_autocorr_findbeta_10x/raw/summary.dat`.

