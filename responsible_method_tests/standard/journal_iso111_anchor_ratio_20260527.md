# Journal

## 2026-05-27: Iso111 Standard Pointwise Landscape Rechecked With Anchor-Ratio Normalization

### Goal

Test whether the standard raw pointwise target score for `iso111` is being misled by a common-mode amplitude offset against the finite `144 x 144` target lattice.

The concrete symptom was the counterintuitive raw ranking

- `(r1, r2) = (0.8, 0.8)` scoring better than `(1.0, 1.0)`

even though the isotropic point should be the most natural control candidate for the `1-1-1` family.

### Raw-Score Background

The earlier raw per-point diagnostics already suggested that the issue was not primarily a side-orientation bug.

- For `(0.8, 0.8)`, the raw target residuals had near-zero mean offset:
  - `mean_delta = -4.5536e-4`
  - `mean_abs_delta = 2.1449e-3`
  - `rms_z = 0.2626`
  - `pos_frac = 0.46875`
- For `(1.0, 1.0)`, the raw target residuals overshot the target at every sampled point:
  - `mean_delta = 6.4008e-3`
  - `mean_abs_delta = 6.4008e-3`
  - `rms_z = 0.6994`
  - `pos_frac = 1.0`

That pattern is consistent with a finite-target common-mode offset: the raw score rewards whichever candidate best cancels that offset, even if its geometry is less physically plausible.

### Implementation Change

To test that hypothesis directly, the standard pointwise scorer was extended to support normalization-aware scoring modes.

- `scripts/fit_standard_acute456_all_base_points.py` now accepts
  - `normalization_mode in {raw, anchor_ratio, l8_ratio}`
  - `anchor_m`, `anchor_n`
- `scripts/plot_standard_pointwise_score_landscapes.py` now forwards those options and writes mode-suffixed outputs.

For the main check here, the new mode was:

- `normalization_mode = anchor_ratio`
- anchor point `(m, n) = (0, -1)`

In this mode, each point is compared through the ratio `G(p) / G(anchor)`, and the anchor point itself is removed from the score so the result is not dominated by a trivial exact match at the normalization point.

### Command Used

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py \
    --family iso111 \
    --score-mode target \
    --normalization-mode anchor_ratio
```

### Artifacts Written

- figure:
  - `responsible_method_tests/standard/standard_pointwise_score_landscapes_anchor_ratio.png`
- table:
  - `responsible_method_tests/standard/iso111_pointwise_score_landscape_anchor_ratio.tsv`

### Main Quantitative Result

The specific confusing inversion is removed by anchor-ratio normalization.

| candidate | raw orbit RMS z | anchor-ratio orbit RMS z |
| --- | ---: | ---: |
| `(0.8, 0.8)` | `0.2575` | `0.5507` |
| `(1.0, 1.0)` | `0.7089` | `0.3510` |

So after removing the common-mode amplitude degree of freedom, the isotropic point becomes clearly better than `(0.8, 0.8)`.

### Full Anchor-Ratio Landscape Ranking

Top five candidates in the `iso111` anchor-ratio target landscape:

| rank | `(r1, r2)` | orbit RMS z |
| --- | --- | ---: |
| 1 | `(1.1, 1.1)` | `0.1728` |
| 2 | `(1.0, 0.9)` | `0.1815` |
| 3 | `(0.9, 1.0)` | `0.2083` |
| 4 | `(0.9, 0.9)` | `0.2310` |
| 5 | `(1.1, 1.0)` | `0.2552` |

Two control ranks of interest:

- `(1.0, 1.0)` ranks `11 / 25` with orbit RMS z `0.3510`
- `(0.8, 0.8)` ranks `17 / 25` with orbit RMS z `0.5507`

The anchor-ratio table contains `62` scored points per candidate rather than `63` because the anchor point is excluded.

### Interpretation

- The raw preference for `(0.8, 0.8)` over `(1.0, 1.0)` was real, but it was not robust to normalization.
- This strongly supports the interpretation that the raw standard target score is contaminated by a common-mode finite-target amplitude bias.
- The isotropic point is recovered as better than `(0.8, 0.8)` once that bias is removed.
- However, `(1.0, 1.0)` is still not the global minimum of the anchor-ratio landscape.
- The remaining displacement of the minimum toward `(1.1, 1.1)` means the raw amplitude bias was not the whole story. Some combination of residual shape mismatch, anchor dependence, or finite-target effects remains in the standard pointwise score.

### Validation

- `python -m py_compile responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py` passed
- File-level error check on both edited scripts reported no errors
- The figure write emitted a `tight_layout` warning, but the output PNG and TSV were written successfully

### Current Takeaway

Anchor-ratio normalization answers the narrow question cleanly:

- yes, the raw `iso111` landscape was overstating `(0.8, 0.8)` because of normalization-sensitive bias
- no, normalization alone does not fully restore `(1.0, 1.0)` as the unique best point in the standard target landscape

### Next Useful Check

- rerun the same `iso111` landscape with `normalization_mode = l8_ratio`
- compare the top few candidates point-by-point to see whether the residual preference for `(1.1, 1.1)` is anchor-specific or persists under a second normalization choice

### Resume Later

The controlling open question is now narrow:

- does the remaining shift of the standard `iso111` minimum toward `(1.1, 1.1)` survive a second normalization choice, or is it mainly an artifact of the specific anchor-ratio construction?

The files to start from next time are:

- detailed note:
  - `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`
- top-level pointer:
  - `responsible_method_tests/journal.md`
- edited scorer:
  - `responsible_method_tests/scripts/fit_standard_acute456_all_base_points.py`
- edited landscape builder:
  - `responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py`
- raw iso111 target table:
  - `responsible_method_tests/standard/iso111_pointwise_score_landscape.tsv`
- anchor-ratio iso111 target table:
  - `responsible_method_tests/standard/iso111_pointwise_score_landscape_anchor_ratio.tsv`
- anchor-ratio figure:
  - `responsible_method_tests/standard/standard_pointwise_score_landscapes_anchor_ratio.png`

Current anchor-ratio facts worth preserving:

- normalization used anchor `(m, n) = (0, -1)`
- the anchor point is excluded from scoring, so each candidate uses `62` scored points rather than `63`
- best anchor-ratio candidate: `(1.1, 1.1)` with orbit RMS z `0.1727771692`
- isotropic control: `(1.0, 1.0)` with orbit RMS z `0.3510409533`
- confusing raw control: `(0.8, 0.8)` with orbit RMS z `0.5507415728`

Most direct next command:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/plot_standard_pointwise_score_landscapes.py \
    --family iso111 \
    --score-mode target \
    --normalization-mode l8_ratio
```

If a quick two-candidate check is enough before a full sweep, use:

```bash
source /workspaces/newQFE/.venv/bin/activate
PYTHONPATH=$PWD/responsible_method_tests/scripts python - <<'PY'
import os
from fit_standard_acute456_all_base_points import _build_summary_rows, _compute_orbit_reduced_score
from plot_standard_acute456_center_fss import STANDARD_ROOT

family_root = os.path.join(STANDARD_ROOT, 'data', 'iso111', 'untwisted')
twisted_dat = os.path.join(STANDARD_ROOT, 'data', 'iso111', 'twisted', 'reference', 'Lx144_Ly144_Tx0_Ty0', 'two_point_all_to_all.dat')
twisted_lattice = (144, 144, 0, 0)

for candidate in ['r1_0p800000__r2_0p800000', 'r1_1p000000__r2_1p000000']:
    rows = _build_summary_rows(
        untwisted_dir=os.path.join(family_root, candidate),
        twisted_dat=twisted_dat,
        twisted_lattice=twisted_lattice,
        untwisted_embedding_cycles=(0, 1),
        twisted_embedding_cycles=(0, 1),
        include_origin=False,
        normalization_mode='l8_ratio',
    )
    orbit = _compute_orbit_reduced_score(rows, include_origin=False, z_key='target_z')
    print(candidate, orbit['rms_z'])
PY
```

Before expanding scope beyond `iso111`, first answer just this:

- does `l8_ratio` also rank `(1.0, 1.0)` above `(0.8, 0.8)`?
- does the full `l8_ratio` landscape still prefer `(1.1, 1.1)` or move back toward isotropic?

## 2026-05-27: Iso111 Quarter-Point Boundary FSS Jackknife Sweep In Progress

### Goal

Regenerate the `iso111` quarter-point boundary FSS sweep with true jackknife error bars, rather than the fast sweep-mode summary fallback.

This is specifically for the quarter-point plots and z-score table produced by

- `responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py --fraction 0.25`

using the new opt-in sweep flag

- `--force-jackknife-sweep`

which allows the sweep to regenerate missing single-displacement sample files instead of forcing `allow_regeneration=False`.

### Why This Matters

The already-written quarter sweep under

- `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep/`

is useful, but its error bars are not uniformly true jackknife.

- For correlator panels, cached single-displacement jackknife samples are used when present, but otherwise the sweep falls back to the stored `two_point_all_to_all.dat` summary errors.
- For ratio panels, the fast sweep uses propagated errors from those sigmas rather than full ratio jackknife.

The forced-jackknife rerun is intended to remove that ambiguity.

### Command To Run Or Resume

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py \
    --fraction 0.25 \
    --sweep-output-dir responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife \
    --force-jackknife-sweep
```

### Output Directory

Write the true-jackknife rerun to a separate directory so it does not overwrite the existing fast sweep:

- jackknife sweep target:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/`
- existing non-jackknife sweep to keep for comparison:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep/`

Expected final summary table for the jackknife rerun:

- `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/iso111_boundary_quarter_fss_zscores.tsv`

### Progress Snapshot When This Note Was Written

Quarter-point cache accounting at the time of the status check:

- total required single-displacement sample files: `453`
- already present: `24`
- still missing: `429`

At the time of the latest check, the forced-jackknife sweep had not yet written the first PNG/TSV pair in the new output directory, but the run was active and the simulator was working on the twisted-reference quarter-point sample

- `(m, n) = (36, 0)` on the `144 x 144` twisted target

which is expected to be one of the slower sample generations.

### Runtime Expectation

Best current estimate from the live run snapshot:

- rough remaining time: about `1.5 hours`
- reasonable range: `1 to 2 hours`

This estimate is only a planning number, not a guarantee. The dominant uncertainty is how much time the larger lattices and twisted-reference quarter samples take relative to the already-finished untwisted ones.

### Safe Restart Rule

If the run is interrupted, rerun the exact command above.

That is safe because the expensive part is cached in

- `responsible_method_tests/standard/_jackknife_samples/`

so regenerated single-displacement sample files do not need to be recomputed once they already exist. A rerun should therefore resume effectively from the accumulated cache, even though the plot script itself rewrites the final PNG/TSV outputs.

### Important Files To Check First Next Time

- driver script:
  - `responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py`
- current fast quarter sweep:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep/`
- target jackknife quarter sweep:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/`
- cache root:
  - `responsible_method_tests/standard/_jackknife_samples/`

### Resume Question

Before doing any new analysis, first check just these two things:

- has `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/iso111_boundary_quarter_fss_zscores.tsv` been written?
- if not, how many of the `453` quarter-point sample files are now present in `_jackknife_samples`?

### Latest Status Before Pausing

Later checks during the same run showed that the job was still healthy but still had not reached the first plot write.

- current completed cache files: `37 / 453`
- current missing cache files: `416 / 453`
- current output files in `iso111_boundary_quarter_fss_sweep_jackknife/`: `0`

At that checkpoint the active simulator subprocess was working on the twisted-reference quarter-point sample

- `(m, n) = (36, -36)` on the `144 x 144` twisted target

and the updated practical ETA was revised upward to roughly

- about `2.5 hours` remaining
- with a reasonable planning range of `2 to 3 hours`

So if this session is resumed later and the jackknife output directory is still empty, do not assume the job stalled just because no PNG/TSV has appeared yet; it may still be progressing through the expensive cached sample-generation stage.

## 2026-05-28: Standard Boundary Dual-Target Follow-up And fit_method Trial

### Goal

Extend the standard boundary analysis beyond the pointwise normalization check.

This session had three concrete aims:

- finish the `iso111` quarter-point boundary jackknife sweep and preserve its clean `chi2` summary
- generate dual-target boundary comparison products for `acute456` using both the large and smaller twisted references
- test a new `fit_method` score in which the twisted ladder and untwisted ladder are both fit first, then compared through their continuum intercepts

### Infrastructure Added Or Changed

- `responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py` now
  - accepts secondary twisted target overlays via `--secondary-twisted-{dat,lattice,label}`
  - uses dataset-aware sweep naming instead of hard-coded `iso111_*` summary stems
  - uses a tight nearest-wrapped-point fallback for the secondary target only, which is needed for the smaller `acute456` twisted lattice `(66,66,33,11)` because it does not realize the exact midpoint keys `(0.5,0.0)`, `(0.0,0.5)`, `(0.5,0.5)` exactly
  - writes separate sweep summary tables `..._large_target_zscores.tsv` and `..._small_target_zscores.tsv` when a secondary target is present
- new script `responsible_method_tests/scripts/score_standard_iso111_boundary_fit_method.py`
  - stages a standard `iso111` twisted reference ladder at `L = 32, 48, 64, 96, 144`
  - supports `fit_family in {continuum, blind_power}`
  - default `continuum` family uses `workflow_common.fit_observable_continuum_power`
  - `blind_power` preserves the original `A + B x^omega` fit on `x = 1 / sqrt(V)`
  - can render diagnostic plots for selected coupling tags

### Artifacts Written

- dual-target boundary figures:
  - `responsible_method_tests/standard/iso111_boundary_midpoint_fss_dual_target.png`
  - `responsible_method_tests/standard/acute456_boundary_midpoint_fss_dual_target.png`
- `acute456` dual-target sweep directories:
  - `responsible_method_tests/standard/acute456_boundary_midpoint_fss_sweep/`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/`
- `iso111` fit-vs-fit score tables:
  - `responsible_method_tests/standard/iso111_boundary_midpoint_fit_method_zscores.tsv`
  - `responsible_method_tests/standard/iso111_boundary_quarter_fit_method_zscores.tsv`
  - `responsible_method_tests/standard/iso111_boundary_midpoint_fit_method_blind_power_zscores.tsv`
  - `responsible_method_tests/standard/iso111_boundary_quarter_fit_method_blind_power_zscores.tsv`
- fit diagnostics:
  - `responsible_method_tests/standard/fit_method_plots/`
  - `responsible_method_tests/standard/fit_method_plots_blind_power/`

### Old Boundary chi2 And Direction-Ratio Result Worth Preserving

The most reliable old-style boundary benchmark produced in this session is the completed quarter-point jackknife sweep:

- summary table:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/iso111_boundary_quarter_fss_zscores.tsv`

Key ranking from that clean quarter-point jackknife table:

- best coupling: `(r1, r2) = (1.0, 1.1)` with `chi2_sum = 0.6071`
- isotropic control `(1.0, 1.0)` ranks `6 / 25` with `chi2_sum = 3.9524`
  - correlator contribution: `corr_chi2 = 3.1630`
  - direction-ratio contribution: `ratio_chi2 = 0.7895`

This matters because the ratio channels are not dominant, but they are also not negligible. For the isotropic control they contribute about `20%` of the total penalty. So the older direct boundary `chi2` score still contains real direction-ratio structure rather than behaving like pure noise.

### acute456 Dual-Target Sweep Findings

The exact `acute456` target couplings are `(r1, r2) = (4.702783, 7.353910)`.

Using the large twisted target summary:

- midpoint sweep:
  - best coupling: `(5.173061, 8.824692)` with `chi2_sum = 0.4806`
  - exact target rank: `5 / 25` with `chi2_sum = 1.0418`
  - target split: `corr_chi2 = 0.8675`, `ratio_chi2 = 0.1743`
- quarter sweep:
  - best coupling: `(3.762226, 5.883128)` with `chi2_sum = 0.1527`
  - exact target rank: `5 / 25` with `chi2_sum = 1.7487`
  - target split: `corr_chi2 = 1.3850`, `ratio_chi2 = 0.3637`

Using the smaller twisted target summary:

- midpoint sweep:
  - best coupling: `(5.173061, 7.353910)` with `chi2_sum = 5.8357`
  - exact target rank: `9 / 25` with `chi2_sum = 12.3344`
  - target split: `corr_chi2 = 11.5643`, `ratio_chi2 = 0.7702`
- quarter sweep:
  - best coupling: `(5.173061, 7.353910)` with `chi2_sum = 8.7707`
  - exact target rank: `2 / 25` with `chi2_sum = 9.4755`
  - target split: `corr_chi2 = 8.7466`, `ratio_chi2 = 0.7289`

These `acute456` dual-target sweeps again show that the correlator terms dominate the total `chi2`, but the ratio terms remain structured and can move the ordering at the `O(0.1-1)` level. That is consistent with the broader impression that the direction-ratio data is not empty.

### fit_method Prototype: What Was Tested

The new `fit_method` score does not compare an untwisted ladder directly against one finite twisted target point.

Instead it:

- fits the full twisted reference ladder at `L = 32, 48, 64, 96, 144`
- fits each untwisted candidate ladder separately
- compares the fitted continuum intercepts panel-by-panel through the usual `z` / `chi2` aggregation

Two fit families were tested:

- `continuum` (default): `workflow_common.fit_observable_continuum_power`, typically with `taylor2`
- `blind_power`: the original `A + B x^omega` fit in `x = 1 / sqrt(V)`

### fit_method Quantitative Results

#### Continuum-Fit Family

- midpoint:
  - best coupling: `(1.1, 0.9)` with `chi2_sum = 72.0744`
  - isotropic control `(1.0, 1.0)` rank: `5 / 25` with `chi2_sum = 82.9620`
  - score split for the best point: `corr_chi2 = 71.9796`, `ratio_chi2 = 0.0948`
- quarter:
  - best coupling: `(0.8, 1.0)` with `chi2_sum = 93.8253`
  - isotropic control rank: `10 / 25` with `chi2_sum = 118.9606`
  - score split for the best point: `corr_chi2 = 93.6373`, `ratio_chi2 = 0.1880`

#### blind_power Family

- midpoint:
  - best coupling: `(0.8, 0.8)` with `chi2_sum = 0.2561`
  - isotropic control rank: `14 / 25` with `chi2_sum = 1.7433`
  - score split for the best point: `corr_chi2 = 0.0926`, `ratio_chi2 = 0.1634`
- quarter:
  - best coupling: `(0.9, 0.8)` with `chi2_sum = 0.6309`
  - isotropic control rank: `21 / 25` with `chi2_sum = 8.5710`
  - score split for the best point: `corr_chi2 = 0.6308`, `ratio_chi2 = 0.0000`

### Why I Do Not Trust fit_method Yet

The diagnostic plots make the main defect visible.

- The raw correlator panels can often be made to track each other reasonably well.
- The ratio-panel extrapolations are much less healthy.
- In the `blind_power` plots, the ratio intercept uncertainties can blow up to absurd scales.
- In some cases the dashed twisted ratio fit is visibly pathological even when the finite-volume points themselves look harmless.

So even though `fit_method` is mathematically well-defined, it does not currently behave like a trustworthy geometry score. The comparison is too sensitive to extrapolation pathology, especially in the ratio channels, and the final ranking ends up being driven almost entirely by correlator intercept fits instead of by a balanced use of the direction-ratio information.

### Current Session Takeaway

My current view after this session is:

- I do **not** trust the new `fit_method` in its present form.
- The old direct boundary `chi2` summaries, including the direction-ratio panels, still look more promising.
- Those older scores do not trivially recover the exact target, but they show structured rankings and nonzero ratio contributions that look more like signal than pure numerical junk.

Stated plainly: the old boundary `chi2` plus direction-ratio data may have some “there there”, while the new fit-vs-fit score does not currently earn much confidence.

### Best Files To Start From Next Time

- this detailed note:
  - `responsible_method_tests/standard/journal_iso111_anchor_ratio_20260527.md`
- boundary scorer / plotter:
  - `responsible_method_tests/scripts/plot_standard_iso111_boundary_midpoint_fss.py`
- fit-vs-fit prototype:
  - `responsible_method_tests/scripts/score_standard_iso111_boundary_fit_method.py`
- clean old-style quarter jackknife summary:
  - `responsible_method_tests/standard/iso111_boundary_quarter_fss_sweep_jackknife/iso111_boundary_quarter_fss_zscores.tsv`
- `acute456` dual-target summaries:
  - `responsible_method_tests/standard/acute456_boundary_midpoint_fss_sweep/acute456_boundary_midpoint_fss_large_target_zscores.tsv`
  - `responsible_method_tests/standard/acute456_boundary_midpoint_fss_sweep/acute456_boundary_midpoint_fss_small_target_zscores.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_target_zscores.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_small_target_zscores.tsv`
- fit diagnostics:
  - `responsible_method_tests/standard/fit_method_plots/`
  - `responsible_method_tests/standard/fit_method_plots_blind_power/`

## 2026-05-28: Acute456 Small-vs-Large Target Interpretation And Minimal Large-L Campaign

### Current Interpretation

The most plausible current reading is that the `acute456` smaller twisted target is benefiting from a much shorter extrapolation lever arm, even if its twisted/untwisted cutoff mismatch is not literally zero.

For the square untwisted ladder, the fit variable is

- `x = 1 / sqrt(V) = 1 / L`

The present untwisted data stop at `L = 64`, so

- `x64 = 0.015625000`
- `x66 = 0.015151515`, so `delta(66-64) = -0.000473485`
- `x144 = 0.006944444`, so `delta(144-64) = -0.008680556`

Measured against the last observed training step `48 -> 64`, this means

- `|64 -> 66|` is only about `0.09` of one last-step spacing
- `|64 -> 144|` is about `1.67` full last-step spacings

So the smaller target is almost a local holdout, while the `144 x 144` target is a genuine long-range extrapolation beyond the current window.

### Working Hypothesis

The present large-target miss is more likely extrapolation-error dominated than cutoff-mismatch dominated.

In other words:

- the smaller target may look better not because `66 x 66` is more physically faithful than `144 x 144`
- but because `64 -> 66` is such a short bridge that fit error is tiny compared with the long `64 -> 144` bridge

If that is right, then the exact `acute456` couplings should improve materially against the large twisted target once the untwisted branch is extended to large `L` and the fit is only asked to propagate over a short local window near `144`.

### Cheapest Decisive Test

The cheapest first pass does **not** need any new twisted jobs.

The large twisted reference at `144 x 144` already exists, so the natural minimal test is to add only large untwisted sizes and see whether the exact target catches up once the comparison is local in `L`.

### Minimal Large-L Campaign

Scope the first pass narrowly:

- dataset: `acute456`
- observable surface: quarter-point boundary score first
- do **not** rerun the full `5 x 5` grid yet

Use this 3-point coupling shortlist:

- exact target: `(4.702783, 7.353910)`
- current quarter large-target winner: `(3.762226, 5.883128)`
- current quarter small-target winner: `(5.173061, 7.353910)`

Add these new untwisted square sizes for each shortlisted point:

- `L = 96, 112, 128, 144`

That is only:

- `3 couplings x 4 new sizes = 12` new untwisted jobs

and it moves the untwisted ladder into the same large-`L` regime as the existing twisted `144 x 144` target.

### Analysis Sequence

1. Direct same-size `144 -> 144` comparison, with no fit at all.

   - If the exact target improves sharply here, that is already strong evidence that long-range extrapolation from `64` was the main problem.

2. Local holdout fit `96,112,128 -> 144`.

   - This is the closest large-target analogue of the current smaller-target logic, where `64 -> 66` is only a tiny endpoint propagation.

3. Stability cross-check using two training windows:

   - local: `96,112,128 -> 144`
   - slightly wider: `64,96,112,128 -> 144`

   If the ranking changes strongly between those two windows, the fit family is still too fragile and the problem has not yet been localized enough.

4. Only if the exact target closes the gap under this local-large-L scheme, widen the follow-up to a small neighborhood around the exact target.

   - Do **not** go back to the full `25`-point grid before the 3-point test is understood.

### Optional Controls After The 12-Job Pass

- add an untwisted `L = 66` run at the exact target, to make the smaller-target comparison a literal same-size control rather than an almost-same-size control
- if the direct `144 -> 144` check still does not help the exact target, then the issue is probably not mainly the long extrapolation; at that point revisit cutoff mismatch, observable weighting, or a more local residual score instead

### Concrete Next Question

The next scientifically useful question is now very narrow:

- does the exact `acute456` target beat the current large-target quarter winner once both are judged with local untwisted data near `L = 144`, instead of by extrapolation from `L = 64`?

## 2026-05-28: Acute456 Quarter-Point Ranking Experiments, Reweighting Note, And Run Decision

### What Was Added After The Large-L Note

Several exploratory ranking/aggregation diagnostics were generated for the `acute456` quarter-point boundary summaries, mainly to test whether a more consensus-style score would stabilize the large-target versus small-target story.

Useful artifacts written during that pass:

- shared-scale quarter `chi2` heatmap comparison:
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_vs_small_target_chi2_heatmaps_shared_scale.png`
- full 7-metric rank aggregation:
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_target_rank_choice_points.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_small_target_rank_choice_points.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_target_rank_choice_squared_points.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_small_target_rank_choice_squared_points.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_rank_choice_points_heatmaps.png`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_rank_choice_squared_points_heatmaps.png`
- ratio-only ranking:
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_target_ratio_z_rank_choice.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_small_target_ratio_z_rank_choice.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_ratio_z_rank_choice_heatmaps.png`
- correlator-only ranking:
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_large_target_corr_z_rank_choice.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_fss_small_target_corr_z_rank_choice.tsv`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_corr_z_rank_choice_heatmaps.png`
  - `responsible_method_tests/standard/acute456_boundary_quarter_fss_sweep/acute456_boundary_quarter_corr_vs_ratio_z_rank_choice.tsv`

### Full 7-Metric Rank Aggregation

The full rank score used the five primitive quarter-point `chi2` panels plus `corr_chi2` and `ratio_chi2`.

Linear rank sum:

- large target winner: `(3.762226, 5.883128)` with rank-point sum `16`
- exact target `(4.702783, 7.353910)`: rank-point sum `60`, rank `6`
- small target winner: `(4.702783, 8.089301)` with rank-point sum `40`
- exact target for the small target: rank-point sum `44`, rank `2`

Squared-rank sum:

- large target winner stays `(3.762226, 5.883128)` with squared-rank sum `54`
- exact target on the large target remains behind with squared-rank sum `672`, rank `6`
- small target winner becomes the exact target `(4.702783, 7.353910)` with squared-rank sum `334`
- runner-up is `(4.702783, 8.089301)` with squared-rank sum `354`

So the exact target is recovered only for the **small** twisted target and only under the harsher mixed squared-rank rule.

### Why That Rank Rule Should Not Be Trusted Too Quickly

This ranking idea is not nonsense; it is basically a Borda/rank-aggregation heuristic with a convex penalty when the ranks are squared.

But it should still be treated as exploratory rather than as the new scientific score because it:

- throws away magnitude information from the underlying `chi2` values
- is sensitive to local rank swaps
- double-counts evidence in the full 7-metric version because `corr_chi2` and `ratio_chi2` are already derived from the primitive panels

So the small-target exact-target win is suggestive, but not yet something to build the next campaign around by itself.

### Subsystem Split: Ratio-Only Versus Correlator-Only

To see whether the mixed-rank result was really coming from one subsystem, the aggregation was rerun separately on the ratio channels and on the correlator channels, ranking by absolute `z` so smaller `|z|` is better.

Ratio-only, using only `(v/4)/(u/4)_z` and `(w/4)/(u/4)_z`:

- large target winner: `(4.232505, 6.618519)` with linear sum `5`, squared sum `13`
- exact target on the large target: linear rank `9`, squared-rank `10`
- small target winner: `(4.702783, 8.089301)` with linear sum `6`, squared sum `18`
- exact target on the small target: linear rank `7`, squared-rank `7`

Correlator-only, using only `v/4_z`, `u/4_z`, and `w/4_z`:

- large target winner: `(3.762226, 5.883128)` with linear sum `6`, squared sum `18`
- exact target on the large target: linear rank `7`, squared-rank `8`
- small target winner: `(5.173061, 7.353910)` with linear sum `14`, squared sum `74`
- exact target on the small target: linear rank `2`, squared-rank `2`

This is the most important interpretive result from the ranking pass:

- ratio-only does **not** recover the exact target
- correlator-only does **not** recover the exact target either
- the mixed squared-rank small-target win is therefore a genuinely hybrid effect, not just a restatement of one subsystem

### Reweighting Note

The idea of reweighting the untwisted branch to the twisted target was considered and then tabled.

Current view:

- same-size geometry reweighting such as `untwisted 144 -> twisted 144` may be viable in principle
- this is **not** a replacement for the current `64 -> 144` size extrapolation problem, because reweighting does not change the lattice size or configuration space
- a real test would require per-configuration or per-block observable data plus enough bond/action information to evaluate `Delta S = S_twisted - S_untwisted` on the same sampled configuration

That is potentially interesting later, but it is not the next cleanest move.

### Current Decision On What To Do Next

The ranking experiments were useful for diagnosis, but they should not drive the next expensive run.

The operational next step should be a cleaner large-`L` signal test.

Recommended first pass:

- dataset: `acute456`
- surface: quarter-point boundary score
- couplings:
  - exact target `(4.702783, 7.353910)`
  - current large-target winner `(3.762226, 5.883128)`
  - current small-target winner `(5.173061, 7.353910)`
- new untwisted sizes: `96, 112, 128, 144`
- initial statistics: start at `n_traj = 40000`

That is still just:

- `3 couplings x 4 sizes = 12` new untwisted jobs

and it is the right compromise between getting a cleaner signal and not overcommitting before the extrapolation hypothesis is actually tested.

### Preferred Analysis Order For That 12-Job Pass

1. direct same-size `144 -> 144` comparison with no fit
2. local holdout fit `96,112,128 -> 144`
3. slightly wider fit `64,96,112,128 -> 144`
4. only if those local-large-L checks help the exact target, widen the coupling neighborhood

### Resume Summary

If resuming later, the main status is now:

- the rank experiments produced interesting but non-decisive diagnostics
- the mixed squared-rank rule rescues the exact target only for the small target, and only heuristically
- ratio-only and correlator-only splits disagree and do not by themselves recover the exact target
- the best next test is **not** more scoring rules; it is the 12-job large-`L`, high-stat untwisted shortlist described above
