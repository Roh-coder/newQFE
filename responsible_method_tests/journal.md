# Journal

## 2026-05-21: Twisted-Reference 4x Trial

### Goal

Test the user hypothesis that the remaining twisted-versus-untwisted mismatch is mainly caused by the two continuum manifolds being sampled on different lattice point sets. The rerun changes the comparison workflow to mimic the intended coupling-to-geometry pipeline:

- build a denser twisted reference family
- fit the twisted continuum manifold first
- interpolate that twisted continuum manifold onto the untwisted continuum points
- use only the twisted-on-untwisted residuals to decide distinguishability

### Configuration

- Run config: `configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json`
- Comparison mode: `python scripts/align_manifolds.py --comparison-mode twisted_reference`
- Twisted family rule: base cell chosen to be the nearest integer multiple of the original pickGeo cell whose area is about `4 x` the untwisted `12x12` base, then scaled by `{1,2,3,4,5}`
- Untwisted family rule: keep the original size ladder `(12,24,36,48,60)`

### Run Status

- Completed benchmark: `geometry_acute_111`
- Completed comparison artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- Completed summary artifact: `results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/SMOKE_SUMMARY.md`
- Interrupted benchmark: `geometry_acute_122`
- Interruption point: twisted scale 5, lattice `(160,100,40,20)`
- Interruption mode: manual `KeyboardInterrupt`, so no trustworthy campaign-wide conclusion exists yet for the 10-geometry set

### Measured Result For `geometry_acute_111`

This is the cleanest control point because both branches are isotropic and target the same equilateral torus.

| quantity | previous symmetric FSS run | 4x twisted-reference run |
| --- | ---: | ---: |
| run tag | `raw_manifold_fss_pickgeo10_fss_20260521` | `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` |
| comparison mode | `symmetric` | `twisted_reference` |
| twisted comparable points | `143` | `575` |
| untwisted comparable points | `143` | `143` |
| twisted alpha | `0.3047` | `0.4128` |
| twisted chi2/dof | `0.1843` | `0.2111` |
| untwisted alpha | `0.2863` | `0.4573` |
| untwisted chi2/dof | `0.0505` | `0.2859` |
| twisted-on-untwisted RMS z | `2.2363` | `3.2160` |
| untwisted-on-twisted RMS z | `2.1648` | `3.9717` |
| common-grid relative RMS | `0.0696` | `0.0907` |
| verdict | `distinguishable` | `distinguishable` |

### Interpretation

- The individual twisted and untwisted continuum fits still look healthy. Both chi2/dof values remain comfortably below `1` in the 4x trial.
- The asymmetric reference-style comparison does **not** reduce the disagreement for the equilateral control. The primary metric, twisted-on-untwisted RMS z, rises from `2.2363` to `3.2160`.
- The denser twisted family increases the twisted continuum point count from `143` to `575`, so this is not a case where the result is being driven by a lack of reference points.
- For `geometry_acute_111`, the hypothesis "the manifolds only look different because the point sets do not sit on top of each other" is therefore disfavored. Even when the twisted manifold is treated as the dense reference and sampled only at the untwisted continuum points, the mismatch stays well above the `Z = 2` marginal threshold.

### Operational Note

The 4x twisted-reference workflow is substantially more expensive than the previous symmetric run.

- `geometry_acute_111` twisted scale 5 already reached lattice `(120,120,0,0)` and took about `51 s` for one production job.
- `geometry_acute_122` was interrupted during twisted scale 5 at `(160,100,40,20)`.
- At `n_workers = 1`, the full 10-geometry campaign is practical but materially slower than the earlier FSS run.

### Current Scientific Takeaway

From the one completed control benchmark, denser twisted sampling by itself does not rescue indistinguishability. If the broader 10-geometry campaign is resumed, it should be treated as a stronger test of the continuum-equivalence claim, not as a formality expected to automatically collapse the manifolds.

### Next Useful Step

- Resume `raw_manifold_fss_pickgeo10_twisted_reference4x_20260521` from the current partial tree, ideally benchmark-by-benchmark, so runtime is easier to control.
- After at least a few additional geometries complete, regenerate the comparison plots with `--comparison-mode twisted_reference` and build a dedicated campaign summary for that run.

### Resume Later

Current saved state before stopping:

- run root: `responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/`
- completed benchmark manifest: `manifest_geometry_acute_111.json`
- completed comparison artifact: `geometry_acute_111_twisted_vs_untwisted_common_grid.json`
- next blocked slice when interrupted: `geometry_acute_122`, twisted scale 5, lattice `(160,100,40,20)`
- config already has `"resume": true`, so rerunning the same generation command should reuse the completed payloads under the current run root

Resume command:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/generate_pointwise_manifold_dataset.py \
	--config responsible_method_tests/configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json
```

After enough benchmarks finish, regenerate the asymmetric comparison artifacts with:

```bash
source /workspaces/newQFE/.venv/bin/activate
for manifest in responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_*.json; do
	python responsible_method_tests/scripts/align_manifolds.py \
		--benchmark-manifest "$manifest" \
		--comparison-mode twisted_reference
done
```

Then build the campaign summary with:

```bash
source /workspaces/newQFE/.venv/bin/activate
python responsible_method_tests/scripts/summarize_benchmark.py \
	--comparison-mode twisted_reference \
	--output responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/PICKGEO10_TWISTED_REFERENCE4X_SUMMARY.md \
	--benchmark-manifest responsible_method_tests/results/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521/manifest_geometry_acute_111.json
```

and add the remaining `--benchmark-manifest` arguments once those benchmark manifests exist.