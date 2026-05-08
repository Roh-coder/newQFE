# Journal

## 2026-05-08 — Continuum Boundary Optimization Smoke Run and Production Plan

### Context
Implemented and executed an end-to-end pipeline for continuum-calibrated boundary-correlator optimization over `(r1, r2)`:
- Script: `run_boundary_continuum_opt.py`
- Output tag: `boundary_opt_smoke`
- Output root: `results/boundary_opt_smoke/`

The pipeline performs:
1. MC precompute for test and reference lattices on size ladders.
2. Boundary-direction sampling at `t_k = k/8`, `k=0..7`, for all 3 cycles.
3. Continuum extrapolation per channel using quadratic fit in `1/L`.
4. Score construction:
   `S(r1,r2) = sum_{c=0..2} sum_{k=0..7} [G_test_inf(c,k)-G_ref_inf(c,k)]^2`
   (uniform weighting in this run).
5. Export of reusable `.dat` files.

---

### What was run
Smoke command (reduced scope for turnaround):
- `r-grid`: `r1,r2 in {0.9, 1.1}` (2x2 points)
- test sizes: `L = [8,16,24]`
- reference mode: `continuum`
- reference sizes: `L = [8,16,24]`
- reference geometries:
  - `L=8  -> (Lx,Ly,Tx,Ty)=(8,8,2,2)`
  - `L=16 -> (16,16,4,4)`
  - `L=24 -> (24,24,6,6)`
- `n_traj=300`, `n_traj_coarse=80`, `n_traj_fine=120`
- workers: `2`

Run status:
- jobs: `24/24 ok`
- errors: `0`
- precompute wall time: `~103.9 s`

---

### Key results
From `results/boundary_opt_smoke/dat/score_map.dat`:
- `(0.9, 0.9) -> 8.726826633`
- `(0.9, 1.1) -> 4.487569639`
- `(1.1, 0.9) -> 0.6815159369`  **minimum**
- `(1.1, 1.1) -> 3.39490116`

Minimum from `results/boundary_opt_smoke/dat/score_minimum.dat`:
- `r1_min = 1.1`, `r2_min = 0.9`, `score_min = 0.6815159369`

Data completeness:
- `n_channels_used = 24` at each point (3 cycles x 8 fractions), confirming full channel coverage.

Continuum fit usage:
- `n_used = 3` points per channel in this smoke run (expected, given 3 sizes).

---

### Interpretation
1. **Pipeline correctness is established** at smoke scale:
   - full precompute-compute-export loop ran successfully,
   - score map and minimum were produced,
   - all requested `.dat` products were emitted.

2. **Physics conclusion is preliminary**:
   - the run is intentionally low-stat and coarse-grid,
   - only 3 sizes were used for quadratic-in-`1/L` fits,
   - therefore uncertainty on channel intercepts is broad.

3. **Directionality signal exists**:
   - the minimum is not flat across the 4 points,
   - suggests structure is already emerging and worth refining with production settings.

---

### Reusable data generated
All reusable text datasets are under:
- `results/boundary_opt_smoke/dat/`

Files:
- `raw_test_channels.dat`
- `raw_ref_channels.dat`
- `continuum_test.dat`
- `continuum_ref.dat`
- `score_map.dat`
- `score_minimum.dat`

Metadata:
- `results/boundary_opt_smoke/manifest_boundary_opt.json`

Visualization:
- `results/boundary_opt_smoke/score_heatmap.png`

---

### Production plan (synthesized)

#### Phase 1 — Baseline production grid
- Test size ladder: `L_test = [8,16,24,32,40]`.
- Keep reference in `continuum` mode.
- Use matching reference size ladder initially (`[8,16,24,32,40]`) with user-defined twists.
- Increase statistics substantially (`n_traj >= 20000` minimum; higher if surface remains noisy).
- Start with moderate `(r1,r2)` mesh to map basin topology.

#### Phase 2 — Basin refinement
- Identify low-score basin from Phase 1.
- Re-run local subgrid around minimum with finer `r-step`.
- Increase `n_traj` only in local region for compute efficiency.

#### Phase 3 — Reference strategy cross-check
- Re-run selected points with `reference_mode=large` using a very large twisted reference lattice.
- Compare minimum location and basin shape against continuum-reference result.
- Use agreement as confidence check for final reported optimum.

#### Phase 4 — Final selection and handoff
- Report final `(r1,r2)` with neighboring score contrast.
- Archive manifest + all `.dat` files as authoritative reusable products.
- Preserve exact CLI command lines used for reproducibility.

---

### Operational notes for next run
- Keep output tags explicit (for example `boundary_opt_prod_v1`, `boundary_opt_prod_refined`).
- Do not overwrite smoke results; treat smoke as baseline benchmark.
- Continue using `.dat` as canonical exported data format.
