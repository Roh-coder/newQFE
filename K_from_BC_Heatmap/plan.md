# Continuum Boundary-Correlator Optimization Plan

## 1. Objective
Build a production pipeline that scans coupling-ratio space `(r1, r2)` and computes a continuum-level mismatch score between:
- a test lattice family (untwisted, multiple sizes), and
- a user-defined twisted reference lattice family (user controls size and twist per level).

The output is a smooth score manifold over `(r1, r2)` with a clear minimum for geometry-to-coupling calibration.

## 2. Parameter-Space and Geometry Definition

### 2.1 Coupling scan
- Scan a rectangular grid in `(r1, r2)`.
- Use fixed `k3` (default `1.0`) and set:
	- `k1 = r1`
	- `k2 = r2`

### 2.2 Test lattices (fixed production set)
- Test sizes:
	- `L_test = [8, 16, 24, 32, 40]`
- Test geometry default:
	- untwisted square: `(Lx, Ly, Tx, Ty) = (L, L, 0, 0)`

### 2.3 Reference lattices (user-configurable)
Allow user input for a list of reference geometry levels:
- `(Lx_ref_i, Ly_ref_i, Tx_ref_i, Ty_ref_i)` for each level `i`
- Optional independent reference size list `L_ref`

Two operation modes:
1. **Reference continuum mode**: run multiple reference sizes and extrapolate to continuum.
2. **Reference large-L mode**: run a single very large twisted reference lattice and use it as a near-continuum proxy.

## 3. Observable Construction

### 3.1 Boundary-direction correlators
For each lattice realization, sample connected correlators along all 3 boundary cycles at 8 fractional positions:
- `t_k = k/8`, `k = 0..7`

Channel index set:
- cycle index `c in {0,1,2}`
- fraction index `k in {0,...,7}`

Stored channels:
- `G_test(c, k, L)`
- `G_ref(c, k, L)` (or `G_ref_large(c, k)` in large-L mode)

### 3.2 Critical-point policy
At each lattice size and `(r1, r2)` point:
- locate pseudo-critical `beta_c(L)` via the existing Gram-Charlier scan, then
- run production MC at `beta_c(L)` to measure correlators.

## 4. Continuum Extrapolation

For each channel `(c, k)`, extrapolate in `x = 1/L` using quadratic form:

`G(L) = a + b x + c x^2`

Continuum value is the intercept:
- `G_inf = a`

Apply this to:
- test channels: `G_test_inf(c,k)`
- reference channels (continuum mode): `G_ref_inf(c,k)`

In reference large-L mode:
- set `G_ref_inf(c,k) := G_ref_large(c,k)`

## 5. Optimization Score Definition

Residual per channel:

`R(c,k) = G_test_inf(c,k) - G_ref_inf(c,k)`

Point score at each `(r1, r2)`:

`S(r1,r2) = sum_{c=0..2} sum_{k=0..7} w(c,k) * R(c,k)^2`

Default:
- uniform weights `w(c,k)=1`

Optional production weighting:
- inverse-variance weights `w(c,k)=1/sigma(c,k)^2` when propagated continuum uncertainties are available.

## 6. Production Workflow

### Phase A: Precompute data grid
For every `(r1, r2)` point:
1. Run all test sizes in `L_test`.
2. Run reference set (either multi-size continuum mode or one-shot large-L mode).
3. Cache all raw outputs in a resumable directory layout.

### Phase B: Build continuum datasets
1. For each channel `(c,k)`, fit test vs `1/L` quadratically and store `G_test_inf(c,k)`.
2. Build reference continuum channels similarly, or ingest large-L proxy channels.
3. Save continuum tensors and uncertainty fields.

### Phase C: Compute score manifold
1. Evaluate `S(r1,r2)` on the full grid.
2. Save `score_map` and minimum location.
3. Generate heatmaps and residual-diagnostic plots around low-score basin.

### Phase D: Refinement near best basin
1. Locally refine `(r1, r2)` grid near minimum.
2. Increase MC stats only in refined region.
3. Recompute and report final optimum with confidence intervals.

## 7. Data Products
For each production tag `results/<tag>/`:
- `manifest.json`: full run configuration
- `grid/...`: cached per-size per-geometry payloads
- `continuum_test.npz`: `G_test_inf(c,k,r1,r2)` and uncertainty
- `continuum_ref.npz`: `G_ref_inf(c,k,r1,r2)` and uncertainty (or large-L snapshots)
- `score_map.npz`: `S(r1,r2)` and metadata
- `plots/score_heatmap.png`: score manifold with marked minimum
- `plots/residual_channels_*.png`: channel-level residual structure

## 8. Operational Notes
- Keep caching and resumability as first-class behavior to support long runs.
- Maintain fixed MC policy across the scan to avoid introducing procedural bias into the manifold.
- Prefer reference continuum mode for highest fidelity; use large-L reference mode when budget is constrained.
- Use local refinement around the coarse minimum to improve final parameter precision at modest additional cost.

## 9. Success Criteria
The production run is complete when:
1. a stable minimum in `S(r1,r2)` is identified,
2. local refinement preserves the same basin,
3. final best-fit `(r1, r2)` and uncertainty are reported with full run metadata.
