# Cost-landscape precompute / replay

Goal: separate the *expensive* part of optimization (MC at each (r1, r2))
from the *cheap* part (cost evaluation given correlators).  Build a dense
test-data cache once, then replay any cost function and visualise its
landscape in seconds.

## Why

Single-seed CMA-ES debugging on a real-MC problem is slow: each evaluation
costs 30–100 s and we need ~30 of them per run to draw any conclusion.
That budget is much better spent on a *one-time* dense grid that all
future cost-function experiments can reuse.

## Architecture

```
                          (one-time, expensive)
   precompute_landscape.py  ──►  results/_landscape/<test_geom>/
       │                              ├── manifest.json
       │                              └── grid/r1_X.XXX_r2_Y.YYY.pkl
       │                                     {beta_c, beta_c_sigma,
       │                                      n_traj_prod, test_data dict}
       ▼
   landscape_cost.py        ──►  results/_landscape/<test_geom>/
       (cheap, replayable)            cost_<mode>_<refTag>.{npz,png}
       │
       └─ for chosen ref_data + cost_mode:
              for each grid point: l2_cost(ref, pt, ...)
              → 2-D array → heatmap PNG
```

## Test-data file format

Each `grid/r1_X.XXX_r2_Y.YYY.pkl` holds:

```
{
  "r1": float, "r2": float, "k3": 1.0,
  "Lx": int, "Ly": int, "Tx": int, "Ty": int,
  "beta_c": float, "beta_c_sigma": float,
  "n_traj_prod": int, "wall_s": float,
  "test_data": {(m, n): {"d": int, "corr": float, "err": float,
                         "conn": float, "conn_err": float}, ...},
}
```

`test_data` is exactly what `mc_engine.load_all_to_all` returns and what
`cost.l2_cost` already consumes — so replay needs zero adapters.

## 4-5-6 problem (initial run)

- test_geom: `(Lx, Ly, Tx, Ty) = (16, 16, 0, 0)`
- ref_geom : `(13, 16, -3, 3)` (existing cache at
  `results/_reference_Lx13_Ly16_Tx-3_Ty3/`)
- truth   : `(r1, r2) = (5.0652, 7.7429)`, `β_c ≈ 0.0628`
- grid    : `r1, r2 ∈ {0.5, 1.0, ..., 10.0}` (step 0.5 → 20×20 = 400 pts)
  for v0 we use step 1.0 → 100 pts to fit in ≈10 min wall on 6 workers
- per-point MC: `n_traj_prod=5000`, `n_traj_scan_coarse=2000`,
  `n_traj_scan_fine=4000` (lighter than CMA evals — we want shape not
  publication-grade error bars).

## Plan

1. **precompute_landscape.py** — multiprocessing pool over the (r1, r2)
   grid.  Each worker: `find_beta_c` (warm-start chained per worker
   using a cell-local seed bracket) then `run_simulator` at β_c, then
   `load_all_to_all`, then pickle.  Skips points whose pkl already
   exists (resumable).  Writes `manifest.json` listing all completed
   points.
2. **landscape_cost.py** — load manifest + reference cache.  Loop over
   grid; for each point call `cost.l2_cost(ref_data, test_data, ...)`
   with the requested `cost_mode`.  Save `(R1, R2, COST, SIGMA)` to
   `.npz` and render a `pcolormesh` heatmap with truth overlay.
3. **evaluate** — visually compare landscapes for `huber_log`,
   `affine_fit`, `test_native`.  Look for: location of global minimum,
   width of basin around truth, presence of the spurious small-(r1, r2)
   trough we observed when starting CMA from (1, 1).

## Cost / time budget

- v0 grid (10×10 = 100 pts) at ≤ 30 s/pt with 6 workers ≈ 8–15 min.
- v1 refinement to 20×20 = 400 pts: ≈ 35–60 min.
- All subsequent cost-function experiments are O(seconds).

## Reuse

Any future cost mode added to `cost.l2_cost` can be heat-mapped against
the same cache without re-running MC.  The cache also serves as a
ground-truth oracle for synthetic-noise stress tests.

## Findings — first 4-5-6 landscape (2026-05-02)

Grid: r1, r2 ∈ {0.5, …, 10.0} step 0.5 (400 points), test=(16,16,0,0),
ref=(13,16,−3,3), n_traj_prod=5000. Truth: (5.07, 7.74). Wall ≈ 60 min
on 6 workers (one-time).

Renders:
- `results/_landscape/Lx16_Ly16_Tx0_Ty0/cost_huber_log_p2_ref-_reference_Lx13_Ly16_Tx-3_Ty3.png`
- `results/_landscape/Lx16_Ly16_Tx0_Ty0/cost_affine_fit_p2_ref-_reference_Lx13_Ly16_Tx-3_Ty3.png`
- `results/_landscape/Lx16_Ly16_Tx0_Ty0/cost_test_native_p2_ref-_reference_Lx13_Ly16_Tx-3_Ty3.png`

**Headline:** all three cost modes place their global minimum at the
bottom-left corner (0.5, 0.5), **not at truth**. Truth sits in a
generic mid-cost region with no resolvable basin. The cost surface is
a bowl tilted toward small geometries, with MC speckle on top.

Per-mode notes:
| mode | argmin | min log10 cost | basin around truth |
|---|---|---|---|
| `huber_log` | (0.5, 0.5) | ≈ −2.5 | none; surface dominated by MC speckle |
| `affine_fit` | (0.5, 0.5) | ≈ −3.8 | none; surface smoother but same tilt |
| `test_native` | (0.5, 0.5) | similar | same |

The right-hand β_c(r1, r2) panel is smooth and physical: β_c grows
monotonically toward the small-(r1, r2) corner, confirming the MC and
β-cascade work correctly. The pathology is purely in the cost
construction, not in the physics.

**Why the corner wins (root cause).** All three costs build a residual
on test sample points whose period vectors scale with (r1, r2). Small
r ⇒ short sample distances ⇒ both G_test and G_ref are near 1 ⇒ log
residuals are trivially small ⇒ cost is small. log G is asymmetric
near 1 (slow growth at short distance, rapid at long), so MC noise
shrinks the cost more on small geometries. Huber clipping caps the
tail but does not cure the tilt.

**Why CMA from (1, 1) failed (now obvious).** The initial cloud sat
inside the artifact basin. CMA correctly found the global min — which
is *not* truth. From (5, 6) the initial cloud straddles truth and the
artifact basin is far away in the ranking, so CMA gets pulled toward
truth before saturating into the corner.

**Implications / fixes to try (cheapest first):**
1. **Lower-bound the search box** at (r_min, r_min) ≈ (3, 3). Excludes
   the artifact corner entirely with no other side effects.
2. **Drop short-distance samples** from the cost: require
   |n1| + |n2| ≥ 2 on the test lattice so the cost cannot be gamed by
   short-distance saturation.
3. **β-ceiling penalty:** add a hard penalty if β_c brackets the upper
   scan bound or never finds an interior peak after 2 sweeps; small
   geometries pin against this in practice.
4. **Subtract leading short-distance behaviour** before residual:
   compare `log G + m·|x|` for a fitted m, so the bowl-tilt is removed.
5. **Geometry-volume penalty:** add `+λ·log(1/(r1·r2))` to the cost so
   the corner is no longer a free lunch.

Once a candidate fix is implemented, re-render the same landscape
(seconds) before running any new CMA. Any cost that does *not* show a
basin at truth on this grid will not work in optimization, period.

## Findings — 14-kernel sweep (2026-05-02)

`landscape_dozen.py` evaluates 14 cost kernels on the same precompute
cache (no MC, ~14 s total wall). Output:
`results/_landscape/Lx16_Ly16_Tx0_Ty0/{dozen.png,dozen.npz}`.

Kernels tested (loss shape × domain × variance weighting × structural):

| family | kernels |
|---|---|
| raw L_p | `l1_diff`, `l2_diff`, `l4_diff` |
| log L_p | `l1_log`, `l2_log` |
| robust | `huber_log` |
| variance-weighted | `chi2`, `chi2_log` |
| relative | `relative` |
| short-distance handling | `drop_short_log`, `effmass`, `slope_loglog` |
| projective | `affine`, `cosine` |

**Result: 14 / 14 kernels argmin at (0.5, 0.5). Distance from argmin
to truth = 8.56 r-units for every single one.** No basin at truth in
any kernel. The bowl tilt toward the small-geometry corner is
structural — it is a property of "compare correlators at sample sites
that scale with the test geometry's r1, r2", and no choice of loss
shape, domain, weighting, or short-distance handling cures it.

**Implication:** the cost design is the problem, not the kernel
choice. The fix must come from one (or a combination) of:

1. **Bounding the search box** away from the corner (CMA-side fix; no
   cost change required, immediately usable).
2. **A geometry-volume regulariser** added to *any* kernel:
   `cost + λ · log(1 / (r1 · r2))` — explicitly penalises the corner.
3. **A β-ceiling penalty** added to *any* kernel: detects pin-against
   upper β-scan bound and adds a fixed cost. Diagnoses + repels the
   corner because small geometries pin β_c.
4. **Reformulation** to use *physical* sample distances (rescale the
   test sample positions back to a reference physical scale before
   comparing) so the cost is no longer a function of r through the
   sample positions.

Options 1, 2, 3 are mechanical and can be added without touching the
kernels. Option 4 is the principled fix but requires reworking the
sample-extraction layer in `cost.py`.

**Workflow rule (now confirmed):** any new cost mode must be replayed
on this landscape grid before any CMA test; if its argmin is anywhere
near (0.5, 0.5) the kernel cannot succeed regardless of CMA tuning.

---

## 8.  Matching-cost replay on cached landscape (2026-05-04)

**Script:** [`landscape_matching_456.py`](landscape_matching_456.py)
(sha256 `408d6afc…`).
**Output dir:** `results/_landscape/_match_456/`.
**Commit:** `d66b5dd`.

### Setup

Replays the **correct matching cost** (algo.md §4) on the existing
cached landscapes (`results/_landscape/Lx{8,16}_Ly{8,16}_Tx0_Ty0/grid/`,
each 400 cells precomputed at `n_traj_prod = 5000`, step = 0.5) on a
sub-sampled **8 × 8 grid** at `r1, r2 ∈ {2, 3, …, 9}` containing the
analytic truth (5.07, 7.74) at d ≈ 0.27 from grid point (5, 8). Two
test sizes (`L = 8, 16`) are processed and a 1/L extrapolation
`cost(L) = c_∞ + a/L` is applied per cell.

The observable is `Z_σ`-invariant:

$$\mathcal{L}(t)=\log|G(t)|-\langle\log|G|\rangle_{t\in 1/8..7/8},$$

paired by **physical side** (algo.md §2.2):
`ref c0 ↔ test c0` (side 5),  `ref c1 ↔ test c2` (side 6),
`ref c2 ↔ test c1` (side 4). Cost = $\sum_{\text{sides}}\sum_k(\Delta\mathcal{L})^2$
over 7 t-samples × 3 sides = 21 residuals.

### Result (figure: `match_cost_grid.png`)

> **Pearson(cost, distance-to-truth):**
> `L=8: ρ = −0.04`,  `L=16: ρ = +0.01`,  `L→∞: ρ = −0.01`.

The cost field is statistically **indistinguishable from noise** with
respect to the truth basin at the current precompute statistics.
Top-5 lowest-cost cells per L are scattered across the plane:
- L=8 argmin (2,4) at d = 4.84 (truth-distance 4.84)
- L=16 argmin (9,8) at d = 3.94
- L→∞ argmin (2,3) at d = 5.65

The five lowest cells per panel include both near-truth points
(e.g. L=16 (5,6) d=1.74, L=8 (4,7) d=1.30) **and** far-from-truth
points; you cannot distinguish them at this signal-to-noise.

### Why it failed (and it is **not** the cost design)

The matching-cost design itself is validated to **0.28 %** on the
twisted↔untwisted reference comparison
([`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)).
What fails here is the **statistics × grid-resolution product** of the
existing precompute:

1. **`n_traj_prod = 5000`** per cell is appropriate for the previous
   landscape diagnostic (where the *shape* — small-r corner divergence
   — dominated by orders of magnitude). The matching cost has **no
   such structural feature**: a correct match yields cost ≈ 0 and
   nearby mismatches yield cost ≈ MC noise. The dynamic range we are
   trying to resolve is comparable to the per-cell statistical error.
2. **Two FSS sizes** (L = 8, 16) is the bare minimum for a 1/L
   extrapolation — it cannot diagnose whether 1/L is even the right
   functional form. 35 / 64 cells of `c_∞_raw` came out **negative**
   and were clipped to 0; that is FSS-noise on cells already at
   matching minimum, not signal.
3. **Step = 1.0** in r is *coarser* than the precompute step (0.5);
   sub-sampling the existing grid wastes the cells closer to truth
   (e.g. (5, 7.5)).

### Recommended production setup

To resolve the truth basin with the matching cost, the precompute
needs:

| knob              | current | recommended (production)               |
|-------------------|---------|----------------------------------------|
| `n_traj_prod`     | 5 000   | **≥ 50 000** (10× ⇒ noise ÷ √10 ≈ 3×)  |
| FSS sizes         | 2       | **≥ 4** (e.g. L ∈ {8, 16, 24, 32})     |
| grid step         | 0.5     | 0.5 near truth, 1.0 in periphery       |
| obs               | log-mean ✓ | log-mean (Z_σ-invariant) ✓          |
| sample density    | 7 t-samples × 3 sides = 21 | retain |
| cost              | un-weighted `Σ(Δ𝓛)²` ✓ | un-weighted ✓                  |

A back-of-envelope: matching cost target sensitivity is
`Δcost ~ 21 · (Δ𝓛)²`, with `Δ𝓛 ~ 1/√n_traj`. At `n_traj = 5 000`
the noise floor is `cost ≳ 21 · (0.014)² ≈ 4 × 10⁻³`, which is
exactly the median cost we observe (median 0.012-0.014). At
`n_traj = 50 000` the floor drops to `~ 4 × 10⁻⁴`, which is also
roughly the truth-cell value at L = 16 (0.002) — i.e. the basin
becomes resolvable with one additional decade of stats.

### Conclusion

**The matching cost is the right observable but the existing
diagnostic landscape is not a fit-for-purpose precompute for it.** The
two artefacts (cost design ↔ corner-divergence test on it; matching
cost ↔ truth-basin test on it) require different statistics tiers.
The cached landscape was built for the former; the latter needs a
**dedicated higher-stats precompute** before the matching cost can
seed a CMA-ES production run.

The script and its negative-result figures are kept under
`results/_landscape/_match_456/` so future work has the noise-floor
calibration to compare against.
