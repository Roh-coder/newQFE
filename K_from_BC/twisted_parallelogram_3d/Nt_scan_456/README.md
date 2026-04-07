# Nt Extrusion Scan on the 4-5-6 Triangle Torus (s=1)

This directory contains data and scripts for a study of how the 3D Ising
model on the 4-5-6 triangular torus changes as the temporal extent Nt is
extruded from thin-slab (near-2D) to thick-slab (3D bulk), with the spatial
geometry held fixed at the s=1 base lattice.

## Physics Motivation

The 4-5-6 torus has a unique modular structure: the spatial geometry is fixed
but the temporal direction adds a third dimension.  As Nt grows from 1 toward
infinity:

- At **Nt = 1** the system is a pure 2D Ising model, deep in the **disordered
  phase** at K = 0.161 (below the 2D triangular critical value K_2D ≈ 0.2747).
  The spatial correlator decays exponentially; there is a finite correlation
  length ξ_2D(K).

- As **Nt increases** from 1 the system gains temporal fluctuations; the
  effective 3D coupling strengthens.  Somewhere around Nt ~ O(10) the system
  crosses over to 3D-like behavior.

- At **Nt → ∞** the system approaches the **isotropic 3D Ising critical
  point** at K ≈ 0.161, where the spatial correlator decays as a power-law
  with exponent 2Δ_σ ≈ 1.036.

The scan therefore maps out the full 2D-to-3D crossover on a topologically
non-trivial spatial manifold.

## Lattice Parameters

The spatial geometry is frozen at s=1 throughout:

| Parameter | Value  | Notes                              |
|-----------|--------|------------------------------------|
| Lx        | 39     | 4-5-6 torus base                  |
| Ly        | 48     |                                    |
| Tx        | 9      | twist                              |
| Ty        | −9     | twist                              |
| Ncell     | 1791   | sites per time slice               |
| K_spatial | 0.161  | fixed; sits at the 3D Ising K_c    |
| K_t       | 0.161  | isotropic; same as K_spatial       |

Nt values to scan:

    Nt ∈ { 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32 }

Known reference critical couplings at this spatial geometry:

| Nt | K_c (approx) | Notes                                             |
|----|-------------|---------------------------------------------------|
| 1  | 0.2747      | 2D triangular K_c = (1/4) ln 3                   |
| 2  | ~0.140      | from existing time-scan data in thinkDoubleTwist/ |
| 4  | ~0.161      | approximate 3D isotropic K_c used here            |
| ∞  | ~0.161      | 3D Ising fixed point                              |

**Note:** At K = 0.161, the system is *not* at criticality for small Nt.
For Nt ≤ ~3 it is in the disordered phase; the gap is finite and decreases
as Nt increases. This non-critical approach is the main observable.

---

## Plan

### Phase 1 — Production runs (Nt scan at fixed K = 0.161, s=1)

**Goal**: Measure spatial and temporal two-point functions for each Nt.

**Command template**:

    ising_tri_twisted_parallelogram \
      --L_x=39 --L_y=48 --T_x=9 --T_y=-9 \
      --N_t=<NT> \
      --k1=0.161 --k2=0.161 --k3=0.161 --k_t=0.161 \
      --n_therm=500 --n_traj=5000 --n_skip=10 \
      --n_wolff=3 --n_metropolis=5 \
      --data_dir=<OUTDIR>

**Observables** (from binary output):
1. `two_point_all_to_all.dat` — equal-time spatial G(r) averaged over time slices
2. `two_point_time.dat`       — space-averaged temporal G(Δt)
3. `two_point_typed.dat`      — directional correlators along e1, e2, e2-e1

**Output directories**: `Lx39_Ly48_Tx9_Ty-9_Nt{NT}_k0.161.../`

**Statistics**: 5000 trajectories for all Nt (spatial lattice is always 1791
sites regardless of Nt, so cost is ∝ Nt).  For Nt ≥ 16, reduce to 2000
trajectories if wall time is limiting.

**Run script**: `run_nt_scan.sh`

---

### Phase 2 — Temporal gap vs Nt

**Goal**: Extract the temporal mass gap m(Nt) for each Nt and track how it
approaches zero as Nt → ∞ (i.e., as the system approaches 3D criticality).

**Method**:
- Fit G_time(Δt) = A · cosh[m(Nt − Δt/2)] for each Nt ≥ 4
  (smaller Nt have too few Δt points for a clean fit; use direct ratio instead)
- Extract m_eff(Nt) = arccosh[(G(Δt−1)+G(Δt+1))/(2G(Δt))] at the midpoint
- Plot m(Nt) vs 1/Nt and vs Nt on a log scale

**Expected behavior**:
- For small Nt (~1–3): m(Nt) large and decreasing rapidly — system is off-critical
- For Nt ~ 4–8: m(Nt) falls toward zero as the system approaches 3D K_c
- For large Nt: m(Nt) ∝ 1/Nt^z or saturates if K ≠ K_c exactly

**New script**: `plot_nt_scan_gap.py`
- Reads `two_point_time.dat` for all Nt values
- Plots m(Nt) vs Nt (log-log and linear scales)
- Fits power-law or exponential to locate the effective crossover Nt*
- Output: `gap_vs_nt.png`

---

### Phase 3 — Spatial correlator evolution with Nt

**Goal**: Watch the spatial two-point function G_spatial(r) transform from
an exponentially decaying (massive, 2D-disordered) form at small Nt to a
power-law (critical, 3D) form at large Nt.

**Method**:
- Overlay G_spatial(r) · r^(2Δ) vs r/Lx for all Nt on one plot
  (scaling collapse would occur at the 3D critical Nt)
- For small Nt, fit G(r) ∝ e^{−r/ξ(Nt)} to extract ξ(Nt)
- For large Nt, check approach to power-law with Δ = 0.5182

**References**:
- 2D modular Ising correlator from `thinkDoubleTwist/out_time_scan_456_critical/Nt1/`
  as the Nt=1 baseline
- 3D image-sum CFT prediction from `compare_456_3d_twopoint.py` (Phase 4 of
  the 3D_456_test plan) as the Nt→∞ target

**New script**: `plot_nt_scan_spatial.py`
- 2 panels:
  - Left: G(r) vs r with exponential fits at small Nt
  - Right: G(r)·r^(2Δ) vs r/Lx with power-law overlay at large Nt
- Output: `spatial_corr_vs_nt.png`

---

### Phase 4 — 2D-to-3D crossover characterization

**Goal**: Quantify the crossover scale Nt* at which the system transitions from
2D-like to 3D-like behavior, and compare it to theoretical expectations.

**Observables to track vs Nt**:
1. Temporal mass gap m(Nt) from Phase 2
2. Spatial correlation length ξ(Nt) from Phase 3 (for small Nt)
3. Scaling collapse residual ε(Nt) = deviation of G(r)·r^(2Δ) from flat (for all Nt)
4. Susceptibility χ(Nt) = sum_r G_spatial(r) per time slice

**Crossover criterion**: Nt* is where ξ(Nt) ~ Nt (temporal and spatial scales
become comparable) — the transition from "pancake" to "isotropic" geometry.

**Theoretical expectation**: From dimensional reduction arguments,
Nt* ~ ξ_2D(K) where ξ_2D is the 2D correlation length at K = 0.161.
This can be estimated from the Nt=1 run.

**New script**: `plot_nt_scan_crossover.py`
- 4-panel summary plot:
  - m(Nt) vs Nt
  - ξ(Nt) vs Nt with ξ = Nt line overlaid
  - Scaling collapse quality ε(Nt) vs Nt
  - χ(Nt) vs Nt
- Output: `crossover_vs_nt.png`

---

## Directory Structure (as populated)

    Nt_scan_456/
    ├── README.md                          ← this file
    ├── run_nt_scan.sh                     ← Phase 1 run script (to write)
    ├── plot_nt_scan_gap.py                ← Phase 2 gap vs Nt (to write)
    ├── plot_nt_scan_spatial.py            ← Phase 3 spatial correlator (to write)
    ├── plot_nt_scan_crossover.py          ← Phase 4 crossover summary (to write)
    ├── Lx39_Ly48_Tx9_Ty-9_Nt1_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt2_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt3_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt4_k0.161.../   ← already done in 3D_456_test
    ├── Lx39_Ly48_Tx9_Ty-9_Nt6_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt8_k0.161.../   ← already done in 3D_456_test (s=2 spatial differs)
    ├── Lx39_Ly48_Tx9_Ty-9_Nt10_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt12_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt16_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt20_k0.161.../
    ├── Lx39_Ly48_Tx9_Ty-9_Nt24_k0.161.../
    └── Lx39_Ly48_Tx9_Ty-9_Nt32_k0.161.../

---

## Immediate Next Step

Write `run_nt_scan.sh` and launch Phase 1.  The Nt=4 data from `3D_456_test`
at the same spatial geometry can be reused directly (symlink or copy).
