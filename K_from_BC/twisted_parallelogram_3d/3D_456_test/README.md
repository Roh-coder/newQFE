# 3D Ising on the 4-5-6 Triangle Torus

This directory contains data and scripts for a 3D finite-size scaling study
of the Ising model on the 4-5-6 triangular torus, extruded in the time direction.

## Background

The 2D (Nt=1) baseline data lives in:

    thinkDoubleTwist/out_time_scan_456_critical/Lx39_Ly48_Tx9_Ty-9_Nt1_k0.260.../

The spatial geometry is the 4-5-6 triangle torus:

| Parameter | Base value (s=1) |
|-----------|-----------------|
| Lx        | 39              |
| Ly        | 48              |
| Tx        | 9               |
| Ty        | -9              |
| Ncell     | 1791            |

The binary `K_from_BC/twisted_parallelogram_3d/ising_tri_twisted_parallelogram`
already supports 3D via `--N_t` and `--k_t`. The standard scaling series uses
integer scale factors s = 1, 2, 3, ... applied to the full geometry:

    (Lx, Ly, Tx, Ty, Nt) = s × (39, 48, 9, -9, Nt0)

with Nt0 = 4 so that the temporal extent matches the isotropic 3D critical
coupling K = 0.161 (used as a fixed approximation throughout).

Known approximate critical couplings (isotropic k1=k2=k3=kt=K):

| Nt  | K (approx) | Notes                             |
|-----|-----------|-----------------------------------|
| 1   | 0.2747    | 2D triangular Kc = (1/4)ln(3)    |
| 2   | 0.140     | from out_time_scan_456_critical   |
| 4   | 0.161     | approximate value used for production runs |

---

## Plan

### Phase 1 — Production correlator runs

**Goal**: High-statistics spatial and temporal two-point functions at K = 0.161.

**Method**: Direct calls to the existing binary at each scale factor s:

    ising_tri_twisted_parallelogram \
      --L_x=39s --L_y=48s --T_x=9s --T_y=-9s \
      --N_t=4s --k1=0.161 --k2=0.161 --k3=0.161 --k_t=0.161 \
      --n_therm=500 --n_traj=2000 --n_skip=10

**Observables collected per run**:
1. **Equal-time spatial two-point function** — `two_point_all_to_all.dat`:  
   G_spatial(r) = <sigma(x,t) sigma(0,t)>_c averaged over all t slices.
2. **Space-averaged temporal correlator** — `two_point_time.dat`:  
   G_time(tau) = (1/Ncell) sum_x <sigma(x, tau) sigma(x, 0)>_c, averaged over the spatial lattice.

**Output directories**: `out_456_3d/Lx{39s}_Ly{48s}_Tx{9s}_Ty{-9s}_Nt{4s}_k0.161.../`

**Statistics**: ~1000–5000 trajectories (larger s = fewer traj due to cost).
Suggested:  s=1 → 5000, s=2 → 2000, s=3 → 1000.

### Phase 2 — Two-point function analysis

**Goal**: Compare spatial and temporal connected correlators across scales s
at fixed K=0.161, Nt/Lx = 4, and verify approach to 3D Ising scaling.

**2a — Equal-time spatial correlator**

**Existing tool**: `thinkDoubleTwist/plot_spatial_two_point_vs_nt.py`
- Reads `two_point_all_to_all.dat`; overlays by Nt value.

**New tool needed**: `plot_456_3d_scaling.py`
- Overlay G_spatial(r) × r^{2*Delta_sigma} vs r/L for each s (scaling collapse plot)
- Use 3D Ising spin scaling dimension: Delta_sigma ≈ 0.5182
- Expected: curves at different s collapse onto a single master curve at K_c

**2b — Space-averaged temporal correlator**

**New tool needed**: `plot_456_3d_time_corr.py`
- Plot G_time(tau) vs tau/Nt for each s
- Fit to cosh form to extract effective mass m_eff(s) and check approach to
  the continuum dispersion relation
- Overlay temporal and spatial effective masses as consistency check

### Phase 3 — Comparison to 2D analytical modular Ising solution

**Goal**: Compare the equal-time spatial correlators from the 3D MC runs to
the exact two-point function from the analytical modular Ising solution on the
same 4-5-6 torus geometry in 2D, as a baseline and cross-check.

**Background**: The 2D modular Ising CFT on the 4-5-6 triangle torus has an
analytical solution for the spin-spin correlator ⟨σ(x) σ(0)⟩ that encodes
the full modular structure of the boundary conditions (twist sectors Tx=9,
Ty=-9). The 3D equal-time spatial correlator at large Nt should approach this
2D result in the limit where the temporal direction decouples.

**New tool needed**: `compare_456_3d_vs_2d_modular.py`
- Read the existing 2D modular Ising correlator data (or compute analytically
  from `ModularIsingFEM/`) for the 4-5-6 torus boundary conditions
- Overlay with the equal-time spatial G_spatial(r) from Phase 1 MC runs at
  each scale s, normalized appropriately
- Also compare to the Nt=1 baseline data in
  `thinkDoubleTwist/out_time_scan_456_critical/Lx39_Ly48_Tx9_Ty-9_Nt1_k0.260.../`
- Quantify the deviation of the 3D correlator from the 2D analytical result
  as a function of Nt/L (finite temporal-extent effect)

### Phase 4 — Comparison to 3D Ising CFT prediction

**Goal**: The two-point function on the 4-5-6 torus inherits non-trivial
modular/geometric structure. Compare measured G(r) to the CFT prediction on
this specific torus geometry.

**New tool needed**: `compare_456_3d_twopoint.py`
- Compute the free-field / conformal block prediction for <sigma sigma> on
  the 4-5-6 torus at the 3D Ising fixed point
- Overlay with MC data from Phase 2
- This is novel physics — no existing comparison script

---

## Directory structure (as populated)

    out_456_3d/
    ├── README.md                        ← this file
    ├── plot_456_3d_scaling.py              ← Phase 2a scaling collapse (to write)
    ├── plot_456_3d_time_corr.py            ← Phase 2b temporal correlator (to write)
    ├── compare_456_3d_vs_2d_modular.py     ← Phase 3 vs 2D analytical modular Ising (to write)
    ├── compare_456_3d_twopoint.py          ← Phase 4 CFT comparison (to write)
    ├── Lx39_Ly48_Tx9_Ty-9_Nt4_k0.161.../   ← s=1 production run
    ├── Lx78_Ly96_Tx18_Ty-18_Nt8_k0.161.../ ← s=2 production run
    └── Lx117_Ly144_Tx27_Ty-27_Nt12_k0.161./ ← s=3 production run

---

## Immediate next step

Write a production run script for Phase 1 (s=1,2,3 at K=0.161), using the
already-built binary at `K_from_BC/twisted_parallelogram_3d/ising_tri_twisted_parallelogram`.
