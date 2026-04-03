# K_from_TwoPoint

## Overview

Determine the **critical couplings $(k_1, k_2, k_3)$ of an anisotropic
triangular-lattice model** by matching its all-to-all two-point function to
that of an equilateral twisted-boundary reference lattice.

**No analytical critical condition is used at any step.**  The critical
$\beta_c$ for every trial coupling triple is found purely numerically from a
**susceptibility peak scan**.  This makes the method model-independent and
directly extensible to higher dimensions, different lattice actions, or models
where no exact solution is known.  The known 2D Ising duality relation is used
only for post-hoc validation.

### Core physics idea

The CFT at a second-order phase transition is universal.  On a fixed physical
torus, **any** lattice discretisation of the same universality class must
produce the same dimensionless connected two-point function $G_{\rm conn}(m,n)$
at criticality, indexed by lattice displacement $(m,n)$.  We exploit this by:

1. Building a reference $G_{\rm conn}^{\rm ref}(m,n)$ from an equilateral
   lattice at its known critical point.
2. Tuning anisotropic couplings $(k_1, k_2, k_3)$ until the same $G(m,n)$
   is reproduced, **finding $\beta_c$ numerically** at each step.

---

## Physical Setup

### Reference lattice (equilateral)

A twisted-parallelogram triangular lattice with parameters
$(L_x, L_y, T_x, T_y)$ and **isotropic** coupling

$$K_c^{\rm eq} = \tfrac{1}{4}\ln 3 \approx 0.27465307$$

(for the 2D Ising model; for other models, the reference $\beta_c$ would also
be found via susceptibility scan).  The physical torus period vectors are

$$\mathbf{V}_1 = L_x\,\mathbf{e}_1 + T_y\,\mathbf{e}_2, \qquad
\mathbf{V}_2 = T_x\,\mathbf{e}_1 - L_y\,\mathbf{e}_2,$$

where $\mathbf{e}_1 = (1, 0)$ and $\mathbf{e}_2 = (\tfrac{1}{2}, \tfrac{\sqrt{3}}{2})$.

Running `thinkDoubleTwist/ising_tri_twisted_parallelogram` produces a file
`two_point_all_to_all.dat` with columns:

```
# d  m  n  corr  err  corr_conn  err_conn
```

Each row gives the **displaced two-point function** for a single lattice
displacement $(m, n)$, **averaged over all $V$ source sites** using
translation invariance:

$$G(m, n) = \frac{1}{V} \sum_{i=0}^{V-1}
\langle \sigma_i \, \sigma_{i + (m,n)} \rangle.$$

There are exactly $V = N_{\rm cell}$ unique displacements on the torus
(one per site in the fundamental domain), so the file has $V$ data rows.
No $V \times V$ pair matrix is ever stored — only the $V$ displaced
correlators.  The translational average reduces the statistical variance
by a factor of $V$ compared to a single source-site measurement.

### Anisotropic target lattice

The same graph topology $(L_x, L_y, T_x, T_y)$ but with independent couplings
$(k_1, k_2, k_3)$ on the three triangular bond directions.  The critical
surface in the 4-dimensional space $(k_1, k_2, k_3, \beta)$ is **not assumed
known**.  Instead, for each trial $(k_1, k_2, k_3)$, we locate $\beta_c$
numerically.

---

## Finding $\beta_c$ Numerically: Susceptibility Peak

For a given $(k_1, k_2, k_3)$ the critical $\beta$ is found by scanning over
$\beta$ and locating the peak in the magnetic susceptibility:

$$\chi = \frac{1}{V}\left(\langle M^2 \rangle - \langle |M| \rangle^2\right),
\qquad M = \sum_i \sigma_i.$$

On a finite lattice, $\chi(\beta)$ has a rounded peak at $\beta_{\rm peak}(L)$
which converges to $\beta_c$ as $L \to \infty$.  The procedure for a single
$(k_1, k_2, k_3)$ is:

1. **Coarse $\beta$ sweep** — run short simulations at 5–8 $\beta$ values
   bracketing the expected transition.  Identify the $\beta$ interval
   containing the susceptibility maximum.
2. **Refine** — bisect / golden-section on that interval with progressively
   more statistics until $\beta_{\rm peak}$ is determined to the desired
   precision.
3. **Production run** — run a long simulation at $\beta_{\rm peak}$ to
   collect the all-to-all two-point function.

This is the same procedure one would use in 3D or for any model.  The Binder
cumulant $U_4 = 1 - \langle m^4\rangle/(3\langle m^2\rangle^2)$ can
supplement the susceptibility: its crossing between two lattice sizes pins
$\beta_c$ more precisely, but at the cost of needing multiple-$L$ runs.

---

## Matching Condition

At criticality, matching displacement-by-displacement:

$$\chi^2(k_1, k_2, k_3) = \sum_{m,n}
\left(\frac{G_{\rm conn}^{\rm aniso}(m, n;\, \beta_c[k]) \;-\;
            G_{\rm conn}^{\rm ref}(m, n)}
           {\delta G_{\rm combined}(m, n)}\right)^2,$$

where $\beta_c[k]$ is the susceptibility-peak $\beta$ for couplings
$(k_1,k_2,k_3)$, and $\delta G_{\rm combined}$ combines jackknife errors from
both simulations.  The matching is only over displacements with
$r_{\rm min} < r < L/3$ to avoid contact terms and wrap-around noise.

---

## Coupling-Space Optimisation (Adaptive Grid Search)

The outer loop searches coupling space using an **adaptive 4×4 grid**.
This avoids the sensitivity of derivative-free simplex methods to stochastic
noise: every level evaluates 16 independent points, so the minimum is
identified robustly before deciding on the next action.

### Parameterisation

We fix $k_3 = 1$ and optimise over:

$$\mathbf{p} = (k_1/k_3,\; k_2/k_3) \equiv (r_1,\; r_2).$$

The overall scale is absorbed into $\beta$, which is retuned numerically at
each grid point.  This reduces the outer optimisation to **two dimensions**.

### Algorithm

At each level a 4×4 grid of $(r_1, r_2)$ values is evaluated over the square
$[r_{1,c} \pm \Delta,\; r_{2,c} \pm \Delta]$ centred on $(r_{1,c}, r_{2,c})$
with half-span $\Delta$:

1. Run the inner susceptibility scan + production run at each of the
   $4 \times 4 = 16$ grid points.
2. Find the grid point $(r_1^*, r_2^*)$ with the smallest $\chi^2 / N_{\rm dof}$.
3. Decide:
   - **Interior minimum** (best point not on any edge): set
     $(r_{1,c}, r_{2,c}) \leftarrow (r_1^*, r_2^*)$ and halve $\Delta \leftarrow \Delta / 2$.
   - **Border minimum** (best point on an edge): translate the grid centre
     to $(r_1^*, r_2^*)$ keeping $\Delta$ unchanged.
4. Repeat until $\chi^2 / N_{\rm dof} \leq 1$ or the level limit is reached.

The border-translate rule ensures the algorithm will chase the minimum even
if the initial grid is far from the critical surface.  The interior-refine
rule zooms in once the minimum is bracketed.

### Warm start

The susceptibility scan at level $\ell+1$ is warm-started from the best
$\beta_c^{(\ell)}$; if the coupling change is small the peak barely moves,
so the inner loop converges in 2–3 evaluations instead of a full coarse sweep.

### Convergence criterion

Stop when $\chi^2 / N_{\rm dof} \lesssim 1$, meaning the anisotropic
correlator is statistically indistinguishable from the reference.

### Live monitor (heatmap)

A four-panel PNG is updated after every point evaluation:

| Panel | Content |
|-------|---------|
| Top-left | Scatter plot of **all evaluated $(r_1, r_2)$ points** across all grid levels, coloured by $\chi^2/N_{\rm dof}$ (`RdYlGn_r` colourmap); the current grid extent shown as a dashed rectangle; best point marked with a star. |
| Top-right | $\chi^2/N_{\rm dof}$ vs evaluation number (semi-log). |
| Bottom-left | $\beta_c$ vs evaluation number. |
| Bottom-right | Latest susceptibility scan $\chi(\beta)$ **with error bars**; vertical line at $\beta_c$. |

---

## Planned Workflow

### Step 0 — Build the simulator

The existing `thinkDoubleTwist/ising_tri_twisted_parallelogram.cc` supports
anisotropic couplings `--k1 --k2 --k3`.  No new C++ code is needed for data
generation.

```bash
g++ -std=c++14 -O3 -Wall -Wno-sign-compare -I include \
    thinkDoubleTwist/ising_tri_twisted_parallelogram.cc \
    -o K_from_TwoPoint/bin/ising_tri_twisted_parallelogram
```

### Step 1 — Generate the equilateral reference (self-consistent)

The reference lattice uses $k_1 = k_2 = k_3 = 1$ (isotropic, equilateral)
and locates $\beta_c$ **from the susceptibility peak on the same lattice
size**, exactly as the optimizer does for anisotropic couplings.  This
ensures that finite-size shifts in $\beta_{\rm peak}$ cancel between
reference and target.

Default test geometry: $L_x = L_y = 32$, $T_x = T_y = 0$
($N_{\rm cell} = 1024$).  Procedure:

1. **Susceptibility scan** at $k = (1, 1, 1)$ over $\beta \in [0.24, 0.30]$
   to locate $\beta_{\rm peak}$.
2. **Production run** at $\beta_{\rm peak}$ with high statistics.

```bash
# Step 1a: find beta_peak (example values for 32x32)
for beta in 0.24 0.25 0.26 0.265 0.27 0.2747 0.28; do
  K_from_TwoPoint/bin/ising_tri_twisted_parallelogram \
      --L_x 32 --L_y 32 --T_x 0 --T_y 0 \
      --k1 1.0 --k2 1.0 --k3 1.0 --beta $beta \
      --n_traj 100000 --n_skip 20 \
      --data_dir /tmp/ref_scan
done
# Read m_susc from stdout of each run; pick beta with largest chi.

# Step 1b: production run at beta_peak
K_from_TwoPoint/bin/ising_tri_twisted_parallelogram \
    --L_x 32 --L_y 32 --T_x 0 --T_y 0 \
    --k1 1.0 --k2 1.0 --k3 1.0 --beta <beta_peak> \
    --n_traj 500000 --n_skip 20 \
    --data_dir K_from_TwoPoint/ref_equilateral_32x32
```

For the 2D Ising model the known exact $\beta_c = K_c / k = 0.27465$ can be
used as a cross-check, but **it is not required** — the susceptibility scan
is the primary method.  This makes the workflow portable to models without
a known critical point.

### Step 2 — Adaptive grid optimisation

`K_from_TwoPoint/optimise_couplings.py`:

```
Input:  reference two_point_all_to_all.dat
        initial grid centre (r1, r2) = (k1/k3, k2/k3)
        initial half-span Delta
        geometry (Lx, Ly, Tx, Ty)

For each grid level:
    Evaluate 4x4 grid over [r1_c ± Delta, r2_c ± Delta]

    For each grid point (r1, r2):
        k1 = r1,  k2 = r2,  k3 = 1.0

        ┌─ Inner loop: find beta_c ──────────────────────────┐
        │ 1. Coarse sweep: run simulator at 7 beta values    │
        │    measuring chi ± err = <M^2>/V.  Bracket peak.  │
        │ 2. Refine: golden-section on bracket (3 steps).   │
        │ 3. Record beta_c; warm-start next evaluation.     │
        └────────────────────────────────────────────────────┘

        Production run at (k1, k2, k3, beta_c)
        → two_point_all_to_all.dat

        Compute chi2 vs reference
        Update live heatmap PNG

    Find best grid point:
        if interior → refine (Delta /= 2, centre on best)
        if border   → translate (centre on best, keep Delta)

    Break if chi2/ndof <= 1.0

Output: best (r1, r2, beta_c), chi2_ndof, all evaluated points
        fit_result.json, grid_level_NN.json per level
        optimisation_live.png / optimisation_final.png
```

Example invocation:

```bash
python K_from_TwoPoint/optimise_couplings.py \
    --exe K_from_TwoPoint/bin/ising_tri_twisted_parallelogram \
    --Lx 32 --Ly 32 --Tx 0 --Ty 0 \
    --ref K_from_TwoPoint/ref_equilateral_32x32/.../two_point_all_to_all.dat \
    --r1_init 1.2 --r2_init 1.0 --half_span 0.3 \
    --beta_init 0.262 \
    --n_traj_prod 100000 \
    --max_iter 6 \
    --output_dir K_from_TwoPoint/results/aniso_test
```

### Step 3 — Validation against exact solution

For the 2D triangular Ising model the critical surface is known analytically:

$$e^{-2\beta k_1} + e^{-2\beta k_2} + e^{-2\beta k_3} = 1$$

(equivalent to $\sinh(2\beta k_i)\sinh(2\beta k_j)$ relations).  We compare
the numerically determined $(\beta_c, k_1, k_2, k_3)$ against this surface as
a **check only** — it is never used during the optimisation.

Also:
- Run at 2× and 4× lattice size; check finite-size scaling collapse.
- Plot the matched $G_{\rm conn}(m,n)$ side-by-side.

### Step 4 — Apply to other models / 3D

Because the method never invokes an analytical critical condition, the same
code works for:
- **3D Ising** (with `--N_t` time extrusion, susceptibility peak in $\beta$)
- **Potts models** (change the simulator, keep the same outer loop)
- **$\phi^4$ theory** (continuous spins, susceptibility still peaks at $\beta_c$)

Only the simulator binary and the susceptibility observable need to change.

---

## Data Format Conventions

**`two_point_all_to_all.dat`**

```
# d  m  n  corr  err  corr_conn  err_conn
```

- `d` — displacement index ($0, 1, \ldots, V-1$).
- `(m, n)` — integer lattice displacement in the $(e_1, e_2)$ basis.
- `corr` — displaced two-point function
  $G(m,n) = \frac{1}{V}\sum_i \langle\sigma_i \sigma_{i+(m,n)}\rangle$,
  averaged over all $V$ source sites (translation invariance) and over
  jackknife blocks.
- `corr_conn` — connected correlator
  $G_{\rm conn}(m,n) = G(m,n) - \langle m \rangle^2$,
  with jackknife error propagation.
- `err`, `err_conn` — jackknife standard errors.
- The file has exactly $V = N_{\rm cell}$ rows (one per unique displacement
  on the torus).  **No $V \times V$ matrix is stored** — only the $V$
  translation-averaged displaced correlators.

The equilateral physical distance for row $(m, n)$ is:

$$r_{\rm eq}(m, n) = \sqrt{m^2 + mn + n^2}$$

(with $a = 1$).

---

## Key Parameters and Notation

| Symbol | Code flag | Meaning |
|--------|-----------|---------|
| $L_x, L_y$ | `--L_x --L_y` | Parallelogram periods |
| $T_x, T_y$ | `--T_x --T_y` | Twist offsets |
| $k_1$ | `--k1` | Coupling along $+e_1$ |
| $k_2$ | `--k2` | Coupling along $+e_2$ |
| $k_3$ | `--k3` | Coupling along $+(e_2-e_1)$ |
| $\beta$ | `--beta` | Inverse temperature (retuned at each step) |
| $r_1 = k_1/k_3$ | — | Coupling ratio 1 (outer-loop parameter) |
| $r_2 = k_2/k_3$ | — | Coupling ratio 2 (outer-loop parameter) |
| $\chi$ | — | Magnetic susceptibility (used to locate $\beta_c$) |
| $K_c^{\rm eq}$ | — | $\tfrac{1}{4}\ln 3 \approx 0.27465$ (2D Ising, validation only) |
| $\Delta_\sigma$ | — | Ising spin dimension $= \tfrac{1}{8}$ |

---

## Files in this directory

| File | Purpose |
|------|---------|
| `README.md` | This planning document |
| `bin/` | Compiled simulator binary |
| `optimise_couplings.py` | Adaptive 4×4 grid + inner susceptibility-peak driver |
| `plot_match.py` | Side-by-side plot of matched $G(m,n)$ curves |
| `ref_equilateral_32x32/` | Reference data for the 32×32 equilateral torus |
| `results/` | Fit outputs and plots |

---

## Open Questions / Design Decisions

1. **Which displacement pairs to include in $\chi^2$?**  Small separations
   have large correlator values and small statistical errors but are most
   contaminated by lattice artefacts; large separations are noisier.  A
   practical cut of $r_{\rm min} \lesssim r \lesssim L/3$ is a starting point.

2. **Amplitude normalisation.**  Should we normalise each correlator by
   $G(r_{\rm min})$ before differencing, or include an overall amplitude $A$
   as a free fit parameter?  An unconstrained $A$ absorbs finite-size
   corrections but widens the degeneracy direction.

3. **Grid density vs compute cost.**  A 4×4 grid costs 16 inner-loop
   evaluations per level.  A 3×3 grid halves the cost (9 points) but gives
   less information about the landscape and makes border/interior distinction
   harder.  If the inner-loop susceptibility scan dominates, coarser grids
   with more refinement levels may be preferable.  Histogram reweighting
   across nearby $\beta$ values could also reduce inner-loop cost by 2–3×.

4. **Statistical errors on the fit.**  Jackknife on the equilateral reference
   gives $\delta G$ used in $\chi^2$.  For the anisotropic side, fit
   uncertainty can be propagated by repeating the minimisation for each
   jackknife block of the reference data.

5. **Geometry choice for a first test.**  Start with a mildly anisotropic case
   close to equilateral (e.g., $k_1 = 1.2, k_2 = k_3 = 1.0$), where the
   answer is near the equilateral point and convergence should be fast.
   Then validate against the known exact critical surface.

6. **Finite-size effects on $\beta_c$.**  The susceptibility peak drifts with
   $L$.  For the two-point matching we need both lattices at the same effective
   criticality; using the same $L$ and the same peak-finding procedure ensures
   the finite-size shift cancels to leading order.
