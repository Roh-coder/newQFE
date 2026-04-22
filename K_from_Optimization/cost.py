#!/usr/bin/env python3
"""
cost.py — Cost functions for the curve-collapse matching problem.

The goal is to find anisotropic couplings (k1, k2, k3) on an untwisted test
lattice whose two-point correlator G_test(t) matches the reference correlator
G_ref(t) along each of the three fundamental torus directions.

Base cost — plain L²
--------------------
For each boundary direction d ∈ {v, u, w}, the correlator is sampled along a
straight path parameterised by t ∈ [0, 1] (arc-length fraction).  The base
cost is the integrated squared difference summed over directions:

    C_d   = ∫₀¹ [G_test^d(t) - G_ref^d(t)]² dt
    C(r1,r2) = Σ_d C_d               ("l2_sum", current default)

This is the plain L² norm of the residual curves.  It has units of [G]² and
goes to zero when the curves collapse onto each other.

Anisotropy-penalising variants
-------------------------------
The plain sum allows the optimizer to zero one direction while ignoring the
others (trading off anisotropic mismatch for a lower total).  Four variants
penalise this:

  A. Variance-normalised L² (weighted sum)
       C = Σ_d  C_d / ‹σ_d²›
     Directions where the signal is tight (low noise) count more; noisy
     directions are down-weighted automatically.  No free parameter.

  B. L⁴ (biquadratic / quartic) — current large-deviation penaliser
       C = Σ_d  ∫₀¹ [G_test - G_ref]⁴ dt
     Raises each per-point deviation to the fourth power before integrating.
     One badly-matching point hurts enormously more than in L².  Gradient is
     steeper near the optimum so convergence is faster once nearby.
     Propagated uncertainty:
       σ_C² ≈ Σ_d  ∫₀¹ (4 Δ³)² (σ_ref² + σ_test²) dt

  C. Max-direction + regularisation (L∞-inspired)
       C = max(C_v, C_u, C_w) + λ·(C_v + C_u + C_w),    λ ≈ 0.1
     Forces the optimizer to reduce the worst direction first.  λ keeps the
     surface smooth and avoids flat gradients when two directions tie.

  D. L² + anisotropy penalty (sum of squared per-direction costs)
       C = Σ_d C_d  +  μ · Σ_d C_d²,    μ ≈ 1 / ‹C_d›
     The second term is zero only when all C_d are equal; it grows when one
     direction dominates.  Cheapest to compute — just add one line to the
     existing summation.  Recommended first experiment.

Uncertainty in C
----------------
Both G_test and G_ref are MC averages with per-point standard errors σ_test(t)
and σ_ref(t).  By linear error propagation, the variance of the integrand at
each t is:

    Var[(G_test - G_ref)²] ≈ 4 (G_test - G_ref)² (σ_test² + σ_ref²)

so the variance of C is approximately:

    σ_C² ≈ Σ_d  ∫₀¹ 4 (G_test - G_ref)² (σ_test² + σ_ref²) dt

This is also computed and returned so the caller can check the SNR = C / σ_C.

SNR criterion
-------------
SNR < 1  → C is consistent with zero within one standard deviation.  The
           curves are already matching within the current statistical precision.
           Increasing n_traj (ideally doubling it) is required before the
           optimizer can make further meaningful progress.

1 ≤ SNR < 3 → Minimum is marginally resolved.  One more refinement step may
               be possible but convergence claims are uncertain.

SNR ≥ 3  → Minimum is clearly distinguished from its neighbours.  Safe to
            continue refining.

Public API
----------
  boundary_paths(Lx, Ly, Tx, Ty) -> list of three (dm, dn) direction vectors
  l2_cost(ref_data, test_data, ...) -> (cost, sigma_cost, per_dir_costs)
  snr(cost, sigma_cost) -> float
  snr_status(cost, sigma_cost) -> str  ("ok" | "marginal" | "need_more_stats")
"""

import math

import numpy as np
from scipy.interpolate import LinearNDInterpolator


# ---------------------------------------------------------------------------
# Torus geometry helpers
# ---------------------------------------------------------------------------

def boundary_paths(Lx, Ly, Tx, Ty):
    """
    Return the three fundamental boundary direction vectors of the torus.

    The torus is spanned by periodicity vectors
        v = (Lx, Ty)   in triangular-lattice (m, n) coordinates
        u = (Tx, -Ly)

    The three boundary directions used for correlator matching are:
        v  = (Lx, Ty)           — first periodicity vector
        u  = (Tx, -Ly)          — second periodicity vector
        w  = -(v + u)           — diagonal (closes the triangle)

    Returns a list of three (dm, dn) tuples.
    """
    return [
        (Lx,        Ty),
        (Tx,        -Ly),
        (-Lx - Tx,  Ly - Ty),
    ]


# ---------------------------------------------------------------------------
# Interpolation helpers (torus-tiling, same as mc_engine.py)
# ---------------------------------------------------------------------------

_SQRT3_2 = math.sqrt(3.0) / 2.0


def _triangular_xy(m, n):
    """Convert triangular-lattice (m, n) → Cartesian (x, y)."""
    return m + 0.5 * n, _SQRT3_2 * n


def _tile_interp(data, Lx, Ly, Tx, Ty, field="conn", copies=2):
    """
    Build a LinearNDInterpolator for *field* by tiling the all-to-all data
    over neighbouring copies of the torus.  This allows clean interpolation
    anywhere along a boundary path without hitting the fundamental-domain
    boundary.
    """
    keys = list(data.keys())
    m_arr = np.array([k[0] for k in keys], dtype=float)
    n_arr = np.array([k[1] for k in keys], dtype=float)
    c_arr = np.array([data[k][field] for k in keys], dtype=float)

    x0, y0 = _triangular_xy(m_arr, n_arr)

    # Periodicity vectors in Cartesian coordinates
    vx = Lx + 0.5 * Ty;  vy = _SQRT3_2 * Ty
    ux = Tx - 0.5 * Ly;  uy = -_SQRT3_2 * Ly

    x_all, y_all, c_all = [x0], [y0], [c_arr]
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            if a == 0 and b == 0:
                continue
            x_all.append(x0 + a * vx + b * ux)
            y_all.append(y0 + a * vy + b * uy)
            c_all.append(c_arr)

    x_all = np.concatenate(x_all)
    y_all = np.concatenate(y_all)
    c_all = np.concatenate(c_all)

    pts = np.column_stack([x_all, y_all])
    _, uid = np.unique(np.round(pts, 8), axis=0, return_index=True)
    return LinearNDInterpolator(pts[uid], c_all[uid])


# ---------------------------------------------------------------------------
# Core cost function
# ---------------------------------------------------------------------------

def l2_cost(ref_data, test_data,
            test_Lx, test_Ly, test_Tx, test_Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            n_samples=400, copies=2, mu=0.0):
    """
    Compute the plain L² cost between the reference and test correlators
    along the three boundary directions.

    Parameters
    ----------
    ref_data : dict
        All-to-all correlator data for the reference lattice, keyed by (m, n).
        Each value must have "conn" and "conn_err" fields.
    test_data : dict
        All-to-all correlator data for the test lattice, same format.
    test_Lx, test_Ly, test_Tx, test_Ty : int
        Geometry of the test lattice.  For the untwisted case Tx=Ty=0.
    ref_Lx, ref_Ly, ref_Tx, ref_Ty : int
        Geometry of the reference lattice.
    n_samples : int
        Number of t ∈ [0, 1] points used for the trapezoid-rule integration.
    copies : int
        Number of torus copies to use when building the interpolants.  2 is
        sufficient for all standard torus geometries.
    mu : float
        Anisotropy-penalty weight (option D).  Set mu > 0 to add the term
        μ · Σ_d C_d² to the cost.  This is zero when all per-direction
        costs are equal and grows when one direction dominates.  A reasonable
        starting value is mu = 1 / mean(C_d) so the penalty is O(1) at the
        initial point.  Default 0.0 reproduces the plain L² cost.

    Returns
    -------
    cost : float
        Total L² cost C = Σ_d ∫ (G_test - G_ref)² dt  (summed over 3 directions).
    sigma_cost : float
        Propagated one-sigma uncertainty in cost from MC errors.
    per_dir : list of float
        Per-direction contributions [C_v, C_u, C_w].
    per_dir_sigma : list of float
        Per-direction uncertainty σ_C for each direction.
    """
    # Build interpolants for the connected correlator and its error
    iref      = _tile_interp(ref_data,  ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty,  "conn",     copies)
    itest     = _tile_interp(test_data, test_Lx, test_Ly, test_Tx, test_Ty, "conn",     copies)
    iref_err  = _tile_interp(ref_data,  ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty,  "conn_err", copies)
    itest_err = _tile_interp(test_data, test_Lx, test_Ly, test_Tx, test_Ty, "conn_err", copies)

    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    t = np.linspace(0.0, 1.0, n_samples)

    per_dir       = []
    per_dir_sigma = []
    dir_labels    = ["v", "u", "w"]

    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        # Cartesian endpoints of each path
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        tex, tey = tdm + 0.5 * tdn, _SQRT3_2 * tdn

        pts_ref  = np.column_stack([t * rex, t * rey])
        pts_test = np.column_stack([t * tex, t * tey])

        G_ref  = np.asarray(iref(pts_ref),       dtype=float)
        G_test = np.asarray(itest(pts_test),      dtype=float)
        e_ref  = np.abs(np.asarray(iref_err(pts_ref),  dtype=float))
        e_test = np.abs(np.asarray(itest_err(pts_test), dtype=float))

        mask = (np.isfinite(G_ref) & np.isfinite(G_test)
                & np.isfinite(e_ref) & np.isfinite(e_test))

        if mask.sum() < 4:
            # Not enough valid interpolation points along this path;
            # skip the direction and contribute 0 (conservative).
            per_dir.append(0.0)
            per_dir_sigma.append(0.0)
            continue

        diff     = G_test[mask] - G_ref[mask]
        var_prop = e_ref[mask]**2 + e_test[mask]**2   # variance of G_test - G_ref
        t_m      = t[mask]
        dt       = t_m[1:] - t_m[:-1]

        # L² cost: trapezoid integral of diff²
        integrand_c = diff ** 2
        C_d = float(np.sum(0.5 * (integrand_c[:-1] + integrand_c[1:]) * dt))

        # Propagated variance: Var(diff²) ≈ (2 diff)² · Var(diff) = 4 diff² · (σ_ref² + σ_test²)
        integrand_v = 4.0 * diff**2 * var_prop
        V_d = float(np.sum(0.5 * (integrand_v[:-1] + integrand_v[1:]) * dt))

        per_dir.append(C_d)
        per_dir_sigma.append(math.sqrt(max(V_d, 0.0)))

    # Total cost and its uncertainty (add in quadrature across directions)
    cost       = sum(per_dir)
    sigma_cost = math.sqrt(sum(s**2 for s in per_dir_sigma))

    # Option D: anisotropy penalty  C += μ · Σ_d C_d²
    # Propagated variance adds  μ² · Σ_d (2 C_d · σ_{C_d})²
    if mu != 0.0:
        penalty = mu * sum(c**2 for c in per_dir)
        var_penalty = mu**2 * sum((2.0 * c * s)**2
                                  for c, s in zip(per_dir, per_dir_sigma))
        cost       += penalty
        sigma_cost  = math.sqrt(sigma_cost**2 + var_penalty)

    # Print per-direction breakdown for diagnostics
    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(dir_labels, per_dir, per_dir_sigma))
    label = f"[L²+μΣC²(μ={mu:.3g})]" if mu != 0.0 else "[L²]"
    print(f"    {label} per-dir: {parts}  total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


# ---------------------------------------------------------------------------
# SNR helpers
# ---------------------------------------------------------------------------

def snr(cost, sigma_cost):
    """
    Return the signal-to-noise ratio C / σ_C.

    This measures how many standard deviations the cost is above zero.
    A value < 1 means the curves are already consistent with matching
    within the current statistical precision.
    """
    if sigma_cost <= 0.0:
        return float("inf")
    return cost / sigma_cost


def snr_status(cost, sigma_cost):
    """
    Classify the current statistics adequacy based on SNR = C / σ_C.

    Returns one of:
        "need_more_stats"  — SNR < 1: cost consistent with zero; increase n_traj
        "marginal"         — 1 ≤ SNR < 3: proceed cautiously
        "ok"               — SNR ≥ 3: optimizer can safely refine further
    """
    s = snr(cost, sigma_cost)
    if s < 1.0:
        return "need_more_stats"
    if s < 3.0:
        return "marginal"
    return "ok"
