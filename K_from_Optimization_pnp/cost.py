#!/usr/bin/env python3
"""
cost.py — Plain L² cost function for the curve-collapse matching problem.

For each boundary direction d ∈ {v, u, w}, the correlator is sampled along a
straight path parameterised by t ∈ [0, 1] (arc-length fraction).  The cost is
the integrated squared difference summed over all three directions:

    C_d   = ∫₀¹ [G_test^d(t) - G_ref^d(t)]² dt
    C(r1,r2) = Σ_d C_d

The integral is approximated by the trapezoid rule on n_samples uniformly
spaced points.

Uncertainty in C
----------------
By linear error propagation:

    σ_C² ≈ Σ_d  ∫₀¹ 4 (G_test - G_ref)² (σ_test² + σ_ref²) dt

so that SNR = C / σ_C tells the optimizer when its current best is
statistically distinguishable from the reference.

Public API
----------
  boundary_paths(Lx, Ly, Tx, Ty) -> list of three (dm, dn) direction vectors
  l2_cost(ref_data, test_data, ...) -> (cost, sigma_cost, per_dir, per_dir_sigma)
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
    """Three fundamental boundary direction vectors of the torus."""
    return [
        (Lx,        Ty),
        (Tx,        -Ly),
        (-Lx - Tx,  Ly - Ty),
    ]


_SQRT3_2 = math.sqrt(3.0) / 2.0


def _triangular_xy(m, n):
    return m + 0.5 * n, _SQRT3_2 * n


def _tile_interp(data, Lx, Ly, Tx, Ty, field, copies=2):
    """LinearNDInterpolator over the data tiled across torus copies."""
    keys = list(data.keys())
    m_arr = np.array([k[0] for k in keys], dtype=float)
    n_arr = np.array([k[1] for k in keys], dtype=float)
    c_arr = np.array([data[k][field] for k in keys], dtype=float)

    x0, y0 = _triangular_xy(m_arr, n_arr)
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
            n_samples=400, copies=2):
    """
    Plain L² cost between reference and test correlators along the three
    boundary directions.

    Returns
    -------
    cost           : float — Σ_d ∫(G_test - G_ref)² dt
    sigma_cost     : float — propagated 1-σ uncertainty in cost
    per_dir        : list[float] — [C_v, C_u, C_w]
    per_dir_sigma  : list[float] — per-direction σ
    """
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
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        tex, tey = tdm + 0.5 * tdn, _SQRT3_2 * tdn
        pts_ref  = np.column_stack([t * rex, t * rey])
        pts_test = np.column_stack([t * tex, t * tey])

        G_ref  = np.asarray(iref(pts_ref),       dtype=float)
        G_test = np.asarray(itest(pts_test),     dtype=float)
        e_ref  = np.abs(np.asarray(iref_err(pts_ref),   dtype=float))
        e_test = np.abs(np.asarray(itest_err(pts_test), dtype=float))

        mask = (np.isfinite(G_ref) & np.isfinite(G_test)
                & np.isfinite(e_ref) & np.isfinite(e_test))

        if mask.sum() < 4:
            per_dir.append(0.0)
            per_dir_sigma.append(0.0)
            continue

        diff     = G_test[mask] - G_ref[mask]
        var_prop = e_ref[mask]**2 + e_test[mask]**2
        t_m      = t[mask]
        dt       = t_m[1:] - t_m[:-1]

        integrand_c = diff ** 2
        C_d = float(np.sum(0.5 * (integrand_c[:-1] + integrand_c[1:]) * dt))

        integrand_v = 4.0 * diff**2 * var_prop
        V_d = float(np.sum(0.5 * (integrand_v[:-1] + integrand_v[1:]) * dt))

        per_dir.append(C_d)
        per_dir_sigma.append(math.sqrt(max(V_d, 0.0)))

    cost       = sum(per_dir)
    sigma_cost = math.sqrt(sum(s**2 for s in per_dir_sigma))

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(dir_labels, per_dir, per_dir_sigma))
    print(f"    [L²] per-dir: {parts}  total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


# ---------------------------------------------------------------------------
# SNR helpers
# ---------------------------------------------------------------------------

def snr(cost, sigma_cost):
    if sigma_cost <= 0.0:
        return float("inf")
    return cost / sigma_cost


def snr_status(cost, sigma_cost):
    s = snr(cost, sigma_cost)
    if s < 1.0:
        return "need_more_stats"
    if s < 3.0:
        return "marginal"
    return "ok"
