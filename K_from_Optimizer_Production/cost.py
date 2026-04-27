#!/usr/bin/env python3
"""
cost.py — Anisotropy-penalising p-mean cost (p=4) over the three torus
boundary directions.

For each of the three fundamental torus boundary directions d ∈ {v, u, w},
the correlator is sampled along the corresponding lattice period vector,
parameterised by the arc-length fraction t ∈ [0, 1] *of that period*.
The per-direction L² mismatch is

    C_d = ∫₀¹ [G_test^d(t) - G_ref^d(t)]² dt

and these are aggregated with a p-power mean (p = 4):

    C(r1,r2) = ( (1/3) Σ_d C_d^p )^(1/p)

This punishes anisotropy: a single direction with large C_d dominates
the total, so the optimizer cannot trade a small gain in one direction
for a large drift in another.  At perfect isotropy (C_v = C_u = C_w)
the p-mean equals the common per-direction value.

Curve collapse across differing ref / test geometries
-----------------------------------------------------
The reference and test lattices need NOT share (Lx, Ly, Tx, Ty).  Each
direction's path endpoint is its *own* lattice period vector, so e.g.

    ref  = (13, 16,  3, -3)  →  periods  (13,-3), (3,-16), (-16,19)
    test = (15, 15,  0,  0)  →  periods  (15, 0), (0,-15), (-15,15)

are each parameterised over the unit interval and compared point-by-point
in t.  This is the fractional-arclength curve-collapse hypothesis: at the
critical point, the boundary-directed two-point function along cycle i of
torus A should collapse onto the same function of t as cycle i of torus B.

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

    # Torus periodicity: t=0 and t=1 are the same lattice site, so we sample
    # only t ∈ [0, 1) and use the periodic trapezoid rule (= uniform sum /
    # n_samples).  Including t=1 pushes LinearNDInterpolator onto the
    # convex-hull boundary of the tiled point cloud, which gives a spurious
    # spike in the residuals at t ≈ 1.
    t = np.linspace(0.0, 1.0, n_samples, endpoint=False)

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

        # Periodic trapezoid rule over t ∈ [0, 1) with n_samples equal
        # intervals of width 1/n_samples — which reduces to a mean.
        # Both the cost integral and its variance-propagation integral
        # use the same ∫…dt convention as before, so SNR thresholds
        # retain their previous calibration.
        C_d = float(np.mean(diff ** 2))
        V_d = float(np.mean(4.0 * diff ** 2 * var_prop))

        per_dir.append(C_d)
        per_dir_sigma.append(math.sqrt(max(V_d, 0.0)))

    # Aggregate across directions with a p-power mean (p=4) instead of a
    # plain sum.  This penalises anisotropy: a single direction with a
    # large C_d dominates the total, so the optimizer cannot trade a
    # small gain in one direction for a large drift in another.  At
    # perfect isotropy (C_v = C_u = C_w) the p-mean equals each C_d.
    #
    #     cost = ( (1/N) Σ_d C_d^p )^(1/p),  p = 4, N = 3
    #
    # Error propagation:
    #     ∂cost/∂C_d = (1/N) C_d^(p-1) · cost^(1-p)
    #     σ_cost²    = Σ_d (∂cost/∂C_d)² σ_{C_d}²
    P     = 4
    N_dir = len(per_dir)
    if N_dir == 0 or all(c == 0.0 for c in per_dir):
        cost       = 0.0
        sigma_cost = 0.0
    else:
        mean_p = sum(c ** P for c in per_dir) / N_dir
        cost   = mean_p ** (1.0 / P)
        if cost > 0.0:
            grads = [(c ** (P - 1)) / (N_dir * cost ** (P - 1))
                     for c in per_dir]
            sigma_cost = math.sqrt(sum((g * s) ** 2
                                       for g, s in zip(grads, per_dir_sigma)))
        else:
            sigma_cost = 0.0

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(dir_labels, per_dir, per_dir_sigma))
    print(f"    [L⁴-mean] per-dir: {parts}  "
          f"total={cost:.4e}±{sigma_cost:.4e}")

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
