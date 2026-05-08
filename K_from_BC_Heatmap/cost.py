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
from math import gcd

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

def _direction_lattice_steps(dm, dn):
    """Lattice points along the boundary direction (dm, dn) of a torus.

    The line t·(dm, dn) for t ∈ [0, 1] passes through ``g`` lattice
    sites, where ``g = gcd(|dm|, |dn|)`` (with gcd(x, 0) = |x|).  Sites
    are at (k·dm/g, k·dn/g) for k = 0, …, g − 1.

    Returns
    -------
    ks : np.ndarray  — integer indices k = 0, 1, …, g − 1
    ms : np.ndarray  — integer site coordinate m = k·dm/g
    ns : np.ndarray  — integer site coordinate n = k·dn/g
    """
    g = gcd(abs(int(dm)), abs(int(dn)))
    if g == 0:
        return np.array([0]), np.array([0]), np.array([0])
    ks = np.arange(g, dtype=int)
    ms = ks * (int(dm) // g)
    ns = ks * (int(dn) // g)
    return ks, ms, ns


def _wrap_to_test_torus(m, n, Lx, Ly, Tx, Ty):
    """Reduce (m, n) into a canonical torus key by undoing torus shifts.

    The two-point file stores (m, n) with the convention used by the
    simulator.  For lookup we try the original key, then the eight
    neighbouring torus copies, returning the first one that exists in
    ``data``.  Centralised here so callers don't reinvent it.
    """
    # Generated by callers; this helper is intentionally tiny.
    raise NotImplementedError  # placeholder, see _lookup_test_value


def _lookup_test_value(data, m, n, Lx, Ly, Tx, Ty, copies=2):
    """Return data[(m, n)] with torus periodicity, or None if absent."""
    vx_m, vx_n = Lx, Ty
    ux_m, ux_n = Tx, -Ly
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            key = (int(m - a * vx_m - b * ux_m),
                   int(n - a * vx_n - b * ux_n))
            if key in data:
                return data[key]
    return None


def _l2_cost_test_native(ref_data, test_data,
                         test_Lx, test_Ly, test_Tx, test_Ty,
                         ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                         copies=2, power=2):
    """Test-native cost: ref interpolated, test sampled at lattice sites.

    For each of the three boundary directions of the *test* lattice we
    enumerate the integer lattice sites along that period vector and
    look up the test correlator there directly (no interpolation).
    The reference correlator is interpolated at the matching physical
    (x, y) positions of the *reference* boundary direction.

    Per-direction cost:  C_d = mean_k (G_test_k − G_ref_k)^power
    Aggregate cost:      C   = mean_d C_d
    """
    if power not in (2, 4):
        raise ValueError(f"cost_power must be 2 or 4, got {power}")

    iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn", copies)
    iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn_err", copies)

    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    per_dir       = []
    per_dir_sigma = []
    dir_labels    = ["v", "u", "w"]

    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            per_dir.append(0.0)
            per_dir_sigma.append(0.0)
            continue

        # Test values at the integer lattice sites along the test path.
        G_test_list, e_test_list, t_list = [], [], []
        N = len(ks)
        for k, mk, nk in zip(ks, ms, ns):
            entry = _lookup_test_value(test_data, mk, nk,
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is None:
                continue
            G_test_list.append(entry["conn"])
            e_test_list.append(entry["conn_err"])
            t_list.append(k / N)
        if len(G_test_list) < 2:
            per_dir.append(0.0)
            per_dir_sigma.append(0.0)
            continue

        G_test = np.asarray(G_test_list, dtype=float)
        e_test = np.abs(np.asarray(e_test_list, dtype=float))
        t_arr  = np.asarray(t_list, dtype=float)

        # Reference, interpolated at the matching fractional positions
        # along the reference period vector.
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts_ref  = np.column_stack([t_arr * rex, t_arr * rey])
        G_ref = np.asarray(iref(pts_ref), dtype=float)
        e_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))

        mask = (np.isfinite(G_ref) & np.isfinite(G_test)
                & np.isfinite(e_ref) & np.isfinite(e_test))
        if mask.sum() < 2:
            per_dir.append(0.0)
            per_dir_sigma.append(0.0)
            continue

        diff     = G_test[mask] - G_ref[mask]
        var_prop = e_ref[mask]**2 + e_test[mask]**2

        if power == 2:
            C_d = float(np.mean(diff ** 2))
            # σ²(C_d) = (1/N²) Σ_k (2·diff_k)² · var_k
            V_d = float(np.sum((2.0 * diff) ** 2 * var_prop)
                        / (mask.sum() ** 2))
        else:  # power == 4
            C_d = float(np.mean(diff ** 4))
            V_d = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                        / (mask.sum() ** 2))

        per_dir.append(C_d)
        per_dir_sigma.append(math.sqrt(max(V_d, 0.0)))

    # Plain mean over the three directions.
    N_dir = len(per_dir)
    if N_dir == 0:
        cost, sigma_cost = 0.0, 0.0
    else:
        cost       = float(sum(per_dir) / N_dir)
        sigma_cost = math.sqrt(sum(s * s for s in per_dir_sigma)) / N_dir

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(dir_labels, per_dir, per_dir_sigma))
    print(f"    [test-native p={power}] per-dir: {parts}  "
          f"total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


def _l2_cost_huber_log(ref_data, test_data,
                       test_Lx, test_Ly, test_Tx, test_Ty,
                       ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                       copies=2, delta=0.5, eps=1e-12):
    """Huber-on-log-residual cost (test-native sampling).

    Same boundary-direction sampling as ``_l2_cost_test_native`` (ref
    interpolated, test read at lattice sites along each period vector),
    but the per-bin residual is computed in log-space and aggregated with
    the Huber kernel:

        r_k   = log G_test_k − log G_ref_k
        ρ(r)  = ½ r²            if |r| ≤ δ
              = δ(|r| − ½δ)     otherwise
        C_d   = mean_k ρ(r_k)

    Aggregate cost is the plain mean over the three boundary directions.

    The 10-cost stress test (results/_native_triage/stress_costs10_prod)
    showed this kernel had the best mean noise robustness across
    σ ∈ [0, 0.10] of any of 14 candidate cost functions.
    """
    iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn", copies)
    iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn_err", copies)

    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    per_dir, per_dir_sigma = [], []
    dir_labels = ["v", "u", "w"]

    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            per_dir.append(0.0); per_dir_sigma.append(0.0); continue

        G_test_list, e_test_list, t_list = [], [], []
        N = len(ks)
        for k, mk, nk in zip(ks, ms, ns):
            entry = _lookup_test_value(test_data, mk, nk,
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is None:
                continue
            G_test_list.append(entry["conn"])
            e_test_list.append(entry["conn_err"])
            t_list.append(k / N)
        if len(G_test_list) < 2:
            per_dir.append(0.0); per_dir_sigma.append(0.0); continue

        G_test = np.asarray(G_test_list, dtype=float)
        e_test = np.abs(np.asarray(e_test_list, dtype=float))
        t_arr  = np.asarray(t_list, dtype=float)

        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts_ref  = np.column_stack([t_arr * rex, t_arr * rey])
        G_ref = np.asarray(iref(pts_ref), dtype=float)
        e_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))

        # Log-residuals require strictly-positive correlator values.
        mask = (np.isfinite(G_ref) & np.isfinite(G_test)
                & (G_ref > eps) & (G_test > eps)
                & np.isfinite(e_ref) & np.isfinite(e_test))
        if mask.sum() < 2:
            per_dir.append(0.0); per_dir_sigma.append(0.0); continue

        Gt_ = G_test[mask]; Gr_ = G_ref[mask]
        et_ = e_test[mask]; er_ = e_ref[mask]
        r   = np.log(Gt_) - np.log(Gr_)
        a   = np.abs(r)
        # Huber kernel
        rho      = np.where(a <= delta, 0.5 * r * r,
                             delta * (a - 0.5 * delta))
        # gradient ψ(r) for variance propagation: r if |r|≤δ else δ·sign(r)
        psi      = np.where(a <= delta, r, np.sign(r) * delta)
        # var(log G) ≈ (e/G)²  (delta method)
        var_logt = (et_ / Gt_) ** 2
        var_logr = (er_ / Gr_) ** 2

        Nk  = mask.sum()
        C_d = float(np.mean(rho))
        # σ²(C_d) = (1/N²) Σ_k ψ_k² · (var_log_t + var_log_r)
        V_d = float(np.sum(psi ** 2 * (var_logt + var_logr)) / (Nk ** 2))

        per_dir.append(C_d)
        per_dir_sigma.append(math.sqrt(max(V_d, 0.0)))

    N_dir = len(per_dir)
    if N_dir == 0:
        cost, sigma_cost = 0.0, 0.0
    else:
        cost       = float(sum(per_dir) / N_dir)
        sigma_cost = math.sqrt(sum(s * s for s in per_dir_sigma)) / N_dir

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(dir_labels, per_dir, per_dir_sigma))
    print(f"    [huber_log δ={delta:.3f}] per-dir: {parts}  "
          f"total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


def _l2_cost_spectral(ref_data, test_data,
                      k_max=1.6, n_r=6, n_theta=12, k_min=None):
    """All-to-all spectral cost: NDFT both correlators on a polar momentum
    grid, compare modulus-squared residuals.  No interpolation anywhere.

    Per-mode estimate (vectorised NDFT):
        G̃(k) = (1/|Ω|) Σ_{(m,n)∈Ω} G(m,n) exp(-i k·r(m,n))
    with r(m,n) = (m + n/2, n·√3/2) the triangular-lattice physical xy.

    Cost:
        C = mean_α |G̃_test(k_α) − G̃_ref(k_α)|²

    Variance is propagated assuming uncorrelated ``conn_err`` per (m, n)
    entry — same approximation already used by the other cost modes.

    The "per-direction" return values bin the modes into the three angular
    wedges aligned with the boundary directions {v, u, w} for diagnostics.
    """
    if k_min is None:
        k_min = k_max / max(n_r, 1)
    rs  = np.linspace(k_min, k_max, n_r)
    ths = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    R, T = np.meshgrid(rs, ths, indexing="ij")
    kx = (R * np.cos(T)).ravel()
    ky = (R * np.sin(T)).ravel()

    def _ndft(data):
        keys = list(data.keys())
        if not keys:
            return None, 0.0
        m = np.array([k[0] for k in keys], dtype=float)
        n = np.array([k[1] for k in keys], dtype=float)
        G = np.array([data[k]["conn"]     for k in keys], dtype=float)
        e = np.array([data[k]["conn_err"] for k in keys], dtype=float)
        x = m + 0.5 * n
        y = _SQRT3_2 * n
        # phase[α, j] = exp(-i (kx[α]·x[j] + ky[α]·y[j]))
        phase = np.exp(-1j * (np.outer(kx, x) + np.outer(ky, y)))
        Nj = float(len(G))
        Gt  = phase @ G / Nj
        # var per mode α: Σ_j |phase|² · e_j² / Nj² = Σ e_j² / Nj² (|e^iφ|=1)
        var = float(np.sum(e * e)) / (Nj * Nj)
        return Gt, var

    Gt, var_t = _ndft(test_data)
    Gr, var_r = _ndft(ref_data)
    if Gt is None or Gr is None:
        return 0.0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]

    diff = Gt - Gr
    sq   = diff.real ** 2 + diff.imag ** 2
    cost = float(np.mean(sq))
    # Rough scalar: σ²(C) ≈ (4/Nα) Σ_α |Δ|² · (var_t + var_r) / Nα,
    # but with var_t, var_r mode-independent this collapses to:
    sigma_cost = math.sqrt(max(4.0 * cost * (var_t + var_r), 0.0))

    # Per-direction wedge bins (just for printout / panel).
    bp = boundary_paths(1, 1, 0, 0)  # template directions in (m, n)
    # Use canonical triangular-lattice direction angles:
    angles_dir = []
    for dm, dn in [(1, 0), (0, -1), (-1, 1)]:
        x_, y_ = _triangular_xy(dm, dn)
        angles_dir.append(math.atan2(y_, x_))
    mode_ang = np.arctan2(ky, kx)
    per_dir, per_dir_sigma = [], []
    for a0 in angles_dir:
        # closest-wedge mask: angle within ±60° of a0 (mod π for reflection)
        d = np.angle(np.exp(1j * (mode_ang - a0)))
        mask = np.abs(d) < (math.pi / 3.0)
        if mask.any():
            per_dir.append(float(np.mean(sq[mask])))
            per_dir_sigma.append(math.sqrt(
                max(4.0 * per_dir[-1] * (var_t + var_r), 0.0)))
        else:
            per_dir.append(0.0); per_dir_sigma.append(0.0)

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(["v", "u", "w"],
                                         per_dir, per_dir_sigma))
    print(f"    [spectral n_r={n_r} n_θ={n_theta} k∈[{k_min:.2f},{k_max:.2f}]] "
          f"per-wedge: {parts}  total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


def _l2_cost_affine_fit(ref_data, test_data,
                        test_Lx, test_Ly, test_Tx, test_Ty,
                        ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                        copies=2, mode="test_native"):
    """Affine-invariant collapse residual.

    Stack G_ref, G_test on a common set of sites and solve

        (a*, b*) = argmin_{a,b}  || a · G_t + b · 1  -  G_r ||²

    Cost = (1/N) || a* G_t + b* 1 - G_r ||².  Invariant under any
    overall rescaling or constant shift of the test correlator —
    projects out β_c-jitter normalisation drift and residual ⟨m⟩²
    offset.  Validated in cost_zoo_test.py: ~1500× MC-noise
    suppression, |Z|≈12-160 vs native 0.04-30 across 6 mismatch tests.

    mode :
      "all_to_all"  — use the intersection of (m,n) keys present in
                      both dicts (works for matching geometries).
      "test_native" — only sample G_t at test-lattice sites along the
                      three test boundary directions; lookup G_r via
                      the same tiled interpolant used by other modes,
                      so this works cross-geometry.

    Per-direction values returned for diagnostics; in 'all_to_all'
    mode the directions are angular wedges of the (m,n) cloud.
    """
    if mode == "all_to_all":
        common = sorted(set(ref_data.keys()) & set(test_data.keys()))
        if len(common) < 4:
            return 0.0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        Gr = np.array([ref_data[k]["conn"]      for k in common], float)
        Gt = np.array([test_data[k]["conn"]     for k in common], float)
        er = np.array([ref_data[k]["conn_err"]  for k in common], float)
        et = np.array([test_data[k]["conn_err"] for k in common], float)
        # Per-direction wedges (angles of (m, n) → physical xy).
        m = np.array([k[0] for k in common], float)
        n = np.array([k[1] for k in common], float)
        ang = np.arctan2(_SQRT3_2 * n, m + 0.5 * n)
    elif mode == "test_native":
        iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                                "conn", copies)
        iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                                "conn_err", copies)
        ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
        test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)
        Gr_l, Gt_l, er_l, et_l, ang_l = [], [], [], [], []
        for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
            ks, ms, ns = _direction_lattice_steps(tdm, tdn)
            if len(ks) < 2:
                continue
            N = len(ks)
            for k, mk, nk in zip(ks, ms, ns):
                entry = _lookup_test_value(test_data, mk, nk,
                                           test_Lx, test_Ly, test_Tx, test_Ty,
                                           copies=copies)
                if entry is None:
                    continue
                t = k / N
                rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
                Gr_v = float(iref(np.array([[t * rex, t * rey]]))[0])
                er_v = float(np.abs(iref_err(np.array([[t * rex, t * rey]]))[0]))
                if not (np.isfinite(Gr_v) and np.isfinite(er_v)):
                    continue
                Gt_l.append(entry["conn"]); et_l.append(entry["conn_err"])
                Gr_l.append(Gr_v); er_l.append(er_v)
                ang_l.append(math.atan2(rey, rex))
        if len(Gr_l) < 4:
            return 0.0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        Gr = np.array(Gr_l, float); Gt = np.array(Gt_l, float)
        er = np.array(er_l, float); et = np.array(et_l, float)
        ang = np.array(ang_l, float)
    else:
        raise ValueError(f"unknown affine_fit mode={mode!r}")

    # Closed-form 2-parameter LSQ.
    A = np.column_stack([Gt, np.ones_like(Gt)])
    sol, *_ = np.linalg.lstsq(A, Gr, rcond=None)
    a_star, b_star = float(sol[0]), float(sol[1])
    res = Gr - (a_star * Gt + b_star)
    cost = float(np.mean(res ** 2))

    # Variance: σ²_residual_k = a²·σ²_test_k + σ²_ref_k (b is a global
    # constant; the LSQ fit reduces a fluctuation in 'a' that is small
    # for large N, so we treat (a*, b*) as fixed when propagating.)
    var_k = (a_star ** 2) * (et ** 2) + (er ** 2) + 1e-30
    Nk = len(res)
    var_cost = float(np.sum((2.0 * res) ** 2 * var_k) / (Nk ** 2))
    sigma_cost = math.sqrt(max(var_cost, 0.0))

    # Per-direction wedges (angular bins of |a*Gt+b - Gr|² aligned with
    # the three triangular-lattice canonical directions).
    angles_dir = []
    for dm, dn in [(1, 0), (0, -1), (-1, 1)]:
        x_, y_ = _triangular_xy(dm, dn)
        angles_dir.append(math.atan2(y_, x_))
    per_dir, per_dir_sigma = [], []
    sq = res ** 2
    for a0 in angles_dir:
        d = np.angle(np.exp(1j * (ang - a0)))
        mask = np.abs(d) < (math.pi / 3.0)
        if mask.any():
            per_dir.append(float(np.mean(sq[mask])))
            per_dir_sigma.append(math.sqrt(
                float(np.sum((2.0 * res[mask]) ** 2 * var_k[mask])
                      / (mask.sum() ** 2))))
        else:
            per_dir.append(0.0); per_dir_sigma.append(0.0)

    parts = "  ".join(f"{l}={c:.4e}±{s:.4e}"
                      for l, c, s in zip(["v", "u", "w"],
                                         per_dir, per_dir_sigma))
    print(f"    [affine_fit mode={mode} a={a_star:+.4f} b={b_star:+.2e} "
          f"N={Nk}] per-wedge: {parts}  total={cost:.4e}±{sigma_cost:.4e}")

    return cost, sigma_cost, per_dir, per_dir_sigma


def l2_cost(ref_data, test_data,
            test_Lx, test_Ly, test_Tx, test_Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            n_samples=400, copies=2,
            cost_mode="l4mean_both_interp",
            cost_power=2,
            spectral_kwargs=None,
            affine_kwargs=None,
            huber_kwargs=None):
    """
    Cost between reference and test correlators along the three boundary
    directions.

    Parameters
    ----------
    cost_mode : str
        - ``"l4mean_both_interp"`` (default, legacy):
              Both reference and test correlators are interpolated and
              sampled on the same uniform grid of ``n_samples`` points
              per direction; per-direction L² values are aggregated with
              an L⁴ power-mean.
        - ``"test_native"``:
              Only the reference correlator is interpolated.  The test
              correlator is read directly at the test lattice sites that
              lie along each boundary direction (no interpolation).  The
              cost is a plain mean over directions of the mean of
              residuals raised to ``cost_power`` (2 → squared, 4 → quart).
    cost_power : int
        Exponent applied to per-site residuals in ``"test_native"`` mode.
        Must be 2 or 4.  Ignored in ``"l4mean_both_interp"`` mode.

    Returns
    -------
    cost           : float
    sigma_cost     : float
    per_dir        : list[float]
    per_dir_sigma  : list[float]
    """
    if cost_mode == "test_native":
        return _l2_cost_test_native(
            ref_data, test_data,
            test_Lx, test_Ly, test_Tx, test_Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            copies=copies, power=int(cost_power),
        )
    if cost_mode == "spectral":
        kw = dict(spectral_kwargs or {})
        return _l2_cost_spectral(ref_data, test_data, **kw)
    if cost_mode == "affine_fit":
        kw = dict(affine_kwargs or {})
        return _l2_cost_affine_fit(
            ref_data, test_data,
            test_Lx, test_Ly, test_Tx, test_Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            copies=copies, **kw)
    if cost_mode == "huber_log":
        kw = dict(huber_kwargs or {})
        return _l2_cost_huber_log(
            ref_data, test_data,
            test_Lx, test_Ly, test_Tx, test_Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            copies=copies, **kw)
    if cost_mode != "l4mean_both_interp":
        raise ValueError(f"unknown cost_mode={cost_mode!r}")

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
