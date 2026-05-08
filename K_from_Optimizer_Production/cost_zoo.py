#!/usr/bin/env python3
"""
cost_zoo.py — alternative cost functions targeting the failure modes
diagnosed in the realmc_456_from11 run.

DIAGNOSIS RECAP
---------------
``cost._l2_cost_test_native`` reparameterizes both ref and test paths to
t∈[0, 1] regardless of physical length, then takes the plain mean of
three directions.  This admits a misleading minimum near (r1, r2)=(1, 1)
where the test correlator's *shape* (in normalized t) accidentally
matches the ref's, even though the physical lattice spacing of the two
periods differs by ~16% and the (k1, k2) coupling is nowhere near truth.
The "w" diagonal direction also dominates the cost by 1-2 decades, so
the plain mean reduces to "match the longest diagonal."

EACH ZOO COST ADDRESSES A SPECIFIC LEAK
---------------------------------------
1. ``_zoo_relative_native``     — divide each per-dir cost by mean(G_ref²) on
                                  the same path; dimensionless, kills the
                                  magnitude-driven dominance of one direction.

2. ``_zoo_log_native``          — compare log|G_test| vs log|G_ref|; this is
                                  effective-mass language, scale-invariant
                                  under any overall (β-driven) rescaling.

3. ``_zoo_arclength_native``    — sample both ref and test correlators at
                                  the SAME physical lattice distances
                                  s_k = k * Δs (k = 1, …, K) along each
                                  direction, with K bounded by the SHORTER
                                  physical period.  This is the direct fix
                                  for the t = k/N reparameterization bug.

4. ``_zoo_effective_mass``      — for each direction compute
                                  m_eff_k = -log(G_{k+1}/G_k) and compare
                                  m_eff curves between ref (interpolated)
                                  and test.  Most physically meaningful;
                                  insensitive to overall amplitude.

5. ``_zoo_pmean_native_p4``     — p=4 anisotropy-penalising mean of the
                                  three test_native per-dir costs (existing
                                  l4mean_both_interp uses p=4 too but with
                                  different sampling); penalises any single
                                  direction running away.

All five share the signature of ``cost._l2_cost_test_native``:
    (ref_data, test_data,
     test_Lx, test_Ly, test_Tx, test_Ty,
     ref_Lx, ref_Ly, ref_Tx, ref_Ty,
     copies=2, power=2)  →  (cost, sigma, per_dir, per_dir_sigma)

USAGE
-----
Monkey-patch into a CMA-ES run:

    import cost as cost_mod, cost_zoo
    cost_mod._l2_cost_test_native = cost_zoo._zoo_arclength_native

then create an Evaluator with cost_mode="test_native".  See
``cost_zoo_driver.py`` for an end-to-end runner.
"""
from __future__ import annotations

import math
from math import gcd

import numpy as np

# Reuse helpers from the production cost module.
from cost import (
    _SQRT3_2,
    _tile_interp,
    _lookup_test_value,
    _direction_lattice_steps,
    boundary_paths,
)

DIR_LABELS = ["v", "u", "w"]


# ---------------------------------------------------------------------------
# Shared sampler: returns (G_test, G_ref, e_test, e_ref, t_arr) per direction
# using the existing test_native logic (k/N reparameterization).
# ---------------------------------------------------------------------------
def _sample_per_direction_native(ref_data, test_data,
                                 test_Lx, test_Ly, test_Tx, test_Ty,
                                 ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                                 copies=2):
    iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn", copies)
    iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn_err", copies)
    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    out = []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            out.append(None); continue

        Gt_l, et_l, t_l = [], [], []
        N = len(ks)
        for k, mk, nk in zip(ks, ms, ns):
            entry = _lookup_test_value(test_data, mk, nk,
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is None:
                continue
            Gt_l.append(entry["conn"])
            et_l.append(entry["conn_err"])
            t_l.append(k / N)
        if len(Gt_l) < 2:
            out.append(None); continue

        Gt = np.asarray(Gt_l, float)
        et = np.abs(np.asarray(et_l, float))
        t  = np.asarray(t_l, float)

        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts_ref = np.column_stack([t * rex, t * rey])
        Gr = np.asarray(iref(pts_ref), float)
        er = np.abs(np.asarray(iref_err(pts_ref), float))

        mask = (np.isfinite(Gt) & np.isfinite(Gr)
                & np.isfinite(et) & np.isfinite(er))
        if mask.sum() < 2:
            out.append(None); continue
        out.append((Gt[mask], Gr[mask], et[mask], er[mask], t[mask]))
    return out


def _aggregate(per_dir, per_sig, label):
    """Plain mean over non-empty directions; identical to T3-patched native."""
    valid = [(c, s) for c, s in zip(per_dir, per_sig)
             if not (c == 0.0 and s == 0.0)]
    if not valid:
        return float("nan"), float("nan")
    N = len(valid)
    cost = sum(c for c, _ in valid) / N
    sig  = math.sqrt(sum(s * s for _, s in valid)) / N
    parts = "  ".join(f"{l}={c:.3e}±{s:.3e}"
                      for l, c, s in zip(DIR_LABELS, per_dir, per_sig))
    print(f"    [zoo {label}] {parts}  total={cost:.3e}±{sig:.3e}")
    return cost, sig


# ===========================================================================
# 1. RELATIVE NATIVE — per-direction cost normalized by mean(G_ref²)
# ===========================================================================
def _zoo_relative_native(ref_data, test_data,
                         test_Lx, test_Ly, test_Tx, test_Ty,
                         ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                         copies=2, power=2):
    samples = _sample_per_direction_native(
        ref_data, test_data,
        test_Lx, test_Ly, test_Tx, test_Ty,
        ref_Lx, ref_Ly, ref_Tx, ref_Ty, copies=copies)
    per_dir, per_sig = [], []
    for s in samples:
        if s is None:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er, _ = s
        diff = Gt - Gr
        var_prop = er ** 2 + et ** 2
        denom = float(np.mean(Gr ** 2)) + 1e-30
        if power == 2:
            C = float(np.mean(diff ** 2)) / denom
            V = float(np.sum((2.0 * diff) ** 2 * var_prop)
                      / (len(diff) ** 2)) / (denom ** 2)
        else:
            C = float(np.mean(diff ** 4)) / denom
            V = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                      / (len(diff) ** 2)) / (denom ** 2)
        per_dir.append(C); per_sig.append(math.sqrt(max(V, 0.0)))
    cost, sig = _aggregate(per_dir, per_sig, f"relative p={power}")
    return cost, sig, per_dir, per_sig


# ===========================================================================
# 2. LOG NATIVE — compare log|G| (effective-mass language)
# ===========================================================================
def _zoo_log_native(ref_data, test_data,
                    test_Lx, test_Ly, test_Tx, test_Ty,
                    ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                    copies=2, power=2):
    EPS = 1e-12
    samples = _sample_per_direction_native(
        ref_data, test_data,
        test_Lx, test_Ly, test_Tx, test_Ty,
        ref_Lx, ref_Ly, ref_Tx, ref_Ty, copies=copies)
    per_dir, per_sig = [], []
    for s in samples:
        if s is None:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er, _ = s
        m = (Gt > EPS) & (Gr > EPS)
        if m.sum() < 2:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er = Gt[m], Gr[m], et[m], er[m]
        lt, lr = np.log(Gt), np.log(Gr)
        diff = lt - lr
        # σ(log G) ≈ σ(G)/G
        var_prop = (et / Gt) ** 2 + (er / Gr) ** 2
        if power == 2:
            C = float(np.mean(diff ** 2))
            V = float(np.sum((2.0 * diff) ** 2 * var_prop)
                      / (len(diff) ** 2))
        else:
            C = float(np.mean(diff ** 4))
            V = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                      / (len(diff) ** 2))
        per_dir.append(C); per_sig.append(math.sqrt(max(V, 0.0)))
    cost, sig = _aggregate(per_dir, per_sig, f"log p={power}")
    return cost, sig, per_dir, per_sig


# ===========================================================================
# 3. ARCLENGTH NATIVE — sample at SAME physical distances on both lattices
# ===========================================================================
def _zoo_arclength_native(ref_data, test_data,
                          test_Lx, test_Ly, test_Tx, test_Ty,
                          ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                          copies=2, power=2):
    """Both correlators sampled at physical distances s = 1, 2, …, K
    along each direction, where K = floor(min(L_test_d, L_ref_d)).
    Test sampled by torus lookup at the nearest lattice point along the
    direction, ref via tiled interpolation at the matching (x, y).
    """
    iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn", copies)
    iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn_err", copies)
    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    per_dir, per_sig = [], []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        # Physical lengths along this direction, in lattice units.
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        tex, tey = tdm + 0.5 * tdn, _SQRT3_2 * tdn
        Lref = math.hypot(rex, rey)
        Ltest = math.hypot(tex, tey)
        Lmin = min(Lref, Ltest)
        if Lmin < 2.0:
            per_dir.append(0.0); per_sig.append(0.0); continue

        K = int(math.floor(Lmin))      # samples at s = 1, 2, …, K-1
        if K < 3:
            per_dir.append(0.0); per_sig.append(0.0); continue
        ss = np.arange(1, K, dtype=float)        # physical distances

        # Reference at (s/Lref) * (rex, rey)
        u_ref = ss / Lref
        Gr = np.asarray(iref(np.column_stack([u_ref * rex, u_ref * rey])),
                        float)
        er = np.abs(np.asarray(iref_err(np.column_stack(
            [u_ref * rex, u_ref * rey])), float))

        # Test: nearest lattice site along the test direction at distance s
        # In test_native gcd-stride mode the available steps are at
        # s_k = k * (Ltest / N_test_steps), k = 0..N_test_steps-1.
        ks_t, ms_t, ns_t = _direction_lattice_steps(tdm, tdn)
        if len(ks_t) < 2:
            per_dir.append(0.0); per_sig.append(0.0); continue
        N_t = len(ks_t)
        step_phys = Ltest / N_t  # spacing of available test sites in physical units
        Gt_l, et_l = [], []
        for s in ss:
            k_idx = int(round(s / step_phys))
            if k_idx <= 0 or k_idx >= N_t:
                Gt_l.append(np.nan); et_l.append(np.nan); continue
            mk, nk = ms_t[k_idx], ns_t[k_idx]
            entry = _lookup_test_value(test_data, mk, nk,
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is None:
                Gt_l.append(np.nan); et_l.append(np.nan); continue
            Gt_l.append(entry["conn"]); et_l.append(entry["conn_err"])
        Gt = np.asarray(Gt_l, float); et = np.abs(np.asarray(et_l, float))

        mask = (np.isfinite(Gt) & np.isfinite(Gr)
                & np.isfinite(et) & np.isfinite(er))
        if mask.sum() < 2:
            per_dir.append(0.0); per_sig.append(0.0); continue
        diff = Gt[mask] - Gr[mask]
        var_prop = er[mask] ** 2 + et[mask] ** 2
        if power == 2:
            C = float(np.mean(diff ** 2))
            V = float(np.sum((2.0 * diff) ** 2 * var_prop)
                      / (mask.sum() ** 2))
        else:
            C = float(np.mean(diff ** 4))
            V = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                      / (mask.sum() ** 2))
        per_dir.append(C); per_sig.append(math.sqrt(max(V, 0.0)))
    cost, sig = _aggregate(per_dir, per_sig, f"arclength p={power}")
    return cost, sig, per_dir, per_sig


# ===========================================================================
# 4. EFFECTIVE-MASS — diff of m_eff(s) = -log(G_{k+1}/G_k)
# ===========================================================================
def _zoo_effective_mass(ref_data, test_data,
                        test_Lx, test_Ly, test_Tx, test_Ty,
                        ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                        copies=2, power=2):
    EPS = 1e-12
    samples = _sample_per_direction_native(
        ref_data, test_data,
        test_Lx, test_Ly, test_Tx, test_Ty,
        ref_Lx, ref_Ly, ref_Tx, ref_Ty, copies=copies)
    per_dir, per_sig = [], []
    for s in samples:
        if s is None:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er, _ = s
        m = (Gt > EPS) & (Gr > EPS)
        if m.sum() < 3:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er = Gt[m], Gr[m], et[m], er[m]
        # m_eff at midpoint k+1/2: -log(G_{k+1}/G_k)
        meff_t = -np.log(Gt[1:] / Gt[:-1])
        meff_r = -np.log(Gr[1:] / Gr[:-1])
        diff = meff_t - meff_r
        # σ²(meff) = (σ_{k+1}/G_{k+1})² + (σ_k/G_k)²
        var_t = (et[1:] / Gt[1:]) ** 2 + (et[:-1] / Gt[:-1]) ** 2
        var_r = (er[1:] / Gr[1:]) ** 2 + (er[:-1] / Gr[:-1]) ** 2
        var_prop = var_t + var_r
        if power == 2:
            C = float(np.mean(diff ** 2))
            V = float(np.sum((2.0 * diff) ** 2 * var_prop)
                      / (len(diff) ** 2))
        else:
            C = float(np.mean(diff ** 4))
            V = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                      / (len(diff) ** 2))
        per_dir.append(C); per_sig.append(math.sqrt(max(V, 0.0)))
    cost, sig = _aggregate(per_dir, per_sig, f"effmass p={power}")
    return cost, sig, per_dir, per_sig


# ===========================================================================
# 5. P4-MEAN NATIVE — anisotropy-penalising mean of test_native per-dir
# ===========================================================================
def _zoo_pmean_native_p4(ref_data, test_data,
                         test_Lx, test_Ly, test_Tx, test_Ty,
                         ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                         copies=2, power=2, p_aggregate=4):
    samples = _sample_per_direction_native(
        ref_data, test_data,
        test_Lx, test_Ly, test_Tx, test_Ty,
        ref_Lx, ref_Ly, ref_Tx, ref_Ty, copies=copies)
    per_dir, per_sig = [], []
    for s in samples:
        if s is None:
            per_dir.append(0.0); per_sig.append(0.0); continue
        Gt, Gr, et, er, _ = s
        diff = Gt - Gr
        var_prop = er ** 2 + et ** 2
        if power == 2:
            C = float(np.mean(diff ** 2))
            V = float(np.sum((2.0 * diff) ** 2 * var_prop)
                      / (len(diff) ** 2))
        else:
            C = float(np.mean(diff ** 4))
            V = float(np.sum((4.0 * diff ** 3) ** 2 * var_prop)
                      / (len(diff) ** 2))
        per_dir.append(C); per_sig.append(math.sqrt(max(V, 0.0)))
    valid = [(c, s) for c, s in zip(per_dir, per_sig)
             if not (c == 0.0 and s == 0.0)]
    if not valid:
        return float("nan"), float("nan"), per_dir, per_sig
    p = int(p_aggregate)
    N = len(valid)
    inner = sum(c ** p for c, _ in valid) / N
    cost = inner ** (1.0 / p)
    # σ_cost via chain rule on the p-mean
    # ∂cost/∂C_d = (1/N) C_d^{p-1} * inner^{(1/p) - 1}
    sig2 = 0.0
    for (c, s) in valid:
        if c <= 0:
            continue
        deriv = (c ** (p - 1)) * (inner ** (1.0 / p - 1.0)) / N
        sig2 += (deriv * s) ** 2
    sig = math.sqrt(sig2)
    parts = "  ".join(f"{l}={c:.3e}±{s:.3e}"
                      for l, c, s in zip(DIR_LABELS, per_dir, per_sig))
    print(f"    [zoo pmean p_agg={p}] {parts}  total={cost:.3e}±{sig:.3e}")
    return cost, sig, per_dir, per_sig


# ===========================================================================
# Registry
# ===========================================================================
ZOO = {
    "relative":   _zoo_relative_native,
    "log":        _zoo_log_native,
    "arclength":  _zoo_arclength_native,
    "effmass":    _zoo_effective_mass,
    "pmean4":     _zoo_pmean_native_p4,
}


def install(name: str):
    """Monkey-patch cost._l2_cost_test_native with a zoo variant."""
    if name not in ZOO:
        raise KeyError(f"unknown zoo cost {name!r}; choices: {sorted(ZOO)}")
    import cost as cost_mod
    cost_mod._l2_cost_test_native = ZOO[name]
    print(f"  cost_zoo: installed '{name}' as cost._l2_cost_test_native")
