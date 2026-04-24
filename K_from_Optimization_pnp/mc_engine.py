#!/usr/bin/env python3
"""
mc_engine.py — Monte Carlo engine for the standalone critical-surface pipeline.

Provides:
  - run_simulator()           : invoke the C++ Ising simulator
  - find_beta_c()             : 3-pass Gram-Charlier susceptibility-peak finder
  - load_all_to_all()         : load two_point_all_to_all.dat
  - boundary_slices_normed()  : Z² cost (variance-normalised boundary slices)
  - boundary_slices_quartic() : Z⁴ cost (quartic variant)
  - extract_boundary_slices() : per-direction slices for plotting

Adapted from K_from_TwoPoint_standalone/optimise_couplings.py,
stripped to the essentials needed by the pipeline.
"""

import math
import os
import subprocess
import sys
import time
import warnings

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import curve_fit, minimize_scalar


# ===================================================================
# Parsing helpers
# ===================================================================

def parse_stdout_value(stdout: str, key: str) -> float:
    for line in stdout.splitlines():
        if line.startswith(key):
            return float(line.split()[1])
    raise ValueError(f"Key '{key}' not found in simulator output")


def parse_stdout_value_with_err(stdout: str, key: str):
    for line in stdout.splitlines():
        if line.startswith(key):
            parts = line.split()
            return float(parts[1]), float(parts[2])
    raise ValueError(f"Key '{key}' not found in simulator output")


def load_all_to_all(path: str):
    """Load two_point_all_to_all.dat → dict keyed by (m, n)."""
    data = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            d, m, n = int(parts[0]), int(parts[1]), int(parts[2])
            data[(m, n)] = {
                "d": d,
                "corr": float(parts[3]), "err": float(parts[4]),
                "conn": float(parts[5]), "conn_err": float(parts[6]),
            }
    return data


# ===================================================================
# Simulator runner
# ===================================================================

def run_simulator(exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta,
                  n_traj=20000, n_skip=20, n_therm=2000,
                  seed=0, data_dir="/tmp/k_from_cs",
                  wall_time=0.0):
    """Run the C++ simulator and return (stdout_text, subdir_path)."""
    if seed == 0:
        seed = int(time.time() * 1000) & 0xFFFFFFFF
    os.makedirs(data_dir, exist_ok=True)
    cmd = [
        exe,
        "--L_x", str(Lx), "--L_y", str(Ly),
        "--T_x", str(Tx), "--T_y", str(Ty),
        "--k1", f"{k1:.12f}", "--k2", f"{k2:.12f}", "--k3", f"{k3:.12f}",
        "--beta", f"{beta:.12f}",
        "--n_traj", str(n_traj), "--n_skip", str(n_skip),
        "--n_therm", str(n_therm),
        "--seed", str(seed),
        "--data_dir", data_dir,
    ]
    if wall_time > 0:
        cmd += ["--wall_time", f"{wall_time:.1f}"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except FileNotFoundError:
        print(f"\nERROR: simulator not found: {exe}", file=sys.stderr)
        raise
    if result.returncode != 0:
        print(f"ERROR: simulator failed:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError("simulator failed")

    subdir = None
    for line in result.stdout.splitlines():
        if line.startswith("opening file:"):
            full_path = line.split(":", 1)[1].strip()
            subdir = os.path.dirname(full_path)
            break
    return result.stdout, subdir


# ===================================================================
# Beta_c finder (3-pass Gram-Charlier)
# ===================================================================

def _gram_charlier(b, mu, sigma, skew, kurt, A):
    z = (b - mu) / sigma
    H3 = z**3 - 3.0 * z
    H4 = z**4 - 6.0 * z**2 + 3.0
    return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)


def _gc_fit(betas_arr, chis_arr, beta_init):
    """Fit Gram-Charlier to data; return (popt_tuple, mode) or (None, beta_init)."""
    n = len(betas_arr)
    b_range = max(betas_arr) - min(betas_arr)
    if b_range < 1e-12:
        return None, beta_init
    use_gc = (n >= 6)
    try:
        if use_gc:
            p0 = [beta_init, b_range / 2, 0.0, 0.0, float(max(chis_arr))]
            bounds = ([min(betas_arr), 1e-8, -5.0, -3.0, 0.0],
                      [max(betas_arr), b_range * 2, 5.0, 10.0, np.inf])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(_gram_charlier, np.array(betas_arr),
                                    np.array(chis_arr),
                                    p0=p0, bounds=bounds, maxfev=8000)
        else:
            def _gauss3(b, mu, sig, A):
                return A * np.exp(-0.5 * ((b - mu) / sig) ** 2)
            p0 = [beta_init, b_range / 2, float(max(chis_arr))]
            bounds = ([min(betas_arr), 1e-8, 0.0],
                      [max(betas_arr), b_range * 2, np.inf])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pg, _ = curve_fit(_gauss3, np.array(betas_arr),
                                  np.array(chis_arr),
                                  p0=p0, bounds=bounds, maxfev=8000)
            popt = (pg[0], pg[1], 0.0, 0.0, pg[2])
        mu_fit, sigma_fit = popt[0], popt[1]
        mode_lo = max(min(betas_arr), mu_fit - 3.0 * sigma_fit)
        mode_hi = min(max(betas_arr), mu_fit + 3.0 * sigma_fit)
        res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                              bounds=(mode_lo, mode_hi), method='bounded')
        mode = float(np.clip(res.x, min(betas_arr), max(betas_arr)))
        return tuple(popt), mode
    except Exception:
        return None, beta_init


def _gc_jackknife_sigma(scan_betas, scan_chis, beta_init):
    """Leave-one-out jackknife estimate of σ(β_c) over the scan points.

    For each i, refit the GC curve on the remaining N-1 points and record
    the resulting peak β_c^{(i)}.  Then
        σ² = (N-1)/N * Σ (β_c^{(i)} - <β_c>)².

    Cost is N extra `_gc_fit` calls (~1–5 ms each on ~30 points), so far
    below MC time.  Returns 0.0 if fewer than 6 jackknife samples succeed
    — the GC fit needs at least ~6 points to be well-posed.
    """
    n = len(scan_betas)
    if n < 7:
        return 0.0
    samples = []
    for i in range(n):
        bi = scan_betas[:i] + scan_betas[i + 1:]
        ci = scan_chis[:i] + scan_chis[i + 1:]
        _, beta_loo = _gc_fit(bi, ci, beta_init)
        if beta_loo is None or not np.isfinite(beta_loo):
            continue
        samples.append(float(beta_loo))
    if len(samples) < 6:
        return 0.0
    arr = np.asarray(samples, dtype=float)
    mean = arr.mean()
    var = (len(arr) - 1) / len(arr) * np.sum((arr - mean) ** 2)
    return float(np.sqrt(var))


def find_beta_c(exe, Lx, Ly, Tx, Ty, k1, k2, k3,
                beta_lo, beta_hi,
                n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
                n_traj_coarse=10000, n_traj_fine=30000,
                data_dir="/tmp/k_from_cs_scan",
                max_shifts=4, progress_cb=None, jackknife=False):
    """
    Locate susceptibility peak in [beta_lo, beta_hi] via 3-pass GC scan.
    Returns (beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis, scan_chi_errs).
    ``beta_c_sigma`` is the leave-one-out jackknife estimate when
    ``jackknife=True``, otherwise 0.0.

    If *progress_cb* is not None it is called after each GC pass with a dict:
        {"pass_num": 0-3, "all_betas": [...], "all_chis": [...],
         "pass_ids": [...], "gc_params": tuple|None, "beta_estimate": float}
    """
    scan_betas, scan_chis, scan_chi_errs = [], [], []
    scan_pass_ids = []
    ran = {}
    _current_pass = [0]
    _last_gc_params = [None]
    _last_beta_est = [None]
    _pass_counter = [0, 0]  # [done_in_pass, total_in_pass]
    _traj_done = [0]

    def _notify(pass_num, gc_params=None, beta_est=None, progress=None):
        if gc_params is not None:
            _last_gc_params[0] = gc_params
        if beta_est is not None:
            _last_beta_est[0] = beta_est
        if progress_cb is not None:
            progress_cb({
                "pass_num": pass_num,
                "all_betas": list(scan_betas),
                "all_chis": list(scan_chis),
                "all_chi_errs": list(scan_chi_errs),
                "pass_ids": list(scan_pass_ids),
                "gc_params": _last_gc_params[0],
                "beta_estimate": _last_beta_est[0],
                "scan_progress": progress,
                "traj_done": _traj_done[0],
            })

    def _run_beta(b, n_traj):
        key = round(b, 10)
        if key in ran:
            return ran[key]
        stdout, _ = run_simulator(exe, Lx, Ly, Tx, Ty, k1, k2, k3, b,
                                  n_traj=n_traj, data_dir=data_dir)
        chi, chi_err = parse_stdout_value_with_err(stdout, "m_susc:")
        scan_betas.append(b)
        scan_chis.append(chi)
        scan_chi_errs.append(chi_err)
        scan_pass_ids.append(_current_pass[0])
        ran[key] = (chi, chi_err)
        _pass_counter[0] += 1
        _traj_done[0] += n_traj
        _notify(_current_pass[0],
                progress={"pts_done": _pass_counter[0],
                          "pts_total": _pass_counter[1]})
        return chi, chi_err

    span = beta_hi - beta_lo

    # --- Coarse sweeps with window translation ---
    # If the susceptibility maximum lands at the edge of the current
    # bracket (or the chi values are monotonic), translate the window
    # by half its span in the appropriate direction and sweep n_coarse
    # *additional* points there.  Accumulated points across all shifts
    # feed the GC fit, so each translation strictly adds information.
    for shift_iter in range(max_shifts + 1):
        _pass_counter[:] = [0, n_coarse]
        for b in np.linspace(beta_lo, beta_hi, n_coarse):
            _run_beta(float(b), n_traj_coarse)

        order = np.argsort(scan_betas)
        sbetas = np.array(scan_betas)[order]
        schis = np.array(scan_chis)[order]
        idx_max = int(np.argmax(schis))
        # Edge detection: argmax within first/last 10% of accumulated points.
        n_pts = len(sbetas)
        edge_band = max(1, n_pts // 10)
        at_lo = (idx_max < edge_band)
        at_hi = (idx_max >= n_pts - edge_band)
        if not at_lo and not at_hi:
            # Found an interior peak in the accumulated sweep.
            print(f"[scan] interior peak found after {shift_iter + 1} sweep(s) "
                  f"({n_pts} pts in [{sbetas[0]:.4f}, {sbetas[-1]:.4f}])")
            break
        if shift_iter < max_shifts:
            # Translate the window by half its span in the direction of
            # the apparent peak.  This adds a fresh chunk of phase space
            # rather than backtracking over points we already have.
            shift = 0.5 * span
            if at_lo:
                beta_hi = beta_lo
                beta_lo = max(0.01, beta_lo - shift)
                print(f"[scan] no interior peak; translating window LEFT to "
                      f"[{beta_lo:.4f}, {beta_hi:.4f}] "
                      f"(shift {shift:.4f}, sweep {shift_iter + 2}/{max_shifts + 1})")
            else:
                beta_lo = beta_hi
                beta_hi = beta_hi + shift
                print(f"[scan] no interior peak; translating window RIGHT to "
                      f"[{beta_lo:.4f}, {beta_hi:.4f}] "
                      f"(shift {shift:.4f}, sweep {shift_iter + 2}/{max_shifts + 1})")
            sys.stdout.flush()
        else:
            print(f"[scan] WARNING: peak still at edge after {max_shifts + 1} "
                  f"window translations; proceeding with edge estimate.")

    # Re-sort
    order = np.argsort(scan_betas)
    sbetas = np.array(scan_betas)[order]
    idx_max = int(np.argmax(np.array(scan_chis)[order]))

    # Initial GC fit for centering
    beta_parabola = float(sbetas[idx_max])
    try:
        b_range = max(scan_betas) - min(scan_betas)
        ord_top = np.argsort(scan_chis)[::-1]
        top_n = min(7, len(scan_chis))
        betas_top = np.array(scan_betas)[ord_top[:top_n]]
        chis_top = np.array(scan_chis)[ord_top[:top_n]]
        p0 = [beta_parabola, b_range / 4, 0.0, 0.0, float(max(chis_top))]
        bounds = ([min(scan_betas), 1e-6, -5.0, -3.0, 0.0],
                  [max(scan_betas), b_range, 5.0, 10.0, np.inf])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(_gram_charlier, betas_top, chis_top,
                                p0=p0, bounds=bounds, maxfev=8000)
        beta_parabola = float(np.clip(popt[0], min(scan_betas), max(scan_betas)))
    except Exception:
        pass

    _notify(0, None, beta_parabola)
    _current_pass[0] = 1

    # ---------------------------------------------------------------
    # Helper: choose the half-width of the next pass.
    #
    # Higher passes scan around the current peak estimate over a window
    # whose half-width is a multiple of the *square root of the variance*
    # of the current Gram-Charlier fit (i.e. the fitted sigma).  This is
    # the natural physical width of the susceptibility peak — it matches
    # the data instead of being a geometric fraction of the original
    # bracket — and it shrinks automatically as the fit tightens.
    #
    # Fall back to a fraction of the original bracket if we have no
    # usable GC fit yet (e.g. pass 0 produced too few good points).
    # ---------------------------------------------------------------
    coarse_step = (max(scan_betas) - min(scan_betas)) / (2 * n_coarse)
    PASS_NSIGMA = {1: 2.0, 2: 1.0, 3: 0.5}

    def _half_width(pass_num, last_gc_params):
        if last_gc_params is not None:
            sigma_fit = float(last_gc_params[1])
            if np.isfinite(sigma_fit) and sigma_fit > 0:
                hw = PASS_NSIGMA[pass_num] * sigma_fit
                # Floor: never collapse below one MC step's worth so we
                # don't sample a degenerate window when sigma_fit is tiny
                # because of an under-constrained fit.
                return max(hw, coarse_step / (2 ** (pass_num + 1)))
        # Fallback: geometric halving of the original coarse step.
        return coarse_step / (2 ** (pass_num - 1))

    # --- Pass 1: ± n_sigma · sigma_fit around initial estimate -------
    # Pass 0 has not yet produced gc1_params; use the initial GC fit
    # captured above (popt) for the width estimate.
    last_gc = popt if 'popt' in locals() else None
    half = _half_width(1, last_gc)
    fine_betas = np.linspace(beta_parabola - half, beta_parabola + half, n_refine + 2)
    _pass_counter[:] = [0, len(fine_betas)]
    for b in fine_betas:
        _run_beta(float(b), n_traj_fine)

    gc1_params, beta_pass1 = _gc_fit(scan_betas, scan_chis, beta_parabola)
    _notify(1, gc1_params, beta_pass1)
    _current_pass[0] = 2

    # --- Pass 2: tighter bracket from gc1 sigma ----------------------
    half2 = _half_width(2, gc1_params)
    fine2_betas = np.linspace(beta_pass1 - half2, beta_pass1 + half2, n_refine2 + 2)
    _pass_counter[:] = [0, len(fine2_betas)]
    for b in fine2_betas:
        _run_beta(float(b), n_traj_fine)

    gc2_params, beta_c = _gc_fit(scan_betas, scan_chis, beta_pass1)
    _notify(2, gc2_params, beta_c)
    _current_pass[0] = 3

    # --- Pass 3: very tight bracket from gc2 sigma -------------------
    half3 = _half_width(3, gc2_params)
    fine3_betas = np.linspace(beta_c - half3, beta_c + half3, n_refine3 + 2)
    _pass_counter[:] = [0, len(fine3_betas)]
    for b in fine3_betas:
        _run_beta(float(b), n_traj_fine)

    gc3_params, beta_c_final = _gc_fit(scan_betas, scan_chis, beta_c)
    _notify(3, gc3_params, beta_c_final)

    chi_peak = max(scan_chis)

    # --- Optional jackknife uncertainty on β_c ------------------------
    if jackknife:
        beta_c_sigma = _gc_jackknife_sigma(scan_betas, scan_chis, beta_c_final)
    else:
        beta_c_sigma = 0.0

    return (beta_c_final, beta_c_sigma, chi_peak,
            scan_betas, scan_chis, scan_chi_errs)


# ===================================================================
# Boundary-slice cost functions
# ===================================================================

def _a2a_to_xy_arrays(data):
    keys = list(data.keys())
    m_arr = np.array([k[0] for k in keys], dtype=float)
    n_arr = np.array([k[1] for k in keys], dtype=float)
    c_arr = np.array([data[k]["conn"] for k in keys], dtype=float)
    x_arr = m_arr + 0.5 * n_arr
    y_arr = (math.sqrt(3.0) / 2.0) * n_arr
    return x_arr, y_arr, c_arr


def _tile_interp_field(data, Lx, Ly, Tx, Ty, copies=2, field="conn"):
    """Build a LinearNDInterpolator of *field* by tiling the torus."""
    keys = list(data.keys())
    m_arr = np.array([k[0] for k in keys], dtype=float)
    n_arr = np.array([k[1] for k in keys], dtype=float)
    c_arr = np.array([data[k][field] for k in keys], dtype=float)
    x0 = m_arr + 0.5 * n_arr
    sqrt3_2 = math.sqrt(3.0) / 2.0
    y0 = sqrt3_2 * n_arr
    vx = Lx + 0.5 * Ty;  vy = sqrt3_2 * Ty
    ux = Tx - 0.5 * Ly;  uy = -sqrt3_2 * Ly
    x_list, y_list, c_list = [x0], [y0], [c_arr]
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            if a == 0 and b == 0:
                continue
            x_list.append(x0 + a * vx + b * ux)
            y_list.append(y0 + a * vy + b * uy)
            c_list.append(c_arr)
    x_all = np.concatenate(x_list)
    y_all = np.concatenate(y_list)
    c_all = np.concatenate(c_list)
    pts = np.column_stack([x_all, y_all])
    _, uid = np.unique(np.round(pts, 8), axis=0, return_index=True)
    return LinearNDInterpolator(pts[uid], c_all[uid])


def _tile_interp_conn(data, Lx, Ly, Tx, Ty, copies=2):
    return _tile_interp_field(data, Lx, Ly, Tx, Ty, copies, "conn")


def _boundary_paths(Lx, Ly, Tx, Ty):
    return [
        (Lx, Ty),
        (Tx, -Ly),
        (-Lx - Tx, Ly - Ty),
    ]


def boundary_slices_normed(ref_data, test_data,
                           Lx, Ly, Tx, Ty,
                           ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                           copies=2, n_samples=400):
    """
    Variance-normalised Z² boundary-slice cost.
    Returns (total_cost, n_valid_paths, total_err).
    """
    rLx, rLy = ref_Lx, ref_Ly
    rTx, rTy = ref_Tx, ref_Ty

    iref = _tile_interp_conn(ref_data, rLx, rLy, rTx, rTy, copies)
    itest = _tile_interp_conn(test_data, Lx, Ly, Tx, Ty, copies)
    iref_err = _tile_interp_field(ref_data, rLx, rLy, rTx, rTy, copies, "conn_err")
    itest_err = _tile_interp_field(test_data, Lx, Ly, Tx, Ty, copies, "conn_err")

    ref_paths = _boundary_paths(rLx, rLy, rTx, rTy)
    test_paths = _boundary_paths(Lx, Ly, Tx, Ty)

    sqrt3_2 = math.sqrt(3.0) / 2.0
    n_valid = 0
    dir_labels = ["v", "u", "w"]
    per_dir = []
    t = np.linspace(0.0, 1.0, n_samples)
    for idx, ((rdm, rdn), (tdm, tdn)) in enumerate(zip(ref_paths, test_paths)):
        rex = rdm + 0.5 * rdn;  rey = sqrt3_2 * rdn
        pts_ref = np.column_stack([t * rex, t * rey])
        tex = tdm + 0.5 * tdn;  tey = sqrt3_2 * tdn
        pts_test = np.column_stack([t * tex, t * tey])

        cc_ref = np.array(iref(pts_ref), dtype=float)
        cc_test = np.array(itest(pts_test), dtype=float)
        ee_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))
        ee_test = np.abs(np.asarray(itest_err(pts_test), dtype=float))
        mask = (np.isfinite(cc_ref) & np.isfinite(cc_test)
                & np.isfinite(ee_ref) & np.isfinite(ee_test))
        if mask.sum() < 4:
            per_dir.append(0.0)
            continue
        diff = cc_test[mask] - cc_ref[mask]
        var = np.maximum(ee_ref[mask]**2 + ee_test[mask]**2, 1e-30)
        tm = t[mask]
        dt = tm[1:] - tm[:-1]
        integrand = 0.5 * (diff[:-1]**2 / var[:-1] + diff[1:]**2 / var[1:])
        per_dir.append(float(np.sum(integrand * dt)))
        n_valid += 1

    total = sum(c ** 2 for c in per_dir)
    if per_dir:
        parts = "  ".join(f"{l}={c:.6e}" for l, c in zip(dir_labels, per_dir))
        print(f"    [Z²] per-direction: {parts}")
    total_err = math.sqrt(2.0 * abs(total)) if total > 0 else 0.0
    return total, max(n_valid, 1), total_err


def boundary_slices_quartic(ref_data, test_data,
                            Lx, Ly, Tx, Ty,
                            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            copies=2, n_samples=400):
    """
    Variance-normalised Z⁴ boundary-slice cost.
    Returns (total_cost, n_valid_paths, total_err).
    """
    rLx, rLy = ref_Lx, ref_Ly
    rTx, rTy = ref_Tx, ref_Ty

    iref = _tile_interp_conn(ref_data, rLx, rLy, rTx, rTy, copies)
    itest = _tile_interp_conn(test_data, Lx, Ly, Tx, Ty, copies)
    iref_err = _tile_interp_field(ref_data, rLx, rLy, rTx, rTy, copies, "conn_err")
    itest_err = _tile_interp_field(test_data, Lx, Ly, Tx, Ty, copies, "conn_err")

    ref_paths = _boundary_paths(rLx, rLy, rTx, rTy)
    test_paths = _boundary_paths(Lx, Ly, Tx, Ty)

    sqrt3_2 = math.sqrt(3.0) / 2.0
    n_valid = 0
    dir_labels = ["v", "u", "w"]
    per_dir = []
    t = np.linspace(0.0, 1.0, n_samples)
    for idx, ((rdm, rdn), (tdm, tdn)) in enumerate(zip(ref_paths, test_paths)):
        rex = rdm + 0.5 * rdn;  rey = sqrt3_2 * rdn
        pts_ref = np.column_stack([t * rex, t * rey])
        tex = tdm + 0.5 * tdn;  tey = sqrt3_2 * tdn
        pts_test = np.column_stack([t * tex, t * tey])

        cc_ref = np.array(iref(pts_ref), dtype=float)
        cc_test = np.array(itest(pts_test), dtype=float)
        ee_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))
        ee_test = np.abs(np.asarray(itest_err(pts_test), dtype=float))
        mask = (np.isfinite(cc_ref) & np.isfinite(cc_test)
                & np.isfinite(ee_ref) & np.isfinite(ee_test))
        if mask.sum() < 4:
            per_dir.append(0.0)
            continue
        diff = cc_test[mask] - cc_ref[mask]
        var = np.maximum(ee_ref[mask]**2 + ee_test[mask]**2, 1e-30)
        z = diff / np.sqrt(var)
        z4 = z ** 4
        tm = t[mask]
        dt = tm[1:] - tm[:-1]
        integrand = 0.5 * (z4[:-1] + z4[1:])
        per_dir.append(float(np.sum(integrand * dt)))
        n_valid += 1

    total = sum(c ** 2 for c in per_dir)
    if per_dir:
        parts = "  ".join(f"{l}={c:.6e}" for l, c in zip(dir_labels, per_dir))
        print(f"    [Z⁴] per-direction: {parts}")
    total_err = math.sqrt(96.0 * abs(total) / max(n_samples, 1)) if total > 0 else 0.0
    return total, max(n_valid, 1), total_err


def extract_boundary_slices(ref_data, test_data, Lx, Ly, Tx, Ty,
                            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            copies=2, n_samples=200):
    """Extract per-direction correlator slices for plotting."""
    rLx, rLy = ref_Lx, ref_Ly
    rTx, rTy = ref_Tx, ref_Ty

    iref = _tile_interp_conn(ref_data, rLx, rLy, rTx, rTy, copies)
    itest = _tile_interp_conn(test_data, Lx, Ly, Tx, Ty, copies)
    iref_err = _tile_interp_field(ref_data, rLx, rLy, rTx, rTy, copies, "conn_err")
    itest_err = _tile_interp_field(test_data, Lx, Ly, Tx, Ty, copies, "conn_err")

    ref_paths = _boundary_paths(rLx, rLy, rTx, rTy)
    test_paths = _boundary_paths(Lx, Ly, Tx, Ty)

    sqrt3_2 = math.sqrt(3.0) / 2.0
    labels = ["v", "u", "w"]
    slices = []
    t = np.linspace(0.0, 1.0, n_samples)
    for i, ((rdm, rdn), (tdm, tdn)) in enumerate(zip(ref_paths, test_paths)):
        rex = rdm + 0.5 * rdn;  rey = sqrt3_2 * rdn
        pts_ref = np.column_stack([t * rex, t * rey])
        tex = tdm + 0.5 * tdn;  tey = sqrt3_2 * tdn
        pts_test = np.column_stack([t * tex, t * tey])

        cc_ref = np.asarray(iref(pts_ref), dtype=float)
        cc_test = np.asarray(itest(pts_test), dtype=float)
        ee_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))
        ee_test = np.abs(np.asarray(itest_err(pts_test), dtype=float))

        mask = (np.isfinite(cc_ref) & np.isfinite(cc_test) &
                np.isfinite(ee_ref) & np.isfinite(ee_test))
        slices.append({
            "t": t[mask],
            "g_ref": cc_ref[mask],
            "g_test": cc_test[mask],
            "diff": cc_test[mask] - cc_ref[mask],
            "diff_err": np.sqrt(ee_ref[mask]**2 + ee_test[mask]**2),
            "label": labels[i],
        })
    return slices
