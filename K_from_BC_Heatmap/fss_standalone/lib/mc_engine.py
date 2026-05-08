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
# Beta_c finder (Speedup 4: Ferrenberg–Swendsen reweighting)
# ===================================================================

def find_beta_c_reweight(exe, Lx, Ly, Tx, Ty, k1, k2, k3,
                         beta_lo, beta_hi,
                         n_traj_parent=40000,
                         n_grid=201,
                         n_eff_floor=0.10,
                         max_recenters=3,
                         data_dir="/tmp/k_from_cs_rw",
                         progress_cb=None,
                         jackknife=False):
    """
    Locate the susceptibility peak in [beta_lo, beta_hi] by running ONE
    long parent MC at the bracket midpoint and Ferrenberg–Swendsen
    reweighting χ(β) onto a dense β grid.

    Returns the same 6-tuple as :func:`find_beta_c`:
        (beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis, scan_chi_errs)

    where ``scan_betas`` / ``scan_chis`` are the dense reweighted grid
    (so downstream code that plots the scan still works) and
    ``scan_chi_errs`` is the jackknife σ on χ(β) at each grid point.

    If the reweighted peak lands too close to the edge of the validity
    window — defined by N_eff(β) ≥ n_eff_floor · N — the parent MC is
    re-spawned at the new estimate and the reweight repeats, up to
    ``max_recenters`` times.  This keeps the optimum away from the
    Boltzmann reweighting horizon without resorting to the legacy
    3-pass scan.
    """
    import reweight as rw_mod  # local import: avoids circulars at module load

    beta_parent = 0.5 * (beta_lo + beta_hi)
    beta_c_final = beta_parent
    beta_c_sigma = 0.0
    rw = None
    chis_grid = sigmas_grid = neff_grid = None
    grid = np.linspace(beta_lo, beta_hi, n_grid)

    for recenter_iter in range(max_recenters + 1):
        # --- Spawn one parent ensemble at beta_parent with traces dump ---
        # Note: the C++ binary supports --dump_traces 1 (S4 patch).
        seed = int(time.time() * 1000) & 0xFFFFFFFF
        cmd = [
            exe,
            "--L_x", str(Lx), "--L_y", str(Ly),
            "--T_x", str(Tx), "--T_y", str(Ty),
            "--k1", f"{k1:.12f}", "--k2", f"{k2:.12f}", "--k3", f"{k3:.12f}",
            "--beta", f"{beta_parent:.12f}",
            "--n_traj", str(n_traj_parent), "--n_skip", "20", "--n_therm", "2000",
            "--seed", str(seed),
            "--data_dir", data_dir,
            "--dump_traces", "1",
        ]
        os.makedirs(data_dir, exist_ok=True)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
        if result.returncode != 0:
            raise RuntimeError(
                f"reweight: parent simulator failed:\n{result.stderr}"
            )
        # Find the subdir the simulator wrote to.
        subdir = None
        for line in result.stdout.splitlines():
            if line.startswith("opening file:"):
                subdir = os.path.dirname(line.split(":", 1)[1].strip())
                break
        if subdir is None:
            raise RuntimeError("reweight: could not find simulator output subdir")
        traces_path = rw_mod.find_traces_file(subdir)
        if traces_path is None:
            raise RuntimeError(
                f"reweight: no traces_*.dat in {subdir} — rebuild simulator?"
            )

        # --- Reweight χ(β) over the dense grid ---
        rw = rw_mod.Reweighter(traces_path, n_eff_floor=n_eff_floor)
        chis_grid, sigmas_grid, neff_grid = rw.chi_grid(grid)

        valid = np.isfinite(chis_grid)
        if not valid.any():
            raise RuntimeError(
                f"reweight: N_eff < {n_eff_floor} on the entire bracket "
                f"[{beta_lo:.4f}, {beta_hi:.4f}]; widen beta_seed."
            )

        # --- GC fit on the reweighted (β, χ) curve ---
        valid_betas = grid[valid].tolist()
        valid_chis = chis_grid[valid].tolist()
        gc_params, beta_c_est = _gc_fit(valid_betas, valid_chis, beta_parent)

        if progress_cb is not None:
            progress_cb({
                "pass_num": recenter_iter,
                "all_betas": list(grid),
                "all_chis": list(chis_grid),
                "all_chi_errs": list(sigmas_grid),
                "pass_ids": [recenter_iter] * len(grid),
                "gc_params": gc_params,
                "beta_estimate": beta_c_est,
                "scan_progress": {"pts_done": int(valid.sum()),
                                  "pts_total": n_grid},
                "traj_done": n_traj_parent * (recenter_iter + 1),
                "n_eff_min": float(np.nanmin(neff_grid / max(rw.n_samples, 1))),
            })

        # --- Validity gates ---------------------------------------
        # (A) Edge gate: did the χ(β) maximum (or GC peak) land at the
        #     boundary of the current bracket?  If so the *true* peak
        #     is almost certainly outside the window — translate the
        #     bracket in that direction by half its span (à la
        #     find_beta_c) regardless of what N_eff says, since N_eff
        #     is happy at the parent centre but tells us nothing about
        #     phase space we have not yet visited.
        # (B) N_eff gate: even if the peak is interior, require
        #     N_eff ≥ floor at the peak and at ±1 GC-σ around it.
        # ---------------------------------------------------------
        idx_argmax = int(np.nanargmax(np.where(valid, chis_grid, -np.inf)))
        edge_band = max(1, n_grid // 10)
        chi_at_lo = (idx_argmax < edge_band)
        chi_at_hi = (idx_argmax >= n_grid - edge_band)
        gc_at_lo  = (gc_params is not None and
                     beta_c_est <= grid[edge_band - 1] + 1e-12)
        gc_at_hi  = (gc_params is not None and
                     beta_c_est >= grid[n_grid - edge_band] - 1e-12)
        at_lo_edge = chi_at_lo or gc_at_lo
        at_hi_edge = chi_at_hi or gc_at_hi

        if at_lo_edge or at_hi_edge:
            if recenter_iter == max_recenters:
                print(f"[rw] WARNING: peak still at bracket edge after "
                      f"{max_recenters} translations; using edge estimate "
                      f"β={beta_c_est:.5f}.")
                beta_c_final = beta_c_est
                break
            # Translate bracket by half its span in the direction of
            # the peak.  Keep the same half-width so we sweep fresh
            # phase space symmetrically about the new parent.
            span = beta_hi - beta_lo
            shift = 0.5 * span
            if at_lo_edge:
                beta_hi = beta_lo
                beta_lo = max(0.01, beta_lo - shift)
                direction = "LEFT"
            else:
                beta_lo = beta_hi
                beta_hi = beta_hi + shift
                direction = "RIGHT"
            beta_parent = 0.5 * (beta_lo + beta_hi)
            grid = np.linspace(beta_lo, beta_hi, n_grid)
            print(f"[rw] translate {direction} → "
                  f"β∈[{beta_lo:.4f},{beta_hi:.4f}] "
                  f"(parent={beta_parent:.4f}, "
                  f"recenter {recenter_iter + 1}/{max_recenters})")
            continue

        if gc_params is None:
            beta_c_final = beta_c_est
            break
        sigma_fit = float(gc_params[1])
        peak_lo = beta_c_est - sigma_fit
        peak_hi = beta_c_est + sigma_fit
        # Find nearest-grid N_eff at the peak / ±1σ
        def _neff_frac(b):
            i = int(np.clip(np.searchsorted(grid, b), 0, n_grid - 1))
            return float(neff_grid[i] / max(rw.n_samples, 1))
        ok_centre = _neff_frac(beta_c_est) >= n_eff_floor
        ok_lo = _neff_frac(peak_lo) >= n_eff_floor
        ok_hi = _neff_frac(peak_hi) >= n_eff_floor
        if ok_centre and ok_lo and ok_hi:
            beta_c_final = beta_c_est
            break
        if recenter_iter == max_recenters:
            print(f"[rw] WARNING: peak outside N_eff window after "
                  f"{max_recenters} recenters; using current estimate.")
            beta_c_final = beta_c_est
            break
        # Recenter parent on the new estimate; tighten bracket to ±2σ.
        new_half = max(2.0 * sigma_fit, (beta_hi - beta_lo) * 0.25)
        beta_parent = float(beta_c_est)
        beta_lo = max(0.01, beta_parent - new_half)
        beta_hi = beta_parent + new_half
        grid = np.linspace(beta_lo, beta_hi, n_grid)
        print(f"[rw] recentre {recenter_iter + 1}: β_parent={beta_parent:.5f} "
              f"window±{new_half:.4f}")

    # --- Optional jackknife on β_c using the final reweighter ---
    if jackknife and rw is not None:
        # Re-evaluate β_c by leaving out blocks of the parent histogram
        # and refitting the GC peak.  Reuses Reweighter._n_blocks.
        B = rw._n_blocks
        block = rw._n // B
        peaks = []
        valid_grid = grid[np.isfinite(chis_grid)]
        for b in range(B):
            lo = b * block
            hi = (b + 1) * block if b < B - 1 else rw._n
            mask = np.ones(rw._n, dtype=bool)
            mask[lo:hi] = False
            sub = {
                "E_per_site": rw._E_total[mask] / rw._n_sites,
                "abs_m": rw._abs_m[mask],
                "m2": rw._m2[mask],
                "n_sites": rw._n_sites,
                "beta_parent": rw._beta,
                "K": rw._K,
                "path": "<jackknife block>",
            }
            try:
                rw_b = rw_mod.Reweighter(sub, n_eff_floor=n_eff_floor)
            except ValueError:
                continue
            cb, _, _ = rw_b.chi_grid(valid_grid)
            ok = np.isfinite(cb)
            if ok.sum() < 6:
                continue
            _, peak_b = _gc_fit(valid_grid[ok].tolist(),
                                cb[ok].tolist(), beta_c_final)
            if peak_b is not None and np.isfinite(peak_b):
                peaks.append(float(peak_b))
        if len(peaks) >= 6:
            arr = np.asarray(peaks)
            beta_c_sigma = float(
                np.sqrt((len(arr) - 1) / len(arr) * np.sum((arr - arr.mean()) ** 2))
            )

    chi_peak = float(np.nanmax(chis_grid))

    return (float(beta_c_final), float(beta_c_sigma), chi_peak,
            grid.tolist(),
            [float(c) if np.isfinite(c) else float("nan") for c in chis_grid],
            [float(s) if np.isfinite(s) else float("nan") for s in sigmas_grid])


# ===================================================================
# Multi-donor reweighting
# ===================================================================

def _spawn_parent_traces(exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_parent,
                         n_traj, data_dir, *, n_skip=20, n_therm=2000,
                         seed=None, timeout=1200):
    """Run one MC parent at *beta_parent* with traces dump.

    Returns the path to ``traces_*.dat`` in the simulator's output
    subdirectory.  Raises RuntimeError on simulator failure.
    """
    import reweight as rw_mod  # local import (avoid circulars)
    if seed is None:
        seed = int(time.time() * 1_000_000) & 0xFFFFFFFF
    cmd = [
        exe,
        "--L_x", str(Lx), "--L_y", str(Ly),
        "--T_x", str(Tx), "--T_y", str(Ty),
        "--k1", f"{k1:.12f}", "--k2", f"{k2:.12f}", "--k3", f"{k3:.12f}",
        "--beta", f"{beta_parent:.12f}",
        "--n_traj", str(n_traj), "--n_skip", str(n_skip),
        "--n_therm", str(n_therm),
        "--seed", str(seed),
        "--data_dir", data_dir,
        "--dump_traces", "1",
    ]
    os.makedirs(data_dir, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(
            f"multidonor: simulator failed at β={beta_parent:.5f}:\n"
            f"{result.stderr}"
        )
    subdir = None
    for line in result.stdout.splitlines():
        if line.startswith("opening file:"):
            subdir = os.path.dirname(line.split(":", 1)[1].strip())
            break
    if subdir is None:
        raise RuntimeError(
            f"multidonor: could not parse simulator subdir at β={beta_parent:.5f}"
        )
    traces_path = rw_mod.find_traces_file(subdir)
    if traces_path is None:
        raise RuntimeError(
            f"multidonor: no traces_*.dat in {subdir} (rebuild simulator?)"
        )
    return traces_path


def find_beta_c_multidonor(exe, Lx, Ly, Tx, Ty, k1, k2, k3,
                           beta_lo, beta_hi,
                           n_traj_parent=40000,
                           n_grid=201,
                           n_eff_floor=0.10,
                           data_dir="/tmp/k_from_cs_md",
                           progress_cb=None,
                           jackknife=False,
                           *,
                           var_E=None,
                           n_donors=None,
                           donor_overlap_alpha=0.75,
                           pilot_n_traj=2000,
                           min_donors=1,
                           max_donors=16):
    """Locate the susceptibility peak in [beta_lo, beta_hi] using
    multiple parent histograms tiled across the bracket.

    The donor spacing is set by the action-variance heuristic

        |Δβ|_safe ≈ sqrt( -ln(n_eff_floor) / Var(E_extensive) )

    so that adjacent parents overlap with N_eff/N ≥ ``n_eff_floor``.
    Donors are placed at ``donor_overlap_alpha · |Δβ|_safe`` apart
    (alpha < 1 enforces overlap).

    Total MC budget ``n_traj_parent`` is split evenly across donors.

    Parameters
    ----------
    var_E : float or None
        If provided, skip the pilot run.  This is the extensive-energy
        variance ``Var(E_total)`` near the bracket midpoint, typically
        cached from the reference build.
    n_donors : int or None
        If given, override the auto-determined donor count.
    donor_overlap_alpha : float
        Donor-spacing multiplier (in units of |Δβ|_safe).  0.75 gives
        ~50%% N_eff overlap between neighbours.
    pilot_n_traj : int
        Trajectories for the pilot run when var_E is unknown.

    Returns the same 6-tuple as :func:`find_beta_c_reweight`:
        (beta_c, beta_c_sigma, chi_peak,
         scan_betas, scan_chis, scan_chi_errs)
    """
    import reweight as rw_mod  # local import (circular-import guard)

    beta_mid = 0.5 * (beta_lo + beta_hi)
    grid = np.linspace(beta_lo, beta_hi, n_grid)

    # ---------------------------------------------------------------
    # 1. Estimate Var(E_extensive) — cheap pilot if not cached.
    # ---------------------------------------------------------------
    if var_E is None:
        pilot_path = _spawn_parent_traces(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_mid,
            n_traj=pilot_n_traj,
            data_dir=os.path.join(data_dir, "pilot"),
        )
        pilot_rw = rw_mod.Reweighter(pilot_path, n_eff_floor=n_eff_floor)
        var_E = float(np.var(pilot_rw._E_total))
    var_E = max(var_E, 1e-30)

    # ---------------------------------------------------------------
    # 2. Determine donor count and spacing.
    # ---------------------------------------------------------------
    delta_safe = math.sqrt(-math.log(max(n_eff_floor, 1e-12)) / var_E)
    span = float(beta_hi - beta_lo)
    if n_donors is None:
        # span / spacing + 1, with overlap factor alpha
        spacing = max(donor_overlap_alpha * delta_safe, 1e-12)
        n_auto = int(math.ceil(span / spacing)) + 1
        n_donors = int(np.clip(n_auto, min_donors, max_donors))
    n_donors = max(1, int(n_donors))

    if n_donors == 1:
        donor_betas = np.array([beta_mid])
    else:
        donor_betas = np.linspace(beta_lo, beta_hi, n_donors)

    n_traj_per = max(500, int(n_traj_parent // n_donors))

    if progress_cb is not None:
        progress_cb({
            "pass_num": 0,
            "all_betas": [], "all_chis": [], "all_chi_errs": [],
            "pass_ids": [],
            "gc_params": None, "beta_estimate": beta_mid,
            "scan_progress": {"pts_done": 0, "pts_total": n_grid},
            "traj_done": 0,
            "n_eff_min": 0.0,
            "multidonor": {
                "n_donors": int(n_donors),
                "delta_safe": float(delta_safe),
                "var_E": float(var_E),
                "n_traj_per": int(n_traj_per),
                "donor_betas": [float(b) for b in donor_betas],
            },
        })

    # ---------------------------------------------------------------
    # 3. Spawn donors and build a CombinedReweighter.
    # ---------------------------------------------------------------
    rws = []
    traj_done = 0
    for i, b_par in enumerate(donor_betas):
        sub = os.path.join(data_dir, f"donor_{i:02d}")
        traces_path = _spawn_parent_traces(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, float(b_par),
            n_traj=n_traj_per, data_dir=sub,
        )
        rws.append(rw_mod.Reweighter(traces_path, n_eff_floor=n_eff_floor))
        traj_done += n_traj_per
        if progress_cb is not None:
            progress_cb({
                "pass_num": 1,
                "all_betas": [], "all_chis": [], "all_chi_errs": [],
                "pass_ids": [],
                "gc_params": None, "beta_estimate": beta_mid,
                "scan_progress": {"pts_done": i + 1, "pts_total": n_donors},
                "traj_done": traj_done,
                "n_eff_min": 0.0,
            })

    combo = rw_mod.CombinedReweighter(rws, n_eff_floor=n_eff_floor)
    chis_grid, sigmas_grid, neff_grid = combo.chi_grid(grid)

    valid = np.isfinite(chis_grid)
    if not valid.any():
        raise RuntimeError(
            f"multidonor: no donor passed N_eff floor on bracket "
            f"[{beta_lo:.4f}, {beta_hi:.4f}] (n_donors={n_donors}); "
            f"increase n_traj_parent or widen donor_overlap_alpha."
        )

    # ---------------------------------------------------------------
    # 4. GC fit on the combined (β, χ) curve.
    # ---------------------------------------------------------------
    valid_betas = grid[valid].tolist()
    valid_chis = chis_grid[valid].tolist()
    gc_params, beta_c_est = _gc_fit(valid_betas, valid_chis, beta_mid)
    beta_c_final = float(beta_c_est) if beta_c_est is not None else beta_mid

    if progress_cb is not None:
        progress_cb({
            "pass_num": 2,
            "all_betas": list(grid),
            "all_chis": list(chis_grid),
            "all_chi_errs": list(sigmas_grid),
            "pass_ids": [2] * len(grid),
            "gc_params": gc_params,
            "beta_estimate": beta_c_final,
            "scan_progress": {"pts_done": int(valid.sum()),
                              "pts_total": n_grid},
            "traj_done": traj_done,
            "n_eff_min": float(np.nanmin([
                neff_grid[i] / max(rws[0]._n, 1) for i in range(len(grid))
            ])),
        })

    # ---------------------------------------------------------------
    # 5. Optional jackknife on β_c (per-donor leave-one-out).
    # ---------------------------------------------------------------
    beta_c_sigma = 0.0
    if jackknife and len(rws) >= 3:
        peaks = []
        for j in range(len(rws)):
            sub_rws = [r for k, r in enumerate(rws) if k != j]
            sub_combo = rw_mod.CombinedReweighter(sub_rws,
                                                  n_eff_floor=n_eff_floor)
            cb, _, _ = sub_combo.chi_grid(grid)
            ok = np.isfinite(cb)
            if ok.sum() < 6:
                continue
            _, peak_b = _gc_fit(grid[ok].tolist(), cb[ok].tolist(),
                                beta_c_final)
            if peak_b is not None and np.isfinite(peak_b):
                peaks.append(float(peak_b))
        if len(peaks) >= 3:
            arr = np.asarray(peaks)
            B = len(arr)
            beta_c_sigma = float(
                math.sqrt((B - 1) / B * float(np.sum((arr - arr.mean()) ** 2)))
            )

    chi_peak = float(np.nanmax(chis_grid))

    return (float(beta_c_final), float(beta_c_sigma), chi_peak,
            grid.tolist(),
            [float(c) if np.isfinite(c) else float("nan") for c in chis_grid],
            [float(s) if np.isfinite(s) else float("nan") for s in sigmas_grid])


def find_beta_c_multidonor_2pass(exe, Lx, Ly, Tx, Ty, k1, k2, k3,
                                 beta_lo, beta_hi,
                                 n_traj_parent=40000,
                                 n_grid=201,
                                 n_eff_floor=0.10,
                                 data_dir="/tmp/k_from_cs_md2",
                                 progress_cb=None,
                                 jackknife=False,
                                 *,
                                 var_E=None,
                                 pass1_n_donors=4,
                                 pass1_budget_frac=0.30,
                                 donor_overlap_alpha=0.75,
                                 pilot_n_traj=2000,
                                 max_donors=16):
    """Two-pass multi-donor reweighting.

    Pass 1
        ``pass1_n_donors`` donors evenly spaced across [beta_lo, beta_hi]
        using ``pass1_budget_frac`` of the total trajectory budget.
        A GC fit on the combined chi(beta) locates the approximate peak
        beta* and identifies which gap [p1_betas[i], p1_betas[i+1]]
        contains it.

    Pass 2
        Dense donors are placed over the *interval-based* window:
        the gap that contains beta* PLUS one neighbouring gap on each
        side.  Concretely, if beta* falls in gap i (0-based), the
        Pass-2 bracket is

            p2_lo = p1_betas[max(0, i-1)]
            p2_hi = p1_betas[min(N-1, i+2)]

        so Pass-2 always spans at most 3 consecutive Pass-1 gaps
        (fewer at the bracket edges).  Donor spacing within that window
        is set by the action-variance heuristic
        |Delta_beta|_safe = sqrt(-ln(n_eff_floor) / Var(E)).
        All Pass-1 donors are kept in the combined reweighter so the
        chi(beta) wings outside the refinement window stay covered.

    The final GC fit is restricted to the Pass-2 window for sharper
    peak resolution.

    Falls back gracefully if the Pass-1 GC fit fails: uses the raw
    chi(beta) argmax to choose the gap.

    Returns the same 6-tuple as :func:`find_beta_c_reweight`:
        (beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis,
         scan_chi_errs)

    Next steps / planned augmentations
    ------------------------------------
    1. Adaptive pass1_n_donors: set N so that the Pass-1 gap width
       matches the action-variance delta_safe from a cheap pilot run
       (``pilot_n_traj`` trajectories at beta_mid), giving guaranteed
       single-gap resolution without user tuning.
    2. Budget recycling: after the Pass-2 window is chosen, redistribute
       any over-allocated Pass-1 budget (from gaps outside the window)
       to the Pass-2 donors, keeping the total trajectory count fixed.
    3. Three-pass extension: after Pass-2 the function already knows a
       much tighter beta_c estimate; add an optional Pass-3 that places
       a handful of donors within a single delta_safe of beta_c for
       sub-percent beta_c uncertainty at large lattice sizes.
    4. Wire into prod_runtime.py: expose via a ``multidonor_2pass``
       CONFIG flag and ``--multidonor-2pass`` CLI switch in run.py,
       slotting in as the primary reweighting path ahead of the current
       single-pass multidonor.
    5. Test-script extension: add the 2-pass method as a third curve in
       test_multidonor_vs_3pass.py so all three methods appear on the
       same plot with a shared budget.
    """
    import reweight as rw_mod

    beta_mid = 0.5 * (beta_lo + beta_hi)
    span = float(beta_hi - beta_lo)
    grid = np.linspace(beta_lo, beta_hi, n_grid)

    # Budget split.
    pass1_budget_frac = float(np.clip(pass1_budget_frac, 0.05, 0.9))
    p1_total = max(2000, int(n_traj_parent * pass1_budget_frac))
    p2_total = max(2000, n_traj_parent - p1_total)

    # ---------------------------------------------------------------
    # PASS 1 — sparse coarse tile, with bracket translation
    #
    # If the χ(β) maximum lands at the edge of [beta_lo, beta_hi] the
    # true peak is outside the window.  Translate by half-span in the
    # appropriate direction (à la find_beta_c) and resweep, accumulating
    # donors so wing coverage strictly grows.
    # ---------------------------------------------------------------
    pass1_n_donors = max(2, int(pass1_n_donors))
    rws_p1 = []
    p1_beta_list = []         # ALL pass-1 donor betas across translations
    p1_per   = max(500, p1_total // pass1_n_donors)
    traj_done = 0

    # Trajectory budget is fixed; allow up to MAX_TRANSLATIONS extra
    # half-span hops without inflating the per-donor sample count.
    MAX_TRANSLATIONS = 6
    for shift_iter in range(MAX_TRANSLATIONS + 1):
        new_betas = np.linspace(beta_lo, beta_hi, pass1_n_donors)
        for i, b_par in enumerate(new_betas):
            sub = os.path.join(data_dir,
                               f"p1_shift{shift_iter:02d}_donor_{i:02d}")
            traces_path = _spawn_parent_traces(
                exe, Lx, Ly, Tx, Ty, k1, k2, k3, float(b_par),
                n_traj=p1_per, data_dir=sub,
            )
            rws_p1.append(rw_mod.Reweighter(traces_path,
                                            n_eff_floor=n_eff_floor))
            p1_beta_list.append(float(b_par))
            traj_done += p1_per
            if progress_cb is not None:
                progress_cb({
                    "pass_num": 1, "all_betas": [], "all_chis": [],
                    "all_chi_errs": [], "pass_ids": [],
                    "gc_params": None, "beta_estimate": beta_mid,
                    "scan_progress": {"pts_done": i + 1,
                                      "pts_total": pass1_n_donors},
                    "traj_done": traj_done, "n_eff_min": 0.0,
                })

        # Re-grid so the χ check spans every donor accumulated so far.
        full_lo = min(min(p1_beta_list), beta_lo)
        full_hi = max(max(p1_beta_list), beta_hi)
        grid = np.linspace(full_lo, full_hi, n_grid)
        combo1 = rw_mod.CombinedReweighter(rws_p1, n_eff_floor=n_eff_floor)
        chis1, sigmas1, _ = combo1.chi_grid(grid)
        valid1 = np.isfinite(chis1)
        if not valid1.any():
            raise RuntimeError(
                f"multidonor_2pass: pass-1 has no valid points on "
                f"[{full_lo:.4f}, {full_hi:.4f}]; widen bracket or "
                f"increase pass1_budget_frac."
            )

        # Edge detection on the χ(β) profile.
        chis_v = np.where(valid1, chis1, -np.inf)
        idx_max = int(np.nanargmax(chis_v))
        edge_band = max(1, n_grid // 10)
        at_lo = (idx_max < edge_band)
        at_hi = (idx_max >= n_grid - edge_band)
        if not at_lo and not at_hi:
            break  # interior peak — proceed to Pass 2
        if shift_iter == MAX_TRANSLATIONS:
            print(f"[md2] WARNING: pass-1 peak still at edge of "
                  f"[{full_lo:.4f}, {full_hi:.4f}] after "
                  f"{MAX_TRANSLATIONS} translations; using edge.")
            break
        shift = 0.5 * span
        if at_lo:
            beta_hi = beta_lo
            beta_lo = max(0.01, beta_lo - shift)
            direction = "LEFT"
        else:
            beta_lo = beta_hi
            beta_hi = beta_hi + shift
            direction = "RIGHT"
        beta_mid = 0.5 * (beta_lo + beta_hi)
        print(f"[md2] pass-1 peak at edge — translate {direction} → "
              f"adding donors in [{beta_lo:.4f}, {beta_hi:.4f}] "
              f"(shift {shift_iter + 1}/{MAX_TRANSLATIONS})")

    # Final pass-1 grid spans every donor accumulated.
    p1_betas = np.array(sorted(set(p1_beta_list)))
    pass1_n_donors = len(p1_betas)
    full_lo = float(grid[0])
    full_hi = float(grid[-1])
    span = full_hi - full_lo
    beta_mid = 0.5 * (full_lo + full_hi)
    gc_p1, beta_est_p1 = _gc_fit(grid[valid1].tolist(),
                                 chis1[valid1].tolist(), beta_mid)

    # ---------------------------------------------------------------
    # Pass-2 window: the Pass-1 gap that contains β* plus one
    # neighbour on each side.
    #
    # With pass1_n_donors = N donors there are N-1 gaps.  If β* falls
    # in gap i (0-based), the window is
    #   p2_lo = p1_betas[max(0, i-1)]
    #   p2_hi = p1_betas[min(N-1, i+2)]
    # so Pass-2 spans at most 3 consecutive Pass-1 gaps.
    # ---------------------------------------------------------------
    if gc_p1 is not None and beta_est_p1 is not None and np.isfinite(beta_est_p1):
        beta_star = float(beta_est_p1)
    else:
        # Pass 1 GC failed — fall back to raw chi argmax.
        idx = int(np.nanargmax(np.where(valid1, chis1, -np.inf)))
        beta_star = float(grid[idx])

    p1_arr = np.asarray(p1_betas)
    # right_idx: index of first Pass-1 donor strictly right of β*
    right_idx = int(np.searchsorted(p1_arr, beta_star, side="right"))
    # gap_i: index of interval [p1_arr[gap_i], p1_arr[gap_i+1]]
    gap_i = int(np.clip(right_idx - 1, 0, len(p1_arr) - 2))
    lo_idx = max(0, gap_i - 1)
    hi_idx = min(len(p1_arr) - 1, gap_i + 2)
    p2_lo = float(p1_arr[lo_idx])
    p2_hi = float(p1_arr[hi_idx])

    # Action-variance heuristic governs donor *spacing* within the window.
    if var_E is None:
        nearest = min(rws_p1, key=lambda r: abs(r._beta - beta_star))
        var_E = float(np.var(nearest._E_total))
    var_E = max(var_E, 1e-30)
    delta_safe = math.sqrt(-math.log(max(n_eff_floor, 1e-12)) / var_E)

    # Extend the Pass-2 window by one delta_safe on each side so that
    # Pass-2 donors reach into the Pass-1 wing coverage, eliminating the
    # N_eff gap (and resulting chi shear) at the p2_lo / p2_hi boundaries.
    overlap = donor_overlap_alpha * delta_safe
    p2_lo = max(float(beta_lo), p2_lo - overlap)
    p2_hi = min(float(beta_hi), p2_hi + overlap)

    spacing = max(donor_overlap_alpha * delta_safe, 1e-12)
    n_p2 = int(math.ceil((p2_hi - p2_lo) / spacing)) + 1
    n_p2 = int(np.clip(n_p2, 2, max_donors))
    p2_betas = np.linspace(p2_lo, p2_hi, n_p2)
    p2_per   = max(500, p2_total // n_p2)

    # ---------------------------------------------------------------
    # PASS 2 — dense refinement around β*
    # ---------------------------------------------------------------
    rws_p2 = []
    for i, b_par in enumerate(p2_betas):
        sub = os.path.join(data_dir, f"p2_donor_{i:02d}")
        traces_path = _spawn_parent_traces(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, float(b_par),
            n_traj=p2_per, data_dir=sub,
        )
        rws_p2.append(rw_mod.Reweighter(traces_path, n_eff_floor=n_eff_floor))
        traj_done += p2_per
        if progress_cb is not None:
            progress_cb({
                "pass_num": 2, "all_betas": [], "all_chis": [],
                "all_chi_errs": [], "pass_ids": [],
                "gc_params": gc_p1, "beta_estimate": beta_star,
                "scan_progress": {"pts_done": i + 1, "pts_total": n_p2},
                "traj_done": traj_done, "n_eff_min": 0.0,
            })

    # Combine ALL donors (both passes) — pass-1 covers the wings,
    # pass-2 nails the peak.
    rws_all = rws_p1 + rws_p2
    combo = rw_mod.CombinedReweighter(rws_all, n_eff_floor=n_eff_floor)
    chis_grid, sigmas_grid, neff_grid = combo.chi_grid(grid)
    valid = np.isfinite(chis_grid)
    if not valid.any():
        raise RuntimeError(
            "multidonor_2pass: combined coverage failed N_eff floor; "
            "increase budget or n_eff_floor."
        )

    # Final GC fit on the combined curve, restricted to the pass-2
    # window when possible (gives a sharper peak fit).
    fit_mask = valid & (grid >= p2_lo) & (grid <= p2_hi)
    if fit_mask.sum() >= 6:
        gc_params, beta_c_est = _gc_fit(
            grid[fit_mask].tolist(),
            chis_grid[fit_mask].tolist(),
            beta_star,
        )
    else:
        gc_params, beta_c_est = _gc_fit(
            grid[valid].tolist(),
            chis_grid[valid].tolist(),
            beta_star,
        )
    beta_c_final = float(beta_c_est) if beta_c_est is not None else beta_star

    # Build pass_ids for the combined grid: points inside the Pass-2
    # refinement window get id=1; Pass-1 wings get id=0.
    _pass_ids = [1 if p2_lo <= float(b) <= p2_hi else 0 for b in grid]
    if progress_cb is not None:
        progress_cb({
            "pass_num": 3,
            "all_betas": list(grid),
            "all_chis": list(chis_grid),
            "all_chi_errs": list(sigmas_grid),
            "pass_ids": _pass_ids,
            "gc_params": gc_params, "beta_estimate": beta_c_final,
            "p2_lo": float(p2_lo), "p2_hi": float(p2_hi),
            "scan_progress": {"pts_done": int(valid.sum()),
                              "pts_total": n_grid},
            "traj_done": traj_done, "n_eff_min": 0.0,
        })

    # Per-donor jackknife on β_c (leave-one-out across pass-2 donors).
    beta_c_sigma = 0.0
    if jackknife and len(rws_p2) >= 3:
        peaks = []
        for j in range(len(rws_p2)):
            sub_p2 = [r for k, r in enumerate(rws_p2) if k != j]
            sub_combo = rw_mod.CombinedReweighter(rws_p1 + sub_p2,
                                                  n_eff_floor=n_eff_floor)
            cb, _, _ = sub_combo.chi_grid(grid)
            ok = np.isfinite(cb) & (grid >= p2_lo) & (grid <= p2_hi)
            if ok.sum() < 6:
                continue
            _, peak_b = _gc_fit(grid[ok].tolist(), cb[ok].tolist(),
                                beta_c_final)
            if peak_b is not None and np.isfinite(peak_b):
                peaks.append(float(peak_b))
        if len(peaks) >= 3:
            arr = np.asarray(peaks)
            B = len(arr)
            beta_c_sigma = float(
                math.sqrt((B - 1) / B * float(np.sum((arr - arr.mean()) ** 2)))
            )

    chi_peak = float(np.nanmax(chis_grid))

    return (float(beta_c_final), float(beta_c_sigma), chi_peak,
            grid.tolist(),
            [float(c) if np.isfinite(c) else float("nan") for c in chis_grid],
            [float(s) if np.isfinite(s) else float("nan") for s in sigmas_grid])


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
