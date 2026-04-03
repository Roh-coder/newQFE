#!/usr/bin/env python3
"""
optimise_couplings.py — Find critical couplings by matching all-to-all
two-point functions between an anisotropic lattice and an equilateral reference.

Two-level optimisation:
  Outer loop: adaptive grid search over (r1, r2) = (k1/k3, k2/k3)
    - Evaluate chi2 on a 4x4 grid
    - If the minimum is in the interior, refine the grid by 2x around it
    - If the minimum is on the border, translate the grid in that direction
  Inner loop: susceptibility-peak scan over beta to find beta_c

A live matplotlib figure updates after each refinement level showing:
  - Panel 1: grid points in (r1, r2) space, coloured by chi2
  - Panel 2: chi2 at best grid point vs refinement level
  - Panel 3: beta_c at best point vs refinement level
  - Panel 4: susceptibility vs beta for the latest inner-loop scan
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time

import matplotlib
# Only force non-interactive Agg when not running inside IPython/Spyder.
# Inside IPython the configured backend (inline, qt5, etc.) is used instead,
# and the figure is pushed to the console after every update.
try:
    from IPython import get_ipython as _gip
    _in_ipython = (_gip() is not None)
except ImportError:
    _in_ipython = False
if not _in_ipython:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_stdout_value(stdout: str, key: str) -> float:
    """Extract a numeric value from simulator stdout like 'm_susc: 1.23e-02 4.56e-03'."""
    for line in stdout.splitlines():
        if line.startswith(key):
            parts = line.split()
            return float(parts[1])
    raise ValueError(f"Key '{key}' not found in simulator output")


def parse_stdout_value_with_err(stdout: str, key: str):
    """Extract mean and error from simulator stdout like 'm_susc: 1.23e-02 4.56e-03'."""
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
            corr = float(parts[3])
            err = float(parts[4])
            conn = float(parts[5])
            conn_err = float(parts[6])
            data[(m, n)] = {
                "d": d, "corr": corr, "err": err,
                "conn": conn, "conn_err": conn_err,
            }
    return data


# ---------------------------------------------------------------------------
# Simulator runner
# ---------------------------------------------------------------------------

def run_simulator(exe: str, Lx: int, Ly: int, Tx: int, Ty: int,
                  k1: float, k2: float, k3: float, beta: float,
                  n_traj: int = 20000, n_skip: int = 20, n_therm: int = 2000,
                  seed: int = 0, data_dir: str = "/tmp/k_from_tp",
                  wall_time: float = 0.0):
    """Run the simulator and return (stdout, subdir_path)."""
    if seed == 0:
        seed = int(time.time() * 1000) & 0xFFFFFFFF

    # Ensure data_dir exists before spawning the process (the C++ MakeDir
    # only handles one level; Python makedirs handles arbitrary depth).
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
        import platform
        print(f"\nERROR: simulator not found: {exe}", file=sys.stderr)
        if platform.system() == "Windows":
            print("  The bundled binary is Linux-only and cannot run on Windows.", file=sys.stderr)
            print("  Please compile the simulator on Windows using MSYS2/MinGW:", file=sys.stderr)
            print("    1. Install MSYS2 from https://www.msys2.org/", file=sys.stderr)
            print("    2. In the MSYS2 terminal run:", file=sys.stderr)
            print("       pacman -S mingw-w64-x86_64-gcc make", file=sys.stderr)
            print("       cd /path/to/K_from_TwoPoint_standalone", file=sys.stderr)
            print("       make", file=sys.stderr)
            print("  Alternatively use WSL (Windows Subsystem for Linux).", file=sys.stderr)
        raise
    if result.returncode != 0:
        print(f"ERROR: simulator failed:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError("simulator failed")

    # Parse the subdir from stdout ("opening file: <path>")
    subdir = None
    for line in result.stdout.splitlines():
        if line.startswith("opening file:"):
            full_path = line.split(":", 1)[1].strip()
            subdir = os.path.dirname(full_path)
            break

    return result.stdout, subdir


# ---------------------------------------------------------------------------
# Inner loop: find beta_c via susceptibility peak
# ---------------------------------------------------------------------------

def find_beta_c(exe: str, Lx: int, Ly: int, Tx: int, Ty: int,
                k1: float, k2: float, k3: float,
                beta_lo: float, beta_hi: float,
                n_coarse: int = 7, n_refine: int = 4, n_refine2: int = 4,
                n_traj_coarse: int = 10000, n_traj_fine: int = 30000,
                data_dir: str = "/tmp/k_from_tp_scan",
                max_shifts: int = 4):
    """
    Locate the susceptibility peak in [beta_lo, beta_hi].

    After each coarse sweep the sorted scan is inspected.  If the maximum
    is at either boundary the window is shifted one full span in that
    direction and new points are added.  This repeats up to max_shifts
    times so the peak is always interior before the golden-section refine.

    Two GC refinement passes are performed:
      Pass 1 (n_refine points, n_traj_fine stats): bracket from coarse GC fit.
      Pass 2 (n_refine2 points, n_traj_fine stats): tighter bracket from pass-1 mode.

    Returns (beta_c, chi_peak, scan_betas, scan_chis, scan_chi_errs,
             (fit_params_1, fit_params_2))
    where each fit_params is (mu, sigma, skew, kurt, A) or None on failure.
    """
    scan_betas: list = []
    scan_chis:  list = []
    scan_chi_errs: list = []
    # Cache already-run betas so we never duplicate a run
    ran: dict = {}  # beta -> (chi, chi_err)

    def _run_beta(b):
        key = round(b, 10)
        if key in ran:
            return ran[key]
        stdout, _ = run_simulator(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, b,
            n_traj=n_traj_coarse, data_dir=data_dir)
        chi, chi_err = parse_stdout_value_with_err(stdout, "m_susc:")
        scan_betas.append(b)
        scan_chis.append(chi)
        scan_chi_errs.append(chi_err)
        ran[key] = (chi, chi_err)
        return chi, chi_err

    span = beta_hi - beta_lo

    # --- Coarse sweeps with automatic window shifting ---
    for shift_iter in range(max_shifts + 1):
        for b in np.linspace(beta_lo, beta_hi, n_coarse):
            _run_beta(b)

        # Sort all accumulated points and find peak
        order = np.argsort(scan_betas)
        sbetas = np.array(scan_betas)[order]
        schis  = np.array(scan_chis)[order]
        idx_max = int(np.argmax(schis))

        at_lo = (idx_max == 0)
        at_hi = (idx_max == len(sbetas) - 1)

        if not at_lo and not at_hi:
            break  # peak is interior — done

        if shift_iter < max_shifts:
            if at_lo:
                new_hi = sbetas[1]               # keep one overlapping point
                new_lo = max(0.01, beta_lo - span)
                print(f"  [scan] peak at lower edge, shifting window to "
                      f"[{new_lo:.4f}, {new_hi:.4f}]")
                beta_lo, beta_hi = new_lo, new_hi
            else:
                new_lo = sbetas[-2]              # keep one overlapping point
                new_hi = beta_hi + span
                print(f"  [scan] peak at upper edge, shifting window to "
                      f"[{new_lo:.4f}, {new_hi:.4f}]")
                beta_lo, beta_hi = new_lo, new_hi
            sys.stdout.flush()

    # Re-sort after all shifts
    order = np.argsort(scan_betas)
    sbetas = np.array(scan_betas)[order]
    idx_max = int(np.argmax(np.array(scan_chis)[order]))

    # --- Refine: Gram-Charlier type A fit (non-zero cumulants κ₁–κ₄, zero higher)
    #     χ(β) = A·exp(-z²/2)·[1 + (s/6)H₃(z) + (κ/24)H₄(z)]
    #     where z = (β-μ)/σ, H₃=z³-3z, H₄=z⁴-6z²+3  (tails → 0)
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    order2 = np.argsort(scan_chis)[::-1]
    top_n = min(7, len(scan_chis))
    betas_top = np.array(scan_betas)[order2[:top_n]]
    chis_top  = np.array(scan_chis)[order2[:top_n]]

    b_range = max(scan_betas) - min(scan_betas)
    beta_parabola = float(sbetas[idx_max])  # fallback
    try:
        p0 = [beta_parabola, b_range / 4, 0.0, 0.0,
              float(max(chis_top))]
        bounds = ([min(scan_betas), 1e-6, -5.0, -3.0, 0.0],
                  [max(scan_betas), b_range,  5.0, 10.0, np.inf])
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(_gram_charlier, betas_top, chis_top,
                                p0=p0, bounds=bounds, maxfev=8000)
        beta_parabola = float(np.clip(popt[0], min(scan_betas), max(scan_betas)))
    except Exception:
        pass

    # Run n_refine fine-stats points symmetrically around the parabola estimate
    half = (max(scan_betas) - min(scan_betas)) / (2 * n_coarse)
    fine_betas = np.linspace(beta_parabola - half, beta_parabola + half, n_refine + 2)
    fine_chis  = []
    for b in fine_betas:
        stdout_f, _ = run_simulator(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, float(b),
            n_traj=n_traj_fine, data_dir=data_dir)
        chi_f, chi_f_err = parse_stdout_value_with_err(stdout_f, "m_susc:")
        scan_betas.append(float(b))
        scan_chis.append(chi_f)
        scan_chi_errs.append(chi_f_err)
        fine_chis.append(chi_f)

    def _gc_fit(betas_arr, chis_arr, beta_init):
        """Fit Gram-Charlier to (betas_arr, chis_arr); return (popt, mode) or (None, beta_init)."""
        n = len(betas_arr)
        b_range = max(betas_arr) - min(betas_arr)
        if b_range < 1e-12:
            return None, beta_init
        use_gc = (n >= 6)
        try:
            import warnings
            if use_gc:
                p0f = [beta_init, b_range / 2, 0.0, 0.0, float(max(chis_arr))]
                bounds_f = ([min(betas_arr), 1e-8, -5.0, -3.0, 0.0],
                            [max(betas_arr), b_range * 2, 5.0, 10.0, np.inf])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    popt, _ = curve_fit(_gram_charlier, betas_arr, chis_arr,
                                        p0=p0f, bounds=bounds_f, maxfev=8000)
            else:
                def _gauss3(b, mu, sigma, A):
                    return A * np.exp(-0.5 * ((b - mu) / sigma) ** 2)
                p0g = [beta_init, b_range / 2, float(max(chis_arr))]
                bounds_g = ([min(betas_arr), 1e-8, 0.0],
                            [max(betas_arr), b_range * 2, np.inf])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pg, _ = curve_fit(_gauss3, betas_arr, chis_arr,
                                      p0=p0g, bounds=bounds_g, maxfev=8000)
                popt = (pg[0], pg[1], 0.0, 0.0, pg[2])
            res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                                  bounds=(min(betas_arr), max(betas_arr)),
                                  method='bounded')
            mode = float(np.clip(res.x, min(betas_arr), max(betas_arr)))
            return tuple(popt), mode
        except Exception:
            return None, beta_init

    # --- Pass 1 GC fit over the first fine points ---
    fit_params_1, beta_pass1 = _gc_fit(list(fine_betas), fine_chis, beta_parabola)
    if fit_params_1 is None:
        beta_pass1 = fine_betas[int(np.argmax(fine_chis))]

    # --- Pass 2: tighter bracket centred on pass-1 mode, half the width ---
    half2 = half / 2.0
    fine2_betas = np.linspace(beta_pass1 - half2, beta_pass1 + half2, n_refine2 + 2)
    fine2_chis  = []
    for b in fine2_betas:
        stdout_f, _ = run_simulator(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, float(b),
            n_traj=n_traj_fine, data_dir=data_dir)
        chi_f, chi_f_err = parse_stdout_value_with_err(stdout_f, "m_susc:")
        scan_betas.append(float(b))
        scan_chis.append(chi_f)
        scan_chi_errs.append(chi_f_err)
        fine2_chis.append(chi_f)

    fit_params_2, beta_c = _gc_fit(list(fine2_betas), fine2_chis, beta_pass1)
    if fit_params_2 is None:
        beta_c = fine2_betas[int(np.argmax(fine2_chis))]

    fit_params = (fit_params_1, fit_params_2)
    parabola_coeffs = fit_params  # keep variable name for compatibility

    chi_peak = max(scan_chis)

    return beta_c, chi_peak, scan_betas, scan_chis, scan_chi_errs, parabola_coeffs


# ---------------------------------------------------------------------------
# Chi-squared between two all-to-all datasets
# ---------------------------------------------------------------------------

def _min_image_r(m: int, n: int, Lx: int, Ly: int) -> float:
    """Physical equilateral distance using minimum-image convention."""
    m_mi = m - Lx * round(m / Lx)
    n_mi = n - Ly * round(n / Ly)
    return math.sqrt(m_mi * m_mi + m_mi * n_mi + n_mi * n_mi)


def _to_grid_conn(data: dict, Lx: int, Ly: int) -> np.ndarray:
    """Fill a (Lx, Ly) array with connected correlator values, zeroing (0,0)."""
    g = np.zeros((Lx, Ly))
    for (m, n), v in data.items():
        g[m % Lx, n % Ly] = v["conn"]
    g[0, 0] = 0.0  # drop singular self-correlator
    return g


def _to_grid_err(data: dict, Lx: int, Ly: int) -> np.ndarray:
    """Fill a (Lx, Ly) array with conn_err values, zeroing (0,0)."""
    g = np.zeros((Lx, Ly))
    for (m, n), v in data.items():
        g[m % Lx, n % Ly] = v["conn_err"]
    g[0, 0] = 0.0
    return g


def _fem_bilinear_integral(h: np.ndarray) -> float:
    """
    Exact integral of (bilinear manifold of h)^2 over the full periodic domain.

    For each (Lx x Ly) cell [i,i+1] x [j,j+1] with corner values a,b,c,d at
    (i,j),(i+1,j),(i,j+1),(i+1,j+1) (periodic), the exact unit-cell integral is:

        I = a^2 + u^2/3 + v^2/3 + w^2/9
            + a*u + a*v + a*w/2 + u*v/2 + u*w/3 + v*w/3

    where u = b-a,  v = c-a,  w = d-b-c+a.
    Summing over all Lx*Ly cells gives the total integral (each cell has area 1).
    """
    Lx, Ly = h.shape
    ip = (np.arange(Lx) + 1) % Lx   # i+1 mod Lx
    jp = (np.arange(Ly) + 1) % Ly   # j+1 mod Ly
    a = h                            # shape (Lx, Ly)
    b = h[ip, :]                     # h[i+1, j]
    c = h[:, jp]                     # h[i,   j+1]
    d = h[np.ix_(ip, jp)]            # h[i+1, j+1]
    u = b - a
    v = c - a
    w = d - b - c + a
    cell_integrals = (a*a + u*u/3.0 + v*v/3.0 + w*w/9.0
                     + a*u + a*v + a*w/2.0 + u*v/2.0 + u*w/3.0 + v*w/3.0)
    return float(cell_integrals.sum())


def compute_chi2(ref_data: dict, test_data: dict, r_min: float = 0.0,
                 r_max_frac: float = 0.33, L_eff: float = 12.0,
                 Lx: int = 0, Ly: int = 0,
                 cost: str = "log_ratio",
                 test_data_shift: dict = None, delta_beta: float = 0.0,
                 beta_diff: float = 0.0):
    """
    Compute chi2 between reference and test all-to-all correlators.

    cost options
    ------------
    'log_ratio'  (default)
        Amplitude-free via GLS weighted-mean subtraction.
        chi2 = sum (ld_i - mu_w)^2 / sigma_i^2,  ndof = N-1.

    'beta_deriv'  (option 4 — first-order beta-mismatch correction)
        Requires test_data_shift (run at beta_c_test + delta_beta),
        beta_diff = beta_c_test - beta_c_ref.
        Corrects each log-diff: ld_corr = ld - beta_diff * dlogG/dbeta
        where dlogG/dbeta ≈ (log G_shift - log G_test) / delta_beta.
        Then applies GLS mean subtraction and L2 chi2.

    'pair_ratio'  (option 2)
        All i<j pairs of log-diffs: amplitude cancels exactly.
        chi2 = sum_{i<j} (ld_i - ld_j)^2 / (sigma_i^2 + sigma_j^2),
        ndof = N*(N-1)/2.

    'residuals'   (no amplitude correction)
        chi2 = sum ld_i^2 / sigma_i^2,  ndof = N.

    'fem_integral'
        Both G_conn fields are placed on an Lx×Ly grid and treated as
        piecewise bilinear manifolds (periodic BC).  The cost is the exact
        analytical integral of (G_test - G_ref)^2 over the full domain,
        computed per-cell from the 4 corner values.  Physical amplitudes
        are preserved.  Returned as the raw integral (ndof=1), so the
        displayed chi2/ndof equals the total integrated difference.
    """
    # ---- FEM integral cost (completely different code path) ----
    if cost == "fem_integral":
        if Lx <= 0 or Ly <= 0:
            raise ValueError("fem_integral cost requires Lx and Ly > 0")
        g_ref  = _to_grid_conn(ref_data,  Lx, Ly)
        g_test = _to_grid_conn(test_data, Lx, Ly)
        e_ref  = _to_grid_err(ref_data,   Lx, Ly)
        e_test = _to_grid_err(test_data,  Lx, Ly)
        h = g_test - g_ref
        sigma_h = np.sqrt(e_ref**2 + e_test**2)  # pointwise uncertainty on h
        total     = _fem_bilinear_integral(h)
        total_err = 2.0 * np.sqrt(_fem_bilinear_integral(h * sigma_h))
        return total, 1, total_err
    r_max = r_max_frac * L_eff

    # ---- single pass: collect all per-displacement quantities ----
    log_diffs    = []   # log(G_test) - log(G_ref)
    rel_errs     = []   # combined relative errors
    d_log_g_dbeta = []  # (log G_shift - log G_test) / delta_beta  for beta_deriv

    for (m, n), ref in ref_data.items():
        if (m, n) not in test_data:
            continue
        if m == 0 and n == 0:
            continue
        r = (_min_image_r(m, n, Lx, Ly) if Lx > 0 and Ly > 0
             else math.sqrt(m*m + m*n + n*n))
        if r < r_min or r > r_max:
            continue

        test = test_data[(m, n)]
        g_ref  = ref["conn"]
        g_test = test["conn"]
        if g_ref <= 0 or g_test <= 0:
            continue

        rel_err = math.sqrt((ref["conn_err"] / g_ref)**2 +
                            (test["conn_err"] / g_test)**2)
        if rel_err < 1e-15:
            continue

        log_diffs.append(math.log(g_test) - math.log(g_ref))
        rel_errs.append(rel_err)

        # derivative estimate for beta_deriv
        if (cost == "beta_deriv" and test_data_shift and
                (m, n) in test_data_shift and delta_beta != 0.0):
            g_sh = test_data_shift[(m, n)]["conn"]
            if g_sh > 0:
                d_log_g_dbeta.append((math.log(g_sh) - math.log(g_test)) / delta_beta)
            else:
                d_log_g_dbeta.append(0.0)
        else:
            d_log_g_dbeta.append(None)

    ndof = len(log_diffs)
    if ndof == 0:
        return 0.0, 0, 0.0

    weights = [1.0 / re**2 for re in rel_errs]

    # ---- cost-function branches ----

    if cost == "residuals":
        chi2 = sum(w * ld**2 for w, ld in zip(weights, log_diffs))
        return chi2, ndof, 0.0

    if cost == "pair_ratio":
        chi2 = 0.0
        n_pairs = 0
        for i in range(ndof):
            for j in range(i + 1, ndof):
                pair_ld  = log_diffs[i] - log_diffs[j]
                pair_var = rel_errs[i]**2 + rel_errs[j]**2
                chi2 += pair_ld**2 / pair_var
                n_pairs += 1
        return chi2, max(n_pairs, 1), 0.0

    if cost == "beta_deriv" and delta_beta != 0.0:
        # Apply correction and fall through to GLS-mean chi2
        log_diffs = [
            ld - beta_diff * (dg if dg is not None else 0.0)
            for ld, dg in zip(log_diffs, d_log_g_dbeta)
        ]

    # log_ratio (default) and corrected beta_deriv:
    # GLS weighted mean, then L2
    sum_w = sum(weights)
    mu_w  = sum(w * ld for w, ld in zip(weights, log_diffs)) / sum_w
    chi2  = sum(w * (ld - mu_w)**2 for w, ld in zip(weights, log_diffs))
    return chi2, max(ndof - 1, 1), 0.0


# ---------------------------------------------------------------------------
# Live-updating plot
# ---------------------------------------------------------------------------

class LivePlotter:
    def __init__(self, output_dir: str, cost: str = "log_ratio"):
        self.output_dir = output_dir
        self.cost = cost
        # All evaluated points across all levels: (r1, r2, chi2_ndof)
        self.all_points = []
        # Current grid axes (set by start_level)
        self.current_grid_r1 = None
        self.current_grid_r2 = None
        self.current_level = 0
        # Latest susceptibility scan
        self.last_scan_betas = []
        self.last_scan_chis = []
        self.last_scan_chi_errs = []
        self.last_beta_c = None
        self.last_parabola = None  # (fit_params_1, fit_params_2) tuple
        # Current in-progress grid point (r1, r2) or None
        self.current_point = None
        # Per-evaluation history for panels 2 & 3
        self.history_chi2_ndof = []
        self.history_chi2_err = []
        self.history_beta_c = []
        # Build figure with a fixed-width colorbar column next to ax1 so
        # adding/removing colorbars never resizes the scatter panel.
        from matplotlib.gridspec import GridSpec
        self.fig = plt.figure(figsize=(14, 10))
        gs = GridSpec(2, 3, figure=self.fig,
                      width_ratios=[1, 0.04, 1],
                      wspace=0.45, hspace=0.38)
        self.ax1     = self.fig.add_subplot(gs[0, 0])   # heatmap
        self.cbar_ax = self.fig.add_subplot(gs[0, 1])   # dedicated colorbar
        self.ax2     = self.fig.add_subplot(gs[0, 2])   # chi2 history
        self.ax3     = self.fig.add_subplot(gs[1, 0])   # beta_c history
        self.ax4     = self.fig.add_subplot(gs[1, 2])   # susceptibility scan
        # keep .axes for flat iteration in _update
        self.axes = [self.ax1, self.ax2, self.ax3, self.ax4]
        self.fig.suptitle("K_from_TwoPoint: Adaptive Grid Search", fontsize=14)

    def start_level(self, level, r1_vals, r2_vals):
        self.current_level = level
        self.current_grid_r1 = np.array(r1_vals)
        self.current_grid_r2 = np.array(r2_vals)

    def set_current_point(self, r1, r2):
        """Call before evaluate_point to highlight the in-progress grid point."""
        self.current_point = (r1, r2)

    def record_point(self, r1, r2, beta_c, chi2, ndof,
                     scan_betas, scan_chis, scan_chi_errs=None, parabola=None,
                     chi2_err=0.0):
        chi2_ndof = chi2 / max(ndof, 1)
        self.all_points.append((r1, r2, chi2_ndof))
        self.last_scan_betas = scan_betas
        self.last_scan_chis = scan_chis
        self.last_scan_chi_errs = scan_chi_errs or []
        self.last_beta_c = beta_c
        self.last_parabola = parabola  # (fit_params_1, fit_params_2)
        self.history_chi2_ndof.append(chi2_ndof)
        self.history_chi2_err.append(chi2_err)
        self.history_beta_c.append(beta_c)
        self.current_point = None  # clear after recording
        self._update()

    def _update(self):
        for ax in self.axes:
            ax.clear()
        self.cbar_ax.clear()

        n = len(self.history_chi2_ndof)
        iters = list(range(1, n + 1))

        # Panel 1: heatmap coloured by chi2/ndof, fixed cbar_ax so ax1 never resizes
        ax1 = self.ax1
        if self.all_points:
            r1s = np.array([p[0] for p in self.all_points])
            r2s = np.array([p[1] for p in self.all_points])
            vals = np.array([p[2] for p in self.all_points])
            vmin, vmax = vals.min(), max(vals.max(), vals.min() + 1e-12)
            sc = ax1.scatter(r1s, r2s, c=vals, cmap="RdYlGn_r",
                             vmin=vmin, vmax=vmax, s=80, zorder=3,
                             edgecolors="k", linewidths=0.5)
            _score_label = ("Integrated difference" if self.cost == "fem_integral"
                            else r"$\chi^2/N_{\rm dof}$")
            self.fig.colorbar(sc, cax=self.cbar_ax, label=_score_label)
            best_idx = int(np.argmin(vals))
            ax1.plot(r1s[best_idx], r2s[best_idx], "*", color="gold",
                     markersize=16, zorder=5, label="best")
            if self.current_point is not None:
                ax1.plot(self.current_point[0], self.current_point[1],
                         "P", color="cyan", markersize=14, zorder=6,
                         label="in progress", markeredgecolor="k", markeredgewidth=0.8)
            ax1.legend(fontsize=9)
        if self.current_grid_r1 is not None and len(self.current_grid_r1) > 1:
            r1lo, r1hi = self.current_grid_r1[0], self.current_grid_r1[-1]
            r2lo, r2hi = self.current_grid_r2[0], self.current_grid_r2[-1]
            from matplotlib.patches import Rectangle
            rect = Rectangle((r1lo, r2lo), r1hi - r1lo, r2hi - r2lo,
                              fill=False, edgecolor="royalblue", lw=1.8,
                              linestyle="--", zorder=4)
            ax1.add_patch(rect)
        ax1.set_xlabel("r1 = k1/k3")
        ax1.set_ylabel("r2 = k2/k3")
        ax1.set_title(f"Grid heatmap (level {self.current_level})")
        ax1.grid(True, alpha=0.2)

        # Panel 2: chi2/ndof vs evaluation number
        ax2 = self.ax2
        if n:
            errs = np.array(self.history_chi2_err)
            vals = np.array(self.history_chi2_ndof)
            has_errs = errs.any()
            if has_errs:
                ax2.errorbar(iters, vals, yerr=errs, fmt="o", color="crimson",
                             ms=4, capsize=3, lw=1.0, zorder=3)
            else:
                ax2.semilogy(iters, vals, "o", color="crimson",
                             ms=4, zorder=3)
        if self.cost != "fem_integral":
            ax2.axhline(1.0, ls="--", color="gray", alpha=0.5,
                        label=r"$\chi^2/N_{\rm dof}=1$")
        ax2.set_xlabel("Evaluation #")
        _y_label = ("Integrated difference" if self.cost == "fem_integral"
                    else r"$\chi^2 / N_{\rm dof}$")
        ax2.set_ylabel(_y_label)
        ax2.set_title("Fit quality")
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        # Panel 3: beta_c vs evaluation number
        ax3 = self.ax3
        if n:
            ax3.plot(iters, self.history_beta_c, "o", color="darkorange", ms=4)
        ax3.set_xlabel("Evaluation #")
        ax3.set_ylabel(r"$\beta_c$")
        ax3.set_title(r"Critical $\beta$ (susceptibility peak)")
        ax3.grid(True, alpha=0.3)

        # Panel 4: latest susceptibility scan with error bars
        ax4 = self.ax4
        if self.last_scan_betas:
            order = np.argsort(self.last_scan_betas)
            sb = np.array(self.last_scan_betas)[order]
            sc_arr = np.array(self.last_scan_chis)[order]
            if (self.last_scan_chi_errs and
                    len(self.last_scan_chi_errs) == len(self.last_scan_betas)):
                se = np.array(self.last_scan_chi_errs)[order]
                ax4.errorbar(sb, sc_arr, yerr=se, fmt="o-", color="teal",
                             ms=4, capsize=3, lw=1.2)
            else:
                ax4.plot(sb, sc_arr, "o-", color="teal", ms=4)
            if self.last_parabola is not None:
                # last_parabola is (fit_params_1, fit_params_2)
                try:
                    b_lo_p = min(self.last_scan_betas)
                    b_hi_p = max(self.last_scan_betas)
                    b_fit = np.linspace(b_lo_p, b_hi_p, 300)

                    def _gc_curve(b_arr, fp):
                        mu, sigma, skew, kurt, A = fp
                        z = (b_arr - mu) / sigma
                        H3 = z**3 - 3.0 * z
                        H4 = z**4 - 6.0 * z**2 + 3.0
                        return A * np.exp(-0.5*z**2) * (1.0 + (skew/6.0)*H3 + (kurt/24.0)*H4)

                    fp1, fp2 = self.last_parabola
                    if fp1 is not None:
                        ax4.plot(b_fit, _gc_curve(b_fit, fp1), "--",
                                 color="orange", lw=1.6, zorder=2,
                                 label="GC fit 1 (coarse)")
                    if fp2 is not None:
                        ax4.plot(b_fit, _gc_curve(b_fit, fp2), "-",
                                 color="red", lw=2.0, zorder=3,
                                 label="GC fit 2 (fine)")
                except Exception:
                    pass
            if self.last_beta_c is not None:
                ax4.axvline(self.last_beta_c, ls="--", color="red", alpha=0.7,
                            label=rf"$\beta_c$ = {self.last_beta_c:.6f}")
            ax4.legend(fontsize=9)
        ax4.set_xlabel(r"$\beta$")
        ax4.set_ylabel(r"$\chi$ (susceptibility)")
        ax4.set_title("Latest susceptibility scan")
        ax4.grid(True, alpha=0.3)

        self.fig.suptitle("K_from_TwoPoint: Adaptive Grid Search",
                           fontsize=14, y=0.995)
        path = os.path.join(self.output_dir, "optimisation_live.png")
        self.fig.savefig(path, dpi=120)
        print(f"  [plot] saved {path}")
        # Display inline in Spyder / IPython console
        try:
            from IPython.display import clear_output, display
            clear_output(wait=True)
            display(self.fig)
        except Exception:
            try:
                self.fig.canvas.draw()
                self.fig.canvas.flush_events()
            except Exception:
                pass

    def save_final(self):
        path = os.path.join(self.output_dir, "optimisation_final.png")
        self.fig.savefig(path, dpi=150)
        print(f"  [plot] saved final: {path}")
        plt.close(self.fig)


# ---------------------------------------------------------------------------
# Point evaluation
# ---------------------------------------------------------------------------

def evaluate_point(exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_guess,
                   n_traj_prod, data_dir, ref_data, L_eff, label, plotter,
                   n_traj_scan_coarse=30000, n_traj_scan_fine=60000,
                   r_min=0.0, r_max_frac=0.33,
                   cost="log_ratio", beta_ref=None):
    """Find beta_c, run production, compute chi2 for a single (k1, k2, k3) point."""
    r1 = k1 / k3
    r2 = k2 / k3

    margin = max(0.15 * beta_guess, 0.02)
    beta_lo = max(0.01, beta_guess - margin)
    beta_hi = beta_guess + margin

    print(f"  [{label}] r1={r1:.6f} r2={r2:.6f}: scanning beta [{beta_lo:.4f}, {beta_hi:.4f}]")
    sys.stdout.flush()

    scan_dir = os.path.join(data_dir, f"scan_{label}")
    beta_c, chi_peak, scan_betas, scan_chis, scan_chi_errs, parabola = find_beta_c(
        exe, Lx, Ly, Tx, Ty, k1, k2, k3,
        beta_lo, beta_hi,
        n_coarse=7, n_refine=3, n_refine2=3,
        n_traj_coarse=n_traj_scan_coarse, n_traj_fine=n_traj_scan_fine,
        data_dir=scan_dir,
    )
    print(f"  [{label}] beta_c = {beta_c:.8f} (chi_peak = {chi_peak:.4e})")
    sys.stdout.flush()

    prod_dir = os.path.join(data_dir, f"prod_{label}")
    stdout, subdir = run_simulator(
        exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_c,
        n_traj=n_traj_prod, n_therm=3000,
        data_dir=prod_dir,
    )
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    test_data = load_all_to_all(a2a_path)

    # beta_deriv: extra short run at beta_c + delta to estimate d log G / d beta
    test_data_shift = None
    delta_beta = 0.0
    beta_diff  = 0.0
    if cost == "beta_deriv" and beta_ref is not None:
        delta_beta = max(abs(beta_c - beta_ref) * 0.5, 5e-4)
        beta_shift = beta_c + delta_beta
        n_shift = max(n_traj_prod // 2, 500)
        shift_dir = os.path.join(data_dir, f"shift_{label}")
        print(f"  [{label}] beta_deriv: extra run at beta={beta_shift:.6f} "
              f"(delta={delta_beta:.4e}, n={n_shift})")
        sys.stdout.flush()
        _, shift_subdir = run_simulator(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_shift,
            n_traj=n_shift, n_therm=1000, data_dir=shift_dir,
        )
        shift_path = os.path.join(shift_subdir, "two_point_all_to_all.dat")
        test_data_shift = load_all_to_all(shift_path)
        beta_diff = beta_c - beta_ref

    chi2, ndof, chi2_err = compute_chi2(
        ref_data, test_data,
        r_min=r_min, r_max_frac=r_max_frac, L_eff=L_eff, Lx=Lx, Ly=Ly,
        cost=cost,
        test_data_shift=test_data_shift,
        delta_beta=delta_beta,
        beta_diff=beta_diff,
    )
    chi2_ndof = chi2 / max(ndof, 1)
    err_str = f" ± {chi2_err:.4f}" if chi2_err > 0 else ""
    print(f"  [{label}] chi2/ndof = {chi2_ndof:.4f}{err_str} (ndof={ndof}, cost={cost})")
    sys.stdout.flush()

    plotter.record_point(r1, r2, beta_c, chi2, ndof,
                         scan_betas, scan_chis, scan_chi_errs, parabola=parabola,
                         chi2_err=chi2_err)

    return beta_c, chi2, ndof, chi2_ndof


# ---------------------------------------------------------------------------
# Adaptive 4x4 grid search
# ---------------------------------------------------------------------------

def grid_search(exe, Lx, Ly, Tx, Ty, ref_data, L_eff,
                r1_center, r2_center, half_span,
                beta_guess, n_traj_prod, data_dir, output_dir,
                n_grid=4, max_levels=10, plotter=None,
                n_traj_scan_coarse=30000, n_traj_scan_fine=60000,
                r_min=1.5, r_max_frac=0.33,
                cost="log_ratio", beta_ref=None):
    """
    Adaptive grid search over (r1, r2) = (k1/k3, k2/k3).

    Rules per level:
      - Evaluate all n_grid x n_grid points.
      - Interior minimum: refine (halve half_span) centred on the best point.
      - Border minimum:   translate the grid so the best point becomes the
                         new centre (half_span unchanged).
    Stops when chi2/ndof <= 1.0 or max_levels reached.
    """
    all_results = []   # (r1, r2, beta_c, chi2, ndof, chi2_ndof)
    best_beta_c = beta_guess
    grid_results = {}

    for level in range(max_levels):
        r1_lo = max(r1_center - half_span, 0.1)
        r1_hi = r1_center + half_span
        r2_lo = max(r2_center - half_span, 0.1)
        r2_hi = r2_center + half_span

        r1_vals = np.linspace(r1_lo, r1_hi, n_grid)
        r2_vals = np.linspace(r2_lo, r2_hi, n_grid)

        print(f"\n{'='*60}")
        print(f"Grid level {level}: centre=({r1_center:.6f}, {r2_center:.6f}),"
              f" half_span={half_span:.6f}")
        print(f"  r1 in [{r1_lo:.6f}, {r1_hi:.6f}]")
        print(f"  r2 in [{r2_lo:.6f}, {r2_hi:.6f}]")
        print(f"{'='*60}")
        sys.stdout.flush()

        if plotter:
            plotter.start_level(level, r1_vals, r2_vals)

        grid_results = {}
        for i, r1 in enumerate(r1_vals):
            for j, r2 in enumerate(r2_vals):
                k3 = 1.0
                k1, k2 = r1 * k3, r2 * k3
                label = f"L{level}_i{i}_j{j}"
                level_dir = os.path.join(data_dir, f"level{level:02d}")
                os.makedirs(level_dir, exist_ok=True)
                if plotter:
                    plotter.set_current_point(r1, r2)
                    plotter._update()
                beta_c, chi2, ndof, chi2_ndof = evaluate_point(
                    exe, Lx, Ly, Tx, Ty, k1, k2, k3, best_beta_c,
                    n_traj_prod, level_dir, ref_data, L_eff, label, plotter,
                    n_traj_scan_coarse=n_traj_scan_coarse,
                    n_traj_scan_fine=n_traj_scan_fine,
                    r_min=r_min, r_max_frac=r_max_frac,
                    cost=cost, beta_ref=beta_ref,
                )
                grid_results[(i, j)] = (r1, r2, beta_c, chi2, chi2_ndof)
                all_results.append((r1, r2, beta_c, chi2, ndof, chi2_ndof))

        # Best grid point
        best_ij = min(grid_results, key=lambda k: grid_results[k][4])
        bi, bj = best_ij
        r1_best, r2_best, beta_c_best, chi2_best, chi2_ndof_best = grid_results[best_ij]

        print(f"\n  Level {level} best: ({bi},{bj})  r1={r1_best:.8f}  "
              f"r2={r2_best:.8f}  chi2/ndof={chi2_ndof_best:.4f}")
        sys.stdout.flush()
        best_beta_c = beta_c_best

        # Save per-level summary
        level_file = os.path.join(output_dir, f"grid_level_{level:02d}.json")
        with open(level_file, "w") as f:
            json.dump({
                "level": level,
                "r1_center": r1_center, "r2_center": r2_center,
                "half_span": half_span,
                "r1_vals": r1_vals.tolist(), "r2_vals": r2_vals.tolist(),
                "best_ij": [bi, bj],
                "best_r1": r1_best, "best_r2": r2_best,
                "best_chi2_ndof": chi2_ndof_best,
                "best_beta_c": beta_c_best,
            }, f, indent=2)

        is_border = (bi == 0 or bi == n_grid - 1 or
                     bj == 0 or bj == n_grid - 1)

        if not is_border and chi2_ndof_best <= 1.0:
            print(f"  Converged at level {level}: chi2/ndof = {chi2_ndof_best:.4f}")
            break

        if is_border:
            r1_center = r1_best
            r2_center = r2_best
            print(f"  Border minimum → translating centre to "
                  f"({r1_center:.6f}, {r2_center:.6f})")
        else:
            r1_center = r1_best
            r2_center = r2_best
            half_span /= 2.0
            print(f"  Interior minimum → refining, new half_span={half_span:.6f}")
        sys.stdout.flush()

    return all_results, grid_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Find critical couplings by two-point matching")
    parser.add_argument("--exe", default="bin/ising_tri_twisted_parallelogram")
    parser.add_argument("--Lx", type=int, default=12)
    parser.add_argument("--Ly", type=int, default=12)
    parser.add_argument("--Tx", type=int, default=0)
    parser.add_argument("--Ty", type=int, default=0)
    parser.add_argument("--ref", default=None,
                        help="Path to reference two_point_all_to_all.dat "
                             "(or ref_metadata.json from --gen_ref)")
    parser.add_argument("--r1_init", type=float, default=1.0,
                        help="Initial grid centre r1 = k1/k3")
    parser.add_argument("--r2_init", type=float, default=1.0,
                        help="Initial grid centre r2 = k2/k3")
    parser.add_argument("--half_span", type=float, default=0.2,
                        help="Initial half-width of the grid in (r1,r2) space")
    parser.add_argument("--beta_init", type=float, default=0.28,
                        help="Initial guess for beta_c")
    parser.add_argument("--n_traj_prod", type=int, default=50000,
                        help="Trajectories for production runs")
    parser.add_argument("--n_traj_scan_coarse", type=int, default=30000,
                        help="Trajectories per beta point for coarse susceptibility scan")
    parser.add_argument("--n_traj_scan_fine", type=int, default=60000,
                        help="Trajectories per beta point for fine (refine) susceptibility scan")
    parser.add_argument("--r_min", type=float, default=0.0,
                        help="Minimum displacement distance included in chi2 sum")
    parser.add_argument("--r_max_frac", type=float, default=0.33,
                        help="Maximum displacement as fraction of L_eff")
    parser.add_argument("--cost", default="log_ratio",
                        choices=["log_ratio", "beta_deriv", "pair_ratio", "residuals", "fem_integral"],
                        help=("Cost function: "
                              "'log_ratio' (GLS-mean amplitude-free, default); "
                              "'beta_deriv' (first-order beta-mismatch correction); "
                              "'pair_ratio' (all log-ratio pairs, exact amplitude cancellation); "
                              "'residuals' (plain L2, no amplitude correction)"))
    parser.add_argument("--beta_ref", type=float, default=None,
                        help="Reference beta_c (used by beta_deriv cost). "
                             "Loaded from ref_metadata.json if not given.")
    parser.add_argument("--max_iter", type=int, default=10,
                        help="Maximum grid refinement/translation levels")
    parser.add_argument("--n_grid", type=int, default=5,
                        help="Grid points per side (default 5 → 5x5 grid per level)")
    parser.add_argument("--output_dir", default="results/default_run")
    parser.add_argument("--gen_ref", action="store_true",
                        help="Generate reference data at (r1,r2)=(1,1) using "
                             "the same susceptibility-peak finder, then exit.")
    parser.add_argument("--ref_n_traj", type=int, default=500000,
                        help="Trajectories for reference production run")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Self-consistent reference generation mode ---
    if args.gen_ref:
        print("Generating self-consistent reference at k=(1,1,1)...")
        sys.stdout.flush()
        ref_dir = os.path.join(args.output_dir, "ref_scan")
        os.makedirs(ref_dir, exist_ok=True)
        margin = max(0.15 * args.beta_init, 0.03)
        beta_lo = max(0.01, args.beta_init - margin)
        beta_hi = args.beta_init + margin
        print(f"  Susceptibility scan: beta in [{beta_lo:.4f}, {beta_hi:.4f}]")
        sys.stdout.flush()
        beta_c, chi_peak, scan_betas, scan_chis, scan_chi_errs, _ = find_beta_c(
            args.exe, args.Lx, args.Ly, args.Tx, args.Ty,
            1.0, 1.0, 1.0,
            beta_lo, beta_hi,
            n_coarse=11, n_refine=5, n_refine2=5,
            n_traj_coarse=100000, n_traj_fine=200000,
            data_dir=ref_dir,
        )
        print(f"  beta_c = {beta_c:.8f} (chi_peak = {chi_peak:.6e})")
        print(f"  Running production at beta_c with {args.ref_n_traj} trajectories...")
        sys.stdout.flush()
        prod_dir = os.path.join(args.output_dir, "ref_production")
        stdout, subdir = run_simulator(
            args.exe, args.Lx, args.Ly, args.Tx, args.Ty,
            1.0, 1.0, 1.0, beta_c,
            n_traj=args.ref_n_traj, n_therm=5000,
            seed=42, data_dir=prod_dir,
        )
        a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
        print(f"  Reference data written to: {a2a_path}")
        print(f"  beta_c = {beta_c:.10f}")
        meta = {"beta_c": beta_c, "chi_peak": chi_peak,
                "Lx": args.Lx, "Ly": args.Ly, "Tx": args.Tx, "Ty": args.Ty,
                "a2a_path": a2a_path}
        meta_path = os.path.join(args.output_dir, "ref_metadata.json")
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)
        print(f"  Metadata saved to: {meta_path}")
        return

    # Load reference
    if args.ref is None:
        meta_path = os.path.join(args.output_dir, "ref_metadata.json")
        if os.path.exists(meta_path):
            with open(meta_path) as f:
                meta = json.load(f)
            args.ref = meta["a2a_path"]
            if args.beta_ref is None and "beta_c" in meta:
                args.beta_ref = meta["beta_c"]
            print(f"Loaded reference path from {meta_path}")
        else:
            parser.error("--ref is required (or run --gen_ref first)")
    # If --ref points to a .json metadata file, resolve the a2a_path from it
    if args.ref.endswith(".json"):
        with open(args.ref) as f:
            meta = json.load(f)
        meta_dir = os.path.dirname(os.path.abspath(args.ref))
        a2a = meta["a2a_path"]
        # Make a2a_path absolute if it is relative (resolve relative to package root
        # first, then relative to the metadata file's directory as a fallback)
        if not os.path.isabs(a2a):
            candidate = a2a  # relative to cwd (package root)
            if not os.path.isfile(candidate):
                candidate = os.path.join(meta_dir, a2a)
            a2a = candidate
        args.ref = a2a
        if args.beta_ref is None and "beta_c" in meta:
            args.beta_ref = meta["beta_c"]
        print(f"Resolved reference from metadata: {args.ref}")
    # Try to read beta_ref from metadata alongside the .dat file if still None
    if args.beta_ref is None:
        meta_guess = os.path.join(os.path.dirname(os.path.dirname(args.ref)),
                                  "ref_metadata.json")
        if os.path.exists(meta_guess):
            with open(meta_guess) as f:
                _m = json.load(f)
            args.beta_ref = _m.get("beta_c")
    if args.cost == "beta_deriv" and args.beta_ref is None:
        parser.error("--cost beta_deriv requires --beta_ref (or a ref_metadata.json "
                     "with a beta_c field)")
    print(f"Loading reference: {args.ref}")
    sys.stdout.flush()
    ref_data = load_all_to_all(args.ref)
    print(f"  {len(ref_data)} displacement entries loaded")

    L_eff = min(args.Lx, args.Ly)
    plotter = LivePlotter(args.output_dir, cost=args.cost)

    data_dir = os.path.join(args.output_dir, "runs")
    os.makedirs(data_dir, exist_ok=True)

    print(f"\nStarting adaptive 4x4 grid search")
    print(f"  centre    = ({args.r1_init}, {args.r2_init})")
    print(f"  half_span = {args.half_span}")
    print(f"  max levels = {args.max_iter}")
    print(f"  cost function = {args.cost}")
    if args.beta_ref is not None:
        print(f"  beta_ref = {args.beta_ref:.8f}")
    sys.stdout.flush()

    all_results, final_grid = grid_search(
        exe=args.exe,
        Lx=args.Lx, Ly=args.Ly, Tx=args.Tx, Ty=args.Ty,
        ref_data=ref_data, L_eff=L_eff,
        r1_center=args.r1_init,
        r2_center=args.r2_init,
        half_span=args.half_span,
        beta_guess=args.beta_init,
        n_traj_prod=args.n_traj_prod,
        data_dir=data_dir,
        output_dir=args.output_dir,
        n_grid=args.n_grid,
        max_levels=args.max_iter,
        plotter=plotter,
        n_traj_scan_coarse=args.n_traj_scan_coarse,
        n_traj_scan_fine=args.n_traj_scan_fine,
        r_min=args.r_min,
        r_max_frac=args.r_max_frac,
        cost=args.cost,
        beta_ref=args.beta_ref,
    )

    if not all_results:
        print("No results collected.")
        return

    best = min(all_results, key=lambda x: x[5])
    r1b, r2b, bc, chi2b, ndofb, chi2_ndofb = best

    print(f"\n{'='*60}")
    print("OPTIMISATION COMPLETE")
    print(f"{'='*60}")
    print(f"  r1 = {r1b:.8f},  r2 = {r2b:.8f},  k3 = 1.0")
    print(f"  beta_c = {bc:.8f}")
    print(f"  chi2/ndof = {chi2_ndofb:.4f}")
    p1 = math.exp(-2 * bc * r1b)
    p2 = math.exp(-2 * bc * r2b)
    p3 = math.exp(-2 * bc * 1.0)
    print(f"\n  Analytical check: exp(-2β·k1)+exp(-2β·k2)+exp(-2β·k3) = {p1+p2+p3:.6f}"
          f"  (should be 1.0)")

    plotter.save_final()

    result_file = os.path.join(args.output_dir, "fit_result.json")
    with open(result_file, "w") as f:
        json.dump({
            "best_r1": r1b, "best_r2": r2b, "best_k3": 1.0,
            "best_beta_c": bc,
            "best_chi2": chi2b, "best_ndof": ndofb,
            "best_chi2_ndof": chi2_ndofb,
            "n_eval": len(all_results),
            "all_results": [
                {"r1": x[0], "r2": x[1], "beta_c": x[2],
                 "chi2": x[3], "ndof": x[4], "chi2_ndof": x[5]}
                for x in all_results
            ],
        }, f, indent=2)
    print(f"Saved: {result_file}")


if __name__ == "__main__":
    main()
