#!/usr/bin/env python3
"""
run_fss_betac.py  —  Finite-Size Scaling β_c test
==================================================
For several lattice sizes L, each β_c finder method locates the
susceptibility peak N_TRIALS times.  The results are:

  1. β_c(L) vs L for each method (with error bars)
  2. Finite-size scaling fit:  β_c(L) = β_c(∞) + a·L^{-1/ν}  (ν=1 for 2D Ising)
  3. Extrapolated β_c(∞) compared to the known exact value ln(3)/4

This comparison does NOT rely on any "exact" reference — the only ground truth
is the analytically known β_c(∞) = ln(3)/4 ≈ 0.274653 for the equilateral
triangular 2D Ising model.
"""

import os, sys, json, time, math, copy, warnings
import numpy as np

# ===========================================================================
# CONFIG
# ===========================================================================
# Lattice sizes to scan.  Keep Lx=Ly (square equilateral).
# Larger sizes are slower (MC cost scales as L²) — adjust N_TRAJ if needed.
LATTICE_SIZES = [4, 6, 8, 12, 16]

# Trials per (method, L) pair.  More trials → tighter error bars on β_c(L).
N_TRIALS = 10

# MC statistics: scale with volume to keep signal roughly comparable.
# These are the BASE values at L=8; they scale as (L/8)² for larger L.
BASE_N_TRAJ_COARSE = 20000
BASE_N_TRAJ_FINE   = 40000

# Coupling ratios — equilateral point.
k1, k2, k3 = 1.0, 1.0, 1.0

# β_c guess for scan window.  The true value is ln(3)/4 ≈ 0.2747 at L=∞,
# but finite-size peaks are shifted lower for small L.  The window is set
# wide enough (±40%) to capture peaks from L=4 onwards.
BETA_GUESS = 0.24
BETA_MARGIN_FRAC = 0.40     # ±40% window
BETA_MARGIN_MIN  = 0.06     # at least this wide

EXE = "bin/ising_tri_twisted_parallelogram"
OUTPUT_DIR = "results/fss_betac"
# ===========================================================================

EXACT_BETA_C_INF = math.log(3) / 4  # 0.274653072167...

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
os.chdir(_HERE)

import optimise_couplings as oc
from scipy.optimize import curve_fit, minimize_scalar

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Scan window — auto-adjust per L
# ---------------------------------------------------------------------------
def scan_window(L):
    """Return (beta_lo, beta_hi) for a given lattice size.
    Smaller lattices have peaks shifted further from β_c(∞),
    so we use a wider window."""
    margin = max(BETA_MARGIN_FRAC * BETA_GUESS, BETA_MARGIN_MIN)
    # Shift centre toward expected peak: smaller L → lower peak
    # Rough empirical shift: Δβ ≈ -0.04·(8/L)
    centre = BETA_GUESS + 0.01 * (L / 8 - 1)
    centre = min(centre, EXACT_BETA_C_INF - 0.005)  # don't overshoot true β_c
    return max(0.01, centre - margin), centre + margin


def traj_for_L(L, base):
    """Scale trajectory count with volume, capped to avoid excessive runtime."""
    scale = max(1, (L / 8) ** 2)
    # Cap at 4× base to keep large-L runs tractable
    return int(min(base * scale, base * 4))


# ---------------------------------------------------------------------------
# Method definitions (same logic as run_betac_diagnostic.py)
# ---------------------------------------------------------------------------

def _gc_refit(sb, sc, beta_init):
    """Fit Gram-Charlier to scan data and return mode."""
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    sb_arr, sc_arr = np.array(sb), np.array(sc)
    b_rng = sb_arr.max() - sb_arr.min()
    p0 = [beta_init, b_rng / 2, 0.0, 0.0, float(sc_arr.max())]
    bounds = ([sb_arr.min(), 1e-8, -5.0, -3.0, 0.0],
              [sb_arr.max(), b_rng * 2, 5.0, 10.0, np.inf])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        popt, _ = curve_fit(_gram_charlier, sb_arr, sc_arr,
                            p0=p0, bounds=bounds, maxfev=8000)
    mu_fit, sigma_fit = popt[0], popt[1]
    mode_lo = max(sb_arr.min(), mu_fit - 3.0 * sigma_fit)
    mode_hi = min(sb_arr.max(), mu_fit + 3.0 * sigma_fit)
    res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                          bounds=(mode_lo, mode_hi), method='bounded')
    return float(np.clip(res.x, sb_arr.min(), sb_arr.max()))


def find_beta_gc2(L, trial_dir, n_coarse, n_fine):
    """Standard 2-pass Gram-Charlier."""
    blo, bhi = scan_window(L)
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, L, L, 0, 0, k1, k2, k3,
        blo, bhi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=n_coarse, n_traj_fine=n_fine,
        data_dir=scan_dir,
    )
    return beta_c, sb, sc, se


def find_beta_gc3(L, trial_dir, n_coarse, n_fine):
    """3-pass Gram-Charlier: standard 2-pass + extra tight refinement."""
    blo, bhi = scan_window(L)
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c_2, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, L, L, 0, 0, k1, k2, k3,
        blo, bhi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=n_coarse, n_traj_fine=n_fine,
        data_dir=scan_dir,
    )
    # Pass 3: very tight bracket
    b_range = max(sb) - min(sb)
    half3 = b_range / (2 * 11) / 4
    fine3_betas = np.linspace(beta_c_2 - half3, beta_c_2 + half3, 7)
    for b in fine3_betas:
        scan_dir3 = os.path.join(trial_dir, "scan_p3")
        stdout_f, _ = oc.run_simulator(
            EXE, L, L, 0, 0, k1, k2, k3, float(b),
            n_traj=n_fine, data_dir=scan_dir3)
        chi_f, chi_f_err = oc.parse_stdout_value_with_err(stdout_f, "m_susc:")
        sb.append(float(b))
        sc.append(chi_f)
        se.append(chi_f_err)

    try:
        beta_c_3 = _gc_refit(sb, sc, beta_c_2)
    except Exception:
        beta_c_3 = beta_c_2
    return beta_c_3, sb, sc, se


def find_beta_gc2_histat(L, trial_dir, n_coarse, n_fine):
    """2-pass GC with 2× trajectories."""
    blo, bhi = scan_window(L)
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, L, L, 0, 0, k1, k2, k3,
        blo, bhi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=n_coarse * 2, n_traj_fine=n_fine * 2,
        data_dir=scan_dir,
    )
    return beta_c, sb, sc, se


def find_beta_gc2_boot(L, trial_dir, n_coarse, n_fine):
    """2-pass GC + bootstrap median."""
    blo, bhi = scan_window(L)
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c_orig, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, L, L, 0, 0, k1, k2, k3,
        blo, bhi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=n_coarse, n_traj_fine=n_fine,
        data_dir=scan_dir,
    )
    # Bootstrap
    sb_arr, sc_arr = np.array(sb), np.array(sc)
    boot_modes = [beta_c_orig]
    rng = np.random.default_rng()
    for _ in range(10):
        idx = rng.choice(len(sb_arr), size=len(sb_arr), replace=True)
        try:
            mode = _gc_refit(sb_arr[idx], sc_arr[idx], beta_c_orig)
            boot_modes.append(mode)
        except Exception:
            boot_modes.append(beta_c_orig)
    return float(np.median(boot_modes)), sb, sc, se


METHODS = {
    "gc2":        ("GC 2-pass (status quo)", find_beta_gc2),
    "gc3":        ("GC 3-pass", find_beta_gc3),
    "gc2_histat": ("GC 2-pass hi-stat (2×)", find_beta_gc2_histat),
    "gc2_boot":   ("GC2 + bootstrap", find_beta_gc2_boot),
}


# ---------------------------------------------------------------------------
# Resume support
# ---------------------------------------------------------------------------
RESULTS_FILE = os.path.join(OUTPUT_DIR, "fss_results.json")

results_flat = []
done_keys = set()

if os.path.exists(RESULTS_FILE):
    with open(RESULTS_FILE) as f:
        results_flat = json.load(f)
    for r in results_flat:
        done_keys.add((r["method"], r["L"], r["trial"]))
    print(f"Resuming: {len(results_flat)} entries loaded, "
          f"{len(done_keys)} (method,L,trial) combos done")


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
total_combos = len(METHODS) * len(LATTICE_SIZES) * N_TRIALS
done_count = len(done_keys)
print(f"FSS β_c test: {len(METHODS)} methods × {len(LATTICE_SIZES)} sizes × "
      f"{N_TRIALS} trials = {total_combos} total")
print(f"  Lattice sizes: {LATTICE_SIZES}")
print(f"  Methods: {list(METHODS.keys())}")
print(f"  Exact β_c(∞) = ln(3)/4 = {EXACT_BETA_C_INF:.12f}")
print(f"  Already done: {done_count}/{total_combos}")
print()

for L in LATTICE_SIZES:
    n_coarse = traj_for_L(L, BASE_N_TRAJ_COARSE)
    n_fine   = traj_for_L(L, BASE_N_TRAJ_FINE)
    print(f"=== L = {L} ===  (coarse: {n_coarse}, fine: {n_fine})")

    for trial in range(N_TRIALS):
        for method_name, (desc, finder_fn) in METHODS.items():
            if (method_name, L, trial) in done_keys:
                continue

            trial_dir = os.path.join(OUTPUT_DIR, "runs",
                                     f"L{L}", method_name, f"trial{trial:03d}")
            os.makedirs(trial_dir, exist_ok=True)

            t0 = time.time()
            print(f"  L={L:2d} [{trial+1:2d}/{N_TRIALS}] {method_name:15s} ...",
                  end="", flush=True)

            try:
                beta_c, sb, sc, se = finder_fn(L, trial_dir, n_coarse, n_fine)
                elapsed = time.time() - t0

                result = {
                    "method": method_name,
                    "L": L,
                    "trial": trial,
                    "beta_c": beta_c,
                    "n_scan_points": len(sb),
                    "elapsed_s": round(elapsed, 1),
                }
                results_flat.append(result)
                done_keys.add((method_name, L, trial))

                bias = beta_c - EXACT_BETA_C_INF
                print(f"  β_c={beta_c:.8f}  Δ={bias:+.6f}  ({elapsed:.0f}s)")
                sys.stdout.flush()

            except Exception as e:
                elapsed = time.time() - t0
                print(f"  FAILED: {e}  ({elapsed:.0f}s)")
                sys.stdout.flush()
                continue

        # Save after each trial across all methods
        with open(RESULTS_FILE, "w") as f:
            json.dump(results_flat, f, indent=2)


# ---------------------------------------------------------------------------
# Analysis and plotting
# ---------------------------------------------------------------------------
print("\n" + "=" * 90)
print("FINITE-SIZE SCALING ANALYSIS")
print("=" * 90)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit as scipy_curve_fit

method_colors = {
    "gc2": "#3498db", "gc3": "#e74c3c",
    "gc2_histat": "#e67e22", "gc2_boot": "#9b59b6",
}
method_labels = {
    "gc2": "GC 2-pass", "gc3": "GC 3-pass",
    "gc2_histat": "GC2 hi-stat (2×)", "gc2_boot": "GC2 + bootstrap",
}


# Collect per-method, per-L statistics
fss_data = {}  # method -> {L: [beta_c values]}
for r in results_flat:
    m, L = r["method"], r["L"]
    fss_data.setdefault(m, {}).setdefault(L, []).append(r["beta_c"])

# Print summary table
print(f"\n{'Method':15s}  {'L':>4s}  {'N':>3s}  {'β_c mean':>12s}  {'std':>10s}  "
      f"{'Δ from exact':>12s}  {'std/√N':>10s}")
print("-" * 80)
for m in METHODS:
    if m not in fss_data:
        continue
    for L in sorted(fss_data[m].keys()):
        vals = np.array(fss_data[m][L])
        mn = vals.mean()
        sd = vals.std(ddof=1) if len(vals) > 1 else 0
        se = sd / np.sqrt(len(vals)) if len(vals) > 1 else 0
        delta = mn - EXACT_BETA_C_INF
        print(f"{m:15s}  {L:4d}  {len(vals):3d}  {mn:.10f}  {sd:.8f}  "
              f"{delta:+.10f}  {se:.8f}")
    print()


# ---------------------------------------------------------------------------
# FSS fits:  β_c(L) = β_c(∞) + a / L^p
# We try two models:
#   (A) Fixed ν=1 (2D Ising):  β_c(L) = β_c_inf + a/L
#   (B) Free exponent:         β_c(L) = β_c_inf + a/L^p
# ---------------------------------------------------------------------------

def fss_model_fixed(L, beta_inf, a):
    """FSS with fixed exponent 1/ν = 1 (2D Ising)."""
    return beta_inf + a / np.asarray(L, dtype=float)

def fss_model_free(L, beta_inf, a, p):
    """FSS with free exponent."""
    return beta_inf + a / np.asarray(L, dtype=float) ** p

fit_results = {}
print("\nFSS EXTRAPOLATION TO L → ∞")
print("=" * 90)
print(f"{'Method':15s}  {'Model':>10s}  {'β_c(∞)':>14s}  {'a':>10s}  {'p':>6s}  "
      f"{'Δ from exact':>14s}  {'χ²/ndof':>10s}")
print("-" * 90)

for m in METHODS:
    if m not in fss_data:
        continue
    Ls = sorted(fss_data[m].keys())
    if len(Ls) < 3:
        print(f"{m:15s}  -- not enough sizes ({len(Ls)}) for FSS fit --")
        continue

    L_arr = np.array(Ls, dtype=float)
    bc_means = np.array([np.mean(fss_data[m][L]) for L in Ls])
    bc_stds = np.array([np.std(fss_data[m][L], ddof=1) / np.sqrt(len(fss_data[m][L]))
                        if len(fss_data[m][L]) > 1 else 1e-4 for L in Ls])

    # Model A: fixed exponent (1/ν = 1)
    try:
        popt_a, pcov_a = scipy_curve_fit(fss_model_fixed, L_arr, bc_means,
                                         sigma=bc_stds, absolute_sigma=True,
                                         p0=[EXACT_BETA_C_INF, -1.0])
        resid = (bc_means - fss_model_fixed(L_arr, *popt_a)) / bc_stds
        chi2_ndof_a = np.sum(resid**2) / max(1, len(Ls) - 2)
        perr_a = np.sqrt(np.diag(pcov_a))
        delta_a = popt_a[0] - EXACT_BETA_C_INF
        print(f"{m:15s}  {'fixed p=1':>10s}  {popt_a[0]:.10f}  {popt_a[1]:+.6f}  {'1.00':>6s}  "
              f"{delta_a:+.10f}  {chi2_ndof_a:.4f}")
        fit_results[(m, "fixed")] = {
            "beta_inf": popt_a[0], "beta_inf_err": perr_a[0],
            "a": popt_a[1], "p": 1.0, "chi2_ndof": chi2_ndof_a,
        }
    except Exception as e:
        print(f"{m:15s}  {'fixed p=1':>10s}  FAILED: {e}")

    # Model B: free exponent
    try:
        popt_b, pcov_b = scipy_curve_fit(fss_model_free, L_arr, bc_means,
                                         sigma=bc_stds, absolute_sigma=True,
                                         p0=[EXACT_BETA_C_INF, -1.0, 1.0],
                                         bounds=([0.20, -10.0, 0.1],
                                                 [0.35, 10.0, 5.0]))
        resid = (bc_means - fss_model_free(L_arr, *popt_b)) / bc_stds
        chi2_ndof_b = np.sum(resid**2) / max(1, len(Ls) - 3)
        perr_b = np.sqrt(np.diag(pcov_b))
        delta_b = popt_b[0] - EXACT_BETA_C_INF
        print(f"{m:15s}  {'free p':>10s}  {popt_b[0]:.10f}  {popt_b[1]:+.6f}  "
              f"{popt_b[2]:.3f}  {delta_b:+.10f}  {chi2_ndof_b:.4f}")
        fit_results[(m, "free")] = {
            "beta_inf": popt_b[0], "beta_inf_err": perr_b[0],
            "a": popt_b[1], "p": popt_b[2], "chi2_ndof": chi2_ndof_b,
        }
    except Exception as e:
        print(f"{m:15s}  {'free p':>10s}  FAILED: {e}")

    print()


# ---------------------------------------------------------------------------
# Save summary
# ---------------------------------------------------------------------------
summary = {
    "exact_beta_c_inf": EXACT_BETA_C_INF,
    "lattice_sizes": LATTICE_SIZES,
    "n_trials": N_TRIALS,
    "methods": {},
}
for m in METHODS:
    if m not in fss_data:
        continue
    m_summary = {"per_L": {}}
    for L in sorted(fss_data[m].keys()):
        vals = np.array(fss_data[m][L])
        m_summary["per_L"][str(L)] = {
            "mean": float(vals.mean()),
            "std": float(vals.std(ddof=1)) if len(vals) > 1 else 0,
            "n": len(vals),
        }
    for model in ["fixed", "free"]:
        key = (m, model)
        if key in fit_results:
            m_summary[f"fss_{model}"] = fit_results[key]
    summary["methods"][m] = m_summary

with open(os.path.join(OUTPUT_DIR, "fss_summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=str)


# ---------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(20, 14))
fig.suptitle(
    f"Finite-Size Scaling β_c Test — {N_TRIALS} trials/point, "
    f"equilateral triangular Ising",
    fontsize=14, fontweight="bold",
)
gs = GridSpec(2, 3, hspace=0.32, wspace=0.30,
              left=0.06, right=0.97, top=0.92, bottom=0.06)


# Panel 1: β_c(L) vs L
ax1 = fig.add_subplot(gs[0, 0])
for m in METHODS:
    if m not in fss_data:
        continue
    Ls = sorted(fss_data[m].keys())
    means = [np.mean(fss_data[m][L]) for L in Ls]
    sems = [np.std(fss_data[m][L], ddof=1) / np.sqrt(len(fss_data[m][L]))
            if len(fss_data[m][L]) > 1 else 0 for L in Ls]
    ax1.errorbar(Ls, means, yerr=sems, fmt="o-", color=method_colors[m],
                 label=method_labels[m], markersize=5, capsize=3, linewidth=1.5)
ax1.axhline(EXACT_BETA_C_INF, ls="--", color="k", lw=1.5, alpha=0.7,
            label=f"ln(3)/4 = {EXACT_BETA_C_INF:.6f}")
ax1.set_xlabel("L (lattice size)")
ax1.set_ylabel("β_c(L)")
ax1.set_title("β_c(L) vs lattice size")
ax1.legend(fontsize=7, loc="lower right")


# Panel 2: β_c(L) vs 1/L with FSS fits
ax2 = fig.add_subplot(gs[0, 1])
inv_L_fine = np.linspace(0, 1.0 / min(LATTICE_SIZES) + 0.02, 100)
for m in METHODS:
    if m not in fss_data:
        continue
    Ls = sorted(fss_data[m].keys())
    inv_Ls = [1.0 / L for L in Ls]
    means = [np.mean(fss_data[m][L]) for L in Ls]
    sems = [np.std(fss_data[m][L], ddof=1) / np.sqrt(len(fss_data[m][L]))
            if len(fss_data[m][L]) > 1 else 0 for L in Ls]
    ax2.errorbar(inv_Ls, means, yerr=sems, fmt="o", color=method_colors[m],
                 label=method_labels[m], markersize=6, capsize=3)
    # Plot fixed-exponent fit if available
    key = (m, "fixed")
    if key in fit_results:
        fr = fit_results[key]
        ax2.plot(inv_L_fine, fr["beta_inf"] + fr["a"] * inv_L_fine,
                 "--", color=method_colors[m], alpha=0.6, lw=1.5)
ax2.axhline(EXACT_BETA_C_INF, ls="-", color="k", lw=1.5, alpha=0.5)
ax2.axvline(0, ls=":", color="gray", lw=0.8, alpha=0.5)
ax2.set_xlabel("1/L")
ax2.set_ylabel("β_c(L)")
ax2.set_title("FSS extrapolation (fixed ν=1)")
ax2.legend(fontsize=7)
# Mark the L→∞ extrapolations
for m in METHODS:
    key = (m, "fixed")
    if key in fit_results:
        fr = fit_results[key]
        ax2.plot(0, fr["beta_inf"], "s", color=method_colors[m],
                 markersize=10, markeredgecolor="k", zorder=5)


# Panel 3: Extrapolated β_c(∞) comparison (bar chart)
ax3 = fig.add_subplot(gs[0, 2])
bar_methods = [m for m in METHODS if (m, "fixed") in fit_results]
extraps = [fit_results[(m, "fixed")]["beta_inf"] for m in bar_methods]
extrap_errs = [fit_results[(m, "fixed")]["beta_inf_err"] for m in bar_methods]
x = np.arange(len(bar_methods))
bars = ax3.bar(x, extraps, yerr=extrap_errs, capsize=5, width=0.5,
               color=[method_colors[m] for m in bar_methods], alpha=0.7,
               edgecolor="k", linewidth=0.5)
ax3.axhline(EXACT_BETA_C_INF, ls="--", color="k", lw=2,
            label=f"Exact = {EXACT_BETA_C_INF:.8f}")
ax3.set_xticks(x)
ax3.set_xticklabels([method_labels[m] for m in bar_methods], fontsize=8,
                     rotation=20, ha="right")
ax3.set_ylabel("Extrapolated β_c(∞)")
ax3.set_title("β_c(∞) from FSS  (closer to line = better)")
ax3.legend(fontsize=8)
for i, (v, e) in enumerate(zip(extraps, extrap_errs)):
    delta = v - EXACT_BETA_C_INF
    ax3.text(i, v + e + 0.0003, f"Δ={delta:+.5f}", ha="center",
             va="bottom", fontsize=7, fontweight="bold")


# Panel 4: β_c scatter (std) vs L for each method
ax4 = fig.add_subplot(gs[1, 0])
for m in METHODS:
    if m not in fss_data:
        continue
    Ls = sorted(fss_data[m].keys())
    stds = [np.std(fss_data[m][L], ddof=1) if len(fss_data[m][L]) > 1 else 0
            for L in Ls]
    ax4.plot(Ls, stds, "o-", color=method_colors[m], label=method_labels[m],
             markersize=5, linewidth=1.5)
ax4.set_xlabel("L")
ax4.set_ylabel("β_c std (over trials)")
ax4.set_title("Finder precision vs lattice size")
ax4.legend(fontsize=7)
ax4.set_yscale("log")


# Panel 5: Residuals from FSS fit (fixed ν=1)
ax5 = fig.add_subplot(gs[1, 1])
for m in METHODS:
    if m not in fss_data or (m, "fixed") not in fit_results:
        continue
    fr = fit_results[(m, "fixed")]
    Ls = sorted(fss_data[m].keys())
    L_arr = np.array(Ls, dtype=float)
    means = np.array([np.mean(fss_data[m][L]) for L in Ls])
    sems = np.array([np.std(fss_data[m][L], ddof=1) / np.sqrt(len(fss_data[m][L]))
                     if len(fss_data[m][L]) > 1 else 1e-4 for L in Ls])
    predicted = fss_model_fixed(L_arr, fr["beta_inf"], fr["a"])
    resid = (means - predicted) / sems
    ax5.plot(Ls, resid, "o-", color=method_colors[m], label=method_labels[m],
             markersize=6, linewidth=1.5)
ax5.axhline(0, ls="-", color="k", lw=0.8)
ax5.axhline(2, ls=":", color="gray", lw=0.8, alpha=0.5)
ax5.axhline(-2, ls=":", color="gray", lw=0.8, alpha=0.5)
ax5.set_xlabel("L")
ax5.set_ylabel("Normalised residual (σ)")
ax5.set_title("FSS fit residuals (fixed ν=1)")
ax5.legend(fontsize=7)


# Panel 6: Time cost vs L
ax6 = fig.add_subplot(gs[1, 2])
for m in METHODS:
    if m not in fss_data:
        continue
    time_data = {}
    for r in results_flat:
        if r["method"] == m:
            time_data.setdefault(r["L"], []).append(r["elapsed_s"])
    Ls = sorted(time_data.keys())
    mean_t = [np.mean(time_data[L]) for L in Ls]
    ax6.plot(Ls, mean_t, "o-", color=method_colors[m], label=method_labels[m],
             markersize=5, linewidth=1.5)
ax6.set_xlabel("L")
ax6.set_ylabel("Mean wall time per trial (s)")
ax6.set_title("Computational cost")
ax6.legend(fontsize=7)
ax6.set_yscale("log")


out_path = os.path.join(OUTPUT_DIR, "fss_betac_plots.png")
fig.savefig(out_path, dpi=150)
print(f"\nPlot saved: {out_path}")
plt.close()

print(f"\nResults saved: {RESULTS_FILE}")
print(f"Summary saved: {os.path.join(OUTPUT_DIR, 'fss_summary.json')}")
print("\nDone.")
