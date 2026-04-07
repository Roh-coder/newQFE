#!/usr/bin/env python3
"""
run_betac_diagnostic.py
=======================
Diagnose how much of the cost-function noise comes from the beta_c finder
vs. pure MC production noise.

Tests 5 methods side-by-side, all on the SAME lattice point (r1=r2=1, 8x8):

  1. exact     — use the known exact beta_c (bypass finder entirely)
  2. gc2       — current 2-pass GC (status quo)
  3. gc3       — 3-pass GC (extra refinement layer)
  4. gc2_histat — 2-pass GC with 2x scan trajectories
  5. gc2_boot  — 2-pass GC + bootstrap median (10 resamples of scan data)

For each method we run N_TRIALS independent evaluations, computing all
three cost functions (raw quartic, normed Z², normed Z⁴) on the same
production data.  This gives a perfectly controlled comparison.
"""

import os, sys, json, time, math, copy
import numpy as np

# ===========================================================================
# CONFIG
# ===========================================================================
Lx, Ly   = 8, 8
Tx, Ty   = 0, 0
r1, r2   = 1.0, 1.0
# The "exact" control uses the reference's own beta_c — the finite-size
# susceptibility peak at 8x8.  This skips the finder entirely and measures
# the noise floor from MC production alone.
EXACT_BETA_C = None   # will be set from reference metadata

N_TRIALS = 30       # per method (total = 5 × 30 = 150 evaluations)

N_TRAJ_PROD        = 50000
N_TRAJ_SCAN_COARSE = 50000
N_TRAJ_SCAN_FINE   = 100000

BETA_GUESS = 0.24
R_MIN      = 0.0
R_MAX_FRAC = 0.33

REF_META = "ref_data/ref_8x8/ref_metadata.json"
EXE = "bin/ising_tri_twisted_parallelogram"
OUTPUT_DIR = "results/betac_diagnostic"

COSTS = [
    "boundary_slices",
    "boundary_slices_normed",
    "boundary_slices_normed_quartic",
]
# ===========================================================================

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
os.chdir(_HERE)

import optimise_couplings as oc

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load reference
with open(REF_META) as f:
    meta = json.load(f)
ref_data = oc.load_all_to_all(meta["a2a_path"])
beta_ref = meta["beta_c"]
ref_Lx, ref_Ly = meta["Lx"], meta["Ly"]
ref_Tx, ref_Ty = meta.get("Tx", 0), meta.get("Ty", 0)
L_eff = min(Lx, Ly)

# Use the reference's finite-size beta_c as the "exact" control
EXACT_BETA_C = meta["beta_c"]
print(f"  Reference beta_c (8×8 finite-size): {EXACT_BETA_C:.10f}")

k3 = 1.0
k1_val, k2_val = r1 * k3, r2 * k3

margin = max(0.20 * BETA_GUESS, 0.04)
BETA_LO = max(0.01, BETA_GUESS - margin)
BETA_HI = BETA_GUESS + margin


# ---------------------------------------------------------------------------
# Beta_c finder variants
# ---------------------------------------------------------------------------

def find_beta_exact(trial_dir):
    """Control: use exact beta_c, no scan at all."""
    return EXACT_BETA_C, [], [], []


def find_beta_gc2(trial_dir):
    """Status quo: 2-pass GC, current stats."""
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3,
        BETA_LO, BETA_HI,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE, n_traj_fine=N_TRAJ_SCAN_FINE,
        data_dir=scan_dir,
    )
    return beta_c, sb, sc, se


def find_beta_gc3(trial_dir):
    """3-pass GC: run a third tighter refinement pass."""
    scan_dir = os.path.join(trial_dir, "scan")
    # Pass 1+2 via the standard finder
    beta_c_2, chi_peak, sb, sc, se, parabola = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3,
        BETA_LO, BETA_HI,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE, n_traj_fine=N_TRAJ_SCAN_FINE,
        data_dir=scan_dir,
    )
    # Pass 3: even tighter bracket, 5 more fine-stats points
    b_range = max(sb) - min(sb)
    half3 = b_range / (2 * 11) / 4  # quarter of the pass-2 half-width
    fine3_betas = np.linspace(beta_c_2 - half3, beta_c_2 + half3, 7)
    for b in fine3_betas:
        scan_dir3 = os.path.join(trial_dir, "scan_p3")
        stdout_f, _ = oc.run_simulator(
            EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3, float(b),
            n_traj=N_TRAJ_SCAN_FINE, data_dir=scan_dir3)
        chi_f, chi_f_err = oc.parse_stdout_value_with_err(stdout_f, "m_susc:")
        sb.append(float(b))
        sc.append(chi_f)
        se.append(chi_f_err)

    # Final GC fit on ALL data
    from scipy.optimize import curve_fit, minimize_scalar
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    sb_arr, sc_arr = np.array(sb), np.array(sc)
    b_rng = sb_arr.max() - sb_arr.min()
    try:
        import warnings
        p0 = [beta_c_2, b_rng / 2, 0.0, 0.0, float(sc_arr.max())]
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
        beta_c_3 = float(np.clip(res.x, sb_arr.min(), sb_arr.max()))
    except Exception:
        beta_c_3 = beta_c_2  # fallback to 2-pass result

    return beta_c_3, sb, sc, se


def find_beta_gc2_histat(trial_dir):
    """2-pass GC with 2x scan trajectories."""
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3,
        BETA_LO, BETA_HI,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE * 2,
        n_traj_fine=N_TRAJ_SCAN_FINE * 2,
        data_dir=scan_dir,
    )
    return beta_c, sb, sc, se


def find_beta_gc2_bootstrap(trial_dir):
    """2-pass GC + bootstrap: run the standard scan, then resample scan data
    10 times with replacement, fit GC each time, take the median mode."""
    scan_dir = os.path.join(trial_dir, "scan")
    beta_c_orig, chi_peak, sb, sc, se, _ = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3,
        BETA_LO, BETA_HI,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE, n_traj_fine=N_TRAJ_SCAN_FINE,
        data_dir=scan_dir,
    )

    from scipy.optimize import curve_fit, minimize_scalar
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    sb_arr, sc_arr = np.array(sb), np.array(sc)
    n_pts = len(sb_arr)
    b_rng = sb_arr.max() - sb_arr.min()

    boot_modes = [beta_c_orig]
    N_BOOT = 10
    rng = np.random.default_rng()
    for _ in range(N_BOOT):
        idx = rng.choice(n_pts, size=n_pts, replace=True)
        b_boot = sb_arr[idx]
        c_boot = sc_arr[idx]
        # Sort for stable fitting
        order = np.argsort(b_boot)
        b_boot, c_boot = b_boot[order], c_boot[order]
        try:
            import warnings
            p0 = [beta_c_orig, b_rng / 2, 0.0, 0.0, float(c_boot.max())]
            bounds = ([b_boot.min(), 1e-8, -5.0, -3.0, 0.0],
                      [b_boot.max(), b_rng * 2, 5.0, 10.0, np.inf])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(_gram_charlier, b_boot, c_boot,
                                    p0=p0, bounds=bounds, maxfev=8000)
            mu_fit, sigma_fit = popt[0], popt[1]
            mode_lo = max(b_boot.min(), mu_fit - 3.0 * sigma_fit)
            mode_hi = min(b_boot.max(), mu_fit + 3.0 * sigma_fit)
            res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                                  bounds=(mode_lo, mode_hi), method='bounded')
            boot_modes.append(float(np.clip(res.x, b_boot.min(), b_boot.max())))
        except Exception:
            boot_modes.append(beta_c_orig)

    beta_c_median = float(np.median(boot_modes))
    return beta_c_median, sb, sc, se


METHODS = {
    "exact":       ("Exact β_c (control)",    find_beta_exact),
    "gc2":         ("2-pass GC (current)",     find_beta_gc2),
    "gc3":         ("3-pass GC",               find_beta_gc3),
    "gc2_histat":  ("2-pass GC, 2× stats",    find_beta_gc2_histat),
    "gc2_boot":    ("2-pass GC + bootstrap",   find_beta_gc2_bootstrap),
}

# ---------------------------------------------------------------------------
# Run production + evaluate costs
# ---------------------------------------------------------------------------

def evaluate_with_method(method_name, finder_fn, trial_idx):
    """Run one trial with the given beta_c finder method."""
    trial_dir = os.path.join(OUTPUT_DIR, "mc_data", method_name, f"trial_{trial_idx:03d}")
    os.makedirs(trial_dir, exist_ok=True)

    t0 = time.time()
    beta_c, sb, sc, se = finder_fn(trial_dir)
    t_scan = time.time() - t0

    # Production run at the found beta_c
    prod_dir = os.path.join(trial_dir, "prod")
    stdout, subdir = oc.run_simulator(
        EXE, Lx, Ly, Tx, Ty, k1_val, k2_val, k3, beta_c,
        n_traj=N_TRAJ_PROD, n_therm=3000,
        data_dir=prod_dir,
    )
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    test_data = oc.load_all_to_all(a2a_path)

    # Evaluate all three costs
    result = {"trial": trial_idx, "method": method_name, "beta_c": beta_c,
              "scan_time_s": round(t_scan, 1)}
    for cost_name in COSTS:
        chi2, ndof, chi2_err = oc.compute_chi2(
            ref_data, test_data,
            r_min=R_MIN, r_max_frac=R_MAX_FRAC, L_eff=L_eff,
            cost=cost_name,
            Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
            ref_Lx=ref_Lx, ref_Ly=ref_Ly,
            ref_Tx=ref_Tx, ref_Ty=ref_Ty,
        )
        result[cost_name] = {
            "chi2": chi2, "ndof": ndof,
            "chi2_ndof": chi2 / max(ndof, 1),
        }

    result["total_time_s"] = round(time.time() - t0, 1)
    return result


# ---------------------------------------------------------------------------
# Main loop: interleave methods for fairness
# ---------------------------------------------------------------------------

all_results = {m: [] for m in METHODS}
results_flat = []

# Resume from existing results if available
_resume_file = os.path.join(OUTPUT_DIR, "diagnostic_results.json")
_start_trial = 0
if os.path.exists(_resume_file):
    with open(_resume_file) as _f:
        _existing = json.load(_f)
    if _existing:
        results_flat = _existing
        for _r in _existing:
            m = _r["method"]
            if m in all_results:
                all_results[m].append(_r)
        _done_trials = set(_r["trial"] for _r in _existing)
        _start_trial = max(_done_trials) + 1
        print(f"Resuming from trial {_start_trial} ({len(_existing)} entries loaded)")

print(f"Beta_c diagnostic: {N_TRIALS} trials × {len(METHODS)} methods")
print(f"  Lattice: {Lx}×{Ly}  prod: {N_TRAJ_PROD}  scan: {N_TRAJ_SCAN_COARSE}/{N_TRAJ_SCAN_FINE}")
print(f"  Methods: {list(METHODS.keys())}")
print(f"  Exact beta_c (ref): {EXACT_BETA_C:.10f}")
print()

for trial in range(_start_trial, N_TRIALS):
    for method_name, (desc, finder_fn) in METHODS.items():
        print(f"  [{trial+1:2d}/{N_TRIALS}] {method_name:15s} ...", end="", flush=True)
        result = evaluate_with_method(method_name, finder_fn, trial)
        all_results[method_name].append(result)
        results_flat.append(result)

        raw = result["boundary_slices"]["chi2_ndof"]
        nq = result["boundary_slices_normed_quartic"]["chi2_ndof"]
        print(f"  beta_c={result['beta_c']:.8f}  raw={raw:.4e}  Z⁴={nq:.4e}  "
              f"({result['total_time_s']}s)")
        sys.stdout.flush()

    # Save incrementally after each round
    with open(os.path.join(OUTPUT_DIR, "diagnostic_results.json"), "w") as f:
        json.dump(results_flat, f, indent=2)

    # Print interim summary every 5 trials
    if (trial + 1) % 5 == 0 or trial == N_TRIALS - 1:
        print(f"\n{'='*90}")
        print(f"  Interim summary after {trial+1} trials:")
        print(f"  {'Method':20s}  {'β_c mean':>12s}  {'β_c std':>10s}  "
              f"{'raw CoV':>8s}  {'Z² CoV':>8s}  {'Z⁴ CoV':>8s}  {'scan_t':>7s}")
        print(f"  {'-'*85}")
        for method_name, (desc, _) in METHODS.items():
            res = all_results[method_name]
            if len(res) < 2:
                continue
            betas = np.array([r["beta_c"] for r in res])
            raw_v = np.array([r["boundary_slices"]["chi2_ndof"] for r in res])
            z2_v  = np.array([r["boundary_slices_normed"]["chi2_ndof"] for r in res])
            z4_v  = np.array([r["boundary_slices_normed_quartic"]["chi2_ndof"] for r in res])
            scan_t = np.mean([r["scan_time_s"] for r in res])
            def cov(v):
                return f"{v.std()/v.mean():.3f}" if v.mean() != 0 else "inf"
            print(f"  {desc:20s}  {betas.mean():12.8f}  {betas.std():10.8f}  "
                  f"{cov(raw_v):>8s}  {cov(z2_v):>8s}  {cov(z4_v):>8s}  {scan_t:6.0f}s")
        print(f"{'='*90}\n")


# =================== Final Summary & Plots ===================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print(f"\n{'='*90}")
print(f"FINAL RESULTS: {N_TRIALS} trials per method")
print(f"{'='*90}")

summary = {}
for method_name, (desc, _) in METHODS.items():
    res = all_results[method_name]
    betas = np.array([r["beta_c"] for r in res])
    s = {"desc": desc, "beta_c_mean": betas.mean(), "beta_c_std": betas.std()}
    for cost_name in COSTS:
        vals = np.array([r[cost_name]["chi2_ndof"] for r in res])
        cov = vals.std() / vals.mean() if vals.mean() != 0 else float("inf")
        s[cost_name] = {"mean": vals.mean(), "std": vals.std(), "cov": cov}
    summary[method_name] = s
    print(f"\n  {desc}:")
    print(f"    beta_c:  mean={betas.mean():.8f}  std={betas.std():.8f}")
    for cost_name in COSTS:
        v = s[cost_name]
        print(f"    {cost_name:35s}  mean={v['mean']:.4e}  CoV={v['cov']:.3f}")

with open(os.path.join(OUTPUT_DIR, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=str)

# --- Plot: beta_c distributions ---
fig, axes = plt.subplots(1, len(METHODS), figsize=(4*len(METHODS), 4), sharey=True)
for ax, (method_name, (desc, _)) in zip(axes, METHODS.items()):
    betas = np.array([r["beta_c"] for r in all_results[method_name]])
    ax.hist(betas, bins=15, color="steelblue", edgecolor="white")
    ax.axvline(EXACT_BETA_C, color="red", ls="--", lw=1.5, label=f"exact={EXACT_BETA_C:.6f}")
    ax.axvline(betas.mean(), color="orange", ls="-", lw=1.5, label=f"mean={betas.mean():.6f}")
    ax.set_xlabel(r"$\beta_c$")
    ax.set_title(f"{desc}\nstd={betas.std():.6f}", fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
axes[0].set_ylabel("Count")
fig.suptitle(f"β_c finder comparison ({N_TRIALS} trials, {Lx}×{Ly})", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "betac_distributions.png"), dpi=150)
print(f"\nPlot saved: {OUTPUT_DIR}/betac_distributions.png")

# --- Plot: CoV comparison bar chart ---
fig2, ax2 = plt.subplots(figsize=(10, 5))
x = np.arange(len(METHODS))
width = 0.25
colours = ["steelblue", "darkorange", "seagreen"]
short_costs = ["raw quartic", "normed Z²", "normed Z⁴"]
for i, (cost_name, colour, sname) in enumerate(zip(COSTS, colours, short_costs)):
    covs = [summary[m][cost_name]["cov"] for m in METHODS]
    ax2.bar(x + i*width, covs, width, label=sname, color=colour, edgecolor="white")
ax2.set_xticks(x + width)
ax2.set_xticklabels([METHODS[m][0] for m in METHODS], fontsize=8, rotation=15)
ax2.set_ylabel("Coefficient of Variation (CoV)")
ax2.set_title(f"Cost stability by β_c method ({N_TRIALS} trials)")
ax2.legend()
ax2.grid(True, alpha=0.3, axis="y")
ax2.axhline(0.1, ls=":", color="gray", alpha=0.5, label="target CoV = 0.1")
fig2.tight_layout()
fig2.savefig(os.path.join(OUTPUT_DIR, "cov_comparison.png"), dpi=150)
print(f"Plot saved: {OUTPUT_DIR}/cov_comparison.png")

# --- Plot: cost distributions (log scale), one row per cost, one col per method ---
fig3, axes3 = plt.subplots(len(COSTS), len(METHODS),
                            figsize=(4*len(METHODS), 3.5*len(COSTS)),
                            sharex='row', sharey='row')
for ci, (cost_name, sname) in enumerate(zip(COSTS, short_costs)):
    for mi, method_name in enumerate(METHODS):
        ax = axes3[ci, mi]
        vals = np.array([r[cost_name]["chi2_ndof"] for r in all_results[method_name]])
        lv = np.log10(np.maximum(vals, 1e-30))
        ax.hist(lv, bins=15, color=colours[ci], edgecolor="white", alpha=0.8)
        cov = vals.std() / vals.mean() if vals.mean() != 0 else float("inf")
        ax.set_title(f"{METHODS[method_name][0]}\nCoV={cov:.3f}", fontsize=8)
        if ci == len(COSTS) - 1:
            ax.set_xlabel("log₁₀(cost)")
        if mi == 0:
            ax.set_ylabel(f"{sname}\nCount")
fig3.suptitle(f"Cost distributions by method ({N_TRIALS} trials)", fontsize=11)
fig3.tight_layout()
fig3.savefig(os.path.join(OUTPUT_DIR, "cost_distributions.png"), dpi=150)
print(f"Plot saved: {OUTPUT_DIR}/cost_distributions.png")
