#!/usr/bin/env python3
"""
run_stability_normed_compare.py
===============================
Side-by-side stability comparison of three cost functions:

    1. boundary_slices           — raw quartic ∫(diff⁴)dt
    2. boundary_slices_normed    — normalised quadratic ∫(Z²)dt
    3. boundary_slices_normed_quartic — normalised quartic ∫(Z⁴)dt

Each trial generates ONE MC dataset (same beta_c scan + production run)
and evaluates all three costs on it, so the comparison is perfectly paired.
"""

import os, sys, json, time
import numpy as np

# ===========================================================================
# CONFIG
# ===========================================================================
Lx, Ly   = 8, 8
Tx, Ty   = 0, 0
r1, r2   = 1.0, 1.0      # equilateral point
BETA_GUESS = 0.24

N_TRIALS = 100

N_TRAJ_PROD        = 50000
N_TRAJ_SCAN_COARSE = 50000
N_TRAJ_SCAN_FINE   = 100000

R_MIN      = 0.0
R_MAX_FRAC = 0.33

REF_META = "ref_data/ref_8x8/ref_metadata.json"
EXE = "bin/ising_tri_twisted_parallelogram"
OUTPUT_DIR = "results/stability_normed_compare"

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

# Null plotter
class NullPlotter:
    def __init__(self):
        self.history_chi2_ndof = []
        self.history_chi2_err = []
        self.history_beta_c = []
    def record_point(self, *a, **kw): pass
    def set_current_point(self, *a, **kw): pass
    def start_level(self, *a, **kw): pass
    def _update(self, *a, **kw): pass

plotter = NullPlotter()

k3 = 1.0
k1, k2 = r1 * k3, r2 * k3

data_dir = os.path.join(OUTPUT_DIR, "mc_data")
os.makedirs(data_dir, exist_ok=True)

results = []
print(f"Normed cost comparison: {N_TRIALS} evaluations at (r1={r1}, r2={r2})")
print(f"  Lattice: {Lx}x{Ly}  stats: {N_TRAJ_PROD}/{N_TRAJ_SCAN_COARSE}/{N_TRAJ_SCAN_FINE}")
print(f"  costs: {COSTS}")
print()

for trial in range(N_TRIALS):
    t0 = time.time()
    label = f"trial_{trial:03d}"
    trial_dir = os.path.join(data_dir, label)

    # Run the full pipeline with the first cost to generate MC data
    # (beta_c scan + production). Subsequent costs reuse the SAME data.
    beta_c, chi2_raw, ndof_raw, chi2_ndof_raw = oc.evaluate_point(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3, BETA_GUESS,
        N_TRAJ_PROD, trial_dir, ref_data, L_eff, label, plotter,
        n_traj_scan_coarse=N_TRAJ_SCAN_COARSE,
        n_traj_scan_fine=N_TRAJ_SCAN_FINE,
        r_min=R_MIN, r_max_frac=R_MAX_FRAC,
        cost="boundary_slices", beta_ref=beta_ref,
        ref_Lx=ref_Lx, ref_Ly=ref_Ly,
        ref_Tx=ref_Tx, ref_Ty=ref_Ty,
    )

    # Load the production data that was just generated
    # evaluate_point stores production in {trial_dir}/prod_{label}/{subdir}/two_point_all_to_all.dat
    prod_parent = os.path.join(trial_dir, f"prod_{label}")
    # find the a2a file inside the simulator's output subdirectory
    prod_a2a = None
    for root, dirs, files in os.walk(prod_parent):
        if "two_point_all_to_all.dat" in files:
            prod_a2a = os.path.join(root, "two_point_all_to_all.dat")
            break
    if prod_a2a is None:
        print(f"  WARNING: could not find a2a data for trial {trial}, skipping")
        continue
    test_data = oc.load_all_to_all(prod_a2a)

    # Evaluate the other two costs on the SAME test_data
    trial_result = {
        "trial": trial,
        "beta_c": beta_c,
    }
    for cost_name in COSTS:
        if cost_name == "boundary_slices":
            trial_result[cost_name] = {
                "chi2": chi2_raw, "ndof": ndof_raw, "chi2_ndof": chi2_ndof_raw,
            }
        else:
            chi2_alt, ndof_alt, chi2_err_alt = oc.compute_chi2(
                ref_data, test_data,
                r_min=R_MIN, r_max_frac=R_MAX_FRAC, L_eff=L_eff,
                cost=cost_name,
                Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
                ref_Lx=ref_Lx, ref_Ly=ref_Ly,
                ref_Tx=ref_Tx, ref_Ty=ref_Ty,
            )
            trial_result[cost_name] = {
                "chi2": chi2_alt, "ndof": ndof_alt,
                "chi2_ndof": chi2_alt / max(ndof_alt, 1),
            }

    elapsed = time.time() - t0
    trial_result["elapsed_s"] = round(elapsed, 1)
    results.append(trial_result)

    parts = "  ".join(
        f"{c}={trial_result[c]['chi2_ndof']:.4e}" for c in COSTS
    )
    print(f"  [{trial+1:3d}/{N_TRIALS}]  beta_c={beta_c:.8f}  {parts}  ({elapsed:.0f}s)")
    sys.stdout.flush()

    # Save incrementally
    with open(os.path.join(OUTPUT_DIR, "compare_results.json"), "w") as f:
        json.dump(results, f, indent=2)

# =================== Summary & Plots ===================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print(f"\n{'='*70}")
print(f"Comparison complete: {N_TRIALS} trials at (r1={r1}, r2={r2}), {Lx}x{Ly}")
for cost_name in COSTS:
    vals = np.array([r[cost_name]["chi2_ndof"] for r in results])
    cov = vals.std() / vals.mean() if vals.mean() != 0 else float("inf")
    print(f"  {cost_name:35s}  mean={vals.mean():.4e}  std={vals.std():.4e}  "
          f"CoV={cov:.3f}  min={vals.min():.4e}  max={vals.max():.4e}")
print(f"{'='*70}")

# --- Histogram comparison ---
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
colours = ["steelblue", "darkorange", "seagreen"]
short_names = ["raw quartic", "normed Z²", "normed Z⁴"]

for ax, cost_name, colour, sname in zip(axes, COSTS, colours, short_names):
    vals = np.array([r[cost_name]["chi2_ndof"] for r in results])
    cov = vals.std() / vals.mean() if vals.mean() != 0 else float("inf")
    ax.hist(np.log10(np.maximum(vals, 1e-30)), bins=25, color=colour,
            edgecolor="white", alpha=0.85)
    ax.set_xlabel("log₁₀(chi2/ndof)")
    ax.set_ylabel("Count")
    ax.set_title(f"{sname}\nCoV={cov:.3f}")

fig.suptitle(f"Cost stability: {N_TRIALS} trials, (r1={r1}, r2={r2}), "
             f"{Lx}x{Ly}, {N_TRAJ_PROD} traj", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "normed_compare_histograms.png"), dpi=150)
print(f"\nPlot saved: {OUTPUT_DIR}/normed_compare_histograms.png")

# --- Paired scatter: normed costs vs raw ---
fig2, axes2 = plt.subplots(1, 2, figsize=(11, 5))
raw_vals = np.array([r["boundary_slices"]["chi2_ndof"] for r in results])

for ax, cost_name, colour, sname in zip(
    axes2,
    ["boundary_slices_normed", "boundary_slices_normed_quartic"],
    ["darkorange", "seagreen"],
    ["normed Z²", "normed Z⁴"],
):
    alt_vals = np.array([r[cost_name]["chi2_ndof"] for r in results])
    ax.scatter(np.log10(np.maximum(raw_vals, 1e-30)),
               np.log10(np.maximum(alt_vals, 1e-30)),
               c=colour, s=20, alpha=0.6, edgecolors="k", linewidths=0.3)
    ax.set_xlabel("log₁₀(raw quartic)")
    ax.set_ylabel(f"log₁₀({sname})")
    ax.set_title(f"Paired: raw quartic vs {sname}")
    ax.grid(True, alpha=0.3)
    # Correlation
    corr = np.corrcoef(np.log10(np.maximum(raw_vals, 1e-30)),
                        np.log10(np.maximum(alt_vals, 1e-30)))[0, 1]
    ax.text(0.05, 0.95, f"r={corr:.3f}", transform=ax.transAxes,
            va="top", fontsize=10, bbox=dict(boxstyle="round", fc="white", alpha=0.8))

fig2.suptitle("Paired cost correlation (same MC data)", fontsize=11)
fig2.tight_layout()
fig2.savefig(os.path.join(OUTPUT_DIR, "normed_compare_scatter.png"), dpi=150)
print(f"Plot saved: {OUTPUT_DIR}/normed_compare_scatter.png")
