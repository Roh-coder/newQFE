#!/usr/bin/env python3
"""
run_stability_test.py
=====================
Repeatedly evaluate the SAME grid point to measure how stable the optimizer
score is under MC noise.  Each trial runs the full pipeline:

    beta_c scan  →  production run  →  boundary_slices cost

Results are written to a JSON file and a histogram is plotted.
"""

import os, sys, json, time
import numpy as np

# ===========================================================================
# CONFIG — edit here
# ===========================================================================
Lx, Ly   = 8, 8          # test lattice
Tx, Ty   = 0, 0
r1, r2   = 1.0, 1.0      # the equilateral point (should give cost ≈ 0)
BETA_GUESS = 0.24         # close to true beta_c ≈ 0.239 for 8×8

N_TRIALS = 100

# MC statistics (match the 10× run)
N_TRAJ_PROD        = 50000
N_TRAJ_SCAN_COARSE = 50000
N_TRAJ_SCAN_FINE   = 100000

COST       = "boundary_slices"
R_MIN      = 0.0
R_MAX_FRAC = 0.33

# Reference
REF_META = "ref_data/ref_8x8/ref_metadata.json"

# Paths
EXE = "bin/ising_tri_twisted_parallelogram"
OUTPUT_DIR = "results/stability_test_8x8"
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

# Minimal plotter that just collects data without rendering
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
print(f"Stability test: {N_TRIALS} evaluations at (r1={r1}, r2={r2})")
print(f"  Lattice: {Lx}x{Ly}  stats: {N_TRAJ_PROD}/{N_TRAJ_SCAN_COARSE}/{N_TRAJ_SCAN_FINE}")
print(f"  cost={COST}  ref={REF_META}")
print()

for trial in range(N_TRIALS):
    t0 = time.time()
    label = f"trial_{trial:03d}"
    trial_dir = os.path.join(data_dir, label)

    beta_c, chi2, ndof, chi2_ndof = oc.evaluate_point(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3, BETA_GUESS,
        N_TRAJ_PROD, trial_dir, ref_data, L_eff, label, plotter,
        n_traj_scan_coarse=N_TRAJ_SCAN_COARSE,
        n_traj_scan_fine=N_TRAJ_SCAN_FINE,
        r_min=R_MIN, r_max_frac=R_MAX_FRAC,
        cost=COST, beta_ref=beta_ref,
        ref_Lx=ref_Lx, ref_Ly=ref_Ly,
        ref_Tx=ref_Tx, ref_Ty=ref_Ty,
    )
    elapsed = time.time() - t0
    results.append({
        "trial": trial,
        "beta_c": beta_c,
        "chi2": chi2,
        "ndof": ndof,
        "chi2_ndof": chi2_ndof,
        "elapsed_s": round(elapsed, 1),
    })
    print(f"  [{trial+1:3d}/{N_TRIALS}]  beta_c={beta_c:.8f}  "
          f"chi2/ndof={chi2_ndof:.4e}  ({elapsed:.0f}s)")
    sys.stdout.flush()

    # Save incrementally
    out_path = os.path.join(OUTPUT_DIR, "stability_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

# Summary statistics
betas = np.array([r["beta_c"] for r in results])
costs = np.array([r["chi2_ndof"] for r in results])

print(f"\n{'='*60}")
print(f"Stability test complete: {N_TRIALS} trials")
print(f"  beta_c:    mean={betas.mean():.8f}  std={betas.std():.8f}  "
      f"min={betas.min():.8f}  max={betas.max():.8f}")
print(f"  chi2/ndof: mean={costs.mean():.4e}  std={costs.std():.4e}  "
      f"min={costs.min():.4e}  max={costs.max():.4e}")
print(f"  CoV(cost): {costs.std()/costs.mean():.2f}")
print(f"{'='*60}")

# Plot histograms
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

ax = axes[0]
ax.hist(betas, bins=20, color="steelblue", edgecolor="white")
ax.axvline(betas.mean(), color="red", ls="--", label=f"mean={betas.mean():.6f}")
ax.set_xlabel(r"$\beta_c$")
ax.set_ylabel("Count")
ax.set_title(r"$\beta_c$ distribution")
ax.legend(fontsize=8)

ax = axes[1]
ax.hist(costs, bins=20, color="darkorange", edgecolor="white")
ax.axvline(costs.mean(), color="red", ls="--", label=f"mean={costs.mean():.2e}")
ax.set_xlabel("chi2/ndof (boundary_slices)")
ax.set_ylabel("Count")
ax.set_title("Cost distribution")
ax.legend(fontsize=8)

ax = axes[2]
ax.hist(np.log10(costs), bins=20, color="seagreen", edgecolor="white")
ax.axvline(np.log10(costs.mean()), color="red", ls="--",
           label=f"mean={np.log10(costs.mean()):.2f}")
ax.set_xlabel("log10(chi2/ndof)")
ax.set_ylabel("Count")
ax.set_title("Cost distribution (log scale)")
ax.legend(fontsize=8)

fig.suptitle(f"Stability: {N_TRIALS} evals at (r1={r1}, r2={r2}), "
             f"{Lx}x{Ly}, {N_TRAJ_PROD} traj", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stability_histograms.png"), dpi=150)
print(f"\nPlot saved: {OUTPUT_DIR}/stability_histograms.png")
