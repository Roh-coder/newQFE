#!/usr/bin/env python3
"""
test_interpolation.py — Compare P1 linear FEM vs 6-parameter quadratic Taylor
approximation to the β_c surface.

Fits β_c(r₁,r₂) ≈ a₀ + a₁r₁ + a₂r₂ + a₃r₁² + a₄r₂² + a₅r₁r₂
to the existing 25 grid points, then generates 25 interstitial test points,
runs actual MC β_c scans at each, and compares truth vs predictions.
"""

import json
import os
import sys
import time

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import mc_engine as mc

# ── Configuration ─────────────────────────────────────────────────
SURFACE = "results/quick_test/betac_surfaces/surface_8x8_T0x0.json"
EXE = "bin/ising_tri_twisted_parallelogram"
Lx, Ly = 8, 8
N_TRAJ_COARSE = 5000
N_TRAJ_FINE = 10000
N_TEST = 25

# ── Load existing surface ─────────────────────────────────────────
with open(SURFACE) as f:
    raw = json.load(f)

r1_data = np.array([p["r1"] for p in raw["points"]])
r2_data = np.array([p["r2"] for p in raw["points"]])
bc_data = np.array([p["beta_c"] for p in raw["points"]])

print(f"Loaded {len(bc_data)} surface points from {SURFACE}")
print(f"  r1 ∈ [{r1_data.min():.2f}, {r1_data.max():.2f}]")
print(f"  r2 ∈ [{r2_data.min():.2f}, {r2_data.max():.2f}]")

# ── Fit 1: P1 linear FEM (Delaunay + piecewise linear) ───────────
pts = np.column_stack([r1_data, r2_data])
tri = Delaunay(pts)
interp_linear = LinearNDInterpolator(tri, bc_data)

# ── Fit 2: 6-parameter quadratic Taylor ──────────────────────────
#    β_c = a0 + a1*r1 + a2*r2 + a3*r1² + a4*r2² + a5*r1*r2
def design_matrix(r1, r2):
    return np.column_stack([
        np.ones_like(r1), r1, r2, r1**2, r2**2, r1 * r2
    ])

A = design_matrix(r1_data, r2_data)
coeffs, residuals, rank, sv = np.linalg.lstsq(A, bc_data, rcond=None)
a0, a1, a2, a3, a4, a5 = coeffs

print(f"\nQuadratic fit: β_c = {a0:.6f} + {a1:.6f}·r₁ + {a2:.6f}·r₂"
      f" + {a3:.6f}·r₁² + {a4:.6f}·r₂² + {a5:.6f}·r₁r₂")

# Residuals on training data
bc_quad_train = A @ coeffs
bc_lin_train = np.array([float(interp_linear(r1_data[i], r2_data[i]))
                         for i in range(len(r1_data))])
print(f"  Training RMS (linear):    {np.sqrt(np.mean((bc_lin_train - bc_data)**2)):.2e}")
print(f"  Training RMS (quadratic): {np.sqrt(np.mean((bc_quad_train - bc_data)**2)):.2e}")

def predict_quad(r1, r2):
    return a0 + a1*r1 + a2*r2 + a3*r1**2 + a4*r2**2 + a5*r1*r2

# ── Generate 25 interstitial test points ──────────────────────────
# Place them at centres of the 4×4 cells formed by the 5×5 grid,
# plus 9 extra random interior points to fill out to 25.
r1_vals = np.array(raw["r1_vals"])
r2_vals = np.array(raw["r2_vals"])

test_r1, test_r2 = [], []
# 16 cell centres (midpoints of the 4×4 grid cells)
for i in range(len(r1_vals) - 1):
    for j in range(len(r2_vals) - 1):
        test_r1.append(0.5 * (r1_vals[i] + r1_vals[i+1]))
        test_r2.append(0.5 * (r2_vals[j] + r2_vals[j+1]))

# 9 more at random interior positions (reproducible seed)
rng = np.random.RandomState(42)
for _ in range(N_TEST - 16):
    test_r1.append(rng.uniform(r1_vals[0] + 0.05, r1_vals[-1] - 0.05))
    test_r2.append(rng.uniform(r2_vals[0] + 0.05, r2_vals[-1] - 0.05))

test_r1 = np.array(test_r1[:N_TEST])
test_r2 = np.array(test_r2[:N_TEST])

print(f"\n{'='*70}")
print(f"Running MC β_c scans at {N_TEST} interstitial test points")
print(f"{'='*70}\n")

# ── Run MC scans at test points ───────────────────────────────────
mc_beta_c = np.zeros(N_TEST)
scan_dir = os.path.join("results", "quick_test", "interpolation_test")

for idx in range(N_TEST):
    r1, r2 = float(test_r1[idx]), float(test_r2[idx])
    k3 = 1.0
    k1, k2 = r1 * k3, r2 * k3

    # Use quadratic prediction as starting guess (it's cheap)
    beta_guess = predict_quad(r1, r2)
    margin = max(0.20 * beta_guess, 0.04)
    beta_lo = max(0.01, beta_guess - margin)
    beta_hi = beta_guess + margin

    pt_dir = os.path.join(scan_dir, f"pt_{idx:02d}")
    print(f"  [{idx+1:2d}/{N_TEST}] r1={r1:.4f} r2={r2:.4f}  "
          f"(guess={beta_guess:.6f})", end="", flush=True)
    t0 = time.time()
    try:
        bc, chi_peak, _, _, _ = mc.find_beta_c(
            EXE, Lx, Ly, 0, 0, k1, k2, k3,
            beta_lo, beta_hi,
            n_traj_coarse=N_TRAJ_COARSE, n_traj_fine=N_TRAJ_FINE,
            data_dir=pt_dir)
    except Exception as e:
        print(f"  FAILED: {e}")
        bc = np.nan
    elapsed = time.time() - t0
    mc_beta_c[idx] = bc
    print(f"  → β_c={bc:.8f} ({elapsed:.0f}s)")

# ── Predict with both methods ─────────────────────────────────────
pred_linear = np.array([float(interp_linear(test_r1[i], test_r2[i]))
                        for i in range(N_TEST)])
pred_quad = predict_quad(test_r1, test_r2)

# ── Comparison table ──────────────────────────────────────────────
err_lin = pred_linear - mc_beta_c
err_quad = pred_quad - mc_beta_c

print(f"\n{'='*90}")
print(f"{'#':>3}  {'r1':>6}  {'r2':>6}  {'MC β_c':>10}  "
      f"{'Linear':>10}  {'err_lin':>10}  {'Quadratic':>10}  {'err_quad':>10}  {'winner':>8}")
print(f"{'-'*90}")
lin_wins = 0
quad_wins = 0
for i in range(N_TEST):
    w = "QUAD" if abs(err_quad[i]) < abs(err_lin[i]) else "LIN"
    if abs(err_quad[i]) < abs(err_lin[i]):
        quad_wins += 1
    else:
        lin_wins += 1
    print(f"{i+1:3d}  {test_r1[i]:6.3f}  {test_r2[i]:6.3f}  "
          f"{mc_beta_c[i]:10.7f}  {pred_linear[i]:10.7f}  {err_lin[i]:+10.7f}  "
          f"{pred_quad[i]:10.7f}  {err_quad[i]:+10.7f}  {w:>8}")

mask = np.isfinite(mc_beta_c)
rms_lin = np.sqrt(np.mean(err_lin[mask]**2))
rms_quad = np.sqrt(np.mean(err_quad[mask]**2))
mae_lin = np.mean(np.abs(err_lin[mask]))
mae_quad = np.mean(np.abs(err_quad[mask]))
max_lin = np.max(np.abs(err_lin[mask]))
max_quad = np.max(np.abs(err_quad[mask]))

print(f"{'-'*90}")
print(f"\n  Summary ({N_TEST} test points):")
print(f"  {'':20s}  {'P1 Linear':>12s}  {'Quadratic':>12s}")
print(f"  {'RMS error':20s}  {rms_lin:12.7f}  {rms_quad:12.7f}")
print(f"  {'MAE':20s}  {mae_lin:12.7f}  {mae_quad:12.7f}")
print(f"  {'Max |error|':20s}  {max_lin:12.7f}  {max_quad:12.7f}")
print(f"  {'Wins':20s}  {lin_wins:12d}  {quad_wins:12d}")
print(f"\n  Quadratic RMS / Linear RMS = {rms_quad/rms_lin:.3f}")
if rms_quad < rms_lin:
    print(f"  → Quadratic is {(1 - rms_quad/rms_lin)*100:.1f}% better by RMS")
else:
    print(f"  → Linear is {(1 - rms_lin/rms_quad)*100:.1f}% better by RMS")

# ── Save results ──────────────────────────────────────────────────
results = {
    "quadratic_coeffs": {"a0": a0, "a1": a1, "a2": a2,
                         "a3": a3, "a4": a4, "a5": a5},
    "rms_linear": rms_lin, "rms_quadratic": rms_quad,
    "mae_linear": mae_lin, "mae_quadratic": mae_quad,
    "max_linear": max_lin, "max_quadratic": max_quad,
    "linear_wins": lin_wins, "quadratic_wins": quad_wins,
    "test_points": [
        {"r1": float(test_r1[i]), "r2": float(test_r2[i]),
         "mc_beta_c": float(mc_beta_c[i]),
         "pred_linear": float(pred_linear[i]),
         "pred_quadratic": float(pred_quad[i])}
        for i in range(N_TEST)
    ],
}
out_json = os.path.join("results", "quick_test", "interpolation_comparison.json")
with open(out_json, "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  Results saved: {out_json}")

# ── Visualization ─────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(17, 5))

# Panel 1: error scatter
ax = axes[0]
ax.scatter(np.abs(err_lin), np.abs(err_quad), c=mc_beta_c, cmap="viridis",
           s=50, edgecolors="k", linewidths=0.5)
lim = max(np.max(np.abs(err_lin[mask])), np.max(np.abs(err_quad[mask]))) * 1.15
ax.plot([0, lim], [0, lim], "k--", alpha=0.4, label="equal error")
ax.set_xlabel("|error|  P1 Linear")
ax.set_ylabel("|error|  Quadratic Taylor")
ax.set_title("Error comparison (below diagonal = quadratic wins)")
ax.legend(fontsize=9)
ax.set_aspect("equal")
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)

# Panel 2: residuals vs r1
ax = axes[1]
ax.scatter(test_r1, err_lin, marker="s", label="P1 Linear", alpha=0.7, s=30)
ax.scatter(test_r1, err_quad, marker="o", label="Quadratic", alpha=0.7, s=30)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("r₁")
ax.set_ylabel("Prediction − MC truth")
ax.set_title("Residuals vs r₁")
ax.legend(fontsize=9)

# Panel 3: residuals vs r2
ax = axes[2]
ax.scatter(test_r2, err_lin, marker="s", label="P1 Linear", alpha=0.7, s=30)
ax.scatter(test_r2, err_quad, marker="o", label="Quadratic", alpha=0.7, s=30)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("r₂")
ax.set_ylabel("Prediction − MC truth")
ax.set_title("Residuals vs r₂")
ax.legend(fontsize=9)

fig.suptitle(f"P1 Linear vs Quadratic Taylor — {N_TEST} interstitial test points  "
             f"(RMS: lin={rms_lin:.5f}, quad={rms_quad:.5f})", fontsize=12, y=1.02)
fig.tight_layout()
out_png = os.path.join("results", "quick_test", "interpolation_comparison.png")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"  Plot saved: {out_png}")
plt.close(fig)

print("\nDone.")
