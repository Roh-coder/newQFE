#!/usr/bin/env python3
"""
plot_surface_3d.py — 3D visualisation of the β_c critical surface.

Panel 1: Raw MC data points as scatter in (r₁, r₂, β_c) space.
Panel 2: Interpolated FEM manifold over the same domain.
"""

import json
import sys
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

# ── Load surface data ─────────────────────────────────────────────
SURFACE = os.path.join(os.path.dirname(__file__),
                       "results", "quick_test", "betac_surfaces",
                       "surface_8x8_T0x0.json")
if len(sys.argv) > 1:
    SURFACE = sys.argv[1]

with open(SURFACE) as f:
    raw = json.load(f)

r1 = np.array([p["r1"] for p in raw["points"]])
r2 = np.array([p["r2"] for p in raw["points"]])
bc = np.array([p["beta_c"] for p in raw["points"]])

# ── Build interpolant ─────────────────────────────────────────────
pts = np.column_stack([r1, r2])
tri = Delaunay(pts)
interp = LinearNDInterpolator(tri, bc)

r1_fine = np.linspace(r1.min(), r1.max(), 80)
r2_fine = np.linspace(r2.min(), r2.max(), 80)
R1, R2 = np.meshgrid(r1_fine, r2_fine)
BC = interp(R1.ravel(), R2.ravel()).reshape(R1.shape)

# ── Figure ────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 7))

# --- Panel 1: scatter points ---
ax1 = fig.add_subplot(121, projection="3d")
sc = ax1.scatter(r1, r2, bc, c=bc, cmap="viridis", s=80,
                 edgecolors="k", linewidths=0.6, depthshade=True,
                 zorder=5)
# Draw Delaunay edges on the surface for visual context
for simplex in tri.simplices:
    idx = list(simplex) + [simplex[0]]
    ax1.plot(r1[idx], r2[idx], bc[idx], "k-", alpha=0.25, lw=0.5)
ax1.set_xlabel("r₁ = k₁/k₃", labelpad=8)
ax1.set_ylabel("r₂ = k₂/k₃", labelpad=8)
ax1.set_zlabel("β_c", labelpad=6)
ax1.set_title("MC scan data points", fontsize=13, pad=12)
ax1.view_init(elev=28, azim=-55)
fig.colorbar(sc, ax=ax1, shrink=0.55, pad=0.10, label="β_c")

# --- Panel 2: interpolated manifold ---
ax2 = fig.add_subplot(122, projection="3d")
surf = ax2.plot_surface(R1, R2, BC, cmap="viridis", alpha=0.85,
                        edgecolor="none", rstride=1, cstride=1,
                        antialiased=True)
# Overlay the data points
ax2.scatter(r1, r2, bc, c="red", s=40, edgecolors="k",
            linewidths=0.5, zorder=10, label="MC data")
ax2.set_xlabel("r₁ = k₁/k₃", labelpad=8)
ax2.set_ylabel("r₂ = k₂/k₃", labelpad=8)
ax2.set_zlabel("β_c", labelpad=6)
ax2.set_title("P1 FEM interpolated manifold", fontsize=13, pad=12)
ax2.view_init(elev=28, azim=-55)
ax2.legend(loc="upper right", fontsize=9)
fig.colorbar(surf, ax=ax2, shrink=0.55, pad=0.10, label="β_c")

Lx, Ly = raw["Lx"], raw["Ly"]
fig.suptitle(f"Critical surface  β_c(r₁, r₂)   —   {Lx}×{Ly} triangular Ising",
             fontsize=14, y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.94])

out = os.path.splitext(SURFACE)[0] + "_3d.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved: {out}")
plt.close(fig)
