#!/usr/bin/env python3
"""
plot_fss_correlator_111.py — FSS two-point correlator approach to continuum.

For the 1-1-1 equilateral geometry (Lx=Ly=L, Tx=Ty=0, r1=r2=k3=1 at β_c),
plot G_conn(t) = G_conn(m, 0) vs t = m/L along the e1=(1,0) direction for
five lattice sizes: L ∈ {8, 12, 16, 24, 32}.

On each curve, highlight the 7 "cost sample points" shared across all L:
  t_k = k/8  for  k = 1, 2, 3, 4, 5, 6, 7

These are the fractional positions used inside landscape_l2_ladder.py
(N_samp=8). At the truth (r1,r2)=(1,1) these sample points fall on
different lattice sites for each L:
  L=8  → m = 1,2,3,4,5,6,7        (every site)
  L=12 → m = 1.5,3,4.5,...        (interpolated via _tile_interp)
  L=16 → m = 2,4,6,8,10,12,14     (every-other site)
  L=24 → m = 3,6,9,12,15,18,21    (every-third)
  L=32 → m = 4,8,12,16,20,24,28   (every-fourth)

A CFT continuum curve G(t) ∝ [sin(πt)]^{-1/4} (2D Ising Δ_σ=1/8) is
overlaid after fitting the amplitude to the L=32 data.

Data sources (all at ~20k trajectories):
  Test L=8,12,16 : results/_ladder_111_line20k/test/grid/r1_1.000_r2_1.000_L{L}.pkl
  Ref  L=16,24,32: results/_ladder_111/ref/L{L}/two_point_all_to_all.dat

Output: results/_ladder_111_line20k/fss_correlator_111.png
"""
from __future__ import annotations

import os
import sys
import pickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────
_BASE20K = os.path.join(_HERE, "results", "_ladder_111_line20k")
_BASE2K  = os.path.join(_HERE, "results", "_ladder_111")
_OUT_DIR = _BASE20K
os.makedirs(_OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Load data
# ──────────────────────────────────────────────────────────────────────────────

def load_test(L: int) -> dict:
    """Load test PKL at truth (r1=1, r2=1) for given L from line20k."""
    path = os.path.join(_BASE20K, "test", "grid",
                        f"r1_1.000_r2_1.000_L{L}.pkl")
    with open(path, "rb") as f:
        return pickle.load(f)

def load_ref(L: int) -> dict:
    """Load reference all-to-all dat for given L from _ladder_111/ref."""
    path = os.path.join(_BASE2K, "ref", f"L{L}", "two_point_all_to_all.dat")
    return mc_engine.load_all_to_all(path)

def extract_e1_curve(data: dict, L: int):
    """
    Extract G_conn(m,0) along e1=(1,0) for m=1..L-1, using symmetry
    G(L-m) = G(m) to fold the full period onto t ∈ (0,1).
    Returns (t, G, G_err) sorted by t, covering the full period.
    """
    t_vals, G_vals, G_err_vals = [], [], []
    for m in range(1, L):
        # use the nearer half: m and L-m give the same G; keep both for display
        key = (m, 0)
        if key in data:
            entry = data[key]
            t_vals.append(m / L)
            G_vals.append(entry["conn"])
            G_err_vals.append(entry["conn_err"])
    order = np.argsort(t_vals)
    return (np.array(t_vals)[order],
            np.array(G_vals)[order],
            np.array(G_err_vals)[order])

# ──────────────────────────────────────────────────────────────────────────────
# Datasets: (label, L, t, G, G_err)
# ──────────────────────────────────────────────────────────────────────────────
datasets = []

# Use only L that are multiples of N_SAMP=8 so that t_k=k/8 hits integer
# lattice sites m = k*(L/8) exactly — no interpolation needed.
#
#   L= 8: m = 1,  2,  3,  4,  5,  6,  7
#   L=16: m = 2,  4,  6,  8,  10, 12, 14
#   L=24: m = 3,  6,  9,  12, 15, 18, 21
#   L=32: m = 4,  8,  12, 16, 20, 24, 28

# L=8 — from 20k test PKL
d8 = load_test(8)
t, G, Ge = extract_e1_curve(d8["test_data"], 8)
datasets.append(("L= 8  (test, 20k)", 8, t, G, Ge))

# L=16 — from 20k test PKL (isotropic truth; same as ref up to FSS)
d16 = load_test(16)
t, G, Ge = extract_e1_curve(d16["test_data"], 16)
datasets.append(("L=16  (test, 20k)", 16, t, G, Ge))

# L=24,32 — from _ladder_111/ref (20k traj each)
for L in [24, 32]:
    dat = load_ref(L)
    t, G, Ge = extract_e1_curve(dat, L)
    datasets.append((f"L={L}  (ref,  20k)", L, t, G, Ge))

# ──────────────────────────────────────────────────────────────────────────────
# Continuum fit:  G_cont(t) = A * (sin(π t))^{-1/4}
# Use L=32 midrange to fit amplitude
# ──────────────────────────────────────────────────────────────────────────────
t32, G32, Ge32 = datasets[-1][2], datasets[-1][3], datasets[-1][4]
# half-period only for the fit, midrange t to avoid short-distance artifacts
mask_fit = (t32 >= 0.15) & (t32 <= 0.50)

def _cfr_cont(t, A):
    return A / np.sin(np.pi * t) ** 0.25

popt, _ = curve_fit(_cfr_cont, t32[mask_fit], G32[mask_fit],
                    p0=[G32[mask_fit].mean()])
A_cont = float(popt[0])
print(f"CFT amplitude fit (L=32): A = {A_cont:.4f}")

# ──────────────────────────────────────────────────────────────────────────────
# Common cost sample-point t-values (N_samp = 8; full period, drop t=0,1)
# k = 1 .. N_samp-1  gives  t_k = k/8
# ──────────────────────────────────────────────────────────────────────────────
N_SAMP = 8
t_sample_full = np.array([k / N_SAMP for k in range(1, N_SAMP)])
# t_sample_full = [1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8]
# Corresponding lattice sites for each L:
#   L= 8: m = 1, 2, 3, 4, 5, 6, 7
#   L=16: m = 2, 4, 6, 8, 10, 12, 14
#   L=24: m = 3, 6, 9, 12, 15, 18, 21
#   L=32: m = 4, 8, 12, 16, 20, 24, 28

COLORS  = ["#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#984ea3"]
MARKERS = ["o", "s", "^", "D", "v"]

# ──────────────────────────────────────────────────────────────────────────────
# Figure: one clean axes, full period t ∈ (0, 1)
# ──────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 1, figsize=(11, 6.2))
ax = axes

# CFT continuum curve — full period
t_cont = np.linspace(0.018, 0.982, 900)
G_cont = _cfr_cont(t_cont, A_cont)
ax.plot(t_cont, G_cont, color="k", lw=2.2, ls="--", zorder=1,
        label=r"Continuum: $A\,/\sin(\pi t)^{1/4}$  (CFT, $\Delta_\sigma=1/8$, fit to $L=32$)")

# Vertical guide lines at each cost sample position
for ts in t_sample_full:
    ax.axvline(ts, color="gray", lw=0.7, ls=":", alpha=0.5, zorder=0)

for (label, L, t, G, Ge), color, marker in zip(datasets, COLORS, MARKERS):
    # thin background curve through every available lattice site
    ax.plot(t, G, lw=1.1, color=color, alpha=0.40, zorder=2)
    ax.fill_between(t, G - Ge, G + Ge, alpha=0.11, color=color)

    # highlighted sample points at t_k = k/8
    # Because L is a multiple of 8, m_k = k*(L//8) is exactly an integer
    # lattice site — no interpolation needed.
    G_samp, Ge_samp = [], []
    step = L // N_SAMP          # e.g. L=8→1, L=16→2, L=24→3, L=32→4
    for k in range(1, N_SAMP):
        m_k = k * step          # exact integer site
        idx = np.searchsorted(t, m_k / L - 1e-9)
        idx = min(idx, len(t) - 1)
        G_samp.append(float(G[idx]))
        Ge_samp.append(float(Ge[idx]))

    ax.errorbar(t_sample_full, np.array(G_samp), yerr=np.array(Ge_samp),
                fmt=marker, color=color, ms=10, lw=1.6, capsize=4,
                zorder=5, label=label,
                path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])

    # annotate each marker with its integer lattice site m
    step = L // N_SAMP
    for k, (ts, gs) in enumerate(zip(t_sample_full, G_samp), start=1):
        m_k = k * step
        ax.text(ts, gs + 0.012, f"$m\\!=\\!{m_k}$",
                ha="center", va="bottom", fontsize=6.5, color=color,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

# ── axis decoration ──────────────────────────────────────────────────────────
ax.set_xlabel(r"Fractional cycle position  $t = m / L$", fontsize=13)
ax.set_ylabel(r"$G_{\mathrm{conn}}(t)$", fontsize=13)
ax.set_title(
    r"FSS of two-point correlator along $e_1$-cycle — 1-1-1 equilateral geometry"
    "\n"
    r"$(L_x=L_y=L,\; T_x=T_y=0,\; r_1=r_2=k_3=1)$ at $\beta_c(L)$,  all 20 000 traj.  "
    r"Large markers = 7 shared cost-sample positions $t_k=k/8$",
    fontsize=11,
)

ax.set_xlim(-0.01, 1.01)
ylim_top = max(d[3].max() for d in datasets) * 1.06
ax.set_ylim(0, ylim_top)
ax.legend(loc="upper center", fontsize=9.5, framealpha=0.92, ncol=3,
          bbox_to_anchor=(0.5, 1.0))
ax.tick_params(labelsize=11)

# ── in-plot annotation: site-index table ────────────────────────────────────
# Shows which lattice site m corresponds to each t_k = k/8 for each L
table_lines = ["Lattice sites  $m$  at shared sample positions $t_k = k/8$:"]
for (_, L, *__) in datasets:
    sites_str = "  ".join(
        f"{int(round(k * L / N_SAMP)):>2d}" for k in range(1, N_SAMP)
    )
    table_lines.append(f"$L={L:2d}$:  m = " + sites_str.replace("  ", r"$,$  "))
table_lines.insert(1, r"$t_k$:       " + "  ".join(
    f"{k/N_SAMP:.3f}" for k in range(1, N_SAMP)
))

ax.text(0.50, 0.03, "\n".join(table_lines),
        transform=ax.transAxes, fontsize=7.8,
        va="bottom", ha="center",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.85, ec="gray"),
        family="monospace", zorder=10)

fig.tight_layout()
out_path = os.path.join(_OUT_DIR, "fss_correlator_111.png")
fig.savefig(out_path, dpi=160, bbox_inches="tight")
print(f"Saved: {out_path}")
