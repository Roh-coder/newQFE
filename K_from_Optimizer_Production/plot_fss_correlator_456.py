#!/usr/bin/env python3
"""
plot_fss_correlator_456.py — FSS two-point correlator approach to continuum
for the 4-5-6 triangle geometry.

Three panels, one per boundary cycle (u=cycle0, v=cycle1, w=cycle2).

On each panel two FSS series are overlaid:
  - TWISTED REFERENCE (4-5-6 torus):   α=1 (13,16,-3,3)
                                         α=2 (26,32,-6,6)
                                         α=3 (39,48,-9,9)
  - UNTWISTED TEST at exact truth couplings (r1,r2,k3=1):
                                         L=8, 16, 24, 32

x-axis: fractional cycle position  t = k / N_SAMP  k=1..N_SAMP-1
y-axis: G_conn(t) along the cycle (evaluated via tiled linear interpolation)

Large markers: N_SAMP=8 sample positions t_k = k/8 (k=1..7).
Vertical dotted lines at t_k mark where every curve is evaluated identically.
CFT: G(t) = A / sin(π t)^{1/4}  (2D Ising Δ_σ = 1/8).

Data from:  results/_fss_456/{ref/a{α},test/L{L}}/two_point_all_to_all.dat
Output:     results/_fss_456/fss_correlator_456.png
"""
from __future__ import annotations

import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import boundary_paths, _SQRT3_2, _tile_interp

# ──────────────────────────────────────────────────────────────────────────────
# Paths / truth couplings
# ──────────────────────────────────────────────────────────────────────────────
_RES  = os.path.join(_HERE, "results", "_fss_456")
_OUT  = _RES
os.makedirs(_OUT, exist_ok=True)

_R1 = (2*math.log(5) - math.log(7)) / (2*math.log(3) - math.log(7))   # ≈5.0652
_R2 = math.log(7) / (2*math.log(3) - math.log(7))                      # ≈7.7429

CYCLE_LABELS = ["cycle 0  (u)", "cycle 1  (v)", "cycle 2  (w)"]

# ──────────────────────────────────────────────────────────────────────────────
# Helper: sample one boundary cycle from an all-to-all dataset
# ──────────────────────────────────────────────────────────────────────────────
N_SAMP = 8
T_SAMPLE = np.array([k / N_SAMP for k in range(1, N_SAMP)])  # [1/8 … 7/8]

def _sample_cycle(dat, Lx, Ly, Tx, Ty, cycle_idx, n_t=200):
    """
    Returns:
      (t_dense, G_dense, Ge_dense)  — full curve for plotting (n_t points)
      (T_SAMPLE, G_samp, Ge_samp)   — values at the 7 shared sample positions
    """
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    dm, dn = paths[cycle_idx]
    ex = dm + 0.5 * dn
    ey = _SQRT3_2 * dn

    iG  = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn",     copies=2)
    iGe = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn_err", copies=2)

    def _eval(t_arr):
        xy = np.column_stack([t_arr * ex, t_arr * ey])
        return iG(xy).ravel(), iGe(xy).ravel()

    t_d = np.linspace(1 / (n_t + 1), 1 - 1 / (n_t + 1), n_t)
    G_d, Ge_d = _eval(t_d)
    G_s, Ge_s = _eval(T_SAMPLE)

    return t_d, G_d, Ge_d, T_SAMPLE, G_s, Ge_s


# ──────────────────────────────────────────────────────────────────────────────
# Load datasets
# ──────────────────────────────────────────────────────────────────────────────

def _try_load(path):
    """Load all-to-all dat, or return None with a warning if missing."""
    if not os.path.exists(path):
        print(f"  WARNING: missing {path}")
        return None
    return mc_engine.load_all_to_all(path)

# --- twisted reference series (isotropic, 4-5-6 torus at α=1,2,3) ---
REF_ALPHAS = [1, 2, 3]
REF_GEOM   = {a: (13*a, 16*a, -3*a, 3*a) for a in REF_ALPHAS}
REF_COLORS  = ["#aec6e8", "#4292c6", "#08306b"]   # light→dark blue
REF_MARKERS = ["o", "s", "D"]

ref_data = {}
for a in REF_ALPHAS:
    path = os.path.join(_RES, "ref", f"a{a}", "two_point_all_to_all.dat")
    dat = _try_load(path)
    if dat is not None:
        print(f"  ref α={a}  ({len(dat)} keys)")
        ref_data[a] = dat

# --- untwisted test at exact truth couplings ---
TEST_SIZES  = [8, 16, 24, 32]
TEST_GEOM   = {L: (L, L, 0, 0) for L in TEST_SIZES}
TEST_COLORS = ["#fdd0a2", "#fd8d3c", "#d94801", "#7f2704"]  # light→dark orange
TEST_MARKERS = ["^", "v", "<", ">"]

test_data = {}
for L in TEST_SIZES:
    path = os.path.join(_RES, "test", f"L{L}", "two_point_all_to_all.dat")
    dat = _try_load(path)
    if dat is not None:
        print(f"  test L={L}  ({len(dat)} keys)")
        test_data[L] = dat

# ──────────────────────────────────────────────────────────────────────────────
# CFT continuum curve
# ──────────────────────────────────────────────────────────────────────────────
def _cft(t, A):
    return A / np.sin(np.pi * t) ** 0.25

# ──────────────────────────────────────────────────────────────────────────────
# Cycle physical lengths
# ──────────────────────────────────────────────────────────────────────────────
_ref3_geom = REF_GEOM[3]
paths3 = boundary_paths(*_ref3_geom)
phys_lengths = []
for dm, dn in paths3:
    ex = dm + 0.5 * dn; ey = _SQRT3_2 * dn
    phys_lengths.append(math.sqrt(ex**2 + ey**2))

# ──────────────────────────────────────────────────────────────────────────────
# Figure: 3 panels (one per cycle), horizontal
# ──────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), sharey=False)

for ci, ax in enumerate(axes):
    # ── fit CFT amplitude to largest available data per cycle ─────────────
    # Use α=3 ref if available, else α=2, else α=1
    A_ref_fit = np.nan
    for a_big in [3, 2, 1]:
        if a_big in ref_data:
            try:
                iGb = _tile_interp(ref_data[a_big], *REF_GEOM[a_big], "conn", 2)
                dm_b, dn_b = boundary_paths(*REF_GEOM[a_big])[ci]
                ex_b = dm_b + 0.5 * dn_b; ey_b = _SQRT3_2 * dn_b
                t_mid = np.linspace(0.20, 0.80, 100)
                G_mid = iGb(np.column_stack([t_mid * ex_b, t_mid * ey_b])).ravel()
                popt, _ = curve_fit(_cft, t_mid, G_mid, p0=[G_mid.mean()])
                A_ref_fit = float(popt[0])
                print(f"  Cycle {ci} CFT A (α={a_big}): {A_ref_fit:.4f}")
                break
            except Exception as e:
                print(f"  cycle {ci} fit failed (α={a_big}): {e}")

    t_cont = np.linspace(0.015, 0.985, 800)
    if np.isfinite(A_ref_fit):
        ax.plot(t_cont, _cft(t_cont, A_ref_fit),
                color="k", lw=2.2, ls=":", zorder=1,
                label=r"CFT  $A/\sin(\pi t)^{1/4}$" + f"  (A={A_ref_fit:.3f})")

    # vertical guide lines
    for ts in T_SAMPLE:
        ax.axvline(ts, color="gray", lw=0.6, ls=":", alpha=0.45, zorder=0)

    # ── twisted reference series ───────────────────────────────────────────
    for a, color, marker in zip(REF_ALPHAS, REF_COLORS, REF_MARKERS):
        if a not in ref_data:
            continue
        Lx, Ly, Tx, Ty = REF_GEOM[a]
        t_d, G_d, Ge_d, t_s, G_s, Ge_s = _sample_cycle(
            ref_data[a], Lx, Ly, Tx, Ty, ci)
        label = rf"Ref $\alpha={a}$  ({Lx}$\times${Ly}, iso)"
        ax.plot(t_d, G_d, lw=1.1, color=color, alpha=0.55, zorder=2)
        ax.fill_between(t_d, G_d - Ge_d, G_d + Ge_d, alpha=0.12, color=color)
        ax.errorbar(t_s, G_s, yerr=Ge_s,
                    fmt=marker, color=color, ms=9, lw=1.5, capsize=3.5,
                    zorder=5, label=label,
                    path_effects=[pe.withStroke(linewidth=2.8, foreground="white")])

    # ── untwisted test series ──────────────────────────────────────────────
    for L, color, marker in zip(TEST_SIZES, TEST_COLORS, TEST_MARKERS):
        if L not in test_data:
            continue
        Lx, Ly, Tx, Ty = TEST_GEOM[L]
        t_d, G_d, Ge_d, t_s, G_s, Ge_s = _sample_cycle(
            test_data[L], Lx, Ly, Tx, Ty, ci)
        label = rf"Test $L={L}$  (truth, sq)"
        ax.plot(t_d, G_d, lw=1.1, color=color, alpha=0.55, ls="--", zorder=3)
        ax.fill_between(t_d, G_d - Ge_d, G_d + Ge_d, alpha=0.12, color=color)
        ax.errorbar(t_s, G_s, yerr=Ge_s,
                    fmt=marker, color=color, ms=8, lw=1.5, capsize=3.5,
                    zorder=5, label=label,
                    path_effects=[pe.withStroke(linewidth=2.8, foreground="white")])

    # ── panel decoration ───────────────────────────────────────────────────
    ax.set_xlabel(r"$t = $ fractional cycle position", fontsize=11)
    if ci == 0:
        ax.set_ylabel(r"$G_{\mathrm{conn}}(t)$", fontsize=12)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=10)

    dm_r3, dn_r3 = paths3[ci]
    ax.set_title(
        f"$\\mathbf{{Cycle\\ {ci}}}$   "
        + CYCLE_LABELS[ci].split()[1]
        + f"\n$|p_{{\\alpha=3}}| \\approx {phys_lengths[ci]:.1f}$  "
        + f"  (period $({dm_r3},{dn_r3})$ at $\\alpha=3$)",
        fontsize=10,
    )

    ax.legend(loc="upper center", fontsize=7.5, framealpha=0.90,
              ncol=2, bbox_to_anchor=(0.5, 1.0))

# ── figure-level title ────────────────────────────────────────────────────────
fig.suptitle(
    "FSS two-point correlator per boundary cycle — 4-5-6 geometry   "
    rf"(truth: $r_1\approx{_R1:.3f}$, $r_2\approx{_R2:.3f}$, $\beta_c^\infty\approx0.0628$)"
    "\n"
    r"Blue: twisted iso ref  $\alpha=1,2,3$.   "
    r"Orange: untwisted truth test  $L=8,16,24,32$.   "
    r"Shape converges to CFT $1/\sin(\pi t)^{1/4}$ in both families.",
    fontsize=10.5, y=1.01,
)

fig.tight_layout()
out_path = os.path.join(_OUT, "fss_correlator_456.png")
fig.savefig(out_path, dpi=160, bbox_inches="tight")
print(f"\nSaved: {out_path}")
