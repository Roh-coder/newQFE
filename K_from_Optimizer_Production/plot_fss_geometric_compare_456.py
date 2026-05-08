#!/usr/bin/env python3
"""
plot_fss_geometric_compare_456.py

Geometry-invariant comparison of test (untwisted, anisotropic at truth couplings)
versus ref (twisted iso) using the two-point function in the continuum limit.

PROCEDURE (eliminates Z_sigma operator-normalization ambiguity):

  1. For each cycle c and each fractional position t_k = k/8 (k=1..7), extract
     G_inf^test(c, t_k) and G_inf^ref(c, t_k) by an  A * exp(b / L)  fit
     in L_phys  (same as plot_fss_continuum_limit_456.py).
  2. Per-cycle CFT amplitude:  A_c^F = best fit of  A / sin(pi t)^{1/4}
     to the 7 G_inf^F(c, t_k) values  (F ∈ {test, ref}).
  3. SHAPE check:  plot  G_inf^F(t) * sin(pi t)^{1/4}  vs t for each cycle.
     Should be a flat line at value A_c^F  if the CFT prediction holds.
  4. GEOMETRIC RATIOS:  R_c^F = A_c^F / A_0^F.  These are independent of
     Z_sigma (cancels in the ratio).  If test and ref realise the same
     continuum geometry, R_c^test ≈ R_c^ref.

Output:  results/_fss_456/fss_geometric_compare_456.png
"""
from __future__ import annotations

import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# Reuse all heavy machinery (data loading, interpolants, FSS fit) from the
# main script — it imports cleanly because it's data-driven (no __main__ guard
# for plotting, but its module body runs fine; we suppress its plot here).
import matplotlib as _mpl
_mpl.rcParams['figure.max_open_warning'] = 0

# Import without re-running its plotting if possible.  Easiest: import as module.
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_fss_main", os.path.join(_HERE, "plot_fss_continuum_limit_456.py"))
_mod = importlib.util.module_from_spec(_spec)
# Suppress its figure save by stubbing savefig before exec
import matplotlib.figure as _figmod
_orig_savefig = _figmod.Figure.savefig
_figmod.Figure.savefig = lambda *a, **k: None
try:
    _spec.loader.exec_module(_mod)
finally:
    _figmod.Figure.savefig = _orig_savefig

T_SAMPLE      = _mod.T_SAMPLE
fits          = _mod.fits           # fits[c][ki] dict with Ginf_test, Ginf_ref
REF_BASE_LEN  = _mod._REF_BASE_LEN  # cycle base lengths

# ─────────────────────────────────────────────────────────────────────────────
# CFT amplitude per cycle per family
# ─────────────────────────────────────────────────────────────────────────────
def _cft(t, A):  # 2D Ising,  Δ_σ = 1/8  →  exponent 2Δ = 1/4
    return A / np.sin(np.pi * t) ** 0.25

def fit_amplitude(G, Ge):
    """Return (A, A_err) from CFT fit to G_inf(t_k)."""
    fin = np.isfinite(G) & np.isfinite(Ge) & (Ge > 0)
    if fin.sum() < 2:
        return np.nan, np.nan
    try:
        popt, pcov = curve_fit(_cft, T_SAMPLE[fin], G[fin],
                               p0=[float(np.nanmean(G[fin]))],
                               sigma=Ge[fin], absolute_sigma=True)
        return float(popt[0]), float(np.sqrt(abs(pcov[0, 0])))
    except Exception:
        return np.nan, np.nan

# Collect per-cycle G_inf and amplitude
Gt, Gte, Gr, Gre = {}, {}, {}, {}
A_test, Ae_test, A_ref, Ae_ref = {}, {}, {}, {}
for c in range(3):
    Gt[c]  = np.array([fits[c][ki]['Ginf_test']     for ki in range(len(T_SAMPLE))])
    Gte[c] = np.array([fits[c][ki]['Ginf_test_err'] for ki in range(len(T_SAMPLE))])
    Gr[c]  = np.array([fits[c][ki]['Ginf_ref']      for ki in range(len(T_SAMPLE))])
    Gre[c] = np.array([fits[c][ki]['Ginf_ref_err']  for ki in range(len(T_SAMPLE))])
    A_test[c],  Ae_test[c]  = fit_amplitude(Gt[c],  Gte[c])
    A_ref[c],   Ae_ref[c]   = fit_amplitude(Gr[c],  Gre[c])

print("\n=== Per-cycle CFT amplitudes (continuum) ===")
print(f"  cycle | A_test          | A_ref           | A_test/A_ref")
for c in range(3):
    print(f"    {c}   | {A_test[c]:.4f}±{Ae_test[c]:.4f}  | "
          f"{A_ref[c]:.4f}±{Ae_ref[c]:.4f}  | {A_test[c]/A_ref[c]:+.4f}")

# Geometry-invariant ratios A_c / A_0
R_test = {c: A_test[c] / A_test[0] for c in range(3)}
Re_test = {c: abs(R_test[c]) * math.sqrt((Ae_test[c]/A_test[c])**2 +
                                          (Ae_test[0]/A_test[0])**2)
           for c in range(3)}
R_ref  = {c: A_ref[c]  / A_ref[0]  for c in range(3)}
Re_ref = {c: abs(R_ref[c]) * math.sqrt((Ae_ref[c]/A_ref[c])**2 +
                                        (Ae_ref[0]/A_ref[0])**2)
          for c in range(3)}

print("\n=== Geometry-invariant ratios A_c / A_0  (Z_sigma cancels) ===")
print("  cycle |  R_test          |  R_ref           |  R_test/R_ref")
for c in range(3):
    print(f"    {c}   |  {R_test[c]:.4f}±{Re_test[c]:.4f}  |  "
          f"{R_ref[c]:.4f}±{Re_ref[c]:.4f}  |  {R_test[c]/R_ref[c]:.4f}")

# Length ratios for context
L0 = REF_BASE_LEN[0]
print("\n=== Cycle-length ratios |p_c| / |p_0|  (ref α=1 baseline) ===")
for c in range(3):
    print(f"    cycle {c}:  {REF_BASE_LEN[c]/L0:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure: 2 rows × 3 columns
#   Row 0 — shape check  G_inf(t)*sin(πt)^{1/4}  vs t   (should be flat at A_c)
#   Row 1 — geometry-invariant amplitude ratios A_c/A_0  test vs ref
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

t_dense = np.linspace(0.02, 0.98, 400)
shape_factor = np.sin(np.pi * t_dense) ** 0.25

CYCLE_PTITLES = [
    rf"Cycle 0   $|p^{{\alpha=1}}|={REF_BASE_LEN[0]:.2f}$",
    rf"Cycle 1   $|p^{{\alpha=1}}|={REF_BASE_LEN[1]:.2f}$",
    rf"Cycle 2   $|p^{{\alpha=1}}|={REF_BASE_LEN[2]:.2f}$",
]

# Row 0 — shape check
for c in range(3):
    ax = axes[0][c]
    sf = np.sin(np.pi * T_SAMPLE) ** 0.25

    # test points
    yt  = Gt[c]  * sf
    yte = Gte[c] * sf
    ax.errorbar(T_SAMPLE, yt, yerr=yte, fmt='o', color='#d94801', ms=9,
                lw=1.6, capsize=4, label=rf"test  $A_{{{c}}}={A_test[c]:.4f}$",
                zorder=4)
    ax.axhline(A_test[c], color='#d94801', lw=1.4, ls='-', alpha=0.7)

    # ref points
    yr  = Gr[c]  * sf
    yre = Gre[c] * sf
    ax.errorbar(T_SAMPLE, yr, yerr=yre, fmt='s', color='#08306b', ms=9,
                lw=1.6, capsize=4, mfc='white', mew=1.7,
                label=rf"ref   $A_{{{c}}}={A_ref[c]:.4f}$", zorder=4)
    ax.axhline(A_ref[c], color='#08306b', lw=1.4, ls='--', alpha=0.7)

    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$t$", fontsize=11)
    ax.set_ylabel(r"$G_\infty(t)\,\sin(\pi t)^{1/4}$  (should be flat)",
                  fontsize=10)
    ax.set_title(CYCLE_PTITLES[c]
                 + "\nShape check: flat ⇒ CFT model holds, plateau = $A_c$",
                 fontsize=10)
    ax.legend(fontsize=9, loc='lower center')
    ax.grid(alpha=0.3)

# Row 1 — geometry-invariant ratios
ax = axes[1][0]
cs = np.arange(3)
xs = cs - 0.10
xs2 = cs + 0.10
ax.errorbar(xs,  [R_test[c] for c in cs], yerr=[Re_test[c] for c in cs],
            fmt='o', color='#d94801', ms=12, lw=2, capsize=5,
            label="test  $A_c/A_0$", zorder=5)
ax.errorbar(xs2, [R_ref[c]  for c in cs], yerr=[Re_ref[c]  for c in cs],
            fmt='s', color='#08306b', ms=12, lw=2, capsize=5,
            mfc='white', mew=1.8, label="ref   $A_c/A_0$", zorder=5)
# also show the cycle-length ratios
ax.scatter(cs, [REF_BASE_LEN[c]/REF_BASE_LEN[0] for c in cs],
           marker='*', s=250, c='gray', zorder=3,
           label=r"$|p_c|/|p_0|$  (ref geom)")
ax.scatter(cs, [(REF_BASE_LEN[0]/REF_BASE_LEN[c]) for c in cs],
           marker='x', s=120, c='k', zorder=3, lw=2,
           label=r"$|p_0|/|p_c|$  (inverse)")
ax.axhline(1.0, color='k', lw=0.7, ls=':', alpha=0.5)
ax.set_xticks(cs);  ax.set_xticklabels([f"c={c}" for c in cs])
ax.set_ylabel(r"Amplitude ratio  $A_c / A_0$", fontsize=11)
ax.set_title("Geometry-invariant ratios "
             r"($Z_\sigma$ cancels)"
             "\nIf same continuum geometry: test ≈ ref",
             fontsize=10)
ax.legend(fontsize=9, loc='best')
ax.grid(alpha=0.3)

# Row 1, panel 1 — ratio of ratios test/ref
ax = axes[1][1]
ratio_of_ratio = [R_test[c]/R_ref[c] for c in cs]
err_ror = [abs(ratio_of_ratio[c]) *
           math.sqrt((Re_test[c]/R_test[c])**2 + (Re_ref[c]/R_ref[c])**2)
           for c in cs]
ax.errorbar(cs, ratio_of_ratio, yerr=err_ror,
            fmt='D', color='#54278f', ms=14, lw=2, capsize=5)
ax.axhline(1.0, color='g', lw=1.4, ls='--', alpha=0.8,
           label="agreement = 1")
ax.set_xticks(cs);  ax.set_xticklabels([f"c={c}" for c in cs])
ax.set_ylabel(r"$R_c^{\rm test} / R_c^{\rm ref}$", fontsize=11)
ax.set_title("Test/Ref of geometry-invariant ratios"
             "\n(Deviation from 1 = geometry mismatch)",
             fontsize=10)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
for c in cs:
    ax.text(c, ratio_of_ratio[c],
            f"  {ratio_of_ratio[c]:.3f}",
            fontsize=10, va='center')

# Row 1, panel 2 — joint CFT-residual scatter plot
ax = axes[1][2]
# overlay all G_inf*sin(πt)^{1/4} normalized by per-family-per-cycle amplitude
for c in range(3):
    sf = np.sin(np.pi * T_SAMPLE) ** 0.25
    if np.isfinite(A_test[c]):
        ax.errorbar(T_SAMPLE + 0.005*c, Gt[c]*sf/A_test[c],
                    yerr=Gte[c]*sf/A_test[c],
                    fmt='o', ms=7, lw=1.2, capsize=3,
                    color=plt.cm.Oranges(0.3 + 0.25*c),
                    label=f"test c={c}")
    if np.isfinite(A_ref[c]):
        ax.errorbar(T_SAMPLE - 0.005*c, Gr[c]*sf/A_ref[c],
                    yerr=Gre[c]*sf/A_ref[c],
                    fmt='s', ms=7, lw=1.2, capsize=3, mfc='white', mew=1.3,
                    color=plt.cm.Blues(0.4 + 0.2*c),
                    label=f"ref  c={c}")
ax.axhline(1.0, color='k', lw=1, ls='--', alpha=0.7)
ax.set_xlim(0, 1);  ax.set_ylim(0.85, 1.15)
ax.set_xlabel("t", fontsize=11)
ax.set_ylabel(r"$G_\infty(t)\sin(\pi t)^{1/4} / A_c$", fontsize=11)
ax.set_title("Normalized shape (all cycles, both families)\n"
             "Spread around 1 = non-CFT residuals",
             fontsize=10)
ax.legend(fontsize=7, loc='lower center', ncol=3)
ax.grid(alpha=0.3)

fig.suptitle(
    "Geometry-invariant comparison of two-point function — 4-5-6 geometry"
    "\n"
    r"Continuum amplitudes $A_c$ from CFT fit to $G_\infty(t_k)=A_c/\sin(\pi t_k)^{1/4}$;  "
    r"ratios $A_c/A_0$ are $Z_\sigma$-independent.",
    fontsize=12, y=1.02)

fig.tight_layout()
out_path = os.path.join(_mod._OUT, "fss_geometric_compare_456.png")
fig.savefig(out_path, dpi=160, bbox_inches='tight')
print(f"\nSaved: {out_path}")
