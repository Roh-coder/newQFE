#!/usr/bin/env python3
"""
plot_fss_scaling_function_456.py

PROPER continuum comparison via the universal CFT scaling function.

The two-point function on a torus at criticality satisfies

    G(x; L) = L^{-2Δ} F(x/L)            as L → ∞,

where F is the *universal scaling function* and 2Δ_σ = 1/4 for 2D Ising.
Therefore the right FSS observable to extrapolate is

    Y(t; L) ≡ G(t·L_phys; L) · L_phys^{1/4}     →     F(t)    as L → ∞.

Y has a finite continuum limit independent of the lattice realization
(up to the operator normalization Z_σ²).  Comparing Y_∞^test(t) to
Y_∞^ref(t)  is the right way to test whether the test and ref lattices
realise the same continuum geometry.

PROCEDURE
  1. For each (cycle c, position t_k):
       Y(t_k; L) = G(t_k; L) · L_phys^{1/4}
       Fit  Y = A * exp(b / L_phys)  → continuum amplitude  A_F^test(c, t_k),
                                                            A_F^ref(c, t_k).
  2. CFT amplitude per cycle:  fit  A_F(c, t_k) = A_c / sin(π t_k)^{1/4}
       → A_c^test, A_c^ref.
  3. GEOMETRY-INVARIANT RATIOS  R_c = A_c / A_0  (Z_σ cancels):
       compare R_c^test vs R_c^ref.
  4. Plot:
       (row 0) Y(t; L) vs 1/L_phys with FSS fits  per cycle (color = t_k)
       (row 1) shape check  Y_∞(t) · sin(πt)^{1/4}  flat at  A_c
       (row 2) ratio  R_c^test / R_c^ref  vs cycle-length ratio  |p_c|/|p_0|

Output: results/_fss_456/fss_scaling_function_456.png
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

# ── Reuse data-loading machinery from the main FSS script (suppress its plot)
import importlib.util, matplotlib.figure as _figmod
_orig_savefig = _figmod.Figure.savefig
_figmod.Figure.savefig = lambda *a, **k: None
_spec = importlib.util.spec_from_file_location(
    "_fss_main", os.path.join(_HERE, "plot_fss_continuum_limit_456.py"))
_mod = importlib.util.module_from_spec(_spec)
try:
    _spec.loader.exec_module(_mod)
finally:
    _figmod.Figure.savefig = _orig_savefig

T_SAMPLE      = _mod.T_SAMPLE
REF_BASE_LEN  = _mod._REF_BASE_LEN
REF_ALPHAS    = _mod.REF_ALPHAS
TEST_SIZES    = _mod.TEST_SIZES
REF_GEOM      = _mod.REF_GEOM
TEST_GEOM     = _mod.TEST_GEOM
ref_data      = _mod.ref_data
test_data     = _mod.test_data
_sample_at_t  = _mod._sample_at_t
N_FIT         = _mod.N_FIT
TWO_DELTA     = 0.25     # 2Δ_σ for 2D Ising

# ─────────────────────────────────────────────────────────────────────────────
# Sample G at all t_k for each (geom, dataset) and rescale to Y = G·L_phys^{2Δ}
# ─────────────────────────────────────────────────────────────────────────────
print("\nSampling Y(t_k; L) = G * L_phys^{1/4} ...")

# storage:  Y[c]['test'][L] = (LP_phys, Y, Ye) with Y = G·LP^{1/4}
Y = {c: {'test': {}, 'ref': {}} for c in range(3)}

for c in range(3):
    # TEST family
    for L in sorted(test_data):
        Lx, Ly, Tx, Ty = TEST_GEOM[L]
        # test cycle has bare length L (square torus), so L_phys = L
        LP = float(L)
        G_arr, Ge_arr = _sample_at_t(test_data[L], Lx, Ly, Tx, Ty, c, T_SAMPLE)
        Y_arr  = G_arr  * LP**TWO_DELTA
        Ye_arr = Ge_arr * LP**TWO_DELTA
        Y[c]['test'][L] = (LP, Y_arr, Ye_arr)

    # REF family
    for a in sorted(ref_data):
        Lx, Ly, Tx, Ty = REF_GEOM[a]
        LP = float(a) * REF_BASE_LEN[c]   # cycle-c physical length
        G_arr, Ge_arr = _sample_at_t(ref_data[a], Lx, Ly, Tx, Ty, c, T_SAMPLE)
        Y_arr  = G_arr  * LP**TWO_DELTA
        Ye_arr = Ge_arr * LP**TWO_DELTA
        Y[c]['ref'][a] = (LP, Y_arr, Ye_arr)

# ─────────────────────────────────────────────────────────────────────────────
# FSS fit  Y(t_k; L) = A_F(t_k) * exp(b / L_phys)  →  A_F(t_k) at L→∞
# ─────────────────────────────────────────────────────────────────────────────
def _fit_Y(LP_arr, Y_arr, Ye_arr):
    LP_arr = np.asarray(LP_arr, float)
    Y_arr  = np.asarray(Y_arr,  float)
    Ye_arr = np.asarray(Ye_arr, float)
    order  = np.argsort(LP_arr)[::-1][:N_FIT]
    LP, Y_, Ye_ = LP_arr[order], Y_arr[order], Ye_arr[order]
    n = len(LP)
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan
    def _model(L_, A, b):
        return A * np.exp(b / L_)
    A0 = float(Y_[0])
    b0 = float(np.log(max(Y_[-1], 1e-12)/max(A0, 1e-12)) * LP[-1])
    try:
        popt, pcov = curve_fit(_model, LP, Y_, p0=[A0, b0],
                               sigma=Ye_, absolute_sigma=True, maxfev=10_000)
        A, b = popt
        Ae   = float(np.sqrt(abs(pcov[0, 0])))
        chi2 = float(np.sum(((Y_ - _model(LP, A, b))/Ye_)**2))
        dof  = max(n - 2, 0)
        return float(A), Ae, float(b), chi2/dof if dof > 0 else np.nan
    except Exception:
        return float(np.mean(Y_)), float(np.std(Y_)), 0.0, np.nan

# Collect FSS-extrapolated A_F(c, t_k) for each family
A_F = {c: {'test': np.zeros(len(T_SAMPLE)),
           'ref':  np.zeros(len(T_SAMPLE))} for c in range(3)}
A_F_err = {c: {'test': np.zeros(len(T_SAMPLE)),
               'ref':  np.zeros(len(T_SAMPLE))} for c in range(3)}
b_F = {c: {'test': np.zeros(len(T_SAMPLE)),
           'ref':  np.zeros(len(T_SAMPLE))} for c in range(3)}

for c in range(3):
    for ki, tk in enumerate(T_SAMPLE):
        # TEST
        LP_arr  = np.array([Y[c]['test'][L][0]      for L in sorted(test_data)])
        Y_arr   = np.array([Y[c]['test'][L][1][ki]  for L in sorted(test_data)])
        Ye_arr  = np.array([Y[c]['test'][L][2][ki]  for L in sorted(test_data)])
        A, Ae, b, chi2 = _fit_Y(LP_arr, Y_arr, Ye_arr)
        A_F[c]['test'][ki]     = A
        A_F_err[c]['test'][ki] = Ae
        b_F[c]['test'][ki]     = b
        # REF
        LP_arr  = np.array([Y[c]['ref'][a][0]      for a in sorted(ref_data)])
        Y_arr   = np.array([Y[c]['ref'][a][1][ki]  for a in sorted(ref_data)])
        Ye_arr  = np.array([Y[c]['ref'][a][2][ki]  for a in sorted(ref_data)])
        A, Ae, b, chi2 = _fit_Y(LP_arr, Y_arr, Ye_arr)
        A_F[c]['ref'][ki]     = A
        A_F_err[c]['ref'][ki] = Ae
        b_F[c]['ref'][ki]     = b

# ─────────────────────────────────────────────────────────────────────────────
# CFT amplitude per cycle:   A_F(c, t_k) = A_c / sin(πt_k)^{1/4}
# ─────────────────────────────────────────────────────────────────────────────
def _cft(t, A):
    return A / np.sin(np.pi * t) ** TWO_DELTA

def fit_amplitude(A_arr, Ae_arr):
    fin = np.isfinite(A_arr) & np.isfinite(Ae_arr) & (Ae_arr > 0)
    if fin.sum() < 2:
        return np.nan, np.nan
    try:
        popt, pcov = curve_fit(_cft, T_SAMPLE[fin], A_arr[fin],
                               p0=[float(np.nanmean(A_arr[fin]))],
                               sigma=Ae_arr[fin], absolute_sigma=True)
        return float(popt[0]), float(np.sqrt(abs(pcov[0, 0])))
    except Exception:
        return float(np.nanmean(A_arr[fin])), float(np.nanstd(A_arr[fin]))

A_cyc = {c: {} for c in range(3)};  Ae_cyc = {c: {} for c in range(3)}
for c in range(3):
    for fam in ('test', 'ref'):
        A_cyc[c][fam], Ae_cyc[c][fam] = fit_amplitude(A_F[c][fam], A_F_err[c][fam])

print("\n=== CFT scaling-function amplitudes  A_c  (universal continuum) ===")
print("  cycle |  A_c^test         |  A_c^ref          | A_test/A_ref")
for c in range(3):
    print(f"    {c}   |  {A_cyc[c]['test']:.4f}±{Ae_cyc[c]['test']:.4f} |  "
          f"{A_cyc[c]['ref']:.4f}±{Ae_cyc[c]['ref']:.4f} |  "
          f"{A_cyc[c]['test']/A_cyc[c]['ref']:.4f}")

# Geometry-invariant ratios
R_test = {c: A_cyc[c]['test']/A_cyc[0]['test'] for c in range(3)}
R_ref  = {c: A_cyc[c]['ref'] /A_cyc[0]['ref']  for c in range(3)}
LRatio = {c: REF_BASE_LEN[c]/REF_BASE_LEN[0]   for c in range(3)}

print("\n=== Geometry-invariant ratios  A_c / A_0  (Z_σ-independent) ===")
print("  cycle |  R_test    |  R_ref     |  R_test/R_ref  |  |p_c|/|p_0|")
for c in range(3):
    print(f"    {c}   |  {R_test[c]:.4f}   |  {R_ref[c]:.4f}   |  "
          f"{R_test[c]/R_ref[c]:.4f}        |  {LRatio[c]:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure: 3 rows × 3 columns
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 3, figsize=(19, 14))
cmap   = plt.get_cmap('plasma')
tcols  = [cmap(0.10 + 0.85 * (k+1)/(len(T_SAMPLE)+1)) for k in range(len(T_SAMPLE))]

CYCLE_PTITLES = [
    rf"Cycle 0   $|p^{{\alpha=1}}|={REF_BASE_LEN[0]:.2f}$",
    rf"Cycle 1   $|p^{{\alpha=1}}|={REF_BASE_LEN[1]:.2f}$",
    rf"Cycle 2   $|p^{{\alpha=1}}|={REF_BASE_LEN[2]:.2f}$",
]

# ── ROW 0: FSS of Y(t_k; L) vs 1/L_phys per cycle
xd = np.linspace(0.001, 0.14, 200)
for c in range(3):
    ax = axes[0][c]
    for ki in range(len(T_SAMPLE)):
        col = tcols[ki]
        # test
        LP = np.array([Y[c]['test'][L][0]     for L in sorted(test_data)])
        Yk = np.array([Y[c]['test'][L][1][ki] for L in sorted(test_data)])
        Ye = np.array([Y[c]['test'][L][2][ki] for L in sorted(test_data)])
        ax.errorbar(1/LP, Yk, yerr=Ye, fmt='o', color=col, ms=6, lw=1.2,
                    capsize=3, alpha=0.9, zorder=3)
        # ref
        LP = np.array([Y[c]['ref'][a][0]      for a in sorted(ref_data)])
        Yk = np.array([Y[c]['ref'][a][1][ki]  for a in sorted(ref_data)])
        Ye = np.array([Y[c]['ref'][a][2][ki]  for a in sorted(ref_data)])
        ax.errorbar(1/LP, Yk, yerr=Ye, fmt='s', color=col, ms=6, lw=1.2,
                    capsize=3, alpha=0.9, mfc='white', mew=1.4, zorder=3)
        # fit curves
        At, bt = A_F[c]['test'][ki], b_F[c]['test'][ki]
        Ar, br = A_F[c]['ref'][ki],  b_F[c]['ref'][ki]
        if np.isfinite(At):
            yd_t = At * np.exp(bt * xd)
            ax.plot(xd, yd_t, color=col, lw=0.9, alpha=0.55, ls='-')
            ax.plot(0, At, marker='|', ms=10, color=col, mew=2)
        if np.isfinite(Ar):
            yd_r = Ar * np.exp(br * xd)
            ax.plot(xd, yd_r, color=col, lw=0.9, alpha=0.55, ls='--')
            ax.plot(0, Ar, marker='|', ms=10, color=col, mew=2)
    ax.axvline(0, color='k', lw=0.7, ls=':', alpha=0.5)
    ax.set_xlim(-0.005, 0.135);  ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$1/L_{\rm phys}$", fontsize=11)
    ax.set_ylabel(r"$Y(t_k; L) = G \cdot L_{\rm phys}^{1/4}$", fontsize=11)
    ax.set_title(CYCLE_PTITLES[c]
                 + r"  — universal scaling FSS"
                 + f"\n(circle=test, open sq=ref; N_FIT={N_FIT})",
                 fontsize=10)

# ── ROW 1: shape check  A_F(t_k) · sin(πt_k)^{1/4}  flat at A_c
sf = np.sin(np.pi * T_SAMPLE) ** TWO_DELTA
for c in range(3):
    ax = axes[1][c]
    yt = A_F[c]['test'] * sf;  yte = A_F_err[c]['test'] * sf
    yr = A_F[c]['ref']  * sf;  yre = A_F_err[c]['ref']  * sf
    ax.errorbar(T_SAMPLE, yt, yerr=yte, fmt='o', color='#d94801', ms=10,
                lw=1.5, capsize=4, label=rf"test  $A_{{{c}}}={A_cyc[c]['test']:.4f}$",
                zorder=4)
    ax.axhline(A_cyc[c]['test'], color='#d94801', lw=1.4, ls='-', alpha=0.7)
    ax.errorbar(T_SAMPLE, yr, yerr=yre, fmt='s', color='#08306b', ms=10,
                lw=1.5, capsize=4, mfc='white', mew=1.6,
                label=rf"ref   $A_{{{c}}}={A_cyc[c]['ref']:.4f}$", zorder=4)
    ax.axhline(A_cyc[c]['ref'], color='#08306b', lw=1.4, ls='--', alpha=0.7)
    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$t$", fontsize=11)
    ax.set_ylabel(r"$A_F(t)\,\sin(\pi t)^{1/4}$  (flat $\Rightarrow$ CFT)", fontsize=10)
    ax.set_title(CYCLE_PTITLES[c] + "\nUniversal-amplitude shape check",
                 fontsize=10)
    ax.legend(fontsize=9, loc='lower center')
    ax.grid(alpha=0.3)

# ── ROW 2: ratio comparison and shape overlay
# panel 0:   A_c^test vs A_c^ref  bar comparison
ax = axes[2][0]
cs = np.arange(3)
ax.bar(cs - 0.18, [A_cyc[c]['test'] for c in cs], width=0.34,
       yerr=[Ae_cyc[c]['test'] for c in cs], capsize=4,
       color='#d94801', alpha=0.85, label='test')
ax.bar(cs + 0.18, [A_cyc[c]['ref']  for c in cs], width=0.34,
       yerr=[Ae_cyc[c]['ref']  for c in cs], capsize=4,
       color='#08306b', alpha=0.85, label='ref',
       edgecolor='k')
ax.set_xticks(cs);  ax.set_xticklabels([f"c={c}" for c in cs])
ax.set_ylabel(r"Universal CFT amplitude  $A_c$", fontsize=11)
ax.set_title("Per-cycle universal amplitudes\n"
             r"(differ by $Z_\sigma^2$ between families)",
             fontsize=10)
ax.legend(fontsize=9);  ax.grid(alpha=0.3, axis='y')

# panel 1: geometry-invariant ratios A_c/A_0
ax = axes[2][1]
xs = cs
ax.errorbar(xs - 0.06, [R_test[c] for c in cs], fmt='o',
            color='#d94801', ms=12, lw=2, label=r"$R^{\rm test}_c=A_c/A_0$",
            zorder=5)
ax.errorbar(xs + 0.06, [R_ref[c]  for c in cs], fmt='s',
            color='#08306b', ms=12, lw=2, mfc='white', mew=1.7,
            label=r"$R^{\rm ref}_c=A_c/A_0$", zorder=5)
ax.scatter(cs, [LRatio[c] for c in cs], marker='*', s=260, c='gray',
           label=r"$|p_c|/|p_0|$ (geom)", zorder=4)
ax.axhline(1.0, color='k', lw=0.7, ls=':', alpha=0.5)
ax.set_xticks(cs);  ax.set_xticklabels([f"c={c}" for c in cs])
ax.set_ylabel(r"Amplitude ratio $A_c/A_0$", fontsize=11)
ax.set_title("Geometry-invariant ratios\n"
             r"agreement $\Rightarrow$ same continuum geometry",
             fontsize=10)
ax.legend(fontsize=9, loc='best');  ax.grid(alpha=0.3)

# panel 2:  R_test/R_ref vs cycle-length ratio
ax = axes[2][2]
ratio_of_ratio = [R_test[c]/R_ref[c] for c in cs]
ax.errorbar(cs, ratio_of_ratio, fmt='D', color='#54278f', ms=14, lw=2,
            label="data:  $R_c^{\\rm test}/R_c^{\\rm ref}$", zorder=5)
ax.scatter(cs, [LRatio[c] for c in cs], marker='*', s=260, c='gray',
           label=r"$|p_c|/|p_0|$  (cycle-length geom)", zorder=4)
ax.axhline(1.0, color='g', lw=1.2, ls='--', alpha=0.7,
           label="perfect agreement = 1")
for c in cs:
    ax.text(c, ratio_of_ratio[c], f"  {ratio_of_ratio[c]:.3f}",
            fontsize=10, va='center')
ax.set_xticks(cs);  ax.set_xticklabels([f"c={c}" for c in cs])
ax.set_ylabel(r"$R_c^{\rm test}/R_c^{\rm ref}$", fontsize=11)
ax.set_title("Test/Ref ratio of amplitude ratios\n"
             "(ideal: equals 1 if same continuum geometry)",
             fontsize=10)
ax.legend(fontsize=9, loc='best');  ax.grid(alpha=0.3)

fig.suptitle(
    "Universal CFT scaling-function comparison — 4-5-6 geometry"
    "\n"
    r"FSS observable: $Y(t;L) = G(t;L)\cdot L_{\rm phys}^{1/4}$,  "
    r"continuum limit $Y \to F(t) = A/\sin(\pi t)^{1/4}$.",
    fontsize=12, y=1.01)

fig.tight_layout()
out_path = os.path.join(_mod._OUT, "fss_scaling_function_456.png")
fig.savefig(out_path, dpi=160, bbox_inches='tight')
print(f"\nSaved: {out_path}")
