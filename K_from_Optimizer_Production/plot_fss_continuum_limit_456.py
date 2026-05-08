#!/usr/bin/env python3
"""
plot_fss_continuum_limit_456.py

Empirical FSS continuum-limit extraction for the 4-5-6 geometry.

At each of N_SAMP-1 = 7 sample positions t_k = k/8 along each boundary cycle,
fit  G(t_k; L_phys) = G_inf(t_k) + a(t_k) * L_phys^{-omega}  independently for:
  - UNTWISTED TEST  (L,L,0,0) at truth couplings, L ∈ {8,16,24,32,48,64}
  - TWISTED REF     (13α,16α,-3α,3α) iso, α ∈ {1,2,3,4}

omega is a FREE parameter (empirically determined).  Fit uses the N_FIT=4 largest
available sizes per family (one dof).  Smaller sizes are shown but not fitted.

Figure layout: 3 rows × 3 columns
  Row 0 — FSS fit (G vs 1/L_phys); circles=test, open squares=ref
  Row 1 — exponent b(t_k) per family per cycle
  Row 2 — continuum comparison G_inf test vs ref vs CFT A/sin(πt)^{1/4}

Output: results/_fss_456/fss_continuum_limit_456.png
"""
from __future__ import annotations

import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines  import Line2D
from scipy.optimize    import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine
from cost import boundary_paths, _SQRT3_2, _tile_interp

# ─────────────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────────────
_DATA = os.path.join(_HERE, "results", "_fss_456")
_OUT  = _DATA

N_SAMP   = 8
T_SAMPLE = np.array([k / N_SAMP for k in range(1, N_SAMP)])   # 1/8 … 7/8

N_FIT = 4   # use N_FIT largest available sizes per family for the nonlinear fit

REF_ALPHAS  = [1, 2, 3, 4]
REF_GEOM    = {a: (13*a, 16*a, -3*a, 3*a) for a in REF_ALPHAS}
TEST_SIZES  = [8, 16, 24, 32, 48, 64]
TEST_GEOM   = {L: (L, L, 0, 0) for L in TEST_SIZES}

# physical base lengths per cycle for the ref family (α=1)
_REF_BASE_LEN = {}   # computed below
for c in range(3):
    dm, dn = boundary_paths(13, 16, -3, 3)[c]
    ex = dm + 0.5*dn;  ey = _SQRT3_2 * dn
    _REF_BASE_LEN[c] = math.sqrt(ex**2 + ey**2)

# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────
def _load(path):
    if not os.path.exists(path):
        print(f"  missing: {path}"); return None
    return mc_engine.load_all_to_all(path)

ref_data  = {a: _load(os.path.join(_DATA,"ref",f"a{a}","two_point_all_to_all.dat"))
             for a in REF_ALPHAS}
test_data = {L: _load(os.path.join(_DATA,"test",f"L{L}","two_point_all_to_all.dat"))
             for L in TEST_SIZES}
ref_data  = {k:v for k,v in ref_data.items()  if v is not None}
test_data = {k:v for k,v in test_data.items() if v is not None}
print(f"Ref  alphas available : {sorted(ref_data.keys())}")
print(f"Test sizes  available : {sorted(test_data.keys())}")

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation helper: build interpolants once per (dataset, geom), cache them.
# Then evaluate at all t_k in one call.
# ─────────────────────────────────────────────────────────────────────────────
_interp_cache = {}   # (id(dat), Lx, Ly, Tx, Ty) -> (iG, iGe, paths)

def _get_interps(dat, Lx, Ly, Tx, Ty, copies=2):
    key = (id(dat), Lx, Ly, Tx, Ty)
    if key not in _interp_cache:
        iG  = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn",     copies)
        iGe = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn_err", copies)
        paths = boundary_paths(Lx, Ly, Tx, Ty)
        _interp_cache[key] = (iG, iGe, paths)
    return _interp_cache[key]

def _sample_at_t(dat, Lx, Ly, Tx, Ty, c, t_arr, copies=2):
    iG, iGe, paths = _get_interps(dat, Lx, Ly, Tx, Ty, copies)
    dm, dn = paths[c]
    ex = dm + 0.5*dn;  ey = _SQRT3_2 * dn
    t  = np.asarray(t_arr)
    xy = np.column_stack([t * ex, t * ey])
    return np.asarray(iG(xy), float), np.asarray(iGe(xy), float)

# ─────────────────────────────────────────────────────────────────────────────
# FSS fit:  G(L) = A * exp(b / L)
# Continuum limit L→∞: G→ A.  Two free params (A, b), WLS.
# Returns (A, A_err, 0.0, b, chi2_dof)  — A is G∞, b in slot 4 for middle panel.
# ─────────────────────────────────────────────────────────────────────────────

def _power_fit(L_arr, G_arr, Ge_arr):
    """Fit G = A * exp(b / L)  by WLS.  Returns (A, A_err, 0.0, b, chi2_dof)."""
    L_arr  = np.asarray(L_arr,  float)
    G_arr  = np.asarray(G_arr,  float)
    Ge_arr = np.asarray(Ge_arr, float)

    order = np.argsort(L_arr)[::-1][:N_FIT]
    L  = L_arr[order];  G = G_arr[order];  Ge = Ge_arr[order]
    n  = len(L)
    if n < 2:
        return np.nan, np.nan, 0.0, 0.0, np.nan

    def _model(L_, A, b):
        return A * np.exp(b / L_)

    # seeds
    A0 = float(G[0])                             # G at largest L ≈ A
    b0 = float(np.log(max(G[-1], 1e-12) / max(A0, 1e-12)) * L[-1])  # from smallest L

    try:
        popt, pcov = curve_fit(
            _model, L, G, p0=[A0, b0],
            sigma=Ge, absolute_sigma=True,
            maxfev=10_000,
        )
        A, b = popt
        A_err = float(np.sqrt(abs(pcov[0, 0])))
        resid = G - _model(L, A, b)
        chi2  = float(np.sum((resid / Ge)**2))
        dof   = max(n - 2, 0)
        return float(A), A_err, 0.0, float(b), chi2/dof if dof > 0 else np.nan
    except Exception:
        return float(np.mean(G)), float(np.std(G)), 0.0, 0.0, np.nan

# ─────────────────────────────────────────────────────────────────────────────
# CFT model
# ─────────────────────────────────────────────────────────────────────────────
def _cft(t, A):
    return A / np.sin(np.pi * t) ** 0.25

# Colormap for t_k values (warm = large t, cool = small t)
_CMAP = plt.cm.plasma
_tk_colors = [_CMAP(0.1 + 0.8 * k / (N_SAMP - 1)) for k in range(N_SAMP - 1)]

# ─────────────────────────────────────────────────────────────────────────────
# Pre-compute all G values and fits per (cycle, t_k)
# ─────────────────────────────────────────────────────────────────────────────
# Shape: [3 cycles][7 t_k] → dict with 'test_LP','test_G','test_Ge', ... 'Ginf_test',..
fits = [[None]*len(T_SAMPLE) for _ in range(3)]

print("Pre-sampling all datasets (building interpolants once per dataset)...")
for c in range(3):
    # build interpolants and sample ALL t_k at once → O(n_datasets) not O(n_datasets * n_tk)
    test_G_all  = {}
    test_Ge_all = {}
    for L in TEST_SIZES:
        if L not in test_data: continue
        gv, ge = _sample_at_t(test_data[L], *TEST_GEOM[L], c, T_SAMPLE)
        test_G_all[L]  = gv
        test_Ge_all[L] = ge
        print(f"  cycle {c} test L={L} sampled")

    ref_G_all  = {}
    ref_Ge_all = {}
    for a_r in REF_ALPHAS:
        if a_r not in ref_data: continue
        gv, ge = _sample_at_t(ref_data[a_r], *REF_GEOM[a_r], c, T_SAMPLE)
        ref_G_all[a_r]  = gv
        ref_Ge_all[a_r] = ge
        print(f"  cycle {c} ref α={a_r} sampled")

    for ki, tk in enumerate(T_SAMPLE):
        entry = {}

        # --- TEST family (cycle-dependent L_phys: rescale by cycle-length ratio) ---
        LP_t, G_t, Ge_t = [], [], []
        for L in TEST_SIZES:
            if L not in test_G_all: continue
            LP_t.append(float(L) * _REF_BASE_LEN[c] / _REF_BASE_LEN[0])
            G_t.append(float(test_G_all[L][ki]))
            Ge_t.append(float(test_Ge_all[L][ki]))
        entry['test_LP']  = np.array(LP_t)
        entry['test_G']   = np.array(G_t)
        entry['test_Ge']  = np.array(Ge_t)
        if len(LP_t) >= 2:
            gi, gi_e, a, b, chi2 = _power_fit(LP_t, G_t, Ge_t)
            entry['Ginf_test']     = gi
            entry['Ginf_test_err'] = gi_e
            entry['a_test']        = a
            entry['b_test']        = b
            entry['chi2_test']     = chi2
        else:
            entry['Ginf_test'] = entry['Ginf_test_err'] = np.nan
            entry['a_test'] = 0.0;  entry['b_test'] = 0.0;  entry['chi2_test'] = np.nan

        # --- REF family (use physical cycle length as L_phys) ---
        LP_r, G_r, Ge_r = [], [], []
        for a_r in REF_ALPHAS:
            if a_r not in ref_G_all: continue
            LP_r.append(float(a_r) * _REF_BASE_LEN[c])
            G_r.append(float(ref_G_all[a_r][ki]))
            Ge_r.append(float(ref_Ge_all[a_r][ki]))
        entry['ref_LP']  = np.array(LP_r)
        entry['ref_G']   = np.array(G_r)
        entry['ref_Ge']  = np.array(Ge_r)
        if len(LP_r) >= 2:
            gi, gi_e, a, b, chi2 = _power_fit(LP_r, G_r, Ge_r)
            entry['Ginf_ref']     = gi
            entry['Ginf_ref_err'] = gi_e
            entry['a_ref']        = a
            entry['b_ref']        = b
            entry['chi2_ref']     = chi2
        else:
            entry['Ginf_ref'] = entry['Ginf_ref_err'] = np.nan
            entry['a_ref'] = 0.0;  entry['b_ref'] = 0.0;  entry['chi2_ref'] = np.nan

        fits[c][ki] = entry

# ─────────────────────────────────────────────────────────────────────────────
# Figure: 3 rows × 3 columns
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 3, figsize=(19, 14))

CYCLE_PTITLES = [
    r"Cycle 0   $|p^{\alpha=1}|=14.73$",
    r"Cycle 1   $|p^{\alpha=1}|=17.69$",
    r"Cycle 2   $|p^{\alpha=1}|=11.79$",
]

x_fit_dense = np.linspace(1.0/200, 0.14, 200)   # 1/L_phys (avoid 0 for power law)

for c in range(3):
    ax_fss  = axes[0][c]
    ax_om   = axes[1][c]
    ax_cont = axes[2][c]

    # ── ROW 0: FSS fit panel ────────────────────────────────────────────────
    Ginf_test_vals, Ginf_test_errs = [], []
    Ginf_ref_vals,  Ginf_ref_errs  = [], []
    b_test_vals,    b_ref_vals     = [], []

    for ki, tk in enumerate(T_SAMPLE):
        e      = fits[c][ki]
        color  = _tk_colors[ki]

        # test data points (circles) — all sizes, lighter for non-fitted ones
        n_avail_t = len(e['test_LP'])
        for idx, (lp, gv, ge) in enumerate(zip(e['test_LP'], e['test_G'], e['test_Ge'])):
            in_fit = (n_avail_t - 1 - idx) < N_FIT   # largest N_FIT
            ax_fss.errorbar(1.0/lp, gv, yerr=ge,
                            fmt='o', color=color,
                            ms=7 if in_fit else 4,
                            alpha=1.0 if in_fit else 0.35,
                            lw=1.4, capsize=3, zorder=5)

        # ref data points (squares)
        n_avail_r = len(e['ref_LP'])
        for idx, (lp, gv, ge) in enumerate(zip(e['ref_LP'], e['ref_G'], e['ref_Ge'])):
            in_fit = (n_avail_r - 1 - idx) < N_FIT
            ax_fss.errorbar(1.0/lp, gv, yerr=ge,
                            fmt='s', color=color,
                            ms=7 if in_fit else 4,
                            alpha=1.0 if in_fit else 0.35,
                            lw=1.4, capsize=3,
                            zorder=5, mfc='white', mew=1.6)

        # fit curves
        gi_t = e['Ginf_test'];  a_t = e.get('a_test', 0.0);  b_t = e.get('b_test', 0.0)
        gi_r = e['Ginf_ref'];   a_r = e.get('a_ref',  0.0);  b_r = e.get('b_ref',  0.0)

        if np.isfinite(gi_t):
            # G = A * exp(b/L) = A * exp(b*x),  x=1/L;  gi_t=A, b_t=b
            y_fit = gi_t * np.exp(b_t * x_fit_dense)
            ax_fss.plot(x_fit_dense, y_fit, color=color, lw=1.0, alpha=0.55, ls='-', zorder=2)
            ax_fss.plot(0, gi_t, marker='|', ms=12, mew=2, color=color, zorder=6)
        if np.isfinite(gi_r):
            y_fit = gi_r * np.exp(b_r * x_fit_dense)
            ax_fss.plot(x_fit_dense, y_fit, color=color, lw=1.0, alpha=0.55, ls='--', zorder=2)
            ax_fss.plot(0, gi_r, marker='|', ms=12, mew=2, color=color, zorder=6)

        Ginf_test_vals.append(gi_t);  Ginf_test_errs.append(e.get('Ginf_test_err', np.nan))
        Ginf_ref_vals.append(gi_r);   Ginf_ref_errs.append(e.get('Ginf_ref_err', np.nan))
        b_test_vals.append(b_t);      b_ref_vals.append(b_r)

    ax_fss.axvline(0, color='k', lw=0.8, ls=':', alpha=0.5)
    ax_fss.axvspan(-0.005, 0.002, alpha=0.06, color='k')
    ax_fss.set_xlim(-0.004, 0.135)
    ax_fss.set_xlabel(r"$1/L_{\rm phys}$", fontsize=11)
    ax_fss.set_ylabel(r"$G_{\rm conn}(t_k;\, L_{\rm phys})$", fontsize=11)
    ax_fss.set_title(
        CYCLE_PTITLES[c]
        + r"  — FSS: $G = A\,e^{b/L}$  (NL-WLS)"
        + f"\n(N_FIT={N_FIT} largest sizes; dim = excluded points)",
        fontsize=9.5)
    ax_fss.set_ylim(bottom=0)
    tick_L = [64, 48, 32, 24, 16, 8]
    tick_x = [1.0/L for L in tick_L]
    ax_fss.set_xticks(tick_x)
    ax_fss.set_xticklabels([f"$1/{L}$" for L in tick_L], fontsize=8)
    if c == 0:
        handles = [
            Line2D([0],[0], marker='o', color='gray', ms=7, lw=0,
                   label='Untwisted test  (circle)'),
            Line2D([0],[0], marker='s', color='gray', ms=7, lw=0,
                   mfc='white', mew=1.6, label='Twisted ref  (square, open)'),
            Line2D([0],[0], color='gray', lw=1.4, ls='-',  label='Test fit'),
            Line2D([0],[0], color='gray', lw=1.4, ls='--', label='Ref fit'),
            Line2D([0],[0], color='gray', lw=2, marker='|', ms=10,
                   label='$G_\\infty$ intercept'),
        ]
        ax_fss.legend(handles=handles, fontsize=8.0, loc='upper right', framealpha=0.90)

    # ── ROW 1: exponent coefficient b  (b<0: approach from above, b>0: from below) ──
    b_t = np.array(b_test_vals, float)
    b_r = np.array(b_ref_vals,  float)
    fin_t = np.isfinite(b_t);  fin_r = np.isfinite(b_r)

    ax_om.scatter(T_SAMPLE[fin_t], b_t[fin_t],
                  c=[_tk_colors[i] for i in np.where(fin_t)[0]],
                  s=80, marker='o', zorder=5, label='Test')
    ax_om.scatter(T_SAMPLE[fin_r], b_r[fin_r],
                  c=[_tk_colors[i] for i in np.where(fin_r)[0]],
                  s=80, marker='s', zorder=5, facecolors='none', linewidths=1.8,
                  label='Ref')
    ax_om.axhline(0, color='k', lw=1.0, ls=':', alpha=0.6)
    ax_om.set_xlim(0, 1)
    ax_om.set_xlabel(r"$t_k$", fontsize=11)
    ax_om.set_ylabel(r"Exponent $b$  in  $A\,e^{b/L}$", fontsize=11)
    ax_om.set_title(CYCLE_PTITLES[c]
                    + r"  — $b$ coefficient ($b{<}0$: approach from above, $b{>}0$: from below)",
                    fontsize=9.5)
    if c == 0:
        ax_om.legend(fontsize=9, loc='upper right', framealpha=0.9)

    # ── ROW 2: Continuum comparison ─────────────────────────────────────────
    Ginf_test_vals = np.array(Ginf_test_vals)
    Ginf_test_errs = np.array(Ginf_test_errs)
    Ginf_ref_vals  = np.array(Ginf_ref_vals)
    Ginf_ref_errs  = np.array(Ginf_ref_errs)

    # Fit CFT amplitude separately to test and ref G_inf
    t_cont = np.linspace(0.02, 0.98, 500)
    A_test = A_ref = np.nan
    for (Gv, Ge_v, col, lbl) in [
        (Ginf_test_vals, Ginf_test_errs, '#d94801', 'test'),
        (Ginf_ref_vals,  Ginf_ref_errs,  '#08306b', 'ref'),
    ]:
        fin = np.isfinite(Gv) & np.isfinite(Ge_v) & (Ge_v > 0)
        if fin.sum() >= 2:
            try:
                popt, _ = curve_fit(_cft, T_SAMPLE[fin], Gv[fin],
                                    p0=[Gv[fin].mean()],
                                    sigma=Ge_v[fin], absolute_sigma=True)
                A = float(popt[0])
                if lbl == 'test': A_test = A
                else:             A_ref  = A
                ls = '-' if lbl == 'test' else '--'
                ax_cont.plot(t_cont, _cft(t_cont, A), color=col, lw=1.6,
                             ls=ls, alpha=0.7, zorder=1,
                             label=rf"CFT fit ({lbl}): $A={A:.4f}$")
            except Exception:
                pass

    # vertical guides
    for ts in T_SAMPLE:
        ax_cont.axvline(ts, color='gray', lw=0.6, ls=':', alpha=0.4, zorder=0)

    # test continuum values
    finite_t = np.isfinite(Ginf_test_vals)
    ax_cont.errorbar(T_SAMPLE[finite_t],
                     Ginf_test_vals[finite_t], yerr=Ginf_test_errs[finite_t],
                     fmt='o', color='#d94801', ms=10, lw=1.8, capsize=4,
                     zorder=5, label=r"$G_\infty$ test ($L\to\infty$)",
                     path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])

    # ref continuum values
    finite_r = np.isfinite(Ginf_ref_vals)
    ax_cont.errorbar(T_SAMPLE[finite_r],
                     Ginf_ref_vals[finite_r], yerr=Ginf_ref_errs[finite_r],
                     fmt='s', color='#08306b', ms=10, lw=1.8, capsize=4,
                     zorder=5, label=r"$G_\infty$ ref ($\alpha\to\infty$)",
                     mfc='white', mew=1.8,
                     path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])

    # % deviation (test − ref)/ref at each t_k
    for ti, tk in enumerate(T_SAMPLE):
        gt = Ginf_test_vals[ti];  gr = Ginf_ref_vals[ti]
        if np.isfinite(gt) and np.isfinite(gr) and gr > 0:
            dev = (gt - gr) / gr * 100
            ax_cont.annotate(f"{dev:+.1f}%",
                             xy=(tk, (gt+gr)/2),
                             fontsize=6.5, ha='center', va='center',
                             color='purple',
                             path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])

    ax_cont.set_xlabel(r"Fractional cycle position $t$", fontsize=11)
    ax_cont.set_ylabel(r"$G_\infty(t)$  [continuum limit]", fontsize=11)
    ax_cont.set_title(
        CYCLE_PTITLES[c]
        + "\nContinuum comparison: test vs ref (% = (test−ref)/ref)",
        fontsize=10)
    ax_cont.set_xlim(0.0, 1.0)
    ax_cont.set_ylim(bottom=0)
    ax_cont.tick_params(labelsize=10)
    ax_cont.legend(loc='upper center', fontsize=8.5, framealpha=0.92, ncol=2)

    print(f"\nCycle {c}:  A_test={A_test:.4f}  A_ref={A_ref:.4f}")
    for ki, tk in enumerate(T_SAMPLE):
        e = fits[c][ki]
        gt = e['Ginf_test'];  gte = e.get('Ginf_test_err') or 0
        gr = e['Ginf_ref'];   gre = e.get('Ginf_ref_err')  or 0
        bt = e.get('b_test', np.nan);  br = e.get('b_ref', np.nan)
        c2t = e.get('chi2_test', np.nan);  c2r = e.get('chi2_ref', np.nan)
        dev = (gt-gr)/gr*100 if np.isfinite(gt) and np.isfinite(gr) and gr>0 else np.nan
        c2ts = f"{c2t:.2f}" if np.isfinite(c2t) else "nan"
        c2rs = f"{c2r:.2f}" if np.isfinite(c2r) else "nan"
        print(f"  t={tk:.3f}: "
              f"G\u221e_test={gt:.4f}\u00b1{gte:.4f} b={bt:.1f} chi2={c2ts}  "
              f"G\u221e_ref={gr:.4f}\u00b1{gre:.4f} b={br:.1f} chi2={c2rs}  "
              f"dev={dev:+.1f}%")

# ── Colourbar for t_k values ──────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap=_CMAP,
                            norm=plt.Normalize(vmin=T_SAMPLE[0], vmax=T_SAMPLE[-1]))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.35, pad=0.01,
                    aspect=20, location='right')
cbar.set_label(r"Sample position $t_k = k/8$", fontsize=11)
cbar.set_ticks(T_SAMPLE)
cbar.set_ticklabels([f"{t:.3f}" for t in T_SAMPLE], fontsize=8)

fig.suptitle(
    r"FSS continuum limit — 4-5-6 geometry:  $G = A\,e^{b/L}$  (NL-WLS, $N_{\rm fit}=$"
    + str(N_FIT) + r" largest sizes)"
    "\n"
    r"Circles = untwisted test  ($L\in\{" + ",".join(str(L) for L in sorted(test_data)) + r"\}$);  "
    r"Open squares = twisted iso ref  ($\alpha\in\{" + ",".join(str(a) for a in sorted(ref_data)) + r"\}$).  "
    r"Purple \% = $(G_\infty^{\rm test} - G_\infty^{\rm ref})/G_\infty^{\rm ref}$.",
    fontsize=11, y=1.01,
)

fig.tight_layout(rect=[0, 0, 0.96, 1.0])
out_path = os.path.join(_OUT, "fss_continuum_limit_456.png")
fig.savefig(out_path, dpi=160, bbox_inches='tight')
print(f"\nSaved: {out_path}")
