#!/usr/bin/env python3
"""
match_456_twisted_vs_untwisted.py

Hypothesis test: does the TWISTED 4-5-6 lattice (13α,16α,-3α,3α) at isotropic
critical k1=k2=k3=1 match the UNTWISTED (L,L,0,0) lattice at the truth
anisotropic critical couplings  k = (r1, r2, 1) ≈ (5.065, 7.743, 1) ?

Key insight that fixes the previous "no-match" plots:
the cycle-index → physical-side correspondence is DIFFERENT for the two
lattices, because the bond-coupling assignment K1↔angle B, K2↔A, K3↔C
permutes which torus cycle ends up corresponding to which side of the
4-5-6 triangle.

  REF  cycles 0,1,2  -> sides 5, 6, 4
  TEST cycles 0,1,2  -> sides 5, 4, 6

Pairing by SIDE rather than cycle index gives:
  side 5 : ref c0   <->  test c0
  side 6 : ref c1   <->  test c2
  side 4 : ref c2   <->  test c1

Right observables (Z_σ-invariant) :
  R(t)  = G(t) / G(1/2)
The CFT prediction at criticality is the same on both lattices (universal),
so R_test(t) on side s should equal R_ref(t) on the same side s in the
continuum limit.

Right range :
  t ∈ [1/8, 7/8] (avoid lattice contact terms near t=0,1)
  Largest available sizes: TEST L≥32, REF α≥3

Output: results/_fss_456/match_twisted_vs_untwisted.png
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
import mc_engine
from cost import boundary_paths, _SQRT3_2, _tile_interp

# ─────────────────────────────────────────────────────────────────────────────
_DATA = os.path.join(_HERE, "results", "_fss_456")
_OUT  = os.path.join(_DATA, "match_twisted_vs_untwisted.png")

REF_ALPHAS = [1, 2, 3, 4]
REF_GEOM   = {a: (13*a, 16*a, -3*a, 3*a) for a in REF_ALPHAS}
TEST_SIZES = [8, 16, 24, 32, 48, 64]
TEST_GEOM  = {L: (L, L, 0, 0) for L in TEST_SIZES}

# Cycle index per lattice  ->  physical side label (4, 5, or 6)
REF_CYCLE_SIDE  = {0: 5, 1: 6, 2: 4}
TEST_CYCLE_SIDE = {0: 5, 1: 4, 2: 6}

# Side  -> (ref_cycle, test_cycle)
SIDE_PAIR = {
    s: (next(c for c, ss in REF_CYCLE_SIDE.items()  if ss == s),
        next(c for c, ss in TEST_CYCLE_SIDE.items() if ss == s))
    for s in (4, 5, 6)
}

T_DENSE = np.linspace(1/16, 15/16, 31)   # dense for plotting
T_FIT   = np.array([k/8 for k in range(1, 8)])  # fit on coarse grid (matches data sampling)

# ─────────────────────────────────────────────────────────────────────────────
# Loaders / interp cache
# ─────────────────────────────────────────────────────────────────────────────
def _load(p):
    if not os.path.exists(p):
        return None
    return mc_engine.load_all_to_all(p)

ref_data  = {a: _load(os.path.join(_DATA, "ref",  f"a{a}", "two_point_all_to_all.dat"))
             for a in REF_ALPHAS}
test_data = {L: _load(os.path.join(_DATA, "test", f"L{L}", "two_point_all_to_all.dat"))
             for L in TEST_SIZES}
ref_data  = {k: v for k, v in ref_data.items()  if v is not None}
test_data = {k: v for k, v in test_data.items() if v is not None}
print("Ref alphas :", sorted(ref_data.keys()))
print("Test sizes :", sorted(test_data.keys()))

_cache = {}
def _interps(dat, Lx, Ly, Tx, Ty, copies=2):
    key = (id(dat), Lx, Ly, Tx, Ty)
    if key not in _cache:
        iG  = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn",     copies)
        iGe = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn_err", copies)
        _cache[key] = (iG, iGe, boundary_paths(Lx, Ly, Tx, Ty))
    return _cache[key]

def sample_along_cycle(dat, Lx, Ly, Tx, Ty, cycle_idx, t_arr):
    iG, iGe, paths = _interps(dat, Lx, Ly, Tx, Ty)
    dm, dn = paths[cycle_idx]
    ex = dm + 0.5*dn;  ey = _SQRT3_2 * dn
    t  = np.asarray(t_arr, float)
    xy = np.column_stack([t * ex, t * ey])
    return np.asarray(iG(xy), float), np.asarray(iGe(xy), float)

# ─────────────────────────────────────────────────────────────────────────────
# CFT short-distance form  G(t) ~ A / sin(π t)^{1/4}
# ─────────────────────────────────────────────────────────────────────────────
def _cft_amp(t, G, Ge):
    """Fit A in G = A · sin(πt)^(-1/4) by WLS over t in [1/8, 7/8]."""
    s = np.sin(np.pi * np.asarray(t, float))
    f = s ** (-0.25)
    w = 1.0 / np.maximum(np.asarray(Ge, float), 1e-12)**2
    A = float(np.sum(w * f * G) / np.sum(w * f * f))
    return A

# ─────────────────────────────────────────────────────────────────────────────
# Build matched curves (largest sizes only)
# ─────────────────────────────────────────────────────────────────────────────
LARGEST_TEST = max(TEST_SIZES)            # 64
LARGEST_REF  = max(REF_ALPHAS)            # 4

print(f"\nUsing largest TEST L = {LARGEST_TEST}, largest REF α = {LARGEST_REF}")

# Per-side data
matched = {}
for side, (rc, tc) in SIDE_PAIR.items():
    Gt, Gt_e = sample_along_cycle(test_data[LARGEST_TEST], *TEST_GEOM[LARGEST_TEST], tc, T_FIT)
    Gr, Gr_e = sample_along_cycle(ref_data[LARGEST_REF],   *REF_GEOM[LARGEST_REF],   rc, T_FIT)

    # Z_σ-invariant ratio: divide by value at t=1/2 (interpolated separately)
    Gt_half, _ = sample_along_cycle(test_data[LARGEST_TEST], *TEST_GEOM[LARGEST_TEST], tc, [0.5])
    Gr_half, _ = sample_along_cycle(ref_data[LARGEST_REF],   *REF_GEOM[LARGEST_REF],   rc, [0.5])

    R_t = Gt / Gt_half[0]
    R_r = Gr / Gr_half[0]
    R_t_e = Gt_e / abs(Gt_half[0])
    R_r_e = Gr_e / abs(Gr_half[0])

    # CFT amplitudes (Z_σ × geometric factor) — ratio of amps gives Z_t/Z_r per side
    A_t = _cft_amp(T_FIT, Gt, Gt_e)
    A_r = _cft_amp(T_FIT, Gr, Gr_e)

    matched[side] = dict(
        rc=rc, tc=tc,
        Gt=Gt, Gt_e=Gt_e, Gr=Gr, Gr_e=Gr_e,
        Rt=R_t, Rt_e=R_t_e, Rr=R_r, Rr_e=R_r_e,
        A_t=A_t, A_r=A_r,
    )
    print(f"  side {side}: ref_c{rc} A_r={A_r:.4f}   test_c{tc} A_t={A_t:.4f}   "
          f"A_t/A_r={A_t/A_r:.4f}")

# Z-amplitude ratios — should be equal across sides if continuum-CFT-equivalent
amp_ratios = np.array([matched[s]['A_t'] / matched[s]['A_r'] for s in (4, 5, 6)])
mean_ratio = float(np.mean(amp_ratios))
std_ratio  = float(np.std(amp_ratios))
print(f"\nA_test/A_ref ratios per side : {amp_ratios.round(4).tolist()}")
print(f"  mean   = {mean_ratio:.4f}")
print(f"  spread = ±{std_ratio:.4f}  ({100*std_ratio/mean_ratio:.2f}% of mean)")
print(f"  → if << 1% the two lattices realise the SAME continuum CFT.")

# ─────────────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
side_order = [4, 5, 6]

t_smooth = np.linspace(1/16, 15/16, 200)
cft_shape = np.sin(np.pi * t_smooth) ** (-0.25)

for j, side in enumerate(side_order):
    m = matched[side]

    # ── ROW 0: rescaled CFT-shape comparison G·sin(πt)^{1/4} vs t ──────────
    ax = axes[0][j]
    sft = np.sin(np.pi * T_FIT) ** 0.25
    yt  = m['Gt']  * sft
    ye_t= m['Gt_e']* sft
    yr  = m['Gr']  * sft
    ye_r= m['Gr_e']* sft
    ax.errorbar(T_FIT, yt, yerr=ye_t, fmt='o', color='C0',
                label=f"TEST L={LARGEST_TEST} (untwisted, k anisotropic)", capsize=3)
    ax.errorbar(T_FIT, yr, yerr=ye_r, fmt='s', mfc='none', color='C3',
                label=f"REF α={LARGEST_REF} (twisted, k=1)",            capsize=3)
    ax.axhline(m['A_t'], color='C0', ls=':', alpha=0.7, label=f"A_test={m['A_t']:.4f}")
    ax.axhline(m['A_r'], color='C3', ls=':', alpha=0.7, label=f"A_ref ={m['A_r']:.4f}")
    ax.set_title(f"Side {side}   (ref c{m['rc']}  ↔  test c{m['tc']})", fontsize=12)
    ax.set_ylabel(r"$G(t)\cdot\sin(\pi t)^{1/4}$  (CFT-flat)")
    ax.legend(fontsize=8, loc='best')
    ax.grid(alpha=0.3)

    # ── ROW 1: Z-invariant ratio R(t) = G(t)/G(1/2) vs t ───────────────────
    ax = axes[1][j]
    ax.errorbar(T_FIT, m['Rt'], yerr=m['Rt_e'], fmt='o', color='C0',
                label="TEST  (untwisted)", capsize=3)
    ax.errorbar(T_FIT, m['Rr'], yerr=m['Rr_e'], fmt='s', mfc='none', color='C3',
                label="REF   (twisted)",   capsize=3)
    cft_R = (np.sin(np.pi * t_smooth) / np.sin(np.pi * 0.5)) ** (-0.25)
    ax.plot(t_smooth, cft_R, 'k-', lw=1.5, alpha=0.6,
            label=r"CFT $1/\sin(\pi t)^{1/4}$ (norm.)")

    # residual band
    resid = (m['Rt'] - m['Rr'])
    rmsd  = float(np.sqrt(np.mean(resid**2)))
    ax.set_xlabel("t  (fraction along cycle)")
    ax.set_ylabel(r"$R(t)=G(t)/G(1/2)$  (Z-invariant)")
    ax.set_title(f"RMSD(test−ref) = {rmsd:.4f}", fontsize=11)
    ax.legend(fontsize=8, loc='best')
    ax.grid(alpha=0.3)

fig.suptitle(
    f"Matching test: untwisted (L={LARGEST_TEST}, k anisotropic) vs twisted (α={LARGEST_REF}, k=1)\n"
    f"A_test/A_ref per side = {amp_ratios.round(4).tolist()}   "
    f"mean={mean_ratio:.4f} ±{std_ratio:.4f} ({100*std_ratio/mean_ratio:.2f}%)",
    fontsize=12,
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(_OUT, dpi=140)
print(f"\nWrote {_OUT}")
