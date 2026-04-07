#!/usr/bin/env python3
"""Compare all-to-all spatial correlator (all Nt) to the 2D modular Ising CFT.

Layout: one row per Nt value.
  Column 0 : MC data          — 2×2 tiled, piecewise-linear (Delaunay) interpolation, log10 scale
  Column 1 : Modular Ising 2D — same tiling, same scale; alpha fit per Nt
  Column 2 : Ratio MC / theory (linear scale, centered on 1)

The 2D modular Ising CFT prediction from arXiv:2209.15546 eq.(57):
  G(nu|tau) ~ |theta1'(0|tau)/theta1(pi*nu|tau)|^{1/4}
              * [|theta2(pi*nu/2)| + |theta3(pi*nu/2)| + |theta4(pi*nu/2)|]
              / [|theta2(0)| + |theta3(0)| + |theta4(0)|]

The theory depends only on (m,n) displacement and the geometric tau — not on Nt.
It is computed once and reused for all rows; only the normalization alpha is
refitted per Nt by weighted least squares.

Geometry (fixed throughout):
  Lx=39  Ly=48  Tx=9  Ty=-9   =>  tau = u/v with Im(tau) > 0
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mcolors
from pathlib import Path
from mpmath import mp
import mpmath
import re, glob, sys

# ── parameters ────────────────────────────────────────────────────────────────
LX, LY, TX, TY = 39, 48, 9, -9
MP_DPS   = 25          # mpmath precision; increase for higher accuracy
LEVELS   = 40          # tricontourf contour levels
CMAP_SIG = "inferno"   # for log10 correlator panels
CMAP_RAT = "RdBu_r"    # for ratio panel
OUTFILE  = "compare_modular_ising_vs_nt_cutoff_scaled.png"
DATADIR  = Path(".")

# ── helpers ───────────────────────────────────────────────────────────────────

def lattice_to_xy(m, n):
    return m + 0.5*n, (np.sqrt(3)/2)*n


def find_files(datadir):
    """Return sorted list of (Nt, path) for all available two_point_all_to_all.dat."""
    result = []
    # nested layout
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            result.append((int(mt.group(1)), p))
    # flat / symlink layout
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            nt = int(mt.group(1))
            if not any(x[0] == nt for x in result):
                result.append((nt, p))
    return sorted(result)


def load_spatial(path):
    ms, ns, cc, ec = [], [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            ms.append(int(p[1])); ns.append(int(p[2]))
            cc.append(float(p[5])); ec.append(float(p[6]))
    m = np.array(ms, float); n = np.array(ns, float)
    return m, n, np.array(cc), np.array(ec)


def modular_tau():
    """Return (tau, b_sign) for the 4-5-6 torus BC."""
    omega = 0.5 + 0.5j*np.sqrt(3)
    v = complex(LX + TY*omega)
    u = complex(TX - LY*omega)
    tau = u / v
    b_sign = 1
    if tau.imag < 0:
        tau = -tau
        b_sign = -1
    return tau, b_sign


def to_torus_uv(m, n, tau, b_sign):
    """Map lattice (m,n) to torus coordinates (a,b) with nu = a + b*tau."""
    ncell = float(LX*LY + TX*TY)
    a = (LY*m + TX*n) / ncell
    b = b_sign * (TY*m - LX*n) / ncell
    a = np.mod(a, 1.0)
    b = np.mod(b, 1.0)
    return a, b


def compute_theory(m, n, tau, b_sign):
    """Evaluate the modular Ising shape function at each (m,n).
    Returns g_th array (same length as m), with nan at singularities."""
    mp.dps = MP_DPS
    tau_mp = mpmath.mpc(tau.real, tau.imag)
    q = mpmath.exp(mpmath.pi * 1j * tau_mp)
    theta1p0 = mpmath.diff(lambda z: mpmath.jtheta(1, z, q), mpmath.mpf("0"))

    a, b = to_torus_uv(m, n, tau, b_sign)
    g_th = np.zeros(len(m))
    for i in range(len(m)):
        nu = complex(float(a[i] + b[i]*tau.real), float(b[i]*tau.imag))
        if abs(nu.real) < 1e-13 and abs(nu.imag) < 1e-13:
            g_th[i] = np.nan
            continue
        z  = mpmath.pi * mpmath.mpc(nu.real, nu.imag)
        z2 = 0.5 * z
        th1 = mpmath.jtheta(1, z, q)
        if abs(th1) == 0:
            g_th[i] = np.nan
            continue
        pref = abs(theta1p0 / th1) ** mpmath.mpf("0.25")
        num = (abs(mpmath.jtheta(1, z2, q)) + abs(mpmath.jtheta(2, z2, q)) +
               abs(mpmath.jtheta(3, z2, q)) + abs(mpmath.jtheta(4, z2, q)))
        den = (abs(mpmath.jtheta(2, 0, q)) + abs(mpmath.jtheta(3, 0, q)) +
               abs(mpmath.jtheta(4, 0, q)))
        if den == 0:
            g_th[i] = np.nan
        else:
            val = float(pref * (num / den))
            g_th[i] = val if np.isfinite(val) else np.nan
        if (i+1) % 200 == 0:
            print(f"  theory {i+1}/{len(m)} ...", flush=True)
    return g_th


def fit_alpha(m, n, g_mc, g_th, eg_mc,
             rmin_frac=0.05, rmax_frac=0.10):
    """Weighted least-squares alpha restricted to points whose physical distance
    r = sqrt(m²+m·n+n²) satisfies rmin_frac*L <= r <= rmax_frac*L, where
    L = sqrt(LX*LY + TX*TY) is the characteristic lattice linear scale.
    Matches theories in the 5–10% near-cutoff band rather than the far field."""
    L = np.sqrt(float(LX * LY + TX * TY))
    r = np.sqrt(m * m + m * n + n * n)
    rmin = rmin_frac * L
    rmax = rmax_frac * L
    mask = ((eg_mc > 0) & np.isfinite(g_mc) & np.isfinite(g_th) & (g_th > 0)
            & (r >= rmin) & (r <= rmax))
    if not mask.any():
        mask = ((eg_mc > 0) & np.isfinite(g_mc) & np.isfinite(g_th)
                & (g_th > 0) & (r <= 0.25 * L))
    if not mask.any():
        return 1.0
    w = 1.0 / eg_mc[mask]**2
    return float((w * g_mc[mask] * g_th[mask]).sum() / (w * g_th[mask]**2).sum())


def make_tiled(m, n, z, Lx=LX, Ly=LY, Tx=TX, Ty=TY):
    """Return 2×2 tiled (x_all, y_all, z_all), masking long Delaunay edges."""
    vm, vn = Lx, Ty
    um, un = Tx, -Ly
    xs, ys, zs = [], [], []
    for ia in range(2):
        for ib in range(2):
            mi = m + ia*vm + ib*um
            ni = n + ia*vn + ib*un
            xi, yi = lattice_to_xy(mi, ni)
            xs.append(xi); ys.append(yi); zs.append(z)
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    z_all = np.concatenate(zs)
    tri = mtri.Triangulation(x_all, y_all)
    # mask long edges (artefact triangles spanning tile boundaries)
    verts = np.stack([x_all[tri.triangles], y_all[tri.triangles]], axis=-1)
    e01 = verts[:,1]-verts[:,0]; e12 = verts[:,2]-verts[:,1]; e20 = verts[:,0]-verts[:,2]
    max_e2 = np.maximum.reduce([(e01**2).sum(1), (e12**2).sum(1), (e20**2).sum(1)])
    tri.set_mask(max_e2 > 3.5**2)
    return tri, z_all


def draw_cells(ax, Lx=LX, Ly=LY, Tx=TX, Ty=TY):
    """Outline the 4 tiled fundamental cells in white/black."""
    vm, vn = Lx, Ty
    um, un = Tx, -Ly
    corners_mn = np.array([[0,0],[vm,vn],[vm+um,vn+un],[um,un],[0,0]], float)
    for ia in range(2):
        for ib in range(2):
            shift = np.array([ia*vm+ib*um, ia*vn+ib*un], float)
            poly = corners_mn + shift
            xb, yb = lattice_to_xy(poly[:,0], poly[:,1])
            ax.plot(xb, yb, color="white", lw=1.3, alpha=0.9, zorder=5)
            ax.plot(xb, yb, color="black", lw=0.4, alpha=0.5, zorder=5)


# ── main ──────────────────────────────────────────────────────────────────────

files = find_files(DATADIR)
if not files:
    sys.exit("No two_point_all_to_all.dat files found.")

print(f"Found {len(files)} Nt values: {[nt for nt,_ in files]}")

# Load first file to get (m,n) reference grid (same for all Nt)
m_ref, n_ref, _, _ = load_spatial(files[0][1])

# Compute modular Ising theory once on the reference grid
tau, b_sign = modular_tau()
print(f"tau = {tau.real:.6f} + {tau.imag:.6f}i  (b_sign={b_sign})")
print(f"Computing modular Ising theory for {len(m_ref)} sites at dps={MP_DPS}...")
g_theory = compute_theory(m_ref, n_ref, tau, b_sign)
print("Theory computation done.")

# ── figure: N_nt rows × 3 columns ────────────────────────────────────────────
N = len(files)
fig, axes = plt.subplots(N, 3, figsize=(18, 5.5*N),
                         gridspec_kw={"wspace": 0.05, "hspace": 0.12})
if N == 1:
    axes = axes[None, :]  # ensure 2D

# Compute shared log colour limits across all MC + theory panels
all_logvals = []
for nt, path in files:
    _, _, cc, _ = load_spatial(path)
    pos = cc > 0
    if pos.any():
        all_logvals.append(np.log10(cc[pos]))
g_th_pos = g_theory[np.isfinite(g_theory) & (g_theory > 0)]
if g_th_pos.size:
    all_logvals.append(np.log10(g_th_pos))
all_logvals = np.concatenate(all_logvals)
vmin_log = float(np.percentile(all_logvals, 2))
vmax_log = float(np.percentile(all_logvals, 99))
floor = vmin_log - 1.0

cmap_sig = plt.get_cmap(CMAP_SIG)
cmap_rat = plt.get_cmap(CMAP_RAT)

for row, (nt, path) in enumerate(files):
    ax_mc, ax_th, ax_rat = axes[row]

    m, n, g_mc, eg_mc = load_spatial(path)

    # Fit normalization alpha using r∈[5%,10%]·L band
    alpha = fit_alpha(m, n, g_mc, g_theory, eg_mc)
    g_th_scaled = alpha * g_theory

    # log10 fields (floor for non-positive)
    log_mc = np.where(g_mc > 0, np.log10(np.where(g_mc > 0, g_mc, 1e-30)), floor)
    log_th = np.where(np.isfinite(g_th_scaled) & (g_th_scaled > 0),
                      np.log10(np.where(g_th_scaled > 0, g_th_scaled, 1e-30)), floor)

    # 2×2 tiled triangulations
    tri_mc,  z_mc  = make_tiled(m, n, log_mc)
    tri_th,  z_th  = make_tiled(m, n, log_th)

    # ratio MC/theory (only where both finite and positive); replace NaN with median
    ratio_raw = np.where((g_mc > 0) & np.isfinite(g_th_scaled) & (g_th_scaled > 0),
                         g_mc / g_th_scaled, np.nan)
    ratio_med = float(np.nanmedian(ratio_raw[np.isfinite(ratio_raw)])) if np.isfinite(ratio_raw).any() else 1.0
    ratio = np.where(np.isfinite(ratio_raw), ratio_raw, ratio_med)
    tri_rat, z_rat = make_tiled(m, n, ratio)
    # Mask ratio triangulation triangles that touch originally-NaN nodes
    nan_nodes = ~np.isfinite(ratio_raw)
    # tile the nan mask same way as make_tiled (4 copies concatenated)
    nan_tiled = np.concatenate([nan_nodes]*4)
    bad_tris = nan_tiled[tri_rat.triangles].any(axis=1) | (tri_rat.mask if tri_rat.mask is not None else np.zeros(len(tri_rat.triangles), bool))
    tri_rat.set_mask(bad_tris)

    # MC panel
    cf0 = ax_mc.tricontourf(tri_mc, z_mc, levels=LEVELS,
                            vmin=vmin_log, vmax=vmax_log, cmap=cmap_sig, extend="both")
    draw_cells(ax_mc)
    ax_mc.set_aspect("equal"); ax_mc.set_title(f"Nt={nt} — MC data", fontsize=9)
    ax_mc.set_xlabel("x", fontsize=8); ax_mc.set_ylabel("y", fontsize=8)
    ax_mc.tick_params(labelsize=7)

    # Theory panel
    cf1 = ax_th.tricontourf(tri_th, z_th, levels=LEVELS,
                            vmin=vmin_log, vmax=vmax_log, cmap=cmap_sig, extend="both")
    draw_cells(ax_th)
    ax_th.set_aspect("equal")
    ax_th.set_title(f"Nt={nt} — Modular Ising CFT (α={alpha:.3e}, r∈[5%,10%]·L fit)", fontsize=9)
    ax_th.set_xlabel("x", fontsize=8); ax_th.tick_params(labelsize=7)
    ax_th.set_yticklabels([])

    # Ratio panel (log scale around 1)
    ratio_finite = z_rat[np.isfinite(z_rat)]
    if ratio_finite.size > 0:
        rmed = float(np.nanmedian(ratio_finite))
        rdev = max(float(np.nanstd(ratio_finite)), 0.05)
        r_lo = max(0.01, rmed - 3*rdev)
        r_hi = rmed + 3*rdev
    else:
        r_lo, r_hi = 0.5, 2.0
    cf2 = ax_rat.tricontourf(tri_rat, z_rat, levels=LEVELS,
                             vmin=r_lo, vmax=r_hi, cmap=cmap_rat, extend="both")
    draw_cells(ax_rat)
    ax_rat.set_aspect("equal")
    ax_rat.set_title(f"Nt={nt} — Ratio MC/Theory", fontsize=9)
    ax_rat.set_xlabel("x", fontsize=8); ax_rat.tick_params(labelsize=7)
    ax_rat.set_yticklabels([])

    # per-row colorbars
    plt.colorbar(cf0, ax=ax_mc, shrink=0.88, pad=0.02).set_label(r"$\log_{10}G$", fontsize=7)
    plt.colorbar(cf1, ax=ax_th, shrink=0.88, pad=0.02).set_label(r"$\log_{10}G$", fontsize=7)
    plt.colorbar(cf2, ax=ax_rat, shrink=0.88, pad=0.02).set_label("MC/theory", fontsize=7)

    print(f"  Nt={nt:>2} done  alpha={alpha:.4e}")

fig.suptitle(
    "4-5-6 torus (Lx=39, Ly=48, Tx=9, Ty=-9), K=0.161  —  "
    "MC spatial correlator vs. 2D modular Ising CFT\n"
    r"2×2 tiled Delaunay interpolation  |  α fit: r∈[5%,10%]·L  |  "
    rf"$\tau = {tau.real:.4f}+{tau.imag:.4f}i$",
    fontsize=11, y=1.002
)

outpath = DATADIR / OUTFILE
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nSaved: {outpath}")
