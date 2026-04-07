#!/usr/bin/env python3
"""Per-Nt comparison: MC spatial correlator vs 2D modular Ising CFT.

One PNG per Nt, each with 4 panels:
  Panel 0 : MC data           — log10 connected correlator, 2×2 tiled Delaunay
  Panel 1 : Modular Ising 2D  — same colour scale; alpha fit on near-cutoff points
  Panel 2 : Ratio MC / theory — linear, centred on 1
  Panel 3 : Z-score tension   — (G_MC - alpha*G_th) / err_MC

Scaling convention: alpha is fit using only the top 20% of theory values
(near-cutoff / short-distance regime) so that the theories are matched close
to the UV cutoff rather than in the far field.

Geometry: Lx=39  Ly=48  Tx=9  Ty=-9   (4-5-6 torus, K=0.161)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import re
import sys
from pathlib import Path
from mpmath import mp
import mpmath

# ── parameters ────────────────────────────────────────────────────────────────
LX, LY, TX, TY = 39, 48, 9, -9
MP_DPS          = 25
LEVELS          = 40
CMAP_SIG        = "inferno"
CMAP_RAT        = "RdBu_r"
CMAP_ZSCORE     = "coolwarm"
# Alpha is fit using points whose distance from origin is in [ALPHA_RMIN, ALPHA_RMAX]
# as a fraction of the lattice linear scale sqrt(ncell).
ALPHA_RMIN_FRAC = 0.05
ALPHA_RMAX_FRAC = 0.10
DATADIR         = Path(".")
OUTDIR          = Path(".")

# ── lattice / torus helpers ───────────────────────────────────────────────────

def lattice_to_xy(m, n):
    return m + 0.5 * n, (np.sqrt(3) / 2) * n


def modular_tau():
    omega = 0.5 + 0.5j * np.sqrt(3)
    v = complex(LX + TY * omega)
    u = complex(TX - LY * omega)
    tau = u / v
    b_sign = 1
    if tau.imag < 0:
        tau = -tau
        b_sign = -1
    return tau, b_sign


def to_torus_uv(m, n, tau, b_sign):
    ncell = float(LX * LY + TX * TY)
    a = (LY * m + TX * n) / ncell
    b = b_sign * (TY * m - LX * n) / ncell
    a = np.mod(a, 1.0)
    b = np.mod(b, 1.0)
    return a, b


# ── file discovery ────────────────────────────────────────────────────────────

def find_files(datadir):
    """Return sorted (Nt, path) for K=0.161 data only, deduplicating by Nt."""
    result = {}
    for p in sorted(Path(datadir).glob(
            "Lx39_Ly48_Tx9_Ty-9_Nt*_k*/*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt and "_k0.161_" in str(p):
            nt = int(mt.group(1))
            if nt not in result:
                result[nt] = p
    for p in sorted(Path(datadir).glob(
            "Lx39_Ly48_Tx9_Ty-9_Nt*_k*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt and "_k0.161_" in str(p):
            nt = int(mt.group(1))
            if nt not in result:
                result[nt] = p
    return sorted(result.items())


def load_spatial(path):
    ms, ns, cc, ec = [], [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            ms.append(int(p[1]))
            ns.append(int(p[2]))
            cc.append(float(p[5]))
            ec.append(float(p[6]))
    return (np.array(ms, float), np.array(ns, float),
            np.array(cc), np.array(ec))


# ── theory ────────────────────────────────────────────────────────────────────

def compute_theory(m, n, tau, b_sign):
    """Evaluate modular Ising shape function (arXiv:2209.15546 eq.57)."""
    mp.dps = MP_DPS
    tau_mp = mpmath.mpc(tau.real, tau.imag)
    q = mpmath.exp(mpmath.pi * 1j * tau_mp)
    theta1p0 = mpmath.diff(lambda z: mpmath.jtheta(1, z, q), mpmath.mpf("0"))

    a, b = to_torus_uv(m, n, tau, b_sign)
    g_th = np.zeros(len(m))
    for i in range(len(m)):
        nu = complex(float(a[i] + b[i] * tau.real), float(b[i] * tau.imag))
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
        if (i + 1) % 200 == 0:
            print(f"  theory {i+1}/{len(m)} ...", flush=True)
    return g_th


def fit_alpha(m, n, g_mc, g_th, eg_mc,
             rmin_frac=ALPHA_RMIN_FRAC, rmax_frac=ALPHA_RMAX_FRAC):
    """Weighted least-squares alpha restricted to points whose physical distance
    r = sqrt(m²+m·n+n²) satisfies rmin_frac*L <= r <= rmax_frac*L, where
    L = sqrt(LX*LY + TX*TY) is the characteristic lattice linear scale.
    This matches the theories in the near-cutoff band (5–10% of lattice size)
    rather than in the far field."""
    L = np.sqrt(float(LX * LY + TX * TY))
    r = np.sqrt(m * m + m * n + n * n)
    rmin = rmin_frac * L
    rmax = rmax_frac * L
    mask = ((eg_mc > 0) & np.isfinite(g_mc) & np.isfinite(g_th) & (g_th > 0)
            & (r >= rmin) & (r <= rmax))
    if not mask.any():
        # fallback: widen to full near-cutoff quarter
        mask = ((eg_mc > 0) & np.isfinite(g_mc) & np.isfinite(g_th)
                & (g_th > 0) & (r <= 0.25 * L))
    if not mask.any():
        return 1.0
    w = 1.0 / eg_mc[mask] ** 2
    return float((w * g_mc[mask] * g_th[mask]).sum() /
                 (w * g_th[mask] ** 2).sum())


# ── tiling / drawing ──────────────────────────────────────────────────────────

def make_tiled(m, n, z, Lx=LX, Ly=LY, Tx=TX, Ty=TY):
    vm, vn = Lx, Ty
    um, un = Tx, -Ly
    xs, ys, zs = [], [], []
    for ia in range(2):
        for ib in range(2):
            mi = m + ia * vm + ib * um
            ni = n + ia * vn + ib * un
            xi, yi = lattice_to_xy(mi, ni)
            xs.append(xi); ys.append(yi); zs.append(z)
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    z_all = np.concatenate(zs)
    tri = mtri.Triangulation(x_all, y_all)
    verts = np.stack([x_all[tri.triangles], y_all[tri.triangles]], axis=-1)
    e01 = verts[:, 1] - verts[:, 0]
    e12 = verts[:, 2] - verts[:, 1]
    e20 = verts[:, 0] - verts[:, 2]
    max_e2 = np.maximum.reduce(
        [(e01 ** 2).sum(1), (e12 ** 2).sum(1), (e20 ** 2).sum(1)])
    tri.set_mask(max_e2 > 3.5 ** 2)
    return tri, z_all


def draw_cells(ax, Lx=LX, Ly=LY, Tx=TX, Ty=TY):
    vm, vn = Lx, Ty
    um, un = Tx, -Ly
    corners_mn = np.array([[0, 0], [vm, vn], [vm+um, vn+un], [um, un], [0, 0]], float)
    for ia in range(2):
        for ib in range(2):
            shift = np.array([ia*vm + ib*um, ia*vn + ib*un], float)
            poly = corners_mn + shift
            xb, yb = lattice_to_xy(poly[:, 0], poly[:, 1])
            ax.plot(xb, yb, color="white", lw=1.3, alpha=0.9, zorder=5)
            ax.plot(xb, yb, color="black", lw=0.4, alpha=0.5, zorder=5)


# ── main ──────────────────────────────────────────────────────────────────────

files = find_files(DATADIR)
if not files:
    sys.exit("No two_point_all_to_all.dat files found.")

print(f"Found {len(files)} Nt values: {[nt for nt, _ in files]}")

# Load reference grid for theory computation
m_ref, n_ref, _, _ = load_spatial(files[0][1])

tau, b_sign = modular_tau()
print(f"tau = {tau.real:.6f} + {tau.imag:.6f}i  (b_sign={b_sign})")
print(f"Computing modular Ising theory for {len(m_ref)} sites at dps={MP_DPS}...")
g_theory = compute_theory(m_ref, n_ref, tau, b_sign)
print("Theory computation done.")

# Shared log colour limits (across all Nt MC + theory values)
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
floor    = vmin_log - 1.0

cmap_sig    = plt.get_cmap(CMAP_SIG)
cmap_rat    = plt.get_cmap(CMAP_RAT)
cmap_zscore = plt.get_cmap(CMAP_ZSCORE)

for nt, path in files:
    print(f"\nProcessing Nt={nt} ...", flush=True)

    m, n, g_mc, eg_mc = load_spatial(path)

    # Near-cutoff band alpha fit (5–10% of lattice scale)
    alpha = fit_alpha(m, n, g_mc, g_theory, eg_mc)
    g_th_scaled = alpha * g_theory
    print(f"  alpha = {alpha:.4e}")

    # ── log10 fields ──────────────────────────────────────────────────────────
    log_mc = np.where(g_mc > 0,
                      np.log10(np.where(g_mc > 0, g_mc, 1e-30)), floor)
    log_th = np.where(np.isfinite(g_th_scaled) & (g_th_scaled > 0),
                      np.log10(np.where(g_th_scaled > 0, g_th_scaled, 1e-30)),
                      floor)

    # ── ratio ─────────────────────────────────────────────────────────────────
    ratio_raw = np.where(
        (g_mc > 0) & np.isfinite(g_th_scaled) & (g_th_scaled > 0),
        g_mc / g_th_scaled, np.nan)
    ratio_med = (float(np.nanmedian(ratio_raw[np.isfinite(ratio_raw)]))
                 if np.isfinite(ratio_raw).any() else 1.0)
    ratio = np.where(np.isfinite(ratio_raw), ratio_raw, ratio_med)

    # ── z-score (MC - theory) / err_MC ────────────────────────────────────────
    zscore_raw = np.where(
        (eg_mc > 0) & np.isfinite(g_th_scaled),
        (g_mc - g_th_scaled) / eg_mc, np.nan)
    zscore_med = (float(np.nanmedian(zscore_raw[np.isfinite(zscore_raw)]))
                  if np.isfinite(zscore_raw).any() else 0.0)
    zscore = np.where(np.isfinite(zscore_raw), zscore_raw, zscore_med)

    # ── triangulations ────────────────────────────────────────────────────────
    tri_mc,  z_mc  = make_tiled(m, n, log_mc)
    tri_th,  z_th  = make_tiled(m, n, log_th)

    # ratio: mask triangles touching originally-NaN nodes
    nan_nodes = ~np.isfinite(ratio_raw)
    nan_tiled  = np.concatenate([nan_nodes] * 4)
    tri_rat, z_rat = make_tiled(m, n, ratio)
    bad_rat = (nan_tiled[tri_rat.triangles].any(axis=1) |
               (tri_rat.mask if tri_rat.mask is not None
                else np.zeros(len(tri_rat.triangles), bool)))
    tri_rat.set_mask(bad_rat)

    # z-score: same masking
    nan_z     = ~np.isfinite(zscore_raw)
    nan_z_til = np.concatenate([nan_z] * 4)
    tri_z, z_zs = make_tiled(m, n, zscore)
    bad_z = (nan_z_til[tri_z.triangles].any(axis=1) |
             (tri_z.mask if tri_z.mask is not None
              else np.zeros(len(tri_z.triangles), bool)))
    tri_z.set_mask(bad_z)

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 13),
                             gridspec_kw={"wspace": 0.12, "hspace": 0.22})
    ax_mc, ax_th = axes[0]
    ax_rat, ax_zs = axes[1]

    # Panel 0 — MC data
    cf0 = ax_mc.tricontourf(tri_mc, z_mc, levels=LEVELS,
                            vmin=vmin_log, vmax=vmax_log,
                            cmap=cmap_sig, extend="both")
    draw_cells(ax_mc)
    ax_mc.set_aspect("equal")
    ax_mc.set_title(f"Nt={nt} — MC data", fontsize=10)
    ax_mc.set_xlabel("x", fontsize=9); ax_mc.set_ylabel("y", fontsize=9)
    ax_mc.tick_params(labelsize=8)
    plt.colorbar(cf0, ax=ax_mc, shrink=0.88, pad=0.02).set_label(
        r"$\log_{10}G$", fontsize=8)

    # Panel 1 — Modular Ising theory
    cf1 = ax_th.tricontourf(tri_th, z_th, levels=LEVELS,
                            vmin=vmin_log, vmax=vmax_log,
                            cmap=cmap_sig, extend="both")
    draw_cells(ax_th)
    ax_th.set_aspect("equal")
    ax_th.set_title(
        f"Nt={nt} — Modular Ising CFT (α={alpha:.3e}, r∈[5%,10%]·L fit)",
        fontsize=10)
    ax_th.set_xlabel("x", fontsize=9)
    ax_th.tick_params(labelsize=8)
    plt.colorbar(cf1, ax=ax_th, shrink=0.88, pad=0.02).set_label(
        r"$\log_{10}G$", fontsize=8)

    # Panel 2 — Ratio MC/theory
    ratio_finite = z_rat[np.isfinite(z_rat)]
    if ratio_finite.size > 0:
        rmed = float(np.nanmedian(ratio_finite))
        rdev = max(float(np.nanstd(ratio_finite)), 0.05)
        r_lo = max(0.01, rmed - 3 * rdev)
        r_hi = rmed + 3 * rdev
    else:
        r_lo, r_hi = 0.5, 2.0
    cf2 = ax_rat.tricontourf(tri_rat, z_rat, levels=LEVELS,
                             vmin=r_lo, vmax=r_hi,
                             cmap=cmap_rat, extend="both")
    draw_cells(ax_rat)
    ax_rat.set_aspect("equal")
    ax_rat.set_title(f"Nt={nt} — Ratio MC/Theory", fontsize=10)
    ax_rat.set_xlabel("x", fontsize=9)
    ax_rat.tick_params(labelsize=8)
    plt.colorbar(cf2, ax=ax_rat, shrink=0.88, pad=0.02).set_label(
        "MC / theory", fontsize=8)

    # Panel 3 — Z-score tension
    zs_finite = z_zs[np.isfinite(z_zs)]
    if zs_finite.size > 0:
        zs_absmax = max(float(np.percentile(np.abs(zs_finite), 98)), 0.5)
    else:
        zs_absmax = 3.0
    cf3 = ax_zs.tricontourf(tri_z, z_zs, levels=LEVELS,
                            vmin=-zs_absmax, vmax=zs_absmax,
                            cmap=cmap_zscore, extend="both")
    draw_cells(ax_zs)
    ax_zs.set_aspect("equal")
    ax_zs.set_title(f"Nt={nt} — Z-score tension", fontsize=10)
    ax_zs.set_xlabel("x", fontsize=9)
    ax_zs.tick_params(labelsize=8)
    plt.colorbar(cf3, ax=ax_zs, shrink=0.88, pad=0.02).set_label(
        r"$(G_{\rm MC} - \alpha G_{\rm th})\,/\,\sigma_{\rm MC}$", fontsize=8)

    # Compute chi²/dof for the title summary
    fit_mask = (eg_mc > 0) & np.isfinite(g_th_scaled) & np.isfinite(g_mc)
    if fit_mask.any():
        chi2_dof = float(np.mean(((g_mc[fit_mask] - g_th_scaled[fit_mask])
                                  / eg_mc[fit_mask]) ** 2))
    else:
        chi2_dof = np.nan

    fig.suptitle(
        f"4-5-6 torus (Lx=39, Ly=48, Tx=9, Ty=-9), Nt={nt}, K=0.161  —  "
        "MC spatial correlator vs. 2D modular Ising CFT\n"
        r"2×2 tiled Delaunay interpolation  |  α fit: r∈[5%,10%]·L  |  "
        rf"$\tau = {tau.real:.4f}+{tau.imag:.4f}i$  |  "
        rf"$\chi^2/\mathrm{{dof}} = {chi2_dof:.4f}$",
        fontsize=10, y=1.003
    )

    outname = OUTDIR / f"compare_modular_ising_Nt{nt:02d}.png"
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outname}")

print("\nAll done.")
