#!/usr/bin/env python3
"""Phase 4: 2D-to-3D crossover characterization on the 4-5-6 torus.

Four-panel summary:
  Top-left  : mass gap m(Nt) vs Nt  (from cosh fits)
  Top-right : spatial correlation length xi(Nt) vs Nt
              from exponential fits G(r) ~ exp(-r/xi); overlaid with xi=Nt line
  Bot-left  : susceptibility chi(Nt) = sum_r G_conn(r) per time slice vs Nt
  Bot-right : scaling collapse residual eps(Nt) = std[ G(r)*r^(2*Delta) ] / mean vs r/Lx
              (flat => critical; large => off-critical)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
import re

DATADIR = Path(".")
K = 0.161
LX = 39
DELTA_SIGMA = 0.5182


def find_files(datadir, name):
    results = []
    for p in sorted(Path(datadir).glob(f"Lx39_Ly48_Tx9_Ty-9_Nt*_k*/*/{name}")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            results.append((int(mt.group(1)), p))
    for p in sorted(Path(datadir).glob(f"Lx39_Ly48_Tx9_Ty-9_Nt*_k*/{name}")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            nt = int(mt.group(1))
            if not any(x[0] == nt for x in results):
                results.append((nt, p))
    return sorted(results)


def load_time(path):
    dt, cc, ec = [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            dt.append(int(p[0])); cc.append(float(p[3])); ec.append(float(p[4]))
    return np.array(dt), np.array(cc), np.array(ec)


def load_spatial(path):
    m_arr, n_arr, conn, err_conn = [], [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            m_arr.append(int(p[1])); n_arr.append(int(p[2]))
            conn.append(float(p[5])); err_conn.append(float(p[6]))
    m = np.array(m_arr, float); n = np.array(n_arr, float)
    r = np.sqrt(m**2 + n**2 + m*n)
    return r, np.array(conn), np.array(err_conn)


def bin_radial(r, g, eg, nbins=40):
    r_max = r.max()
    edges = np.linspace(0, r_max + 1e-6, nbins + 1)
    rb, gb, eb = [], [], []
    for i in range(nbins):
        mask = (r >= edges[i]) & (r < edges[i+1]) & (g > 0) & (eg > 0)
        if mask.sum() == 0:
            continue
        w = 1.0 / eg[mask]**2
        wsum = w.sum()
        rb.append(0.5 * (edges[i] + edges[i+1]))
        gb.append((w * g[mask]).sum() / wsum)
        eb.append(1.0 / np.sqrt(wsum))
    return np.array(rb), np.array(gb), np.array(eb)


def fit_exponential(rb, gb, eb):
    """Fit G ~ G0 * exp(-r/xi); return (xi, dxi) or (None, None)."""
    pos = (gb > 0) & (eb > 0)
    if pos.sum() < 3:
        return None, None
    def model(r, G0, xi):
        return G0 * np.exp(-r / xi)
    try:
        popt, pcov = curve_fit(model, rb[pos], gb[pos], p0=[gb[pos][0], 5.0],
                               sigma=eb[pos], absolute_sigma=True, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        if popt[1] <= 0:
            return None, None
        return popt[1], perr[1]
    except Exception:
        return None, None


def fit_gap(dt, g, eg, Nt):
    if Nt == 1:
        return None, None
    if Nt == 2:
        if len(g) >= 2 and g[1] > 0 and g[0] / g[1] >= 1.0:
            ratio = g[0] / g[1]
            m = np.arccosh(ratio)
            sq = np.sqrt(ratio**2 - 1 + 1e-30)
            dm = np.sqrt((eg[0] / (g[1] * sq))**2 + (eg[1] * g[0] / (g[1]**2 * sq))**2)
            return m, dm
        return None, None
    mask = eg > 0
    if mask.sum() < 2:
        return None, None
    def model(dt_arr, A, m):
        return A * np.cosh(m * (dt_arr - Nt / 2.0))
    try:
        p0 = [g[mask][len(g[mask])//2], 0.3]
        popt, pcov = curve_fit(model, dt[mask], g[mask], p0=p0,
                               sigma=eg[mask], absolute_sigma=True, maxfev=8000)
        perr = np.sqrt(np.diag(pcov))
        if popt[1] <= 0:
            return None, None
        return popt[1], perr[1]
    except Exception:
        return None, None


def scaling_residual(rb, gb, eb):
    """Coefficient of variation of G(r)*r^(2*Delta) over r in [0.1*Lx, 0.5*Lx]."""
    rr = rb / LX
    mask = (rr >= 0.1) & (rr <= 0.5) & (gb > 0)
    if mask.sum() < 3:
        return None
    vals = gb[mask] * rb[mask]**(2 * DELTA_SIGMA)
    return np.std(vals) / np.mean(vals)


# ── gather all data ───────────────────────────────────────────────────────────
time_files    = dict(find_files(DATADIR, "two_point_time.dat"))
spatial_files = dict(find_files(DATADIR, "two_point_all_to_all.dat"))
all_nt = sorted(set(time_files) | set(spatial_files))

rows = []
for nt in all_nt:
    row = {"nt": nt}
    # temporal gap
    if nt in time_files:
        dt, g, eg = load_time(time_files[nt])
        row["m"], row["dm"] = fit_gap(dt, g, eg, nt)
    else:
        row["m"] = row["dm"] = None
    # spatial
    if nt in spatial_files:
        r, g, eg = load_spatial(spatial_files[nt])
        rb, gb, eb = bin_radial(r, g, eg)
        row["chi"] = float(g[g > 0].sum())   # crude susceptibility
        row["eps"] = scaling_residual(rb, gb, eb)
        row["xi"], row["dxi"] = fit_exponential(rb, gb, eb)
        row["rb"] = rb; row["gb"] = gb; row["eb"] = eb
    else:
        row["chi"] = row["eps"] = row["xi"] = row["dxi"] = None
    rows.append(row)

# ── colour map ────────────────────────────────────────────────────────────────
cmap = plt.get_cmap("plasma")
nt_vals_all = [r["nt"] for r in rows]
norm = matplotlib.colors.LogNorm(vmin=min(nt_vals_all), vmax=max(nt_vals_all))

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_gap, ax_xi = axes[0]
ax_chi, ax_eps = axes[1]

# ── Top-left: gap ─────────────────────────────────────────────────────────────
nt_g  = np.array([r["nt"] for r in rows if r["m"]  is not None])
m_g   = np.array([r["m"]  for r in rows if r["m"]  is not None])
dm_g  = np.array([r["dm"] for r in rows if r["dm"] is not None])
ax_gap.errorbar(nt_g, m_g, yerr=dm_g, fmt="o-", color="#e05c5c", ms=6, capsize=3, lw=1.5)
if len(nt_g) >= 3:
    sel = nt_g >= 4
    if sel.sum() >= 2:
        alpha, logC = np.polyfit(np.log(nt_g[sel]), np.log(m_g[sel]), 1)
        nt_fit = np.logspace(np.log10(nt_g.min()), np.log10(nt_g.max()), 200)
        ax_gap.plot(nt_fit, np.exp(logC)*nt_fit**alpha, "k--", lw=1, alpha=0.6,
                    label=rf"$\propto N_t^{{{alpha:.2f}}}$")
        ax_gap.legend(fontsize=9)
ax_gap.set_xscale("log"); ax_gap.set_yscale("log")
ax_gap.set_xlabel(r"$N_t$", fontsize=12); ax_gap.set_ylabel("mass gap $m$", fontsize=12)
ax_gap.set_title("Temporal mass gap", fontsize=11)
ax_gap.grid(True, which="both", alpha=0.25)

# ── Top-right: xi vs Nt with xi=Nt guide ──────────────────────────────────────
nt_xi = np.array([r["nt"] for r in rows if r["xi"]  is not None])
xi    = np.array([r["xi"] for r in rows if r["xi"]  is not None])
dxi   = np.array([r["dxi"] for r in rows if r["dxi"] is not None])
ax_xi.errorbar(nt_xi, xi, yerr=dxi, fmt="s-", color="#44aadd", ms=6, capsize=3, lw=1.5,
               label=r"$\xi$ (exp. fit)")
nt_guide = np.linspace(1, max(nt_xi) if len(nt_xi) else 32, 100)
ax_xi.plot(nt_guide, nt_guide, "k--", lw=1, alpha=0.6, label=r"$\xi = N_t$")
ax_xi.set_xlabel(r"$N_t$", fontsize=12); ax_xi.set_ylabel(r"$\xi$ (lattice units)", fontsize=12)
ax_xi.set_title(r"Spatial correlation length $\xi(N_t)$", fontsize=11)
ax_xi.legend(fontsize=9)
ax_xi.grid(True, alpha=0.25)

# ── Bot-left: susceptibility ──────────────────────────────────────────────────
nt_chi = [r["nt"]  for r in rows if r["chi"] is not None]
chi    = [r["chi"] for r in rows if r["chi"] is not None]
ax_chi.plot(nt_chi, chi, "o-", color="#55cc77", ms=6, lw=1.5)
ax_chi.set_xlabel(r"$N_t$", fontsize=12)
ax_chi.set_ylabel(r"$\chi = \sum_r G_{\rm conn}(r)$", fontsize=12)
ax_chi.set_title("Integrated susceptibility", fontsize=11)
ax_chi.set_yscale("log")
ax_chi.grid(True, which="both", alpha=0.25)

# ── Bot-right: scaling collapse residual ──────────────────────────────────────
nt_eps = [r["nt"]  for r in rows if r["eps"] is not None]
eps    = [r["eps"] for r in rows if r["eps"] is not None]
ax_eps.plot(nt_eps, eps, "D-", color="#cc7722", ms=6, lw=1.5)
ax_eps.axhline(0, color="gray", lw=0.8, ls="--")
ax_eps.set_xlabel(r"$N_t$", fontsize=12)
ax_eps.set_ylabel(r"$\epsilon = \sigma / \mu\;[G(r)\cdot r^{2\Delta}]$", fontsize=12)
ax_eps.set_title(r"Scaling collapse quality ($\epsilon\to 0$ at $N_t^*$)", fontsize=11)
ax_eps.grid(True, alpha=0.25)

fig.suptitle(f"4-5-6 torus Nt scan — crossover summary, K={K}, s=1 (Lx=39)", fontsize=13)
fig.tight_layout()
out = DATADIR / "crossover_vs_nt.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved: {out}")

# Summary table
print(f"\n{'Nt':>4}  {'m_gap':>9}  {'xi':>9}  {'chi':>12}  {'eps':>8}")
for r in rows:
    m_s   = f"{r['m']:.4f}"   if r['m']   is not None else "     n/a"
    xi_s  = f"{r['xi']:.2f}"  if r['xi']  is not None else "     n/a"
    chi_s = f"{r['chi']:.4f}" if r['chi'] is not None else "        n/a"
    eps_s = f"{r['eps']:.4f}" if r['eps'] is not None else "     n/a"
    print(f"{r['nt']:>4}  {m_s:>9}  {xi_s:>9}  {chi_s:>12}  {eps_s:>8}")
