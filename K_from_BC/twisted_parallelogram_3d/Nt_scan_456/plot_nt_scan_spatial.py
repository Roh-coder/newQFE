#!/usr/bin/env python3
"""Phase 3: Spatial correlator evolution with Nt on the 4-5-6 torus.

Left panel : G_conn(r) vs r for all Nt, with exponential fits at small Nt.
Right panel: G_conn(r) * r^(2*Delta) vs r/Lx for all Nt.
             At 3D criticality (large Nt) this should flatten => horizontal.
             The 3D Ising spin dimension is Delta_sigma = 0.5182.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
import glob, re

DATADIR = Path(".")
K = 0.161
LX = 39
DELTA_SIGMA = 0.5182


def find_spatial_files(datadir):
    results = []
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            results.append((int(mt.group(1)), p))
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/two_point_all_to_all.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            nt = int(mt.group(1))
            if not any(x[0] == nt for x in results):
                results.append((nt, p))
    return sorted(results)


def load_spatial(path):
    m, n, conn, err_conn = [], [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            m.append(int(p[1])); n.append(int(p[2]))
            conn.append(float(p[5])); err_conn.append(float(p[6]))
    m = np.array(m, dtype=float); n = np.array(n, dtype=float)
    conn = np.array(conn); err_conn = np.array(err_conn)
    r = np.sqrt(m**2 + n**2 + m*n)
    return r, conn, err_conn


def bin_radial(r, g, eg, nbins=40):
    """Bin by distance; return (r_bin, g_bin, eg_bin)."""
    r_max = r.max()
    edges = np.linspace(0, r_max + 1e-6, nbins + 1)
    rb, gb, eb = [], [], []
    for i in range(nbins):
        mask = (r >= edges[i]) & (r < edges[i+1]) & (g > 0)
        if mask.sum() == 0:
            continue
        # inverse-variance weighted mean
        w = 1.0 / np.where(eg[mask] > 0, eg[mask]**2, 1e-20)
        wsum = w.sum()
        g_mean = (w * g[mask]).sum() / wsum
        e_mean = 1.0 / np.sqrt(wsum)
        rb.append(0.5 * (edges[i] + edges[i+1]))
        gb.append(g_mean)
        eb.append(e_mean)
    return np.array(rb), np.array(gb), np.array(eb)


# ── collect ───────────────────────────────────────────────────────────────────
files = find_spatial_files(DATADIR)
if not files:
    raise SystemExit("No two_point_all_to_all.dat files found.")

datasets = []
for nt, path in files:
    r, g, eg = load_spatial(path)
    rb, gb, eb = bin_radial(r, g, eg)
    datasets.append(dict(nt=nt, r=rb, g=gb, eg=eb))

cmap = plt.get_cmap("plasma")
nt_vals = [d["nt"] for d in datasets]
norm = matplotlib.colors.LogNorm(vmin=min(nt_vals), vmax=max(nt_vals))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

for d in datasets:
    col = cmap(norm(d["nt"]))
    r, g, eg = d["r"], d["g"], d["eg"]
    pos = (g > 0) & (eg > 0)
    lbl = f"Nt={d['nt']}"

    # ── left: raw G(r), log-log ───────────────────────────────────────────
    ax1.errorbar(r[pos], g[pos], yerr=eg[pos], fmt=".", color=col,
                 ms=4, lw=0.8, capsize=1.5, label=lbl, alpha=0.85)

    # Exponential fit for small-Nt runs (Nt <= 6) where system is off-critical
    if d["nt"] <= 6 and pos.sum() >= 4:
        def exp_model(r, G0, xi):
            return G0 * np.exp(-r / xi)
        try:
            popt, _ = curve_fit(exp_model, r[pos], g[pos], p0=[g[pos][0], 5.0],
                                sigma=eg[pos], absolute_sigma=True, maxfev=5000)
            r_fit = np.linspace(r[pos].min(), r[pos].max(), 200)
            ax1.plot(r_fit, exp_model(r_fit, *popt), "-", color=col, lw=1.0, alpha=0.5)
        except Exception:
            pass

    # ── right: G(r) * r^(2*Delta) vs r/Lx ───────────────────────────────
    scaled_g  = g[pos] * r[pos]**(2 * DELTA_SIGMA)
    scaled_eg = eg[pos] * r[pos]**(2 * DELTA_SIGMA)
    ax2.errorbar(r[pos] / LX, scaled_g, yerr=scaled_eg, fmt=".", color=col,
                 ms=4, lw=0.8, capsize=1.5, label=lbl, alpha=0.85)

ax1.set_xscale("log"); ax1.set_yscale("log")
ax1.set_xlabel(r"$r$ (lattice units)", fontsize=12)
ax1.set_ylabel(r"$G_{\rm conn}(r)$", fontsize=12)
ax1.set_title(r"Spatial correlator $G_{\rm conn}(r)$", fontsize=11)
ax1.legend(fontsize=7, ncol=2)
ax1.grid(True, which="both", alpha=0.25)

# power-law guide r^{-2*Delta} shape on right panel
ax2.set_xlabel(r"$r / L_x$", fontsize=12)
ax2.set_ylabel(r"$G_{\rm conn}(r) \cdot r^{2\Delta_\sigma}$  ($\Delta_\sigma=0.5182$)", fontsize=12)
ax2.set_title(r"Scaling collapse check: flat $\Rightarrow$ 3D critical", fontsize=11)
ax2.legend(fontsize=7, ncol=2)
ax2.grid(True, which="both", alpha=0.25)

# colourbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.85, pad=0.02)
cbar.set_label(r"$N_t$", fontsize=11)
cbar.set_ticks(nt_vals)
cbar.set_ticklabels([str(x) for x in nt_vals])

fig.suptitle(f"4-5-6 torus Nt scan — spatial correlator evolution, K={K}, s=1", fontsize=12)
fig.tight_layout()
out = DATADIR / "spatial_corr_vs_nt.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved: {out}")
