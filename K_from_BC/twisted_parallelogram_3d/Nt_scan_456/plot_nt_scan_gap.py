#!/usr/bin/env python3
"""Phase 2: Temporal mass gap vs Nt for the 4-5-6 torus Nt-extrusion scan.

For each Nt with a two_point_time.dat file:
  - Fit G_conn(dt) = A * cosh(m * (dt - Nt/2)) to extract the mass gap m(Nt)
  - Also compute the arccosh effective mass at the midpoint as a cross-check

Left panel : G_conn(dt) for every Nt, offset for clarity (log scale)
Right panel: m(Nt) vs Nt  (log-log scale with power-law fit guide)
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


def find_time_files(datadir):
    """Return list of (Nt, path) sorted by Nt, searching both flat and nested layouts."""
    results = []
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/*/two_point_time.dat")):
        mt = re.search(r"Nt(\d+)_k", str(p))
        if mt:
            results.append((int(mt.group(1)), p))
    # Also flat layout (symlinks like Nt4)
    for p in sorted(Path(datadir).glob("Lx39_Ly48_Tx9_Ty-9_Nt*_k*/two_point_time.dat")):
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


def cosh_model(dt, A, m, Nt):
    return A * np.cosh(m * (dt - Nt / 2.0))


def fit_gap(dt, g, eg, Nt):
    """Fit cosh form; return (m, dm) or (None, None)."""
    if Nt == 1:
        return None, None
    if Nt == 2:
        # analytic: G(0)=A*cosh(m), G(1)=A  → m = arccosh(G(0)/G(1))
        if g[1] > 0 and g[0] / g[1] >= 1.0:
            m = np.arccosh(g[0] / g[1])
            # propagate errors
            dg0, dg1 = eg[0], eg[1]
            ratio = g[0] / g[1]
            dm = np.sqrt((dg0 / (g[1] * np.sqrt(ratio**2 - 1 + 1e-30)))**2 +
                         (dg1 * g[0] / (g[1]**2 * np.sqrt(ratio**2 - 1 + 1e-30)))**2)
            return m, dm
        return None, None
    # General: fit excluding dt=0 (err=0); use all other points
    mask = eg > 0
    if mask.sum() < 2:
        return None, None
    def model(dt, A, m):
        return cosh_model(dt, A, m, Nt)
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


def eff_mass_midpoint(dt, g, eg, Nt):
    """Return arccosh effective mass at the midpoint dt ~ Nt//2."""
    mid = Nt // 2
    idx = np.where(dt == mid)[0]
    if not idx.size:
        return None, None
    i = idx[0]
    if i == 0 or i == len(dt) - 1:
        return None, None
    num = g[i-1] + g[i+1]
    den = 2.0 * g[i]
    ratio = num / den
    if ratio <= 1.0:
        return None, None
    m = np.arccosh(ratio)
    dm_dg_prev = 1.0 / (2 * g[i] * np.sqrt(ratio**2 - 1 + 1e-30))
    dm_dg_next = dm_dg_prev
    dm_dg_curr = -ratio / (g[i] * np.sqrt(ratio**2 - 1 + 1e-30))
    sigma = np.sqrt((dm_dg_prev*eg[i-1])**2 + (dm_dg_curr*eg[i])**2 + (dm_dg_next*eg[i+1])**2)
    return m, sigma


# ── collect ───────────────────────────────────────────────────────────────────
files = find_time_files(DATADIR)
if not files:
    raise SystemExit("No two_point_time.dat files found.")

datasets = []
for nt, path in files:
    dt, g, eg = load_time(path)
    m_fit, dm_fit = fit_gap(dt, g, eg, nt)
    m_mid, dm_mid = eff_mass_midpoint(dt, g, eg, nt)
    datasets.append(dict(nt=nt, dt=dt, g=g, eg=eg, m_fit=m_fit, dm_fit=dm_fit,
                         m_mid=m_mid, dm_mid=dm_mid))

# ── colour map: one colour per Nt ────────────────────────────────────────────
cmap = plt.get_cmap("plasma")
nt_vals = [d["nt"] for d in datasets]
norm = matplotlib.colors.LogNorm(vmin=min(nt_vals), vmax=max(nt_vals))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# ── left: G(dt) curves ───────────────────────────────────────────────────────
for d in datasets:
    col = cmap(norm(d["nt"]))
    lbl = f"Nt={d['nt']}"
    mask = d["eg"] > 0
    ax1.errorbar(d["dt"][mask] / d["nt"], d["g"][mask], yerr=d["eg"][mask],
                 fmt="o-", color=col, ms=4, lw=1.2, capsize=2, label=lbl)

ax1.set_yscale("log")
ax1.set_xlabel(r"$\Delta t / N_t$", fontsize=12)
ax1.set_ylabel(r"$G_{\rm conn}(\Delta t)$", fontsize=12)
ax1.set_title("Temporal correlator (all Nt)", fontsize=11)
ax1.legend(fontsize=7, ncol=2, loc="upper right")
ax1.grid(True, which="both", alpha=0.25)

# ── right: m(Nt) vs Nt ───────────────────────────────────────────────────────
nt_fit, m_fit_arr, dm_fit_arr = [], [], []
nt_mid, m_mid_arr, dm_mid_arr = [], [], []
for d in datasets:
    if d["m_fit"] is not None:
        nt_fit.append(d["nt"]); m_fit_arr.append(d["m_fit"]); dm_fit_arr.append(d["dm_fit"])
    if d["m_mid"] is not None:
        nt_mid.append(d["nt"]); m_mid_arr.append(d["m_mid"]); dm_mid_arr.append(d["dm_mid"])

nt_fit = np.array(nt_fit); m_fit_arr = np.array(m_fit_arr); dm_fit_arr = np.array(dm_fit_arr)
nt_mid = np.array(nt_mid); m_mid_arr = np.array(m_mid_arr); dm_mid_arr = np.array(dm_mid_arr)

ax2.errorbar(nt_fit, m_fit_arr, yerr=dm_fit_arr, fmt="o", color="#e05c5c",
             ms=6, capsize=3, lw=1.4, label="cosh fit")
ax2.errorbar(nt_mid, m_mid_arr, yerr=dm_mid_arr, fmt="s", color="#5588dd",
             ms=6, capsize=3, lw=1.4, label=r"$m_{\rm eff}$ at midpoint", alpha=0.8)

# guide: 1/Nt power law (relative normalisation from largest Nt)
if len(nt_fit) >= 2:
    # fit power law m ~ C * Nt^{-alpha} to the fit values at larger Nt
    sel = nt_fit >= 4
    if sel.sum() >= 2:
        lnt = np.log(nt_fit[sel]); lm = np.log(m_fit_arr[sel])
        alpha, logC = np.polyfit(lnt, lm, 1)
        nt_guide = np.logspace(np.log10(nt_fit.min()), np.log10(nt_fit.max()), 200)
        ax2.plot(nt_guide, np.exp(logC) * nt_guide**alpha, "k--", lw=1,
                 alpha=0.6, label=rf"$\propto N_t^{{{alpha:.2f}}}$")

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r"$N_t$", fontsize=12)
ax2.set_ylabel(r"mass gap $m$  (lattice units)", fontsize=12)
ax2.set_title(r"Mass gap $m(N_t)$ — 2D-to-3D crossover", fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, which="both", alpha=0.25)

fig.suptitle(f"4-5-6 torus Nt scan — K={K}, s=1 (Lx=39, Ly=48)", fontsize=13)
fig.tight_layout()
out = DATADIR / "gap_vs_nt.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved: {out}")

# Print table
print(f"\n{'Nt':>4}  {'m_fit':>10}  {'dm_fit':>10}  {'m_mid':>10}  {'dm_mid':>10}")
for d in datasets:
    mf  = f"{d['m_fit']:.5f}"  if d['m_fit']  is not None else "    n/a  "
    dmf = f"{d['dm_fit']:.5f}" if d['dm_fit'] is not None else "    n/a  "
    mm  = f"{d['m_mid']:.5f}"  if d['m_mid']  is not None else "    n/a  "
    dmm = f"{d['dm_mid']:.5f}" if d['dm_mid'] is not None else "    n/a  "
    print(f"{d['nt']:>4}  {mf:>10}  {dmf:>10}  {mm:>10}  {dmm:>10}")
