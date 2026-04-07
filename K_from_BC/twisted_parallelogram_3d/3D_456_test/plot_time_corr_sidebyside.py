#!/usr/bin/env python3
"""Side-by-side temporal correlator plots for s=1 and s=2 on the 4-5-6 torus.

Left panel : G_conn(dt) vs dt, both scales, with error bars (log y-scale).
             Overlay cosh fit.
Right panel: effective mass  m_eff(dt) = acosh[(G(dt-1)+G(dt+1)) / (2 G(dt))]
             evaluated at interior dt values.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
from scipy.optimize import curve_fit
import glob, re

DATADIR = Path(".")
OUTFILE = "time_corr_side_by_side.png"

# ── load ─────────────────────────────────────────────────────────────────────

def load_time(path):
    dt, corr, err, cc, ec = [], [], [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            dt.append(int(p[0]))
            corr.append(float(p[1])); err.append(float(p[2]))
            cc.append(float(p[3]));   ec.append(float(p[4]))
    return (np.array(dt), np.array(corr), np.array(err),
            np.array(cc),  np.array(ec))

def parse_geo(dirname):
    mt = re.search(r"Lx(\d+)_Ly(\d+)_Tx(-?\d+)_Ty(-?\d+)_Nt(\d+)", dirname)
    if not mt: return None
    return {k: int(v) for k, v in zip(["Lx","Ly","Tx","Ty","Nt"], mt.groups())}

# ── cosh fit  G = A * cosh(m*(dt - Nt/2)) ────────────────────────────────────

def cosh_model(dt, A, m, Nt):
    return A * np.cosh(m * (dt - Nt / 2.0))

def fit_cosh(dt, g, eg, Nt):
    def model(dt, A, m):
        return cosh_model(dt, A, m, Nt)
    try:
        p0 = [g[Nt//2], 0.5]
        popt, pcov = curve_fit(model, dt, g, p0=p0, sigma=eg, absolute_sigma=True,
                               maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception as e:
        print(f"  fit failed: {e}")
        return None, None

# ── effective mass at interior points ────────────────────────────────────────

def eff_mass(dt, g, eg, Nt):
    """acosh[(G(dt-1)+G(dt+1)) / (2 G(dt))] at interior dt=1..Nt-2."""
    pts_dt, pts_m, pts_em = [], [], []
    for i in range(1, len(dt)-1):
        num = g[i-1] + g[i+1]
        den = 2.0 * g[i]
        ratio = num / den
        if ratio <= 1.0:
            continue
        m = np.arccosh(ratio)
        # error propagation (leading order)
        dm_dg_prev  =  1.0 / (2.0 * g[i] * np.sqrt(ratio**2 - 1))
        dm_dg_next  =  dm_dg_prev
        dm_dg_curr  = -ratio / (g[i] * np.sqrt(ratio**2 - 1))
        sigma_m = np.sqrt((dm_dg_prev*eg[i-1])**2 +
                          (dm_dg_curr*eg[i]  )**2 +
                          (dm_dg_next*eg[i+1])**2)
        pts_dt.append(dt[i])
        pts_m.append(m)
        pts_em.append(sigma_m)
    return np.array(pts_dt), np.array(pts_m), np.array(pts_em)

# ── collect datasets ──────────────────────────────────────────────────────────

files = sorted(glob.glob(str(DATADIR / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_time.dat")))
if not files:
    raise SystemExit("No two_point_time.dat found.")

datasets = []
for fpath in files:
    geo = parse_geo(Path(fpath).parent.name)
    if geo is None: continue
    dt, corr, err, cc, ec = load_time(fpath)
    s = geo["Lx"] // 39
    datasets.append(dict(s=s, geo=geo, dt=dt, cc=cc, ec=ec))

datasets.sort(key=lambda d: d["s"])

# ── figure ────────────────────────────────────────────────────────────────────

colors = {1: "#e6821e", 2: "#3a86ff", 3: "#5cb85c"}
markers = {1: "o", 2: "s", 3: "^"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

for d in datasets:
    s, geo, dt, g, eg = d["s"], d["geo"], d["dt"], d["cc"], d["ec"]
    Nt = geo["Nt"]
    col = colors.get(s, "gray")
    mk  = markers.get(s, "D")
    lbl = f"s={s}  (Nt={Nt})"

    # ── left: correlator ──────────────────────────────────────────────────
    ax1.errorbar(dt, g, yerr=eg, fmt=mk, color=col, capsize=3,
                 ms=6, lw=1.3, label=lbl, zorder=4)

    # cosh fit
    popt, perr = fit_cosh(dt, g, eg, Nt)
    if popt is not None:
        dt_fine = np.linspace(0, Nt-1, 300)
        gfit = cosh_model(dt_fine, *popt, Nt)
        ax1.plot(dt_fine, gfit, color=col, lw=1.2, ls="--", alpha=0.7)
        A, m = popt
        ax1.text(0.97, 0.97 - 0.12*(s-1),
                 rf"s={s}: $m_{{fit}}={m:.4f}\pm{perr[1]:.4f}$",
                 transform=ax1.transAxes, ha="right", va="top",
                 fontsize=9, color=col)

    # ── right: effective mass ─────────────────────────────────────────────
    em_dt, em_m, em_e = eff_mass(dt, g, eg, Nt)
    if len(em_dt):
        ax2.errorbar(em_dt, em_m, yerr=em_e, fmt=mk, color=col, capsize=3,
                     ms=6, lw=1.3, label=lbl, zorder=4)

# ── left panel formatting ──────────────────────────────────────────────────
ax1.set_yscale("log")
ax1.set_xlabel(r"$\Delta t$", fontsize=12)
ax1.set_ylabel(r"$G_{\rm conn}(\Delta t)$", fontsize=12)
ax1.set_title("Temporal correlator", fontsize=12)
ax1.legend(fontsize=9)
ax1.set_xticks(range(0, max(d["geo"]["Nt"] for d in datasets)))
ax1.grid(True, which="both", alpha=0.3)

# ── right panel formatting ─────────────────────────────────────────────────
ax2.set_xlabel(r"$\Delta t$", fontsize=12)
ax2.set_ylabel(r"$m_{\rm eff}$", fontsize=12)
ax2.set_title(r"Effective mass  $m_{\rm eff} = \mathrm{arccosh}\!\left[\frac{G(\Delta t-1)+G(\Delta t+1)}{2\,G(\Delta t)}\right]$",
              fontsize=10)
ax2.legend(fontsize=9)
ax2.set_xticks(range(0, max(d["geo"]["Nt"] for d in datasets)))
ax2.grid(True, alpha=0.3)

fig.suptitle("4-5-6 torus temporal correlator — K=0.161", fontsize=13)
fig.tight_layout()

outpath = DATADIR / OUTFILE
fig.savefig(outpath, dpi=180, bbox_inches="tight")
print(f"Saved: {outpath}")
