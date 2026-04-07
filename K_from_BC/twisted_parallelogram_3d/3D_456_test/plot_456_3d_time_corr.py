#!/usr/bin/env python3
"""Phase 2b: Temporal correlator analysis and effective mass extraction.

Reads two_point_time.dat from each scale factor s=1,2,3 directory and:
  - Plots G_time(tau) vs tau/Nt for each s
  - Fits to A*cosh(m*(tau - Nt/2)) to extract effective mass m_eff(s)
  - Overlays effective mass vs 1/L to check approach to continuum

File format (two_point_time.dat):
  # dt corr err corr_conn err_conn
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Physical lattice spacing (triangular lattice unit): a=1 throughout


def load_time_corr(path: Path):
    dts, corr, err, corr_conn, err_conn = [], [], [], [], []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            dts.append(int(parts[0]))
            corr.append(float(parts[1]))
            err.append(float(parts[2]))
            corr_conn.append(float(parts[3]))
            err_conn.append(float(parts[4]))
    return (np.array(dts), np.array(corr), np.array(err),
            np.array(corr_conn), np.array(err_conn))


def cosh_model(tau, A, m, Nt):
    return A * np.cosh(m * (tau - Nt / 2.0))


def fit_cosh(dts, g, ge, Nt):
    """Fit G(tau) = A * cosh(m*(tau - Nt/2)). Returns (A, m, m_err)."""
    # Use only tau in [1, Nt-1] and skip tau=0 (self-correlator)
    mid = Nt // 2
    mask = (dts >= 1) & (dts <= Nt - 1)
    x = dts[mask].astype(float)
    y = g[mask]
    ye = ge[mask]
    if len(x) < 3:
        return None, None, None
    try:
        p0 = [y[len(y) // 2], 0.5, Nt]
        popt, pcov = curve_fit(
            lambda tau, A, m: cosh_model(tau, A, m, Nt),
            x, y, p0=p0[:2], sigma=ye, absolute_sigma=True, maxfev=5000
        )
        m_err = np.sqrt(pcov[1, 1])
        return popt[0], popt[1], m_err
    except Exception:
        return None, None, None


def parse_geometry(dirname: str):
    m = re.search(r"Lx(\d+).*Nt(\d+)", dirname)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def main():
    ap = argparse.ArgumentParser(description="Temporal correlator analysis for 3D Ising 4-5-6 torus")
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output-corr", default="plot_456_3d_time_corr.png",
                    help="Output figure path for G_time(tau) overlay")
    ap.add_argument("--output-mass", default="plot_456_3d_eff_mass.png",
                    help="Output figure path for effective mass vs 1/L")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_time.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No two_point_time.dat files found matching: {pattern}")

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(files)))
    fig_corr, ax_corr = plt.subplots(figsize=(8, 5))
    fig_mass, ax_mass = plt.subplots(figsize=(6, 5))

    masses = []
    mass_errs = []
    inv_L = []

    for fpath, color in zip(files, colors):
        fpath = Path(fpath)
        dirname = fpath.parent.name
        Lx, Nt = parse_geometry(dirname)
        if Lx is None:
            continue
        s = Lx // 39

        dts, corr, err, corr_conn, err_conn = load_time_corr(fpath)

        tau_over_Nt = dts / float(Nt)
        label = f"s={s}  (Lx={Lx}, Nt={Nt})"

        # Plot connected temporal correlator vs tau/Nt
        mask = corr_conn > 0
        ax_corr.errorbar(tau_over_Nt[mask], corr_conn[mask], yerr=err_conn[mask],
                         fmt="o-", markersize=3, color=color, label=label, lw=1)

        # Fit cosh to connected correlator
        A, m, m_err = fit_cosh(dts, corr_conn, err_conn, Nt)
        if m is not None and m > 0:
            masses.append(m)
            mass_errs.append(m_err if m_err is not None else 0.0)
            inv_L.append(1.0 / float(Lx))
            # Overlay cosh fit
            tau_fine = np.linspace(0, Nt, 200)
            fit_curve = A * np.cosh(m * (tau_fine - Nt / 2.0))
            ax_corr.plot(tau_fine / Nt, fit_curve, "--", color=color, lw=0.8,
                         label=f"cosh fit m={m:.4f}±{m_err:.4f}" if m_err else f"cosh fit m={m:.4f}")

    ax_corr.set_xlabel(r"$\tau / N_t$")
    ax_corr.set_ylabel(r"$G_{\rm time}(\tau)$ (connected)")
    ax_corr.set_title("Space-averaged temporal correlator: 4-5-6 torus, K=0.161")
    ax_corr.set_yscale("log")
    ax_corr.legend(fontsize=7)
    fig_corr.tight_layout()
    out_corr = datadir / args.output_corr if not Path(args.output_corr).is_absolute() else Path(args.output_corr)
    fig_corr.savefig(out_corr, dpi=150)
    print(f"Saved: {out_corr}")

    # Effective mass vs 1/L
    if masses:
        inv_L = np.array(inv_L)
        masses = np.array(masses)
        mass_errs = np.array(mass_errs)
        ax_mass.errorbar(inv_L, masses, yerr=mass_errs, fmt="o", color="steelblue",
                         markersize=6, capsize=4, label="m_eff(s)")
        # Linear extrapolation to 1/L → 0 (continuum) if enough points
        if len(inv_L) >= 2:
            p = np.polyfit(inv_L, masses, 1)
            x_fit = np.linspace(0, inv_L.max() * 1.1, 100)
            ax_mass.plot(x_fit, np.polyval(p, x_fit), "--", color="gray",
                         label=f"linear extrap → {p[1]:.4f}")
        ax_mass.axvline(0, color="k", lw=0.5, ls=":")
        ax_mass.set_xlabel(r"$1/L_x$")
        ax_mass.set_ylabel(r"$m_{\rm eff}$")
        ax_mass.set_title("Effective temporal mass vs 1/L: 4-5-6 torus, K=0.161")
        ax_mass.legend()
        fig_mass.tight_layout()
        out_mass = datadir / args.output_mass if not Path(args.output_mass).is_absolute() else Path(args.output_mass)
        fig_mass.savefig(out_mass, dpi=150)
        print(f"Saved: {out_mass}")
    else:
        print("WARNING: no successful cosh fits; skipping effective mass plot.")


if __name__ == "__main__":
    main()
