#!/usr/bin/env python3
"""Phase 3: Compare 3D MC spatial correlator to 2D analytical modular Ising solution.

For each scale factor s=1,2,3, reads the equal-time spatial two_point_all_to_all.dat
and overlays it against the analytical 2D Ising CFT two-point function on the
4-5-6 torus with boundary conditions Lx=39s, Ly=48s, Tx=9s, Ty=-9s.

The 2D modular Ising CFT spin-spin two-point function on the torus is:

  G(nu | tau) ~ |theta1'(0|tau) / theta1(pi*nu | tau)|^(1/4)
                 * [|theta1(pi*nu/2)| + |theta2(pi*nu/2)| + |theta3(pi*nu/2)|
                    + |theta4(pi*nu/2)|]
                 / [|theta2(0)| + |theta3(0)| + |theta4(0)|]

normalized to the MC data by weighted least squares.

References: arXiv:2209.15546 eq.(57); Di Francesco, Saleur, Zuber.

Usage:
  python compare_456_3d_vs_2d_modular.py --datadir . --output compare_2d_modular.pdf
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from mpmath import mp, jtheta, diff, pi, mpc, mpf, e as mp_e
except ImportError:
    sys.exit("ERROR: mpmath is required. Install with: pip install mpmath")

# ── 2D modular Ising CFT utilities ──────────────────────────────────────────

def modular_tau(lx: int, ly: int, tx: int, ty: int):
    """Compute the modular parameter tau for the twisted parallelogram."""
    omega = 0.5 + 0.5j * np.sqrt(3.0)
    v = complex(lx + ty * omega)
    u = complex(tx - ly * omega)
    tau = u / v
    b_sign = 1
    if tau.imag < 0:
        tau = -tau
        b_sign = -1
    return tau, b_sign


def torus_coords(m, n, lx, ly, tx, ty, b_sign):
    """Map triangular-lattice displacement (m,n) to torus coordinates (a,b)."""
    ncell = float(lx * ly + tx * ty)
    a = (ly * m + tx * n) / ncell
    b = b_sign * (ty * m - lx * n) / ncell
    return np.mod(a, 1.0), np.mod(b, 1.0)


def ising_torus_shape(nu_arr, tau, mp_dps=40):
    """Evaluate the 2D Ising CFT two-point shape at an array of nu = a + b*tau values.

    Returns an array of real positive values (unnormalized).
    """
    mp.dps = mp_dps
    q = mp_e ** (pi * 1j * mpc(tau.real, tau.imag))
    theta1p0 = diff(lambda zz: jtheta(1, zz, q), mpf("0"))

    g_th = np.zeros(len(nu_arr))
    for i, nu in enumerate(nu_arr):
        z = pi * mpc(nu.real, nu.imag)
        z2 = 0.5 * z
        th1 = jtheta(1, z, q)
        if abs(th1) == 0:
            g_th[i] = np.nan
            continue
        pref = abs(theta1p0 / th1) ** mpf("0.25")
        num = (abs(jtheta(1, z2, q)) + abs(jtheta(2, z2, q))
               + abs(jtheta(3, z2, q)) + abs(jtheta(4, z2, q)))
        den = (abs(jtheta(2, mpf("0"), q)) + abs(jtheta(3, mpf("0"), q))
               + abs(jtheta(4, mpf("0"), q)))
        g_th[i] = float(pref * num / den) if den != 0 else np.nan
    return g_th


# ── Data loading ─────────────────────────────────────────────────────────────

def load_connected(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            rows.append((int(parts[1]), int(parts[2]), float(parts[5]), float(parts[6])))
    if not rows:
        raise SystemExit(f"No data in {path}")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def parse_geometry(dirname: str):
    m = re.search(r"Lx(\d+)_Ly(\d+)_Tx(-?\d+)_Ty(-?\d+)_Nt(\d+)", dirname)
    if not m:
        return None
    return {k: int(v) for k, v in zip(["Lx", "Ly", "Tx", "Ty", "Nt"], m.groups())}


# ── Plotting ──────────────────────────────────────────────────────────────────

def radial_bin_weighted(r, vals, errs, nbins=50):
    rmax = np.nanmax(r[np.isfinite(r)])
    edges = np.linspace(0, rmax, nbins + 1)
    cx, cy, cye = [], [], []
    for i in range(nbins):
        mask = (r >= edges[i]) & (r < edges[i + 1]) & np.isfinite(vals) & (errs > 0)
        if mask.sum() == 0:
            continue
        w = 1.0 / errs[mask] ** 2
        cx.append(0.5 * (edges[i] + edges[i + 1]))
        cy.append(np.sum(w * vals[mask]) / np.sum(w))
        cye.append(1.0 / np.sqrt(np.sum(w)))
    return np.array(cx), np.array(cy), np.array(cye)


def main():
    ap = argparse.ArgumentParser(
        description="Compare 3D MC spatial correlator to 2D analytical modular Ising"
    )
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="compare_2d_modular.png", help="Output figure path")
    ap.add_argument("--mp-dps", type=int, default=40, help="mpmath decimal places")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_all_to_all.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No data files found matching: {pattern}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_abs = axes[0]   # absolute values
    ax_ratio = axes[1] # ratio MC / 2D theory

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(files)))

    for fpath, color in zip(files, colors):
        fpath = Path(fpath)
        geo = parse_geometry(fpath.parent.name)
        if geo is None:
            print(f"Skipping (unrecognized dir): {fpath.parent.name}")
            continue
        lx, ly, tx, ty, nt = geo["Lx"], geo["Ly"], geo["Tx"], geo["Ty"], geo["Nt"]
        s = lx // 39
        label = f"s={s}  Nt={nt}"

        print(f"Processing {fpath.parent.name} ...", flush=True)
        m_arr, n_arr, g_mc, g_err = load_connected(fpath)

        # Triangular lattice physical distance
        r = np.sqrt(m_arr ** 2 + n_arr ** 2 + m_arr * n_arr)

        # 2D CFT computation (at base-s=1 modular parameter — tau is scale-invariant)
        tau, b_sign = modular_tau(lx, ly, tx, ty)
        a, b = torus_coords(m_arr, n_arr, lx, ly, tx, ty, b_sign)
        nu_arr = np.array([complex(float(ai + bi * tau.real), float(bi * tau.imag))
                           for ai, bi in zip(a, b)])

        print(f"  Computing 2D CFT shape for {len(nu_arr)} points ...", flush=True)
        g_theory = ising_torus_shape(nu_arr, tau, mp_dps=args.mp_dps)

        # Fit normalization alpha via weighted least squares: g_mc ≈ alpha * g_theory
        valid = np.isfinite(g_theory) & np.isfinite(g_mc) & (g_err > 0) & (g_theory > 0)
        if valid.sum() < 2:
            print(f"  WARNING: too few valid points for normalization fit; skipping.")
            continue
        w = 1.0 / g_err[valid] ** 2
        alpha = np.sum(w * g_mc[valid] * g_theory[valid]) / np.sum(w * g_theory[valid] ** 2)
        g_theory_scaled = alpha * g_theory

        # Bin by physical distance for cleaner plots
        rc, gc_mc_bin, gc_mc_err = radial_bin_weighted(r, g_mc, g_err)
        _, gc_th_bin, _ = radial_bin_weighted(r[valid], g_theory_scaled[valid],
                                               g_err[valid], nbins=50)

        # Rebin theory on same radius grid as MC
        rc2, gc_th2, _ = radial_bin_weighted(r, g_theory_scaled, np.ones_like(r) * 1e-10)

        ax_abs.errorbar(rc, gc_mc_bin, yerr=gc_mc_err, fmt="o", markersize=3,
                        color=color, label=f"MC {label}", lw=1)
        ax_abs.plot(rc2, gc_th2, "--", color=color, lw=1.2, label=f"2D theory {label}")

        # Ratio MC / 2D theory (binned)
        ratio_r, ratio_g, ratio_e = radial_bin_weighted(
            r[valid],
            g_mc[valid] / g_theory_scaled[valid],
            g_err[valid] / g_theory_scaled[valid]
        )
        ax_ratio.errorbar(ratio_r / float(lx), ratio_g, yerr=ratio_e,
                          fmt="o-", markersize=3, color=color, label=label, lw=1)

    ax_abs.set_yscale("log")
    ax_abs.set_xlabel(r"$r$ (lattice units)")
    ax_abs.set_ylabel(r"$G(r)$")
    ax_abs.set_title("MC vs 2D modular Ising (spatial, equal-time)")
    ax_abs.legend(fontsize=7)

    ax_ratio.axhline(1.0, color="k", lw=0.8, ls="--", label="ratio=1 (perfect 2D)")
    ax_ratio.set_xlabel(r"$r / L_x$")
    ax_ratio.set_ylabel(r"$G_{\rm MC}(r) / G_{\rm 2D\,CFT}(r)$")
    ax_ratio.set_title("Ratio MC / 2D theory (deviation due to 3D)")
    ax_ratio.legend(fontsize=7)

    fig.suptitle("3D Ising on 4-5-6 torus at K=0.161 vs 2D modular Ising CFT", fontsize=11)
    fig.tight_layout()
    outpath = datadir / args.output if not Path(args.output).is_absolute() else Path(args.output)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
