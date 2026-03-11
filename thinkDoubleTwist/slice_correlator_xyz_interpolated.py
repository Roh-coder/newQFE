#!/usr/bin/env python3
"""Interpolate all-to-all correlator and slice along x, y, and x+y directions.

Input row format:
  d m n corr err corr_conn err_conn

Directions:
  x = (1, 0)
  y = (0, 1)
  xy = (1, 1)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def lattice_to_xy(m: np.ndarray, n: np.ndarray):
    # Triangular basis: e1=(1,0), e2=(1/2,sqrt(3)/2)
    return m + 0.5 * n, (np.sqrt(3.0) * 0.5) * n


def load_full(path: Path):
    m, n = [], []
    c, cc = [], []
    e, ec = [], []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            _d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            m.append(float(m_s))
            n.append(float(n_s))
            c.append(float(c_s))
            e.append(float(e_s))
            cc.append(float(cc_s))
            ec.append(float(ec_s))
    if not m:
        raise SystemExit("No all-to-all rows found.")
    return np.array(m), np.array(n), np.array(c), np.array(e), np.array(cc), np.array(ec)


def tile_periodic(m, n, c, e, cc, ec, lx, ly):
    # Build 3x3 periodic image tiles in (x, y) shifts to avoid boundary artifacts.
    m_out, n_out, c_out, e_out, cc_out, ec_out = [], [], [], [], [], []
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            dm = a * lx
            dn = b * ly
            m_out.append(m + dm)
            n_out.append(n + dn)
            c_out.append(c)
            e_out.append(e)
            cc_out.append(cc)
            ec_out.append(ec)
    return (
        np.concatenate(m_out),
        np.concatenate(n_out),
        np.concatenate(c_out),
        np.concatenate(e_out),
        np.concatenate(cc_out),
        np.concatenate(ec_out),
    )


def make_interpolators(m, n, c, e, cc, ec):
    x, y = lattice_to_xy(m, n)
    tri = mtri.Triangulation(x, y)
    c_interp = mtri.LinearTriInterpolator(tri, c)
    e_interp = mtri.LinearTriInterpolator(tri, e)
    cc_interp = mtri.LinearTriInterpolator(tri, cc)
    ec_interp = mtri.LinearTriInterpolator(tri, ec)
    return c_interp, e_interp, cc_interp, ec_interp


def sample_slice(interp, dm, dn, n_samples):
    t = np.linspace(0.0, 1.0, n_samples)
    m = t * dm
    n = t * dn
    x, y = lattice_to_xy(m, n)
    z = interp(x, y)
    z = np.asarray(z.filled(np.nan) if hasattr(z, "filled") else z)
    return t, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--samples", type=int, default=400, help="Samples per slice")
    ap.add_argument("--output", required=True, help="Output .png path")
    args = ap.parse_args()

    lx, ly = args.L_x, args.L_y

    m, n, c, e, cc, ec = load_full(Path(args.full))
    mt, nt, ct, et, cct, ect = tile_periodic(m, n, c, e, cc, ec, lx, ly)
    c_interp, e_interp, cc_interp, ec_interp = make_interpolators(mt, nt, ct, et, cct, ect)

    dirs = {
        "x": (1, 0),
        "y": (0, 1),
        "xy": (1, 1),
    }
    colors = {"x": "red", "y": "blue", "xy": "green"}

    fig, axes = plt.subplots(2, 1, figsize=(8.3, 7.1), dpi=180, sharex=True)

    for name in ("x", "y", "xy"):
        dm, dn = dirs[name]
        t, z = sample_slice(c_interp, dm, dn, args.samples)
        _, ez = sample_slice(e_interp, dm, dn, args.samples)
        mask = np.isfinite(z)
        axes[0].plot(t[mask], z[mask], color=colors[name], lw=1.8, label=name)
        axes[0].fill_between(t[mask], z[mask] - ez[mask], z[mask] + ez[mask], color=colors[name], alpha=0.22)

    axes[0].set_title("Interpolated Raw Correlator Slices Along x/y/xy")
    axes[0].set_ylabel("C")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    for name in ("x", "y", "xy"):
        dm, dn = dirs[name]
        t, z = sample_slice(cc_interp, dm, dn, args.samples)
        _, ez = sample_slice(ec_interp, dm, dn, args.samples)
        mask = np.isfinite(z)
        axes[1].plot(t[mask], z[mask], color=colors[name], lw=1.8, label=name)
        axes[1].fill_between(t[mask], z[mask] - ez[mask], z[mask] + ez[mask], color=colors[name], alpha=0.22)

    axes[1].set_title("Interpolated Connected Correlator Slices Along x/y/xy")
    axes[1].set_xlabel("fraction of full period (0 to 1)")
    axes[1].set_ylabel("C_conn")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
