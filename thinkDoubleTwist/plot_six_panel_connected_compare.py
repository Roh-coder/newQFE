#!/usr/bin/env python3
"""Create a six-panel connected-correlator comparison plot.

Layout (2x3):
  Row 1: connected correlator curves
  Row 2: corresponding error curves

Columns:
  1) Untwisted geometry (direct site values, wrapped over full period)
  2) Twisted geometry (interpolated manifold slices along v/u/w)
  3) Difference (twisted - untwisted) with propagated errors

Difference error propagation:
  sigma_diff = sqrt(sigma_twisted^2 + sigma_untwisted^2)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def load_full_rows(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append(
                {
                    "d": int(d_s),
                    "m": int(m_s),
                    "n": int(n_s),
                    "corr": float(c_s),
                    "err": float(e_s),
                    "corr_conn": float(cc_s),
                    "err_conn": float(ec_s),
                }
            )
    if not rows:
        raise SystemExit(f"No data rows found in {path}")
    return rows


def extract_untwisted_direction(rows, dm, dn, period):
    key_to_row = {(r["m"] % period, r["n"] % period): r for r in rows}
    rs, cc, ec = [], [], []
    for r in range(period):
        m = (dm * r) % period
        n = (dn * r) % period
        row = key_to_row.get((m, n))
        if row is not None:
            rs.append(r)
            cc.append(row["corr_conn"])
            ec.append(row["err_conn"])

    if not rs:
        return np.array([]), np.array([]), np.array([])

    # Add N+1 endpoint to close the full period.
    rs.append(period)
    cc.append(cc[0])
    ec.append(ec[0])
    return np.array(rs, dtype=float), np.array(cc, dtype=float), np.array(ec, dtype=float)


def lattice_to_xy(m: np.ndarray, n: np.ndarray):
    return m + 0.5 * n, (np.sqrt(3.0) * 0.5) * n


def load_twisted_arrays(path: Path):
    m, n = [], []
    cc, ec = [], []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            _d_s, m_s, n_s, _c_s, _e_s, cc_s, ec_s = s.split()
            m.append(float(m_s))
            n.append(float(n_s))
            cc.append(float(cc_s))
            ec.append(float(ec_s))
    if not m:
        raise SystemExit(f"No data rows found in {path}")
    return np.array(m), np.array(n), np.array(cc), np.array(ec)


def tile_periodic(m, n, cc, ec, lx, ly, tx, ty):
    m_out, n_out, cc_out, ec_out = [], [], [], []
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            dm = a * lx + b * tx
            dn = a * ty - b * ly
            m_out.append(m + dm)
            n_out.append(n + dn)
            cc_out.append(cc)
            ec_out.append(ec)
    return np.concatenate(m_out), np.concatenate(n_out), np.concatenate(cc_out), np.concatenate(ec_out)


def make_interpolators(m, n, cc, ec):
    x, y = lattice_to_xy(m, n)
    tri = mtri.Triangulation(x, y)
    cc_interp = mtri.LinearTriInterpolator(tri, cc)
    ec_interp = mtri.LinearTriInterpolator(tri, ec)
    return cc_interp, ec_interp


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
    ap.add_argument("--untwisted", required=True, help="Untwisted all-to-all .dat")
    ap.add_argument("--twisted", required=True, help="Twisted all-to-all .dat")
    ap.add_argument("--u_L", type=int, required=True, help="Untwisted period (Lx=Ly)")
    ap.add_argument("--t_Lx", type=int, required=True)
    ap.add_argument("--t_Ly", type=int, required=True)
    ap.add_argument("--t_Tx", type=int, required=True)
    ap.add_argument("--t_Ty", type=int, required=True)
    ap.add_argument("--samples", type=int, default=601, help="Number of points from 0..1 (use odd for midpoint)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    u_rows = load_full_rows(Path(args.untwisted))

    tw_m, tw_n, tw_cc, tw_ec = load_twisted_arrays(Path(args.twisted))
    tw_mt, tw_nt, tw_cct, tw_ect = tile_periodic(tw_m, tw_n, tw_cc, tw_ec, args.t_Lx, args.t_Ly, args.t_Tx, args.t_Ty)
    tw_cc_interp, tw_ec_interp = make_interpolators(tw_mt, tw_nt, tw_cct, tw_ect)

    # Direction correspondence used in earlier comparisons.
    dir_keys = ("d1", "d2", "d3")
    untw_dirs = {"d1": (1, 0), "d2": (0, 1), "d3": (-1, 1)}
    tw_dirs = {
        "d1": (args.t_Lx, args.t_Ty),
        "d2": (args.t_Tx, -args.t_Ly),
        "d3": (-(args.t_Lx + args.t_Tx), args.t_Ly - args.t_Ty),
    }
    untw_labels = {"d1": "e1", "d2": "e2", "d3": "e2-e1"}
    tw_labels = {"d1": "v", "d2": "u", "d3": "w"}
    diff_labels = {"d1": "v-e1", "d2": "u-e2", "d3": "w-(e2-e1)"}
    colors = {"d1": "red", "d2": "blue", "d3": "black"}

    untw_data = {}
    tw_data = {}
    diff_data = {}

    for name in dir_keys:
        dm_u, dn_u = untw_dirs[name]
        rs, cc_u, ec_u = extract_untwisted_direction(u_rows, dm_u, dn_u, args.u_L)
        if len(rs) == 0:
            continue
        t_u = rs / float(args.u_L)
        cc_u_direct = cc_u
        ec_u_direct = ec_u

        dm_t, dn_t = tw_dirs[name]
        # Twisted panels are sampled densely from interpolation.
        t_t, cc_t = sample_slice(tw_cc_interp, dm_t, dn_t, args.samples)
        _t_t, ec_t = sample_slice(tw_ec_interp, dm_t, dn_t, args.samples)
        mask = np.isfinite(cc_t) & np.isfinite(ec_t)
        cc_t_use = np.where(mask, cc_t, np.nan)
        ec_t_use = np.where(mask, ec_t, np.nan)

        # For difference, evaluate twisted interpolant exactly at untwisted site fractions.
        # This keeps untwisted fully site-based with no interpolation.
        t_d = t_u
        m_d = t_d * dm_t
        n_d = t_d * dn_t
        x_d, y_d = lattice_to_xy(m_d, n_d)
        cc_td = tw_cc_interp(x_d, y_d)
        ec_td = tw_ec_interp(x_d, y_d)
        cc_td = np.asarray(cc_td.filled(np.nan) if hasattr(cc_td, "filled") else cc_td, dtype=float)
        ec_td = np.asarray(ec_td.filled(np.nan) if hasattr(ec_td, "filled") else ec_td, dtype=float)

        d = cc_td - cc_u_direct
        ed = np.sqrt(ec_td * ec_td + ec_u_direct * ec_u_direct)

        untw_data[name] = (t_u, cc_u_direct, ec_u_direct)
        tw_data[name] = (t_t, cc_t_use, ec_t_use)
        diff_data[name] = (t_d, d, ed)

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.2), dpi=200, sharex="col")

    u_geom = f"L_x={args.u_L}, L_y={args.u_L}, T_x=0, T_y=0"
    t_geom = f"L_x={args.t_Lx}, L_y={args.t_Ly}, T_x={args.t_Tx}, T_y={args.t_Ty}"

    # Top row: connected correlators.
    for name in dir_keys:
        if name not in untw_data:
            continue

        t, y, e = untw_data[name]
        axes[0, 0].errorbar(
            t,
            y,
            yerr=e,
            fmt="-o",
            ms=2.8,
            lw=1.2,
            capsize=2,
            color=colors[name],
            label=untw_labels[name],
        )

        t, y, e = tw_data[name]
        m = np.isfinite(y) & np.isfinite(e)
        axes[0, 1].plot(t[m], y[m], color=colors[name], lw=1.8, label=tw_labels[name])
        axes[0, 1].fill_between(t[m], y[m] - e[m], y[m] + e[m], color=colors[name], alpha=0.20)

        t, y, e = diff_data[name]
        m = np.isfinite(y) & np.isfinite(e)
        axes[0, 2].plot(t[m], y[m], color=colors[name], lw=1.8, label=diff_labels[name])
        axes[0, 2].fill_between(t[m], y[m] - e[m], y[m] + e[m], color=colors[name], alpha=0.20)

    axes[0, 0].set_title(f"Untwisted Sites ({u_geom})")
    axes[0, 1].set_title(f"Twisted Interp ({t_geom})")
    axes[0, 2].set_title("Diff (T - U)")
    axes[0, 2].axhline(0.0, color="gray", lw=1.0, alpha=0.6)

    axes[0, 0].set_ylabel("C_conn")
    for j in range(3):
        axes[0, j].grid(alpha=0.25)
        axes[0, j].legend(frameon=False, fontsize=9)

    # Bottom row: error magnitudes.
    for name in dir_keys:
        if name not in untw_data:
            continue

        t, _y, e = untw_data[name]
        axes[1, 0].plot(t, e, "-o", ms=2.8, lw=1.2, color=colors[name], label=untw_labels[name])

        t, _y, e = tw_data[name]
        m = np.isfinite(e)
        axes[1, 1].plot(t[m], e[m], color=colors[name], lw=1.8, label=tw_labels[name])

        t, _y, e = diff_data[name]
        m = np.isfinite(e)
        axes[1, 2].plot(t[m], e[m], color=colors[name], lw=1.8, label=diff_labels[name])

    axes[1, 0].set_title("U Error")
    axes[1, 1].set_title("T Error")
    axes[1, 2].set_title("Diff Error")

    axes[1, 0].set_xlabel("fraction of full period/vector (0 to 1)")
    axes[1, 1].set_xlabel("fraction of full period/vector (0 to 1)")
    axes[1, 2].set_xlabel("fraction of full period/vector (0 to 1)")

    axes[1, 0].set_ylabel("error(C_conn)")
    for j in range(3):
        axes[1, j].grid(alpha=0.25)
        axes[1, j].legend(frameon=False, fontsize=9)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
