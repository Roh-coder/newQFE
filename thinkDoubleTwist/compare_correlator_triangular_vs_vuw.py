#!/usr/bin/env python3
"""Compare correlators along triangular directions (untwisted) and v/u/w (twisted).

Inputs:
  --untwisted: all-to-all .dat for untwisted (Lx=Ly, Tx=Ty=0)
  --twisted: all-to-all .dat for twisted (Lx,Ly,Tx,Ty)
  --L: period for untwisted
  --L_x, --L_y, --T_x, --T_y: for twisted
  --output: output .png

Plots both raw and connected correlators with error bars.
"""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def load_full(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append({
                "d": int(d_s),
                "m": int(m_s),
                "n": int(n_s),
                "corr": float(c_s),
                "err": float(e_s),
                "corr_conn": float(cc_s),
                "err_conn": float(ec_s),
            })
    return rows

def extract_direction(rows, dm, dn, period):
    key_to_row = {(r["m"], r["n"]): r for r in rows}
    rs, corr, err, corr_conn, err_conn = [], [], [], [], []
    for r in range(period):
        m = dm * r
        n = dn * r
        row = key_to_row.get((m, n))
        if row is not None:
            rs.append(r)
            corr.append(row["corr"])
            err.append(row["err"])
            corr_conn.append(row["corr_conn"])
            err_conn.append(row["err_conn"])
    return np.array(rs), np.array(corr), np.array(err), np.array(corr_conn), np.array(err_conn)

def period_for_direction(dm, dn, lx, ly, tx, ty):
    # Find period for direction in twisted cell
    ncell = lx * ly + tx * ty
    def coset_key(m, n):
        def mod(a, n):
            r = a % n
            return r + n if r < 0 else r
        k1 = mod(-ly * m - tx * n, ncell)
        k2 = mod(-ty * m + lx * n, ncell)
        return (k1, k2)
    k0 = coset_key(0, 0)
    r = 1
    while True:
        k = coset_key(dm * r, dn * r)
        if k == k0:
            return r
        if r > ncell:
            raise RuntimeError("Could not find period within ncell steps")
        r += 1

def extract_direction_twisted(rows, dm, dn, lx, ly, tx, ty):
    ncell = lx * ly + tx * ty
    def coset_key(m, n):
        def mod(a, n):
            r = a % n
            return r + n if r < 0 else r
        k1 = mod(-ly * m - tx * n, ncell)
        k2 = mod(-ty * m + lx * n, ncell)
        return (k1, k2)
    key_to_row = {}
    for r in rows:
        key_to_row[coset_key(r["m"], r["n"])] = r
    period = period_for_direction(dm, dn, lx, ly, tx, ty)
    rs, corr, err, corr_conn, err_conn = [], [], [], [], []
    for r in range(period):
        m = dm * r
        n = dn * r
        row = key_to_row.get(coset_key(m, n))
        if row is not None:
            rs.append(r)
            corr.append(row["corr"])
            err.append(row["err"])
            corr_conn.append(row["corr_conn"])
            err_conn.append(row["err_conn"])
    return np.array(rs), np.array(corr), np.array(err), np.array(corr_conn), np.array(err_conn)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--untwisted", required=True)
    ap.add_argument("--twisted", required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    # Untwisted directions
    untwisted_rows = load_full(Path(args.untwisted))
    tri_dirs = {
        "e1 (base)": (1, 0),
        "e2 (left)": (0, 1),
        "e2-e1 (right)": (-1, 1),
    }
    tri_colors = {"e1 (base)": "red", "e2 (left)": "blue", "e2-e1 (right)": "green"}

    # Twisted v/u/w
    twisted_rows = load_full(Path(args.twisted))
    v = (args.L_x, args.T_y)
    u = (args.T_x, -args.L_y)
    w = (-(args.L_x + args.T_x), args.L_y - args.T_y)
    vuw_dirs = {"v": v, "u": u, "w": w}
    vuw_colors = {"v": "orange", "u": "purple", "w": "black"}

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=180, sharex=False)

    # Raw correlators
    ax = axes[0]
    for name, (dm, dn) in tri_dirs.items():
        rs, c, e, cc, ec = extract_direction(untwisted_rows, dm, dn, args.L)
        ax.errorbar(rs, c, yerr=e, fmt="-o", ms=2.5, lw=1.1, color=tri_colors[name], capsize=2, label=f"untwisted {name}")
    for name, (dm, dn) in vuw_dirs.items():
        rs, c, e, cc, ec = extract_direction_twisted(twisted_rows, dm, dn, args.L_x, args.L_y, args.T_x, args.T_y)
        ax.errorbar(rs, c, yerr=e, fmt="--s", ms=2.5, lw=1.1, color=vuw_colors[name], capsize=2, label=f"twisted {name}")
    ax.set_title("Raw Correlator: Triangular vs v/u/w Directions")
    ax.set_ylabel("C(r)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)

    # Connected correlators
    ax = axes[1]
    for name, (dm, dn) in tri_dirs.items():
        rs, c, e, cc, ec = extract_direction(untwisted_rows, dm, dn, args.L)
        ax.errorbar(rs, cc, yerr=ec, fmt="-o", ms=2.5, lw=1.1, color=tri_colors[name], capsize=2, label=f"untwisted {name}")
    for name, (dm, dn) in vuw_dirs.items():
        rs, c, e, cc, ec = extract_direction_twisted(twisted_rows, dm, dn, args.L_x, args.L_y, args.T_x, args.T_y)
        ax.errorbar(rs, cc, yerr=ec, fmt="--s", ms=2.5, lw=1.1, color=vuw_colors[name], capsize=2, label=f"twisted {name}")
    ax.set_title("Connected Correlator: Triangular vs v/u/w Directions")
    ax.set_xlabel("step r along direction")
    ax.set_ylabel("C_conn(r)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)

if __name__ == "__main__":
    main()
