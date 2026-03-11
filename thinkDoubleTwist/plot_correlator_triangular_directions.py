#!/usr/bin/env python3
"""Plot correlator along the three base triangular directions for untwisted lattice.

Input row format:
  d m n corr err corr_conn err_conn

Directions:
  e1 = (1, 0)   (base)
  e2 = (0, 1)   (left)
  e2-e1 = (-1, 1) (right)

Assumes L_x = L_y = period, T_x = T_y = 0.
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
    # Map (m, n) to row for fast lookup
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L", type=int, required=True, help="Lattice period (L_x = L_y)")
    ap.add_argument("--output", required=True, help="Output .png path")
    args = ap.parse_args()

    period = args.L
    rows = load_full(Path(args.full))
    dirs = {
        "e1 (base)": (1, 0),
        "e2 (left)": (0, 1),
        "e2-e1 (right)": (-1, 1),
    }
    colors = {"e1 (base)": "red", "e2 (left)": "blue", "e2-e1 (right)": "green"}

    fig, axes = plt.subplots(2, 1, figsize=(8.3, 7.1), dpi=180, sharex=True)

    for name, (dm, dn) in dirs.items():
        rs, c, e, cc, ec = extract_direction(rows, dm, dn, period)
        axes[0].errorbar(rs, c, yerr=e, fmt="-o", ms=2.8, lw=1.1, color=colors[name], capsize=2, label=name)
    axes[0].set_title("Raw Correlator Along Triangular Directions")
    axes[0].set_ylabel("C(r)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    for name, (dm, dn) in dirs.items():
        rs, c, e, cc, ec = extract_direction(rows, dm, dn, period)
        axes[1].errorbar(rs, cc, yerr=ec, fmt="-o", ms=2.8, lw=1.1, color=colors[name], capsize=2, label=name)
    axes[1].set_title("Connected Correlator Along Triangular Directions")
    axes[1].set_xlabel("step r along direction")
    axes[1].set_ylabel("C_conn(r)")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)

if __name__ == "__main__":
    main()
