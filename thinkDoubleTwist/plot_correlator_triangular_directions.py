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
    # Map wrapped (m, n) to row for periodic lookup on untwisted torus.
    key_to_row = {(r["m"] % period, r["n"] % period): r for r in rows}
    rs, corr, err, corr_conn, err_conn = [], [], [], [], []
    for r in range(period):
        m = (dm * r) % period
        n = (dn * r) % period
        row = key_to_row.get((m, n))
        if row is not None:
            rs.append(r)
            corr.append(row["corr"])
            err.append(row["err"])
            corr_conn.append(row["corr_conn"])
            err_conn.append(row["err_conn"])
    return np.array(rs), np.array(corr), np.array(err), np.array(corr_conn), np.array(err_conn)


def append_period_endpoint(rs, ys, es, period):
    # Add the N+1th point (r=period) by repeating r=0 to close one full orbit.
    if len(rs) == 0:
        return rs, ys, es
    if rs[-1] == period:
        return rs, ys, es
    return (
        np.append(rs, period),
        np.append(ys, ys[0]),
        np.append(es, es[0]),
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L", type=int, required=True, help="Lattice period (L_x = L_y)")
    ap.add_argument("--output", required=True, help="Output .png path")
    ap.add_argument(
        "--connected-only",
        action="store_true",
        help="Plot only connected correlator with a lower panel showing direct errors",
    )
    args = ap.parse_args()

    period = args.L
    rows = load_full(Path(args.full))
    dirs = {
        "e1 (base)": (1, 0),
        "e2 (left)": (0, 1),
        "e2-e1 (right)": (-1, 1),
    }
    colors = {"e1 (base)": "red", "e2 (left)": "blue", "e2-e1 (right)": "green"}
    geom = f"(L_x={period}, L_y={period}, T_x=0, T_y=0)"

    fig, axes = plt.subplots(2, 1, figsize=(8.3, 7.1), dpi=180, sharex=True)

    if args.connected_only:
        for name, (dm, dn) in dirs.items():
            rs, _c, _e, cc, ec = extract_direction(rows, dm, dn, period)
            rs, cc, ec = append_period_endpoint(rs, cc, ec, period)
            axes[0].errorbar(rs, cc, yerr=ec, fmt="-o", ms=2.8, lw=1.2, color=colors[name], capsize=2, label=name)
            axes[1].plot(rs, ec, "-o", ms=2.8, lw=1.2, color=colors[name], label=name)

        axes[0].set_title(f"Connected Correlator Along Triangular Directions (Direct Site Values) {geom}")
        axes[0].set_ylabel("C_conn(r)")
        axes[0].grid(alpha=0.25)
        axes[0].legend(frameon=False)

        axes[1].set_title(f"Direct Error Along Full Period in Triangular Directions {geom}")
        axes[1].set_xlabel("step r along direction")
        axes[1].set_ylabel("error(C_conn)")
        axes[1].grid(alpha=0.25)
        axes[1].legend(frameon=False)
    else:
        for name, (dm, dn) in dirs.items():
            rs, c, e, cc, ec = extract_direction(rows, dm, dn, period)
            rs, c, e = append_period_endpoint(rs, c, e, period)
            axes[0].errorbar(rs, c, yerr=e, fmt="-o", ms=2.8, lw=1.1, color=colors[name], capsize=2, label=name)
        axes[0].set_title(f"Raw Correlator Along Triangular Directions {geom}")
        axes[0].set_ylabel("C(r)")
        axes[0].grid(alpha=0.25)
        axes[0].legend(frameon=False)

        for name, (dm, dn) in dirs.items():
            rs, c, e, cc, ec = extract_direction(rows, dm, dn, period)
            rs, cc, ec = append_period_endpoint(rs, cc, ec, period)
            axes[1].errorbar(rs, cc, yerr=ec, fmt="-o", ms=2.8, lw=1.1, color=colors[name], capsize=2, label=name)
        axes[1].set_title(f"Connected Correlator Along Triangular Directions {geom}")
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
