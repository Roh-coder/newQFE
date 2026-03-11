#!/usr/bin/env python3
"""Plot correlator along v, u, w directions from full all-to-all output.

Input rows format:
  d m n corr err corr_conn err_conn

Directions are defined in lattice coordinates:
  v = (L_x, T_y)
  u = (T_x, -L_y)
  w = -(u + v) = (-(L_x + T_x), L_y - T_y)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def mod(a: int, n: int) -> int:
    r = a % n
    return r + n if r < 0 else r


def make_key(k1: int, k2: int) -> int:
    return ((k1 & 0xFFFFFFFF) << 32) | (k2 & 0xFFFFFFFF)


def coset_key(m: int, n: int, lx: int, ly: int, tx: int, ty: int, ncell: int) -> int:
    k1 = mod(-ly * m - tx * n, ncell)
    k2 = mod(-ty * m + lx * n, ncell)
    return make_key(k1, k2)


def load_full(path: Path):
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
    return rows


def period_for_direction(dm: int, dn: int, lx: int, ly: int, tx: int, ty: int, ncell: int) -> int:
    k0 = coset_key(0, 0, lx, ly, tx, ty, ncell)
    r = 1
    while True:
        k = coset_key(r * dm, r * dn, lx, ly, tx, ty, ncell)
        if k == k0:
            return r
        if r > ncell:
            raise RuntimeError("Could not find period within ncell steps")
        r += 1


def primitive_dir(dm: int, dn: int):
    g = abs(__import__("math").gcd(dm, dn))
    if g == 0:
        return dm, dn, 1
    return dm // g, dn // g, g


def extract_line(dm: int, dn: int, key_to_row: dict, lx: int, ly: int, tx: int, ty: int, ncell: int):
    p = period_for_direction(dm, dn, lx, ly, tx, ty, ncell)
    rs = []
    corr = []
    err = []
    corr_conn = []
    err_conn = []

    for r in range(p):
        k = coset_key(r * dm, r * dn, lx, ly, tx, ty, ncell)
        row = key_to_row[k]
        rs.append(r)
        corr.append(row["corr"])
        err.append(row["err"])
        corr_conn.append(row["corr_conn"])
        err_conn.append(row["err_conn"])

    return rs, corr, err, corr_conn, err_conn, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    ap.add_argument("--output", required=True, help="Output image path")
    args = ap.parse_args()

    lx, ly, tx, ty = args.L_x, args.L_y, args.T_x, args.T_y
    ncell = lx * ly + tx * ty
    if ncell <= 0:
        raise SystemExit("Need L_x*L_y + T_x*T_y > 0")

    rows = load_full(Path(args.full))
    if not rows:
        raise SystemExit("No data rows found")

    key_to_row = {}
    for row in rows:
        k = coset_key(row["m"], row["n"], lx, ly, tx, ty, ncell)
        key_to_row[k] = row

    dirs = {
        "v": (lx, ty),
        "u": (tx, -ly),
        "w": (-(lx + tx), ly - ty),
    }

    extracted = {}
    scale = {}
    for name, (dm, dn) in dirs.items():
        dmp, dnp, g = primitive_dir(dm, dn)
        scale[name] = g
        extracted[name] = extract_line(dmp, dnp, key_to_row, lx, ly, tx, ty, ncell)

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.2), dpi=180, sharex=False)

    colors = {"v": "red", "u": "blue", "w": "black"}

    ax = axes[0]
    for name in ("v", "u", "w"):
        rs, corr, err, _cc, _ec, p = extracted[name]
        ax.errorbar(rs, corr, yerr=err, fmt="-o", ms=2.8, lw=1.1,
                    color=colors[name], capsize=2,
                    label=f"{name} primitive step (full={scale[name]}x, period {p})")
    ax.set_title("Raw Correlator Along v/u/w")
    ax.set_ylabel("C(r)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[1]
    for name in ("v", "u", "w"):
        rs, _c, _e, cc, ec, p = extracted[name]
        ax.errorbar(rs, cc, yerr=ec, fmt="-o", ms=2.8, lw=1.1,
                    color=colors[name], capsize=2,
                    label=f"{name} primitive step (full={scale[name]}x, period {p})")
    ax.set_title("Connected Correlator Along v/u/w")
    ax.set_xlabel("step r along direction")
    ax.set_ylabel("C_conn(r)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
