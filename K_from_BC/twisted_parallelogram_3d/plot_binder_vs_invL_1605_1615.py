#!/usr/bin/env python3
"""
plot_binder_vs_invL_1605_1615.py
=================================
Plot Binder cumulant U4 vs 1/L from the scan
  out_binder_1605_1615_nk15_ntraj500_Nt4L_L4to64/summary.tsv

Can be run while the scan is still in progress; it will plot whatever
rows exist so far.

Usage:
    python plot_binder_vs_invL_1605_1615.py [--data-dir DIR]
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.join(
    SCRIPT_DIR,
    "out_binder_1605_1615_nk15_ntraj500_Nt4L_L4to64",
)


def read_summary(path):
    """
    Parse summary.tsv.
    Deduplicates by (L, dat_file), then inverse-variance-averages
    any multiple independent runs at the same (L, K).
    Returns {L: [(K, U4, U4_err), ...]} sorted by K.
    """
    seen = {}
    with open(path) as fh:
        for raw in fh:
            raw = raw.strip()
            if not raw or raw.startswith("L\t") or raw.startswith("#"):
                continue
            p = raw.split("\t")
            if len(p) < 7:
                continue
            try:
                L   = int(p[0])
                K   = float(p[1])
                u4  = float(p[2])
                u4e = float(p[3])
                dat = p[6]
            except (ValueError, IndexError):
                continue
            key = (L, dat)
            if key not in seen:
                seen[key] = (K, u4, u4e)

    by_LK = defaultdict(list)
    for (L, _dat), (K, u4, u4e) in seen.items():
        by_LK[(L, K)].append((u4, u4e))

    records = defaultdict(list)
    for (L, K), pts in by_LK.items():
        if len(pts) == 1:
            u4, u4e = pts[0]
        else:
            weights = np.array([1.0 / e**2 for _, e in pts if e > 0])
            vals    = np.array([v           for v, _ in pts])
            w_total = weights.sum()
            u4      = float((weights * vals).sum() / w_total)
            u4e     = float(1.0 / np.sqrt(w_total))
        records[L].append((K, u4, u4e))

    for L in records:
        records[L].sort()
    return dict(records)


def plot_vs_inv_L(records, out_path):
    """U4 vs 1/L, one line per K value."""
    by_K = defaultdict(list)
    for L, pts in records.items():
        for K, u4, u4e in pts:
            by_K[K].append((L, u4, u4e))
    for K in by_K:
        by_K[K].sort()

    K_vals = sorted(by_K.keys())
    n_K    = len(K_vals)

    cmap   = plt.cm.turbo
    colors = [cmap(i / max(n_K - 1, 1)) for i in range(n_K)]

    fig, ax = plt.subplots(figsize=(13, 6))

    for K, color in zip(K_vals, colors):
        pts   = by_K[K]
        Larr  = np.array([p[0] for p in pts], dtype=float)
        U4arr = np.array([p[1] for p in pts])
        U4err = np.abs(np.array([p[2] for p in pts]))
        x     = 1.0 / Larr

        ax.errorbar(x, U4arr, yerr=U4err,
                    fmt="o", color=color, ecolor=color,
                    markersize=4, markeredgewidth=0.8,
                    elinewidth=0.9, capsize=2,
                    linestyle="-", linewidth=1.1,
                    label=f"K={K:.7f}",
                    zorder=2)

    Ls      = sorted(records.keys())
    x_ticks = [1.0 / L for L in Ls]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"1/{L}\n(L={L})" for L in Ls], fontsize=8)

    ax.legend(fontsize=7, loc="upper left",
              ncol=2, framealpha=0.9, title=r"$K$", title_fontsize=8)

    ax.set_xlabel(r"$1/L$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title(
        r"$U_4$ vs $1/L$ — equilateral 3D Ising ($N_t=4L$, $T_x=T_y=0$, "
        r"$K\in[0.1605,0.1615]$)",
        fontsize=11)
    ax.set_xlim(left=0.0, right=max(x_ticks) * 1.18)
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DIR,
                    help="directory containing summary.tsv")
    args = ap.parse_args()

    tsv_path = os.path.join(args.data_dir, "summary.tsv")
    if not os.path.exists(tsv_path):
        print(f"ERROR: {tsv_path} not found.", file=sys.stderr)
        sys.exit(1)

    records = read_summary(tsv_path)
    total   = sum(len(v) for v in records.values())
    print(f"Loaded {total} unique (L, K) points from {tsv_path}")
    for L in sorted(records):
        print(f"  L={L}: {len(records[L])} K points, "
              f"K=[{records[L][0][0]:.7f}, {records[L][-1][0]:.7f}]")

    out_path = os.path.join(args.data_dir, "binder_vs_inv_L.png")
    plot_vs_inv_L(records, out_path)
    print("Done.")


if __name__ == "__main__":
    main()
