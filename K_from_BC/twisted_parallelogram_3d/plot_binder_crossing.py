#!/usr/bin/env python3
"""
plot_binder_crossing.py  (ribbon edition)
=========================================
Read out_binder_3d/summary.tsv (written live by run_binder_scan.py) or fall
back to scanning .dat files, then produce a Binder-cumulant crossing plot
with shaded ±1σ ribbons per lattice size L.

Safe to run while run_binder_scan.py is still going — just re-run to refresh.

Usage:
    python plot_binder_crossing.py               # use summary.tsv
    python plot_binder_crossing.py --rescan      # re-parse .dat files
    python plot_binder_crossing.py --data-dir X
    python plot_binder_crossing.py --out foo.png
"""

import os
import sys
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
SCRIPT_DIR       = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.join(SCRIPT_DIR, "out_binder_3d")
KC_ESTIMATE      = 0.165
BINDER_3D_ISING  = 0.4655     # universal Binder value for 3D Ising, periodic BC
K_PLOT_MIN       = 0.16
K_PLOT_MAX       = 0.162

COLORS  = {8: "#1f77b4", 16: "#ff7f0e", 24: "#2ca02c", 32: "#d62728"}
MARKERS = {8: "o",        16: "s",        24: "^",        32: "D"}
# ---------------------------------------------------------------------------


def read_summary(path):
    """Read summary.tsv → {L: sorted [(K, U4, U4_err), ...]}."""
    from collections import defaultdict
    records = defaultdict(list)
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("L\t") or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                try:
                    L   = int(parts[0])
                    K   = float(parts[1])
                    u4  = float(parts[2])
                    u4e = float(parts[3])
                    records[L].append((K, u4, u4e))
                except ValueError:
                    pass
    except OSError:
        pass
    for L in records:
        # deduplicate: keep last entry per K
        seen = {}
        for item in records[L]:
            seen[item[0]] = item
        records[L] = sorted(seen.values())
    return dict(records)


def scan_dat_files(data_dir):
    """Fallback: scan all .dat files under data_dir."""
    from collections import defaultdict
    records = defaultdict(list)
    pat = os.path.join(data_dir, "*", "*.dat")
    for path in glob.glob(pat):
        d = {}
        try:
            with open(path) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            d[parts[0]] = [float(x) for x in parts[1:]]
                        except ValueError:
                            pass
        except OSError:
            continue
        try:
            L   = int(d["L_x"][0])
            K   = float(d["K1"][0])
            u4  = float(d["U4"][0])
            u4e = float(d["U4"][1])
        except (KeyError, IndexError):
            continue
        records[L].append((K, u4, u4e))
    for L in records:
        records[L].sort()
    return dict(records)


def find_crossings(K1, U1, K2, U2):
    """Pairwise crossing locations via dense linear interpolation."""
    K_lo = max(K1[0], K2[0])
    K_hi = min(K1[-1], K2[-1])
    if K_lo >= K_hi:
        return []
    K_shared = np.linspace(K_lo, K_hi, 3000)
    U1i = np.interp(K_shared, K1, U1)
    U2i = np.interp(K_shared, K2, U2)
    diff = U1i - U2i
    crossings = []
    for i in range(len(diff) - 1):
        if diff[i] * diff[i + 1] <= 0 and diff[i] != diff[i + 1]:
            t  = -diff[i] / (diff[i + 1] - diff[i])
            Kc = K_shared[i] + t * (K_shared[i + 1] - K_shared[i])
            Uc = float(0.5 * (
                np.interp(Kc, K1, U1) +
                np.interp(Kc, K2, U2)
            ))
            crossings.append((float(Kc), Uc))
    return crossings


def make_plot(records, out_path):
    sizes = sorted(records.keys())
    if not sizes:
        print("No data to plot.")
        return

    K_fine = np.linspace(0.065, 0.265, 1200)

    fig, ax = plt.subplots(figsize=(9, 6))

    # --- draw scatter points with error bars only ---
    for L in sizes:
        pts    = records[L]
        K_arr  = np.array([p[0] for p in pts])
        U4_arr = np.array([p[1] for p in pts])
        U4_err = np.abs(np.array([p[2] for p in pts]))
        c = COLORS.get(L, f"C{sizes.index(L)}")
        m = MARKERS.get(L, "o")

        # Keep data sorted by K so lines/markers are plotted at the right x positions.
        order = np.argsort(K_arr)
        K_arr = K_arr[order]
        U4_arr = U4_arr[order]
        U4_err = U4_err[order]

        # Plot measured points with 1-sigma error bars; no ribbons/lines.
        ax.errorbar(
            K_arr,
            U4_arr,
            yerr=U4_err,
            fmt=m,
            color=c,
            ecolor=c,
            markersize=5,
            markeredgewidth=1.1,
            elinewidth=0.9,
            capsize=2,
            linestyle="none",
            zorder=2,
            label=f"L={L}",
        )

        # Add linear interpolation curve through measured points.
        if len(K_arr) >= 2:
            k_line = np.linspace(float(K_arr[0]), float(K_arr[-1]), 200)
            u4_line = np.interp(k_line, K_arr, U4_arr)
            err_line = np.interp(k_line, K_arr, U4_err)
            ax.plot(k_line, u4_line, color=c, linewidth=1.1, alpha=0.9, zorder=1)
            ax.fill_between(
                k_line,
                u4_line - err_line,
                u4_line + err_line,
                color=c,
                alpha=0.14,
                zorder=0.8,
            )

    ax.set_xscale("linear")
    ax.set_xlim(left=K_PLOT_MIN, right=K_PLOT_MAX)
    ax.set_xticks(np.linspace(K_PLOT_MIN, K_PLOT_MAX, 5))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(r"$K = K_1 = K_2 = K_3 = K_t$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title("Binder cumulant scatter with error bars", fontsize=11)
    ax.set_ylim(-0.05, 1.10)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.12)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    ap.add_argument("--out", default=None,
                    help="output PNG (default: <data-dir>/binder_crossing.png)")
    ap.add_argument("--rescan", action="store_true",
                    help="ignore summary.tsv, re-parse all .dat files")
    args = ap.parse_args()

    data_dir     = args.data_dir
    out_path     = args.out or os.path.join(data_dir, "binder_crossing.png")
    summary_path = os.path.join(data_dir, "summary.tsv")

    if not args.rescan and os.path.exists(summary_path):
        records = read_summary(summary_path)
        total = sum(len(v) for v in records.values())
        print(f"Read {total} point(s) from {summary_path}")
    else:
        records = scan_dat_files(data_dir)
        total = sum(len(v) for v in records.values())
        print(f"Scanned {total} point(s) from .dat files in {data_dir}")

    if not records:
        print("No data yet. Run run_binder_scan.py first.")
        sys.exit(0)

    make_plot(records, out_path)
    print("Done.")


if __name__ == "__main__":
    main()
