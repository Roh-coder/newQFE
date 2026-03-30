#!/usr/bin/env python3
"""
plot_binder_crossing_new.py
===========================
Produce Binder-cumulant crossing plots from summary.tsv.

Handles the aliasing issue: deduplicates rows by (L, dat_file) so each
independent Monte Carlo run contributes exactly one point.
If multiple independent runs share the same K and L, their U4 values are
combined via inverse-variance weighted average.

Outputs two PNGs next to the data directory:
    binder_crossing_all.png   — all four L values with ribbon bands
    binder_crossing_zoom.png  — zoomed into the L=16..32 crossing region

Usage:
    python plot_binder_crossing_new.py [--data-dir DIR]
"""

import os
import sys
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
import matplotlib.cm

# ── aesthetics ──────────────────────────────────────────────────────────────
COLORS  = {8: "#1f77b4", 16: "#ff7f0e", 24: "#2ca02c", 32: "#d62728"}
MARKERS = {8: "o",        16: "s",        24: "^",        32: "D"}
BINDER_3D_ISING = 0.4655   # 3D Ising universal crossing value (periodic BC)


# ── I/O ─────────────────────────────────────────────────────────────────────

def read_summary(path):
    """
    Parse summary.tsv → {L: [(K, U4, U4_err), ...]} sorted by K.
    Deduplicates by (L, dat_file), then inverse-variance-averages
    any multiple runs at the same (L, K).
    """
    # collect unique (L, dat) → (K, U4, U4_err)
    seen = {}           # (L, dat) -> (K, U4, U4_err)
    with open(path) as fh:
        for raw in fh:
            raw = raw.strip()
            if not raw or raw.startswith("L\t") or raw.startswith("#"):
                continue
            p = raw.split("\t")
            if len(p) < 7:
                continue
            try:
                L    = int(p[0])
                K    = float(p[1])
                u4   = float(p[2])
                u4e  = float(p[3])
                dat  = p[6]
            except (ValueError, IndexError):
                continue
            key = (L, dat)
            if key not in seen:
                seen[key] = (K, u4, u4e)

    # group by (L, K) and inverse-variance average
    by_LK = defaultdict(list)
    for (L, _dat), (K, u4, u4e) in seen.items():
        by_LK[(L, K)].append((u4, u4e))

    records = defaultdict(list)
    for (L, K), pts in by_LK.items():
        if len(pts) == 1:
            u4, u4e = pts[0]
        else:
            # inverse-variance weighted average
            weights = np.array([1.0 / e**2 for _, e in pts if e > 0])
            vals    = np.array([v           for v, _ in pts])
            w_total = weights.sum()
            u4      = float((weights * vals).sum() / w_total)
            u4e     = float(1.0 / np.sqrt(w_total))
        records[L].append((K, u4, u4e))

    for L in records:
        records[L].sort()
    return dict(records)


# ── crossing finder ──────────────────────────────────────────────────────────

def find_crossings(K1, U1, K2, U2):
    """Linear-interpolation pairwise crossing between two sorted (K, U) curves."""
    K_lo = max(K1[0], K2[0])
    K_hi = min(K1[-1], K2[-1])
    if K_lo >= K_hi:
        return []
    K_shared = np.linspace(K_lo, K_hi, 5000)
    U1i = np.interp(K_shared, K1, U1)
    U2i = np.interp(K_shared, K2, U2)
    diff = U1i - U2i
    crossings = []
    for i in range(len(diff) - 1):
        if diff[i] * diff[i + 1] <= 0 and diff[i] != diff[i + 1]:
            t  = -diff[i] / (diff[i + 1] - diff[i])
            Kc = K_shared[i] + t * (K_shared[i + 1] - K_shared[i])
            Uc = 0.5 * (np.interp(Kc, K1, U1) + np.interp(Kc, K2, U2))
            crossings.append((float(Kc), float(Uc)))
    return crossings


# ── plot helpers ─────────────────────────────────────────────────────────────

def draw_curves(ax, records, sizes=None, label_prefix="L="):
    """Draw error-bar scatter + interpolation ribbon for each L."""
    if sizes is None:
        sizes = sorted(records.keys())
    for L in sizes:
        if L not in records:
            continue
        pts    = records[L]
        Karr   = np.array([p[0] for p in pts])
        U4arr  = np.array([p[1] for p in pts])
        U4err  = np.abs(np.array([p[2] for p in pts]))
        order  = np.argsort(Karr)
        Karr, U4arr, U4err = Karr[order], U4arr[order], U4err[order]

        c = COLORS.get(L, f"C{sizes.index(L)}")
        m = MARKERS.get(L, "o")

        ax.errorbar(Karr, U4arr, yerr=U4err,
                    fmt=m, color=c, ecolor=c,
                    markersize=5, markeredgewidth=1.1,
                    elinewidth=0.9, capsize=2,
                    linestyle="none", zorder=3,
                    label=f"{label_prefix}{L}")

        if len(Karr) >= 2:
            Kline = np.linspace(Karr[0], Karr[-1], 400)
            Uline = np.interp(Kline, Karr, U4arr)
            Eline = np.interp(Kline, Karr, U4err)
            ax.plot(Kline, Uline, color=c, linewidth=1.3, alpha=0.9, zorder=2)
            ax.fill_between(Kline, Uline - Eline, Uline + Eline,
                            color=c, alpha=0.15, zorder=1)


def annotate_crossings(ax, records, pairs):
    """Find and annotate crossings between pairs of L values."""
    crossing_Ks = []
    # Collect all crossings first, then stagger annotations vertically
    all_crosses = []
    for La, Lb in pairs:
        if La not in records or Lb not in records:
            continue
        pts_a = records[La]; pts_b = records[Lb]
        Ka = np.array([p[0] for p in pts_a])
        Ua = np.array([p[1] for p in pts_a])
        Kb = np.array([p[0] for p in pts_b])
        Ub = np.array([p[1] for p in pts_b])
        order_a = np.argsort(Ka); Ka, Ua = Ka[order_a], Ua[order_a]
        order_b = np.argsort(Kb); Kb, Ub = Kb[order_b], Ub[order_b]
        crosses = find_crossings(Ka, Ua, Kb, Ub)
        for Kc, Uc in crosses:
            crossing_Ks.append(Kc)
            all_crosses.append((Kc, Uc, La, Lb))

    # Draw vertical lines
    for Kc, Uc, La, Lb in all_crosses:
        ax.axvline(Kc, color=COLORS.get(Lb, "gray"), linestyle="--",
                   linewidth=0.9, alpha=0.6, zorder=0)

    # Stagger text annotations along the y-axis to avoid overlap
    ylim = ax.get_ylim()
    yrange = ylim[1] - ylim[0]
    y_text_positions = np.linspace(ylim[0] + 0.08 * yrange,
                                   ylim[0] + 0.30 * yrange,
                                   max(1, len(all_crosses)))
    for i, (Kc, Uc, La, Lb) in enumerate(all_crosses):
        ax.text(Kc + 0.00008, y_text_positions[i],
                f"$K_c={Kc:.5f}$\n(L={La},{Lb})",
                fontsize=7.5, color=COLORS.get(Lb, "gray"),
                va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          alpha=0.75, edgecolor="none"))

    return crossing_Ks


# ── main plots ───────────────────────────────────────────────────────────────

def plot_all(records, out_path):
    """Full-range plot with all L values — data, error bars and ribbons only."""
    sizes = sorted(records.keys())
    fig, ax = plt.subplots(figsize=(9, 6))

    draw_curves(ax, records, sizes)

    ax.set_xlabel(r"$K = K_1 = K_2 = K_3 = K_t$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title("Binder cumulant crossing — equilateral 3D ($N_t = 4L$)", fontsize=11)
    ax.set_xlim(0.1598, 0.1622)
    ax.set_ylim(-0.05, 1.05)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def plot_zoom(records, out_path):
    """
    Zoom into the L=16..32 crossing region — data, error bars and ribbons only.
    """
    sizes   = sorted(records.keys())
    zoom_Ls = [L for L in sizes if L >= 16]

    # Find approximate crossing range from L=24,32
    Kc_estimates = []
    for La, Lb in [(24, 32), (16, 32), (16, 24)]:
        if La not in records or Lb not in records:
            continue
        pts_a = records[La]; pts_b = records[Lb]
        Ka = np.array([p[0] for p in pts_a]); Ua = np.array([p[1] for p in pts_a])
        Kb = np.array([p[0] for p in pts_b]); Ub = np.array([p[1] for p in pts_b])
        Ka, Ua = Ka[np.argsort(Ka)], Ua[np.argsort(Ka)]
        Kb, Ub = Kb[np.argsort(Kb)], Ub[np.argsort(Kb)]
        for Kc, _ in find_crossings(Ka, Ua, Kb, Ub):
            Kc_estimates.append(Kc)

    if Kc_estimates:
        Kc_mean = np.mean(Kc_estimates)
        half    = 0.0025
        Klo, Khi = Kc_mean - half, Kc_mean + half
    else:
        Klo, Khi = 0.1605, 0.1620

    fig, ax = plt.subplots(figsize=(9, 6))

    draw_curves(ax, records, zoom_Ls)

    ax.set_xlim(Klo, Khi)
    all_U = []
    for L in zoom_Ls:
        if L not in records: continue
        for K, u, _ in records[L]:
            if Klo <= K <= Khi:
                all_U.append(u)
    if all_U:
        ax.set_ylim(max(0.0, min(all_U) - 0.05), min(1.0, max(all_U) + 0.1))

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    ax.set_xlabel(r"$K = K_1 = K_2 = K_3 = K_t$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title("Binder cumulant crossing — equilateral 3D ($N_t = 4L$, zoom)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def plot_crossing_summary(records, out_path):
    """
    Simple Binder crossing plot — data, error bars and ribbons only.
    """
    sizes = sorted(records.keys())
    fig, ax = plt.subplots(figsize=(9, 6))

    draw_curves(ax, records, sizes)

    ax.set_xlabel(r"$K = K_1 = K_2 = K_3 = K_t$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title(
        r"Binder cumulant crossing — equilateral 3D Ising ($N_t = 4L$, $T_x=T_y=0$)",
        fontsize=11)
    ax.set_xlim(0.1598, 0.1622)
    ax.set_ylim(-0.02, 1.05)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def plot_vs_inv_L(records, out_path):
    """
    U4 vs 1/L, one line per K value.
    Lines for K < Kc fan down toward 0; lines for K > Kc fan up toward 1.
    At Kc all lines meet — the crossing is visible as the bundle of lines
    converging around a common 1/L value.
    """
    # Invert the records dict: {K: [(L, U4, U4_err), ...]}
    by_K = defaultdict(list)
    for L, pts in records.items():
        for K, u4, u4e in pts:
            by_K[K].append((L, u4, u4e))
    for K in by_K:
        by_K[K].sort()   # sort by L

    K_vals = sorted(by_K.keys())
    n_K    = len(K_vals)

    # One distinct color per K using a cycled colormap
    cmap   = plt.cm.turbo
    colors = [cmap(i / (n_K - 1)) for i in range(n_K)]

    fig, ax = plt.subplots(figsize=(11, 6))

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

    # x-axis ticks at the actual 1/L values; label both 1/L and L
    Ls      = sorted(records.keys())
    x_ticks = [1.0 / L for L in Ls]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"1/{L}\n(L={L})" for L in Ls], fontsize=9)

    ax.legend(fontsize=7.5, loc="upper left",
              ncol=2, framealpha=0.9, title=r"$K$", title_fontsize=8)

    ax.set_xlabel(r"$1/L$", fontsize=12)
    ax.set_ylabel(r"Binder cumulant $U_4$", fontsize=12)
    ax.set_title(
        r"$U_4$ vs $1/L$ — equilateral 3D Ising ($N_t=4L$, $T_x=T_y=0$)",
        fontsize=11)
    ax.set_xlim(left=0.0, right=max(x_ticks) * 1.15)
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


# ── entry point ──────────────────────────────────────────────────────────────

def main():
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    DEFAULT_DIR = os.path.join(
        SCRIPT_DIR,
        "out_binder_refine_160_162_nk15_ntraj10000_Nt4L",
    )

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DIR,
                    help="directory containing summary.tsv")
    args = ap.parse_args()

    data_dir = args.data_dir
    tsv_path = os.path.join(data_dir, "summary.tsv")

    if not os.path.exists(tsv_path):
        print(f"ERROR: {tsv_path} not found.", file=sys.stderr)
        sys.exit(1)

    records = read_summary(tsv_path)
    total   = sum(len(v) for v in records.values())
    print(f"Loaded {total} unique (L, K) points from {tsv_path}")
    for L in sorted(records):
        print(f"  L={L}: {len(records[L])} points, "
              f"K=[{records[L][0][0]:.5f}, {records[L][-1][0]:.5f}]")

    out_all    = os.path.join(data_dir, "binder_crossing_all.png")
    out_zoom   = os.path.join(data_dir, "binder_crossing_zoom.png")
    out_summ   = os.path.join(data_dir, "binder_crossing_summary.png")
    out_invL   = os.path.join(data_dir, "binder_vs_inv_L.png")

    plot_all(records, out_all)
    plot_zoom(records, out_zoom)
    plot_crossing_summary(records, out_summ)
    plot_vs_inv_L(records, out_invL)
    print("Done.")


if __name__ == "__main__":
    main()
