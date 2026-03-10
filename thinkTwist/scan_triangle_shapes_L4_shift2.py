#!/usr/bin/env python3
"""Brute-force scan of triangle shapes for L=4, shift=2.

Goal: find local triangle shape that makes twisted-BC two-point functions look
closest to square-target two-point functions.

Model per direction d in {x, y, y-x}:
  C_d(r) = exp(- d_eff_d(r) / xi )
  d_eff_d(r) = folded_distance(r, p_d) * ell_d

Twisted periods use BC:
  (i+L, j) ~ (i, j)
  (i, j+L) ~ (i+shift, j)

Square target uses period L in all directions and local lengths:
  ell_x=1, ell_y=1, ell_y-x=sqrt(2)
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def canonicalize(x: int, y: int, L: int, shift: int) -> tuple[int, int]:
    y0 = y % L
    b = (y - y0) // L
    x0 = (x - b * shift) % L
    return x0, y0


def period(step: tuple[int, int], L: int, shift: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, L, shift)
    for n in range(1, L * L + 1):
        x, y = canonicalize(x + sx, y + sy, L, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError("No period found")


def folded_distance(r: np.ndarray, p: int) -> np.ndarray:
    rem = np.mod(r, p)
    return np.minimum(rem, p - rem)


def correlator(r: np.ndarray, p: int, ell: float, xi: float) -> np.ndarray:
    return np.exp(-(folded_distance(r, p) * ell) / xi)


def mae(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.abs(a - b)))


def shape_lengths(lx: float, ly: float, theta_deg: float) -> tuple[float, float, float]:
    theta = math.radians(theta_deg)
    lr2 = lx * lx + ly * ly - 2.0 * lx * ly * math.cos(theta)
    if lr2 <= 0:
        return lx, ly, float("nan")
    return lx, ly, math.sqrt(lr2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--shift", type=int, default=2)
    parser.add_argument("--xi", type=float, default=2.0)
    parser.add_argument("--rmax", type=int, default=32)
    parser.add_argument("--lx", type=float, default=1.0)
    parser.add_argument("--ly-min", type=float, default=0.2)
    parser.add_argument("--ly-max", type=float, default=2.0)
    parser.add_argument("--ly-step", type=float, default=0.01)
    parser.add_argument("--th-min", type=float, default=20.0)
    parser.add_argument("--th-max", type=float, default=160.0)
    parser.add_argument("--th-step", type=float, default=1.0)
    parser.add_argument("--out-dir", default="thinkTwist/L4_shift2_shape_scan")
    args = parser.parse_args()

    L = args.L
    shift = args.shift
    xi = args.xi

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    p_x = period((1, 0), L, shift)
    p_y = period((0, 1), L, shift)
    p_d = period((-1, 1), L, shift)

    r = np.arange(0, args.rmax + 1)

    c_ref_x = correlator(r, L, 1.0, xi)
    c_ref_y = correlator(r, L, 1.0, xi)
    c_ref_d = correlator(r, L, math.sqrt(2.0), xi)

    rows: list[dict[str, float]] = []
    best = None

    ly_vals = np.arange(args.ly_min, args.ly_max + 0.5 * args.ly_step, args.ly_step)
    th_vals = np.arange(args.th_min, args.th_max + 0.5 * args.th_step, args.th_step)

    for ly in ly_vals:
        for th in th_vals:
            lx, ly_use, lr = shape_lengths(args.lx, float(ly), float(th))
            if not np.isfinite(lr):
                continue

            c_x = correlator(r, p_x, lx, xi)
            c_y = correlator(r, p_y, ly_use, xi)
            c_d = correlator(r, p_d, lr, xi)

            err_x = mae(c_x, c_ref_x)
            err_y = mae(c_y, c_ref_y)
            err_d = mae(c_d, c_ref_d)
            err_tot = (err_x + err_y + err_d) / 3.0

            row = {
                "lx": lx,
                "ly": ly_use,
                "theta_deg": float(th),
                "lr": lr,
                "err_x": err_x,
                "err_y": err_y,
                "err_y_minus_x": err_d,
                "err_total": err_tot,
            }
            rows.append(row)
            if best is None or err_tot < best["err_total"]:
                best = row

    if best is None:
        raise RuntimeError("No shapes scanned")

    # Save full scan table.
    csv_path = out_dir / "shape_scan.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "lx",
                "ly",
                "theta_deg",
                "lr",
                "err_x",
                "err_y",
                "err_y_minus_x",
                "err_total",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    # Save concise report.
    report = out_dir / "best_shape_report.txt"
    with report.open("w") as f:
        f.write("L=4 shift=2 triangle-shape scan for two-point matching\n")
        f.write("===================================================\n\n")
        f.write(f"L={L}, shift={shift}, xi={xi}, rmax={args.rmax}\n")
        f.write(f"twisted periods: p_x={p_x}, p_y={p_y}, p_y-x={p_d}\n\n")
        f.write("Best shape (minimum mean MAE vs square-target two-point curves):\n")
        f.write(f"  lx={best['lx']:.6f}\n")
        f.write(f"  ly={best['ly']:.6f}\n")
        f.write(f"  theta={best['theta_deg']:.6f} deg\n")
        f.write(f"  lr=|y-x|={best['lr']:.6f}\n")
        f.write(f"  err_x={best['err_x']:.6f}\n")
        f.write(f"  err_y={best['err_y']:.6f}\n")
        f.write(f"  err_y-x={best['err_y_minus_x']:.6f}\n")
        f.write(f"  err_total={best['err_total']:.6f}\n")

    # Plot: target vs raw-equilateral vs best-scanned.
    raw_lx, raw_ly, raw_lr = 1.0, 1.0, math.sqrt(2.0)

    c_raw_x = correlator(r, p_x, raw_lx, xi)
    c_raw_y = correlator(r, p_y, raw_ly, xi)
    c_raw_d = correlator(r, p_d, raw_lr, xi)

    c_best_x = correlator(r, p_x, best["lx"], xi)
    c_best_y = correlator(r, p_y, best["ly"], xi)
    c_best_d = correlator(r, p_d, best["lr"], xi)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    cfg = [
        ("x", c_ref_x, c_raw_x, c_best_x, "#1f77b4", best["err_x"]),
        ("y", c_ref_y, c_raw_y, c_best_y, "#d62728", best["err_y"]),
        ("y-x", c_ref_d, c_raw_d, c_best_d, "#2ca02c", best["err_y_minus_x"]),
    ]

    for ax, (name, cref, craw, cbest, col, eb) in zip(axes, cfg):
        eraw = mae(craw, cref)
        ax.plot(r, cref, "k--", lw=2.0, label="square target")
        ax.plot(r, craw, color="#9e9e9e", lw=2.0, label="twisted raw equilateral")
        ax.plot(r, cbest, color=col, lw=2.3, label="twisted best-scanned shape")
        ax.set_ylabel("C(r)")
        ax.set_ylim(-0.02, 1.05)
        ax.grid(alpha=0.25)
        ax.set_title(f"dir {name}: MAE raw={eraw:.4f}, MAE best={eb:.4f}")
        ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("r")
    fig.suptitle(
        f"L=4, shift=2: brute-force triangle-shape fit for two-point functions\n"
        f"best ly={best['ly']:.3f}, theta={best['theta_deg']:.2f} deg, lr={best['lr']:.3f}",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    fig_path = out_dir / "two_point_best_shape_comparison.png"
    fig.savefig(fig_path, dpi=180)
    fig.savefig(out_dir / "two_point_best_shape_comparison.svg")
    plt.close(fig)

    # Heatmap of objective over (ly, theta).
    ly_unique = sorted({row["ly"] for row in rows})
    th_unique = sorted({row["theta_deg"] for row in rows})
    z = np.full((len(th_unique), len(ly_unique)), np.nan)
    ly_to_i = {v: i for i, v in enumerate(ly_unique)}
    th_to_i = {v: i for i, v in enumerate(th_unique)}
    for row in rows:
        z[th_to_i[row["theta_deg"]], ly_to_i[row["ly"]]] = row["err_total"]

    hfig, hax = plt.subplots(figsize=(9, 5.5))
    im = hax.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[min(ly_unique), max(ly_unique), min(th_unique), max(th_unique)],
        cmap="viridis",
    )
    hax.scatter([best["ly"]], [best["theta_deg"]], c="red", s=40, label="best")
    hax.set_xlabel("ly")
    hax.set_ylabel("theta (deg)")
    hax.set_title("Mean MAE vs square-target two-point curves")
    hax.legend(loc="upper right", fontsize=8)
    cbar = hfig.colorbar(im, ax=hax)
    cbar.set_label("err_total")
    hfig.tight_layout()
    hfig.savefig(out_dir / "shape_scan_heatmap.png", dpi=180)
    hfig.savefig(out_dir / "shape_scan_heatmap.svg")
    plt.close(hfig)

    print(f"Wrote {csv_path}")
    print(f"Wrote {report}")
    print(f"Wrote {fig_path}")
    print(f"Best: ly={best['ly']:.6f}, theta={best['theta_deg']:.6f}, lr={best['lr']:.6f}, err_total={best['err_total']:.6f}")


if __name__ == "__main__":
    main()
