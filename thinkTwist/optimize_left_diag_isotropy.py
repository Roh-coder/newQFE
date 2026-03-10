#!/usr/bin/env python3
"""Optimize left/diag lengths for directional isotropy of two-point functions.

Model:
    C_d(n) = exp(-(n_fold_d(n)/p_d) * ell_d_over_a)

with n_fold_d(n) = min(n mod p_d, p_d-(n mod p_d)).

Optimization target:
  x length fixed to 1.0.
  Choose left and diag lengths to make C_left(n), C_diag(n) as close as possible
  to C_x(n) over one x-period n=0..p_x-1.
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


def n_fold(n: np.ndarray, p: int) -> np.ndarray:
    rem = np.mod(n, p)
    return np.minimum(rem, p - rem)


def corr(n: np.ndarray, p: int, ell_over_a: float) -> np.ndarray:
    frac_sep = n_fold(n, p) / float(p)
    return np.exp(-(frac_sep * ell_over_a))


def mse(a: np.ndarray, b: np.ndarray) -> float:
    d = a - b
    return float(np.mean(d * d))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=16)
    parser.add_argument("--shift", type=int, default=8)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument("--left-min", type=float, default=0.05)
    parser.add_argument("--left-max", type=float, default=4.0)
    parser.add_argument("--left-step", type=float, default=0.002)
    parser.add_argument("--diag-min", type=float, default=0.05)
    parser.add_argument("--diag-max", type=float, default=4.0)
    parser.add_argument("--diag-step", type=float, default=0.002)
    parser.add_argument("--out-dir", default="thinkTwist/isotropy_optimization")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    px = period((1, 0), args.L, args.shift)
    pl = period((0, 1), args.L, args.shift)
    pd = period((-1, 1), args.L, args.shift)

    # One x-period, as requested.
    n = np.arange(0, px)
    cx = corr(n, px, args.x_length)

    left_vals = np.arange(args.left_min, args.left_max + 0.5 * args.left_step, args.left_step)
    diag_vals = np.arange(args.diag_min, args.diag_max + 0.5 * args.diag_step, args.diag_step)

    best = None

    # Cache per-direction correlator arrays for speed.
    left_cache = {float(v): corr(n, pl, float(v)) for v in left_vals}
    diag_cache = {float(v): corr(n, pd, float(v)) for v in diag_vals}

    for lv in left_vals:
        cl = left_cache[float(lv)]
        err_l = mse(cl, cx)
        for dv in diag_vals:
            cd = diag_cache[float(dv)]
            err_d = mse(cd, cx)
            obj = 0.5 * (err_l + err_d)
            if best is None or obj < best["obj"]:
                best = {
                    "left": float(lv),
                    "diag": float(dv),
                    "obj": float(obj),
                    "err_left": float(err_l),
                    "err_diag": float(err_d),
                }

    assert best is not None

    # Save report.
    report = out_dir / f"best_lengths_L{args.L}_shift{args.shift}.txt"
    with report.open("w") as f:
        f.write("Directional isotropy optimization (one x-period)\n")
        f.write("=============================================\n\n")
        f.write(f"L={args.L}, shift={args.shift}\n")
        f.write(f"periods: p_x={px}, p_left={pl}, p_diag={pd}\n")
        f.write(f"x_length/a fixed = {args.x_length:.6f}\n")
        f.write(f"n-range = 0..{px-1} (one x-period)\n\n")
        f.write("Objective: 0.5*(MSE(C_left,C_x)+MSE(C_diag,C_x))\n")
        f.write("with C_d(n)=exp(-(n_fold_d(n)/p_d)*ell_d/a).\n\n")
        f.write("Best lengths:\n")
        f.write(f"  left_length/a = {best['left']:.6f}\n")
        f.write(f"  diag_length/a = {best['diag']:.6f}\n")
        f.write(f"  err_left (MSE vs x) = {best['err_left']:.8f}\n")
        f.write(f"  err_diag (MSE vs x) = {best['err_diag']:.8f}\n")
        f.write(f"  objective = {best['obj']:.8f}\n")

    # Save a local neighborhood table around optimum for transparency.
    csv_path = out_dir / f"best_neighborhood_L{args.L}_shift{args.shift}.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["left_length_over_a", "diag_length_over_a", "err_left", "err_diag", "objective"])
        for lv in np.arange(max(args.left_min, best["left"] - 0.05), min(args.left_max, best["left"] + 0.05) + 1e-12, 0.01):
            cl = corr(n, pl, float(lv))
            err_l = mse(cl, cx)
            for dv in np.arange(max(args.diag_min, best["diag"] - 0.05), min(args.diag_max, best["diag"] + 0.05) + 1e-12, 0.01):
                cd = corr(n, pd, float(dv))
                err_d = mse(cd, cx)
                obj = 0.5 * (err_l + err_d)
                w.writerow([f"{lv:.6f}", f"{dv:.6f}", f"{err_l:.10f}", f"{err_d:.10f}", f"{obj:.10f}"])

    # Plot curves with optimum.
    cl_best = corr(n, pl, best["left"])
    cd_best = corr(n, pd, best["diag"])

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.plot(n, cx, "-o", color="#1f77b4", label=f"x (p={px}, ell/a={args.x_length:.3f})")
    ax.plot(n, cl_best, "-o", color="#d62728", label=f"left (p={pl}, ell/a={best['left']:.3f})")
    ax.plot(n, cd_best, "-o", color="#2ca02c", label=f"diag (p={pd}, ell/a={best['diag']:.3f})")
    ax.set_xlabel("n (within one x-period)")
    ax.set_ylabel("C_d(n)")
    ax.set_ylim(-0.02, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(
        f"Optimized directional isotropy: L={args.L}, shift={args.shift}\n"
        f"objective={best['obj']:.7f}, err_left={best['err_left']:.7f}, err_diag={best['err_diag']:.7f}"
    )
    fig.tight_layout()

    png = out_dir / f"optimized_isotropy_curves_L{args.L}_shift{args.shift}.png"
    svg = out_dir / f"optimized_isotropy_curves_L{args.L}_shift{args.shift}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    print(f"Wrote {report}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {png}")
    print(
        "Best "
        f"left={best['left']:.6f} diag={best['diag']:.6f} "
        f"obj={best['obj']:.8f}"
    )


if __name__ == "__main__":
    main()
