#!/usr/bin/env python3
"""Solve compensating triangle shape for twisted BC over multiple simple cases.

For BC
  (i + L, j) ~ (i, j)
  (i, j + L) ~ (i + shift, j)

and directions x=(1,0), y=(0,1), r=(y-x)=(-1,1), we compute periods
p_x, p_y, p_r under the quotient and choose local lengths (scale |x|=1):

  |y|   = p_x / p_y
  |r|   = sqrt(2) * p_x / p_r

so that half-orbit physical scales match square-continuum targets.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class CaseResult:
    L: int
    shift: int
    px: int
    py: int
    pr: int
    lx: float
    ly: float
    lr: float
    feasible: bool
    lr_min: float
    lr_max: float
    cos_theta: float
    theta_deg: float
    theta_projected_deg: float


def is_equivalent_to_zero(x: int, y: int, L: int, shift: int) -> bool:
    if y % L != 0:
        return False
    b = y // L
    rem_x = x - b * shift
    return rem_x % L == 0


def canonicalize(x: int, y: int, L: int, shift: int) -> tuple[int, int]:
    """Map integer lattice coords to a canonical representative in the quotient.

    Quotient generators are g1=(L,0), g2=(shift,L). We reduce y modulo L and
    compensate x by the removed g2 component, then reduce x modulo L.
    """

    y0 = y % L
    b = (y - y0) // L
    x0 = (x - b * shift) % L
    return x0, y0


def period(step: tuple[int, int], L: int, shift: int) -> int:
    """Exact period on the finite quotient lattice.

    Since the quotient has L^2 states, a period is guaranteed in <= L^2 steps.
    """

    sx, sy = step
    x, y = 0, 0
    x, y = canonicalize(x, y, L, shift)
    for n in range(1, L * L + 1):
        x, y = canonicalize(x + sx, y + sy, L, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period found for step={step}, L={L}, shift={shift}")


def solve_case(L: int, shift: int) -> CaseResult:
    if L <= 0:
        raise ValueError("L must be positive")
    if shift < 0:
        raise ValueError("shift must be >= 0")

    px = period((1, 0), L, shift)
    py = period((0, 1), L, shift)
    pr = period((-1, 1), L, shift)

    lx = 1.0
    ly = px / py
    lr = math.sqrt(2.0) * px / pr

    lr_min = abs(lx - ly)
    lr_max = lx + ly
    feasible = (lr_min <= lr <= lr_max)

    # From |y-x|^2 = |x|^2 + |y|^2 - 2|x||y|cos(theta)
    denom = 2.0 * lx * ly
    if denom == 0.0:
        raise RuntimeError("Degenerate geometry: zero denominator in angle solve")

    cos_raw = (lx * lx + ly * ly - lr * lr) / denom

    if feasible:
        cos_theta = cos_raw
        theta_deg = math.degrees(math.acos(cos_theta))
    else:
        cos_theta = float("nan")
        theta_deg = float("nan")

    # Projected boundary angle after clipping (what the old code effectively showed).
    cos_proj = max(-1.0, min(1.0, cos_raw))
    theta_projected_deg = math.degrees(math.acos(cos_proj))

    return CaseResult(
        L,
        shift,
        px,
        py,
        pr,
        lx,
        ly,
        lr,
        feasible,
        lr_min,
        lr_max,
        cos_theta,
        theta_deg,
        theta_projected_deg,
    )


def parse_case(s: str) -> tuple[int, int]:
    try:
        l_str, sh_str = s.split(":", 1)
        return int(l_str), int(sh_str)
    except Exception as exc:
        raise argparse.ArgumentTypeError(
            f"Case '{s}' must be formatted as L:shift (e.g., 64:32)"
        ) from exc


def make_sweep_plots(results: list[CaseResult], out_prefix: Path) -> tuple[Path, Path]:
    shifts = np.array([r.shift for r in results], dtype=float)
    theta = np.array([r.theta_deg if r.feasible else np.nan for r in results], dtype=float)
    theta_proj = np.array([r.theta_projected_deg for r in results], dtype=float)
    ly = np.array([r.ly for r in results], dtype=float)
    lr = np.array([r.lr for r in results], dtype=float)
    py = np.array([r.py for r in results], dtype=float)
    pr = np.array([r.pr for r in results], dtype=float)
    L = results[0].L

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(shifts, theta, "o-", color="#1f77b4", label="theta_xy feasible (deg)")
    axes[0].plot(shifts, theta_proj, "--", color="#999999", alpha=0.7, label="theta projected (deg)")
    axes[0].set_ylabel("triangle angle theta_xy (deg)")
    axes[0].set_title(f"Compensating shape vs twist shift (L={L})")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].plot(shifts, ly, "o-", color="#d62728", label="|y| / |x|")
    axes[1].plot(shifts, lr, "o-", color="#2ca02c", label="|y-x| / |x|")
    axes[1].plot(shifts, py / L, "--", color="#d62728", alpha=0.35, label="p_y / L")
    axes[1].plot(shifts, pr / L, "--", color="#2ca02c", alpha=0.35, label="p_(y-x) / L")
    axes[1].set_xlabel("shift")
    axes[1].set_ylabel("relative length")
    axes[1].grid(alpha=0.25)
    axes[1].legend(ncol=2, fontsize=9)

    fig.tight_layout()

    png_path = out_prefix.with_suffix(".png")
    svg_path = out_prefix.with_suffix(".svg")
    fig.savefig(png_path, dpi=180)
    fig.savefig(svg_path)
    plt.close(fig)
    return png_path, svg_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        type=parse_case,
        action="append",
        dest="cases",
        help="Case as L:shift. Repeat this flag for multiple cases.",
    )
    parser.add_argument(
        "--out-prefix",
        default="thinkTwist/twisted_shape_cases",
        help="Prefix for generated report files (without extension)",
    )
    parser.add_argument(
        "--sweep-L",
        type=int,
        default=None,
        help="If set, ignore --case and sweep shift=0..L-1 for this L.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable PNG/SVG plot generation in sweep mode.",
    )
    args = parser.parse_args()

    # A small default suite of simple benchmark cases.
    if args.sweep_L is not None:
        L = args.sweep_L
        if L <= 0:
            raise ValueError("--sweep-L must be positive")
        cases = [(L, shift) for shift in range(L)]
    else:
        cases = args.cases or [
            (8, 0),
            (8, 2),
            (8, 4),
            (16, 0),
            (16, 8),
            (64, 0),
            (64, 32),
        ]

    results = [solve_case(L, shift) for (L, shift) in cases]

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    txt_path = out_prefix.with_suffix(".txt")
    csv_path = out_prefix.with_suffix(".csv")

    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "L",
                "shift",
                "p_x",
                "p_y",
                "p_y_minus_x",
                "|x|",
                "|y|",
                "|y-x|",
                "feasible",
                "|y-x|_min",
                "|y-x|_max",
                "cos(theta_xy)",
                "theta_xy_deg",
                "theta_projected_deg",
            ]
        )
        for r in results:
            w.writerow(
                [
                    r.L,
                    r.shift,
                    r.px,
                    r.py,
                    r.pr,
                    f"{r.lx:.8f}",
                    f"{r.ly:.8f}",
                    f"{r.lr:.8f}",
                    str(r.feasible),
                    f"{r.lr_min:.8f}",
                    f"{r.lr_max:.8f}",
                    "nan" if math.isnan(r.cos_theta) else f"{r.cos_theta:.8f}",
                    "nan" if math.isnan(r.theta_deg) else f"{r.theta_deg:.8f}",
                    f"{r.theta_projected_deg:.8f}",
                ]
            )

    lines = []
    lines.append("Compensating Triangle Cases")
    lines.append("==========================")
    lines.append("")
    lines.append("Columns: L, shift, periods (px, py, p(y-x)), and solved shape with |x|=1")
    lines.append("feasible=False means exact Euclidean triangle realization is impossible for that case")
    lines.append("theta_projected is angle after clipping to boundary (legacy artifact diagnostic)")
    lines.append("")
    lines.append(
        "{:>4} {:>6} {:>4} {:>4} {:>7} {:>8} {:>10} {:>10} {:>9} {:>12} {:>12}".format(
            "L", "shift", "px", "py", "p(y-x)", "|x|", "|y|", "|y-x|", "feasible", "theta", "theta_proj"
        )
    )
    lines.append("-" * 108)

    for r in results:
        theta_text = "nan" if math.isnan(r.theta_deg) else f"{r.theta_deg:.6f}"
        lines.append(
            "{:>4d} {:>6d} {:>4d} {:>4d} {:>7d} {:>8.3f} {:>10.6f} {:>10.6f} {:>9} {:>12} {:>12.6f}".format(
                r.L,
                r.shift,
                r.px,
                r.py,
                r.pr,
                r.lx,
                r.ly,
                r.lr,
                str(r.feasible),
                theta_text,
                r.theta_projected_deg,
            )
        )

    feasible_count = sum(1 for r in results if r.feasible)
    lines.append("")
    lines.append(f"Feasible cases: {feasible_count} / {len(results)}")

    txt_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {txt_path}")
    print(f"Wrote {csv_path}")

    if args.sweep_L is not None and not args.no_plots:
        png_path, svg_path = make_sweep_plots(results, out_prefix)
        print(f"Wrote {png_path}")
        print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
