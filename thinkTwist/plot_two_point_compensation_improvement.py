#!/usr/bin/env python3
"""Show two-point improvement after compensation.

For each feasible case (Nx, Ny, Ns), compare per-direction correlators:
- Square-theory target
- Twisted BC with raw equilateral local lengths
- Twisted BC with compensating local lengths

Model: C(r) = exp(-d_eff(r)/xi), with folded periodic distance.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_case(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace(" ", "").split("x")
    if len(parts) != 3:
        raise ValueError(f"Invalid case '{text}', expected Nx x Ny x Ns")
    return int(parts[0]), int(parts[1]), int(parts[2])


def parse_cases(text: str) -> list[tuple[int, int, int]]:
    out: list[tuple[int, int, int]] = []
    for chunk in text.split(","):
        c = chunk.strip()
        if c:
            out.append(parse_case(c))
    return out


def periods(nx: int, ny: int, ns: int) -> tuple[int, int, int]:
    px = nx
    py = nx * ny // math.gcd(nx, ns)
    pr = nx * ny // math.gcd(nx, ny + ns)
    return px, py, pr


def folded_distance(r: np.ndarray, p: int) -> np.ndarray:
    rem = np.mod(r, p)
    return np.minimum(rem, p - rem)


def compensating_lengths(px: int, py: int, pr: int) -> tuple[bool, tuple[float, float, float]]:
    lx = 1.0
    ly = px / py
    lr = math.sqrt(2.0) * px / pr
    feasible = abs(lx - ly) <= lr <= (lx + ly)
    return feasible, (lx, ly, lr)


def correlator(r: np.ndarray, p: int, ell: float, xi: float) -> np.ndarray:
    d = folded_distance(r, p) * ell
    return np.exp(-d / xi)


def draw_case(nx: int, ny: int, ns: int, xi: float, out_dir: Path) -> tuple[bool, str]:
    px, py, pr = periods(nx, ny, ns)
    feasible, (lx, ly, lr) = compensating_lengths(px, py, pr)
    if not feasible:
        return False, f"({nx},{ny},{ns}) infeasible"

    # Periods for directions under twisted BC.
    p_tw = {
        "x": px,
        "y": py,
        "y-x": pr,
    }

    # Square-theory target periods and local lengths.
    # Use p_ref=px for all directions so x,y,y-x are isotropically periodic in the target.
    p_ref = px
    ell_ref = {
        "x": 1.0,
        "y": 1.0,
        "y-x": math.sqrt(2.0),
    }

    # Raw equilateral local lengths on triangular lattice.
    ell_raw = {
        "x": 1.0,
        "y": 1.0,
        "y-x": math.sqrt(2.0),
    }

    # Compensating local lengths from periodicities.
    ell_cmp = {
        "x": lx,
        "y": ly,
        "y-x": lr,
    }

    pmax = max(px, py, pr, p_ref)
    r = np.arange(0, 4 * pmax + 1)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    directions = [
        ("x", "#1f77b4"),
        ("y", "#d62728"),
        ("y-x", "#2ca02c"),
    ]

    improvements = []

    for ax, (name, color) in zip(axes, directions):
        c_ref = correlator(r, p_ref, ell_ref[name], xi)
        c_raw = correlator(r, p_tw[name], ell_raw[name], xi)
        c_cmp = correlator(r, p_tw[name], ell_cmp[name], xi)

        mae_raw = float(np.mean(np.abs(c_raw - c_ref)))
        mae_cmp = float(np.mean(np.abs(c_cmp - c_ref)))
        improvements.append((name, mae_raw, mae_cmp))

        ax.plot(r, c_ref, "k--", lw=2.0, label="square-theory target")
        ax.plot(r, c_raw, color="#9e9e9e", lw=2.0, label="twisted + raw equilateral")
        ax.plot(r, c_cmp, color=color, lw=2.3, label="twisted + compensated")

        ax.set_ylabel("C(r)")
        ax.set_ylim(-0.02, 1.05)
        ax.grid(alpha=0.25)
        ax.set_title(
            f"dir {name}: p_tw={p_tw[name]}, p_ref={p_ref}, "
            f"MAE raw={mae_raw:.4f}, MAE comp={mae_cmp:.4f}"
        )
        ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("trajectory separation r")

    n_better = sum(1 for _, a, b in improvements if b < a)
    fig.suptitle(
        f"Two-point compensation test: Nx={nx}, Ny={ny}, Ns={ns}, xi={xi:.3g}"
        f"\ncompensation better in {n_better}/3 directions (lower MAE to square target)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    stem = f"two_point_compensation_Nx{nx}_Ny{ny}_Ns{ns}"
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    return True, (
        f"({nx},{ny},{ns}) better_dirs={n_better}/3 "
        f"raw_mae={[round(x[1],4) for x in improvements]} "
        f"comp_mae={[round(x[2],4) for x in improvements]}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        default="4x4x2,6x5x4,7x7x7,8x8x4",
        help="Comma-separated Nx x Ny x Ns cases",
    )
    parser.add_argument("--xi", type=float, default=2.0, help="Decay length")
    parser.add_argument(
        "--out-dir",
        default="thinkTwist/two_point_compensation",
        help="Output directory",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    done: list[str] = []
    skipped: list[str] = []

    for nx, ny, ns in parse_cases(args.cases):
        ok, msg = draw_case(nx, ny, ns, args.xi, out_dir)
        if ok:
            done.append(msg)
        else:
            skipped.append(msg)

    readme = out_dir / "README.txt"
    lines = [
        "Two-Point Compensation Improvement",
        "===============================",
        "",
        f"cases = {args.cases}",
        f"xi = {args.xi}",
        f"generated = {len(done)}",
    ]
    lines.extend(done)
    lines.append("")
    lines.append(f"skipped = {len(skipped)}")
    lines.extend(skipped)
    readme.write_text("\n".join(lines) + "\n")

    print(f"Wrote {readme}")
    for d in done:
        print(d)
    for s in skipped:
        print(f"SKIP {s}")


if __name__ == "__main__":
    main()
