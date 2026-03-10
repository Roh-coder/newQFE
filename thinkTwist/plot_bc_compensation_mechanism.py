#!/usr/bin/env python3
"""Visualize exactly how compensating triangle lengths offset twisted-BC periodicities.

For each case (Nx, Ny, Ns), this script shows:
1) BC-induced periods p_x, p_y, p_y-x.
2) Distortion with raw equilateral local lengths.
3) Compensation by choosing local lengths inversely to periods.
4) Effective vectors after period-weighting, matching square-theory ratios.
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


def compensating_lengths(px: int, py: int, pr: int) -> tuple[float, float, float, float]:
    lx = 1.0
    ly = px / py
    lr = math.sqrt(2.0) * px / pr
    low = abs(lx - ly)
    high = lx + ly
    if not (low <= lr <= high):
        raise ValueError("Infeasible compensating triangle")
    return lx, ly, lr, math.degrees(math.acos((lx * lx + ly * ly - lr * lr) / (2.0 * lx * ly)))


def draw_case(nx: int, ny: int, ns: int, out_dir: Path) -> tuple[bool, str]:
    px, py, pr = periods(nx, ny, ns)

    try:
        lx, ly, lr, theta_deg = compensating_lengths(px, py, pr)
    except ValueError:
        return False, f"({nx},{ny},{ns}) infeasible"

    # Equilateral reference local lengths.
    lx_eq, ly_eq, lr_eq = 1.0, 1.0, math.sqrt(2.0)

    # Normalized full-orbit scales relative to x-orbit scale.
    # target square-theory ratio: 1, 1, sqrt(2)
    sx_eq = 1.0
    sy_eq = (py * ly_eq) / (px * lx_eq)
    sr_eq = (pr * lr_eq) / (px * lx_eq)

    sx_cp = 1.0
    sy_cp = (py * ly) / (px * lx)
    sr_cp = (pr * lr) / (px * lx)

    # Build compensating local vectors.
    ex = np.array([lx, 0.0])
    th = math.radians(theta_deg)
    ey = np.array([ly * math.cos(th), ly * math.sin(th)])

    x_eff = px * ex
    y_eff = py * ey

    # Panel figure.
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    # Panel A: Equilateral lattice + BC labels.
    ax = axes[0, 0]
    ex_lat = np.array([1.0, 0.0])
    ey_lat = np.array([0.5, math.sqrt(3.0) / 2.0])

    O = np.array([0.0, 0.0])
    A = nx * ex_lat
    B = ny * ey_lat
    C = A + B

    ax.plot([O[0], A[0]], [O[1], A[1]], color="#111111", lw=2)
    ax.plot([B[0], C[0]], [B[1], C[1]], color="#111111", lw=2)
    ax.plot([O[0], B[0]], [O[1], B[1]], color="#111111", lw=2)
    ax.plot([A[0], C[0]], [A[1], C[1]], color="#111111", lw=2)

    for i in range(nx + 1):
        for j in range(ny + 1):
            p = i * ex_lat + j * ey_lat
            ax.scatter(p[0], p[1], s=6, color="#222222", zorder=3)
            if i + 1 <= nx:
                q = (i + 1) * ex_lat + j * ey_lat
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#cccccc", lw=0.6)
            if j + 1 <= ny:
                q = i * ex_lat + (j + 1) * ey_lat
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#cccccc", lw=0.6)
            if i - 1 >= 0 and j + 1 <= ny:
                q = (i - 1) * ex_lat + (j + 1) * ey_lat
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#dddddd", lw=0.5)

    # BC arrows.
    m_left = 0.5 * (O + B)
    m_right = 0.5 * (A + C)
    ax.annotate("", xy=m_right, xytext=m_left, arrowprops=dict(arrowstyle="<->", lw=1.6, color="#1f77b4"))
    ax.text(*(0.5 * (m_left + m_right)), "(i+Nx,j)~(i,j)", color="#1f77b4", fontsize=8, ha="center", va="bottom")

    shift = ns % nx
    for i in sorted(set([0, nx // 2, nx - 1])):
        p_top = i * ex_lat + ny * ey_lat
        p_bot = ((i + shift) % nx) * ex_lat
        ax.annotate("", xy=p_bot, xytext=p_top, arrowprops=dict(arrowstyle="->", lw=1.3, color="#d62728"))

    ax.text(B[0] + 0.5 * A[0], B[1] + 0.3, f"(i,Ny)~(i+Ns,0), Ns={shift}", color="#d62728", fontsize=8, ha="center")
    ax.text(O[0], O[1] - 0.45, f"p_x={px}, p_y={py}, p_y-x={pr}", fontsize=9, ha="left")

    ax.set_title("A) Equilateral Lattice With Twisted BC")
    ax.set_aspect("equal")
    ax.axis("off")

    # Panel B: Raw equilateral orbit-scale distortion.
    ax = axes[0, 1]
    labels = ["x", "y", "y-x"]
    vals_eq = [sx_eq, sy_eq, sr_eq]
    target = [1.0, 1.0, math.sqrt(2.0)]
    xloc = np.arange(3)
    w = 0.36
    ax.bar(xloc - w / 2, vals_eq, width=w, color="#9e9e9e", label="raw equilateral")
    ax.bar(xloc + w / 2, target, width=w, color="#4caf50", alpha=0.7, label="square target")
    ax.set_xticks(xloc, labels)
    ax.set_ylabel("normalized full-orbit scale")
    ax.set_title("B) BC Distortion Before Compensation")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    # Panel C: Compensating side-length choice.
    ax = axes[1, 0]
    ax.bar([0, 1, 2], [lx_eq, ly_eq, lr_eq], width=0.35, color="#bbbbbb", label="equilateral local lengths")
    ax.bar([0.4, 1.4, 2.4], [lx, ly, lr], width=0.35, color="#1976d2", label="compensating local lengths")
    ax.set_xticks([0.2, 1.2, 2.2], ["|x|", "|y|", "|y-x|"])
    ax.set_title("C) Length Adjustment That Compensates BC")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    txt = (
        f"|x|=1\n"
        f"|y|=p_x/p_y={ly:.4f}\n"
        f"|y-x|=sqrt(2)*p_x/p_y-x={lr:.4f}\n"
        f"theta={theta_deg:.2f} deg"
    )
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="#f5f5f5", ec="#cccccc"))

    # Panel D: Post-compensation effective vectors.
    ax = axes[1, 1]
    ax.plot([0, x_eff[0]], [0, x_eff[1]], color="#1f77b4", lw=3, label="p_x x")
    ax.plot([0, y_eff[0]], [0, y_eff[1]], color="#d62728", lw=3, label="p_y y")
    ax.plot([x_eff[0], y_eff[0]], [x_eff[1], y_eff[1]], color="#2ca02c", lw=3, label="p_y y - p_x x")
    ax.scatter([0, x_eff[0], y_eff[0]], [0, x_eff[1], y_eff[1]], color="black", s=20)
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_title("D) After Compensation: Square-Theory Ratios")
    ax.legend(fontsize=8, loc="upper right")

    ang = math.degrees(math.acos(max(-1.0, min(1.0, np.dot(x_eff, y_eff) / (np.linalg.norm(x_eff) * np.linalg.norm(y_eff))))))
    ratio = np.linalg.norm(y_eff - x_eff) / np.linalg.norm(x_eff)
    txt2 = (
        f"|p_x x| : |p_y y| = 1 : {np.linalg.norm(y_eff)/np.linalg.norm(x_eff):.6f}\n"
        f"angle = {ang:.6f} deg\n"
        f"|p_y y - p_x x|/|p_x x| = {ratio:.6f}\n"
        f"target = 1, 90 deg, sqrt(2)"
    )
    ax.text(0.03, 0.97, txt2, transform=ax.transAxes, ha="left", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="#f5f5f5", ec="#cccccc"))

    fig.suptitle(
        f"How Triangle Compensation Offsets BC Distortion: Nx={nx}, Ny={ny}, Ns={ns}",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = f"bc_compensation_Nx{nx}_Ny{ny}_Ns{ns}"
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    return True, png.name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        default="4x4x2,6x5x4,7x7x7",
        help="Comma-separated Nx x Ny x Ns cases",
    )
    parser.add_argument(
        "--out-dir",
        default="thinkTwist/bc_compensation_diagrams",
        help="Output directory",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = parse_cases(args.cases)
    made: list[str] = []
    skipped: list[str] = []

    for nx, ny, ns in cases:
        ok, msg = draw_case(nx, ny, ns, out_dir)
        if ok:
            made.append(msg)
        else:
            skipped.append(msg)

    readme = out_dir / "README.txt"
    lines = [
        "BC Compensation Mechanism Diagrams",
        "==================================",
        "",
        f"Requested cases: {args.cases}",
        f"Generated: {len(made)}",
    ]
    for m in made:
        lines.append(m)
    lines.append("")
    lines.append(f"Skipped (infeasible): {len(skipped)}")
    for s in skipped:
        lines.append(s)

    readme.write_text("\n".join(lines) + "\n")

    print(f"Wrote {readme}")
    for m in made:
        print(m)
    for s in skipped:
        print(f"SKIP {s}")


if __name__ == "__main__":
    main()
