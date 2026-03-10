#!/usr/bin/env python3
"""Plot feasible twisted-lattice cases: BC schematic + compensating triangle.

For each (Nx, Ny, Ns) case, this script:
1) Computes periods p_x, p_y, p_y_minus_x.
2) Computes square-target compensating side lengths (scale s=1):
     l_x = 1
     l_y = p_x / p_y
     l_r = sqrt(2) * p_x / p_y_minus_x
3) Keeps only geometrically feasible cases.
4) Plots:
   - Left: finite lattice patch in deformed basis with BC annotations.
   - Right: the compensating triangle (x, y, y-x edges).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def p_x(nx: int, ny: int, ns: int) -> int:
    _ = ny, ns
    return nx


def p_y(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ns)


def p_diag(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ny + ns)


def parse_case(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace(" ", "").split("x")
    if len(parts) != 3:
        raise ValueError(f"Case '{text}' must be Nx x Ny x Ns")
    nx, ny, ns = (int(parts[0]), int(parts[1]), int(parts[2]))
    return nx, ny, ns


def parse_cases(text: str) -> list[tuple[int, int, int]]:
    out: list[tuple[int, int, int]] = []
    for chunk in text.split(","):
        c = chunk.strip()
        if c:
            out.append(parse_case(c))
    return out


def case_geometry(nx: int, ny: int, ns: int) -> dict[str, float]:
    px = p_x(nx, ny, ns)
    py = p_y(nx, ny, ns)
    pr = p_diag(nx, ny, ns)

    lx = 1.0
    ly = px / py
    lr = math.sqrt(2.0) * px / pr

    lower = abs(lx - ly)
    upper = lx + ly
    feasible = lower <= lr <= upper

    cos_theta = (lx * lx + ly * ly - lr * lr) / (2.0 * lx * ly)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.degrees(math.acos(cos_theta))

    return {
        "px": float(px),
        "py": float(py),
        "pr": float(pr),
        "lx": lx,
        "ly": ly,
        "lr": lr,
        "lower": lower,
        "upper": upper,
        "feasible": 1.0 if feasible else 0.0,
        "theta": theta,
    }


def draw_case(nx: int, ny: int, ns: int, out_dir: Path) -> tuple[bool, Path | None, str]:
    g = case_geometry(nx, ny, ns)
    feasible = bool(g["feasible"])

    if not feasible:
        return False, None, (
            f"({nx},{ny},{ns}) infeasible: l_r={g['lr']:.4f} not in "
            f"[{g['lower']:.4f}, {g['upper']:.4f}]"
        )

    # Left panel uses a fixed equilateral triangular lattice for clear BC display.
    ex_lat = np.array([1.0, 0.0])
    ey_lat = np.array([0.5, math.sqrt(3.0) / 2.0])

    # Right panel uses the compensating local triangle from the solved side lengths.
    lx = g["lx"]
    ly = g["ly"]
    theta = math.radians(g["theta"])
    ex_def = np.array([lx, 0.0])
    ey_def = np.array([ly * math.cos(theta), ly * math.sin(theta)])

    O = np.array([0.0, 0.0])
    A = nx * ex_lat
    B = ny * ey_lat
    C = A + B

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    ax = axes[0]

    # Draw patch boundary.
    ax.plot([O[0], A[0]], [O[1], A[1]], color="#111111", lw=2.0)  # bottom
    ax.plot([B[0], C[0]], [B[1], C[1]], color="#111111", lw=2.0)  # top
    ax.plot([O[0], B[0]], [O[1], B[1]], color="#111111", lw=2.0)  # left
    ax.plot([A[0], C[0]], [A[1], C[1]], color="#111111", lw=2.0)  # right

    # Draw lattice nodes and local bonds.
    pts = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            p = i * ex_lat + j * ey_lat
            pts[(i, j)] = p
            ax.scatter(p[0], p[1], s=10, color="#222222", zorder=3)

    for i in range(nx + 1):
        for j in range(ny + 1):
            p = pts[(i, j)]
            if i + 1 <= nx:
                q = pts[(i + 1, j)]
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#bbbbbb", lw=0.8, zorder=1)
            if j + 1 <= ny:
                q = pts[(i, j + 1)]
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#bbbbbb", lw=0.8, zorder=1)
            if i - 1 >= 0 and j + 1 <= ny:
                q = pts[(i - 1, j + 1)]
                ax.plot([p[0], q[0]], [p[1], q[1]], color="#d9d9d9", lw=0.7, zorder=1)

    # Boundary-condition annotations.
    mid_left = 0.5 * (O + B)
    mid_right = 0.5 * (A + C)
    ax.annotate(
        "",
        xy=(mid_right[0], mid_right[1]),
        xytext=(mid_left[0], mid_left[1]),
        arrowprops=dict(arrowstyle="<->", color="#1f77b4", lw=1.8),
    )
    ax.text(
        *(0.5 * (mid_left + mid_right)),
        "(i+Nx,j)~(i,j)",
        color="#1f77b4",
        fontsize=9,
        ha="center",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8),
    )

    shift_mod = ns % nx
    # Draw a few explicit top->bottom twisted identifications: (i,Ny) ~ (i+Ns,0).
    sample_i = [0, nx // 2, nx - 1] if nx > 1 else [0]
    sample_i = sorted(set(sample_i))
    for i in sample_i:
        p_top = i * ex_lat + ny * ey_lat
        p_bot = ((i + shift_mod) % nx) * ex_lat
        ax.annotate(
            "",
            xy=(p_bot[0], p_bot[1]),
            xytext=(p_top[0], p_top[1]),
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.5, alpha=0.9),
        )

    t_top = B + 0.5 * A
    ax.text(
        t_top[0],
        t_top[1] + 0.08 * (np.linalg.norm(B - O) + 1e-6),
        f"(i,Ny)~(i+Ns,0), Ns={shift_mod}",
        color="#d62728",
        fontsize=9,
        ha="center",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8),
    )

    ax.set_title(f"Equilateral Triangular Lattice + BC: Nx={nx}, Ny={ny}, Ns={ns}")
    ax.set_aspect("equal")
    ax.axis("off")

    # Triangle panel.
    ax2 = axes[1]
    X = ex_def
    Y = ey_def
    ax2.plot([0, X[0]], [0, X[1]], color="#1f77b4", lw=3, label="x edge")
    ax2.plot([0, Y[0]], [0, Y[1]], color="#d62728", lw=3, label="y edge")
    ax2.plot([X[0], Y[0]], [X[1], Y[1]], color="#2ca02c", lw=3, label="y-x edge")
    ax2.scatter([0, X[0], Y[0]], [0, X[1], Y[1]], color="black", s=30)

    ax2.text(0.5 * X[0], 0.5 * X[1] - 0.02, f"|x|={g['lx']:.3f}", color="#1f77b4")
    ax2.text(0.5 * Y[0], 0.5 * Y[1] + 0.02, f"|y|={g['ly']:.3f}", color="#d62728")
    ax2.text(
        0.5 * (X[0] + Y[0]) + 0.02,
        0.5 * (X[1] + Y[1]),
        f"|y-x|={g['lr']:.3f}",
        color="#2ca02c",
    )

    txt = (
        f"p_x={int(g['px'])}, p_y={int(g['py'])}, p_y-x={int(g['pr'])}\n"
        f"theta_xy={g['theta']:.2f} deg"
    )
    ax2.text(
        0.03,
        0.97,
        txt,
        transform=ax2.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f5f5f5", ec="#cccccc"),
    )

    ax2.set_title("Compensating Triangle (Square-Target)")
    ax2.set_aspect("equal")
    ax2.grid(alpha=0.2)
    ax2.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "Feasible twisted case: lattice BC and local triangle deformation",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    stem = f"case_Nx{nx}_Ny{ny}_Ns{ns}"
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    return True, png, f"({nx},{ny},{ns}) theta={g['theta']:.2f} deg"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        default="3x4x2,4x4x2,6x5x4,6x6x3,7x7x7,8x8x4",
        help="Comma-separated case list as Nx x Ny x Ns",
    )
    parser.add_argument(
        "--out-dir",
        default="thinkTwist/feasible_case_plots",
        help="Output directory for plots",
    )
    args = parser.parse_args()

    cases = parse_cases(args.cases)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    made: list[str] = []
    skipped: list[str] = []

    for nx, ny, ns in cases:
        ok, out_path, msg = draw_case(nx, ny, ns, out_dir)
        if ok:
            assert out_path is not None
            made.append(f"{msg} -> {out_path}")
        else:
            skipped.append(msg)

    summary = out_dir / "README_cases.txt"
    lines = [
        "Feasible Case Plots",
        "===================",
        "",
        f"Requested cases: {args.cases}",
        "",
        f"Generated: {len(made)}",
    ]
    lines.extend(made)
    lines.append("")
    lines.append(f"Skipped (infeasible): {len(skipped)}")
    lines.extend(skipped)
    summary.write_text("\n".join(lines) + "\n")

    print(f"Wrote {summary}")
    for line in made:
        print(line)
    for line in skipped:
        print(f"SKIP {line}")


if __name__ == "__main__":
    main()
