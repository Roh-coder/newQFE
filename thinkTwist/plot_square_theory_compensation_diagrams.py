#!/usr/bin/env python3
"""Create diagrams showing compensating triangles induce square-theory effective geometry.

For each feasible (Nx, Ny, Ns):
1) Build compensating local triangle (x, y, y-x) using periodicity rules.
2) Form period-weighted effective vectors:
     X_eff = (p_x/2) * x
     Y_eff = (p_y/2) * y
     D_eff = Y_eff - X_eff = (p_y_minus_x/2) * (y-x)
3) Plot and report that:
     |X_eff| = |Y_eff|,
     angle(X_eff, Y_eff) = 90 deg,
     |D_eff| = sqrt(2) * |X_eff|.
These are the square-theory directional relations.
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
        raise ValueError(f"Case '{text}' must be Nx x Ny x Ns")
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


def compensating_triangle(nx: int, ny: int, ns: int) -> tuple[bool, dict[str, float]]:
    px, py, pr = periods(nx, ny, ns)

    lx = 1.0
    ly = px / py
    lr = math.sqrt(2.0) * px / pr

    lower = abs(lx - ly)
    upper = lx + ly
    feasible = lower <= lr <= upper

    if not feasible:
        return False, {
            "px": float(px),
            "py": float(py),
            "pr": float(pr),
            "lx": lx,
            "ly": ly,
            "lr": lr,
            "lower": lower,
            "upper": upper,
        }

    cos_theta = (lx * lx + ly * ly - lr * lr) / (2.0 * lx * ly)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.acos(cos_theta)

    x_vec = np.array([lx, 0.0])
    y_vec = np.array([ly * math.cos(theta), ly * math.sin(theta)])
    r_vec = y_vec - x_vec

    x_eff = 0.5 * px * x_vec
    y_eff = 0.5 * py * y_vec
    d_eff = y_eff - x_eff

    nx_eff = float(np.linalg.norm(x_eff))
    ny_eff = float(np.linalg.norm(y_eff))
    nd_eff = float(np.linalg.norm(d_eff))

    dot = float(np.dot(x_eff, y_eff))
    denom = max(nx_eff * ny_eff, 1e-15)
    cos_eff = max(-1.0, min(1.0, dot / denom))
    angle_eff = math.degrees(math.acos(cos_eff))

    return True, {
        "px": float(px),
        "py": float(py),
        "pr": float(pr),
        "lx": lx,
        "ly": ly,
        "lr": lr,
        "theta_deg": math.degrees(theta),
        "x_vec": x_vec,
        "y_vec": y_vec,
        "r_vec": r_vec,
        "x_eff": x_eff,
        "y_eff": y_eff,
        "d_eff": d_eff,
        "nx_eff": nx_eff,
        "ny_eff": ny_eff,
        "nd_eff": nd_eff,
        "angle_eff": angle_eff,
    }


def draw_case(nx: int, ny: int, ns: int, out_dir: Path) -> tuple[bool, str]:
    ok, data = compensating_triangle(nx, ny, ns)
    if not ok:
        return False, (
            f"({nx},{ny},{ns}) infeasible: l_r={data['lr']:.4f} not in "
            f"[{data['lower']:.4f}, {data['upper']:.4f}]"
        )

    px = int(data["px"])
    py = int(data["py"])
    pr = int(data["pr"])

    x_vec = data["x_vec"]
    y_vec = data["y_vec"]
    r_vec = data["r_vec"]
    x_eff = data["x_eff"]
    y_eff = data["y_eff"]
    d_eff = data["d_eff"]

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

    # Panel 1: local compensating triangle.
    ax = axes[0]
    O = np.array([0.0, 0.0])
    X = x_vec
    Y = y_vec
    ax.plot([O[0], X[0]], [O[1], X[1]], color="#1f77b4", lw=3, label="x")
    ax.plot([O[0], Y[0]], [O[1], Y[1]], color="#d62728", lw=3, label="y")
    ax.plot([X[0], Y[0]], [X[1], Y[1]], color="#2ca02c", lw=3, label="y-x")
    ax.scatter([0, X[0], Y[0]], [0, X[1], Y[1]], color="black", s=24)
    ax.set_title("Local Compensating Triangle")
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="upper right")

    txt_local = (
        f"|x|={data['lx']:.4f}\n"
        f"|y|={data['ly']:.4f}\n"
        f"|y-x|={data['lr']:.4f}\n"
        f"theta={data['theta_deg']:.2f} deg"
    )
    ax.text(
        0.03,
        0.97,
        txt_local,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="#f7f7f7", ec="#cccccc"),
    )

    # Panel 2: period-weighting map.
    ax2 = axes[1]
    ax2.arrow(0, 0, x_vec[0], x_vec[1], width=0.01, color="#1f77b4", alpha=0.8)
    ax2.arrow(0, 0, y_vec[0], y_vec[1], width=0.01, color="#d62728", alpha=0.8)
    ax2.arrow(0, 0, x_eff[0], x_eff[1], width=0.02, color="#1f77b4")
    ax2.arrow(0, 0, y_eff[0], y_eff[1], width=0.02, color="#d62728")
    ax2.set_title("Period Weighting")
    ax2.set_aspect("equal")
    ax2.grid(alpha=0.2)

    txt_map = (
        f"p_x={px}, p_y={py}, p_y-x={pr}\n"
        f"X_eff=(p_x/2)x\n"
        f"Y_eff=(p_y/2)y"
    )
    ax2.text(
        0.03,
        0.97,
        txt_map,
        transform=ax2.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="#f7f7f7", ec="#cccccc"),
    )

    # Panel 3: square-theory effective geometry.
    ax3 = axes[2]
    ax3.plot([0, x_eff[0]], [0, x_eff[1]], color="#1f77b4", lw=3, label="X_eff")
    ax3.plot([0, y_eff[0]], [0, y_eff[1]], color="#d62728", lw=3, label="Y_eff")
    ax3.plot([x_eff[0], y_eff[0]], [x_eff[1], y_eff[1]], color="#2ca02c", lw=3, label="Y_eff-X_eff")
    ax3.scatter([0, x_eff[0], y_eff[0]], [0, x_eff[1], y_eff[1]], color="black", s=24)
    ax3.set_title("Effective Square-Theory Geometry")
    ax3.set_aspect("equal")
    ax3.grid(alpha=0.2)
    ax3.legend(fontsize=8, loc="upper right")

    ratio = data["nd_eff"] / max(data["nx_eff"], 1e-15)
    txt_eff = (
        f"|X_eff|={data['nx_eff']:.4f}\n"
        f"|Y_eff|={data['ny_eff']:.4f}\n"
        f"angle={data['angle_eff']:.4f} deg\n"
        f"|Y_eff-X_eff|/|X_eff|={ratio:.6f}"
    )
    ax3.text(
        0.03,
        0.97,
        txt_eff,
        transform=ax3.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="#f7f7f7", ec="#cccccc"),
    )

    fig.suptitle(
        f"Case Nx={nx}, Ny={ny}, Ns={ns}: compensating triangle -> square-theory effective vectors",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    stem = f"square_theory_case_Nx{nx}_Ny{ny}_Ns{ns}"
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)

    return True, f"({nx},{ny},{ns}) -> {png.name}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases",
        default="3x4x2,4x4x2,6x5x4,6x6x3,7x7x7,8x8x4",
        help="Comma-separated Nx x Ny x Ns cases",
    )
    parser.add_argument(
        "--out-dir",
        default="thinkTwist/square_theory_diagrams",
        help="Output directory",
    )
    args = parser.parse_args()

    cases = parse_cases(args.cases)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    done: list[str] = []
    skipped: list[str] = []

    for nx, ny, ns in cases:
        ok, msg = draw_case(nx, ny, ns, out_dir)
        if ok:
            done.append(msg)
        else:
            skipped.append(msg)

    report = out_dir / "README.txt"
    lines = [
        "Compensating Triangle -> Square Theory Diagrams",
        "==============================================",
        "",
        f"Requested cases: {args.cases}",
        f"Generated: {len(done)}",
    ]
    lines.extend(done)
    lines.append("")
    lines.append(f"Skipped (infeasible): {len(skipped)}")
    lines.extend(skipped)
    report.write_text("\n".join(lines) + "\n")

    print(f"Wrote {report}")
    for d in done:
        print(d)
    for s in skipped:
        print(f"SKIP {s}")


if __name__ == "__main__":
    main()
