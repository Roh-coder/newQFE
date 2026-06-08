#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize the L=4 base displacement grid and the v/u/w directions used in the reweighting plots."
    )
    parser.add_argument(
        "--output",
        default=str(RESPONSIBLE_ROOT / "reweighting" / "base_grid_boundary_directions_L4.png"),
    )
    return parser.parse_args()


def _draw_direction(ax: plt.Axes, end: tuple[float, float], *, color: str, label: str) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=(0.0, 0.0),
        arrowprops={"arrowstyle": "->", "lw": 2.2, "color": color},
        zorder=2,
    )
    ax.text(end[0] + 0.08, end[1] + 0.08, label, color=color, fontsize=12, weight="bold")


def main() -> None:
    args = parse_args()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.4, 8.0))

    # Draw the 4x4 base cell in displacement coordinates (m, n).
    for coord in range(5):
        ax.axhline(coord, color="#d4d4d8", lw=0.9, zorder=0)
        ax.axvline(coord, color="#d4d4d8", lw=0.9, zorder=0)

    lattice_points_x: list[float] = []
    lattice_points_y: list[float] = []
    for m_coord in range(4):
        for n_coord in range(4):
            lattice_points_x.append(float(m_coord))
            lattice_points_y.append(float(n_coord))
    ax.scatter(lattice_points_x, lattice_points_y, s=18, color="#9ca3af", zorder=1)

    selected_points = {
        "origin": {"xy": (0.0, 0.0), "color": "#111827"},
        "q_v = (1,0) = (L/4,0)": {"xy": (1.0, 0.0), "color": "#7c2d12"},
        "q_u = (0,1) = (0,L/4)": {"xy": (0.0, 1.0), "color": "#7c3aed"},
        "q_w = (1,1) = (L/4,L/4)": {"xy": (1.0, 1.0), "color": "#ca8a04"},
        "unused = (1,3) = (L/4,3L/4)": {"xy": (1.0, 3.0), "color": "#dc2626"},
        "mid_v = (2,0) = (L/2,0)": {"xy": (2.0, 0.0), "color": "#1d4ed8"},
        "mid_u = (0,2) = (0,L/2)": {"xy": (0.0, 2.0), "color": "#047857"},
        "mid_w = (2,2) = (L/2,L/2)": {"xy": (2.0, 2.0), "color": "#b45309"},
    }

    for label, meta in selected_points.items():
        x_coord, y_coord = meta["xy"]
        color = meta["color"]
        marker_size = 90 if label.startswith("unused") else 60
        ax.scatter([x_coord], [y_coord], s=marker_size, color=color, edgecolor="white", linewidth=0.9, zorder=3)
        text_dx = 0.12 if x_coord < 2.0 else 0.15
        text_dy = -0.18 if label.startswith("unused") else (0.12 if y_coord < 2.0 else 0.15)
        ax.text(x_coord + text_dx, y_coord + text_dy, label, fontsize=10, color=color)

    _draw_direction(ax, (2.0, 0.0), color="#1d4ed8", label="v direction")
    _draw_direction(ax, (0.0, 2.0), color="#047857", label="u direction")
    _draw_direction(ax, (2.0, 2.0), color="#b45309", label="w direction")

    ax.text(
        2.55,
        0.55,
        "Use (L/4, 3L/4)\n as the unused interior normalization point",
        fontsize=11,
        color="#92400e",
        bbox={"facecolor": "#fef3c7", "edgecolor": "#f59e0b", "boxstyle": "round,pad=0.35"},
    )

    fig.suptitle("4x4 Base Displacement Grid (L = 4)", fontsize=18, y=0.97)
    fig.text(
        0.5,
        0.935,
        "Selected sample points are at quarter and half displacements along v, u, and w; q_w lies on the w diagonal, while (L/4,3L/4) is the unused point.",
        ha="center",
        va="top",
        fontsize=10.5,
        color="#444444",
    )
    ax.set_xlabel("m displacement")
    ax.set_ylabel("n displacement")
    ax.set_xlim(-0.45, 4.1)
    ax.set_ylim(-0.45, 4.1)
    ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_yticks([0, 1, 2, 3, 4])
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.90])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()