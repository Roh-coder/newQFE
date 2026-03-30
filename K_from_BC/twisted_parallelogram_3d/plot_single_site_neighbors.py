#!/usr/bin/env python3
"""Plot a single site and its nearest neighbors for the triangular x time lattice.

This diagram is intended for presentations:
- 6 in-plane nearest neighbors on the triangular lattice
- 2 temporal nearest neighbors (+t and -t)
"""

from __future__ import annotations

import argparse
import math
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


Vec3 = Tuple[float, float, float]


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(c: float, v: Vec3) -> Vec3:
    return (c * v[0], c * v[1], c * v[2])


def build_geometry(z_scale: float = 1.10) -> Dict[str, Vec3]:
    # 2D triangular basis vectors embedded in the x-y plane.
    e1 = (1.0, 0.0, 0.0)
    e2 = (0.5, math.sqrt(3.0) / 2.0, 0.0)
    e2me1 = (e2[0] - e1[0], e2[1] - e1[1], 0.0)

    # Time direction orthogonal to the triangular plane in this schematic.
    et = (0.0, 0.0, z_scale)

    points = {
        "0": (0.0, 0.0, 0.0),
        "+e1": e1,
        "-e1": scale(-1.0, e1),
        "+e2": e2,
        "-e2": scale(-1.0, e2),
        "+(e2-e1)": e2me1,
        "-(e2-e1)": scale(-1.0, e2me1),
        "+t": et,
        "-t": scale(-1.0, et),
    }
    return points


def draw(out_path: str, title: str) -> None:
    pts = build_geometry()

    fig = plt.figure(figsize=(10, 5.2))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)

    center = pts["0"]
    planar_labels: List[str] = [
        "+e1",
        "-e1",
        "+e2",
        "-e2",
        "+(e2-e1)",
        "-(e2-e1)",
    ]
    temporal_labels: List[str] = ["+t", "-t"]

    planar_color = "#1f77b4"
    temporal_color = "#d62728"
    center_color = "#111111"

    # 3D panel.
    for lab in planar_labels:
        p = pts[lab]
        ax3d.plot([center[0], p[0]], [center[1], p[1]], [center[2], p[2]],
                  color=planar_color, linewidth=2.0, alpha=0.9)
    for lab in temporal_labels:
        p = pts[lab]
        ax3d.plot([center[0], p[0]], [center[1], p[1]], [center[2], p[2]],
                  color=temporal_color, linewidth=2.4, alpha=0.95)

    ax3d.scatter([center[0]], [center[1]], [center[2]], s=90, c=center_color, zorder=10)

    for lab in planar_labels:
        p = pts[lab]
        ax3d.scatter([p[0]], [p[1]], [p[2]], s=65, c=planar_color)
        ax3d.text(p[0] * 1.10, p[1] * 1.10, p[2], lab, fontsize=9, color=planar_color)

    for lab in temporal_labels:
        p = pts[lab]
        ax3d.scatter([p[0]], [p[1]], [p[2]], s=70, c=temporal_color)
        ax3d.text(p[0] * 1.07, p[1] * 1.07, p[2] * 1.03, lab, fontsize=9, color=temporal_color)

    ax3d.text(0.02, 0.02, 0.03, "site 0", fontsize=9, color=center_color)
    ax3d.set_title("3D view", fontsize=11)
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("t")
    ax3d.view_init(elev=20, azim=-55)
    ax3d.set_box_aspect((1.3, 1.2, 1.0))

    # 2D projection panel.
    for lab in planar_labels:
        p = pts[lab]
        ax2d.plot([center[0], p[0]], [center[1], p[1]], color=planar_color, linewidth=2.0)
    ax2d.scatter([center[0]], [center[1]], s=90, c=center_color, zorder=5)
    for lab in planar_labels:
        p = pts[lab]
        ax2d.scatter([p[0]], [p[1]], s=65, c=planar_color)
        ax2d.text(p[0] * 1.14, p[1] * 1.14, lab, fontsize=9, color=planar_color,
                  ha="center", va="center")

    # Show temporal neighbors in projection with dashed arrows to indicate out-of-plane links.
    for sgn, lab in [(1.0, "+t"), (-1.0, "-t")]:
        ax2d.annotate(
            "",
            xy=(0.0, 0.0),
            xytext=(0.35 * sgn, 0.35 * sgn),
            arrowprops=dict(arrowstyle="->", linestyle="--", lw=1.8, color=temporal_color),
        )
        ax2d.text(0.42 * sgn, 0.42 * sgn, lab, color=temporal_color, fontsize=9,
                  ha="center", va="center")

    ax2d.set_title("Triangular plane projection", fontsize=11)
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.set_xlabel("x")
    ax2d.set_ylabel("y")
    ax2d.grid(alpha=0.25)

    fig.suptitle(title, fontsize=13)
    fig.text(
        0.5,
        0.01,
        "Nearest neighbors of one site: 6 planar (+/-e1, +/-e2, +/-[e2-e1]) and 2 temporal (+/-t).",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(out_path, dpi=220)
    print(f"Saved: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default="single_site_neighbors.png",
        help="Output image path.",
    )
    parser.add_argument(
        "--title",
        default="Triangular x Time Ising Lattice: Single Site Neighborhood",
        help="Figure title.",
    )
    args = parser.parse_args()
    draw(args.out, args.title)


if __name__ == "__main__":
    main()
