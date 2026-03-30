#!/usr/bin/env python3
"""Plot selected twisted-parallelogram triangles against a 4-5-6 reference.

Each candidate (Lx, Ly, Tx, Ty) defines lattice vectors
  v = (Lx, Ty), u = (Tx, -Ly), w = -(v + u)
in integer triangular-lattice coordinates (m, n).

We convert (m, n) to Euclidean coordinates using basis
  e1 = (1, 0), e2 = (1/2, sqrt(3)/2),
then draw candidate triangles normalized so the shortest side is 4,
overlaid with an ideal 4-5-6 triangle.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class Candidate:
    lx: int
    ly: int
    tx: int
    ty: int


def tri_len_sq(m: int, n: int) -> int:
    return m * m + n * n + m * n


def mn_to_xy(m: float, n: float) -> np.ndarray:
    return np.array([m + 0.5 * n, 0.5 * math.sqrt(3.0) * n], dtype=float)


def reference_triangle() -> np.ndarray:
    # Fix AB=4, AC=5, BC=6 with A=(0,0), B=(4,0), C above x-axis.
    ax, ay = 0.0, 0.0
    bx, by = 4.0, 0.0
    cx = (5.0 * 5.0 - 6.0 * 6.0 + 4.0 * 4.0) / (2.0 * 4.0)
    cy = math.sqrt(max(0.0, 5.0 * 5.0 - cx * cx))
    return np.array([[ax, ay], [bx, by], [cx, cy]], dtype=float)


def candidate_vertices(c: Candidate) -> np.ndarray:
    # Triangle vertices: A=0, B=v, C=v+u.
    v_mn = (c.lx, c.ty)
    u_mn = (c.tx, -c.ly)
    a = np.array([0.0, 0.0], dtype=float)
    b = mn_to_xy(*v_mn)
    cpt = mn_to_xy(v_mn[0] + u_mn[0], v_mn[1] + u_mn[1])
    return np.array([a, b, cpt], dtype=float)


def edge_lengths(verts: np.ndarray) -> tuple[float, float, float]:
    a, b, c = verts
    ab = float(np.linalg.norm(b - a))
    bc = float(np.linalg.norm(c - b))
    ca = float(np.linalg.norm(a - c))
    return ab, bc, ca


def _transform_from_order(verts: np.ndarray, order: tuple[int, int, int]) -> np.ndarray:
    # Given vertex order (A, B, C), map A->(0,0), AB->+x and |AB|->4.
    out = verts[list(order)].copy()
    out -= out[0]

    ab = out[1] - out[0]
    ab_len = float(np.linalg.norm(ab))
    if ab_len <= 0.0:
        return out

    out *= 4.0 / ab_len
    ab = out[1] - out[0]
    angle = math.atan2(ab[1], ab[0])
    c = math.cos(-angle)
    s = math.sin(-angle)
    rot = np.array([[c, -s], [s, c]], dtype=float)
    out = out @ rot.T

    # Keep third vertex above the base for consistent up/down orientation.
    if out[2, 1] < 0.0:
        out[:, 1] *= -1.0
    return out


def normalize_and_orient(verts: np.ndarray, ref: np.ndarray) -> np.ndarray:
    # Align shortest side with reference base (4) and longest side with reference long edge (6).
    # This removes mirrored edge assignments between candidate and 4-5-6 overlays.
    edges = [
        (0, 1, float(np.linalg.norm(verts[1] - verts[0]))),
        (1, 2, float(np.linalg.norm(verts[2] - verts[1]))),
        (2, 0, float(np.linalg.norm(verts[0] - verts[2]))),
    ]
    shortest = min(edges, key=lambda e: e[2])
    longest = max(edges, key=lambda e: e[2])

    short_set = {shortest[0], shortest[1]}
    long_set = {longest[0], longest[1]}
    shared = list(short_set & long_set)

    if len(shared) != 1:
        # Fallback for pathological near-degenerate cases.
        return _transform_from_order(verts, (0, 1, 2))

    b = shared[0]
    a = list(short_set - {b})[0]
    c_idx = list(long_set - {b})[0]

    # Two possible base directions; pick the one closest to the reference vertex layout.
    cand1 = _transform_from_order(verts, (a, b, c_idx))
    cand2 = _transform_from_order(verts, (b, a, c_idx))

    err1 = float(np.sum((cand1 - ref) ** 2))
    err2 = float(np.sum((cand2 - ref) ** 2))
    return cand1 if err1 <= err2 else cand2


def score_vs_456(c: Candidate) -> float:
    lv2 = tri_len_sq(c.lx, c.ty)
    lu2 = tri_len_sq(c.tx, -c.ly)
    lw2 = tri_len_sq(-c.lx - c.tx, c.ly - c.ty)
    target = (4, 5, 6)
    scales = [
        lv2 / float(target[0] * target[0]),
        lu2 / float(target[1] * target[1]),
        lw2 / float(target[2] * target[2]),
    ]
    mean = sum(scales) / 3.0
    return (max(scales) - min(scales)) / mean if mean else float("inf")


def draw_triangle(ax: plt.Axes, verts: np.ndarray, color: str, label: str, ls: str = "-") -> None:
    loop = np.vstack([verts, verts[0]])
    ax.plot(loop[:, 0], loop[:, 1], color=color, lw=2.0, ls=ls, label=label)


def make_plot(candidates: list[Candidate], out_path: str) -> None:
    ref = reference_triangle()

    cols = 4
    rows = 2
    fig, axes = plt.subplots(rows, cols, figsize=(16, 8), constrained_layout=True)
    flat = axes.ravel()

    for i, c in enumerate(candidates):
        ax = flat[i]
        cand = normalize_and_orient(candidate_vertices(c), ref)

        draw_triangle(ax, ref, "#666666", "ideal 4-5-6", ls="--")
        draw_triangle(ax, cand, "#1f77b4", "candidate")

        ab, bc, ca = edge_lengths(cand)
        score = score_vs_456(c)
        ax.set_title(f"Lx={c.lx}, Ly={c.ly}, Tx={c.tx}, Ty={c.ty}", fontsize=10)
        ax.text(
            0.02,
            0.98,
            f"norm edges: {ab:.3f}, {bc:.3f}, {ca:.3f}\nscore={score:.6f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc", alpha=0.9),
        )
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.25)

    if len(candidates) < len(flat):
        for j in range(len(candidates), len(flat)):
            flat[j].axis("off")

    flat[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("Candidate Twisted-Parallelogram Triangles vs Ideal 4-5-6", fontsize=14)
    fig.savefig(out_path, dpi=180)
    print(f"Wrote {out_path}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot selected candidate triangles vs ideal 4-5-6.")
    ap.add_argument(
        "--out",
        default="thinkDoubleTwist/triangle_candidates_vs_456.png",
        help="Output image path.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    # From previous search output (fixed v:u:w order ~ 4:5:6).
    candidates = [
        Candidate(13, 16, 3, -3),
        Candidate(10, 3, -13, -13),
        Candidate(3, 13, 16, 10),
        Candidate(23, 19, -10, -16),
        Candidate(16, 29, 19, 7),
        Candidate(26, 32, 6, -6),
    ]

    make_plot(candidates, args.out)


if __name__ == "__main__":
    main()
