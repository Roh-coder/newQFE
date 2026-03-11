#!/usr/bin/env python3
"""Draw an equilateral triangular lattice patch tiling the plane."""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt


def hex_norm(i: int, j: int) -> int:
    return max(abs(i), abs(j), abs(i + j))


def make_hex_sites(radius: int):
    sites = set()
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            if hex_norm(i, j) <= radius:
                sites.add((i, j))
    return sites


def xy(i: int, j: int, a: float):
    # Triangular-lattice basis: e1=(a,0), e2=(a/2, sqrt(3)a/2)
    return (a * (i + 0.5 * j), a * (math.sqrt(3.0) * 0.5 * j))


def draw_lattice(
    nx: int,
    ny: int,
    a: float,
    node_size: float,
    lw: float,
    l_x: int,
    l_y: int,
    t_x: int,
    t_y: int,
):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    params_label = ""
    v_label = ""
    u_label = ""
    w_label = ""

    # Build a hexagonal patch of triangular-lattice sites in axial (i,j) coords.
    # Radius is enlarged when needed so the full v/u/w construction fits.
    lx_req = max(1, int(l_x))
    ly_req = max(1, int(l_y))
    tx_req = int(t_x)
    ty_req = int(t_y)
    offsets = [
        (0, 0),
        (lx_req, ty_req),
        (tx_req, -ly_req),
        (lx_req + tx_req, ty_req - ly_req),
    ]
    needed_radius = max(hex_norm(i, j) for i, j in offsets) + 2
    base_radius = max(2, (max(nx, ny) + 1) // 2)
    radius = max(base_radius, needed_radius)
    sites = make_hex_sites(radius)

    # Draw each undirected edge once using +e1, +e2, and +(e2-e1).
    for i, j in sites:
        x0, y0 = xy(i, j, a)
        for ni, nj in ((i + 1, j), (i, j + 1), (i - 1, j + 1)):
            if (ni, nj) in sites:
                x1, y1 = xy(ni, nj, a)
                ax.plot([x0, x1], [y0, y1], color="gray", lw=lw)

    # Draw lattice sites.
    xs, ys = [], []
    for i, j in sites:
        x, y = xy(i, j, a)
        xs.append(x)
        ys.append(y)

    ax.scatter(xs, ys, s=node_size, color="tab:blue", zorder=3)

    # Draw vectors from a centered interior anchor in the hex patch.
    if radius >= 2:
        # Four-parameter form:
        # v = L_x e1 + T_y e2
        # u = T_x e1 - L_y e2
        lx_use = max(1, min(int(l_x), max(1, nx - 2)))
        ly_use = max(1, min(int(l_y), max(1, ny - 2)))
        tx_use = int(t_x)
        ty_use = int(t_y)

        # Centered anchor; radius was chosen to ensure this construction fits.
        i_mark = 0
        j_mark = 0

        xr, yr = xy(i_mark, j_mark, a)
        ax.scatter([xr], [yr], s=3.0 * node_size, color="red", zorder=5)

        xv_end, yv_end = xy(i_mark + lx_use, j_mark + ty_use, a)
        xu_end, yu_end = xy(i_mark + tx_use, j_mark - ly_use, a)

        # Opposite corner for the completed parallelogram.
        xw_end, yw_end = xy(i_mark + lx_use + tx_use, j_mark + ty_use - ly_use, a)

        ax.annotate(
            "",
            xy=(xv_end, yv_end),
            xytext=(xr, yr),
            arrowprops=dict(arrowstyle="->", lw=2.4, color="red"),
            zorder=6,
        )
        ax.annotate(
            "",
            xy=(xu_end, yu_end),
            xytext=(xr, yr),
            arrowprops=dict(arrowstyle="->", lw=2.4, color="blue"),
            zorder=6,
        )

        # Translated copies to complete the parallelogram.
        ax.annotate(
            "",
            xy=(xw_end, yw_end),
            xytext=(xv_end, yv_end),
            arrowprops=dict(arrowstyle="->", lw=2.0, color="blue", alpha=0.9),
            zorder=6,
        )
        ax.annotate(
            "",
            xy=(xw_end, yw_end),
            xytext=(xu_end, yu_end),
            arrowprops=dict(arrowstyle="->", lw=2.0, color="red", alpha=0.9),
            zorder=6,
        )

        # Black vector on the other parallelogram diagonal with flipped orientation.
        ax.annotate(
            "",
            xy=(xr, yr),
            xytext=(xw_end, yw_end),
            arrowprops=dict(arrowstyle="->", lw=2.2, color="black"),
            zorder=7,
        )
        params_label = f"L_x={lx_use},  L_y={ly_use},  T_x={tx_use},  T_y={ty_use}"
        v_label = f"v = {lx_use}e_1 + {ty_use}e_2"
        u_label = f"u = {tx_use}e_1 - {ly_use}e_2"
        w_e1 = -(lx_use + tx_use)
        w_e2 = (ly_use - ty_use)
        w_label = f"w = -(v + u) = {w_e1}e_1 + {w_e2}e_2"

    ax.set_aspect("equal")
    ax.set_axis_off()
    y_text = 0.03
    text_box = dict(facecolor="white", edgecolor="none", alpha=0.9, pad=0.2)
    if v_label:
        # Reserve extra room below the lattice so labels never overlap geometry.
        fig.subplots_adjust(bottom=0.28)
        y_text = 0.19
        if params_label:
            fig.text(
                0.14,
                y_text,
                params_label,
                ha="left",
                va="bottom",
                fontsize=10,
                color="black",
                bbox=text_box,
            )
        fig.text(
            0.14,
            max(0.01, y_text - 0.04),
            v_label,
            ha="left",
            va="bottom",
            fontsize=10,
            color="red",
            bbox=text_box,
        )
    if u_label:
        fig.text(
            0.14,
            max(0.01, y_text - 0.08),
            u_label,
            ha="left",
            va="bottom",
            fontsize=10,
            color="blue",
            bbox=text_box,
        )
    if w_label:
        fig.text(
            0.14,
            max(0.01, y_text - 0.12),
            w_label,
            ha="left",
            va="bottom",
            fontsize=10,
            color="black",
            bbox=text_box,
        )
    fig.tight_layout(rect=(0.0, 0.16, 1.0, 1.0))
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nx", type=int, default=14, help="Hex-size control (used with ny to set base hex radius)")
    p.add_argument("--ny", type=int, default=10, help="Hex-size control (used with nx to set base hex radius)")
    p.add_argument("--a", type=float, default=1.0, help="Lattice spacing")
    p.add_argument("--node-size", type=float, default=10.0)
    p.add_argument("--line-width", type=float, default=0.8)
    p.add_argument("--L_x", type=int, default=4,
                   help="Regular lattice length in x direction (coefficient of e1 in v)")
    p.add_argument("--L_y", type=int, default=3,
                   help="Regular lattice length in left/Y direction (magnitude in u)")
    p.add_argument("--T_x", type=int, default=1,
                   help="Twist in x direction (coefficient of e1 in u)")
    p.add_argument("--T_y", type=int, default=1,
                   help="Twist in left/Y direction (coefficient of e2 in v)")
    p.add_argument("--output", default="thinkDoubleTwist/equilateral_lattice.png")
    args = p.parse_args()

    fig = draw_lattice(
        args.nx,
        args.ny,
        args.a,
        args.node_size,
        args.line_width,
        args.L_x,
        args.L_y,
        args.T_x,
        args.T_y,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
