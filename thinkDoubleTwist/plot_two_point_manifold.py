#!/usr/bin/env python3
"""Plot directional two-point data as a 3D manifold-like surface."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def load_typed(path: Path):
    data = {}
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            t_s, r_s, c_s, e_s = s.split()
            t = int(t_s)
            r = int(r_s)
            c = float(c_s)
            e = float(e_s)
            data.setdefault(t, []).append((r, c, e))

    for t in data:
        data[t].sort(key=lambda x: x[0])
    return data


def lattice_to_xy(m: float, n: float):
    # Triangular basis e1=(1,0), e2=(1/2,sqrt(3)/2)
    return (m + 0.5 * n, (np.sqrt(3.0) / 2.0) * n)


def reduce_to_cell(m: float, n: float, Lx: int, Ly: int, Tx: int, Ty: int):
    # Solve [m,n]^T = a*v + b*u with v=(Lx,Ty), u=(Tx,-Ly), then wrap a,b to [0,1).
    N = float(Lx * Ly + Tx * Ty)
    a = (Ly * m + Tx * n) / N
    b = (Ty * m - Lx * n) / N
    af = a - np.floor(a)
    bf = b - np.floor(b)
    mr = af * Lx + bf * Tx
    nr = af * Ty - bf * Ly
    return mr, nr


def channel_points(typed, type_map, Lx, Ly, Tx, Ty):
    # Lattice-coordinate directions for type channels.
    # 0/4: +e1 -> (1,0), 1/5: +e2 -> (0,1), 2/6: +(e2-e1) -> (-1,1)
    dirs = {
        type_map[0]: (1, 0),
        type_map[1]: (0, 1),
        type_map[2]: (-1, 1),
    }

    xs, ys, zs = [], [], []
    for t in (type_map[0], type_map[1], type_map[2]):
        if t not in typed or not typed[t]:
            continue
        dm, dn = dirs[t]
        for r, c, _e in typed[t]:
            m = dm * r
            n = dn * r
            mr, nr = reduce_to_cell(m, n, Lx, Ly, Tx, Ty)
            x, y = lattice_to_xy(mr, nr)
            xs.append(x)
            ys.append(y)
            zs.append(c)

    if not xs:
        return np.array([]), np.array([]), np.array([])
    return np.array(xs), np.array(ys), np.array(zs)


def draw_cell_boundary(ax, Lx, Ly, Tx, Ty, z_ref):
    corners = [(0.0, 0.0), (Lx, Ty), (Lx + Tx, Ty - Ly), (Tx, -Ly), (0.0, 0.0)]
    pts = [lattice_to_xy(m, n) for m, n in corners]
    xb = [p[0] for p in pts]
    yb = [p[1] for p in pts]
    zb = [z_ref] * len(pts)
    ax.plot(xb, yb, zb, color="black", lw=1.5, alpha=0.9)


def plot_surface(ax, x, y, z, title, Lx, Ly, Tx, Ty):
    # Show all raw points first (scatter), then fill a manifold through them.
    ax.scatter(x, y, z, color="#ff3b30", s=28, alpha=0.95, depthshade=False)
    tri = mtri.Triangulation(x, y)
    surf = ax.plot_trisurf(tri, z, cmap="viridis", linewidth=0.10, antialiased=True, alpha=0.42)
    draw_cell_boundary(ax, Lx, Ly, Tx, Ty, float(np.min(z)))
    ax.set_title(title)
    ax.set_xlabel("x displacement")
    ax.set_ylabel("y displacement")
    ax.set_zlabel("C(r)")
    return surf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--typed", required=True, help="Typed two-point .dat file")
    ap.add_argument("--output", required=True, help="Output image path")
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    args = ap.parse_args()

    typed = load_typed(Path(args.typed))

    fig = plt.figure(figsize=(12, 5.2), dpi=180)

    # Left: raw manifold (types 0,1,2)
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    x, y, z = channel_points(typed, (0, 1, 2), args.L_x, args.L_y, args.T_x, args.T_y)
    if len(x) > 0:
        plot_surface(ax1, x, y, z, "Raw Two-Point Manifold (Parallelogram Cell)",
                     args.L_x, args.L_y, args.T_x, args.T_y)
    else:
        ax1.set_title("Raw Two-Point Manifold (missing types 0/1/2)")

    # Right: connected manifold (types 4,5,6), if present.
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    xc, yc, zc = channel_points(typed, (4, 5, 6), args.L_x, args.L_y, args.T_x, args.T_y)
    if len(xc) > 0:
        plot_surface(ax2, xc, yc, zc, "Connected Two-Point Manifold (Parallelogram Cell)",
                     args.L_x, args.L_y, args.T_x, args.T_y)
    else:
        ax2.set_title("Connected Two-Point Manifold (missing types 4/5/6)")

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
