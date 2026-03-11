#!/usr/bin/env python3
"""Plot full all-to-all two-point data as a spike-centered 3D manifold."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def lattice_to_xy(m: float, n: float):
    # Triangular basis e1=(1,0), e2=(1/2,sqrt(3)/2)
    return (m + 0.5 * n, (np.sqrt(3.0) / 2.0) * n)


def to_uv(m: float, n: float, Lx: int, Ly: int, Tx: int, Ty: int):
    # [m,n]^T = a*(Lx,Ty) + b*(Tx,-Ly)
    N = float(Lx * Ly + Tx * Ty)
    a = (Ly * m + Tx * n) / N
    b = (Ty * m - Lx * n) / N
    return a, b


def from_uv(a: float, b: float, Lx: int, Ly: int, Tx: int, Ty: int):
    m = a * Lx + b * Tx
    n = a * Ty - b * Ly
    return m, n


def load_full(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append(
                (
                    int(d_s),
                    int(m_s),
                    int(n_s),
                    float(c_s),
                    float(e_s),
                    float(cc_s),
                    float(ec_s),
                )
            )
    return rows


def centered_arrays(rows, Lx, Ly, Tx, Ty, connected=True):
    # Center around the largest raw correlator spike.
    _d0, m0, n0, _c0, _e0, _cc0, _ec0 = max(rows, key=lambda r: r[3])

    xs, ys, zs = [], [], []
    for _d, m, n, c, _e, cc, _ec in rows:
        dm = m - m0
        dn = n - n0

        # Half-wrap in both periods so distances are centered around spike.
        a, b = to_uv(dm, dn, Lx, Ly, Tx, Ty)
        a -= np.round(a)
        b -= np.round(b)
        mw, nw = from_uv(a, b, Lx, Ly, Tx, Ty)

        x, y = lattice_to_xy(mw, nw)
        xs.append(x)
        ys.append(y)
        zs.append(cc if connected else c)

    return np.array(xs), np.array(ys), np.array(zs)


def draw_centered_cell_boundary(ax, Lx, Ly, Tx, Ty, z_ref):
    # Corners of centered fundamental cell: +/- v/2 +/- u/2
    corners = [
        (-(Lx + Tx) * 0.5, (-Ty + Ly) * 0.5),
        ((Lx - Tx) * 0.5, (Ty + Ly) * 0.5),
        ((Lx + Tx) * 0.5, (Ty - Ly) * 0.5),
        ((-Lx + Tx) * 0.5, (-Ty - Ly) * 0.5),
        (-(Lx + Tx) * 0.5, (-Ty + Ly) * 0.5),
    ]
    pts = [lattice_to_xy(m, n) for m, n in corners]
    xb = [p[0] for p in pts]
    yb = [p[1] for p in pts]
    zb = [z_ref] * len(pts)
    ax.plot(xb, yb, zb, color="black", lw=1.4, alpha=0.9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True)
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--connected", action="store_true", default=True)
    args = ap.parse_args()

    rows = load_full(Path(args.full))
    if not rows:
        raise SystemExit("No data rows found in full all-to-all file.")

    x, y, z = centered_arrays(
        rows,
        Lx=args.L_x,
        Ly=args.L_y,
        Tx=args.T_x,
        Ty=args.T_y,
        connected=args.connected,
    )

    fig = plt.figure(figsize=(7.2, 5.8), dpi=220)
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    tri = mtri.Triangulation(x, y)
    surf = ax.plot_trisurf(
        tri,
        z,
        cmap="viridis",
        linewidth=0.08,
        antialiased=True,
        alpha=0.96,
    )
    fig.colorbar(surf, ax=ax, shrink=0.63, pad=0.08, label="C_conn")

    draw_centered_cell_boundary(ax, args.L_x, args.L_y, args.T_x, args.T_y, float(np.min(z)))

    ax.set_title("Spike-Centered Connected Manifold (All-to-All)")
    ax.set_xlabel("x displacement")
    ax.set_ylabel("y displacement")
    ax.set_zlabel("C_conn")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
