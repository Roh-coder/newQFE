#!/usr/bin/env python3
"""Piecewise-linear interpolated heatmap of the all-to-all connected correlator,
tiled 2×2 on the 4-5-6 twisted-parallelogram torus.

Interpolation method: Delaunay triangulation + tricontourf (piecewise linear).
Tiling: 4 copies of the fundamental cell (ia ∈ {0,1}, ib ∈ {0,1}).

Period vectors in (m,n) lattice coords (from the BC code convention):
  v = (Lx,  Ty)    – boundary gluing along x
  u = (Tx, -Ly)    – boundary gluing along y

Physical embedding (triangular lattice):
  x = m + 0.5*n,   y = (√3/2)*n
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


# ── helpers ──────────────────────────────────────────────────────────────────

def lattice_to_xy(m: np.ndarray, n: np.ndarray):
    return m + 0.5 * n, (np.sqrt(3.0) / 2.0) * n


def load_data(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            rows.append((int(p[1]), int(p[2]), float(p[4]), float(p[5])))
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]   # m, n, corr, corr_conn


def parse_geometry(dirname: str):
    mt = re.search(r"Lx(\d+)_Ly(\d+)_Tx(-?\d+)_Ty(-?\d+)_Nt(\d+)", dirname)
    if not mt:
        return None
    keys = ("Lx", "Ly", "Tx", "Ty", "Nt")
    return {k: int(v) for k, v in zip(keys, mt.groups())}


def tile_2x2(m, n, z, Lx, Ly, Tx, Ty, floor_log10: float):
    """Return tiled (x_all, y_all, log10_z_all) for a 2×2 grid of copies."""
    vm, vn = Lx, Ty
    um, un = Tx, -Ly

    # Replace non-positive with a floor value before log
    z_safe = np.where(z > 0, z, np.nan)
    logz = np.where(np.isfinite(z_safe), np.log10(z_safe), floor_log10)

    xs, ys, zs = [], [], []
    for ia in range(2):
        for ib in range(2):
            mi = m + ia * vm + ib * um
            ni = n + ia * vn + ib * un
            xi, yi = lattice_to_xy(mi, ni)
            xs.append(xi)
            ys.append(yi)
            zs.append(logz)

    return np.concatenate(xs), np.concatenate(ys), np.concatenate(zs)


def cell_outline(Lx, Ly, Tx, Ty):
    """Corner (m,n) coords of the fundamental-cell parallelogram in (m,n)."""
    vm, vn = Lx, Ty
    um, un = Tx, -Ly
    corners_mn = np.array([[0, 0], [vm, vn], [vm + um, vn + un], [um, un], [0, 0]], float)
    return corners_mn


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="heatmap_tiled.png")
    ap.add_argument("--levels", type=int, default=48, help="tricontourf contour levels")
    ap.add_argument("--vmin", type=float, default=None, help="log10 colour scale min")
    ap.add_argument("--vmax", type=float, default=None, help="log10 colour scale max")
    ap.add_argument("--cmap", default="inferno")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    files = sorted(glob.glob(str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_all_to_all.dat")))
    if not files:
        raise SystemExit(f"No data files found under {datadir}")

    # ── collect all datasets and shared colour limits ────────────────────────
    datasets = []
    all_logz = []
    for fpath in files:
        fpath = Path(fpath)
        geo = parse_geometry(fpath.parent.name)
        if geo is None:
            continue
        m, n, _, conn = load_data(fpath)
        pos = conn[conn > 0]
        if pos.size:
            all_logz.append(np.log10(pos))
        datasets.append((geo, m, n, conn, fpath))

    if not datasets:
        raise SystemExit("No valid datasets.")

    all_logz = np.concatenate(all_logz)
    floor = float(np.percentile(all_logz, 1)) - 1.0
    vmin = args.vmin if args.vmin is not None else float(np.percentile(all_logz, 2))
    vmax = args.vmax if args.vmax is not None else float(np.percentile(all_logz, 99))

    cmap = plt.get_cmap(args.cmap)

    n_panels = len(datasets)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 7.5))
    if n_panels == 1:
        axes = [axes]

    for ax, (geo, m, n, conn, fpath) in zip(axes, datasets):
        Lx, Ly, Tx, Ty, Nt = geo["Lx"], geo["Ly"], geo["Tx"], geo["Ty"], geo["Nt"]
        s = Lx // 39

        x_all, y_all, logz_all = tile_2x2(m, n, conn, Lx, Ly, Tx, Ty, floor)

        # ── piecewise-linear interpolation via Delaunay triangulation ────────
        tri = mtri.Triangulation(x_all, y_all)

        # Remove triangles whose edges span more than ~3× the mean edge length
        # (avoids long "ghost" triangles at tile boundaries)
        verts = np.stack([x_all[tri.triangles], y_all[tri.triangles]], axis=-1)
        # edge vectors
        e01 = verts[:, 1] - verts[:, 0]
        e12 = verts[:, 2] - verts[:, 1]
        e20 = verts[:, 0] - verts[:, 2]
        max_edge2 = np.maximum.reduce([
            (e01**2).sum(1), (e12**2).sum(1), (e20**2).sum(1)
        ])
        threshold = (3.5 * (Lx / 39))**2    # ~3.5 lattice spacings at s=1
        tri.set_mask(max_edge2 > threshold)

        cntr = ax.tricontourf(tri, logz_all, levels=args.levels,
                              vmin=vmin, vmax=vmax, cmap=cmap, extend="both")

        # ── draw parallelogram outlines for each of the 4 copies ─────────────
        vm, vn = Lx, Ty
        um, un = Tx, -Ly
        base_mn = cell_outline(Lx, Ly, Tx, Ty)
        for ia in range(2):
            for ib in range(2):
                shift = np.array([ia * vm + ib * um, ia * vn + ib * un], float)
                poly_mn = base_mn + shift
                xb, yb = lattice_to_xy(poly_mn[:, 0], poly_mn[:, 1])
                ax.plot(xb, yb, color="white", lw=1.4, alpha=0.85, zorder=5)
                ax.plot(xb, yb, color="black", lw=0.5,  alpha=0.55, zorder=5)

        # mark the shared origin
        ax.scatter([0], [0], c="cyan", s=30, zorder=10, label="origin")

        ax.set_aspect("equal")
        ax.set_title(f"s={s}  (Lx={Lx}, Nt={Nt})  ×2×2 tiling", fontsize=11)
        ax.set_xlabel("x (lattice units)", fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel("y (lattice units)", fontsize=9)
        ax.tick_params(labelsize=8)

    # ── shared colorbar ───────────────────────────────────────────────────────
    fig.subplots_adjust(right=0.87, wspace=0.06)
    cbar_ax = fig.add_axes([0.89, 0.12, 0.018, 0.74])
    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cbar_ax, extend="both")
    cb.set_label(r"$\log_{10}\,G_{\rm conn}(r)$", fontsize=10)

    fig.suptitle(
        "Equal-time connected spatial correlator — 4-5-6 torus, K=0.161\n"
        "Piecewise-linear interpolation (Delaunay), 2×2 tiled",
        fontsize=12, y=1.01)

    out = (datadir / args.output
           if not Path(args.output).is_absolute() else Path(args.output))
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
