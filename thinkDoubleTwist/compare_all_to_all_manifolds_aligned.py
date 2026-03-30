#!/usr/bin/env python3
"""Align and compare two all-to-all triangular-lattice two-point manifolds.

The two datasets can come from different twisted/untwisted boundary conditions.
Each is mapped into normalized cell coordinates (a,b) where:
  [m,n]^T = a * v + b * u,
  v = (Lx, Ty), u = (Tx, -Ly).

Workflow:
1) Load all-to-all rows: d m n corr err corr_conn err_conn
2) Convert each row to (a,b) and center around the peak (max channel).
3) Interpolate both manifolds on a common centered grid in (a,b).
4) Plot side-by-side surfaces and a difference surface.
5) Plot aligned boundary-direction line cuts for v, u, and w=v+u.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def load_all_to_all(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append(
                {
                    "d": int(d_s),
                    "m": int(m_s),
                    "n": int(n_s),
                    "corr": float(c_s),
                    "err": float(e_s),
                    "corr_conn": float(cc_s),
                    "err_conn": float(ec_s),
                }
            )
    if not rows:
        raise SystemExit(f"No data rows found in {path}")
    return rows


def to_ab(m: np.ndarray, n: np.ndarray, lx: int, ly: int, tx: int, ty: int):
    ncell = float(lx * ly + tx * ty)
    a = (ly * m + tx * n) / ncell
    b = (ty * m - lx * n) / ncell
    return a, b


def wrap_centered(x: np.ndarray):
    return x - np.floor(x + 0.5)


def prepare_points(rows, lx: int, ly: int, tx: int, ty: int, channel: str):
    m = np.array([r["m"] for r in rows], dtype=float)
    n = np.array([r["n"] for r in rows], dtype=float)
    z = np.array([r[channel] for r in rows], dtype=float)

    a, b = to_ab(m, n, lx, ly, tx, ty)

    # Align manifolds by centering at the maximum of the chosen channel.
    i0 = int(np.argmax(z))
    a0 = a[i0]
    b0 = b[i0]

    ac = wrap_centered(a - a0)
    bc = wrap_centered(b - b0)
    return ac, bc, z


def tile_3x3(a: np.ndarray, b: np.ndarray, z: np.ndarray):
    aa, bb, zz = [], [], []
    for da in (-1.0, 0.0, 1.0):
        for db in (-1.0, 0.0, 1.0):
            aa.append(a + da)
            bb.append(b + db)
            zz.append(z)
    return np.concatenate(aa), np.concatenate(bb), np.concatenate(zz)


def make_interpolator(a: np.ndarray, b: np.ndarray, z: np.ndarray):
    tri = mtri.Triangulation(a, b)
    return mtri.LinearTriInterpolator(tri, z)


def eval_interp(interp, a_grid: np.ndarray, b_grid: np.ndarray):
    val = interp(a_grid, b_grid)
    if hasattr(val, "filled"):
        val = val.filled(np.nan)
    return np.asarray(val, dtype=float)


def sample_segment(interp, a1: float, b1: float, n: int):
    t = np.linspace(0.0, 1.0, n)
    a = t * a1
    b = t * b1
    z = eval_interp(interp, a, b)
    return t, z


def plot_surfaces(a_grid, b_grid, z_u, z_t, z_d, out_path: Path, channel_label: str):
    fig = plt.figure(figsize=(14, 4.8), dpi=200)

    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax1.plot_surface(a_grid, b_grid, z_u, cmap="viridis", linewidth=0, antialiased=True, alpha=0.95)
    ax1.set_title("Untwisted Aligned")
    ax1.set_xlabel("a (v direction)")
    ax1.set_ylabel("b (u direction)")
    ax1.set_zlabel(channel_label)

    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.plot_surface(a_grid, b_grid, z_t, cmap="viridis", linewidth=0, antialiased=True, alpha=0.95)
    ax2.set_title("Twisted Aligned")
    ax2.set_xlabel("a (v direction)")
    ax2.set_ylabel("b (u direction)")
    ax2.set_zlabel(channel_label)

    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    ax3.plot_surface(a_grid, b_grid, z_d, cmap="coolwarm", linewidth=0, antialiased=True, alpha=0.95)
    ax3.set_title("Difference (Twisted - Untwisted)")
    ax3.set_xlabel("a (v direction)")
    ax3.set_ylabel("b (u direction)")
    ax3.set_zlabel("Delta")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_line_cuts(interp_u, interp_t, out_path: Path, channel_label: str):
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), dpi=200, sharex=True)

    dirs = {
        "v (a-axis)": (0.5, 0.0),
        "u (b-axis)": (0.0, 0.5),
        "w (a=b)": (0.5, 0.5),
    }
    colors = {
        "v (a-axis)": "red",
        "u (b-axis)": "blue",
        "w (a=b)": "black",
    }

    for name, (a1, b1) in dirs.items():
        t, zu = sample_segment(interp_u, a1, b1, 300)
        _t, zt = sample_segment(interp_t, a1, b1, 300)
        mask = np.isfinite(zu) & np.isfinite(zt)

        axes[0].plot(t[mask], zu[mask], color=colors[name], lw=1.5, ls="-", label=f"untwisted {name}")
        axes[0].plot(t[mask], zt[mask], color=colors[name], lw=1.5, ls="--", label=f"twisted {name}")

        axes[1].plot(t[mask], zt[mask] - zu[mask], color=colors[name], lw=1.7, label=name)

    axes[0].set_title(f"Aligned Boundary Direction Cuts ({channel_label})")
    axes[0].set_ylabel(channel_label)
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8, ncol=2)

    axes[1].set_title("Difference Along Aligned Cuts (Twisted - Untwisted)")
    axes[1].set_xlabel("fraction from center to boundary")
    axes[1].set_ylabel("Delta")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--untwisted", required=True, help="Untwisted all-to-all .dat")
    ap.add_argument("--twisted", required=True, help="Twisted all-to-all .dat")

    ap.add_argument("--u_Lx", type=int, required=True)
    ap.add_argument("--u_Ly", type=int, required=True)
    ap.add_argument("--u_Tx", type=int, required=True)
    ap.add_argument("--u_Ty", type=int, required=True)

    ap.add_argument("--t_Lx", type=int, required=True)
    ap.add_argument("--t_Ly", type=int, required=True)
    ap.add_argument("--t_Tx", type=int, required=True)
    ap.add_argument("--t_Ty", type=int, required=True)

    ap.add_argument("--channel", choices=["corr", "corr_conn"], default="corr_conn")
    ap.add_argument("--grid", type=int, default=180, help="Grid points per axis for aligned surfaces")
    ap.add_argument("--surface_out", required=True, help="Output png for 3-surface figure")
    ap.add_argument("--cuts_out", required=True, help="Output png for aligned direction cuts")
    args = ap.parse_args()

    untwisted_rows = load_all_to_all(Path(args.untwisted))
    twisted_rows = load_all_to_all(Path(args.twisted))

    au, bu, zu = prepare_points(untwisted_rows, args.u_Lx, args.u_Ly, args.u_Tx, args.u_Ty, args.channel)
    at, bt, zt = prepare_points(twisted_rows, args.t_Lx, args.t_Ly, args.t_Tx, args.t_Ty, args.channel)

    au3, bu3, zu3 = tile_3x3(au, bu, zu)
    at3, bt3, zt3 = tile_3x3(at, bt, zt)

    interp_u = make_interpolator(au3, bu3, zu3)
    interp_t = make_interpolator(at3, bt3, zt3)

    grid = args.grid
    avals = np.linspace(-0.5, 0.5, grid)
    bvals = np.linspace(-0.5, 0.5, grid)
    a_grid, b_grid = np.meshgrid(avals, bvals, indexing="xy")

    z_u = eval_interp(interp_u, a_grid, b_grid)
    z_t = eval_interp(interp_t, a_grid, b_grid)
    z_d = z_t - z_u

    channel_label = "C_conn" if args.channel == "corr_conn" else "C"
    plot_surfaces(a_grid, b_grid, z_u, z_t, z_d, Path(args.surface_out), channel_label)
    plot_line_cuts(interp_u, interp_t, Path(args.cuts_out), channel_label)

    finite = np.isfinite(z_d)
    if np.any(finite):
        rms = float(np.sqrt(np.nanmean(z_d[finite] ** 2)))
        mad = float(np.nanmean(np.abs(z_d[finite])))
        print(f"surface_out={args.surface_out}")
        print(f"cuts_out={args.cuts_out}")
        print(f"delta_rms={rms:.8e}")
        print(f"delta_mean_abs={mad:.8e}")
    else:
        print("No finite overlap region found in aligned grid.")


if __name__ == "__main__":
    main()
