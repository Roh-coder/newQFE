#!/usr/bin/env python3
"""Phase 4: Compare 3D MC spatial correlator to 3D Ising CFT prediction.

The 3D Ising CFT two-point function of the spin operator sigma on an infinite
flat space takes the power-law form:

  G_flat(r) = C / r^(2 * Delta_sigma)

with Delta_sigma ≈ 0.5182 (Kos, Poland, Simmons-Duffin, Vichi 2016).

On a finite torus with geometry set by (Lx, Ly, Tx, Ty, Nt) we include the
leading periodic image sum (method-of-images) as a proxy for finite-volume
corrections. The full sum over the lattice of torus images of the power-law
kernel approximates the CFT two-point function on the flat torus T^2 × S^1:

  G_torus(x) ~ sum_{n in Z^3} C / |x - n_1 v1 - n_2 v2 - n_3 t_hat|^(2*Delta)

truncated at N_images shells. This is not an exact CFT result but captures the
leading finite-L corrections and serves as a reference for deviation from
free-field / mean-field behaviour.

The isotropic 3D flat-space correlator is also drawn as a straight line on a
log-log plot as the "bulk" reference.

For a proper 3D Ising CFT calculation on the torus one would need the full
conformal block decomposition; this script provides the simpler power-law
plus image-sum benchmark for comparison.

Usage:
  python compare_456_3d_twopoint.py --datadir . --output compare_3d_cft.pdf
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DELTA_SIGMA = 0.5182  # 3D Ising critical exponent (bootstrap, 2016)

# ── Lattice geometry helpers ──────────────────────────────────────────────────

def triangular_basis():
    """Real-space unit vectors of the triangular lattice (e1, e2)."""
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.5, np.sqrt(3.0) / 2.0])
    return e1, e2


def lattice_to_xy(m, n):
    e1, e2 = triangular_basis()
    return m[:, None] * e1 + n[:, None] * e2  # shape (N, 2)


def torus_period_vectors(lx, ly, tx, ty, nt, kt=1.0):
    """Compute the three real-space period vectors of the 3D torus.

    v1 = Lx * e1 + Ty * e2  (spatial)
    v2 = Tx * e1 - Ly * e2  (spatial)  [sign convention from twisted parallelogram]
    v3 = (0, 0, Nt * kt_spacing)  (temporal, using kt=1 lattice spacing)

    Returns v1, v2 as 2-vectors and nt_len as the temporal period.
    """
    e1, e2 = triangular_basis()
    v1 = lx * e1 + ty * e2
    v2 = tx * e1 - ly * e2  # note: stored as (Tx, -Ly) in the code
    return v1, v2, float(nt)  # temporal period in units of kt spacing


# ── 3D Ising CFT power-law in infinite space ─────────────────────────────────

def g_flat(r, delta=DELTA_SIGMA):
    """G(r) = 1 / r^(2*delta) (normalized C=1)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(r > 0, r ** (-2.0 * delta), np.nan)


# ── Method-of-images torus sum ────────────────────────────────────────────────

def g_torus_image_sum(xy, v1, v2, nt_period, n_images=3, delta=DELTA_SIGMA):
    """Evaluate power-law two-point function on the 3D torus via image sum.

    xy:    (N, 2) physical positions of lattice sites
    v1,v2: 2-vectors: spatial periodicity vectors
    nt_period: float temporal period (in temporal lattice units)
    n_images: number of image shells in each direction (total (2N+1)^3 images)

    Returns array of length N with G_torus(xy) (origin excluded).
    """
    N = len(xy)
    g = np.zeros(N)
    for n1 in range(-n_images, n_images + 1):
        for n2 in range(-n_images, n_images + 1):
            for n3 in range(-n_images, n_images + 1):
                shift_xy = n1 * v1 + n2 * v2  # 2D spatial shift
                shift_t = n3 * nt_period       # temporal shift
                dx = xy[:, 0] - shift_xy[0]
                dy = xy[:, 1] - shift_xy[1]
                dt = shift_t
                r2 = dx * dx + dy * dy + dt * dt
                # Skip the (0,0,0) image at the origin itself
                if n1 == 0 and n2 == 0 and n3 == 0:
                    mask_nonzero = r2 > 0
                    g[mask_nonzero] += r2[mask_nonzero] ** (-delta)
                else:
                    g += r2 ** (-delta)
    return g  # unnormalized


# ── Data loading ──────────────────────────────────────────────────────────────

def load_connected(path: Path):
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            rows.append((int(parts[1]), int(parts[2]), float(parts[5]), float(parts[6])))
    if not rows:
        raise SystemExit(f"No data in {path}")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def parse_geometry(dirname: str):
    m = re.search(r"Lx(\d+)_Ly(\d+)_Tx(-?\d+)_Ty(-?\d+)_Nt(\d+)", dirname)
    if not m:
        return None
    return {k: int(v) for k, v in zip(["Lx", "Ly", "Tx", "Ty", "Nt"], m.groups())}


def radial_bin_weighted(r, vals, errs, nbins=50):
    rmax = np.nanmax(r[np.isfinite(r)])
    edges = np.linspace(0, rmax, nbins + 1)
    cx, cy, cye = [], [], []
    for i in range(nbins):
        mask = (r >= edges[i]) & (r < edges[i + 1]) & np.isfinite(vals) & (errs > 0)
        if mask.sum() == 0:
            continue
        w = 1.0 / errs[mask] ** 2
        cx.append(0.5 * (edges[i] + edges[i + 1]))
        cy.append(np.sum(w * vals[mask]) / np.sum(w))
        cye.append(1.0 / np.sqrt(np.sum(w)))
    return np.array(cx), np.array(cy), np.array(cye)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Compare 3D MC spatial correlator to 3D Ising CFT power-law prediction"
    )
    ap.add_argument("--datadir", default=".", help="Directory containing run subdirectories")
    ap.add_argument("--output", default="compare_3d_cft.png", help="Output figure path")
    ap.add_argument("--delta", type=float, default=DELTA_SIGMA,
                    help=f"3D Ising spin scaling dimension (default: {DELTA_SIGMA})")
    ap.add_argument("--n-images", type=int, default=2,
                    help="Number of image shells in each direction for torus sum (default: 2)")
    args = ap.parse_args()

    datadir = Path(args.datadir)
    pattern = str(datadir / "Lx*_Ly*_Tx*_Ty*_Nt*_k*/two_point_all_to_all.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No data files found: {pattern}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_log = axes[0]   # log-log: G(r) vs r
    ax_scaled = axes[1]  # scaled: G(r)*r^(2*Delta) vs r/L (collapse check)

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(files)))

    for fpath, color in zip(files, colors):
        fpath = Path(fpath)
        geo = parse_geometry(fpath.parent.name)
        if geo is None:
            print(f"Skipping: {fpath.parent.name}")
            continue
        lx, ly, tx, ty, nt = geo["Lx"], geo["Ly"], geo["Tx"], geo["Ty"], geo["Nt"]
        s = lx // 39
        label = f"s={s}  (Lx={lx}, Nt={nt})"

        print(f"Processing {fpath.parent.name} ...", flush=True)
        m_arr, n_arr, g_mc, g_err = load_connected(fpath)

        # Physical positions (equal-time, spatial only)
        xy = lattice_to_xy(m_arr.astype(int), n_arr.astype(int))
        r = np.sqrt(m_arr ** 2 + n_arr ** 2 + m_arr * n_arr)

        # 3D torus period vectors
        v1, v2, nt_period = torus_period_vectors(lx, ly, tx, ty, nt)

        # Method-of-images torus sum (equal-time: set temporal component of image to 0)
        # For equal-time spatial correlator, we sum only spatial images (n3=0 slice),
        # giving an effective 2D torus sum.
        g_img = np.zeros(len(m_arr))
        for n1 in range(-args.n_images, args.n_images + 1):
            for n2 in range(-args.n_images, args.n_images + 1):
                shift = n1 * v1 + n2 * v2
                dx = xy[:, 0] - shift[0]
                dy = xy[:, 1] - shift[1]
                r2 = dx * dx + dy * dy
                if n1 == 0 and n2 == 0:
                    mask = r2 > 0
                    g_img[mask] += r2[mask] ** (-args.delta)
                else:
                    g_img += r2 ** (-args.delta)

        # Fit normalization to MC data
        valid = (g_img > 0) & np.isfinite(g_mc) & (g_err > 0)
        if valid.sum() < 2:
            print("  WARNING: too few valid points for fit.")
            continue
        w = 1.0 / g_err[valid] ** 2
        alpha = np.sum(w * g_mc[valid] * g_img[valid]) / np.sum(w * g_img[valid] ** 2)
        g_theory = alpha * g_img

        # Bin
        rc, gc_mc_bin, gc_mc_e = radial_bin_weighted(r, g_mc, g_err)
        rc2, gc_th_bin, _ = radial_bin_weighted(r, g_theory, np.ones_like(r) * 1e-10)

        # Plot log-log
        ax_log.errorbar(rc, gc_mc_bin, yerr=gc_mc_e, fmt="o", markersize=3,
                        color=color, label=f"MC {label}")
        ax_log.plot(rc2, gc_th_bin, "--", color=color, lw=1.2,
                    label=f"CFT image sum {label}")

        # Scaling collapse: G(r) * r^(2*Delta) vs r/L
        ax_scaled.errorbar(rc / lx, gc_mc_bin * rc ** (2 * args.delta),
                           yerr=gc_mc_e * rc ** (2 * args.delta),
                           fmt="o", markersize=3, color=color, label=f"MC {label}")
        ax_scaled.plot(rc2 / lx, gc_th_bin * rc2 ** (2 * args.delta),
                       "--", color=color, lw=1.2)

    # Flat-space bulk power law reference on log-log
    r_ref = np.logspace(np.log10(0.5), np.log10(30), 100)
    ax_log.plot(r_ref, g_flat(r_ref, args.delta) * 0.05, "k:", lw=1,
                label=rf"$r^{{-2\Delta}}$ bulk, $\Delta$={args.delta}")

    ax_log.set_xscale("log")
    ax_log.set_yscale("log")
    ax_log.set_xlabel(r"$r$ (lattice units)")
    ax_log.set_ylabel(r"$G_{\rm spatial}(r)$")
    ax_log.set_title("Spatial correlator vs 3D Ising CFT image sum")
    ax_log.legend(fontsize=7)

    ax_scaled.set_xlabel(r"$r / L_x$")
    ax_scaled.set_ylabel(rf"$G(r) \cdot r^{{2\Delta_\sigma}}$,  $\Delta_\sigma={args.delta}$")
    ax_scaled.set_title("Scaling collapse check: MC vs 3D CFT image sum")
    ax_scaled.legend(fontsize=7)

    fig.suptitle("3D Ising on 4-5-6 torus at K=0.161 vs 3D CFT power-law + images", fontsize=11)
    fig.tight_layout()
    outpath = datadir / args.output if not Path(args.output).is_absolute() else Path(args.output)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
