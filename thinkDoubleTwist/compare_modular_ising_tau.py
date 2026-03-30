#!/usr/bin/env python3
"""Compare MC all-to-all correlator to modular-torus Ising CFT expectation.

Given twisted-parallelogram BC vectors
  v = (L_x, T_y), u = (T_x, -L_y)
in triangular-lattice coordinates, define complex embedding with
  e1 = 1, e2 = 1/2 + i*sqrt(3)/2.
Then tau is chosen from geometry as tau = u / v (with Im(tau) > 0 convention).

For each lattice displacement, map to torus coordinate nu = a + b*tau,
where (a,b) solves p = a v + b u. Compare the connected MC correlator against
critical Ising CFT torus spin-spin shape:

G(nu|tau) ~ 0.5 * |theta1'(0|tau)/theta1(pi*nu|tau)|^(1/4)
          * [ |theta2(pi*nu/2)/theta2(0)|^(1/2)
            + |theta3(pi*nu/2)/theta3(0)|^(1/2)
            + |theta4(pi*nu/2)/theta4(0)|^(1/2) ]

The overall normalization is fit by weighted least squares to MC data.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from mpmath import mp


def lattice_to_xy(m: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = m + 0.5 * n
    y = (np.sqrt(3.0) / 2.0) * n
    return x, y


def rotate_xy(x: np.ndarray, y: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    c = np.cos(theta)
    s = np.sin(theta)
    return c * x - s * y, s * x + c * y


def load_full(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            d_s, m_s, n_s, c_s, e_s, cc_s, ec_s = s.split()
            rows.append((int(m_s), int(n_s), float(cc_s), float(ec_s)))
    if not rows:
        raise SystemExit("No all-to-all rows found.")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def to_uv(m: np.ndarray, n: np.ndarray, lx: int, ly: int, tx: int, ty: int) -> tuple[np.ndarray, np.ndarray]:
    ncell = float(lx * ly + tx * ty)
    a = (ly * m + tx * n) / ncell
    b = (ty * m - lx * n) / ncell
    return a, b


def modular_tau(lx: int, ly: int, tx: int, ty: int) -> tuple[complex, int]:
    omega = 0.5 + 0.5j * np.sqrt(3.0)
    v = complex(lx + ty * omega)
    u = complex(tx - ly * omega)
    tau = u / v
    b_sign = 1
    if tau.imag < 0:
        tau = -tau
        b_sign = -1
    return tau, b_sign


def ising_torus_shape(nu: complex, tau: complex, theta1p0: mp.mpf, q: mp.mpc) -> float:
    # Eq. (57) of arXiv:2209.15546:
    # <sigma(0)sigma(z)> = |theta1'(0|tau)/theta1(z|tau)|^(1/4)
    #                      * (sum_{nu=1}^4 |theta_nu(z/2|tau)|)
    #                        /(sum_{nu=2}^4 |theta_nu(0|tau)|)
    # mpmath jtheta uses argument with periods (pi, pi*tau), so we pass z->pi*nu.
    z = mp.pi * mp.mpc(nu.real, nu.imag)
    z2 = 0.5 * z

    th1 = mp.jtheta(1, z, q)
    if abs(th1) == 0:
        return float("nan")

    pref = abs(theta1p0 / th1) ** mp.mpf("0.25")
    num = (
        abs(mp.jtheta(1, z2, q))
        + abs(mp.jtheta(2, z2, q))
        + abs(mp.jtheta(3, z2, q))
        + abs(mp.jtheta(4, z2, q))
    )
    den = (
        abs(mp.jtheta(2, mp.mpf("0.0"), q))
        + abs(mp.jtheta(3, mp.mpf("0.0"), q))
        + abs(mp.jtheta(4, mp.mpf("0.0"), q))
    )
    if den == 0:
        return float("nan")

    g = pref * (num / den)
    return float(g)


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare MC correlator to modular Ising-torus expectation.")
    ap.add_argument("--full", required=True, help="Path to full all-to-all .dat")
    ap.add_argument("--L_x", type=int, required=True)
    ap.add_argument("--L_y", type=int, required=True)
    ap.add_argument("--T_x", type=int, required=True)
    ap.add_argument("--T_y", type=int, required=True)
    ap.add_argument("--output", required=True, help="Output PNG path")
    ap.add_argument("--report", required=True, help="Output text report path")
    ap.add_argument("--mp-dps", type=int, default=50, help="mpmath precision")
    ap.add_argument(
        "--fit-max-fraction",
        type=float,
        default=1.0,
        help=(
            "Use only points with torus-distance from singularity <= this fraction "
            "(in wrapped [0, 0.5] units) for the alpha fit."
        ),
    )
    ap.add_argument(
        "--alpha-override",
        type=float,
        default=None,
        help="If provided, use this alpha normalization instead of weighted least-squares fit.",
    )
    ap.add_argument(
        "--enforce-short-distance-peak",
        action="store_true",
        help=(
            "If set, renormalize theory so max(model) at displacements within "
            "one lattice spacing from origin is above --peak-threshold."
        ),
    )
    ap.add_argument(
        "--peak-threshold",
        type=float,
        default=1.05,
        help="Target lower bound for short-distance theory peak when enforcing peak normalization.",
    )
    ap.add_argument("--align-v-x", action="store_true", help="Rotate x/y so v aligns to +x")
    ap.add_argument("--tile-a-min", type=int, default=0, help="Minimum tile index along v")
    ap.add_argument("--tile-a-max", type=int, default=1, help="Maximum tile index along v")
    ap.add_argument("--tile-b-min", type=int, default=0, help="Minimum tile index along u")
    ap.add_argument("--tile-b-max", type=int, default=1, help="Maximum tile index along u")
    args = ap.parse_args()

    mp.dps = args.mp_dps

    m, n, c_mc, ec_mc = load_full(Path(args.full))

    tau, b_sign = modular_tau(args.L_x, args.L_y, args.T_x, args.T_y)

    a, b = to_uv(m, n, args.L_x, args.L_y, args.T_x, args.T_y)
    b = b_sign * b

    # Wrap into principal torus domain.
    a = np.mod(a, 1.0)
    b = np.mod(b, 1.0)

    q = mp.e ** (mp.pi * 1j * mp.mpc(tau.real, tau.imag))
    theta1p0 = mp.diff(lambda zz: mp.jtheta(1, zz, q), mp.mpf("0.0"))

    g_th = np.zeros_like(c_mc)
    mask_good = np.ones_like(c_mc, dtype=bool)
    for i in range(len(c_mc)):
        nu = complex(float(a[i] + b[i] * tau.real), float(b[i] * tau.imag))
        # Skip exact origin singularity.
        if abs(nu.real) < 1.0e-14 and abs(nu.imag) < 1.0e-14:
            g_th[i] = np.nan
            mask_good[i] = False
            continue
        val = ising_torus_shape(nu, tau, theta1p0, q)
        if not np.isfinite(val):
            g_th[i] = np.nan
            mask_good[i] = False
        else:
            g_th[i] = val

    # Weighted normalization fit: c_mc ~= alpha * g_th.
    w = np.zeros_like(c_mc)
    pos_err = ec_mc > 0
    w[pos_err] = 1.0 / (ec_mc[pos_err] ** 2)

    # Restrict fit region to a neighborhood of the singularity in torus coords.
    # a,b are in [0,1); distance uses wrapped coords in [0,0.5].
    a_wrap = np.minimum(a, 1.0 - a)
    b_wrap = np.minimum(b, 1.0 - b)
    fit_radius = np.sqrt(a_wrap * a_wrap + b_wrap * b_wrap)
    fit_region_mask = fit_radius <= float(args.fit_max_fraction)

    fit_mask = mask_good & np.isfinite(c_mc) & np.isfinite(g_th) & (w > 0) & fit_region_mask
    if not np.any(fit_mask):
        raise SystemExit("No valid points for weighted fit.")

    num = float(np.sum(w[fit_mask] * c_mc[fit_mask] * g_th[fit_mask]))
    den = float(np.sum(w[fit_mask] * g_th[fit_mask] * g_th[fit_mask]))
    alpha_fit = num / den if den != 0.0 else 0.0
    alpha = float(args.alpha_override) if args.alpha_override is not None else alpha_fit

    c_th = alpha * g_th

    # Optional post-fit normalization: make sure the short-distance theory
    # peak (within one lattice spacing) is above a requested threshold.
    peak_scale = 1.0
    peak_before = float("nan")
    peak_after = float("nan")
    if args.enforce_short_distance_peak:
        x0, y0 = lattice_to_xy(m, n)
        r0 = np.sqrt(x0 * x0 + y0 * y0)
        short_mask = (r0 > 1.0e-12) & (r0 <= 1.0 + 1.0e-12) & np.isfinite(c_th)
        if np.any(short_mask):
            peak_before = float(np.max(c_th[short_mask]))
            if peak_before <= args.peak_threshold:
                if peak_before > 0.0:
                    peak_scale = args.peak_threshold / peak_before
                else:
                    peak_scale = 1.0
                alpha *= peak_scale
                c_th *= peak_scale
            peak_after = float(np.max(c_th[short_mask]))
    residual = c_mc - c_th
    pull = np.full_like(c_mc, np.nan)
    good_pull = ec_mc > 0
    pull[good_pull] = residual[good_pull] / ec_mc[good_pull]

    chi2 = float(np.sum(w[fit_mask] * residual[fit_mask] ** 2))
    dof = int(np.sum(fit_mask) - 1)
    chi2_dof = chi2 / dof if dof > 0 else float("nan")
    rmse = float(np.sqrt(np.mean(residual[fit_mask] ** 2)))

    # Coordinates for plotting on tiled copies of the fundamental domain.
    v_m, v_n = args.L_x, args.T_y
    u_m, u_n = args.T_x, -args.L_y

    m_tiles = []
    n_tiles = []
    c_mc_tiles = []
    c_th_tiles = []
    residual_tiles = []
    ec_tiles = []
    pull_tiles = []
    for ia in range(args.tile_a_min, args.tile_a_max + 1):
        for ib in range(args.tile_b_min, args.tile_b_max + 1):
            m_tiles.append(m + ia * v_m + ib * u_m)
            n_tiles.append(n + ia * v_n + ib * u_n)
            c_mc_tiles.append(c_mc)
            c_th_tiles.append(c_th)
            residual_tiles.append(residual)
            ec_tiles.append(ec_mc)
            pull_tiles.append(pull)

    m_plot = np.concatenate(m_tiles)
    n_plot = np.concatenate(n_tiles)
    c_mc_plot = np.concatenate(c_mc_tiles)
    c_th_plot = np.concatenate(c_th_tiles)
    residual_plot = np.concatenate(residual_tiles)
    ec_plot = np.concatenate(ec_tiles)
    pull_plot = np.concatenate(pull_tiles)

    x, y = lattice_to_xy(m_plot, n_plot)
    if args.align_v_x:
        vx, vy = lattice_to_xy(np.array([args.L_x], dtype=float), np.array([args.T_y], dtype=float))
        theta = -math.atan2(float(vy[0]), float(vx[0]))
        x, y = rotate_xy(x, y, theta)

    tri_data = mtri.Triangulation(x, y)
    tri_theory = mtri.Triangulation(x, y)
    bad = ~(np.isfinite(c_th_plot) & np.isfinite(residual_plot))
    if np.any(bad):
        t = tri_theory.triangles
        tri_mask = bad[t[:, 0]] | bad[t[:, 1]] | bad[t[:, 2]]
        tri_theory.set_mask(tri_mask)

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.6), dpi=220, constrained_layout=True)
    axes = axes.ravel()

    finite_both = np.isfinite(c_mc_plot) & np.isfinite(c_th_plot)
    if np.any(finite_both):
        vmin_ct = float(min(np.nanmin(c_mc_plot[finite_both]), np.nanmin(c_th_plot[finite_both])))
        vmax_ct = float(max(np.nanmax(c_mc_plot[finite_both]), np.nanmax(c_th_plot[finite_both])))
    else:
        vmin_ct = float(np.nanmin(c_mc_plot))
        vmax_ct = float(np.nanmax(c_mc_plot))

    p0 = axes[0].tricontourf(tri_data, c_mc_plot, levels=24, cmap="viridis", vmin=vmin_ct, vmax=vmax_ct)
    fig.colorbar(p0, ax=axes[0], shrink=0.86, pad=0.02)
    axes[0].set_title("MC connected correlator")

    p1 = axes[1].tricontourf(tri_theory, c_th_plot, levels=24, cmap="viridis", vmin=vmin_ct, vmax=vmax_ct)
    fig.colorbar(p1, ax=axes[1], shrink=0.86, pad=0.02)
    axes[1].set_title("Modular Ising prediction (fit scale)")

    vmax = np.nanmax(np.abs(residual_plot[np.isfinite(residual_plot)]))
    p2 = axes[2].tricontourf(tri_theory, residual_plot, levels=24, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    fig.colorbar(p2, ax=axes[2], shrink=0.86, pad=0.02)
    axes[2].set_title("Residual: MC - theory")

    # Uncertainty-aware panel: pull = residual / sigma_MC.
    finite_pull = np.isfinite(pull_plot)
    if np.any(finite_pull):
        pmax = np.nanmax(np.abs(pull_plot[finite_pull]))
        pmax = max(1.0, min(pmax, 8.0))
        p3 = axes[3].tricontourf(tri_theory, pull_plot, levels=24, cmap="RdBu_r", vmin=-pmax, vmax=pmax)
        fig.colorbar(p3, ax=axes[3], shrink=0.86, pad=0.02)
    else:
        p3 = axes[3].tricontourf(tri_theory, np.zeros_like(c_th_plot), levels=24, cmap="RdBu_r", vmin=-1.0, vmax=1.0)
        fig.colorbar(p3, ax=axes[3], shrink=0.86, pad=0.02)
    axes[3].set_title("z-score tension: (MC - theory) / sigma_MC")

    for ax in axes:
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle(
        f"Modular-space Ising comparison (tiles a={args.tile_a_min}..{args.tile_a_max}, b={args.tile_b_min}..{args.tile_b_max}): "
        f"tau={tau.real:.6f}+{tau.imag:.6f}i, alpha={alpha:.6e}, chi2/dof={chi2_dof:.3f}",
        fontsize=11,
    )

    out_png = Path(args.output)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)

    report = Path(args.report)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "Modular Ising vs MC correlator comparison",
                f"full_data = {args.full}",
                f"geometry = (Lx={args.L_x}, Ly={args.L_y}, Tx={args.T_x}, Ty={args.T_y})",
                f"tau = {tau.real:.15f} + {tau.imag:.15f} i",
                f"fit_points = {int(np.sum(fit_mask))}",
                f"fit_max_fraction = {float(args.fit_max_fraction):.16e}",
                f"alpha = {alpha:.16e}",
                f"alpha_fit = {alpha_fit:.16e}",
                f"alpha_override = {float(args.alpha_override):.16e}" if args.alpha_override is not None else "alpha_override = none",
                f"peak_norm_enabled = {int(args.enforce_short_distance_peak)}",
                f"peak_threshold = {float(args.peak_threshold):.16e}",
                f"peak_scale = {peak_scale:.16e}",
                f"short_distance_peak_before = {peak_before:.16e}",
                f"short_distance_peak_after = {peak_after:.16e}",
                f"chi2 = {chi2:.16e}",
                f"dof = {dof}",
                f"chi2_per_dof = {chi2_dof:.16e}",
                f"rmse = {rmse:.16e}",
                f"mc_sigma_mean = {float(np.nanmean(ec_mc)):.16e}",
                f"mc_sigma_median = {float(np.nanmedian(ec_mc)):.16e}",
                f"z_score_tension_rms = {float(np.sqrt(np.nanmean(pull[fit_mask] ** 2))):.16e}",
                "",
            ]
        )
    )

    print(f"Wrote {out_png}")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
