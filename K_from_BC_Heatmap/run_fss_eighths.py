#!/usr/bin/env python3
"""
run_fss_eighths.py — minimal FSS pipeline for a single (r1, r2) point.

For one chosen (r1, r2):
  1. Run MC at every L in --sizes for both the "test" (untwisted) and
     "ref" (twisted) geometries, caching pickles under
     results/<tag>/grid/L<L>/{test,ref}/.
  2. Sample the connected two-point function G_conn at fractional cycle
     positions  t_k = k/8,  k = 1..7  along each of the three torus
     boundary directions.
  3. Fit G(t_k; L) → G_∞(t_k) with three continuum extrapolations
     (linear in 1/L, quadratic in 1/L, exponential a + b·exp(−L/ξ))
     and overlay all three on a single per-cycle panel.

Usage (defaults: r1=r2=1.0, sizes=8 12 16 20, n_traj=20000):
    python3 run_fss_eighths.py --tag fss_eighths --n-workers 4
"""
from __future__ import annotations

import argparse
import math
import os
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine                                         # noqa: E402
from cost import boundary_paths, _tile_interp, _SQRT3_2  # noqa: E402
from precompute_grid import _point_path, _run_one        # noqa: E402


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def _make_geom(L, lx_mult, ly_mult, tx_frac, ty_frac):
    return (int(round(lx_mult * L)), int(round(ly_mult * L)),
            int(round(tx_frac * L)), int(round(ty_frac * L)))


# ---------------------------------------------------------------------------
# Stage 1 — MC precompute (resumable; one pickle per (L, kind))
# ---------------------------------------------------------------------------

def stage_mc(args):
    exe = args.exe or os.path.join(_HERE, "bin",
                                   "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        print(f"ERROR: simulator not found: {exe}", file=sys.stderr)
        sys.exit(1)
    out_root = os.path.join(_HERE, "results", args.tag)
    scratch_root = os.path.join(out_root, "_mc_scratch")
    os.makedirs(scratch_root, exist_ok=True)

    geoms = {}
    jobs = []
    for L in args.sizes:
        gt = _make_geom(L, 1.0, 1.0,
                        args.test_Tx_frac, args.test_Ty_frac)
        gr = _make_geom(L, 1.0, 1.0,
                        args.ref_Tx_frac, args.ref_Ty_frac)
        geoms[L] = {"test": gt, "ref": gr}
        for kind, geom in (("test", gt), ("ref", gr)):
            Lx, Ly, Tx, Ty = geom
            grid_dir = os.path.join(out_root, "grid", f"L{L}", kind)
            os.makedirs(grid_dir, exist_ok=True)
            out_pkl = _point_path(grid_dir, args.r1, args.r2)
            label = f"L{L}_{kind}"
            jobs.append((exe, label, Lx, Ly, Tx, Ty, args.r1, args.r2,
                          args.n_traj,
                          args.n_traj_coarse, args.n_traj_fine,
                          args.beta_seed[0], args.beta_seed[1],
                          scratch_root, out_pkl))

    n_total = len(jobs)
    n_cached = sum(1 for j in jobs if os.path.exists(j[-1]))
    print(f"[mc] {n_total} jobs ({n_cached} cached, "
          f"{n_total - n_cached} to run)  workers={args.n_workers}")
    for L in args.sizes:
        print(f"[mc]   L={L:3d}  test={geoms[L]['test']}  "
              f"ref={geoms[L]['ref']}")

    if n_total - n_cached > 0:
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futs = {ex.submit(_run_one, j): j for j in jobs}
            done = 0
            for fut in as_completed(futs):
                done += 1
                label, status, info = fut.result()
                if status == "ok":
                    print(f"[mc] [{done}/{n_total}] ok   {label}  "
                          f"beta_c={info['beta_c']:.5f}  "
                          f"wall={info['wall']:.1f}s", flush=True)
                elif status == "skip":
                    pass
                else:
                    print(f"[mc] [{done}/{n_total}] ERR  {label}  "
                          f"{info.get('msg','?')}", flush=True)
        print(f"[mc] done in {time.time()-t0:.1f}s")
    shutil.rmtree(scratch_root, ignore_errors=True)
    return geoms, out_root


# ---------------------------------------------------------------------------
# Stage 2 — sample correlator at t_k = k/8 along each cycle
# ---------------------------------------------------------------------------

# Eighth-spaced fractional positions (interior only).
T_FRACS = np.arange(1, 8) / 8.0   # 0.125 .. 0.875


def _sample_eighths(pkl_path):
    """Return arrays of shape (3, len(T_FRACS)): G and sigma_G per cycle."""
    with open(pkl_path, "rb") as f:
        pkl = pickle.load(f)
    data = pkl["data"]
    Lx, Ly, Tx, Ty = pkl["Lx"], pkl["Ly"], pkl["Tx"], pkl["Ty"]

    iG   = _tile_interp(data, Lx, Ly, Tx, Ty, "conn",     copies=2)
    iE   = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)

    paths = boundary_paths(Lx, Ly, Tx, Ty)
    G  = np.zeros((3, len(T_FRACS)))
    sG = np.zeros((3, len(T_FRACS)))
    for c, (dm, dn) in enumerate(paths):
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        pts = np.column_stack([T_FRACS * ex, T_FRACS * ey])
        G[c]  = np.asarray(iG(pts), dtype=float)
        sG[c] = np.abs(np.asarray(iE(pts), dtype=float))
    return G, sG, pkl


# ---------------------------------------------------------------------------
# Continuum fits
# ---------------------------------------------------------------------------

def _wls_poly(invL, y, sigma, deg):
    """Weighted least squares to a polynomial in 1/L of given degree.
    Returns (intercept, sigma_intercept, popt, cov, ok).
    """
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    x, yv, sv = invL[mask], y[mask], sigma[mask]
    if len(x) < deg + 1:
        return np.nan, np.nan, None, None, False
    X = np.vander(x, deg + 1, increasing=True)
    w = 1.0 / sv ** 2
    XtW = X.T * w
    try:
        cov = np.linalg.inv(XtW @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, None, None, False
    beta = cov @ (XtW @ yv)
    return float(beta[0]), float(np.sqrt(max(cov[0, 0], 0.0))), beta, cov, True


def _fit_power(L, y, sigma):
    """Fit  y = a + b * (1/L)**A   (NLS, ≥3 pts).

    Free exponent A bounded to (0.1, 6) to keep the model meaningful
    and the intercept identifiable.  The intercept is also clamped to
    a generous neighbourhood of the data range to suppress occasional
    runaways near the optimiser's bounds.
    """
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    Lm, ym, sm = L[mask], y[mask], sigma[mask]
    if len(Lm) < 3:
        return np.nan, np.nan, None, False
    invLm = 1.0 / Lm
    a0 = float(ym[np.argmax(Lm)])
    b0 = float(ym[np.argmin(Lm)] - a0)
    A0 = 1.5

    def model(LL, a, b, A):
        return a + b * (1.0 / LL) ** A

    y_min, y_max = float(ym.min()), float(ym.max())
    y_pad = max(y_max - y_min, 1e-6)
    bounds = (
        [y_min - 5 * y_pad, -100 * y_pad, 0.1],
        [y_max + 5 * y_pad,  100 * y_pad, 6.0],
    )
    try:
        popt, pcov = curve_fit(model, Lm, ym, p0=[a0, b0, A0],
                                sigma=sm, absolute_sigma=True,
                                bounds=bounds, maxfev=8000)
        a_fit = float(popt[0])
        sa_fit = float(np.sqrt(max(pcov[0, 0], 0.0)))
        if not np.isfinite(a_fit) or not np.isfinite(sa_fit):
            return np.nan, np.nan, None, False
        if abs(a_fit - 0.5 * (y_min + y_max)) > 10 * y_pad:
            return np.nan, np.nan, None, False
        return a_fit, sa_fit, popt, True
    except Exception:  # noqa: BLE001
        return np.nan, np.nan, None, False


# ---------------------------------------------------------------------------
# Stage 3 — plot
# ---------------------------------------------------------------------------

_FAMILIES = ("test", "ref")
_FAM_TITLE = {"test": "test (untwisted)", "ref": "ref (twisted)"}
_CYCLE_TITLE = ("cycle 0  (Lx, Ty)",
                "cycle 1  (Tx, −Ly)",
                "cycle 2  (−Lx−Tx, Ly−Ty)")


def stage_plot(args, geoms, out_root):
    sizes = sorted(args.sizes)
    Larr = np.asarray(sizes, dtype=float)
    invL = 1.0 / Larr

    # Gather samples: G[fam][cycle, k, L] and sigma analogously.
    G_data  = {fam: np.full((3, len(T_FRACS), len(sizes)), np.nan)
               for fam in _FAMILIES}
    sG_data = {fam: np.full((3, len(T_FRACS), len(sizes)), np.nan)
               for fam in _FAMILIES}
    for li, L in enumerate(sizes):
        for fam in _FAMILIES:
            pkl = _point_path(
                os.path.join(out_root, "grid", f"L{L}", fam),
                args.r1, args.r2)
            if not os.path.exists(pkl):
                print(f"[plot] WARNING missing pickle: {pkl}")
                continue
            G, sG, _ = _sample_eighths(pkl)
            G_data[fam][:, :, li]  = G
            sG_data[fam][:, :, li] = sG

    # Prepare a fine 1/L grid for plotting fit curves all the way to 0.
    invL_dense = np.linspace(0.0, invL.max() * 1.05, 200)
    L_dense    = np.where(invL_dense > 0, 1.0 / np.maximum(invL_dense, 1e-12),
                           np.inf)

    # Color per t_k fraction (viridis gradient).
    cmap = plt.get_cmap("viridis")
    n_k  = len(T_FRACS)
    frac_colors = [cmap(ki / max(n_k - 1, 1)) for ki in range(n_k)]

    # Marker shape per lattice size.
    _MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
    size_markers = {L: _MARKERS[li % len(_MARKERS)] for li, L in enumerate(sizes)}

    # Fixed colors for the three fit types.
    FIT_COLORS = {"linear": "#e05c2b", "quadratic": "#2b7fe0", "power": "#2bb55e"}
    FIT_LS     = {"linear": "-",       "quadratic": "--",       "power": ":"}
    FIT_MARK   = {"linear": "s",       "quadratic": "^",        "power": "D"}

    fig, axes = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    fig.suptitle(
        f"FSS continuum extrapolation @ (r1, r2) = "
        f"({args.r1:g}, {args.r2:g})  —  "
        f"sample positions t_k = k/8, k = 1..7  along each boundary cycle\n"
        f"sizes L = {sizes}   "
        f"test geometry frac=({args.test_Tx_frac},{args.test_Ty_frac})   "
        f"ref geometry frac=({args.ref_Tx_frac},{args.ref_Ty_frac})   "
        f"n_traj={args.n_traj}",
        fontsize=11)

    # Track summary of per-fit intercepts for printing.
    summary_rows = []

    for fi, fam in enumerate(_FAMILIES):
        for c in range(3):
            ax = axes[fi, c]
            for ki in range(n_k):
                fc = frac_colors[ki]
                y  = G_data[fam][c, ki, :]
                sy = sG_data[fam][c, ki, :]

                # Data points — shape by L, color by t_k.
                for li, L in enumerate(sizes):
                    if np.isfinite(y[li]):
                        ax.errorbar(invL[li], y[li], yerr=sy[li],
                                    fmt=size_markers[L], ms=5.5,
                                    color=fc, ecolor=fc, capsize=2,
                                    mec="k", mew=0.5, zorder=4)

                # ---- fits — color by fit type, trace follows t_k color ----
                a_lin, sa_lin, b_lin, _, ok_lin = _wls_poly(invL, y, sy, deg=1)
                a_qd,  sa_qd,  b_qd,  _, ok_qd  = _wls_poly(invL, y, sy, deg=2)
                a_pw,  sa_pw,  p_pw,     ok_pw  = _fit_power(Larr, y, sy)

                if ok_lin:
                    yfit = b_lin[0] + b_lin[1] * invL_dense
                    ax.plot(invL_dense, yfit,
                            ls=FIT_LS["linear"], color=FIT_COLORS["linear"],
                            lw=1.0, alpha=0.75)
                    ax.errorbar([0.0], [a_lin], yerr=[sa_lin],
                                fmt=FIT_MARK["linear"], ms=6,
                                color=FIT_COLORS["linear"],
                                mfc=fc, mec=FIT_COLORS["linear"], mew=1.2,
                                zorder=5)
                if ok_qd:
                    yfit = (b_qd[0] + b_qd[1] * invL_dense
                            + b_qd[2] * invL_dense ** 2)
                    ax.plot(invL_dense, yfit,
                            ls=FIT_LS["quadratic"], color=FIT_COLORS["quadratic"],
                            lw=1.0, alpha=0.75)
                    ax.errorbar([0.0], [a_qd], yerr=[sa_qd],
                                fmt=FIT_MARK["quadratic"], ms=6,
                                color=FIT_COLORS["quadratic"],
                                mfc=fc, mec=FIT_COLORS["quadratic"], mew=1.2,
                                zorder=5)
                if ok_pw:
                    a, b, A = p_pw
                    yfit = a + b * invL_dense ** A
                    ax.plot(invL_dense, yfit,
                            ls=FIT_LS["power"], color=FIT_COLORS["power"],
                            lw=1.2, alpha=0.85)
                    ax.errorbar([0.0], [a_pw], yerr=[sa_pw],
                                fmt=FIT_MARK["power"], ms=6,
                                color=FIT_COLORS["power"],
                                mfc=fc, mec=FIT_COLORS["power"], mew=1.2,
                                zorder=5)

                summary_rows.append((fam, c, ki, T_FRACS[ki],
                                     a_lin, sa_lin,
                                     a_qd,  sa_qd,
                                     a_pw,  sa_pw,
                                     (p_pw[2] if ok_pw else np.nan)))

            ax.axvline(0.0, color="k", lw=0.6, ls=":", alpha=0.5)
            ax.set_xlim(-0.005, invL.max() * 1.07)
            # Y-limits driven by data ± errors, not by fit curves, so
            # an ill-conditioned exp fit cannot blow up the scale.
            yvals = G_data[fam][c, :, :]
            evals = sG_data[fam][c, :, :]
            finite = np.isfinite(yvals) & np.isfinite(evals)
            if finite.any():
                lo = float(np.nanmin((yvals - evals)[finite]))
                hi = float(np.nanmax((yvals + evals)[finite]))
                pad = max(0.05 * (hi - lo), 1e-3)
                ax.set_ylim(lo - pad, hi + pad)
            ax.set_xlabel(r"$1/L$")
            if c == 0:
                ax.set_ylabel(rf"$G_{{\rm conn}}(t_k)$  [{_FAM_TITLE[fam]}]")
            geom_str = geoms[sizes[-1]][fam]
            ax.set_title(
                f"{_FAM_TITLE[fam]} — {_CYCLE_TITLE[c]}\n"
                f"L={sizes[-1]} geom={geom_str}",
                fontsize=9)

    # Legend: fit types + lattice size shapes.
    from matplotlib.lines import Line2D
    fit_handles = [
        Line2D([0], [0], color=FIT_COLORS["linear"],      ls="-",  lw=1.8,
               marker=FIT_MARK["linear"],      mfc="w", mec=FIT_COLORS["linear"],
               mew=1.2, ms=7, label="linear  a + b/L"),
        Line2D([0], [0], color=FIT_COLORS["quadratic"],   ls="--", lw=1.8,
               marker=FIT_MARK["quadratic"],   mfc="w", mec=FIT_COLORS["quadratic"],
               mew=1.2, ms=7, label="quadratic  a + b/L + c/L²"),
        Line2D([0], [0], color=FIT_COLORS["power"], ls=":",  lw=1.8,
               marker=FIT_MARK["power"], mfc="w", mec=FIT_COLORS["power"],
               mew=1.2, ms=7, label="power  a + b·(1/L)^A"),
    ]
    size_handles = [
        Line2D([0], [0], color="0.4", ls="", marker=size_markers[L],
               ms=6, mec="k", mew=0.6, label=f"L = {L}")
        for L in sizes
    ]
    leg_ax = axes[0, -1]
    leg1 = leg_ax.legend(handles=fit_handles, fontsize=8, loc="upper left",
                         framealpha=0.9, title="fit type", title_fontsize=8)
    leg_ax.add_artist(leg1)
    leg_ax.legend(handles=size_handles, fontsize=8, loc="upper right",
                  framealpha=0.9, title="lattice size", title_fontsize=8)

    # Colorbar for t_k.
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=T_FRACS[0],
                                                   vmax=T_FRACS[-1]))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.018, pad=0.02)
    cb.set_label(r"sample fraction  $t_k = k/8$")
    cb.set_ticks(T_FRACS)
    cb.set_ticklabels([f"{t:.3f}" for t in T_FRACS])

    fig.tight_layout(rect=[0, 0, 0.96, 0.94])

    plot_dir = os.path.join(out_root, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    out_png = os.path.join(plot_dir,
        f"fss_eighths_r1_{args.r1:.3f}_r2_{args.r2:.3f}.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"[plot] {out_png}")

    # Print compact summary of intercepts.
    print()
    print(f"{'fam':4} {'c':>1} {'k':>1}  {'t_k':>5}   "
          f"{'lin':>10} ± {'σ':>8}   "
          f"{'quad':>10} ± {'σ':>8}   "
          f"{'pow':>10} ± {'σ':>8}   {'A':>5}")
    for row in summary_rows:
        fam, c, ki, tk, a_l, s_l, a_q, s_q, a_p, s_p, A_p = row
        A_str = f"{A_p:5.2f}" if np.isfinite(A_p) else "  -- "
        print(f"{fam:4} {c:>1} {ki+1:>1}  {tk:>5.3f}   "
              f"{a_l:>10.5f} ± {s_l:>8.5f}   "
              f"{a_q:>10.5f} ± {s_q:>8.5f}   "
              f"{a_p:>10.5f} ± {s_p:>8.5f}   {A_str}")
    return out_png


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag",   default="fss_eighths")
    ap.add_argument("--r1",    type=float, default=1.0)
    ap.add_argument("--r2",    type=float, default=1.0)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[8, 12, 16, 20])
    ap.add_argument("--test-Tx-frac", type=float, default=0.0)
    ap.add_argument("--test-Ty-frac", type=float, default=0.0)
    ap.add_argument("--ref-Tx-frac",  type=float, default=0.25)
    ap.add_argument("--ref-Ty-frac",  type=float, default=0.25)
    ap.add_argument("--n-traj",        type=int, default=20000)
    ap.add_argument("--n-traj-coarse", type=int, default=2000)
    ap.add_argument("--n-traj-fine",   type=int, default=4000)
    ap.add_argument("--beta-seed", type=float, nargs=2,
                    default=(0.05, 0.40))
    ap.add_argument("--n-workers", type=int, default=4)
    ap.add_argument("--exe", default=None)
    ap.add_argument("--stages", nargs="+",
                    default=["mc", "plot"],
                    choices=["mc", "plot"])
    args = ap.parse_args()

    geoms = {L: {"test": _make_geom(L, 1.0, 1.0,
                                     args.test_Tx_frac, args.test_Ty_frac),
                  "ref":  _make_geom(L, 1.0, 1.0,
                                     args.ref_Tx_frac, args.ref_Ty_frac)}
             for L in args.sizes}
    out_root = os.path.join(_HERE, "results", args.tag)

    if "mc" in args.stages:
        geoms, out_root = stage_mc(args)
    if "plot" in args.stages:
        stage_plot(args, geoms, out_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
