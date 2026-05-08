#!/usr/bin/env python3
"""
run_fss.py — minimal FSS driver.

For each lattice in LATTICES (each entry is (L_x, L_y, T_x, T_y)):
  1. (optional) locate beta_c via a Gram-Charlier susceptibility scan
  2. run the C++ Ising simulator at beta_c (or BETA_FIXED)
  3. load the connected two-point function

Then plot:
  - the susceptibility scan + GC fit per lattice (beta_c finder visualizer)
  - G_conn(t) along each of the 3 boundary cycles, t in [0, 1] (one panel
    per cycle, all lattice sizes overlaid)
  - FSS extrapolation of the boundary midpoint G_conn(t = 1/2) vs 1/L for
    each cycle, with a linear (a + b/L) and quadratic (a + b/L + c/L^2) fit
    — this is the geometry/coupling-optimization continuum target.

Edit the CONFIG block below and click Run.
"""
# ===========================================================================
# CONFIG
# ===========================================================================
# Each entry: (L_x, L_y, T_x, T_y).  Untwisted = T_x = T_y = 0.
# For a twisted parallelogram set T_x, T_y to integer twist offsets;
# e.g. (16, 16, 4, 4) is a 1/4 twist on a 16x16 lattice.
LATTICES    = [
    (8,  8,  0, 0),
    (12, 12, 0, 0),
    (16, 16, 0, 0),
    (24, 24, 0, 0),
    (32, 32, 0, 0),
    (48, 48, 0, 0),
    (64, 64, 0, 0),
]
K1, K2, K3  = 1.0, 1.0, 1.0         # nearest-neighbour couplings (equilateral)
DELTA_SIGMA = 0.125                 # spin scaling dimension (2D Ising = 1/8)

# beta_c selection ---------------------------------------------------------
FIND_BETA_C = True                  # True → run scan; False → use BETA_FIXED
BETA_FIXED  = 0.27465               # used only when FIND_BETA_C is False
SCAN_LO     = 0.05                  # scan bracket for beta_c finder
SCAN_HI     = 0.40

# beta_c finder MC stats (cheap, used only during the susceptibility scan)
SCAN_N_TRAJ_COARSE = 4000           # traj/point in the initial coarse pass
SCAN_N_TRAJ_FINE   = 8000           # traj/point in the 3 refinement passes
SCAN_N_COARSE      = 13             # # of beta points in the coarse pass
SCAN_N_REFINE      = 7              # # of beta points in each refine pass

# production MC ------------------------------------------------------------
N_TRAJ      = 80000
N_SKIP      = 20
N_THERM     = 4000

TAG         = "hi"                  # output dir name → results/<tag>/

# boundary-cycle sampling --------------------------------------------------
N_T_SAMPLES = 33                    # # of points per cycle in [0, 1]
# ===========================================================================

import os
import sys
import pickle
import subprocess
import time

import numpy as np
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "lib"))
import mc_engine  # noqa: E402
from cost import boundary_paths, _tile_interp, _SQRT3_2  # noqa: E402


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _exe_path():
    base = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    cands = (base + ".exe", base) if os.name == "nt" else (base, base + ".exe")
    for c in cands:
        if os.path.exists(c):
            return c
    sys.exit(f"ERROR: simulator binary not found in {os.path.join(_HERE, 'bin')}.\n"
             f"Build it with `make` from this directory.")


def _label(Lx, Ly, Tx, Ty):
    return (f"L{Lx}x{Ly}" if Tx == 0 and Ty == 0
            else f"L{Lx}x{Ly}_T{Tx}x{Ty}")


def run_one_lattice(Lx, Ly, Tx, Ty, exe, out_pkl):
    """Find beta_c (optionally), run production MC, load two-point data.
    Returns dict with Lx, Ly, Tx, Ty, beta_c, scan_betas, scan_chis, rows."""
    lab = _label(Lx, Ly, Tx, Ty)
    if os.path.exists(out_pkl):
        with open(out_pkl, "rb") as f:
            return pickle.load(f)

    scratch = os.path.join(_HERE, "results", TAG, "_scratch", lab)
    os.makedirs(scratch, exist_ok=True)

    # ---- (1) beta_c
    if FIND_BETA_C:
        print(f"[scan] {lab:>14}  finding beta_c in [{SCAN_LO}, {SCAN_HI}] ...",
              flush=True)
        t0 = time.time()
        beta_c, _, chi_peak, sb, sc, _ = mc_engine.find_beta_c(
            exe, Lx, Ly, Tx, Ty, K1, K2, K3, SCAN_LO, SCAN_HI,
            n_coarse=SCAN_N_COARSE,
            n_refine=SCAN_N_REFINE,
            n_refine2=SCAN_N_REFINE,
            n_refine3=SCAN_N_REFINE,
            n_traj_coarse=SCAN_N_TRAJ_COARSE,
            n_traj_fine=SCAN_N_TRAJ_FINE,
            data_dir=os.path.join(scratch, "scan"))
        print(f"[scan] {lab:>14}  beta_c = {beta_c:.5f}  "
              f"chi_peak = {chi_peak:.2f}  ({time.time()-t0:.1f}s)", flush=True)
    else:
        beta_c = BETA_FIXED
        sb, sc = [], []
        print(f"[scan] {lab:>14}  using BETA_FIXED = {beta_c:.5f}", flush=True)

    # ---- (2) production MC at beta_c
    print(f"[mc]   {lab:>14}  running production at beta = {beta_c:.5f} ...",
          flush=True)
    t0 = time.time()
    _, subdir = mc_engine.run_simulator(
        exe, Lx, Ly, Tx, Ty, K1, K2, K3, beta_c,
        n_traj=N_TRAJ, n_skip=N_SKIP, n_therm=N_THERM,
        data_dir=os.path.join(scratch, "prod"))
    if subdir is None:
        sys.exit(f"could not locate simulator output for {lab}")

    # ---- (3) load two-point file
    rows = []
    with open(os.path.join(subdir, "two_point_all_to_all.dat")) as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.split()
            if len(p) >= 7:
                rows.append((int(p[1]), int(p[2]),
                             float(p[5]), float(p[6])))  # (m, n, conn, conn_err)
    print(f"[mc]   {lab:>14}  done in {time.time()-t0:.1f}s  "
          f"({len(rows)} two-point entries)", flush=True)

    # Convert rows → dict expected by _tile_interp.
    data = {(m, n): {"conn": g, "conn_err": e} for (m, n, g, e) in rows}
    payload = {"Lx": Lx, "Ly": Ly, "Tx": Tx, "Ty": Ty,
               "beta_c": float(beta_c),
               "scan_betas": list(sb), "scan_chis": list(sc),
               "data": data}
    with open(out_pkl, "wb") as f:
        pickle.dump(payload, f)
    return payload


def sample_boundary_cycles(payload, t_fracs):
    """Sample G_conn along the 3 boundary cycles at fractions t in [0, 1].
    Returns G[c, k], sG[c, k] arrays of shape (3, len(t_fracs))."""
    Lx, Ly = payload["Lx"], payload["Ly"]
    Tx, Ty = payload["Tx"], payload["Ty"]
    data = payload["data"]
    iG = _tile_interp(data, Lx, Ly, Tx, Ty, "conn",     copies=2)
    iE = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    G  = np.zeros((3, len(t_fracs)))
    sG = np.zeros_like(G)
    for c, (dm, dn) in enumerate(paths):
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        pts = np.column_stack([t_fracs * ex, t_fracs * ey])
        G[c]  = np.asarray(iG(pts), dtype=float).ravel()
        sG[c] = np.abs(np.asarray(iE(pts), dtype=float).ravel())
    return G, sG


def _wls_poly(invL, y, sigma, deg):
    """Weighted least-squares polynomial fit; return (intercept, sigma, coeffs)."""
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    x, yv, sv = invL[mask], y[mask], sigma[mask]
    if len(x) < deg + 1:
        return np.nan, np.nan, None
    X = np.vander(x, deg + 1, increasing=True)
    w = 1.0 / sv ** 2
    try:
        cov = np.linalg.inv((X.T * w) @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, None
    beta = cov @ ((X.T * w) @ yv)
    return float(beta[0]), float(np.sqrt(max(cov[0, 0], 0.0))), beta


def _power_law_fit(Larr, y, sigma):
    """Fit G = A * L^{-B} via WLS in log-log space.

    Error propagation: sigma_{log G} = sigma_G / G  (first order).
    Returns (A, sigma_A, B, sigma_B) or (nan, nan, nan, nan) on failure.
    """
    mask = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0) & (y > 0)
    L, yv, sv = Larr[mask], y[mask], sigma[mask]
    if len(L) < 2:
        return np.nan, np.nan, np.nan, np.nan
    logL = np.log(L)
    logY = np.log(yv)
    slogY = sv / yv                    # sigma_{log G} ≈ sigma_G / G
    w = 1.0 / slogY ** 2
    # Design matrix: [1, logL] → [logA, -B]
    X = np.column_stack([np.ones_like(logL), logL])
    try:
        cov = np.linalg.inv((X.T * w) @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan
    coeffs = cov @ ((X.T * w) @ logY)
    logA, neg_B = float(coeffs[0]), float(coeffs[1])
    A  = np.exp(logA)
    B  = -neg_B
    slogA = float(np.sqrt(max(cov[0, 0], 0.0)))
    sB    = float(np.sqrt(max(cov[1, 1], 0.0)))
    sA    = A * slogA                  # delta A = A * delta(log A)
    return A, sA, B, sB


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

def plot_beta_c_scan(payloads, out_png):
    """Visualise the beta_c finder: chi(beta) scan points + GC fit + peak."""
    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.get_cmap("viridis")
    for i, p in enumerate(payloads):
        if not p["scan_betas"]:
            continue
        col = cmap(i / max(len(payloads) - 1, 1))
        sb = np.array(p["scan_betas"])
        sc = np.array(p["scan_chis"])
        order = np.argsort(sb)
        sb, sc = sb[order], sc[order]
        ax.plot(sb, sc, "o", ms=4, color=col, alpha=0.6,
                label=f"{_label(p['Lx'], p['Ly'], p['Tx'], p['Ty'])}  "
                      f"(β_c = {p['beta_c']:.4f})")
        # Smooth GC curve through the scan points (best-effort).
        try:
            from mc_engine import _gc_fit, _gram_charlier
            params, _ = _gc_fit(list(sb), list(sc), p["beta_c"])
            if params is not None:
                bb = np.linspace(sb.min(), sb.max(), 400)
                ax.plot(bb, _gram_charlier(bb, *params),
                        "-", lw=1.0, color=col, alpha=0.8)
        except Exception:
            pass
        ax.axvline(p["beta_c"], color=col, ls=":", lw=0.8, alpha=0.8)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"susceptibility  $\chi(\beta)$")
    ax.set_title("β_c finder — susceptibility scan + Gram-Charlier fit")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"[plot] saved {out_png}")
    plt.show()


def plot_boundary_cycles(payloads, out_png):
    """G_conn vs t along each of the 3 boundary cycles, all lattices overlaid."""
    t = np.linspace(0.0, 1.0, N_T_SAMPLES)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    cmap = plt.get_cmap("viridis")
    for i, p in enumerate(payloads):
        col = cmap(i / max(len(payloads) - 1, 1))
        G, sG = sample_boundary_cycles(p, t)
        for c in range(3):
            axes[c].errorbar(t, G[c], yerr=sG[c], fmt="o-", ms=3, lw=1,
                             capsize=2, color=col,
                             label=_label(p["Lx"], p["Ly"], p["Tx"], p["Ty"]))
    for c, ax in enumerate(axes):
        ax.set_xlabel(r"$t$ along cycle")
        ax.set_title(f"boundary cycle {c}")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"$G_{\rm conn}(t)$")
    axes[-1].legend(fontsize=8, loc="best")
    fig.suptitle(f"Boundary correlators   (k1, k2, k3) = ({K1:g}, {K2:g}, {K3:g})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=150)
    print(f"[plot] saved {out_png}")
    plt.show()


def plot_midpoint_fss(payloads, out_png):
    """Theory-agnostic FSS extrapolation of the boundary midpoint.

    Fits G_c(L/2) = A * L^{-B} in log-log space for each of the three
    boundary cycles.  No scaling dimension is assumed.

    The optimizer target is amplitude isotropy: all three A values should
    agree when the lattice geometry correctly reproduces the continuum CFT.
    The fitted exponent B is a bonus check — for 2D Ising it should converge
    to 2*Delta = 0.25, but the optimizer does not need to know that.

    A single-panel comparison plot with all three cycles overlaid (log-log)
    is produced alongside the per-cycle detail panels.
    """
    t_mid = np.array([0.5])
    Larr  = np.array([p["Lx"] for p in payloads], dtype=float)
    invL  = 1.0 / Larr
    G_mid  = np.full((3, len(payloads)), np.nan)
    sG_mid = np.full_like(G_mid, np.nan)
    for i, p in enumerate(payloads):
        G, sG = sample_boundary_cycles(p, t_mid)
        G_mid[:, i]  = G[:, 0]
        sG_mid[:, i] = sG[:, 0]

    cyc_color = ["#1f77b4", "#d62728", "#2ca02c"]
    L_dense   = np.geomspace(Larr.min() * 0.85, Larr.max() * 1.15, 300)
    summary   = []

    # ---- 4-panel layout: 3 detail + 1 log-log comparison ----
    fig = plt.figure(figsize=(18, 5))
    detail_axes = [fig.add_subplot(1, 4, c + 1) for c in range(3)]
    ax_ll = fig.add_subplot(1, 4, 4)

    for c in range(3):
        col = cyc_color[c]
        ax  = detail_axes[c]

        A, sA, B, sB = _power_law_fit(Larr, G_mid[c], sG_mid[c])
        summary.append((A, sA, B, sB))

        ax.errorbar(invL, G_mid[c], yerr=sG_mid[c], fmt="o", ms=6,
                    capsize=3, color=col, mec="k", mew=0.4,
                    label="data", zorder=4)
        if np.isfinite(A):
            ax.plot(1.0 / L_dense, A * L_dense ** (-B), "-", color="k", lw=1.5,
                    label=fr"$A L^{{-B}}$"
                          fr"  $A={A:.4f}\pm{sA:.4f}$"
                          fr"  $B={B:.4f}\pm{sB:.4f}$")
            # mark the L→∞ intercept at 1/L=0 — meaningful only as a visual ref
            ax.axvline(0.0, color="k", lw=0.6, ls=":", alpha=0.4)
        ax.set_xlabel(r"$1/L$")
        ax.set_ylabel(r"$G_c(L/2)$")
        ax.set_title(f"cycle {c}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper left")

        # log-log comparison panel
        ax_ll.errorbar(Larr, G_mid[c], yerr=sG_mid[c], fmt="o", ms=5,
                       capsize=2, color=col, mec="k", mew=0.3,
                       label=f"cycle {c}" + (fr"  $B={B:.3f}$" if np.isfinite(B) else ""),
                       zorder=4)
        if np.isfinite(A):
            ax_ll.plot(L_dense, A * L_dense ** (-B), "-", color=col,
                       lw=1.2, alpha=0.7)

    ax_ll.set_xscale("log"); ax_ll.set_yscale("log")
    ax_ll.set_xlabel(r"$L$")
    ax_ll.set_ylabel(r"$G_c(L/2)$")
    ax_ll.set_title("log-log: all cycles")
    ax_ll.grid(alpha=0.3, which="both")
    ax_ll.legend(fontsize=8)

    detail_axes[0].set_ylabel(r"$G_c(L/2)$")
    fig.suptitle("Boundary-midpoint FSS  —  $G_c(L/2) = A\\,L^{-B}$  (theory-agnostic)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=150)
    print(f"[plot] saved {out_png}")

    print("\n=== Power-law fit  G_c(L/2) = A * L^{-B} ===")
    print(f"     {'cycle':>6}   {'A':>18}   {'B':>18}   note")
    for c, (A, sA, B, sB) in enumerate(summary):
        note = f"B/2 = {B/2:.4f}" if np.isfinite(B) else ""
        print(f"     {c:>6d}   {A:>8.5f} ± {sA:<7.5f}   "
              f"{B:>8.5f} ± {sB:<7.5f}   {note}")
    A_arr = np.array([s[0] for s in summary])
    B_arr = np.array([s[2] for s in summary])
    if np.all(np.isfinite(A_arr)):
        print(f"\n     amplitude anisotropy  max|A_c - <A>|/<A> = "
              f"{np.max(np.abs(A_arr - A_arr.mean()))/abs(A_arr.mean()):.3e}  "
              f"← optimizer drives this to 0")
        print(f"     mean fitted exponent  <B> = {B_arr.mean():.4f}  "
              f"(continuum value = 2*Delta, geometry-independent)")
    plt.show()


def plot_fss_scaling(payloads, out_png,
                     t_fracs=None,
                     L_select=(8, 16, 24, 32, 48, 64)):
    """Every-eighth correlator FSS quality check.

    For each boundary cycle and each t-fraction t_k = k/8 (k=1..7):
      - fits  G_c(t_k, L) = A(t_k) * L^{-B(t_k)}  in log-log space
      - top row:    raw G vs 1/L with power-law curve overlaid
      - bottom row: L^{B(t_k)} * G vs 1/L  — should be flat at A(t_k)
                    if the pure power-law holds.  Curvature = corrections.

    L_select filters payloads to a specific set of sizes (default: 8,16,24,32,48,64).
    """
    if t_fracs is None:
        t_fracs = np.arange(1, 8) / 8.0   # 1/8, 2/8, …, 7/8

    plist = sorted([p for p in payloads if p["Lx"] in set(L_select)],
                   key=lambda p: p["Lx"])
    if not plist:
        print("[plot] plot_fss_scaling: no payloads match L_select, skipping.")
        return
    Larr = np.array([p["Lx"] for p in plist], dtype=float)
    invL = 1.0 / Larr

    # G_all[c, t_idx, L_idx]
    G_all  = np.full((3, len(t_fracs), len(plist)), np.nan)
    sG_all = np.full_like(G_all, np.nan)
    for i, p in enumerate(plist):
        G, sG = sample_boundary_cycles(p, t_fracs)
        G_all[:, :, i]  = G
        sG_all[:, :, i] = sG

    cmap    = plt.get_cmap("plasma")
    L_dense = np.geomspace(Larr.min() * 0.9, Larr.max() * 1.1, 300)
    cyc_names = ["cycle 0", "cycle 1", "cycle 2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 8))

    for c in range(3):
        ax_r = axes[0, c]
        ax_s = axes[1, c]

        for ti, t in enumerate(t_fracs):
            col = cmap(ti / max(len(t_fracs) - 1, 1))
            y   = G_all[c, ti]
            sy  = sG_all[c, ti]
            A, sA, B, sB = _power_law_fit(Larr, y, sy)
            lbl = f"$t={t:.3f}$"

            # top: raw G vs 1/L + power-law curve
            ax_r.errorbar(invL, y, yerr=sy, fmt="o", ms=4, capsize=2,
                          color=col, mec="k", mew=0.3, zorder=4)
            if np.isfinite(A) and np.isfinite(B):
                ax_r.plot(1.0 / L_dense, A * L_dense**(-B),
                          "-", color=col, lw=1.0, alpha=0.8, label=lbl)

            # bottom: L^B * G vs 1/L — flat at A means scaling holds
            if np.isfinite(B):
                scaled  = y  * Larr**B
                sscaled = sy * Larr**B
                ax_s.errorbar(invL, scaled, yerr=sscaled, fmt="o", ms=4,
                              capsize=2, color=col, mec="k", mew=0.3,
                              zorder=4, label=lbl)
                if np.isfinite(A):
                    ax_s.axhline(A, color=col, lw=0.8, ls="--", alpha=0.6)

        ax_r.set_title(f"{cyc_names[c]}: $G_c(t,L)$ vs $1/L$")
        ax_r.set_ylabel(r"$G_c(t,\,L)$")
        ax_r.grid(alpha=0.3)
        ax_r.legend(fontsize=6, ncol=2, loc="upper left")

        ax_s.set_title(f"{cyc_names[c]}: $L^B\\,G_c$ vs $1/L$  (flat = good)")
        ax_s.set_xlabel(r"$1 / L$")
        ax_s.set_ylabel(r"$L^{B(t)}\,G_c(t,\,L)$")
        ax_s.grid(alpha=0.3)
        ax_s.legend(fontsize=6, ncol=2, loc="best")

    L_str = ", ".join(str(int(l)) for l in sorted(L_select))
    fig.suptitle(
        f"Power-law scaling quality  $G_c(t,L) = A(t)\\,L^{{-B(t)}}$\n"
        f"t = k/8, k=1…7 · sizes L ∈ {{{L_str}}}\n"
        f"Bottom: $L^B G$ vs $1/L$ — flat dashes at $A(t)$ if correction-free",
        fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig(out_png, dpi=150)
    print(f"[plot] saved {out_png}")
    plt.show()


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    exe = _exe_path()
    out_root = os.path.join(_HERE, "results", TAG)
    os.makedirs(out_root, exist_ok=True)

    payloads = [run_one_lattice(Lx, Ly, Tx, Ty, exe,
                                os.path.join(out_root,
                                             f"{_label(Lx, Ly, Tx, Ty)}.pkl"))
                for (Lx, Ly, Tx, Ty) in LATTICES]

    if FIND_BETA_C and any(p["scan_betas"] for p in payloads):
        plot_beta_c_scan(payloads, os.path.join(out_root, "beta_c_scan.png"))
    plot_boundary_cycles(payloads, os.path.join(out_root, "boundary_cycles.png"))
    plot_midpoint_fss(payloads, os.path.join(out_root, "midpoint_fss.png"))
    plot_fss_scaling(payloads, os.path.join(out_root, "fss_scaling_quality.png"))


if __name__ == "__main__":
    main()
