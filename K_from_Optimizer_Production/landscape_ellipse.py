"""
landscape_ellipse.py — fit ellipses to the half-decay contour of the
2D correlator G(x, y) and use ellipse parameters as the cost.

For each lattice (test or reference) we build a 2D interpolator on
the connected correlator and find r_half(θ) = the radius along
direction θ where G drops to G_max / 2. Because the connected
correlator is symmetric, G(x,y) = G(-x,-y), the half-decay contour
is centred at the origin (sample θ ∈ [0, π)). A centred ellipse

    M11 x² + 2 M12 x y + M22 y² = 1

inverts to

    1 / r_half(θ)² = M11 cos²θ + 2 M12 cosθ sinθ + M22 sin²θ

which is linear in (M11, M12, M22). One least-squares solve per
lattice gives M; eigenvalues of M give the principal half-axes
(a, b) = (1/√λ_min, 1/√λ_max) and the eigenvector direction gives
the orientation angle φ.

Six cost kernels:
  1. ellipse_axes      (Δa)² + (Δb)² + sin²(Δφ)
  2. ellipse_axes_log  (Δlog a)² + (Δlog b)² + sin²(Δφ)
  3. ellipse_M_frob    Frobenius || M_test − M_ref ||²  (matrix L2)
  4. ellipse_M_log     Frobenius || log(M_test) − log(M_ref) ||²
  5. ellipse_area_eps  (Δlog area)² + (Δlog aspect)² + sin²(Δφ)
                       — area/aspect/orientation factored
  6. ellipse_residual  Predict r_half_test at each θ from M_ref,
                       compare in log to actual r_half_test.

Loads precompute pkl grid; writes one combined PNG.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine  # noqa: E402
from cost import _tile_interp  # noqa: E402

_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
    "_reference_Lx39_Ly48_Tx-9_Ty9": (39, 48, -9, 9),
    "_reference_Lx16_Ly16_Tx0_Ty0":  (16, 16, 0, 0),
}
_TRUTH_456 = (5.0652, 7.7429)


# --------------------- ellipse extraction ---------------------------- #
def _gmax_estimate(interp, dx_probe=1.0, n_probe=12):
    """Estimate G(0+) by averaging interp at small radii in many directions."""
    th = np.linspace(0, 2 * np.pi, n_probe, endpoint=False)
    pts = np.column_stack([dx_probe * np.cos(th), dx_probe * np.sin(th)])
    vals = np.asarray(interp(pts), dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) == 0:
        return np.nan
    return float(np.mean(vals))


def _r_half_at_theta(interp, theta, target,
                     r_min=0.5, r_max=14.0, n_steps=120):
    """Find smallest r > r_min where interp(r cosθ, r sinθ) ≤ target.
    Linear interp in (r, log G) between bracketing samples.
    Returns NaN if no crossing in [r_min, r_max].
    """
    if target <= 0:
        return np.nan
    rs = np.linspace(r_min, r_max, n_steps)
    pts = np.column_stack([rs * np.cos(theta), rs * np.sin(theta)])
    g = np.asarray(interp(pts), dtype=float)
    # mask non-finite or non-positive
    ok = np.isfinite(g) & (g > 0)
    if ok.sum() < 2:
        return np.nan
    # find first sample where g <= target
    below = np.where(ok & (g <= target))[0]
    if len(below) == 0:
        return np.nan
    j = below[0]
    if j == 0:
        return rs[0]
    if not ok[j - 1]:
        return rs[j]
    g0, g1 = g[j - 1], g[j]
    r0, r1 = rs[j - 1], rs[j]
    if g0 == g1 or g0 <= 0 or g1 <= 0:
        return 0.5 * (r0 + r1)
    lt = np.log(target)
    alpha = (lt - np.log(g0)) / (np.log(g1) - np.log(g0))
    return r0 + alpha * (r1 - r0)


def _fit_centred_ellipse(interp, r_min=0.5, r_max=14.0, n_theta=24,
                         half_target=None):
    """Return dict(M, a, b, phi, ok) by fitting a centred ellipse to the
    half-decay contour. half_target overrides G_max/2 if supplied."""
    if half_target is None:
        gmax = _gmax_estimate(interp, dx_probe=1.0)
        if not np.isfinite(gmax):
            return None
        target = 0.5 * gmax
    else:
        target = float(half_target)

    thetas = np.linspace(0.0, np.pi, n_theta, endpoint=False)
    rh = np.array([_r_half_at_theta(interp, th, target,
                                    r_min=r_min, r_max=r_max)
                   for th in thetas])
    ok = np.isfinite(rh) & (rh > 0)
    if ok.sum() < 5:
        return None
    th = thetas[ok]
    rh = rh[ok]
    y = 1.0 / (rh ** 2)
    c = np.cos(th); s = np.sin(th)
    A = np.column_stack([c * c, 2.0 * c * s, s * s])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    M11, M12, M22 = coef
    M = np.array([[M11, M12], [M12, M22]], dtype=float)
    # eigen-decompose: ensure SPD-ish
    w, v = np.linalg.eigh(M)
    if np.any(w <= 0):
        return None
    # axes: r_half along principal axis = 1/sqrt(eigenvalue)
    # smaller eigenvalue → larger axis → 'a' (semi-major)
    order = np.argsort(w)
    a = 1.0 / np.sqrt(w[order[0]])
    b = 1.0 / np.sqrt(w[order[1]])
    # orientation of semi-major axis
    vec = v[:, order[0]]
    phi = np.arctan2(vec[1], vec[0])
    # canonicalise to [0, π)
    phi = phi % np.pi
    return dict(M=M, a=float(a), b=float(b), phi=float(phi),
                target=target, n_used=int(ok.sum()))


# --------------------- cost kernels ---------------------------------- #
def _angle_diff(p1, p2):
    """Smallest angle between two orientations in [0, π)."""
    d = (p1 - p2 + np.pi / 2) % np.pi - np.pi / 2
    return float(d)


def k_ellipse_axes(test, ref):
    if test is None or ref is None: return np.nan
    da = test["a"] - ref["a"]; db = test["b"] - ref["b"]
    dp = _angle_diff(test["phi"], ref["phi"])
    return da * da + db * db + np.sin(dp) ** 2


def k_ellipse_axes_log(test, ref):
    if test is None or ref is None: return np.nan
    da = np.log(test["a"]) - np.log(ref["a"])
    db = np.log(test["b"]) - np.log(ref["b"])
    dp = _angle_diff(test["phi"], ref["phi"])
    return da * da + db * db + np.sin(dp) ** 2


def k_ellipse_M_frob(test, ref):
    if test is None or ref is None: return np.nan
    D = test["M"] - ref["M"]
    return float(np.sum(D * D))


def k_ellipse_M_log(test, ref):
    if test is None or ref is None: return np.nan
    # matrix log via eigendecomp (M is SPD by construction)
    def mlog(X):
        w, v = np.linalg.eigh(X)
        return v @ np.diag(np.log(w)) @ v.T
    D = mlog(test["M"]) - mlog(ref["M"])
    return float(np.sum(D * D))


def k_ellipse_area_aspect(test, ref):
    if test is None or ref is None: return np.nan
    area_t = np.pi * test["a"] * test["b"]
    area_r = np.pi * ref["a"]  * ref["b"]
    asp_t  = test["a"] / test["b"]
    asp_r  = ref["a"]  / ref["b"]
    da = np.log(area_t) - np.log(area_r)
    dr = np.log(asp_t)  - np.log(asp_r)
    dp = _angle_diff(test["phi"], ref["phi"])
    return da * da + dr * dr + np.sin(dp) ** 2


def k_ellipse_residual(test_interp, ref_M, half_target,
                       r_min=0.5, r_max=14.0, n_theta=24):
    """Sample r_half_test at each θ via interp; predict r_half from
    ref_M at the same θ; compare in log."""
    if test_interp is None or ref_M is None:
        return np.nan
    thetas = np.linspace(0.0, np.pi, n_theta, endpoint=False)
    rh_t = np.array([_r_half_at_theta(test_interp, th, half_target,
                                      r_min=r_min, r_max=r_max)
                     for th in thetas])
    ok = np.isfinite(rh_t) & (rh_t > 0)
    if ok.sum() < 4:
        return np.nan
    th = thetas[ok]
    c = np.cos(th); s = np.sin(th)
    M11, M22 = ref_M[0, 0], ref_M[1, 1]; M12 = ref_M[0, 1]
    inv_r2 = M11 * c * c + 2 * M12 * c * s + M22 * s * s
    if np.any(inv_r2 <= 0):
        return np.nan
    rh_r_pred = 1.0 / np.sqrt(inv_r2)
    d = np.log(rh_t[ok]) - np.log(rh_r_pred)
    return float(np.mean(d * d))


# --------------------- driver ---------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag",   default="_reference_Lx39_Ly48_Tx-9_Ty9")
    ap.add_argument("--out",       default=None)
    ap.add_argument("--n-theta",   type=int, default=24)
    ap.add_argument("--r-min",     type=float, default=0.5)
    ap.add_argument("--r-max",     type=float, default=14.0)
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    ref_geom = _REF_GEOMS[args.ref_tag]
    ref_dir  = os.path.join(_HERE, "results", args.ref_tag)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    rL, rL2, rT, rT2 = ref_geom
    iref = _tile_interp(ref_data, rL, rL2, rT, rT2, "conn", 2)
    ref_fit = _fit_centred_ellipse(iref, r_min=args.r_min, r_max=args.r_max,
                                   n_theta=args.n_theta)
    if ref_fit is None:
        print("[error] ref ellipse fit failed"); return 1
    print(f"[ref] M=\n{ref_fit['M']}\n   a={ref_fit['a']:.3f}  "
          f"b={ref_fit['b']:.3f}  phi={np.degrees(ref_fit['phi']):.1f}°  "
          f"target={ref_fit['target']:.3f}  n_used={ref_fit['n_used']}")

    print(f"[load+fit] {n*n} grid pts ...")
    t0 = time.time()
    fits = {}
    interps = {}
    skipped = 0
    for r1 in rs:
        for r2 in rs:
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl):
                continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            it = _tile_interp(pt["test_data"], Lx, Ly, Tx, Ty, "conn", 2)
            f_ = _fit_centred_ellipse(it, r_min=args.r_min,
                                      r_max=args.r_max, n_theta=args.n_theta)
            if f_ is None:
                skipped += 1
                continue
            fits[(float(r1), float(r2))] = f_
            interps[(float(r1), float(r2))] = it
    print(f"[load+fit] {len(fits)} pts (skipped {skipped}) in "
          f"{time.time()-t0:.1f}s")

    KERNELS = [
        ("ellipse_axes",     lambda t: k_ellipse_axes(t, ref_fit)),
        ("ellipse_axes_log", lambda t: k_ellipse_axes_log(t, ref_fit)),
        ("ellipse_M_frob",   lambda t: k_ellipse_M_frob(t, ref_fit)),
        ("ellipse_M_log",    lambda t: k_ellipse_M_log(t, ref_fit)),
        ("ellipse_area_asp", lambda t: k_ellipse_area_aspect(t, ref_fit)),
    ]

    results = {}
    print("[eval] kernels:")
    for name, fn in KERNELS:
        Z = np.full((n, n), np.nan)
        t0 = time.time()
        for i, r1 in enumerate(rs):
            for j, r2 in enumerate(rs):
                f_ = fits.get((float(r1), float(r2)))
                try:
                    Z[i, j] = fn(f_)
                except Exception:
                    pass
        if np.all(np.isnan(Z)):
            print(f"  {name:20s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        results[name] = dict(Z=Z, argmin=(rs[ij[0]], rs[ij[1]]),
                             d_truth=d_truth,
                             cost_truth=float(Z[truth_i, truth_j]),
                             cost_argmin=float(Z[ij]))
        print(f"  {name:20s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:5.2f}  cost(truth)={Z[truth_i,truth_j]:.3e}  "
              f"cost(min)={Z[ij]:.3e}  wall={time.time()-t0:.1f}s")

    # residual cost — needs interp + ref_fit
    Z = np.full((n, n), np.nan)
    t0 = time.time()
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            it = interps.get((float(r1), float(r2)))
            try:
                Z[i, j] = k_ellipse_residual(it, ref_fit["M"],
                                             ref_fit["target"],
                                             r_min=args.r_min,
                                             r_max=args.r_max,
                                             n_theta=args.n_theta)
            except Exception:
                pass
    if not np.all(np.isnan(Z)):
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        results["ellipse_residual"] = dict(
            Z=Z, argmin=(rs[ij[0]], rs[ij[1]]),
            d_truth=d_truth, cost_truth=float(Z[truth_i, truth_j]),
            cost_argmin=float(Z[ij]))
        print(f"  {'ellipse_residual':20s} argmin=({rs[ij[0]]:5.2f},"
              f"{rs[ij[1]]:5.2f})  d_truth={d_truth:5.2f}  "
              f"cost(truth)={Z[truth_i,truth_j]:.3e}  "
              f"cost(min)={Z[ij]:.3e}  wall={time.time()-t0:.1f}s")

    # render
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    K = len(results)
    cols = 3
    rows = (K + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 3.8 * rows))
    axes = axes.flatten()
    for ax, (name, info) in zip(axes, results.items()):
        Z = np.log10(np.maximum(info["Z"], 1e-15))
        vmin, vmax = np.nanpercentile(Z, [2, 98])
        im = ax.pcolormesh(R1, R2, Z, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=11, mec="k", mew=1)
        am = info["argmin"]
        ax.plot(am[0], am[1], "w*", ms=12, mec="k", mew=1)
        ax.set_title(f"{name}\nargmin=({am[0]:.1f},{am[1]:.1f}) "
                     f"d={info['d_truth']:.1f}", fontsize=10)
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_aspect("equal")
    for ax in axes[K:]: ax.axis("off")
    fig.suptitle(f"Ellipse-fit cost kernels  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)\n"
                 f"ref ellipse: a={ref_fit['a']:.2f} b={ref_fit['b']:.2f} "
                 f"φ={np.degrees(ref_fit['phi']):.1f}°",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = args.out or os.path.join(root, "ellipse.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n[done] figure → {out}")

    npz_out = (os.path.splitext(out)[0]) + ".npz"
    np.savez(npz_out, R1=R1, R2=R2, rs=rs,
             ref_M=ref_fit["M"], ref_a=ref_fit["a"],
             ref_b=ref_fit["b"], ref_phi=ref_fit["phi"],
             **{k: v["Z"] for k, v in results.items()})
    print(f"[done] data   → {npz_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
