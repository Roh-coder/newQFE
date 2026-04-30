#!/usr/bin/env python3
"""
paraboloid_interp_test.py — conceptual test of the "interpolate-the-dense-
reference, sample-the-test-exactly" idea, using a smooth analytic function
(a paraboloid) so the only error source is interpolation itself.

Setup
-----
  f(x, y) = x² + y²

  Domain: circular disk of radius R = 1.

  Reference: an N_ref × N_ref Cartesian grid restricted to the disk, then
             rotated by 30° about the z-axis.  Values are exact f(x',y').

  Test:     a fixed 10 × 10 Cartesian grid restricted to the disk, NOT
            rotated.  Values are exact f(x,y).

  We build a LinearNDInterpolator over the (rotated) reference points
  (mirroring the cost.py choice) and evaluate it at every test-grid
  position to obtain G_ref(x_test, y_test).  The residual at each test
  point is

      r_i = f(x_i, y_i)  –  G_ref_interp(x_i, y_i)

  We sweep N_ref ∈ {5, 10, 20, 50, 100, 200} and report

      • RMS residual  (= sqrt(mean(r²)))
      • max residual
      • residual at the disk centre  (where interp is best)
      • residual at the disk edge    (where interp is worst)

  We also compare scipy's LinearND (k=1, what cost.py uses) against
  CloughTocher2D (k=3, smoother) so we can see whether linear interp is
  the bottleneck.

  A figure paraboloid_interp_test.png shows residual maps for each
  reference size.

Conclusion the test is designed to surface
------------------------------------------
  Because f is smooth and analytic, a LinearND interp residual scales as
  O(h²) where h is the reference grid spacing.  For h = 2R/N_ref we
  expect RMS residual ∝ 1/N_ref².  If the actual residual flattens out
  for large N_ref, something other than linear-interp error dominates
  (e.g. boundary / hull artifacts).  Comparing linear vs cubic shows
  whether LinearND is the limiting factor for the cost function.

  This isolates the *geometric* contribution to the residual mismatch
  and lets us calibrate how big the reference correlator grid needs to
  be before interpolation noise drops below the MC noise (~ 1e-3 in the
  4-5-6 problem).

Usage:  python paraboloid_interp_test.py [--ref-sizes 5 10 20 50 100 200]
                                          [--n-test 10] [--show]
"""
from __future__ import annotations

import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator


def make_disk_grid(N: int, R: float = 1.0) -> np.ndarray:
    """Cartesian N×N grid clipped to the disk of radius R; returns (M,2)."""
    xs = np.linspace(-R, R, N)
    ys = np.linspace(-R, R, N)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([X.ravel(), Y.ravel()])
    pts = pts[np.hypot(pts[:, 0], pts[:, 1]) <= R + 1e-12]
    return pts


def rotate(pts: np.ndarray, angle_deg: float) -> np.ndarray:
    a = math.radians(angle_deg)
    c, s = math.cos(a), math.sin(a)
    R = np.array([[c, -s], [s, c]])
    return pts @ R.T


def paraboloid(p: np.ndarray) -> np.ndarray:
    return p[..., 0] ** 2 + p[..., 1] ** 2


def run_one(N_ref: int, test_pts: np.ndarray, R: float, angle_deg: float):
    ref_pts_unrot = make_disk_grid(N_ref, R)
    ref_pts       = rotate(ref_pts_unrot, angle_deg)
    ref_vals      = paraboloid(ref_pts_unrot)   # f is rotation-invariant,
                                                # use either pre- or post-rot
    # (paraboloid is invariant under z-rotation, so f(rotated)=f(unrotated))

    interp_lin   = LinearNDInterpolator(ref_pts, ref_vals)
    interp_cub   = CloughTocher2DInterpolator(ref_pts, ref_vals)

    g_lin = interp_lin(test_pts)
    g_cub = interp_cub(test_pts)
    g_true = paraboloid(test_pts)

    res_lin = g_lin - g_true
    res_cub = g_cub - g_true
    return ref_pts, res_lin, res_cub


def summary(name: str, res: np.ndarray, test_pts: np.ndarray):
    finite = np.isfinite(res)
    if finite.sum() == 0:
        return dict(rms=float("nan"), max=float("nan"),
                    n_finite=0, n_total=len(res))
    r = res[finite]
    rho = np.hypot(test_pts[finite, 0], test_pts[finite, 1])
    centre_mask = rho < 0.3
    edge_mask   = rho > 0.85
    return dict(
        rms=float(np.sqrt(np.mean(r ** 2))),
        max=float(np.max(np.abs(r))),
        rms_centre=float(np.sqrt(np.mean(r[centre_mask] ** 2)))
                    if centre_mask.sum() else float("nan"),
        rms_edge  =float(np.sqrt(np.mean(r[edge_mask]   ** 2)))
                    if edge_mask.sum() else float("nan"),
        n_finite=int(finite.sum()),
        n_total=int(len(res)),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--ref-sizes", nargs="+", type=int,
                   default=[5, 10, 20, 50, 100, 200])
    p.add_argument("--n-test", type=int, default=10)
    p.add_argument("--R", type=float, default=1.0)
    p.add_argument("--angle-deg", type=float, default=30.0)
    p.add_argument("--out-png", type=str, default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "paraboloid_interp_test.png"))
    args = p.parse_args()

    test_pts = make_disk_grid(args.n_test, args.R)
    print(f"Test grid: {args.n_test}×{args.n_test} → {len(test_pts)} pts in disk")
    print(f"Reference rotation: {args.angle_deg}°  domain radius R={args.R}")
    print(f"f(x,y) = x² + y²,  exact at all sample points\n")

    print(f"{'N_ref':>6} {'h_ref':>8}  "
          f"{'RMS_lin':>10} {'RMS_cub':>10}  "
          f"{'RMS_lin/h²':>11} {'RMS_lin_ctr':>11} {'RMS_lin_edge':>12}  "
          f"{'n_fin/n':>8}")
    print("-" * 105)

    all_results = []
    for N in args.ref_sizes:
        h = 2 * args.R / N
        ref_pts, res_lin, res_cub = run_one(N, test_pts, args.R, args.angle_deg)
        s_lin = summary("lin", res_lin, test_pts)
        s_cub = summary("cub", res_cub, test_pts)
        all_results.append(dict(
            N=N, h=h, ref_pts=ref_pts,
            res_lin=res_lin, res_cub=res_cub,
            s_lin=s_lin, s_cub=s_cub,
        ))
        ratio = s_lin["rms"] / (h ** 2)  # should be roughly constant if O(h²)
        print(f"{N:>6} {h:>8.4f}  "
              f"{s_lin['rms']:>10.2e} {s_cub['rms']:>10.2e}  "
              f"{ratio:>11.3f} {s_lin['rms_centre']:>11.2e} "
              f"{s_lin['rms_edge']:>12.2e}  "
              f"{s_lin['n_finite']}/{s_lin['n_total']}")

    # Plot
    n = len(all_results)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 4 * rows),
                             squeeze=False)
    vmax = max(np.nanmax(np.abs(r["res_lin"])) for r in all_results)
    for ax, r in zip(axes.flat, all_results):
        sc = ax.scatter(test_pts[:, 0], test_pts[:, 1],
                        c=r["res_lin"], s=40, cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax)
        ax.scatter(r["ref_pts"][:, 0], r["ref_pts"][:, 1],
                   s=2, c="k", alpha=0.15)
        ax.set_aspect("equal")
        ax.set_title(f"N_ref={r['N']}  h={r['h']:.3f}\n"
                     f"RMS_lin={r['s_lin']['rms']:.2e}  "
                     f"max={r['s_lin']['max']:.2e}")
        ax.set_xlim(-args.R * 1.05, args.R * 1.05)
        ax.set_ylim(-args.R * 1.05, args.R * 1.05)
        plt.colorbar(sc, ax=ax, fraction=0.046)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle(f"LinearND residual on f=x²+y²  (test {args.n_test}², "
                 f"ref rotated {args.angle_deg}°)", y=1.02)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=120, bbox_inches="tight")
    print(f"\nFigure saved to {args.out_png}")


if __name__ == "__main__":
    main()
