"""
landscape_l2_ladder.py — replay-only L2 residual cost on the FSS ladder.

For each (L_test, L_ref) ladder cell α and each test point (r1, r2):

  J_α(r1, r2) = Σ_s ( G_test_α(r1, r2; x_s) − G_ref_α(x_s) )^2

where x_s are fixed sample points along the three boundary cycles
(t = k/N_samp for k=1..N_samp-1, all three cycles -> 3*(N_samp-1)
samples). Both G_test and G_ref are evaluated by tiled linear
interpolation of their respective all-to-all data.

Then continuum-extrapolate J_α(θ) → J_∞(θ) by a per-θ weighted
least-squares fit J(L) = J_∞ + a/L (3 ladder rungs; we use the
test-side L as the FSS variable since the ref scales with it).

Outputs:
  results/<tag>/l2_ladder.png      4-panel heatmap
  results/<tag>/l2_ladder.npz      grids + extrapolation arrays
  Also prints a per-point table when --points is set (sparse mode).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import _tile_interp, _SQRT3_2, boundary_paths


# ---------------------------------------------------------------------
# Sample points: t = k/N_samp on each cycle of the TEST torus, mapped
# to physical (x,y) coordinates. Same fractional parameterisation is
# used to sample the REFERENCE torus (so the comparison is at matched
# fractional positions along matched cycles, not matched absolute
# physical coordinates — the right convention for FSS where the box
# itself scales).
# ---------------------------------------------------------------------
def _sample_xy_on_torus(Lx, Ly, Tx, Ty, N_samp):
    """Return (x, y) arrays of shape (3*(N_samp-1),) in physical units
    on the given torus. Drops t=0 and t=1."""
    xs, ys = [], []
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    for (dm, dn) in paths:
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        for k in range(1, N_samp):
            t = k / N_samp
            xs.append(t * ex); ys.append(t * ey)
    return np.asarray(xs), np.asarray(ys)


def _sample_data(data, Lx, Ly, Tx, Ty, N_samp, copies=2):
    """Sample G_conn at the cycle parameterisation t=k/N_samp on this
    torus. Returns array of length 3*(N_samp-1)."""
    iref = _tile_interp(data, Lx, Ly, Tx, Ty, "conn", copies)
    xs, ys = _sample_xy_on_torus(Lx, Ly, Tx, Ty, N_samp)
    pts = np.column_stack([xs, ys])
    return np.asarray(iref(pts), dtype=float)


def _wls_intercept(L, J, sigma=None):
    """Weighted least squares fit J = a + b/L; return (intercept, sigma).
    If sigma is None, do unweighted OLS."""
    L = np.asarray(L, float); J = np.asarray(J, float)
    invL = 1.0 / L
    if sigma is None:
        w = np.ones_like(J)
    else:
        sigma = np.asarray(sigma, float)
        w = 1.0 / np.maximum(sigma**2, 1e-30)
    # design matrix [1, 1/L]
    X = np.column_stack([np.ones_like(L), invL])
    WX = X * w[:, None]
    A = X.T @ WX
    b = X.T @ (w * J)
    try:
        beta = np.linalg.solve(A, b)
        cov = np.linalg.inv(A)
        return float(beta[0]), float(np.sqrt(cov[0, 0]))
    except np.linalg.LinAlgError:
        return float(np.mean(J)), float(np.std(J))


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="_ladder_111")
    ap.add_argument("--test-sizes", type=int, nargs="+",
                    default=[8, 12, 16])
    ap.add_argument("--ref-sizes", type=int, nargs="+",
                    default=[16, 24, 32])
    ap.add_argument("--N-samp", type=int, default=8,
                    help="cycle samples per direction (drops t=0,1)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if len(args.test_sizes) != len(args.ref_sizes):
        print("ERROR: ladder needs len(test_sizes)==len(ref_sizes)",
              file=sys.stderr); return 2
    rungs = list(zip(args.test_sizes, args.ref_sizes))
    n_rungs = len(rungs)

    root = os.path.join(_HERE, "results", args.tag)
    test_root = os.path.join(root, "test")
    ref_root  = os.path.join(root, "ref")

    # Pre-load reference samples per ref size.
    ref_vec = {}
    for L_ref in args.ref_sizes:
        path = os.path.join(ref_root, f"L{L_ref}",
                            "two_point_all_to_all.dat")
        if not os.path.exists(path):
            print(f"ERROR: missing reference {path}", file=sys.stderr)
            return 1
        rd = mc_engine.load_all_to_all(path)
        ref_vec[L_ref] = _sample_data(rd, L_ref, L_ref, 0, 0, args.N_samp)
        print(f"[ref] L={L_ref}  samples={len(ref_vec[L_ref])}  "
              f"|G|_med={np.nanmedian(np.abs(ref_vec[L_ref])):.4g}")

    # Discover test points.
    grid_dir = os.path.join(test_root, "grid")
    files = sorted(os.listdir(grid_dir))
    pts = {}  # (r1, r2) -> {L: pkl_path}
    for fn in files:
        if not fn.endswith(".pkl"): continue
        # r1_X.XXX_r2_Y.YYY_LZZ.pkl
        try:
            base = fn[:-4]
            parts = base.split("_")
            r1 = float(parts[1]); r2 = float(parts[3])
            L  = int(parts[4][1:])
        except Exception:
            continue
        pts.setdefault((r1, r2), {})[L] = os.path.join(grid_dir, fn)

    pt_list = sorted(pts.keys())
    print(f"[test] discovered {len(pt_list)} test points across "
          f"sizes {args.test_sizes}")

    # Per-rung J_α.
    Js = {α: {} for α in range(n_rungs)}  # α -> {(r1,r2): J}
    for (r1, r2), files_for_pt in pts.items():
        for α, (L_test, L_ref) in enumerate(rungs):
            if L_test not in files_for_pt:
                continue
            with open(files_for_pt[L_test], "rb") as f:
                pt = pickle.load(f)
            G_t = _sample_data(pt["test_data"], L_test, L_test, 0, 0,
                               args.N_samp)
            G_r = ref_vec[L_ref]
            mask = np.isfinite(G_t) & np.isfinite(G_r)
            if not mask.any():
                continue
            res = G_t[mask] - G_r[mask]
            Js[α][(r1, r2)] = float(np.sum(res * res))

    # Continuum extrapolation per (r1, r2).
    L_arr = np.array(args.test_sizes, float)
    J_inf = {}
    J_inf_sigma = {}
    for pt_xy in pt_list:
        Js_pt = []
        for α in range(n_rungs):
            if pt_xy not in Js[α]:
                Js_pt = None; break
            Js_pt.append(Js[α][pt_xy])
        if Js_pt is None or len(Js_pt) < 2:
            continue
        Js_pt = np.array(Js_pt, float)
        if len(Js_pt) >= 2:
            inter, sig = _wls_intercept(L_arr[:len(Js_pt)], Js_pt)
            J_inf[pt_xy] = inter
            J_inf_sigma[pt_xy] = sig

    # Print per-point table.
    print("\n[ladder cost table]")
    hdr = "  ".join([f"J(L_t={Lt},L_r={Lr})" for Lt, Lr in rungs])
    print(f"  (r1, r2)         {hdr}    J_inf      σ_inf")
    for pt_xy in pt_list:
        row = []
        for α in range(n_rungs):
            v = Js[α].get(pt_xy, np.nan)
            row.append(f"{v:14.5g}" if np.isfinite(v) else "         nan ")
        ji = J_inf.get(pt_xy, np.nan)
        si = J_inf_sigma.get(pt_xy, np.nan)
        print(f"  ({pt_xy[0]:.3f},{pt_xy[1]:.3f})  " + "  ".join(row) +
              f"  {ji:10.4g} ±{si:8.3g}")

    # Plot: scatter for sparse points, heatmap for grid.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r1s = sorted({p[0] for p in pt_list})
    r2s = sorted({p[1] for p in pt_list})
    is_grid = (len(pt_list) == len(r1s) * len(r2s)) and len(pt_list) >= 9

    ncols = n_rungs + 1
    fig, axes = plt.subplots(1, ncols, figsize=(4.8 * ncols, 5.0))
    titles = ([f"Rung {α+1}: L_t={Lt}, L_r={Lr}"
               for α, (Lt, Lr) in enumerate(rungs)]
              + ["J_∞  (FSS extrap.)"])
    panels = [Js[α] for α in range(n_rungs)] + [J_inf]

    # Shared colour scale across all rungs (not J_∞ which can go negative).
    all_finite = [v for Z in [Js[α] for α in range(n_rungs)]
                  for v in Z.values() if np.isfinite(v) and v > 0]
    global_vmin = np.log10(np.percentile(all_finite, 2))  if all_finite else -3
    global_vmax = np.log10(np.percentile(all_finite, 98)) if all_finite else  0

    def _panel(ax, name, Z, shared_scale=True):
        xs = np.array([p[0] for p in pt_list], float)
        ys = np.array([p[1] for p in pt_list], float)
        zs = np.array([Z.get(p, np.nan) for p in pt_list], float)
        valid = np.isfinite(zs)
        if not valid.any():
            ax.axis("off"); return
        if is_grid:
            R1, R2 = np.meshgrid(r1s, r2s, indexing="ij")
            G = np.full(R1.shape, np.nan)
            for p, v in Z.items():
                ii = r1s.index(p[0]); jj = r2s.index(p[1])
                G[ii, jj] = v
            Zp = np.log10(np.maximum(G, 1e-12))
            if shared_scale:
                vmin, vmax = global_vmin, global_vmax
            else:
                vmin, vmax = np.nanpercentile(Zp, [2, 98])
            im = ax.pcolormesh(R1, R2, Zp, shading="auto",
                               cmap="viridis", vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                         label="log₁₀ J")
        else:
            zlog = np.log10(np.maximum(zs[valid], 1e-12))
            vmin = global_vmin if shared_scale else zlog.min()
            vmax = global_vmax if shared_scale else zlog.max()
            sc = ax.scatter(xs[valid], ys[valid], c=zlog,
                            vmin=vmin, vmax=vmax,
                            s=500, cmap="viridis",
                            edgecolor="k", linewidth=0.6)
            for x_, y_, z_ in zip(xs[valid], ys[valid], zs[valid]):
                ax.annotate(f"{z_:.2g}", (x_, y_), fontsize=8,
                            ha="center", va="center", color="white")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04,
                         label="log₁₀ J")
        # truth (1,1)
        ax.plot(1.0, 1.0, "rX", ms=13, mec="k", mew=1.2, zorder=5,
                label="truth (1,1)")
        # argmin
        finite_pts = [(p, Z[p]) for p in pt_list
                      if p in Z and np.isfinite(Z[p])]
        if finite_pts:
            am = min(finite_pts, key=lambda z: z[1])[0]
            d_t = math.hypot(am[0] - 1.0, am[1] - 1.0)
            ax.plot(am[0], am[1], "w*", ms=14, mec="k", mew=1, zorder=6,
                    label=f"argmin d={d_t:.2f}")
            ax.legend(fontsize=7, loc="upper right")
            title_str = f"{name}\nargmin=({am[0]:.2f},{am[1]:.2f})  d={d_t:.3f}"
        else:
            title_str = name
        ax.set_title(title_str, fontsize=10)
        ax.set_xlabel("r₁"); ax.set_ylabel("r₂"); ax.set_aspect("equal")
        ax.set_xlim(min(r1s) - 0.15, max(r1s) + 0.15)
        ax.set_ylim(min(r2s) - 0.15, max(r2s) + 0.15)

    for ax, name, Z, is_last in zip(axes, titles, panels,
                                    [False]*n_rungs + [True]):
        _panel(ax, name, Z, shared_scale=not is_last)

    fig.suptitle(
        f"FSS ladder  L2 cost  1-1-1 triangle\n"
        f"test L∈{args.test_sizes}  ref L∈{args.ref_sizes}  "
        f"truth=(1,1)  n_traj≈2k/pt",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out_png = args.out or os.path.join(root, "l2_ladder.png")
    fig.savefig(out_png, dpi=130)
    print(f"\n→ {out_png}")

    # NPZ dump
    np.savez(os.path.join(root, "l2_ladder.npz"),
             rungs=np.array(rungs),
             pt_list=np.array(pt_list),
             J_per_rung=np.array(
                 [[Js[α].get(p, np.nan) for p in pt_list]
                  for α in range(n_rungs)]),
             J_inf=np.array([J_inf.get(p, np.nan) for p in pt_list]),
             J_inf_sigma=np.array(
                 [J_inf_sigma.get(p, np.nan) for p in pt_list]))
    print(f"→ {os.path.join(root, 'l2_ladder.npz')}")


if __name__ == "__main__":
    main()
