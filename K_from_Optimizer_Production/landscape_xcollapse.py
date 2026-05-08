"""
landscape_xcollapse.py — cross-geometry collapse cost.

Brower-Owen (Sec 5.2): at criticality, G(r) along each lattice axis,
when r is rescaled by the physical edge length ell_axis, lies on a
SINGLE UNIVERSAL CURVE — the c=1/2 CFT correlator. This is geometry-
blind in the short-distance regime where bulk power-law dominates.

Cost: take 3 axes from REF (fixed at true couplings) + 3 axes from
TEST (varying r1,r2). Rescale each curve to t = r/ell_axis. Measure
collapse residual on the short-distance window. Twisted reference is
fine — it just contributes its own (already correct) curves which
test must match in the universal regime.

Kernels:
  1. xcollapse_var      pointwise variance of log G across all 6 axes
  2. xcollapse_pairs    pairwise (log G_i - log G_j)^2
  3. xcollapse_test_vs_ref   only test-vs-ref pairings (not within-lattice)
  4. xcollapse_short    var restricted to short-distance window (t<0.3)
  5. xcollapse_powerlaw fit slope -2*Delta to all 6, cost = scatter of
                         fitted Delta values
  6. xcollapse_residual  fit ALL 6 to single power law, sum sq residual
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
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)  # noqa: E402

_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
    "_reference_Lx39_Ly48_Tx-9_Ty9": (39, 48, -9, 9),
    "_reference_Lx16_Ly16_Tx0_Ty0":  (16, 16, 0, 0),
}
_TRUTH_456 = (5.0652, 7.7429)


# --------------------- per-axis sample extraction --------------------- #
def _test_axes(test_data, Lx, Ly, Tx, Ty, copies=2):
    """Return list of 3 dicts: per cycle axis on test lattice."""
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    out = []
    for (dm, dn) in paths:
        ks, ms, ns = _direction_lattice_steps(dm, dn)
        N = len(ks)
        if N < 3: out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        p_len = float(np.hypot(ex, ey))
        Gt = np.full(N, np.nan); et = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       Lx, Ly, Tx, Ty, copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]; et[i] = abs(entry["conn_err"])
        r_phys = t_arr * p_len
        out.append(dict(r=r_phys, G=Gt, e=et, p_len=p_len, src="test"))
    return out


def _ref_axes(ref_data, Lx, Ly, Tx, Ty, n_pts=24, r_max_frac=0.5, copies=2):
    """Return list of 3 dicts: per cycle axis on REF lattice via interp."""
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    iref     = _tile_interp(ref_data, Lx, Ly, Tx, Ty, "conn",     copies)
    iref_err = _tile_interp(ref_data, Lx, Ly, Tx, Ty, "conn_err", copies)
    out = []
    for (dm, dn) in paths:
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        p_len = float(np.hypot(ex, ey))
        # sample evenly along the cycle in [0, r_max_frac * p_len]
        t_arr = np.linspace(1.0 / n_pts, r_max_frac, n_pts)
        pts = np.column_stack([t_arr * ex, t_arr * ey])
        G  = np.asarray(iref(pts), dtype=float)
        e  = np.abs(np.asarray(iref_err(pts), dtype=float))
        r_phys = t_arr * p_len
        out.append(dict(r=r_phys, G=G, e=e, p_len=p_len, src="ref"))
    return out


def _gather_rescaled(axes, eps=1e-12):
    """Per axis: t=r/p_len, finite & positive G, sorted."""
    out = []
    for ax in axes:
        if ax is None: out.append(None); continue
        r = ax["r"]; G = ax["G"]
        m = np.isfinite(G) & np.isfinite(r) & (G > eps)
        if m.sum() < 3: out.append(None); continue
        t = r[m] / ax["p_len"]
        Gv = G[m]
        order = np.argsort(t)
        out.append(dict(t=t[order], G=Gv[order], src=ax["src"]))
    return out


def _common_grid(rescaled, n_grid=20, lo=0.05, hi=0.35):
    """Interpolate each curve onto a common t-grid (log G vs t)."""
    tg = np.linspace(lo, hi, n_grid)
    cols = []; srcs = []
    for d in rescaled:
        if d is None: cols.append(None); srcs.append(None); continue
        if d["t"][-1] < hi or d["t"][0] > lo:
            cols.append(None); srcs.append(None); continue
        lg = np.interp(tg, d["t"], np.log(d["G"]))
        cols.append(lg); srcs.append(d["src"])
    valid = [(c, s) for c, s in zip(cols, srcs) if c is not None]
    if len(valid) < 3:
        return None, None, None
    M = np.vstack([c for c, s in valid])
    S = [s for c, s in valid]
    return tg, M, S


# --------------------- cost kernels ---------------------------------- #
def k_xcollapse_var(rescaled):
    tg, M, S = _common_grid(rescaled)
    if M is None: return np.nan
    return float(np.mean(np.var(M, axis=0)))


def k_xcollapse_pairs(rescaled):
    tg, M, S = _common_grid(rescaled)
    if M is None: return np.nan
    n = M.shape[0]
    s = 0.0; cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = M[i] - M[j]; s += float(np.mean(d * d)); cnt += 1
    return s / max(cnt, 1)


def k_xcollapse_test_vs_ref(rescaled):
    """Only test↔ref pairs; ignore within-test and within-ref."""
    tg, M, S = _common_grid(rescaled)
    if M is None: return np.nan
    s = 0.0; cnt = 0
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            if S[i] == S[j]: continue
            d = M[i] - M[j]; s += float(np.mean(d * d)); cnt += 1
    if cnt == 0: return np.nan
    return s / cnt


def k_xcollapse_short(rescaled):
    tg, M, S = _common_grid(rescaled, n_grid=10, lo=0.05, hi=0.20)
    if M is None: return np.nan
    return float(np.mean(np.var(M, axis=0)))


def k_xcollapse_powerlaw(rescaled):
    """Fit each axis to log G = -2*Delta * log t + const over short-t.
    Cost = variance of fitted Delta values across axes."""
    deltas = []
    for d in rescaled:
        if d is None: continue
        m = (d["t"] >= 0.05) & (d["t"] <= 0.30)
        if m.sum() < 3: continue
        lt = np.log(d["t"][m]); lg = np.log(d["G"][m])
        slope = np.polyfit(lt, lg, 1)[0]
        deltas.append(-0.5 * slope)
    if len(deltas) < 3: return np.nan
    return float(np.var(deltas))


def k_xcollapse_residual(rescaled):
    """Fit ALL axes simultaneously to a single power law in t (one slope,
    one intercept). Cost = sum sq residual / (n_axes * n_pts)."""
    tg, M, S = _common_grid(rescaled, n_grid=12, lo=0.05, hi=0.25)
    if M is None: return np.nan
    lt = np.log(tg)
    # design matrix: y = a*lt + b; pool all axes equally
    n_ax, n_pt = M.shape
    X = np.column_stack([np.tile(lt, n_ax), np.ones(n_ax * n_pt)])
    y = M.flatten()
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    return float(np.mean(res * res))


KERNELS = [
    ("xcollapse_var",       k_xcollapse_var),
    ("xcollapse_pairs",     k_xcollapse_pairs),
    ("xcollapse_test_ref",  k_xcollapse_test_vs_ref),
    ("xcollapse_short",     k_xcollapse_short),
    ("xcollapse_powerlaw",  k_xcollapse_powerlaw),
    ("xcollapse_residual",  k_xcollapse_residual),
]


# --------------------- driver ---------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag",   default="_reference_Lx39_Ly48_Tx-9_Ty9")
    ap.add_argument("--out",       default=None)
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    # extract reference axes once
    ref_geom = _REF_GEOMS[args.ref_tag]
    ref_dir  = os.path.join(_HERE, "results", args.ref_tag)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    rax = _ref_axes(ref_data, *ref_geom)
    print(f"[ref] axes p_len = "
          f"{[a['p_len'] if a else None for a in rax]}")

    print(f"[load+extract] {n*n} grid pts ...")
    t0 = time.time()
    samples = {}
    for r1 in rs:
        for r2 in rs:
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            tax = _test_axes(pt["test_data"], Lx, Ly, Tx, Ty)
            # combine ref + test axes
            samples[(float(r1), float(r2))] = _gather_rescaled(rax + tax)
    print(f"[load+extract] {len(samples)} pts in {time.time()-t0:.1f}s")

    results = {}
    print("[eval] kernels:")
    for name, fn in KERNELS:
        Z = np.full((n, n), np.nan)
        t0 = time.time()
        for i, r1 in enumerate(rs):
            for j, r2 in enumerate(rs):
                d = samples.get((float(r1), float(r2)))
                if d is None: continue
                try:
                    Z[i, j] = fn(d)
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
        ax.set_xlabel("r1"); ax.set_ylabel("r2"); ax.set_aspect("equal")
    for ax in axes[K:]: ax.axis("off")
    fig.suptitle(f"Cross-geometry collapse  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or os.path.join(root, "xcollapse.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n[done] figure → {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
