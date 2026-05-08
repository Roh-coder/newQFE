"""
landscape_dozen.py — evaluate ~12 cost kernels against the precompute grid.

Loads results/_landscape/<tag>/grid/*.pkl, extracts per-direction
(G_test, σ_test, G_ref, σ_ref) arrays once per grid point (using the
same boundary_paths/_tile_interp/_lookup_test_value plumbing as
cost._l2_cost_test_native), then evaluates many cost kernels on those
arrays without any further MC.

Output: one combined PNG (4×3 grid of heatmaps) + per-kernel summary.
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
import cost as cost_module  # noqa: E402
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)  # noqa: E402

_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
    "_reference_Lx39_Ly48_Tx-9_Ty9": (39, 48, -9, 9),
    "_reference_Lx16_Ly16_Tx0_Ty0":  (16, 16, 0, 0),
    "_reference_Lx12_Ly12_Tx0_Ty0":  (12, 12, 0, 0),
    "_reference_Lx10_Ly10_Tx0_Ty0":  (10, 10, 0, 0),
    "_reference_Lx8_Ly8_Tx0_Ty0":    (8, 8, 0, 0),
}
_TRUTH_456 = (5.0652, 7.7429)


# ------------------------- sample extraction -------------------------- #
def _ref_samples(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                 test_Lx, test_Ly, test_Tx, test_Ty, copies=2):
    """Return list of 3 dicts with per-direction (G_ref, e_ref, t)."""
    iref     = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn", copies)
    iref_err = _tile_interp(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                            "conn_err", copies)
    ref_paths  = boundary_paths(ref_Lx,  ref_Ly,  ref_Tx,  ref_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    out = []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 2:
            out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts_ref  = np.column_stack([t_arr * rex, t_arr * rey])
        G_ref = np.asarray(iref(pts_ref), dtype=float)
        e_ref = np.abs(np.asarray(iref_err(pts_ref), dtype=float))
        out.append(dict(t=t_arr, G_ref=G_ref, e_ref=e_ref,
                        ms=np.asarray(ms), ns=np.asarray(ns),
                        rex=rex, rey=rey, N=N))
    return out


def _test_samples(test_data, ref_pack,
                  test_Lx, test_Ly, test_Tx, test_Ty, copies=2):
    """For each direction, look up test correlator at integer lattice sites."""
    dirs = []
    for d in ref_pack:
        if d is None:
            dirs.append(None); continue
        Gt = np.full(d["N"], np.nan)
        et = np.full(d["N"], np.nan)
        for i, (mk, nk) in enumerate(zip(d["ms"], d["ns"])):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
                et[i] = abs(entry["conn_err"])
        dirs.append(dict(G_test=Gt, e_test=et,
                         G_ref=d["G_ref"], e_ref=d["e_ref"],
                         t=d["t"]))
    return dirs


# ------------------------- cost kernels ------------------------------- #
def _mask(d):
    return (np.isfinite(d["G_test"]) & np.isfinite(d["G_ref"])
            & np.isfinite(d["e_test"]) & np.isfinite(d["e_ref"]))


def _per_dir_mean(per_dir):
    arr = [c for c in per_dir if c is not None and np.isfinite(c)]
    return float(np.mean(arr)) if arr else np.nan


def k_l2_diff(dirs):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 2: pd.append(None); continue
        r = d["G_test"][m] - d["G_ref"][m]
        pd.append(float(np.mean(r*r)))
    return _per_dir_mean(pd)


def k_l4_diff(dirs):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 2: pd.append(None); continue
        r = d["G_test"][m] - d["G_ref"][m]
        pd.append(float(np.mean(r**4)))
    return _per_dir_mean(pd)


def k_l1_diff(dirs):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 2: pd.append(None); continue
        r = d["G_test"][m] - d["G_ref"][m]
        pd.append(float(np.mean(np.abs(r))))
    return _per_dir_mean(pd)


def k_l2_log(dirs, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 2: pd.append(None); continue
        r = np.log(d["G_test"][m]) - np.log(d["G_ref"][m])
        pd.append(float(np.mean(r*r)))
    return _per_dir_mean(pd)


def k_l1_log(dirs, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 2: pd.append(None); continue
        r = np.log(d["G_test"][m]) - np.log(d["G_ref"][m])
        pd.append(float(np.mean(np.abs(r))))
    return _per_dir_mean(pd)


def k_huber_log(dirs, delta=0.5, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 2: pd.append(None); continue
        r = np.log(d["G_test"][m]) - np.log(d["G_ref"][m])
        a = np.abs(r)
        rho = np.where(a <= delta, 0.5*r*r, delta*(a - 0.5*delta))
        pd.append(float(np.mean(rho)))
    return _per_dir_mean(pd)


def k_chi2(dirs, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 2: pd.append(None); continue
        r = d["G_test"][m] - d["G_ref"][m]
        v = d["e_test"][m]**2 + d["e_ref"][m]**2 + eps
        pd.append(float(np.mean(r*r / v)))
    return _per_dir_mean(pd)


def k_chi2_log(dirs, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 2: pd.append(None); continue
        r  = np.log(d["G_test"][m]) - np.log(d["G_ref"][m])
        vt = (d["e_test"][m] / d["G_test"][m])**2
        vr = (d["e_ref"][m]  / d["G_ref"][m]) **2
        pd.append(float(np.mean(r*r / (vt + vr + eps))))
    return _per_dir_mean(pd)


def k_relative(dirs, eps=1e-12):
    """((Gt-Gr)/Gr)^2 — penalise relative misfit, ignore tiny values' weight."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (np.abs(d["G_ref"]) > eps)
        if m.sum() < 2: pd.append(None); continue
        r = (d["G_test"][m] - d["G_ref"][m]) / d["G_ref"][m]
        pd.append(float(np.mean(r*r)))
    return _per_dir_mean(pd)


def k_drop_short_log(dirs, drop=2, eps=1e-12):
    """l2 on log G but drop the first `drop` samples (short distances)."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        idx = np.where(m)[0]
        if len(idx) <= drop + 1: pd.append(None); continue
        idx = idx[drop:]
        r = np.log(d["G_test"][idx]) - np.log(d["G_ref"][idx])
        pd.append(float(np.mean(r*r)))
    return _per_dir_mean(pd)


def k_effmass(dirs, eps=1e-12):
    """squared diff of effective masses: m_k = -log(G[k+1]/G[k]).
    Short-distance saturation contributes only via the slope, so the
    constant offset (G≈1) doesn't dominate.
    """
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 3: pd.append(None); continue
        Gt = d["G_test"][m]; Gr = d["G_ref"][m]
        mt = -np.diff(np.log(Gt))
        mr = -np.diff(np.log(Gr))
        pd.append(float(np.mean((mt - mr)**2)))
    return _per_dir_mean(pd)


def k_slope_loglog(dirs, eps=1e-12):
    """Diff of fitted slope (log G vs t) per direction."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d) & (d["G_test"] > eps) & (d["G_ref"] > eps)
        if m.sum() < 3: pd.append(None); continue
        t = d["t"][m]
        st = np.polyfit(t, np.log(d["G_test"][m]), 1)[0]
        sr = np.polyfit(t, np.log(d["G_ref"][m]),  1)[0]
        pd.append(float((st - sr)**2))
    return _per_dir_mean(pd)


def k_affine(dirs, eps=1e-12):
    """Best a,b minimising ||a*Gt + b - Gr||^2 — residual / N."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 3: pd.append(None); continue
        Gt = d["G_test"][m]; Gr = d["G_ref"][m]
        A = np.column_stack([Gt, np.ones_like(Gt)])
        coef, *_ = np.linalg.lstsq(A, Gr, rcond=None)
        res = Gr - A @ coef
        pd.append(float(np.mean(res*res)))
    return _per_dir_mean(pd)


def k_cosine(dirs):
    """1 - cosine similarity of G_test and G_ref per direction."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = _mask(d)
        if m.sum() < 2: pd.append(None); continue
        a = d["G_test"][m]; b = d["G_ref"][m]
        denom = np.linalg.norm(a) * np.linalg.norm(b)
        if denom <= 0: pd.append(None); continue
        pd.append(float(1.0 - np.dot(a, b) / denom))
    return _per_dir_mean(pd)


KERNELS = [
    ("l2_diff",         k_l2_diff),
    ("l4_diff",         k_l4_diff),
    ("l1_diff",         k_l1_diff),
    ("l2_log",          k_l2_log),
    ("l1_log",          k_l1_log),
    ("huber_log",       k_huber_log),
    ("chi2",            k_chi2),
    ("chi2_log",        k_chi2_log),
    ("relative",        k_relative),
    ("drop_short_log",  k_drop_short_log),
    ("effmass",         k_effmass),
    ("slope_loglog",    k_slope_loglog),
    ("affine",          k_affine),
    ("cosine",          k_cosine),
]


# ------------------------- driver ------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag",   default="_reference_Lx13_Ly16_Tx-3_Ty3")
    ap.add_argument("--out",       default=None,
                    help="output PNG (default: results/_landscape/<tag>/dozen.png)")
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
    ref_pack = _ref_samples(ref_data, rL, rL2, rT, rT2,
                            Lx, Ly, Tx, Ty)

    # Pre-load all test pkls and extract sample arrays.
    print(f"[load] {n*n} grid pts ...")
    t0 = time.time()
    samples = {}
    for r1 in rs:
        for r2 in rs:
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            samples[(float(r1), float(r2))] = _test_samples(
                pt["test_data"], ref_pack, Lx, Ly, Tx, Ty)
    print(f"[load] {len(samples)} pts in {time.time()-t0:.1f}s")

    # Evaluate every kernel.
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
        # argmin + truth cell
        if np.all(np.isnan(Z)):
            print(f"  {name:18s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        # distance from argmin to truth cell (in r-units)
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        results[name] = dict(Z=Z, argmin=(rs[ij[0]], rs[ij[1]]),
                             d_truth=d_truth,
                             cost_truth=float(Z[truth_i, truth_j]),
                             cost_argmin=float(Z[ij]))
        print(f"  {name:18s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:5.2f}  cost(truth)={Z[truth_i,truth_j]:.3e}  "
              f"cost(min)={Z[ij]:.3e}  wall={time.time()-t0:.1f}s")

    # Render combined heatmap figure.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    K = len(results)
    cols = 4
    rows = (K + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.2*cols, 3.6*rows))
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
        ax.set_title(f"{name}\nargmin=({am[0]:.1f},{am[1]:.1f}) d={info['d_truth']:.1f}",
                     fontsize=10)
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_aspect("equal")
    for ax in axes[K:]: ax.axis("off")
    fig.suptitle(f"Cost-kernel landscapes  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = args.out or os.path.join(root, "dozen.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n[done] figure → {out}")

    # Save .npz of all kernels
    npz_out = os.path.join(root, "dozen.npz")
    np.savez(npz_out, R1=R1, R2=R2, rs=rs,
             **{k: v["Z"] for k, v in results.items()})
    print(f"[done] data   → {npz_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
