"""
landscape_halflife.py — cost kernels keyed on the half-decay length scale.

Each cycle has a correlator G(r) sampled at physical distances r=t*|p|.
We define r_half = first r where G(r) = G(r_min)/2 (linear interp in
log G vs r, per cycle). r_half is a physical length scale; comparing it
between test and reference is geometry/period invariant.

Six kernels:
  1. half_len_diff     — sum_dir (r_half_test - r_half_ref)^2
  2. half_len_logdiff  — sum_dir (log r_half_test - log r_half_ref)^2
  3. half_ratio        — match ratios r_half[u]:r_half[v]:r_half[w]
                          (pure anisotropy fingerprint, scale-free)
  4. window_l2_log     — l2-log of G restricted to r where G_ref/G_ref(rmin)
                          ∈ [0.3, 0.9] (the decay window only)
  5. window_window     — same but mask both sides (G_test AND G_ref in band)
  6. rescaled_l2_log   — rescale test r-axis so r_half_test → r_half_ref,
                          then l2-log on common g-grid (interp in log g)

Loads the precompute pkl grid (no extra MC); writes one combined PNG.
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


# ------------------------- sample extraction -------------------------- #
def _ref_samples(ref_data, ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                 test_Lx, test_Ly, test_Tx, test_Ty, copies=2):
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
        # physical radial distances r = t * |p_ref|
        p_ref_len = float(np.hypot(rex, rey))
        r_ref = t_arr * p_ref_len
        out.append(dict(t=t_arr, G_ref=G_ref, e_ref=e_ref,
                        ms=np.asarray(ms), ns=np.asarray(ns),
                        rex=rex, rey=rey, N=N,
                        r_ref=r_ref, p_ref_len=p_ref_len))
    return out


def _test_samples(test_data, ref_pack,
                  test_Lx, test_Ly, test_Tx, test_Ty, copies=2):
    """Look up test correlator at integer (m,n) sites; r_test = |t * p_test|.
    p_test is the cycle vector for the TEST geometry."""
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)
    dirs = []
    for d, (tdm, tdn) in zip(ref_pack, test_paths):
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
        # physical test r at each sample
        tex, tey = tdm + 0.5 * tdn, _SQRT3_2 * tdn
        p_test_len = float(np.hypot(tex, tey))
        r_test = d["t"] * p_test_len
        dirs.append(dict(G_test=Gt, e_test=et,
                         G_ref=d["G_ref"], e_ref=d["e_ref"],
                         t=d["t"],
                         r_test=r_test, r_ref=d["r_ref"],
                         p_test_len=p_test_len, p_ref_len=d["p_ref_len"]))
    return dirs


# ------------------------- helpers ----------------------------------- #
def _r_half(r, G):
    """Half-decay length: first r where G drops to G[0]/2.
    Linear interp in (r, log G). Returns NaN if not bracketed.
    `r` and `G` must already be finite, sorted by r ascending,
    and G must be > 0."""
    if len(G) < 2:
        return np.nan
    g0 = G[0]
    target = 0.5 * g0
    if g0 <= 0:
        return np.nan
    # find first index where G drops below target
    below = np.where(G <= target)[0]
    if len(below) == 0:
        return np.nan
    j = below[0]
    if j == 0:
        return r[0]
    # interp in log G vs r between j-1 and j
    r0, r1_ = r[j-1], r[j]
    lg0, lg1 = np.log(G[j-1]), np.log(G[j])
    lt = np.log(target)
    if lg1 == lg0:
        return 0.5 * (r0 + r1_)
    alpha = (lt - lg0) / (lg1 - lg0)
    return r0 + alpha * (r1_ - r0)


def _clean(d, key_r, key_G, eps=1e-12):
    """Return (r, G) sorted, finite, positive."""
    r = d[key_r]; G = d[key_G]
    m = np.isfinite(r) & np.isfinite(G) & (G > eps)
    r = r[m]; G = G[m]
    order = np.argsort(r)
    return r[order], G[order]


def _per_dir_mean(per_dir):
    arr = [c for c in per_dir if c is not None and np.isfinite(c)]
    return float(np.mean(arr)) if arr else np.nan


# ------------------------- cost kernels ------------------------------- #
def k_half_len_diff(dirs):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        rt, Gt = _clean(d, "r_test", "G_test")
        rr, Gr = _clean(d, "r_ref",  "G_ref")
        ht = _r_half(rt, Gt); hr = _r_half(rr, Gr)
        if not (np.isfinite(ht) and np.isfinite(hr)):
            pd.append(None); continue
        pd.append(float((ht - hr) ** 2))
    return _per_dir_mean(pd)


def k_half_len_logdiff(dirs):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        rt, Gt = _clean(d, "r_test", "G_test")
        rr, Gr = _clean(d, "r_ref",  "G_ref")
        ht = _r_half(rt, Gt); hr = _r_half(rr, Gr)
        if not (np.isfinite(ht) and np.isfinite(hr) and ht > 0 and hr > 0):
            pd.append(None); continue
        pd.append(float((np.log(ht) - np.log(hr)) ** 2))
    return _per_dir_mean(pd)


def k_half_ratio(dirs):
    """Compare anisotropy fingerprints r_half_u : r_half_v : r_half_w
    Normalised to mean=1, then sum sq diff. Scale-free."""
    ht = []; hr = []
    for d in dirs:
        if d is None: continue
        rt, Gt = _clean(d, "r_test", "G_test")
        rr, Gr = _clean(d, "r_ref",  "G_ref")
        a = _r_half(rt, Gt); b = _r_half(rr, Gr)
        if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0:
            ht.append(a); hr.append(b)
    if len(ht) < 2:
        return np.nan
    ht = np.asarray(ht); hr = np.asarray(hr)
    ht = ht / ht.mean()
    hr = hr / hr.mean()
    return float(np.mean((ht - hr) ** 2))


def k_window_l2_log(dirs, lo=0.3, hi=0.9, eps=1e-12):
    """l2-log of G restricted to samples where G_ref(r)/G_ref(r_min) ∈ [lo, hi].
    Window is anchored to ref's own decay — physical, no period dependence."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        Gr = d["G_ref"]; Gt = d["G_test"]
        m = (np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps))
        if m.sum() < 2: pd.append(None); continue
        # anchor: largest G_ref within mask (~ smallest r)
        anchor = float(np.max(Gr[m]))
        ratio = Gr / max(anchor, eps)
        m2 = m & (ratio >= lo) & (ratio <= hi)
        if m2.sum() < 2: pd.append(None); continue
        r = np.log(Gt[m2]) - np.log(Gr[m2])
        pd.append(float(np.mean(r * r)))
    return _per_dir_mean(pd)


def k_window_both(dirs, lo=0.3, hi=0.9, eps=1e-12):
    """Same window mask but applied to BOTH G_ref and G_test, intersected."""
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        Gr = d["G_ref"]; Gt = d["G_test"]
        m = (np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps))
        if m.sum() < 2: pd.append(None); continue
        ar = float(np.max(Gr[m])); at = float(np.max(Gt[m]))
        rr = Gr / max(ar, eps); rt = Gt / max(at, eps)
        m2 = m & (rr >= lo) & (rr <= hi) & (rt >= lo) & (rt <= hi)
        if m2.sum() < 2: pd.append(None); continue
        r = np.log(Gt[m2]) - np.log(Gr[m2])
        pd.append(float(np.mean(r * r)))
    return _per_dir_mean(pd)


def k_rescaled_l2_log(dirs, n_grid=12, lo=0.2, hi=0.95, eps=1e-12):
    """Rescale test r-axis so r_half_test → r_half_ref. Then resample both
    correlators on a common g-target grid (g = G/G[0]) by linear interp in
    (log G vs r). Compare the resulting r values: cost = sum_dir
    mean((r_test_at_g - r_ref_at_g)^2 / r_ref_at_g^2). Pure shape-after-scale.
    """
    pd = []
    targets = np.linspace(lo, hi, n_grid)
    for d in dirs:
        if d is None: pd.append(None); continue
        rt, Gt = _clean(d, "r_test", "G_test")
        rr, Gr = _clean(d, "r_ref",  "G_ref")
        if len(Gt) < 3 or len(Gr) < 3:
            pd.append(None); continue
        # normalise by first sample
        gt = Gt / Gt[0]; gr = Gr / Gr[0]
        # we want r(g) — invert. Need monotone decreasing g.
        # take cumulative running min so g is non-increasing.
        gt_mono = np.minimum.accumulate(gt)
        gr_mono = np.minimum.accumulate(gr)
        # r at target g via linear interp in log g
        # numpy.interp expects increasing x; flip.
        try:
            r_t_at = np.interp(np.log(targets[::-1]),
                               np.log(gt_mono[::-1]), rt[::-1])[::-1]
            r_r_at = np.interp(np.log(targets[::-1]),
                               np.log(gr_mono[::-1]), rr[::-1])[::-1]
        except Exception:
            pd.append(None); continue
        ok = np.isfinite(r_t_at) & np.isfinite(r_r_at) & (r_r_at > 0)
        if ok.sum() < 3:
            pd.append(None); continue
        # rescale test by half-life ratio
        ht = _r_half(rt, Gt); hr = _r_half(rr, Gr)
        if np.isfinite(ht) and np.isfinite(hr) and ht > 0:
            r_t_at = r_t_at * (hr / ht)
        diff = (r_t_at[ok] - r_r_at[ok]) / r_r_at[ok]
        pd.append(float(np.mean(diff * diff)))
    return _per_dir_mean(pd)


KERNELS = [
    ("half_len_diff",     k_half_len_diff),
    ("half_len_logdiff",  k_half_len_logdiff),
    ("half_ratio",        k_half_ratio),
    ("window_l2_log",     k_window_l2_log),
    ("window_both",       k_window_both),
    ("rescaled_l2_log",   k_rescaled_l2_log),
]


# ------------------------- driver ------------------------------------- #
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

    ref_geom = _REF_GEOMS[args.ref_tag]
    ref_dir  = os.path.join(_HERE, "results", args.ref_tag)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    rL, rL2, rT, rT2 = ref_geom
    ref_pack = _ref_samples(ref_data, rL, rL2, rT, rT2,
                            Lx, Ly, Tx, Ty)

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
            print(f"  {name:18s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        results[name] = dict(Z=Z, argmin=(rs[ij[0]], rs[ij[1]]),
                             d_truth=d_truth,
                             cost_truth=float(Z[truth_i, truth_j]),
                             cost_argmin=float(Z[ij]))
        print(f"  {name:18s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:5.2f}  cost(truth)={Z[truth_i,truth_j]:.3e}  "
              f"cost(min)={Z[ij]:.3e}  wall={time.time()-t0:.1f}s")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    K = len(results)
    cols = 3
    rows = (K + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6*cols, 3.8*rows))
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
    fig.suptitle(f"Half-decay-length cost kernels  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or os.path.join(root, "halflife.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n[done] figure → {out}")

    npz_out = (os.path.splitext(out)[0]) + ".npz"
    np.savez(npz_out, R1=R1, R2=R2, rs=rs,
             **{k: v["Z"] for k, v in results.items()})
    print(f"[done] data   → {npz_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
