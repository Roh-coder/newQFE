"""
landscape_offset.py — landscape of "fit per-direction amplitude c, look
at residual-after-offset and the spread of c across directions".

Diagnostic plot evidence: c=[c_u, c_v, c_w] is near-uniform at small/
isotropic couplings (spread~0.03–0.09) and clearly anisotropic near
truth (spread~0.20). Meanwhile the shape-residual after offset is
~uniformly small everywhere (0.06–0.10). So spread(c) carries the
anisotropy fingerprint.

Kernels:
  shape_only      mean rms residual after per-dir multiplicative offset
  spread_c        std(c) across the 3 directions  (LARGER near truth)
  spread_c_minus  -spread_c (so that argmin matches truth)
  shape_over_spread  shape_only / spread_c   (small near truth)
  shape_over_spread2 shape_only^2 / spread_c
  shape_minus_spread  shape_only - lambda*spread_c
  spread_pure     same as spread_c (computed separately)

Cost is over a grid of (r1,r2) precompute pkls.
"""
from __future__ import annotations

import argparse, json, os, pickle, sys, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)

_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
    "_reference_Lx39_Ly48_Tx-9_Ty9": (39, 48, -9, 9),
    "_reference_Lx16_Ly16_Tx0_Ty0":  (16, 16, 0, 0),
}
_TRUTH_456 = (5.0652, 7.7429)


def _ref_pack(ref_data, rL, rL2, rT, rT2, tL, tL2, tT, tT2, copies=2):
    iref = _tile_interp(ref_data, rL, rL2, rT, rT2, "conn", copies)
    rp = boundary_paths(rL, rL2, rT, rT2)
    tp = boundary_paths(tL, tL2, tT, tT2)
    out = []
    for (rdm, rdn), (tdm, tdn) in zip(rp, tp):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 3: out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5*rdn, _SQRT3_2*rdn
        pts_ref = np.column_stack([t_arr*rex, t_arr*rey])
        Gr = np.asarray(iref(pts_ref), dtype=float)
        out.append(dict(t=t_arr, Gr=Gr, ms=ms, ns=ns, N=N))
    return out


def _per_pt_obs(test_data, ref_pack, tL, tL2, tT, tT2, copies=2,
                drop_first=1, drop_last=1, eps=1e-12):
    """Return (rms_total, spread_c, c_arr) for this (r1,r2) point.
    drop_first: skip the t≈0 sample(s) which carry the contact-term
                singularity.
    drop_last:  skip the t≈1 sample(s) which sit next to the
                wraparound (period-boundary) singularity."""
    cs = []; rms_per = []
    for d in ref_pack:
        if d is None: continue
        Gr = d["Gr"]; ms = d["ms"]; ns = d["ns"]
        N = d["N"]
        Gt = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       tL, tL2, tT, tT2, copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
        m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps)
        # drop the singular endpoints t≈0 and t≈1
        idx = np.where(m)[0]
        if len(idx) <= drop_first + drop_last + 1: continue
        if drop_last > 0:
            idx = idx[drop_first:-drop_last]
        else:
            idx = idx[drop_first:]
        log_diff = np.log(Gt[idx]) - np.log(Gr[idx])
        log_c = float(np.mean(log_diff))
        cs.append(np.exp(log_c))
        res = log_diff - log_c
        rms_per.append(float(np.sqrt(np.mean(res * res))))
    if len(cs) < 2:
        return np.nan, np.nan, []
    cs = np.asarray(cs)
    rms_total = float(np.sqrt(np.mean(np.asarray(rms_per) ** 2)))
    spread_c = float(np.std(cs))
    return rms_total, spread_c, cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag",   default="_reference_Lx39_Ly48_Tx-9_Ty9")
    ap.add_argument("--out",       default=None)
    ap.add_argument("--lam",       type=float, default=0.5,
                    help="lambda for shape − lam·spread combo")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    rg = _REF_GEOMS[args.ref_tag]
    ref_dir = os.path.join(_HERE, "results", args.ref_tag)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    ref_pack = _ref_pack(ref_data, *rg, Lx, Ly, Tx, Ty)

    print(f"[load+eval] {n*n} grid pts ...")
    t0 = time.time()
    SHAPE  = np.full((n, n), np.nan)
    SPREAD = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            sh, sp, _ = _per_pt_obs(pt["test_data"], ref_pack,
                                    Lx, Ly, Tx, Ty)
            SHAPE[i, j] = sh
            SPREAD[i, j] = sp
    print(f"[load+eval] done in {time.time()-t0:.1f}s")

    # cost panels
    panels = {
        "shape_only":            SHAPE,
        "spread_c":              -SPREAD,        # negate so argmin near truth
        "shape/spread":          SHAPE / np.maximum(SPREAD, 1e-9),
        "shape^2/spread":        SHAPE**2 / np.maximum(SPREAD, 1e-9),
        f"shape - {args.lam}*spread": SHAPE - args.lam * SPREAD,
        "shape/spread^2":        SHAPE / np.maximum(SPREAD**2, 1e-9),
    }

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    K = len(panels)
    cols = 3; rows = (K + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6*cols, 3.8*rows))
    axes = axes.flatten()
    print("\n[summary]")
    for ax, (name, Z) in zip(axes, panels.items()):
        if np.all(np.isnan(Z)):
            ax.axis("off"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        am = (rs[ij[0]], rs[ij[1]])
        d_t = float(np.hypot(am[0]-_TRUTH_456[0], am[1]-_TRUTH_456[1]))
        # render (use raw Z, not log, since some are negative)
        Zp = Z.copy()
        vmin, vmax = np.nanpercentile(Zp, [2, 98])
        im = ax.pcolormesh(R1, R2, Zp, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=11, mec="k", mew=1)
        ax.plot(am[0], am[1], "w*", ms=12, mec="k", mew=1)
        ax.set_title(f"{name}\nargmin=({am[0]:.1f},{am[1]:.1f}) d={d_t:.1f}",
                     fontsize=10)
        ax.set_xlabel("r1"); ax.set_ylabel("r2"); ax.set_aspect("equal")
        print(f"  {name:30s} argmin=({am[0]:5.2f},{am[1]:5.2f})  d_truth={d_t:.2f}")
    for ax in axes[K:]: ax.axis("off")
    fig.suptitle(f"Offset-removed cost  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or os.path.join(root, "offset_cost.png")
    fig.savefig(out, dpi=130)
    print("→", out)


if __name__ == "__main__":
    main()
