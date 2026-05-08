"""
landscape_zoom_bowl.py — focused heatmap near 4-5-6 truth, comparing
coarse-5k vs fine-20k grids for the shape/spread cost.

Produces a 3-panel figure:
  Left:   coarse grid L=16 5k  (existing Lx16_Ly16_Tx0_Ty0, step=0.5)
  Middle: fine grid   L=16 20k (new Lx16_Ly16_Tx0_Ty0_zoom25, step=0.25)
          — both zoomed to the truth region [r1_lo,r1_hi]×[r2_lo,r2_hi]
  Right:  1D horizontal slices at r2≈truth through both grids

Output: results/_landscape/zoom_bowl.png
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)

_TRUTH_456 = (5.0652, 7.7429)
_REF_TAG   = "_reference_Lx39_Ly48_Tx-9_Ty9"
_REF_GEOM  = (39, 48, -9, 9)

# Zoom window for display
_R1_LO, _R1_HI = 2.5, 8.5
_R2_LO, _R2_HI = 2.5, 8.5


# ---- identical helpers as landscape_offset.py ----

def _ref_pack(ref_data, rL, rL2, rT, rT2, tL, tL2, tT, tT2, copies=2):
    iref = _tile_interp(ref_data, rL, rL2, rT, rT2, "conn", copies)
    rp = boundary_paths(rL, rL2, rT, rT2)
    tp = boundary_paths(tL, tL2, tT, tT2)
    out = []
    for (rdm, rdn), (tdm, tdn) in zip(rp, tp):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 3:
            out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5*rdn, _SQRT3_2*rdn
        pts_ref = np.column_stack([t_arr*rex, t_arr*rey])
        Gr = np.asarray(iref(pts_ref), dtype=float)
        out.append(dict(t=t_arr, Gr=Gr, ms=ms, ns=ns, N=N))
    return out


def _per_pt_obs(test_data, ref_pack, tL, tL2, tT, tT2,
                copies=2, drop_first=1, drop_last=1, eps=1e-12):
    cs = []; rms_per = []
    for d in ref_pack:
        if d is None: continue
        Gr = d["Gr"]; ms = d["ms"]; ns = d["ns"]; N = d["N"]
        Gt = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       tL, tL2, tT, tT2, copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
        m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps)
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
        return np.nan, np.nan
    cs = np.asarray(cs)
    return float(np.sqrt(np.mean(np.asarray(rms_per)**2))), float(np.std(cs))


def _load_grid(landscape_tag, ref_data):
    root = os.path.join(_HERE, "results", "_landscape", landscape_tag)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    n_traj = manifest["n_traj_prod"]
    r_step = manifest["r_step"]

    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9, r_step)
    rp = _ref_pack(ref_data, *_REF_GEOM, Lx, Ly, Tx, Ty)

    n = len(rs)
    SHAPE  = np.full((n, n), np.nan)
    SPREAD = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid", f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl):
                continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            sh, sp = _per_pt_obs(pt["test_data"], rp, Lx, Ly, Tx, Ty)
            SHAPE[i, j]  = sh
            SPREAD[i, j] = sp

    RATIO = SHAPE / np.maximum(SPREAD, 1e-9)
    return rs, Lx, n_traj, r_step, SHAPE, SPREAD, RATIO


def _argmin_str(Z, rs):
    ij = np.unravel_index(np.nanargmin(Z), Z.shape)
    am = (rs[ij[0]], rs[ij[1]])
    d = float(np.hypot(am[0]-_TRUTH_456[0], am[1]-_TRUTH_456[1]))
    return am, d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coarse-tag", default="Lx16_Ly16_Tx0_Ty0",
                    help="existing coarse-grid landscape tag")
    ap.add_argument("--fine-tag",   default="Lx16_Ly16_Tx0_Ty0_zoom25",
                    help="new fine-grid landscape tag")
    ap.add_argument("--cost", choices=["shape_only", "shape/spread"],
                    default="shape/spread")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ref_dat = mc_engine.load_all_to_all(
        os.path.join(_HERE, "results", _REF_TAG, "two_point_all_to_all.dat"))

    datasets = {}
    for tag, label in [(args.coarse_tag, "coarse"), (args.fine_tag, "fine")]:
        mf = os.path.join(_HERE, "results", "_landscape", tag, "manifest.json")
        if not os.path.exists(mf):
            print(f"[warn] {tag} not found — skipping")
            continue
        print(f"[load] {tag} ...")
        rs, Lx, n_traj, r_step, SH, SP, RA = _load_grid(tag, ref_dat)
        datasets[label] = dict(rs=rs, Lx=Lx, n_traj=n_traj,
                               r_step=r_step, shape=SH, spread=SP, ratio=RA)
        key = "ratio" if args.cost == "shape/spread" else "shape"
        am, d = _argmin_str(datasets[label][key], rs)
        print(f"  {args.cost}  argmin=({am[0]:.2f},{am[1]:.2f})  d={d:.2f}")

    if not datasets:
        print("ERROR: no datasets found"); return 1

    cost_key = "ratio" if args.cost == "shape/spread" else "shape"

    # ---- figure ----
    n_panels = len(datasets)
    fig, axes = plt.subplots(1, n_panels + 1,
                             figsize=(5.5*(n_panels+1), 5.5),
                             gridspec_kw={"width_ratios": [1]*n_panels + [1.1]})

    # shared colour scale across both grids (interpolated to the same window)
    all_z = []
    for d in datasets.values():
        Z = d[cost_key]
        rs = d["rs"]
        mask_r1 = (rs >= _R1_LO) & (rs <= _R1_HI)
        mask_r2 = (rs >= _R2_LO) & (rs <= _R2_HI)
        sub = Z[np.ix_(mask_r1, mask_r2)]
        all_z.append(sub[np.isfinite(sub)])
    if all_z:
        flat = np.concatenate(all_z)
        vmin = float(np.percentile(flat, 2))
        vmax = float(np.percentile(flat, 98))
    else:
        vmin, vmax = 0, 1

    panel_titles = {"coarse": "coarse (step=0.5)", "fine": "fine (step=0.25)"}

    ax_1d = axes[-1]
    colours = {"coarse": "tab:blue", "fine": "tab:orange"}
    linestyles = {"coarse": "-", "fine": "--"}

    for ax_idx, (label, data) in enumerate(datasets.items()):
        ax = axes[ax_idx]
        rs  = data["rs"]
        Z   = data[cost_key]
        Lx  = data["Lx"]
        ntraj = data["n_traj"]
        step  = data["r_step"]

        # zoom mask
        mask_r1 = (rs >= _R1_LO) & (rs <= _R1_HI)
        mask_r2 = (rs >= _R2_LO) & (rs <= _R2_HI)
        rs_r1 = rs[mask_r1]
        rs_r2 = rs[mask_r2]
        Z_zoom = Z[np.ix_(mask_r1, mask_r2)]

        R1z, R2z = np.meshgrid(rs_r1, rs_r2, indexing="ij")
        im = ax.pcolormesh(R1z, R2z, Z_zoom, shading="auto",
                           cmap="viridis", vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=13, mec="k",
                mew=1.2, zorder=5, label="truth")
        am, d = _argmin_str(Z_zoom, rs_r1)  # argmin within zoom window
        ax.plot(am[0], am[1], "w*", ms=14, mec="k", mew=1, zorder=6)

        ax.set_xlim(_R1_LO, _R1_HI)
        ax.set_ylim(_R2_LO, _R2_HI)
        ax.set_xlabel("r₁", fontsize=12)
        ax.set_ylabel("r₂", fontsize=12)
        ax.set_title(
            f"{args.cost}   L={Lx}  step={step}  ({ntraj//1000}k traj)\n"
            f"argmin=({am[0]:.2f},{am[1]:.2f})  d={d:.2f}",
            fontsize=10)

        # 1D slice at r2 closest to truth
        truth_r2_idx = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        col = colours.get(label, "grey")
        ls  = linestyles.get(label, "-")
        sl  = Z[:, truth_r2_idx]
        # normalise
        sl_n = sl - np.nanmin(sl)
        mx = np.nanmax(sl_n)
        if mx > 0:
            sl_n /= mx
        ax_1d.plot(rs, sl_n, color=col, ls=ls, lw=2,
                   label=f"{label} L={Lx} {ntraj//1000}k")

    ax_1d.axvline(_TRUTH_456[0], color="red", ls=":", lw=1.5, label="truth r₁")
    ax_1d.set_xlabel("r₁", fontsize=12)
    ax_1d.set_ylabel("normalised cost", fontsize=12)
    ax_1d.set_title(f"{args.cost} — 1D slice at r₂≈{_TRUTH_456[1]:.1f}", fontsize=10)
    ax_1d.set_xlim(_R1_LO, _R1_HI)
    ax_1d.legend(fontsize=9)

    fig.suptitle(
        "4-5-6  shape/spread cost:  coarse-5k  vs  fine-20k  (L=16, ref=39×48)\n"
        f"Truth at ({_TRUTH_456[0]:.2f}, {_TRUTH_456[1]:.2f})",
        fontsize=12, fontweight="bold")

    out = args.out or os.path.join(
        _HERE, "results", "_landscape", "zoom_bowl.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"→ {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
