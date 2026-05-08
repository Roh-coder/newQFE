"""
landscape_fss_compare.py — side-by-side FSS comparison of the 4-5-6 landscape.

Shows shape_only and shape/spread at two test sizes (L=8, L=16) against the
same 39x48 reference, in a single 2x2 figure, demonstrating whether the
cost basin tilts toward truth as L increases.

Output: results/_landscape/fss_compare.png
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
import matplotlib.colors as mcolors

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)

_TRUTH_456 = (5.0652, 7.7429)
_REF_TAG   = "_reference_Lx39_Ly48_Tx-9_Ty9"
_REF_GEOM  = (39, 48, -9, 9)


# ---------- identical helpers from landscape_offset.py ----------

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
    rms_total = float(np.sqrt(np.mean(np.asarray(rms_per)**2)))
    spread_c  = float(np.std(cs))
    return rms_total, spread_c


# ---------- compute grids ----------

def _compute(landscape_tag, ref_data, rs):
    root = os.path.join(_HERE, "results", "_landscape", landscape_tag)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    n_traj = manifest["n_traj_prod"]

    rp = _ref_pack(ref_data, *_REF_GEOM, Lx, Ly, Tx, Ty)

    n = len(rs)
    SHAPE  = np.full((n, n), np.nan)
    SPREAD = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl):
                continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            sh, sp = _per_pt_obs(pt["test_data"], rp, Lx, Ly, Tx, Ty)
            SHAPE[i, j]  = sh
            SPREAD[i, j] = sp

    RATIO = SHAPE / np.maximum(SPREAD, 1e-9)
    return Lx, n_traj, SHAPE, SPREAD, RATIO


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ref_dat = mc_engine.load_all_to_all(
        os.path.join(_HERE, "results", _REF_TAG, "two_point_all_to_all.dat"))

    # Use the same r grid as the existing runs (both share 0.5..10 step 0.5)
    rs = np.arange(0.5, 10.0 + 1e-9, 0.5)
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")

    sizes = ["Lx8_Ly8_Tx0_Ty0", "Lx16_Ly16_Tx0_Ty0"]
    results = {}
    for tag in sizes:
        print(f"[compute] {tag} ...")
        Lx, n_traj, SH, SP, RA = _compute(tag, ref_dat, rs)
        results[tag] = dict(Lx=Lx, n_traj=n_traj,
                            shape=SH, spread=SP, ratio=RA)
        am_sh = np.unravel_index(np.nanargmin(SH), SH.shape)
        am_ra = np.unravel_index(np.nanargmin(RA), RA.shape)
        print(f"  shape_only  argmin=({rs[am_sh[0]]:.1f},{rs[am_sh[1]]:.1f})  "
              f"d={np.hypot(rs[am_sh[0]]-_TRUTH_456[0], rs[am_sh[1]]-_TRUTH_456[1]):.2f}")
        print(f"  shape/spread argmin=({rs[am_ra[0]]:.1f},{rs[am_ra[1]]:.1f})  "
              f"d={np.hypot(rs[am_ra[0]]-_TRUTH_456[0], rs[am_ra[1]]-_TRUTH_456[1]):.2f}")

    # ---- figure: 3 rows × 2 cols ----
    #   row 0: shape_only   (L=8 | L=16)
    #   row 1: shape/spread (L=8 | L=16)
    #   row 2: 1D slice through truth along r2=r1 diagonal
    fig = plt.figure(figsize=(11, 15))
    gs  = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.30)

    cost_rows = [("shape_only", "shape"), ("shape/spread", "ratio")]

    # Compute shared colour limits so L=8 and L=16 panels are directly comparable
    for row_idx, (cost_name, key) in enumerate(cost_rows):
        Zarr = [results[t][key] for t in sizes]
        # shared percentile clamp across both sizes
        allvals = np.concatenate([z[np.isfinite(z)] for z in Zarr])
        vmin = float(np.percentile(allvals, 2))
        vmax = float(np.percentile(allvals, 98))

        for col_idx, tag in enumerate(sizes):
            r    = results[tag]
            Z    = r[key]
            Lx   = r["Lx"]
            ntraj = r["n_traj"]
            ax   = fig.add_subplot(gs[row_idx, col_idx])

            im = ax.pcolormesh(R1, R2, Z, shading="auto", cmap="viridis",
                               vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            # Truth marker
            ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX",
                    ms=13, mec="k", mew=1.2, label="truth", zorder=5)

            # Argmin marker
            ij = np.unravel_index(np.nanargmin(Z), Z.shape)
            am = (rs[ij[0]], rs[ij[1]])
            d_t = float(np.hypot(am[0]-_TRUTH_456[0], am[1]-_TRUTH_456[1]))
            ax.plot(am[0], am[1], "w*", ms=14, mec="k", mew=1, zorder=6)

            ax.set_xlabel("r₁", fontsize=11)
            ax.set_ylabel("r₂", fontsize=11)
            ax.set_title(
                f"{cost_name}   L={Lx}  ({ntraj//1000}k traj)\n"
                f"argmin=({am[0]:.1f},{am[1]:.1f})  d={d_t:.2f}",
                fontsize=10)

    # row 2: 1D cross-section along the diagonal r1=r2 (where truth would lie
    # if it were isotropic) and along the r2~1.52*r1 truth-parallel ridge
    ax_1d_sh  = fig.add_subplot(gs[2, 0])
    ax_1d_ra  = fig.add_subplot(gs[2, 1])

    colours = {"Lx8_Ly8_Tx0_Ty0": "tab:blue", "Lx16_Ly16_Tx0_Ty0": "tab:orange"}
    linestyles = {"Lx8_Ly8_Tx0_Ty0": "-", "Lx16_Ly16_Tx0_Ty0": "--"}

    # Slice: fix r2 = truth_r2 and vary r1
    truth_r2_idx = int(np.argmin(np.abs(rs - _TRUTH_456[1])))

    for tag in sizes:
        r     = results[tag]
        Lx    = r["Lx"]
        ntraj = r["n_traj"]
        col   = colours[tag]
        ls    = linestyles[tag]
        lbl   = f"L={Lx} ({ntraj//1000}k)"

        sh_slice = r["shape"][  :, truth_r2_idx]
        ra_slice = r["ratio"][  :, truth_r2_idx]

        # normalise to [0,1] so the two sizes can be overlaid
        def _norm(v):
            v = v.copy()
            v -= np.nanmin(v)
            mx = np.nanmax(v)
            return v / mx if mx > 0 else v

        ax_1d_sh.plot(rs, _norm(sh_slice), color=col, ls=ls, lw=2, label=lbl)
        ax_1d_ra.plot(rs, _norm(ra_slice), color=col, ls=ls, lw=2, label=lbl)

    for ax, title in [(ax_1d_sh, "shape_only  (slice at r₂=r₂_truth, normalised)"),
                      (ax_1d_ra, "shape/spread (slice at r₂=r₂_truth, normalised)")]:
        ax.axvline(_TRUTH_456[0], color="red", ls=":", lw=1.5, label="truth r₁")
        ax.set_xlabel("r₁", fontsize=11)
        ax.set_ylabel("normalised cost", fontsize=11)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=9)
        ax.set_xlim(rs[0], rs[-1])

    fig.suptitle(
        "FSS comparison:  4-5-6 test lattice vs ref (39,48,−9,9)\n"
        "Does the basin approach truth as L increases?",
        fontsize=12, fontweight="bold", y=0.995)

    out_dir = os.path.join(_HERE, "results", "_landscape")
    out_png = args.out or os.path.join(out_dir, "fss_compare.png")
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"→ {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
