"""
landscape_scale_aware.py — replay scale-aware cost variants on the
precompute grid. Diagnoses whether the small-geometry trough goes away
once absolute scale information is reintroduced.

Variants:
  A) huber_log alone  (baseline; bowl tilted to corner)
  B) huber_log + λ · |log(r1 r2 / R0^2)|   (volume regularizer)
  C) pure volume cost  (r1 r2 - R0^2)^2  — to show the regularizer alone
  D) β_c match cost   (β_c(r1,r2) - β_c_target)^2  — diagnostic using
     known truth β_c=0.0628; shows the achievable basin if we had that
     handle.

Saves results/_landscape/<tag>/scale_aware.png and .npz.
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

import mc_engine
import landscape_dozen as ld

_REF_GEOMS = {"_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3)}
_TRUTH_456 = (5.0652, 7.7429)
_BETA_C_TRUTH = 0.0628
_R0_TRUTH = float(np.sqrt(_TRUTH_456[0] * _TRUTH_456[1]))  # geom-mean of truth


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag", default="_reference_Lx13_Ly16_Tx-3_Ty3")
    ap.add_argument("--lambda-vol", type=float, default=0.005,
                    help="weight λ on |log(r1 r2 / R0²)| penalty")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    manifest = json.load(open(os.path.join(root, "manifest.json")))
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    rL, rL2, rT, rT2 = _REF_GEOMS[args.ref_tag]
    ref_data = mc_engine.load_all_to_all(
        os.path.join(_HERE, "results", args.ref_tag,
                     "two_point_all_to_all.dat"))
    ref_pack = ld._ref_samples(ref_data, rL, rL2, rT, rT2, Lx, Ly, Tx, Ty)

    print("[load] grid pts ...")
    t0 = time.time()
    samples = {}
    betac_grid = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl):
                continue
            pt = pickle.load(open(pkl, "rb"))
            samples[(float(r1), float(r2))] = ld._test_samples(
                pt["test_data"], ref_pack, Lx, Ly, Tx, Ty)
            betac_grid[i, j] = pt["beta_c"]
    print(f"[load] {len(samples)} pts in {time.time()-t0:.1f}s")

    # A: huber_log
    Z_huber = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            d = samples.get((float(r1), float(r2)))
            if d is None: continue
            Z_huber[i, j] = ld.k_huber_log(d)

    # Volume penalty (geom-mean form): pen = |log(r1 r2 / R0²)|
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    PEN = np.abs(np.log(np.maximum(R1 * R2, 1e-12) / (_R0_TRUTH ** 2)))

    # B: huber_log + λ · pen
    Z_combo = Z_huber + args.lambda_vol * PEN

    # C: pure volume cost = (r1 r2 - R0²)²
    Z_vol = (R1 * R2 - _R0_TRUTH ** 2) ** 2

    # D: β_c match
    Z_betac = (betac_grid - _BETA_C_TRUTH) ** 2

    panels = [
        ("A) huber_log (baseline)",        Z_huber,  True),
        (f"B) huber_log + {args.lambda_vol}·|log(r1·r2/R0²)|",
                                            Z_combo,  True),
        ("C) pure volume cost (r1·r2−R0²)²", Z_vol,   True),
        ("D) (β_c − 0.0628)²  diagnostic",   Z_betac, True),
    ]

    print("\n[summary]")
    for title, Z, _ in panels:
        if np.all(np.isnan(Z)):
            print(f"  {title:50s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        print(f"  {title:50s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:.2f}  cost(truth)={Z[truth_i,truth_j]:.3e}")

    # Render
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (title, Z, do_log) in zip(axes, panels):
        Zp = np.log10(np.maximum(Z, 1e-15)) if do_log else Z
        vmin, vmax = np.nanpercentile(Zp, [2, 98])
        im = ax.pcolormesh(R1, R2, Zp, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label="log10 cost" if do_log else "cost")
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=14, mec="k", mew=1.2,
                label="truth")
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        ax.plot(rs[ij[0]], rs[ij[1]], "w*", ms=14, mec="k", mew=1.2,
                label=f"argmin")
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_title(title, fontsize=11)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"Scale-aware cost diagnostics  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)  R0=√(r1·r2)≈{_R0_TRUTH:.2f}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_png = os.path.join(root, "scale_aware.png")
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"\n[done] figure → {out_png}")

    np.savez(os.path.join(root, "scale_aware.npz"),
             R1=R1, R2=R2, rs=rs,
             Z_huber=Z_huber, Z_combo=Z_combo, Z_vol=Z_vol, Z_betac=Z_betac,
             PEN=PEN, betac_grid=betac_grid,
             lambda_vol=args.lambda_vol, R0=_R0_TRUTH,
             beta_c_truth=_BETA_C_TRUTH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
