"""
landscape_cost.py — replay any cost mode against a precomputed grid.

Loads a precompute_landscape.py output for a chosen test_geom, plus a
reference correlator cache, and computes cost(r1, r2) for every grid
point. Saves an .npz and renders a heatmap PNG.

Examples
--------
python landscape_cost.py --landscape Lx16_Ly16_Tx0_Ty0 \
    --ref-tag _reference_Lx13_Ly16_Tx-3_Ty3 --cost-mode huber_log
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

# Reference geometries that ship with this repo (tag -> (Lx,Ly,Tx,Ty)).
_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
}

# Truth for the canonical 4-5-6 problem (only used as a marker).
_TRUTH_456 = (5.0652, 7.7429)


def _load_landscape(tag: str):
    root = os.path.join(_HERE, "results", "_landscape", tag)
    mpath = os.path.join(root, "manifest.json")
    if not os.path.exists(mpath):
        raise FileNotFoundError(f"manifest missing: {mpath}")
    with open(mpath) as f:
        manifest = json.load(f)
    return root, manifest


def _load_ref(ref_tag: str):
    if ref_tag in _REF_GEOMS:
        ref_geom = _REF_GEOMS[ref_tag]
    else:
        raise ValueError(f"unknown ref tag '{ref_tag}'; known: {list(_REF_GEOMS)}")
    ref_dir = os.path.join(_HERE, "results", ref_tag)
    a2a = os.path.join(ref_dir, "two_point_all_to_all.dat")
    if not os.path.exists(a2a):
        raise FileNotFoundError(f"reference cache missing: {a2a}")
    return ref_geom, mc_engine.load_all_to_all(a2a)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", required=True,
                    help="tag inside results/_landscape/")
    ap.add_argument("--ref-tag", default="_reference_Lx13_Ly16_Tx-3_Ty3")
    ap.add_argument("--cost-mode", default="huber_log",
                    choices=["l4mean_both_interp", "test_native",
                             "spectral", "affine_fit", "huber_log"])
    ap.add_argument("--cost-power", type=int, default=2)
    ap.add_argument("--vmax", type=float, default=None,
                    help="cap heatmap colour at this cost value (default: 99th pct)")
    ap.add_argument("--log", action="store_true",
                    help="render heatmap on log10 scale")
    args = ap.parse_args()

    root, manifest = _load_landscape(args.landscape)
    ref_geom, ref_data = _load_ref(args.ref_tag)
    Lx, Ly, Tx, Ty = manifest["geom"]
    ref_Lx, ref_Ly, ref_Tx, ref_Ty = ref_geom

    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    COST = np.full((n, n), np.nan)
    SIGMA = np.full((n, n), np.nan)
    BETAC = np.full((n, n), np.nan)

    t0 = time.time()
    n_loaded = n_missing = 0
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl):
                n_missing += 1
                continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            try:
                c, s, _pd, _pds = cost_module.l2_cost(
                    ref_data, pt["test_data"],
                    Lx, Ly, Tx, Ty,
                    ref_Lx, ref_Ly, ref_Tx, ref_Ty,
                    cost_mode=args.cost_mode,
                    cost_power=args.cost_power,
                )
                COST[i, j] = c
                SIGMA[i, j] = s
                BETAC[i, j] = pt["beta_c"]
                n_loaded += 1
            except Exception as e:  # noqa: BLE001
                print(f"  [skip] r=({r1:.2f},{r2:.2f}): {e}")

    wall = time.time() - t0
    print(f"[replay] {n_loaded}/{n*n} pts evaluated ({n_missing} missing)  "
          f"cost={args.cost_mode}  wall={wall:.1f}s")

    out_tag = f"{args.cost_mode}_p{args.cost_power}_ref-{args.ref_tag}"
    npz = os.path.join(root, f"cost_{out_tag}.npz")
    np.savez(npz, R1=R1, R2=R2, COST=COST, SIGMA=SIGMA, BETAC=BETAC,
             cost_mode=args.cost_mode, cost_power=args.cost_power,
             ref_tag=args.ref_tag)

    # ---- heatmap ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Cost heatmap.
    Z = COST.copy()
    valid = np.isfinite(Z)
    if not valid.any():
        print("ERROR: no finite cost values"); return 1

    if args.log:
        Z = np.log10(np.maximum(Z, 1e-12))
        cbar_label = f"log10 cost ({args.cost_mode})"
    else:
        cbar_label = f"cost ({args.cost_mode})"

    vmax = args.vmax if args.vmax is not None else np.nanpercentile(Z, 99)
    vmin = np.nanmin(Z)
    ax = axes[0]
    im = ax.pcolormesh(R1, R2, Z, shading="auto",
                       cmap="viridis", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label=cbar_label)
    # Mark argmin and truth.
    flat_min = np.nanargmin(COST)
    i_m, j_m = np.unravel_index(flat_min, COST.shape)
    ax.plot(R1[i_m, j_m], R2[i_m, j_m], "w*", ms=18, mec="k", mew=1.2,
            label=f"argmin ({R1[i_m,j_m]:.2f},{R2[i_m,j_m]:.2f})")
    ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=14, mec="k", mew=1.2,
            label=f"4-5-6 truth ({_TRUTH_456[0]:.2f},{_TRUTH_456[1]:.2f})")
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title(f"cost landscape — {args.cost_mode}  "
                 f"test=({Lx},{Ly},{Tx},{Ty})  ref={args.ref_tag}")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_aspect("equal")

    # β_c heatmap (sanity-check the MC).
    ax2 = axes[1]
    im2 = ax2.pcolormesh(R1, R2, BETAC, shading="auto", cmap="magma")
    fig.colorbar(im2, ax=ax2, label="β_c")
    ax2.plot(_TRUTH_456[0], _TRUTH_456[1], "wX", ms=14, mec="k", mew=1.2)
    ax2.set_xlabel("r1"); ax2.set_ylabel("r2")
    ax2.set_title("β_c(r1, r2) on the precompute grid")
    ax2.set_aspect("equal")

    fig.tight_layout()
    png = os.path.join(root, f"cost_{out_tag}.png")
    fig.savefig(png, dpi=140)
    plt.close(fig)
    print(f"[replay] npz → {npz}")
    print(f"[replay] png → {png}")

    # Quick numeric summary.
    print("\n=== summary ===")
    print(f"argmin(cost) = ({R1[i_m,j_m]:.3f}, {R2[i_m,j_m]:.3f}) "
          f"cost={COST[i_m,j_m]:.4e}  β_c={BETAC[i_m,j_m]:.5f}")
    print(f"4-5-6 truth = ({_TRUTH_456[0]:.3f}, {_TRUTH_456[1]:.3f})")
    # Cost at the cell nearest truth.
    di = np.argmin(np.abs(rs - _TRUTH_456[0]))
    dj = np.argmin(np.abs(rs - _TRUTH_456[1]))
    if np.isfinite(COST[di, dj]):
        print(f"cost at nearest-truth cell ({rs[di]:.2f},{rs[dj]:.2f}) "
              f"= {COST[di, dj]:.4e}  β_c={BETAC[di, dj]:.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
