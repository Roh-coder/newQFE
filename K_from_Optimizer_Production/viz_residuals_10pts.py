"""
viz_residuals_10pts.py — visualize per-axis G_test vs G_ref residuals
at 10 selected (r1, r2) points.

For each (r1,r2) we extract the 3 cycle correlators on test (16,16,0,0)
and the matched-fractional-position interpolation on ref (39,48,-9,9),
then plot:
  top row:  log G_ref (line) and log G_test (markers) vs t in [0,0.5)
  bottom:   residual log G_test - log G_ref vs t

3 directions per (r1,r2) point, color-coded.
"""
from __future__ import annotations

import os, sys, pickle, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)

REF_TAG  = "_reference_Lx39_Ly48_Tx-9_Ty9"
REF_GEOM = (39, 48, -9, 9)
TEST_GEOM = (16, 16, 0, 0)
LANDSCAPE = "Lx16_Ly16_Tx0_Ty0"
TRUTH = (5.0652, 7.7429)

# 10 points chosen to span the landscape
POINTS = [
    (5.0, 7.5, "near truth"),
    (5.0, 8.0, "near truth +"),
    (4.5, 4.5, "powerlaw min"),
    (2.5, 2.5, "bowl min"),
    (1.0, 1.0, "isotropic"),
    (0.5, 0.5, "small corner"),
    (8.0, 8.0, "diag large"),
    (3.0, 6.0, "moderate"),
    (7.0, 3.0, "skew"),
    (10.0, 10.0, "far corner"),
]

DIR_COLORS = ["tab:blue", "tab:orange", "tab:green"]
DIR_NAMES  = ["u", "v", "w"]


def extract(test_data, iref):
    """Return list of 3 dicts: per direction (t, G_ref, G_test, p_ref_len, p_test_len)."""
    ref_paths  = boundary_paths(*REF_GEOM)
    test_paths = boundary_paths(*TEST_GEOM)
    out = []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 3: out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5*rdn, _SQRT3_2*rdn
        pts_ref = np.column_stack([t_arr*rex, t_arr*rey])
        Gr = np.asarray(iref(pts_ref), dtype=float)
        Gt = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       *TEST_GEOM, copies=2)
            if entry is not None:
                Gt[i] = entry["conn"]
        p_ref_len  = float(np.hypot(rex, rey))
        tex, tey   = tdm + 0.5*tdn, _SQRT3_2*tdn
        p_test_len = float(np.hypot(tex, tey))
        out.append(dict(t=t_arr, Gr=Gr, Gt=Gt,
                        pr=p_ref_len, pt=p_test_len))
    return out


def main():
    root = os.path.join(_HERE, "results", "_landscape", LANDSCAPE)
    ref_dir = os.path.join(_HERE, "results", REF_TAG)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    iref = _tile_interp(ref_data, *REF_GEOM, "conn", 2)

    n_pts = len(POINTS)
    fig, axes = plt.subplots(2, n_pts, figsize=(2.6*n_pts, 5.5),
                             sharex=True)
    for col, (r1, r2, label) in enumerate(POINTS):
        pkl = os.path.join(root, "grid", f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
        if not os.path.exists(pkl):
            for row in range(2):
                axes[row, col].set_title(f"({r1},{r2})\nMISSING", fontsize=8)
                axes[row, col].axis("off")
            continue
        with open(pkl, "rb") as f:
            pt = pickle.load(f)
        dirs = extract(pt["test_data"], iref)

        ax_top = axes[0, col]; ax_bot = axes[1, col]
        for k, d in enumerate(dirs):
            if d is None: continue
            t = d["t"]
            Gr = d["Gr"]; Gt = d["Gt"]
            m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > 0) & (Gt > 0)
            if m.sum() < 2: continue
            tt = t[m]
            color = DIR_COLORS[k]
            ax_top.semilogy(tt, Gr[m], "-", color=color, lw=1.2,
                            label=f"{DIR_NAMES[k]} ref" if col == 0 else None)
            ax_top.semilogy(tt, Gt[m], "o", color=color, ms=3.5,
                            label=f"{DIR_NAMES[k]} test" if col == 0 else None)
            res = np.log(Gt[m]) - np.log(Gr[m])
            ax_bot.plot(tt, res, "-o", color=color, ms=3.5, lw=1.0)

        is_truth = (abs(r1 - TRUTH[0]) < 0.6 and abs(r2 - TRUTH[1]) < 0.6)
        title_color = "red" if is_truth else "black"
        ax_top.set_title(f"({r1},{r2})\n{label}", fontsize=9,
                         color=title_color)
        ax_top.set_xlim(0, 0.5)
        ax_bot.axhline(0, color="k", lw=0.6)
        ax_bot.set_xlim(0, 0.5)
        ax_bot.set_ylim(-1.0, 1.0)
        if col == 0:
            ax_top.set_ylabel("G(t)")
            ax_bot.set_ylabel("log G_test − log G_ref")
        ax_bot.set_xlabel("t = k/N")

    axes[0, 0].legend(loc="lower left", fontsize=6, ncol=2)
    fig.suptitle(f"Per-axis correlator residuals  test={TEST_GEOM}  ref={REF_GEOM}  "
                 f"truth=(5.07, 7.74)  [u=blue, v=orange, w=green]",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(root, "residuals_10pts.png")
    fig.savefig(out, dpi=130)
    print("→", out)


if __name__ == "__main__":
    main()
