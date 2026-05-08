"""
viz_residuals_normalized.py — Same 10 (r1,r2) points but with the
per-direction multiplicative offset removed.

For each (r1,r2) and each direction, fit a single scalar c such that
log G_test - log G_ref - log c = 0 in least-squares sense
(i.e. c = exp(mean(log G_test - log G_ref))). Then plot:
  top:    log G_test/c (markers) overlaid on log G_ref (line)
  bottom: residual after offset removed
"""
from __future__ import annotations

import os, sys, pickle
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
        out.append(dict(t=t_arr, Gr=Gr, Gt=Gt))
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
    # collect per-direction c values to print
    c_table = []
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
        cs = []
        rms = []
        for k, d in enumerate(dirs):
            if d is None: cs.append(np.nan); rms.append(np.nan); continue
            t = d["t"]; Gr = d["Gr"]; Gt = d["Gt"]
            m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > 0) & (Gt > 0)
            if m.sum() < 2:
                cs.append(np.nan); rms.append(np.nan); continue
            tt = t[m]
            log_diff = np.log(Gt[m]) - np.log(Gr[m])
            log_c = float(np.mean(log_diff))
            c = float(np.exp(log_c))
            cs.append(c)
            res_after = log_diff - log_c
            rms.append(float(np.sqrt(np.mean(res_after * res_after))))
            color = DIR_COLORS[k]
            ax_top.semilogy(tt, Gr[m], "-", color=color, lw=1.2,
                            label=f"{DIR_NAMES[k]} ref" if col == 0 else None)
            ax_top.semilogy(tt, Gt[m] / c, "o", color=color, ms=3.5,
                            label=f"{DIR_NAMES[k]} test/c" if col == 0 else None)
            ax_bot.plot(tt, res_after, "-o", color=color, ms=3.5, lw=1.0)
        c_table.append((r1, r2, label, cs, rms))

        is_truth = (abs(r1 - TRUTH[0]) < 0.6 and abs(r2 - TRUTH[1]) < 0.6)
        title_color = "red" if is_truth else "black"
        ttl = (f"({r1},{r2})\n{label}\n"
               f"c=[{cs[0]:.2f},{cs[1]:.2f},{cs[2]:.2f}]")
        ax_top.set_title(ttl, fontsize=8, color=title_color)
        ax_top.set_xlim(0, 0.5)
        ax_bot.axhline(0, color="k", lw=0.6)
        ax_bot.set_xlim(0, 0.5)
        ax_bot.set_ylim(-0.25, 0.25)
        if col == 0:
            ax_top.set_ylabel("G(t)  [test rescaled by c]")
            ax_bot.set_ylabel("log(G_test/c) − log G_ref")
        ax_bot.set_xlabel("t = k/N")

    axes[0, 0].legend(loc="lower left", fontsize=6, ncol=2)
    fig.suptitle(f"Per-axis residuals AFTER multiplicative offset removed  "
                 f"test={TEST_GEOM}  ref={REF_GEOM}  truth=(5.07, 7.74)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(root, "residuals_10pts_normalized.png")
    fig.savefig(out, dpi=130)
    print("→", out)

    # print c and shape-residual table
    print(f"\n{'pt':<22} {'c_u':>6} {'c_v':>6} {'c_w':>6}  "
          f"{'rms_u':>7} {'rms_v':>7} {'rms_w':>7}  spread(c)")
    for r1, r2, label, cs, rms in c_table:
        spread = float(np.std(cs))
        rms_total = float(np.sqrt(np.mean(np.asarray(rms)**2)))
        print(f"({r1:.1f},{r2:.1f}) {label:<10}  "
              f"{cs[0]:6.3f} {cs[1]:6.3f} {cs[2]:6.3f}  "
              f"{rms[0]:7.4f} {rms[1]:7.4f} {rms[2]:7.4f}  "
              f"{spread:.3f}   rms_tot={rms_total:.4f}")


if __name__ == "__main__":
    main()
