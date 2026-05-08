"""
verify_cycles.py — sanity-check that u,v,w cycle vectors are
constructed correctly for both test and ref lattices and that the
resulting G(t) along each cycle picks up the singularities at t=0 and
t=1 (the periodic-torus return points).

Prints the (m,n) period vectors, their physical (xy) coords, lengths,
angles, and the modular tau extracted from each. Plots G(t) for
t in [-0.05, 1.05] along each cycle on each lattice (interp on ref,
integer-site lookup on test).
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
DIR_NAMES = ["u", "v", "w"]
DIR_COLORS = ["tab:blue", "tab:orange", "tab:green"]


def describe(geom, label):
    Lx, Ly, Tx, Ty = geom
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    print(f"\n=== {label}  geom=(Lx={Lx}, Ly={Ly}, Tx={Tx}, Ty={Ty}) ===")
    info = []
    for k, (dm, dn) in enumerate(paths):
        x = dm + 0.5 * dn
        y = _SQRT3_2 * dn
        plen = float(np.hypot(x, y))
        angle = float(np.degrees(np.arctan2(y, x)))
        # number of integer sites along the cycle: g = gcd(|dm|,|dn|)
        from math import gcd
        g = gcd(abs(int(dm)), abs(int(dn))) or max(abs(int(dm)), abs(int(dn)))
        info.append((dm, dn, x, y, plen, angle, g))
        print(f"  {DIR_NAMES[k]}: (m,n)=({dm:+4d},{dn:+4d})   "
              f"xy=({x:+8.3f}, {y:+8.3f})   |p|={plen:7.3f}   "
              f"angle={angle:+7.2f}°   #sites along cycle = {g}")
    sum_mn = (sum(p[0] for p in info), sum(p[1] for p in info))
    print(f"  sum of (m,n) across u+v+w = {sum_mn}  (must be (0,0))")
    # ratio of cycle lengths (sorted)
    lens = sorted([i[4] for i in info])
    print(f"  sorted lens = {[f'{x:.3f}' for x in lens]}   "
          f"ratios vs shortest = {[f'{x/lens[0]:.4f}' for x in lens]}")
    # modular tau via Brower-Owen Eq. 56:
    #   |tau| = ell_2 / ell_1 (ell_1 <= ell_2 <= ell_3)
    #   arg(tau) = arccos(-e1*.e2*) but with our convention let's
    #   instead just print the (Lx,Ly,Tx,Ty)→tau via the standard
    #   parallelogram parametrisation: tau = (Tx + 0.5 Ly + i sqrt3/2 Ly)/Lx
    Lx_, Ly_, Tx_, Ty_ = geom
    tau_re = (Tx_ + 0.5 * Ly_) / Lx_
    tau_im = _SQRT3_2 * Ly_ / Lx_
    print(f"  parallelogram tau = {tau_re:.4f} + {tau_im:.4f} i   "
          f"|tau|={np.hypot(tau_re, tau_im):.4f}")
    return info


def trace_G_test(test_data, geom, n_sites=None):
    """Walk integer sites along each cycle ON THE TEST LATTICE for
    t = k/g, k = 0,..,g-1. Optionally extend a bit past t=1 by walking
    further (which should wrap back via torus equivalence)."""
    Lx, Ly, Tx, Ty = geom
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    out = []
    for (dm, dn) in paths:
        from math import gcd
        g = gcd(abs(int(dm)), abs(int(dn))) or max(abs(int(dm)), abs(int(dn)))
        # integer sites at k = 0..2g-1 (i.e. covering [0,2] in t)
        nk = 2 * g
        ts, Gs = [], []
        for k in range(nk):
            mk = k * dm // g
            nkn = k * dn // g
            entry = _lookup_test_value(test_data, int(mk), int(nkn),
                                       Lx, Ly, Tx, Ty, copies=2)
            ts.append(k / g)
            Gs.append(entry["conn"] if entry is not None else np.nan)
        out.append((np.asarray(ts), np.asarray(Gs)))
    return out


def trace_G_ref(iref, geom, n_pts=121, t_lo=-0.05, t_hi=1.05):
    """Sample the ref interpolator densely along each cycle, including
    a bit past 0 and 1 to see the singularity wrap behaviour."""
    paths = boundary_paths(*geom)
    out = []
    ts = np.linspace(t_lo, t_hi, n_pts)
    for (dm, dn) in paths:
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        pts = np.column_stack([ts * ex, ts * ey])
        Gs = np.asarray(iref(pts), dtype=float)
        out.append((ts, Gs))
    return out


def main():
    info_ref  = describe(REF_GEOM,  "REFERENCE (39,48,-9,9)")
    info_test = describe(TEST_GEOM, "TEST      (16,16,0,0)")

    # load data + interpolator
    ref_dir = os.path.join(_HERE, "results", REF_TAG)
    ref_data = mc_engine.load_all_to_all(
        os.path.join(ref_dir, "two_point_all_to_all.dat"))
    iref = _tile_interp(ref_data, *REF_GEOM, "conn", 2)

    # use a precompute pkl near truth for the test correlator
    landscape = os.path.join(_HERE, "results", "_landscape", LANDSCAPE)
    pkl = os.path.join(landscape, "grid", "r1_5.000_r2_7.500.pkl")
    if not os.path.exists(pkl):
        pkl = os.path.join(landscape, "grid", "r1_5.000_r2_8.000.pkl")
    print(f"\n[load] test data from {pkl}")
    with open(pkl, "rb") as f:
        pt = pickle.load(f)
    test_traces = trace_G_test(pt["test_data"], TEST_GEOM)
    ref_traces  = trace_G_ref(iref, REF_GEOM)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for k in range(3):
        # ref panel (top)
        axt = axes[0, k]
        ts_r, Gs_r = ref_traces[k]
        axt.plot(ts_r, Gs_r, "-", color=DIR_COLORS[k], lw=1.5)
        for tline in [0.0, 1.0]:
            axt.axvline(tline, color="k", lw=0.6, ls="--")
        axt.set_title(f"REF cycle {DIR_NAMES[k]}\n"
                      f"|p|={info_ref[k][4]:.2f}  ang={info_ref[k][5]:.1f}°  "
                      f"sites={info_ref[k][6]}")
        axt.set_yscale("log"); axt.set_ylim(0.05, 1.5)
        if k == 0: axt.set_ylabel("G_ref(t)  [interp]")

        # test panel (bottom)
        axb = axes[1, k]
        ts_t, Gs_t = test_traces[k]
        axb.plot(ts_t, Gs_t, "o-", color=DIR_COLORS[k], lw=1.0, ms=4)
        for tline in [0.0, 1.0]:
            axb.axvline(tline, color="k", lw=0.6, ls="--")
        axb.set_title(f"TEST cycle {DIR_NAMES[k]}\n"
                      f"|p|={info_test[k][4]:.2f}  ang={info_test[k][5]:.1f}°  "
                      f"sites={info_test[k][6]}")
        axb.set_xlabel("t (period fractions)")
        axb.set_yscale("log"); axb.set_ylim(0.05, 1.5)
        if k == 0: axb.set_ylabel("G_test(t)  [integer sites]")
    fig.suptitle(f"u/v/w cycles: G(t) over t ∈ [-0.05, 2.0]\n"
                 f"vertical dashed lines = period boundaries (singularities)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(landscape, "verify_cycles.png")
    fig.savefig(out, dpi=130)
    print(f"\n→ {out}")


if __name__ == "__main__":
    main()
