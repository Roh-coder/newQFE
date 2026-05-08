#!/usr/bin/env python3
"""
zoo_ising_tau_scan.py — discriminate cost variants by scanning a 2-D
family of test geometries against a fixed Ising-CFT reference, using
the Jacobi-theta torus correlator as ground truth (no MC).

Setup
-----
  ref:  geometry (Lx_ref, Ly_ref, Tx_ref, Ty_ref) → periods (ω₁_r, ω₂_r), τ_r
  test: geometry (Lx_t,   Ly_t,   Tx_t,   Ty_t)   → periods (ω₁_t, ω₂_t)

Then for each (a, b) on a 2-D scan we override the test periods to
    ω₁_test_eff = a * ω₁_t
    ω₂_test_eff = b * ω₂_t  (real scaling along ω₂ only)
which makes τ_test_eff = (b/a) · τ_t, so the scan continuously varies
the test modular parameter.  The "TRUTH" point is (a*, b*) such that
τ_test_eff = τ_ref, i.e. b*/a* = τ_ref / τ_test (as complex ratios; we
project to the real axis since this 2-D scan only varies real
scales — the truth set is a 1-D curve, not a point, in general).

For each (a, b) we build the Ising data dict, evaluate every cost
variant against the ref Ising data, and record the cost.  A good cost
is one whose minimum on the scan grid lies on or near the truth curve.

Usage:  python zoo_ising_tau_scan.py [--n N] [--span S]
"""
from __future__ import annotations

import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cft_torus as cft
import cost as cost_mod
import cost_zoo

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=15,
                    help="grid resolution (n × n)")
    ap.add_argument("--span", type=float, default=0.6,
                    help="±span around the truth scaling")
    ap.add_argument("--ref",  type=int, nargs=4, default=[13, 16, -3, 3])
    ap.add_argument("--test", type=int, nargs=4, default=[16, 16,  0, 0])
    args = ap.parse_args()

    Rx, Ry, RTx, RTy = args.ref
    Tx, Ty, TTx, TTy = args.test
    print(f"  ref geom  = ({Rx}, {Ry}, {RTx}, {RTy})")
    print(f"  test geom = ({Tx}, {Ty}, {TTx}, {TTy})")

    w1_r, w2_r = cft.torus_periods(Rx, Ry, RTx, RTy)
    w1_t, w2_t = cft.torus_periods(Tx, Ty, TTx, TTy)
    tau_r = w2_r / w1_r
    tau_t = w2_t / w1_t
    if tau_r.imag < 0: tau_r = complex(tau_r.real, -tau_r.imag)
    if tau_t.imag < 0: tau_t = complex(tau_t.real, -tau_t.imag)

    print(f"  ω₁_r={w1_r}  ω₂_r={w2_r}  τ_r={tau_r}")
    print(f"  ω₁_t={w1_t}  ω₂_t={w2_t}  τ_t={tau_t}")
    print(f"  τ_r / τ_t = {tau_r / tau_t}")
    # Truth curve: b/a = (τ_r / τ_t).  In our real-only scan we project
    # to the real axis, so truth on the (a, b) grid is the line
    # b/a = Re(τ_r/τ_t) — only an approximation when τ_r/τ_t is complex.
    truth_ratio = (tau_r / tau_t).real
    print(f"  TRUTH ratio b/a = Re(τ_r/τ_t) = {truth_ratio:.4f}")

    # Build ref Ising correlator data once.
    print("  building ref Ising correlator…")
    ref_data = cft.make_ising_data_dict(Rx, Ry, RTx, RTy)

    # Scan grid centred on (a, b) = (1, truth_ratio).
    a_centre, b_centre = 1.0, truth_ratio
    a_lo, a_hi = a_centre * (1 - args.span), a_centre * (1 + args.span)
    b_lo, b_hi = b_centre * (1 - args.span), b_centre * (1 + args.span)
    a_arr = np.linspace(a_lo, a_hi, args.n)
    b_arr = np.linspace(b_lo, b_hi, args.n)
    print(f"  scan: a ∈ [{a_lo:.3f}, {a_hi:.3f}],  "
          f"b ∈ [{b_lo:.3f}, {b_hi:.3f}],  {args.n}×{args.n}")

    # Cost variants under test.
    variants = {
        "production":  cost_mod._l2_cost_test_native,  # current production
        **{f"zoo_{k}": v for k, v in cost_zoo.ZOO.items()},
    }
    print(f"  variants: {list(variants.keys())}")

    grids = {name: np.full((args.n, args.n), np.nan) for name in variants}
    nan_count = {name: 0 for name in variants}

    t0 = time.perf_counter()
    n_total = args.n * args.n
    cell = 0
    # silence cost-function chatter
    import io
    _real_stdout = sys.stdout
    null = io.StringIO()
    for i, a in enumerate(a_arr):
        for j, b in enumerate(b_arr):
            cell += 1
            if cell % max(1, n_total // 20) == 0:
                _real_stdout.write(f"    cell {cell}/{n_total}  "
                                   f"({100*cell/n_total:.0f}%)\n")
                _real_stdout.flush()
            # build test data with overridden periods
            try:
                test_data = cft.make_ising_data_dict(
                    Tx, Ty, TTx, TTy,
                    omega1_override=a * w1_t,
                    omega2_override=b * w2_t)
            except Exception as e:
                continue
            for name, fn in variants.items():
                sys.stdout = null
                try:
                    c, s, pd, ps = fn(ref_data, test_data,
                                      Tx, Ty, TTx, TTy,
                                      Rx, Ry, RTx, RTy,
                                      copies=2, power=2)
                    if math.isfinite(c):
                        grids[name][i, j] = c
                    else:
                        nan_count[name] += 1
                except Exception:
                    nan_count[name] += 1
                finally:
                    sys.stdout = _real_stdout
    wall = time.perf_counter() - t0
    print(f"  scan done in {wall:.1f}s")

    # Score each variant: where is the min, how close to truth-line?
    print()
    print("  variant          min(a,b)         min_cost      "
          "dist_to_truth   nan_cells")
    print("  " + "-" * 76)
    summary = dict(ref=args.ref, test=args.test,
                   tau_r=[tau_r.real, tau_r.imag],
                   tau_t=[tau_t.real, tau_t.imag],
                   truth_ratio=truth_ratio,
                   span=args.span, n=args.n,
                   a_arr=a_arr.tolist(), b_arr=b_arr.tolist(),
                   variants={})
    for name in variants:
        g = grids[name]
        if np.all(np.isnan(g)):
            print(f"  {name:<15}  ALL NaN")
            summary["variants"][name] = dict(all_nan=True)
            continue
        i_min, j_min = np.unravel_index(np.nanargmin(g), g.shape)
        a_min, b_min = float(a_arr[i_min]), float(b_arr[j_min])
        c_min = float(g[i_min, j_min])
        # distance from (a_min, b_min) to the truth-LINE b = truth_ratio·a
        # in the (a, b) plane (perpendicular distance from a point to a line).
        # line: truth_ratio·a − b = 0 → unit normal (truth_ratio, -1)/√(1+r²)
        dist = abs(truth_ratio * a_min - b_min) / math.sqrt(
            truth_ratio ** 2 + 1.0)
        print(f"  {name:<15}  ({a_min:.3f},{b_min:.3f})  {c_min:.4e}    "
              f"{dist:.4f}     {nan_count[name]}/{n_total}")
        summary["variants"][name] = dict(
            a_min=a_min, b_min=b_min, c_min=c_min,
            dist_to_truth=dist,
            nan_count=nan_count[name],
            grid=g.tolist())

    with open(os.path.join(OUT, "zoo_ising_scan_summary.json"), "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    # Plot all variants on a 2 × ⌈N/2⌉ grid.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ncols = 3
        nrows = (len(variants) + ncols - 1) // ncols
        fig, axs = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 4 * nrows),
                                squeeze=False)
        for k, (name, fn) in enumerate(variants.items()):
            r, c = k // ncols, k % ncols
            ax = axs[r, c]
            g = grids[name]
            if np.all(np.isnan(g)):
                ax.set_title(f"{name}: all NaN"); continue
            # log scale where positive, with NaN cells masked
            with np.errstate(invalid="ignore"):
                logg = np.where(g > 0, np.log10(np.abs(g) + 1e-30), np.nan)
            im = ax.imshow(logg.T, origin="lower",
                           extent=[a_arr[0], a_arr[-1], b_arr[0], b_arr[-1]],
                           aspect="auto", cmap="viridis")
            plt.colorbar(im, ax=ax, label="log₁₀ |cost|")
            # truth line
            aa = np.linspace(a_arr[0], a_arr[-1], 100)
            ax.plot(aa, truth_ratio * aa, "--", color="red", lw=2,
                    label="truth: b=τ_r/τ_t · a")
            # min marker
            i_min, j_min = np.unravel_index(np.nanargmin(g), g.shape)
            ax.scatter([a_arr[i_min]], [b_arr[j_min]], marker="X",
                       s=140, color="orange", edgecolor="black",
                       zorder=10, label="cost min")
            ax.set_xlabel("a (scale ω₁)")
            ax.set_ylabel("b (scale ω₂)")
            ax.set_title(name)
            ax.legend(fontsize=7, loc="lower right")
        # turn off any extra panels
        for k in range(len(variants), nrows * ncols):
            axs[k // ncols, k % ncols].axis("off")
        plt.tight_layout()
        png = os.path.join(OUT, "zoo_ising_scan.png")
        plt.savefig(png, dpi=120); plt.close(fig)
        print(f"  plot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"  plot failed: {e}")


if __name__ == "__main__":
    main()
