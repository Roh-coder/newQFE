#!/usr/bin/env python3
"""
zoo_ising_boundary_scan.py — discriminate cost variants with NO MC noise.

Pick a fixed REFERENCE torus geometry and compute its analytic Ising-CFT
spin-spin two-point function (Jacobi theta) along its 3 boundary
directions.  Then sweep a one-parameter family of TEST geometries (same
discrete (Lx, Ly, Tx, Ty), but rescale ω₂ by a real factor s) so the
test modular parameter τ_test = s·τ_ref sweeps through the truth
(s = 1) and the optimizer's "trap" region (s ≪ 1 or ≫ 1).

For each test geometry we sample the analytic Ising correlator along
the test's 3 boundary directions and feed (ref_data, test_data) to
every cost variant.  Output: per-cost curve cost(s) with the truth
marked at s = 1 — a good cost is one whose minimum sits at s = 1 with
a clean, monotone basin.

Usage:
    python zoo_ising_boundary_scan.py [--n N] [--smin S] [--smax S]
                                       [--ref Lx Ly Tx Ty]
                                       [--test Lx Ly Tx Ty]

Notes
-----
The 5 zoo variants live in ``cost_zoo`` and the 1 production variant is
``cost._l2_cost_test_native``.  All take the same signature.

Boundary-direction sampling is exactly what every test_native-family
cost already does internally (via _sample_per_direction_native or the
production sibling).  We pass them whole-torus data dicts; the costs
themselves restrict to the boundary paths.  This means the per-direction
N (number of sample points) is determined by the geometry's gcds — not a
free parameter — and matches what runs in production.
"""
from __future__ import annotations

import argparse, io, json, math, os, sys, time
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
    ap.add_argument("--n",    type=int,   default=41,
                    help="number of s values in the sweep")
    ap.add_argument("--smin", type=float, default=0.4)
    ap.add_argument("--smax", type=float, default=2.5)
    ap.add_argument("--ref",  type=int, nargs=4, default=[13, 16, -3, 3])
    ap.add_argument("--test", type=int, nargs=4, default=[16, 16,  0, 0])
    args = ap.parse_args()

    Rx, Ry, RTx, RTy = args.ref
    Tx, Ty, TTx, TTy = args.test
    print(f"  ref  geom: ({Rx}, {Ry}, {RTx}, {RTy})")
    print(f"  test geom: ({Tx}, {Ty}, {TTx}, {TTy})")

    w1_r, w2_r = cft.torus_periods(Rx, Ry, RTx, RTy)
    w1_t, w2_t = cft.torus_periods(Tx, Ty, TTx, TTy)
    tau_r = w2_r / w1_r
    tau_t = w2_t / w1_t
    if tau_r.imag < 0: tau_r = complex(tau_r.real, -tau_r.imag)
    if tau_t.imag < 0: tau_t = complex(tau_t.real, -tau_t.imag)
    print(f"  τ_ref  = {tau_r:.4f}")
    print(f"  τ_test = {tau_t:.4f}  (will be rescaled by s)")

    # The "truth" point is when test data uses the SAME analytic
    # correlator at the SAME modular point as the ref → set test
    # periods to (w1_r, w2_r).  When s = 1 we match ref's τ exactly.
    # For s != 1 we scale ω₂_test by s so τ_test_eff = s · τ_ref.

    # Build ref Ising data ONCE, using ref's actual modular parameter.
    print("  building ref Ising correlator once…")
    ref_data = cft.make_ising_data_dict(Rx, Ry, RTx, RTy)

    s_arr = np.linspace(args.smin, args.smax, args.n)

    variants = {
        "production":  cost_mod._l2_cost_test_native,
        **{f"zoo_{k}": v for k, v in cost_zoo.ZOO.items()},
    }
    print(f"  variants: {list(variants.keys())}")
    print(f"  scan: s ∈ [{args.smin}, {args.smax}],  N = {args.n}")

    curves = {n: np.full(args.n, np.nan) for n in variants}
    nan_count = {n: 0 for n in variants}

    # silence cost-function chatter
    real_stdout = sys.stdout
    null = io.StringIO()
    t0 = time.perf_counter()
    for k, s in enumerate(s_arr):
        # Test data: same DISCRETE (Lx, Ly, Tx, Ty) so the test boundary
        # paths and gcds are unchanged, but we override the *physical*
        # periods so the analytic correlator lives on a different torus.
        # Ref periods at s=1 are (w1_r, w2_r); for general s we use
        # ω₁_eff = w1_r,  ω₂_eff = s·w2_r.  This puts τ_test_eff = s·τ_ref.
        try:
            test_data = cft.make_ising_data_dict(
                Tx, Ty, TTx, TTy,
                omega1_override=w1_r,
                omega2_override=s * w2_r)
        except Exception as e:
            real_stdout.write(f"    s={s:.3f}: build failed: {e}\n")
            continue

        for name, fn in variants.items():
            sys.stdout = null
            try:
                c, sigma, pd, ps = fn(ref_data, test_data,
                                      Tx, Ty, TTx, TTy,
                                      Rx, Ry, RTx, RTy,
                                      copies=2, power=2)
                if math.isfinite(c):
                    curves[name][k] = c
                else:
                    nan_count[name] += 1
            except Exception:
                nan_count[name] += 1
            finally:
                sys.stdout = real_stdout

        if (k + 1) % max(1, args.n // 10) == 0:
            real_stdout.write(f"    {k+1}/{args.n}  s={s:.3f}\n")
            real_stdout.flush()
    wall = time.perf_counter() - t0
    print(f"  scan done in {wall:.1f}s")

    # Score each variant: where is the min, AND how steeply does cost
    # respond to anisotropy (i.e. how sensitive)?
    print()
    hdr = (f"  {'variant':<15}  {'s_min':>6}  {'|s−1|':>6}  "
           f"{'curv@1':>9}  {'asym(2)':>8}  {'asym(½)':>8}  "
           f"{'mono':>5}  {'sens_score':>10}")
    print(hdr)
    print("  " + "-" * (len(hdr)))
    summary = dict(ref=args.ref, test=args.test,
                   tau_r=[tau_r.real, tau_r.imag],
                   tau_t=[tau_t.real, tau_t.imag],
                   s=s_arr.tolist(), variants={})

    # Helper: interpolate cost at a target s
    def _at(c, s_target):
        m = np.isfinite(c)
        if m.sum() < 2: return float("nan")
        return float(np.interp(s_target, s_arr[m], c[m]))

    rows = []
    for name in variants:
        c = curves[name]
        if np.all(np.isnan(c)):
            print(f"  {name:<15}  ALL NaN")
            summary["variants"][name] = dict(all_nan=True); continue
        i_min = int(np.nanargmin(c))
        s_min = float(s_arr[i_min])
        c_min = float(c[i_min])
        c_at_1 = _at(c, 1.0)
        c_at_05 = _at(c, 0.5)
        c_at_2  = _at(c, 2.0)

        # Anisotropy sensitivity: ratios cost(off-truth)/cost(at-truth)
        # → larger means the cost rises more sharply when anisotropy is
        # wrong.  An ideal cost has both ratios > 1 (basin around s=1).
        asym2  = c_at_2  / c_at_1 if c_at_1 > 0 else float("nan")
        asym05 = c_at_05 / c_at_1 if c_at_1 > 0 else float("nan")

        # Log-curvature at s=1 estimated from the 3 neighbours of s=1.
        i1 = int(np.argmin(np.abs(s_arr - 1.0)))
        if 1 <= i1 < args.n - 1:
            cL, c0, cR = c[i1 - 1], c[i1], c[i1 + 1]
            ds = (s_arr[i1 + 1] - s_arr[i1 - 1]) / 2.0
            with np.errstate(invalid="ignore", divide="ignore"):
                if np.all(np.array([cL, c0, cR]) > 0):
                    curv = float((np.log(cL) + np.log(cR) - 2 * np.log(c0))
                                 / (ds ** 2))
                else:
                    curv = float("nan")
        else:
            curv = float("nan")

        # Monotone-down score
        i_truth = int(np.argmin(np.abs(s_arr - 1.0)))
        with np.errstate(invalid="ignore"):
            ldiff = np.diff(c[:i_truth + 1][np.isfinite(c[:i_truth + 1])])
            rdiff = np.diff(c[i_truth:][np.isfinite(c[i_truth:])])
            ls = float(np.mean(ldiff < 0)) if len(ldiff) else 1.0
            rs = float(np.mean(rdiff > 0)) if len(rdiff) else 1.0
        mono = 0.5 * (ls + rs)

        # Composite sensitivity score:
        #   s = log10(min(asym2, asym05))  (worst-direction signal-to-truth)
        # Higher is better.  Negative or nan = cost goes the wrong way at
        # one of the test points → unusable for optimization.
        worst_asym = min(asym2, asym05) if (np.isfinite(asym2)
                                            and np.isfinite(asym05)) else float("nan")
        sens = math.log10(worst_asym) if (worst_asym
                                          and worst_asym > 0) else float("nan")

        print(f"  {name:<15}  {s_min:>6.3f}  {abs(s_min-1):>6.3f}  "
              f"{curv:>9.2e}  {asym2:>8.2f}  {asym05:>8.2f}  "
              f"{mono:>5.2f}  {sens:>10.3f}")
        summary["variants"][name] = dict(
            s_min=s_min, c_min=c_min,
            err_s=abs(s_min - 1.0),
            curv_log_at_1=curv,
            asym_s_2=asym2, asym_s_05=asym05,
            sensitivity_score=sens,
            monotone_score=mono,
            curve=c.tolist())
        rows.append((name, sens, abs(s_min - 1)))

    # Final ranking by sensitivity score (descending), with err_s tie-break
    print()
    print("  RANKING (by sensitivity score; higher → cost rises faster off truth):")
    rows_sorted = sorted([r for r in rows if math.isfinite(r[1])],
                         key=lambda r: (-r[1], r[2]))
    for i, (name, s, errs) in enumerate(rows_sorted, 1):
        print(f"    {i}. {name:<15}  sens={s:+.3f}  err_s={errs:.3f}")
    rows_bad = [r for r in rows if not math.isfinite(r[1])]
    for name, _, _ in rows_bad:
        print(f"    --. {name:<15}  (sens=NaN — not usable)")

    with open(os.path.join(OUT, "zoo_ising_boundary_scan_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    # Plot all variants on a single panel.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax = axs[0]
        colours = plt.cm.tab10(np.linspace(0, 1, len(variants)))
        for (name, c), col in zip(curves.items(), colours):
            mask = np.isfinite(c)
            if mask.sum() < 2: continue
            # normalise each curve to its min so we can compare shapes
            c_norm = c / np.nanmin(c[c > 0]) if np.any(c > 0) else c
            ax.plot(s_arr[mask], c_norm[mask], "-o", ms=4, color=col,
                    label=name, lw=1.5)
            i_min = int(np.nanargmin(c))
            ax.scatter([s_arr[i_min]], [c_norm[i_min]], marker="X",
                       s=100, color=col, edgecolor="black", zorder=10)
        ax.axvline(1.0, color="red", lw=2, ls="--", label="TRUTH (s=1)")
        ax.set_ylabel("cost / cost_min  (per-variant normalised)")
        ax.set_yscale("log")
        ax.set_title(f"Cost vs ω₂ scaling s  (ref={tuple(args.ref)}  "
                     f"test={tuple(args.test)})")
        ax.legend(loc="best", fontsize=8, ncol=2)
        ax.grid(alpha=0.3, which="both")

        ax = axs[1]
        for (name, c), col in zip(curves.items(), colours):
            mask = np.isfinite(c)
            if mask.sum() < 2: continue
            ax.plot(s_arr[mask], c[mask], "-o", ms=4, color=col,
                    label=name, lw=1.5)
        ax.axvline(1.0, color="red", lw=2, ls="--", label="TRUTH (s=1)")
        ax.set_xlabel("s  (test ω₂ rescaling factor)")
        ax.set_ylabel("absolute cost")
        ax.set_yscale("log")
        ax.legend(loc="best", fontsize=8, ncol=2)
        ax.grid(alpha=0.3, which="both")

        plt.tight_layout()
        png = os.path.join(OUT, "zoo_ising_boundary_scan.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")


if __name__ == "__main__":
    main()
