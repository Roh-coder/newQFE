#!/usr/bin/env python3
"""
cost_zoo_analytic_test.py — noise-free A/B of every cost variant on
the 4-5-6 problem using the analytic torus correlator from cft_torus.py.

Setup
-----
ref geometry  = (13, 16, -3, 3) → ω₁_ref, ω₂_ref → τ_ref
test geometry = (16, 16,  0, 0) → ω₁_test, ω₂_test → τ_test

Scan parameter (r1, r2) ∈ ℝ² acts as a COMPLEX stretch factor
(s := r1 + i r2) applied to ω₂_test:

    ω₂_test_eff(r1, r2) = (r1 + i r2) · ω₂_test
    ⇒ τ_test_eff(r1, r2) = (r1 + i r2) · τ_test

The TRUE point (r1*, r2*) is the complex value of  τ_ref / τ_test.
At truth the test torus has the same modular parameter as the ref
torus, so the curve-collapse hypothesis predicts cost = 0.

We scan a grid of (r1, r2), generate test data analytically at every
point, and evaluate every cost variant (production test_native +
spectral + 5 zoo variants).  Plot 7 heatmaps; mark truth.

Run:
    python cost_zoo_analytic_test.py [--n 21] [--m_sq 0.01]
"""
from __future__ import annotations

import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=21,
                    help="grid resolution per axis (n × n cost surfaces)")
    ap.add_argument("--m_sq", type=float, default=0.01,
                    help="IR mass² for free-field correlator (small = critical)")
    ap.add_argument("--n_cut", type=int, default=18)
    ap.add_argument("--r1_range", type=float, nargs=2, default=[0.3, 2.0])
    ap.add_argument("--r2_range", type=float, nargs=2, default=[-0.5, 1.3])
    args = ap.parse_args()

    import cft_torus, cost as cost_mod, cost_zoo

    ref_geom  = (13, 16, -3,  3)
    test_geom = (16, 16,  0,  0)
    w1r, w2r  = cft_torus.torus_periods(*ref_geom)
    w1t, w2t  = cft_torus.torus_periods(*test_geom)
    tau_ref   = cft_torus.modular_tau(*ref_geom)
    tau_test  = cft_torus.modular_tau(*test_geom)
    s_truth   = tau_ref / tau_test           # complex stretch factor at truth
    r1_true, r2_true = s_truth.real, s_truth.imag

    print("=" * 76)
    print(f"cost_zoo_analytic_test — noise-free cost-surface diagnostic")
    print("=" * 76)
    print(f"  ref  = {ref_geom}  τ_ref  = {tau_ref:.4f}")
    print(f"  test = {test_geom}  τ_test = {tau_test:.4f}")
    print(f"  TRUTH (analytic coords): (r1*, r2*) = "
          f"({r1_true:.4f}, {r2_true:.4f})  [s = τ_ref/τ_test]")
    print(f"  grid: {args.n}×{args.n}  r1∈[{args.r1_range[0]}, "
          f"{args.r1_range[1]}]  r2∈[{args.r2_range[0]}, {args.r2_range[1]}]")
    print(f"  m² = {args.m_sq}  n_cut = {args.n_cut}")

    # Reference data (computed once)
    print("\n  generating ref data ...")
    t0 = time.perf_counter()
    ref_data = cft_torus.make_data_dict(*ref_geom,
                                        m_sq=args.m_sq, n_cut=args.n_cut)
    print(f"    done ({time.perf_counter()-t0:.2f}s, {len(ref_data)} sites)")

    # Cost functions to test
    variants = {
        "test_native":  cost_mod._l2_cost_test_native,   # production
        **{f"zoo_{k}": v for k, v in cost_zoo.ZOO.items()},
    }

    r1s = np.linspace(*args.r1_range, args.n)
    r2s = np.linspace(*args.r2_range, args.n)
    surfaces = {name: np.full((args.n, args.n), np.nan) for name in variants}

    n_total = args.n * args.n
    print(f"\n  scanning {n_total} grid points × {len(variants)} variants ...")
    t0 = time.perf_counter()
    for i, r1 in enumerate(r1s):
        for j, r2 in enumerate(r2s):
            s = complex(r1, r2)
            # Generate test data with overridden ω₂_test_eff
            test_data = cft_torus.make_data_dict(
                *test_geom,
                omega1_override=w1t,
                omega2_override=s * w2t,
                m_sq=args.m_sq, n_cut=args.n_cut)
            # Suppress the per-eval prints from cost functions
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                for name, fn in variants.items():
                    try:
                        c, _, _, _ = fn(ref_data, test_data,
                                        *test_geom, *ref_geom)
                        surfaces[name][i, j] = c if math.isfinite(c) else np.nan
                    except Exception:
                        pass
        elapsed = time.perf_counter() - t0
        eta = elapsed / (i + 1) * (args.n - i - 1)
        print(f"    row {i+1}/{args.n}  elapsed {elapsed:.0f}s  eta {eta:.0f}s")

    # Argmin per surface
    print("\n  argmin per cost surface (analytic coords; truth = "
          f"({r1_true:.3f}, {r2_true:.3f})):")
    summary = dict(
        truth=dict(r1=r1_true, r2=r2_true,
                   tau_ref=str(tau_ref), tau_test=str(tau_test)),
        grid=dict(n=args.n, r1_range=args.r1_range, r2_range=args.r2_range,
                  m_sq=args.m_sq, n_cut=args.n_cut),
        variants={},
    )
    for name, S in surfaces.items():
        if np.all(np.isnan(S)):
            print(f"    {name:>22}: ALL NaN")
            summary["variants"][name] = dict(status="all_nan")
            continue
        idx = np.unravel_index(np.nanargmin(S), S.shape)
        r1_min, r2_min = r1s[idx[0]], r2s[idx[1]]
        d = math.hypot(r1_min - r1_true, r2_min - r2_true)
        c_min = float(S[idx]); c_at_truth = float(np.nan)
        # interpolate cost at truth
        if (args.r1_range[0] <= r1_true <= args.r1_range[1]
                and args.r2_range[0] <= r2_true <= args.r2_range[1]):
            i_t = int(round((r1_true - args.r1_range[0]) /
                            (args.r1_range[1] - args.r1_range[0]) * (args.n-1)))
            j_t = int(round((r2_true - args.r2_range[0]) /
                            (args.r2_range[1] - args.r2_range[0]) * (args.n-1)))
            i_t = max(0, min(args.n-1, i_t))
            j_t = max(0, min(args.n-1, j_t))
            c_at_truth = float(S[i_t, j_t])
        verdict = "✓" if d < 0.15 else ("~" if d < 0.4 else "✗")
        print(f"    {name:>22}: argmin=({r1_min:.3f}, {r2_min:.3f})  "
              f"dist={d:.3f}  cost_min={c_min:.3e}  cost@truth={c_at_truth:.3e}  "
              f"{verdict}")
        summary["variants"][name] = dict(
            argmin=[float(r1_min), float(r2_min)],
            dist_to_truth=d,
            cost_min=c_min,
            cost_at_truth=c_at_truth,
        )

    # Plot 7 heatmaps
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n_var = len(variants)
        cols = 4
        rows = (n_var + cols - 1) // cols
        fig, axs = plt.subplots(rows, cols, figsize=(4.8 * cols, 4.4 * rows))
        axs = np.atleast_2d(axs).ravel()
        R1, R2 = np.meshgrid(r1s, r2s, indexing="ij")
        for ax, (name, S) in zip(axs, surfaces.items()):
            if np.all(np.isnan(S)):
                ax.set_title(f"{name}\n(all NaN)"); ax.axis("off"); continue
            # log-scale heatmap; clip negatives
            S_pos = np.where(np.isfinite(S) & (S > 0), S, np.nan)
            with np.errstate(divide="ignore"):
                logS = np.log10(S_pos)
            pcm = ax.pcolormesh(R1, R2, logS, cmap="viridis", shading="auto")
            plt.colorbar(pcm, ax=ax, label="log₁₀ cost")
            idx = np.unravel_index(np.nanargmin(S), S.shape)
            r1_min, r2_min = r1s[idx[0]], r2s[idx[1]]
            ax.scatter([r1_true], [r2_true], marker="*", s=300, color="red",
                       edgecolor="black", zorder=10, label="TRUTH")
            ax.scatter([r1_min], [r2_min], marker="X", s=120, color="orange",
                       edgecolor="black", zorder=9, label="argmin")
            d = math.hypot(r1_min - r1_true, r2_min - r2_true)
            ax.set_title(f"{name}\nargmin=({r1_min:.2f},{r2_min:.2f})  d={d:.2f}")
            ax.set_xlabel("r1 (= Re of stretch s)")
            ax.set_ylabel("r2 (= Im of stretch s)")
            ax.legend(loc="upper right", fontsize=7)
        for ax in axs[n_var:]:
            ax.axis("off")
        plt.tight_layout()
        png = os.path.join(OUT, "cost_zoo_analytic.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"\n  plot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"  plot failed: {e}")

    with open(os.path.join(OUT, "cost_zoo_analytic_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  summary → results/_native_triage/cost_zoo_analytic_summary.json")


if __name__ == "__main__":
    main()
