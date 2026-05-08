#!/usr/bin/env python3
"""
zoo_ising_tau_grid.py — 2-D scan over modular parameter τ.

For a fixed REFERENCE torus geometry we hold ω₁_test = ω₁_ref and let
ω₂_test sweep over a 2-D grid in the complex plane:

    τ_test_eff = ω₂_test / ω₁_ref

so the scan IS a scan over τ.  Truth is at τ_test_eff = τ_ref.

For each grid cell we build the analytic Ising-CFT data dict (no MC
noise) and compute every cost variant against the ref data.  We then
analyze the resulting 2-D cost landscape:

  * dist (min → truth) in (Re τ, Im τ) plane
  * |basin Hessian| at truth (steepness)
  * basin aspect ratio (anisotropy of the cost landscape)
  * # spurious local minima
  * dynamic range = max(cost)/min(cost) on the grid

and emit a 6-panel heatmap with truth marked.

Usage:  python zoo_ising_tau_grid.py [--n N]
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


def _local_minima_count(g: np.ndarray) -> int:
    n = 0
    H, W = g.shape
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            v = g[i, j]
            if not np.isfinite(v): continue
            nbrs = [g[i-1, j], g[i+1, j], g[i, j-1], g[i, j+1],
                    g[i-1, j-1], g[i-1, j+1], g[i+1, j-1], g[i+1, j+1]]
            nbrs = [x for x in nbrs if np.isfinite(x)]
            if len(nbrs) >= 4 and v < min(nbrs):
                n += 1
    return n


def _hessian_at(g: np.ndarray, ix: int, iy: int,
                dx: float, dy: float):
    """Return (H_xx, H_yy, H_xy) of log(g) at grid point (ix, iy)."""
    H, W = g.shape
    if not (1 <= ix < H - 1 and 1 <= iy < W - 1):
        return float("nan"), float("nan"), float("nan")
    vc = g[ix, iy]; vL = g[ix-1, iy]; vR = g[ix+1, iy]
    vD = g[ix, iy-1]; vU = g[ix, iy+1]
    vDL = g[ix-1, iy-1]; vDR = g[ix+1, iy-1]
    vUL = g[ix-1, iy+1]; vUR = g[ix+1, iy+1]
    vals = [vc, vL, vR, vD, vU, vDL, vDR, vUL, vUR]
    if any(not (np.isfinite(v) and v > 0) for v in vals):
        return float("nan"), float("nan"), float("nan")
    lc = math.log(vc)
    Hxx = (math.log(vL) + math.log(vR) - 2 * lc) / (dx * dx)
    Hyy = (math.log(vD) + math.log(vU) - 2 * lc) / (dy * dy)
    Hxy = (math.log(vUR) + math.log(vDL)
           - math.log(vUL) - math.log(vDR)) / (4 * dx * dy)
    return Hxx, Hyy, Hxy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",   type=int, default=25)
    ap.add_argument("--ref", type=int, nargs=4, default=[13, 16, -3, 3])
    ap.add_argument("--test",type=int, nargs=4, default=[16, 16,  0,  0])
    # τ-window centred on τ_ref
    ap.add_argument("--span-re", type=float, default=1.2)
    ap.add_argument("--span-im", type=float, default=1.0)
    args = ap.parse_args()

    Rx, Ry, RTx, RTy = args.ref
    Tx, Ty, TTx, TTy = args.test
    print(f"  ref geom : ({Rx}, {Ry}, {RTx}, {RTy})")
    print(f"  test geom: ({Tx}, {Ty}, {TTx}, {TTy})")

    w1_r, w2_r = cft.torus_periods(Rx, Ry, RTx, RTy)
    tau_r = w2_r / w1_r
    if tau_r.imag < 0: tau_r = complex(tau_r.real, -tau_r.imag)
    print(f"  τ_ref = {tau_r:.4f}")

    # Grid: τ = τ_re + i τ_im around (Re τ_ref, Im τ_ref).
    re_lo = tau_r.real - args.span_re / 2
    re_hi = tau_r.real + args.span_re / 2
    im_lo = max(0.05, tau_r.imag - args.span_im / 2)
    im_hi = tau_r.imag + args.span_im / 2
    re_arr = np.linspace(re_lo, re_hi, args.n)
    im_arr = np.linspace(im_lo, im_hi, args.n)
    print(f"  scan: Re τ ∈ [{re_lo:.3f}, {re_hi:.3f}],  "
          f"Im τ ∈ [{im_lo:.3f}, {im_hi:.3f}],  {args.n}×{args.n}")

    # Ref Ising data once.
    print("  building ref Ising correlator…")
    ref_data = cft.make_ising_data_dict(Rx, Ry, RTx, RTy)

    variants = {
        "production":  cost_mod._l2_cost_test_native,
        **{f"zoo_{k}": v for k, v in cost_zoo.ZOO.items()},
    }
    print(f"  variants: {list(variants.keys())}")

    grids = {n: np.full((args.n, args.n), np.nan) for n in variants}
    nan_count = {n: 0 for n in variants}

    real_stdout = sys.stdout
    null = io.StringIO()
    t0 = time.perf_counter()
    n_total = args.n * args.n
    cell = 0
    for i, re in enumerate(re_arr):
        for j, im in enumerate(im_arr):
            cell += 1
            if cell % max(1, n_total // 20) == 0:
                real_stdout.write(f"    cell {cell}/{n_total}  "
                                  f"({100 * cell / n_total:.0f}%)\n")
                real_stdout.flush()
            # Set τ_test = re + i·im → ω₂_test = τ · ω₁_ref
            tau = complex(re, im)
            try:
                test_data = cft.make_ising_data_dict(
                    Tx, Ty, TTx, TTy,
                    omega1_override=w1_r,
                    omega2_override=tau * w1_r)
            except Exception:
                continue
            for name, fn in variants.items():
                sys.stdout = null
                try:
                    c, *_ = fn(ref_data, test_data,
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
                    sys.stdout = real_stdout
    wall = time.perf_counter() - t0
    print(f"  scan done in {wall:.1f}s  ({wall / n_total * 1000:.1f} ms/cell)")

    # Index of truth on the grid (nearest cell).
    i_truth = int(np.argmin(np.abs(re_arr - tau_r.real)))
    j_truth = int(np.argmin(np.abs(im_arr - tau_r.imag)))
    truth_re = re_arr[i_truth]; truth_im = im_arr[j_truth]
    print(f"\n  truth grid cell: (re={truth_re:.3f}, im={truth_im:.3f})  "
          f"index ({i_truth}, {j_truth})")

    dre = re_arr[1] - re_arr[0]
    dim = im_arr[1] - im_arr[0]

    print()
    hdr = (f"  {'variant':<14}  {'min_τ':>16}  {'dist→truth':>10}  "
           f"{'cost@truth':>11}  "
           f"{'|H|':>9}  {'aniso':>6}  {'#minima':>8}  "
           f"{'dyn':>9}  {'nan':>5}")
    print(hdr)
    print("  " + "-" * len(hdr))

    summary = dict(ref=args.ref, test=args.test,
                   tau_r=[tau_r.real, tau_r.imag],
                   re=re_arr.tolist(), im=im_arr.tolist(),
                   variants={})

    rows = []
    for name in variants:
        g = grids[name]
        if np.all(np.isnan(g)):
            print(f"  {name:<14}  ALL NaN")
            summary["variants"][name] = dict(all_nan=True); continue

        i_min, j_min = np.unravel_index(np.nanargmin(g), g.shape)
        re_min = float(re_arr[i_min]); im_min = float(im_arr[j_min])
        c_min = float(g[i_min, j_min])
        c_at_truth = float(g[i_truth, j_truth])
        dist = math.hypot(re_min - tau_r.real, im_min - tau_r.imag)

        # Hessian of log(g) at truth
        Hxx, Hyy, Hxy = _hessian_at(g, i_truth, j_truth, dre, dim)
        if all(map(np.isfinite, [Hxx, Hyy, Hxy])):
            # eigvals of [[Hxx, Hxy],[Hxy, Hyy]]
            tr = Hxx + Hyy
            det = Hxx * Hyy - Hxy * Hxy
            disc = max(tr * tr / 4 - det, 0.0)
            ev_max = tr / 2 + math.sqrt(disc)
            ev_min = tr / 2 - math.sqrt(disc)
            H_norm = math.sqrt(Hxx ** 2 + Hyy ** 2 + 2 * Hxy ** 2)
            aniso = (abs(ev_max) / abs(ev_min)
                     if ev_min and ev_max * ev_min > 0 else float("nan"))
        else:
            H_norm = float("nan"); aniso = float("nan")

        n_minima = _local_minima_count(g)
        with np.errstate(invalid="ignore", divide="ignore"):
            finite = g[np.isfinite(g) & (g > 0)]
            dyn = float(np.max(finite) / np.min(finite)) if len(finite) else float("nan")

        print(f"  {name:<14}  ({re_min:+.2f},{im_min:+.2f})  "
              f"{dist:>10.3f}  {c_at_truth:>11.3e}  "
              f"{H_norm:>9.2e}  {aniso:>6.2f}  {n_minima:>8d}  "
              f"{dyn:>9.2e}  {nan_count[name]:>5d}")

        summary["variants"][name] = dict(
            min_re=re_min, min_im=im_min, c_min=c_min,
            c_at_truth=c_at_truth,
            dist_to_truth=dist,
            hessian_norm=H_norm, basin_aniso=aniso,
            n_local_minima=n_minima,
            dynamic_range=dyn,
            nan_count=nan_count[name],
            grid=g.tolist())
        rows.append((name, dist, H_norm, n_minima, c_at_truth / c_min
                     if c_min > 0 else float("nan")))

    # Composite ranking (lower dist + higher curvature + fewer spurious minima)
    print()
    print("  RANKING (composite: dist↓, |H|↑, #minima↓):")
    rows_ok = [r for r in rows if all(np.isfinite(x) for x in r[1:])]
    rows_ok.sort(key=lambda r: (r[1],   # dist to truth
                                -r[2],  # -|H| (higher is better)
                                r[3]))  # # minima
    for i, (n, d, H_, nm, ratio) in enumerate(rows_ok, 1):
        print(f"    {i}. {n:<14}  dist={d:.3f}  |H|={H_:.2e}  "
              f"#min={nm}  cost(truth)/cost(min)={ratio:.3f}")
    bad = [r for r in rows if not all(np.isfinite(x) for x in r[1:])]
    for n, d, H_, nm, ratio in bad:
        print(f"    --. {n:<14}  dist={d:.3f}  (Hessian/aniso ill-defined)")

    # Save summary
    with open(os.path.join(OUT, "zoo_ising_tau_grid_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    # Plot 6-panel heatmap of log10(cost).
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ncols = 3
        nrows = (len(variants) + ncols - 1) // ncols
        fig, axs = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.4 * nrows),
                                squeeze=False)
        for k, name in enumerate(variants):
            r, c = k // ncols, k % ncols
            ax = axs[r, c]
            g = grids[name]
            if np.all(np.isnan(g)):
                ax.set_title(f"{name}: ALL NaN"); ax.axis("off"); continue
            with np.errstate(invalid="ignore", divide="ignore"):
                logg = np.where(g > 0, np.log10(g), np.nan)
            # vmin/vmax from valid percentiles
            valid = logg[np.isfinite(logg)]
            if len(valid) > 0:
                vmin, vmax = np.percentile(valid, [2, 98])
            else:
                vmin, vmax = -1, 1
            im = ax.imshow(logg.T, origin="lower",
                           extent=[re_arr[0], re_arr[-1],
                                   im_arr[0], im_arr[-1]],
                           aspect="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
            plt.colorbar(im, ax=ax, label="log₁₀ cost")
            # truth & min markers
            ax.scatter([tau_r.real], [tau_r.imag], marker="*",
                       s=300, color="red", edgecolor="white", lw=1.5,
                       zorder=10, label="truth")
            i_min, j_min = np.unravel_index(np.nanargmin(g), g.shape)
            ax.scatter([re_arr[i_min]], [im_arr[j_min]], marker="X",
                       s=140, color="orange", edgecolor="black",
                       zorder=9, label="cost min")
            # contour lines
            try:
                CS = ax.contour(re_arr, im_arr, logg.T,
                                levels=8, colors="white",
                                linewidths=0.6, alpha=0.6)
            except Exception:
                pass
            ax.set_xlabel("Re τ"); ax.set_ylabel("Im τ")
            d = summary["variants"][name].get("dist_to_truth", float("nan"))
            H_ = summary["variants"][name].get("hessian_norm", float("nan"))
            nm = summary["variants"][name].get("n_local_minima", -1)
            ax.set_title(f"{name}\ndist={d:.2f}  |H|={H_:.1e}  "
                         f"#min={nm}", fontsize=9)
            ax.legend(loc="upper right", fontsize=7)
        for k in range(len(variants), nrows * ncols):
            axs[k // ncols, k % ncols].axis("off")
        plt.tight_layout()
        png = os.path.join(OUT, "zoo_ising_tau_grid.png")
        plt.savefig(png, dpi=130); plt.close(fig)
        print(f"\nplot → {os.path.relpath(png, HERE)}")
    except Exception as e:
        print(f"plot failed: {e}")


if __name__ == "__main__":
    main()
