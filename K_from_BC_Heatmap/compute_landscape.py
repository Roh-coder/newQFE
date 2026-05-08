#!/usr/bin/env python3
"""
compute_landscape.py — replay cost(r1, r2; L) over the precomputed grid
and extrapolate to the continuum (L → ∞) per (r1, r2).

For each (r1, r2):
  - Load test (untwisted) and ref (twisted) pickles for every available L.
  - Compute the cost on each L with the requested ``cost_mode``.
  - Fit cost(L) vs 1/L (linear by default; quadratic optional) to get
    cost_inf(r1, r2) and an extrapolation uncertainty.
  - Store both finite-L costs and the continuum extrapolation in a single
    .npz, plus a JSON sidecar with the fit per point.

The output is what ``plot_landscape.py`` consumes.

Examples
--------
python compute_landscape.py --tag default --cost-mode huber_log
python compute_landscape.py --tag default --cost-mode test_native --fit quadratic
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

import numpy as np
from scipy.optimize import curve_fit

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import cost as cost_module  # noqa: E402


_COST_MODES = ["l4mean_both_interp", "test_native",
               "spectral", "affine_fit", "huber_log"]


def _point_path(grid_dir: str, r1: float, r2: float) -> str:
    return os.path.join(grid_dir, f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")


def _fit_continuum(Ls, costs, sigmas, mode="linear"):
    """Fit cost(L) → cost_inf, extrapolating to L → ∞.

    Modes
    -----
    linear      cost(L) = a + b/L                          (WLS, 2 params)
    quadratic   cost(L) = a + b/L + c/L²                   (WLS, 3 params, ≥3 pts)
    exponential cost(L) = a + b·exp(−L/ξ)                  (NLS, 3 params, ≥3 pts)

    Returns (cost_inf, sigma_inf, slope, ok).  ``ok`` is False when there
    are too few finite-cost points to fit; in that case we fall back to
    the largest-L value.  For exponential, ``slope`` stores ξ (scale).
    """
    Ls     = np.asarray(Ls,     dtype=float)
    costs  = np.asarray(costs,  dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    mask = np.isfinite(costs) & np.isfinite(sigmas) & (sigmas >= 0.0)
    Ls, costs, sigmas = Ls[mask], costs[mask], sigmas[mask]
    n = len(Ls)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), False
    if n == 1:
        return float(costs[0]), float(sigmas[0]), float("nan"), False

    # ---- exponential: cost = a + b·exp(−L/ξ) ----
    if mode == "exponential" and n >= 3:
        def _exp_model(L, a, b, xi):
            return a + b * np.exp(-L / np.abs(xi))

        # Initial guesses: a from largest L, b from range, ξ from mid-L
        a0 = float(costs[np.argmax(Ls)])
        b0 = float(costs[np.argmin(Ls)] - a0)
        xi0 = float(np.mean(Ls) / 2.0)
        try:
            popt, pcov = curve_fit(
                _exp_model, Ls, costs,
                p0=[a0, b0, xi0],
                sigma=np.maximum(sigmas, 1e-30),
                absolute_sigma=True,
                maxfev=4000,
            )
            cost_inf = float(popt[0])
            sigma_inf = float(np.sqrt(max(pcov[0, 0], 0.0)))
            slope = float(np.abs(popt[2]))   # ξ stored in slope slot
            return cost_inf, sigma_inf, slope, True
        except Exception:  # noqa: BLE001
            pass  # fall through to linear below

    # ---- WLS: linear or quadratic ----
    x = 1.0 / Ls
    w = 1.0 / np.maximum(sigmas, 1e-30) ** 2
    if mode == "quadratic" and n >= 3:
        # cost = a + b*x + c*x²
        X = np.column_stack([np.ones_like(x), x, x ** 2])
    else:
        # cost = a + b*x   (also fallback for failed exponential)
        X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    XtW = X.T @ W
    try:
        cov = np.linalg.inv(XtW @ X)
    except np.linalg.LinAlgError:
        return float(costs[-1]), float(sigmas[-1]), float("nan"), False
    beta = cov @ XtW @ costs
    cost_inf = float(beta[0])
    sigma_inf = float(np.sqrt(max(cov[0, 0], 0.0)))
    slope = float(beta[1]) if len(beta) > 1 else float("nan")
    return cost_inf, sigma_inf, slope, True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--cost-mode", default="huber_log", choices=_COST_MODES)
    ap.add_argument("--cost-power", type=int, default=2)
    ap.add_argument("--fit", default="linear",
                    choices=["linear", "quadratic", "exponential", "all"],
                    help="continuum extrapolation type; 'all' computes all three "
                         "in one pass and stores them with suffixed keys")
    ap.add_argument("--out-name", default=None,
                    help="output basename (default: cost_<mode> or cost_<mode>_all_fits)")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", args.tag)
    mpath = os.path.join(root, "manifest.json")
    if not os.path.exists(mpath):
        print(f"ERROR: manifest missing: {mpath}", file=sys.stderr)
        return 1
    with open(mpath) as f:
        manifest = json.load(f)
    sizes = list(manifest["sizes"])
    twist_frac = manifest["twist_frac"]

    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")

    COST_L  = {L: np.full((n, n), np.nan) for L in sizes}
    SIGMA_L = {L: np.full((n, n), np.nan) for L in sizes}
    BETAC_T = {L: np.full((n, n), np.nan) for L in sizes}
    BETAC_R = {L: np.full((n, n), np.nan) for L in sizes}

    _FIT_MODES = ["linear", "quadratic", "exponential"] if args.fit == "all" else [args.fit]
    COST_INF  = {m: np.full((n, n), np.nan) for m in _FIT_MODES}
    SIGMA_INF = {m: np.full((n, n), np.nan) for m in _FIT_MODES}
    SLOPE_INF = {m: np.full((n, n), np.nan) for m in _FIT_MODES}
    OK_INF    = {m: np.zeros((n, n), dtype=bool) for m in _FIT_MODES}

    fit_records = []  # JSON sidecar
    t0 = time.time()
    n_loaded = n_partial = n_missing = 0

    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            costs_per_L  = []
            sigmas_per_L = []
            Ls_used      = []
            for L in sizes:
                test_pkl = _point_path(
                    os.path.join(root, "grid", f"L{L}", "test"), r1, r2)
                ref_pkl = _point_path(
                    os.path.join(root, "grid", f"L{L}", "ref"), r1, r2)
                if not (os.path.exists(test_pkl) and os.path.exists(ref_pkl)):
                    continue
                with open(test_pkl, "rb") as f:
                    pt = pickle.load(f)
                with open(ref_pkl, "rb") as f:
                    pr = pickle.load(f)
                BETAC_T[L][i, j] = pt["beta_c"]
                BETAC_R[L][i, j] = pr["beta_c"]
                try:
                    c, s, _pd, _pds = cost_module.l2_cost(
                        pr["data"], pt["data"],
                        pt["Lx"], pt["Ly"], pt["Tx"], pt["Ty"],
                        pr["Lx"], pr["Ly"], pr["Tx"], pr["Ty"],
                        cost_mode=args.cost_mode,
                        cost_power=args.cost_power,
                    )
                except Exception as e:  # noqa: BLE001
                    print(f"  [skip] L={L} r=({r1:.2f},{r2:.2f}): {e}")
                    continue
                COST_L[L][i, j]  = c
                SIGMA_L[L][i, j] = s
                costs_per_L.append(c)
                sigmas_per_L.append(s)
                Ls_used.append(L)
            if not costs_per_L:
                n_missing += 1
                continue
            if len(costs_per_L) < len(sizes):
                n_partial += 1
            else:
                n_loaded += 1
            rec = {
                "r1": float(r1), "r2": float(r2),
                "Ls": [int(L) for L in Ls_used],
                "costs": [float(c) for c in costs_per_L],
                "sigmas": [float(s) for s in sigmas_per_L],
            }
            for m in _FIT_MODES:
                cinf, sinf, slope, ok = _fit_continuum(
                    Ls_used, costs_per_L, sigmas_per_L, mode=m)
                COST_INF[m][i, j]  = cinf
                SIGMA_INF[m][i, j] = sinf
                SLOPE_INF[m][i, j] = slope
                OK_INF[m][i, j]    = ok
                rec[m] = {"cost_inf": float(cinf), "sigma_inf": float(sinf),
                          "slope": float(slope), "ok": bool(ok)}
            # backward-compat alias for single-fit runs
            if len(_FIT_MODES) == 1:
                m0 = _FIT_MODES[0]
                rec.update({"cost_inf": rec[m0]["cost_inf"],
                            "sigma_inf": rec[m0]["sigma_inf"],
                            "slope": rec[m0]["slope"],
                            "ok": rec[m0]["ok"]})
            fit_records.append(rec)

    wall = time.time() - t0
    print(f"[replay] full={n_loaded}  partial={n_partial}  missing={n_missing}  "
          f"cost={args.cost_mode}  fit={args.fit}  wall={wall:.1f}s")

    # Build npz payload.  For `--fit all` we store cost_inf_<mode> etc.
    # For a single fit we also keep the legacy un-suffixed keys for
    # backward compatibility with plot_landscape.py.
    out_name = args.out_name or (
        f"cost_{args.cost_mode}_all_fits" if args.fit == "all"
        else f"cost_{args.cost_mode}"
    )
    npz_payload = dict(
        rs=rs, R1=R1, R2=R2,
        sizes=np.array(sizes, dtype=int),
        twist_frac=twist_frac,
        **{f"cost_L{L}": COST_L[L] for L in sizes},
        **{f"sigma_L{L}": SIGMA_L[L] for L in sizes},
        **{f"betac_test_L{L}": BETAC_T[L] for L in sizes},
        **{f"betac_ref_L{L}": BETAC_R[L] for L in sizes},
    )
    for m in _FIT_MODES:
        npz_payload[f"cost_inf_{m}"]  = COST_INF[m]
        npz_payload[f"sigma_inf_{m}"] = SIGMA_INF[m]
        npz_payload[f"slope_inf_{m}"] = SLOPE_INF[m]
        npz_payload[f"ok_inf_{m}"]    = OK_INF[m]
    if len(_FIT_MODES) == 1:
        # legacy un-suffixed aliases
        m0 = _FIT_MODES[0]
        npz_payload.update(cost_inf=COST_INF[m0], sigma_inf=SIGMA_INF[m0],
                           slope_inf=SLOPE_INF[m0], ok_inf=OK_INF[m0])

    npz = os.path.join(root, f"{out_name}.npz")
    np.savez(npz, **npz_payload)
    json_path = os.path.join(root, f"{out_name}_fits.json")
    with open(json_path, "w") as f:
        json.dump({
            "tag": args.tag, "cost_mode": args.cost_mode,
            "cost_power": args.cost_power, "fit": args.fit,
            "fit_modes": _FIT_MODES,
            "sizes": sizes, "twist_frac": twist_frac,
            "points": fit_records,
        }, f, indent=2)
    print(f"[done] {npz}")
    print(f"[done] {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
