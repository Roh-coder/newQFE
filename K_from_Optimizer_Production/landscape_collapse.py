"""
landscape_collapse.py — Brower-Owen multi-axis curve-collapse cost.

Reference: "Ising Model on the Affine Plane" (Brower & Owen, Sec 5.2,
Fig. 13). At critical couplings sinh(2K_i)=ell*_i/ell_i, the spin-spin
correlator measured along all 6 lattice axes collapses onto a single
universal curve when distances are rescaled by the appropriate physical
edge length. Off-criticality the curves diverge.

Cost: extract G(r) along the 6 axes (X, Y, Z and the three opposite
diagonals XY, XZ, YZ are equivalent up to symmetry on a periodic
parallelogram, so we use the 3 independent cycle directions u, v, w
already in the precompute), rescale physical distance by the geometric
edge length, then measure the residual scatter of the 3 curves around
their pointwise mean / fit.

Cost kernels (no reference needed!):
  1. collapse_var      pointwise variance across rescaled axes
  2. collapse_log      same but on log G
  3. collapse_band     median-absolute deviation around pointwise median
  4. collapse_pairs    sum of pairwise (l2)^2 between rescaled curves
  5. collapse_smooth   collapse + smoothness penalty (curvature of pointwise mean)

Also runs an "anti-collapse cross-check" against the paper's prediction:
  6. axis_aniso        var of asymptotic decay slopes across axes (small
                       at criticality, large off-criticality)
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

from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _SQRT3_2)  # noqa: E402

_TRUTH_456 = (5.0652, 7.7429)


# --------------------- per-axis sample extraction --------------------- #
def _axis_samples(test_data, test_Lx, test_Ly, test_Tx, test_Ty, copies=2):
    """Return list of 3 dicts: per cycle axis, (r_phys, G, e_G).

    r_phys = t * |cycle_vector| = physical distance along that axis on
    the test lattice.
    """
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)
    out = []
    for (tdm, tdn) in test_paths:
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 3:
            out.append(None); continue
        t_arr = np.asarray([k / N for k in ks], dtype=float)
        tex, tey = tdm + 0.5 * tdn, _SQRT3_2 * tdn
        p_len = float(np.hypot(tex, tey))
        Gt = np.full(N, np.nan); et = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       test_Lx, test_Ly, test_Tx, test_Ty,
                                       copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
                et[i] = abs(entry["conn_err"])
        r_phys = t_arr * p_len
        out.append(dict(r=r_phys, G=Gt, e=et, p_len=p_len, N=N))
    return out


# --------------------- collapse helpers ------------------------------ #
def _gather_rescaled(axes, eps=1e-12):
    """Rescale each axis's r by 1 / max(G)... no wait, scale by p_len.

    We want a SHARED parameter across axes such that G_axis(r/scale_axis)
    collapses. Brower-Owen use scale_axis = ell_axis (the physical edge
    length). On our test (16,16,0,0) the cycle lengths are
    |u|=|v|=16, |w|=22.6. After rescaling, all axes use t in [0,1).

    Returns dict per axis: t (in [0,1]), G(t), e(t), each filtered for
    finite & positive G.
    """
    out = []
    for ax in axes:
        if ax is None: out.append(None); continue
        r = ax["r"]; G = ax["G"]; e = ax["e"]
        m = np.isfinite(G) & np.isfinite(r) & (G > eps)
        if m.sum() < 3:
            out.append(None); continue
        t = r[m] / ax["p_len"]   # fractional position along cycle
        Gv = G[m]; ev = e[m]
        order = np.argsort(t)
        out.append(dict(t=t[order], G=Gv[order], e=ev[order]))
    return out


def _common_grid(rescaled, n_grid=24, lo=0.05, hi=0.45, eps=1e-12):
    """Interpolate each axis curve onto a common t-grid (log G vs t)."""
    tg = np.linspace(lo, hi, n_grid)
    cols = []
    for d in rescaled:
        if d is None: cols.append(None); continue
        if d["t"][-1] < hi or d["t"][0] > lo:
            cols.append(None); continue
        # interp log G at tg
        lg = np.interp(tg, d["t"], np.log(d["G"]))
        cols.append(lg)
    valid = [c for c in cols if c is not None]
    if len(valid) < 2:
        return None, None
    M = np.vstack(valid)  # (n_axes, n_grid)
    return tg, M


# --------------------- cost kernels ---------------------------------- #
def k_collapse_var(rescaled):
    tg, M = _common_grid(rescaled)
    if M is None: return np.nan
    # variance of log G across axes at each t, summed
    return float(np.mean(np.var(M, axis=0)))


def k_collapse_log(rescaled):
    """Same kernel but normalised by mean magnitude (stabilises tilt)."""
    tg, M = _common_grid(rescaled)
    if M is None: return np.nan
    mu = M.mean(axis=0)
    dev = M - mu[None, :]
    return float(np.mean(dev * dev))


def k_collapse_band(rescaled):
    """Median-abs-deviation around pointwise median of log G."""
    tg, M = _common_grid(rescaled)
    if M is None: return np.nan
    med = np.median(M, axis=0)
    return float(np.mean(np.abs(M - med[None, :])))


def k_collapse_pairs(rescaled):
    """Sum of pairwise (log G_i - log G_j)^2 across axes."""
    tg, M = _common_grid(rescaled)
    if M is None: return np.nan
    n = M.shape[0]
    s = 0.0; cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = M[i] - M[j]
            s += float(np.mean(d * d)); cnt += 1
    return s / max(cnt, 1)


def k_collapse_smooth(rescaled, lam=0.2):
    """Collapse + curvature penalty on the pointwise mean (rules out
    degenerate flat solutions where all G≈const)."""
    tg, M = _common_grid(rescaled)
    if M is None: return np.nan
    mu = M.mean(axis=0)
    dev = M - mu[None, :]
    var_term = float(np.mean(dev * dev))
    # second derivative of mu along t
    d2 = np.diff(mu, n=2)
    curv = float(np.mean(d2 * d2))
    # we WANT some curvature (decay), so penalise FLAT means by adding
    # 1/curv? Easier: simply require curvature above threshold by adding
    # lambda * exp(-curv/scale). For now just cost = var.
    return var_term + lam * np.exp(-curv * 1e3)


def k_axis_aniso(rescaled, t_lo=0.10, t_hi=0.40):
    """Variance of fitted decay-slope across axes in window [t_lo, t_hi].

    Slope = d(log G)/dt fitted by linregress per axis. At criticality
    all axes have the same slope (= conformal scaling on R^2). Off
    criticality, slopes differ.
    """
    slopes = []
    for d in rescaled:
        if d is None: continue
        m = (d["t"] >= t_lo) & (d["t"] <= t_hi)
        if m.sum() < 3: continue
        t = d["t"][m]; lg = np.log(d["G"][m])
        s = np.polyfit(t, lg, 1)[0]
        slopes.append(s)
    if len(slopes) < 2:
        return np.nan
    return float(np.var(slopes))


KERNELS = [
    ("collapse_var",     k_collapse_var),
    ("collapse_log",     k_collapse_log),
    ("collapse_band",    k_collapse_band),
    ("collapse_pairs",   k_collapse_pairs),
    ("collapse_smooth",  k_collapse_smooth),
    ("axis_aniso",       k_axis_aniso),
]


# --------------------- driver ---------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--out",       default=None)
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    print(f"[load+extract] {n*n} grid pts ...")
    t0 = time.time()
    samples = {}
    for r1 in rs:
        for r2 in rs:
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            axes = _axis_samples(pt["test_data"], Lx, Ly, Tx, Ty)
            samples[(float(r1), float(r2))] = _gather_rescaled(axes)
    print(f"[load+extract] {len(samples)} pts in {time.time()-t0:.1f}s")

    results = {}
    print("[eval] kernels:")
    for name, fn in KERNELS:
        Z = np.full((n, n), np.nan)
        t0 = time.time()
        for i, r1 in enumerate(rs):
            for j, r2 in enumerate(rs):
                d = samples.get((float(r1), float(r2)))
                if d is None: continue
                try:
                    Z[i, j] = fn(d)
                except Exception:
                    pass
        if np.all(np.isnan(Z)):
            print(f"  {name:18s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        truth_i = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        truth_j = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        results[name] = dict(Z=Z, argmin=(rs[ij[0]], rs[ij[1]]),
                             d_truth=d_truth,
                             cost_truth=float(Z[truth_i, truth_j]),
                             cost_argmin=float(Z[ij]))
        print(f"  {name:18s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:5.2f}  cost(truth)={Z[truth_i,truth_j]:.3e}  "
              f"cost(min)={Z[ij]:.3e}  wall={time.time()-t0:.1f}s")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    K = len(results)
    cols = 3
    rows = (K + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 3.8 * rows))
    axes = axes.flatten()
    for ax, (name, info) in zip(axes, results.items()):
        Z = np.log10(np.maximum(info["Z"], 1e-15))
        vmin, vmax = np.nanpercentile(Z, [2, 98])
        im = ax.pcolormesh(R1, R2, Z, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=11, mec="k", mew=1)
        am = info["argmin"]
        ax.plot(am[0], am[1], "w*", ms=12, mec="k", mew=1)
        ax.set_title(f"{name}\nargmin=({am[0]:.1f},{am[1]:.1f}) "
                     f"d={info['d_truth']:.1f}", fontsize=10)
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_aspect("equal")
    for ax in axes[K:]: ax.axis("off")
    fig.suptitle(f"Brower-Owen multi-axis collapse cost  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"truth=(5.07,7.74)  [reference-FREE]", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or os.path.join(root, "collapse.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n[done] figure → {out}")

    npz_out = (os.path.splitext(out)[0]) + ".npz"
    np.savez(npz_out, R1=R1, R2=R2, rs=rs,
             **{k: v["Z"] for k, v in results.items()})
    print(f"[done] data   → {npz_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
