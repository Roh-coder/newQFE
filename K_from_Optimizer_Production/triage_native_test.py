#!/usr/bin/env python3
"""
triage_native_test.py — compare CURRENT test_native cost-mode behaviour
against the PROPOSED reliability fixes from the workflow triage.

Read-only: imports cost.py / mc_engine.py from this directory but does not
modify them.  Reuses cached two_point_all_to_all.dat replicates from
results/_betac_stability_test/{equil_12, aniso_r2_12, aniso_r5_12} (10
RW + 10 GC runs per case), all at L=12, T=0.

Tests
-----
T1  Determinism vs key order:
      For one (ref,test) pair, shuffle dict-insertion order of test_data
      and recompute the current cost.  Report max cost-spread over 8
      shuffles relative to printed cost.

T2  Convex-hull dropouts (copies=2 vs copies=4):
      For a contrived twisted geometry (Lx=12,Ly=12,Tx=3,Ty=-3) on the
      same equil_12 data, count how many of 400 sample points the legacy
      'l4mean_both_interp' path drops as NaN.

T3  Silent missing-direction:
      Run test_native cost with a contrived test geometry where one of
      the three boundary directions has gcd=1 (only one lattice site).
      Compare current behaviour (silent C_d=0) against a "valid-dir mean"
      that drops collapsed directions.

T4  Reference-error interpolation bias:
      For one (ref,test) pair, recompute test_native σ_cost using
      (a) current LinearNDInterpolator on conn_err,
      (b) "nearest-vertex max" — for each query point take the max
          conn_err across the surrounding 3 Delaunay vertices.
      Report (b)/(a) ratio per direction.

T5  Discrimination (Z) — current vs proposed:
      "Proposed" cost = inverse-variance weighted mean over directions
      using honest per-direction σ from (T4-b).  Compute pairwise costs
      on the 10×10 RW-replicate grid for null (equil_12 vs equil_12) and
      signal (equil_12 vs aniso_r2_12, equil_12 vs aniso_r5_12).
      Z = (mean_signal - mean_null) / std_null .

Outputs
-------
  results/_native_triage/triage_summary.json
  results/_native_triage/triage.log     (full stdout)
"""
from __future__ import annotations

import io
import contextlib
import glob
import json
import math
import os
import random
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cost as cost_mod                 # noqa: E402
import mc_engine                        # noqa: E402

OUT_DIR = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT_DIR, exist_ok=True)

STABILITY_ROOT = os.path.join(
    HERE, "results", "_betac_stability_test")

CASES = {
    "equil":    "equil_12",
    "aniso_r2": "aniso_r2_12",
    "aniso_r5": "aniso_r5_12",
}
LXY = (12, 12, 0, 0)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _replicate_files(case_dir: str, kind: str = "rw") -> List[str]:
    pat = os.path.join(STABILITY_ROOT, case_dir,
                       f"{kind}_t*", "*", "two_point_all_to_all.dat")
    return sorted(glob.glob(pat))


def load_replicates() -> Dict[str, List[dict]]:
    out = {}
    for tag, sub in CASES.items():
        files = _replicate_files(sub, kind="rw")
        if not files:
            files = _replicate_files(sub, kind="3pass")
        out[tag] = [mc_engine.load_all_to_all(p) for p in files]
        print(f"[load] {tag}: {len(files)} replicates  "
              f"(first={os.path.basename(os.path.dirname(files[0]))})")
    return out


# ---------------------------------------------------------------------------
# Quiet wrapper around cost_mod.l2_cost (it prints unconditionally)
# ---------------------------------------------------------------------------

def quiet_cost(ref_data, test_data, geom_test, geom_ref, **kw):
    Lx, Ly, Tx, Ty = geom_test
    rLx, rLy, rTx, rTy = geom_ref
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        c, s, pd, ps = cost_mod.l2_cost(
            ref_data, test_data,
            Lx, Ly, Tx, Ty,
            rLx, rLy, rTx, rTy,
            **kw)
    return c, s, pd, ps


# ===========================================================================
# T1 — determinism vs dict key order
# ===========================================================================

def test_T1_determinism(ref, test) -> dict:
    print("\n=== T1: determinism vs key order ===")
    base, _, _, _ = quiet_cost(ref, test, LXY, LXY,
                               cost_mode="test_native", cost_power=2)
    costs = [base]
    rng = random.Random(0xC0FFEE)
    keys = list(test.keys())
    for trial in range(8):
        rng.shuffle(keys)
        test_shuf = {k: test[k] for k in keys}
        c, _, _, _ = quiet_cost(ref, test_shuf, LXY, LXY,
                                cost_mode="test_native", cost_power=2)
        costs.append(c)
    spread = max(costs) - min(costs)
    rel = spread / abs(base) if base else 0.0
    print(f"  base cost          = {base:.6e}")
    print(f"  min..max over 8 sh = {min(costs):.6e} .. {max(costs):.6e}")
    print(f"  abs spread         = {spread:.3e}")
    print(f"  rel spread         = {rel:.3e}")
    return dict(base=base, min=min(costs), max=max(costs),
                spread=spread, rel_spread=rel)


# ===========================================================================
# T2 — convex hull dropouts in legacy interp path
# ===========================================================================

def _hull_dropouts(ref_data, test_data, geom_test, geom_ref,
                   copies: int, n_samples: int = 400):
    Lx, Ly, Tx, Ty = geom_test
    rLx, rLy, rTx, rTy = geom_ref
    iref  = cost_mod._tile_interp(ref_data,  rLx, rLy, rTx, rTy, "conn",  copies)
    itest = cost_mod._tile_interp(test_data, Lx, Ly, Tx, Ty, "conn",  copies)
    ref_paths  = cost_mod.boundary_paths(rLx, rLy, rTx, rTy)
    test_paths = cost_mod.boundary_paths(Lx, Ly, Tx, Ty)
    t = np.linspace(0.0, 1.0, n_samples, endpoint=False)
    drops_ref, drops_test = 0, 0
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        rex, rey = rdm + 0.5*rdn, math.sqrt(3)/2*rdn
        tex, tey = tdm + 0.5*tdn, math.sqrt(3)/2*tdn
        Gr = iref(np.column_stack([t*rex, t*rey]))
        Gt = itest(np.column_stack([t*tex, t*tey]))
        drops_ref  += int(np.isnan(Gr).sum())
        drops_test += int(np.isnan(Gt).sum())
    return drops_ref, drops_test


def test_T2_hull_dropouts(ref, test) -> dict:
    print("\n=== T2: convex-hull dropouts (copies=2 vs 4) ===")
    twisted = (12, 12, 3, -3)
    res = {}
    for c in (2, 4):
        dr, dt = _hull_dropouts(ref, test, twisted, twisted, copies=c)
        print(f"  copies={c}  ref-NaN={dr}/1200   test-NaN={dt}/1200")
        res[f"copies_{c}"] = dict(ref_drops=dr, test_drops=dt)
    return res


# ===========================================================================
# T3 — silent missing-direction
# ===========================================================================

def test_T3_silent_zero(ref, test) -> dict:
    print("\n=== T3: silent missing-direction (gcd=1 collapse) ===")
    # Geometry where direction (Lx, Ty) = (12, 1) has gcd=1 → 1 site → collapses.
    Lx, Ly, Tx, Ty = 12, 12, 0, 1
    paths = cost_mod.boundary_paths(Lx, Ly, Tx, Ty)
    gcds = [math.gcd(abs(dm), abs(dn)) for dm, dn in paths]
    print(f"  contrived test geom (Lx,Ly,Tx,Ty)={Lx,Ly,Tx,Ty}  paths={paths}  gcds={gcds}")
    # We use the equil_12 reference at its native geometry; the "test" data
    # is loaded but at LXY — so the look-up at (12,1) etc. will return None
    # for at least the gcd=1 direction and still find sites along the gcd=12
    # direction.  The point is: does the cost see any direction collapse?
    c, s, pd, ps = quiet_cost(ref, test, (Lx, Ly, Tx, Ty), LXY,
                              cost_mode="test_native", cost_power=2)
    n_zero = sum(1 for v in pd if v == 0.0 and ps[ pd.index(v) ] == 0.0)
    print(f"  per_dir = {[f'{v:.3e}' for v in pd]}")
    print(f"  per_sig = {[f'{v:.3e}' for v in ps]}")
    print(f"  cost    = {c:.4e} ± {s:.2e}")
    print(f"  collapsed dirs counted = {n_zero}/3   "
          "(current code averages them as 0 — biased low)")
    # Proposed behaviour: drop collapsed dirs from the mean.
    valid = [(v, sv) for v, sv in zip(pd, ps) if not (v == 0.0 and sv == 0.0)]
    if valid:
        c_prop = sum(v for v, _ in valid) / len(valid)
        s_prop = math.sqrt(sum(sv**2 for _, sv in valid)) / len(valid)
    else:
        c_prop, s_prop = float("nan"), float("nan")
    print(f"  proposed (drop-collapsed): cost={c_prop:.4e} ± {s_prop:.2e}  "
          f"(N_valid={len(valid)})")
    bias = (c_prop - c) / c if c else float("nan")
    print(f"  current/proposed cost ratio = {c / c_prop:.3f}  "
          f"(current is {bias*100:.1f}% lower)")
    return dict(per_dir=pd, per_sig=ps,
                cost_current=c, sigma_current=s,
                cost_proposed=c_prop, sigma_proposed=s_prop,
                n_collapsed=n_zero,
                bias_frac=(c - c_prop) / c_prop if c_prop else None)


# ===========================================================================
# T4 — reference-error interpolation bias
# ===========================================================================

def _nearest_vertex_max_err(ref_data, geom, query_xy, copies=2):
    """For each (x,y) query, return the max conn_err of the surrounding
    Delaunay simplex vertices in the tiled point cloud."""
    rLx, rLy, rTx, rTy = geom
    keys = list(ref_data.keys())
    m = np.array([k[0] for k in keys], float)
    n = np.array([k[1] for k in keys], float)
    e_arr = np.array([ref_data[k]["conn_err"] for k in keys], float)
    x0, y0 = m + 0.5*n, math.sqrt(3)/2 * n
    vx, vy =  rLx + 0.5*rTy,            math.sqrt(3)/2 * rTy
    ux, uy =  rTx - 0.5*rLy,           -math.sqrt(3)/2 * rLy
    xs, ys, es = [x0], [y0], [e_arr]
    for a in range(-copies, copies+1):
        for b in range(-copies, copies+1):
            if a == 0 and b == 0:
                continue
            xs.append(x0 + a*vx + b*ux)
            ys.append(y0 + a*vy + b*uy)
            es.append(e_arr)
    pts = np.column_stack([np.concatenate(xs), np.concatenate(ys)])
    es  = np.concatenate(es)
    _, uid = np.unique(np.round(pts, 8), axis=0, return_index=True)
    pts = pts[uid]; es = es[uid]
    tri = Delaunay(pts)
    out = np.full(len(query_xy), np.nan)
    simp = tri.find_simplex(query_xy)
    for i, s in enumerate(simp):
        if s < 0:
            continue
        verts = tri.simplices[s]
        out[i] = es[verts].max()
    return out


def test_T4_err_interp_bias(ref, test) -> dict:
    print("\n=== T4: ref-σ interpolation vs nearest-vertex-max ===")
    Lx, Ly, Tx, Ty = LXY
    iref_err = cost_mod._tile_interp(ref, Lx, Ly, Tx, Ty, "conn_err", copies=2)
    ref_paths  = cost_mod.boundary_paths(*LXY)
    test_paths = cost_mod.boundary_paths(*LXY)
    res_per_dir = []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = cost_mod._direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            continue
        N = len(ks)
        rex, rey = rdm + 0.5*rdn, math.sqrt(3)/2*rdn
        t = ks / N
        pts = np.column_stack([t*rex, t*rey])
        e_interp = np.abs(np.asarray(iref_err(pts), float))
        e_max    = _nearest_vertex_max_err(ref, LXY, pts, copies=2)
        ratio = np.nanmean(e_max / e_interp) if np.any(np.isfinite(e_interp)) else float("nan")
        res_per_dir.append(dict(
            dir=(tdm, tdn),
            mean_interp=float(np.nanmean(e_interp)),
            mean_vertex_max=float(np.nanmean(e_max)),
            ratio_max_over_interp=float(ratio),
        ))
        print(f"  dir={tdm,tdn}: mean σ_interp={np.nanmean(e_interp):.3e}  "
              f"mean σ_vertex_max={np.nanmean(e_max):.3e}  "
              f"ratio={ratio:.3f}")
    return dict(per_dir=res_per_dir)


# ===========================================================================
# T5 — discrimination Z (current vs proposed)
# ===========================================================================

def _proposed_native(ref, test, geom):
    """Proposed test-native cost: drops collapsed dirs, inverse-variance
    weighted mean across surviving dirs, σ_ref from nearest-vertex max."""
    Lx, Ly, Tx, Ty = geom
    iref     = cost_mod._tile_interp(ref, Lx, Ly, Tx, Ty, "conn",     copies=2)
    ref_paths  = cost_mod.boundary_paths(Lx, Ly, Tx, Ty)
    test_paths = cost_mod.boundary_paths(Lx, Ly, Tx, Ty)
    pds, sds = [], []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = cost_mod._direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            continue
        N = len(ks)
        rex, rey = rdm + 0.5*rdn, math.sqrt(3)/2*rdn
        Gt_l, et_l, t_l = [], [], []
        for k, mk, nk in zip(ks, ms, ns):
            ent = cost_mod._lookup_test_value(test, mk, nk, Lx, Ly, Tx, Ty, copies=2)
            if ent is None:
                continue
            Gt_l.append(ent["conn"]); et_l.append(ent["conn_err"]); t_l.append(k/N)
        if len(Gt_l) < 2:
            continue
        Gt = np.array(Gt_l); et = np.abs(np.array(et_l)); t_arr = np.array(t_l)
        pts = np.column_stack([t_arr*rex, t_arr*rey])
        Gr = np.asarray(iref(pts), float)
        er = _nearest_vertex_max_err(ref, geom, pts, copies=2)
        mask = np.isfinite(Gr) & np.isfinite(Gt) & np.isfinite(er) & np.isfinite(et)
        if mask.sum() < 2:
            continue
        diff = Gt[mask] - Gr[mask]
        var  = er[mask]**2 + et[mask]**2
        Cd = float(np.mean(diff**2))
        Vd = float(np.sum((2*diff)**2 * var) / mask.sum()**2)
        pds.append(Cd); sds.append(math.sqrt(max(Vd, 0.0)))
    if not pds:
        return float("nan"), float("nan")
    # Inverse-variance weighted mean.
    w = np.array([1.0 / (s*s) if s > 0 else 0.0 for s in sds])
    if w.sum() == 0:
        return float(np.mean(pds)), float(np.sqrt(np.sum(np.array(sds)**2)) / len(pds))
    pds_a = np.array(pds)
    cost = float(np.sum(w * pds_a) / w.sum())
    sig  = float(math.sqrt(1.0 / w.sum()))
    return cost, sig


def test_T5_discrimination(reps: Dict[str, List[dict]]) -> dict:
    print("\n=== T5: discrimination Z (current vs proposed) ===")
    eq = reps["equil"]
    a2 = reps["aniso_r2"]
    a5 = reps["aniso_r5"]

    def _matrix(refs, tests, fn):
        out = []
        for i, r in enumerate(refs):
            for j, t in enumerate(tests):
                if refs is tests and i == j:
                    continue
                out.append(fn(r, t))
        return np.array(out)

    def cur(r, t):
        c, _, _, _ = quiet_cost(r, t, LXY, LXY,
                                cost_mode="test_native", cost_power=2)
        return c

    def prop(r, t):
        c, _ = _proposed_native(r, t, LXY)
        return c

    summary = {}
    for label, fn in [("current", cur), ("proposed_ivw", prop)]:
        null   = _matrix(eq, eq, fn)
        sig_r2 = _matrix(eq, a2, fn)
        sig_r5 = _matrix(eq, a5, fn)
        mu_n, sd_n = float(np.mean(null)), float(np.std(null, ddof=1))
        mu_2, sd_2 = float(np.mean(sig_r2)), float(np.std(sig_r2, ddof=1))
        mu_5, sd_5 = float(np.mean(sig_r5)), float(np.std(sig_r5, ddof=1))
        z2 = (mu_2 - mu_n) / sd_n if sd_n else float("inf")
        z5 = (mu_5 - mu_n) / sd_n if sd_n else float("inf")
        summary[label] = dict(
            null=dict(mean=mu_n, std=sd_n, n=len(null)),
            sig_r2=dict(mean=mu_2, std=sd_2, n=len(sig_r2), z=z2),
            sig_r5=dict(mean=mu_5, std=sd_5, n=len(sig_r5), z=z5),
        )
        print(f"  [{label:13s}] null  μ={mu_n:.3e}±{sd_n:.1e}  "
              f"r2 μ={mu_2:.3e} (Z={z2:6.2f})  "
              f"r5 μ={mu_5:.3e} (Z={z5:6.2f})")
    return summary


# ===========================================================================
# Main
# ===========================================================================

def main():
    log_path = os.path.join(OUT_DIR, "triage.log")
    json_path = os.path.join(OUT_DIR, "triage_summary.json")
    log_buf = io.StringIO()

    class Tee:
        def __init__(self, *streams): self.s = streams
        def write(self, x):
            for s in self.s:
                s.write(x)
        def flush(self):
            for s in self.s:
                s.flush()

    sys.stdout = Tee(sys.__stdout__, log_buf)

    print("=" * 72)
    print("triage_native_test — read-only diagnostic of test_native cost")
    print("=" * 72)

    reps = load_replicates()
    ref0 = reps["equil"][0]
    test0 = reps["equil"][1]      # different replicate, same lattice
    test_aniso = reps["aniso_r2"][0]

    out = {}
    out["T1_determinism"]      = test_T1_determinism(ref0, test0)
    out["T2_hull_dropouts"]    = test_T2_hull_dropouts(ref0, test0)
    out["T3_silent_zero"]      = test_T3_silent_zero(ref0, test0)
    out["T4_err_interp_bias"]  = test_T4_err_interp_bias(ref0, test_aniso)
    out["T5_discrimination"]   = test_T5_discrimination(reps)

    print("\n" + "=" * 72)
    print("summary →", os.path.relpath(json_path, HERE))
    print("=" * 72)

    sys.stdout = sys.__stdout__
    with open(log_path, "w") as f:
        f.write(log_buf.getvalue())
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2, default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
