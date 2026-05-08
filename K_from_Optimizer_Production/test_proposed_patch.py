#!/usr/bin/env python3
"""
test_proposed_patch.py — apply the recommended T3+T4 patch as a SHADOW
function (no edits to cost.py) and test it head-to-head against the
current ``test_native`` cost.

Patch under test
----------------
P3: drop collapsed (gcd==1 or N_used<2) directions from the per-direction
    average instead of treating them as C_d=0, σ_d=0.
P4: per-point reference σ uses the maximum conn_err over the surrounding
    Delaunay simplex vertices (tiled point cloud) instead of a linear
    interpolation of conn_err.

Test matrix
-----------
A. T3 stress: vary geometry to force 0, 1, or 2 directions to collapse,
   confirm current bias is exactly (N_valid/3) and proposed cost is the
   true single-direction cost.

B. T2/Hull stress: scan twist (Tx,Ty) over (12,12,Tx,Ty) for
   |Tx|,|Ty|≤6, count NaN drop-outs at copies=2 (current default) for
   the legacy l4mean_both_interp path.  Identify the smallest twist for
   which copies=2 fails.

C. Discrimination Z (test_native vs proposed) on the cached
   _betac_stability_test data, isotropic L=12 case (5 RW replicates per
   group: equil, aniso_r2, aniso_r5).  Repeats the T5 measurement of
   the previous run with σ_ref from vertex-max.

D. SNR-status flip count: how many evaluator decisions
   (ok/marginal/need_more_stats) change between current σ_cost and
   proposed σ_cost across the 5×5 cross-replicate grid in each group.

Outputs
-------
  results/_native_triage/patch_test_summary.json
  results/_native_triage/patch_test.log
"""
from __future__ import annotations

import contextlib
import glob
import io
import json
import math
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cost as cost_mod          # noqa: E402
import mc_engine                 # noqa: E402

OUT_DIR = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT_DIR, exist_ok=True)

STABILITY_ROOT = os.path.join(HERE, "results", "_betac_stability_test")
LXY = (12, 12, 0, 0)
SQRT3_2 = math.sqrt(3.0) / 2.0


# ---------------------------------------------------------------------------
# Loaders / helpers
# ---------------------------------------------------------------------------

def _load_replicates(case_dir: str) -> List[dict]:
    pat = os.path.join(STABILITY_ROOT, case_dir, "rw_t*", "*",
                       "two_point_all_to_all.dat")
    files = sorted(glob.glob(pat))
    return [mc_engine.load_all_to_all(p) for p in files]


def _quiet(fn, *a, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return fn(*a, **kw)


def _vertex_max_err(ref_data, geom, query_xy, copies=2,
                    field="conn_err"):
    """Per-query maximum of `field` over the surrounding Delaunay simplex
    vertices in the tiled point cloud.  Returns NaN where the query is
    outside the tiled hull."""
    rLx, rLy, rTx, rTy = geom
    keys = list(ref_data.keys())
    m = np.array([k[0] for k in keys], float)
    n = np.array([k[1] for k in keys], float)
    e = np.array([ref_data[k][field] for k in keys], float)
    x0, y0 = m + 0.5 * n, SQRT3_2 * n
    vx, vy = rLx + 0.5 * rTy, SQRT3_2 * rTy
    ux, uy = rTx - 0.5 * rLy, -SQRT3_2 * rLy
    xs, ys, es = [x0], [y0], [e]
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            if a == 0 and b == 0:
                continue
            xs.append(x0 + a * vx + b * ux)
            ys.append(y0 + a * vy + b * uy)
            es.append(e)
    pts = np.column_stack([np.concatenate(xs), np.concatenate(ys)])
    es = np.concatenate(es)
    _, uid = np.unique(np.round(pts, 8), axis=0, return_index=True)
    pts, es = pts[uid], es[uid]
    tri = Delaunay(pts)
    out = np.full(len(query_xy), np.nan)
    simp = tri.find_simplex(query_xy)
    for i, s in enumerate(simp):
        if s >= 0:
            out[i] = es[tri.simplices[s]].max()
    return out


# ---------------------------------------------------------------------------
# Proposed test_native cost (shadow, with P3 + P4 patch)
# ---------------------------------------------------------------------------

def proposed_test_native(ref_data, test_data,
                         test_geom, ref_geom,
                         copies=2, power=2,
                         use_vertex_max_err=True,
                         drop_collapsed=True):
    Lx, Ly, Tx, Ty = test_geom
    rLx, rLy, rTx, rTy = ref_geom
    iref = cost_mod._tile_interp(ref_data, rLx, rLy, rTx, rTy, "conn", copies)
    iref_err_interp = cost_mod._tile_interp(ref_data, rLx, rLy, rTx, rTy,
                                            "conn_err", copies)
    ref_paths = cost_mod.boundary_paths(rLx, rLy, rTx, rTy)
    test_paths = cost_mod.boundary_paths(Lx, Ly, Tx, Ty)
    per_dir, per_sig, valid_flag = [], [], []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = cost_mod._direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            per_dir.append(0.0); per_sig.append(0.0); valid_flag.append(False)
            continue
        N = len(ks)
        Gt_l, et_l, t_l = [], [], []
        for k, mk, nk in zip(ks, ms, ns):
            ent = cost_mod._lookup_test_value(test_data, mk, nk,
                                              Lx, Ly, Tx, Ty, copies=copies)
            if ent is None:
                continue
            Gt_l.append(ent["conn"]); et_l.append(ent["conn_err"]); t_l.append(k / N)
        if len(Gt_l) < 2:
            per_dir.append(0.0); per_sig.append(0.0); valid_flag.append(False)
            continue
        Gt = np.array(Gt_l); et = np.abs(np.array(et_l)); t_arr = np.array(t_l)
        rex, rey = rdm + 0.5 * rdn, SQRT3_2 * rdn
        pts = np.column_stack([t_arr * rex, t_arr * rey])
        Gr = np.asarray(iref(pts), float)
        if use_vertex_max_err:
            er = _vertex_max_err(ref_data, ref_geom, pts, copies=copies)
        else:
            er = np.abs(np.asarray(iref_err_interp(pts), float))
        mask = np.isfinite(Gr) & np.isfinite(Gt) & np.isfinite(er) & np.isfinite(et)
        if mask.sum() < 2:
            per_dir.append(0.0); per_sig.append(0.0); valid_flag.append(False)
            continue
        diff = Gt[mask] - Gr[mask]
        var = er[mask] ** 2 + et[mask] ** 2
        if power == 2:
            Cd = float(np.mean(diff ** 2))
            Vd = float(np.sum((2 * diff) ** 2 * var) / mask.sum() ** 2)
        else:
            Cd = float(np.mean(diff ** 4))
            Vd = float(np.sum((4 * diff ** 3) ** 2 * var) / mask.sum() ** 2)
        per_dir.append(Cd); per_sig.append(math.sqrt(max(Vd, 0.0)))
        valid_flag.append(True)
    if drop_collapsed:
        keep = [(c, s) for c, s, ok in zip(per_dir, per_sig, valid_flag) if ok]
    else:
        keep = list(zip(per_dir, per_sig))
    if not keep:
        return float("nan"), float("nan"), per_dir, per_sig, valid_flag
    cost = sum(c for c, _ in keep) / len(keep)
    sig = math.sqrt(sum(s * s for _, s in keep)) / len(keep)
    return cost, sig, per_dir, per_sig, valid_flag


def current_test_native(ref_data, test_data, test_geom, ref_geom):
    c, s, pd, ps = _quiet(cost_mod.l2_cost,
                          ref_data, test_data,
                          *test_geom, *ref_geom,
                          cost_mode="test_native", cost_power=2)
    return c, s, pd, ps


# ===========================================================================
# A. T3 stress — collapse 0, 1, or 2 directions
# ===========================================================================

def test_A_collapse_bias(ref, test) -> dict:
    print("\n=== A: collapse-bias stress ===")
    cases = [
        ("0 collapse  (12,12,0,0)",  (12, 12, 0, 0)),
        ("1 collapse  (12,11,0,0)",  (12, 11, 0, 0)),
        ("2 collapse  (12,12,0,1)",  (12, 12, 0, 1)),
    ]
    out = []
    for label, geom in cases:
        paths = cost_mod.boundary_paths(*geom)
        gcds = [math.gcd(abs(dm), abs(dn)) for dm, dn in paths]
        n_collapse = sum(1 for g in gcds if g <= 1)
        c_cur, s_cur, _, _ = current_test_native(ref, test, geom, LXY)
        c_pr, s_pr, pd, ps, vf = proposed_test_native(
            ref, test, geom, LXY, drop_collapsed=True, use_vertex_max_err=True)
        n_valid = sum(vf)
        ratio = c_cur / c_pr if c_pr and not math.isnan(c_pr) else float("nan")
        expected = n_valid / 3 if n_valid else float("nan")
        ok = abs(ratio - expected) < 1e-3 if not math.isnan(ratio) else False
        print(f"  {label}  gcds={gcds}  N_valid={n_valid}/3")
        print(f"    current  cost={c_cur:.4e} ± {s_cur:.2e}")
        print(f"    proposed cost={c_pr:.4e} ± {s_pr:.2e}")
        print(f"    ratio cur/prop = {ratio:.4f}  expected N_valid/3 = {expected:.4f}  "
              f"({'OK' if ok else 'MISMATCH'})")
        out.append(dict(label=label, geom=geom, gcds=gcds,
                        n_valid=n_valid,
                        cost_current=c_cur, sigma_current=s_cur,
                        cost_proposed=c_pr, sigma_proposed=s_pr,
                        ratio=ratio, expected=expected, ok=ok))
    return out


# ===========================================================================
# B. Hull-dropout twist scan (legacy path with copies=2)
# ===========================================================================

def test_B_hull_twist_scan(ref, test) -> dict:
    print("\n=== B: hull dropouts vs twist (copies=2 vs 4) ===")
    n_samp = 400
    out = []
    for Tx in (-6, -4, -2, 0, 2, 4, 6):
        for Ty in (-6, -4, -2, 0, 2, 4, 6):
            geom = (12, 12, Tx, Ty)
            try:
                paths = cost_mod.boundary_paths(*geom)
                # Skip degenerate.
                if any(dm == 0 and dn == 0 for dm, dn in paths):
                    continue
                row = dict(geom=geom)
                for copies in (2, 4):
                    iref = cost_mod._tile_interp(ref, *geom, "conn", copies)
                    itest = cost_mod._tile_interp(test, *geom, "conn", copies)
                    drops_r, drops_t = 0, 0
                    t = np.linspace(0.0, 1.0, n_samp, endpoint=False)
                    for (rdm, rdn), (tdm, tdn) in zip(paths, paths):
                        rex, rey = rdm + 0.5 * rdn, SQRT3_2 * rdn
                        Gr = iref(np.column_stack([t * rex, t * rey]))
                        Gt = itest(np.column_stack([t * rex, t * rey]))
                        drops_r += int(np.isnan(Gr).sum())
                        drops_t += int(np.isnan(Gt).sum())
                    row[f"drops_copies_{copies}"] = (drops_r, drops_t)
                if (row["drops_copies_2"] != (0, 0)
                        or row["drops_copies_4"] != (0, 0)):
                    print(f"  {geom}  c2={row['drops_copies_2']}  "
                          f"c4={row['drops_copies_4']}")
                out.append(row)
            except Exception as e:
                out.append(dict(geom=geom, error=str(e)))
    n_bad_2 = sum(1 for r in out
                  if "drops_copies_2" in r and r["drops_copies_2"] != (0, 0))
    n_bad_4 = sum(1 for r in out
                  if "drops_copies_4" in r and r["drops_copies_4"] != (0, 0))
    print(f"  scan summary: {len(out)} geometries  "
          f"with drops at copies=2: {n_bad_2}  copies=4: {n_bad_4}")
    return dict(rows=out, n_bad_copies_2=n_bad_2, n_bad_copies_4=n_bad_4)


# ===========================================================================
# C / D. Discrimination + SNR-status flips on the cached stability data
# ===========================================================================

CASES = {"equil": "equil_12", "aniso_r2": "aniso_r2_12", "aniso_r5": "aniso_r5_12"}


def test_C_D_discrimination_and_snr(reps) -> dict:
    print("\n=== C/D: discrimination + SNR-status flips ===")
    eq, a2, a5 = reps["equil"], reps["aniso_r2"], reps["aniso_r5"]

    def _mat(refs, tests, fn):
        out = []
        for i, r in enumerate(refs):
            for j, t in enumerate(tests):
                if refs is tests and i == j:
                    continue
                out.append(fn(r, t))
        return out

    def cur(r, t):
        c, s, _, _ = current_test_native(r, t, LXY, LXY)
        return c, s

    def prop(r, t):
        c, s, _, _, _ = proposed_test_native(r, t, LXY, LXY,
                                             drop_collapsed=True,
                                             use_vertex_max_err=True)
        return c, s

    def _summarise(label, fn):
        null = _mat(eq, eq, fn)
        s_r2 = _mat(eq, a2, fn)
        s_r5 = _mat(eq, a5, fn)
        c_null = np.array([x[0] for x in null])
        c_r2 = np.array([x[0] for x in s_r2])
        c_r5 = np.array([x[0] for x in s_r5])
        sig_null = np.array([x[1] for x in null])
        sig_r2 = np.array([x[1] for x in s_r2])
        sig_r5 = np.array([x[1] for x in s_r5])
        sd_n = float(np.std(c_null, ddof=1))
        z2 = (np.mean(c_r2) - np.mean(c_null)) / sd_n if sd_n else float("inf")
        z5 = (np.mean(c_r5) - np.mean(c_null)) / sd_n if sd_n else float("inf")
        # SNR per pair
        def _status_counts(c_arr, s_arr):
            stats = {"ok": 0, "marginal": 0, "need_more_stats": 0}
            for c, s in zip(c_arr, s_arr):
                stats[cost_mod.snr_status(c, s)] += 1
            return stats
        st = {
            "null": _status_counts(c_null, sig_null),
            "r2":   _status_counts(c_r2,  sig_r2),
            "r5":   _status_counts(c_r5,  sig_r5),
        }
        print(f"  [{label:9s}] Z(r2)={z2:6.2f}   Z(r5)={z5:6.2f}   "
              f"null σ={sd_n:.2e}")
        print(f"            status counts:  null={st['null']}  "
              f"r2={st['r2']}  r5={st['r5']}")
        return dict(z_r2=z2, z_r5=z5, std_null=sd_n,
                    status_counts=st,
                    means=dict(null=float(np.mean(c_null)),
                               r2=float(np.mean(c_r2)),
                               r5=float(np.mean(c_r5))),
                    sigma_means=dict(null=float(np.mean(sig_null)),
                                     r2=float(np.mean(sig_r2)),
                                     r5=float(np.mean(sig_r5))))

    summary = {"current": _summarise("current", cur),
               "proposed": _summarise("proposed", prop)}
    # Per-pair status flips
    flips = {"null": 0, "r2": 0, "r5": 0}
    for grp, refs, tests in [("null", eq, eq), ("r2", eq, a2), ("r5", eq, a5)]:
        for i, r in enumerate(refs):
            for j, t in enumerate(tests):
                if refs is tests and i == j:
                    continue
                cc, sc = cur(r, t)
                pc, ps = prop(r, t)
                if cost_mod.snr_status(cc, sc) != cost_mod.snr_status(pc, ps):
                    flips[grp] += 1
    summary["snr_status_flips"] = flips
    print(f"  per-pair SNR-status flips: {flips}")
    return summary


# ===========================================================================
# Main
# ===========================================================================

def main():
    log_buf = io.StringIO()

    class Tee:
        def __init__(self, *s): self.s = s
        def write(self, x):
            for s in self.s:
                s.write(x)
        def flush(self):
            for s in self.s:
                s.flush()

    sys.stdout = Tee(sys.__stdout__, log_buf)

    print("=" * 72)
    print("test_proposed_patch — shadow-test of T3+T4 patch (no code edits)")
    print("=" * 72)

    reps = {tag: _load_replicates(sub) for tag, sub in CASES.items()}
    for tag, lst in reps.items():
        print(f"  loaded {tag}: {len(lst)} replicates")

    out = {}
    out["A_collapse_bias"]    = test_A_collapse_bias(reps["equil"][0],
                                                     reps["equil"][1])
    out["B_hull_twist_scan"]  = test_B_hull_twist_scan(reps["equil"][0],
                                                       reps["equil"][1])
    out["CD_discrimination"]  = test_C_D_discrimination_and_snr(reps)

    sys.stdout = sys.__stdout__
    log_path = os.path.join(OUT_DIR, "patch_test.log")
    json_path = os.path.join(OUT_DIR, "patch_test_summary.json")
    with open(log_path, "w") as f:
        f.write(log_buf.getvalue())
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    print(log_buf.getvalue())
    print(f"\nwrote {os.path.relpath(log_path, HERE)}")
    print(f"wrote {os.path.relpath(json_path, HERE)}")


if __name__ == "__main__":
    main()
