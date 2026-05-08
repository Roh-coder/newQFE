#!/usr/bin/env python3
"""
precompute_grid.py — MC precompute for cost landscapes over (r1, r2).

For every (r1, r2) on a user-defined grid, run the Ising simulator at the
finite-size pseudo-critical β_c on TWO independently configurable geometries:

  - "test" lattice: shape set by --test-{Lx,Ly}-mult, --test-{Tx,Ty}-frac
  - "ref"  lattice: shape set by --ref-{Lx,Ly}-mult,  --ref-{Tx,Ty}-frac

Geometry parameterisation
--------------------------
Each lattice dimension is expressed as a multiple (Lx, Ly) or fraction
(Tx, Ty) of the FSS scale L:

    Lx = round(Lx_mult * L)      Ly = round(Ly_mult * L)
    Tx = round(Tx_frac * L)      Ty = round(Ty_frac * L)

  Default test:  (mult=1, 1, frac=0,    0)    → untwisted square torus
  Default ref:   (mult=1, 1, frac=0.25, 0.25) → symmetric quarter-twist

Repeat for several lattice sizes L ∈ {L_1, L_2, …}.  Each (r1, r2, L)
point produces one pickle per geometry under

    results/<tag>/grid/L<L>/{test,ref}/r1_X.XXX_r2_Y.YYY.pkl

Resumable: existing pkls are reused.  Designed to be run offline at high
statistics; the matching ``compute_landscape.py`` does a continuum
extrapolation of the cost in 1/L per (r1, r2).

Sanity check: set all ref fracs equal to all test fracs — cost should
vanish inside statistical noise (test = ref geometry).
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine  # noqa: E402


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _make_geom(L, lx_mult, ly_mult, tx_frac, ty_frac):
    """Return (Lx, Ly, Tx, Ty) for scale L and fractional parameters."""
    return (int(round(lx_mult * L)), int(round(ly_mult * L)),
            int(round(tx_frac * L)), int(round(ty_frac * L)))


# Legacy shims — kept so external code that imports them still works.
def _geom_test(L: int) -> tuple[int, int, int, int]:
    return L, L, 0, 0


def _geom_ref(L: int, twist_frac: float) -> tuple[int, int, int, int]:
    T = int(round(twist_frac * L))
    return L, L, T, T


def _point_path(grid_dir: str, r1: float, r2: float) -> str:
    return os.path.join(grid_dir, f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")


# ---------------------------------------------------------------------------
# Per-point worker
# ---------------------------------------------------------------------------

def _run_one(args):
    """Worker: one (geom, L, r1, r2). Returns (label, status, info)."""
    (exe, label, Lx, Ly, Tx, Ty, r1, r2,
     n_traj_prod, n_traj_scan_coarse, n_traj_scan_fine,
     beta_seed_lo, beta_seed_hi, scratch_root, out_pkl) = args

    if os.path.exists(out_pkl):
        return (label, "skip", {})

    scratch = os.path.join(scratch_root, label)
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_c, beta_c_sigma, chi_peak, *_ = mc_engine.find_beta_c(
            exe, Lx, Ly, Tx, Ty, r1, r2, 1.0,
            beta_seed_lo, beta_seed_hi,
            n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
            n_traj_coarse=n_traj_scan_coarse,
            n_traj_fine=n_traj_scan_fine,
            max_shifts=4, jackknife=False,
            data_dir=os.path.join(scratch, "scan"),
        )
        prod_dir = os.path.join(scratch, "prod")
        _, subdir = mc_engine.run_simulator(
            exe, Lx, Ly, Tx, Ty, r1, r2, 1.0, beta_c,
            n_traj=n_traj_prod, data_dir=prod_dir,
        )
        if subdir is None:
            return (label, "err", {"msg": "no MC subdir"})
        a2a = os.path.join(subdir, "two_point_all_to_all.dat")
        test_data = mc_engine.load_all_to_all(a2a)
        wall = time.time() - t0
        payload = {
            "r1": float(r1), "r2": float(r2), "k3": 1.0,
            "Lx": int(Lx), "Ly": int(Ly), "Tx": int(Tx), "Ty": int(Ty),
            "beta_c": float(beta_c),
            "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "n_traj_prod": int(n_traj_prod),
            "wall_s": float(wall),
            "data": test_data,
        }
        tmp = out_pkl + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(payload, f)
        os.replace(tmp, out_pkl)
        shutil.rmtree(scratch, ignore_errors=True)
        return (label, "ok",
                {"beta_c": beta_c, "wall": wall})
    except Exception as e:  # noqa: BLE001
        return (label, "err", {"msg": str(e)})


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[8, 12, 16],
                    help="lattice sizes L for continuum extrapolation")
    # test geometry
    ap.add_argument("--test-Lx-mult", type=float, default=1.0, metavar="M",
                    help="test Lx = round(M*L)  [default 1.0]")
    ap.add_argument("--test-Ly-mult", type=float, default=1.0, metavar="M",
                    help="test Ly = round(M*L)  [default 1.0]")
    ap.add_argument("--test-Tx-frac", type=float, default=0.0, metavar="F",
                    help="test Tx = round(F*L)  [default 0.0 = untwisted]")
    ap.add_argument("--test-Ty-frac", type=float, default=0.0, metavar="F",
                    help="test Ty = round(F*L)  [default 0.0]")
    # ref geometry
    ap.add_argument("--ref-Lx-mult", type=float, default=1.0, metavar="M",
                    help="ref  Lx = round(M*L)  [default 1.0]")
    ap.add_argument("--ref-Ly-mult", type=float, default=1.0, metavar="M",
                    help="ref  Ly = round(M*L)  [default 1.0]")
    ap.add_argument("--ref-Tx-frac", type=float, default=0.25, metavar="F",
                    help="ref  Tx = round(F*L)  [default 0.25]")
    ap.add_argument("--ref-Ty-frac", type=float, default=0.25, metavar="F",
                    help="ref  Ty = round(F*L)  [default 0.25]")
    ap.add_argument("--r-min", type=float, default=0.5)
    ap.add_argument("--r-max", type=float, default=3.0)
    ap.add_argument("--r-step", type=float, default=0.5)
    ap.add_argument("--n-traj", type=int, default=20000,
                    help="production MC trajectories per (geom, L, r) point")
    ap.add_argument("--n-traj-coarse", type=int, default=2000)
    ap.add_argument("--n-traj-fine", type=int, default=4000)
    ap.add_argument("--beta-seed", type=float, nargs=2,
                    default=(0.05, 0.40),
                    help="initial β bracket")
    ap.add_argument("--n-workers", type=int, default=4)
    ap.add_argument("--tag", type=str, default="default",
                    help="results/<tag>/  (default: 'default')")
    ap.add_argument("--exe", type=str, default=None,
                    help="path to ising_tri_twisted_parallelogram "
                         "(default: bin/ in this folder)")
    args = ap.parse_args()

    exe = args.exe or os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        print(f"ERROR: simulator not found: {exe}", file=sys.stderr)
        return 1

    out_root = os.path.join(_HERE, "results", args.tag)
    scratch_root = os.path.join(out_root, "_mc_scratch")
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(scratch_root, exist_ok=True)

    rs = np.arange(args.r_min, args.r_max + 1e-9, args.r_step)
    points = [(float(r1), float(r2)) for r1 in rs for r2 in rs]

    jobs = []
    for L in args.sizes:
        geom_test = _make_geom(L, args.test_Lx_mult, args.test_Ly_mult,
                                args.test_Tx_frac, args.test_Ty_frac)
        geom_ref  = _make_geom(L, args.ref_Lx_mult,  args.ref_Ly_mult,
                                args.ref_Tx_frac,  args.ref_Ty_frac)
        print(f"[plan] L={L:3d}  test={geom_test}  ref={geom_ref}")
        for kind, geom in (("test", geom_test), ("ref", geom_ref)):
            Lx, Ly, Tx, Ty = geom
            grid_dir = os.path.join(out_root, "grid", f"L{L}", kind)
            os.makedirs(grid_dir, exist_ok=True)
            for r1, r2 in points:
                out_pkl = _point_path(grid_dir, r1, r2)
                label = f"L{L}_{kind}_r1{r1:.3f}_r2{r2:.3f}"
                jobs.append((exe, label, Lx, Ly, Tx, Ty, r1, r2,
                             args.n_traj, args.n_traj_coarse, args.n_traj_fine,
                             args.beta_seed[0], args.beta_seed[1],
                             scratch_root, out_pkl))

    n_total = len(jobs)
    n_existing = sum(1 for j in jobs if os.path.exists(j[-1]))
    print(f"[plan] sizes={args.sizes}")
    print(f"[plan] r∈[{args.r_min},{args.r_max}] step={args.r_step} "
          f"→ {len(points)} (r1,r2) pts × {len(args.sizes)} L × 2 geoms = "
          f"{n_total} jobs")
    print(f"[plan] {n_existing} cached → {n_total - n_existing} to compute  "
          f"workers={args.n_workers}")
    print(f"[plan] n_traj_prod={args.n_traj}  scan_coarse={args.n_traj_coarse}  "
          f"scan_fine={args.n_traj_fine}")
    print(f"[plan] out: {out_root}")

    t_wall = time.time()
    n_ok = n_skip = n_err = 0
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs = {ex.submit(_run_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label, status, info = fut.result()
            if status == "ok":
                n_ok += 1
                print(f"[{i:4d}/{n_total}] ok  {label}  "
                      f"β_c={info['beta_c']:.5f}  wall={info['wall']:.1f}s")
            elif status == "skip":
                n_skip += 1
            else:
                n_err += 1
                print(f"[{i:4d}/{n_total}] ERR {label}  "
                      f"{info.get('msg', '?')}")
    wall = time.time() - t_wall

    manifest = {
        "tag": args.tag,
        "sizes": args.sizes,
        "test_Lx_mult": args.test_Lx_mult, "test_Ly_mult": args.test_Ly_mult,
        "test_Tx_frac": args.test_Tx_frac, "test_Ty_frac": args.test_Ty_frac,
        "ref_Lx_mult":  args.ref_Lx_mult,  "ref_Ly_mult":  args.ref_Ly_mult,
        "ref_Tx_frac":  args.ref_Tx_frac,  "ref_Ty_frac":  args.ref_Ty_frac,
        "r_min": args.r_min, "r_max": args.r_max, "r_step": args.r_step,
        "n_traj_prod": args.n_traj,
        "n_traj_scan_coarse": args.n_traj_coarse,
        "n_traj_scan_fine": args.n_traj_fine,
        "beta_seed": list(args.beta_seed),
        "wall_s": round(wall, 1),
        "n_total": n_total, "n_ok_this_run": n_ok,
        "n_skip_this_run": n_skip, "n_err_this_run": n_err,
        "test_geoms": {str(L): list(_make_geom(L, args.test_Lx_mult,
                                               args.test_Ly_mult,
                                               args.test_Tx_frac,
                                               args.test_Ty_frac))
                       for L in args.sizes},
        "ref_geoms":  {str(L): list(_make_geom(L, args.ref_Lx_mult,
                                               args.ref_Ly_mult,
                                               args.ref_Tx_frac,
                                               args.ref_Ty_frac))
                       for L in args.sizes},
    }
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\n[done] wall={wall:.1f}s  ok={n_ok} skip={n_skip} err={n_err}")
    print(f"[done] manifest → {os.path.join(out_root, 'manifest.json')}")
    shutil.rmtree(scratch_root, ignore_errors=True)
    return 0 if n_err == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
