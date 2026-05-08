"""
precompute_landscape.py — dense (r1, r2) MC grid for offline cost replay.

For each grid point on test_geom (Lx, Ly, Tx, Ty):
  1. Find β_c via mc_engine.find_beta_c (light scan budget).
  2. Run production MC at β_c.
  3. Pickle {r1, r2, β_c, β_c_sigma, n_traj_prod, test_data} into
     results/_landscape/<test_geom>/grid/r1_X.XXX_r2_Y.YYY.pkl
Skips points whose pkl already exists (resumable).

Designed to be the one-time expensive precursor to landscape_cost.py.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine  # noqa: E402


def _point_path(grid_dir: str, r1: float, r2: float) -> str:
    return os.path.join(grid_dir, f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")


def _run_one(args):
    """Worker: compute one grid point. Returns (r1, r2, status, info)."""
    (exe, Lx, Ly, Tx, Ty, r1, r2,
     n_traj_prod, n_traj_scan_coarse, n_traj_scan_fine,
     beta_seed_lo, beta_seed_hi, scratch_root, out_pkl) = args

    if os.path.exists(out_pkl):
        return (r1, r2, "skip", {})

    label = f"r1_{r1:.3f}_r2_{r2:.3f}"
    scratch = os.path.join(scratch_root, label)
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_c, beta_c_sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
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
            return (r1, r2, "err", {"msg": "no MC subdir"})
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
            "test_data": test_data,
        }
        tmp = out_pkl + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(payload, f)
        os.replace(tmp, out_pkl)
        # Clean MC scratch to save disk.
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
        return (r1, r2, "ok",
                {"beta_c": beta_c, "beta_c_sigma": beta_c_sigma,
                 "wall": wall})
    except Exception as e:  # noqa: BLE001
        return (r1, r2, "err", {"msg": str(e)})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Lx", type=int, default=16)
    ap.add_argument("--Ly", type=int, default=16)
    ap.add_argument("--Tx", type=int, default=0)
    ap.add_argument("--Ty", type=int, default=0)
    ap.add_argument("--r-min", type=float, default=0.5)
    ap.add_argument("--r-max", type=float, default=10.0)
    ap.add_argument("--r-step", type=float, default=1.0)
    ap.add_argument("--n-traj", type=int, default=5000,
                    help="production MC trajectories per point")
    ap.add_argument("--n-traj-coarse", type=int, default=2000)
    ap.add_argument("--n-traj-fine", type=int, default=4000)
    ap.add_argument("--beta-seed", type=float, nargs=2,
                    default=(0.05, 0.20),
                    help="initial β bracket (broad — covers small + large r)")
    ap.add_argument("--n-workers", type=int, default=6)
    ap.add_argument("--out-tag", type=str, default=None,
                    help="results/_landscape/<tag>/  (default: derived)")
    args = ap.parse_args()

    Lx, Ly, Tx, Ty = args.Lx, args.Ly, args.Tx, args.Ty
    tag = args.out_tag or f"Lx{Lx}_Ly{Ly}_Tx{Tx}_Ty{Ty}"
    out_root = os.path.join(_HERE, "results", "_landscape", tag)
    grid_dir = os.path.join(out_root, "grid")
    scratch_root = os.path.join(out_root, "_mc_scratch")
    os.makedirs(grid_dir, exist_ok=True)
    os.makedirs(scratch_root, exist_ok=True)

    exe = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        print(f"ERROR: simulator not found: {exe}", file=sys.stderr)
        return 1

    rs = np.arange(args.r_min, args.r_max + 1e-9, args.r_step)
    points = [(float(r1), float(r2)) for r1 in rs for r2 in rs]

    # Build job list (skip already-done).
    jobs = []
    for r1, r2 in points:
        out_pkl = _point_path(grid_dir, r1, r2)
        jobs.append((exe, Lx, Ly, Tx, Ty, r1, r2,
                     args.n_traj, args.n_traj_coarse, args.n_traj_fine,
                     args.beta_seed[0], args.beta_seed[1],
                     scratch_root, out_pkl))
    n_total = len(jobs)
    n_existing = sum(1 for j in jobs if os.path.exists(j[-1]))
    print(f"[plan] {n_total} grid points  ({n_existing} already cached)  "
          f"→ {n_total - n_existing} to compute  workers={args.n_workers}")
    print(f"[plan] geom=({Lx},{Ly},{Tx},{Ty})  r∈[{args.r_min},{args.r_max}] step={args.r_step}")
    print(f"[plan] n_traj_prod={args.n_traj}  scan_coarse={args.n_traj_coarse}  "
          f"scan_fine={args.n_traj_fine}")
    print(f"[plan] out: {out_root}")

    t_wall = time.time()
    n_ok = n_skip = n_err = 0
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs = {ex.submit(_run_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r1, r2, status, info = fut.result()
            if status == "ok":
                n_ok += 1
                print(f"[{i:4d}/{n_total}] ok  r=({r1:5.2f},{r2:5.2f})  "
                      f"β_c={info['beta_c']:.5f}  wall={info['wall']:.1f}s")
            elif status == "skip":
                n_skip += 1
            else:
                n_err += 1
                print(f"[{i:4d}/{n_total}] ERR r=({r1:5.2f},{r2:5.2f})  "
                      f"{info.get('msg', '?')}")
    wall = time.time() - t_wall

    # Manifest
    manifest = {
        "tag": tag, "geom": [Lx, Ly, Tx, Ty],
        "r_min": args.r_min, "r_max": args.r_max, "r_step": args.r_step,
        "n_traj_prod": args.n_traj,
        "n_traj_scan_coarse": args.n_traj_coarse,
        "n_traj_scan_fine": args.n_traj_fine,
        "beta_seed": list(args.beta_seed),
        "wall_s": round(wall, 1),
        "n_total": n_total, "n_ok_this_run": n_ok,
        "n_skip_this_run": n_skip, "n_err_this_run": n_err,
        "points": [],
    }
    for r1, r2 in points:
        p = _point_path(grid_dir, r1, r2)
        if os.path.exists(p):
            manifest["points"].append({"r1": r1, "r2": r2,
                                        "path": os.path.relpath(p, out_root)})
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\n[done] wall={wall:.1f}s  ok={n_ok} skip={n_skip} err={n_err}")
    print(f"[done] manifest → {os.path.join(out_root, 'manifest.json')}")
    # Clean MC scratch
    import shutil
    shutil.rmtree(scratch_root, ignore_errors=True)
    return 0 if n_err == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
