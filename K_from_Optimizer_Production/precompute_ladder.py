"""
precompute_ladder.py — multi-L MC for the FSS ladder (plan §9.10).

Two modes:

  --mode ref     Run untwisted-square reference MC at sizes (16,24,32)
                 at the isotropic point (r1=r2=k3=1), one run per size.
                 Output: results/_ladder_111/ref/L{L}/two_point_all_to_all.dat

  --mode test    Sweep test (r1, r2) grid at sizes (8,12,16). For each
                 (r1, r2, L_test) find β_c then run production MC.
                 Output: results/_ladder_111/test/grid/r1_X_r2_Y_L{L}.pkl

Resumable: skips outputs that already exist.
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


# ---------------------------------------------------------------------
# Single-job workers
# ---------------------------------------------------------------------

def _run_ref_one(args):
    (exe, L, n_traj, n_traj_coarse, n_traj_fine,
     beta_seed_lo, beta_seed_hi, scratch_root, out_dir) = args
    out_dat = os.path.join(out_dir, "two_point_all_to_all.dat")
    if os.path.exists(out_dat):
        return ("skip", L, {})
    os.makedirs(out_dir, exist_ok=True)
    scratch = os.path.join(scratch_root, f"ref_L{L}")
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_c, beta_c_sigma, chi_peak, *_ = mc_engine.find_beta_c(
            exe, L, L, 0, 0, 1.0, 1.0, 1.0,
            beta_seed_lo, beta_seed_hi,
            n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
            n_traj_coarse=n_traj_coarse, n_traj_fine=n_traj_fine,
            max_shifts=4, jackknife=False,
            data_dir=os.path.join(scratch, "scan"),
        )
        prod_dir = os.path.join(scratch, "prod")
        _, subdir = mc_engine.run_simulator(
            exe, L, L, 0, 0, 1.0, 1.0, 1.0, beta_c,
            n_traj=n_traj, data_dir=prod_dir,
        )
        if subdir is None:
            return ("err", L, {"msg": "no MC subdir"})
        a2a = os.path.join(subdir, "two_point_all_to_all.dat")
        shutil.copy(a2a, out_dat)
        meta = {
            "L": int(L), "Tx": 0, "Ty": 0,
            "k1": 1.0, "k2": 1.0, "k3": 1.0,
            "beta_c": float(beta_c),
            "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "n_traj_prod": int(n_traj),
            "wall_s": time.time() - t0,
        }
        with open(os.path.join(out_dir, "meta.json"), "w") as f:
            json.dump(meta, f, indent=2)
        shutil.rmtree(scratch, ignore_errors=True)
        return ("ok", L, meta)
    except Exception as e:  # noqa: BLE001
        return ("err", L, {"msg": str(e)})


def _run_test_one(args):
    (exe, L, r1, r2, n_traj, n_traj_coarse, n_traj_fine,
     beta_seed_lo, beta_seed_hi, scratch_root, out_pkl) = args
    if os.path.exists(out_pkl):
        return ("skip", L, r1, r2, {})
    label = f"L{L}_r1_{r1:.3f}_r2_{r2:.3f}"
    scratch = os.path.join(scratch_root, label)
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_c, beta_c_sigma, chi_peak, *_ = mc_engine.find_beta_c(
            exe, L, L, 0, 0, r1, r2, 1.0,
            beta_seed_lo, beta_seed_hi,
            n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
            n_traj_coarse=n_traj_coarse, n_traj_fine=n_traj_fine,
            max_shifts=4, jackknife=False,
            data_dir=os.path.join(scratch, "scan"),
        )
        prod_dir = os.path.join(scratch, "prod")
        _, subdir = mc_engine.run_simulator(
            exe, L, L, 0, 0, r1, r2, 1.0, beta_c,
            n_traj=n_traj, data_dir=prod_dir,
        )
        if subdir is None:
            return ("err", L, r1, r2, {"msg": "no MC subdir"})
        a2a = os.path.join(subdir, "two_point_all_to_all.dat")
        test_data = mc_engine.load_all_to_all(a2a)
        payload = {
            "r1": float(r1), "r2": float(r2), "k3": 1.0,
            "Lx": int(L), "Ly": int(L), "Tx": 0, "Ty": 0,
            "beta_c": float(beta_c),
            "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "n_traj_prod": int(n_traj),
            "wall_s": time.time() - t0,
            "test_data": test_data,
        }
        tmp = out_pkl + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(payload, f)
        os.replace(tmp, out_pkl)
        shutil.rmtree(scratch, ignore_errors=True)
        return ("ok", L, r1, r2,
                {"beta_c": beta_c, "wall": time.time() - t0})
    except Exception as e:  # noqa: BLE001
        return ("err", L, r1, r2, {"msg": str(e)})


# ---------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------

def _do_ref(args):
    out_root = os.path.join(_HERE, "results", args.tag, "ref")
    scratch = os.path.join(_HERE, "results", args.tag, "_scratch_ref")
    os.makedirs(scratch, exist_ok=True)
    exe = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    jobs = []
    for L in args.ref_sizes:
        out_dir = os.path.join(out_root, f"L{L}")
        jobs.append((exe, int(L), args.n_traj_ref,
                     args.n_traj_coarse, args.n_traj_fine,
                     args.beta_seed[0], args.beta_seed[1],
                     scratch, out_dir))
    print(f"[ref] sizes={args.ref_sizes}  workers={args.n_workers}")
    n_ok = n_skip = n_err = 0
    t_wall = time.time()
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs = [ex.submit(_run_ref_one, j) for j in jobs]
        for fut in as_completed(futs):
            status, L, info = fut.result()
            if status == "ok":
                n_ok += 1
                print(f"[ref] ok  L={L:3d}  β_c={info['beta_c']:.5f}  "
                      f"wall={info['wall_s']:.1f}s")
            elif status == "skip":
                n_skip += 1
                print(f"[ref] skip L={L}")
            else:
                n_err += 1
                print(f"[ref] ERR L={L}: {info.get('msg','?')}")
    shutil.rmtree(scratch, ignore_errors=True)
    manifest = {
        "ref_sizes": args.ref_sizes,
        "n_traj_prod": args.n_traj_ref,
        "wall_s": round(time.time() - t_wall, 1),
        "n_ok": n_ok, "n_skip": n_skip, "n_err": n_err,
    }
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    return 0 if n_err == 0 else 2


def _do_test(args):
    out_root = os.path.join(_HERE, "results", args.tag, "test")
    grid_dir = os.path.join(out_root, "grid")
    scratch = os.path.join(_HERE, "results", args.tag, "_scratch_test")
    os.makedirs(grid_dir, exist_ok=True)
    os.makedirs(scratch, exist_ok=True)
    exe = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if args.points:
        if len(args.points) % 2:
            print("ERROR: --points needs an even number of values "
                  "(r1 r2 r1 r2 ...)", file=sys.stderr)
            return 1
        points = [(float(args.points[2*i]), float(args.points[2*i+1]))
                  for i in range(len(args.points) // 2)]
    else:
        rs = np.arange(args.r_min, args.r_max + 1e-9, args.r_step)
        points = [(float(r1), float(r2)) for r1 in rs for r2 in rs]
    jobs = []
    for L in args.test_sizes:
        for r1, r2 in points:
            out_pkl = os.path.join(
                grid_dir, f"r1_{r1:.3f}_r2_{r2:.3f}_L{L}.pkl")
            jobs.append((exe, int(L), r1, r2, args.n_traj,
                         args.n_traj_coarse, args.n_traj_fine,
                         args.beta_seed[0], args.beta_seed[1],
                         scratch, out_pkl))
    n_total = len(jobs)
    n_existing = sum(1 for j in jobs if os.path.exists(j[-1]))
    print(f"[test] sizes={args.test_sizes}  grid={len(points)}  "
          f"jobs={n_total} ({n_existing} cached)  workers={args.n_workers}")
    n_ok = n_skip = n_err = 0
    t_wall = time.time()
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs = {ex.submit(_run_test_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            status, L, r1, r2, info = fut.result()
            if status == "ok":
                n_ok += 1
                print(f"[{i:4d}/{n_total}] ok  L={L:2d} r=({r1:5.2f},{r2:5.2f})  "
                      f"β_c={info['beta_c']:.5f}  wall={info['wall']:.1f}s")
            elif status == "skip":
                n_skip += 1
            else:
                n_err += 1
                print(f"[{i:4d}/{n_total}] ERR L={L} r=({r1:5.2f},{r2:5.2f}): "
                      f"{info.get('msg','?')}")
    shutil.rmtree(scratch, ignore_errors=True)
    manifest = {
        "test_sizes": args.test_sizes,
        "r_min": args.r_min, "r_max": args.r_max, "r_step": args.r_step,
        "n_traj_prod": args.n_traj,
        "wall_s": round(time.time() - t_wall, 1),
        "n_total": n_total, "n_ok_this_run": n_ok,
        "n_skip_this_run": n_skip, "n_err_this_run": n_err,
    }
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    return 0 if n_err == 0 else 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["ref", "test"], required=True)
    ap.add_argument("--tag", default="_ladder_111",
                    help="results subdirectory")
    ap.add_argument("--ref-sizes", type=int, nargs="+",
                    default=[16, 24, 32])
    ap.add_argument("--test-sizes", type=int, nargs="+",
                    default=[8, 12, 16])
    ap.add_argument("--r-min", type=float, default=0.5)
    ap.add_argument("--r-max", type=float, default=2.5)
    ap.add_argument("--r-step", type=float, default=0.25)
    ap.add_argument("--points", type=float, nargs="+", default=None,
                    help="explicit list of test points: r1 r2 r1 r2 ...; "
                         "overrides r-min/r-max/r-step grid")
    ap.add_argument("--n-traj", type=int, default=2000,
                    help="production traj per test point")
    ap.add_argument("--n-traj-ref", type=int, default=20000,
                    help="production traj per reference run")
    ap.add_argument("--n-traj-coarse", type=int, default=2000)
    ap.add_argument("--n-traj-fine", type=int, default=4000)
    ap.add_argument("--beta-seed", type=float, nargs=2,
                    default=(0.15, 0.40),
                    help="initial β bracket; for 1-1-1 truth β_c≈0.275")
    ap.add_argument("--n-workers", type=int, default=6)
    args = ap.parse_args()
    return _do_ref(args) if args.mode == "ref" else _do_test(args)


if __name__ == "__main__":
    raise SystemExit(main())
