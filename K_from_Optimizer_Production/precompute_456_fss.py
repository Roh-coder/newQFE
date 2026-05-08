#!/usr/bin/env python3
"""
precompute_456_fss.py  —  Run MC for the 4-5-6 FSS correlator study.

Produces two families:

  REF — twisted (13α,16α,-3α,3α), k1=k2=k3=1, α ∈ {1,2,3}
        α=1 and α=3 already exist; only α=2 is new.
        Output: results/_fss_456/ref/a{α}/two_point_all_to_all.dat

  TEST — untwisted (L,L,0,0) at exact truth couplings
         r1 = (2ln5-ln7)/(2ln3-ln7) ≈ 5.06523
         r2 = ln7/(2ln3-ln7)        ≈ 7.74293
         L ∈ {8, 16, 24, 32}  (all multiples of 8 for clean sample sites)
         Output: results/_fss_456/test/L{L}/two_point_all_to_all.dat

Runs missing pieces only (resumable).  Parallelises across 4 workers.
"""
from __future__ import annotations

import json
import math
import os
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine

# ──────────────────────────────────────────────────────────────────────────────
# Truth couplings  (exact analytic)
# ──────────────────────────────────────────────────────────────────────────────
_R1   = (2*math.log(5) - math.log(7)) / (2*math.log(3) - math.log(7))   # ≈ 5.06523
_R2   = math.log(7) / (2*math.log(3) - math.log(7))                      # ≈ 7.74293
_K3   = 1.0
_BINF = math.log(3.0 / math.sqrt(7.0)) / 2.0                             # ≈ 0.06283

_OUT = os.path.join(_HERE, "results", "_fss_456")
_EXE = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")


# ──────────────────────────────────────────────────────────────────────────────
# Job definitions
# ──────────────────────────────────────────────────────────────────────────────

def _jobs():
    jobs = []

    # ---- REFERENCE: twisted (13α,16α,-3α,3α) --------------------------------
    for alpha in [1, 2, 3]:
        Lx, Ly = 13*alpha, 16*alpha
        Tx, Ty = -3*alpha, 3*alpha
        out_dir = os.path.join(_OUT, "ref", f"a{alpha}")
        out_dat = os.path.join(out_dir, "two_point_all_to_all.dat")
        # α=1 already at results/_reference_Lx13_Ly16_Tx-3_Ty3
        # α=3 already at results/_reference_Lx39_Ly48_Tx-9_Ty9
        src = {
            1: os.path.join(_HERE, "results",
                            "_reference_Lx13_Ly16_Tx-3_Ty3",
                            "two_point_all_to_all.dat"),
            3: os.path.join(_HERE, "results",
                            "_reference_Lx39_Ly48_Tx-9_Ty9",
                            "two_point_all_to_all.dat"),
        }
        jobs.append(("ref", alpha, Lx, Ly, Tx, Ty,
                     1.0, 1.0, 1.0,            # k1,k2,k3
                     0.24, 0.29,               # β_c seeds
                     50000,                    # n_traj_prod
                     out_dir, out_dat,
                     src.get(alpha)))          # pre-existing source (or None)

    # ---- TEST: untwisted (L,L,0,0) at truth couplings -----------------------
    for L in [8, 16, 24, 32]:
        out_dir = os.path.join(_OUT, "test", f"L{L}")
        out_dat = os.path.join(out_dir, "two_point_all_to_all.dat")
        jobs.append(("test", L, L, L, 0, 0,
                     _R1, _R2, _K3,
                     0.03, 0.12,               # β_c seeds
                     100000 if L <= 16 else 50000,
                     out_dir, out_dat,
                     None))
    return jobs


# ──────────────────────────────────────────────────────────────────────────────
# Workers
# ──────────────────────────────────────────────────────────────────────────────

def _symlink_or_copy(src, dst_dir, out_dat):
    """Copy or symlink pre-existing all-to-all file into the output dir."""
    os.makedirs(dst_dir, exist_ok=True)
    if os.path.exists(out_dat):
        return ("skip", {})
    if not os.path.exists(src):
        return ("err", {"msg": f"src missing: {src}"})
    shutil.copy2(src, out_dat)
    return ("ok", {"src": src})


def _run_one(args):
    (kind, label, Lx, Ly, Tx, Ty,
     k1, k2, k3,
     beta_lo, beta_hi, n_traj_prod,
     out_dir, out_dat, preexisting_src) = args

    if os.path.exists(out_dat):
        return (kind, label, "skip", {})

    # Pre-existing data: just copy
    if preexisting_src is not None:
        status, info = _symlink_or_copy(preexisting_src, out_dir, out_dat)
        # also copy meta.json if present
        src_meta = os.path.join(os.path.dirname(preexisting_src),
                                "reference_meta.json")
        if os.path.exists(src_meta):
            shutil.copy2(src_meta,
                         os.path.join(out_dir, "meta.json"))
        return (kind, label, status, info)

    os.makedirs(out_dir, exist_ok=True)
    scratch = os.path.join(out_dir, "_scratch")
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()

    try:
        beta_c, beta_c_sigma, chi_peak, *_ = mc_engine.find_beta_c(
            _EXE, Lx, Ly, Tx, Ty, k1, k2, k3,
            beta_lo, beta_hi,
            n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
            n_traj_coarse=3000, n_traj_fine=8000,
            max_shifts=6, jackknife=False,
            data_dir=os.path.join(scratch, "scan"),
        )
        _, subdir = mc_engine.run_simulator(
            _EXE, Lx, Ly, Tx, Ty, k1, k2, k3, beta_c,
            n_traj=n_traj_prod,
            data_dir=os.path.join(scratch, "prod"),
        )
        if subdir is None:
            return (kind, label, "err", {"msg": "no MC subdir"})
        src_dat = os.path.join(subdir, "two_point_all_to_all.dat")
        shutil.copy2(src_dat, out_dat)
        meta = {
            "kind": kind, "label": str(label),
            "Lx": int(Lx), "Ly": int(Ly), "Tx": int(Tx), "Ty": int(Ty),
            "k1": float(k1), "k2": float(k2), "k3": float(k3),
            "beta_c": float(beta_c), "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "n_traj_prod": int(n_traj_prod),
            "wall_s": time.time() - t0,
        }
        with open(os.path.join(out_dir, "meta.json"), "w") as f:
            json.dump(meta, f, indent=2)
        shutil.rmtree(scratch, ignore_errors=True)
        return (kind, label, "ok",
                {"beta_c": beta_c, "wall": time.time() - t0})
    except Exception as e:
        return (kind, label, "err", {"msg": str(e)})


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=4)
    args = ap.parse_args()

    jobs = _jobs()
    print(f"Total jobs: {len(jobs)}")
    t_total = time.time()

    with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
        futs = {pool.submit(_run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            kind, label, status, info = fut.result()
            if status == "ok":
                bc = info.get("beta_c", "?")
                wall = info.get("wall", 0)
                print(f"  OK  {kind}/{label}  beta_c={bc:.5f}  wall={wall:.1f}s")
            elif status == "skip":
                print(f"  --  {kind}/{label}  (already exists)")
            else:
                print(f"  ERR {kind}/{label}  {info}")

    print(f"\nDone.  Total wall: {time.time()-t_total:.1f}s")
    print(f"Output: {_OUT}")


if __name__ == "__main__":
    main()
