#!/usr/bin/env python3
"""
precompute_456_fss_extend.py  —  Additional MC data for free-omega FSS fits.

Adds to results/_fss_456/:
  REF  α=4: (52,64,-12,12), iso,  50k traj
  TEST L=48: (48,48,0,0),   truth couplings,  50k traj
  TEST L=64: (64,64,0,0),   truth couplings,  30k traj  (larger lattice)

With these additions each family has 4 data points, giving 1 dof for a
3-parameter free-omega fit:  G(L) = G_inf + a * L^(-omega).
"""
from __future__ import annotations
import json, math, os, shutil, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import mc_engine

_R1   = (2*math.log(5) - math.log(7)) / (2*math.log(3) - math.log(7))
_R2   = math.log(7) / (2*math.log(3) - math.log(7))
_OUT  = os.path.join(_HERE, "results", "_fss_456")
_EXE  = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")

JOBS = [
    # (kind, label, Lx, Ly, Tx, Ty, k1, k2, k3, blo, bhi, n_traj)
    ("ref",  "a4",  52, 64, -12, 12, 1.0, 1.0, 1.0,  0.24, 0.29, 50_000),
    ("test", "L48", 48, 48,   0,  0, _R1, _R2, 1.0,  0.03, 0.12, 50_000),
    ("test", "L64", 64, 64,   0,  0, _R1, _R2, 1.0,  0.03, 0.12, 30_000),
]


def _run_one(args):
    kind, label, Lx, Ly, Tx, Ty, k1, k2, k3, blo, bhi, n_traj = args
    out_dir = os.path.join(_OUT, kind, label)
    out_dat = os.path.join(out_dir, "two_point_all_to_all.dat")
    if os.path.exists(out_dat):
        return kind, label, "skip", {}
    os.makedirs(out_dir, exist_ok=True)
    scratch = os.path.join(out_dir, "_scratch")
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_c, beta_c_sigma, chi_peak, *_ = mc_engine.find_beta_c(
            _EXE, Lx, Ly, Tx, Ty, k1, k2, k3, blo, bhi,
            n_coarse=11, n_refine=5, n_refine2=5, n_refine3=5,
            n_traj_coarse=3000, n_traj_fine=8000,
            max_shifts=6, jackknife=False,
            data_dir=os.path.join(scratch, "scan"),
        )
        _, subdir = mc_engine.run_simulator(
            _EXE, Lx, Ly, Tx, Ty, k1, k2, k3, beta_c,
            n_traj=n_traj,
            data_dir=os.path.join(scratch, "prod"),
        )
        if subdir is None:
            return kind, label, "err", {"msg": "no subdir"}
        shutil.copy2(os.path.join(subdir, "two_point_all_to_all.dat"), out_dat)
        meta = dict(kind=kind, label=label, Lx=int(Lx), Ly=int(Ly),
                    Tx=int(Tx), Ty=int(Ty),
                    k1=float(k1), k2=float(k2), k3=float(k3),
                    beta_c=float(beta_c), beta_c_sigma=float(beta_c_sigma),
                    chi_peak=float(chi_peak), n_traj_prod=int(n_traj),
                    wall_s=time.time()-t0)
        with open(os.path.join(out_dir, "meta.json"), "w") as f:
            json.dump(meta, f, indent=2)
        shutil.rmtree(scratch, ignore_errors=True)
        return kind, label, "ok", {"beta_c": beta_c, "wall": time.time()-t0}
    except Exception as e:
        return kind, label, "err", {"msg": str(e)}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=3)
    args = ap.parse_args()
    print(f"Jobs: {[j[1] for j in JOBS]}  workers={args.n_workers}")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
        futs = {pool.submit(_run_one, j): j for j in JOBS}
        for fut in as_completed(futs):
            kind, label, status, info = fut.result()
            if status == "ok":
                print(f"  OK  {kind}/{label}  beta_c={info['beta_c']:.5f}  wall={info['wall']:.0f}s")
            elif status == "skip":
                print(f"  --  {kind}/{label}  (already exists)")
            else:
                print(f"  ERR {kind}/{label}  {info}")
    print(f"Total wall: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
