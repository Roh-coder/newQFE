#!/usr/bin/env python3
"""
run_binder_scan.py
==================
Binder-crossing scan for the 3D Ising model on the triangular-prism lattice.

Geometry: Lx=Ly=Nt=L, Tx=Ty=0 (equilateral triangular torus x time), isotropic K.
11 log-spaced K values centred on Kc=0.165, covering [0.065, 0.265].
Four lattice sizes: L = 8, 16, 24, 32.

Job ordering: largest L first; within each L, K nearest Kc first (worst
critical slowing-down first) so the hardest jobs start immediately.

Outputs:
  out_binder_3d/progress.log   -- timestamped per-job START/DONE lines
  out_binder_3d/summary.tsv    -- U4±err, m_susc±err per (L,K) as they finish
  out_binder_3d/log_L??_iK??.txt  -- full stdout from each run

Monitor with:  tail -f out_binder_3d/progress.log
Re-plot with:  python plot_binder_crossing.py
Skips existing completed .dat files automatically.
"""

import os
import math
import glob
import subprocess
import concurrent.futures
import sys
import time
import threading
import argparse

# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BINARY       = os.path.join(SCRIPT_DIR, "ising_tri_twisted_parallelogram")
DATA_DIR     = os.path.join(SCRIPT_DIR, "out_binder_3d")
PROGRESS_LOG = ""
SUMMARY_TSV  = ""

KC_CENTER    = 0.165
K_MIN, K_MAX, N_K = 0.065, 0.265, 11


def build_k_values(k_min, k_max, n_k, mode):
    if n_k < 2:
        raise ValueError("n_k must be >= 2")
    if k_min <= 0.0 or k_max <= 0.0 or k_max <= k_min:
        raise ValueError("Require 0 < k_min < k_max")
    if mode == "log":
        return [
            math.exp(math.log(k_min) + i * (math.log(k_max) - math.log(k_min)) / (n_k - 1))
            for i in range(n_k)
        ]
    if mode == "linear":
        return [k_min + i * (k_max - k_min) / (n_k - 1) for i in range(n_k)]
    raise ValueError(f"Unknown k spacing mode: {mode}")

# Quick-start defaults for this batch; can be overridden via CLI.
DEFAULT_N_TRAJ = 1000
DEFAULT_N_THERM = 300
DEFAULT_N_SKIP = 10

# Ordered largest L first.
L_VALUES = [32, 24, 16, 8]

MAX_WORKERS = 8

# ---------------------------------------------------------------------------
_log_lock     = threading.Lock()
_summary_lock = threading.Lock()
_t0           = None
_done_count   = 0
_done_lock    = threading.Lock()
_total_jobs   = 0


def _log(msg):
    ts   = time.strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with _log_lock:
        with open(PROGRESS_LOG, "a") as f:
            f.write(line + "\n")


def _append_summary(L, K, u4, u4e, ms, mse, dat_path):
    rel  = os.path.relpath(dat_path, os.path.dirname(DATA_DIR))
    line = f"{L}\t{K:.10f}\t{u4:.12e}\t{u4e:.12e}\t{ms:.12e}\t{mse:.12e}\t{rel}\n"
    with _summary_lock:
        with open(SUMMARY_TSV, "a") as f:
            f.write(line)


def _dat_exists(L, Nt, K):
    run_id = (f"Lx{L}_Ly{L}_Tx0_Ty0_Nt{Nt}"
              f"_k{K:.3f}_{K:.3f}_{K:.3f}_kt{K:.3f}")
    hits = glob.glob(os.path.join(DATA_DIR, run_id, "*.dat"))
    # Multiple nearby K values can share the same rounded run_id directory.
    # Match by exact K1 stored in file content to avoid cross-K collisions.
    tol = 5e-11
    candidates = []
    for h in sorted(hits):
        base = os.path.basename(h)
        if base.startswith(run_id + "_") and base.endswith(".dat"):
            candidates.append(h)

    for h in candidates:
        d = _parse_dat(h)
        try:
            k1 = float(d["K1"][0])
        except (KeyError, IndexError, ValueError):
            continue
        if abs(k1 - K) <= tol:
            return h

    return None


def _parse_dat(path):
    d = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        d[parts[0]] = [float(x) for x in parts[1:]]
                    except ValueError:
                        pass
    except OSError:
        pass
    return d


def _record_dat(L, K, dat):
    d = _parse_dat(dat)
    try:
        u4  = d["U4"][0];       u4e = d["U4"][1]
        ms  = d["m_susc"][0];   mse = d["m_susc"][1]
        _append_summary(L, K, u4, u4e, ms, mse, dat)
        return u4, u4e
    except (KeyError, IndexError):
        return None, None


def run_one(args):
    global _done_count
    L, Nt, iK, K, n_traj, n_therm, n_skip = args

    existing = _dat_exists(L, Nt, K)
    if existing is not None:
        _log(f"SKIP  L={L:2d}  Nt={Nt:3d}  K={K:.6f}  -> {os.path.basename(existing)}")
        _record_dat(L, K, existing)
        with _done_lock:
            _done_count += 1
        return (L, K, 0, existing)

    seed  = ((L * 997 + iK * 31337) ^ 0xDEADBEEF) & 0x7FFFFFFF
    Kstr  = f"{K:.10f}"
    cmd   = [
        BINARY,
        f"--L_x={L}",   f"--L_y={L}",
        "--T_x=0",      "--T_y=0",
        f"--N_t={Nt}",
        f"--k1={Kstr}",  f"--k2={Kstr}",
        f"--k3={Kstr}",  f"--k_t={Kstr}",
        f"--n_therm={n_therm}",
        f"--n_traj={n_traj}",
        f"--n_skip={n_skip}",
        "--n_wolff=3",
        "--n_metropolis=5",
        f"--seed={seed}",
        f"--data_dir={DATA_DIR}",
    ]

    _log(f"START L={L:2d}  Nt={Nt:3d}  K={K:.6f}  seed=0x{seed:08X}  "
         f"n_sites={L*L*Nt}  n_traj={n_traj}")
    t_start  = time.time()
    log_path = os.path.join(DATA_DIR, f"log_L{L:02d}_iK{iK:02d}.txt")
    with open(log_path, "w") as lf:
        lf.write("CMD: " + " ".join(cmd) + "\n\n")
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t_start
    status  = "OK" if result.returncode == 0 else f"FAIL(rc={result.returncode})"

    u4_str = ""
    dat = _dat_exists(L, Nt, K)
    if dat and result.returncode == 0:
        u4, u4e = _record_dat(L, K, dat)
        if u4 is not None:
            u4_str = f"  U4={u4:.4f}±{u4e:.4f}"

    with _done_lock:
        _done_count += 1
        dc = _done_count

    elapsed_total = time.time() - _t0
    remaining     = _total_jobs - dc
    eta_s = (elapsed_total / dc) * remaining if dc > 0 else 0
    _log(f"DONE  L={L:2d}  Nt={Nt:3d}  K={K:.6f}  {status}  "
         f"t={elapsed:.0f}s{u4_str}")
    _log(f"PROGRESS {dc}/{_total_jobs}  remaining={remaining}  "
         f"ETA ~{eta_s/60:.1f} min")
    return (L, K, result.returncode, log_path)


def main():
    global _t0, _total_jobs, DATA_DIR, PROGRESS_LOG, SUMMARY_TSV

    ap = argparse.ArgumentParser(description="Run Binder scan with live progress logs.")
    ap.add_argument("--data-dir", default=DATA_DIR,
                    help="output directory for logs, summary, and run subdirectories")
    ap.add_argument("--n-traj", type=int, default=DEFAULT_N_TRAJ,
                    help="production trajectories per job")
    ap.add_argument("--n-therm", type=int, default=DEFAULT_N_THERM,
                    help="thermalization trajectories per job")
    ap.add_argument("--n-skip", type=int, default=DEFAULT_N_SKIP,
                    help="measurement skip interval")
    ap.add_argument("--k-min", type=float, default=K_MIN,
                    help="minimum K value for scan")
    ap.add_argument("--k-max", type=float, default=K_MAX,
                    help="maximum K value for scan")
    ap.add_argument("--n-k", type=int, default=N_K,
                    help="number of K points")
    ap.add_argument("--k-spacing", choices=["log", "linear"], default="log",
                    help="K grid spacing")
    ap.add_argument("--nt-factor", type=int, default=1,
                    help="set Nt = nt_factor * L (default: 1)")
    args = ap.parse_args()

    if args.nt_factor < 1:
        raise ValueError("--nt-factor must be >= 1")

    K_VALUES = build_k_values(args.k_min, args.k_max, args.n_k, args.k_spacing)

    DATA_DIR = os.path.abspath(args.data_dir)
    PROGRESS_LOG = os.path.join(DATA_DIR, "progress.log")
    SUMMARY_TSV = os.path.join(DATA_DIR, "summary.tsv")

    if not os.path.isfile(BINARY):
        print(f"ERROR: binary not found: {BINARY}")
        sys.exit(1)

    os.makedirs(DATA_DIR, exist_ok=True)

    # Write summary header once (don't clobber existing)
    if not os.path.exists(SUMMARY_TSV):
        with open(SUMMARY_TSV, "w") as f:
            f.write("L\tK\tU4\tU4_err\tm_susc\tm_susc_err\tdat\n")

    # Build and sort jobs: L descending, then |log(K/Kc)| ascending (Kc-nearest first)
    all_jobs = [
        (L, args.nt_factor * L, iK, K, args.n_traj, args.n_therm, args.n_skip)
        for L in L_VALUES
        for iK, K in enumerate(K_VALUES)
    ]
    all_jobs.sort(key=lambda j: (-j[0], abs(math.log(j[3] / KC_CENTER))))

    _total_jobs = len(all_jobs)
    _t0         = time.time()

    _log(f"Binary    : {BINARY}")
    _log(f"Data dir  : {DATA_DIR}")
    _log(f"n_traj    : {args.n_traj}")
    _log(f"n_therm   : {args.n_therm}")
    _log(f"n_skip    : {args.n_skip}")
    _log(f"nt_factor : {args.nt_factor}  (Nt = nt_factor * L)")
    _log(f"k_spacing : {args.k_spacing}")
    _log(f"K values  : " + " ".join(f"{k:.4f}" for k in K_VALUES))
    _log(f"Sizes     : {L_VALUES}")
    _log(f"Total jobs: {_total_jobs}   max_workers={MAX_WORKERS}")
    _log(f"First 8   : " + "  ".join(f"L={j[0]}/Nt={j[1]}/K={j[3]:.3f}" for j in all_jobs[:8]))
    _log(f"Monitor   : tail -f {PROGRESS_LOG}")
    _log("")

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        # Submit in sorted order so the first MAX_WORKERS picked up are the costliest
        futs = [ex.submit(run_one, job) for job in all_jobs]
        for fut in concurrent.futures.as_completed(futs):
            fut.result()   # surface exceptions if any

    _log("All jobs complete.")
    _log(f"Summary : {SUMMARY_TSV}")
    _log(f"Plot    : python plot_binder_crossing.py")


if __name__ == "__main__":
    main()
