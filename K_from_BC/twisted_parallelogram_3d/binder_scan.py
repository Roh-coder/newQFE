#!/usr/bin/env python3
"""
Binder cumulant scan for the 3D 4-5-6 triangle × N_t Ising model.

11 log-spaced K values centred on K_c ≈ 0.165, range [0.065, 0.265].
N_t = 8, 16, 24, 32 with spatial geometry scale=1 (Lx=13, Ly=16, Tx=3, Ty=-3).

Usage:
    python binder_scan.py          # run sims + plot
    python binder_scan.py --plot   # re-plot from existing data
"""

import argparse
import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(SCRIPT_DIR, "ising_tri_twisted_parallelogram")
DATA_DIR = os.path.join(SCRIPT_DIR, "binder_scan_out")

# Spatial geometry — 4-5-6 triangle, scale=1 → N_cell = 13*16 + 3*3 = 217
LX, LY, TX, TY = 13, 16, 3, -3

# 11 log-spaced K values centred on 0.165 with range ±0.1
K_VALUES = np.logspace(np.log10(0.065), np.log10(0.265), 11)
NT_VALUES = [8, 16, 24, 32]

# MC parameters
N_THERM = 4000
N_TRAJ = 40000
N_SKIP = 5
N_WOLFF = 3
N_METROPOLIS = 5
SEED = 31415926

MAX_WORKERS = 4  # parallel sims

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_id_str(nt: int, k: float) -> str:
    """Match the C sprintf format: Lx%d_Ly%d_Tx%d_Ty%d_Nt%d_k%.3f_%.3f_%.3f_kt%.3f"""
    return (f"Lx{LX}_Ly{LY}_Tx{TX}_Ty{TY}_Nt{nt}"
            f"_k{k:.3f}_{k:.3f}_{k:.3f}_kt{k:.3f}")


def find_dat(nt: int, k: float) -> str | None:
    """Glob for the .dat file produced by the binary for this (nt, k) pair."""
    pattern = os.path.join(DATA_DIR, run_id_str(nt, k), "*.dat")
    hits = glob.glob(pattern)
    return hits[0] if hits else None


def run_sim(nt: int, k: float) -> bool:
    if find_dat(nt, k) is not None:
        print(f"  [skip] Nt={nt:2d}  K={k:.4f}")
        return True

    print(f"  [run ] Nt={nt:2d}  K={k:.4f}")
    cmd = [
        BINARY,
        "--L_x", str(LX), "--L_y", str(LY),
        "--T_x", str(TX), "--T_y", str(TY),
        "--N_t", str(nt),
        "--k1",  f"{k:.8f}",
        "--k2",  f"{k:.8f}",
        "--k3",  f"{k:.8f}",
        "--k_t", f"{k:.8f}",
        "--beta", "1.0",
        "--seed", str(SEED),
        "--n_therm", str(N_THERM),
        "--n_traj",  str(N_TRAJ),
        "--n_skip",  str(N_SKIP),
        "--n_wolff", str(N_WOLFF),
        "--n_metropolis", str(N_METROPOLIS),
        "--data_dir", DATA_DIR,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  [FAIL] Nt={nt}  K={k:.4f}\n{res.stderr[:400]}")
        return False
    return True


def parse_dat(path: str) -> dict:
    data = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                key = parts[0]
                try:
                    data[key] = (float(parts[1]), float(parts[2]))
                except ValueError:
                    pass
    return data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true",
                        help="Skip simulations; just re-plot existing data.")
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    args = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)

    # ------------------------------------------------------------------
    # Run simulations
    # ------------------------------------------------------------------
    if not args.plot:
        jobs = [(nt, k) for nt in NT_VALUES for k in K_VALUES]
        n_needed = sum(1 for nt, k in jobs if find_dat(nt, k) is None)
        print(f"Submitting {len(jobs)} runs ({n_needed} not yet done) "
              f"with up to {args.workers} workers …")
        print(f"K values: {' '.join(f'{k:.4f}' for k in K_VALUES)}\n")

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_sim, nt, k): (nt, k) for nt, k in jobs}
            for fut in as_completed(futures):
                nt, k = futures[fut]
                if not fut.result():
                    print(f"  WARNING: failed  Nt={nt}  K={k:.4f}")
        print()

    # ------------------------------------------------------------------
    # Collect results
    # ------------------------------------------------------------------
    results = {}
    for nt in NT_VALUES:
        ks, u4s, u4_errs = [], [], []
        for k in K_VALUES:
            path = find_dat(nt, k)
            if path is None:
                print(f"  Missing: Nt={nt}  K={k:.4f}")
                continue
            d = parse_dat(path)
            if "U4" not in d:
                print(f"  No U4 in {path}")
                continue
            u4_mean, u4_err = d["U4"]
            ks.append(k)
            u4s.append(u4_mean)
            u4_errs.append(u4_err)
        results[nt] = (np.array(ks), np.array(u4s), np.array(u4_errs))

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(NT_VALUES)))

    for i, nt in enumerate(NT_VALUES):
        ks, u4s, errs = results[nt]
        if len(ks) == 0:
            continue
        ax.errorbar(ks, u4s, yerr=errs,
                    marker="o", ms=5, linestyle="-", linewidth=1.5,
                    color=colors[i], label=f"$N_t = {nt}$", capsize=3, capthick=1)

    ax.axvline(0.165, color="gray", linestyle="--", linewidth=1, alpha=0.7,
               label=r"$K_c^{\rm guess}=0.165$")

    ax.set_xscale("log")
    ax.set_xlabel(r"$K$", fontsize=14)
    ax.set_ylabel(r"$U_4 = \tfrac{3}{2}\!\left(1 - \tfrac{\langle m^4\rangle}"
                  r"{3\langle m^2\rangle^2}\right)$", fontsize=13)
    ax.set_title(
        r"Binder cumulant — 4-5-6 triangle $\times\, N_t$ (3D Ising)" "\n"
        rf"Lx={LX}, Ly={LY}, Tx={TX}, Ty={TY},  $N_{{\rm cell}}=217$",
        fontsize=12)
    ax.legend(fontsize=11, loc="upper left")
    ax.grid(True, alpha=0.25, which="both")

    out = os.path.join(SCRIPT_DIR, "binder_crossing_scan.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")

    # Also an inset zoom around the crossing
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    k_lo, k_hi = 0.12, 0.22
    for i, nt in enumerate(NT_VALUES):
        ks, u4s, errs = results[nt]
        if len(ks) == 0:
            continue
        mask = (ks >= k_lo) & (ks <= k_hi)
        if mask.sum() < 2:
            continue
        ax2.errorbar(ks[mask], u4s[mask], yerr=errs[mask],
                     marker="o", ms=6, linestyle="-", linewidth=2,
                     color=colors[i], label=f"$N_t = {nt}$", capsize=4, capthick=1.5)

    ax2.axvline(0.165, color="gray", linestyle="--", linewidth=1.2, alpha=0.7,
                label=r"$K_c^{\rm guess}=0.165$")
    ax2.set_xlabel(r"$K$", fontsize=14)
    ax2.set_ylabel(r"$U_4$", fontsize=14)
    ax2.set_title(
        r"Binder cumulant — crossing region zoom" "\n"
        rf"Lx={LX}, Ly={LY}, Tx={TX}, Ty={TY},  $N_{{\rm cell}}=217$",
        fontsize=12)
    ax2.legend(fontsize=11, loc="upper left")
    ax2.grid(True, alpha=0.3)

    out2 = os.path.join(SCRIPT_DIR, "binder_crossing_zoom.png")
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved: {out2}")
    plt.close("all")


if __name__ == "__main__":
    main()
