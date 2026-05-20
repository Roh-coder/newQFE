#!/usr/bin/env python3
"""Standalone FSS / autocorrelation diagnostic.

Runs the simulator at the symmetric isotropic triangular-Ising critical point
(K1=K2=K3=1, beta = ln(3)/4) for several L values with --dump_traces 1.
Then computes:
  * integrated autocorrelation time tau_int for |m| and E (Madras-Sokal
    automatic windowing),
  * thermalization traces,
  * FSS of the connected two-point function at fractional separations
    t = k/8 along the three lattice cycles.

Outputs go to results/diag_fss_autocorr/{plots,raw}.
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from lib import mc_engine  # noqa: E402
from workflow_common import sample_directional_channels  # noqa: E402

BETA_C_SYM = math.log(3.0) / 4.0  # exact symmetric isotropic triangular Ising

DEFAULT_SIZES = (8, 16, 24, 32, 40, 48, 56)
EXE = os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram")


def run_one(L, n_therm, n_traj, n_skip, scratch_root, seed,
            find_beta=False, beta_finder_cfg=None, beta_fixed=BETA_C_SYM):
    data_dir = os.path.join(scratch_root, f"L{L}")
    os.makedirs(data_dir, exist_ok=True)

    # ---- Determine beta to use for production -------------------------------
    beta_used = float(beta_fixed)
    beta_sigma = 0.0
    chi_peak = float("nan")
    beta_finder_wall = 0.0
    if find_beta:
        bf = beta_finder_cfg or {}
        scan_dir = os.path.join(data_dir, "scan")
        os.makedirs(scan_dir, exist_ok=True)
        t_bf = time.time()
        beta_used, beta_sigma, chi_peak, _sb, _sc, _se = mc_engine.find_beta_c(
            EXE, L, L, 0, 0, 1.0, 1.0, 1.0,
            float(bf.get("beta_lo", 0.10)),
            float(bf.get("beta_hi", 0.45)),
            n_coarse=int(bf.get("n_coarse", 11)),
            n_refine=int(bf.get("n_refine", 5)),
            n_refine2=int(bf.get("n_refine2", 5)),
            n_refine3=int(bf.get("n_refine3", 5)),
            n_traj_coarse=int(bf.get("n_traj_coarse", 120)),
            n_traj_fine=int(bf.get("n_traj_fine", 240)),
            data_dir=scan_dir,
            max_shifts=int(bf.get("max_shifts", 4)),
            jackknife=bool(bf.get("jackknife", False)),
        )
        beta_finder_wall = time.time() - t_bf

    # ---- Production run with traces ----------------------------------------
    prod_dir = os.path.join(data_dir, "prod")
    os.makedirs(prod_dir, exist_ok=True)
    cmd = [
        EXE,
        "--L_x", str(L), "--L_y", str(L),
        "--T_x", "0", "--T_y", "0",
        "--k1", "1.000000000000",
        "--k2", "1.000000000000",
        "--k3", "1.000000000000",
        "--beta", f"{beta_used:.12f}",
        "--n_traj", str(n_traj),
        "--n_skip", str(n_skip),
        "--n_therm", str(n_therm),
        "--seed", str(seed),
        "--data_dir", prod_dir,
        "--dump_traces", "1",
    ]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return {"L": L, "ok": False, "stderr": res.stderr[-2000:]}
    a2a_files = sorted(glob(os.path.join(prod_dir, "*", "two_point_all_to_all.dat")),
                       key=os.path.getmtime)
    trace_files = sorted(glob(os.path.join(prod_dir, "*", "traces_*.dat")),
                         key=os.path.getmtime)
    if not a2a_files or not trace_files:
        return {"L": L, "ok": False, "stderr": "missing outputs"}
    return {
        "L": L,
        "ok": True,
        "wall_s": time.time() - t0,
        "beta_finder_wall_s": beta_finder_wall,
        "beta_used": beta_used,
        "beta_sigma": beta_sigma,
        "chi_peak": chi_peak,
        "all_to_all": a2a_files[-1],
        "traces": trace_files[-1],
    }


def integrated_autocorr(x):
    """Madras-Sokal automatic-windowing tau_int for a 1D series."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 8:
        return np.nan, np.nan
    x = x - x.mean()
    var0 = np.dot(x, x) / n
    if var0 <= 0:
        return np.nan, np.nan
    # Compute autocorrelation via FFT.
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, n=nfft)
    acf = np.fft.irfft(f * np.conj(f), n=nfft)[:n] / (n * var0)
    # Automatic windowing (Sokal): find smallest M with M >= c*tau(M).
    c = 6.0
    tau = 0.5
    M = 1
    while M < n:
        tau = 0.5 + np.sum(acf[1 : M + 1])
        if M >= c * tau:
            break
        M += 1
    # crude error estimate (Wolff)
    err = tau * np.sqrt((2.0 * (2 * M + 1)) / max(n, 1))
    return float(tau), float(err)


def parse_traces(path):
    """Return n_sites, beta, and arrays (E_per_site, abs_m, m2)."""
    n_sites = None
    beta = None
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                # header lines may carry n_sites, beta
                if "n_sites" in s:
                    toks = s.split()
                    for i, t in enumerate(toks):
                        if t == "n_sites":
                            n_sites = int(toks[i + 1])
                        elif t == "beta":
                            beta = float(toks[i + 1])
                continue
            parts = s.split()
            if len(parts) >= 3:
                rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
    arr = np.asarray(rows)
    return n_sites, beta, arr  # arr[:,0]=E_per_site, [:,1]=|m|, [:,2]=m^2


def make_payload(L, all_to_all_file):
    data = mc_engine.load_all_to_all(all_to_all_file)
    return {
        "Lx": L,
        "Ly": L,
        "Tx": 0,
        "Ty": 0,
        "data": data,
        "all_to_all_file": all_to_all_file,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*", default=list(DEFAULT_SIZES))
    ap.add_argument("--n-therm", type=int, default=4000)
    ap.add_argument("--n-traj", type=int, default=8000)
    ap.add_argument("--n-skip", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out-root", default=os.path.join(HERE, "results", "diag_fss_autocorr"))
    ap.add_argument("--find-beta", action="store_true",
                    help="Locate beta_c per L via 3-pass susceptibility scan "
                         "instead of using the analytic value.")
    ap.add_argument("--beta-from-largest-size", action="store_true",
                    help="With --find-beta, locate beta_c once at max(sizes) "
                         "and reuse that same beta for every smaller L.")
    ap.add_argument("--beta-lo", type=float, default=0.10)
    ap.add_argument("--beta-hi", type=float, default=0.45)
    ap.add_argument("--bf-n-coarse", type=int, default=11)
    ap.add_argument("--bf-n-refine", type=int, default=5)
    ap.add_argument("--bf-n-refine2", type=int, default=5)
    ap.add_argument("--bf-n-refine3", type=int, default=5)
    ap.add_argument("--bf-n-traj-coarse", type=int, default=120)
    ap.add_argument("--bf-n-traj-fine", type=int, default=240)
    args = ap.parse_args()

    if args.beta_from_largest_size and not args.find_beta:
        ap.error("--beta-from-largest-size requires --find-beta")

    bf_cfg = {
        "beta_lo": args.beta_lo,
        "beta_hi": args.beta_hi,
        "n_coarse": args.bf_n_coarse,
        "n_refine": args.bf_n_refine,
        "n_refine2": args.bf_n_refine2,
        "n_refine3": args.bf_n_refine3,
        "n_traj_coarse": args.bf_n_traj_coarse,
        "n_traj_fine": args.bf_n_traj_fine,
    }

    out_root = args.out_root
    plots_dir = os.path.join(out_root, "plots")
    raw_dir = os.path.join(out_root, "raw")
    scratch_root = os.path.join(out_root, "_scratch")
    for d in (plots_dir, raw_dir, scratch_root):
        os.makedirs(d, exist_ok=True)

    print(f"[diag] beta_c (sym triang Ising) = {BETA_C_SYM:.10f}")
    print(f"[diag] sizes = {args.sizes}, n_therm={args.n_therm}, "
          f"n_traj={args.n_traj}, n_skip={args.n_skip}, workers={args.workers}")
    if args.find_beta:
        beta_mode = "largest-size shared" if args.beta_from_largest_size else "per-size"
        print(f"[diag] find_beta=ON ({beta_mode}), scan range=[{args.beta_lo},{args.beta_hi}], "
              f"bf_traj=({args.bf_n_traj_coarse}/{args.bf_n_traj_fine}), "
              f"bf_n=({args.bf_n_coarse}/{args.bf_n_refine}/{args.bf_n_refine2}/{args.bf_n_refine3})")
    else:
        print("[diag] find_beta=OFF (using analytic beta_c)")

    shared_beta_info = None
    shared_beta_L = None
    beta_plot_label = r"$\beta_c=\ln 3/4$"
    if args.find_beta and args.beta_from_largest_size:
        shared_beta_L = max(args.sizes)
        shared_scan_dir = os.path.join(scratch_root, f"_shared_beta_L{shared_beta_L}")
        os.makedirs(shared_scan_dir, exist_ok=True)
        print(f"[diag] locating shared beta_c at largest size L={shared_beta_L}")
        t_bf = time.time()
        beta_used, beta_sigma, chi_peak, _sb, _sc, _se = mc_engine.find_beta_c(
            EXE, shared_beta_L, shared_beta_L, 0, 0, 1.0, 1.0, 1.0,
            float(bf_cfg.get("beta_lo", 0.10)),
            float(bf_cfg.get("beta_hi", 0.45)),
            n_coarse=int(bf_cfg.get("n_coarse", 11)),
            n_refine=int(bf_cfg.get("n_refine", 5)),
            n_refine2=int(bf_cfg.get("n_refine2", 5)),
            n_refine3=int(bf_cfg.get("n_refine3", 5)),
            n_traj_coarse=int(bf_cfg.get("n_traj_coarse", 120)),
            n_traj_fine=int(bf_cfg.get("n_traj_fine", 240)),
            data_dir=shared_scan_dir,
            max_shifts=int(bf_cfg.get("max_shifts", 4)),
            jackknife=bool(bf_cfg.get("jackknife", False)),
        )
        shared_beta_info = {
            "beta_used": beta_used,
            "beta_sigma": beta_sigma,
            "chi_peak": chi_peak,
            "beta_finder_wall_s": time.time() - t_bf,
        }
        beta_plot_label = (rf"$\beta_c(L={shared_beta_L})={beta_used:.6f}$"
                           if np.isfinite(beta_used)
                           else rf"shared $\beta_c(L={shared_beta_L})$")
        print(f"[diag] shared beta_c from L={shared_beta_L}: "
              f"{beta_used:.6f}+-{beta_sigma:.6f} "
              f"(bf_wall={shared_beta_info['beta_finder_wall_s']:.1f}s)")
    elif args.find_beta:
        beta_plot_label = r"$\beta_c(L)$ from per-size susceptibility scans"

    # 1. Run simulations in parallel.
    sim_results = {}
    seed_base = int(time.time()) & 0xFFFF
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(run_one, L, args.n_therm, args.n_traj, args.n_skip,
                      scratch_root, seed_base + i + 1,
                      args.find_beta and not args.beta_from_largest_size,
                      bf_cfg,
                      shared_beta_info["beta_used"] if shared_beta_info is not None else BETA_C_SYM): L
            for i, L in enumerate(args.sizes)
        }
        for fut in as_completed(futs):
            r = fut.result()
            L = r["L"]
            if shared_beta_info is not None:
                r["beta_used"] = float(shared_beta_info["beta_used"])
                r["beta_sigma"] = float(shared_beta_info["beta_sigma"])
                r["chi_peak"] = float(shared_beta_info["chi_peak"])
                r["beta_finder_wall_s"] = (
                    float(shared_beta_info["beta_finder_wall_s"])
                    if L == shared_beta_L else 0.0
                )
            sim_results[L] = r
            if r["ok"]:
                if args.find_beta and shared_beta_info is not None:
                    bf_msg = (f" beta_c={r['beta_used']:.6f}+-{r['beta_sigma']:.6f}"
                              f" (shared from L={shared_beta_L}"
                              f", bf_wall={shared_beta_info['beta_finder_wall_s']:.1f}s)")
                else:
                    bf_msg = (f" beta_c={r['beta_used']:.6f}+-{r['beta_sigma']:.6f}"
                              f" (bf_wall={r['beta_finder_wall_s']:.1f}s)"
                              if args.find_beta else "")
                print(f"[diag] L={L:>3}  ok  prod_wall={r['wall_s']:.1f}s{bf_msg}")
            else:
                print(f"[diag] L={L:>3}  ERR  {r.get('stderr','')[:200]}")

    sizes_ok = sorted(L for L, r in sim_results.items() if r["ok"])
    if not sizes_ok:
        print("[diag] no successful runs"); return

    # 2. Per-L parse traces + autocorrelation; collect FSS data.
    fractions = np.array([k / 8.0 for k in range(1, 8)])
    summary_rows = []
    traces_by_L = {}
    G_by_L = {}  # G[L] shape (3, 7)
    for L in sizes_ok:
        r = sim_results[L]
        n_sites, beta_run, traces = parse_traces(r["traces"])
        E_series = traces[:, 0]
        absm_series = traces[:, 1]
        m2_series = traces[:, 2]
        tau_E, tau_E_err = integrated_autocorr(E_series)
        tau_m, tau_m_err = integrated_autocorr(absm_series)
        # Thermalization assessment: split n_therm trajectories already removed
        # by simulator; here we look at the production series stability via
        # a 5-block mean.
        nb = 5
        bs = len(absm_series) // nb
        block_means = [absm_series[i * bs : (i + 1) * bs].mean() for i in range(nb)]
        therm_drift = (block_means[-1] - block_means[0]) / np.std(absm_series, ddof=1)
        traces_by_L[L] = {
            "E": E_series, "absm": absm_series, "m2": m2_series,
            "tau_E": tau_E, "tau_E_err": tau_E_err,
            "tau_m": tau_m, "tau_m_err": tau_m_err,
            "therm_drift_sigmas": float(therm_drift),
        }

        payload = make_payload(L, r["all_to_all"])
        G, sG = sample_directional_channels(payload, fractions)
        G_by_L[L] = (G, sG)

        summary_rows.append((L, len(absm_series), tau_E, tau_m,
                             float(np.mean(absm_series)),
                             float(np.std(absm_series, ddof=1)),
                             float(therm_drift)))
        print(f"[diag] L={L:>3}  n_meas={len(absm_series):>5}  "
              f"tau_int(|m|)={tau_m:5.2f}+-{tau_m_err:.2f}  "
              f"tau_int(E)={tau_E:5.2f}+-{tau_E_err:.2f}  "
              f"therm_drift={therm_drift:+.2f} sigma")

    # 3. Save summary table.
    sumpath = os.path.join(raw_dir, "summary.dat")
    with open(sumpath, "w") as f:
        f.write("# L  n_meas  tau_int_E  tau_int_|m|  mean_|m|  sd_|m|  "
                "therm_drift_sigmas  beta_used  beta_sigma\n")
        for L, row in zip(sizes_ok, summary_rows):
            r = sim_results[L]
            extra = (float(r.get("beta_used", float("nan"))),
                     float(r.get("beta_sigma", 0.0)))
            full = list(row) + list(extra)
            f.write(" ".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in full) + "\n")
    print(f"[diag] wrote {sumpath}")

    # 4. Plots ----------------------------------------------------------------
    # 4a. Thermalization (|m| vs trajectory index) overlaid per L.
    fig, ax = plt.subplots(figsize=(9, 5))
    for L in sizes_ok:
        s = traces_by_L[L]["absm"]
        ax.plot(np.arange(len(s)), s, lw=0.6, alpha=0.7, label=f"L={L}")
    ax.set_xlabel("measured trajectory index (post-thermalization)")
    ax.set_ylabel("|m|")
    ax.set_title("Thermalization traces of |m|")
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "thermalization_absm.png"), dpi=120)
    plt.close(fig)

    # 4b. tau_int vs L
    Ls = np.array(sizes_ok)
    tau_m = np.array([traces_by_L[L]["tau_m"] for L in sizes_ok])
    tau_m_e = np.array([traces_by_L[L]["tau_m_err"] for L in sizes_ok])
    tau_E = np.array([traces_by_L[L]["tau_E"] for L in sizes_ok])
    tau_E_e = np.array([traces_by_L[L]["tau_E_err"] for L in sizes_ok])
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(Ls, tau_m, yerr=tau_m_e, fmt="o-", label=r"$\tau_{\rm int}(|m|)$")
    ax.errorbar(Ls, tau_E, yerr=tau_E_e, fmt="s-", label=r"$\tau_{\rm int}(E)$")
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\tau_{\rm int}$  (in units of measured sweeps, n_skip="
                  f"{args.n_skip})")
    ax.set_title(rf"Integrated autocorrelation time at {beta_plot_label}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "tau_int_vs_L.png"), dpi=120)
    plt.close(fig)

    # 4c. FSS panels: G vs 1/L for each (cycle, k). 3 rows x 7 cols.
    fig, axes = plt.subplots(3, 7, figsize=(20, 8.5), sharex=True)
    inv_L = 1.0 / Ls
    for cyc in range(3):
        for ki, k in enumerate(range(1, 8)):
            ax = axes[cyc, ki]
            ys = []
            es = []
            for L in sizes_ok:
                G, sG = G_by_L[L]
                ys.append(G[cyc, ki])
                es.append(sG[cyc, ki])
            ys = np.array(ys); es = np.array(es)
            ax.errorbar(inv_L, ys, yerr=es, fmt="o-", lw=1)
            ax.set_title(f"cycle={cyc} k={k} t={k/8:.3f}", fontsize=9)
            ax.grid(alpha=0.3)
            if cyc == 2:
                ax.set_xlabel("1/L")
            if ki == 0:
                ax.set_ylabel("G(t)")
    fig.suptitle(rf"FSS of two-point function (sym triang Ising, {beta_plot_label},"
                 f" K=(1,1,1)); n_therm={args.n_therm}, n_traj={args.n_traj}, "
                 f"n_skip={args.n_skip}")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(plots_dir, "fss_two_point.png"), dpi=120)
    plt.close(fig)

    # 4d. Save per-L raw G, sG.
    with open(os.path.join(raw_dir, "fss_G.dat"), "w") as f:
        f.write("# L cycle k t G sG\n")
        for L in sizes_ok:
            G, sG = G_by_L[L]
            for cyc in range(3):
                for ki, k in enumerate(range(1, 8)):
                    f.write(f"{L} {cyc} {k} {k/8:.6f} {G[cyc,ki]:.10g} {sG[cyc,ki]:.10g}\n")

    print(f"[diag] plots: {plots_dir}")
    print(f"[diag] raw  : {raw_dir}")


if __name__ == "__main__":
    main()
