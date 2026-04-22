"""
run_benchmark.py — drive one or more optimizers against the same reference run.

Workflow:
  1. Build (or locate) the reference correlator: an 8x8 equilateral
     (k1=k2=k3=1) untwisted run at its measured beta_c.  Cached in
     `results/_reference/`.
  2. For each optimizer method, instantiate an Evaluator wired to a fresh
     OptimizerPlotter, then call `run_<method>(evaluator, ...)`.
  3. Write a top-level `results/summary.json` comparing all methods.

Usage:
    python run_benchmark.py --method all
    python run_benchmark.py --method nelder_mead --max-evals 30
    python run_benchmark.py --method gp --no-vis     # skip plotting
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import mc_engine
from evaluator import Evaluator
import optimizers as opt_mod


HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EXE = os.path.join(HERE, "bin", "ising_tri_twisted_parallelogram")
# Default results directory — tests should override with --results-dir to
# keep runs isolated without touching each other's data.
DEFAULT_RESULTS = os.path.join(HERE, "results")


# ---------------------------------------------------------------------------
#  Reference run
# ---------------------------------------------------------------------------

def build_reference(exe, L=8, n_traj=60000,
                    beta_lo=0.20, beta_hi=0.32,
                    out_dir=None, results_root=None):
    """Locate beta_c for an LxL equilateral, run a long MC there, cache it."""
    if out_dir is None:
        root = results_root or DEFAULT_RESULTS
        out_dir = os.path.join(root, "_reference")
    os.makedirs(out_dir, exist_ok=True)
    cache_a2a = os.path.join(out_dir, "two_point_all_to_all.dat")
    cache_meta = os.path.join(out_dir, "reference_meta.json")
    cache_scan = os.path.join(out_dir, "reference_betac_scan.json")
    if os.path.exists(cache_a2a) and os.path.exists(cache_meta):
        with open(cache_meta) as f:
            meta = json.load(f)
        print(f"[ref] using cached reference at {cache_a2a} "
              f"(β_c={meta['beta_c']:.6f})")
        data = mc_engine.load_all_to_all(cache_a2a)
        return data, meta

    print(f"[ref] building {L}x{L} equilateral reference …")
    scratch = os.path.join(out_dir, "_scratch")
    beta_c, chi_peak, ref_sb, ref_sc, ref_sce = mc_engine.find_beta_c(
        exe, L, L, 0, 0, 1.0, 1.0, 1.0,
        beta_lo, beta_hi,
        n_traj_coarse=4000, n_traj_fine=10000,
        data_dir=os.path.join(scratch, "scan"),
    )
    print(f"[ref]  β_c = {beta_c:.6f}  χ_peak={chi_peak:.3g}")

    # Persist the β-scan so the criticality-check plot can show it.
    with open(cache_scan, "w") as f:
        json.dump({
            "beta_c": float(beta_c),
            "chi_peak": float(chi_peak),
            "scan_betas": [float(b) for b in ref_sb],
            "scan_chis": [float(c) for c in ref_sc],
            "scan_chi_errs": [float(e) for e in ref_sce],
        }, f, indent=2)

    print(f"[ref] running production MC: n_traj={n_traj}")
    _, subdir = mc_engine.run_simulator(
        exe, L, L, 0, 0, 1.0, 1.0, 1.0, beta_c,
        n_traj=n_traj, data_dir=os.path.join(scratch, "prod"),
    )
    src_a2a = os.path.join(subdir, "two_point_all_to_all.dat")

    import shutil
    shutil.copy(src_a2a, cache_a2a)
    meta = {"L": L, "beta_c": float(beta_c), "n_traj": n_traj,
            "geometry": [L, L, 0, 0]}
    with open(cache_meta, "w") as f:
        json.dump(meta, f, indent=2)
    shutil.rmtree(scratch, ignore_errors=True)

    data = mc_engine.load_all_to_all(cache_a2a)
    return data, meta


# ---------------------------------------------------------------------------
#  Run one optimizer
# ---------------------------------------------------------------------------

def run_one(method, evaluator, *, max_evals, x0=(1.0, 1.0), snr_floor=1.0,
            method_kwargs=None):
    fn = opt_mod.ALL_METHODS[method]
    print(f"\n========== {method} ==========")
    t0 = time.time()
    kwargs = dict(x0=x0, max_evals=max_evals, snr_floor=snr_floor)
    if method_kwargs:
        kwargs.update(method_kwargs)
    summary = fn(evaluator, **kwargs)
    summary["wall_total_s"] = time.time() - t0
    print(f"[{method}] done: {summary.get('status')}  "
          f"n={summary.get('n_evals')}  "
          f"best_cost={summary.get('best_cost')}")
    return summary


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=None,
                    metavar="DIR",
                    help="root directory for all outputs of this run "
                         "(default: results/).  Each run should use its own "
                         "directory so results are never overwritten. "
                         "The shared reference is always read from/written to "
                         "results/_reference/ regardless of this setting.")
    ap.add_argument("--method", default="nelder_mead",
                    help="one of " + ", ".join(opt_mod.ALL_METHODS) + ", or 'all'")
    ap.add_argument("--exe", default=DEFAULT_EXE)
    ap.add_argument("--ref-L", type=int, default=8)
    ap.add_argument("--test-L", type=int, default=None,
                    help="defaults to ref-L (self-consistency test)")
    ap.add_argument("--ref-n-traj", type=int, default=60000)
    ap.add_argument("--n-traj-prod", type=int, default=30000)
    ap.add_argument("--n-traj-scan-coarse", type=int, default=4000)
    ap.add_argument("--n-traj-scan-fine", type=int, default=10000)
    ap.add_argument("--max-evals", type=int, default=30)
    ap.add_argument("--snr-floor", type=float, default=0.0,
                    help="early-stop when running-best SNR drops below this. "
                         "0 (default) disables the floor (use for benchmarks); "
                         "1.0 is the production setting.")
    ap.add_argument("--x0", type=float, nargs=2, default=(1.0, 1.0),
                    metavar=("R1", "R2"))
    ap.add_argument("--beta-seed", type=float, nargs=2, default=(0.20, 0.32),
                    metavar=("LO", "HI"),
                    help="initial β bracket for the very first evaluator call. "
                         "Widen this when --x0 is far from (1, 1) so the "
                         "shifting GC scan can find the peak.")
    # --- Nelder-Mead knobs ---
    ap.add_argument("--nm-sigma0", type=float, default=0.1,
                    help="NM initial simplex edge length. "
                         "Use ~0.4–0.5×‖x0−truth‖ for far-start runs "
                         "(default 0.1).")
    ap.add_argument("--nm-xatol", type=float, default=0.005,
                    help="NM position convergence tolerance (default 0.005).")
    ap.add_argument("--nm-shrink", type=float, default=0.75,
                    help="NM shrink coefficient σ applied to all non-best "
                         "vertices when contraction fails (default 0.75; "
                         "scipy uses 0.5).")
    ap.add_argument("--noise-stop-snr", type=float, default=0.0,
                    help="NM early stop: halt if the last 5 consecutive evals "
                         "all have SNR below this value (0 = disabled). "
                         "Use ~1.5 to stop when the curve collapse is "
                         "noise-dominated and further steps are meaningless.")
    ap.add_argument("--nm-indist-stop-snr", type=float, default=0.0,
                    help="NM early stop: halt when the running best has "
                         "SNR < this value, meaning the test correlator is "
                         "statistically indistinguishable from the reference. "
                         "Use 1.0 (default disabled = 0).")
    ap.add_argument("--nm-snr-target", type=float, default=0.0,
                    help="NM adaptive stats: boost n_traj_prod by 1.5x whenever "
                         "SNR < this value (default 0 = disabled).")
    ap.add_argument("--nm-snr-max-traj-factor", type=int, default=4,
                    help="NM adaptive stats: cap boost at this multiple of "
                         "starting n_traj_prod (default 4).")
    # --- BFGS-hand knobs ---
    ap.add_argument("--bfgs-fd-eps", type=float, default=0.05,
                    help="BFGS central-FD step size (default 0.05). "
                         "Increase if MC noise sigma_cost > fd_eps × |grad|.")
    ap.add_argument("--bfgs-sigma0", type=float, default=0.1,
                    help="BFGS initial step length scale in (r1,r2) space "
                         "(default 0.1). Bump up for far-start runs.")
    ap.add_argument("--bfgs-max-step", type=float, default=0.4,
                    help="BFGS maximum step length in (r1,r2) space (default 0.4).")
    ap.add_argument("--bfgs-noise-stop-snr", type=float, default=0.0,
                    help="BFGS noise-dominated stop threshold (same semantics "
                         "as --noise-stop-snr for NM, default 0 = off).")
    ap.add_argument("--bfgs-indist-stop-snr", type=float, default=0.0,
                    help="BFGS indistinguishable-from-reference stop (default 0 = off).")
    ap.add_argument("--bfgs-snr-target", type=float, default=3.0,
                    help="BFGS adaptive stats: boost n_traj_prod by 1.5x whenever "
                         "SNR < this value (default 3.0; 0 disables).")
    ap.add_argument("--bfgs-snr-max-traj-factor", type=int, default=4,
                    help="BFGS adaptive stats: cap boost at this multiple of "
                         "starting n_traj_prod (default 4).")
    # --- Visualizer knobs ---
    ap.add_argument("--mu", type=float, default=0.0,
                    help="Anisotropy-penalty weight (option D): cost += mu*sum(C_d^2). "
                         "Default 0 = plain L². Typical starting value: 1/mean(C_d).")
    ap.add_argument("--save-every", type=int, default=5,
                    help="write one optimizer PNG frame every N evaluations "
                         "(default 5; use 1 to recover old every-eval behaviour).")
    ap.add_argument("--no-vis", action="store_true",
                    help="disable PNG frame writing (faster)")
    args = ap.parse_args()

    if not os.path.exists(args.exe):
        sys.exit(f"simulator not found at {args.exe}.  Run `make` first.")

    RESULTS = args.results_dir or DEFAULT_RESULTS
    # Reference directory: shared per lattice size so different L runs don't
    # collide.  L=8 keeps the canonical name for backward compatibility.
    if args.ref_L == 8:
        REF_DIR = os.path.join(DEFAULT_RESULTS, "_reference")
    else:
        REF_DIR = os.path.join(DEFAULT_RESULTS, f"_reference_L{args.ref_L}")

    test_L = args.test_L or args.ref_L
    os.makedirs(RESULTS, exist_ok=True)

    # 1. Reference
    ref_data, ref_meta = build_reference(
        args.exe, L=args.ref_L, n_traj=args.ref_n_traj,
        out_dir=REF_DIR,
    )
    ref_geom = tuple(ref_meta["geometry"])
    test_geom = (test_L, test_L, 0, 0)

    # 2. Method list
    if args.method == "all":
        methods = list(opt_mod.ALL_METHODS)
    else:
        if args.method not in opt_mod.ALL_METHODS:
            sys.exit(f"unknown method: {args.method}")
        methods = [args.method]

    # 3. Run each method
    all_summaries = []
    for m in methods:
        m_dir = os.path.join(RESULTS, m)
        os.makedirs(m_dir, exist_ok=True)
        frames_dir = os.path.join(m_dir, "frames")

        opt_plot = None
        beta_plot_dir = None
        if not args.no_vis:
            from visualization import OptimizerPlotter
            opt_plot = OptimizerPlotter(
                frames_dir, m, ref_data,
                test_geom[0], test_geom[1], test_geom[2], test_geom[3],
                save_every=args.save_every,
            )
            beta_plot_dir = frames_dir

        ev = Evaluator(
            args.exe, ref_data, ref_geom, test_geom, m_dir,
            n_traj_prod=args.n_traj_prod,
            n_traj_scan_coarse=args.n_traj_scan_coarse,
            n_traj_scan_fine=args.n_traj_scan_fine,
            beta_seed=tuple(args.beta_seed),
            optimizer_plot=opt_plot,
            beta_plot_dir=beta_plot_dir,
            mu=args.mu,
        )

        nm_kwargs = {"sigma0": args.nm_sigma0,
                     "xatol": args.nm_xatol,
                     "shrink": args.nm_shrink,
                     "noise_stop_snr": args.noise_stop_snr,
                     "indist_stop_snr": args.nm_indist_stop_snr,
                     "snr_target": args.nm_snr_target,
                     "snr_max_traj_factor": args.nm_snr_max_traj_factor}
        bfgs_hand_kwargs = {"fd_eps": args.bfgs_fd_eps,
                            "sigma0": args.bfgs_sigma0,
                            "max_step": args.bfgs_max_step,
                            "noise_stop_snr": args.bfgs_noise_stop_snr,
                            "indist_stop_snr": args.bfgs_indist_stop_snr,
                            "snr_target": args.bfgs_snr_target,
                            "snr_max_traj_factor": args.bfgs_snr_max_traj_factor}
        _method_kwargs = {
            "nelder_mead": nm_kwargs,
            "bfgs_hand":   bfgs_hand_kwargs,
        }
        summary = run_one(m, ev, max_evals=args.max_evals,
                          x0=tuple(args.x0), snr_floor=args.snr_floor,
                          method_kwargs=_method_kwargs.get(m, {}))
        with open(os.path.join(m_dir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        all_summaries.append(summary)

    # 4. Top-level summary
    out = os.path.join(RESULTS, "summary.json")
    with open(out, "w") as f:
        json.dump({"reference": ref_meta,
                   "test_geom": list(test_geom),
                   "x0": list(args.x0),
                   "methods": all_summaries}, f, indent=2)
    print(f"\n[done] wrote {out}")
    print("       per-method dirs in", RESULTS)


if __name__ == "__main__":
    main()
