#!/usr/bin/env python3
"""
run.py — Front-end driver for the K_from_Optimization standalone pipeline.

EDIT THE `CONFIG` BLOCK BELOW to set up your run.

Pipeline
--------
1. Build (or load from cache) the reference correlator G_ref by:
     a. running a 3-pass Gram-Charlier susceptibility-peak scan on the
        ref lattice with k1=k2=k3=1 to find β_c
     b. running a long production MC at β_c on the ref lattice
2. Run Nelder-Mead in the (r1, r2) plane on the test lattice.  Each
   evaluation:
     a. finds β_c on the test lattice for couplings (r1, r2, 1)
     b. runs a production MC at that β_c
     c. computes the L² cost = Σ_d ∫(G_test_d - G_ref_d)² dt
3. Write summary.json and the visualization frames.
"""
from __future__ import annotations

import json
import os
import sys
import time

# ============================================================================
# CONFIG — edit these and run `python run.py`
# ============================================================================

CONFIG = {
    # --- Lattice geometries ---
    # Reference lattice: equilateral (k1=k2=k3=1), Lref × Lref, no twist.
    "ref_L":          24,

    # Test lattice geometry (the one we're optimising couplings for).
    # For an untwisted square the Lx=Ly and Tx=Ty=0.
    "test_Lx":        24,
    "test_Ly":        24,
    "test_Tx":        0,
    "test_Ty":        0,

    # --- Reference statistics (one-time cost) ---
    "ref_n_traj":               60_000,   # production trajectories at ref β_c
    "ref_scan_n_traj_coarse":   4_000,    # per-point during ref β_c scan
    "ref_scan_n_traj_fine":    10_000,    # per-point during ref refinement passes
    "ref_beta_seed":            (0.20, 0.32),  # initial bracket for β_c scan

    # --- Test-lattice statistics (per optimizer evaluation) ---
    "n_traj_prod":              30_000,   # production trajectories at test β_c
    "n_traj_scan_coarse":       4_000,    # per-point coarse β_c scan
    "n_traj_scan_fine":        10_000,    # per-point refinement
    "beta_seed":                (0.20, 0.32),  # initial bracket (warm-started after eval 1)

    # --- Nelder-Mead optimizer knobs ---
    "x0":               (1.0, 1.0),  # starting point in (r1, r2)
    "max_evals":        30,
    "nm_sigma0":        0.10,        # initial simplex side length
    "nm_xatol":         0.005,       # convergence tolerance on simplex diameter
    "nm_fatol":         1e-6,        # convergence tolerance on cost spread
    "nm_shrink":        0.75,        # shrink coeff (0.75 is more noise-tolerant than scipy's 0.5)

    # Adaptive statistics: when SNR < snr_target, multiply n_traj_prod by 1.5
    # (capped at snr_max_traj_factor × starting value).  Set snr_target=0 to
    # disable.  Recommended: snr_target=2.0–3.0, snr_max_traj_factor=4.
    "snr_target":               2.0,
    "snr_max_traj_factor":      4,

    # Stopping rules (set 0 to disable each).
    "snr_floor":                0.0,    # stop if running best SNR < this
    "indist_stop_snr":          1.0,    # stop if best is statistically indistinguishable

    # --- Output ---
    "results_dir":      "results",     # everything goes here
    "save_every":       1,             # write optimizer PNG every N evals (1 = every eval)
    "no_vis":           False,         # disable PNG frames entirely (faster)
    "dashboard":        True,          # live terminal dashboard (rich); disables PNG vis when on

    # --- Simulator binary ---
    "exe":              "bin/ising_tri_twisted_parallelogram",
}

# ============================================================================
# End of CONFIG
# ============================================================================


HERE = os.path.dirname(os.path.abspath(__file__))


def build_reference(cfg, ref_dir):
    """Locate β_c for an L×L equilateral lattice, run a long MC there, cache it."""
    import mc_engine

    cache_a2a  = os.path.join(ref_dir, "two_point_all_to_all.dat")
    cache_meta = os.path.join(ref_dir, "reference_meta.json")
    cache_scan = os.path.join(ref_dir, "reference_betac_scan.json")

    if os.path.exists(cache_a2a) and os.path.exists(cache_meta):
        with open(cache_meta) as f:
            meta = json.load(f)
        print(f"[ref] using cached reference at {cache_a2a} "
              f"(β_c={meta['beta_c']:.6f})")
        data = mc_engine.load_all_to_all(cache_a2a)
        return data, meta

    L = cfg["ref_L"]
    print(f"[ref] building {L}×{L} equilateral reference …")
    scratch = os.path.join(ref_dir, "_scratch")
    beta_lo, beta_hi = cfg["ref_beta_seed"]
    beta_c, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
        cfg["exe"], L, L, 0, 0, 1.0, 1.0, 1.0,
        beta_lo, beta_hi,
        n_traj_coarse=cfg["ref_scan_n_traj_coarse"],
        n_traj_fine=cfg["ref_scan_n_traj_fine"],
        data_dir=os.path.join(scratch, "scan"),
    )
    print(f"[ref]  β_c = {beta_c:.6f}  χ_peak={chi_peak:.3g}")

    with open(cache_scan, "w") as f:
        json.dump({
            "beta_c": float(beta_c),
            "chi_peak": float(chi_peak),
            "scan_betas": [float(b) for b in sb],
            "scan_chis":  [float(c) for c in sc],
            "scan_chi_errs": [float(e) for e in sce],
        }, f, indent=2)

    print(f"[ref] running production MC: n_traj={cfg['ref_n_traj']}")
    _, subdir = mc_engine.run_simulator(
        cfg["exe"], L, L, 0, 0, 1.0, 1.0, 1.0, beta_c,
        n_traj=cfg["ref_n_traj"],
        data_dir=os.path.join(scratch, "prod"),
    )
    src_a2a = os.path.join(subdir, "two_point_all_to_all.dat")
    import shutil
    shutil.copy(src_a2a, cache_a2a)
    meta = {"L": L, "beta_c": float(beta_c), "n_traj": cfg["ref_n_traj"],
            "geometry": [L, L, 0, 0]}
    with open(cache_meta, "w") as f:
        json.dump(meta, f, indent=2)
    shutil.rmtree(scratch, ignore_errors=True)

    data = mc_engine.load_all_to_all(cache_a2a)
    return data, meta


def main():
    cfg = CONFIG
    os.chdir(HERE)

    if not os.path.exists(cfg["exe"]):
        sys.exit(f"simulator not found at {cfg['exe']!r} — run `make` first.")

    results_dir = cfg["results_dir"]
    os.makedirs(results_dir, exist_ok=True)

    # Reference cache: one folder per L so changing ref_L doesn't collide.
    ref_dir = os.path.join(results_dir, f"_reference_L{cfg['ref_L']}")
    os.makedirs(ref_dir, exist_ok=True)

    # 1. Reference correlator
    ref_data, ref_meta = build_reference(cfg, ref_dir)
    ref_geom  = tuple(ref_meta["geometry"])
    test_geom = (cfg["test_Lx"], cfg["test_Ly"], cfg["test_Tx"], cfg["test_Ty"])

    # 2. Set up evaluator + visualization
    nm_dir = os.path.join(results_dir, "nelder_mead")
    os.makedirs(nm_dir, exist_ok=True)
    frames_dir = os.path.join(nm_dir, "frames")

    # Dashboard and PNG visualizer can run together — dashboard handles
    # the live β_c-scan progress (so we leave beta_plot_dir off to avoid
    # double-rendering), while OptimizerPlotter still saves the per-eval
    # PNG snapshots to disk.
    use_dashboard = bool(cfg.get("dashboard", False))
    use_png_vis = not cfg["no_vis"]

    opt_plot = None
    beta_plot_dir = None
    if use_png_vis:
        from visualization import OptimizerPlotter
        opt_plot = OptimizerPlotter(
            frames_dir, "nelder_mead", ref_data,
            test_geom[0], test_geom[1], test_geom[2], test_geom[3],
            save_every=cfg["save_every"],
        )
        # Only enable the per-eval β_c PNG plotter when the dashboard is OFF
        # (otherwise the dashboard already shows the live scan progress).
        if not use_dashboard:
            beta_plot_dir = frames_dir

    dash = None
    if use_dashboard:
        from dashboard import Dashboard
        dash = Dashboard(ref_meta, test_geom, cfg["max_evals"],
                         method_name="nelder_mead")

    from evaluator import Evaluator
    ev = Evaluator(
        cfg["exe"], ref_data, ref_geom, test_geom, nm_dir,
        n_traj_prod=cfg["n_traj_prod"],
        n_traj_scan_coarse=cfg["n_traj_scan_coarse"],
        n_traj_scan_fine=cfg["n_traj_scan_fine"],
        beta_seed=tuple(cfg["beta_seed"]),
        optimizer_plot=opt_plot,
        beta_plot_dir=beta_plot_dir,
        dashboard=dash,
    )

    # 3. Run Nelder-Mead
    print(f"\n========== nelder_mead ==========")
    t0 = time.time()
    from optimizer import run_nelder_mead

    def _run():
        return run_nelder_mead(
            ev,
            x0=tuple(cfg["x0"]),
            max_evals=cfg["max_evals"],
            sigma0=cfg["nm_sigma0"],
            xatol=cfg["nm_xatol"],
            fatol=cfg["nm_fatol"],
            shrink=cfg["nm_shrink"],
            snr_floor=cfg["snr_floor"],
            indist_stop_snr=cfg["indist_stop_snr"],
            snr_target=cfg["snr_target"],
            snr_max_traj_factor=cfg["snr_max_traj_factor"],
        )

    if dash is not None:
        with dash:
            summary = _run()
    else:
        summary = _run()
    summary["wall_total_s"] = time.time() - t0
    summary["config"] = {k: (list(v) if isinstance(v, tuple) else v)
                         for k, v in cfg.items()}

    if use_png_vis:
        from visualization import flush_render_queue
        flush_render_queue()

    print(f"\n[nm] done: {summary.get('status')}  "
          f"n={summary.get('n_evals')}  "
          f"best=(r1={summary.get('best_r1'):.4f}, r2={summary.get('best_r2'):.4f})  "
          f"cost={summary.get('best_cost'):.4e}  SNR={summary.get('best_snr'):.2f}")

    summary_path = os.path.join(nm_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[nm] summary written to {summary_path}")
    if opt_plot is not None:
        print(f"[nm] frames in {frames_dir}/")


if __name__ == "__main__":
    main()
