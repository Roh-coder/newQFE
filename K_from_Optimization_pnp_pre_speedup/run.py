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
    # Reference lattice: equilateral couplings (k1=k2=k3=1) on an arbitrary
    # twisted parallelogram (Lx, Ly, Tx, Ty).  For an untwisted square set
    # Lx=Ly and Tx=Ty=0.
    "ref_Lx":         24,
    "ref_Ly":         24,
    "ref_Tx":         0,
    "ref_Ty":         0,

    # Test lattice geometry (the one we're optimising couplings for).
    "test_Lx":        24,
    "test_Ly":        24,
    "test_Tx":        0,
    "test_Ty":        0,

    # --- Reference statistics (one-time cost) ---
    "ref_n_traj":               1_000_000,   # production trajectories at ref β_c
    "ref_scan_n_traj_coarse":   10_000,    # per-point during ref β_c scan
    "ref_scan_n_traj_fine":    20_000,    # per-point during ref refinement passes
    "ref_beta_seed":            (0.20, 0.32),  # initial bracket for β_c scan
    # β_c scan-method knobs for the REFERENCE pass (3-pass Gram-Charlier).
    # Pass 0 sweeps `n_coarse` points across `[beta_lo, beta_hi]` (translating
    # the window up to `max_shifts` times if no interior peak is found).
    # Passes 1–3 then sample `n_refine[1/2/3] + 2` points around the running
    # estimate, with the bracket shrinking each pass to a multiple of the
    # GC-fit σ.  Larger numbers = more MC = a more carefully resolved peak.
    "ref_scan_n_coarse":        11,
    "ref_scan_n_refine":        5,
    "ref_scan_n_refine2":       5,
    "ref_scan_n_refine3":       5,
    "ref_scan_max_shifts":      4,

    # --- Test-lattice statistics (per optimizer evaluation) ---
    "n_traj_prod":              50_000,   # production trajectories at test β_c
    "n_traj_scan_coarse":       10_000,    # per-point coarse β_c scan
    "n_traj_scan_fine":        20_000,    # per-point refinement
    "beta_seed":                (0.20, 0.32),  # initial bracket (warm-started after eval 1)
    # Per-evaluation β_c scan-method knobs (mirror the ref_scan_* keys).
    # Per-eval scans typically use fewer refinement points than the
    # reference because the optimizer already needs many evaluations.
    "scan_n_coarse":            11,
    "scan_n_refine":            5,
    "scan_n_refine2":           5,
    "scan_n_refine3":           5,
    "scan_max_shifts":          4,
    # Leave-one-out jackknife on the GC scan points to estimate σ(β_c).
    # Adds ~N extra GC fits per scan (fast: O(10–200 ms) vs minutes of MC).
    # Applies to BOTH the reference build and per-eval scans.
    "scan_jackknife":           False,

    # --- Optimizer choice ---
    # "nelder_mead": hand-coded NM simplex (default; deterministic).
    # "cmaes":       Covariance Matrix Adaptation ES (requires `pip install cma`).
    #                More robust on noisy/flat cost surfaces; uses a
    #                population of `cma_popsize` evals per generation.
    "optimizer":        "nelder_mead",

    # --- Nelder-Mead optimizer knobs ---
    "x0":               (0.5, 1.0),  # starting point in (r1, r2)
    "max_evals":        30,
    "nm_sigma0":        0.10,        # initial simplex side length
    "nm_xatol":         0.005,       # convergence tolerance on simplex diameter
    "nm_fatol":         1e-6,        # convergence tolerance on cost spread
    "nm_shrink":        0.75,        # shrink coeff (0.75 is more noise-tolerant than scipy's 0.5)

    # --- CMA-ES optimizer knobs (used only if optimizer == "cmaes") ---
    "cma_sigma0":       0.10,        # initial step size
    "cma_popsize":      6,           # λ: evals per generation (4+3*log(N) ≈ 6 for N=2)
    "cma_max_gens":     0,           # hard cap on generations (0 = no cap, only max_evals limits)
    "cma_tolx":         0.005,       # convergence tol on solution spread
    "cma_tolfun":       1e-6,        # convergence tol on cost spread
    "cma_seed":         None,        # int for reproducibility, None = random

    # Adaptive statistics: when SNR < snr_target, multiply n_traj_prod by 1.5
    # (capped at snr_max_traj_factor × starting value).  Set snr_target=0 to
    # disable.  Recommended: snr_target=2.0–3.0, snr_max_traj_factor=4.
    "snr_target":               2.0,
    "snr_max_traj_factor":      4,

    # Stopping rules (set 0 to disable each).
    "snr_floor":                0.0,    # stop if running best SNR < this
    "indist_stop_snr":          1.0,    # stop if best is statistically indistinguishable

    # --- Output ---
    "results_dir":      "results",     # root directory for all outputs
    # Per-job output subdirectory under results_dir.  Leave as None to
    # auto-generate a tag based on the test geometry, wall-clock timestamp,
    # and process id — this lets multiple concurrent runs coexist without
    # clobbering each other's eval_log.jsonl / summary.json / frames.  Set
    # to a string to pin a specific name.
    "run_name":         None,
    "save_every":       1,             # write optimizer PNG every N evals (1 = every eval)
    "no_vis":           False,         # disable PNG frames entirely (faster)
    "dashboard":        True,          # live terminal dashboard (rich); disables PNG vis when on

    # --- Simulator binary ---
    "exe":              "bin/ising_tri_twisted_parallelogram",

    # ------------------------------------------------------------------
    # Speedup flags (Phases 1–3 of the implementation strategy in the
    # README).  Defaults are the LEGACY values so the workflow that
    # produced the recorded baselines is preserved bit-for-bit when all
    # three flags are left untouched.  Flip them on individually to
    # benchmark each speedup against the legacy path.
    # ------------------------------------------------------------------

    # Speedup 2 — persistent β_c interpolation cache.
    # When True, every β_c the optimizer ever computes is appended to
    #   results/_betac_cache_<test_geom>/samples.jsonl
    # and a Delaunay+LinearND interpolant is rebuilt lazily.  Subsequent
    # queries within tol_r of an existing point and inside the convex
    # hull skip the full 3-pass MC scan and use the interpolant.
    "betac_cache":            False,
    "betac_tol_r":            0.05,    # max hull-distance tol for a hit
    "betac_tol_beta":         2e-3,    # require LOO interp residual < this
    "betac_min_neighbours":   4,       # skip lookups until cache has this many

    # Speedup 3 — parallel CMA-ES population.
    # n_workers=1 preserves the legacy serial loop.  n_workers>1 fans
    # the λ candidates of each generation out across a ProcessPool.
    # Has no effect on Nelder-Mead (which is fundamentally serial).
    "n_workers":              1,
    "master_seed":            None,    # int for reproducibility, None=random

    # Speedup 1 — persistent C++ back-end via line-delimited JSON RPC.
    # "subprocess" = legacy fork/exec per simulator call (default).
    # "rpc"        = single long-lived simulator child driven over
    #                stdin/stdout (Phase 3.1 ships only the production
    #                verb; find_beta_c and cost still go through Python).
    "backend":                "subprocess",
}

# ============================================================================
# End of CONFIG
# ============================================================================


HERE = os.path.dirname(os.path.abspath(__file__))


def ensure_binary(exe: str) -> None:
    """Compile the simulator with a C++ compiler if missing or out-of-date.

    Tries ``$CXX`` first, then ``g++``, ``clang++``, and ``c++`` from PATH.
    On Windows, also appends ``.exe`` to the output if not present.
    """
    import shutil as _shutil
    import subprocess

    src = os.path.join(HERE, "src", "ising_tri_twisted_parallelogram.cc")
    inc_dir = os.path.join(HERE, "include")
    headers = []
    if os.path.isdir(inc_dir):
        headers = [os.path.join(inc_dir, f) for f in os.listdir(inc_dir)
                   if f.endswith(".h")]

    # On Windows, ensure .exe extension
    exe_abs = exe if os.path.isabs(exe) else os.path.join(HERE, exe)
    if os.name == "nt" and not exe_abs.lower().endswith(".exe"):
        exe_abs += ".exe"

    needs_build = not os.path.exists(exe_abs)
    if not needs_build:
        try:
            exe_mtime = os.path.getmtime(exe_abs)
            for dep in [src, *headers]:
                if os.path.exists(dep) and os.path.getmtime(dep) > exe_mtime:
                    needs_build = True
                    break
        except OSError:
            needs_build = True

    if not needs_build:
        return

    os.makedirs(os.path.dirname(exe_abs) or ".", exist_ok=True)

    # Find a usable C++ compiler.
    candidates = []
    if os.environ.get("CXX"):
        candidates.append(os.environ["CXX"])
    candidates += ["g++", "clang++", "c++"]
    cxx = next((c for c in candidates if _shutil.which(c)), None)
    if cxx is None:
        sys.exit(
            "[build] No C++ compiler found on PATH (tried: "
            + ", ".join(candidates) + ").\n"
            "  • Linux/macOS: install g++ (e.g. `sudo apt install build-essential` "
            "or `xcode-select --install`)\n"
            "  • Windows: install MSYS2/MinGW-w64 g++, or run inside WSL.\n"
            "Or set the environment variable CXX to your compiler."
        )

    cxxflags = [
        "-std=c++14", "-O3", "-Wall",
        "-Wno-sign-compare", "-Wno-unknown-warning-option",
        "-Wno-deprecated-declarations",
    ]
    cmd = [cxx, *cxxflags, f"-I{inc_dir}", src, "-o", exe_abs]

    print(f"[build] {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[build] {cxx} failed (exit {e.returncode}); "
                 "fix the build error and re-run.")
    if not os.path.exists(exe_abs):
        sys.exit(f"[build] compile succeeded but {exe_abs!r} still missing.")
    print(f"[build] built {exe_abs}")


def build_reference(cfg, ref_dir):
    """Locate β_c for the reference lattice, run a long MC there, cache it."""
    import mc_engine

    cache_a2a  = os.path.join(ref_dir, "two_point_all_to_all.dat")
    cache_meta = os.path.join(ref_dir, "reference_meta.json")
    cache_scan = os.path.join(ref_dir, "reference_betac_scan.json")

    Lx = int(cfg["ref_Lx"]); Ly = int(cfg["ref_Ly"])
    Tx = int(cfg["ref_Tx"]); Ty = int(cfg["ref_Ty"])

    if os.path.exists(cache_a2a) and os.path.exists(cache_meta):
        with open(cache_meta) as f:
            meta = json.load(f)
        cached_geom = tuple(meta.get("geometry", []))
        if cached_geom == (Lx, Ly, Tx, Ty):
            print(f"[ref] using cached reference at {cache_a2a} "
                  f"(β_c={meta['beta_c']:.6f})")
            data = mc_engine.load_all_to_all(cache_a2a)
            return data, meta
        else:
            print(f"[ref] cache geometry {cached_geom} != requested "
                  f"{(Lx, Ly, Tx, Ty)}; rebuilding")

    print(f"[ref] building reference Lx={Lx} Ly={Ly} Tx={Tx} Ty={Ty} "
          f"(equilateral couplings) …")
    scratch = os.path.join(ref_dir, "_scratch")
    beta_lo, beta_hi = cfg["ref_beta_seed"]
    beta_c, beta_c_sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
        cfg["exe"], Lx, Ly, Tx, Ty, 1.0, 1.0, 1.0,
        beta_lo, beta_hi,
        n_coarse=cfg["ref_scan_n_coarse"],
        n_refine=cfg["ref_scan_n_refine"],
        n_refine2=cfg["ref_scan_n_refine2"],
        n_refine3=cfg["ref_scan_n_refine3"],
        n_traj_coarse=cfg["ref_scan_n_traj_coarse"],
        n_traj_fine=cfg["ref_scan_n_traj_fine"],
        max_shifts=cfg["ref_scan_max_shifts"],
        jackknife=cfg.get("scan_jackknife", False),
        data_dir=os.path.join(scratch, "scan"),
    )
    if beta_c_sigma > 0:
        print(f"[ref]  β_c = {beta_c:.6f} ± {beta_c_sigma:.2e} (jackknife)  "
              f"χ_peak={chi_peak:.3g}")
    else:
        print(f"[ref]  β_c = {beta_c:.6f}  χ_peak={chi_peak:.3g}")

    with open(cache_scan, "w") as f:
        json.dump({
            "beta_c": float(beta_c),
            "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "scan_betas": [float(b) for b in sb],
            "scan_chis":  [float(c) for c in sc],
            "scan_chi_errs": [float(e) for e in sce],
        }, f, indent=2)

    print(f"[ref] running production MC: n_traj={cfg['ref_n_traj']}")
    _, subdir = mc_engine.run_simulator(
        cfg["exe"], Lx, Ly, Tx, Ty, 1.0, 1.0, 1.0, beta_c,
        n_traj=cfg["ref_n_traj"],
        data_dir=os.path.join(scratch, "prod"),
    )
    src_a2a = os.path.join(subdir, "two_point_all_to_all.dat")
    import shutil
    shutil.copy(src_a2a, cache_a2a)
    meta = {"beta_c": float(beta_c), "beta_c_sigma": float(beta_c_sigma),
            "n_traj": cfg["ref_n_traj"],
            "geometry": [Lx, Ly, Tx, Ty]}
    with open(cache_meta, "w") as f:
        json.dump(meta, f, indent=2)
    shutil.rmtree(scratch, ignore_errors=True)

    data = mc_engine.load_all_to_all(cache_a2a)
    return data, meta


def main():
    cfg = CONFIG
    os.chdir(HERE)

    ensure_binary(cfg["exe"])

    # On Windows the compiled binary has a .exe suffix — propagate that
    # to the rest of the pipeline.
    if os.name == "nt" and not cfg["exe"].lower().endswith(".exe"):
        cfg["exe"] = cfg["exe"] + ".exe"

    results_dir = cfg["results_dir"]
    os.makedirs(results_dir, exist_ok=True)

    # Reference cache: one folder per geometry so changing ref params doesn't collide.
    ref_tag = (f"Lx{cfg['ref_Lx']}_Ly{cfg['ref_Ly']}"
               f"_Tx{cfg['ref_Tx']}_Ty{cfg['ref_Ty']}")
    ref_dir = os.path.join(results_dir, f"_reference_{ref_tag}")
    os.makedirs(ref_dir, exist_ok=True)

    # 1. Reference correlator
    ref_data, ref_meta = build_reference(cfg, ref_dir)
    ref_geom  = tuple(ref_meta["geometry"])
    test_geom = (cfg["test_Lx"], cfg["test_Ly"], cfg["test_Tx"], cfg["test_Ty"])

    # 2. Set up evaluator + visualization
    # Per-job output directory (avoid clobbering when multiple runs share
    # the same results/ tree).  Keep the tag SHORT — on Windows the total
    # path down to an MC output file must fit in MAX_PATH (260 chars), and
    # the simulator writes files like
    #   <run>/mc_scratch/eval####_r1_x.xxxx_r2_x.xxxx/scan/
    #         <lx>x<ly>_t<tx>x<ty>_k<...>_<seed>.dat
    # which already adds ~125 chars.  Default = "j<HHMMSS>_p<pid>" (~15 chars).
    run_name = cfg.get("run_name")
    if not run_name:
        stamp = time.strftime("%H%M%S")
        run_name = f"j{stamp}_p{os.getpid()}"
    nm_dir = os.path.join(results_dir, run_name)
    os.makedirs(nm_dir, exist_ok=True)
    print(f"[run] output dir: {nm_dir}")
    # Drop a tiny meta file so you can identify which test/ref geometry
    # this short-named run directory belonged to.
    with open(os.path.join(nm_dir, "run_meta.json"), "w") as _f:
        json.dump({
            "run_name": run_name,
            "ref_geometry":  [cfg["ref_Lx"],  cfg["ref_Ly"],
                              cfg["ref_Tx"],  cfg["ref_Ty"]],
            "test_geometry": [cfg["test_Lx"], cfg["test_Ly"],
                              cfg["test_Tx"], cfg["test_Ty"]],
            "started":   time.strftime("%Y-%m-%d %H:%M:%S"),
            "pid":       os.getpid(),
        }, _f, indent=2)
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
            frames_dir, cfg.get("optimizer", "nelder_mead"), ref_data,
            test_geom[0], test_geom[1], test_geom[2], test_geom[3],
            save_every=cfg["save_every"],
            ref_Lx=ref_geom[0], ref_Ly=ref_geom[1],
            ref_Tx=ref_geom[2], ref_Ty=ref_geom[3],
        )
        # The β_c scan is now embedded as the 5th panel of the optimizer
        # PNG, so we don't need a separate betac_scan.png file.

    dash = None
    if use_dashboard:
        from dashboard import Dashboard
        dash = Dashboard(ref_meta, test_geom, cfg["max_evals"],
                         method_name=cfg.get("optimizer", "nelder_mead"))

    from evaluator import Evaluator

    # Speedup 2 — persistent β_c interpolation cache.
    bc_cache = None
    if cfg.get("betac_cache", False):
        from betac_cache import BetacCache
        bc_cache = BetacCache(
            test_geom, root=results_dir,
            tol_r=cfg.get("betac_tol_r", 0.05),
            tol_beta=cfg.get("betac_tol_beta", 2e-3),
            min_neighbours=cfg.get("betac_min_neighbours", 4),
        )
        print(f"[cache] β_c cache enabled at {bc_cache.dir} "
              f"(currently {len(bc_cache)} samples)")

    ev = Evaluator(
        cfg["exe"], ref_data, ref_geom, test_geom, nm_dir,
        n_traj_prod=cfg["n_traj_prod"],
        n_traj_scan_coarse=cfg["n_traj_scan_coarse"],
        n_traj_scan_fine=cfg["n_traj_scan_fine"],
        scan_n_coarse=cfg["scan_n_coarse"],
        scan_n_refine=cfg["scan_n_refine"],
        scan_n_refine2=cfg["scan_n_refine2"],
        scan_n_refine3=cfg["scan_n_refine3"],
        scan_max_shifts=cfg["scan_max_shifts"],
        scan_jackknife=cfg.get("scan_jackknife", False),
        beta_seed=tuple(cfg["beta_seed"]),
        optimizer_plot=opt_plot,
        beta_plot_dir=beta_plot_dir,
        dashboard=dash,
        betac_cache=bc_cache,
    )

    # 3. Run optimizer (Nelder-Mead or CMA-ES)
    optimizer_name = cfg.get("optimizer", "nelder_mead")
    print(f"\n========== {optimizer_name} ==========")
    t0 = time.time()

    if optimizer_name == "nelder_mead":
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
    elif optimizer_name == "cmaes":
        from optimizer import run_cmaes

        def _run():
            return run_cmaes(
                ev,
                x0=tuple(cfg["x0"]),
                max_evals=cfg["max_evals"],
                max_gens=cfg.get("cma_max_gens", 0),
                sigma0=cfg["cma_sigma0"],
                popsize=cfg["cma_popsize"],
                tolx=cfg["cma_tolx"],
                tolfun=cfg["cma_tolfun"],
                snr_floor=cfg["snr_floor"],
                indist_stop_snr=cfg["indist_stop_snr"],
                snr_target=cfg["snr_target"],
                snr_max_traj_factor=cfg["snr_max_traj_factor"],
                seed=cfg.get("cma_seed"),
            )
    else:
        raise ValueError(
            f"Unknown optimizer {optimizer_name!r}; "
            "expected 'nelder_mead' or 'cmaes'."
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

    tag = "nm" if optimizer_name == "nelder_mead" else "cma"
    print(f"\n[{tag}] done: {summary.get('status')}  "
          f"n={summary.get('n_evals')}  "
          f"best=(r1={summary.get('best_r1'):.4f}, r2={summary.get('best_r2'):.4f})  "
          f"cost={summary.get('best_cost'):.4e}  SNR={summary.get('best_snr'):.2f}")

    summary_path = os.path.join(nm_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[{tag}] summary written to {summary_path}")
    if opt_plot is not None:
        print(f"[{tag}] frames in {frames_dir}/")


if __name__ == "__main__":
    main()
