#!/usr/bin/env python3
"""
run.py — K_from_Optimizer Production: plug-and-play CMA-ES + reweight + fallback.

Usage
-----
Edit the ``CONFIG`` block below, then::

    python run.py                   # run with CONFIG defaults
    python run.py --max-evals 30    # override any CONFIG key from the CLI
    python run.py --live            # pop a live PNG viewer window
    python run.py --save-frames     # also keep one numbered PNG per eval

All CONFIG keys are exposed as ``--<key-with-dashes>``; tuple-valued keys
take 2+ space-separated numbers (e.g. ``--x0 0.9 1.1``).

Pipeline
--------
1. Build (or load from cache) the reference correlator (long MC at the
   reference β_c).
2. Drive CMA-ES on the test lattice; each candidate (r1, r2) goes through
   a Ferrenberg–Swendsen β_c reweight.  If reweighting fails (typically
   because the dense-grid peak is outside the validity window), the
   evaluator transparently falls back to the legacy 3-pass Gram–Charlier
   scan.
3. Optionally fan out the λ candidates of each generation across N
   worker processes.
4. Emit a 4-panel PNG per evaluation (or every ``save_every``).  With
   ``--live`` a Tk window reloads the rolling frame in real time; with
   ``--save-frames`` every rendered PNG is also archived to
   ``<run>/frames_history/``.

Standalone — no dependencies on sister directories.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import threading
import time
from typing import Optional

# Wire local imports first.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


# ============================================================================
# CONFIG — edit these and run `python run.py`
# ============================================================================

CONFIG = {
    # --- Lattice geometries -----------------------------------------------
    "ref_Lx":  12, "ref_Ly":  12, "ref_Tx":  0, "ref_Ty":  0,
    "test_Lx": 12, "test_Ly": 12, "test_Tx": 0, "test_Ty": 0,

    # --- Reference statistics (one-time per geometry; cached on disk) -----
    "ref_n_traj":             100_000,
    "ref_scan_n_traj_coarse":   4_000,
    "ref_scan_n_traj_fine":    10_000,
    "ref_beta_seed":           (0.20, 0.32),
    "ref_scan_n_coarse":       11,
    "ref_scan_n_refine":        5,
    "ref_scan_n_refine2":       5,
    "ref_scan_n_refine3":       5,
    "ref_scan_max_shifts":      4,

    # --- Per-evaluation MC ------------------------------------------------
    "n_traj_prod":             20_000,
    "n_traj_scan_coarse":       4_000,   # used only on fallback
    "n_traj_scan_fine":        10_000,   # used only on fallback
    "beta_seed":               (0.20, 0.32),
    "scan_n_coarse":           11,
    "scan_n_refine":            5,
    "scan_n_refine2":           5,
    "scan_n_refine3":           5,
    "scan_max_shifts":          4,
    "scan_jackknife":          True,

    # --- Optimizer (cmaes recommended for production) ---------------------
    "optimizer":     "cmaes",
    "x0":            (1.0, 1.0),
    "max_evals":     40,
    "cma_sigma0":    0.10,
    "cma_popsize":   6,
    "cma_max_gens":  0,
    "cma_tolx":      0.005,
    "cma_tolfun":    1e-6,
    "cma_seed":      None,

    # --- Adaptive statistics ----------------------------------------------
    "snr_target":           2.0,
    "snr_max_traj_factor":  4,
    "snr_floor":            0.0,
    "indist_stop_snr":      1.0,

    # --- Speedup S4: Ferrenberg–Swendsen β reweighting --------------------
    # In production this directory ALWAYS uses reweight=True; on failure
    # it transparently falls back to the 3-pass GC scan.
    "reweight":              True,
    "reweight_n_traj":      40_000,
    "reweight_n_grid":         201,
    "reweight_n_eff_floor":   0.10,
    "reweight_max_recenters":    3,

    # --- Speedup S3: parallel CMA-ES population ---------------------------
    # n_workers=0 (default) → auto: uses min(os.cpu_count(), cma_popsize).
    # n_workers=1 → serial (no subprocess overhead).
    # n_workers=N → exactly N workers.
    "n_workers":   0,
    "master_seed": None,

    # --- Output / visualization ------------------------------------------
    "results_dir":    "results",
    "run_name":       None,         # auto = "j<HHMMSS>_p<pid>"
    "save_every":     1,            # write rolling PNG every N evals
    "no_vis":         False,        # True → no PNGs at all (fastest)
    "dashboard":      True,         # rich terminal dashboard
    "save_frames":    False,        # archive every rolling PNG to frames_history/
    "live":           False,        # pop a Tk window watching the rolling PNG
    "live_poll_ms":   500,          # viewer refresh interval (ms)

    # --- Simulator binary -------------------------------------------------
    "exe": ("bin/ising_tri_twisted_parallelogram.exe"
            if sys.platform == "win32"
            else "bin/ising_tri_twisted_parallelogram"),
}

# ============================================================================


# ---------------------------------------------------------------------------
# Reference build + binary compile (lifted from upstream run.py, trimmed).
# ---------------------------------------------------------------------------

def ensure_binary(exe: str) -> None:
    """Compile the C++ simulator if missing, using $CXX or g++/clang++/c++."""
    # Resolve paths relative to this file so the function works regardless
    # of the caller's cwd.
    exe_abs = exe if os.path.isabs(exe) else os.path.join(_HERE, exe)
    if os.path.exists(exe_abs):
        return
    out_dir = os.path.dirname(exe_abs)
    os.makedirs(out_dir, exist_ok=True)
    src_abs = os.path.join(_HERE, "src", "ising_tri_twisted_parallelogram.cc")
    inc_abs = os.path.join(_HERE, "include")
    if not os.path.exists(src_abs):
        raise FileNotFoundError(f"Cannot build {exe_abs}: {src_abs} not found")
    candidates = [os.environ.get("CXX")] + ["g++", "clang++", "c++"]
    import subprocess
    for cxx in filter(None, candidates):
        if not shutil.which(cxx):
            continue
        cmd = [cxx, "-O3", "-std=c++17", f"-I{inc_abs}", src_abs, "-o", exe_abs]
        print(f"[build] {' '.join(cmd)}", flush=True)
        if subprocess.run(cmd).returncode == 0:
            print(f"[build] OK → {exe_abs}", flush=True)
            return
    raise RuntimeError(f"No working C++ compiler found to build {exe_abs}")


def build_reference(cfg, ref_dir):
    """Locate β_c on the reference lattice and run a long MC there.  Cached."""
    import mc_engine

    cache_a2a  = os.path.join(ref_dir, "two_point_all_to_all.dat")
    cache_meta = os.path.join(ref_dir, "reference_meta.json")
    cache_scan = os.path.join(ref_dir, "reference_betac_scan.json")

    Lx = int(cfg["ref_Lx"]); Ly = int(cfg["ref_Ly"])
    Tx = int(cfg["ref_Tx"]); Ty = int(cfg["ref_Ty"])

    if os.path.exists(cache_a2a) and os.path.exists(cache_meta):
        with open(cache_meta) as f:
            meta = json.load(f)
        if tuple(meta.get("geometry", [])) == (Lx, Ly, Tx, Ty):
            print(f"[ref] using cached reference at {cache_a2a} "
                  f"(β_c={meta['beta_c']:.6f})")
            return mc_engine.load_all_to_all(cache_a2a), meta

    print(f"[ref] building reference Lx={Lx} Ly={Ly} Tx={Tx} Ty={Ty} …")
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
    print(f"[ref]  β_c = {beta_c:.6f}"
          f"{(' ± ' + format(beta_c_sigma, '.2e')) if beta_c_sigma > 0 else ''}"
          f"  χ_peak={chi_peak:.3g}")

    with open(cache_scan, "w") as f:
        json.dump({
            "beta_c": float(beta_c), "beta_c_sigma": float(beta_c_sigma),
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
    shutil.copy(os.path.join(subdir, "two_point_all_to_all.dat"), cache_a2a)
    meta = {"beta_c": float(beta_c), "beta_c_sigma": float(beta_c_sigma),
            "n_traj": cfg["ref_n_traj"],
            "geometry": [Lx, Ly, Tx, Ty]}
    with open(cache_meta, "w") as f:
        json.dump(meta, f, indent=2)
    shutil.rmtree(scratch, ignore_errors=True)

    return mc_engine.load_all_to_all(cache_a2a), meta


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------

def _coerce(orig, val):
    """Coerce a string CLI value back into the type stored in CONFIG."""
    if orig is None or isinstance(orig, str):
        return val
    if isinstance(orig, bool):
        return val.lower() in ("1", "true", "yes", "on") if isinstance(val, str) else bool(val)
    if isinstance(orig, int):
        return int(val)
    if isinstance(orig, float):
        return float(val)
    return val


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__.split("\n\n", 1)[0],
    )
    p.add_argument("--config", type=str, default=None,
                   help="Optional JSON file whose keys are merged on top of CONFIG.")
    p.add_argument("--results-dir", type=str, default=None)
    p.add_argument("--run-name", type=str, default=None)
    p.add_argument("--max-evals", type=int, default=None)
    p.add_argument("--n-workers", type=int, default=None,
                   help="Workers for parallel CMA-ES. 0=auto (default): "
                        "min(cpu_count, popsize). 1=serial.")
    p.add_argument("--cma-popsize", type=int, default=None)
    p.add_argument("--cma-max-gens", type=int, default=None,
                   help="Stop after this many CMA-ES generations (0=unlimited).")
    p.add_argument("--cma-sigma0", type=float, default=None)
    p.add_argument("--n-traj-prod", type=int, default=None)
    p.add_argument("--reweight-n-traj", type=int, default=None)
    p.add_argument("--save-every", type=int, default=None)
    p.add_argument("--ref-Lx", type=int, default=None)
    p.add_argument("--ref-Ly", type=int, default=None)
    p.add_argument("--test-Lx", type=int, default=None)
    p.add_argument("--test-Ly", type=int, default=None)
    p.add_argument("--x0", type=float, nargs=2, default=None,
                   metavar=("R1", "R2"))
    p.add_argument("--no-reweight", action="store_true",
                   help="Force the legacy 3-pass GC scan on every eval.")
    p.add_argument("--indist-stop-snr", type=float, default=None,
                   help="Stop when best SNR < this (default 1.0). Set 0 to disable.")
    p.add_argument("--no-vis", action="store_true",
                   help="Disable PNG rendering entirely (fastest).")
    p.add_argument("--no-dashboard", action="store_true",
                   help="Disable the rich terminal dashboard.")
    p.add_argument("--live", action="store_true",
                   help="Pop a Tk window that auto-reloads the rolling PNG.")
    p.add_argument("--save-frames", action="store_true",
                   help="Archive every rolling PNG into <run>/frames_history/.")
    p.add_argument("--live-poll-ms", type=int, default=None)
    p.add_argument("--print-config", action="store_true",
                   help="Print the resolved CONFIG and exit.")
    return p


def apply_cli(cfg: dict, args: argparse.Namespace) -> dict:
    """Overlay CLI overrides on top of CONFIG (after optional --config JSON)."""
    cfg = dict(cfg)
    if args.config:
        with open(args.config) as f:
            patch = json.load(f)
        for k in ("ref_beta_seed", "beta_seed", "x0"):
            if k in patch and isinstance(patch[k], list):
                patch[k] = tuple(patch[k])
        cfg.update(patch)

    mapping = {
        "results_dir":    args.results_dir,
        "run_name":       args.run_name,
        "max_evals":      args.max_evals,
        "n_workers":      args.n_workers,
        "cma_popsize":    args.cma_popsize,
        "cma_max_gens":   args.cma_max_gens,
        "cma_sigma0":     args.cma_sigma0,
        "n_traj_prod":    args.n_traj_prod,
        "reweight_n_traj": args.reweight_n_traj,
        "save_every":     args.save_every,
        "ref_Lx":         args.ref_Lx,
        "ref_Ly":         args.ref_Ly,
        "test_Lx":        args.test_Lx,
        "test_Ly":        args.test_Ly,
        "live_poll_ms":   args.live_poll_ms,
        "indist_stop_snr": args.indist_stop_snr,
    }
    for k, v in mapping.items():
        if v is not None:
            cfg[k] = _coerce(cfg.get(k), v)

    if args.x0 is not None:
        cfg["x0"] = tuple(args.x0)
    if args.no_reweight:
        cfg["reweight"] = False
    if args.no_vis:
        cfg["no_vis"] = True
    if args.no_dashboard:
        cfg["dashboard"] = False
    if args.live:
        cfg["live"] = True
    if args.save_frames:
        cfg["save_frames"] = True
    return cfg


# ---------------------------------------------------------------------------
# Optimizer driver
# ---------------------------------------------------------------------------

def _make_evaluator_kwargs(cfg, ref_data, ref_geom, test_geom, nm_dir):
    return dict(
        exe=cfg["exe"], ref_data=ref_data,
        ref_geom=ref_geom, test_geom=test_geom, output_dir=nm_dir,
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
        reweight=cfg.get("reweight", True),
        reweight_n_traj=cfg["reweight_n_traj"],
        reweight_n_grid=cfg["reweight_n_grid"],
        reweight_n_eff_floor=cfg["reweight_n_eff_floor"],
        reweight_max_recenters=cfg["reweight_max_recenters"],
    )


def _fb_kwargs(cfg, log_path):
    return dict(
        scan_kwargs=dict(
            n_coarse=cfg["scan_n_coarse"],
            n_refine=cfg["scan_n_refine"],
            n_refine2=cfg["scan_n_refine2"],
            n_refine3=cfg["scan_n_refine3"],
            n_traj_coarse=cfg["n_traj_scan_coarse"],
            n_traj_fine=cfg["n_traj_scan_fine"],
            max_shifts=cfg["scan_max_shifts"],
        ),
        log_path=log_path,
        verbose=True,
    )


def run_optimizer(cfg: dict) -> dict:
    """Run one optimizer end-to-end; returns the summary dict."""
    import prod_runtime
    from evaluator    import Evaluator
    from optimizer    import run_cmaes, run_nelder_mead

    os.chdir(_HERE)
    ensure_binary(cfg["exe"])

    results_dir = cfg["results_dir"]
    os.makedirs(results_dir, exist_ok=True)

    ref_tag = (f"Lx{cfg['ref_Lx']}_Ly{cfg['ref_Ly']}"
               f"_Tx{cfg['ref_Tx']}_Ty{cfg['ref_Ty']}")
    ref_dir = os.path.join(results_dir, f"_reference_{ref_tag}")
    os.makedirs(ref_dir, exist_ok=True)

    ref_data, ref_meta = build_reference(cfg, ref_dir)
    ref_geom  = tuple(ref_meta["geometry"])
    test_geom = (cfg["test_Lx"], cfg["test_Ly"], cfg["test_Tx"], cfg["test_Ty"])

    run_name = cfg.get("run_name") or f"j{time.strftime('%H%M%S')}_p{os.getpid()}"
    nm_dir = os.path.join(results_dir, run_name)
    os.makedirs(nm_dir, exist_ok=True)
    print(f"[run] output dir: {nm_dir}")
    with open(os.path.join(nm_dir, "run_meta.json"), "w") as f:
        json.dump({
            "run_name": run_name,
            "ref_geometry":  list(ref_geom),
            "test_geometry": list(test_geom),
            "started":   time.strftime("%Y-%m-%d %H:%M:%S"),
            "pid":       os.getpid(),
            "config":    {k: (list(v) if isinstance(v, tuple) else v)
                          for k, v in cfg.items()},
        }, f, indent=2)
    frames_dir = os.path.join(nm_dir, "frames")

    # Install reweight → 3-pass fallback in the main process.
    fb_log = os.path.join(nm_dir, "fallback_log.jsonl")
    prod_runtime.install_reweight_fallback(
        scan_kwargs=_fb_kwargs(cfg, fb_log)["scan_kwargs"],
        log_path=fb_log, verbose=True,
    )

    # Build the optimizer-side plotter / dashboard.
    use_png_vis  = not cfg.get("no_vis", False)
    use_dash     = bool(cfg.get("dashboard", True))
    opt_plot = None
    if use_png_vis:
        from visualization import OptimizerPlotter
        opt_plot = OptimizerPlotter(
            frames_dir, cfg.get("optimizer", "cmaes"), ref_data,
            test_geom[0], test_geom[1], test_geom[2], test_geom[3],
            save_every=cfg["save_every"],
            ref_Lx=ref_geom[0], ref_Ly=ref_geom[1],
            ref_Tx=ref_geom[2], ref_Ty=ref_geom[3],
        )
    dash = None
    if use_dash:
        from dashboard import Dashboard
        dash = Dashboard(ref_meta, test_geom, cfg["max_evals"],
                         method_name=cfg.get("optimizer", "cmaes"))

    # Build Evaluator (serial) or ProdGenerationPool (parallel).
    _requested = int(cfg.get("n_workers", 0))
    if _requested == 0:  # auto
        _cpus = os.cpu_count() or 1
        _pop  = int(cfg.get("cma_popsize", 6))
        n_workers = min(_cpus, _pop)
        print(f"[run] n_workers=auto → {n_workers} (cpus={_cpus}, popsize={_pop})")
    else:
        n_workers = max(1, _requested)
    optimizer_name = cfg.get("optimizer", "cmaes")

    if optimizer_name == "cmaes" and n_workers > 1:
        # Parallel CMA-ES.  Workers are headless; opt_plot/dash run in main.
        ev_kwargs = _make_evaluator_kwargs(cfg, ref_data, ref_geom,
                                           test_geom, nm_dir)
        pool = prod_runtime.ProdGenerationPool(
            n_workers=n_workers,
            evaluator_kwargs=ev_kwargs,
            fb_kwargs=_fb_kwargs(cfg, fb_log),
            master_seed=cfg.get("master_seed"),
        )
        # We still need a "host" Evaluator-like object for opt_plot / dash
        # forwarding from worker results.  Use a small shim.
        ev = _PoolEvaluatorShim(opt_plot, dash, n_traj_prod=cfg["n_traj_prod"])
        # The CMA-ES path looks for `current_simplex`/`current_gaussian` on
        # the evaluator — the shim already exposes both.
        run_kwargs = dict(
            x0=tuple(cfg["x0"]), max_evals=cfg["max_evals"],
            max_gens=cfg.get("cma_max_gens", 0),
            sigma0=cfg["cma_sigma0"], popsize=cfg["cma_popsize"],
            tolx=cfg["cma_tolx"], tolfun=cfg["cma_tolfun"],
            snr_floor=cfg["snr_floor"],
            indist_stop_snr=cfg["indist_stop_snr"],
            snr_target=cfg["snr_target"],
            snr_max_traj_factor=cfg["snr_max_traj_factor"],
            seed=cfg.get("cma_seed"),
            pool=_PoolForwarder(pool, opt_plot, dash, ev_shim=ev),
        )
        runner = lambda: run_cmaes(ev, **run_kwargs)
    else:
        ev = Evaluator(
            **_make_evaluator_kwargs(cfg, ref_data, ref_geom, test_geom, nm_dir),
            optimizer_plot=opt_plot,
            beta_plot_dir=None,
            dashboard=dash,
            betac_cache=None,
        )
        if optimizer_name == "cmaes":
            runner = lambda: run_cmaes(
                ev, x0=tuple(cfg["x0"]), max_evals=cfg["max_evals"],
                max_gens=cfg.get("cma_max_gens", 0),
                sigma0=cfg["cma_sigma0"], popsize=cfg["cma_popsize"],
                tolx=cfg["cma_tolx"], tolfun=cfg["cma_tolfun"],
                snr_floor=cfg["snr_floor"],
                indist_stop_snr=cfg["indist_stop_snr"],
                snr_target=cfg["snr_target"],
                snr_max_traj_factor=cfg["snr_max_traj_factor"],
                seed=cfg.get("cma_seed"),
            )
        elif optimizer_name == "nelder_mead":
            runner = lambda: run_nelder_mead(
                ev, x0=tuple(cfg["x0"]), max_evals=cfg["max_evals"],
                sigma0=cfg.get("nm_sigma0", 0.10),
                xatol=cfg.get("nm_xatol", 0.005),
                fatol=cfg.get("nm_fatol", 1e-6),
                shrink=cfg.get("nm_shrink", 0.75),
                snr_floor=cfg["snr_floor"],
                indist_stop_snr=cfg["indist_stop_snr"],
                snr_target=cfg["snr_target"],
                snr_max_traj_factor=cfg["snr_max_traj_factor"],
            )
        else:
            raise ValueError(f"Unknown optimizer {optimizer_name!r}")

    # Optional frame-archival thread.
    history_dir = os.path.join(nm_dir, "frames_history")
    archiver_stop = threading.Event()
    archiver_thread: Optional[threading.Thread] = None
    if cfg.get("save_frames", False) and use_png_vis:
        latest = os.path.join(frames_dir,
                              f"optimizer_{cfg.get('optimizer','cmaes')}.png")

        def _archiver():
            while not archiver_stop.is_set():
                prod_runtime.snapshot_frames(latest, history_dir)
                archiver_stop.wait(0.5)
            prod_runtime.snapshot_frames(latest, history_dir)  # final snap

        archiver_thread = threading.Thread(target=_archiver, daemon=True)
        archiver_thread.start()
        print(f"[run] archiving frames → {history_dir}/")

    print(f"\n========== {optimizer_name} (n_workers={n_workers}) ==========")
    t0 = time.time()
    if dash is not None:
        with dash:
            summary = runner()
    else:
        summary = runner()
    summary["wall_total_s"] = time.time() - t0

    if use_png_vis:
        from visualization import flush_render_queue
        flush_render_queue()
    if archiver_thread is not None:
        archiver_stop.set()
        archiver_thread.join(timeout=2.0)

    # Append fallback stats and persist summary.
    summary["fallback_stats"] = prod_runtime.fallback_stats()
    summary["config"] = {k: (list(v) if isinstance(v, tuple) else v)
                         for k, v in cfg.items()}
    summary_path = os.path.join(nm_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[run] done: {summary.get('status')}  "
          f"n={summary.get('n_evals')}  "
          f"best=(r1={summary.get('best_r1', float('nan')):.4f}, "
          f"r2={summary.get('best_r2', float('nan')):.4f})  "
          f"cost={summary.get('best_cost', float('nan')):.4e}")
    print(f"[run] fallback stats: {summary['fallback_stats']}")
    print(f"[run] summary → {summary_path}")
    return summary


# ---------------------------------------------------------------------------
# Glue for parallel CMA-ES: shim that exposes the attributes run_cmaes
# touches on the evaluator (plot/dash forwarding lives on the pool side).
# ---------------------------------------------------------------------------

class _PoolEvaluatorShim:
    """Stand-in evaluator exposing the attrs run_cmaes mutates."""
    def __init__(self, opt_plot, dash, *, n_traj_prod: int):
        self.optimizer_plot   = opt_plot
        self.dashboard        = dash
        self.n_traj_prod      = int(n_traj_prod)
        self.current_simplex  = None
        self.current_gaussian = None
        self.current_gp_surface = None


class _DictAttr:
    """Wrap a dict so dashboard.update_eval can use attribute access."""
    __slots__ = ("__dict__",)
    def __init__(self, d): self.__dict__.update(d)


class _PoolForwarder:
    """Wrap GenerationPool so each batch result also fires plot/dash callbacks."""
    def __init__(self, pool, opt_plot, dash, ev_shim=None):
        self._pool = pool
        self._opt  = opt_plot
        self._dash = dash
        self._ev   = ev_shim   # _PoolEvaluatorShim → current_gaussian

    def map_generation(self, points, eval_id_base):
        results = self._pool.map_generation(points, eval_id_base)
        # Snapshot the CMA-ES Gaussian from the shim (set by run_cmaes before
        # calling map_generation so it's the same for all points in this gen).
        gaussian = (self._ev.current_gaussian
                    if self._ev is not None else None)
        # Forward each EvalResult-shaped dict to the plotter / dashboard
        # in eval_id order so the visualizer's monotone buffer is happy.
        # Suppress mid-batch renders: only render on the last result so the
        # whole generation costs exactly one matplotlib figure build.
        n = len(results)
        if self._opt is not None and n > 1:
            _orig_save_every = self._opt.save_every
            self._opt.save_every = n + 1   # suppress renders for all but last
        else:
            _orig_save_every = None
        for i, d in enumerate(results):
            is_last = (i == n - 1)
            if _orig_save_every is not None and is_last:
                self._opt.save_every = _orig_save_every  # restore before final update
            test_data = d.pop("_test_data", None)
            if self._opt is not None:
                # Populate the susceptibility scan panel from EvalResult data.
                if d.get("scan_betas"):
                    self._opt.update_scan_from_result(
                        d["scan_betas"], d["scan_chis"],
                        d.get("scan_chi_errs", []), d["beta_c"])
                self._opt.update(
                    r1=d["r1"], r2=d["r2"], cost=d["cost"],
                    sigma_cost=d.get("sigma_cost", 0.0),
                    beta_c=d["beta_c"], eval_id=d["eval_id"],
                    test_data=test_data,
                    gaussian=gaussian,
                )
            if self._dash is not None:
                try:
                    self._dash.update_eval(_DictAttr(d))
                except Exception:
                    pass
        return results


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    cfg  = apply_cli(CONFIG, args)

    if args.print_config:
        print(json.dumps({k: (list(v) if isinstance(v, tuple) else v)
                          for k, v in cfg.items()}, indent=2))
        return 0

    if cfg.get("live") and not cfg.get("no_vis"):
        # Run the optimizer in a background thread; Tk in main thread.
        from live_viewer import LiveViewer
        method = cfg.get("optimizer", "cmaes")
        run_name = cfg.get("run_name") or f"j{time.strftime('%H%M%S')}_p{os.getpid()}"
        cfg["run_name"] = run_name   # pin so the viewer can find the PNG
        latest_png = os.path.join(cfg["results_dir"], run_name,
                                  "frames", f"optimizer_{method}.png")
        viewer = LiveViewer(latest_png, poll_ms=cfg.get("live_poll_ms", 500),
                            title=f"K_from_Optimizer — {run_name}")
        viewer.start_optimizer_thread(lambda: run_optimizer(cfg))
        viewer.mainloop()
        return 0

    run_optimizer(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
