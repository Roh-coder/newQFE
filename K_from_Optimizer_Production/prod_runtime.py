"""
prod_runtime.py — production glue for K_from_Optimizer.

Two responsibilities:

1. ``install_reweight_fallback(...)``
       Monkey-patches ``mc_engine.find_beta_c_reweight`` so that any
       failure of the single-donor Ferrenberg–Swendsen path falls back
       to the legacy 3-pass Gram–Charlier susceptibility-peak scan.
       Logs each call (status, geometry, β_c, error) to a JSONL file so
       the user can audit the fallback rate after a run.

2. ``ProdGenerationPool``
       Drop-in replacement for ``parallel.GenerationPool`` whose worker
       initializer also installs the reweight fallback inside each
       worker process.  When ``n_workers <= 1`` it degenerates to the
       serial path (the patch already applied in the main process is
       inherited by the in-process Evaluator).

Both pieces are intentionally additive: the upstream Speedup_Upgrade
modules are imported unmodified.  Bug fixes that land there propagate
here automatically.
"""
from __future__ import annotations

import json
import os
import pickle
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# 0.  Ensure THIS directory is importable so workers (spawn-mode) and any
#     `python -m` invocations both find the local modules.
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import mc_engine  # noqa: E402


# ---------------------------------------------------------------------------
# 1.  Reweight → 3-pass fallback
# ---------------------------------------------------------------------------
#
# State lives at module-level so it is inherited by forked workers; in
# spawned (Windows / macOS-default) workers we reinstall via the
# ``initargs`` of ProdGenerationPool below.
# ---------------------------------------------------------------------------

_FB_LOCK = threading.Lock()

_FB_CFG: dict = {
    "log_path":     None,    # Optional[str]
    "scan_kwargs":  {},      # passed to mc_engine.find_beta_c on fallback
    "enabled":      False,
    "verbose":      True,
    # Two-pass multi-donor reweighting (highest priority when enabled).
    "multidonor_2pass":  False,
    "md2_kwargs":        {},  # passed to mc_engine.find_beta_c_multidonor_2pass
    # Single-pass multi-donor reweighting (second priority when enabled).
    "multidonor":   False,
    "md_kwargs":    {},      # passed to mc_engine.find_beta_c_multidonor
}

_FB_STATS: dict = {
    "multidonor_2pass_ok": 0,
    "multidonor_ok":  0,
    "reweight_ok":   0,
    "fallback_used": 0,
    "fallback_failed": 0,
}

# Hold a reference to the unpatched original so we can fall back / restore.
_ORIG_REWEIGHT = mc_engine.find_beta_c_reweight
_ORIG_LEGACY   = mc_engine.find_beta_c


def _audit(record: dict) -> None:
    path = _FB_CFG["log_path"]
    if not path:
        return
    record = {"ts": time.time(), **record}
    with _FB_LOCK:
        with open(path, "a") as f:
            f.write(json.dumps(record) + "\n")


def _patched_reweight(exe, Lx, Ly, Tx, Ty, k1, k2, k3,
                      beta_lo, beta_hi, *,
                      n_traj_parent=40000, n_grid=201,
                      n_eff_floor=0.10, max_recenters=3,
                      data_dir="/tmp/k_from_cs_rw",
                      progress_cb=None, jackknife=False):
    """Try 2-pass multi-donor → 1-pass multi-donor → single-donor reweight → 3-pass GC scan.

    The 2-pass multi-donor path is taken only when ``install_reweight_fallback``
    was called with ``multidonor_2pass=True``.
    The 1-pass multi-donor path is taken when ``multidonor=True`` (and 2-pass
    is not enabled or has failed).
    Otherwise behaves as the original single-donor reweight wrapper.
    """
    # ----- -1. Two-pass multi-donor reweighting (highest priority) -----
    if _FB_CFG.get("multidonor_2pass"):
        try:
            md2_kw = dict(_FB_CFG.get("md2_kwargs", {}))
            md2_data_dir = data_dir.rstrip("/").rstrip("\\") + "_md2"
            result = mc_engine.find_beta_c_multidonor_2pass(
                exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_lo, beta_hi,
                n_traj_parent=n_traj_parent, n_grid=n_grid,
                n_eff_floor=n_eff_floor,
                data_dir=md2_data_dir,
                progress_cb=progress_cb, jackknife=jackknife,
                **md2_kw,
            )
            with _FB_LOCK:
                _FB_STATS["multidonor_2pass_ok"] += 1
            _audit({
                "status": "multidonor_2pass",
                "geom":   [Lx, Ly, Tx, Ty],
                "k":      [k1, k2, k3],
                "beta":   [beta_lo, beta_hi],
                "beta_c": float(result[0]),
                "sigma":  float(result[1]),
            })
            return result
        except Exception as exc_md2:
            if _FB_CFG["verbose"]:
                print(f"[prod] multidonor_2pass FAILED "
                      f"({type(exc_md2).__name__}: {exc_md2}); "
                      f"falling back to 1-pass multidonor/reweight…")
            _audit({
                "status": "multidonor_2pass_failed",
                "geom":   [Lx, Ly, Tx, Ty],
                "k":      [k1, k2, k3],
                "beta":   [beta_lo, beta_hi],
                "error":  f"{type(exc_md2).__name__}: {exc_md2}",
            })
            # Fall through to 1-pass multi-donor or single-donor path below.

    # ----- 0. Single-pass multi-donor reweighting (secondary, if enabled) -----
    if _FB_CFG.get("multidonor"):
        try:
            md_kw = dict(_FB_CFG.get("md_kwargs", {}))
            # Use a sibling scratch dir so multidonor + reweight + fallback
            # don't collide on simulator output paths.
            md_data_dir = data_dir.rstrip("/").rstrip("\\") + "_md"
            result = mc_engine.find_beta_c_multidonor(
                exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_lo, beta_hi,
                n_traj_parent=n_traj_parent, n_grid=n_grid,
                n_eff_floor=n_eff_floor,
                data_dir=md_data_dir,
                progress_cb=progress_cb, jackknife=jackknife,
                **md_kw,
            )
            with _FB_LOCK:
                _FB_STATS["multidonor_ok"] += 1
            _audit({
                "status": "multidonor",
                "geom":   [Lx, Ly, Tx, Ty],
                "k":      [k1, k2, k3],
                "beta":   [beta_lo, beta_hi],
                "beta_c": float(result[0]),
                "sigma":  float(result[1]),
            })
            return result
        except Exception as exc_md:
            if _FB_CFG["verbose"]:
                print(f"[prod] multidonor FAILED "
                      f"({type(exc_md).__name__}: {exc_md}); "
                      f"falling back to single-donor reweight…")
            _audit({
                "status": "multidonor_failed",
                "geom":   [Lx, Ly, Tx, Ty],
                "k":      [k1, k2, k3],
                "beta":   [beta_lo, beta_hi],
                "error":  f"{type(exc_md).__name__}: {exc_md}",
            })
            # Fall through to single-donor path below.

    # ----- 1. Single-donor reweight (original path) -----
    try:
        result = _ORIG_REWEIGHT(
            exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_lo, beta_hi,
            n_traj_parent=n_traj_parent, n_grid=n_grid,
            n_eff_floor=n_eff_floor, max_recenters=max_recenters,
            data_dir=data_dir, progress_cb=progress_cb,
            jackknife=jackknife,
        )
        with _FB_LOCK:
            _FB_STATS["reweight_ok"] += 1
        _audit({
            "status": "reweight",
            "geom":   [Lx, Ly, Tx, Ty],
            "k":      [k1, k2, k3],
            "beta":   [beta_lo, beta_hi],
            "beta_c": float(result[0]),
            "sigma":  float(result[1]),
        })
        return result
    except Exception as exc:
        if _FB_CFG["verbose"]:
            print(f"[prod] reweight FAILED ({type(exc).__name__}: {exc}); "
                  f"falling back to 3-pass GC scan…")
        with _FB_LOCK:
            _FB_STATS["fallback_used"] += 1

        sk = dict(_FB_CFG["scan_kwargs"])
        # Use a fresh scratch sibling so reweight + fallback don't share
        # output files (otherwise the C++ binary's `opening file:` parse
        # in find_beta_c is ambiguous).
        fb_data_dir = data_dir.rstrip("/").rstrip("\\") + "_fb"
        try:
            result = _ORIG_LEGACY(
                exe, Lx, Ly, Tx, Ty, k1, k2, k3, beta_lo, beta_hi,
                data_dir=fb_data_dir,
                progress_cb=progress_cb, jackknife=jackknife,
                **sk,
            )
        except Exception as exc2:
            with _FB_LOCK:
                _FB_STATS["fallback_failed"] += 1
            _audit({
                "status": "fallback_failed",
                "geom":   [Lx, Ly, Tx, Ty],
                "k":      [k1, k2, k3],
                "beta":   [beta_lo, beta_hi],
                "error_reweight": f"{type(exc).__name__}: {exc}",
                "error_legacy":   f"{type(exc2).__name__}: {exc2}",
            })
            raise
        _audit({
            "status": "fallback_3pass",
            "geom":   [Lx, Ly, Tx, Ty],
            "k":      [k1, k2, k3],
            "beta":   [beta_lo, beta_hi],
            "beta_c": float(result[0]),
            "sigma":  float(result[1]),
            "error_reweight": f"{type(exc).__name__}: {exc}",
        })
        return result


def install_reweight_fallback(scan_kwargs: dict,
                              log_path: Optional[str] = None,
                              verbose: bool = True,
                              *,
                              multidonor: bool = False,
                              md_kwargs: Optional[dict] = None,
                              multidonor_2pass: bool = False,
                              md2_kwargs: Optional[dict] = None) -> None:
    """Activate the reweight→3-pass fallback in the current process.

    ``scan_kwargs`` must contain the legacy ``find_beta_c`` knobs —
    ``n_coarse``, ``n_refine``, ``n_refine2``, ``n_refine3``,
    ``n_traj_coarse``, ``n_traj_fine``, ``max_shifts`` — to use on the
    fallback path.  Identical defaults to the production CONFIG keys
    ``scan_n_*`` / ``n_traj_scan_*`` / ``scan_max_shifts``.

    When ``multidonor_2pass=True``, the patched function tries
    :func:`mc_engine.find_beta_c_multidonor_2pass` first (two-pass
    interval refinement); on failure it falls through to single-pass
    multidonor (if ``multidonor=True``), then single-donor reweight,
    then the legacy 3-pass GC scan.
    ``md2_kwargs`` is forwarded to ``find_beta_c_multidonor_2pass``
    (e.g. ``pass1_n_donors``, ``pass1_budget_frac``,
    ``donor_overlap_alpha``, ``pilot_n_traj``).

    When ``multidonor=True``, the patched function tries
    :func:`mc_engine.find_beta_c_multidonor` (tiles the bracket
    with multiple parent histograms); on failure it falls through to the
    single-donor reweight, then to the legacy 3-pass GC scan.
    ``md_kwargs`` is forwarded to ``find_beta_c_multidonor`` (e.g.
    ``var_E``, ``n_donors``, ``donor_overlap_alpha``, ``pilot_n_traj``).
    """
    _FB_CFG["scan_kwargs"]      = dict(scan_kwargs)
    _FB_CFG["log_path"]         = str(log_path) if log_path else None
    _FB_CFG["verbose"]          = bool(verbose)
    _FB_CFG["enabled"]          = True
    _FB_CFG["multidonor_2pass"] = bool(multidonor_2pass)
    _FB_CFG["md2_kwargs"]       = dict(md2_kwargs) if md2_kwargs else {}
    _FB_CFG["multidonor"]       = bool(multidonor)
    _FB_CFG["md_kwargs"]        = dict(md_kwargs) if md_kwargs else {}
    mc_engine.find_beta_c_reweight = _patched_reweight


def fallback_stats() -> dict:
    with _FB_LOCK:
        return dict(_FB_STATS)


# ---------------------------------------------------------------------------
# 2.  Parallel pool with fallback installed in each worker.
# ---------------------------------------------------------------------------

_W_EVAL = None    # type: ignore[assignment]
_W_RNG  = None    # type: ignore[assignment]


def _prod_worker_init(eval_kwargs: dict, seed_state: bytes,
                      fb_kwargs: dict) -> None:
    """Per-worker startup: install fallback, then build a headless Evaluator."""
    global _W_EVAL, _W_RNG

    # Be defensive when the worker was launched via "spawn" (Win/macOS):
    # the parent's sys.path is replayed by ProcessPoolExecutor, but make
    # sure the production dir is at the front.
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)

    install_reweight_fallback(
        scan_kwargs=fb_kwargs.get("scan_kwargs", {}),
        log_path=fb_kwargs.get("log_path"),
        verbose=fb_kwargs.get("verbose", False),
        multidonor_2pass=fb_kwargs.get("multidonor_2pass", False),
        md2_kwargs=fb_kwargs.get("md2_kwargs"),
        multidonor=fb_kwargs.get("multidonor", False),
        md_kwargs=fb_kwargs.get("md_kwargs"),
    )

    ss = pickle.loads(seed_state)
    _W_RNG = np.random.default_rng(ss)

    from evaluator   import Evaluator
    from betac_cache import BetacCache

    cache_kwargs = eval_kwargs.pop("_betac_cache", None)
    cache = BetacCache(**cache_kwargs) if cache_kwargs is not None else None

    _W_EVAL = Evaluator(
        **eval_kwargs,
        optimizer_plot=None,
        beta_plot_dir=None,
        dashboard=None,
        betac_cache=cache,
    )


def _prod_worker_eval(point_with_id) -> dict:
    eid, r1, r2 = point_with_id
    import shutil as _sh

    _W_EVAL._eval_id = eid - 1
    # Defer scratch cleanup so we can read the production a2a file for
    # visualization (residual curves) in the main process.
    old_keep = _W_EVAL.keep_mc_subdirs
    _W_EVAL.keep_mc_subdirs = True
    try:
        res = _W_EVAL(float(r1), float(r2))
    finally:
        _W_EVAL.keep_mc_subdirs = old_keep

    d = asdict(res)

    # Load the test correlator from the prod MC output before wiping scratch.
    label   = f"ev{eid:04d}"  # must match evaluator.py
    scratch = os.path.join(_W_EVAL._mc_root, label)
    try:
        prod_dir = os.path.join(scratch, "prod")
        for root, _, files in os.walk(prod_dir):
            if "two_point_all_to_all.dat" in files:
                d["_test_data"] = mc_engine.load_all_to_all(
                    os.path.join(root, "two_point_all_to_all.dat"))
                break
    except Exception:
        pass
    # Now clean up (unless the caller asked to keep subdirs).
    if not old_keep:
        _sh.rmtree(scratch, ignore_errors=True)

    return d


class ProdGenerationPool:
    """ProcessPool fan-out for CMA-ES generations with reweight fallback.

    Mirrors ``parallel.GenerationPool`` from the upstream package but:
      * installs ``install_reweight_fallback`` inside each worker, and
      * inherits the same fallback state in the main process for the
        ``n_workers <= 1`` serial path.
    """

    def __init__(self, n_workers: int, evaluator_kwargs: dict,
                 fb_kwargs: dict,
                 master_seed: Optional[int] = None):
        self.n_workers   = max(1, int(n_workers))
        self.master_ss   = np.random.SeedSequence(master_seed)
        self._pool: Optional[ProcessPoolExecutor] = None
        self._spawned_seeds: list = []
        self._eval_kwargs = dict(evaluator_kwargs)
        self._fb_kwargs   = dict(fb_kwargs)
        self._serial_eval = None

        if self.n_workers > 1:
            child_ss = self.master_ss.spawn(self.n_workers)
            self._spawned_seeds = [pickle.dumps(s) for s in child_ss]
            self._pool = ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=_prod_worker_init,
                initargs=(self._eval_kwargs, self._spawned_seeds[0],
                          self._fb_kwargs),
            )
        else:
            from evaluator   import Evaluator
            from betac_cache import BetacCache
            cache_kwargs = self._eval_kwargs.pop("_betac_cache", None)
            cache = (BetacCache(**cache_kwargs)
                     if cache_kwargs is not None else None)
            self._serial_eval = Evaluator(
                **self._eval_kwargs,
                optimizer_plot=None, beta_plot_dir=None, dashboard=None,
                betac_cache=cache,
            )

    def attach_callbacks(self, *, optimizer_plot=None,
                         beta_plot_dir=None, dashboard=None) -> None:
        if self._serial_eval is None:
            return
        self._serial_eval.optimizer_plot = optimizer_plot
        self._serial_eval.beta_plot_dir  = beta_plot_dir
        self._serial_eval.dashboard      = dashboard

    @property
    def serial_evaluator(self):
        return self._serial_eval

    def map_generation(self, points, eval_id_base: int):
        tasks = [(eval_id_base + k, float(r1), float(r2))
                 for k, (r1, r2) in enumerate(points)]
        if self.n_workers <= 1 or self._pool is None:
            results = []
            for eid, r1, r2 in tasks:
                self._serial_eval._eval_id = eid - 1
                res = self._serial_eval(r1, r2)
                results.append(asdict(res))
            return results

        futures = {self._pool.submit(_prod_worker_eval, t): t[0] for t in tasks}
        out = {}
        for fut in as_completed(futures):
            eid = futures[fut]
            try:
                out[eid] = fut.result()
            except Exception as e:
                raise RuntimeError(f"worker failed on eval_id={eid}") from e
        return [out[eid] for eid, _, _ in tasks]

    def shutdown(self) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=True, cancel_futures=True)
            self._pool = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.shutdown()


# ---------------------------------------------------------------------------
# 3.  Helper: snapshot the latest optimizer frame for presentation.
# ---------------------------------------------------------------------------

def snapshot_frames(latest_png: str, history_dir: str,
                    eval_id: Optional[int] = None) -> Optional[str]:
    """Copy the rolling ``optimizer_<method>.png`` to a numbered history PNG.

    Returns the new path on success, ``None`` if the source doesn't exist.
    Idempotent — if the latest mtime matches the previous snapshot, it's
    a no-op.  Safe to call after every evaluation.
    """
    if not os.path.exists(latest_png):
        return None
    os.makedirs(history_dir, exist_ok=True)
    mtime = os.path.getmtime(latest_png)

    state_path = os.path.join(history_dir, ".last_mtime")
    if os.path.exists(state_path):
        try:
            if abs(float(open(state_path).read().strip()) - mtime) < 1e-6:
                return None
        except Exception:
            pass

    if eval_id is None:
        # Auto-number: count existing snapshots.
        n = sum(1 for f in os.listdir(history_dir)
                if f.startswith("frame_") and f.endswith(".png"))
        eval_id = n + 1
    dst = os.path.join(history_dir, f"frame_{eval_id:04d}.png")
    import shutil
    shutil.copy2(latest_png, dst)
    with open(state_path, "w") as f:
        f.write(f"{mtime}\n")
    return dst
