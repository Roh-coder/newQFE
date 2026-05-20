#!/usr/bin/env python3
"""
run_boundary_opt.py

Standalone boundary-continuum optimizer that can be zipped, moved, and run
offline. Edit CONFIG below and press Run in Spyder, or call it from the
terminal with an optional JSON config override:

    python run_boundary_opt.py --config configs/smoke_continuum.json

This runner is self-contained with respect to project code:
- local C++ source in src/ is used for the Monte Carlo executable,
- local Python helpers in lib/ are used for beta_c finding and geometry,
- all outputs go to this folder's results/ tree.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import pickle
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "lib"))

import mc_engine  # noqa: E402
from cost import _SQRT3_2, _tile_interp, boundary_paths  # noqa: E402


# Optional JSON override file. In Spyder you can set this to a config under
# configs/ instead of editing the hard-coded defaults below.
CONFIG_FILE = None


# ---------------------------------------------------------------------------
# IN-SCRIPT QUICK SWITCHES (Spyder-friendly)
# ---------------------------------------------------------------------------
# Set these for a one-line reference strategy switch without editing nested
# JSON fields in CONFIG.
#
# Allowed values for REFERENCE_MODE_OVERRIDE:
#   - None        -> keep mode from CONFIG / JSON override file
#   - "continuum" -> use reference_family.sizes for continuum extrapolation
#   - "large"     -> use one reference lattice at reference_family.large_geom
#
# LARGE-RUN RECIPE (edit these three lines only):
#   REFERENCE_MODE_OVERRIDE = "large"
#   REFERENCE_LARGE_GEOM_OVERRIDE = [80, 80, 0, 0]
#   RUN_TAG = "prod_large_ref_80"
#
# CONTINUUM RECIPE:
#   REFERENCE_MODE_OVERRIDE = "continuum"
#   REFERENCE_LARGE_GEOM_OVERRIDE = None
REFERENCE_MODE_OVERRIDE = None
REFERENCE_LARGE_GEOM_OVERRIDE = None


# ---------------------------------------------------------------------------
# HARD-CODED DEFAULTS (edit these in Spyder)
# ---------------------------------------------------------------------------
# Keep these as plain variables to make the workflow easy to configure without
# digging through nested JSON structures.

# Run / output
RUN_TAG = "smoke_default"
RESULTS_DIR = None
RESUME = True
PROGRESS_LOG_NAME = "progress_log.jsonl"

# Execution
EXE = None
AUTO_BUILD = True
FORCE_REBUILD = False
BUILD_COMMAND = ["make"]
BUILD_TIMEOUT_S = 600
N_WORKERS = 2

# (r1, r2) grid
R_MIN = 0.9
R_MAX = 1.1
R_STEP = 0.2
R1_VALUES = None
R2_VALUES = None
K3 = 1.0

# Test family
TEST_SIZES = [4, 6, 8]
TEST_LX_MULT = 1.0
TEST_LY_MULT = 1.0
TEST_TX_FRAC = 0.0
TEST_TY_FRAC = 0.0
TEST_GEOM_MAP = None
# Optional explicit test lattice table with full parameters per size.
# Format: [[Lx, Ly, Tx, Ty], ...]
# Internal fit-size L is inferred as max(|Lx|, |Ly|).
# Example:
# TEST_LATTICES = [
#     [8, 8, 0, 0],
#     [16, 16, 0, 0],
#     [24, 24, 0, 0],
# ]
TEST_LATTICES = None

# Reference family
REFERENCE_MODE = "continuum"  # continuum | large
REFERENCE_SIZES = [4, 6, 8]
REFERENCE_LX_MULT = 1.0
REFERENCE_LY_MULT = 1.0
REFERENCE_TX_FRAC = 0.25
REFERENCE_TY_FRAC = 0.25
REFERENCE_GEOM_MAP = None
# Optional explicit reference lattice table with full parameters per size.
# Format: [[Lx, Ly, Tx, Ty], ...]
# Internal fit-size L is inferred as max(|Lx|, |Ly|).
# If this table has exactly one row, the script auto-switches to large mode
# and uses that geometry as reference_family.large_geom.
# Example:
# REFERENCE_LATTICES = [
#     [80, 80, 0, 0],
# ]
REFERENCE_LATTICES = None
REFERENCE_LARGE_GEOM = [12, 12, 3, 3]
REFERENCE_FIXED_R1 = 1.0
REFERENCE_FIXED_R2 = 1.0

# Production MC
MC_N_TRAJ = 200
MC_N_SKIP = 20
MC_N_THERM = 500

# beta_c finder
BETA_LO = 0.05
BETA_HI = 0.40
BETA_N_COARSE = 11
BETA_N_REFINE = 5
BETA_N_REFINE2 = 5
BETA_N_REFINE3 = 5
BETA_N_TRAJ_COARSE = 40
BETA_N_TRAJ_FINE = 80
BETA_MAX_SHIFTS = 4
BETA_JACKKNIFE = False

# Analysis
ANALYSIS_K_VALUES = [1, 2, 3, 4, 5, 6, 7]
ANALYSIS_K_DENOMINATOR = 8
ANALYSIS_WEIGHTED = False
ANALYSIS_FIT_MODE = "quadratic"  # linear | quadratic | auto

# Monitor / plots
MONITOR_SHOW_LIVE = True
MONITOR_SHOW_FINAL = True
MONITOR_REFRESH_EVERY_JOBS = 1
MONITOR_REFRESH_EVERY_SECONDS = 0.5
MONITOR_SAVE_SNAPSHOT = True
MONITOR_DASHBOARD_SIZE = [14, 10]
MONITOR_FINAL_SIZE = [13.0, 5.8]
MONITOR_DPI = 150


def _build_default_config() -> dict[str, Any]:
    return {
        "run": {
            "tag": RUN_TAG,
        },
        "paths": {
            "results_dir": RESULTS_DIR,
            "resume": RESUME,
            "progress_log_name": PROGRESS_LOG_NAME,
        },
        "execution": {
            "exe": EXE,
            "auto_build": AUTO_BUILD,
            "force_rebuild": FORCE_REBUILD,
            "build_command": list(BUILD_COMMAND),
            "build_timeout_s": BUILD_TIMEOUT_S,
            "n_workers": N_WORKERS,
        },
        "grid": {
            "r_min": R_MIN,
            "r_max": R_MAX,
            "r_step": R_STEP,
            "r1_values": R1_VALUES,
            "r2_values": R2_VALUES,
            "k3": K3,
        },
        "test_family": {
            "sizes": list(TEST_SIZES),
            "geom_defaults": {
                "Lx_mult": TEST_LX_MULT,
                "Ly_mult": TEST_LY_MULT,
                "Tx_frac": TEST_TX_FRAC,
                "Ty_frac": TEST_TY_FRAC,
            },
            "geom_map": copy.deepcopy(TEST_GEOM_MAP),
            "lattices": copy.deepcopy(TEST_LATTICES),
        },
        "reference_family": {
            "mode": REFERENCE_MODE,
            "sizes": list(REFERENCE_SIZES),
            "geom_defaults": {
                "Lx_mult": REFERENCE_LX_MULT,
                "Ly_mult": REFERENCE_LY_MULT,
                "Tx_frac": REFERENCE_TX_FRAC,
                "Ty_frac": REFERENCE_TY_FRAC,
            },
            "geom_map": copy.deepcopy(REFERENCE_GEOM_MAP),
            "lattices": copy.deepcopy(REFERENCE_LATTICES),
            "large_geom": list(REFERENCE_LARGE_GEOM),
            "fixed_r1": REFERENCE_FIXED_R1,
            "fixed_r2": REFERENCE_FIXED_R2,
        },
        "mc": {
            "n_traj": MC_N_TRAJ,
            "n_skip": MC_N_SKIP,
            "n_therm": MC_N_THERM,
        },
        "beta_c_finder": {
            "beta_lo": BETA_LO,
            "beta_hi": BETA_HI,
            "n_coarse": BETA_N_COARSE,
            "n_refine": BETA_N_REFINE,
            "n_refine2": BETA_N_REFINE2,
            "n_refine3": BETA_N_REFINE3,
            "n_traj_coarse": BETA_N_TRAJ_COARSE,
            "n_traj_fine": BETA_N_TRAJ_FINE,
            "max_shifts": BETA_MAX_SHIFTS,
            "jackknife": BETA_JACKKNIFE,
        },
        "analysis": {
            "k_values": list(ANALYSIS_K_VALUES),
            "k_denominator": ANALYSIS_K_DENOMINATOR,
            "weighted": ANALYSIS_WEIGHTED,
            "fit_mode": ANALYSIS_FIT_MODE,
        },
        "monitor": {
            "show_live_monitor": MONITOR_SHOW_LIVE,
            "show_final_plots": MONITOR_SHOW_FINAL,
            "refresh_every_jobs": MONITOR_REFRESH_EVERY_JOBS,
            "refresh_every_seconds": MONITOR_REFRESH_EVERY_SECONDS,
            "save_monitor_snapshot": MONITOR_SAVE_SNAPSHOT,
            "dashboard_size": list(MONITOR_DASHBOARD_SIZE),
            "final_size": list(MONITOR_FINAL_SIZE),
            "dpi": MONITOR_DPI,
        },
    }


CONFIG = _build_default_config()


_T0_GLOBAL = time.time()


def _json_default(obj: Any):
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _deep_update(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def _load_json_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _elapsed(t0: float | None = None) -> str:
    seconds = time.time() - (t0 if t0 is not None else _T0_GLOBAL)
    return str(timedelta(seconds=int(seconds)))


def log(message: str) -> None:
    print(f"[{_ts()}  +{_elapsed()}]  {message}", flush=True)


def banner(title: str) -> None:
    bar = "=" * 72
    print(f"\n{bar}\n[{_ts()}  +{_elapsed()}]  {title}\n{bar}", flush=True)


class EventLog:
    def __init__(self, path: str):
        self.path = path
        os.makedirs(os.path.dirname(path), exist_ok=True)

    def write(self, event: str, **payload: Any) -> None:
        record = {
            "timestamp": _ts(),
            "elapsed_s": round(time.time() - _T0_GLOBAL, 3),
            "event": event,
            **payload,
        }
        with open(self.path, "a", encoding="utf-8") as handle:
            json.dump(record, handle, default=_json_default)
            handle.write("\n")


class LiveMonitor:
    def __init__(self, cfg: dict[str, Any], out_root: str,
                 r1_values: np.ndarray, r2_values: np.ndarray,
                 group_labels: list[str]):
        self.cfg = cfg
        self.out_root = out_root
        self.r1_values = r1_values
        self.r2_values = r2_values
        self.group_labels = group_labels
        mon = cfg["monitor"]
        self.enabled = bool(mon.get("show_live_monitor", True))
        self.refresh_every_jobs = max(int(mon.get("refresh_every_jobs", 1)), 1)
        self.refresh_every_seconds = float(mon.get("refresh_every_seconds", 0.5))
        self.save_snapshot_enabled = bool(mon.get("save_monitor_snapshot", True))
        self.dpi = int(mon.get("dpi", 150))
        self.fig = None
        self.axes = None
        self.last_refresh = 0.0
        if self.enabled:
            plt.ion()
            figsize = mon.get("dashboard_size", [14, 10])
            self.fig, self.axes = plt.subplots(2, 2, figsize=figsize)
            self.fig.suptitle(f"Boundary optimization monitor  (tag={cfg['run']['tag']})")

    def _heatmap(self, ax, Z: np.ndarray, title: str, cmap: str) -> None:
        ax.clear()
        extent = [self.r1_values.min(), self.r1_values.max(),
                  self.r2_values.min(), self.r2_values.max()]
        if np.any(np.isfinite(Z)):
            ax.imshow(Z, origin="lower", aspect="auto", extent=extent, cmap=cmap)
            jj, ii = np.unravel_index(np.nanargmin(Z), Z.shape)
            ax.plot(self.r1_values[ii], self.r2_values[jj], marker="*", markersize=12,
                    markerfacecolor="white", markeredgecolor="k")
        else:
            ax.imshow(np.zeros_like(Z), origin="lower", aspect="auto", extent=extent,
                      cmap=cmap, vmin=0.0, vmax=1.0)
        ax.set_xlabel("r1")
        ax.set_ylabel("r2")
        ax.set_title(title)

    def refresh(self, state: dict[str, Any], force: bool = False) -> None:
        if not self.enabled:
            return
        now = time.time()
        if not force:
            if state["event_index"] % self.refresh_every_jobs != 0:
                if now - self.last_refresh < self.refresh_every_seconds:
                    return
        ax_summary, ax_groups, ax_score, ax_z = self.axes.ravel()

        ax_summary.clear()
        ax_summary.axis("off")
        summary = [
            f"stage: {state['stage']}",
            f"jobs: {state['jobs_done']}/{state['jobs_total']}  ok={state['jobs_ok']}  cached={state['jobs_skip']}  err={state['jobs_err']}",
            f"analysis cells: {state['analysis_done']}/{state['analysis_total']}",
            f"elapsed: {_elapsed(state['t_start'])}",
        ]
        if state.get("eta"):
            summary.append(f"eta: {state['eta']}")
        if state.get("last_message"):
            summary.append(f"last: {state['last_message']}")
        ax_summary.text(0.02, 0.98, "\n".join(summary), va="top", ha="left",
                        family="monospace", fontsize=10)

        ax_groups.clear()
        fracs = []
        for label in self.group_labels:
            done, total = state["group_counts"].get(label, (0, 1))
            frac = float(done) / float(total) if total else 0.0
            fracs.append(frac)
        xpos = np.arange(len(self.group_labels))
        ax_groups.bar(xpos, fracs, color="#2f6c8f")
        ax_groups.set_xticks(xpos)
        ax_groups.set_xticklabels(self.group_labels, rotation=45, ha="right")
        ax_groups.set_ylim(0.0, 1.0)
        ax_groups.set_ylabel("completion")
        ax_groups.set_title("Per-group completion")
        ax_groups.grid(alpha=0.25, axis="y")

        self._heatmap(ax_score, state["score_grid"], "Partial score heatmap", "viridis")
        self._heatmap(ax_z, state["zscore_grid"], "Partial RMS z-score heatmap", "magma")

        self.fig.tight_layout()
        self.fig.canvas.draw_idle()
        plt.pause(0.001)
        self.last_refresh = now

    def save_snapshot(self) -> None:
        if self.fig is None or not self.save_snapshot_enabled:
            return
        path = os.path.join(self.out_root, "monitor_snapshot.png")
        self.fig.savefig(path, dpi=self.dpi)


def _default_exe_path() -> str:
    base = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    if os.name == "nt":
        return base + ".exe"
    return base


def _newest_source_mtime() -> float:
    newest = os.path.getmtime(os.path.join(_HERE, "Makefile"))
    for root, _, files in os.walk(os.path.join(_HERE, "src")):
        for name in files:
            newest = max(newest, os.path.getmtime(os.path.join(root, name)))
    for root, _, files in os.walk(os.path.join(_HERE, "include")):
        for name in files:
            newest = max(newest, os.path.getmtime(os.path.join(root, name)))
    return newest


def _ensure_simulator(cfg: dict[str, Any], events: EventLog) -> str:
    exec_cfg = cfg["execution"]
    exe = exec_cfg.get("exe") or _default_exe_path()
    if not os.path.isabs(exe):
        exe = os.path.join(_HERE, exe)
    exe = os.path.normpath(exe)

    auto_build = bool(exec_cfg.get("auto_build", True))
    force_rebuild = bool(exec_cfg.get("force_rebuild", False))
    needs_build = force_rebuild or (not os.path.exists(exe))
    if auto_build and os.path.exists(exe):
        needs_build = needs_build or (_newest_source_mtime() > os.path.getmtime(exe))

    if auto_build and needs_build:
        cmd = exec_cfg.get("build_command") or ["make"]
        timeout_s = int(exec_cfg.get("build_timeout_s", 600))
        log(f"building local simulator with: {cmd}")
        events.write("build_start", command=cmd, timeout_s=timeout_s)
        result = subprocess.run(
            cmd,
            cwd=_HERE,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            shell=isinstance(cmd, str),
        )
        events.write(
            "build_complete",
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "Local simulator build failed.\n"
                f"stdout:\n{result.stdout}\n\n"
                f"stderr:\n{result.stderr}"
            )

    if not os.path.exists(exe):
        raise FileNotFoundError(
            f"Simulator executable not found: {exe}. "
            "Set execution.auto_build=True or supply execution.exe."
        )
    return exe


def _parse_cli(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=str, default=None, help="JSON config override")
    ap.add_argument("--tag", type=str, default=None, help="Override run.tag")
    ap.add_argument(
        "--reference-mode",
        choices=("continuum", "large"),
        default=None,
        help="Override reference_family.mode",
    )
    ap.add_argument(
        "--reference-large-geom",
        type=int,
        nargs=4,
        metavar=("LX", "LY", "TX", "TY"),
        default=None,
        help="Override reference_family.large_geom with a single large reference lattice",
    )
    ap.add_argument("--no-monitor", action="store_true", help="Disable live monitor")
    ap.add_argument("--no-show-final", action="store_true", help="Disable final interactive plots")
    return ap.parse_args(argv)


def _resolve_axis(values, lo: float, hi: float, step: float) -> np.ndarray:
    if values is not None:
        arr = np.array([float(v) for v in values], dtype=float)
    else:
        arr = np.arange(float(lo), float(hi) + 1.0e-12, float(step), dtype=float)
    arr = np.unique(np.round(arr, 12))
    if arr.size == 0:
        raise ValueError("Empty r-grid")
    return arr


def _normalise_geom_map(raw: dict[Any, Any] | None) -> dict[int, tuple[int, int, int, int]] | None:
    if raw is None:
        return None
    out: dict[int, tuple[int, int, int, int]] = {}
    for key, value in raw.items():
        if len(value) != 4:
            raise ValueError(f"Geometry entry for {key} must have four integers")
        out[int(key)] = tuple(int(v) for v in value)
    return out


def _normalise_lattice_table(raw: Any, field_name: str) -> list[tuple[int, int, int, int, int]] | None:
    if raw is None:
        return None

    def _infer_size(Lx: int, Ly: int) -> int:
        # Use the larger linear extent as the effective size scale for 1/L fits.
        L = max(abs(int(Lx)), abs(int(Ly)))
        if L <= 0:
            raise ValueError(f"{field_name} has non-positive inferred L from (Lx, Ly)=({Lx}, {Ly})")
        return L

    out: list[tuple[int, int, int, int, int]] = []
    for i, row in enumerate(raw):
        if len(row) != 4:
            raise ValueError(f"{field_name}[{i}] must be [Lx, Ly, Tx, Ty]")
        Lx, Ly, Tx, Ty = [int(v) for v in row]
        L = _infer_size(Lx, Ly)
        out.append((L, Lx, Ly, Tx, Ty))
    return out


def _apply_explicit_lattice_tables(cfg: dict[str, Any]) -> None:
    test_table = _normalise_lattice_table(cfg["test_family"].get("lattices"), "test_family.lattices")
    if test_table is not None:
        sizes = [L for L, *_ in test_table]
        if len(set(sizes)) != len(sizes):
            raise ValueError(
                "test_family.lattices has duplicate inferred L entries; "
                "inferred L = max(|Lx|, |Ly|)"
            )
        cfg["test_family"]["sizes"] = sizes
        cfg["test_family"]["geom_map"] = {
            int(L): [int(Lx), int(Ly), int(Tx), int(Ty)]
            for L, Lx, Ly, Tx, Ty in test_table
        }

    ref_table = _normalise_lattice_table(cfg["reference_family"].get("lattices"), "reference_family.lattices")
    if ref_table is not None:
        sizes = [L for L, *_ in ref_table]
        if len(set(sizes)) != len(sizes):
            raise ValueError(
                "reference_family.lattices has duplicate inferred L entries; "
                "inferred L = max(|Lx|, |Ly|)"
            )
        cfg["reference_family"]["sizes"] = sizes
        cfg["reference_family"]["geom_map"] = {
            int(L): [int(Lx), int(Ly), int(Tx), int(Ty)]
            for L, Lx, Ly, Tx, Ty in ref_table
        }


def _geom_from_defaults(L: int, defaults: dict[str, Any]) -> tuple[int, int, int, int]:
    return (
        int(round(float(defaults.get("Lx_mult", 1.0)) * L)),
        int(round(float(defaults.get("Ly_mult", 1.0)) * L)),
        int(round(float(defaults.get("Tx_frac", 0.0)) * L)),
        int(round(float(defaults.get("Ty_frac", 0.0)) * L)),
    )


def _build_geom_map(sizes: list[int], family_cfg: dict[str, Any]) -> dict[int, tuple[int, int, int, int]]:
    geom_map = _normalise_geom_map(family_cfg.get("geom_map"))
    if geom_map is not None:
        missing = [L for L in sizes if L not in geom_map]
        if missing:
            raise ValueError(f"Missing explicit geometry entries for sizes: {missing}")
        return {int(L): tuple(geom_map[int(L)]) for L in sizes}
    defaults = family_cfg.get("geom_defaults") or {}
    return {int(L): _geom_from_defaults(int(L), defaults) for L in sizes}


def _fit_channel_in_invL(Larr: np.ndarray, y: np.ndarray, sigma: np.ndarray,
                         fit_mode: str) -> tuple[float, float, float, float, int]:
    x = 1.0 / np.asarray(Larr, dtype=float)
    y = np.asarray(y, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    n_used = int(np.count_nonzero(mask))
    if n_used == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    if fit_mode == "auto":
        degree = 2 if n_used >= 3 else 1 if n_used >= 2 else None
    elif fit_mode == "linear":
        degree = 1 if n_used >= 2 else None
    elif fit_mode == "quadratic":
        degree = 2 if n_used >= 3 else None
    else:
        raise ValueError(f"Unsupported fit_mode: {fit_mode}")

    if degree is None:
        return np.nan, np.nan, np.nan, np.nan, n_used

    xv, yv, sv = x[mask], y[mask], sigma[mask]
    X = np.vander(xv, degree + 1, increasing=True)
    w = 1.0 / (sv * sv)
    XtW = X.T * w
    try:
        cov = np.linalg.inv(XtW @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan, n_used
    beta = cov @ (XtW @ yv)
    a = float(beta[0])
    sa = float(np.sqrt(max(cov[0, 0], 0.0)))
    b = float(beta[1]) if degree >= 1 else 0.0
    c = float(beta[2]) if degree >= 2 else 0.0
    return a, sa, b, c, n_used


def _tile_sample(payload: dict[str, Any], t_fracs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    Lx, Ly, Tx, Ty = payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"]
    data = payload["data"]
    iG = _tile_interp(data, Lx, Ly, Tx, Ty, "conn", copies=2)
    iE = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)
    paths = boundary_paths(Lx, Ly, Tx, Ty)
    G = np.full((3, len(t_fracs)), np.nan)
    sG = np.full_like(G, np.nan)
    for cycle, (dm, dn) in enumerate(paths):
        ex = dm + 0.5 * dn
        ey = _SQRT3_2 * dn
        pts = np.column_stack([t_fracs * ex, t_fracs * ey])
        G[cycle] = np.asarray(iG(pts), dtype=float).ravel()
        sG[cycle] = np.abs(np.asarray(iE(pts), dtype=float).ravel())
    return G, sG


def _write_dat(path: str, header_lines: list[str], columns: list[str], rows: list[list[Any]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for header in header_lines:
            handle.write(f"# {header}\n")
        handle.write("# columns: " + " ".join(columns) + "\n")
        for row in rows:
            out = []
            for value in row:
                if isinstance(value, (int, np.integer)):
                    out.append(str(int(value)))
                elif isinstance(value, (float, np.floating)):
                    out.append(f"{float(value):.10g}")
                else:
                    out.append(str(value))
            handle.write(" ".join(out) + "\n")


def _point_path(grid_dir: str, r1: float, r2: float) -> str:
    return os.path.join(grid_dir, f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")


def _run_one_job(job: dict[str, Any]) -> dict[str, Any]:
    out_pkl = job["out_pkl"]
    if job["resume"] and os.path.exists(out_pkl):
        return {
            "label": job["label"],
            "status": "skip",
            "group_label": job["group_label"],
            "kind": job["kind"],
            "r1": job["r1"],
            "r2": job["r2"],
        }

    scratch = os.path.join(job["scratch_root"], job["label"])
    os.makedirs(scratch, exist_ok=True)
    t0 = time.time()
    try:
        beta_cfg = job["beta_c_finder"]
        beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis, scan_chi_errs = mc_engine.find_beta_c(
            job["exe"],
            job["Lx"], job["Ly"], job["Tx"], job["Ty"],
            job["r1"], job["r2"], job["k3"],
            beta_cfg["beta_lo"], beta_cfg["beta_hi"],
            n_coarse=beta_cfg["n_coarse"],
            n_refine=beta_cfg["n_refine"],
            n_refine2=beta_cfg["n_refine2"],
            n_refine3=beta_cfg["n_refine3"],
            n_traj_coarse=beta_cfg["n_traj_coarse"],
            n_traj_fine=beta_cfg["n_traj_fine"],
            data_dir=os.path.join(scratch, "scan"),
            max_shifts=beta_cfg["max_shifts"],
            jackknife=beta_cfg["jackknife"],
        )
        _, subdir = mc_engine.run_simulator(
            job["exe"],
            job["Lx"], job["Ly"], job["Tx"], job["Ty"],
            job["r1"], job["r2"], job["k3"], beta_c,
            n_traj=job["mc"]["n_traj"],
            n_skip=job["mc"]["n_skip"],
            n_therm=job["mc"]["n_therm"],
            data_dir=os.path.join(scratch, "prod"),
        )
        if subdir is None:
            raise RuntimeError("simulator produced no output directory")
        a2a = os.path.join(subdir, "two_point_all_to_all.dat")
        payload = {
            "kind": job["kind"],
            "r1": float(job["r1"]),
            "r2": float(job["r2"]),
            "k3": float(job["k3"]),
            "L": int(job["L"]),
            "Lx": int(job["Lx"]),
            "Ly": int(job["Ly"]),
            "Tx": int(job["Tx"]),
            "Ty": int(job["Ty"]),
            "beta_c": float(beta_c),
            "beta_c_sigma": float(beta_c_sigma),
            "chi_peak": float(chi_peak),
            "scan_betas": list(scan_betas),
            "scan_chis": list(scan_chis),
            "scan_chi_errs": list(scan_chi_errs),
            "n_traj_prod": int(job["mc"]["n_traj"]),
            "n_skip": int(job["mc"]["n_skip"]),
            "n_therm": int(job["mc"]["n_therm"]),
            "wall_s": float(time.time() - t0),
            "data": mc_engine.load_all_to_all(a2a),
        }
        tmp = out_pkl + ".tmp"
        with open(tmp, "wb") as handle:
            pickle.dump(payload, handle)
        os.replace(tmp, out_pkl)
        return {
            "label": job["label"],
            "status": "ok",
            "group_label": job["group_label"],
            "kind": job["kind"],
            "r1": job["r1"],
            "r2": job["r2"],
            "beta_c": float(beta_c),
            "wall": float(time.time() - t0),
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "label": job["label"],
            "status": "err",
            "group_label": job["group_label"],
            "kind": job["kind"],
            "r1": job["r1"],
            "r2": job["r2"],
            "error": str(exc),
        }


def _score_heatmap_plot(r1_values: np.ndarray, r2_values: np.ndarray,
                        score_grid: np.ndarray, out_png: str, title: str,
                        dpi: int, show: bool) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(
        score_grid,
        origin="lower",
        aspect="auto",
        extent=[r1_values.min(), r1_values.max(), r2_values.min(), r2_values.max()],
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("score S(r1,r2)")
    if np.any(np.isfinite(score_grid)):
        jj, ii = np.unravel_index(np.nanargmin(score_grid), score_grid.shape)
        ax.plot(r1_values[ii], r2_values[jj], marker="*", markersize=14,
                markerfacecolor="white", markeredgecolor="k")
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi)
    if show:
        fig.show()
        plt.pause(0.001)
    else:
        plt.close(fig)


def _score_and_zscore_heatmap_plot(r1_values: np.ndarray, r2_values: np.ndarray,
                                   score_grid: np.ndarray, zscore_grid: np.ndarray,
                                   out_png: str, title: str,
                                   figsize: list[float], dpi: int, show: bool) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)
    extent = [r1_values.min(), r1_values.max(), r2_values.min(), r2_values.max()]

    im1 = ax1.imshow(score_grid, origin="lower", aspect="auto", extent=extent)
    cb1 = fig.colorbar(im1, ax=ax1)
    cb1.set_label("score S(r1,r2)")
    if np.any(np.isfinite(score_grid)):
        jj, ii = np.unravel_index(np.nanargmin(score_grid), score_grid.shape)
        ax1.plot(r1_values[ii], r2_values[jj], marker="*", markersize=12,
                 markerfacecolor="white", markeredgecolor="k")
    ax1.set_title("Score heatmap")
    ax1.set_xlabel("r1")
    ax1.set_ylabel("r2")

    im2 = ax2.imshow(zscore_grid, origin="lower", aspect="auto", extent=extent)
    cb2 = fig.colorbar(im2, ax=ax2)
    cb2.set_label("RMS z-score")
    if np.any(np.isfinite(zscore_grid)):
        jj, ii = np.unravel_index(np.nanargmin(zscore_grid), zscore_grid.shape)
        ax2.plot(r1_values[ii], r2_values[jj], marker="*", markersize=12,
                 markerfacecolor="white", markeredgecolor="k")
    ax2.set_title("Test-vs-reference RMS z-score")
    ax2.set_xlabel("r1")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi)
    if show:
        fig.show()
        plt.pause(0.001)
    else:
        plt.close(fig)


def _build_jobs(cfg: dict[str, Any], exe: str, out_root: str,
                r1_values: np.ndarray, r2_values: np.ndarray,
                test_geom_map: dict[int, tuple[int, int, int, int]],
                ref_geom_map: dict[int, tuple[int, int, int, int]] | None):
    jobs: list[dict[str, Any]] = []
    points = [(float(r1), float(r2)) for r2 in r2_values for r1 in r1_values]
    group_totals: dict[str, int] = {}
    scratch_root = os.path.join(out_root, "_mc_scratch")
    os.makedirs(scratch_root, exist_ok=True)
    grid_cfg = cfg["grid"]

    for L in cfg["test_family"]["sizes"]:
        group_label = f"test L{L}"
        group_totals[group_label] = len(points)
        Lx, Ly, Tx, Ty = test_geom_map[L]
        grid_dir = os.path.join(out_root, "grid", f"L{L}", "test")
        os.makedirs(grid_dir, exist_ok=True)
        for r1, r2 in points:
            jobs.append({
                "label": f"test_L{L}_r1{r1:.3f}_r2{r2:.3f}",
                "kind": "test",
                "group_label": group_label,
                "L": int(L),
                "Lx": int(Lx), "Ly": int(Ly), "Tx": int(Tx), "Ty": int(Ty),
                "r1": float(r1), "r2": float(r2), "k3": float(grid_cfg["k3"]),
                "out_pkl": _point_path(grid_dir, r1, r2),
                "scratch_root": scratch_root,
                "resume": bool(cfg["paths"]["resume"]),
                "mc": cfg["mc"],
                "beta_c_finder": cfg["beta_c_finder"],
                "exe": exe,
            })

    ref_cfg = cfg["reference_family"]
    ref_r1 = float(ref_cfg["fixed_r1"])
    ref_r2 = float(ref_cfg["fixed_r2"])
    if ref_cfg["mode"] == "continuum":
        for L in ref_cfg["sizes"]:
            group_label = f"ref L{L}"
            group_totals[group_label] = 1
            Lx, Ly, Tx, Ty = ref_geom_map[L]
            grid_dir = os.path.join(out_root, "grid", f"L{L}", "ref")
            os.makedirs(grid_dir, exist_ok=True)
            jobs.append({
                "label": f"ref_L{L}_fixed_r1{ref_r1:.3f}_r2{ref_r2:.3f}",
                "kind": "ref",
                "group_label": group_label,
                "L": int(L),
                "Lx": int(Lx), "Ly": int(Ly), "Tx": int(Tx), "Ty": int(Ty),
                "r1": ref_r1, "r2": ref_r2, "k3": float(grid_cfg["k3"]),
                "out_pkl": _point_path(grid_dir, ref_r1, ref_r2),
                "scratch_root": scratch_root,
                "resume": bool(cfg["paths"]["resume"]),
                "mc": cfg["mc"],
                "beta_c_finder": cfg["beta_c_finder"],
                "exe": exe,
            })
    else:
        group_label = "ref large"
        group_totals[group_label] = 1
        Lx, Ly, Tx, Ty = [int(v) for v in ref_cfg["large_geom"]]
        grid_dir = os.path.join(out_root, "grid", "Lref_large", "ref")
        os.makedirs(grid_dir, exist_ok=True)
        jobs.append({
            "label": f"ref_large_fixed_r1{ref_r1:.3f}_r2{ref_r2:.3f}",
            "kind": "ref",
            "group_label": group_label,
            "L": -1,
            "Lx": int(Lx), "Ly": int(Ly), "Tx": int(Tx), "Ty": int(Ty),
            "r1": ref_r1, "r2": ref_r2, "k3": float(grid_cfg["k3"]),
            "out_pkl": _point_path(grid_dir, ref_r1, ref_r2),
            "scratch_root": scratch_root,
            "resume": bool(cfg["paths"]["resume"]),
            "mc": cfg["mc"],
            "beta_c_finder": cfg["beta_c_finder"],
            "exe": exe,
        })
    return jobs, points, group_totals


def _run_jobs(jobs: list[dict[str, Any]], cfg: dict[str, Any],
              state: dict[str, Any], monitor: LiveMonitor,
              events: EventLog) -> int:
    n_workers = int(cfg["execution"]["n_workers"])
    state["stage"] = "precompute"
    monitor.refresh(state, force=True)

    def _handle_result(result: dict[str, Any]) -> None:
        state["jobs_done"] += 1
        state["event_index"] += 1
        label = result["group_label"]
        done, total = state["group_counts"][label]
        state["group_counts"][label] = (done + 1, total)
        status = result["status"]
        if status == "ok":
            state["jobs_ok"] += 1
            state["last_message"] = (
                f"ok {result['label']} beta_c={result['beta_c']:.5f} wall={result['wall']:.1f}s"
            )
        elif status == "skip":
            state["jobs_skip"] += 1
            state["last_message"] = f"cached {result['label']}"
        else:
            state["jobs_err"] += 1
            state["last_message"] = f"ERR {result['label']} {result.get('error', '?')}"
        remaining = state["jobs_total"] - state["jobs_done"]
        if state["jobs_done"] > 0:
            mean = (time.time() - state["t_start"]) / state["jobs_done"]
            state["eta"] = str(timedelta(seconds=int(mean * remaining))) if remaining > 0 else "0:00:00"
        log(state["last_message"])
        events.write("job_complete", **result)
        monitor.refresh(state)

    if n_workers <= 1:
        for job in jobs:
            _handle_result(_run_one_job(job))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_run_one_job, job): job for job in jobs}
            for future in as_completed(futures):
                _handle_result(future.result())
    monitor.refresh(state, force=True)
    return state["jobs_err"]


def _load_payload(path: str) -> dict[str, Any]:
    with open(path, "rb") as handle:
        return pickle.load(handle)


def _resolve_config(config: dict[str, Any] | None = None,
                    config_path: str | None = None,
                    cli_args=None) -> dict[str, Any]:
    cfg = copy.deepcopy(CONFIG)
    chosen_config = config_path or getattr(cli_args, "config", None) or CONFIG_FILE
    if chosen_config:
        if not os.path.isabs(chosen_config):
            chosen_config = os.path.join(_HERE, chosen_config)
        _deep_update(cfg, _load_json_config(chosen_config))
    if config is not None:
        _deep_update(cfg, config)

    # Apply in-script quick switches before CLI overrides.
    if REFERENCE_MODE_OVERRIDE is not None:
        cfg["reference_family"]["mode"] = str(REFERENCE_MODE_OVERRIDE)
    if REFERENCE_LARGE_GEOM_OVERRIDE is not None:
        cfg["reference_family"]["large_geom"] = [int(v) for v in REFERENCE_LARGE_GEOM_OVERRIDE]

    if cli_args is not None:
        if cli_args.tag:
            cfg["run"]["tag"] = cli_args.tag
        if cli_args.reference_mode:
            cfg["reference_family"]["mode"] = cli_args.reference_mode
        if cli_args.reference_large_geom is not None:
            cfg["reference_family"]["large_geom"] = [int(v) for v in cli_args.reference_large_geom]
        if cli_args.no_monitor:
            cfg["monitor"]["show_live_monitor"] = False
        if cli_args.no_show_final:
            cfg["monitor"]["show_final_plots"] = False

    _apply_explicit_lattice_tables(cfg)

    if cfg["paths"]["results_dir"] in (None, ""):
        cfg["paths"]["results_dir"] = os.path.join(_HERE, "results")
    elif not os.path.isabs(cfg["paths"]["results_dir"]):
        cfg["paths"]["results_dir"] = os.path.join(_HERE, cfg["paths"]["results_dir"])

    cfg["test_family"]["sizes"] = [int(v) for v in cfg["test_family"]["sizes"]]
    if cfg["reference_family"].get("sizes") is None:
        cfg["reference_family"]["sizes"] = list(cfg["test_family"]["sizes"])
    cfg["reference_family"]["sizes"] = [int(v) for v in cfg["reference_family"]["sizes"]]
    cfg["reference_family"]["large_geom"] = [int(v) for v in cfg["reference_family"]["large_geom"]]

    # Auto-switch to large reference mode if exactly one reference lattice is defined.
    if len(cfg["reference_family"]["sizes"]) == 1:
        only_L = int(cfg["reference_family"]["sizes"][0])
        ref_geom_map = _normalise_geom_map(cfg["reference_family"].get("geom_map"))
        if ref_geom_map is not None:
            if only_L not in ref_geom_map:
                raise ValueError(
                    f"reference_family.sizes has {only_L}, but reference_family.geom_map is missing that entry"
                )
            geom = list(ref_geom_map[only_L])
        else:
            geom = list(_geom_from_defaults(only_L, cfg["reference_family"].get("geom_defaults") or {}))
        cfg["reference_family"]["mode"] = "large"
        cfg["reference_family"]["large_geom"] = [int(v) for v in geom]

    cfg["analysis"]["k_values"] = [int(v) for v in cfg["analysis"]["k_values"]]
    fit_mode = cfg["analysis"]["fit_mode"]
    if fit_mode not in {"linear", "quadratic", "auto"}:
        raise ValueError("analysis.fit_mode must be linear, quadratic, or auto")
    if cfg["reference_family"]["mode"] not in {"continuum", "large"}:
        raise ValueError("reference_family.mode must be continuum or large")
    return cfg


def main(config: dict[str, Any] | None = None,
         config_path: str | None = None,
         cli_args=None) -> dict[str, Any]:
    if cli_args is None:
        cli_args = _parse_cli()
    cfg = _resolve_config(config=config, config_path=config_path, cli_args=cli_args)

    r1_values = _resolve_axis(cfg["grid"].get("r1_values"),
                              cfg["grid"]["r_min"], cfg["grid"]["r_max"], cfg["grid"]["r_step"])
    r2_values = _resolve_axis(cfg["grid"].get("r2_values"),
                              cfg["grid"]["r_min"], cfg["grid"]["r_max"], cfg["grid"]["r_step"])
    test_geom_map = _build_geom_map(cfg["test_family"]["sizes"], cfg["test_family"])
    ref_geom_map = None
    if cfg["reference_family"]["mode"] == "continuum":
        ref_geom_map = _build_geom_map(cfg["reference_family"]["sizes"], cfg["reference_family"])

    out_root = os.path.join(cfg["paths"]["results_dir"], cfg["run"]["tag"])
    dat_dir = os.path.join(out_root, "dat")
    os.makedirs(dat_dir, exist_ok=True)
    events = EventLog(os.path.join(out_root, cfg["paths"]["progress_log_name"]))
    events.write("run_start", config=cfg)

    banner(f"boundary_opt_standalone  tag={cfg['run']['tag']}")
    exe = _ensure_simulator(cfg, events)
    jobs, points, group_totals = _build_jobs(cfg, exe, out_root, r1_values, r2_values,
                                             test_geom_map, ref_geom_map)
    group_labels = list(group_totals)
    state = {
        "stage": "build",
        "jobs_done": 0,
        "jobs_total": len(jobs),
        "jobs_ok": 0,
        "jobs_skip": 0,
        "jobs_err": 0,
        "analysis_done": 0,
        "analysis_total": len(points),
        "group_counts": {label: (0, total) for label, total in group_totals.items()},
        "score_grid": np.full((len(r2_values), len(r1_values)), np.nan),
        "zscore_grid": np.full((len(r2_values), len(r1_values)), np.nan),
        "event_index": 0,
        "last_message": "starting",
        "eta": None,
        "t_start": time.time(),
    }
    monitor = LiveMonitor(cfg, out_root, r1_values, r2_values, group_labels)
    monitor.refresh(state, force=True)

    log(f"using simulator: {exe}")
    log(f"results: {out_root}")
    log(f"test sizes: {cfg['test_family']['sizes']}")
    log(f"reference mode: {cfg['reference_family']['mode']}")
    if cfg["reference_family"]["mode"] == "continuum":
        log(f"ref sizes: {cfg['reference_family']['sizes']}")
    else:
        log(f"ref large geom: {cfg['reference_family']['large_geom']}")
        log("note: reference_family.sizes is ignored in large mode")
    log(f"grid: {len(r1_values)} x {len(r2_values)} points")

    err = _run_jobs(jobs, cfg, state, monitor, events)
    if err > 0:
        monitor.save_snapshot()
        events.write("run_abort", reason="precompute_errors", count=err)
        raise RuntimeError(f"Precompute encountered {err} job error(s); aborting analysis")

    k_values = np.array(cfg["analysis"]["k_values"], dtype=int)
    t_fracs = k_values.astype(float) / float(cfg["analysis"]["k_denominator"])
    fit_mode = cfg["analysis"]["fit_mode"]
    weighted = bool(cfg["analysis"]["weighted"])

    raw_test_rows: list[list[Any]] = []
    raw_ref_rows: list[list[Any]] = []
    cont_test_rows: list[list[Any]] = []
    cont_ref_rows: list[list[Any]] = []
    score_rows: list[list[Any]] = []
    zscore_rows: list[list[Any]] = []

    r1_index = {round(float(v), 12): i for i, v in enumerate(r1_values)}
    r2_index = {round(float(v), 12): i for i, v in enumerate(r2_values)}

    ref_cache: dict[Any, tuple[dict[str, Any], np.ndarray, np.ndarray]] = {}
    ref_cfg = cfg["reference_family"]
    ref_fixed_r1 = float(ref_cfg["fixed_r1"])
    ref_fixed_r2 = float(ref_cfg["fixed_r2"])
    if ref_cfg["mode"] == "continuum":
        for L in ref_cfg["sizes"]:
            pkl = _point_path(os.path.join(out_root, "grid", f"L{L}", "ref"), ref_fixed_r1, ref_fixed_r2)
            payload = _load_payload(pkl)
            ref_cache[int(L)] = (payload, *_tile_sample(payload, t_fracs))
    else:
        pkl = _point_path(os.path.join(out_root, "grid", "Lref_large", "ref"), ref_fixed_r1, ref_fixed_r2)
        payload = _load_payload(pkl)
        ref_cache["large"] = (payload, *_tile_sample(payload, t_fracs))

    state["stage"] = "analysis"
    state["last_message"] = "starting continuum analysis"
    monitor.refresh(state, force=True)

    for r2 in r2_values:
        for r1 in r1_values:
            test_by_L: dict[int, tuple[dict[str, Any], np.ndarray, np.ndarray]] = {}
            for L in cfg["test_family"]["sizes"]:
                pkl = _point_path(os.path.join(out_root, "grid", f"L{L}", "test"), float(r1), float(r2))
                payload = _load_payload(pkl)
                G, sG = _tile_sample(payload, t_fracs)
                test_by_L[int(L)] = (payload, G, sG)
                for cycle in range(3):
                    for ik, kval in enumerate(k_values):
                        raw_test_rows.append([
                            r1, r2, L, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                            payload["beta_c"], cycle, int(kval), float(t_fracs[ik]),
                            G[cycle, ik], sG[cycle, ik],
                        ])

            if ref_cfg["mode"] == "continuum":
                for L in ref_cfg["sizes"]:
                    payload, G, sG = ref_cache[int(L)]
                    for cycle in range(3):
                        for ik, kval in enumerate(k_values):
                            raw_ref_rows.append([
                                r1, r2, L, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                                payload["beta_c"], cycle, int(kval), float(t_fracs[ik]),
                                G[cycle, ik], sG[cycle, ik],
                            ])
            else:
                payload, G, sG = ref_cache["large"]
                for cycle in range(3):
                    for ik, kval in enumerate(k_values):
                        raw_ref_rows.append([
                            r1, r2, -1, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                            payload["beta_c"], cycle, int(kval), float(t_fracs[ik]),
                            G[cycle, ik], sG[cycle, ik],
                        ])

            Gt_inf = np.full((3, len(t_fracs)), np.nan)
            sGt_inf = np.full_like(Gt_inf, np.nan)
            Gr_inf = np.full_like(Gt_inf, np.nan)
            sGr_inf = np.full_like(Gt_inf, np.nan)

            for cycle in range(3):
                for ik, kval in enumerate(k_values):
                    Lt = sorted(test_by_L)
                    yt = np.array([test_by_L[L][1][cycle, ik] for L in Lt], dtype=float)
                    st = np.array([test_by_L[L][2][cycle, ik] for L in Lt], dtype=float)
                    at, sat, bt, ct, nt = _fit_channel_in_invL(np.array(Lt, dtype=float), yt, st, fit_mode)
                    Gt_inf[cycle, ik] = at
                    sGt_inf[cycle, ik] = sat
                    cont_test_rows.append([r1, r2, cycle, int(kval), float(t_fracs[ik]), at, sat, bt, ct, nt])

                    if ref_cfg["mode"] == "continuum":
                        Lr = sorted(k for k in ref_cache if isinstance(k, int))
                        yr = np.array([ref_cache[L][1][cycle, ik] for L in Lr], dtype=float)
                        sr = np.array([ref_cache[L][2][cycle, ik] for L in Lr], dtype=float)
                        ar, sar, br, cr, nr = _fit_channel_in_invL(np.array(Lr, dtype=float), yr, sr, fit_mode)
                        cont_ref_rows.append([r1, r2, cycle, int(kval), float(t_fracs[ik]), ar, sar, br, cr, nr])
                    else:
                        ar = float(ref_cache["large"][1][cycle, ik])
                        sar = float(ref_cache["large"][2][cycle, ik])
                        br = cr = 0.0
                        nr = 1
                        cont_ref_rows.append([r1, r2, cycle, int(kval), float(t_fracs[ik]), ar, sar, br, cr, nr])
                    Gr_inf[cycle, ik] = ar
                    sGr_inf[cycle, ik] = sar

            score = 0.0
            n_used = 0
            z2_sum = 0.0
            n_z = 0
            for cycle in range(3):
                for ik in range(len(t_fracs)):
                    gt = Gt_inf[cycle, ik]
                    st = sGt_inf[cycle, ik]
                    gr = Gr_inf[cycle, ik]
                    sr = sGr_inf[cycle, ik]
                    if not (np.isfinite(gt) and np.isfinite(gr)):
                        continue
                    diff = gt - gr
                    weight = 1.0
                    var = np.nan
                    if np.isfinite(st) and np.isfinite(sr):
                        var = st * st + sr * sr
                    if weighted and np.isfinite(var) and var > 0:
                        weight = 1.0 / var
                    score += weight * diff * diff
                    n_used += 1
                    if np.isfinite(var) and var > 0:
                        z = diff / np.sqrt(var)
                        z2_sum += z * z
                        n_z += 1

            if n_used == 0:
                score = np.nan
            z_rms = np.sqrt(z2_sum / n_z) if n_z > 0 else np.nan
            score_rows.append([r1, r2, score, n_used])
            zscore_rows.append([r1, r2, z_rms, n_z])

            ii = r1_index[round(float(r1), 12)]
            jj = r2_index[round(float(r2), 12)]
            state["score_grid"][jj, ii] = score
            state["zscore_grid"][jj, ii] = z_rms
            state["analysis_done"] += 1
            state["event_index"] += 1
            state["last_message"] = f"analyzed r1={r1:.3f} r2={r2:.3f}"
            events.write("analysis_cell_complete", r1=float(r1), r2=float(r2), score=score, zscore_rms=z_rms)
            monitor.refresh(state)

    header = [
        "Standalone boundary optimization outputs",
        f"tag={cfg['run']['tag']}",
        f"reference_mode={ref_cfg['mode']}",
        f"test_sizes={cfg['test_family']['sizes']}",
        f"ref_sizes={ref_cfg['sizes'] if ref_cfg['mode'] == 'continuum' else 'large_only'}",
        f"fit_mode={fit_mode}",
        f"weighted={weighted}",
        f"r1_values={r1_values.tolist()}",
        f"r2_values={r2_values.tolist()}",
        f"channels: cycle c=0,1,2 and t_k=k/{cfg['analysis']['k_denominator']} for k={cfg['analysis']['k_values']}",
    ]

    _write_dat(
        os.path.join(dat_dir, "raw_test_channels.dat"),
        header,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_test_rows,
    )
    _write_dat(
        os.path.join(dat_dir, "raw_ref_channels.dat"),
        header,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_ref_rows,
    )
    _write_dat(
        os.path.join(dat_dir, "continuum_test.dat"),
        header,
        ["r1", "r2", "cycle", "k", "t", "G_inf", "sigma_inf", "b_1overL", "c_1overL2", "n_used"],
        cont_test_rows,
    )
    _write_dat(
        os.path.join(dat_dir, "continuum_ref.dat"),
        header,
        ["r1", "r2", "cycle", "k", "t", "G_inf", "sigma_inf", "b_1overL", "c_1overL2", "n_used"],
        cont_ref_rows,
    )
    _write_dat(
        os.path.join(dat_dir, "score_map.dat"),
        header,
        ["r1", "r2", "score", "n_channels_used"],
        score_rows,
    )
    _write_dat(
        os.path.join(dat_dir, "zscore_map.dat"),
        header,
        ["r1", "r2", "zscore_rms", "n_channels_used"],
        zscore_rows,
    )

    show_final = bool(cfg["monitor"].get("show_final_plots", True))
    dpi = int(cfg["monitor"].get("dpi", 150))
    final_size = cfg["monitor"].get("final_size", [13.0, 5.8])
    best_summary = None
    if np.any(np.isfinite(state["score_grid"])):
        jj, ii = np.unravel_index(np.nanargmin(state["score_grid"]), state["score_grid"].shape)
        best_summary = {
            "r1_min": float(r1_values[ii]),
            "r2_min": float(r2_values[jj]),
            "score_min": float(state["score_grid"][jj, ii]),
        }
        _write_dat(
            os.path.join(dat_dir, "score_minimum.dat"),
            header,
            ["r1_min", "r2_min", "score_min"],
            [[best_summary["r1_min"], best_summary["r2_min"], best_summary["score_min"]]],
        )
        _score_heatmap_plot(
            r1_values, r2_values, state["score_grid"],
            os.path.join(out_root, "score_heatmap.png"),
            f"Boundary continuum score map  (tag={cfg['run']['tag']})",
            dpi=dpi,
            show=show_final,
        )
        _score_and_zscore_heatmap_plot(
            r1_values, r2_values, state["score_grid"], state["zscore_grid"],
            os.path.join(out_root, "score_and_zscore_heatmaps.png"),
            f"Boundary continuum diagnostics  (tag={cfg['run']['tag']})",
            figsize=final_size,
            dpi=dpi,
            show=show_final,
        )
        log(f"min score at r1={best_summary['r1_min']:.4f}, r2={best_summary['r2_min']:.4f}, S={best_summary['score_min']:.6g}")
    else:
        log("warning: score map has no finite entries")

    manifest = {
        "config": cfg,
        "derived": {
            "r1_values": r1_values.tolist(),
            "r2_values": r2_values.tolist(),
            "test_geom_map": {str(k): list(v) for k, v in test_geom_map.items()},
            "ref_geom_map": ({str(k): list(v) for k, v in ref_geom_map.items()} if ref_geom_map is not None else None),
            "exe": exe,
        },
        "summary": {
            "jobs_total": state["jobs_total"],
            "jobs_ok": state["jobs_ok"],
            "jobs_skip": state["jobs_skip"],
            "jobs_err": state["jobs_err"],
            "analysis_total": state["analysis_total"],
            "analysis_done": state["analysis_done"],
            "best": best_summary,
            "wall_s": round(time.time() - state["t_start"], 3),
        },
    }
    with open(os.path.join(out_root, "manifest_boundary_opt.json"), "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, default=_json_default)

    monitor.refresh(state, force=True)
    monitor.save_snapshot()
    events.write("run_complete", best=best_summary, wall_s=time.time() - state["t_start"])

    return {
        "config": cfg,
        "out_root": out_root,
        "best": best_summary,
        "score_grid": state["score_grid"],
        "zscore_grid": state["zscore_grid"],
    }


if __name__ == "__main__":
    RESULT = main()
    if RESULT.get("best") is not None:
        log(
            "finished successfully: "
            f"best=(r1={RESULT['best']['r1_min']:.4f}, r2={RESULT['best']['r2_min']:.4f})"
        )