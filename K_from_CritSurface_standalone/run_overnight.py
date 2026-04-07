#!/usr/bin/env python3
"""
run_overnight.py — All-in-one standalone script: reference generation,
β_c surface scan, and coupling optimisation with live progress monitor.

All configuration is in the CONFIG dictionary at the top. Edit it, then run:

    python run_overnight.py                     # use defaults in CONFIG
    python run_overnight.py config.json         # load from JSON file
    python run_overnight.py --dry-run           # print config and exit

Every phase saves incrementally so the run can be interrupted and resumed.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import mc_engine as mc
from betac_surface import BetacSurface

import matplotlib
try:
    import matplotlib.pyplot as plt
    plt.ion()  # interactive mode — enables live updates in Spyder/Qt
    _test = plt.figure(); plt.close(_test)
    # Only count as truly interactive if the backend is NOT Agg
    _backend = matplotlib.get_backend().lower()
    _INTERACTIVE = _backend not in ('agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template')
except Exception:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _INTERACTIVE = False

# =========================================================================
# DEFAULT CONFIGURATION — edit here or supply a JSON file
# =========================================================================
CONFIG = {
    # -----------------------------------------------------------------
    # Reference lattice  (equilateral triangles, possibly twisted)
    # k1 = k2 = k3 = 1.0 always.  Physics target to match.
    # -----------------------------------------------------------------
    "ref_Lx": 20,
    "ref_Ly": 20,
    "ref_Tx": 0,            # twist in x (integer lattice units)
    "ref_Ty": 0,            # twist in y
    "ref_beta_c": None,     # set to skip the β_c scan (e.g. 0.27465307)
    "ref_n_traj": 200000,   # production trajectories for reference
    "ref_n_traj_scan_coarse": 20000,
    "ref_n_traj_scan_fine":   40000,

    # -----------------------------------------------------------------
    # Test lattice  (untwisted, anisotropic: k1=r1, k2=r2, k3=1)
    # -----------------------------------------------------------------
    "test_Lx": 20,
    "test_Ly": 20,

    # -----------------------------------------------------------------
    # β_c surface scan region  (Phase 1)
    # -----------------------------------------------------------------
    "surface_r1_lo": 0.5,
    "surface_r1_hi": 2.0,
    "surface_r2_lo": 0.5,
    "surface_r2_hi": 2.0,
    "surface_n_r1": 10,     # grid points in r₁ direction
    "surface_n_r2": 10,     # grid points in r₂ direction  (total = n_r1 × n_r2)
    "surface_n_traj_coarse": 20000,
    "surface_n_traj_fine":   40000,

    # -----------------------------------------------------------------
    # Optimisation  (Phase 2) — adaptive grid search
    # -----------------------------------------------------------------
    "opt_r1_init": 1.0,     # initial grid centre
    "opt_r2_init": 1.0,
    "opt_half_span": 0.40,  # half-width of search grid in r-space
    "opt_n_grid": 5,        # points per side per level (5×5 = 25)
    "opt_max_levels": 4,    # refinement levels
    "opt_top_n": 3,         # re-score top N with Z⁴
    "opt_n_traj": 30000,    # MC trajectories per optimisation point

    # -----------------------------------------------------------------
    # General
    # -----------------------------------------------------------------
    "beta_init": 0.27,      # starting guess for β_c scans
    "output_dir": "results/overnight",
    "reuse_test_scans": True, # reuse existing scan data from disk if available

    # -----------------------------------------------------------------
    # Dashboard / video
    # -----------------------------------------------------------------
    "save_frames": False,   # save every dashboard redraw as numbered PNG
    "make_video": False,    # compile frames into mp4 (implies save_frames)
}

# =========================================================================
# Simulator path
# =========================================================================
EXE = ("bin/ising_tri_twisted_parallelogram.exe"
       if sys.platform == "win32"
       else "bin/ising_tri_twisted_parallelogram")


# =========================================================================
# Config handling
# =========================================================================

def load_config(path=None):
    """Return CONFIG dict, optionally overridden by a JSON file."""
    cfg = dict(CONFIG)
    if path and os.path.isfile(path):
        with open(path) as f:
            overrides = json.load(f)
        cfg.update(overrides)
        print(f"Config loaded from {path} ({len(overrides)} overrides)")
    return cfg


def print_config(cfg):
    print("\n" + "=" * 70)
    print("  CONFIGURATION")
    print("=" * 70)
    sections = [
        ("Reference lattice", ["ref_Lx", "ref_Ly", "ref_Tx", "ref_Ty",
                                "ref_beta_c", "ref_n_traj",
                                "ref_n_traj_scan_coarse", "ref_n_traj_scan_fine"]),
        ("Test lattice", ["test_Lx", "test_Ly"]),
        ("Surface scan", ["surface_r1_lo", "surface_r1_hi",
                          "surface_r2_lo", "surface_r2_hi",
                          "surface_n_r1", "surface_n_r2",
                          "surface_n_traj_coarse", "surface_n_traj_fine"]),
        ("Optimisation", ["opt_r1_init", "opt_r2_init", "opt_half_span",
                          "opt_n_grid", "opt_max_levels", "opt_top_n",
                          "opt_n_traj"]),
        ("General", ["beta_init", "output_dir", "reuse_test_scans"]),
    ]
    for title, keys in sections:
        print(f"\n  [{title}]")
        for k in keys:
            print(f"    {k:30s} = {cfg[k]}")
    n_surface = cfg["surface_n_r1"] * cfg["surface_n_r2"]
    n_opt = cfg["opt_n_grid"] ** 2 * cfg["opt_max_levels"]
    print(f"\n  Surface points: {n_surface}")
    print(f"  Max optimisation evals: {n_opt}")
    print("=" * 70 + "\n")


# =========================================================================
# Auto-build simulator
# =========================================================================

def build_exe():
    if os.path.isfile(EXE):
        return
    print(f"Simulator not found at '{EXE}' — building ...")
    import subprocess
    for make_cmd in (["mingw32-make"] if sys.platform == "win32" else []) + ["make"]:
        try:
            r = subprocess.run([make_cmd], cwd=_HERE, capture_output=True, text=True)
            if r.returncode == 0 and os.path.isfile(EXE):
                print(f"  Build OK via {make_cmd}")
                return
        except FileNotFoundError:
            pass
    src = os.path.join(_HERE, "src", "ising_tri_twisted_parallelogram.cc")
    inc = os.path.join(_HERE, "include")
    os.makedirs(os.path.join(_HERE, "bin"), exist_ok=True)
    cmd = ["g++", "-std=c++14", "-O3", "-Wall", "-Wno-sign-compare",
           "-Wno-unused-variable", f"-I{inc}", src, "-o",
           os.path.join(_HERE, EXE)]
    r = subprocess.run(cmd, cwd=_HERE, capture_output=True, text=True)
    if r.returncode == 0 and os.path.isfile(EXE):
        print(f"  Build OK via g++")
        return
    print(f"  Build FAILED:\n{r.stderr}")
    sys.exit(1)


# =========================================================================
# Progress monitor + live dashboard (single PNG updated across all phases)
# =========================================================================

from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy.interpolate import LinearNDInterpolator as _LNDInterp
from scipy.spatial import Delaunay as _Delaunay


# Pass colours for the β_c scan panels
_PASS_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c", 3: "#d62728"}
_PASS_LABELS = {0: "coarse", 1: "refine 1", 2: "refine 2", 3: "refine 3"}

# Minimum interval (seconds) between dashboard redraws
_DASHBOARD_MIN_INTERVAL = 2.0


class Monitor:
    """Unified single-PNG dashboard updated live across all three phases.

    Layout (3 × 3):

        [0,0] Ref β_c scan        [0,1] Test β_c scan       [0,2] Surface scatter + Delaunay
        [1,0] Surface manifold    [1,1] Z² heatmap (pixel)  [1,2] Z² convergence
        [2,0] Boundary residuals  [2,1] Z² history (log)    [2,2] Run status text

    Gold ★ = global best,  Cyan ◆ = latest evaluation.
    """

    def __init__(self, output_dir, cfg):
        self.output_dir = output_dir
        self.cfg = cfg
        os.makedirs(output_dir, exist_ok=True)
        self.t0 = time.time()
        self._logpath = os.path.join(output_dir, "run.log")
        self._logf = open(self._logpath, "a")
        self._fig = None
        self._last_draw_time = 0.0

        # ---- Phase 0 state ----
        self._phase = 0
        self._phase_desc = "Phase 0: Reference"
        self._ref_scan = None          # latest progress_cb data

        # ---- Phase 1 state ----
        self._test_scan = None         # latest progress_cb data for current test pt
        self._test_scan_label = None   # (r1, r2) of current test pt
        self._surface_pts = []         # completed (r1, r2, beta_c) points
        self._surface_total = cfg["surface_n_r1"] * cfg["surface_n_r2"]
        self._surface_grid_r1 = np.linspace(cfg["surface_r1_lo"],
                                            cfg["surface_r1_hi"],
                                            cfg["surface_n_r1"])
        self._surface_grid_r2 = np.linspace(cfg["surface_r2_lo"],
                                            cfg["surface_r2_hi"],
                                            cfg["surface_n_r2"])
        # Live interpolator rebuilt as points arrive
        self._live_interp = None
        self._live_tri = None

        # ---- Phase 2 state ----
        self.global_best_z2 = float("inf")
        self.global_best = None
        self.all_results = []
        self._level_results = []
        self._current_point = None
        self._current_level = 0
        self._current_grid_r1 = None
        self._current_grid_r2 = None
        self._history_z2 = []
        self._last_slices = None
        self._best_slices = None

        # Full surface (loaded at Phase 2 start)
        self._surf_r1 = None
        self._surf_r2 = None
        self._surf_bc = None
        self._surf_interp = None
        self._surf_tri = None

        # ---- Trajectory (configuration) counters ----
        # Estimate scan traj per find_beta_c call:
        #   coarse: 11 * n_traj_coarse; fine passes: ~21 * n_traj_fine
        ref_scan_est = (11 * cfg["ref_n_traj_scan_coarse"]
                        + 21 * cfg["ref_n_traj_scan_fine"])
        self._ref_configs_done = 0
        self._ref_configs_total = ref_scan_est + cfg["ref_n_traj"]

        # Per-run test lattice progress (current MC run only)
        self._test_run_label = ""     # e.g. "Scan (0.6,0.8)" or "Opt L1 #5"
        self._test_run_done = 0       # traj done in current run
        self._test_run_total = 0      # traj expected in current run

        # ---- Frame saving ----
        self._save_frames = False
        self._frames_dir = None
        self._frame_count = 0

    # ==================================================================
    # Helpers
    # ==================================================================

    def _elapsed(self):
        s = time.time() - self.t0
        m, s = divmod(int(s), 60)
        h, m = divmod(m, 60)
        return f"{h}h{m:02d}m{s:02d}s"

    def _log(self, msg):
        print(msg)
        self._logf.write(msg + "\n")
        self._logf.flush()

    def _should_redraw(self):
        """Rate-limit dashboard redraws."""
        now = time.time()
        if now - self._last_draw_time < _DASHBOARD_MIN_INTERVAL:
            return False
        self._last_draw_time = now
        return True

    # ==================================================================
    # Scan callbacks (called from within find_beta_c via progress_cb)
    # ==================================================================

    def ref_scan_cb(self, data):
        """progress_cb for Phase 0 reference β_c scan."""
        self._ref_scan = data
        self._ref_configs_done = data.get("traj_done", 0)
        if self._should_redraw():
            self._update_dashboard()

    def make_test_scan_cb(self, r1, r2):
        """Return a closure suitable as progress_cb for a Phase 1 surface scan."""
        self._test_scan_label = (r1, r2)
        self._test_scan = None  # reset
        # Estimate total traj for one find_beta_c call
        scan_est = (11 * self.cfg["surface_n_traj_coarse"]
                    + 21 * self.cfg["surface_n_traj_fine"])
        self._test_run_label = f"Scan ({r1:.3f},{r2:.3f})"
        self._test_run_done = 0
        self._test_run_total = scan_est

        def _cb(data):
            self._test_scan = data
            self._test_run_done = data.get("traj_done", 0)
            if self._should_redraw():
                self._update_dashboard()
        return _cb

    # ==================================================================
    # Phase 0: Reference
    # ==================================================================

    def phase0_start(self, cfg):
        self._phase = 0
        self._phase_desc = "Phase 0: Reference β_c scan"
        self._log(f"\n{'='*70}")
        self._log("  PHASE 0: Reference generation")
        self._log(f"  Lattice: {cfg['ref_Lx']}×{cfg['ref_Ly']}  "
                  f"Tx={cfg['ref_Tx']} Ty={cfg['ref_Ty']}")
        self._log(f"{'='*70}")
        self._update_dashboard()

    def phase0_production_start(self, n_traj):
        """Called just before the reference production run."""
        self._phase_desc = "Phase 0: Ref production"
        self._update_dashboard()

    def phase0_done(self, beta_c, elapsed, n_traj_prod=0):
        self._ref_configs_done += n_traj_prod
        self._log(f"  β_c = {beta_c:.10f}  ({elapsed:.0f}s)")
        self._update_dashboard()

    # ==================================================================
    # Phase 1: Surface scan
    # ==================================================================

    def phase1_start(self, cfg, n_total, n_done):
        self._phase = 1
        self._phase_desc = f"Phase 1: Surface scan (0/{n_total})"
        self._log(f"\n{'='*70}")
        self._log("  PHASE 1: β_c surface scan")
        self._log(f"  Test lattice: {cfg['test_Lx']}×{cfg['test_Ly']}")
        self._log(f"  r₁ ∈ [{cfg['surface_r1_lo']}, {cfg['surface_r1_hi']}]  "
                  f"r₂ ∈ [{cfg['surface_r2_lo']}, {cfg['surface_r2_hi']}]")
        self._log(f"  Grid: {cfg['surface_n_r1']}×{cfg['surface_n_r2']} "
                  f"= {n_total} points  ({n_done} already done)")
        self._log(f"{'='*70}")
        self._log(f"  {'#':>5s}  {'r1':>8s}  {'r2':>8s}  {'β_c':>12s}  "
                  f"{'time':>6s}  {'elapsed':>10s}")
        self._log(f"  {'-'*5}  {'-'*8}  {'-'*8}  {'-'*12}  "
                  f"{'-'*6}  {'-'*10}")
        self._update_dashboard()

    def phase1_point(self, idx, total, r1, r2, beta_c, dt):
        """Log a completed surface point and update live surface data."""
        self._log(f"  {idx:>5d}/{total:<5d} {r1:8.5f}  {r2:8.5f}  "
                  f"{beta_c:12.8f}  {dt:5.0f}s  {self._elapsed()}")
        self._surface_pts.append((r1, r2, beta_c))
        self._phase_desc = f"Phase 1: Surface scan ({idx}/{total})"
        # Rebuild live interpolator when we have enough non-degenerate points
        self._rebuild_live_interp()
        self._update_dashboard()

    def phase1_done(self, n_pts, path):
        self._log(f"\n  Surface complete: {n_pts} points → {path}")
        self._log(f"  Phase 1 elapsed: {self._elapsed()}")
        self._update_dashboard()

    def _rebuild_live_interp(self):
        """Rebuild Delaunay + interpolator from accumulated surface points."""
        pts = self._surface_pts
        if len(pts) < 3:
            return
        r1 = np.array([p[0] for p in pts])
        r2 = np.array([p[1] for p in pts])
        bc = np.array([p[2] for p in pts])
        try:
            tri = _Delaunay(np.column_stack([r1, r2]))
            self._live_tri = tri
            self._live_interp = _LNDInterp(tri, bc)
        except Exception:
            pass  # not enough non-degenerate points yet

    # ==================================================================
    # Phase 2: Optimisation (surface injection + grid search)
    # ==================================================================

    def set_surface(self, surface_obj):
        """Inject completed BetacSurface for Phase 2 dashboard panels."""
        self._surf_r1 = surface_obj.r1_arr
        self._surf_r2 = surface_obj.r2_arr
        self._surf_bc = surface_obj.beta_c_arr
        self._surf_tri = surface_obj._tri
        self._surf_interp = surface_obj._interp

    def phase2_level_start(self, level, r1_vals, r2_vals,
                           r1c, r2c, hs1, hs2):
        self._phase = 2
        self._level_results = []
        self._current_level = level
        self._current_grid_r1 = np.array(r1_vals)
        self._current_grid_r2 = np.array(r2_vals)
        n = len(r1_vals) * len(r2_vals)
        self._phase_desc = (f"Phase 2: Opt level {level}  "
                            f"({r1c:.3f},{r2c:.3f})±({hs1:.3f},{hs2:.3f})")
        self._log(f"\n{'='*70}")
        self._log(f"  LEVEL {level}  centre=({r1c:.4f}, {r2c:.4f})  "
                  f"half-span=({hs1:.4f}, {hs2:.4f})")
        self._log(f"  r₁ ∈ [{r1_vals[0]:.4f}, {r1_vals[-1]:.4f}]  "
                  f"r₂ ∈ [{r2_vals[0]:.4f}, {r2_vals[-1]:.4f}]  ({n} pts)")
        self._log(f"{'='*70}")
        self._log(f"  {'#':>4s}  {'r1':>8s}  {'r2':>8s}  {'β_c':>10s}  "
                  f"{'Z²':>10s}  {'time':>6s}  {'best_Z²':>10s}  status")
        self._log(f"  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*10}  "
                  f"{'-'*10}  {'-'*6}  {'-'*10}  {'-'*8}")

    def phase2_opt_run_start(self, idx, r1, r2, n_traj, level):
        """Called just before an optimizer production run."""
        self._test_run_label = f"Opt L{level} #{idx}"
        self._test_run_done = 0
        self._test_run_total = n_traj
        self._update_dashboard()

    def phase2_opt_run_done(self, n_traj):
        """Called right after an optimizer production run finishes."""
        self._test_run_done = n_traj
        # Don't redraw here - phase2_point will do it

    def phase2_point(self, idx, r1, r2, beta_c, z2, dt, inside, level,
                     slices=None, n_traj=0):
        self._test_run_done = n_traj  # mark current run complete
        is_best = z2 < self.global_best_z2
        if is_best:
            self.global_best_z2 = z2
            self.global_best = {"r1": r1, "r2": r2, "beta_c": beta_c, "z2": z2}
        status = "★ BEST" if is_best else ("(extrap)" if not inside else "")
        rec = {"r1": r1, "r2": r2, "beta_c": beta_c, "z2": z2, "level": level}
        self._level_results.append(rec)
        self.all_results.append(rec)
        self._history_z2.append(z2)
        self._current_point = (r1, r2)
        if slices is not None:
            self._last_slices = slices
            if is_best:
                self._best_slices = slices
        self._log(f"  {idx:>4d}  {r1:8.5f}  {r2:8.5f}  {beta_c:10.7f}  "
                  f"{z2:10.5f}  {dt:5.0f}s  {self.global_best_z2:10.5f}  "
                  f"{status}")
        self._update_dashboard()

    def phase2_level_end(self, level, best, top_results):
        self._log(f"\n  Level {level} best:  r1={best['r1']:.6f}  "
                  f"r2={best['r2']:.6f}  Z²={best['z2']:.5f}  "
                  f"Z⁴={best.get('z4', float('nan')):.5f}")
        if len(top_results) > 1:
            self._log(f"  Top-{len(top_results)} Z⁴ re-scores:")
            for i, tr in enumerate(top_results):
                self._log(f"    #{i+1}  r1={tr['r1']:.5f}  r2={tr['r2']:.5f}"
                          f"  Z²={tr['z2']:.5f}  Z⁴={tr['z4']:.5f}")
        self._log(f"  Elapsed: {self._elapsed()}")

    # ==================================================================
    # Live dashboard — 3×3 single PNG
    # ==================================================================

    def _update_dashboard(self):
        """Redraw all panels and save to dashboard_live.png."""
        try:
            # Close previous figure so Spyder inline backend registers
            # each redraw as a new plot (otherwise Plots pane won't update).
            if self._fig is not None:
                plt.close(self._fig)
            self._fig = plt.figure(figsize=(22, 14))
            fig = self._fig
            gs = GridSpec(3, 3, figure=fig, hspace=0.38, wspace=0.30)

            ax00 = fig.add_subplot(gs[0, 0])
            self._draw_betac_scan(ax00, self._ref_scan,
                                  "Ref β_c scan")

            ax01 = fig.add_subplot(gs[0, 1])
            lbl = ""
            if self._test_scan_label:
                lbl = f" (r₁={self._test_scan_label[0]:.3f}, r₂={self._test_scan_label[1]:.3f})"
            self._draw_betac_scan(ax01, self._test_scan,
                                  f"Test β_c scan{lbl}")

            ax02 = fig.add_subplot(gs[0, 2])
            self._draw_surface_scatter(ax02)

            ax10 = fig.add_subplot(gs[1, 0])
            self._draw_manifold(ax10)

            ax11 = fig.add_subplot(gs[1, 1])
            self._draw_pixel_heatmap(ax11)

            # Status panel spans rows 1–2 in column 2
            ax_status = fig.add_subplot(gs[1:3, 2])
            self._draw_status(ax_status)

            ax20 = fig.add_subplot(gs[2, 0])
            self._draw_boundary_residuals(ax20)

            ax21 = fig.add_subplot(gs[2, 1])
            self._draw_z2_convergence_combined(ax21)

            fig.suptitle(f"K_from_CritSurface — {self._phase_desc}  |  "
                         f"{self._elapsed()}", fontsize=14, y=0.995)

            path = os.path.join(self.output_dir, "dashboard_live.png")
            fig.savefig(path, dpi=110, bbox_inches="tight")

            # Save numbered frame if enabled
            if self._save_frames and self._frames_dir:
                frame_path = os.path.join(
                    self._frames_dir, f"frame_{self._frame_count:05d}.png")
                fig.savefig(frame_path, dpi=110, bbox_inches="tight")
                self._frame_count += 1

            # Update live display (Spyder Plots pane / interactive window)
            if _INTERACTIVE:
                plt.show(block=False)
                plt.pause(0.01)
            else:
                # Fallback: push saved PNG to IPython/Spyder console
                try:
                    from IPython.display import display, Image as IPImage
                    display(IPImage(filename=path))
                except Exception:
                    pass
        except Exception as e:
            self._log(f"    [dashboard ERROR] {e}")

    def enable_frame_saving(self):
        """Enable saving numbered frames for video compilation."""
        self._save_frames = True
        self._frames_dir = os.path.join(self.output_dir, "frames")
        os.makedirs(self._frames_dir, exist_ok=True)
        self._frame_count = 0
        self._log(f"  Frame saving enabled → {self._frames_dir}/")

    def make_video(self, fps=10):
        """Compile saved frames into an mp4 video using ffmpeg."""
        if not self._save_frames or self._frame_count == 0:
            self._log("  No frames to compile.")
            return None
        video_path = os.path.join(self.output_dir, "dashboard_video.mp4")
        pattern = os.path.join(self._frames_dir, "frame_%05d.png")
        cmd = [
            "ffmpeg", "-y",
            "-framerate", str(fps),
            "-i", pattern,
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            video_path,
        ]
        self._log(f"  Compiling {self._frame_count} frames → {video_path}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=300)
            if result.returncode == 0:
                self._log(f"  Video saved: {video_path}")
                return video_path
            else:
                self._log(f"  ffmpeg failed: {result.stderr[:500]}")
                return None
        except FileNotFoundError:
            self._log("  ffmpeg not found — frames saved but video skipped")
            return None
        except Exception as e:
            self._log(f"  Video compilation error: {e}")
            return None

    # ------------------------------------------------------------------
    # Panel: β_c susceptibility scan (Phase 0 ref / Phase 1 test)
    # ------------------------------------------------------------------

    def _draw_betac_scan(self, ax, scan_data, title):
        if scan_data is None:
            ax.set_title(title)
            ax.text(0.5, 0.5, "waiting …", transform=ax.transAxes,
                    ha="center", va="center", fontsize=11, color="gray")
            ax.set_xlabel("β")
            ax.set_ylabel("χ")
            return

        betas = np.array(scan_data["all_betas"])
        chis = np.array(scan_data["all_chis"])
        chi_errs = np.array(scan_data.get("all_chi_errs",
                                          [0.0] * len(betas)))
        pids = np.array(scan_data["pass_ids"])

        if len(betas) > 0:
            for p in sorted(set(pids)):
                mask = pids == p
                ax.errorbar(betas[mask], chis[mask], yerr=chi_errs[mask],
                            fmt="o", color=_PASS_COLORS.get(p, "gray"),
                            ms=4, capsize=2, capthick=0.5, elinewidth=0.5,
                            label=_PASS_LABELS.get(p, f"p{p}"),
                            zorder=3 + p, alpha=0.85)

        # GC fit curve
        gc = scan_data.get("gc_params")
        if gc is not None and len(betas) > 0:
            b_fine = np.linspace(betas.min(), betas.max(), 300)
            chi_fit = mc._gram_charlier(b_fine, *gc)
            ax.plot(b_fine, chi_fit, "k-", lw=1.5, alpha=0.6,
                    label="GC fit")

        be = scan_data.get("beta_estimate")
        if be is not None:
            ax.axvline(be, color="red", ls="--", lw=1.2, alpha=0.7)
            if len(betas) == 0:
                ax.text(0.5, 0.5, f"β_c = {be:.8f}\n(cached)",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=11, color="green")

        pn = scan_data.get("pass_num", "?")
        xlbl = f"β   (est={be:.6f})" if be is not None else "β"
        ax.set_xlabel(xlbl)
        ax.set_ylabel("χ")

        # Progress bar in title
        prog = scan_data.get("scan_progress")
        prog_txt = ""
        if prog and prog["pts_total"] > 0:
            d, t = prog["pts_done"], prog["pts_total"]
            n_blk = 12
            filled = int(min(d / t, 1.0) * n_blk)
            prog_txt = f"  {'█' * filled}{'░' * (n_blk - filled)} {d}/{t}"

        ax.set_title(f"{title} — pass {pn}{prog_txt}")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.2)

    # ------------------------------------------------------------------
    # Panel: Surface scatter + Delaunay (populates live in Phase 1)
    # ------------------------------------------------------------------

    def _draw_surface_scatter(self, ax):
        # Pending grid positions (gray empty circles)
        for r1 in self._surface_grid_r1:
            for r2 in self._surface_grid_r2:
                ax.plot(r1, r2, "o", color="lightgray", ms=5,
                        markeredgecolor="silver", zorder=1)

        pts = self._surface_pts
        if not pts:
            ax.set_title(f"Surface scatter (0/{self._surface_total})")
            ax.set_xlabel("r₁"); ax.set_ylabel("r₂")
            ax.grid(True, alpha=0.2)
            return

        r1a = np.array([p[0] for p in pts])
        r2a = np.array([p[1] for p in pts])
        bca = np.array([p[2] for p in pts])

        sc = ax.scatter(r1a, r2a, c=bca, cmap="viridis", s=45,
                        edgecolors="k", linewidths=0.5, zorder=3)
        plt.colorbar(sc, ax=ax, shrink=0.82, label="β_c")

        # Delaunay edges
        if self._live_tri is not None:
            for simplex in self._live_tri.simplices:
                il = list(simplex) + [simplex[0]]
                ax.plot(r1a[il], r2a[il], "k-", alpha=0.2, lw=0.5)

        # Mark current test point being scanned
        if self._phase == 1 and self._test_scan_label:
            ax.plot(self._test_scan_label[0], self._test_scan_label[1],
                    "D", color="cyan", ms=10, markeredgecolor="k",
                    markeredgewidth=1.0, zorder=6, label="scanning")
            ax.legend(fontsize=7)

        ax.set_title(f"Surface scatter ({len(pts)}/{self._surface_total})")
        ax.set_xlabel("r₁ = k₁/k₃")
        ax.set_ylabel("r₂ = k₂/k₃")
        ax.grid(True, alpha=0.2)

    # ------------------------------------------------------------------
    # Panel: Interpolated β_c manifold
    # ------------------------------------------------------------------

    def _draw_manifold(self, ax):
        # Choose interpolator: full surface if available, else live
        interp = self._surf_interp or self._live_interp
        if interp is None:
            ax.set_title("β_c manifold (building…)")
            ax.text(0.5, 0.5, "need ≥ 3 surface pts",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=10, color="gray")
            return

        if self._surf_r1 is not None:
            r1_lim = (self._surf_r1.min(), self._surf_r1.max())
            r2_lim = (self._surf_r2.min(), self._surf_r2.max())
        else:
            pts = self._surface_pts
            r1a = [p[0] for p in pts]
            r2a = [p[1] for p in pts]
            r1_lim = (min(r1a), max(r1a))
            r2_lim = (min(r2a), max(r2a))

        r1f = np.linspace(r1_lim[0], r1_lim[1], 80)
        r2f = np.linspace(r2_lim[0], r2_lim[1], 80)
        R1, R2 = np.meshgrid(r1f, r2f)
        BC = interp(R1.ravel(), R2.ravel()).reshape(R1.shape)

        im = ax.pcolormesh(R1, R2, BC, cmap="viridis", shading="auto")
        plt.colorbar(im, ax=ax, shrink=0.82, label="β_c")

        # Surface data
        if self._surf_r1 is not None:
            ax.scatter(self._surf_r1, self._surf_r2, c="white", s=12,
                       edgecolors="k", linewidths=0.3, zorder=4)
        elif self._surface_pts:
            r1a = [p[0] for p in self._surface_pts]
            r2a = [p[1] for p in self._surface_pts]
            ax.scatter(r1a, r2a, c="white", s=12, edgecolors="k",
                       linewidths=0.3, zorder=4)

        # Opt evals
        if self.all_results:
            ax.scatter([r["r1"] for r in self.all_results],
                       [r["r2"] for r in self.all_results],
                       c="gray", s=16, marker="s", edgecolors="k",
                       linewidths=0.2, zorder=5, alpha=0.5)

        # Current (cyan ◆) + Best (gold ★)
        if self._current_point:
            ax.plot(self._current_point[0], self._current_point[1],
                    "D", color="cyan", ms=11, markeredgecolor="k",
                    markeredgewidth=0.8, zorder=8)
        if self.global_best:
            ax.plot(self.global_best["r1"], self.global_best["r2"],
                    "*", color="gold", ms=16, markeredgecolor="k",
                    markeredgewidth=0.7, zorder=9)

        ax.set_xlabel("r₁")
        ax.set_ylabel("r₂")
        ax.set_title("β_c(r₁,r₂) manifold")

    # ------------------------------------------------------------------
    # Panel: Pixelated Z² heatmap (Phase 2 current level)
    # ------------------------------------------------------------------

    def _draw_pixel_heatmap(self, ax):
        results = self.all_results
        if not results:
            ax.set_title("Z² heatmap (Phase 2)")
            ax.text(0.5, 0.5, "waiting for opt…", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            return

        # Fixed zoom: full surface scan extent
        r1_lo = self.cfg["surface_r1_lo"]
        r1_hi = self.cfg["surface_r1_hi"]
        r2_lo = self.cfg["surface_r2_lo"]
        r2_hi = self.cfg["surface_r2_hi"]

        # Colour scale from all Z² values
        z2_all = [r["z2"] for r in results]
        import matplotlib.colors as mcolors
        norm = mcolors.Normalize(vmin=min(z2_all), vmax=max(z2_all))
        cmap = plt.get_cmap("viridis_r")

        # Draw each level's points as rectangles sized to that level's grid
        # spacing (later levels paint over earlier ones → smaller pixels on top)
        levels = sorted(set(r["level"] for r in results))
        for lv in levels:
            lv_pts = [r for r in results if r["level"] == lv]
            r1v = sorted(set(r["r1"] for r in lv_pts))
            r2v = sorted(set(r["r2"] for r in lv_pts))
            dr1 = ((r1v[-1] - r1v[0]) / max(len(r1v) - 1, 1)
                   if len(r1v) > 1 else (r1_hi - r1_lo) / self.cfg["opt_n_grid"])
            dr2 = ((r2v[-1] - r2v[0]) / max(len(r2v) - 1, 1)
                   if len(r2v) > 1 else (r2_hi - r2_lo) / self.cfg["opt_n_grid"])
            for r in lv_pts:
                c = cmap(norm(r["z2"]))
                rect = Rectangle((r["r1"] - dr1 / 2, r["r2"] - dr2 / 2),
                                 dr1, dr2, facecolor=c,
                                 edgecolor="k", linewidth=0.3,
                                 zorder=3 + lv)
                ax.add_patch(rect)

        # Colourbar via a ScalarMappable
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, shrink=0.82, label="Z²")

        # Current level search box
        if (self._current_grid_r1 is not None
                and len(self._current_grid_r1) > 1):
            gr1lo, gr1hi = self._current_grid_r1[0], self._current_grid_r1[-1]
            gr2lo, gr2hi = self._current_grid_r2[0], self._current_grid_r2[-1]
            rect = Rectangle((gr1lo, gr2lo),
                              gr1hi - gr1lo, gr2hi - gr2lo,
                              fill=False, edgecolor="royalblue", lw=1.5,
                              ls="--", zorder=10)
            ax.add_patch(rect)

        if self._current_point:
            ax.plot(self._current_point[0], self._current_point[1],
                    "D", color="cyan", ms=11, markeredgecolor="k",
                    markeredgewidth=0.8, zorder=11, label="current")
        if self.global_best:
            ax.plot(self.global_best["r1"], self.global_best["r2"],
                    "*", color="gold", ms=16, markeredgecolor="k",
                    markeredgewidth=0.7, zorder=12, label="best")

        ax.set_xlim(r1_lo, r1_hi)
        ax.set_ylim(r2_lo, r2_hi)
        ax.set_aspect("auto")
        ax.set_xlabel("r₁"); ax.set_ylabel("r₂")
        ax.set_title(f"Z² heatmap (L{self._current_level}, "
                     f"{len(results)} evals)")
        ax.legend(fontsize=7, loc="upper right")

    # ------------------------------------------------------------------
    # Panel: Z² convergence combined (scatter + running best, log scale)
    # ------------------------------------------------------------------

    def _draw_z2_convergence_combined(self, ax):
        n = len(self._history_z2)
        if n == 0:
            ax.set_title("Z² convergence")
            ax.text(0.5, 0.5, "waiting for opt…", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            return
        z2s = np.array(self._history_z2)
        running = np.minimum.accumulate(z2s)
        iters = np.arange(1, n+1)
        ax.semilogy(iters, z2s, "o", ms=3, alpha=0.4, color="steelblue",
                    label="Z² per pt", zorder=3)
        ax.semilogy(iters, running, "-", color="red", lw=2,
                    label="running best", zorder=4)
        levels = sorted(set(r["level"] for r in self.all_results))
        cum = 0
        for lv in levels:
            cnt = sum(1 for r in self.all_results if r["level"] == lv)
            if cum > 0:
                ax.axvline(cum+0.5, color="gray", ls="--", lw=0.7)
                ax.text(cum+1, z2s.max()*0.9, f"L{lv}", fontsize=7, va="top")
            cum += cnt
        ax.set_xlabel("Evaluation #"); ax.set_ylabel("Z² (log)")
        ax.set_title("Z² convergence"); ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # Panel: Boundary-slice residuals
    # ------------------------------------------------------------------

    def _draw_boundary_residuals(self, ax):
        dc = ["#d62728", "#1f77b4", "#2ca02c"]
        dl = ["v", "u", "w"]
        ax.axhline(0, ls="-", color="gray", alpha=0.4, lw=0.8)
        drawn = False
        if self._best_slices:
            for i, s in enumerate(self._best_slices):
                if len(s["t"]) > 1:
                    ax.fill_between(s["t"], s["diff"]-s["diff_err"],
                                    s["diff"]+s["diff_err"],
                                    color=dc[i], alpha=0.08)
                    ax.plot(s["t"], s["diff"], "--", color=dc[i],
                            alpha=0.5, lw=1.2, label=f"best {dl[i]}")
                    drawn = True
        if self._last_slices:
            for i, s in enumerate(self._last_slices):
                if len(s["t"]) > 1:
                    ax.fill_between(s["t"], s["diff"]-s["diff_err"],
                                    s["diff"]+s["diff_err"],
                                    color=dc[i], alpha=0.15)
                    ax.plot(s["t"], s["diff"], "-", color=dc[i],
                            lw=1.6, label=f"latest {dl[i]}")
                    drawn = True
        ax.set_xlabel("t (fraction along path)")
        ax.set_ylabel("G_test − G_ref")
        ax.set_title("Boundary-slice residuals")
        ax.grid(True, alpha=0.2)
        if drawn:
            ax.legend(fontsize=7, ncol=2)
        else:
            ax.text(0.5, 0.5, "waiting for opt…", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")

    # ------------------------------------------------------------------
    # Panel: Status text (double-height, rows 1-2 col 2)
    # ------------------------------------------------------------------

    def _draw_status(self, ax):
        ax.axis("off")

        def _bar(done, total, width=20):
            if total <= 0:
                return f"{done}/??"
            frac = min(done / total, 1.0)
            filled = int(frac * width)
            pct = frac * 100
            return f"{'█' * filled}{'░' * (width - filled)} {pct:4.1f}%"

        def _fmt_configs(n):
            if n >= 1_000_000:
                return f"{n/1e6:.2f}M"
            elif n >= 1_000:
                return f"{n/1e3:.1f}k"
            return str(n)

        lines = [
            f"Phase:   {self._phase_desc}",
            f"Elapsed: {self._elapsed()}",
            "",
            "─── Lattice parameters ───",
            f"Ref:  {self.cfg['ref_Lx']}×{self.cfg['ref_Ly']}  "
            f"Tx={self.cfg['ref_Tx']} Ty={self.cfg['ref_Ty']}  "
            f"k=(1,1,1)",
            f"Test: {self.cfg['test_Lx']}×{self.cfg['test_Ly']}  "
            f"Tx=0 Ty=0  k=(r₁,r₂,1)",
            "",
            "─── MC configurations ───",
            "",
            f"Ref lattice  ({self.cfg['ref_Lx']}×{self.cfg['ref_Ly']}):",
            f"  {_bar(self._ref_configs_done, self._ref_configs_total)}",
            f"  {_fmt_configs(self._ref_configs_done)} / "
            f"~{_fmt_configs(self._ref_configs_total)} configs",
        ]

        # Ref β_c estimate
        if self._ref_scan:
            be = self._ref_scan.get("beta_estimate")
            if be is not None:
                lines.append(f"  β_c = {be:.8f}")

        lines += [
            "",
            f"Test lattice ({self.cfg['test_Lx']}×{self.cfg['test_Ly']}):",
        ]

        # Current MC run progress (per-run, not aggregate)
        if self._test_run_total > 0:
            lines.append(f"  {self._test_run_label}:")
            lines.append(f"  {_bar(self._test_run_done, self._test_run_total)}")
            lines.append(f"  {_fmt_configs(self._test_run_done)} / "
                         f"~{_fmt_configs(self._test_run_total)} traj")
        else:
            lines.append("  (idle)")

        # Surface + current scan details
        s_done = len(self._surface_pts)
        lines.append(f"  Surface: {s_done}/{self._surface_total} pts")
        if self._test_scan and self._test_scan_label:
            r1, r2 = self._test_scan_label
            pn = self._test_scan.get("pass_num", "?")
            prog = self._test_scan.get("scan_progress")
            tl = f"  Scanning ({r1:.3f},{r2:.3f}) p{pn}"
            if prog and prog["pts_total"] > 0:
                d, t = prog["pts_done"], prog["pts_total"]
                tl += f" [{d}/{t}]"
            lines.append(tl)

        lines.append("")

        if self.global_best:
            b = self.global_best
            lines += [
                "─── Best result ───",
                f"  r₁  = {b['r1']:.6f}",
                f"  r₂  = {b['r2']:.6f}",
                f"  β_c = {b['beta_c']:.8f}",
                f"  Z²  = {b['z2']:.6f}",
            ]

        if self.all_results:
            n_grid = 0
            if (self._current_grid_r1 is not None
                    and self._current_grid_r2 is not None):
                n_grid = len(self._current_grid_r1) * len(self._current_grid_r2)
            n_done = len(self._level_results)
            lines.append(f"\nOpt L{self._current_level}: "
                         f"{_bar(n_done, n_grid, 15)}")
            lines.append(f"Total evals: {len(self.all_results)}")

        txt = "\n".join(lines)
        ax.text(0.05, 0.97, txt, transform=ax.transAxes,
                va="top", ha="left", fontsize=9, family="monospace",
                bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow",
                          ec="gray", alpha=0.9))
        ax.set_title("Run status")

    # ==================================================================
    # Final outputs
    # ==================================================================

    def final_summary(self):
        self._phase_desc = "COMPLETE"
        self._log(f"\n{'='*70}")
        self._log("  RUN COMPLETE")
        self._log(f"{'='*70}")
        if self.global_best:
            b = self.global_best
            self._log(f"  Best:  r₁ = {b['r1']:.6f}")
            self._log(f"         r₂ = {b['r2']:.6f}")
            self._log(f"         β_c = {b['beta_c']:.8f}")
            self._log(f"         Z²  = {b['z2']:.6f}")
        self._log(f"  Total points evaluated: {len(self.all_results)}")
        self._log(f"  Total wall time: {self._elapsed()}")
        self._log(f"  Log file: {self._logpath}")
        self._log(f"{'='*70}")
        # Final high-res dashboard
        self._last_draw_time = 0  # force redraw
        self._update_dashboard()
        if self._fig is not None:
            path = os.path.join(self.output_dir, "dashboard_final.png")
            self._fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(self._fig)
            self._log(f"  Dashboard (final) → {path}")
        self._save_convergence_standalone()
        self._logf.close()

    def _save_convergence_standalone(self):
        if not self.all_results:
            return
        z2s = [r["z2"] for r in self.all_results]
        running = np.minimum.accumulate(z2s)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(1, len(z2s)+1), z2s, "o", ms=3, alpha=0.5,
                label="Z² per point")
        ax.plot(range(1, len(z2s)+1), running, "r-", lw=2, label="running best")
        levels = sorted(set(r["level"] for r in self.all_results))
        idx = 0
        for lv in levels:
            n = sum(1 for r in self.all_results if r["level"] == lv)
            if idx > 0:
                ax.axvline(idx+0.5, color="gray", ls="--", lw=0.8)
                ax.text(idx+1, max(z2s)*0.9, f"L{lv}", fontsize=8)
            idx += n
        ax.set_xlabel("Evaluation #"); ax.set_ylabel("Z² cost")
        ax.set_title("Optimiser convergence"); ax.legend(); fig.tight_layout()
        path = os.path.join(self.output_dir, "convergence.png")
        fig.savefig(path, dpi=120); plt.close(fig)
        self._log(f"  Convergence plot → {path}")


# =========================================================================
# Phase 0: Generate reference
# =========================================================================

def _check_lattice(saved, cfg, label, ref=False):
    """Compare lattice dimensions in *saved* dict against current config.

    Returns True if they match, False (with a printed warning) otherwise.
    *saved* should have keys Lx, Ly (and optionally Tx, Ty).
    If *ref* is True, compare against ref_Lx/ref_Ly/ref_Tx/ref_Ty;
    otherwise compare against test_Lx/test_Ly (Tx=0, Ty=0).
    """
    if ref:
        want = (cfg["ref_Lx"], cfg["ref_Ly"], cfg["ref_Tx"], cfg["ref_Ty"])
    else:
        want = (cfg["test_Lx"], cfg["test_Ly"], 0, 0)
    got = (saved.get("Lx"), saved.get("Ly"),
           saved.get("Tx", 0), saved.get("Ty", 0))
    if got != want:
        print(f"  WARNING: {label} lattice mismatch — "
              f"saved {got[0]}×{got[1]} T{got[2]}x{got[3]} "
              f"vs config {want[0]}×{want[1]} T{want[2]}x{want[3]}  "
              f"→ ignoring cached data")
        return False
    return True


def _reconstruct_scan_from_dir(scan_dir, beta_c, *, expect_Lx=None,
                                expect_Ly=None, expect_Tx=None,
                                expect_Ty=None):
    """Rebuild scan data dict by reading .dat files from a beta_scan directory.

    If expect_Lx/Ly/Tx/Ty are given, only .dat files whose header matches
    those dimensions are included (others are silently skipped).
    """
    import glob
    betas, chis, chi_errs = [], [], []
    dat_files = glob.glob(os.path.join(scan_dir, "**", "*.dat"), recursive=True)
    for fpath in dat_files:
        if os.path.basename(fpath) in ("two_point_all_to_all.dat",
                                       "two_point_typed.dat"):
            continue
        b_val = chi_val = chi_err_val = None
        file_Lx = file_Ly = file_Tx = file_Ty = None
        try:
            with open(fpath) as f:
                for line in f:
                    if line.startswith("beta "):
                        b_val = float(line.split()[1])
                    elif line.startswith("m_susc "):
                        parts = line.split()
                        chi_val = float(parts[1])
                        chi_err_val = float(parts[2]) if len(parts) > 2 else 0.0
                    elif line.startswith("L_x "):
                        file_Lx = int(line.split()[1])
                    elif line.startswith("L_y "):
                        file_Ly = int(line.split()[1])
                    elif line.startswith("T_x "):
                        file_Tx = int(line.split()[1])
                    elif line.startswith("T_y "):
                        file_Ty = int(line.split()[1])
        except Exception:
            continue
        # Lattice-dimension filter
        if expect_Lx is not None and file_Lx is not None:
            if (file_Lx != expect_Lx or file_Ly != expect_Ly
                    or (file_Tx or 0) != (expect_Tx or 0)
                    or (file_Ty or 0) != (expect_Ty or 0)):
                continue  # skip file from a different lattice
        if b_val is not None and chi_val is not None:
            betas.append(b_val)
            chis.append(chi_val)
            chi_errs.append(chi_err_val or 0.0)
    if not betas:
        return None
    # Sort by beta
    order = np.argsort(betas)
    betas = list(np.array(betas)[order])
    chis = list(np.array(chis)[order])
    chi_errs = list(np.array(chi_errs)[order])
    # Fit GC curve — use the chi-peak beta as initial guess, not the hint
    beta_peak = betas[int(np.argmax(chis))]
    gc_params = None
    beta_est = beta_peak
    try:
        gc_tuple, mode = mc._gc_fit(betas, chis, beta_peak)
        if gc_tuple is not None:
            gc_params = list(gc_tuple)
            beta_est = mode
    except Exception:
        pass
    return {
        "pass_num": "done",
        "all_betas": betas,
        "all_chis": chis,
        "all_chi_errs": chi_errs,
        "pass_ids": [0] * len(betas),  # all shown as one pass
        "gc_params": gc_params,
        "beta_estimate": beta_est,
        "scan_progress": {"pts_done": len(betas), "pts_total": len(betas)},
    }


def phase0(cfg, mon):
    out = cfg["output_dir"]
    ref_dir = os.path.join(out, "ref_data")
    meta_path = os.path.join(ref_dir, "ref_metadata.json")

    if os.path.isfile(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
        if not _check_lattice(meta, cfg, "Phase 0 ref_metadata", ref=True):
            # Lattice mismatch — regenerate
            mon._log("  Phase 0: cached reference has wrong lattice dims — regenerating")
        else:
            mon._log(f"  Phase 0: reference exists → β_c={meta['beta_c']:.10f}")
            # Restore scan data for the dashboard
            scan_cache = os.path.join(ref_dir, "ref_scan_data.json")
            if os.path.isfile(scan_cache):
                with open(scan_cache) as f:
                    mon._ref_scan = json.load(f)
            else:
                # Reconstruct from .dat files in the scan directory
                scan_dir = os.path.join(ref_dir, "beta_scan")
                reconstructed = _reconstruct_scan_from_dir(
                    scan_dir, meta["beta_c"],
                    expect_Lx=cfg["ref_Lx"], expect_Ly=cfg["ref_Ly"],
                    expect_Tx=cfg["ref_Tx"], expect_Ty=cfg["ref_Ty"])
                if reconstructed is not None:
                    mon._ref_scan = reconstructed
                    mon._log(f"    Reconstructed {len(reconstructed['all_betas'])}"
                             f" scan points from {scan_dir}")
                    # Save cache for next time
                    with open(scan_cache, "w") as f:
                        json.dump(reconstructed, f)
                else:
                    mon._ref_scan = {
                        "pass_num": "done",
                        "all_betas": [], "all_chis": [], "all_chi_errs": [],
                        "pass_ids": [], "gc_params": None,
                        "beta_estimate": meta["beta_c"],
                        "scan_progress": {"pts_done": 1, "pts_total": 1},
                    }
            mon._ref_configs_done = mon._ref_configs_total  # 100%
            return meta

    os.makedirs(ref_dir, exist_ok=True)
    mon.phase0_start(cfg)

    t0 = time.time()
    if cfg["ref_beta_c"] is not None:
        beta_c = cfg["ref_beta_c"]
        mon._log(f"  Using provided β_c = {beta_c:.10f}")
    else:
        mon._log("  Scanning for β_c ...")
        margin = max(0.20 * cfg["beta_init"], 0.04)
        beta_lo = max(0.01, cfg["beta_init"] - margin)
        beta_hi = cfg["beta_init"] + margin
        scan_dir = os.path.join(ref_dir, "beta_scan")
        beta_c, chi_peak, _, _, _ = mc.find_beta_c(
            EXE, cfg["ref_Lx"], cfg["ref_Ly"],
            cfg["ref_Tx"], cfg["ref_Ty"],
            k1=1.0, k2=1.0, k3=1.0,
            beta_lo=beta_lo, beta_hi=beta_hi,
            n_traj_coarse=cfg["ref_n_traj_scan_coarse"],
            n_traj_fine=cfg["ref_n_traj_scan_fine"],
            data_dir=scan_dir,
            progress_cb=mon.ref_scan_cb)
        mon._log(f"  Found β_c = {beta_c:.10f}  (χ={chi_peak:.4f})")

    mon._log(f"  Production run ({cfg['ref_n_traj']} trajectories) ...")
    mon.phase0_production_start(cfg["ref_n_traj"])
    prod_dir = os.path.join(ref_dir, "ref_production")
    stdout, subdir = mc.run_simulator(
        EXE, cfg["ref_Lx"], cfg["ref_Ly"],
        cfg["ref_Tx"], cfg["ref_Ty"],
        k1=1.0, k2=1.0, k3=1.0, beta=beta_c,
        n_traj=cfg["ref_n_traj"], n_therm=5000, data_dir=prod_dir)
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    mon.phase0_done(beta_c, time.time() - t0, n_traj_prod=cfg["ref_n_traj"])

    meta = {
        "beta_c": beta_c,
        "Lx": cfg["ref_Lx"], "Ly": cfg["ref_Ly"],
        "Tx": cfg["ref_Tx"], "Ty": cfg["ref_Ty"],
        "k1": 1.0, "k2": 1.0, "k3": 1.0,
        "n_traj": cfg["ref_n_traj"],
        "a2a_path": os.path.relpath(a2a_path, _HERE),
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    # Cache scan data for dashboard on restart
    if mon._ref_scan is not None:
        scan_cache = os.path.join(ref_dir, "ref_scan_data.json")
        # Convert numpy types to native Python for JSON serialisation
        cache = {}
        for k, v in mon._ref_scan.items():
            if isinstance(v, np.ndarray):
                cache[k] = v.tolist()
            elif isinstance(v, (list, tuple)) and len(v) > 0:
                cache[k] = [float(x) if isinstance(x, (np.floating, np.integer))
                            else x for x in v]
            elif isinstance(v, (np.floating, np.integer)):
                cache[k] = float(v)
            else:
                cache[k] = v
        # gc_params may be a numpy array or tuple of numpy scalars
        gcp = cache.get("gc_params")
        if gcp is not None:
            if isinstance(gcp, np.ndarray):
                cache["gc_params"] = gcp.tolist()
            elif isinstance(gcp, (list, tuple)):
                cache["gc_params"] = [float(x) for x in gcp]
        with open(scan_cache, "w") as f:
            json.dump(cache, f)
    return meta


# =========================================================================
# Phase 1: β_c surface
# =========================================================================

def phase1(cfg, mon):
    out = cfg["output_dir"]
    surface_dir = os.path.join(out, "betac_surfaces")
    os.makedirs(surface_dir, exist_ok=True)
    tag = f"surface_{cfg['test_Lx']}x{cfg['test_Ly']}_T0x0"
    surface_path = os.path.join(surface_dir, f"{tag}.json")

    n_r1 = cfg["surface_n_r1"]
    n_r2 = cfg["surface_n_r2"]
    total = n_r1 * n_r2

    r1_vals = np.linspace(cfg["surface_r1_lo"], cfg["surface_r1_hi"], n_r1)
    r2_vals = np.linspace(cfg["surface_r2_lo"], cfg["surface_r2_hi"], n_r2)

    # Load existing points (incremental resume)
    completed = {}
    if os.path.isfile(surface_path):
        with open(surface_path) as f:
            raw = json.load(f)
        if not _check_lattice(raw, cfg, "Phase 1 surface"):
            raw = {"points": []}  # discard mismatched data
        for p in raw.get("points", []):
            completed[(round(p["r1"], 8), round(p["r2"], 8))] = p["beta_c"]
        if len(completed) >= total:
            mon._log(f"  Phase 1: surface complete → {surface_path}")
            # Restore surface points for the dashboard
            for (r1, r2), bc in sorted(completed.items()):
                mon._surface_pts.append((r1, r2, bc))
            mon._rebuild_live_interp()
            return surface_path

    mon.phase1_start(cfg, total, len(completed))

    scan_base = os.path.join(out, "surface_scans")
    last_beta_c = cfg["beta_init"]
    done = len(completed)

    for r1 in r1_vals:
        for r2 in r2_vals:
            key = (round(r1, 8), round(r2, 8))
            if key in completed:
                last_beta_c = completed[key]
                continue
            done += 1
            k1, k2, k3 = r1, r2, 1.0
            label = f"r1_{r1:.4f}_r2_{r2:.4f}"
            scan_dir = os.path.join(scan_base, label)

            # Try to recover β_c from existing scan data on disk
            if cfg.get("reuse_test_scans", True) and os.path.isdir(scan_dir):
                recovered = _reconstruct_scan_from_dir(
                    scan_dir, last_beta_c,
                    expect_Lx=cfg["test_Lx"], expect_Ly=cfg["test_Ly"],
                    expect_Tx=0, expect_Ty=0)
                if recovered is not None and recovered["beta_estimate"] is not None:
                    beta_c = recovered["beta_estimate"]
                    mon._log(f"  {done:>5d}/{total:<5d} {r1:8.5f}  {r2:8.5f}  "
                             f"{beta_c:12.8f}  (reused)  {mon._elapsed()}")
                    # Populate test scan panel so dashboard shows the data
                    mon._test_scan = recovered
                    mon._test_scan_label = (r1, r2)
                    mon._surface_pts.append((r1, r2, beta_c))
                    mon._rebuild_live_interp()
                    completed[key] = beta_c
                    last_beta_c = beta_c
                    _save_surface(surface_path, cfg, r1_vals, r2_vals, completed)
                    mon._phase_desc = f"Phase 1: Surface scan ({done}/{total})"
                    mon._update_dashboard()
                    continue

            margin = max(0.20 * last_beta_c, 0.04)
            beta_lo = max(0.01, last_beta_c - margin)
            beta_hi = last_beta_c + margin

            test_cb = mon.make_test_scan_cb(r1, r2)
            t0 = time.time()
            try:
                beta_c, chi_peak, _, _, _ = mc.find_beta_c(
                    EXE, cfg["test_Lx"], cfg["test_Ly"], 0, 0, k1, k2, k3,
                    beta_lo, beta_hi,
                    n_traj_coarse=cfg["surface_n_traj_coarse"],
                    n_traj_fine=cfg["surface_n_traj_fine"],
                    data_dir=scan_dir,
                    progress_cb=test_cb)
            except Exception as e:
                mon._log(f"  FAILED at r1={r1:.4f} r2={r2:.4f}: {e}")
                continue
            dt = time.time() - t0
            mon.phase1_point(done, total, r1, r2, beta_c, dt)

            completed[key] = beta_c
            last_beta_c = beta_c
            _save_surface(surface_path, cfg, r1_vals, r2_vals, completed)

    _save_surface(surface_path, cfg, r1_vals, r2_vals, completed)
    mon.phase1_done(len(completed), surface_path)

    # Generate surface plot
    try:
        surf = BetacSurface(surface_path)
        plot_path = os.path.join(surface_dir, f"{tag}.png")
        surf.plot(plot_path)
        mon._log(f"  Surface plot → {plot_path}")
    except Exception as e:
        mon._log(f"  Surface plot failed: {e}")

    return surface_path


def _save_surface(path, cfg, r1_vals, r2_vals, completed):
    points = [{"r1": r1, "r2": r2, "beta_c": bc}
              for (r1, r2), bc in sorted(completed.items())]
    data = {
        "Lx": cfg["test_Lx"], "Ly": cfg["test_Ly"], "Tx": 0, "Ty": 0,
        "r1_vals": [round(float(v), 8) for v in r1_vals],
        "r2_vals": [round(float(v), 8) for v in r2_vals],
        "n_r1": len(r1_vals), "n_r2": len(r2_vals),
        "n_traj_scan_coarse": cfg["surface_n_traj_coarse"],
        "n_traj_scan_fine": cfg["surface_n_traj_fine"],
        "points": points,
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


# =========================================================================
# Phase 2: Optimisation
# =========================================================================

def phase2(cfg, surface_path, ref_meta, mon):
    surface = BetacSurface(surface_path)
    surface.summary()
    mon.set_surface(surface)

    ref_data = mc.load_all_to_all(ref_meta["a2a_path"])
    mon._log(f"\n  Reference loaded: β_c={ref_meta['beta_c']:.10f}")
    mon._log(f"\n{'='*70}")
    mon._log("  PHASE 2: Optimisation")
    mon._log(f"  Grid: {cfg['opt_n_grid']}×{cfg['opt_n_grid']}  "
             f"Levels: {cfg['opt_max_levels']}  "
             f"Top-{cfg['opt_top_n']} Z⁴ re-score")
    mon._log(f"  MC per point: {cfg['opt_n_traj']} trajectories")
    mon._log(f"{'='*70}")

    out = cfg["output_dir"]
    opt_dir = os.path.join(out, "optimisation")
    run_dir = os.path.join(opt_dir, "runs")
    os.makedirs(run_dir, exist_ok=True)

    r1c = cfg["opt_r1_init"]
    r2c = cfg["opt_r2_init"]
    # Cap initial half-span to the surface extent so refinement is
    # meaningful (otherwise clamping makes every level the same range).
    max_hs1 = (cfg["surface_r1_hi"] - cfg["surface_r1_lo"]) / 2.0
    max_hs2 = (cfg["surface_r2_hi"] - cfg["surface_r2_lo"]) / 2.0
    hs1 = min(cfg["opt_half_span"], max_hs1)
    hs2 = min(cfg["opt_half_span"], max_hs2)
    ng = cfg["opt_n_grid"]

    for level in range(cfg["opt_max_levels"]):
        # --- Check for completed level on disk (resume support) ---
        level_json = os.path.join(opt_dir, f"grid_level_{level:02d}.json")
        if os.path.isfile(level_json):
            with open(level_json) as f:
                saved = json.load(f)
            # Validate lattice dimensions if present in saved JSON
            if "test_Lx" in saved:
                want = (cfg["ref_Lx"], cfg["ref_Ly"],
                        cfg["ref_Tx"], cfg["ref_Ty"],
                        cfg["test_Lx"], cfg["test_Ly"])
                got = (saved.get("ref_Lx"), saved.get("ref_Ly"),
                       saved.get("ref_Tx", 0), saved.get("ref_Ty", 0),
                       saved.get("test_Lx"), saved.get("test_Ly"))
                if got != want:
                    mon._log(f"  Level {level}: SKIPPING — lattice mismatch "
                             f"(saved ref {got[0]}×{got[1]} test {got[4]}×{got[5]} "
                             f"vs config ref {want[0]}×{want[1]} test {want[4]}×{want[5]})")
                    break  # can't trust any further levels either
            # Restore monitor state from saved grid
            for gk, gv in saved["grid"].items():
                rec = {"r1": gv["r1"], "r2": gv["r2"],
                       "beta_c": gv["beta_c"], "z2": gv["z2"],
                       "level": level}
                mon.all_results.append(rec)
                mon._history_z2.append(gv["z2"])
                if gv["z2"] < mon.global_best_z2:
                    mon.global_best_z2 = gv["z2"]
                    mon.global_best = {"r1": gv["r1"], "r2": gv["r2"],
                                       "beta_c": gv["beta_c"],
                                       "z2": gv["z2"]}
            # Advance centre / half-span for next level
            r1c = saved["best_r1"]
            r2c = saved["best_r2"]
            # Determine whether the best was on the border (translate)
            # or interior (refine).  Reconstruct from saved half-spans
            # vs config: if hs shrank, it was interior; else border.
            if level == 0:
                hs1 = saved["hs_r1"]
                hs2 = saved["hs_r2"]
            # Interior → halve; border → keep.  We can't perfectly tell,
            # so just read the next level's file if it exists, or infer
            # from saved hs vs current hs.
            next_json = os.path.join(opt_dir,
                                     f"grid_level_{level+1:02d}.json")
            if os.path.isfile(next_json):
                with open(next_json) as f2:
                    nxt = json.load(f2)
                hs1 = nxt["hs_r1"]
                hs2 = nxt["hs_r2"]
                r1c = nxt["r1_center"]
                r2c = nxt["r2_center"]
            else:
                # Heuristic: halve if best wasn't on grid boundary
                hs1 = saved["hs_r1"] / 2.0
                hs2 = saved["hs_r2"] / 2.0
            mon._current_level = level
            mon._phase = 2
            mon._phase_desc = (f"Phase 2: Resumed — loaded L{level}  "
                               f"best Z²={saved['best_z2']:.5f}")
            mon._log(f"  Level {level}: loaded {len(saved['grid'])} pts "
                     f"from {level_json}  (best Z²="
                     f"{saved['best_z2']:.5f})")
            mon._last_draw_time = 0  # force redraw
            mon._update_dashboard()
            continue

        r1_lo = max(r1c - hs1, cfg["surface_r1_lo"])
        r1_hi = min(r1c + hs1, cfg["surface_r1_hi"])
        r2_lo = max(r2c - hs2, cfg["surface_r2_lo"])
        r2_hi = min(r2c + hs2, cfg["surface_r2_hi"])

        r1_vals = np.linspace(r1_lo, r1_hi, ng)
        r2_vals = np.linspace(r2_lo, r2_hi, ng)

        mon.phase2_level_start(level, r1_vals, r2_vals, r1c, r2c, hs1, hs2)

        grid = {}
        idx = 0
        for i, r1 in enumerate(r1_vals):
            for j, r2 in enumerate(r2_vals):
                idx += 1
                k1, k2, k3 = r1, r2, 1.0
                beta_c = surface(r1, r2)
                inside = surface.domain_contains(r1, r2)

                t0 = time.time()
                label = f"L{level}_i{i}_j{j}"
                prod_dir = os.path.join(run_dir, f"level{level:02d}",
                                        f"prod_{label}")
                mon.phase2_opt_run_start(idx, r1, r2,
                                         cfg["opt_n_traj"], level)
                try:
                    stdout, subdir = mc.run_simulator(
                        EXE, cfg["test_Lx"], cfg["test_Ly"], 0, 0,
                        k1, k2, k3, beta_c,
                        n_traj=cfg["opt_n_traj"], n_therm=3000,
                        data_dir=prod_dir)
                except Exception as e:
                    mon._log(f"  {idx:>4d}  FAILED: {e}")
                    continue
                mon.phase2_opt_run_done(cfg["opt_n_traj"])

                a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
                test_data = mc.load_all_to_all(a2a_path)

                chi2, ndof, _ = mc.boundary_slices_normed(
                    ref_data, test_data,
                    cfg["test_Lx"], cfg["test_Ly"], 0, 0,
                    cfg["ref_Lx"], cfg["ref_Ly"],
                    cfg["ref_Tx"], cfg["ref_Ty"])
                z2 = chi2 / max(ndof, 1)
                dt = time.time() - t0

                # Compute boundary slices for dashboard residual panel
                slices = mc.extract_boundary_slices(
                    ref_data, test_data,
                    cfg["test_Lx"], cfg["test_Ly"], 0, 0,
                    cfg["ref_Lx"], cfg["ref_Ly"],
                    cfg["ref_Tx"], cfg["ref_Ty"])

                mon.phase2_point(idx, r1, r2, beta_c, z2, dt, inside, level,
                                 slices=slices, n_traj=cfg["opt_n_traj"])

                grid[(i, j)] = {
                    "r1": r1, "r2": r2, "beta_c": beta_c,
                    "z2": z2, "a2a_path": a2a_path,
                }

        if not grid:
            mon._log("  WARNING: no points evaluated!")
            continue

        # Re-score top N with Z⁴
        ranked = sorted(grid.keys(), key=lambda k: grid[k]["z2"])
        top_keys = ranked[:cfg["opt_top_n"]]
        top_results = []
        for key in top_keys:
            pt = grid[key]
            test_data = mc.load_all_to_all(pt["a2a_path"])
            chi2_q, ndof_q, _ = mc.boundary_slices_quartic(
                ref_data, test_data,
                cfg["test_Lx"], cfg["test_Ly"], 0, 0,
                cfg["ref_Lx"], cfg["ref_Ly"],
                cfg["ref_Tx"], cfg["ref_Ty"])
            pt["z4"] = chi2_q / max(ndof_q, 1)
            top_results.append({"r1": pt["r1"], "r2": pt["r2"],
                                "z2": pt["z2"], "z4": pt["z4"]})

        best_key = min(top_keys,
                       key=lambda k: grid[k].get("z4", grid[k]["z2"]))
        best = grid[best_key]
        mon.phase2_level_end(level, best, top_results)

        # Save level JSON
        with open(os.path.join(opt_dir, f"grid_level_{level:02d}.json"), "w") as f:
            json.dump({
                "level": level,
                "ref_Lx": cfg["ref_Lx"], "ref_Ly": cfg["ref_Ly"],
                "ref_Tx": cfg["ref_Tx"], "ref_Ty": cfg["ref_Ty"],
                "test_Lx": cfg["test_Lx"], "test_Ly": cfg["test_Ly"],
                "r1_center": float(r1c), "r2_center": float(r2c),
                "hs_r1": float(hs1), "hs_r2": float(hs2),
                "best_r1": float(best["r1"]), "best_r2": float(best["r2"]),
                "best_z2": float(best["z2"]),
                "best_z4": float(best.get("z4", best["z2"])),
                "best_beta_c": float(best["beta_c"]),
                "grid": {f"{k[0]}_{k[1]}": {
                    "r1": v["r1"], "r2": v["r2"],
                    "beta_c": v["beta_c"], "z2": v["z2"]}
                    for k, v in grid.items()},
            }, f, indent=2)

        # Translate or refine
        bi, bj = best_key
        r1c = float(best["r1"])
        r2c = float(best["r2"])
        if bi == 0 or bi == ng - 1:
            mon._log(f"  r₁: border → translate to {r1c:.4f}")
        else:
            hs1 /= 2.0
            mon._log(f"  r₁: interior → refine, hs={hs1:.4f}")
        if bj == 0 or bj == ng - 1:
            mon._log(f"  r₂: border → translate to {r2c:.4f}")
        else:
            hs2 /= 2.0
            mon._log(f"  r₂: interior → refine, hs={hs2:.4f}")

    # Save final result
    if mon.global_best:
        fit = {
            "best_r1": mon.global_best["r1"],
            "best_r2": mon.global_best["r2"],
            "best_beta_c": mon.global_best["beta_c"],
            "best_z2": mon.global_best["z2"],
            "config": cfg,
            "all_results": mon.all_results,
        }
        with open(os.path.join(opt_dir, "fit_result.json"), "w") as f:
            json.dump(fit, f, indent=2, default=float)


# =========================================================================
# Main
# =========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Overnight pipeline: reference → β_c surface → optimise")
    parser.add_argument("config_file", nargs="?", default=None,
                        help="JSON file with config overrides")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print config and exit without running")
    parser.add_argument("--save-frames", action="store_true",
                        help="Save every dashboard redraw as a numbered PNG")
    parser.add_argument("--video", action="store_true",
                        help="Save frames and compile an mp4 video (implies --save-frames)")
    args = parser.parse_args()

    cfg = load_config(args.config_file)
    print_config(cfg)

    if args.dry_run:
        print("Dry run — exiting.")
        return

    build_exe()
    os.makedirs(cfg["output_dir"], exist_ok=True)

    # Save config snapshot
    with open(os.path.join(cfg["output_dir"], "config.json"), "w") as f:
        json.dump(cfg, f, indent=2)

    mon = Monitor(cfg["output_dir"], cfg)

    do_frames = args.save_frames or args.video or cfg.get("save_frames") or cfg.get("make_video")
    do_video  = args.video or cfg.get("make_video")

    if do_frames:
        mon.enable_frame_saving()

    ref_meta = phase0(cfg, mon)
    surface_path = phase1(cfg, mon)
    phase2(cfg, surface_path, ref_meta, mon)

    mon.final_summary()

    if do_video:
        mon.make_video()

    print(f"\nAll output in: {cfg['output_dir']}/")


if __name__ == "__main__":
    main()
