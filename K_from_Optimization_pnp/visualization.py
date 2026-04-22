"""
visualization.py — live frame-writing for the optimizer benchmark.

Two plotter classes:

  * BetaScanPlotter  — consumes the dicts emitted by mc_engine.find_beta_c's
                       progress_cb; saves a frame per pass showing the
                       Gram-Charlier fit overlaid on the (beta, chi) data.

  * OptimizerPlotter — called by evaluator.evaluate after each call; saves a
                       4-panel frame: (r1, r2) trajectory, cost history,
                       beta_c history, current curve collapse.

Everything is written to disk as PNG via matplotlib's Agg backend.  No
windows, no GUI — safe for headless / CI / batch use.

Disable by passing --no-vis (callers should accept None and skip).
"""

from __future__ import annotations

import os
import queue
import threading
from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
#  Background render thread — savefig never blocks the optimizer loop.
# ---------------------------------------------------------------------------

_render_queue: queue.Queue = queue.Queue()
_SENTINEL = object()


def _render_worker():
    """Worker: each item is a zero-arg callable that creates, saves, closes
    a figure entirely within this thread.  All matplotlib calls happen here;
    the main thread never touches plt after queuing."""
    while True:
        item = _render_queue.get()
        if item is _SENTINEL:
            break
        try:
            item()  # callable: create + savefig + close
        except Exception as exc:
            print(f"[vis] render error: {exc}")
        _render_queue.task_done()


_render_thread = threading.Thread(target=_render_worker, daemon=True,
                                   name="vis-render")
_render_thread.start()


def _queue_render(fn) -> None:
    """Queue a zero-arg render callable; returns immediately."""
    _render_queue.put(fn)


def flush_render_queue() -> None:
    """Block until all queued frames have been written to disk."""
    _render_queue.join()

from mc_engine import _gram_charlier  # re-use the exact fit formula
from cost import boundary_paths, _tile_interp, _triangular_xy


# ---------------------------------------------------------------------------
#  Beta_c scan plotter
# ---------------------------------------------------------------------------

class BetaScanPlotter:
    """Save one PNG per pass of the GC beta_c scan.

    Use as the ``progress_cb`` argument of ``mc_engine.find_beta_c``.

    Parameters
    ----------
    output_dir : str
        Directory in which to write ``betac_scan_<label>_pass<N>.png``.
    label : str
        Label that identifies which evaluator call this scan belongs to
        (e.g. ``"eval0042_r1_1.0123_r2_0.9876"``).
    every_point : bool, default False
        If True, write a frame after every MC point (slow, makes nice GIFs).
        If False (default), write only at the end of each pass.
    """

    def __init__(self, output_dir: str, label: str,
                 every_point: bool = False):
        self.output_dir = output_dir
        self.label = label
        self.every_point = every_point
        os.makedirs(output_dir, exist_ok=True)
        self._last_pass = -1
        self._frame_idx = 0

    # The progress_cb signature: cb(state_dict)
    def __call__(self, state: dict) -> None:
        pass_num = int(state["pass_num"])
        new_pass = (pass_num != self._last_pass)
        self._last_pass = pass_num
        if not (self.every_point or new_pass):
            return
        self._render(state)

    # ------------------------------------------------------------------
    def _render(self, state: dict) -> None:
        # Copy all data before handing off to background thread.
        betas = np.asarray(state["all_betas"], dtype=float).copy()
        chis  = np.asarray(state["all_chis"],  dtype=float).copy()
        errs  = np.asarray(state["all_chi_errs"], dtype=float).copy()
        pids  = np.asarray(state["pass_ids"],  dtype=int).copy()
        gc       = state.get("gc_params")
        b_est    = state.get("beta_estimate")
        pass_num = int(state["pass_num"])
        label    = self.label
        frame_idx = self._frame_idx
        out_path  = os.path.join(self.output_dir, "betac_scan.png")
        self._frame_idx += 1

        def _do_render():
            fig, ax = plt.subplots(figsize=(7.0, 4.5))
            cmap = plt.get_cmap("viridis")
            n_passes_seen = max(int(pids.max()) + 1, 1) if len(pids) else 1
            for p in range(n_passes_seen):
                mask = pids == p
                if not mask.any():
                    continue
                ax.errorbar(betas[mask], chis[mask], yerr=errs[mask],
                            fmt="o", ms=5, lw=0,
                            ecolor=cmap(p / max(n_passes_seen - 1, 1)),
                            color=cmap(p / max(n_passes_seen - 1, 1)),
                            elinewidth=1.0, capsize=2,
                            label=f"pass {p}  (n={int(mask.sum())})")
            if gc is not None and len(betas) >= 4:
                b_fit = np.linspace(betas.min(), betas.max(), 400)
                try:
                    y_fit = _gram_charlier(b_fit, *gc)
                    ax.plot(b_fit, y_fit, "-", color="crimson", lw=1.5,
                            label="GC fit")
                except Exception:
                    pass
            if b_est is not None and np.isfinite(b_est):
                ax.axvline(b_est, color="black", ls="--", lw=1.0,
                           label=f"β_c ≈ {b_est:.5f}")
            ax.set_xlabel("β")
            ax.set_ylabel("χ (susceptibility)")
            ax.set_title(f"Gram-Charlier scan — {label}  (pass {pass_num})")
            ax.legend(loc="best", fontsize=8, framealpha=0.9)
            ax.grid(alpha=0.25)
            fig.text(0.99, 0.99, f"frame {frame_idx}",
                     ha="right", va="top", fontsize=8, color="gray",
                     transform=fig.transFigure)
            fig.tight_layout()
            fig.savefig(out_path, dpi=110)
            plt.close(fig)

        _queue_render(_do_render)


# ---------------------------------------------------------------------------
#  Optimizer trajectory plotter
# ---------------------------------------------------------------------------

class OptimizerPlotter:
    """Save one PNG per evaluator call summarising the optimizer state.

    Parameters
    ----------
    output_dir : str
        Directory for ``optimizer_<method>_step<N>.png`` frames.
    method : str
        Identifier for the optimizer (e.g. "nelder_mead").
    ref_data : dict
        Reference correlator dataset (output of cost.load_correlator).
        Used for the curve-collapse panel.
    Lx, Ly, Tx, Ty : int
        Lattice dimensions of the reference run (for boundary paths).
    n_copies : int, default 1
        Number of periodic tile copies to use in interpolation.
    """

    def __init__(self, output_dir: str, method: str,
                 ref_data: dict, Lx: int, Ly: int, Tx: int, Ty: int,
                 n_copies: int = 1, save_every: int = 5):
        self.output_dir = output_dir
        self.method = method
        self.ref_data = ref_data
        self.Lx = Lx
        self.Ly = Ly
        self.Tx = Tx
        self.Ty = Ty
        self.n_copies = n_copies
        # Write a frame every `save_every` evaluations (+ always the first).
        # Set save_every=1 to recover the old every-eval behaviour.
        self.save_every = max(1, int(save_every))
        os.makedirs(output_dir, exist_ok=True)

        # Per-iteration history.
        self.r1_hist: list[float] = []
        self.r2_hist: list[float] = []
        self.cost_hist: list[float] = []
        self.sigma_hist: list[float] = []
        self.beta_hist: list[float] = []
        self._step = 0
        # Simplex history: each entry is a list of 3 [r1, r2] vertices.
        # Populated by run_nelder_mead via evaluator.current_simplex.
        self.simplex_hist: list = []
        self.current_simplex: Optional[list] = None

    # ------------------------------------------------------------------
    def update(self, r1: float, r2: float, cost: float, sigma_cost: float,
               beta_c: float, test_data: Optional[dict] = None,
               simplex: Optional[list] = None) -> None:
        """Record one evaluation; write a frame every ``save_every`` steps."""
        self.r1_hist.append(float(r1))
        self.r2_hist.append(float(r2))
        self.cost_hist.append(float(cost))
        self.sigma_hist.append(float(sigma_cost))
        self.beta_hist.append(float(beta_c))
        if simplex is not None:
            self.current_simplex = simplex
            # Only append when the simplex actually changes — otherwise a
            # single NM iteration (reflect+expand or reflect+contract) pushes
            # the same triangle 2-3 times and the ghost trail looks stuck.
            if not self.simplex_hist or simplex != self.simplex_hist[-1]:
                self.simplex_hist.append(simplex)
        else:
            self.current_simplex = None
        self._step += 1
        # Only write a PNG on the first step and every save_every steps after.
        if self._step == 1 or self._step % self.save_every == 0:
            self._render(test_data, float(sigma_cost))

    # ------------------------------------------------------------------
    def _render(self, test_data: Optional[dict], sigma_cost: float = 0.0) -> None:
        # Snapshot all mutable state before handing off to background thread.
        r1 = np.asarray(self.r1_hist).copy()
        r2 = np.asarray(self.r2_hist).copy()
        c  = np.asarray(self.cost_hist).copy()
        s  = np.asarray(self.sigma_hist).copy()
        bc = np.asarray(self.beta_hist).copy()
        simplex_snap = list(self.simplex_hist)
        step_snap    = self._step
        method       = self.method
        out_path     = os.path.join(self.output_dir, f"optimizer_{self.method}.png")
        # For curve-collapse panel: copy everything needed.
        ref_data_snap = self.ref_data
        Lx, Ly, Tx, Ty = self.Lx, self.Ly, self.Tx, self.Ty
        n_copies = self.n_copies
        test_data_snap = test_data
        sigma_cost_snap = sigma_cost

        def _do_render():
            import matplotlib.pyplot as _plt
            import matplotlib.patches as _mpatches
            n = len(c)
            idx = np.arange(1, n + 1)
            fig, axes = _plt.subplots(2, 2, figsize=(11.0, 8.0))

            # ---- Panel 1: (r1, r2) trajectory + NM simplex ----------------
            ax = axes[0, 0]
            ax.plot(r1, r2, "-", color="lightgray", lw=0.8, zorder=1)
            c_for_color = np.where(c > 0, c, np.nanmin(c[c > 0]) if np.any(c > 0) else 1e-12)
            sc = ax.scatter(r1, r2, c=np.log10(c_for_color),
                            cmap="plasma", s=35, edgecolor="k", lw=0.3, zorder=2)
            best_i = int(np.argmin(c))
            ax.plot(r1[best_i], r2[best_i], "o", mfc="none", mec="lime",
                    ms=15, mew=2.0, label=f"best (eval {best_i + 1})", zorder=3)
            ax.plot(r1[-1], r2[-1], "*", color="cyan", ms=14,
                    mec="black", mew=0.6, label="current", zorder=4)
            n_ghost = 4
            simplices_to_draw = simplex_snap[-(n_ghost + 1):]
            n_draw = len(simplices_to_draw)
            for k, tri_verts in enumerate(simplices_to_draw):
                tri = np.array(tri_verts)
                is_current = (k == n_draw - 1)
                alpha_face = 0.25 if is_current else 0.04 + 0.06 * (k / max(n_draw - 1, 1))
                alpha_edge = 0.7 if is_current else 0.15
                lw = 1.4 if is_current else 0.6
                poly = _mpatches.Polygon(
                    tri, closed=True,
                    facecolor="deepskyblue", edgecolor="deepskyblue",
                    alpha=alpha_face, linewidth=lw, zorder=0 if not is_current else 5,
                )
                ax.add_patch(poly)
                if is_current:
                    ax.plot(np.append(tri[:, 0], tri[0, 0]),
                            np.append(tri[:, 1], tri[0, 1]),
                            "-", color="deepskyblue", lw=lw, alpha=alpha_edge, zorder=5)
            ax.set_xlabel("r₁"); ax.set_ylabel("r₂")
            ax.set_title(f"(r₁, r₂) trajectory — {method}, eval {n}")
            ax.legend(loc="best", fontsize=8); ax.grid(alpha=0.25)
            cbar = fig.colorbar(sc, ax=ax, shrink=0.85)
            cbar.set_label("log₁₀(cost)")

            # ---- Panel 2: cost vs evaluation ----------------------------
            ax = axes[0, 1]
            ax.errorbar(idx, c, yerr=s, fmt="o-", ms=4, lw=1.0, color="C0",
                        ecolor="C0", elinewidth=0.7, capsize=2, label="cost ±σ")
            ax.plot(idx, np.minimum.accumulate(c), "-", color="crimson", lw=1.5,
                    label="running best")
            ax.set_yscale("log"); ax.set_xlabel("evaluation #")
            ax.set_ylabel("L² cost"); ax.set_title("Cost history")
            ax.legend(loc="best", fontsize=8); ax.grid(alpha=0.25, which="both")

            # ---- Panel 3: beta_c history --------------------------------
            ax = axes[1, 0]
            ax.plot(idx, bc, "o-", color="C2", ms=4, lw=1.0)
            ax.set_xlabel("evaluation #"); ax.set_ylabel("β_c (this evaluation)")
            ax.set_title("β_c history (from on-the-fly GC scan)"); ax.grid(alpha=0.25)

            # ---- Panel 4: curve collapse --------------------------------
            ax = axes[1, 1]
            _render_curve_collapse(ax, test_data_snap, sigma_cost_snap,
                                   ref_data_snap, Lx, Ly, Tx, Ty, n_copies)

            fig.suptitle(
                f"{method}  |  current best: cost={c[best_i]:.4g}, "
                f"r₁={r1[best_i]:.4f}, r₂={r2[best_i]:.4f}, "
                f"β_c={bc[best_i]:.5f}", fontsize=11)
            fig.text(0.99, 0.99, f"step {step_snap}",
                     ha="right", va="top", fontsize=9, color="gray",
                     transform=fig.transFigure)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            fig.savefig(out_path, dpi=110)
            _plt.close(fig)

        _queue_render(_do_render)

    # ------------------------------------------------------------------
    def _plot_curve_collapse(self, ax, test_data, sigma_cost=0.0):
        _render_curve_collapse(ax, test_data, sigma_cost,
                               self.ref_data, self.Lx, self.Ly,
                               self.Tx, self.Ty, self.n_copies)


def _render_curve_collapse(ax, test_data, sigma_cost,
                           ref_data, Lx, Ly, Tx, Ty, n_copies):
        if test_data is None:
            ax.text(0.5, 0.5, "(no test data this step)",
                    ha="center", va="center", transform=ax.transAxes,
                    color="gray")
            ax.set_axis_off()
            return
        try:
            dirs = boundary_paths(Lx, Ly, Tx, Ty)
        except Exception as exc:
            ax.text(0.5, 0.5, f"(boundary_paths error: {exc})",
                    ha="center", va="center", transform=ax.transAxes,
                    color="red")
            ax.set_axis_off()
            return

        colors = {"v": "C0", "u": "C1", "w": "C2"}
        try:
            g_ref = _tile_interp(ref_data, Lx, Ly,
                                 Tx, Ty, "conn",
                                 copies=n_copies)
            g_test = _tile_interp(test_data, Lx, Ly,
                                  Tx, Ty, "conn",
                                  copies=n_copies)
        except Exception as exc:
            ax.text(0.5, 0.5, f"(interp build error: {exc})",
                    ha="center", va="center", transform=ax.transAxes,
                    color="red")
            ax.set_axis_off()
            return

        names = ["v", "u", "w"]
        n_samples = 64
        t = np.linspace(0.0, 1.0, n_samples)
        ax.axhline(0.0, color="black", lw=0.9, ls="--", alpha=0.5,
                   label="zero (perfect match)")
        for name, (dm, dn) in zip(names, dirs):
            try:
                xs, ys = _triangular_xy(t * dm, t * dn)
                ref_vals  = np.asarray(g_ref(np.column_stack([xs, ys])),  float)
                test_vals = np.asarray(g_test(np.column_stack([xs, ys])), float)
            except Exception as exc:
                ax.text(0.5, 0.5, f"(interp error: {exc})",
                        ha="center", va="center", transform=ax.transAxes,
                        color="red")
                ax.set_axis_off()
                return
            diff = test_vals - ref_vals
            ax.plot(t, diff, "-", color=colors[name], lw=1.3,
                    label=f"Δ{name}")
            # Ribbon: ±per-point noise estimate following the curve.
            # Derived from sigma_cost: if sigma_cost is the MC error on the
            # total integrated L² cost, a rough per-point noise scale per
            # direction is sqrt(sigma_cost / (3 * n_samples)).
            if sigma_cost > 0:
                noise = float(np.sqrt(max(sigma_cost, 0.0) / (3 * n_samples)))
            else:
                # Fallback: use RMS of the residual itself.
                noise = float(np.sqrt(np.nanmean(diff**2)))
            ax.fill_between(t, diff - noise, diff + noise,
                            color=colors[name], alpha=0.18)

        ax.set_xlabel("path parameter  t")
        ax.set_ylabel("G_test − G_ref")
        ax.set_title("Boundary residuals (0 = perfect match)")
        ax.legend(loc="best", fontsize=7, ncol=2)
        ax.grid(alpha=0.25)


# ---------------------------------------------------------------------------
#  Convenience: assemble frames into an animated GIF
# ---------------------------------------------------------------------------

def write_gif_command_hint(frames_dir: str) -> str:
    """Return a shell command (string) the user can copy to build a GIF."""
    return (
        f"# Build GIFs from frame sequences:\n"
        f"cd {frames_dir} && \\\n"
        f"ffmpeg -y -framerate 4 -pattern_type glob "
        f"-i 'optimizer_*_step*.png' -vf scale=900:-1 optimizer.gif && \\\n"
        f"ffmpeg -y -framerate 2 -pattern_type glob "
        f"-i 'betac_scan_*.png' -vf scale=800:-1 betac_scan.gif"
    )
