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
from typing import Optional, Sequence

import matplotlib

# Pick a backend.  If we're inside an IPython kernel (Spyder, Jupyter,
# qtconsole) leave whatever inline / Qt backend is already configured so
# `plt.show()` actually pops figures into the Plots pane.  In a plain
# headless interpreter, fall back to the non-GUI Agg backend.
def _select_backend():
    if os.environ.get("KOPT_FORCE_AGG"):
        matplotlib.use("Agg", force=True)
        return
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            return  # respect the kernel's inline / Qt5Agg / etc.
    except Exception:
        pass
    # Plain Python: use Agg so we don't need a display.
    matplotlib.use("Agg", force=True)


_select_backend()
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402


def _is_interactive_backend() -> bool:
    return matplotlib.get_backend().lower() not in ("agg", "pdf", "svg", "ps")


def _show_or_close(fig) -> None:
    """Display the figure (interactive backend) or close it (Agg)."""
    if _is_interactive_backend():
        try:
            plt.show()
        except Exception:
            pass
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Render queue (kept for backward compatibility — now a synchronous no-op).
# ---------------------------------------------------------------------------

def _queue_render(fn) -> None:
    """Run a zero-arg render callable immediately on the main thread.

    Rendering must happen on the main thread when an interactive backend
    (Qt, inline, …) is in use; doing it synchronously also makes plots
    appear in Spyder/Jupyter in the right order.
    """
    try:
        fn()
    except Exception as exc:
        print(f"[vis] render error: {exc}")


def flush_render_queue() -> None:
    """No-op: rendering is now synchronous."""
    return


# Back-compat: some callers may still call display_inline(path).  When the
# matplotlib backend is interactive we don't need to do anything (the figure
# is already shown via plt.show); keep a stub so old code doesn't break.
def display_inline(png_path: str, *, replace: bool = True) -> None:  # noqa: D401
    return

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
            _show_or_close(fig)

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
        Lattice dimensions of the **test** run (for the test boundary
        paths).
    ref_Lx, ref_Ly, ref_Tx, ref_Ty : int, optional
        Lattice dimensions of the **reference** run.  When omitted they
        default to the test geometry (back-compat for callers that don't
        yet pass them).  When ref_geom != test_geom, each correlator is
        sampled on its own torus's boundary path (matching the cost-
        function convention in cost.py) before computing the residual.
    n_copies : int, default 1
        Number of periodic tile copies to use in interpolation.
    """

    def __init__(self, output_dir: str, method: str,
                 ref_data: dict, Lx: int, Ly: int, Tx: int, Ty: int,
                 n_copies: int = 1, save_every: int = 5,
                 ref_Lx: Optional[int] = None, ref_Ly: Optional[int] = None,
                 ref_Tx: Optional[int] = None, ref_Ty: Optional[int] = None):
        self.output_dir = output_dir
        self.method = method
        self.ref_data = ref_data
        # Test geometry (used for sampling the test correlator).
        self.Lx = Lx
        self.Ly = Ly
        self.Tx = Tx
        self.Ty = Ty
        # Reference geometry (used for sampling the reference correlator).
        # Defaults to the test geometry when not supplied.
        self.ref_Lx = Lx if ref_Lx is None else int(ref_Lx)
        self.ref_Ly = Ly if ref_Ly is None else int(ref_Ly)
        self.ref_Tx = Tx if ref_Tx is None else int(ref_Tx)
        self.ref_Ty = Ty if ref_Ty is None else int(ref_Ty)
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
        # Gaussian history: each entry is a dict {mean, cov, sigma, gen}
        # written by run_cmaes via evaluator.current_gaussian, one per
        # CMA-ES generation.
        self.gaussian_hist: list = []
        self.current_gaussian: Optional[dict] = None
        # GP surface snapshot (BO mode): {r1, r2, mean, std, acq, next, step}
        self.current_gp_surface: Optional[dict] = None

        # Latest β_c-scan state, captured via update_scan().  Used to draw
        # the 5th panel of the optimizer PNG.
        self._scan_snap: Optional[dict] = None

        # ---- Speedup 3: out-of-order eval buffering ----
        # When the parallel CMA-ES path delivers EvalResults in arrival
        # order rather than eval_id order, callers pass eval_id=N to
        # update(); we buffer everything until the next contiguous id is
        # available and only then flush in strict eval_id sequence.  This
        # keeps the recorded PNG history monotone even when generations
        # finish on multiple workers concurrently.
        #
        # If callers never supply eval_id (legacy / serial path), the
        # buffering machinery is bypassed and update() behaves as before.
        self._pending_updates: dict = {}     # eval_id -> kwargs
        self._next_eval_id: int = 1

    def update_scan_from_result(self, scan_betas, scan_chis, scan_chi_errs,
                                beta_c: float) -> None:
        """Populate the β_c scan panel from EvalResult fields.

        Called by the parallel forwarder since workers set optimizer_plot=None
        and never fire the progress_cb; the scan arrays ARE stored in
        EvalResult and can be used to reconstruct the panel.
        A GC fit is attempted on the reweighted grid so the peak curve is drawn.
        """
        if not scan_betas:
            return
        n = len(scan_betas)
        errs = list(scan_chi_errs) if scan_chi_errs else [0.0] * n

        # Fit GC to the reweighted grid so the panel shows the fit curve.
        gc_params = None
        try:
            from mc_engine import _gc_fit as _gcf
            valid = [i for i, c in enumerate(scan_chis) if np.isfinite(c)]
            if len(valid) >= 6:
                vb = [scan_betas[i] for i in valid]
                vc = [scan_chis[i]  for i in valid]
                gc_params, _ = _gcf(vb, vc, float(beta_c))
        except Exception:
            pass

        self._scan_snap = {
            "betas": list(scan_betas),
            "chis":  list(scan_chis),
            "errs":  errs,
            "pids":  [0] * n,
            "gc":    gc_params,
            "b_est": float(beta_c),
            "pass":  0,
        }

    # ------------------------------------------------------------------
    def update_scan(self, state: dict) -> None:
        """Hook usable as ``progress_cb`` for ``mc_engine.find_beta_c``.

        Stores a deep-copied snapshot of the current scan so the next
        call to ``update()`` can render it as the optimizer PNG's
        β_c-scan panel.
        """
        try:
            self._scan_snap = {
                "betas": list(state.get("all_betas", [])),
                "chis":  list(state.get("all_chis", [])),
                "errs":  list(state.get("all_chi_errs", [])),
                "pids":  list(state.get("pass_ids", [])),
                "gc":    state.get("gc_params"),
                "b_est": state.get("beta_estimate"),
                "pass":  int(state.get("pass_num", 0)),
            }
        except Exception as exc:
            print(f"[vis] update_scan error: {exc}")

    # ------------------------------------------------------------------
    def update(self, r1: float, r2: float, cost: float, sigma_cost: float,
               beta_c: float, test_data: Optional[dict] = None,
               simplex: Optional[list] = None,
               gaussian: Optional[dict] = None,
               gp_surface: Optional[dict] = None,
               eval_id: Optional[int] = None) -> None:
        """Record one evaluation; write a frame every ``save_every`` steps.

        When ``eval_id`` is given (parallel-CMA path), updates are
        buffered and only applied in strict eval_id order so the saved
        PNG history stays monotone even when workers finish out of
        sequence.  Legacy serial callers omit ``eval_id`` and get
        immediate, in-order updates.
        """
        if eval_id is not None:
            self._pending_updates[int(eval_id)] = {
                "r1": float(r1), "r2": float(r2),
                "cost": float(cost), "sigma_cost": float(sigma_cost),
                "beta_c": float(beta_c),
                "test_data": test_data,
                "simplex": simplex, "gaussian": gaussian,
                "gp_surface": gp_surface,
            }
            while self._next_eval_id in self._pending_updates:
                kw = self._pending_updates.pop(self._next_eval_id)
                self._next_eval_id += 1
                self._apply_update(**kw)
            return
        # Legacy serial path: apply immediately.  Keep _next_eval_id
        # advancing too, so a later switch to the buffered path stays
        # consistent.
        self._next_eval_id += 1
        self._apply_update(r1, r2, cost, sigma_cost, beta_c,
                           test_data, simplex, gaussian, gp_surface)

    def flush_render_queue(self) -> None:
        """Force any out-of-order pending updates to render.

        Call at end-of-run after the optimizer has terminated to make
        sure no straggler eval got dropped.  In legacy / serial mode this
        is a no-op (``_pending_updates`` is always empty).
        """
        for eid in sorted(self._pending_updates):
            kw = self._pending_updates.pop(eid)
            if eid == self._next_eval_id:
                self._next_eval_id += 1
            self._apply_update(**kw)

    def _apply_update(self, r1, r2, cost, sigma_cost, beta_c,
                      test_data=None, simplex=None, gaussian=None,
                      gp_surface=None):
        """Internal: append one evaluation to history and maybe render."""
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
        if gaussian is not None:
            self.current_gaussian = gaussian
            # One entry per CMA-ES generation: dedupe on the "gen" field if
            # present, otherwise on object identity.
            new_gen = gaussian.get("gen") if isinstance(gaussian, dict) else None
            if (not self.gaussian_hist
                    or self.gaussian_hist[-1].get("gen") != new_gen):
                self.gaussian_hist.append(dict(gaussian))
        else:
            self.current_gaussian = None
        if gp_surface is not None:
            self.current_gp_surface = gp_surface
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
        gaussian_snap = [dict(g) for g in self.gaussian_hist]
        gp_surface_snap = dict(self.current_gp_surface) if self.current_gp_surface is not None else None
        step_snap    = self._step
        method       = self.method
        out_path     = os.path.join(self.output_dir, f"optimizer_{self.method}.png")
        # For curve-collapse panel: copy everything needed.
        ref_data_snap = self.ref_data
        Lx, Ly, Tx, Ty = self.Lx, self.Ly, self.Tx, self.Ty
        ref_geom = (self.ref_Lx, self.ref_Ly, self.ref_Tx, self.ref_Ty)
        n_copies = self.n_copies
        test_data_snap = test_data
        sigma_cost_snap = sigma_cost
        scan_snap = dict(self._scan_snap) if self._scan_snap else None

        def _do_render():
            import matplotlib.pyplot as _plt
            import matplotlib.patches as _mpatches
            n = len(c)
            idx = np.arange(1, n + 1)
            fig = _plt.figure(figsize=(15.5, 8.5))
            gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.2])
            axes = np.array([
                [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])],
                [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])],
            ])
            ax_scan = fig.add_subplot(gs[:, 2])

            # ---- Panel 1: (r1, r2) trajectory + NM simplex ----------------
            ax = axes[0, 0]
            # GP surface background (BO mode): draw mean contourf first.
            if gp_surface_snap is not None:
                try:
                    _R1 = np.asarray(gp_surface_snap["r1"])
                    _R2 = np.asarray(gp_surface_snap["r2"])
                    _mu = np.asarray(gp_surface_snap["mean"])
                    _safe = np.clip(_mu, 1e-12, None)
                    ax.contourf(_R1, _R2, np.log10(_safe),
                                levels=20, cmap="plasma", alpha=0.45, zorder=0)
                    # Acquisition function (LCB) contour lines.
                    _acq = np.asarray(gp_surface_snap["acq"])
                    ax.contour(_R1, _R2, _acq, levels=8,
                               colors="white", linewidths=0.4,
                               linestyles="--", alpha=0.35, zorder=1)
                    # Mark the proposed next point.
                    _nxt = gp_surface_snap.get("next")
                    if _nxt is not None:
                        ax.plot(_nxt[0], _nxt[1], "D", color="white",
                                ms=9, mec="black", mew=0.8,
                                label=f"BO proposal (step {gp_surface_snap.get('step','?')})",
                                zorder=6)
                except Exception:
                    pass
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
            # CMA-ES Gaussian ghost trail: 1σ + 2σ confidence ellipses.
            gaussians_to_draw = gaussian_snap[-(n_ghost + 1):]
            n_g = len(gaussians_to_draw)
            _gauss_title = ""
            for k, g in enumerate(gaussians_to_draw):
                try:
                    mean = np.asarray(g["mean"], dtype=float).reshape(2)
                    cov  = np.asarray(g["cov"], dtype=float).reshape(2, 2)
                    sigma = float(g.get("sigma", 1.0))
                except Exception:
                    continue
                # Eigendecomposition of σ²·C gives ellipse axes.
                C = (sigma ** 2) * cov
                try:
                    evals, evecs = np.linalg.eigh(C)
                except np.linalg.LinAlgError:
                    continue
                evals = np.clip(evals, 0.0, None)
                order = np.argsort(evals)[::-1]
                evals = evals[order]; evecs = evecs[:, order]
                angle = float(np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0])))
                width_1s  = 2.0 * float(np.sqrt(evals[0]))
                height_1s = 2.0 * float(np.sqrt(evals[1]))
                is_current = (k == n_g - 1)
                alpha_face = 0.22 if is_current else 0.04 + 0.05 * (k / max(n_g - 1, 1))
                alpha_edge = 0.85 if is_current else 0.20
                lw = 1.4 if is_current else 0.6
                z_face = 5 if is_current else 0
                # 1σ ellipse (filled).
                ell1 = _mpatches.Ellipse(
                    xy=mean, width=width_1s, height=height_1s, angle=angle,
                    facecolor="orange", edgecolor="darkorange",
                    alpha=alpha_face, linewidth=lw, zorder=z_face,
                )
                ax.add_patch(ell1)
                # 2σ ellipse (outline only) for the current generation.
                if is_current:
                    ell2 = _mpatches.Ellipse(
                        xy=mean, width=2.0 * width_1s, height=2.0 * height_1s,
                        angle=angle, facecolor="none", edgecolor="darkorange",
                        alpha=alpha_edge * 0.7, linewidth=lw * 0.8,
                        linestyle="--", zorder=z_face,
                    )
                    ax.add_patch(ell2)
                    ax.plot(mean[0], mean[1], "x", color="darkorange",
                            ms=8, mew=1.6, zorder=z_face + 1,
                            label=f"CMA mean (gen {g.get('gen', '?')})")
                    # Collect params for the panel title.
                    _ax1_r = float(np.sqrt(evals[0]))
                    _ax2_r = float(np.sqrt(evals[1]))
                    _cond  = _ax1_r / _ax2_r if _ax2_r > 0 else float("inf")
                    _gauss_title = (
                        f"mean ({mean[0]:+.4f}, {mean[1]:+.4f})  "
                        f"σ={sigma:.4f}  axes={_ax1_r:.4f}/{_ax2_r:.4f}  "
                        f"θ={angle:.1f}°  cond={_cond:.2f}"
                    )
            ax.set_xlabel("r₁"); ax.set_ylabel("r₂")
            _best_line = (f"best eval {best_i+1}: r₁={r1[best_i]:+.4f}  r₂={r2[best_i]:+.4f}  "
                          f"cost={c[best_i]:.3e}")
            _title = f"(r₁, r₂) trajectory — {method}, eval {n}\n{_best_line}"
            if _gauss_title:
                _title += f"\n{_gauss_title}"
            ax.set_title(_title, fontsize=8)
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
                                   ref_data_snap, Lx, Ly, Tx, Ty, n_copies,
                                   ref_geom=ref_geom)

            # ---- Panel 5: β_c scan (latest pass) -----------------------
            _render_betac_scan_panel(ax_scan, scan_snap)

            fig.suptitle(
                f"{method}  |  current best: cost={c[best_i]:.4g}, "
                f"r₁={r1[best_i]:.4f}, r₂={r2[best_i]:.4f}, "
                f"β_c={bc[best_i]:.5f}", fontsize=11)
            fig.text(0.99, 0.99, f"step {step_snap}",
                     ha="right", va="top", fontsize=9, color="gray",
                     transform=fig.transFigure)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            fig.savefig(out_path, dpi=110)
            _show_or_close(fig)

        _queue_render(_do_render)

    # ------------------------------------------------------------------
    def _plot_curve_collapse(self, ax, test_data, sigma_cost=0.0):
        _render_curve_collapse(ax, test_data, sigma_cost,
                               self.ref_data, self.Lx, self.Ly,
                               self.Tx, self.Ty, self.n_copies,
                               ref_geom=(self.ref_Lx, self.ref_Ly,
                                         self.ref_Tx, self.ref_Ty))


def _render_betac_scan_panel(ax, scan_snap):
    """Render the latest β_c scan onto a single axis (5th panel of opt PNG)."""
    if not scan_snap or not scan_snap.get("betas"):
        ax.text(0.5, 0.5, "(no β_c scan yet)", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        ax.set_axis_off()
        return
    betas = np.asarray(scan_snap["betas"], dtype=float)
    chis  = np.asarray(scan_snap["chis"],  dtype=float)
    errs  = np.asarray(scan_snap["errs"],  dtype=float)
    pids  = np.asarray(scan_snap["pids"],  dtype=int)
    gc    = scan_snap.get("gc")
    b_est = scan_snap.get("b_est")
    pass_num = int(scan_snap.get("pass", 0))

    cmap = plt.get_cmap("viridis")
    n_passes_seen = max(int(pids.max()) + 1, 1) if len(pids) else 1
    for p in range(n_passes_seen):
        mask = pids == p
        if not mask.any():
            continue
        col = cmap(p / max(n_passes_seen - 1, 1))
        ax.errorbar(betas[mask], chis[mask], yerr=errs[mask],
                    fmt="o", ms=5, lw=0,
                    ecolor=col, color=col,
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
    ax.set_title(f"β_c scan (pass {pass_num})")
    ax.legend(loc="best", fontsize=7, framealpha=0.9)
    ax.grid(alpha=0.25)


def _render_curve_collapse(ax, test_data, sigma_cost,
                           ref_data, Lx, Ly, Tx, Ty, n_copies,
                           ref_geom=None):
        """Plot the boundary residuals  G_test(t) − G_ref(t)  along the
        three torus boundary directions.

        Each correlator is sampled on its OWN torus's boundary path
        (parameterised by t ∈ [0, 1)), then differenced at matching t.
        This matches the cost-function convention in cost.py and is
        required when ref_geom != test_geom (e.g. ref 26×32 twisted vs
        test 24×24 untwisted).  The curve is closed by appending the
        t=0 sample at t=1 since the parameterisation is periodic.
        """
        if test_data is None:
            ax.text(0.5, 0.5, "(no test data this step)",
                    ha="center", va="center", transform=ax.transAxes,
                    color="gray")
            ax.set_axis_off()
            return
        if ref_geom is None:
            ref_geom = (Lx, Ly, Tx, Ty)
        rLx, rLy, rTx, rTy = ref_geom
        try:
            test_dirs = boundary_paths(Lx, Ly, Tx, Ty)
            ref_dirs  = boundary_paths(rLx, rLy, rTx, rTy)
        except Exception as exc:
            ax.text(0.5, 0.5, f"(boundary_paths error: {exc})",
                    ha="center", va="center", transform=ax.transAxes,
                    color="red")
            ax.set_axis_off()
            return

        colors = {"v": "C0", "u": "C1", "w": "C2"}
        try:
            g_ref = _tile_interp(ref_data, rLx, rLy,
                                 rTx, rTy, "conn",
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
        # Periodic sampling: t ∈ [0, 1) avoids the LinearNDInterpolator
        # convex-hull spike at t=1 (see cost.py).  We close the loop for
        # plotting by appending the t=0 value at t=1.
        t = np.linspace(0.0, 1.0, n_samples, endpoint=False)
        t_plot = np.concatenate([t, [1.0]])
        ax.axhline(0.0, color="black", lw=0.9, ls="--", alpha=0.5,
                   label="zero (perfect match)")
        for name, (rdm, rdn), (tdm, tdn) in zip(names, ref_dirs, test_dirs):
            try:
                rxs, rys = _triangular_xy(t * rdm, t * rdn)
                txs, tys = _triangular_xy(t * tdm, t * tdn)
                ref_vals  = np.asarray(g_ref(np.column_stack([rxs, rys])),  float)
                test_vals = np.asarray(g_test(np.column_stack([txs, tys])), float)
            except Exception as exc:
                ax.text(0.5, 0.5, f"(interp error: {exc})",
                        ha="center", va="center", transform=ax.transAxes,
                        color="red")
                ax.set_axis_off()
                return
            diff = test_vals - ref_vals
            diff_plot = np.concatenate([diff, diff[:1]])
            ax.plot(t_plot, diff_plot, "-", color=colors[name], lw=1.3,
                    label=f"Δ{name}")
            # Ribbon: ±per-point noise estimate following the curve.
            # Derived from sigma_cost: if sigma_cost is the MC error on the
            # total integrated L² cost, a rough per-point noise scale per
            # direction is sqrt(sigma_cost / (3 * n_samples)).
            if sigma_cost > 0:
                noise = float(np.sqrt(max(sigma_cost, 0.0) / (3 * n_samples)))
            else:
                # Fallback: use RMS of the residual itself.
                noise = float(np.sqrt(np.nanmean(diff_plot**2)))
            ax.fill_between(t_plot, diff_plot - noise, diff_plot + noise,
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
