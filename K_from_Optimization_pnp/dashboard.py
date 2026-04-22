"""
dashboard.py — Live terminal dashboard for the K_from_Optimization pipeline.

Renders a single, continuously-updated view of:

  Top:    static configuration (reference geometry, β_c, test geometry).
  Middle: current β_c scan (pass index, current bracket, points done,
          best estimate so far, most-recent χ values).
  Bottom: optimizer history table (one row per evaluation) plus a
          running "best so far" line.

Uses `rich`.  All `print()` output from mc_engine/optimizer is preserved
above the dashboard (Live handles scroll-back transparently).

Usage (from run.py)::

    from dashboard import Dashboard
    with Dashboard(ref_meta, test_geom, max_evals) as dash:
        ev = Evaluator(..., dashboard=dash)
        run_nelder_mead(ev, ...)
"""
from __future__ import annotations

import time
from typing import Optional

from rich.console import Console
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bar(done: int, total: int, width: int = 24) -> str:
    if total <= 0:
        return "·" * width
    n = max(0, min(width, int(round(width * done / total))))
    return "█" * n + "░" * (width - n)


def _fmt(x, fmt: str = "{:+.4f}") -> str:
    if x is None:
        return "—"
    try:
        return fmt.format(float(x))
    except Exception:
        return str(x)


# ---------------------------------------------------------------------------
# Dashboard
# ---------------------------------------------------------------------------

class Dashboard:
    """Live terminal dashboard.  Use as a context manager."""

    def __init__(self, ref_meta: dict, test_geom: tuple, max_evals: int,
                 method_name: str = "nelder_mead",
                 max_history_rows: int = 12):
        self.ref_meta = ref_meta
        self.test_geom = tuple(test_geom)
        self.max_evals = int(max_evals)
        self.method_name = method_name
        self.max_history_rows = int(max_history_rows)

        # State
        self._t0 = time.time()
        self._history: list = []           # list of EvalResult
        self._scan_state: Optional[dict] = None
        self._scan_eval_id: Optional[int] = None
        self._scan_r1: Optional[float] = None
        self._scan_r2: Optional[float] = None
        self._scan_bracket: Optional[tuple] = None

        self._console = Console()
        self._layout: Optional[Layout] = None
        self._live: Optional[Live] = None

    # -------- context manager --------

    def __enter__(self):
        self._layout = Layout()
        self._layout.split(
            Layout(name="header", size=4),
            Layout(name="scan",   size=10),
            Layout(name="opt"),
        )
        self._refresh()
        self._live = Live(
            self._layout,
            console=self._console,
            refresh_per_second=4,
            screen=False,
            redirect_stdout=False,
            redirect_stderr=False,
        )
        self._live.__enter__()
        return self

    def __exit__(self, *exc):
        if self._live is not None:
            self._live.__exit__(*exc)
        return False

    # -------- public API --------

    def begin_eval(self, eval_id: int, r1: float, r2: float,
                   bracket: tuple) -> None:
        self._scan_eval_id = int(eval_id)
        self._scan_r1, self._scan_r2 = float(r1), float(r2)
        self._scan_bracket = (float(bracket[0]), float(bracket[1]))
        self._scan_state = {
            "pass_num": 0,
            "all_betas": [], "all_chis": [], "all_chi_errs": [],
            "pass_ids": [], "gc_params": None, "beta_estimate": None,
            "scan_progress": None, "traj_done": 0,
        }
        self._refresh()

    def update_scan(self, info: dict) -> None:
        """Hook used as `progress_cb` for `mc_engine.find_beta_c`."""
        self._scan_state = info
        self._refresh()

    def update_eval(self, result) -> None:
        self._history.append(result)
        # Clear the scan once the evaluation is finished.
        self._scan_state = None
        self._scan_eval_id = None
        self._refresh()

    # -------- rendering --------

    def _refresh(self) -> None:
        self._layout["header"].update(self._render_header())
        self._layout["scan"].update(self._render_scan())
        self._layout["opt"].update(self._render_opt())

    def _render_header(self) -> Panel:
        meta = self.ref_meta
        ref_geom = tuple(meta.get("geometry", (meta.get("L"),) * 2 + (0, 0)))
        Lx, Ly, Tx, Ty = self.test_geom
        wall = time.time() - self._t0
        text = Text()
        text.append("ref:  ", style="bold cyan")
        text.append(f"L={meta.get('L', '?')}×{meta.get('L', '?')}  "
                    f"β_c={meta.get('beta_c', float('nan')):.6f}  "
                    f"n_traj={meta.get('n_traj', '?')}\n")
        text.append("test: ", style="bold cyan")
        text.append(f"{Lx}×{Ly}  Tx={Tx} Ty={Ty}    "
                    f"method={self.method_name}    "
                    f"max_evals={self.max_evals}    "
                    f"wall={wall:6.1f}s")
        return Panel(text, title="K_from_Optimization", border_style="cyan")

    def _render_scan(self) -> Panel:
        s = self._scan_state
        if s is None:
            return Panel(Text("(idle — between evaluations)", style="dim"),
                         title="β_c scan", border_style="grey50")

        eid = self._scan_eval_id
        r1, r2 = self._scan_r1, self._scan_r2
        lo, hi = self._scan_bracket if self._scan_bracket else (None, None)
        pass_num = s.get("pass_num", 0)
        beta_est = s.get("beta_estimate")
        n_pts = len(s.get("all_betas", []))
        traj_done = s.get("traj_done", 0)

        prog = s.get("scan_progress") or {}
        pts_done = prog.get("pts_done", 0)
        pts_total = prog.get("pts_total", 0)

        # Best-so-far estimate from raw χ
        chis = s.get("all_chis", [])
        betas = s.get("all_betas", [])
        if chis:
            i = max(range(len(chis)), key=lambda j: chis[j])
            beta_argmax = betas[i]
            chi_max = chis[i]
        else:
            beta_argmax = None
            chi_max = None

        text = Text()
        text.append(f"eval {eid:03d}  ", style="bold")
        text.append(f"r1={r1:+.4f} r2={r2:+.4f}    "
                    f"bracket=[{_fmt(lo)}, {_fmt(hi)}]\n")

        text.append(f"pass {pass_num}", style="bold yellow")
        text.append(f"   pts={n_pts}   traj={traj_done}\n")

        if pts_total:
            text.append(f"  this pass: {_bar(pts_done, pts_total)} "
                        f"{pts_done}/{pts_total}\n")

        text.append(f"  current β_c estimate : ")
        text.append(_fmt(beta_est, "{:.6f}"), style="bold green")
        text.append("\n  argmax χ              : ")
        text.append(_fmt(beta_argmax, "{:.6f}"))
        text.append("    χ_max = ")
        text.append(_fmt(chi_max, "{:.4g}"))
        return Panel(text, title=f"β_c scan", border_style="yellow")

    def _render_opt(self) -> Panel:
        if not self._history:
            return Panel(Text("(no evaluations yet)", style="dim"),
                         title="optimizer", border_style="green")

        tbl = Table(expand=True, header_style="bold magenta",
                    show_lines=False, padding=(0, 1))
        tbl.add_column("ev", justify="right", style="cyan", width=5)
        tbl.add_column("r1", justify="right", width=7)
        tbl.add_column("r2", justify="right", width=7)
        tbl.add_column("β_c", justify="right", width=9)
        tbl.add_column("cost", justify="right", width=10)
        tbl.add_column("σ_cost", justify="right", width=9)
        tbl.add_column("SNR", justify="right", width=6)
        tbl.add_column("wall(s)", justify="right", width=8)

        best = min(self._history, key=lambda r: r.cost)

        rows = self._history[-self.max_history_rows:]
        for r in rows:
            is_best = (r.eval_id == best.eval_id)
            star = "★" if is_best else " "
            style = "bold green" if is_best else None
            # Color SNR by status
            snr_style = ("green" if r.snr_status == "ok"
                         else "yellow" if r.snr_status == "marginal"
                         else "red")
            tbl.add_row(
                f"{star}{r.eval_id:03d}",
                f"{r.r1:+.4f}",
                f"{r.r2:+.4f}",
                f"{r.beta_c:.6f}",
                f"{r.cost:.3e}",
                f"{r.sigma_cost:.2e}",
                Text(f"{r.snr:.2f}", style=snr_style),
                f"{r.wall_time_s:.1f}",
                style=style,
            )

        # Footer / summary
        total_traj = sum(r.n_traj_prod + r.n_traj_scan_total
                         for r in self._history)
        wall_total = sum(r.wall_time_s for r in self._history)
        footer = Text()
        footer.append(f"\nbest: ", style="bold green")
        footer.append(f"r1={best.r1:+.4f}  r2={best.r2:+.4f}  "
                      f"β_c={best.beta_c:.6f}  cost={best.cost:.3e}  "
                      f"SNR={best.snr:.2f}  (eval {best.eval_id:03d})\n")
        footer.append(f"runs: {len(self._history)}/{self.max_evals}    "
                      f"total_traj={total_traj}    "
                      f"wall_in_evals={wall_total:.1f}s",
                      style="dim")

        from rich.console import Group
        return Panel(Group(tbl, footer),
                     title=f"optimizer ({self.method_name})",
                     border_style="green")
