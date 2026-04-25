"""run_panel.py — Driver for the post-speedup validation panel.

Implements §A–§G of the README "Validation test panel" section.
This is the orchestration layer that runs CMA-ES (and NM as a smoke
test) across the (geometry × start × seed × speedup) grid and writes
``PANEL_<speedup>.md`` summaries plus ``PANEL_LATEST.md``.

Usage
-----
    # Run the full panel for one speedup:
    python tests/validation/run_panel.py --speedup s4

    # Run a quick subset (1 geometry × 2 starts × 1 seed) for a smoke check:
    python tests/validation/run_panel.py --speedup s4 \
        --geoms R-equi --starts near,far_lo --seeds 1

    # Dry run — print the matrix it WOULD execute, no MC:
    python tests/validation/run_panel.py --speedup s4 --dry-run

    # Render the rolled-up Markdown table from existing per-cell logs:
    python tests/validation/run_panel.py --report

Speedup IDs
-----------
    s2     = β_c interpolation cache (cold + warm sub-panels)
    s3     = parallel CMA-ES (n_workers ∈ {1, 4, λ})
    s4     = Ferrenberg–Swendsen reweighting
    s5     = BO surrogate
    legacy = baseline serial / no speedups (used to populate
             tests/baselines/ before any speedup is panelled)
    all    = run every speedup sequentially

Outputs
-------
    results/_validation/<speedup>/<tag>/<start>/seed_<n>/
        eval_log.jsonl  — per-eval cost log (handed off to existing infra)
        summary.json    — run summary including final (r1*, r2*)
    results/_validation/<speedup>/PANEL.md    — per-speedup summary
    results/_validation/PANEL_LATEST.md       — combined latest table

This driver is intentionally subprocess-based: it shells out to
``python run.py`` for each cell so each test is fully isolated (no
shared state across cells).  The trade-off is wall time; the panel is
parallelisable via ``--n-workers`` (which forks the cell loop, NOT
the CMA-ES population — that is managed by the s3 speedup itself).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

HERE = Path(__file__).resolve().parent
PKG = HERE.parent.parent  # K_from_Optimization_pnp/
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

from tests.validation.geometries import (  # noqa: E402
    CMAES_START_DISTRIBUTIONS,
    NM_CONTROL_STARTS, NM_CONTROL_MAX_EVALS,
    PANEL_CMAES_MAX_EVALS, PANEL_SEEDS_PER_CELL,
    get_geometry, list_geometries,
)


# ---------------------------------------------------------------------------
# Config templates — speedup-specific CONFIG overrides
# ---------------------------------------------------------------------------

def _base_cfg(geom: dict, start_x0, start_sigma, seed: int,
              optimizer: str, max_evals: int) -> dict:
    """Build a CONFIG dict for one cell.

    Mirrors the keys read by run.py.  Values reflect the panel's
    "small-but-honest" budget: ~10–15 s per eval at L=12, ~50 evals
    per CMA-ES cell ⇒ ~10 min per cell.  Multiply by 5 geoms × 5
    starts × 3 seeds = 75 cells ⇒ ~10 wall hours, parallelisable.
    """
    return {
        # Geometry
        "ref_Lx": geom["ref_Lx"], "ref_Ly": geom["ref_Ly"],
        "ref_Tx": geom["ref_Tx"], "ref_Ty": geom["ref_Ty"],
        "test_Lx": geom["test_Lx"], "test_Ly": geom["test_Ly"],
        "test_Tx": geom["test_Tx"], "test_Ty": geom["test_Ty"],
        # Reference build
        "ref_n_traj":             40_000,
        "ref_scan_n_traj_coarse":  4_000,
        "ref_scan_n_traj_fine":   10_000,
        "ref_beta_seed":          (0.20, 0.32),
        "ref_scan_n_coarse":      11,
        "ref_scan_n_refine":      5,
        "ref_scan_n_refine2":     5,
        "ref_scan_n_refine3":     5,
        "ref_scan_max_shifts":    4,
        # Per-eval MC stats
        "n_traj_prod":            10_000,
        "n_traj_scan_coarse":      2_000,
        "n_traj_scan_fine":        4_000,
        "scan_n_coarse":          11,
        "scan_n_refine":          5,
        "scan_n_refine2":         5,
        "scan_n_refine3":         5,
        "scan_max_shifts":        4,
        "scan_jackknife":         False,
        "beta_seed":              (0.20, 0.32),
        # Optimizer
        "optimizer":              optimizer,
        "x0":                     tuple(start_x0),
        "max_evals":              max_evals,
        # NM knobs (used iff optimizer == "nelder_mead")
        "nm_sigma0":              float(start_sigma),
        "nm_xatol":               0.005,
        "nm_fatol":               1e-6,
        "nm_shrink":              0.75,
        # CMA knobs (used iff optimizer == "cmaes")
        "cma_max_gens":           0,
        "cma_sigma0":             float(start_sigma),
        "cma_popsize":            6,
        "cma_tolx":               0.005,
        "cma_tolfun":             1e-6,
        "cma_seed":               int(seed),
        # SNR ratchet — hold modest for panel runs
        "snr_floor":              1.0,
        "indist_stop_snr":        0.5,
        "snr_target":             3.0,
        "snr_max_traj_factor":    4,
        # Output / I/O
        "save_every":             5,
        "keep_mc_subdirs":        False,
        # Headless: panel cells are subprocesses, no terminal UI.
        "dashboard":              False,
        "no_vis":                 True,
    }


def _apply_speedup(cfg: dict, speedup: str, n_workers: int = 1) -> dict:
    """Layer speedup-specific overrides on top of the base cfg."""
    cfg = dict(cfg)
    if speedup == "legacy":
        cfg.update(betac_cache=False, n_workers=1, reweight=False)
    elif speedup == "s2":
        cfg.update(betac_cache=True, n_workers=1, reweight=False)
    elif speedup == "s3":
        cfg.update(betac_cache=False, n_workers=int(n_workers), reweight=False)
    elif speedup == "s4":
        cfg.update(betac_cache=False, n_workers=1, reweight=True)
    elif speedup == "s5":
        cfg.update(betac_cache=False, n_workers=1, reweight=False,
                   optimizer="bo")
    else:
        raise ValueError(f"unknown speedup id {speedup!r}")
    return cfg


# ---------------------------------------------------------------------------
# Cell execution
# ---------------------------------------------------------------------------

def _cell_dir(speedup: str, tag: str, start_id: str, seed: int,
              root: Path) -> Path:
    return root / speedup / tag / start_id / f"seed_{seed:03d}"


def _write_cfg(cfg: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Convert tuples to lists for JSON.
    serial = {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()}
    path.write_text(json.dumps(serial, indent=2))


def _run_cell(cfg: dict, cell_dir: Path, dry_run: bool = False,
              run_py: Optional[Path] = None) -> dict:
    """Execute one panel cell.

    Returns a dict with status + summary; writes ``summary.json`` and
    ``eval_log.jsonl`` into ``cell_dir`` via the standard run.py path.
    """
    cell_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cell_dir / "config.json"
    _write_cfg(cfg, cfg_path)
    if dry_run:
        return {"status": "DRY-RUN", "cell_dir": str(cell_dir)}

    # Driver pattern: import run.run() programmatically with this cfg.
    # We use subprocess so process state (matplotlib backend, scratch
    # dirs, OMP_NUM_THREADS) is fully isolated per cell.
    run_py = run_py or (PKG / "run.py")
    if not run_py.exists():
        raise FileNotFoundError(f"run.py not found at {run_py}")
    env = os.environ.copy()
    env["KOPT_PNP_CFG_OVERRIDE"] = str(cfg_path)
    env["KOPT_PNP_OUTPUT_DIR"] = str(cell_dir)
    env["MPLBACKEND"] = "Agg"
    t0 = time.time()
    try:
        result = subprocess.run(
            [sys.executable, str(run_py)],
            cwd=str(PKG), env=env,
            capture_output=True, text=True, timeout=8 * 3600,
        )
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT", "cell_dir": str(cell_dir),
                "wall_s": time.time() - t0}
    wall = time.time() - t0
    if result.returncode != 0:
        (cell_dir / "stderr.txt").write_text(result.stderr or "")
        (cell_dir / "stdout.txt").write_text(result.stdout or "")
        return {"status": "FAILED", "cell_dir": str(cell_dir),
                "wall_s": wall, "returncode": result.returncode}
    summary_path = cell_dir / "summary.json"
    if not summary_path.exists():
        return {"status": "NO-SUMMARY", "cell_dir": str(cell_dir),
                "wall_s": wall}
    summary = json.loads(summary_path.read_text())
    summary["status_panel"] = "OK"
    summary["wall_panel_s"] = wall
    return summary


# ---------------------------------------------------------------------------
# Panel orchestration
# ---------------------------------------------------------------------------

def build_matrix(speedup: str, geoms: list, starts: list,
                 seeds: int, optimizer: str = "cmaes") -> list:
    """Enumerate every (geom, start, seed) cell as a list of (cfg, cell_dir)."""
    matrix = []
    for tag in geoms:
        geom = get_geometry(tag)
        for start_id, x0, sigma0 in CMAES_START_DISTRIBUTIONS:
            if start_id not in starts:
                continue
            for seed in range(1, int(seeds) + 1):
                if optimizer == "nelder_mead":
                    max_evals = NM_CONTROL_MAX_EVALS
                else:
                    max_evals = PANEL_CMAES_MAX_EVALS
                cfg = _base_cfg(geom, x0, sigma0, seed, optimizer, max_evals)
                cfg = _apply_speedup(cfg, speedup)
                matrix.append((tag, start_id, seed, cfg))
    return matrix


def run_panel(speedup: str, *, geoms=None, starts=None, seeds=None,
              n_workers: int = 1, dry_run: bool = False,
              optimizer: str = "cmaes",
              root: Optional[Path] = None) -> list:
    """Execute the full panel for one speedup."""
    geoms = geoms or list_geometries()
    starts = starts or [s[0] for s in CMAES_START_DISTRIBUTIONS]
    seeds = seeds or PANEL_SEEDS_PER_CELL
    root = Path(root) if root else (PKG / "results" / "_validation")
    matrix = build_matrix(speedup, geoms, starts, seeds, optimizer)
    print(f"[panel] speedup={speedup} optimizer={optimizer} "
          f"matrix size={len(matrix)} ({len(geoms)} geoms × "
          f"{len(starts)} starts × {seeds} seeds)")
    if dry_run:
        for (tag, sid, seed, cfg) in matrix:
            cell_dir = _cell_dir(speedup, tag, sid, seed, root)
            print(f"  DRY {tag:<8} {sid:<8} seed={seed} -> {cell_dir}")
        return []
    results = []
    if n_workers <= 1:
        for (tag, sid, seed, cfg) in matrix:
            cell_dir = _cell_dir(speedup, tag, sid, seed, root)
            print(f"[panel] running {tag}/{sid}/seed_{seed} ...")
            res = _run_cell(cfg, cell_dir, dry_run=False)
            res.update(tag=tag, start=sid, seed=seed)
            results.append(res)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futures = {}
            for (tag, sid, seed, cfg) in matrix:
                cell_dir = _cell_dir(speedup, tag, sid, seed, root)
                fut = ex.submit(_run_cell, cfg, cell_dir, False)
                futures[fut] = (tag, sid, seed)
            for fut in as_completed(futures):
                tag, sid, seed = futures[fut]
                res = fut.result()
                res.update(tag=tag, start=sid, seed=seed)
                print(f"[panel] DONE {tag}/{sid}/seed_{seed} "
                      f"({res.get('status_panel', res.get('status'))})")
                results.append(res)
    write_panel_md(speedup, results, root)
    return results


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_panel_md(speedup: str, results: list, root: Path) -> Path:
    """Write the per-speedup PANEL.md summary table."""
    out = root / speedup / "PANEL.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = ["| Speedup | Geom | Start | Seed | r1* | r2* | dist | cost* | n_evals | wall_s | Status |",
            "|---------|------|-------|------|-----|-----|------|-------|---------|--------|--------|"]
    for r in sorted(results, key=lambda r: (r.get("tag", ""),
                                             r.get("start", ""),
                                             r.get("seed", 0))):
        tag = r.get("tag", "?")
        try:
            geom = get_geometry(tag)
            xstar = geom.get("expected_optimum")
        except KeyError:
            xstar = None
        r1 = r.get("best_r1", float("nan"))
        r2 = r.get("best_r2", float("nan"))
        if xstar:
            dist = ((r1 - xstar[0]) ** 2 + (r2 - xstar[1]) ** 2) ** 0.5
            dist_s = f"{dist:.3f}"
        else:
            dist_s = "—"
        rows.append(
            f"| {speedup} | {tag} | {r.get('start', '?')} | {r.get('seed', '?')} "
            f"| {r1:.3f} | {r2:.3f} | {dist_s} | {r.get('best_cost', float('nan')):.3g} "
            f"| {r.get('n_evals', '?')} | {r.get('wall_panel_s', float('nan')):.0f} "
            f"| {r.get('status_panel', r.get('status', '?'))} |"
        )
    out.write_text("\n".join(rows) + "\n")
    print(f"[panel] wrote {out}")
    return out


def write_latest_md(root: Path) -> Path:
    """Roll up every per-speedup PANEL.md into PANEL_LATEST.md."""
    latest = root / "PANEL_LATEST.md"
    root.mkdir(parents=True, exist_ok=True)
    chunks = ["# Validation panel — latest results\n"]
    for sp in sorted((root.iterdir() if root.exists() else [])):
        if not sp.is_dir():
            continue
        md = sp / "PANEL.md"
        if md.exists():
            chunks.append(f"\n## Speedup `{sp.name}`\n")
            chunks.append(md.read_text())
    latest.write_text("\n".join(chunks))
    print(f"[panel] rolled up -> {latest}")
    return latest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main():
    p = argparse.ArgumentParser(description="Validation panel driver")
    p.add_argument("--speedup", choices=("legacy", "s2", "s3", "s4", "s5", "all"),
                   help="Which speedup to validate")
    p.add_argument("--geoms", default=None,
                   help="Comma-separated subset of geometry tags")
    p.add_argument("--starts", default=None,
                   help="Comma-separated subset of start IDs (near,offset,far_lo,far_hi,skew)")
    p.add_argument("--seeds", type=int, default=None,
                   help="Number of independent seeds per cell")
    p.add_argument("--optimizer", default="cmaes",
                   choices=("nelder_mead", "cmaes", "bo"))
    p.add_argument("--n-workers", type=int, default=1,
                   help="Cells to execute in parallel (independent processes)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--report", action="store_true",
                   help="Roll up existing PANEL.md files into PANEL_LATEST.md and exit")
    p.add_argument("--root", default=None,
                   help="Root output directory (default: results/_validation/)")
    args = p.parse_args()

    root = Path(args.root) if args.root else (PKG / "results" / "_validation")

    if args.report:
        write_latest_md(root)
        return

    if not args.speedup:
        p.error("--speedup is required (or pass --report)")

    geoms = args.geoms.split(",") if args.geoms else None
    starts = args.starts.split(",") if args.starts else None

    speedups = (["legacy", "s2", "s3", "s4", "s5"]
                if args.speedup == "all" else [args.speedup])
    for sp in speedups:
        run_panel(sp, geoms=geoms, starts=starts, seeds=args.seeds,
                  optimizer=args.optimizer,
                  n_workers=args.n_workers, dry_run=args.dry_run, root=root)
    write_latest_md(root)


if __name__ == "__main__":
    _main()
