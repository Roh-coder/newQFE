"""Head-to-head time-budget benchmark.

Runs each speedup variant on the SAME geometry, x0, MC master seed,
optimizer seed, then SIGTERMs after WALL_BUDGET seconds.  Reads back
eval_log.jsonl from each cell and reports best cost / r1 / r2 / dist
to (1, 1) / n_evals achieved within the budget.
"""
from __future__ import annotations
import json, os, signal, subprocess, sys, time
from pathlib import Path

PKG = Path("/workspaces/newQFE/K_from_Optimizer_Speedup_Upgrade")
OUT_ROOT = PKG / "results" / "_h2h"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

WALL_BUDGET = 300.0         # seconds per method
GRACE = 10.0                # extra seconds to let the child write final lines

# Common config: 12x12 self-consistency (optimum at (1,1)),
# offset start so there is real work to do.
BASE = {
    "ref_Lx": 12, "ref_Ly": 12, "ref_Tx": 0, "ref_Ty": 0,
    "test_Lx": 12, "test_Ly": 12, "test_Tx": 0, "test_Ty": 0,
    "ref_n_traj": 40000,
    "n_traj_prod": 10000,
    "n_traj_scan_coarse": 2000,
    "n_traj_scan_fine":  4000,
    "x0": [0.85, 1.15],
    "max_evals": 9999,
    "indist_stop_snr": 0.0,
    "snr_floor": 0.0,
    "no_vis": False,
    "dashboard": False,
    "save_every": 1,
    "master_seed": 42,
    "keep_mc_subdirs": False,
}

METHODS = {
    "s4_cmaes":   {**BASE, "optimizer": "cmaes", "cma_seed": 42,
                    "cma_sigma0": 0.15, "cma_popsize": 6, "reweight": True},
}


def run_one(name: str, cfg_patch: dict) -> dict:
    cell_dir = OUT_ROOT / name
    if cell_dir.exists():
        # Wipe per-method dir for a clean run; keeps the shared
        # _reference_* cache outside this dir.
        import shutil; shutil.rmtree(cell_dir)
    cell_dir.mkdir()

    # Wipe S2 cache before legacy / s4 / s5 so they don't accidentally
    # benefit from a previous s2_cold run.
    if not cfg_patch.get("betac_cache", False):
        cache_dir = PKG / "results" / "_betac_cache_Lx12_Ly12_Tx0_Ty0"
        if cache_dir.exists():
            import shutil; shutil.rmtree(cache_dir)

    cfg_path = cell_dir / "_override.json"
    cfg_path.write_text(json.dumps(cfg_patch))

    env = os.environ.copy()
    env["KOPT_PNP_CFG_OVERRIDE"] = str(cfg_path)
    env["KOPT_PNP_OUTPUT_DIR"]   = str(cell_dir)

    log_path = cell_dir / "stdout.log"
    log_fh = open(log_path, "w")
    t0 = time.monotonic()
    proc = subprocess.Popen(
        [sys.executable, "run.py"],
        cwd=str(PKG), env=env,
        stdout=log_fh, stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    print(f"[{name}] launched pid={proc.pid}; budget={WALL_BUDGET:.0f}s")

    deadline = t0 + WALL_BUDGET
    while True:
        rc = proc.poll()
        if rc is not None:
            wall = time.monotonic() - t0
            print(f"[{name}] exited early rc={rc} wall={wall:.1f}s")
            break
        if time.monotonic() >= deadline:
            print(f"[{name}] budget reached, sending SIGTERM")
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                proc.wait(timeout=GRACE)
            except subprocess.TimeoutExpired:
                print(f"[{name}] SIGKILL after grace")
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            wall = time.monotonic() - t0
            break
        time.sleep(0.5)
    log_fh.close()

    return parse_cell(name, cell_dir, wall)


def parse_cell(name, cell_dir, wall):
    log = cell_dir / "eval_log.jsonl"
    out = {"method": name, "wall_s": round(wall, 1),
           "n_evals": 0, "best_cost": None, "best_r1": None,
           "best_r2": None, "dist_from_truth": None,
           "best_eval_id": None, "best_at_wall_s": None}
    if not log.exists():
        out["status"] = "no eval_log"
        return out
    rows = []
    with open(log) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass  # truncated trailing line on SIGTERM
    out["n_evals"] = len(rows)
    if not rows:
        out["status"] = "no evals"
        return out
    # Best by cost
    best = min(rows, key=lambda r: r["cost"])
    out["best_cost"] = best["cost"]
    out["best_r1"]   = best["r1"]
    out["best_r2"]   = best["r2"]
    out["dist_from_truth"] = ((best["r1"] - 1.0)**2 + (best["r2"] - 1.0)**2) ** 0.5
    out["best_eval_id"] = best["eval_id"]
    # Cumulative wall at the best eval (sum of per-eval wall_time_s up to and including it)
    cum = 0.0
    for r in rows:
        cum += r.get("wall_time_s", 0.0)
        if r["eval_id"] == best["eval_id"]:
            out["best_at_wall_s"] = round(cum, 1)
            break
    out["status"] = "OK"
    return out


def render_table(results):
    cols = ["method", "n_evals", "best_cost", "best_r1", "best_r2",
            "dist_from_truth", "best_eval_id", "best_at_wall_s",
            "wall_s", "status"]
    widths = {c: max(len(c), max(len(str(_fmt(r.get(c)))) for r in results))
              for c in cols}
    print()
    print(" | ".join(c.ljust(widths[c]) for c in cols))
    print("-+-".join("-" * widths[c] for c in cols))
    for r in results:
        print(" | ".join(_fmt(r.get(c)).ljust(widths[c]) for c in cols))


def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        if v == 0 or 1e-3 <= abs(v) < 1e4:
            return f"{v:.4g}"
        return f"{v:.3e}"
    return str(v)


def main():
    print(f"Head-to-head: WALL_BUDGET={WALL_BUDGET}s, master_seed=42, "
          f"x0=(0.85, 1.15), geom=12x12 self-consistency")
    results = []
    for name, cfg in METHODS.items():
        res = run_one(name, cfg)
        print(f"[{name}] -> {res}")
        results.append(res)
    render_table(results)
    out_json = OUT_ROOT / "summary.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {out_json}")


if __name__ == "__main__":
    main()
