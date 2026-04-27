#!/usr/bin/env python3
"""
check_overnight.py — summarise the 4 runs of run_overnight.sh.

Reads <root>/<run>/summary.json and (optionally) fallback_log.jsonl for
each of A_legacy, B_reweight, C_reweight_fallback, D_parallel; emits a
single Markdown table at <root>/OVERNIGHT_REPORT.md.

Usage:
    python check_overnight.py [<root>]    # default root = results/_overnight/
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter

RUNS = ("A_legacy", "B_reweight", "C_reweight_fallback", "D_parallel")


def _load_summary(run_dir):
    p = os.path.join(run_dir, "summary.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def _fallback_breakdown(run_dir):
    p = os.path.join(run_dir, "fallback_log.jsonl")
    if not os.path.exists(p):
        return Counter()
    counts = Counter()
    with open(p) as f:
        for line in f:
            try:
                counts[json.loads(line).get("status", "?")] += 1
            except Exception:
                continue
    return counts


def main(argv=None) -> int:
    root = (argv[1] if argv and len(argv) > 1
            else "results/_overnight")
    if not os.path.isdir(root):
        print(f"[check] {root} not found", file=sys.stderr)
        return 1

    rows = []
    baseline_wall = None
    for name in RUNS:
        rd = os.path.join(root, name)
        s  = _load_summary(rd)
        fb = _fallback_breakdown(rd)
        if s is None:
            rows.append((name, None, None, None, None, None, None, fb))
            continue
        wall = float(s.get("wall_total_s", 0.0))
        n_ev = int(s.get("n_evals", 0))
        spe  = wall / max(1, n_ev)
        if name == "A_legacy":
            baseline_wall = wall
        speed = (baseline_wall / wall) if (baseline_wall and wall) else None
        rows.append((name,
                     n_ev, wall, spe,
                     s.get("best_r1"), s.get("best_r2"),
                     s.get("best_cost"), s.get("best_beta_c"),
                     speed, fb))

    out = [f"# Overnight validation report",
           f"Root: `{root}`", ""]
    out.append("| Run | n_eval | wall (s) | s/eval | best (r1, r2) | best β_c | "
               "best cost | speedup vs A | reweight | fallback | failed |")
    out.append("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        if r[1] is None:
            out.append(f"| {r[0]} | — | — | — | — | — | — | — | — | — | — |")
            continue
        name, n_ev, wall, spe, r1, r2, cost, bc, speed, fb = r
        out.append(
            f"| {name} | {n_ev} | {wall:.1f} | {spe:.2f} | "
            f"({r1:+.4f}, {r2:+.4f}) | {bc:.6f} | {cost:.3e} | "
            f"{('%.2fx' % speed) if speed else '—'} | "
            f"{fb.get('reweight', 0)} | {fb.get('fallback_3pass', 0)} | "
            f"{fb.get('fallback_failed', 0)} |"
        )

    out += ["",
            "## Pass criteria (per K_from_Optimizer_Production design)",
            "- **B_reweight**: best (r1, r2) within 0.05 of (1.0, 1.0); "
            "fallback count = 0 expected.",
            "- **C_reweight_fallback**: at least one fallback fired; "
            "best β_c agrees with A within 1e-3.",
            "- **D_parallel**: speedup vs A ≥ 2× on 4 workers (depends on hardware).",
            ""]
    txt = "\n".join(out)
    rep = os.path.join(root, "OVERNIGHT_REPORT.md")
    with open(rep, "w") as f:
        f.write(txt)
    print(txt)
    print(f"\n[check] wrote {rep}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
