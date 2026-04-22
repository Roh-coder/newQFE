#!/usr/bin/env python3
"""
analyze_benchmark.py — produce comparison plots and a summary table from
the per-method eval_log.jsonl files in results/.

Outputs (all in results/):
  - convergence.png           : best cost vs evaluation # for each method
  - trajectories.png          : (r1, r2) scatter, one panel per method
  - distance_history.png      : distance from (1, 1) vs evaluation #
  - benchmark_table.md        : summary table for the README
"""

from __future__ import annotations

import argparse
import json
import os
import math
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
_ap = argparse.ArgumentParser(add_help=False)
_ap.add_argument("--results-dir", default=None)
_args, _ = _ap.parse_known_args()
RESULTS = _args.results_dir or os.path.join(HERE, "results")
TRUE_R = (1.0, 1.0)

METHODS = ["nelder_mead", "powell", "bfgs_fd", "gp", "cma"]
COLORS = {"nelder_mead": "C0", "powell": "C1", "bfgs_fd": "C2",
          "gp": "C3", "cma": "C4"}


def load(name):
    p = os.path.join(RESULTS, name, "eval_log.jsonl")
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return [json.loads(l) for l in f]


def dist(r1, r2):
    return math.hypot(r1 - TRUE_R[0], r2 - TRUE_R[1])


# ---------------------------------------------------------------------------
#  Plot 1: convergence — running best cost vs eval
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.0, 5.0))
for name in METHODS:
    rows = load(name)
    if not rows:
        continue
    costs = np.array([r["cost"] for r in rows])
    running_best = np.minimum.accumulate(costs)
    idx = np.arange(1, len(costs) + 1)
    ax.plot(idx, running_best, "o-", ms=4, lw=1.4,
            color=COLORS[name], label=name)
ax.set_yscale("log")
ax.set_xlabel("evaluation #")
ax.set_ylabel("running best L² cost")
ax.set_title("Optimizer convergence (lower is better)")
ax.grid(alpha=0.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "convergence.png"), dpi=120)
plt.close(fig)


# ---------------------------------------------------------------------------
#  Plot 2: trajectories — (r1, r2) scatter per method
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(1, len(METHODS), figsize=(4.0 * len(METHODS), 4.2),
                         sharex=True, sharey=True)
for ax, name in zip(axes, METHODS):
    rows = load(name)
    if not rows:
        ax.text(0.5, 0.5, "(no data)", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(name)
        continue
    r1 = np.array([r["r1"] for r in rows])
    r2 = np.array([r["r2"] for r in rows])
    cost = np.array([r["cost"] for r in rows])
    sc = ax.scatter(r1, r2, c=np.log10(cost), cmap="plasma",
                    s=40, edgecolor="k", lw=0.3)
    ax.plot(r1, r2, "-", color="lightgray", lw=0.5, zorder=0)
    ax.plot(*TRUE_R, "+", color="lime", ms=14, mew=2.5, label="true (1,1)")
    i_best = int(np.argmin(cost))
    ax.plot(r1[i_best], r2[i_best], "o", mfc="none", mec="cyan",
            ms=14, mew=2, label="best by cost")
    ax.set_title(f"{name}  (n={len(rows)})")
    ax.set_xlabel("r₁")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best")
axes[0].set_ylabel("r₂")
fig.suptitle("Optimizer trajectories in (r₁, r₂) — color = log₁₀(cost)",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(RESULTS, "trajectories.png"), dpi=120)
plt.close(fig)


# ---------------------------------------------------------------------------
#  Plot 3: distance from (1, 1) vs eval (with running min)
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.0, 5.0))
for name in METHODS:
    rows = load(name)
    if not rows:
        continue
    d = np.array([dist(r["r1"], r["r2"]) for r in rows])
    idx = np.arange(1, len(d) + 1)
    ax.plot(idx, d, "o", ms=3, color=COLORS[name], alpha=0.4)
    ax.plot(idx, np.minimum.accumulate(d), "-", lw=1.6,
            color=COLORS[name], label=name)
ax.set_xlabel("evaluation #")
ax.set_ylabel("‖(r₁, r₂) − (1, 1)‖")
ax.set_yscale("log")
ax.set_title("Distance from true couplings (running min, line)")
ax.grid(alpha=0.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "distance_history.png"), dpi=120)
plt.close(fig)


# ---------------------------------------------------------------------------
#  Markdown table
# ---------------------------------------------------------------------------

lines = [
    "| Method | n | best cost | best (r₁, r₂) | best dist | final-3 ⟨dist⟩ | min dist seen | total traj | wall (s) |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
]
for name in METHODS:
    rows = load(name)
    if not rows:
        lines.append(f"| {name} | 0 | — | — | — | — | — | — | — |")
        continue
    n = len(rows)
    costs = np.array([r["cost"] for r in rows])
    r1 = np.array([r["r1"] for r in rows])
    r2 = np.array([r["r2"] for r in rows])
    i_best = int(np.argmin(costs))
    best_cost = costs[i_best]
    best_d = dist(r1[i_best], r2[i_best])
    final = rows[-min(3, n):]
    final_md = np.mean([dist(r["r1"], r["r2"]) for r in final])
    min_d = min(dist(rr["r1"], rr["r2"]) for rr in rows)
    total_traj = sum(r["n_traj_prod"] + r["n_traj_scan_total"] for r in rows)
    wall = sum(r["wall_time_s"] for r in rows)
    lines.append(
        f"| {name} | {n} | {best_cost:.2e} "
        f"| ({r1[i_best]:.3f}, {r2[i_best]:.3f}) "
        f"| {best_d:.3f} | {final_md:.3f} | {min_d:.3f} "
        f"| {total_traj:,} | {wall:.0f} |"
    )

with open(os.path.join(RESULTS, "benchmark_table.md"), "w") as f:
    f.write("\n".join(lines) + "\n")

print("wrote:")
for fn in ("convergence.png", "trajectories.png",
           "distance_history.png", "benchmark_table.md"):
    print(" ", os.path.join(RESULTS, fn))
