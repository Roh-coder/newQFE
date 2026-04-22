#!/usr/bin/env python3
"""
plot_criticality_check.py — confirm the reference and best-test lattices
are both tuned to their susceptibility peaks.

What it does
------------
* Loads the reference β-scan dump  (results/_reference/reference_betac_scan.json)
* For each optimizer's eval_log.jsonl, picks the lowest-cost evaluation and
  uses its persisted scan_betas / scan_chis / scan_chi_errs.
* Re-fits a Gram-Charlier curve to each scan via mc_engine._gc_fit and marks
  the inferred β_c.

Output
------
results/criticality_check.png
    A row of panels: leftmost = reference; one per optimizer method.
    Each panel shows the χ(β) data (errorbars), the GC fit, and a vertical
    line at β_c.  This makes it easy to eyeball whether each lattice is
    really sitting at its critical point.

Usage
-----
    python plot_criticality_check.py
    python plot_criticality_check.py --methods nelder_mead cma
"""

from __future__ import annotations

import argparse
import json
import math
import os
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mc_engine import _gc_fit, _gram_charlier


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ALL_METHODS = ["nelder_mead", "powell", "bfgs_fd", "gp", "cma"]


def _load_jsonl(path: str) -> list:
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def _best_eval(rows: list) -> Optional[dict]:
    rows = [r for r in rows if r.get("scan_betas")]
    if not rows:
        return None
    return min(rows, key=lambda r: r["cost"])


def _draw_panel(ax, betas, chis, errs, beta_c, *, title: str,
                color: str, marker: str = "o"):
    betas = np.asarray(betas, dtype=float)
    chis = np.asarray(chis, dtype=float)
    errs = np.asarray(errs, dtype=float)
    if betas.size == 0:
        ax.text(0.5, 0.5, "(no data)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return

    ax.errorbar(betas, chis, yerr=errs, fmt=marker, ms=4,
                color=color, ecolor=color, alpha=0.8,
                elinewidth=1.0, capsize=2, lw=0,
                label=f"χ(β)  (n={len(betas)})")

    gc_params, gc_mode = _gc_fit(betas.tolist(), chis.tolist(),
                                 float(beta_c) if beta_c is not None
                                 else float(betas[np.argmax(chis)]))
    if gc_params is not None:
        b_fit = np.linspace(betas.min(), betas.max(), 400)
        try:
            y_fit = _gram_charlier(b_fit, *gc_params)
            ax.plot(b_fit, y_fit, "-", color="crimson", lw=1.4,
                    label="GC fit")
        except Exception:
            gc_mode = beta_c

    bc_show = float(beta_c) if beta_c is not None else float(gc_mode)
    if bc_show is not None and np.isfinite(bc_show):
        ax.axvline(bc_show, color="black", ls="--", lw=1.0,
                   label=f"β_c ≈ {bc_show:.5f}")

    # Highlight the maximum-χ point in the data (criticality eyeball test).
    i_max = int(np.argmax(chis))
    ax.plot(betas[i_max], chis[i_max], "*", ms=14, mfc="gold",
            mec="black", mew=0.7, label=f"χ_max @ β={betas[i_max]:.5f}",
            zorder=10)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("β")
    ax.set_ylabel("χ")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best", framealpha=0.9)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--methods", nargs="*", default=ALL_METHODS,
                    help="optimizer subdirs under results/ to include")
    ap.add_argument("--results-dir", default=None,
                    help="results root directory (default: results/).")
    ap.add_argument("--out", default=None,
                    help="output PNG path (default: <results-dir>/criticality_check.png)")
    args = ap.parse_args()
    results_root = args.results_dir or os.path.join(HERE, "results")
    out = args.out or os.path.join(results_root, "criticality_check.png")

    # 1. Reference panel data — always from the canonical shared location.
    ref_path = os.path.join(HERE, "results", "_reference", "reference_betac_scan.json")
    if not os.path.exists(ref_path):
        print(f"WARNING: reference β-scan not found at {ref_path}.")
        print("  This file is written by run_benchmark.py when the reference")
        print("  is built fresh.  Delete results/_reference/ and re-run to")
        print("  regenerate it.")
        ref_payload = None
    else:
        with open(ref_path) as f:
            ref_payload = json.load(f)

    # 2. Per-method best eval.
    panels = []  # list of (title, betas, chis, errs, beta_c, color)
    if ref_payload is not None:
        panels.append((
            "REFERENCE (equilateral, untwisted)",
            ref_payload["scan_betas"],
            ref_payload["scan_chis"],
            ref_payload["scan_chi_errs"],
            ref_payload["beta_c"],
            "tab:blue",
        ))

    method_colors = {"nelder_mead": "C0", "powell": "C1", "bfgs_fd": "C2",
                     "gp": "C3", "cma": "C4"}

    for m in args.methods:
        rows = _load_jsonl(os.path.join(results_root, m, "eval_log.jsonl"))
        best = _best_eval(rows)
        if best is None:
            print(f"  [{m}] no scan data found, skipping")
            continue
        title = (f"{m} — best eval (#{best['eval_id']})\n"
                 f"r₁={best['r1']:.3f}, r₂={best['r2']:.3f}, "
                 f"cost={best['cost']:.2e}")
        panels.append((
            title,
            best["scan_betas"],
            best["scan_chis"],
            best.get("scan_chi_errs") or [0.0] * len(best["scan_betas"]),
            best["beta_c"],
            method_colors.get(m, "gray"),
        ))

    if not panels:
        print("Nothing to plot.")
        return

    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.4), squeeze=False)
    for ax, (title, b, c, e, bc, col) in zip(axes[0], panels):
        _draw_panel(ax, b, c, e, bc, title=title, color=col)

    fig.suptitle("Criticality check — susceptibility peak vs β\n"
                 "(★ = data peak,  -- = β_c estimate from GC fit)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(args.out, dpi=120)
    plt.close(fig)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
