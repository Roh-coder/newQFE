#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from align_manifolds import build_method_pair_comparison
from plot_fss_interpolated_manifolds import _load_json


def _format_float(value: Any, digits: int = 4) -> str:
    number = float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def _resolve_target_twisted_manifest(path: str) -> tuple[str, dict[str, Any]]:
    payload = _load_json(path)
    if "methods" in payload:
        methods = payload.get("methods", {})
        if "twisted" not in methods:
            raise ValueError(f"benchmark manifest has no twisted method: {path}")
        return os.path.abspath(str(methods["twisted"])), payload
    if str(payload.get("method")) != "twisted":
        raise ValueError(f"method manifest is not twisted: {path}")
    return os.path.abspath(path), payload


def _resolve_candidate_untwisted_manifest(path: str) -> tuple[str, dict[str, Any]]:
    payload = _load_json(path)
    if "methods" in payload:
        methods = payload.get("methods", {})
        if "untwisted" not in methods:
            raise ValueError(f"benchmark manifest has no untwisted method: {path}")
        return os.path.abspath(str(methods["untwisted"])), payload
    if str(payload.get("method")) != "untwisted":
        raise ValueError(f"method manifest is not untwisted: {path}")
    return os.path.abspath(path), payload


def _candidate_score(metrics: dict[str, Any]) -> float:
    return float(
        max(
            metrics["twisted_on_untwisted"]["rms_z"],
            metrics["untwisted_on_twisted"]["rms_z"],
        )
    )


def _candidate_label(payload: dict[str, Any], candidate_manifest_path: str) -> str:
    if "benchmark_id" in payload:
        label = str(payload["benchmark_id"])
        couplings = payload.get("couplings")
        if isinstance(couplings, dict) and "r1" in couplings and "r2" in couplings:
            return f"{label} (r1={float(couplings['r1']):.3f}, r2={float(couplings['r2']):.3f})"
        return label
    return os.path.basename(candidate_manifest_path)


def _build_rows(
    *,
    target_twisted_manifest_path: str,
    target_context: dict[str, Any],
    candidate_manifest_paths: list[str],
    grid_size: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target_id = str(target_context.get("benchmark_id", os.path.basename(target_twisted_manifest_path)))
    target_description = str(target_context.get("description", target_context.get("benchmark_description", "")))
    target_run_tag = str(target_context.get("run_tag", ""))

    for candidate_manifest_path in candidate_manifest_paths:
        untwisted_manifest_path, candidate_payload = _resolve_candidate_untwisted_manifest(candidate_manifest_path)
        metrics, _ = build_method_pair_comparison(
            twisted_method_manifest_path=target_twisted_manifest_path,
            untwisted_method_manifest_path=untwisted_manifest_path,
            grid_size=grid_size,
            benchmark_id=target_id,
            description=target_description,
            run_tag=target_run_tag,
        )
        couplings = metrics["untwisted"].get("couplings", {})
        row = {
            "candidate_manifest": untwisted_manifest_path,
            "candidate_label": _candidate_label(candidate_payload, untwisted_manifest_path),
            "r1": float(couplings.get("r1", float("nan"))),
            "r2": float(couplings.get("r2", float("nan"))),
            "score": _candidate_score(metrics),
            "twisted_on_untwisted_rms_z": float(metrics["twisted_on_untwisted"]["rms_z"]),
            "untwisted_on_twisted_rms_z": float(metrics["untwisted_on_twisted"]["rms_z"]),
            "common_grid_rel_rms": float(metrics["common_grid"]["relative_rms"]),
            "common_grid_corr": float(metrics["common_grid"]["correlation"]),
            "interpretation": str(metrics["interpretation"]),
        }
        rows.append(row)

    rows.sort(key=lambda item: (item["score"], item["common_grid_rel_rms"], item["candidate_label"]))
    return rows


def _plot_scan(
    rows: list[dict[str, Any]],
    *,
    output_path: str,
    title: str,
    truth_r1: float | None,
    truth_r2: float | None,
    annotate_points: bool,
) -> str:
    r1 = np.asarray([float(row["r1"]) for row in rows], dtype=float)
    r2 = np.asarray([float(row["r2"]) for row in rows], dtype=float)
    score = np.asarray([float(row["score"]) for row in rows], dtype=float)
    common = np.asarray([float(row["common_grid_rel_rms"]) for row in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(title, fontsize=14, y=0.98)

    scatter0 = axes[0].scatter(r1, r2, c=score, s=120, cmap="viridis_r")
    best = rows[0]
    axes[0].scatter([best["r1"]], [best["r2"]], marker="*", s=260, color="red", edgecolors="black", linewidths=0.8)
    if truth_r1 is not None and truth_r2 is not None:
        axes[0].scatter([truth_r1], [truth_r2], marker="x", s=120, color="black", linewidths=2.0)
    if annotate_points:
        for row in rows:
            axes[0].text(float(row["r1"]) + 0.01, float(row["r2"]) + 0.01, row["candidate_label"].split()[0], fontsize=8)
    axes[0].set_title("max directional RMS z")
    axes[0].set_xlabel("r1 = k1 / k3")
    axes[0].set_ylabel("r2 = k2 / k3")
    axes[0].grid(True, alpha=0.25)
    fig.colorbar(scatter0, ax=axes[0], fraction=0.046, pad=0.04, label="score")

    scatter1 = axes[1].scatter(r1, r2, c=common, s=120, cmap="magma_r")
    axes[1].scatter([best["r1"]], [best["r2"]], marker="*", s=260, color="cyan", edgecolors="black", linewidths=0.8)
    if truth_r1 is not None and truth_r2 is not None:
        axes[1].scatter([truth_r1], [truth_r2], marker="x", s=120, color="black", linewidths=2.0)
    axes[1].set_title("common-grid relative RMS")
    axes[1].set_xlabel("r1 = k1 / k3")
    axes[1].set_ylabel("r2 = k2 / k3")
    axes[1].grid(True, alpha=0.25)
    fig.colorbar(scatter1, ax=axes[1], fraction=0.046, pad=0.04, label="rel RMS")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _write_markdown(rows: list[dict[str, Any]], *, output_path: str, title: str) -> str:
    lines = [f"# {title}", "", "| rank | candidate | r1 | r2 | score | tw->un RMS z | un->tw RMS z | common-grid rel RMS | interpretation |", "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for idx, row in enumerate(rows, start=1):
        lines.append(
            "| "
            + str(idx)
            + " | "
            + row["candidate_label"]
            + " | "
            + _format_float(row["r1"], digits=3)
            + " | "
            + _format_float(row["r2"], digits=3)
            + " | "
            + _format_float(row["score"], digits=4)
            + " | "
            + _format_float(row["twisted_on_untwisted_rms_z"], digits=4)
            + " | "
            + _format_float(row["untwisted_on_twisted_rms_z"], digits=4)
            + " | "
            + _format_float(row["common_grid_rel_rms"], digits=4)
            + " | "
            + row["interpretation"]
            + " |"
        )
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank untwisted coupling candidates against a chosen twisted target torus using continuum-manifold comparison metrics."
    )
    parser.add_argument("--target", required=True, help="Twisted method manifest or benchmark manifest containing the twisted target")
    parser.add_argument(
        "--candidate",
        action="append",
        required=True,
        help="Untwisted method manifest or benchmark manifest containing an untwisted candidate; pass once per candidate",
    )
    parser.add_argument("--output-dir", default=None, help="Optional output directory for scan artifacts")
    parser.add_argument("--grid-size", type=int, default=180, help="Interpolation grid size used by pairwise comparison")
    parser.add_argument("--truth-r1", type=float, default=None, help="Optional known target r1 for plotting")
    parser.add_argument("--truth-r2", type=float, default=None, help="Optional known target r2 for plotting")
    parser.add_argument(
        "--annotate-points",
        action="store_true",
        help="Overlay candidate labels on the score panel",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    target_path = os.path.abspath(args.target)
    target_twisted_manifest_path, target_context = _resolve_target_twisted_manifest(target_path)
    candidate_manifest_paths = [os.path.abspath(path) for path in args.candidate]
    rows = _build_rows(
        target_twisted_manifest_path=target_twisted_manifest_path,
        target_context=target_context,
        candidate_manifest_paths=candidate_manifest_paths,
        grid_size=max(int(args.grid_size), 20),
    )

    target_id = str(target_context.get("benchmark_id", os.path.splitext(os.path.basename(target_twisted_manifest_path))[0]))
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else os.path.dirname(target_path)
    os.makedirs(output_dir, exist_ok=True)
    json_path = os.path.join(output_dir, f"{target_id}_coupling_fit_scan.json")
    md_path = os.path.join(output_dir, f"{target_id}_coupling_fit_scan.md")
    png_path = os.path.join(output_dir, f"{target_id}_coupling_fit_scan.png")

    payload = {
        "target": target_twisted_manifest_path,
        "grid_size": max(int(args.grid_size), 20),
        "best_candidate": rows[0] if rows else None,
        "rows": rows,
    }
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    title = f"{target_id}: local untwisted coupling scan"
    _write_markdown(rows, output_path=md_path, title=title)
    _plot_scan(
        rows,
        output_path=png_path,
        title=title,
        truth_r1=args.truth_r1,
        truth_r2=args.truth_r2,
        annotate_points=args.annotate_points,
    )

    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")
    if rows:
        best = rows[0]
        print(
            f"best_candidate={best['candidate_label']} score={best['score']:.6f} "
            f"r1={best['r1']:.6f} r2={best['r2']:.6f}"
        )


if __name__ == "__main__":
    main()