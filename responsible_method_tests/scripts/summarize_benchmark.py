#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import numpy as np

from align_manifolds import COMPARISON_MODES, build_benchmark_comparison


def _format_float(value: Any, digits: int = 4) -> str:
    number = float(value)
    if not np.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def _label_from_metrics(metrics: dict[str, Any]) -> str:
    benchmark_id = str(metrics["benchmark_id"])
    run_tag = str(metrics.get("run_tag", "")).strip()
    if not run_tag:
        return benchmark_id
    if "high_stats_taylor2_refit" in run_tag:
        return f"{benchmark_id} (taylor2 refit)"
    if "high_stats" in run_tag:
        return f"{benchmark_id} (high stat)"
    if "integer_multiples" in run_tag:
        return f"{benchmark_id} (baseline)"
    return f"{benchmark_id} ({run_tag})"


def _find_results_root(paths: list[str]) -> str:
    common = os.path.commonpath(paths)
    probe = common
    while probe not in ("", os.path.sep):
        if os.path.basename(probe) == "results":
            return probe
        parent = os.path.dirname(probe)
        if parent == probe:
            break
        probe = parent
    return os.path.dirname(paths[0])


def _artifact_link(summary_path: str, benchmark_manifest_path: str, benchmark_id: str, suffix: str) -> str:
    artifact_path = os.path.join(os.path.dirname(benchmark_manifest_path), f"{benchmark_id}_{suffix}")
    if not os.path.exists(artifact_path):
        return "n/a"
    rel_path = os.path.relpath(artifact_path, start=os.path.dirname(summary_path))
    return f"[{os.path.basename(artifact_path)}]({rel_path})"


def _build_markdown(
    summary_path: str,
    benchmark_manifest_paths: list[str],
    metrics_list: list[dict[str, Any]],
    *,
    comparison_mode: str,
) -> str:
    lines: list[str] = []
    lines.append("# Continuum Benchmark Summary")
    lines.append("")
    if comparison_mode == "twisted_reference":
        lines.append(
            "Direct comparison uses the twisted continuum interpolant sampled on untwisted continuum points; the reverse-direction and common-grid columns are diagnostics only."
        )
        lines.append(
            "Interpretation rule: twisted-on-untwisted RMS z <= 1.0 means not distinguished within current continuum point errors, 1.0-2.0 is marginal, and > 2.0 is distinguishable."
        )
    else:
        lines.append("Direct twisted-vs-untwisted comparison uses shared-cell interpolation and directional continuum-point residuals.")
        lines.append("Interpretation rule: max directional RMS z <= 1.0 means not distinguished within current continuum point errors, 1.0-2.0 is marginal, and > 2.0 is distinguishable.")
    lines.append("")
    lines.append("| benchmark | twisted alpha | twisted chi2/dof | untwisted alpha | untwisted chi2/dof | tw->un RMS z | un->tw RMS z | common-grid rel RMS | interpretation | comparison | metrics |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |")

    for benchmark_manifest_path, metrics in zip(benchmark_manifest_paths, metrics_list):
        benchmark_id = str(metrics["benchmark_id"])
        figure_link = _artifact_link(summary_path, benchmark_manifest_path, benchmark_id, "twisted_vs_untwisted_common_grid.png")
        metrics_link = _artifact_link(summary_path, benchmark_manifest_path, benchmark_id, "twisted_vs_untwisted_common_grid.json")
        lines.append(
            "| "
            + _label_from_metrics(metrics)
            + " | "
            + _format_float(metrics["twisted"]["alpha"])
            + " | "
            + _format_float(metrics["twisted"]["chi2_per_dof"])
            + " | "
            + _format_float(metrics["untwisted"]["alpha"])
            + " | "
            + _format_float(metrics["untwisted"]["chi2_per_dof"])
            + " | "
            + _format_float(metrics["twisted_on_untwisted"]["rms_z"])
            + " | "
            + _format_float(metrics["untwisted_on_twisted"]["rms_z"])
            + " | "
            + _format_float(metrics["common_grid"]["relative_rms"])
            + " | "
            + str(metrics["interpretation"])
            + " | "
            + figure_link
            + " | "
            + metrics_link
            + " |"
        )
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a compact markdown summary table for one or more continuum benchmark manifests."
    )
    parser.add_argument(
        "--benchmark-manifest",
        action="append",
        required=True,
        help="Path to manifest_geometry_*.json; pass once per benchmark row",
    )
    parser.add_argument("--output", default=None, help="Optional output markdown path")
    parser.add_argument("--grid-size", type=int, default=180, help="Grid size passed to the shared-cell comparer")
    parser.add_argument(
        "--comparison-mode",
        choices=COMPARISON_MODES,
        default="symmetric",
        help="How to turn pointwise residuals into the distinguishability verdict",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    benchmark_manifest_paths = [os.path.abspath(path) for path in args.benchmark_manifest]
    metrics_list = [
        build_benchmark_comparison(
            path,
            grid_size=max(int(args.grid_size), 20),
            comparison_mode=str(args.comparison_mode),
        )[0]
        for path in benchmark_manifest_paths
    ]

    output_path = args.output
    if output_path is None:
        results_root = _find_results_root(benchmark_manifest_paths)
        output_path = os.path.join(results_root, "CONTINUUM_BENCHMARK_SUMMARY.md")
    output_path = os.path.abspath(output_path)

    markdown = _build_markdown(
        output_path,
        benchmark_manifest_paths,
        metrics_list,
        comparison_mode=str(args.comparison_mode),
    )
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(markdown)
        handle.write("\n")

    print(f"wrote {output_path}")
    for metrics in metrics_list:
        print(
            f"row={_label_from_metrics(metrics)} comparison_mode={metrics['comparison_mode']} interpretation={metrics['interpretation']} "
            f"tw->un_rms_z={metrics['twisted_on_untwisted']['rms_z']:.6f} "
            f"un->tw_rms_z={metrics['untwisted_on_twisted']['rms_z']:.6f}"
        )


if __name__ == "__main__":
    main()