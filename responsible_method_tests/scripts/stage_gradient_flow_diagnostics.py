#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from matplotlib.patches import Ellipse

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_raw_reweighting_fss import _bundle_path, _ratio_observables  # noqa: E402
from plot_reweighting_refined_surfaces import _finite_difference, _score_row_from_single_donor  # noqa: E402
from test_geometry_match_grid_interpolation import CouplingPoint, _block_slices, _parse_selected_bundle, _selected_specs_for_size  # noqa: E402


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stage completed gradient-flow artifacts into the reweighting tree and compute "
            "a jackknife error estimate for the local score gradient vector."
        )
    )
    parser.add_argument("--run-root", required=True, help="Completed gradient-flow run directory.")
    parser.add_argument(
        "--stage-dir",
        default=None,
        help="Destination under the reweighting tree. Defaults to responsible_method_tests/reweighting/<run-root-name>.",
    )
    parser.add_argument(
        "--step-index",
        type=int,
        default=-1,
        help="Trajectory row to diagnose. Defaults to the last attempted step.",
    )
    parser.add_argument(
        "--n-blocks",
        type=int,
        default=16,
        help="Nominal jackknife block count for the saved selected-observable traces.",
    )
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def _resolve_step_row(rows: list[dict[str, str]], step_index: int) -> tuple[int, dict[str, str]]:
    if not rows:
        raise ValueError("trajectory.tsv contains no rows")
    resolved = int(step_index)
    if resolved < 0:
        resolved = len(rows) + resolved
    if resolved < 0 or resolved >= len(rows):
        raise IndexError(f"step_index {step_index} is out of range for {len(rows)} trajectory rows")
    return resolved, rows[resolved]


def _subset_payload(payload: dict[str, Any], lo: int, hi: int) -> dict[str, Any]:
    n_samples = int(np.asarray(payload["mag"], dtype=float).size)
    mask = np.ones(n_samples, dtype=bool)
    mask[int(lo) : int(hi)] = False
    subset = dict(payload)
    subset["corr"] = {key: np.asarray(values, dtype=float)[mask] for key, values in dict(payload["corr"]).items()}
    subset["mag"] = np.asarray(payload["mag"], dtype=float)[mask]
    subset["e1"] = np.asarray(payload["e1"], dtype=float)[mask]
    subset["e2"] = np.asarray(payload["e2"], dtype=float)[mask]
    subset["e3"] = np.asarray(payload["e3"], dtype=float)[mask]
    return subset


def _score_rows_at_point(
    donor_point: CouplingPoint,
    *,
    fit_sizes: list[int],
    payloads: dict[tuple[int, str], dict[str, Any]],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
    gradient_step: float,
    metric_key: str,
) -> dict[str, Any]:
    direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]] = {}
    reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]] = {}

    center_row = _score_row_from_single_donor(
        donor_point,
        donor_point,
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
    )
    r1_minus_row = _score_row_from_single_donor(
        donor_point,
        CouplingPoint(float(donor_point.r1) - float(gradient_step), float(donor_point.r2)),
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
    )
    r1_plus_row = _score_row_from_single_donor(
        donor_point,
        CouplingPoint(float(donor_point.r1) + float(gradient_step), float(donor_point.r2)),
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
    )
    r2_minus_row = _score_row_from_single_donor(
        donor_point,
        CouplingPoint(float(donor_point.r1), float(donor_point.r2) - float(gradient_step)),
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
    )
    r2_plus_row = _score_row_from_single_donor(
        donor_point,
        CouplingPoint(float(donor_point.r1), float(donor_point.r2) + float(gradient_step)),
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
    )

    center_value = float(center_row[metric_key])
    grad_r1 = float(
        _finite_difference(center_value, float(r1_minus_row[metric_key]), float(r1_plus_row[metric_key]), float(gradient_step))
    )
    grad_r2 = float(
        _finite_difference(center_value, float(r2_minus_row[metric_key]), float(r2_plus_row[metric_key]), float(gradient_step))
    )
    return {
        "center_score": center_value,
        "r1_minus_score": float(r1_minus_row[metric_key]),
        "r1_plus_score": float(r1_plus_row[metric_key]),
        "r2_minus_score": float(r2_minus_row[metric_key]),
        "r2_plus_score": float(r2_plus_row[metric_key]),
        "grad_r1": grad_r1,
        "grad_r2": grad_r2,
        "grad_norm": float(math.hypot(grad_r1, grad_r2)),
    }


def _jackknife_sigma(values: np.ndarray) -> float:
    if values.size < 2:
        return float("nan")
    mean_value = float(np.mean(values))
    return float(np.sqrt((values.size - 1.0) / values.size * np.sum(np.square(values - mean_value))))


def _jackknife_covariance(vectors: np.ndarray) -> np.ndarray:
    if vectors.ndim != 2 or vectors.shape[0] < 2:
        return np.full((2, 2), np.nan, dtype=float)
    centered = vectors - np.mean(vectors, axis=0, keepdims=True)
    return (vectors.shape[0] - 1.0) / vectors.shape[0] * (centered.T @ centered)


def _snr(value: float, sigma: float) -> float:
    if not math.isfinite(value) or not math.isfinite(sigma) or sigma <= 0.0:
        return float("nan")
    return float(abs(value) / sigma)


def _vector_snr(covariance: np.ndarray, grad_r1: float, grad_r2: float) -> float:
    if covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        return float("nan")
    det = float(np.linalg.det(covariance))
    if det <= 0.0:
        return float("nan")
    vector = np.asarray([float(grad_r1), float(grad_r2)], dtype=float)
    return float(math.sqrt(vector @ np.linalg.inv(covariance) @ vector))


def _plot_gradient_vector(output_path: Path, diagnostics: dict[str, Any]) -> None:
    grad_r1 = float(diagnostics["gradient"]["grad_r1"])
    grad_r2 = float(diagnostics["gradient"]["grad_r2"])
    sigma_r1 = float(diagnostics["gradient"]["grad_r1_sigma"])
    sigma_r2 = float(diagnostics["gradient"]["grad_r2_sigma"])
    covariance = np.asarray(diagnostics["gradient"]["covariance"], dtype=float)

    fig, ax = plt.subplots(figsize=(6.8, 6.0), constrained_layout=True)
    fig.patch.set_facecolor("#f8f5ef")
    ax.set_facecolor("#fffdfa")
    ax.axhline(0.0, color="#cbd5e1", linewidth=1.0)
    ax.axvline(0.0, color="#cbd5e1", linewidth=1.0)

    ax.annotate(
        "",
        xy=(grad_r1, grad_r2),
        xytext=(0.0, 0.0),
        arrowprops={"arrowstyle": "->", "linewidth": 2.2, "color": "#1d4ed8"},
        zorder=3,
    )
    ax.errorbar(
        [grad_r1],
        [grad_r2],
        xerr=[[sigma_r1], [sigma_r1]] if math.isfinite(sigma_r1) else None,
        yerr=[[sigma_r2], [sigma_r2]] if math.isfinite(sigma_r2) else None,
        fmt="o",
        color="#1d4ed8",
        ecolor="#60a5fa",
        elinewidth=1.4,
        capsize=4,
        markersize=6,
        zorder=4,
    )

    if covariance.shape == (2, 2) and np.all(np.isfinite(covariance)):
        eigvals, eigvecs = np.linalg.eigh(covariance)
        eigvals = np.maximum(eigvals, 0.0)
        if float(np.max(eigvals)) > 0.0:
            angle = float(np.degrees(np.arctan2(eigvecs[1, 1], eigvecs[0, 1])))
            ellipse = Ellipse(
                xy=(grad_r1, grad_r2),
                width=2.0 * float(np.sqrt(eigvals[1])),
                height=2.0 * float(np.sqrt(eigvals[0])),
                angle=angle,
                facecolor="#bfdbfe",
                edgecolor="#1d4ed8",
                linewidth=1.2,
                alpha=0.35,
                zorder=2,
            )
            ax.add_patch(ellipse)

    span_x = max(abs(grad_r1) + (3.0 * sigma_r1 if math.isfinite(sigma_r1) else 0.0), 1.0e-8)
    span_y = max(abs(grad_r2) + (3.0 * sigma_r2 if math.isfinite(sigma_r2) else 0.0), 1.0e-8)
    ax.set_xlim(-1.15 * span_x, 1.15 * span_x)
    ax.set_ylim(-1.15 * span_y, 1.15 * span_y)

    point = diagnostics["point"]
    ax.set_title(
        f"Local reweight gradient at step {int(diagnostics['step_index'])}: "
        f"(r1, r2)=({float(point['r1']):.6f}, {float(point['r2']):.6f})",
        fontsize=12,
        color="#0f172a",
    )
    ax.set_xlabel("d(score)/d(r1)")
    ax.set_ylabel("d(score)/d(r2)")
    ax.grid(True, color="#e2e8f0", linewidth=0.9, alpha=0.9)

    stats = diagnostics["gradient"]
    text = "\n".join(
        [
            f"grad_r1 = {float(stats['grad_r1']):+.3e} +/- {float(stats['grad_r1_sigma']):.2e}",
            f"grad_r2 = {float(stats['grad_r2']):+.3e} +/- {float(stats['grad_r2_sigma']):.2e}",
            f"|grad| = {float(stats['grad_norm']):.3e} +/- {float(stats['grad_norm_sigma']):.2e}",
            f"SNR_r1 = {float(stats['grad_r1_snr']):.2f}",
            f"SNR_r2 = {float(stats['grad_r2_snr']):.2f}",
            f"SNR_|grad| = {float(stats['grad_norm_snr']):.2f}",
        ]
    )
    ax.text(
        0.98,
        0.04,
        text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.5,
        family="monospace",
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.40", "fc": "#fff7ed", "ec": "#fdba74", "lw": 1.0},
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _write_tsv(output_path: Path, diagnostics: dict[str, Any]) -> None:
    gradient = diagnostics["gradient"]
    point = diagnostics["point"]
    headers = [
        "step_index",
        "point_r1",
        "point_r2",
        "metric_key",
        "gradient_step",
        "grad_r1",
        "grad_r1_sigma",
        "grad_r1_snr",
        "grad_r2",
        "grad_r2_sigma",
        "grad_r2_snr",
        "grad_norm",
        "grad_norm_sigma",
        "grad_norm_snr",
        "vector_snr_mahalanobis",
        "center_score",
        "r1_minus_score",
        "r1_plus_score",
        "r2_minus_score",
        "r2_plus_score",
        "jackknife_replicas",
        "source_run_root",
    ]
    row = {
        "step_index": int(diagnostics["step_index"]),
        "point_r1": float(point["r1"]),
        "point_r2": float(point["r2"]),
        "metric_key": str(diagnostics["metric_key"]),
        "gradient_step": float(diagnostics["gradient_step"]),
        "grad_r1": float(gradient["grad_r1"]),
        "grad_r1_sigma": float(gradient["grad_r1_sigma"]),
        "grad_r1_snr": float(gradient["grad_r1_snr"]),
        "grad_r2": float(gradient["grad_r2"]),
        "grad_r2_sigma": float(gradient["grad_r2_sigma"]),
        "grad_r2_snr": float(gradient["grad_r2_snr"]),
        "grad_norm": float(gradient["grad_norm"]),
        "grad_norm_sigma": float(gradient["grad_norm_sigma"]),
        "grad_norm_snr": float(gradient["grad_norm_snr"]),
        "vector_snr_mahalanobis": float(gradient["vector_snr_mahalanobis"]),
        "center_score": float(gradient["center_score"]),
        "r1_minus_score": float(gradient["r1_minus_score"]),
        "r1_plus_score": float(gradient["r1_plus_score"]),
        "r2_minus_score": float(gradient["r2_minus_score"]),
        "r2_plus_score": float(gradient["r2_plus_score"]),
        "jackknife_replicas": int(diagnostics["jackknife_replicas"]),
        "source_run_root": str(diagnostics["source_run_root"]),
    }
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def _copy_artifacts(run_root: Path, stage_dir: Path) -> None:
    for name in ["summary.json", "trajectory.tsv", "trajectory.png", "trajectory_path.png"]:
        source = run_root / name
        if source.is_file():
            shutil.copy2(source, stage_dir / name)


def main() -> None:
    args = parse_args()
    run_root = Path(args.run_root).resolve()
    stage_dir = Path(args.stage_dir).resolve() if args.stage_dir else RESPONSIBLE_ROOT / "reweighting" / run_root.name
    stage_dir.mkdir(parents=True, exist_ok=True)

    summary = _load_json(run_root / "summary.json")
    rows = _load_rows(run_root / "trajectory.tsv")
    resolved_step_index, step_row = _resolve_step_row(rows, int(args.step_index))

    metric_key = str(summary["metric_key"])
    if metric_key.endswith("_scaled"):
        raise ValueError(f"scaled metric {metric_key} is not supported by this diagnostic yet")

    fit_sizes = [int(value) for value in summary["fit_sizes"]]
    gradient_step = float(summary["gradient_step"])
    point = CouplingPoint(float(step_row["current_r1"]), float(step_row["current_r2"]))
    target_point = CouplingPoint(float(summary["target_point"]["r1"]), float(summary["target_point"]["r2"]))
    target_size = int(summary["args"]["target_size"])
    target_x = 1.0 / float(target_size)
    labels = [label for label, _ in _selected_specs_for_size(int(fit_sizes[0]))]

    payloads: dict[tuple[int, str], dict[str, Any]] = {}
    payload_slices: dict[int, list[tuple[int, int]]] = {}
    for size in fit_sizes:
        payload = _parse_selected_bundle(_bundle_path(run_root, int(size), point.tag), labels)
        payloads[(int(size), point.tag)] = payload
        payload_slices[int(size)] = _block_slices(int(np.asarray(payload["mag"], dtype=float).size), int(args.n_blocks))

    target_root = run_root / "_target"
    target_payload = _parse_selected_bundle(_bundle_path(target_root, target_size, target_point.tag), labels)
    target_slices = _block_slices(int(np.asarray(target_payload["mag"], dtype=float).size), int(args.n_blocks))
    target_rows = _ratio_observables(target_payload)

    full_sample = _score_rows_at_point(
        point,
        fit_sizes=fit_sizes,
        payloads=payloads,
        target_rows=target_rows,
        target_x=target_x,
        gradient_step=gradient_step,
        metric_key=metric_key,
    )

    n_replicas = min(len(target_slices), *(len(slices) for slices in payload_slices.values()))
    gradients: list[tuple[float, float, float]] = []
    for rep_index in range(n_replicas):
        subset_payloads: dict[tuple[int, str], dict[str, Any]] = {}
        for size in fit_sizes:
            lo, hi = payload_slices[int(size)][rep_index]
            subset_payloads[(int(size), point.tag)] = _subset_payload(payloads[(int(size), point.tag)], lo, hi)
        target_lo, target_hi = target_slices[rep_index]
        subset_target_rows = _ratio_observables(_subset_payload(target_payload, target_lo, target_hi))
        subset_result = _score_rows_at_point(
            point,
            fit_sizes=fit_sizes,
            payloads=subset_payloads,
            target_rows=subset_target_rows,
            target_x=target_x,
            gradient_step=gradient_step,
            metric_key=metric_key,
        )
        gradients.append(
            (
                float(subset_result["grad_r1"]),
                float(subset_result["grad_r2"]),
                float(subset_result["grad_norm"]),
            )
        )

    gradient_array = np.asarray(gradients, dtype=float)
    grad_r1_sigma = _jackknife_sigma(gradient_array[:, 0]) if gradient_array.size else float("nan")
    grad_r2_sigma = _jackknife_sigma(gradient_array[:, 1]) if gradient_array.size else float("nan")
    grad_norm_sigma = _jackknife_sigma(gradient_array[:, 2]) if gradient_array.size else float("nan")
    covariance = _jackknife_covariance(gradient_array[:, :2]) if gradient_array.size else np.full((2, 2), np.nan, dtype=float)

    diagnostics = {
        "source_run_root": str(run_root),
        "staged_output_dir": str(stage_dir),
        "step_index": int(resolved_step_index),
        "point": {"r1": float(point.r1), "r2": float(point.r2), "point_tag": point.tag},
        "metric_key": metric_key,
        "gradient_step": float(gradient_step),
        "jackknife_replicas": int(n_replicas),
        "trajectory_gradient": {
            "grad_r1": float(step_row["grad_r1"]),
            "grad_r2": float(step_row["grad_r2"]),
            "grad_norm": float(step_row["grad_norm"]),
        },
        "gradient": {
            "center_score": float(full_sample["center_score"]),
            "r1_minus_score": float(full_sample["r1_minus_score"]),
            "r1_plus_score": float(full_sample["r1_plus_score"]),
            "r2_minus_score": float(full_sample["r2_minus_score"]),
            "r2_plus_score": float(full_sample["r2_plus_score"]),
            "grad_r1": float(full_sample["grad_r1"]),
            "grad_r1_sigma": float(grad_r1_sigma),
            "grad_r1_snr": float(_snr(float(full_sample["grad_r1"]), float(grad_r1_sigma))),
            "grad_r2": float(full_sample["grad_r2"]),
            "grad_r2_sigma": float(grad_r2_sigma),
            "grad_r2_snr": float(_snr(float(full_sample["grad_r2"]), float(grad_r2_sigma))),
            "grad_norm": float(full_sample["grad_norm"]),
            "grad_norm_sigma": float(grad_norm_sigma),
            "grad_norm_snr": float(_snr(float(full_sample["grad_norm"]), float(grad_norm_sigma))),
            "vector_snr_mahalanobis": float(
                _vector_snr(covariance, float(full_sample["grad_r1"]), float(full_sample["grad_r2"]))
            ),
            "covariance": covariance.tolist(),
        },
    }

    _copy_artifacts(run_root, stage_dir)
    json_path = stage_dir / f"gradient_vector_step{resolved_step_index}.json"
    tsv_path = stage_dir / f"gradient_vector_step{resolved_step_index}.tsv"
    png_path = stage_dir / f"gradient_vector_step{resolved_step_index}.png"
    json_path.write_text(json.dumps(diagnostics, indent=2, sort_keys=True), encoding="utf-8")
    _write_tsv(tsv_path, diagnostics)
    _plot_gradient_vector(png_path, diagnostics)

    print(json.dumps({
        "staged_output_dir": str(stage_dir),
        "gradient_vector_json": str(json_path),
        "gradient_vector_tsv": str(tsv_path),
        "gradient_vector_png": str(png_path),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()