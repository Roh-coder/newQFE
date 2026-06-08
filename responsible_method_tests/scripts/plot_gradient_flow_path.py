#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a path-focused figure for a completed reweight gradient-flow run "
            "from its saved summary.json and trajectory.tsv artifacts."
        )
    )
    parser.add_argument(
        "--output-root",
        required=True,
        help="Gradient-flow result directory that contains summary.json and trajectory.tsv.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output PNG path. Defaults to <output-root>/trajectory_path.png.",
    )
    return parser.parse_args()


def _load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def _accepted_points(rows: list[dict[str, str]]) -> list[tuple[float, float]]:
    if not rows:
        return []
    points: list[tuple[float, float]] = [
        (float(rows[0]["current_r1"]), float(rows[0]["current_r2"]))
    ]
    for row in rows:
        if str(row.get("accepted", "")).lower() == "true":
            points.append((float(row["proposal_r1"]), float(row["proposal_r2"])))
    return points


def _rejected_point(rows: list[dict[str, str]]) -> tuple[float, float] | None:
    if not rows:
        return None
    last_row = rows[-1]
    if str(last_row.get("accepted", "")).lower() == "true":
        return None
    return (float(last_row["proposal_r1"]), float(last_row["proposal_r2"]))


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _score_range(rows: list[dict[str, str]]) -> tuple[float, float]:
    scores: list[float] = []
    for row in rows:
        for key in ("current_score", "predicted_score_from_current_donor", "direct_score_at_proposal"):
            value = _safe_float(row.get(key))
            if math.isfinite(value):
                scores.append(value)
    if not scores:
        return (0.0, 0.0)
    return (min(scores), max(scores))


def _bounds(points: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    x_min = min(xs)
    x_max = max(xs)
    y_min = min(ys)
    y_max = max(ys)
    x_pad = max(0.01, 0.10 * max(x_max - x_min, 1e-9))
    y_pad = max(0.01, 0.10 * max(y_max - y_min, 1e-9))
    return (x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad)


def _summary_lines(summary: dict[str, Any], rows: list[dict[str, str]]) -> list[str]:
    score_min, score_max = _score_range(rows)
    return [
        f"accepted/attempted: {summary['n_accepted_steps']}/{summary['n_attempted_steps']}",
        f"final score: {float(summary['final_score']):.3e}",
        f"final reason: {summary['final_reason']}",
        f"step size: {float(summary['args']['step_size']):.4f}",
        f"gradient step: {float(summary['gradient_step']):.4f}",
        f"score window: [{score_min:.3e}, {score_max:.3e}]",
    ]


def _has_covariance(rows: list[dict[str, str]]) -> bool:
    covariance_keys = ("gradient_cov_11", "gradient_cov_12", "gradient_cov_22")
    return any(any(row.get(key, "") != "" for key in covariance_keys) for row in rows)


def _covariance_correlation(row: dict[str, str]) -> float:
    explicit = _safe_float(row.get("gradient_cov_corr"))
    if math.isfinite(explicit):
        return explicit
    cov_11 = _safe_float(row.get("gradient_cov_11"))
    cov_12 = _safe_float(row.get("gradient_cov_12"))
    cov_22 = _safe_float(row.get("gradient_cov_22"))
    if not math.isfinite(cov_11) or not math.isfinite(cov_12) or not math.isfinite(cov_22):
        return float("nan")
    if cov_11 <= 0.0 or cov_22 <= 0.0:
        return float("nan")
    denom = math.sqrt(cov_11 * cov_22)
    if denom <= 0.0:
        return float("nan")
    return cov_12 / denom


def render_plot(output_root: Path, output_path: Path) -> Path:
    summary_path = output_root / "summary.json"
    trajectory_path = output_root / "trajectory.tsv"
    if not summary_path.is_file():
        raise FileNotFoundError(f"missing summary file: {summary_path}")
    if not trajectory_path.is_file():
        raise FileNotFoundError(f"missing trajectory file: {trajectory_path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = _load_rows(trajectory_path)
    accepted_points = _accepted_points(rows)
    if not accepted_points:
        raise ValueError("trajectory.tsv does not contain any plottable points")

    target = (
        float(summary["target_point"]["r1"]),
        float(summary["target_point"]["r2"]),
    )
    rejected_point = _rejected_point(rows)
    has_covariance = _has_covariance(rows)
    all_points = list(accepted_points) + [target]
    if rejected_point is not None:
        all_points.append(rejected_point)

    fig = plt.figure(figsize=(12.6, 8.6), constrained_layout=True)
    fig.patch.set_facecolor("#f8f5ef")
    gs = fig.add_gridspec(2, 2)
    trajectory_ax = fig.add_subplot(gs[0, 0])
    score_ax = fig.add_subplot(gs[0, 1])
    covariance_ax = fig.add_subplot(gs[1, 0])
    latest_cov_ax = fig.add_subplot(gs[1, 1])
    trajectory_ax.set_facecolor("#fffdfa")
    score_ax.set_facecolor("#fffdfa")
    covariance_ax.set_facecolor("#fffdfa")
    latest_cov_ax.set_facecolor("#fffdfa")

    x_values = [point[0] for point in accepted_points]
    y_values = [point[1] for point in accepted_points]
    trajectory_ax.plot(x_values, y_values, color="#1d4ed8", linewidth=2.8, alpha=0.9, zorder=2)

    colors = ["#0f172a"] + ["#2563eb"] * max(len(accepted_points) - 2, 0) + ["#1d4ed8"]
    if len(accepted_points) == 1:
        colors = ["#0f172a"]
    elif len(accepted_points) == 2:
        colors = ["#0f172a", "#1d4ed8"]
    trajectory_ax.scatter(
        x_values,
        y_values,
        s=78,
        c=colors,
        edgecolors="#eff6ff",
        linewidths=1.2,
        zorder=3,
    )

    for index, (x_value, y_value) in enumerate(accepted_points):
        trajectory_ax.annotate(
            str(index),
            xy=(x_value, y_value),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            color="#0f172a",
            bbox={"boxstyle": "round,pad=0.18", "fc": "#ffffff", "ec": "#cbd5e1", "lw": 0.8},
            zorder=4,
        )

    trajectory_ax.scatter(
        [target[0]],
        [target[1]],
        marker="x",
        s=190,
        color="#dc2626",
        linewidths=2.8,
        zorder=5,
    )
    trajectory_ax.annotate(
        "target",
        xy=target,
        xytext=(8, -12),
        textcoords="offset points",
        fontsize=10,
        color="#991b1b",
        fontweight="bold",
    )

    if rejected_point is not None:
        last_accepted = accepted_points[-1]
        trajectory_ax.plot(
            [last_accepted[0], rejected_point[0]],
            [last_accepted[1], rejected_point[1]],
            linestyle=(0, (4, 3)),
            color="#dc2626",
            linewidth=1.8,
            alpha=0.85,
            zorder=1,
        )
        trajectory_ax.scatter(
            [rejected_point[0]],
            [rejected_point[1]],
            marker="D",
            s=68,
            facecolors="#fee2e2",
            edgecolors="#dc2626",
            linewidths=1.2,
            zorder=4,
        )
        trajectory_ax.annotate(
            "rejected",
            xy=rejected_point,
            xytext=(8, 6),
            textcoords="offset points",
            fontsize=9,
            color="#991b1b",
        )

    fig.suptitle("Gradient-flow diagnostics toward the iso111 target", fontsize=18, y=0.985, color="#0f172a")

    trajectory_ax.set_title("parameter-space trajectory")
    trajectory_ax.set_xlabel("r1", fontsize=12, color="#0f172a")
    trajectory_ax.set_ylabel("r2", fontsize=12, color="#0f172a")
    trajectory_ax.grid(True, color="#cbd5e1", linewidth=0.9, alpha=0.85)
    for spine in trajectory_ax.spines.values():
        spine.set_color("#94a3b8")
        spine.set_linewidth(1.0)

    x_min, x_max, y_min, y_max = _bounds(all_points)
    trajectory_ax.set_xlim(x_min, x_max)
    trajectory_ax.set_ylim(y_min, y_max)

    summary_box = "\n".join(_summary_lines(summary, rows))
    trajectory_ax.text(
        0.98,
        0.96,
        summary_box,
        transform=trajectory_ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        family="monospace",
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.40", "fc": "#fff7ed", "ec": "#fdba74", "lw": 1.0},
    )

    step_ids = [int(row["step_index"]) for row in rows]
    current_scores = [_safe_float(row.get("current_score")) for row in rows]
    predicted_scores = [_safe_float(row.get("predicted_score_from_current_donor")) for row in rows]
    direct_scores = [_safe_float(row.get("direct_score_at_proposal")) for row in rows]
    score_ax.plot(step_ids, current_scores, "-o", color="#111827", label="current direct score")
    score_ax.plot(step_ids, predicted_scores, "--s", color="#ea580c", label="predicted proposal score")
    score_ax.plot(step_ids, direct_scores, ":^", color="#2563eb", label="direct proposal score")
    score_ax.set_title(f"{summary['metric_key']} by step")
    score_ax.set_xlabel("step index")
    score_ax.set_ylabel(summary["metric_key"])
    score_ax.grid(True, color="#cbd5e1", linewidth=0.9, alpha=0.85)
    score_ax.legend(fontsize=8)

    covariance_ax.set_title("gradient covariance by step")
    covariance_ax.set_xlabel("step index")
    covariance_ax.set_ylabel("covariance")
    covariance_ax.grid(True, color="#cbd5e1", linewidth=0.9, alpha=0.85)
    covariance_ax.axhline(0.0, color="#94a3b8", linewidth=0.9, alpha=0.8)

    latest_cov_ax.set_title("latest gradient covariance")
    latest_cov_ax.set_xticks([0, 1], labels=["g_r1", "g_r2"])
    latest_cov_ax.set_yticks([0, 1], labels=["g_r1", "g_r2"])

    if has_covariance:
        cov_11_values = [_safe_float(row.get("gradient_cov_11")) for row in rows]
        cov_12_values = [_safe_float(row.get("gradient_cov_12")) for row in rows]
        cov_22_values = [_safe_float(row.get("gradient_cov_22")) for row in rows]
        covariance_ax.plot(step_ids, cov_11_values, "-o", color="#7c3aed", label="var(g_r1)")
        covariance_ax.plot(step_ids, cov_22_values, "-s", color="#059669", label="var(g_r2)")
        covariance_ax.plot(step_ids, cov_12_values, "-^", color="#dc2626", label="cov(g_r1,g_r2)")
        covariance_ax.legend(fontsize=8)

        latest_cov_11 = cov_11_values[-1]
        latest_cov_12 = cov_12_values[-1]
        latest_cov_22 = cov_22_values[-1]
        latest_matrix = [
            [latest_cov_11, latest_cov_12],
            [latest_cov_12, latest_cov_22],
        ]
        finite_entries = [abs(value) for row in latest_matrix for value in row if math.isfinite(value)]
        scale = max(finite_entries) if finite_entries else 1.0
        if scale <= 0.0:
            scale = 1.0
        latest_cov_ax.imshow(latest_matrix, cmap="coolwarm", vmin=-scale, vmax=scale)
        for row_index, matrix_row in enumerate(latest_matrix):
            for col_index, value in enumerate(matrix_row):
                latest_cov_ax.text(
                    col_index,
                    row_index,
                    f"{value:.2e}" if math.isfinite(value) else "nan",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="#0f172a",
                    bbox={"boxstyle": "round,pad=0.15", "fc": "#ffffff", "ec": "none", "alpha": 0.80},
                )
        latest_corr = _covariance_correlation(rows[-1])
        latest_cov_ax.text(
            0.5,
            -0.14,
            f"corr={latest_corr:.3f}" if math.isfinite(latest_corr) else "corr=nan",
            transform=latest_cov_ax.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color="#334155",
        )
    else:
        covariance_ax.text(0.5, 0.5, "covariance data not available", ha="center", va="center", fontsize=10, color="#64748b")
        latest_cov_ax.text(0.5, 0.5, "covariance data not available", ha="center", va="center", fontsize=10, color="#64748b")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    output_path = Path(args.output).resolve() if args.output else output_root / "trajectory_path.png"
    written_path = render_plot(output_root, output_path)
    print(json.dumps({"trajectory_path_png": str(written_path)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()