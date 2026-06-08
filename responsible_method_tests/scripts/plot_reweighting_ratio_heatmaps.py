#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_raw_reweighting_fss import _bundle_path, _infer_sizes, _ratio_observables  # noqa: E402
from test_geometry_match_grid_interpolation import (  # noqa: E402
    _fit_blind_power_model,
    _parse_selected_bundle,
    _predict_at_target_x,
    _selected_specs_for_size,
    _token,
)


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent

BALANCE_CHOICES = ("none", "mean", "max")


@dataclass(frozen=True)
class Observable:
    key: str
    title: str
    group: str


OBSERVABLES: tuple[Observable, ...] = (
    Observable(key="mid_v_over_q_w", title="log residual[mid_v / q_w]", group="midpoint"),
    Observable(key="mid_u_over_q_w", title="log residual[mid_u / q_w]", group="midpoint"),
    Observable(key="mid_w_over_q_w", title="log residual[mid_w / q_w]", group="midpoint"),
    Observable(key="q_v_over_q_w", title="log residual[q_v / q_w]", group="quarter"),
    Observable(key="q_u_over_q_w", title="log residual[q_u / q_w]", group="quarter"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot modular-style heatmaps over the reweighting r1,r2 grid using q_w-normalized fitted raw ratios."
    )
    parser.add_argument(
        "--output-root",
        default=str(
            RESPONSIBLE_ROOT / "reweighting" / "iso111_grid5x5_sizes32_28_24_20_16_12_8_4_hi100k_20260606"
        ),
    )
    parser.add_argument(
        "--target-root",
        default=str(RESPONSIBLE_ROOT / "reweighting" / "iso111_target_L64_hi100k_20260606"),
    )
    parser.add_argument("--target-size", type=int, default=64)
    parser.add_argument("--target-r1", type=float, default=1.0)
    parser.add_argument("--target-r2", type=float, default=1.0)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for heatmap outputs. Defaults to <output-root>/heatmaps_target_L64_qw_norm.",
    )
    parser.add_argument("--balance-mode", choices=BALANCE_CHOICES, default="none")
    parser.add_argument("--title-prefix", default="Iso111 q_w-normalized raw-ratio")
    return parser.parse_args()


def _aggregate_metric_keys() -> list[str]:
    return ["midpoint_score", "quarter_score", "total_score"]


def _scaled_metric_key(key: str) -> str:
    return f"{key}_scaled"


def _signed_metric_keys() -> list[str]:
    return [f"{observable.key}_log_residual" for observable in OBSERVABLES]


def _signed_metric_titles() -> list[str]:
    return [observable.title for observable in OBSERVABLES]


def _positive_metric_keys(balance_mode: str) -> list[str]:
    raw_keys = [*[f"{observable.key}_chi2" for observable in OBSERVABLES], *_aggregate_metric_keys()]
    if balance_mode == "none":
        return raw_keys
    return [_scaled_metric_key(key) for key in raw_keys]


def _aggregate_plot_keys(balance_mode: str) -> list[str]:
    keys = _aggregate_metric_keys()
    if balance_mode == "none":
        return keys
    return [_scaled_metric_key(key) for key in keys]


def _aggregate_plot_titles(balance_mode: str) -> list[str]:
    titles = ["midpoint score", "quarter score", "total score"]
    if balance_mode == "none":
        return titles
    return [f"{title} ({balance_mode}-balanced)" for title in titles]


def _positive_metric_titles(balance_mode: str) -> list[str]:
    titles = [
        *[observable.title.replace("log residual", "chi2") for observable in OBSERVABLES],
        *_aggregate_plot_titles(balance_mode),
    ]
    if balance_mode == "none":
        return titles
    adjusted: list[str] = []
    for index, title in enumerate(titles):
        if index < len(OBSERVABLES):
            adjusted.append(f"{title} / grid {balance_mode}")
        else:
            adjusted.append(title)
    return adjusted


def _grid_extent(values: list[float]) -> tuple[float, float]:
    if len(values) == 1:
        return values[0] - 0.5, values[0] + 0.5
    steps = np.diff(np.asarray(values, dtype=float))
    step = float(np.median(steps))
    return float(values[0] - 0.5 * step), float(values[-1] + 0.5 * step)


def _build_grid(rows: list[dict[str, Any]], key: str) -> tuple[list[float], list[float], np.ndarray]:
    r1_values = sorted({float(row["r1"]) for row in rows})
    r2_values = sorted({float(row["r2"]) for row in rows})
    lookup = {(float(row["r1"]), float(row["r2"])): row for row in rows}
    grid = np.full((len(r1_values), len(r2_values)), np.nan, dtype=float)
    for i, r1_value in enumerate(r1_values):
        for j, r2_value in enumerate(r2_values):
            row = lookup.get((r1_value, r2_value))
            if row is None:
                continue
            value = float(row[key])
            if math.isfinite(value):
                grid[i, j] = value
    return r1_values, r2_values, grid


def _best_index(values: np.ndarray, *, use_abs: bool) -> tuple[int, int] | None:
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        return None
    scored = np.abs(values) if use_abs else values
    masked = np.where(finite_mask, scored, np.inf)
    flat_index = int(np.argmin(masked))
    return tuple(int(value) for value in np.unravel_index(flat_index, values.shape))


def _sparse_ticks(values: list[float], *, max_ticks: int = 6) -> list[float]:
    if len(values) <= max_ticks:
        return [float(value) for value in values]
    indices = np.linspace(0, len(values) - 1, max_ticks)
    indices = sorted({int(round(index)) for index in indices})
    return [float(values[index]) for index in indices]


def _format_tick_labels(values: list[float]) -> list[str]:
    return [f"{value:.2f}" for value in values]


def _style_heatmap_axis(
    *,
    axis: plt.Axes,
    r1_values: list[float],
    r2_values: list[float],
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    show_xlabel: bool,
    show_ylabel: bool,
) -> None:
    xticks = _sparse_ticks(r1_values)
    yticks = _sparse_ticks(r2_values)
    axis.set_xticks(xticks)
    axis.set_yticks(yticks)
    axis.set_xticklabels(_format_tick_labels(xticks), fontsize=8)
    axis.set_yticklabels(_format_tick_labels(yticks), fontsize=8)
    axis.set_xlim(x0, x1)
    axis.set_ylim(y0, y1)
    axis.set_xlabel("r1" if show_xlabel else "")
    axis.set_ylabel("r2" if show_ylabel else "")


def _panel_summary(
    *,
    key: str,
    values: np.ndarray,
    r1_values: list[float],
    r2_values: list[float],
    signed: bool,
) -> str:
    best_index = _best_index(values, use_abs=signed)
    if best_index is None:
        return "best: n/a"
    best_r1 = float(r1_values[best_index[0]])
    best_r2 = float(r2_values[best_index[1]])
    best_value = float(values[best_index])
    label = "log resid" if signed else key.replace("_", " ")
    value_text = f"{best_value:+.4f}" if signed else f"{best_value:.4f}"
    return f"best r=({best_r1:.2f}, {best_r2:.2f})\n{label}={value_text}"


def _grid_balance_scale(rows: list[dict[str, Any]], key: str, balance_mode: str) -> float:
    if balance_mode == "none":
        return 1.0
    values = np.asarray([float(row[key]) for row in rows if math.isfinite(float(row[key]))], dtype=float)
    if values.size == 0:
        return 1.0
    scale = float(np.mean(values)) if balance_mode == "mean" else float(np.max(values))
    if not math.isfinite(scale) or scale <= 0.0:
        return 1.0
    return scale


def _add_balanced_scores(rows: list[dict[str, Any]], balance_mode: str) -> dict[str, float]:
    scale_keys = [f"{observable.key}_chi2" for observable in OBSERVABLES]
    scales = {key: _grid_balance_scale(rows, key, balance_mode) for key in scale_keys}
    for row in rows:
        midpoint_scaled = 0.0
        quarter_scaled = 0.0
        for observable in OBSERVABLES:
            chi2_key = f"{observable.key}_chi2"
            scaled_value = float(row[chi2_key]) / float(scales[chi2_key])
            row[_scaled_metric_key(chi2_key)] = scaled_value
            if observable.group == "midpoint":
                midpoint_scaled += scaled_value
            else:
                quarter_scaled += scaled_value
        row[_scaled_metric_key("midpoint_score")] = float(midpoint_scaled)
        row[_scaled_metric_key("quarter_score")] = float(quarter_scaled)
        row[_scaled_metric_key("total_score")] = float(midpoint_scaled + quarter_scaled)
    return scales


def _load_summary(output_root: Path) -> dict[str, Any]:
    return json.loads((output_root / "summary.json").read_text(encoding="utf-8"))


def _family_rows(
    *,
    output_root: Path,
    point_tag: str,
    sizes: list[int],
    labels: list[str],
) -> dict[int, dict[str, dict[str, float]]]:
    rows: dict[int, dict[str, dict[str, float]]] = {}
    for size in sizes:
        payload = _parse_selected_bundle(_bundle_path(output_root, int(size), point_tag), labels)
        rows[int(size)] = _ratio_observables(payload)
    return rows


def _score_row(
    *,
    r1: float,
    r2: float,
    point_tag: str,
    family_rows: dict[int, dict[str, dict[str, float]]],
    fit_sizes: list[int],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
) -> dict[str, Any]:
    x_values = np.asarray([1.0 / float(size) for size in fit_sizes], dtype=float)
    row: dict[str, Any] = {"r1": float(r1), "r2": float(r2), "point_tag": str(point_tag)}
    midpoint_score = 0.0
    quarter_score = 0.0
    for observable in OBSERVABLES:
        key = str(observable.key)
        y_values = np.asarray([float(family_rows[size][key]["value"]) for size in fit_sizes], dtype=float)
        sigma_values = np.asarray([float(family_rows[size][key]["sigma"]) for size in fit_sizes], dtype=float)
        fit_payload = _fit_blind_power_model(x_values, y_values, sigma_values)
        pred_value, pred_sigma = _predict_at_target_x(fit_payload, float(target_x))
        target_value = float(target_rows[key]["value"])
        if pred_value <= 0.0 or target_value <= 0.0 or not math.isfinite(pred_value) or not math.isfinite(target_value):
            residual = float("nan")
            chi2 = float("nan")
        else:
            residual = float(math.log(pred_value) - math.log(target_value))
            chi2 = float(residual * residual)
        row[f"{key}_pred"] = float(pred_value)
        row[f"{key}_pred_sigma"] = float(pred_sigma)
        row[f"{key}_target"] = float(target_value)
        row[f"{key}_log_residual"] = float(residual)
        row[f"{key}_chi2"] = float(chi2)
        if math.isfinite(chi2):
            if observable.group == "midpoint":
                midpoint_score += chi2
            else:
                quarter_score += chi2
    row["midpoint_score"] = float(midpoint_score)
    row["quarter_score"] = float(quarter_score)
    row["total_score"] = float(midpoint_score + quarter_score)
    return row


def _write_table(output_path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "r1",
        "r2",
        *[f"{observable.key}_pred" for observable in OBSERVABLES],
        *[f"{observable.key}_target" for observable in OBSERVABLES],
        *[f"{observable.key}_log_residual" for observable in OBSERVABLES],
        *[f"{observable.key}_chi2" for observable in OBSERVABLES],
        *(_aggregate_metric_keys()),
        *[_scaled_metric_key(f"{observable.key}_chi2") for observable in OBSERVABLES],
        *[_scaled_metric_key(key) for key in _aggregate_metric_keys()],
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([f"{float(row[key]):.10f}" for key in headers])


def _write_summary(output_path: Path, rows: list[dict[str, Any]], balance_mode: str) -> None:
    metric_keys = _positive_metric_keys(balance_mode)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "best_r1", "best_r2", "best_value"])
        for key in metric_keys:
            best_row = min(rows, key=lambda row: float(row[key]))
            writer.writerow([
                key,
                f"{float(best_row['r1']):.10f}",
                f"{float(best_row['r2']):.10f}",
                f"{float(best_row[key]):.10f}",
            ])


def _write_notes(
    output_path: Path,
    *,
    output_root: Path,
    target_root: Path,
    target_size: int,
    target_point: tuple[float, float],
    balance_mode: str,
    chi2_scales: dict[str, float],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("# Reweighting q_w-Normalized Heatmaps\n\n")
        handle.write(f"Main raw family root: {output_root}\n\n")
        handle.write(f"Target raw root: {target_root}\n\n")
        handle.write(f"Target holdout: L={int(target_size)}, (r1,r2)=({float(target_point[0]):.3f},{float(target_point[1]):.3f})\n\n")
        handle.write("Observable set:\n")
        handle.write("- midpoint ratios normalized by q_w: mid_v/q_w, mid_u/q_w, mid_w/q_w\n")
        handle.write("- quarter ratios normalized by q_w: q_v/q_w, q_u/q_w\n\n")
        handle.write("Each grid cell is scored by fitting the raw finite-size ladder with A + B x^omega at x = 1/L, then comparing the predicted value at the target x = 1/64 to the direct L64 target using signed log residuals and chi2 = residual^2.\n\n")
        handle.write("Important caveat: the truly unused point (L/4, 3L/4) is not present in the current saved bundles. These heatmaps therefore use q_w = (L/4, L/4) as the currently available interior normalization point.\n\n")
        handle.write(f"Balance mode: {balance_mode}\n")
        if balance_mode != "none":
            handle.write("Per-observable chi2 scales:\n")
            for observable in OBSERVABLES:
                key = f"{observable.key}_chi2"
                handle.write(f"- {key}: {float(chi2_scales[key]):.8g}\n")


def _plot_signed_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: Path,
    title: str,
    target_point: tuple[float, float],
    balance_mode: str,
) -> None:
    signed_keys = _signed_metric_keys()
    aggregate_keys = _aggregate_plot_keys(balance_mode)

    r1_values: list[float] | None = None
    r2_values: list[float] | None = None
    grids: dict[str, np.ndarray] = {}
    signed_values: list[np.ndarray] = []
    aggregate_values: list[np.ndarray] = []
    for key in [*signed_keys, *aggregate_keys]:
        r1_values, r2_values, grid = _build_grid(rows, key)
        grids[key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            if key in signed_keys:
                signed_values.append(finite)
            else:
                aggregate_values.append(finite)
    if r1_values is None or r2_values is None:
        raise ValueError("failed to build signed heatmap grids")

    signed_limit = float(np.max(np.abs(np.concatenate(signed_values)))) if signed_values else 1.0
    aggregate_limit = float(np.max(np.concatenate(aggregate_values))) if aggregate_values else 1.0
    if signed_limit <= 0.0:
        signed_limit = 1.0
    if aggregate_limit <= 0.0:
        aggregate_limit = 1.0

    signed_norm = TwoSlopeNorm(vcenter=0.0, vmin=-signed_limit, vmax=signed_limit)
    aggregate_norm = Normalize(vmin=0.0, vmax=aggregate_limit)
    signed_cmap = plt.get_cmap("RdBu_r").copy()
    positive_cmap = plt.get_cmap("viridis_r").copy()
    signed_cmap.set_bad(color="#e5e7eb")
    positive_cmap.set_bad(color="#e5e7eb")

    x0, x1 = _grid_extent(r1_values)
    y0, y1 = _grid_extent(r2_values)
    extent = [x0, x1, y0, y1]

    plot_keys = [*signed_keys, *aggregate_keys]
    plot_titles = [*_signed_metric_titles(), *_aggregate_plot_titles(balance_mode)]
    n_panels = len(plot_keys)
    ncols = 4
    nrows = int(math.ceil(float(n_panels) / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18.5, 4.9 * nrows), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    used_axes = axes_flat[:n_panels]
    for axis in axes_flat[n_panels:]:
        axis.set_visible(False)
    signed_image = None
    aggregate_image = None
    signed_axes: list[plt.Axes] = []
    aggregate_axes: list[plt.Axes] = []
    for panel_index, (axis, key, panel_title) in enumerate(zip(used_axes, plot_keys, plot_titles)):
        values = np.asarray(grids[key], dtype=float)
        signed = key in signed_keys
        image = axis.imshow(
            values.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap=(signed_cmap if signed else positive_cmap),
            norm=(signed_norm if signed else aggregate_norm),
        )
        if signed:
            signed_image = image
            signed_axes.append(axis)
        else:
            aggregate_image = image
            aggregate_axes.append(axis)
        axis.scatter([float(target_point[0])], [float(target_point[1])], marker="x", s=90, color="black", linewidths=1.8)
        best_index = _best_index(values, use_abs=signed)
        if best_index is not None:
            best_r1 = r1_values[best_index[0]]
            best_r2 = r2_values[best_index[1]]
            axis.scatter([best_r1], [best_r2], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        axis.set_title(panel_title, fontsize=11)
        row_index, col_index = divmod(panel_index, ncols)
        _style_heatmap_axis(
            axis=axis,
            r1_values=r1_values,
            r2_values=r2_values,
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            show_xlabel=(row_index == nrows - 1),
            show_ylabel=(col_index == 0),
        )
        axis.text(
            0.02,
            0.98,
            _panel_summary(key=key, values=values, r1_values=r1_values, r2_values=r2_values, signed=signed),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#9ca3af"},
        )
    if signed_image is not None:
        fig.colorbar(signed_image, ax=signed_axes, fraction=0.025, pad=0.025, label="signed log residual")
    if aggregate_image is not None:
        label = "aggregate score" if balance_mode == "none" else f"{balance_mode}-balanced aggregate score"
        fig.colorbar(aggregate_image, ax=aggregate_axes, fraction=0.025, pad=0.025, label=label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target (1,1)", fontsize=14)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _plot_positive_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: Path,
    title: str,
    target_point: tuple[float, float],
    balance_mode: str,
) -> None:
    keys = _positive_metric_keys(balance_mode)
    titles = _positive_metric_titles(balance_mode)
    r1_values: list[float] | None = None
    r2_values: list[float] | None = None
    grids: dict[str, np.ndarray] = {}
    finite_values: list[np.ndarray] = []
    for key in keys:
        r1_values, r2_values, grid = _build_grid(rows, key)
        grids[key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            finite_values.append(finite)
    if r1_values is None or r2_values is None:
        raise ValueError("failed to build positive heatmap grids")
    vmax = float(np.max(np.concatenate(finite_values))) if finite_values else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad(color="#e5e7eb")
    x0, x1 = _grid_extent(r1_values)
    y0, y1 = _grid_extent(r2_values)
    extent = [x0, x1, y0, y1]

    n_panels = len(keys)
    ncols = 4
    nrows = int(math.ceil(float(n_panels) / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18.5, 4.9 * nrows), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    used_axes = axes_flat[:n_panels]
    for axis in axes_flat[n_panels:]:
        axis.set_visible(False)
    image = None
    for panel_index, (axis, key, panel_title) in enumerate(zip(used_axes, keys, titles)):
        values = np.asarray(grids[key], dtype=float)
        image = axis.imshow(values.T, origin="lower", extent=extent, aspect="auto", cmap=cmap, norm=norm)
        axis.scatter([float(target_point[0])], [float(target_point[1])], marker="x", s=90, color="black", linewidths=1.8)
        best_index = _best_index(values, use_abs=False)
        if best_index is not None:
            best_r1 = r1_values[best_index[0]]
            best_r2 = r2_values[best_index[1]]
            axis.scatter([best_r1], [best_r2], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        axis.set_title(panel_title, fontsize=11)
        row_index, col_index = divmod(panel_index, ncols)
        _style_heatmap_axis(
            axis=axis,
            r1_values=r1_values,
            r2_values=r2_values,
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            show_xlabel=(row_index == nrows - 1),
            show_ylabel=(col_index == 0),
        )
        axis.text(
            0.02,
            0.98,
            _panel_summary(key=key, values=values, r1_values=r1_values, r2_values=r2_values, signed=False),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#9ca3af"},
        )
    if image is not None:
        label = "positive score" if balance_mode == "none" else f"{balance_mode}-balanced positive score"
        fig.colorbar(image, ax=used_axes, fraction=0.025, pad=0.025, label=label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target (1,1)", fontsize=14)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    target_root = Path(args.target_root).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else output_root / "heatmaps_target_L64_qw_norm"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = _load_summary(output_root)
    fit_sizes = [int(value) for value in summary["fit_sizes"]]
    labels = [label for label, _ in _selected_specs_for_size(fit_sizes[0])]
    target_tag = f"r1_{_token(float(args.target_r1))}__r2_{_token(float(args.target_r2))}"
    target_payload = _parse_selected_bundle(_bundle_path(target_root, int(args.target_size), target_tag), labels)
    target_rows = _ratio_observables(target_payload)
    target_x = 1.0 / float(args.target_size)

    rows: list[dict[str, Any]] = []
    for point_row in list(summary["rows"]):
        point_tag = str(point_row["point_tag"])
        family_rows = _family_rows(output_root=output_root, point_tag=point_tag, sizes=fit_sizes, labels=labels)
        rows.append(
            _score_row(
                r1=float(point_row["r1"]),
                r2=float(point_row["r2"]),
                point_tag=point_tag,
                family_rows=family_rows,
                fit_sizes=fit_sizes,
                target_rows=target_rows,
                target_x=target_x,
            )
        )

    chi2_scales = _add_balanced_scores(rows, str(args.balance_mode))
    base_name = f"target_L{int(args.target_size)}_qw_norm"
    table_path = output_dir / f"{base_name}_landscape.tsv"
    summary_path = output_dir / f"{base_name}_best_points.tsv"
    notes_path = output_dir / f"{base_name}_README.md"
    signed_path = output_dir / f"{base_name}_signed_heatmaps.png"
    positive_path = output_dir / f"{base_name}_score_heatmaps.png"

    _write_table(table_path, rows)
    _write_summary(summary_path, rows, str(args.balance_mode))
    _write_notes(
        notes_path,
        output_root=output_root,
        target_root=target_root,
        target_size=int(args.target_size),
        target_point=(float(args.target_r1), float(args.target_r2)),
        balance_mode=str(args.balance_mode),
        chi2_scales=chi2_scales,
    )
    _plot_signed_heatmaps(
        rows=rows,
        output_path=signed_path,
        title=f"{args.title_prefix} signed heatmaps",
        target_point=(float(args.target_r1), float(args.target_r2)),
        balance_mode=str(args.balance_mode),
    )
    _plot_positive_heatmaps(
        rows=rows,
        output_path=positive_path,
        title=f"{args.title_prefix} score heatmaps",
        target_point=(float(args.target_r1), float(args.target_r2)),
        balance_mode=str(args.balance_mode),
    )

    print(f"wrote {table_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {notes_path}")
    print(f"wrote {signed_path}")
    print(f"wrote {positive_path}")


if __name__ == "__main__":
    main()