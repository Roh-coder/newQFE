#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

import matplotlib
import numpy as np
from mpmath import mp

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm


HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(HERE, ".."))
RMT_SCRIPTS = os.path.join(REPO_ROOT, "responsible_method_tests", "scripts")
if RMT_SCRIPTS not in sys.path:
    sys.path.insert(0, RMT_SCRIPTS)

from generate_pointwise_manifold_dataset import _ising_torus_shape, _modular_tau  # noqa: E402


DEFAULT_TARGET_GEOMETRY = (66, 66, 33, 11)
DEFAULT_OUTPUT_DIR = os.path.join(HERE, "results")
DEFAULT_TAU_REAL_HALFSPAN = 0.08
DEFAULT_TAU_IMAG_HALFSPAN = 0.20
DEFAULT_GRID_SIZE = 61
DEFAULT_ANCHOR_DENOMINATOR = 8.0
DEFAULT_POINT_FRACTION = 0.25
DEFAULT_MAX_TICKS = 7
DEFAULT_BALANCE_MODE = "none"
BALANCE_CHOICES = ("none", "mean", "max", "local")


@dataclass(frozen=True)
class Observable:
    key: str
    title: str
    nu_kind: str
    group: str


OBSERVABLES: tuple[Observable, ...] = (
    Observable(key="v4", title="log residual[v/4 anchor-ratio]", nu_kind="v4", group="anchored"),
    Observable(key="u4", title="log residual[u/4 anchor-ratio]", nu_kind="u4", group="anchored"),
    Observable(key="w4", title="log residual[w/4 anchor-ratio]", nu_kind="w4", group="anchored"),
    Observable(key="v4_over_u4", title="log residual[(v/4)/(u/4)]", nu_kind="v4_over_u4", group="ratio"),
    Observable(key="w4_over_u4", title="log residual[(w/4)/(u/4)]", nu_kind="w4_over_u4", group="ratio"),
)


def _aggregate_metric_keys() -> list[str]:
    return ["corr_score", "ratio_score", "total_score", "vw_score", "vw_diff_score"]


def _raw_positive_metric_keys() -> list[str]:
    return [*[f"{observable.key}_chi2" for observable in OBSERVABLES], *_aggregate_metric_keys()]


def _scaled_metric_key(key: str) -> str:
    return f"{key}_scaled"


def _observable_balance_title(title: str, balance_mode: str) -> str:
    if balance_mode == "mean":
        return f"{title} / grid mean"
    if balance_mode == "max":
        return f"{title} / grid max"
    if balance_mode == "local":
        return f"{title} / local ring"
    return title


def _aggregate_balance_title(title: str, balance_mode: str) -> str:
    if balance_mode == "none":
        return title
    return f"{title} ({balance_mode}-balanced)"


def _aggregate_plot_keys(balance_mode: str) -> list[str]:
    aggregate_keys = _aggregate_metric_keys()
    if balance_mode == "none":
        return aggregate_keys
    return [_scaled_metric_key(key) for key in aggregate_keys]


def _aggregate_plot_titles(balance_mode: str) -> list[str]:
    return [
        _aggregate_balance_title(title, balance_mode)
        for title in ["corr score", "ratio score", "total score", "v+w score", "v-w score"]
    ]


def _positive_plot_keys(balance_mode: str) -> list[str]:
    raw_keys = _raw_positive_metric_keys()
    if balance_mode == "none":
        return raw_keys
    return [_scaled_metric_key(key) for key in raw_keys]


def _positive_plot_titles(point_fraction: float, balance_mode: str) -> list[str]:
    observable_titles = _observable_titles(point_fraction, prefix="chi2")
    if balance_mode != "none":
        observable_titles = [_observable_balance_title(title, balance_mode) for title in observable_titles]
    return [*observable_titles, *_aggregate_plot_titles(balance_mode)]


def _balance_mode_description(balance_mode: str) -> str:
    if balance_mode == "mean":
        return "each observable chi2 map is divided by its grid mean before aggregation"
    if balance_mode == "max":
        return "each observable chi2 map is divided by its grid maximum before aggregation"
    if balance_mode == "local":
        return "each observable chi2 map is divided by the mean on the target-centered 3x3 neighbor ring before aggregation"
    return "raw observable chi2 maps are summed without rescaling"


def _signed_metric_keys() -> list[str]:
    return [*[f"{observable.key}_log_residual" for observable in OBSERVABLES], "vw_diff_log_residual"]


def _signed_metric_titles(point_fraction: float) -> list[str]:
    v_label = _direction_fraction_label("v", point_fraction)
    w_label = _direction_fraction_label("w", point_fraction)
    return [
        *_observable_titles(point_fraction, prefix="log residual"),
        f"log residual[({v_label})-({w_label})]",
    ]


def _build_grid_generic(
    rows: list[dict[str, Any]],
    key: str,
    *,
    x_key: str,
    y_key: str,
) -> tuple[list[float], list[float], np.ndarray]:
    x_values = sorted({float(row[x_key]) for row in rows})
    y_values = sorted({float(row[y_key]) for row in rows})
    lookup = {(float(row[x_key]), float(row[y_key])): row for row in rows}
    grid = np.full((len(x_values), len(y_values)), np.nan, dtype=float)
    for i, x_value in enumerate(x_values):
        for j, y_value in enumerate(y_values):
            row = lookup.get((x_value, y_value))
            if row is None:
                continue
            value = float(row[key])
            if math.isfinite(value):
                grid[i, j] = value
    return x_values, y_values, grid


def _local_ring_scale(
    grid: np.ndarray,
    x_values: list[float],
    y_values: list[float],
    *,
    target_x: float,
    target_y: float,
) -> float:
    if grid.size == 0:
        return 1.0
    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)
    ix = int(np.argmin(np.abs(x_arr - float(target_x))))
    iy = int(np.argmin(np.abs(y_arr - float(target_y))))
    x_lo = max(0, ix - 1)
    x_hi = min(len(x_values), ix + 2)
    y_lo = max(0, iy - 1)
    y_hi = min(len(y_values), iy + 2)
    subgrid = np.asarray(grid[x_lo:x_hi, y_lo:y_hi], dtype=float)
    mask = np.isfinite(subgrid)
    mask[ix - x_lo, iy - y_lo] = False
    ring_values = subgrid[mask]
    positive = ring_values[ring_values > 0.0]
    if positive.size:
        return float(np.mean(positive))
    finite = ring_values[np.isfinite(ring_values)]
    if finite.size:
        fallback = float(np.mean(np.abs(finite)))
        return fallback if fallback > 0.0 else 1.0
    return 1.0


def _grid_balance_scale(rows: list[dict[str, Any]], key: str, balance_mode: str) -> float:
    if balance_mode == "none":
        return 1.0
    values = np.asarray([float(row[key]) for row in rows if math.isfinite(float(row[key]))], dtype=float)
    if values.size == 0:
        return 1.0
    if balance_mode == "mean":
        scale = float(np.mean(values))
    elif balance_mode == "max":
        scale = float(np.max(values))
    else:
        raise ValueError(f"unsupported balance mode: {balance_mode}")
    if not math.isfinite(scale) or scale <= 0.0:
        return 1.0
    return scale


def _add_balanced_scores(
    rows: list[dict[str, Any]],
    balance_mode: str,
    *,
    coord_keys: tuple[str, str] | None = None,
    target_point: tuple[float, float] | None = None,
) -> dict[str, float]:
    scale_metric_keys = [*[f"{observable.key}_chi2" for observable in OBSERVABLES], "vw_diff_score"]
    if balance_mode == "local":
        if coord_keys is None or target_point is None:
            raise ValueError("local balance mode requires coord_keys and target_point")
        metric_scales = {}
        for metric_key in scale_metric_keys:
            x_values, y_values, grid = _build_grid_generic(rows, metric_key, x_key=coord_keys[0], y_key=coord_keys[1])
            metric_scales[metric_key] = _local_ring_scale(
                grid,
                x_values,
                y_values,
                target_x=float(target_point[0]),
                target_y=float(target_point[1]),
            )
    else:
        metric_scales = {metric_key: _grid_balance_scale(rows, metric_key, balance_mode) for metric_key in scale_metric_keys}
    for row in rows:
        corr_scaled = 0.0
        ratio_scaled = 0.0
        for observable in OBSERVABLES:
            chi2_key = f"{observable.key}_chi2"
            scaled_value = float(row[chi2_key]) / float(metric_scales[chi2_key])
            row[_scaled_metric_key(chi2_key)] = scaled_value
            if observable.group == "anchored":
                corr_scaled += scaled_value
            else:
                ratio_scaled += scaled_value
        row[_scaled_metric_key("corr_score")] = float(corr_scaled)
        row[_scaled_metric_key("ratio_score")] = float(ratio_scaled)
        row[_scaled_metric_key("total_score")] = float(corr_scaled + ratio_scaled)
        row[_scaled_metric_key("vw_score")] = float(row[_scaled_metric_key("v4_chi2")] + row[_scaled_metric_key("w4_chi2")])
        row[_scaled_metric_key("vw_diff_score")] = float(row["vw_diff_score"]) / float(metric_scales["vw_diff_score"])
    return metric_scales


def _point_fraction_rational(point_fraction: float) -> Fraction:
    fraction_value = float(point_fraction)
    if not np.isfinite(fraction_value) or fraction_value <= 0.0 or fraction_value >= 1.0:
        raise ValueError(f"point fraction must satisfy 0 < fraction < 1, got {point_fraction}")
    return Fraction(f"{fraction_value:.12g}").limit_denominator(128)


def _direction_fraction_label(direction: str, point_fraction: float) -> str:
    fraction = _point_fraction_rational(point_fraction)
    if fraction.denominator == 1:
        return f"{fraction.numerator}{direction}"
    if fraction.numerator == 1:
        return f"{direction}/{fraction.denominator}"
    return f"{fraction.numerator}{direction}/{fraction.denominator}"


def _observable_titles(point_fraction: float, *, prefix: str) -> list[str]:
    v_label = _direction_fraction_label("v", point_fraction)
    u_label = _direction_fraction_label("u", point_fraction)
    w_label = _direction_fraction_label("w", point_fraction)
    return [
        f"{prefix}[{v_label} anchor-ratio]",
        f"{prefix}[{u_label} anchor-ratio]",
        f"{prefix}[{w_label} anchor-ratio]",
        f"{prefix}[({v_label})/({u_label})]",
        f"{prefix}[({w_label})/({u_label})]",
    ]


def _sample_points(tau: complex, point_fraction: float) -> dict[str, complex]:
    fraction = float(_point_fraction_rational(point_fraction))
    return {
        "v4": complex(fraction, 0.0),
        "u4": fraction * tau,
        "w4": fraction * (1.0 + tau),
    }


def _anchor_point(tau: complex, anchor_denominator: float) -> complex:
    return tau / float(anchor_denominator)


def _q_and_theta1p0(tau: complex, mp_dps: int) -> tuple[mp.mpc, mp.mpf]:
    mp.dps = max(int(mp_dps), 30)
    tau_mp = mp.mpc(float(tau.real), float(tau.imag))
    q = mp.e ** (mp.pi * 1j * tau_mp)
    theta1p0 = mp.diff(lambda zz: mp.jtheta(1, zz, q), mp.mpf("0.0"))
    return q, theta1p0


def _shape_value(nu: complex, tau: complex, q: mp.mpc, theta1p0: mp.mpf) -> float:
    value = _ising_torus_shape(complex(float(nu.real), float(nu.imag)), complex(float(tau.real), float(tau.imag)), theta1p0, q)
    return float(value)


def _anchored_observables(tau: complex, *, anchor_denominator: float, mp_dps: int) -> dict[str, float]:
    q, theta1p0 = _q_and_theta1p0(tau, mp_dps)
def _anchored_observables(
    tau: complex,
    *,
    anchor_denominator: float,
    mp_dps: int,
    point_fraction: float = DEFAULT_POINT_FRACTION,
) -> dict[str, float]:
    q, theta1p0 = _q_and_theta1p0(tau, mp_dps)
    points = _sample_points(tau, point_fraction)
    anchor = _anchor_point(tau, anchor_denominator)
    anchor_value = _shape_value(anchor, tau, q, theta1p0)
    if not np.isfinite(anchor_value) or anchor_value <= 0.0:
        raise ValueError(f"invalid anchor modular value for tau={tau}")

    values: dict[str, float] = {}
    for key, nu in points.items():
        point_value = _shape_value(nu, tau, q, theta1p0)
        if not np.isfinite(point_value) or point_value <= 0.0:
            raise ValueError(f"invalid modular value for {key} at tau={tau}")
        values[key] = float(point_value / anchor_value)

    values["v4_over_u4"] = float(values["v4"] / values["u4"])
    values["w4_over_u4"] = float(values["w4"] / values["u4"])
    return values


def _score_row(
    *,
    tau: complex,
    target_values: dict[str, float],
    anchor_denominator: float,
    mp_dps: int,
    point_fraction: float = DEFAULT_POINT_FRACTION,
) -> dict[str, Any]:
    candidate_values = _anchored_observables(
        tau,
        anchor_denominator=anchor_denominator,
        mp_dps=mp_dps,
        point_fraction=point_fraction,
    )
    row: dict[str, Any] = {
        "tau_real": float(tau.real),
        "tau_imag": float(tau.imag),
    }
    corr_score = 0.0
    ratio_score = 0.0
    residuals: dict[str, float] = {}
    for observable in OBSERVABLES:
        candidate_value = float(candidate_values[observable.key])
        target_value = float(target_values[observable.key])
        residual = math.log(candidate_value) - math.log(target_value)
        residuals[observable.key] = residual
        row[f"{observable.key}_value"] = candidate_value
        row[f"{observable.key}_target"] = target_value
        row[f"{observable.key}_log_residual"] = residual
        row[f"{observable.key}_chi2"] = residual * residual
        if observable.group == "anchored":
            corr_score += residual * residual
        else:
            ratio_score += residual * residual
    vw_diff_log_residual = residuals["v4"] - residuals["w4"]
    row["vw_diff_log_residual"] = float(vw_diff_log_residual)
    row["corr_score"] = float(corr_score)
    row["ratio_score"] = float(ratio_score)
    row["total_score"] = float(corr_score + ratio_score)
    row["vw_score"] = float(row["v4_chi2"] + row["w4_chi2"])
    row["vw_diff_score"] = float(vw_diff_log_residual * vw_diff_log_residual)
    return row


def _grid_extent(values: list[float]) -> tuple[float, float]:
    if len(values) == 1:
        return values[0] - 0.5, values[0] + 0.5
    steps = np.diff(np.asarray(values, dtype=float))
    step = float(np.median(steps))
    return float(values[0] - 0.5 * step), float(values[-1] + 0.5 * step)


def _build_grid(rows: list[dict[str, Any]], key: str) -> tuple[list[float], list[float], np.ndarray]:
    return _build_grid_generic(rows, key, x_key="tau_real", y_key="tau_imag")


def _best_index(values: np.ndarray, *, use_abs: bool) -> tuple[int, int] | None:
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        return None
    scored = np.abs(values) if use_abs else values
    masked = np.where(finite_mask, scored, np.inf)
    flat_index = int(np.argmin(masked))
    return tuple(int(v) for v in np.unravel_index(flat_index, values.shape))


def _sparse_ticks(values: list[float], *, max_ticks: int = DEFAULT_MAX_TICKS) -> list[float]:
    if not values:
        return []
    if len(values) <= max_ticks:
        return [float(value) for value in values]
    indices = np.linspace(0, len(values) - 1, max_ticks)
    indices = sorted({int(round(index)) for index in indices})
    return [float(values[index]) for index in indices]


def _format_tick_labels(values: list[float]) -> list[str]:
    return [f"{value:.3f}" for value in values]


def _panel_summary(
    *,
    key: str,
    values: np.ndarray,
    tau_real_values: list[float],
    tau_imag_values: list[float],
    signed: bool,
) -> str:
    best_index = _best_index(values, use_abs=signed)
    if best_index is None:
        return "best: n/a"
    best_real = float(tau_real_values[best_index[0]])
    best_imag = float(tau_imag_values[best_index[1]])
    best_value = float(values[best_index])
    label = "log resid" if signed else key.replace("_", " ")
    value_text = f"{best_value:+.4f}" if signed else f"{best_value:.4f}"
    return f"best tau=({best_real:.3f}, {best_imag:.3f})\n{label}={value_text}"


def _style_heatmap_axis(
    *,
    axis: plt.Axes,
    tau_real_values: list[float],
    tau_imag_values: list[float],
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    show_xlabel: bool,
    show_ylabel: bool,
) -> None:
    xticks = _sparse_ticks(tau_real_values)
    yticks = _sparse_ticks(tau_imag_values)
    axis.set_xticks(xticks)
    axis.set_yticks(yticks)
    axis.set_xticklabels(_format_tick_labels(xticks), fontsize=8)
    axis.set_yticklabels(_format_tick_labels(yticks), fontsize=8)
    axis.set_xlim(x0, x1)
    axis.set_ylim(y0, y1)
    axis.tick_params(axis="x", rotation=30)
    axis.set_xlabel("tau_real" if show_xlabel else "")
    axis.set_ylabel("tau_imag" if show_ylabel else "")


def _write_table(output_path: str, rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    headers = [
        "tau_real",
        "tau_imag",
        *[f"{observable.key}_value" for observable in OBSERVABLES],
        *[f"{observable.key}_log_residual" for observable in OBSERVABLES],
        "vw_diff_log_residual",
        *[f"{observable.key}_chi2" for observable in OBSERVABLES],
        *_aggregate_metric_keys(),
        *[_scaled_metric_key(f"{observable.key}_chi2") for observable in OBSERVABLES],
        *[_scaled_metric_key(key) for key in _aggregate_metric_keys()],
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([f"{float(row[key]):.10f}" for key in headers])


def _write_summary(output_path: str, rows: list[dict[str, Any]], *, include_scaled_metrics: bool) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    metrics = _raw_positive_metric_keys()
    if include_scaled_metrics:
        metrics.extend([_scaled_metric_key(key) for key in _raw_positive_metric_keys()])
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "best_tau_real", "best_tau_imag", "best_value"])
        for metric in metrics:
            best_row = min(rows, key=lambda row: float(row[metric]))
            writer.writerow(
                [
                    metric,
                    f"{float(best_row['tau_real']):.10f}",
                    f"{float(best_row['tau_imag']):.10f}",
                    f"{float(best_row[metric]):.10f}",
                ]
            )


def _plot_signed_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: str,
    title: str,
    target_tau: complex,
    point_fraction: float,
    balance_mode: str,
) -> None:
    signed_keys = _signed_metric_keys()
    aggregate_keys = _aggregate_plot_keys(balance_mode)

    tau_real_values: list[float] | None = None
    tau_imag_values: list[float] | None = None
    grids: dict[str, np.ndarray] = {}
    signed_values: list[np.ndarray] = []
    aggregate_values: list[np.ndarray] = []
    for key in [*signed_keys, *aggregate_keys]:
        tau_real_values, tau_imag_values, grid = _build_grid(rows, key)
        grids[key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            if key in signed_keys:
                signed_values.append(finite)
            else:
                aggregate_values.append(finite)

    if tau_real_values is None or tau_imag_values is None:
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

    x0, x1 = _grid_extent(tau_real_values)
    y0, y1 = _grid_extent(tau_imag_values)
    extent = [x0, x1, y0, y1]

    fig, axes = plt.subplots(2, 4, figsize=(18.5, 9.8), squeeze=False, constrained_layout=True)
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
    signed_axes = []
    aggregate_axes = []
    plot_keys = [*signed_keys, *aggregate_keys]
    plot_titles = [
        *_signed_metric_titles(point_fraction),
        *_aggregate_plot_titles(balance_mode),
    ]
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
        axis.set_rasterization_zorder(0)
        if signed:
            signed_axes.append(axis)
            signed_image = image
        else:
            aggregate_axes.append(axis)
            aggregate_image = image
        axis.scatter([float(target_tau.real)], [float(target_tau.imag)], marker="x", s=90, color="black", linewidths=1.8)
        best_index = _best_index(values, use_abs=signed)
        if best_index is not None:
            best_real = tau_real_values[best_index[0]]
            best_imag = tau_imag_values[best_index[1]]
            axis.scatter([best_real], [best_imag], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        axis.set_title(panel_title, fontsize=11)
        row_index, col_index = divmod(panel_index, ncols)
        _style_heatmap_axis(
            axis=axis,
            tau_real_values=tau_real_values,
            tau_imag_values=tau_imag_values,
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
            _panel_summary(key=key, values=values, tau_real_values=tau_real_values, tau_imag_values=tau_imag_values, signed=signed),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#9ca3af"},
        )
    if signed_image is not None:
        fig.colorbar(signed_image, ax=signed_axes, fraction=0.025, pad=0.025, label="signed log residual")
    if aggregate_image is not None:
        aggregate_label = "aggregate score" if balance_mode == "none" else f"{balance_mode}-balanced aggregate score"
        fig.colorbar(aggregate_image, ax=aggregate_axes, fraction=0.025, pad=0.025, label=aggregate_label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target tau", fontsize=14)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _plot_positive_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: str,
    title: str,
    target_tau: complex,
    point_fraction: float,
    balance_mode: str,
) -> None:
    keys = _positive_plot_keys(balance_mode)
    tau_real_values: list[float] | None = None
    tau_imag_values: list[float] | None = None
    grids: dict[str, np.ndarray] = {}
    finite_values: list[np.ndarray] = []
    for key in keys:
        tau_real_values, tau_imag_values, grid = _build_grid(rows, key)
        grids[key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            finite_values.append(finite)
    if tau_real_values is None or tau_imag_values is None:
        raise ValueError("failed to build chi2 heatmap grids")
    vmax = float(np.max(np.concatenate(finite_values))) if finite_values else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad(color="#e5e7eb")
    x0, x1 = _grid_extent(tau_real_values)
    y0, y1 = _grid_extent(tau_imag_values)
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
    for panel_index, (axis, key, panel_title) in enumerate(zip(
        used_axes,
        keys,
        _positive_plot_titles(point_fraction, balance_mode),
    )):
        values = np.asarray(grids[key], dtype=float)
        image = axis.imshow(values.T, origin="lower", extent=extent, aspect="auto", cmap=cmap, norm=norm)
        axis.set_rasterization_zorder(0)
        axis.scatter([float(target_tau.real)], [float(target_tau.imag)], marker="x", s=90, color="black", linewidths=1.8)
        best_index = _best_index(values, use_abs=False)
        if best_index is not None:
            best_real = tau_real_values[best_index[0]]
            best_imag = tau_imag_values[best_index[1]]
            axis.scatter([best_real], [best_imag], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        axis.set_title(panel_title, fontsize=11)
        row_index, col_index = divmod(panel_index, ncols)
        _style_heatmap_axis(
            axis=axis,
            tau_real_values=tau_real_values,
            tau_imag_values=tau_imag_values,
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
            _panel_summary(key=key, values=values, tau_real_values=tau_real_values, tau_imag_values=tau_imag_values, signed=False),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#9ca3af"},
        )
    if image is not None:
        colorbar_label = "positive score" if balance_mode == "none" else f"{balance_mode}-balanced positive score"
        fig.colorbar(image, ax=used_axes, fraction=0.025, pad=0.025, label=colorbar_label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target tau", fontsize=14)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _write_notes(
    output_path: str,
    *,
    target_tau: complex,
    anchor_denominator: float,
    target_geometry: tuple[int, int, int, int],
    point_fraction: float,
    balance_mode: str,
    chi2_scales: dict[str, float],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    v_label = _direction_fraction_label("v", point_fraction)
    u_label = _direction_fraction_label("u", point_fraction)
    w_label = _direction_fraction_label("w", point_fraction)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("# Modular anchored-ratio heatmaps\n\n")
        handle.write(f"Target geometry: {target_geometry}\n\n")
        handle.write(f"Target tau: {target_tau.real:.12f} + {target_tau.imag:.12f} i\n\n")
        handle.write(f"Anchor choice: nu_anchor = tau / {anchor_denominator:.3f}\n\n")
        handle.write(f"Point fraction: {str(_point_fraction_rational(point_fraction))}\n\n")
        handle.write("Observable set:\n")
        handle.write(f"- {v_label} anchored to nu_anchor\n")
        handle.write(f"- {u_label} anchored to nu_anchor\n")
        handle.write(f"- {w_label} anchored to nu_anchor\n")
        handle.write(f"- ({v_label})/({u_label})\n")
        handle.write(f"- ({w_label})/({u_label})\n\n")
        handle.write("Additional derived metrics:\n")
        handle.write(f"- v+w score = chi2[{v_label} anchor-ratio] + chi2[{w_label} anchor-ratio]\n")
        handle.write(f"- v-w score = (log residual[{v_label} anchor-ratio] - log residual[{w_label} anchor-ratio])^2\n\n")
        handle.write(f"Aggregate balance mode: {balance_mode}\n")
        handle.write(f"Balance rule: {_balance_mode_description(balance_mode)}\n\n")
        if balance_mode != "none":
            handle.write("Per-observable chi2 scales:\n")
            for observable in OBSERVABLES:
                chi2_key = f"{observable.key}_chi2"
                handle.write(f"- {chi2_key}: {float(chi2_scales[chi2_key]):.8g}\n")
            handle.write(f"- vw_diff_score: {float(chi2_scales['vw_diff_score']):.8g}\n")
            handle.write("\n")
        handle.write("Signed panels use log residuals relative to the target tau. Positive panels use squared log residuals and summed sector scores.\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot tau-space anchored-ratio heatmaps from the modular Ising torus solution.")
    parser.add_argument("--target-geometry", nargs=4, type=int, default=list(DEFAULT_TARGET_GEOMETRY))
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--tau-real-halfspan", type=float, default=DEFAULT_TAU_REAL_HALFSPAN)
    parser.add_argument("--tau-imag-halfspan", type=float, default=DEFAULT_TAU_IMAG_HALFSPAN)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    parser.add_argument("--anchor-denominator", type=float, default=DEFAULT_ANCHOR_DENOMINATOR)
    parser.add_argument("--fraction", type=float, default=DEFAULT_POINT_FRACTION)
    parser.add_argument("--balance-mode", choices=BALANCE_CHOICES, default=DEFAULT_BALANCE_MODE)
    parser.add_argument("--mp-dps", type=int, default=70)
    parser.add_argument("--title-prefix", default="Acute456 modular anchored-ratio")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    point_fraction = float(args.fraction)
    _point_fraction_rational(point_fraction)
    target_geometry = tuple(int(v) for v in args.target_geometry)
    target_tau, _ = _modular_tau(*target_geometry)
    target_tau = complex(float(target_tau.real), float(target_tau.imag))
    output_dir = os.path.abspath(str(args.output_dir))
    os.makedirs(output_dir, exist_ok=True)

    target_values = _anchored_observables(
        target_tau,
        anchor_denominator=float(args.anchor_denominator),
        mp_dps=int(args.mp_dps),
        point_fraction=point_fraction,
    )

    tau_real_values = np.linspace(float(target_tau.real) - float(args.tau_real_halfspan), float(target_tau.real) + float(args.tau_real_halfspan), max(int(args.grid_size), 3))
    tau_imag_values = np.linspace(max(1.0e-3, float(target_tau.imag) - float(args.tau_imag_halfspan)), float(target_tau.imag) + float(args.tau_imag_halfspan), max(int(args.grid_size), 3))
    rows: list[dict[str, Any]] = []
    for tau_real in tau_real_values:
        for tau_imag in tau_imag_values:
            tau = complex(float(tau_real), float(tau_imag))
            rows.append(
                _score_row(
                    tau=tau,
                    target_values=target_values,
                    anchor_denominator=float(args.anchor_denominator),
                    mp_dps=int(args.mp_dps),
                    point_fraction=point_fraction,
                )
            )
    chi2_scales = _add_balanced_scores(
        rows,
        balance_mode=str(args.balance_mode),
        coord_keys=("tau_real", "tau_imag"),
        target_point=(float(target_tau.real), float(target_tau.imag)),
    )

    base_name = f"target_Lx{target_geometry[0]}_Ly{target_geometry[1]}_Tx{target_geometry[2]}_Ty{target_geometry[3]}"
    table_path = os.path.join(output_dir, f"{base_name}_tau_landscape.tsv")
    summary_path = os.path.join(output_dir, f"{base_name}_best_points.tsv")
    notes_path = os.path.join(output_dir, f"{base_name}_README.md")
    signed_heatmap_path = os.path.join(output_dir, f"{base_name}_signed_heatmaps.png")
    chi2_heatmap_path = os.path.join(output_dir, f"{base_name}_chi2_heatmaps.png")

    _write_table(table_path, rows)
    _write_summary(summary_path, rows, include_scaled_metrics=str(args.balance_mode) != "none")
    _write_notes(
        notes_path,
        target_tau=target_tau,
        anchor_denominator=float(args.anchor_denominator),
        target_geometry=target_geometry,
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
        chi2_scales=chi2_scales,
    )
    _plot_signed_heatmaps(
        rows=rows,
        output_path=signed_heatmap_path,
        title=f"{args.title_prefix} signed residual heatmaps",
        target_tau=target_tau,
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
    )
    _plot_positive_heatmaps(
        rows=rows,
        output_path=chi2_heatmap_path,
        title=f"{args.title_prefix} score heatmaps",
        target_tau=target_tau,
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
    )

    print(f"target_tau={target_tau.real:.12f}+{target_tau.imag:.12f}i")
    print(f"wrote {table_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {notes_path}")
    print(f"wrote {signed_heatmap_path}")
    print(f"wrote {chi2_heatmap_path}")


if __name__ == "__main__":
    main()