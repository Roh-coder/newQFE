#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm


HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(HERE, ".."))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
K_FROM_CONTINUUM = os.path.join(REPO_ROOT, "K_from_continuum")
if K_FROM_CONTINUUM not in sys.path:
    sys.path.insert(0, K_FROM_CONTINUUM)

import plot_modular_anchored_ratio_heatmaps as tau_heatmaps  # noqa: E402
from workflow_common import exact_triangular_ising_beta  # noqa: E402


DEFAULT_TARGET_GEOMETRY = (66, 66, 33, 11)
DEFAULT_TARGET_R1 = 4.702782819756
DEFAULT_TARGET_R2 = 7.353910143333
DEFAULT_OUTPUT_DIR = os.path.join(HERE, "results")
DEFAULT_MULTIPLIER_MIN = 0.6
DEFAULT_MULTIPLIER_MAX = 1.4
DEFAULT_GRID_SIZE = 61


def _sinh_rule_angles(r1: float, r2: float) -> tuple[float, float, float, float]:
    couplings = {
        "k1": float(r1),
        "k2": float(r2),
        "k3": 1.0,
    }
    beta_c = float(exact_triangular_ising_beta(couplings))
    alpha = math.atan2(1.0, math.sinh(2.0 * beta_c * couplings["k2"]))
    beta_angle = math.atan2(1.0, math.sinh(2.0 * beta_c * couplings["k1"]))
    gamma = math.atan2(1.0, math.sinh(2.0 * beta_c * couplings["k3"]))
    angle_sum = alpha + beta_angle + gamma
    if not math.isclose(angle_sum, math.pi, rel_tol=1.0e-10, abs_tol=1.0e-10):
        raise ValueError(f"sinh-rule angles do not sum to pi for r1={r1}, r2={r2}: {angle_sum}")
    return alpha, beta_angle, gamma, beta_c


def tau_eff_from_couplings(r1: float, r2: float) -> tuple[complex, float, float, float, float]:
    alpha, beta_angle, gamma, beta_c = _sinh_rule_angles(r1, r2)
    side_ratio = math.sin(alpha) / math.sin(beta_angle)
    tau_eff = complex(side_ratio * math.cos(gamma), side_ratio * math.sin(gamma))
    if not np.isfinite(tau_eff.real) or not np.isfinite(tau_eff.imag) or tau_eff.imag <= 0.0:
        raise ValueError(f"invalid tau_eff for r1={r1}, r2={r2}: {tau_eff}")
    return tau_eff, beta_c, alpha, beta_angle, gamma


def _score_row(
    *,
    r1: float,
    r2: float,
    target_values: dict[str, float],
    anchor_denominator: float,
    mp_dps: int,
    point_fraction: float,
) -> dict[str, Any]:
    tau_eff, beta_c, alpha, beta_angle, gamma = tau_eff_from_couplings(r1, r2)
    row = tau_heatmaps._score_row(
        tau=tau_eff,
        target_values=target_values,
        anchor_denominator=anchor_denominator,
        mp_dps=mp_dps,
        point_fraction=point_fraction,
    )
    row.update(
        {
            "r1": float(r1),
            "r2": float(r2),
            "tau_eff_real": float(tau_eff.real),
            "tau_eff_imag": float(tau_eff.imag),
            "beta_c": float(beta_c),
            "alpha": float(alpha),
            "beta_angle": float(beta_angle),
            "gamma": float(gamma),
        }
    )
    return row


def _build_grid(rows: list[dict[str, Any]], key: str) -> tuple[list[float], list[float], np.ndarray]:
    r1_values = sorted({float(row["r1"]) for row in rows})
    r2_values = sorted({float(row["r2"]) for row in rows})
    lookup = {(float(row["r1"]), float(row["r2"])): row for row in rows}
    grid = np.full((len(r1_values), len(r2_values)), np.nan, dtype=float)
    for i, r1 in enumerate(r1_values):
        for j, r2 in enumerate(r2_values):
            row = lookup.get((r1, r2))
            if row is None:
                continue
            value = float(row[key])
            if math.isfinite(value):
                grid[i, j] = value
    return r1_values, r2_values, grid


def _panel_summary(
    *,
    key: str,
    values: np.ndarray,
    r1_values: list[float],
    r2_values: list[float],
    signed: bool,
) -> str:
    best_index = tau_heatmaps._best_index(values, use_abs=signed)
    if best_index is None:
        return "best: n/a"
    best_r1 = float(r1_values[best_index[0]])
    best_r2 = float(r2_values[best_index[1]])
    best_value = float(values[best_index])
    label = "log resid" if signed else key.replace("_", " ")
    value_text = f"{best_value:+.4f}" if signed else f"{best_value:.4f}"
    return f"best (r1,r2)=({best_r1:.3f}, {best_r2:.3f})\n{label}={value_text}"


def _style_axis(
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
    xticks = tau_heatmaps._sparse_ticks(r1_values)
    yticks = tau_heatmaps._sparse_ticks(r2_values)
    axis.set_xticks(xticks)
    axis.set_yticks(yticks)
    axis.set_xticklabels(tau_heatmaps._format_tick_labels(xticks), fontsize=8)
    axis.set_yticklabels(tau_heatmaps._format_tick_labels(yticks), fontsize=8)
    axis.set_xlim(x0, x1)
    axis.set_ylim(y0, y1)
    axis.tick_params(axis="x", rotation=30)
    axis.set_xlabel("r1" if show_xlabel else "")
    axis.set_ylabel("r2" if show_ylabel else "")


def _plot_signed_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: str,
    title: str,
    target_r1: float,
    target_r2: float,
    point_fraction: float,
    balance_mode: str,
) -> None:
    signed_keys = tau_heatmaps._signed_metric_keys()
    aggregate_keys = tau_heatmaps._aggregate_plot_keys(balance_mode)
    plot_keys = [*signed_keys, *aggregate_keys]
    plot_titles = [
        *tau_heatmaps._signed_metric_titles(point_fraction),
        *tau_heatmaps._aggregate_plot_titles(balance_mode),
    ]

    r1_values: list[float] | None = None
    r2_values: list[float] | None = None
    grids: dict[str, np.ndarray] = {}
    signed_values: list[np.ndarray] = []
    aggregate_values: list[np.ndarray] = []
    for key in plot_keys:
        r1_values, r2_values, grid = _build_grid(rows, key)
        grids[key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            if key in signed_keys:
                signed_values.append(finite)
            else:
                aggregate_values.append(finite)
    if r1_values is None or r2_values is None:
        raise ValueError("failed to build r-grid signed heatmaps")

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

    x0, x1 = tau_heatmaps._grid_extent(r1_values)
    y0, y1 = tau_heatmaps._grid_extent(r2_values)
    extent = [x0, x1, y0, y1]

    n_panels = len(plot_keys)
    ncols = 4
    nrows = int(math.ceil(float(n_panels) / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18.5, 4.9 * nrows), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    used_axes = axes_flat[:n_panels]
    for axis in axes_flat[n_panels:]:
        axis.set_visible(False)
    signed_axes = []
    aggregate_axes = []
    signed_image = None
    aggregate_image = None

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
        axis.scatter([float(target_r1)], [float(target_r2)], marker="x", s=90, color="black", linewidths=1.8)
        best_index = tau_heatmaps._best_index(values, use_abs=signed)
        if best_index is not None:
            best_r1 = r1_values[best_index[0]]
            best_r2 = r2_values[best_index[1]]
            axis.scatter([best_r1], [best_r2], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        row_index, col_index = divmod(panel_index, ncols)
        _style_axis(
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
        axis.set_title(panel_title, fontsize=11)
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
        aggregate_label = "aggregate score" if balance_mode == "none" else f"{balance_mode}-balanced aggregate score"
        fig.colorbar(aggregate_image, ax=aggregate_axes, fraction=0.025, pad=0.025, label=aggregate_label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target (r1, r2)", fontsize=14)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _plot_positive_heatmaps(
    *,
    rows: list[dict[str, Any]],
    output_path: str,
    title: str,
    target_r1: float,
    target_r2: float,
    point_fraction: float,
    balance_mode: str,
) -> None:
    keys = tau_heatmaps._positive_plot_keys(balance_mode)
    titles = tau_heatmaps._positive_plot_titles(point_fraction, balance_mode)

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
        raise ValueError("failed to build r-grid positive heatmaps")

    vmax = float(np.max(np.concatenate(finite_values))) if finite_values else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad(color="#e5e7eb")

    x0, x1 = tau_heatmaps._grid_extent(r1_values)
    y0, y1 = tau_heatmaps._grid_extent(r2_values)
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
        axis.set_rasterization_zorder(0)
        axis.scatter([float(target_r1)], [float(target_r2)], marker="x", s=90, color="black", linewidths=1.8)
        best_index = tau_heatmaps._best_index(values, use_abs=False)
        if best_index is not None:
            best_r1 = r1_values[best_index[0]]
            best_r2 = r2_values[best_index[1]]
            axis.scatter([best_r1], [best_r2], marker="*", s=180, color="white", edgecolors="black", linewidths=0.8, zorder=4)
        row_index, col_index = divmod(panel_index, ncols)
        _style_axis(
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
        axis.set_title(panel_title, fontsize=11)
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
        colorbar_label = "positive score" if balance_mode == "none" else f"{balance_mode}-balanced positive score"
        fig.colorbar(image, ax=used_axes, fraction=0.025, pad=0.025, label=colorbar_label)
    fig.suptitle(f"{title}\nstar = per-panel best, x = target (r1, r2)", fontsize=14)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _write_table(output_path: str, rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    headers = [
        "r1",
        "r2",
        "tau_eff_real",
        "tau_eff_imag",
        "beta_c",
        "alpha",
        "beta_angle",
        "gamma",
        *[f"{observable.key}_value" for observable in tau_heatmaps.OBSERVABLES],
        *[f"{observable.key}_log_residual" for observable in tau_heatmaps.OBSERVABLES],
        "vw_diff_log_residual",
        *[f"{observable.key}_chi2" for observable in tau_heatmaps.OBSERVABLES],
        *tau_heatmaps._aggregate_metric_keys(),
        *[tau_heatmaps._scaled_metric_key(f"{observable.key}_chi2") for observable in tau_heatmaps.OBSERVABLES],
        *[tau_heatmaps._scaled_metric_key(key) for key in tau_heatmaps._aggregate_metric_keys()],
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([f"{float(row[key]):.10f}" for key in headers])


def _write_summary(output_path: str, rows: list[dict[str, Any]], *, include_scaled_metrics: bool) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    metrics = tau_heatmaps._raw_positive_metric_keys()
    if include_scaled_metrics:
        metrics.extend([tau_heatmaps._scaled_metric_key(key) for key in tau_heatmaps._raw_positive_metric_keys()])
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "best_r1", "best_r2", "tau_eff_real", "tau_eff_imag", "best_value"])
        for metric in metrics:
            best_row = min(rows, key=lambda row: float(row[metric]))
            writer.writerow(
                [
                    metric,
                    f"{float(best_row['r1']):.10f}",
                    f"{float(best_row['r2']):.10f}",
                    f"{float(best_row['tau_eff_real']):.10f}",
                    f"{float(best_row['tau_eff_imag']):.10f}",
                    f"{float(best_row[metric]):.10f}",
                ]
            )


def _write_notes(
    output_path: str,
    *,
    target_geometry: tuple[int, int, int, int],
    target_tau: complex,
    target_r1: float,
    target_r2: float,
    r1_min: float,
    r1_max: float,
    r2_min: float,
    r2_max: float,
    anchor_denominator: float,
    point_fraction: float,
    balance_mode: str,
    chi2_scales: dict[str, float],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    v_label = tau_heatmaps._direction_fraction_label("v", point_fraction)
    u_label = tau_heatmaps._direction_fraction_label("u", point_fraction)
    w_label = tau_heatmaps._direction_fraction_label("w", point_fraction)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("# Modular anchored-ratio r-grid heatmaps\n\n")
        handle.write(f"Target geometry: {target_geometry}\n\n")
        handle.write(f"Target tau: {target_tau.real:.12f} + {target_tau.imag:.12f} i\n\n")
        handle.write(f"Target couplings: r1 = {target_r1:.12f}, r2 = {target_r2:.12f}\n\n")
        handle.write(f"r1 range: [{r1_min:.6f}, {r1_max:.6f}]\n")
        handle.write(f"r2 range: [{r2_min:.6f}, {r2_max:.6f}]\n\n")
        handle.write(f"Anchor choice: nu_anchor = tau / {anchor_denominator:.3f}\n\n")
        handle.write(f"Point fraction: {str(tau_heatmaps._point_fraction_rational(point_fraction))}\n")
        handle.write(f"Aggregate balance mode: {balance_mode}\n")
        handle.write(f"Balance rule: {tau_heatmaps._balance_mode_description(balance_mode)}\n")
        handle.write("Observable set:\n")
        handle.write(f"- {v_label} anchored to nu_anchor\n")
        handle.write(f"- {u_label} anchored to nu_anchor\n")
        handle.write(f"- {w_label} anchored to nu_anchor\n")
        handle.write(f"- ({v_label})/({u_label})\n")
        handle.write(f"- ({w_label})/({u_label})\n\n")
        handle.write("Additional derived metrics:\n")
        handle.write(f"- v+w score = chi2[{v_label} anchor-ratio] + chi2[{w_label} anchor-ratio]\n")
        handle.write(f"- v-w score = (log residual[{v_label} anchor-ratio] - log residual[{w_label} anchor-ratio])^2\n\n")
        if balance_mode != "none":
            handle.write("Per-observable chi2 scales:\n")
            for observable in tau_heatmaps.OBSERVABLES:
                chi2_key = f"{observable.key}_chi2"
                handle.write(f"- {chi2_key}: {float(chi2_scales[chi2_key]):.8g}\n")
            handle.write(f"- vw_diff_score: {float(chi2_scales['vw_diff_score']):.8g}\n")
            handle.write("\n")
        handle.write("Sinh-rule mapping used:\n\n")
        handle.write("- sinh(2 beta_c k2) = cot(alpha)\n")
        handle.write("- sinh(2 beta_c k1) = cot(beta)\n")
        handle.write("- sinh(2 beta_c k3) = cot(gamma)\n")
        handle.write("- tau_eff = (sin(alpha) / sin(beta)) * exp(i gamma)\n")


def _resolve_bounds(lower: float | None, upper: float | None, truth: float, multiplier_min: float, multiplier_max: float) -> tuple[float, float]:
    lo = float(lower) if lower is not None else float(multiplier_min) * float(truth)
    hi = float(upper) if upper is not None else float(multiplier_max) * float(truth)
    if lo <= 0.0 or hi <= lo:
        raise ValueError(f"invalid bounds [{lo}, {hi}] for truth={truth}")
    return lo, hi


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot modular anchored-ratio heatmaps in r1,r2 using the exact triangular sinh-rule map to tau_eff.")
    parser.add_argument("--target-geometry", nargs=4, type=int, default=list(DEFAULT_TARGET_GEOMETRY))
    parser.add_argument("--truth-r1", type=float, default=DEFAULT_TARGET_R1)
    parser.add_argument("--truth-r2", type=float, default=DEFAULT_TARGET_R2)
    parser.add_argument("--r1-min", type=float)
    parser.add_argument("--r1-max", type=float)
    parser.add_argument("--r2-min", type=float)
    parser.add_argument("--r2-max", type=float)
    parser.add_argument("--multiplier-min", type=float, default=DEFAULT_MULTIPLIER_MIN)
    parser.add_argument("--multiplier-max", type=float, default=DEFAULT_MULTIPLIER_MAX)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    parser.add_argument("--anchor-denominator", type=float, default=tau_heatmaps.DEFAULT_ANCHOR_DENOMINATOR)
    parser.add_argument("--fraction", type=float, default=tau_heatmaps.DEFAULT_POINT_FRACTION)
    parser.add_argument("--balance-mode", choices=tau_heatmaps.BALANCE_CHOICES, default=tau_heatmaps.DEFAULT_BALANCE_MODE)
    parser.add_argument("--mp-dps", type=int, default=70)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--title-prefix", default="Acute456 small-target modular anchored-ratio via sinh rule")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    point_fraction = float(args.fraction)
    tau_heatmaps._point_fraction_rational(point_fraction)
    target_geometry = tuple(int(v) for v in args.target_geometry)
    target_tau, _ = tau_heatmaps._modular_tau(*target_geometry)
    target_tau = complex(float(target_tau.real), float(target_tau.imag))
    target_values = tau_heatmaps._anchored_observables(
        target_tau,
        anchor_denominator=float(args.anchor_denominator),
        mp_dps=int(args.mp_dps),
        point_fraction=point_fraction,
    )

    r1_min, r1_max = _resolve_bounds(args.r1_min, args.r1_max, float(args.truth_r1), float(args.multiplier_min), float(args.multiplier_max))
    r2_min, r2_max = _resolve_bounds(args.r2_min, args.r2_max, float(args.truth_r2), float(args.multiplier_min), float(args.multiplier_max))
    grid_size = max(int(args.grid_size), 3)
    r1_values = np.linspace(r1_min, r1_max, grid_size)
    r2_values = np.linspace(r2_min, r2_max, grid_size)

    rows: list[dict[str, Any]] = []
    for r1 in r1_values:
        for r2 in r2_values:
            rows.append(
                _score_row(
                    r1=float(r1),
                    r2=float(r2),
                    target_values=target_values,
                    anchor_denominator=float(args.anchor_denominator),
                    mp_dps=int(args.mp_dps),
                    point_fraction=point_fraction,
                )
            )
    chi2_scales = tau_heatmaps._add_balanced_scores(
        rows,
        balance_mode=str(args.balance_mode),
        coord_keys=("r1", "r2"),
        target_point=(float(args.truth_r1), float(args.truth_r2)),
    )

    output_dir = os.path.abspath(str(args.output_dir))
    os.makedirs(output_dir, exist_ok=True)
    base_name = f"target_Lx{target_geometry[0]}_Ly{target_geometry[1]}_Tx{target_geometry[2]}_Ty{target_geometry[3]}_sinh_rule_rgrid"
    table_path = os.path.join(output_dir, f"{base_name}.tsv")
    summary_path = os.path.join(output_dir, f"{base_name}_best_points.tsv")
    notes_path = os.path.join(output_dir, f"{base_name}_README.md")
    signed_heatmap_path = os.path.join(output_dir, f"{base_name}_signed_heatmaps.png")
    chi2_heatmap_path = os.path.join(output_dir, f"{base_name}_chi2_heatmaps.png")

    _write_table(table_path, rows)
    _write_summary(summary_path, rows, include_scaled_metrics=str(args.balance_mode) != "none")
    _write_notes(
        notes_path,
        target_geometry=target_geometry,
        target_tau=target_tau,
        target_r1=float(args.truth_r1),
        target_r2=float(args.truth_r2),
        r1_min=r1_min,
        r1_max=r1_max,
        r2_min=r2_min,
        r2_max=r2_max,
        anchor_denominator=float(args.anchor_denominator),
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
        chi2_scales=chi2_scales,
    )
    _plot_signed_heatmaps(
        rows=rows,
        output_path=signed_heatmap_path,
        title=f"{args.title_prefix} signed residual heatmaps",
        target_r1=float(args.truth_r1),
        target_r2=float(args.truth_r2),
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
    )
    _plot_positive_heatmaps(
        rows=rows,
        output_path=chi2_heatmap_path,
        title=f"{args.title_prefix} score heatmaps",
        target_r1=float(args.truth_r1),
        target_r2=float(args.truth_r2),
        point_fraction=point_fraction,
        balance_mode=str(args.balance_mode),
    )

    print(f"target_tau={target_tau.real:.12f}+{target_tau.imag:.12f}i")
    print(f"target_r1={float(args.truth_r1):.12f}")
    print(f"target_r2={float(args.truth_r2):.12f}")
    print(f"wrote {table_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {notes_path}")
    print(f"wrote {signed_heatmap_path}")
    print(f"wrote {chi2_heatmap_path}")


if __name__ == "__main__":
    main()