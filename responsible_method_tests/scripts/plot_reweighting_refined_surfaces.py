#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_raw_reweighting_fss import _bundle_path, _ratio_observables  # noqa: E402
from plot_reweighting_ratio_heatmaps import (  # noqa: E402
    BALANCE_CHOICES,
    OBSERVABLES,
    _add_balanced_scores,
    _aggregate_plot_keys,
    _aggregate_plot_titles,
    _best_index,
    _build_grid,
    _grid_extent,
    _plot_positive_heatmaps,
    _plot_signed_heatmaps,
    _score_row,
    _style_heatmap_axis,
    _write_summary,
    _write_table,
)
from test_geometry_match_grid_interpolation import (  # noqa: E402
    CouplingPoint,
    _log_weights,
    _parse_selected_bundle,
    _ratio_stat,
    _selected_specs_for_size,
    _token,
)


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a 2x-refined q_w-normalized heatmap surface by reweighting local base-grid donors "
            "to inserted midpoint couplings, then render base-grid gradient fields from direct local reweight stencils."
        )
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
        help=(
            "Directory for refined heatmap outputs. Defaults to "
            "<output-root>/heatmaps_target_L64_qw_norm_refined2x_neighbor_interp."
        ),
    )
    parser.add_argument(
        "--gradient-dir",
        default=None,
        help=(
            "Directory for gradient-field outputs. Defaults to "
            "<output-root>/gradient_fields_target_L64_qw_norm_local_reweight."
        ),
    )
    parser.add_argument(
        "--gradient-step",
        type=float,
        default=None,
        help="Local reweight step for central-difference score gradients. Defaults to one quarter of the base-grid spacing.",
    )
    parser.add_argument("--balance-mode", choices=BALANCE_CHOICES, default="none")
    parser.add_argument("--title-prefix", default="Iso111 q_w-normalized raw-ratio")
    return parser.parse_args()


def _load_summary(output_root: Path) -> dict[str, Any]:
    return json.loads((output_root / "summary.json").read_text(encoding="utf-8"))


def _base_grid_values(summary: dict[str, Any]) -> list[float]:
    values = sorted({float(value) for value in summary["grid_values"]})
    if len(values) < 2:
        raise ValueError("need at least a 2x2 base grid to refine the heatmap surface")
    return values


def _refined_values(base_values: list[float], refine_factor: int = 2) -> list[float]:
    if refine_factor != 2:
        raise ValueError(f"this script currently supports only refine_factor=2, got {refine_factor}")
    refined: list[float] = []
    for index, value in enumerate(base_values[:-1]):
        next_value = float(base_values[index + 1])
        refined.append(float(value))
        refined.append(0.5 * (float(value) + next_value))
    refined.append(float(base_values[-1]))
    return refined


def _axis_bracket(value: float, base_values: list[float], *, tol: float = 1.0e-12) -> tuple[float, float, float]:
    for base_value in base_values:
        if abs(float(base_value) - float(value)) <= tol:
            return float(base_value), float(base_value), 0.0
    for lower, upper in zip(base_values[:-1], base_values[1:]):
        if float(lower) - tol <= float(value) <= float(upper) + tol:
            if abs(float(upper) - float(lower)) <= tol:
                return float(lower), float(upper), 0.0
            fraction = (float(value) - float(lower)) / (float(upper) - float(lower))
            return float(lower), float(upper), float(min(max(fraction, 0.0), 1.0))
    raise ValueError(f"value {value:.8f} lies outside the base grid range {base_values[0]:.8f}..{base_values[-1]:.8f}")


def _neighbor_weights(point: CouplingPoint, base_values: list[float]) -> list[tuple[CouplingPoint, float]]:
    r1_lo, r1_hi, tx = _axis_bracket(float(point.r1), base_values)
    r2_lo, r2_hi, ty = _axis_bracket(float(point.r2), base_values)

    if r1_lo == r1_hi and r2_lo == r2_hi:
        return [(CouplingPoint(float(r1_lo), float(r2_lo)), 1.0)]
    if r1_lo == r1_hi:
        return [
            (CouplingPoint(float(r1_lo), float(r2_lo)), 1.0 - float(ty)),
            (CouplingPoint(float(r1_lo), float(r2_hi)), float(ty)),
        ]
    if r2_lo == r2_hi:
        return [
            (CouplingPoint(float(r1_lo), float(r2_lo)), 1.0 - float(tx)),
            (CouplingPoint(float(r1_hi), float(r2_lo)), float(tx)),
        ]

    return [
        (CouplingPoint(float(r1_lo), float(r2_lo)), (1.0 - float(tx)) * (1.0 - float(ty))),
        (CouplingPoint(float(r1_hi), float(r2_lo)), float(tx) * (1.0 - float(ty))),
        (CouplingPoint(float(r1_lo), float(r2_hi)), (1.0 - float(tx)) * float(ty)),
        (CouplingPoint(float(r1_hi), float(r2_hi)), float(tx) * float(ty)),
    ]


def _is_base_point(point: CouplingPoint, base_values: list[float], *, tol: float = 1.0e-12) -> bool:
    return any(abs(float(point.r1) - float(value)) <= tol for value in base_values) and any(
        abs(float(point.r2) - float(value)) <= tol for value in base_values
    )


def _observable_ratio_specs() -> list[tuple[str, str, str]]:
    return [(str(observable.key), str(observable.key).replace("_over_q_w", ""), "q_w") for observable in OBSERVABLES]


def _reweighted_ratio_observables(payload: dict[str, Any], target_point: CouplingPoint) -> dict[str, dict[str, float]]:
    log_w = _log_weights(payload, target_point)
    corr = dict(payload["corr"])
    mag = np.asarray(payload["mag"], dtype=float)
    rows: dict[str, dict[str, float]] = {}
    for key, numerator, denominator in _observable_ratio_specs():
        rows[key] = _ratio_stat(
            np.asarray(corr[numerator], dtype=float),
            np.asarray(corr[denominator], dtype=float),
            mag,
            log_w=log_w,
        )
    return rows


def _blend_stats(weighted_rows: list[tuple[float, dict[str, dict[str, float]]]]) -> dict[str, dict[str, float]]:
    blended: dict[str, dict[str, float]] = {}
    for key, _, _ in _observable_ratio_specs():
        weights = np.asarray([float(weight) for weight, _ in weighted_rows], dtype=float)
        values = np.asarray([float(rows[key]["value"]) for _, rows in weighted_rows], dtype=float)
        mean_value = float(np.sum(weights * values))
        spread = float(np.sqrt(np.sum(weights * np.square(values - mean_value))))

        sigma_terms: list[float] = []
        neff_values: list[float] = []
        for weight, rows in weighted_rows:
            sigma = float(rows[key].get("sigma", float("nan")))
            if math.isfinite(sigma):
                sigma_terms.append((float(weight) * sigma) ** 2)
            neff = float(rows[key].get("n_eff_fraction", float("nan")))
            if math.isfinite(neff):
                neff_values.append(neff)
        sigma_from_jk = math.sqrt(sum(sigma_terms)) if sigma_terms else float("nan")
        if math.isfinite(sigma_from_jk):
            sigma_value = float(math.sqrt(sigma_from_jk**2 + spread**2))
        else:
            sigma_value = float(spread if spread > 0.0 else float("nan"))
        blended[key] = {
            "value": float(mean_value),
            "sigma": float(sigma_value),
            "n_eff_fraction": float(min(neff_values)) if neff_values else float("nan"),
        }
    return blended


def _load_payloads(output_root: Path, fit_sizes: list[int], base_points: list[CouplingPoint], labels: list[str]) -> dict[tuple[int, str], dict[str, Any]]:
    payloads: dict[tuple[int, str], dict[str, Any]] = {}
    for size in fit_sizes:
        for point in base_points:
            bundle_path = _bundle_path(output_root, int(size), point.tag)
            payloads[(int(size), point.tag)] = _parse_selected_bundle(bundle_path, labels)
    return payloads


def _target_rows(target_root: Path, target_size: int, target_point: CouplingPoint, labels: list[str]) -> dict[str, dict[str, float]]:
    bundle_path = _bundle_path(target_root, int(target_size), target_point.tag)
    payload = _parse_selected_bundle(bundle_path, labels)
    return _ratio_observables(payload)


def _interpolation_label(neighbors: list[tuple[CouplingPoint, float]]) -> str:
    if len(neighbors) == 1:
        return "direct"
    if len(neighbors) == 2:
        return "edge_interp"
    return "cell_interp"


def _write_refined_notes(
    output_path: Path,
    *,
    output_root: Path,
    target_root: Path,
    target_size: int,
    target_point: CouplingPoint,
    balance_mode: str,
    base_values: list[float],
    refined_values: list[float],
    gradient_step: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("# Reweighting q_w-Normalized Refined Surface\n\n")
        handle.write(f"Main raw family root: {output_root}\n\n")
        handle.write(f"Target raw root: {target_root}\n\n")
        handle.write(
            f"Target holdout: L={int(target_size)}, (r1,r2)=({float(target_point.r1):.3f},{float(target_point.r2):.3f})\n\n"
        )
        handle.write(f"Base grid values: {', '.join(f'{float(value):.3f}' for value in base_values)}\n\n")
        handle.write(f"Refined grid values: {', '.join(f'{float(value):.3f}' for value in refined_values)}\n\n")
        handle.write(
            "Inserted midpoint nodes are evaluated by reweighting the enclosing base-grid corner donors to the target coupling, "
            "then blending those donor predictions with bilinear weights inside each base-grid cell. Original 5x5 nodes keep their direct scored values.\n\n"
        )
        handle.write(
            "Gradient fields are rendered separately from a direct local reweight stencil at each original 5x5 node: "
            f"for each arrow, the score is re-evaluated at r +/- delta e_i using that node's own donor payload with delta={float(gradient_step):.6f}. "
            "The gradient figure therefore does not finite-difference the 2x-refined neighbor-interpolated surface.\n\n"
        )
        handle.write(
            "Observable set: midpoint and quarter ratios normalized by q_w = (L/4, L/4): mid_v/q_w, mid_u/q_w, mid_w/q_w, q_v/q_w, q_u/q_w.\n\n"
        )
        handle.write(f"Balance mode: {balance_mode}\n")


def _resolve_gradient_step(base_values: list[float], requested_step: float | None) -> float:
    if requested_step is not None:
        step = float(requested_step)
    elif len(base_values) >= 2:
        spacing = np.diff(np.asarray(base_values, dtype=float))
        step = 0.25 * float(np.min(np.abs(spacing)))
    else:
        step = 0.01
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError(f"gradient step must be positive, got {step}")
    return step


def _score_row_from_single_donor(
    donor_point: CouplingPoint,
    eval_point: CouplingPoint,
    *,
    fit_sizes: list[int],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
    payloads: dict[tuple[int, str], dict[str, Any]],
    direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]],
    reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]],
) -> dict[str, Any]:
    family_rows: dict[int, dict[str, dict[str, float]]] = {}
    same_point = math.isclose(float(donor_point.r1), float(eval_point.r1), abs_tol=1.0e-12) and math.isclose(
        float(donor_point.r2), float(eval_point.r2), abs_tol=1.0e-12
    )
    for size in fit_sizes:
        direct_key = (int(size), donor_point.tag)
        if same_point:
            if direct_key not in direct_family_cache:
                direct_family_cache[direct_key] = _ratio_observables(payloads[direct_key])
            family_rows[int(size)] = direct_family_cache[direct_key]
            continue

        cache_key = (int(size), donor_point.tag, eval_point.tag)
        if cache_key not in reweighted_family_cache:
            reweighted_family_cache[cache_key] = _reweighted_ratio_observables(payloads[direct_key], eval_point)
        family_rows[int(size)] = reweighted_family_cache[cache_key]
    return _score_row(
        r1=float(eval_point.r1),
        r2=float(eval_point.r2),
        point_tag=eval_point.tag,
        family_rows=family_rows,
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=float(target_x),
    )


def _finite_difference(center_value: float, minus_value: float, plus_value: float, step: float) -> float:
    if math.isfinite(minus_value) and math.isfinite(plus_value):
        return (float(plus_value) - float(minus_value)) / (2.0 * float(step))
    if math.isfinite(center_value) and math.isfinite(plus_value):
        return (float(plus_value) - float(center_value)) / float(step)
    if math.isfinite(center_value) and math.isfinite(minus_value):
        return (float(center_value) - float(minus_value)) / float(step)
    return float("nan")


def _normalized_vectors(dx: np.ndarray, dy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    magnitude = np.hypot(dx, dy)
    out_x = np.zeros_like(dx)
    out_y = np.zeros_like(dy)
    mask = magnitude > 0.0
    out_x[mask] = dx[mask] / magnitude[mask]
    out_y[mask] = dy[mask] / magnitude[mask]
    return out_x, out_y


def _plot_gradient_fields(
    rows: list[dict[str, Any]],
    *,
    balance_mode: str,
    base_values: list[float],
    fit_sizes: list[int],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
    payloads: dict[tuple[int, str], dict[str, Any]],
    direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]],
    reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]],
    target_point: CouplingPoint,
    gradient_step: float,
    gradient_dir: Path,
    base_name: str,
    title_prefix: str,
) -> dict[str, dict[str, float]]:
    metric_keys = _aggregate_plot_keys(balance_mode)
    metric_titles = _aggregate_plot_titles(balance_mode)
    base_rows = [
        row
        for row in rows
        if _is_base_point(CouplingPoint(float(row["r1"]), float(row["r2"])), base_values)
    ]
    grid_spacing = float(np.min(np.abs(np.diff(np.asarray(base_values, dtype=float))))) if len(base_values) >= 2 else 0.1
    arrow_length = 0.30 * grid_spacing
    base_x, base_y = np.meshgrid(np.asarray(base_values, dtype=float), np.asarray(base_values, dtype=float), indexing="ij")
    local_score_cache: dict[tuple[str, str], dict[str, Any]] = {}

    def _local_score(donor_point: CouplingPoint, eval_point: CouplingPoint) -> dict[str, Any]:
        cache_key = (donor_point.tag, eval_point.tag)
        if cache_key not in local_score_cache:
            local_score_cache[cache_key] = _score_row_from_single_donor(
                donor_point,
                eval_point,
                fit_sizes=fit_sizes,
                target_rows=target_rows,
                target_x=float(target_x),
                payloads=payloads,
                direct_family_cache=direct_family_cache,
                reweighted_family_cache=reweighted_family_cache,
            )
        return local_score_cache[cache_key]

    fig, axes = plt.subplots(1, len(metric_keys), figsize=(5.6 * len(metric_keys), 5.6), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    fig.suptitle(
        f"{title_prefix} base-grid score gradients from direct local reweight\n"
        f"star = panel minimum, x = target (1,1), arrow color = |grad score|, local step delta = {float(gradient_step):.4f}",
        fontsize=13,
    )

    summary: dict[str, dict[str, float]] = {}
    tsv_headers = ["r1", "r2"]
    tsv_rows: list[dict[str, float]] = []
    quivers: list[Any] = []

    for axis, key, title in zip(axes_flat, metric_keys, metric_titles):
        coarse_r1, coarse_r2, coarse_grid = _build_grid(base_rows, key)
        grad_r1 = np.full_like(coarse_grid, np.nan, dtype=float)
        grad_r2 = np.full_like(coarse_grid, np.nan, dtype=float)
        for i, r1_value in enumerate(coarse_r1):
            for j, r2_value in enumerate(coarse_r2):
                donor_point = CouplingPoint(float(r1_value), float(r2_value))
                center_value = float(coarse_grid[i, j])
                r1_minus = float(
                    _local_score(donor_point, CouplingPoint(float(r1_value) - float(gradient_step), float(r2_value)))[key]
                )
                r1_plus = float(
                    _local_score(donor_point, CouplingPoint(float(r1_value) + float(gradient_step), float(r2_value)))[key]
                )
                r2_minus = float(
                    _local_score(donor_point, CouplingPoint(float(r1_value), float(r2_value) - float(gradient_step)))[key]
                )
                r2_plus = float(
                    _local_score(donor_point, CouplingPoint(float(r1_value), float(r2_value) + float(gradient_step)))[key]
                )
                grad_r1[i, j] = _finite_difference(center_value, r1_minus, r1_plus, float(gradient_step))
                grad_r2[i, j] = _finite_difference(center_value, r2_minus, r2_plus, float(gradient_step))
        grad_mag = np.hypot(grad_r1, grad_r2)
        u, v = _normalized_vectors(grad_r1, grad_r2)
        display_u = arrow_length * u
        display_v = arrow_length * v

        x0, x1 = _grid_extent(coarse_r1)
        y0, y1 = _grid_extent(coarse_r2)
        axis.set_facecolor("white")
        axis.grid(True, color="#d1d5db", linewidth=0.8, alpha=0.8)
        axis.axhline(float(target_point.r2), color="#e5e7eb", linewidth=1.0, zorder=0)
        axis.axvline(float(target_point.r1), color="#e5e7eb", linewidth=1.0, zorder=0)
        axis.quiver(
            base_x,
            base_y,
            display_u.T,
            display_v.T,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            color="white",
            pivot="mid",
            width=0.015,
            headwidth=4.8,
            headlength=6.0,
            headaxislength=5.3,
            zorder=4,
        )
        quiver = axis.quiver(
            base_x,
            base_y,
            display_u.T,
            display_v.T,
            grad_mag.T,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            cmap="magma",
            pivot="mid",
            width=0.009,
            headwidth=4.4,
            headlength=5.6,
            headaxislength=5.0,
            zorder=5,
        )
        quivers.append(quiver)
        axis.scatter(base_x.ravel(), base_y.ravel(), s=18, c="white", edgecolors="#111827", linewidths=0.5, zorder=6)
        axis.scatter([float(target_point.r1)], [float(target_point.r2)], marker="x", s=90, color="#dc2626", linewidths=1.8, zorder=7)
        best_index = _best_index(coarse_grid, use_abs=False)
        best_text = "best: n/a"
        if best_index is not None:
            best_r1 = float(coarse_r1[best_index[0]])
            best_r2 = float(coarse_r2[best_index[1]])
            best_value = float(coarse_grid[best_index])
            axis.scatter([best_r1], [best_r2], marker="*", s=160, color="white", edgecolors="black", linewidths=0.8, zorder=7)
            best_text = f"best = ({best_r1:.2f}, {best_r2:.2f})\nscore = {best_value:.4g}"

        axis.set_title(title)
        _style_heatmap_axis(
            axis=axis,
            r1_values=coarse_r1,
            r2_values=coarse_r2,
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            show_xlabel=True,
            show_ylabel=(axis is axes_flat[0]),
        )
        axis.set_aspect("equal")
        axis.text(
            0.03,
            0.97,
            f"{best_text}\nmean |grad| = {float(np.nanmean(grad_mag)):.4g}\nmax |grad| = {float(np.nanmax(grad_mag)):.4g}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.86, "edgecolor": "#9ca3af"},
        )

        summary[key] = {
            "max_grad_norm": float(np.nanmax(grad_mag)),
            "mean_grad_norm": float(np.nanmean(grad_mag)),
        }

        tsv_headers.extend([f"{key}", f"{key}_dC_dr1", f"{key}_dC_dr2", f"{key}_grad_norm"])
        if not tsv_rows:
            for i, r1_value in enumerate(base_values):
                for j, r2_value in enumerate(base_values):
                    tsv_rows.append({"r1": float(r1_value), "r2": float(r2_value), "_i": float(i), "_j": float(j)})
        row_lookup = {(int(row["_i"]), int(row["_j"])): row for row in tsv_rows}
        for i, r1_value in enumerate(base_values):
            for j, r2_value in enumerate(base_values):
                row = row_lookup[(i, j)]
                row[key] = float(coarse_grid[i, j])
                row[f"{key}_dC_dr1"] = float(grad_r1[i, j])
                row[f"{key}_dC_dr2"] = float(grad_r2[i, j])
                row[f"{key}_grad_norm"] = float(grad_mag[i, j])

    if quivers:
        fig.colorbar(quivers[-1], ax=axes_flat, fraction=0.032, pad=0.02, label="|local score gradient|")

    figure_path = gradient_dir / f"{base_name}_aggregate_gradient_fields.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    tsv_path = gradient_dir / f"{base_name}_aggregate_gradient_fields.tsv"
    unique_headers = ["r1", "r2"] + [header for header in tsv_headers if header not in {"r1", "r2"}]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(unique_headers)
        for row in tsv_rows:
            writer.writerow([f"{float(row[key]):.10f}" for key in unique_headers])

    summary["figure_path"] = {"value": str(figure_path)}
    summary["tsv_path"] = {"value": str(tsv_path)}
    return summary


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    target_root = Path(args.target_root).resolve()
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else output_root / "heatmaps_target_L64_qw_norm_refined2x_neighbor_interp"
    )
    gradient_dir = (
        Path(args.gradient_dir).resolve()
        if args.gradient_dir
        else output_root / "gradient_fields_target_L64_qw_norm_local_reweight"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    gradient_dir.mkdir(parents=True, exist_ok=True)

    summary = _load_summary(output_root)
    fit_sizes = [int(value) for value in summary["fit_sizes"]]
    base_values = _base_grid_values(summary)
    refined_values = _refined_values(base_values, refine_factor=2)
    gradient_step = _resolve_gradient_step(base_values, args.gradient_step)
    base_points = [CouplingPoint(float(r1), float(r2)) for r1 in base_values for r2 in base_values]
    refined_points = [CouplingPoint(float(r1), float(r2)) for r1 in refined_values for r2 in refined_values]
    target_point = CouplingPoint(float(args.target_r1), float(args.target_r2))

    labels = [label for label, _ in _selected_specs_for_size(fit_sizes[0])]
    payloads = _load_payloads(output_root, fit_sizes, base_points, labels)
    target_rows = _target_rows(target_root, int(args.target_size), target_point, labels)
    target_x = 1.0 / float(args.target_size)

    direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]] = {}
    reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]] = {}

    rows: list[dict[str, Any]] = []
    for point in refined_points:
        neighbor_weights = _neighbor_weights(point, base_values)
        family_rows: dict[int, dict[str, dict[str, float]]] = {}
        for size in fit_sizes:
            if _is_base_point(point, base_values):
                direct_key = (int(size), point.tag)
                if direct_key not in direct_family_cache:
                    direct_family_cache[direct_key] = _ratio_observables(payloads[direct_key])
                family_rows[int(size)] = direct_family_cache[direct_key]
                continue

            weighted_rows: list[tuple[float, dict[str, dict[str, float]]]] = []
            for donor, weight in neighbor_weights:
                cache_key = (int(size), donor.tag, point.tag)
                if cache_key not in reweighted_family_cache:
                    reweighted_family_cache[cache_key] = _reweighted_ratio_observables(payloads[(int(size), donor.tag)], point)
                weighted_rows.append((float(weight), reweighted_family_cache[cache_key]))
            family_rows[int(size)] = _blend_stats(weighted_rows)

        row = _score_row(
            r1=float(point.r1),
            r2=float(point.r2),
            point_tag=point.tag,
            family_rows=family_rows,
            fit_sizes=fit_sizes,
            target_rows=target_rows,
            target_x=target_x,
        )
        row["source_mode"] = _interpolation_label(neighbor_weights)
        row["neighbor_tags"] = [donor.tag for donor, _ in neighbor_weights]
        row["neighbor_weights"] = [float(weight) for _, weight in neighbor_weights]
        rows.append(row)

    chi2_scales = _add_balanced_scores(rows, str(args.balance_mode))
    base_name = f"target_L{int(args.target_size)}_qw_norm_refined2x_neighbor_interp"

    _write_table(output_dir / f"{base_name}_landscape.tsv", rows)
    _write_summary(output_dir / f"{base_name}_best_points.tsv", rows, str(args.balance_mode))
    _write_refined_notes(
        output_dir / f"{base_name}_README.md",
        output_root=output_root,
        target_root=target_root,
        target_size=int(args.target_size),
        target_point=target_point,
        balance_mode=str(args.balance_mode),
        base_values=base_values,
        refined_values=refined_values,
        gradient_step=float(gradient_step),
    )
    _plot_signed_heatmaps(
        rows=rows,
        output_path=output_dir / f"{base_name}_signed_heatmaps.png",
        title=f"{args.title_prefix} signed heatmaps (2x refined neighbor reweight)",
        target_point=(float(args.target_r1), float(args.target_r2)),
        balance_mode=str(args.balance_mode),
    )
    _plot_positive_heatmaps(
        rows=rows,
        output_path=output_dir / f"{base_name}_score_heatmaps.png",
        title=f"{args.title_prefix} score heatmaps (2x refined neighbor reweight)",
        target_point=(float(args.target_r1), float(args.target_r2)),
        balance_mode=str(args.balance_mode),
    )

    gradient_base_name = f"target_L{int(args.target_size)}_qw_norm_local_reweight"
    gradient_summary = _plot_gradient_fields(
        rows,
        balance_mode=str(args.balance_mode),
        base_values=base_values,
        fit_sizes=fit_sizes,
        target_rows=target_rows,
        target_x=target_x,
        payloads=payloads,
        direct_family_cache=direct_family_cache,
        reweighted_family_cache=reweighted_family_cache,
        target_point=target_point,
        gradient_step=float(gradient_step),
        gradient_dir=gradient_dir,
        base_name=gradient_base_name,
        title_prefix=str(args.title_prefix),
    )

    best_key = _aggregate_plot_keys(str(args.balance_mode))[-1]
    best_row = min(rows, key=lambda row: float(row[best_key]))
    manifest = {
        "output_root": str(output_root),
        "target_root": str(target_root),
        "fit_sizes": fit_sizes,
        "target_size": int(args.target_size),
        "target_point": {"r1": float(target_point.r1), "r2": float(target_point.r2)},
        "base_values": [float(value) for value in base_values],
        "refined_values": [float(value) for value in refined_values],
        "balance_mode": str(args.balance_mode),
        "gradient_method": "local_reweight_central_difference",
        "gradient_step": float(gradient_step),
        "chi2_scales": {key: float(value) for key, value in chi2_scales.items()},
        "best_total_like_row": {
            "point_tag": str(best_row["point_tag"]),
            "r1": float(best_row["r1"]),
            "r2": float(best_row["r2"]),
            "metric_key": str(best_key),
            "metric_value": float(best_row[best_key]),
            "source_mode": str(best_row["source_mode"]),
            "neighbor_tags": list(best_row["neighbor_tags"]),
            "neighbor_weights": [float(value) for value in best_row["neighbor_weights"]],
        },
        "gradient_summary": gradient_summary,
    }
    (output_dir / f"{base_name}_summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    print(f"wrote {output_dir / f'{base_name}_landscape.tsv'}")
    print(f"wrote {output_dir / f'{base_name}_best_points.tsv'}")
    print(f"wrote {output_dir / f'{base_name}_README.md'}")
    print(f"wrote {output_dir / f'{base_name}_signed_heatmaps.png'}")
    print(f"wrote {output_dir / f'{base_name}_score_heatmaps.png'}")
    print(f"wrote {gradient_dir / f'{gradient_base_name}_aggregate_gradient_fields.png'}")
    print(f"wrote {gradient_dir / f'{gradient_base_name}_aggregate_gradient_fields.tsv'}")


if __name__ == "__main__":
    main()