#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from test_geometry_match_grid_interpolation import (  # noqa: E402
    _block_slices,
    _connected_from_weights,
    _evaluate_power_model_on_x,
    _fit_blind_power_model,
    _panel_specs,
    _parse_selected_bundle,
    _predict_at_target_x,
    _ratio_stat,
    _selected_specs_for_size,
    _token,
)


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render raw direct-only FSS panels from cached selected-observable bundles "
            "with only the finite-size points and the target marker."
        )
    )
    parser.add_argument(
        "--output-root",
        default=str(
            RESPONSIBLE_ROOT
            / "reweighting"
            / "iso111_grid5x5_sizes4_8_12_16_20_24_28_32_hi10k_20260606"
        ),
    )
    parser.add_argument("--target-size", type=int, default=32)
    parser.add_argument("--target-r1", type=float, default=1.0)
    parser.add_argument("--target-r2", type=float, default=1.0)
    parser.add_argument(
        "--point-tags",
        nargs="+",
        default=[f"r1_{_token(1.0)}__r2_{_token(1.0)}"],
        help="Point tags to plot. Defaults to the isotropic target point.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for the raw FSS plots. Defaults to <output-root>/raw_fss_panels.",
    )
    parser.add_argument(
        "--target-root",
        default=None,
        help="Optional separate output root to source the target bundle from.",
    )
    parser.add_argument(
        "--fit-ratios",
        action="store_true",
        help="Overlay A + B x^omega fits on ratio panels and annotate the continuum constant.",
    )
    return parser.parse_args()


def _infer_sizes(output_root: Path) -> list[int]:
    raw_root = output_root / "raw"
    sizes: list[int] = []
    for child in raw_root.iterdir():
        if not child.is_dir() or not child.name.startswith("L"):
            continue
        sizes.append(int(child.name[1:]))
    resolved = sorted(set(sizes))
    if not resolved:
        raise FileNotFoundError(f"no size directories found under {raw_root}")
    return resolved


def _bundle_path(output_root: Path, size: int, point_tag: str) -> Path:
    point_root = output_root / "raw" / f"L{int(size)}" / point_tag
    matches = sorted(point_root.glob("*/selected_observables_bundle.dat"))
    if not matches:
        raise FileNotFoundError(f"missing selected_observables_bundle.dat under {point_root}")
    return matches[0]


def _connected_stat(corr: np.ndarray, mag: np.ndarray, *, n_blocks: int = 16) -> dict[str, float]:
    n_samples = int(mag.size)
    if n_samples < 16:
        raise ValueError("need at least 16 measurements for connected-correlator jackknife")
    weights = np.full(n_samples, 1.0 / float(n_samples), dtype=float)
    value = _connected_from_weights(weights, corr, mag)

    leave_one_out: list[float] = []
    for lo, hi in _block_slices(n_samples, n_blocks):
        mask = np.ones(n_samples, dtype=bool)
        mask[lo:hi] = False
        if np.count_nonzero(mask) < 8:
            continue
        sub_weights = np.full(np.count_nonzero(mask), 1.0 / float(np.count_nonzero(mask)), dtype=float)
        leave_one_out.append(
            _connected_from_weights(sub_weights, np.asarray(corr[mask], dtype=float), np.asarray(mag[mask], dtype=float))
        )
    if len(leave_one_out) < 2:
        sigma = float("nan")
    else:
        jk = np.asarray(leave_one_out, dtype=float)
        sigma = float(np.sqrt((jk.size - 1.0) / jk.size * np.sum((jk - np.mean(jk)) ** 2)))
    return {"value": float(value), "sigma": float(sigma)}


def _connected_observables(payload: dict[str, object]) -> dict[str, dict[str, float]]:
    mag = np.asarray(payload["mag"], dtype=float)
    corr = dict(payload["corr"])
    rows: dict[str, dict[str, float]] = {}
    for label, samples in corr.items():
        rows[str(label)] = _connected_stat(np.asarray(samples, dtype=float), mag)
    return rows


def _ratio_observables(payload: dict[str, object]) -> dict[str, dict[str, float]]:
    mag = np.asarray(payload["mag"], dtype=float)
    corr = dict(payload["corr"])
    rows: dict[str, dict[str, float]] = {}
    for panel in [*_panel_specs(), *_cross_scale_ratio_specs(), *_q_w_normalized_specs()]:
        rows[str(panel["label"])] = _ratio_stat(
            np.asarray(corr[str(panel["numerator"])], dtype=float),
            np.asarray(corr[str(panel["denominator"])], dtype=float),
            mag,
        )
    return rows


def _cross_scale_ratio_specs() -> list[dict[str, str]]:
    return [
        {"label": "mid_v_over_q_v", "numerator": "mid_v", "denominator": "q_v"},
        {"label": "mid_u_over_q_u", "numerator": "mid_u", "denominator": "q_u"},
        {"label": "mid_w_over_q_w", "numerator": "mid_w", "denominator": "q_w"},
    ]


def _q_w_normalized_specs() -> list[dict[str, str]]:
    return [
        {"label": "mid_v_over_q_w", "numerator": "mid_v", "denominator": "q_w"},
        {"label": "mid_u_over_q_w", "numerator": "mid_u", "denominator": "q_w"},
        {"label": "mid_w_over_q_w", "numerator": "mid_w", "denominator": "q_w"},
        {"label": "q_v_over_q_w", "numerator": "q_v", "denominator": "q_w"},
        {"label": "q_u_over_q_w", "numerator": "q_u", "denominator": "q_w"},
    ]


def _correlator_group() -> dict[str, object]:
    return {
        "slug": "correlators_raw",
        "title": "raw connected correlator FSS",
        "y_axis_label": "Connected correlator",
        "panels": [
            {"label": "mid_v", "title": "mid_v", "color": "#1d4ed8"},
            {"label": "mid_u", "title": "mid_u", "color": "#047857"},
            {"label": "mid_w", "title": "mid_w", "color": "#b45309"},
            {"label": "q_v", "title": "q_v", "color": "#7c2d12"},
            {"label": "q_u", "title": "q_u", "color": "#7c3aed"},
            {"label": "q_w", "title": "q_w", "color": "#0891b2"},
        ],
    }


def _midpoint_correlator_group() -> dict[str, object]:
    return {
        "slug": "midpoint_correlators_raw",
        "title": "raw midpoint correlators",
        "y_axis_label": "Connected correlator",
        "panels": [
            {"label": "mid_v", "title": "mid_v", "color": "#1d4ed8"},
            {"label": "mid_u", "title": "mid_u", "color": "#047857"},
            {"label": "mid_w", "title": "mid_w", "color": "#b45309"},
        ],
    }


def _quarter_correlator_group() -> dict[str, object]:
    return {
        "slug": "quarter_correlators_raw",
        "title": "raw quarter-point correlators",
        "y_axis_label": "Connected correlator",
        "panels": [
            {"label": "q_v", "title": "q_v", "color": "#7c2d12"},
            {"label": "q_u", "title": "q_u", "color": "#7c3aed"},
            {"label": "q_w", "title": "q_w", "color": "#0891b2"},
        ],
    }


def _ratio_groups() -> list[dict[str, object]]:
    return [
        {
            "slug": "midpoint_ratios_raw",
            "title": "raw midpoint anchored ratios",
            "y_axis_label": "Correlator ratio",
            "panels": [
                {"label": "mid_v_anchor", "title": "mid_v_anchor", "color": "#1d4ed8"},
                {"label": "mid_u_anchor", "title": "mid_u_anchor", "color": "#047857"},
                {"label": "mid_w_anchor", "title": "mid_w_anchor", "color": "#b45309"},
                {"label": "mid_v_over_u", "title": "mid_v_over_u", "color": "#7c2d12"},
                {"label": "mid_w_over_u", "title": "mid_w_over_u", "color": "#7c3aed"},
            ],
        },
        {
            "slug": "quarter_ratios_raw",
            "title": "raw quarter-point anchored ratios",
            "y_axis_label": "Correlator ratio",
            "panels": [
                {"label": "q_v_anchor", "title": "q_v_anchor", "color": "#1d4ed8"},
                {"label": "q_u_anchor", "title": "q_u_anchor", "color": "#047857"},
                {"label": "q_w_anchor", "title": "q_w_anchor", "color": "#b45309"},
                {"label": "q_v_over_u", "title": "q_v_over_u", "color": "#7c2d12"},
                {"label": "q_w_over_u", "title": "q_w_over_u", "color": "#7c3aed"},
            ],
        },
        {
            "slug": "half_over_quarter_ratios_raw",
            "title": "raw half-point / quarter-point ratios",
            "y_axis_label": "Correlator ratio",
            "panels": [
                {"label": "mid_v_over_q_v", "title": "mid_v_over_q_v", "color": "#1d4ed8"},
                {"label": "mid_u_over_q_u", "title": "mid_u_over_q_u", "color": "#047857"},
                {"label": "mid_w_over_q_w", "title": "mid_w_over_q_w", "color": "#b45309"},
            ],
        },
        {
            "slug": "midpoint_qw_normalized_raw",
            "title": "raw midpoint ratios normalized by q_w",
            "y_axis_label": "Correlator ratio",
            "panels": [
                {"label": "mid_v_over_q_w", "title": "mid_v_over_q_w", "color": "#1d4ed8"},
                {"label": "mid_u_over_q_w", "title": "mid_u_over_q_w", "color": "#047857"},
                {"label": "mid_w_over_q_w", "title": "mid_w_over_q_w", "color": "#b45309"},
            ],
        },
        {
            "slug": "quarter_qw_normalized_raw",
            "title": "raw quarter-point ratios normalized by q_w",
            "y_axis_label": "Correlator ratio",
            "panels": [
                {"label": "q_v_over_q_w", "title": "q_v_over_q_w", "color": "#7c2d12"},
                {"label": "q_u_over_q_w", "title": "q_u_over_q_w", "color": "#7c3aed"},
            ],
        },
    ]


def _layout_for_panel_count(panel_count: int) -> tuple[int, int, tuple[float, float], list[float], float, float, float]:
    if panel_count <= 3:
        cols = max(panel_count, 1)
        return 1, cols, (5.8 * cols, 5.3), [0.0, 0.0, 1.0, 0.82], 0.985, 0.92, 0.885
    return 2, 3, (18.0, 9.8), [0.0, 0.0, 1.0, 0.90], 0.985, 0.945, 0.965


def _pretty_panel_title(raw_title: str) -> str:
    if "_over_" in raw_title:
        left, right = raw_title.split("_over_", 1)
        return f"{left} / {right}"
    if raw_title.endswith("_anchor"):
        return f"{raw_title[:-7]} / anchor"
    return raw_title.replace("_", " ")


def _build_series_details(
    family_rows: dict[int, dict[str, dict[str, float]]],
    *,
    sizes: list[int],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
) -> dict[str, dict[str, object]]:
    x_values = np.asarray([1.0 / float(size) for size in sizes], dtype=float)
    details: dict[str, dict[str, object]] = {}
    all_labels = set(target_rows.keys())
    for size in sizes:
        all_labels.update(family_rows[size].keys())
    for label in sorted(all_labels):
        y_values = np.asarray([float(family_rows[size][label]["value"]) for size in sizes], dtype=float)
        sigma_values = np.asarray([float(family_rows[size][label]["sigma"]) for size in sizes], dtype=float)
        details[label] = {
            "x_values": x_values,
            "y_values": y_values,
            "sigma_values": sigma_values,
            "target": {
                "x": float(target_x),
                "value": float(target_rows[label]["value"]),
                "sigma": float(target_rows[label]["sigma"]),
            },
        }
    return details


def _render_group_sheet(
    output_path: Path,
    *,
    point_tag: str,
    point_label: str,
    sizes: list[int],
    target_size: int,
    group: dict[str, object],
    details: dict[str, dict[str, object]],
    fit_ratios: bool = False,
) -> None:
    panels = list(group["panels"])
    nrows, ncols, figsize, tight_rect, title_y, subtitle_y, legend_y = _layout_for_panel_count(len(panels))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = list(axes.ravel())
    fig.suptitle(f"{point_label} {str(group['title'])}", fontsize=16, y=title_y)
    fit_note = "; fit: A + B x^omega" if fit_ratios else ""
    fig.text(
        0.5,
        subtitle_y,
        f"direct MC only; untwisted sizes={','.join(str(int(size)) for size in sizes)}; target=L={int(target_size)}; no extrapolation{fit_note}",
        ha="center",
        va="top",
        fontsize=10,
        color="#444444",
    )

    legend_handles = None
    legend_labels = None
    for axis, panel in zip(axes_flat, panels):
        label = str(panel["label"])
        panel_detail = details[label]
        x_values = np.asarray(panel_detail["x_values"], dtype=float)
        y_values = np.asarray(panel_detail["y_values"], dtype=float)
        sigma_values = np.asarray(panel_detail["sigma_values"], dtype=float)
        target = dict(panel_detail["target"])
        color = str(panel["color"])
        x_min = min(float(np.min(x_values)), float(target["x"]))
        x_max = max(float(np.max(x_values)), float(target["x"]))

        axis.errorbar(
            x_values,
            y_values,
            yerr=sigma_values,
            fmt="o",
            color=color,
            ecolor=color,
            capsize=3,
            markersize=5,
            markeredgecolor="white",
            markeredgewidth=0.8,
            linewidth=1.0,
            label="direct",
        )
        axis.errorbar(
            [float(target["x"])],
            [float(target["value"])],
            yerr=[float(target["sigma"])],
            fmt="none",
            color="#dc2626",
            ecolor="#dc2626",
            capsize=3,
            linewidth=1.0,
        )
        axis.plot(
            [float(target["x"])],
            [float(target["value"])],
            linestyle="none",
            marker="*",
            color="#dc2626",
            markersize=13,
            markeredgecolor="white",
            markeredgewidth=0.9,
            label=f"target L={int(target_size)}",
            zorder=5,
        )
        if fit_ratios:
            fit_payload = _fit_blind_power_model(x_values, y_values, sigma_values)
            fit_x = np.geomspace(max(x_min * 0.98, 1.0e-6), x_max * 1.02, 300)
            fit_y = _evaluate_power_model_on_x(
                fit_x,
                float(fit_payload["A"]),
                float(fit_payload["B"]),
                float(fit_payload["omega"]),
            )
            pred_target, pred_sigma = _predict_at_target_x(fit_payload, float(target["x"]))
            axis.plot(
                fit_x,
                fit_y,
                color=color,
                linewidth=1.6,
                linestyle="--",
                alpha=0.95,
                label="fit",
                zorder=2,
            )
            axis.text(
                0.03,
                0.97,
                (
                    f"A_inf={float(fit_payload['A']):.6f} +/- {float(fit_payload.get('sigma_A', float('nan'))):.6f}\n"
                    f"omega={float(fit_payload['omega']):.2f}\n"
                    f"pred@L{int(target_size)}={float(pred_target):.6f} +/- {float(pred_sigma):.6f}"
                ),
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8.5,
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
            )
        axis.set_title(_pretty_panel_title(str(panel["title"])), fontsize=14)
        axis.set_xlabel("1 / sqrt(lattice volume)", fontsize=11)
        axis.set_ylabel(str(group["y_axis_label"]), fontsize=11)
        axis.set_xscale("log")
        axis.set_xlim(max(x_min * 0.92, 1.0e-6), x_max * 1.08)
        axis.grid(True, which="both", alpha=0.35)
        axis.tick_params(axis="both", labelsize=10.5)
        if legend_handles is None:
            legend_handles, legend_labels = axis.get_legend_handles_labels()

    for axis in axes_flat[len(panels):]:
        axis.axis("off")

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, legend_y),
            ncol=2,
            frameon=False,
            fontsize=11,
        )
    fig.tight_layout(rect=tight_rect)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _untoken(token: str) -> float:
    if token.startswith("m"):
        return -float(token[1:].replace("p", "."))
    return float(token.replace("p", "."))


def _point_label_from_tag(point_tag: str) -> str:
    r1_part, r2_part = point_tag.split("__", 1)
    r1 = _untoken(r1_part.removeprefix("r1_"))
    r2 = _untoken(r2_part.removeprefix("r2_"))
    return f"Point ({r1:.6f}, {r2:.6f})"


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    sizes = _infer_sizes(output_root)
    target_size = int(args.target_size)
    target_root = Path(args.target_root).resolve() if args.target_root else output_root
    target_sizes = _infer_sizes(target_root)
    if target_size not in target_sizes:
        raise ValueError(f"target size L={target_size} is not present under {target_root / 'raw'}")

    output_dir = Path(args.output_dir).resolve() if args.output_dir else output_root / "raw_fss_panels"
    output_dir.mkdir(parents=True, exist_ok=True)

    labels = [label for label, _ in _selected_specs_for_size(sizes[0])]
    target_tag = f"r1_{_token(float(args.target_r1))}__r2_{_token(float(args.target_r2))}"
    target_payload = _parse_selected_bundle(_bundle_path(target_root, target_size, target_tag), labels)
    target_corr_rows = _connected_observables(target_payload)
    target_ratio_rows = _ratio_observables(target_payload)
    target_x = 1.0 / float(target_size)

    for point_tag in list(args.point_tags):
        direct_payloads = {
            int(size): _parse_selected_bundle(_bundle_path(output_root, int(size), point_tag), labels)
            for size in sizes
        }
        corr_family = {int(size): _connected_observables(payload) for size, payload in direct_payloads.items()}
        ratio_family = {int(size): _ratio_observables(payload) for size, payload in direct_payloads.items()}

        point_label = _point_label_from_tag(point_tag)
        corr_details = _build_series_details(corr_family, sizes=sizes, target_rows=target_corr_rows, target_x=target_x)
        ratio_details = _build_series_details(ratio_family, sizes=sizes, target_rows=target_ratio_rows, target_x=target_x)

        _render_group_sheet(
            output_dir / f"{point_tag}_{str(_correlator_group()['slug'])}.png",
            point_tag=point_tag,
            point_label=point_label,
            sizes=sizes,
            target_size=target_size,
            group=_correlator_group(),
            details=corr_details,
            fit_ratios=False,
        )
        _render_group_sheet(
            output_dir / f"{point_tag}_{str(_midpoint_correlator_group()['slug'])}.png",
            point_tag=point_tag,
            point_label=point_label,
            sizes=sizes,
            target_size=target_size,
            group=_midpoint_correlator_group(),
            details=corr_details,
            fit_ratios=False,
        )
        _render_group_sheet(
            output_dir / f"{point_tag}_{str(_quarter_correlator_group()['slug'])}.png",
            point_tag=point_tag,
            point_label=point_label,
            sizes=sizes,
            target_size=target_size,
            group=_quarter_correlator_group(),
            details=corr_details,
            fit_ratios=False,
        )
        for group in _ratio_groups():
            _render_group_sheet(
                output_dir / f"{point_tag}_{str(group['slug'])}.png",
                point_tag=point_tag,
                point_label=point_label,
                sizes=sizes,
                target_size=target_size,
                group=group,
                details=ratio_details,
                fit_ratios=bool(args.fit_ratios),
            )


if __name__ == "__main__":
    main()