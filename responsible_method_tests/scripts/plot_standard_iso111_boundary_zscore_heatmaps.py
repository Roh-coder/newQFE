#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm


_PANEL_BASE_FIELDS: list[dict[str, Any]] = [
    {"stem": "v4", "z_aliases": ["v4_z", "v/4_z"], "chi2_aliases": ["v4_chi2", "v/4_chi2"], "label": "v/4"},
    {"stem": "u4", "z_aliases": ["u4_z", "u/4_z"], "chi2_aliases": ["u4_chi2", "u/4_chi2"], "label": "u/4"},
    {"stem": "w4", "z_aliases": ["w4_z", "w/4_z"], "chi2_aliases": ["w4_chi2", "w/4_chi2"], "label": "w/4"},
    {
        "stem": "v4_over_u4",
        "z_aliases": ["v4_over_u4_z", "(v/4)/(u/4)_z"],
        "chi2_aliases": ["v4_over_u4_chi2", "(v/4)/(u/4)_chi2"],
        "label": "(v/4)/(u/4)",
    },
    {
        "stem": "w4_over_u4",
        "z_aliases": ["w4_over_u4_z", "(w/4)/(u/4)_z"],
        "chi2_aliases": ["w4_over_u4_chi2", "(w/4)/(u/4)_chi2"],
        "label": "(w/4)/(u/4)",
    },
]

_AGGREGATE_BASE_FIELDS: list[dict[str, Any]] = [
    {
        "stem": "corr",
        "chi2_aliases": ["corr_chi2", "corr_z2"],
        "title": "corr chi^2",
        "label": "corr",
    },
    {
        "stem": "ratio",
        "chi2_aliases": ["ratio_chi2", "ratio_z2"],
        "title": "ratio chi^2",
        "label": "ratio",
    },
    {
        "stem": "total",
        "chi2_aliases": ["chi2_sum", "z2_sum"],
        "title": "total chi^2",
        "label": "total",
    },
]

DEFAULT_INPUT = os.path.join(
    "/workspaces/newQFE",
    "responsible_method_tests",
    "standard",
    "iso111_boundary_quarter_fss_sweep_jackknife",
    "iso111_boundary_quarter_fss_zscores_partial.tsv",
)


def _field_specs(score_mode: str) -> list[dict[str, Any]]:
    if score_mode == "z":
        fields: list[dict[str, Any]] = [
            {
                "key": f"{base['stem']}_z",
                "aliases": list(base["z_aliases"]),
                "title": f"z[{base['label']}]",
                "group": "signed",
            }
            for base in _PANEL_BASE_FIELDS
        ]
        for base in _AGGREGATE_BASE_FIELDS:
            fields.append(
                {
                    "key": f"{base['stem']}_chi2",
                    "aliases": list(base["chi2_aliases"]),
                    "title": str(base["title"]),
                    "group": "aggregate",
                }
            )
        return fields

    if score_mode == "chi2":
        fields = [
            {
                "key": f"{base['stem']}_chi2",
                "aliases": list(base["chi2_aliases"]),
                "square_aliases": list(base["z_aliases"]),
                "title": f"chi^2[{base['label']}]",
                "group": "panel_positive",
            }
            for base in _PANEL_BASE_FIELDS
        ]
        for base in _AGGREGATE_BASE_FIELDS:
            fields.append(
                {
                    "key": f"{base['stem']}_chi2",
                    "aliases": list(base["chi2_aliases"]),
                    "title": str(base["title"]),
                    "group": "aggregate_positive",
                }
            )
        return fields

    raise KeyError(f"unknown score mode: {score_mode}")


def _default_output_path(input_path: str, score_mode: str) -> str:
    root, _ = os.path.splitext(input_path)
    if score_mode == "chi2":
        return f"{root}_chi2_heatmaps.png"
    return f"{root}_heatmaps.png"


def _load_rows(input_path: str) -> list[dict[str, Any]]:
    with open(input_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows


def _grid_extent(values: list[float]) -> tuple[float, float]:
    if len(values) == 1:
        return values[0] - 0.5, values[0] + 0.5
    steps = np.diff(np.asarray(values, dtype=float))
    step = float(np.median(steps))
    return float(values[0] - 0.5 * step), float(values[-1] + 0.5 * step)


def _extract_value(row: dict[str, Any], field: dict[str, Any]) -> float:
    for alias in field.get("aliases", []):
        if alias in row and row[alias] not in (None, ""):
            return float(row[alias])
    for alias in field.get("square_aliases", []):
        if alias in row and row[alias] not in (None, ""):
            value = float(row[alias])
            return float(value ** 2)
    return float("nan")


def _build_grids(rows: list[dict[str, Any]], field_specs: list[dict[str, Any]]) -> dict[str, Any]:
    r1_values = sorted({float(row["r1"]) for row in rows})
    r2_values = sorted({float(row["r2"]) for row in rows})
    lookup = {(float(row["r1"]), float(row["r2"])): row for row in rows}

    grids: dict[str, np.ndarray] = {}
    for field in field_specs:
        grid = np.full((len(r1_values), len(r2_values)), np.nan, dtype=float)
        for i, r1_value in enumerate(r1_values):
            for j, r2_value in enumerate(r2_values):
                row = lookup.get((r1_value, r2_value))
                if row is None:
                    continue
                value = _extract_value(row, field)
                if not math.isfinite(value):
                    continue
                grid[i, j] = value
        grids[str(field["key"])] = grid
    return {
        "r1_values": r1_values,
        "r2_values": r2_values,
        "grids": grids,
        "filled_cells": len(rows),
        "total_cells": len(r1_values) * len(r2_values),
    }


def _finite_values(grids: dict[str, np.ndarray], keys: list[str]) -> np.ndarray:
    collected: list[np.ndarray] = []
    for key in keys:
        values = grids[key]
        finite = values[np.isfinite(values)]
        if finite.size:
            collected.append(finite)
    if not collected:
        return np.asarray([], dtype=float)
    return np.concatenate(collected)


def _best_index(values: np.ndarray, *, use_abs: bool) -> tuple[int, int] | None:
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        return None
    score = np.abs(values) if use_abs else values
    masked_score = np.where(finite_mask, score, np.inf)
    flat_index = int(np.argmin(masked_score))
    return tuple(int(v) for v in np.unravel_index(flat_index, values.shape))


def _annotation_color(value: float, *, group: str, signed_limit: float, aggregate_limit: float) -> str:
    if group == "signed":
        if signed_limit <= 0.0:
            return "black"
        return "white" if abs(value) >= 0.55 * signed_limit else "black"
    if aggregate_limit <= 0.0:
        return "black"
    return "white" if value >= 0.55 * aggregate_limit else "black"


def _format_value(value: float, *, group: str) -> str:
    if group == "signed":
        return f"{value:+.2f}"
    return f"{value:.2f}"


def _plot_heatmaps(
    *,
    input_path: str,
    output_path: str,
    annotate: bool,
    figure_title: str | None,
    score_mode: str,
) -> str:
    rows = _load_rows(input_path)
    field_specs = _field_specs(score_mode)
    grid_payload = _build_grids(rows, field_specs)
    r1_values = list(grid_payload["r1_values"])
    r2_values = list(grid_payload["r2_values"])
    grids = dict(grid_payload["grids"])
    filled_cells = int(grid_payload["filled_cells"])
    total_cells = int(grid_payload["total_cells"])

    signed_keys = [str(field["key"]) for field in field_specs if field["group"] == "signed"]
    panel_positive_keys = [str(field["key"]) for field in field_specs if field["group"] == "panel_positive"]
    aggregate_keys = [
        str(field["key"])
        for field in field_specs
        if field["group"] in {"aggregate", "aggregate_positive"}
    ]
    signed_values = _finite_values(grids, signed_keys)
    panel_positive_values = _finite_values(grids, panel_positive_keys)
    aggregate_values = _finite_values(grids, aggregate_keys)
    signed_limit = float(np.max(np.abs(signed_values))) if signed_values.size else 1.0
    panel_positive_limit = float(np.max(panel_positive_values)) if panel_positive_values.size else 1.0
    aggregate_limit = float(np.max(aggregate_values)) if aggregate_values.size else 1.0
    if signed_limit <= 0.0:
        signed_limit = 1.0
    if panel_positive_limit <= 0.0:
        panel_positive_limit = 1.0
    if aggregate_limit <= 0.0:
        aggregate_limit = 1.0

    signed_norm = TwoSlopeNorm(vcenter=0.0, vmin=-signed_limit, vmax=signed_limit)
    panel_positive_norm = Normalize(vmin=0.0, vmax=panel_positive_limit)
    aggregate_norm = Normalize(vmin=0.0, vmax=aggregate_limit)
    signed_cmap = plt.get_cmap("RdBu_r").copy()
    positive_cmap = plt.get_cmap("viridis_r").copy()
    signed_cmap.set_bad(color="#e5e7eb")
    positive_cmap.set_bad(color="#e5e7eb")

    x0, x1 = _grid_extent(r1_values)
    y0, y1 = _grid_extent(r2_values)
    extent = [x0, x1, y0, y1]

    fig, axes = plt.subplots(2, 4, figsize=(18.5, 9.4), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    signed_axes = []
    panel_positive_axes = []
    aggregate_axes = []
    signed_image = None
    panel_positive_image = None
    aggregate_image = None

    for axis, field in zip(axes_flat, field_specs):
        key = str(field["key"])
        group = str(field["group"])
        values = np.asarray(grids[key], dtype=float)
        if group == "signed":
            cmap = signed_cmap
            norm = signed_norm
        elif group == "panel_positive":
            cmap = positive_cmap
            norm = panel_positive_norm
        else:
            cmap = positive_cmap
            norm = aggregate_norm
        image = axis.imshow(
            values.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap=cmap,
            norm=norm,
        )
        if group == "signed":
            signed_axes.append(axis)
            signed_image = image
        elif group == "panel_positive":
            panel_positive_axes.append(axis)
            panel_positive_image = image
        else:
            aggregate_axes.append(axis)
            aggregate_image = image

        axis.scatter([1.0], [1.0], marker="x", s=90, color="black", linewidths=1.8)
        best_index = _best_index(values, use_abs=(group == "signed"))
        if best_index is not None:
            best_r1 = r1_values[best_index[0]]
            best_r2 = r2_values[best_index[1]]
            axis.scatter(
                [best_r1],
                [best_r2],
                marker="*",
                s=180,
                color="white",
                edgecolors="black",
                linewidths=0.8,
                zorder=4,
            )

        if annotate:
            for i, r1_value in enumerate(r1_values):
                for j, r2_value in enumerate(r2_values):
                    value = values[i, j]
                    if not math.isfinite(float(value)):
                        axis.text(
                            float(r1_value),
                            float(r2_value),
                            "--",
                            ha="center",
                            va="center",
                            fontsize=8,
                            color="#6b7280",
                        )
                        continue
                    axis.text(
                        float(r1_value),
                        float(r2_value),
                        _format_value(float(value), group=group),
                        ha="center",
                        va="center",
                        fontsize=8.2,
                        color=_annotation_color(
                            float(value),
                            group=group,
                            signed_limit=signed_limit,
                            aggregate_limit=(panel_positive_limit if group == "panel_positive" else aggregate_limit),
                        ),
                    )

        axis.set_title(str(field["title"]), fontsize=11)
        axis.set_xlabel("r1 = k1 / k3")
        axis.set_ylabel("r2 = k2 / k3")
        axis.set_xticks(r1_values)
        axis.set_yticks(r2_values)
        axis.set_xlim(x0, x1)
        axis.set_ylim(y0, y1)
        axis.grid(False)

    if signed_image is not None and signed_axes:
        fig.colorbar(
            signed_image,
            ax=signed_axes,
            fraction=0.025,
            pad=0.025,
            label="signed z-score",
        )
    if panel_positive_image is not None and panel_positive_axes:
        fig.colorbar(
            panel_positive_image,
            ax=panel_positive_axes,
            fraction=0.025,
            pad=0.025,
            label="per-panel chi^2",
        )
    if aggregate_image is not None and aggregate_axes:
        fig.colorbar(
            aggregate_image,
            ax=aggregate_axes,
            fraction=0.025,
            pad=0.025,
            label="aggregate chi^2",
        )

    if figure_title:
        title = figure_title
    elif score_mode == "chi2":
        title = "Iso111 quarter-point boundary FSS chi^2 heatmaps"
    else:
        title = "Iso111 quarter-point boundary FSS z-score heatmaps"
    fig.suptitle(
        f"{title}\nfilled cells: {filled_cells}/{total_cells}; x = isotropic (1,1), star = per-panel best",
        fontsize=14,
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot heatmaps from an iso111 boundary z-score TSV.")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input z-score TSV.")
    parser.add_argument("--output", default=None, help="Output PNG path. Defaults to <input stem>_heatmaps.png.")
    parser.add_argument("--title", default=None, help="Optional figure title override.")
    parser.add_argument(
        "--score-mode",
        choices=("chi2", "z"),
        default="chi2",
        help="Plot chi^2 contributions or the original signed z-scores.",
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="Disable per-cell numeric annotations.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = os.path.abspath(str(args.input))
    output_path = os.path.abspath(str(args.output)) if args.output else _default_output_path(input_path, str(args.score_mode))
    written_path = _plot_heatmaps(
        input_path=input_path,
        output_path=output_path,
        annotate=not bool(args.no_annotate),
        figure_title=str(args.title) if args.title else None,
        score_mode=str(args.score_mode),
    )
    print(f"wrote {written_path}")


if __name__ == "__main__":
    main()