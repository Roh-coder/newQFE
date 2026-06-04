#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from typing import Callable

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


DEFAULT_INPUT = os.path.join(
    "/workspaces/newQFE",
    "responsible_method_tests",
    "standard",
    "acute456_boundary_quarter_fss_sweep_anchor_ratio",
    "acute456_boundary_quarter_fss_anchor_ratio_small_target_zscores.tsv",
)
DEFAULT_OUTPUT_DIR = os.path.join(
    "/workspaces/newQFE",
    "responsible_method_tests",
    "standard",
    "mixed_scores",
)


@dataclass(frozen=True)
class Candidate:
    key: str
    title: str
    formula: str
    description: str
    evaluator: Callable[[dict[str, float]], float]


def _row_value(row: dict[str, float], key: str) -> float:
    return float(row[key])


def _mean_sector_chi2(row: dict[str, float], key: str, count: int) -> float:
    return _row_value(row, key) / float(count)


def _corr_mean(row: dict[str, float]) -> float:
    return _mean_sector_chi2(row, "corr_chi2", 3)


def _ratio_mean(row: dict[str, float]) -> float:
    return _mean_sector_chi2(row, "ratio_chi2", 2)


def _safe_sqrt(value: float) -> float:
    return math.sqrt(max(float(value), 0.0))


CANDIDATES: tuple[Candidate, ...] = (
    Candidate(
        key="sector_balanced_rms",
        title="Balanced sectors",
        formula="sqrt(0.50*corr_chi2/3 + 0.50*ratio_chi2/2)",
        description="Equal weight to anchored-correlator and pure-ratio sectors after per-panel normalization.",
        evaluator=lambda row: _safe_sqrt(0.50 * _corr_mean(row) + 0.50 * _ratio_mean(row)),
    ),
    Candidate(
        key="ratio_tilted_rms",
        title="Ratio-tilted",
        formula="sqrt(0.35*corr_chi2/3 + 0.65*ratio_chi2/2)",
        description="Favors direction-shape agreement while still keeping some anchored-correlator information.",
        evaluator=lambda row: _safe_sqrt(0.35 * _corr_mean(row) + 0.65 * _ratio_mean(row)),
    ),
    Candidate(
        key="corr_tilted_rms",
        title="Correlator-tilted",
        formula="sqrt(0.65*corr_chi2/3 + 0.35*ratio_chi2/2)",
        description="Favors the anchored-correlator sector while still penalizing directional-ratio mismatches.",
        evaluator=lambda row: _safe_sqrt(0.65 * _corr_mean(row) + 0.35 * _ratio_mean(row)),
    ),
    Candidate(
        key="anchor_shape_rms",
        title="Anchor + shape",
        formula="sqrt(0.50*u/4_chi2 + 0.25*(v/u)_chi2 + 0.25*(w/u)_chi2)",
        description="A nonredundant-ish amplitude-plus-shape blend: one anchored amplitude channel plus two shape ratios.",
        evaluator=lambda row: _safe_sqrt(
            0.50 * _row_value(row, "u/4_chi2")
            + 0.25 * _row_value(row, "(v/4)/(u/4)_chi2")
            + 0.25 * _row_value(row, "(w/4)/(u/4)_chi2")
        ),
    ),
    Candidate(
        key="sector_guard_rms",
        title="Sector guard",
        formula="max(sqrt(corr_chi2/3), sqrt(ratio_chi2/2))",
        description="A conservative guard metric that only looks good when both sectors look good.",
        evaluator=lambda row: max(_safe_sqrt(_corr_mean(row)), _safe_sqrt(_ratio_mean(row))),
    ),
)


def _load_rows(input_path: str) -> list[dict[str, float]]:
    with open(input_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        loaded: list[dict[str, float]] = []
        for raw_row in reader:
            row: dict[str, float] = {}
            for key, value in raw_row.items():
                if key is None:
                    continue
                row[str(key)] = float(value) if key in raw_row and key not in {"coupling_tag"} else value  # type: ignore[assignment]
            loaded.append(row)
    if not loaded:
        raise ValueError(f"no rows found in {input_path}")
    return loaded


def _compute_scores(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    scored_rows: list[dict[str, float]] = []
    for row in rows:
        scored = dict(row)
        scored["corr_panel_mean_chi2"] = _corr_mean(row)
        scored["ratio_panel_mean_chi2"] = _ratio_mean(row)
        for candidate in CANDIDATES:
            scored[candidate.key] = float(candidate.evaluator(scored))
        scored_rows.append(scored)
    return scored_rows


def _grid_extent(values: list[float]) -> tuple[float, float]:
    if len(values) == 1:
        return values[0] - 0.5, values[0] + 0.5
    steps = np.diff(np.asarray(values, dtype=float))
    step = float(np.median(steps))
    return float(values[0] - 0.5 * step), float(values[-1] + 0.5 * step)


def _build_grid(rows: list[dict[str, float]], key: str) -> tuple[list[float], list[float], np.ndarray]:
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


def _best_index(values: np.ndarray) -> tuple[int, int] | None:
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        return None
    masked = np.where(finite_mask, values, np.inf)
    flat_index = int(np.argmin(masked))
    return tuple(int(v) for v in np.unravel_index(flat_index, values.shape))


def _annotation_color(value: float, vmax: float) -> str:
    if vmax <= 0.0:
        return "black"
    return "white" if value >= 0.58 * vmax else "black"


def _base_name(input_path: str) -> str:
    return os.path.splitext(os.path.basename(input_path))[0]


def _write_score_table(output_path: str, rows: list[dict[str, float]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    headers = [
        "coupling_tag",
        "r1",
        "r2",
        "corr_chi2",
        "ratio_chi2",
        "chi2_sum",
        "corr_panel_mean_chi2",
        "ratio_panel_mean_chi2",
        *[candidate.key for candidate in CANDIDATES],
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow(
                [
                    str(row["coupling_tag"]),
                    f"{float(row['r1']):.10f}",
                    f"{float(row['r2']):.10f}",
                    f"{float(row['corr_chi2']):.10f}",
                    f"{float(row['ratio_chi2']):.10f}",
                    f"{float(row['chi2_sum']):.10f}",
                    f"{float(row['corr_panel_mean_chi2']):.10f}",
                    f"{float(row['ratio_panel_mean_chi2']):.10f}",
                    *[f"{float(row[candidate.key]):.10f}" for candidate in CANDIDATES],
                ]
            )


def _write_best_summary(output_path: str, rows: list[dict[str, float]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "title", "best_r1", "best_r2", "best_score", "formula", "description"])
        for candidate in CANDIDATES:
            best_row = min(rows, key=lambda row: float(row[candidate.key]))
            writer.writerow(
                [
                    candidate.key,
                    candidate.title,
                    f"{float(best_row['r1']):.10f}",
                    f"{float(best_row['r2']):.10f}",
                    f"{float(best_row[candidate.key]):.10f}",
                    candidate.formula,
                    candidate.description,
                ]
            )


def _write_candidate_notes(output_path: str, input_path: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("# Mixed score candidates\n\n")
        handle.write(f"Input TSV: {input_path}\n\n")
        handle.write("All candidates are dimensionless RMS-like costs built from z-based chi^2 summaries, so they avoid mixing raw observables with incompatible units.\n\n")
        handle.write("Definitions:\n\n")
        handle.write("- corr_panel_mean_chi2 = corr_chi2 / 3\n")
        handle.write("- ratio_panel_mean_chi2 = ratio_chi2 / 2\n\n")
        for index, candidate in enumerate(CANDIDATES, start=1):
            handle.write(f"{index}. {candidate.title}: `{candidate.formula}`\n")
            handle.write(f"   {candidate.description}\n")


def _plot_candidates(
    *,
    rows: list[dict[str, float]],
    output_path: str,
    figure_title: str,
    truth_r1: float | None,
    truth_r2: float | None,
) -> None:
    grids: dict[str, np.ndarray] = {}
    r1_values: list[float] | None = None
    r2_values: list[float] | None = None
    finite_values: list[np.ndarray] = []
    for candidate in CANDIDATES:
        cand_r1, cand_r2, grid = _build_grid(rows, candidate.key)
        r1_values = cand_r1
        r2_values = cand_r2
        grids[candidate.key] = grid
        finite = grid[np.isfinite(grid)]
        if finite.size:
            finite_values.append(finite)

    if r1_values is None or r2_values is None:
        raise ValueError("failed to construct candidate grids")

    vmax = float(np.max(np.concatenate(finite_values))) if finite_values else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad(color="#e5e7eb")

    x0, x1 = _grid_extent(r1_values)
    y0, y1 = _grid_extent(r2_values)
    extent = [x0, x1, y0, y1]

    fig, axes = plt.subplots(2, 3, figsize=(18.5, 10.0), squeeze=False, constrained_layout=True)
    axes_flat = list(axes.ravel())
    image = None
    plotted_axes = []
    for axis, candidate in zip(axes_flat, CANDIDATES):
        values = np.asarray(grids[candidate.key], dtype=float)
        image = axis.imshow(
            values.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap=cmap,
            norm=norm,
        )
        plotted_axes.append(axis)

        best_index = _best_index(values)
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

        if truth_r1 is not None and truth_r2 is not None:
            axis.scatter(
                [float(truth_r1)],
                [float(truth_r2)],
                marker="x",
                s=90,
                color="black",
                linewidths=1.8,
                zorder=4,
            )

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
                    f"{float(value):.2f}",
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    color=_annotation_color(float(value), vmax),
                )

        axis.set_title(candidate.title, fontsize=11)
        axis.set_xlabel("r1 = k1 / k3")
        axis.set_ylabel("r2 = k2 / k3")
        axis.set_xticks(r1_values)
        axis.set_yticks(r2_values)
        axis.set_xlim(x0, x1)
        axis.set_ylim(y0, y1)

    text_axis = axes_flat[len(CANDIDATES)]
    text_axis.axis("off")
    lines = [
        "Candidate definitions",
        "",
    ]
    for index, candidate in enumerate(CANDIDATES, start=1):
        lines.append(f"{index}. {candidate.title}")
        lines.append(candidate.formula)
        lines.append("")
    text_axis.text(0.0, 1.0, "\n".join(lines), ha="left", va="top", fontsize=10, family="monospace")

    if image is not None:
        fig.colorbar(image, ax=plotted_axes, fraction=0.025, pad=0.025, label="mixed RMS score")
    fig.suptitle(
        f"{figure_title}\nstar = per-panel best, x = reference target",
        fontsize=14,
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute and plot mixed correlator-plus-ratio score candidates.")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input z-score TSV.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Output directory for tables and heatmaps.")
    parser.add_argument("--title", default=None, help="Optional figure title override.")
    parser.add_argument("--truth-r1", type=float, default=None, help="Optional reference r1 marker.")
    parser.add_argument("--truth-r2", type=float, default=None, help="Optional reference r2 marker.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = os.path.abspath(str(args.input))
    output_dir = os.path.abspath(str(args.output_dir))
    os.makedirs(output_dir, exist_ok=True)

    rows = _load_rows(input_path)
    scored_rows = _compute_scores(rows)
    base_name = _base_name(input_path)

    scores_path = os.path.join(output_dir, f"{base_name}_mixed_score_candidates.tsv")
    summary_path = os.path.join(output_dir, f"{base_name}_mixed_score_best_points.tsv")
    notes_path = os.path.join(output_dir, f"{base_name}_mixed_score_candidates.md")
    heatmap_path = os.path.join(output_dir, f"{base_name}_mixed_score_candidates_heatmaps.png")

    _write_score_table(scores_path, scored_rows)
    _write_best_summary(summary_path, scored_rows)
    _write_candidate_notes(notes_path, input_path)
    _plot_candidates(
        rows=scored_rows,
        output_path=heatmap_path,
        figure_title=str(args.title) if args.title else f"{base_name} mixed-score candidate heatmaps",
        truth_r1=(float(args.truth_r1) if args.truth_r1 is not None else None),
        truth_r2=(float(args.truth_r2) if args.truth_r2 is not None else None),
    )

    print(f"wrote {scores_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {notes_path}")
    print(f"wrote {heatmap_path}")


if __name__ == "__main__":
    main()