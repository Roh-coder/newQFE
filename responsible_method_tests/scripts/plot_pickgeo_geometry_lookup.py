#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Polygon


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "pickgeo10_geometry_lookup.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the normalized triangle shapes in a pickGeo benchmark config.",
    )
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        help="Path to a responsible_method_tests pickGeo config JSON.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Path to the output PNG.",
    )
    parser.add_argument(
        "--table-output",
        default=None,
        help="Optional path for a TSV lookup table. Defaults to <output>.tsv.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_target_ratio(benchmark_id: str) -> tuple[int, int, int]:
    token = benchmark_id.rsplit("_", 1)[-1]
    if len(token) != 3 or not token.isdigit():
        raise ValueError(f"Could not extract 3-digit ratio from benchmark id: {benchmark_id}")
    return tuple(int(ch) for ch in token)


def sorted_cycle_lengths(geometry: tuple[int, int, int, int]) -> tuple[float, float, float]:
    lx, ly, tx, ty = geometry
    omega = complex(0.5, math.sqrt(3.0) / 2.0)
    cycles = (
        lx + ty * omega,
        tx - ly * omega,
        -(lx + tx) + (ly - ty) * omega,
    )
    lengths = sorted(abs(value) for value in cycles)
    return (lengths[0], lengths[1], lengths[2])


def triangle_vertices(a_over_c: float, b_over_c: float) -> list[tuple[float, float]]:
    x_coord = (b_over_c * b_over_c - a_over_c * a_over_c + 1.0) / 2.0
    height_sq = max(b_over_c * b_over_c - x_coord * x_coord, 0.0)
    return [(0.0, 0.0), (1.0, 0.0), (x_coord, math.sqrt(height_sq))]


def shape_type(spread: float, scalene_gap: float) -> str:
    if spread < 1.0e-9:
        return "equilateral"
    if scalene_gap < 1.0e-9:
        return "isosceles"
    return "scalene"


def build_rows(config: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for benchmark in config.get("benchmarks", []):
        benchmark_id = str(benchmark["id"])
        geometry = tuple(int(value) for value in benchmark["modular_target_geometry"])
        target_ratio = parse_target_ratio(benchmark_id)
        target_sorted = sorted(float(value) for value in target_ratio)
        target_norm = [value / target_sorted[-1] for value in target_sorted]

        side_a, side_b, side_c = sorted_cycle_lengths(geometry)
        a_over_c = side_a / side_c
        b_over_c = side_b / side_c
        spread = 1.0 - a_over_c
        scalene_gap = min(side_b - side_a, side_c - side_b) / side_c
        ratio_error = max(
            abs(a_over_c - target_norm[0]),
            abs(b_over_c - target_norm[1]),
            0.0,
        )

        rows.append(
            {
                "benchmark_id": benchmark_id,
                "label": "-".join(str(value) for value in target_ratio),
                "digits": "".join(str(value) for value in target_ratio),
                "geometry": geometry,
                "a": side_a,
                "b": side_b,
                "c": side_c,
                "a_over_c": a_over_c,
                "b_over_c": b_over_c,
                "spread": spread,
                "scalene_gap": scalene_gap,
                "ratio_error": ratio_error,
                "shape_type": shape_type(spread, scalene_gap),
            }
        )
    return rows


def write_lookup_table(rows: list[dict[str, object]], table_output: Path) -> None:
    table_output.parent.mkdir(parents=True, exist_ok=True)
    with table_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "benchmark_id",
                "target_ratio",
                "Lx",
                "Ly",
                "Tx",
                "Ty",
                "short_over_long",
                "mid_over_long",
                "spread",
                "scalene_gap",
                "ratio_error",
                "shape_type",
            ]
        )
        for row in rows:
            lx, ly, tx, ty = row["geometry"]
            writer.writerow(
                [
                    row["benchmark_id"],
                    row["label"],
                    lx,
                    ly,
                    tx,
                    ty,
                    f"{row['a_over_c']:.6f}",
                    f"{row['b_over_c']:.6f}",
                    f"{row['spread']:.6f}",
                    f"{row['scalene_gap']:.6f}",
                    f"{row['ratio_error']:.6f}",
                    row["shape_type"],
                ]
            )


def plot_lookup(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    scalene_max = max(float(row["scalene_gap"]) for row in rows) if rows else 1.0
    norm = colors.Normalize(vmin=0.0, vmax=scalene_max if scalene_max > 0.0 else 1.0)
    cmap = matplotlib.colormaps["viridis"]

    figure = plt.figure(figsize=(16, 11), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, width_ratios=[1.15, 0.85], height_ratios=[1.05, 1.55])

    ax_shape = figure.add_subplot(grid[0, 0])
    ax_table = figure.add_subplot(grid[0, 1])
    triangles_grid = grid[1, :].subgridspec(2, 5, wspace=0.24, hspace=0.38)

    x_values = [float(row["b_over_c"]) for row in rows]
    y_values = [float(row["a_over_c"]) for row in rows]
    colors_values = [float(row["scalene_gap"]) for row in rows]
    sizes = [150.0 + 550.0 * float(row["spread"]) for row in rows]

    scatter = ax_shape.scatter(
        x_values,
        y_values,
        c=colors_values,
        s=sizes,
        cmap=cmap,
        norm=norm,
        edgecolors="black",
        linewidths=0.8,
        alpha=0.95,
        zorder=3,
    )
    ax_shape.plot([0.6, 1.0], [0.6, 1.0], linestyle="--", color="#666666", linewidth=1.0)
    ax_shape.axvline(1.0, linestyle="--", color="#666666", linewidth=1.0)

    for row in rows:
        ax_shape.annotate(
            str(row["digits"]),
            (float(row["b_over_c"]), float(row["a_over_c"])),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=10,
            weight="bold",
        )

    ax_shape.set_xlim(0.82, 1.01)
    ax_shape.set_ylim(0.62, 1.01)
    ax_shape.set_aspect("equal", adjustable="box")
    ax_shape.set_xlabel("middle side / longest side")
    ax_shape.set_ylabel("shortest side / longest side")
    ax_shape.set_title("Normalized shape space")
    ax_shape.text(
        0.823,
        1.002,
        "Closer to y=x or x=1 means closer to an isosceles family.\n"
        "Marker size tracks total anisotropy 1 - A/C.",
        va="top",
        ha="left",
        fontsize=9,
    )

    colorbar = figure.colorbar(scatter, ax=ax_shape, fraction=0.046, pad=0.04)
    colorbar.set_label("scalene gap = min(B-A, C-B) / C")

    ax_table.axis("off")
    table_rows = sorted(
        rows,
        key=lambda row: (-float(row["scalene_gap"]), -float(row["spread"]), str(row["digits"])),
    )
    table_lines = [
        "id   A/C    B/C    spread  scalene  type",
        "---  -----  -----  ------  -------  ----------",
    ]
    for row in table_rows:
        table_lines.append(
            f"{str(row['digits']):>3}  "
            f"{float(row['a_over_c']):.3f}  "
            f"{float(row['b_over_c']):.3f}  "
            f"{float(row['spread']):.3f}  "
            f"{float(row['scalene_gap']):.3f}  "
            f"{str(row['shape_type'])}"
        )
    ax_table.text(
        0.0,
        1.0,
        "Lookup table\n\n" + "\n".join(table_lines) + "\n\n"
        "spread   = 1 - A/C\n"
        "scalene  = min(B-A, C-B) / C\n"
        "A <= B <= C are the sorted side lengths.",
        va="top",
        ha="left",
        family="monospace",
        fontsize=10,
    )

    for axis, row in zip((figure.add_subplot(triangles_grid[index // 5, index % 5]) for index in range(len(rows))), rows):
        color = cmap(norm(float(row["scalene_gap"])))
        vertices = triangle_vertices(float(row["a_over_c"]), float(row["b_over_c"]))
        polygon = Polygon(vertices, closed=True, facecolor=color, edgecolor=color, alpha=0.35, linewidth=2.0)
        axis.add_patch(polygon)
        axis.plot([point[0] for point in vertices + [vertices[0]]], [point[1] for point in vertices + [vertices[0]]], color=color, linewidth=2.0)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(-0.06, 1.06)
        axis.set_ylim(-0.08, 0.92)
        axis.axis("off")
        axis.set_title(str(row["digits"]), fontsize=12, weight="bold", pad=4)
        axis.text(
            0.5,
            -0.09,
            f"A/C={float(row['a_over_c']):.3f}  B/C={float(row['b_over_c']):.3f}\n"
            f"scalene={float(row['scalene_gap']):.3f}",
            ha="center",
            va="top",
            transform=axis.transAxes,
            fontsize=9,
        )

    figure.suptitle(
        "PickGeo10 tested geometries\nShapes normalized by the longest side",
        fontsize=17,
        weight="bold",
    )
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    table_output = Path(args.table_output).expanduser().resolve() if args.table_output else output_path.with_suffix(".tsv")

    config = load_json(config_path)
    rows = build_rows(config)
    if not rows:
        raise ValueError(f"No benchmarks found in config: {config_path}")

    write_lookup_table(rows, table_output)
    plot_lookup(rows, output_path)

    print(f"Wrote plot: {output_path}")
    print(f"Wrote table: {table_output}")


if __name__ == "__main__":
    main()