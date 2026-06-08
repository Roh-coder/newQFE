#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent

DEFAULT_GRID = (
    RESPONSIBLE_ROOT
    / "results"
    / "geometry_match_reweight_interp_iso111_grid5x5_donor3x3_sizes8_12_16_24_target32_20260604"
    / "grid_scores.tsv"
)
DEFAULT_OUTPUT_DIR = DEFAULT_GRID.parent / "gradient_field"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot direct vs histogram-reweighted cost-surface gradient fields from a grid_scores.tsv file."
    )
    parser.add_argument("--grid-tsv", default=str(DEFAULT_GRID))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--title", default="Iso111 cost gradient field from direct vs reweighted surfaces")
    return parser.parse_args()


def _load_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "r1",
            "r2",
            "donor_r1",
            "donor_r2",
            "direct_z2_sum",
            "reweighted_z2_sum",
            "delta_z2_sum",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"missing required columns in {path}")
        for raw in reader:
            rows.append({key: float(raw[key]) for key in required})
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def _build_field(rows: list[dict[str, float]], key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r1_values = np.asarray(sorted({row["r1"] for row in rows}), dtype=float)
    r2_values = np.asarray(sorted({row["r2"] for row in rows}), dtype=float)
    field = np.full((r1_values.size, r2_values.size), np.nan, dtype=float)
    r1_index = {float(value): idx for idx, value in enumerate(r1_values)}
    r2_index = {float(value): idx for idx, value in enumerate(r2_values)}
    for row in rows:
        field[r1_index[float(row["r1"])]][r2_index[float(row["r2"])] ] = float(row[key])
    if np.any(~np.isfinite(field)):
        raise ValueError(f"incomplete rectangular grid for {key}")
    return r1_values, r2_values, field


def _compute_gradient(r1_values: np.ndarray, r2_values: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d_dr1, d_dr2 = np.gradient(field, r1_values, r2_values, edge_order=2)
    magnitude = np.hypot(d_dr1, d_dr2)
    return d_dr1, d_dr2, magnitude


def _normalized_vectors(dx: np.ndarray, dy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mag = np.hypot(dx, dy)
    out_x = np.zeros_like(dx)
    out_y = np.zeros_like(dy)
    mask = mag > 0.0
    out_x[mask] = dx[mask] / mag[mask]
    out_y[mask] = dy[mask] / mag[mask]
    return out_x, out_y


def _write_gradient_tsv(
    path: Path,
    r1_values: np.ndarray,
    r2_values: np.ndarray,
    direct_field: np.ndarray,
    reweight_field: np.ndarray,
    direct_grad_r1: np.ndarray,
    direct_grad_r2: np.ndarray,
    direct_mag: np.ndarray,
    reweight_grad_r1: np.ndarray,
    reweight_grad_r2: np.ndarray,
    reweight_mag: np.ndarray,
) -> None:
    lines = [
        "r1\tr2\tdirect_z2_sum\treweighted_z2_sum\tdelta_z2_sum\tdirect_dC_dr1\tdirect_dC_dr2\tdirect_grad_norm\treweighted_dC_dr1\treweighted_dC_dr2\treweighted_grad_norm\tgrad_delta_dr1\tgrad_delta_dr2\tgrad_delta_norm"
    ]
    for i, r1 in enumerate(r1_values):
        for j, r2 in enumerate(r2_values):
            delta_r1 = float(reweight_grad_r1[i, j] - direct_grad_r1[i, j])
            delta_r2 = float(reweight_grad_r2[i, j] - direct_grad_r2[i, j])
            delta_norm = math.hypot(delta_r1, delta_r2)
            lines.append(
                "\t".join(
                    [
                        f"{float(r1):.6f}",
                        f"{float(r2):.6f}",
                        f"{float(direct_field[i, j]):.16e}",
                        f"{float(reweight_field[i, j]):.16e}",
                        f"{float(reweight_field[i, j] - direct_field[i, j]):.16e}",
                        f"{float(direct_grad_r1[i, j]):.16e}",
                        f"{float(direct_grad_r2[i, j]):.16e}",
                        f"{float(direct_mag[i, j]):.16e}",
                        f"{float(reweight_grad_r1[i, j]):.16e}",
                        f"{float(reweight_grad_r2[i, j]):.16e}",
                        f"{float(reweight_mag[i, j]):.16e}",
                        f"{delta_r1:.16e}",
                        f"{delta_r2:.16e}",
                        f"{delta_norm:.16e}",
                    ]
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plot_panel(
    ax: plt.Axes,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    field: np.ndarray,
    grad_r1: np.ndarray,
    grad_r2: np.ndarray,
    grad_mag: np.ndarray,
    title: str,
) -> None:
    mesh = ax.pcolormesh(x_grid, y_grid, field.T, shading="auto", cmap="magma")
    u, v = _normalized_vectors(grad_r1, grad_r2)
    quiver = ax.quiver(
        x_grid,
        y_grid,
        u.T,
        v.T,
        grad_mag.T,
        cmap="viridis",
        pivot="mid",
        scale=18,
        width=0.008,
    )
    ax.set_title(title)
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    plt.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04, label="cost")
    plt.colorbar(quiver, ax=ax, fraction=0.046, pad=0.10, label="|grad cost|")


def main() -> None:
    args = parse_args()
    grid_path = Path(args.grid_tsv).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = _load_rows(grid_path)
    r1_values, r2_values, direct_field = _build_field(rows, "direct_z2_sum")
    _, _, reweight_field = _build_field(rows, "reweighted_z2_sum")

    direct_grad_r1, direct_grad_r2, direct_mag = _compute_gradient(r1_values, r2_values, direct_field)
    reweight_grad_r1, reweight_grad_r2, reweight_mag = _compute_gradient(r1_values, r2_values, reweight_field)
    delta_field = reweight_field - direct_field
    delta_grad_norm = np.hypot(reweight_grad_r1 - direct_grad_r1, reweight_grad_r2 - direct_grad_r2)

    x_grid, y_grid = np.meshgrid(r1_values, r2_values, indexing="ij")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    fig.suptitle(str(args.title), fontsize=14)

    _plot_panel(
        axes[0, 0],
        x_grid,
        y_grid,
        direct_field,
        direct_grad_r1,
        direct_grad_r2,
        direct_mag,
        "Direct cost gradient field",
    )
    _plot_panel(
        axes[0, 1],
        x_grid,
        y_grid,
        reweight_field,
        reweight_grad_r1,
        reweight_grad_r2,
        reweight_mag,
        "Reweighted cost gradient field",
    )

    delta_mesh = axes[1, 0].pcolormesh(x_grid, y_grid, delta_field.T, shading="auto", cmap="coolwarm")
    axes[1, 0].set_title("Reweighted - direct cost")
    axes[1, 0].set_xlabel("r1")
    axes[1, 0].set_ylabel("r2")
    plt.colorbar(delta_mesh, ax=axes[1, 0], fraction=0.046, pad=0.04, label="delta cost")

    grad_mesh = axes[1, 1].pcolormesh(x_grid, y_grid, delta_grad_norm.T, shading="auto", cmap="cividis")
    axes[1, 1].set_title("Gradient mismatch norm")
    axes[1, 1].set_xlabel("r1")
    axes[1, 1].set_ylabel("r2")
    plt.colorbar(grad_mesh, ax=axes[1, 1], fraction=0.046, pad=0.04, label="|grad_rw - grad_direct|")

    figure_path = output_dir / "cost_gradient_field.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    _write_gradient_tsv(
        output_dir / "cost_gradient_field.tsv",
        r1_values,
        r2_values,
        direct_field,
        reweight_field,
        direct_grad_r1,
        direct_grad_r2,
        direct_mag,
        reweight_grad_r1,
        reweight_grad_r2,
        reweight_mag,
    )

    summary = {
        "grid_tsv": str(grid_path),
        "output_dir": str(output_dir),
        "r1_values": [float(value) for value in r1_values],
        "r2_values": [float(value) for value in r2_values],
        "direct_cost_range": [float(np.min(direct_field)), float(np.max(direct_field))],
        "reweighted_cost_range": [float(np.min(reweight_field)), float(np.max(reweight_field))],
        "delta_cost_range": [float(np.min(delta_field)), float(np.max(delta_field))],
        "max_gradient_mismatch_norm": float(np.max(delta_grad_norm)),
        "mean_gradient_mismatch_norm": float(np.mean(delta_grad_norm)),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()