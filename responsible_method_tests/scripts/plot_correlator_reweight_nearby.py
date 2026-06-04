#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render summary plots for the nearby-point correlator reweighting test."
    )
    parser.add_argument("summary_json", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for output PNGs. Defaults to the summary.json directory.",
    )
    return parser.parse_args()


def _load_summary(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _direct_lookup(direct_rows: list[dict]) -> dict[str, dict]:
    return {str(row["label"]): row for row in direct_rows}


def _matrix(rows: list[dict], labels: list[str], key: str) -> np.ndarray:
    grid = np.full((len(labels), len(labels)), np.nan, dtype=float)
    index = {label: idx for idx, label in enumerate(labels)}
    for row in rows:
        grid[index[str(row["source"])]][index[str(row["target"])] ] = float(row[key])
    return grid


def _matrix_text(ax: plt.Axes, grid: np.ndarray, fmt: str, text_color: str = "black") -> None:
    for iy in range(grid.shape[0]):
        for ix in range(grid.shape[1]):
            value = grid[iy, ix]
            if math.isfinite(float(value)):
                ax.text(ix, iy, format(float(value), fmt), ha="center", va="center", color=text_color, fontsize=9)


def _save_overview(summary: dict, output_path: Path) -> None:
    direct_rows = list(summary["direct"])
    labels = [str(row["label"]) for row in direct_rows]
    xs = np.asarray([float(row["r1"]) for row in direct_rows], dtype=float)
    ys = np.asarray([float(row["r2"]) for row in direct_rows], dtype=float)
    conn = np.asarray([float(row["connected"]) for row in direct_rows], dtype=float)
    sigma = np.asarray([float(row["sigma"]) for row in direct_rows], dtype=float)
    wall = np.asarray([float(row["wall_seconds"]) for row in direct_rows], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.scatter(xs, ys, s=95, c=conn, cmap="viridis", edgecolor="black", linewidth=0.8)
    for row in direct_rows:
        ax.annotate(
            str(row["label"]),
            (float(row["r1"]), float(row["r2"])),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9,
        )
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title("Nearby Points In Parameter Space")
    ax.grid(alpha=0.25)

    ax = axes[1]
    xpos = np.arange(len(labels), dtype=float)
    ax.errorbar(xpos, conn, yerr=sigma, fmt="o", capsize=4, color="#1d4ed8", markersize=7)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Connected correlator")
    ax.set_title("Direct MC Result At Displacement")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[2]
    ax.bar(xpos, wall, color="#0f766e")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Wall seconds")
    ax.set_title("Per-Point Runtime")
    ax.grid(axis="y", alpha=0.25)

    disp = summary["displacement"]
    lattice = summary["lattice"]
    fig.suptitle(
        "Nearby-point correlator reweighting control"
        f" | lattice=({lattice['Lx']},{lattice['Ly']},{lattice['Tx']},{lattice['Ty']})"
        f" | displacement=({disp['m']},{disp['n']})"
        f" | n_traj={summary['n_traj']} n_therm={summary['n_therm']} n_skip={summary['n_skip']}",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _save_matrix(summary: dict, output_path: Path) -> None:
    direct_rows = list(summary["direct"])
    labels = [str(row["label"]) for row in direct_rows]
    reweight_rows = list(summary["reweight"])
    z_grid = _matrix(reweight_rows, labels, "z")
    neff_grid = _matrix(reweight_rows, labels, "n_eff_fraction")
    delta_grid = _matrix(reweight_rows, labels, "delta")

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)

    vmax = max(1.0, float(np.nanmax(np.abs(z_grid))))
    im0 = axes[0].imshow(z_grid, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    axes[0].set_title("Prediction vs Direct: z-score")
    axes[0].set_xticks(np.arange(len(labels)))
    axes[0].set_yticks(np.arange(len(labels)))
    axes[0].set_xticklabels(labels, rotation=18, ha="right")
    axes[0].set_yticklabels(labels)
    axes[0].set_xlabel("target")
    axes[0].set_ylabel("source")
    _matrix_text(axes[0], z_grid, ".2f")
    fig.colorbar(im0, ax=axes[0], shrink=0.9)

    im1 = axes[1].imshow(neff_grid, cmap="YlGn", vmin=float(np.nanmin(neff_grid)), vmax=1.0)
    axes[1].set_title("Effective sample fraction")
    axes[1].set_xticks(np.arange(len(labels)))
    axes[1].set_yticks(np.arange(len(labels)))
    axes[1].set_xticklabels(labels, rotation=18, ha="right")
    axes[1].set_yticklabels(labels)
    axes[1].set_xlabel("target")
    axes[1].set_ylabel("source")
    _matrix_text(axes[1], neff_grid, ".4f")
    fig.colorbar(im1, ax=axes[1], shrink=0.9)

    max_abs_delta = max(float(np.nanmax(np.abs(delta_grid))), 1.0e-12)
    im2 = axes[2].imshow(delta_grid, cmap="PiYG", vmin=-max_abs_delta, vmax=max_abs_delta)
    axes[2].set_title("Prediction - direct")
    axes[2].set_xticks(np.arange(len(labels)))
    axes[2].set_yticks(np.arange(len(labels)))
    axes[2].set_xticklabels(labels, rotation=18, ha="right")
    axes[2].set_yticklabels(labels)
    axes[2].set_xlabel("target")
    axes[2].set_ylabel("source")
    _matrix_text(axes[2], delta_grid, ".4e")
    fig.colorbar(im2, ax=axes[2], shrink=0.9)

    fig.suptitle("Nearby-point correlator reweighting agreement matrices", fontsize=12)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    summary_path = args.summary_json.resolve()
    summary = _load_summary(summary_path)
    output_dir = (args.output_dir.resolve() if args.output_dir else summary_path.parent)
    output_dir.mkdir(parents=True, exist_ok=True)

    _save_overview(summary, output_dir / "correlator_reweight_nearby_overview.png")
    _save_matrix(summary, output_dir / "correlator_reweight_nearby_matrices.png")


if __name__ == "__main__":
    main()