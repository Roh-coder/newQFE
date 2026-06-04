#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render finite-size plots for the nearby-point correlator reweighting ladder."
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


def _sorted_runs(summary: dict) -> list[dict]:
    runs = list(summary.get("runs", []))
    if not runs:
        raise ValueError("expected a multi-size summary with a non-empty 'runs' array")
    return sorted(runs, key=lambda run: int(run["size"]))


def _point_labels(summary: dict, runs: list[dict]) -> list[str]:
    if "points" in summary:
        return [str(row["label"]) for row in summary["points"]]
    return [str(row["label"]) for row in runs[0]["direct"]]


def _direct_row(run: dict, label: str) -> dict:
    for row in run["direct"]:
        if str(row["label"]) == label:
            return row
    raise KeyError(f"missing direct row for label={label}")


def _reweight_row(run: dict, source: str, target: str) -> dict:
    for row in run["reweight"]:
        if str(row["source"]) == source and str(row["target"]) == target:
            return row
    raise KeyError(f"missing reweight row for {source}->{target}")


def _x_values(runs: list[dict]) -> tuple[np.ndarray, list[str]]:
    sizes = [int(run["size"]) for run in runs]
    x = np.asarray(sizes, dtype=float)
    tick_labels = [str(size) for size in sizes]
    return x, tick_labels


def _save_targets(summary: dict, runs: list[dict], output_path: Path) -> None:
    labels = _point_labels(summary, runs)
    x, tick_labels = _x_values(runs)
    colors = {
        labels[0]: "#111827",
        labels[1]: "#1d4ed8",
        labels[2]: "#b91c1c",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)
    ax_direct = axes[0, 0]
    target_axes = {
        labels[0]: axes[0, 1],
        labels[1]: axes[1, 0],
        labels[2]: axes[1, 1],
    }

    for label in labels:
        ys = np.asarray([float(_direct_row(run, label)["connected"]) for run in runs], dtype=float)
        sigmas = np.asarray([float(_direct_row(run, label)["sigma"]) for run in runs], dtype=float)
        ax_direct.errorbar(
            x,
            ys,
            yerr=sigmas,
            marker="o",
            capsize=3,
            linewidth=1.4,
            color=colors[label],
            label=label,
        )
    ax_direct.set_title("Direct MC finite-size trend")
    ax_direct.set_xlabel("L")
    ax_direct.set_ylabel("Connected correlator")
    ax_direct.set_xticks(x)
    ax_direct.set_xticklabels(tick_labels)
    ax_direct.grid(alpha=0.25)
    ax_direct.legend(frameon=False, fontsize=9)

    markers = {labels[0]: "s", labels[1]: "^", labels[2]: "D"}
    for target in labels:
        ax = target_axes[target]
        direct_y = np.asarray([float(_direct_row(run, target)["connected"]) for run in runs], dtype=float)
        direct_sigma = np.asarray([float(_direct_row(run, target)["sigma"]) for run in runs], dtype=float)
        ax.errorbar(
            x,
            direct_y,
            yerr=direct_sigma,
            marker="o",
            capsize=3,
            linewidth=1.6,
            color="#000000",
            label=f"direct {target}",
        )
        for source in labels:
            if source == target:
                continue
            pred_y = np.asarray(
                [float(_reweight_row(run, source, target)["predicted_connected"]) for run in runs],
                dtype=float,
            )
            pred_sigma = np.asarray(
                [float(_reweight_row(run, source, target)["predicted_sigma"]) for run in runs],
                dtype=float,
            )
            ax.errorbar(
                x,
                pred_y,
                yerr=pred_sigma,
                marker=markers[source],
                capsize=3,
                linewidth=1.2,
                linestyle="--",
                color=colors[source],
                label=f"rewt {source}->{target}",
            )
        ax.set_title(f"Target: {target}")
        ax.set_xlabel("L")
        ax.set_ylabel("Connected correlator")
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels)
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)

    disp_fraction = summary.get("disp_fraction")
    disp_label = (
        f"disp_fraction={disp_fraction:.3f}" if isinstance(disp_fraction, (int, float)) else "size-specific displacement"
    )
    fig.suptitle(
        "Finite-size correlator reweighting at three nearby coupling points"
        f" | {disp_label} | n_traj={summary['n_traj']} n_therm={summary['n_therm']} n_skip={summary['n_skip']}",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _save_quality(summary: dict, runs: list[dict], output_path: Path) -> None:
    labels = _point_labels(summary, runs)
    sizes = [int(run["size"]) for run in runs]
    pairs = [(source, target) for source in labels for target in labels if source != target]
    z_grid = np.asarray(
        [[float(_reweight_row(run, source, target)["z"]) for run in runs] for source, target in pairs],
        dtype=float,
    )
    neff_grid = np.asarray(
        [[float(_reweight_row(run, source, target)["n_eff_fraction"]) for run in runs] for source, target in pairs],
        dtype=float,
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.0), constrained_layout=True)

    vmax = max(1.0, float(np.nanmax(np.abs(z_grid))))
    im0 = axes[0].imshow(z_grid, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    axes[0].set_title("z-score by source/target and size")
    axes[0].set_xticks(np.arange(len(sizes)))
    axes[0].set_xticklabels([str(size) for size in sizes])
    axes[0].set_yticks(np.arange(len(pairs)))
    axes[0].set_yticklabels([f"{source}->{target}" for source, target in pairs])
    axes[0].set_xlabel("L")
    axes[0].set_ylabel("reweight pair")
    for iy in range(z_grid.shape[0]):
        for ix in range(z_grid.shape[1]):
            axes[0].text(ix, iy, f"{z_grid[iy, ix]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im0, ax=axes[0], shrink=0.9)

    im1 = axes[1].imshow(neff_grid, aspect="auto", cmap="YlGn", vmin=float(np.nanmin(neff_grid)), vmax=1.0)
    axes[1].set_title("Effective sample fraction by source/target and size")
    axes[1].set_xticks(np.arange(len(sizes)))
    axes[1].set_xticklabels([str(size) for size in sizes])
    axes[1].set_yticks(np.arange(len(pairs)))
    axes[1].set_yticklabels([f"{source}->{target}" for source, target in pairs])
    axes[1].set_xlabel("L")
    axes[1].set_ylabel("reweight pair")
    for iy in range(neff_grid.shape[0]):
        for ix in range(neff_grid.shape[1]):
            axes[1].text(ix, iy, f"{neff_grid[iy, ix]:.4f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im1, ax=axes[1], shrink=0.9)

    fig.suptitle("Finite-size reweighting quality diagnostics", fontsize=12)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    summary_path = args.summary_json.resolve()
    summary = _load_summary(summary_path)
    runs = _sorted_runs(summary)
    output_dir = args.output_dir.resolve() if args.output_dir else summary_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    _save_targets(summary, runs, output_dir / "correlator_reweight_fss_targets.png")
    _save_quality(summary, runs, output_dir / "correlator_reweight_fss_quality.png")


if __name__ == "__main__":
    main()