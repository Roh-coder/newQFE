#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from align_manifolds import build_benchmark_comparison


def _unit_cell_outline(tau_real: float, tau_imag: float) -> tuple[np.ndarray, np.ndarray]:
    x_value = np.asarray([0.0, 1.0, 1.0 + tau_real, tau_real, 0.0], dtype=float)
    y_value = np.asarray([0.0, 0.0, tau_imag, tau_imag, 0.0], dtype=float)
    return x_value, y_value


def _default_output_path(benchmark_manifest_path: str, benchmark_id: str) -> str:
    base_dir = os.path.dirname(os.path.abspath(benchmark_manifest_path))
    return os.path.join(base_dir, f"{benchmark_id}_continuum_difference.png")


def _plot_difference(metrics: dict, plot_data: dict, *, output_path: str) -> None:
    benchmark_id = str(metrics["benchmark_id"])
    tau_real = float(plot_data["tau_real"])
    tau_imag = float(plot_data["tau_imag"])
    outline_x, outline_y = _unit_cell_outline(tau_real, tau_imag)

    X_grid = np.asarray(plot_data["X_grid"], dtype=float)
    Y_grid = np.asarray(plot_data["Y_grid"], dtype=float)
    twisted_grid = np.asarray(plot_data["twisted_grid"], dtype=float)
    untwisted_grid = np.asarray(plot_data["untwisted_grid"], dtype=float)
    diff_grid = twisted_grid - untwisted_grid

    z_on_untwisted = np.asarray(plot_data["z_on_untwisted"], dtype=float)
    z_on_twisted = np.asarray(plot_data["z_on_twisted"], dtype=float)
    untwisted_sample_x = np.asarray(plot_data["untwisted_sample_x"], dtype=float)
    untwisted_sample_y = np.asarray(plot_data["untwisted_sample_y"], dtype=float)
    twisted_sample_x = np.asarray(plot_data["twisted_sample_x"], dtype=float)
    twisted_sample_y = np.asarray(plot_data["twisted_sample_y"], dtype=float)

    finite_common = np.isfinite(diff_grid)
    diff_limit = float(np.nanmax(np.abs(diff_grid[finite_common]))) if np.any(finite_common) else 1.0
    diff_limit = max(diff_limit, 1.0e-6)

    z_limit = 0.0
    if np.any(np.isfinite(z_on_untwisted)):
        z_limit = max(z_limit, float(np.nanmax(np.abs(z_on_untwisted[np.isfinite(z_on_untwisted)]))))
    if np.any(np.isfinite(z_on_twisted)):
        z_limit = max(z_limit, float(np.nanmax(np.abs(z_on_twisted[np.isfinite(z_on_twisted)]))))
    z_limit = max(z_limit, 1.0)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    ax_grid, ax_untwisted, ax_twisted = axes
    fig.suptitle(
        f"{benchmark_id}: post-fit twisted - untwisted differences",
        fontsize=14,
        y=0.98,
    )

    grid_handle = ax_grid.pcolormesh(
        X_grid,
        Y_grid,
        diff_grid,
        shading="auto",
        cmap="coolwarm",
        vmin=-diff_limit,
        vmax=diff_limit,
    )
    ax_grid.plot(outline_x, outline_y, color="black", linewidth=1.2)
    ax_grid.set_title(
        "Common-grid continuum difference\n"
        f"RMS={metrics['common_grid']['rms_abs']:.5f}, max={metrics['common_grid']['max_abs']:.5f}",
        fontsize=11,
    )
    ax_grid.set_xlabel("nu_real")
    ax_grid.set_ylabel("nu_imag")
    ax_grid.set_aspect("equal", adjustable="box")
    ax_grid.grid(True, alpha=0.25)
    fig.colorbar(grid_handle, ax=ax_grid, fraction=0.046, pad=0.04, label="twisted - untwisted")

    untwisted_handle = ax_untwisted.scatter(
        untwisted_sample_x,
        untwisted_sample_y,
        c=z_on_untwisted,
        s=18,
        cmap="coolwarm",
        vmin=-z_limit,
        vmax=z_limit,
    )
    ax_untwisted.plot(outline_x, outline_y, color="black", linewidth=1.2)
    ax_untwisted.set_title(
        "Sampled on untwisted fitted points\n"
        f"RMS z={metrics['twisted_on_untwisted']['rms_z']:.3f}, max |z|={metrics['twisted_on_untwisted']['max_abs_z']:.3f}",
        fontsize=11,
    )
    ax_untwisted.set_xlabel("nu_real")
    ax_untwisted.set_ylabel("nu_imag")
    ax_untwisted.set_aspect("equal", adjustable="box")
    ax_untwisted.grid(True, alpha=0.25)

    twisted_handle = ax_twisted.scatter(
        twisted_sample_x,
        twisted_sample_y,
        c=z_on_twisted,
        s=18,
        cmap="coolwarm",
        vmin=-z_limit,
        vmax=z_limit,
    )
    ax_twisted.plot(outline_x, outline_y, color="black", linewidth=1.2)
    ax_twisted.set_title(
        "Sampled on twisted fitted points\n"
        f"RMS z={metrics['untwisted_on_twisted']['rms_z']:.3f}, max |z|={metrics['untwisted_on_twisted']['max_abs_z']:.3f}",
        fontsize=11,
    )
    ax_twisted.set_xlabel("nu_real")
    ax_twisted.set_ylabel("nu_imag")
    ax_twisted.set_aspect("equal", adjustable="box")
    ax_twisted.grid(True, alpha=0.25)
    fig.colorbar(twisted_handle, ax=[ax_untwisted, ax_twisted], fraction=0.046, pad=0.04, label="cross-method z")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot post-fit differences between twisted and untwisted continuum manifolds."
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--grid-size", type=int, default=180, help="Interpolation grid size in each cell direction")
    parser.add_argument(
        "--comparison-mode",
        choices=["symmetric", "twisted_reference"],
        default="symmetric",
        help="Comparison mode passed through to the continuum comparer",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    benchmark_manifest_path = os.path.abspath(args.benchmark_manifest)
    metrics, plot_data = build_benchmark_comparison(
        benchmark_manifest_path,
        grid_size=max(int(args.grid_size), 20),
        comparison_mode=str(args.comparison_mode),
    )
    output_path = os.path.abspath(args.output) if args.output else _default_output_path(
        benchmark_manifest_path,
        str(metrics["benchmark_id"]),
    )
    _plot_difference(metrics, plot_data, output_path=output_path)
    print(f"wrote {output_path}")
    print(f"common_grid_rms_abs={metrics['common_grid']['rms_abs']:.6f}")
    print(f"twisted_on_untwisted_rms_z={metrics['twisted_on_untwisted']['rms_z']:.6f}")
    print(f"untwisted_on_twisted_rms_z={metrics['untwisted_on_twisted']['rms_z']:.6f}")


if __name__ == "__main__":
    main()