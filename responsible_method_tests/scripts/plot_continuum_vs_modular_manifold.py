#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from plot_fss_interpolated_manifolds import _interpolate_periodic_manifold, _load_dat_rows, _load_json, _unit_cell_outline


def _lookup_continuum(rows: list[dict[str, Any]]) -> dict[int, float]:
    out: dict[int, float] = {}
    for row in rows:
        if int(row["is_origin"]) != 0:
            continue
        value = float(row["A"])
        if np.isfinite(value):
            out[int(row["point_id"])] = value
    return out


def _lookup_modular(rows: list[dict[str, Any]]) -> dict[int, float]:
    out: dict[int, float] = {}
    for row in rows:
        if int(row["is_origin"]) != 0:
            continue
        value = float(row["modular_aligned"])
        if np.isfinite(value):
            out[int(row["point_id"])] = value
    return out


def _load_method_surfaces(method_manifest_path: str, *, grid_size: int) -> dict[str, Any]:
    manifest = _load_json(method_manifest_path)
    continuum_rows = _load_dat_rows(str(manifest["pointwise_continuum"]))
    modular_rows = _load_dat_rows(str(manifest["modular_aligned"]))
    alignment = _load_json(str(manifest["modular_alignment"]))

    point_rows = [row for row in continuum_rows if int(row["is_origin"]) == 0]
    tau_real = float(manifest["target_tau"]["real"])
    tau_imag = float(manifest["target_tau"]["imag"])

    continuum_lookup = _lookup_continuum(continuum_rows)
    modular_lookup = _lookup_modular(modular_rows)

    continuum_xyz = _interpolate_periodic_manifold(
        point_rows,
        continuum_lookup,
        tau_real=tau_real,
        tau_imag=tau_imag,
        grid_size=grid_size,
    )
    modular_xyz = _interpolate_periodic_manifold(
        point_rows,
        modular_lookup,
        tau_real=tau_real,
        tau_imag=tau_imag,
        grid_size=grid_size,
    )

    sample_x = np.asarray([float(row["nu_real"]) for row in point_rows], dtype=float)
    sample_y = np.asarray([float(row["nu_imag"]) for row in point_rows], dtype=float)

    return {
        "benchmark_id": str(manifest["benchmark_id"]),
        "method": str(manifest["method"]),
        "tau_real": tau_real,
        "tau_imag": tau_imag,
        "continuum_xyz": continuum_xyz,
        "modular_xyz": modular_xyz,
        "sample_x": sample_x,
        "sample_y": sample_y,
        "alpha": float(alignment["alpha"]),
        "chi2_per_dof": float(alignment["chi2_per_dof"]),
    }


def _plot_benchmark(benchmark_manifest_path: str, *, output_path: str | None, grid_size: int) -> str:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    method_order = [name for name in ("twisted", "untwisted") if name in benchmark_manifest["methods"]]
    method_data = [
        _load_method_surfaces(str(benchmark_manifest["methods"][method_name]), grid_size=grid_size)
        for method_name in method_order
    ]
    if len(method_data) == 0:
        raise ValueError("benchmark manifest has no methods to plot")

    z_min = np.inf
    z_max = -np.inf
    for data in method_data:
        for _, _, Z in (data["continuum_xyz"], data["modular_xyz"]):
            finite = np.isfinite(Z)
            if np.any(finite):
                z_min = min(z_min, float(np.nanmin(Z[finite])))
                z_max = max(z_max, float(np.nanmax(Z[finite])))

    if not np.isfinite(z_min) or not np.isfinite(z_max):
        raise ValueError("no finite manifold values available for plotting")

    fig = plt.figure(figsize=(14, 10))
    benchmark_id = str(benchmark_manifest["benchmark_id"])
    fig.suptitle(f"{benchmark_id}: continuum manifolds overlaid with the aligned modular manifold", fontsize=15, y=0.98)

    colors = {
        "twisted": "tab:blue",
        "untwisted": "tab:orange",
    }

    for row_index, data in enumerate(method_data):
        method = str(data["method"])
        color = colors.get(method, "tab:green")
        Xc, Yc, Zc = data["continuum_xyz"]
        Xm, Ym, Zm = data["modular_xyz"]

        ax3d = fig.add_subplot(len(method_data), 2, 2 * row_index + 1, projection="3d")
        ax2d = fig.add_subplot(len(method_data), 2, 2 * row_index + 2)

        ax3d.plot_surface(
            Xc,
            Yc,
            Zc,
            rstride=1,
            cstride=1,
            linewidth=0.0,
            antialiased=True,
            shade=False,
            color=color,
            alpha=0.55,
        )
        ax3d.plot_wireframe(
            Xm,
            Ym,
            Zm,
            rstride=max(grid_size // 18, 4),
            cstride=max(grid_size // 18, 4),
            color="black",
            linewidth=0.6,
            alpha=0.9,
        )
        ax3d.set_title(
            f"{method}: continuum surface + modular wireframe\nchi2/dof={data['chi2_per_dof']:.3f}",
            fontsize=11,
        )
        ax3d.set_xlabel("nu_real")
        ax3d.set_ylabel("nu_imag")
        ax3d.set_zlabel("connected correlator")
        ax3d.view_init(elev=28, azim=-60)
        ax3d.set_zlim(z_min, z_max)

        ax2d.contour(Xc, Yc, Zc, levels=10, colors=[color], linewidths=1.4)
        ax2d.contour(Xm, Ym, Zm, levels=10, colors=["black"], linewidths=1.0, linestyles="--")
        outline_x, outline_y = _unit_cell_outline(data["tau_real"], data["tau_imag"])
        ax2d.plot(outline_x, outline_y, color="black", linewidth=1.2)
        ax2d.scatter(data["sample_x"], data["sample_y"], s=8, color="black", alpha=0.3)
        ax2d.set_title(
            f"{method}: contour comparison (alpha={data['alpha']:.4f})",
            fontsize=11,
        )
        ax2d.set_xlabel("nu_real")
        ax2d.set_ylabel("nu_imag")
        ax2d.set_aspect("equal", adjustable="box")
        ax2d.grid(True, alpha=0.25)
        ax2d.legend(
            handles=[
                Line2D([0], [0], color=color, linewidth=2.0, label="continuum manifold"),
                Line2D([0], [0], color="black", linewidth=1.5, linestyle="--", label="aligned modular manifold"),
            ],
            loc="best",
        )

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(benchmark_manifest_path),
            f"{benchmark_id}_continuum_vs_modular_manifolds.png",
        )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"wrote {output_path}")
    print(f"methods={[data['method'] for data in method_data]}")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay the continuum-extrapolated lattice manifolds with the aligned modular manifold for a benchmark."
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--grid-size", type=int, default=140, help="Interpolation grid size in each cell direction")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_benchmark(
        benchmark_manifest_path=os.path.abspath(args.benchmark_manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
        grid_size=max(int(args.grid_size), 20),
    )


if __name__ == "__main__":
    main()