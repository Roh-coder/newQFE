#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


_INT_RE = re.compile(r"^[+-]?\d+$")


def _parse_value(token: str) -> Any:
    if _INT_RE.match(token):
        return int(token)
    try:
        return float(token)
    except ValueError:
        return token


def _load_json(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_dat_rows(path: str) -> list[dict[str, Any]]:
    columns: list[str] | None = None
    rows: list[dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            if columns is None:
                raise ValueError(f"missing '# columns:' header in {path}")
            parts = line.split()
            if len(parts) != len(columns):
                raise ValueError(
                    f"row in {path} has {len(parts)} fields, expected {len(columns)}: {line}"
                )
            rows.append({key: _parse_value(value) for key, value in zip(columns, parts)})
    return rows


def _group_raw_rows_by_scale(raw_rows: list[dict[str, Any]]) -> dict[int, dict[int, dict[str, Any]]]:
    grouped: dict[int, dict[int, dict[str, Any]]] = defaultdict(dict)
    for row in raw_rows:
        grouped[int(row["scale"])][int(row["point_id"])] = row
    return dict(grouped)


def _tile_periodic_points(
    point_rows: list[dict[str, Any]],
    value_lookup: dict[int, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a_values: list[float] = []
    b_values: list[float] = []
    z_values: list[float] = []
    for row in point_rows:
        point_id = int(row["point_id"])
        if point_id not in value_lookup:
            continue
        z = float(value_lookup[point_id])
        if not np.isfinite(z):
            continue
        a0 = float(row["a_wrap"])
        b0 = float(row["b_wrap"])
        for da in (-1.0, 0.0, 1.0):
            for db in (-1.0, 0.0, 1.0):
                a_values.append(a0 + da)
                b_values.append(b0 + db)
                z_values.append(z)
    return np.asarray(a_values), np.asarray(b_values), np.asarray(z_values)


def _interpolate_periodic_manifold(
    point_rows: list[dict[str, Any]],
    value_lookup: dict[int, float],
    *,
    tau_real: float,
    tau_imag: float,
    grid_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a_tiled, b_tiled, z_tiled = _tile_periodic_points(point_rows, value_lookup)
    triangulation = mtri.Triangulation(a_tiled, b_tiled)
    interpolator = mtri.LinearTriInterpolator(triangulation, z_tiled)

    a_grid = np.linspace(0.0, 1.0, grid_size)
    b_grid = np.linspace(0.0, 1.0, grid_size)
    A, B = np.meshgrid(a_grid, b_grid)
    Z = np.asarray(interpolator(A, B), dtype=float)
    X = A + tau_real * B
    Y = tau_imag * B
    return X, Y, Z


def _unit_cell_outline(tau_real: float, tau_imag: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray([0.0, 1.0, 1.0 + tau_real, tau_real, 0.0], dtype=float)
    y = np.asarray([0.0, 0.0, tau_imag, tau_imag, 0.0], dtype=float)
    return x, y


def _build_scale_value_lookup(
    raw_rows_by_scale: dict[int, dict[int, dict[str, Any]]],
) -> dict[int, dict[int, float]]:
    out: dict[int, dict[int, float]] = {}
    for scale, rows in raw_rows_by_scale.items():
        out[scale] = {point_id: float(row["conn"]) for point_id, row in rows.items()}
    return out


def _plot_manifest(
    manifest_path: str,
    *,
    output_path: str | None,
    grid_size: int,
    include_origin: bool,
) -> str:
    manifest = _load_json(manifest_path)
    point_rows = _load_dat_rows(str(manifest["smallest_lattice_points"]))
    raw_rows = _load_dat_rows(str(manifest["pointwise_raw"]))

    if not include_origin:
        point_rows = [row for row in point_rows if int(row["is_origin"]) == 0]

    tau_real = float(manifest["target_tau"]["real"])
    tau_imag = float(manifest["target_tau"]["imag"])
    raw_rows_by_scale = _group_raw_rows_by_scale(raw_rows)
    value_lookup_by_scale = _build_scale_value_lookup(raw_rows_by_scale)
    scales = sorted(value_lookup_by_scale)
    if len(scales) == 0:
        raise ValueError("no scales found in pointwise_raw table")

    interpolated: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    z_min = np.inf
    z_max = -np.inf
    for scale in scales:
        X, Y, Z = _interpolate_periodic_manifold(
            point_rows,
            value_lookup_by_scale[scale],
            tau_real=tau_real,
            tau_imag=tau_imag,
            grid_size=grid_size,
        )
        interpolated[scale] = (X, Y, Z)
        finite = np.isfinite(Z)
        if np.any(finite):
            z_min = min(z_min, float(np.nanmin(Z[finite])))
            z_max = max(z_max, float(np.nanmax(Z[finite])))

    if not np.isfinite(z_min) or not np.isfinite(z_max):
        raise ValueError("interpolated manifold contains no finite values")

    benchmark_id = str(manifest["benchmark_id"])
    method = str(manifest["method"])
    family_sizes = {int(payload["scale"]): int(payload["family_size"]) for payload in manifest["payloads"]}

    fig = plt.figure(figsize=(15, 7))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)
    fig.suptitle(
        f"{benchmark_id} {method}: transparent interpolated FSS manifolds (no extrapolation)",
        fontsize=14,
        y=0.98,
    )

    cmap = plt.get_cmap("viridis")
    legend_handles: list[mpatches.Patch] = []
    n_scales = len(scales)
    outline_x, outline_y = _unit_cell_outline(tau_real, tau_imag)

    for idx, scale in enumerate(scales):
        color = cmap(idx / max(n_scales - 1, 1))
        alpha = 0.16 + 0.10 * idx / max(n_scales - 1, 1)
        X, Y, Z = interpolated[scale]
        ax3d.plot_surface(
            X,
            Y,
            Z,
            rstride=1,
            cstride=1,
            linewidth=0.0,
            antialiased=True,
            shade=False,
            color=color,
            alpha=alpha,
        )
        ax2d.contour(
            X,
            Y,
            Z,
            levels=10,
            colors=[color],
            linewidths=1.0,
            alpha=min(alpha + 0.25, 0.85),
        )
        legend_handles.append(
            mpatches.Patch(color=color, alpha=max(alpha, 0.35), label=f"scale={scale} (L={family_sizes[scale]})")
        )

    sample_x = np.asarray([float(row["nu_real"]) for row in point_rows], dtype=float)
    sample_y = np.asarray([float(row["nu_imag"]) for row in point_rows], dtype=float)
    ax2d.scatter(sample_x, sample_y, s=8, color="black", alpha=0.35, label="sample points")
    ax2d.plot(outline_x, outline_y, color="black", linewidth=1.2)
    ax3d.plot(outline_x, outline_y, zs=[z_min] * len(outline_x), color="black", linewidth=1.0)

    ax3d.set_xlabel("nu_real")
    ax3d.set_ylabel("nu_imag")
    ax3d.set_zlabel("connected correlator")
    ax3d.set_title("Transparent surface overlay", fontsize=11)
    ax3d.view_init(elev=28, azim=-60)
    ax3d.set_zlim(z_min, z_max)

    ax2d.set_xlabel("nu_real")
    ax2d.set_ylabel("nu_imag")
    ax2d.set_title("Contour overlay on the target unit cell", fontsize=11)
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.grid(True, alpha=0.25)

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=min(4, len(legend_handles)),
        frameon=False,
    )
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.90])

    if output_path is None:
        output_path = os.path.join(os.path.dirname(manifest_path), "fss_interpolated_manifolds.png")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"wrote {output_path}")
    print(f"scales={scales}")
    print(f"include_origin={include_origin}")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot transparent interpolated FSS manifolds from a method manifest without continuum extrapolation."
    )
    parser.add_argument("--manifest", required=True, help="Path to manifest_geometry_*_{method}.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--grid-size", type=int, default=140, help="Interpolation grid size in each cell direction")
    parser.add_argument(
        "--include-origin",
        action="store_true",
        help="Include the origin point in the interpolated surface",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_manifest(
        manifest_path=os.path.abspath(args.manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
        grid_size=max(int(args.grid_size), 20),
        include_origin=bool(args.include_origin),
    )


if __name__ == "__main__":
    main()