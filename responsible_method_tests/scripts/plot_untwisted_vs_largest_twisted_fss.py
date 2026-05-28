#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import defaultdict
from typing import Any

import matplotlib
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay, QhullError

matplotlib.use("Agg")
import matplotlib.pyplot as plt


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_REPO_ROOT = os.path.normpath(os.path.join(_PROJECT_ROOT, ".."))
_KFC_ROOT = os.path.join(_REPO_ROOT, "K_from_continuum")
if _KFC_ROOT not in sys.path:
    sys.path.insert(0, _KFC_ROOT)

from workflow_common import evaluate_observable_fit  # noqa: E402


_INT_RE = re.compile(r"^[+-]?\d+$")
_COLUMN_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


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
                if columns is None:
                    maybe_columns = line[1:].strip().split()
                    if maybe_columns and all(_COLUMN_NAME_RE.match(token) for token in maybe_columns):
                        columns = maybe_columns
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


def _boundary_paths(lx: int, ly: int, tx: int, ty: int) -> list[tuple[int, int]]:
    return [
        (lx, ty),
        (tx, -ly),
        (-lx - tx, ly - ty),
    ]


def _to_ab(
    m: int,
    n: int,
    lx: int,
    ly: int,
    tx: int,
    ty: int,
    *,
    embedding_cycles: tuple[int, int],
) -> tuple[float, float]:
    paths = _boundary_paths(lx, ly, tx, ty)
    (dm_a, dn_a), (dm_b, dn_b) = [paths[idx] for idx in embedding_cycles]
    det = float(dm_a * dn_b - dn_a * dm_b)
    if det == 0.0:
        raise ValueError(
            f"embedding cycle basis {embedding_cycles} is degenerate for lattice {(lx, ly, tx, ty)}"
        )
    a_value = (float(dn_b) * float(m) - float(dm_b) * float(n)) / det
    b_value = (-float(dn_a) * float(m) + float(dm_a) * float(n)) / det
    return a_value, b_value


def _wrap_unit(value: float) -> float:
    return value % 1.0


def _unit_cell_outline(tau_real: float, tau_imag: float) -> tuple[np.ndarray, np.ndarray]:
    x_value = np.asarray([0.0, 1.0, 1.0 + tau_real, tau_real, 0.0], dtype=float)
    y_value = np.asarray([0.0, 0.0, tau_imag, tau_imag, 0.0], dtype=float)
    return x_value, y_value


def _select_representative_points(
    continuum_rows: list[dict[str, Any]],
    *,
    n_panels: int,
) -> list[int]:
    candidates = [
        row
        for row in continuum_rows
        if int(row["is_origin"]) == 0 and int(row["m"]) == 0 and 0.0 < float(row["b_wrap"]) <= 0.5 + 1.0e-12
    ]
    if len(candidates) < n_panels:
        candidates = [
            row
            for row in continuum_rows
            if int(row["is_origin"]) == 0 and 0.0 < float(row["b_wrap"]) <= 0.5 + 1.0e-12
        ]
    if len(candidates) == 0:
        raise ValueError("no non-origin points available for plotting")

    candidates.sort(key=lambda row: (float(row["b_wrap"]), float(row["a_wrap"]), int(row["point_id"])))
    idx_values = np.linspace(0, len(candidates) - 1, num=min(n_panels, len(candidates)), dtype=int)

    selected: list[int] = []
    seen: set[int] = set()
    for idx in idx_values:
        point_id = int(candidates[int(idx)]["point_id"])
        if point_id not in seen:
            selected.append(point_id)
            seen.add(point_id)

    for row in candidates:
        if len(selected) >= min(n_panels, len(candidates)):
            break
        point_id = int(row["point_id"])
        if point_id not in seen:
            selected.append(point_id)
            seen.add(point_id)
    return selected


def _group_raw_rows(raw_rows: list[dict[str, Any]]) -> dict[int, list[dict[str, Any]]]:
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        grouped[int(row["point_id"])].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: float(row["scale"]))
    return grouped


def _build_largest_twisted_interpolator(
    twisted_manifest: dict[str, Any],
) -> tuple[Any, Any, np.ndarray, np.ndarray, np.ndarray, int]:
    payloads = list(twisted_manifest["payloads"])
    largest_payload = max(payloads, key=lambda payload: int(payload["scale"]))
    lx, ly, tx, ty = [int(value) for value in largest_payload["lattice"]]
    embedding_cycles = tuple(int(value) for value in twisted_manifest.get("embedding_cycles", [0, 1]))
    largest_rows = _load_dat_rows(str(largest_payload["data_path"]))

    unique_points: dict[tuple[float, float], list[float]] = {}

    for row in largest_rows:
        conn_value = float(row["corr_conn"])
        if not np.isfinite(conn_value):
            continue
        a_raw, b_raw = _to_ab(
            int(row["m"]),
            int(row["n"]),
            lx,
            ly,
            tx,
            ty,
            embedding_cycles=embedding_cycles,
        )
        a_wrap = _wrap_unit(a_raw)
        b_wrap = _wrap_unit(b_raw)
        key = (round(a_wrap, 12), round(b_wrap, 12))
        unique_points.setdefault(key, []).append(conn_value)

    a_tiled: list[float] = []
    b_tiled: list[float] = []
    z_tiled: list[float] = []
    a_core: list[float] = []
    b_core: list[float] = []
    z_core: list[float] = []
    for (a_wrap, b_wrap), values in sorted(unique_points.items()):
        conn_value = float(np.mean(values))
        a_core.append(a_wrap)
        b_core.append(b_wrap)
        z_core.append(conn_value)
        for delta_a in (-1.0, 0.0, 1.0):
            for delta_b in (-1.0, 0.0, 1.0):
                a_tiled.append(a_wrap + delta_a)
                b_tiled.append(b_wrap + delta_b)
                z_tiled.append(conn_value)

    points_tiled = np.column_stack([np.asarray(a_tiled, dtype=float), np.asarray(b_tiled, dtype=float)])
    z_tiled_array = np.asarray(z_tiled, dtype=float)
    try:
        triangulation = Delaunay(points_tiled, qhull_options="QJ Qc Qbb Q12")
        interpolator = LinearNDInterpolator(triangulation, z_tiled_array, fill_value=np.nan)
    except QhullError:
        interpolator = LinearNDInterpolator(points_tiled, z_tiled_array, fill_value=np.nan)
    nearest_interpolator = NearestNDInterpolator(points_tiled, z_tiled_array)
    return (
        interpolator,
        nearest_interpolator,
        np.asarray(a_core, dtype=float),
        np.asarray(b_core, dtype=float),
        np.asarray(z_core, dtype=float),
        int(largest_payload["scale"]),
    )


def _evaluate_interpolator(
    interpolator: Any,
    nearest_interpolator: Any,
    a_value: float,
    b_value: float,
    *,
    a_core: np.ndarray,
    b_core: np.ndarray,
    z_core: np.ndarray,
) -> float:
    interpolated = np.asarray(interpolator(np.asarray([[a_value, b_value]])), dtype=float).reshape(-1)
    if interpolated.size > 0 and np.isfinite(interpolated[0]):
        return float(interpolated[0])

    nearest_value = np.asarray(nearest_interpolator(np.asarray([[a_value, b_value]])), dtype=float).reshape(-1)
    if nearest_value.size > 0 and np.isfinite(nearest_value[0]):
        return float(nearest_value[0])

    best_distance = float("inf")
    best_value = float("nan")
    for delta_a in (-1.0, 0.0, 1.0):
        for delta_b in (-1.0, 0.0, 1.0):
            dist2 = np.square(a_core - (a_value + delta_a)) + np.square(b_core - (b_value + delta_b))
            idx = int(np.argmin(dist2))
            if float(dist2[idx]) < best_distance:
                best_distance = float(dist2[idx])
                best_value = float(z_core[idx])
    return best_value


def _panel_title(row: dict[str, Any], twisted_value: float) -> str:
    return (
        f"point {int(row['point_id'])}: (m,n)=({int(row['m'])},{int(row['n'])})\n"
        f"(a,b)=({float(row['a_wrap']):.3f}, {float(row['b_wrap']):.3f}), twisted Lmax={twisted_value:.4f}"
    )


def _write_table(rows: list[dict[str, Any]], table_path: str) -> None:
    with open(table_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "point_id",
                "m",
                "n",
                "a_wrap",
                "b_wrap",
                "nu_real",
                "nu_imag",
                "A_continuum",
                "sigma_A",
                "twisted_largest_interp",
                "delta_A_minus_twisted",
                "z_A_minus_twisted",
                "is_origin",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    int(row["point_id"]),
                    int(row["m"]),
                    int(row["n"]),
                    f"{float(row['a_wrap']):.8f}",
                    f"{float(row['b_wrap']):.8f}",
                    f"{float(row['nu_real']):.8f}",
                    f"{float(row['nu_imag']):.8f}",
                    f"{float(row['A']):.8f}",
                    f"{float(row['sigma_A']):.8f}",
                    f"{float(row['twisted_largest_interp']):.8f}",
                    f"{float(row['delta_A_minus_twisted']):.8f}",
                    f"{float(row['z_A_minus_twisted']):.8f}",
                    int(row["is_origin"]),
                ]
            )


def _plot_benchmark(
    benchmark_manifest_path: str,
    *,
    output_path: str | None,
    table_output_path: str | None,
    point_ids: list[int] | None,
    n_panels: int,
) -> tuple[str, str, list[int], dict[str, float]]:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    untwisted_manifest = _load_json(str(benchmark_manifest["methods"]["untwisted"]))
    twisted_manifest = _load_json(str(benchmark_manifest["methods"]["twisted"]))

    untwisted_raw_rows = _load_dat_rows(str(untwisted_manifest["pointwise_raw"]))
    untwisted_continuum_rows = _load_dat_rows(str(untwisted_manifest["pointwise_continuum"]))
    untwisted_point_rows = _load_dat_rows(str(untwisted_manifest["smallest_lattice_points"]))

    raw_by_id = _group_raw_rows(untwisted_raw_rows)
    point_by_id = {int(row["point_id"]): row for row in untwisted_point_rows}

    interpolator, nearest_interpolator, a_core, b_core, z_core, largest_twisted_scale = _build_largest_twisted_interpolator(
        twisted_manifest
    )

    comparison_rows: list[dict[str, Any]] = []
    for row in untwisted_continuum_rows:
        point_id = int(row["point_id"])
        twisted_value = _evaluate_interpolator(
            interpolator,
            nearest_interpolator,
            float(row["a_wrap"]),
            float(row["b_wrap"]),
            a_core=a_core,
            b_core=b_core,
            z_core=z_core,
        )
        delta = float(row["A"]) - twisted_value
        sigma_A = float(row["sigma_A"])
        z_value = delta / sigma_A if np.isfinite(sigma_A) and sigma_A > 0.0 else float("nan")
        comparison_rows.append(
            {
                **row,
                "twisted_largest_interp": twisted_value,
                "delta_A_minus_twisted": delta,
                "z_A_minus_twisted": z_value,
            }
        )

    selected_point_ids = point_ids or _select_representative_points(untwisted_continuum_rows, n_panels=n_panels)
    selected_point_ids = [
        point_id
        for point_id in selected_point_ids
        if point_id in raw_by_id and any(int(row["point_id"]) == point_id for row in comparison_rows)
    ]
    if len(selected_point_ids) == 0:
        raise ValueError("no selected points were present in the untwisted raw and continuum tables")

    selected_lookup = {int(row["point_id"]): row for row in comparison_rows}
    tau_real = float(untwisted_manifest["target_tau"]["real"])
    tau_imag = float(untwisted_manifest["target_tau"]["imag"])

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes_flat = list(axes.ravel())
    ax_scatter = axes_flat[0]
    ax_residual = axes_flat[1]
    point_axes = axes_flat[2:]

    benchmark_id = str(benchmark_manifest["benchmark_id"])
    fig.suptitle(
        f"{benchmark_id}: untwisted continuum vs interpolated largest twisted lattice",
        fontsize=15,
        y=0.98,
    )

    non_origin_rows = [row for row in comparison_rows if int(row["is_origin"]) == 0]
    x_values = np.asarray([float(row["twisted_largest_interp"]) for row in non_origin_rows], dtype=float)
    y_values = np.asarray([float(row["A"]) for row in non_origin_rows], dtype=float)
    sigma_values = np.asarray([float(row["sigma_A"]) for row in non_origin_rows], dtype=float)
    ax_scatter.errorbar(x_values, y_values, yerr=sigma_values, fmt="o", ms=3.0, alpha=0.55, capsize=2)
    lim_min = min(float(np.min(x_values)), float(np.min(y_values)))
    lim_max = max(float(np.max(x_values)), float(np.max(y_values)))
    padding = 0.04 * (lim_max - lim_min)
    ax_scatter.plot([lim_min - padding, lim_max + padding], [lim_min - padding, lim_max + padding], linestyle="--", color="black")
    ax_scatter.set_xlim(lim_min - padding, lim_max + padding)
    ax_scatter.set_ylim(lim_min - padding, lim_max + padding)
    ax_scatter.set_xlabel(f"largest twisted interpolation (scale={largest_twisted_scale})")
    ax_scatter.set_ylabel("untwisted continuum A")
    ax_scatter.set_title("All non-origin continuum points")
    ax_scatter.grid(True, alpha=0.3)

    nu_real = np.asarray([float(row["nu_real"]) for row in non_origin_rows], dtype=float)
    nu_imag = np.asarray([float(row["nu_imag"]) for row in non_origin_rows], dtype=float)
    z_values_map = np.asarray([float(row["z_A_minus_twisted"]) for row in non_origin_rows], dtype=float)
    vmax = float(np.nanmax(np.abs(z_values_map))) if len(z_values_map) else 1.0
    scatter = ax_residual.scatter(
        nu_real,
        nu_imag,
        c=z_values_map,
        cmap="coolwarm",
        s=18,
        alpha=0.8,
        vmin=-vmax,
        vmax=vmax,
    )
    outline_x, outline_y = _unit_cell_outline(tau_real, tau_imag)
    ax_residual.plot(outline_x, outline_y, color="black", linewidth=1.0)
    ax_residual.set_aspect("equal", adjustable="box")
    ax_residual.set_xlabel("nu_real")
    ax_residual.set_ylabel("nu_imag")
    ax_residual.set_title("Pointwise z = (A - twisted_interp) / sigma_A")
    ax_residual.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax_residual, fraction=0.046, pad=0.04)
    cbar.set_label("continuum-twisted z")

    for axis, point_id in zip(point_axes, selected_point_ids):
        point_rows = raw_by_id[point_id]
        continuum_row = selected_lookup[point_id]
        inv_scale = np.asarray([1.0 / float(row["scale"]) for row in point_rows], dtype=float)
        y_point = np.asarray([float(row["conn"]) for row in point_rows], dtype=float)
        yerr_point = np.asarray([float(row["conn_err"]) for row in point_rows], dtype=float)

        axis.errorbar(inv_scale, y_point, yerr=yerr_point, fmt="o", capsize=3, label="untwisted MC")

        x_fit = np.linspace(0.0, float(np.max(inv_scale)) * 1.05, 300)
        y_fit = evaluate_observable_fit(
            x_fit,
            float(continuum_row["A"]),
            float(continuum_row["B"]),
            float(continuum_row["C"]),
            str(continuum_row["fit_mode"]),
        )
        axis.plot(x_fit, y_fit, color="tab:orange", label=f"fit: {continuum_row['fit_mode']}")
        axis.errorbar(
            [0.0],
            [float(continuum_row["A"])],
            yerr=[float(continuum_row["sigma_A"])],
            fmt="s",
            color="tab:red",
            capsize=3,
            label="continuum A",
        )
        axis.axhline(
            float(continuum_row["twisted_largest_interp"]),
            color="tab:green",
            linestyle="--",
            linewidth=1.2,
            label="largest twisted interp",
        )

        delta = float(continuum_row["delta_A_minus_twisted"])
        z_value = float(continuum_row["z_A_minus_twisted"])
        axis.set_title(_panel_title(point_by_id[point_id], float(continuum_row["twisted_largest_interp"])), fontsize=10)
        axis.text(
            0.03,
            0.04,
            f"A - twisted = {delta:+.4f}\nz = {z_value:+.2f}",
            transform=axis.transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )
        axis.set_xlabel("1 / scale factor")
        axis.set_ylabel("connected correlator")
        axis.grid(True, alpha=0.3)
        axis.set_xlim(-0.02, float(np.max(inv_scale)) * 1.08)
        axis.legend(fontsize=7)

    for axis in point_axes[len(selected_point_ids):]:
        axis.axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(benchmark_manifest_path),
            f"{benchmark_id}_untwisted_vs_largest_twisted_fss.png",
        )
    if table_output_path is None:
        table_output_path = os.path.join(
            os.path.dirname(benchmark_manifest_path),
            f"{benchmark_id}_untwisted_vs_largest_twisted_points.tsv",
        )

    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    _write_table(comparison_rows, table_output_path)

    valid_z = np.asarray([float(row["z_A_minus_twisted"]) for row in non_origin_rows], dtype=float)
    summary = {
        "n_points": float(len(non_origin_rows)),
        "rms_delta": float(np.sqrt(np.mean(np.square([float(row["delta_A_minus_twisted"]) for row in non_origin_rows])))),
        "rms_z": float(np.sqrt(np.mean(np.square(valid_z)))),
        "max_abs_z": float(np.max(np.abs(valid_z))),
    }
    print(f"wrote {output_path}")
    print(f"wrote {table_output_path}")
    print(f"largest_twisted_scale={largest_twisted_scale}")
    print(f"selected_point_ids={selected_point_ids}")
    print(
        "summary="
        f"n_points:{int(summary['n_points'])} rms_delta:{summary['rms_delta']:.6f} "
        f"rms_z:{summary['rms_z']:.6f} max_abs_z:{summary['max_abs_z']:.6f}"
    )
    return output_path, table_output_path, selected_point_ids, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare untwisted continuum extrapolations against the interpolated largest twisted lattice "
            "and plot selected untwisted pointwise FSS panels."
        )
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--table-output", default=None, help="Optional TSV output path")
    parser.add_argument(
        "--point-id",
        action="append",
        dest="point_ids",
        type=int,
        help="Optional point_id to plot; may be repeated",
    )
    parser.add_argument("--n-panels", type=int, default=4, help="Number of automatically selected FSS panels")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_benchmark(
        benchmark_manifest_path=os.path.abspath(args.benchmark_manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
        table_output_path=os.path.abspath(args.table_output) if args.table_output else None,
        point_ids=list(args.point_ids or []),
        n_panels=max(int(args.n_panels), 1),
    )


if __name__ == "__main__":
    main()