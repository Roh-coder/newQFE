#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from typing import Any

import matplotlib
import numpy as np

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


def _select_representative_points(
    continuum_rows: list[dict[str, Any]],
    *,
    n_panels: int,
) -> list[int]:
    candidates = [
        row for row in continuum_rows
        if int(row["is_origin"]) == 0 and int(row["m"]) == 0 and 0.0 < float(row["b_wrap"]) <= 0.5 + 1.0e-12
    ]
    if len(candidates) < n_panels:
        candidates = [
            row for row in continuum_rows
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


def _build_modular_lookup(rows: list[dict[str, Any]]) -> dict[int, float]:
    lookup: dict[int, float] = {}
    for row in rows:
        value = row.get("modular_aligned")
        if isinstance(value, (int, float)) and np.isfinite(float(value)):
            lookup[int(row["point_id"])] = float(value)
    return lookup


def _panel_title(row: dict[str, Any]) -> str:
    return (
        f"point {int(row['point_id'])}: (m,n)=({int(row['m'])},{int(row['n'])})\n"
        f"(a,b)=({float(row['a_wrap']):.3f}, {float(row['b_wrap']):.3f})"
    )


def _plot_manifest(
    manifest_path: str,
    *,
    output_path: str | None,
    point_ids: list[int] | None,
    n_panels: int,
) -> str:
    manifest = _load_json(manifest_path)
    raw_rows = _load_dat_rows(str(manifest["pointwise_raw"]))
    continuum_rows = _load_dat_rows(str(manifest["pointwise_continuum"]))
    modular_rows = _load_dat_rows(str(manifest["modular_aligned"]))

    continuum_by_id = {int(row["point_id"]): row for row in continuum_rows}
    raw_by_id = _group_raw_rows(raw_rows)
    modular_by_id = _build_modular_lookup(modular_rows)

    selected_point_ids = point_ids or _select_representative_points(continuum_rows, n_panels=n_panels)
    selected_point_ids = [point_id for point_id in selected_point_ids if point_id in continuum_by_id and point_id in raw_by_id]
    if len(selected_point_ids) == 0:
        raise ValueError("no selected points were present in both raw and continuum tables")

    n_cols = 2 if len(selected_point_ids) > 1 else 1
    n_rows = int(np.ceil(len(selected_point_ids) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.8 * n_cols, 4.8 * n_rows), squeeze=False)
    axes_flat = list(axes.ravel())

    benchmark_id = str(manifest["benchmark_id"])
    method = str(manifest["method"])
    fig.suptitle(f"{benchmark_id} {method}: pointwise FSS of connected correlators", fontsize=14)

    for ax, point_id in zip(axes_flat, selected_point_ids):
        continuum_row = continuum_by_id[point_id]
        point_rows = raw_by_id[point_id]
        inv_scale = np.asarray([1.0 / float(row["scale"]) for row in point_rows], dtype=float)
        y = np.asarray([float(row["conn"]) for row in point_rows], dtype=float)
        yerr = np.asarray([float(row["conn_err"]) for row in point_rows], dtype=float)

        ax.errorbar(inv_scale, y, yerr=yerr, fmt="o", capsize=3, label="MC data")

        x_fit = np.linspace(0.0, float(np.max(inv_scale)) * 1.05, 300)
        y_fit = evaluate_observable_fit(
            x_fit,
            float(continuum_row["A"]),
            float(continuum_row["B"]),
            float(continuum_row["C"]),
            str(continuum_row["fit_mode"]),
        )
        ax.plot(x_fit, y_fit, color="tab:orange", label=f"fit: {continuum_row['fit_mode']}")

        A = float(continuum_row["A"])
        sigma_A = float(continuum_row["sigma_A"])
        if np.isfinite(A):
            ax.errorbar([0.0], [A], yerr=[sigma_A] if np.isfinite(sigma_A) else None, fmt="s", color="tab:red", capsize=3, label="continuum A")

        modular_value = modular_by_id.get(point_id)
        if modular_value is not None and np.isfinite(modular_value):
            ax.axhline(modular_value, color="tab:green", linestyle="--", linewidth=1.2, label="aligned modular")

        ax.set_title(_panel_title(continuum_row), fontsize=11)
        ax.set_xlabel("1 / scale factor")
        ax.set_ylabel("connected correlator")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.02, float(np.max(inv_scale)) * 1.08)
        ax.legend(fontsize=8)

    for ax in axes_flat[len(selected_point_ids):]:
        ax.axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])

    if output_path is None:
        output_path = os.path.join(os.path.dirname(manifest_path), "fss_selected_points.png")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"wrote {output_path}")
    print(f"selected_point_ids={selected_point_ids}")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot pointwise FSS panels from a method manifest.")
    parser.add_argument("--manifest", required=True, help="Path to manifest_geometry_*_{method}.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument(
        "--point-id",
        action="append",
        dest="point_ids",
        type=int,
        help="Optional point_id to plot; may be repeated",
    )
    parser.add_argument("--n-panels", type=int, default=4, help="Number of automatically selected panels")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_manifest(
        manifest_path=os.path.abspath(args.manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
        point_ids=list(args.point_ids or []),
        n_panels=max(int(args.n_panels), 1),
    )


if __name__ == "__main__":
    main()