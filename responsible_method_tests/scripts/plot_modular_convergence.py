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


def _compute_scale_metrics(method_manifest_path: str) -> dict[str, Any]:
    manifest = _load_json(method_manifest_path)
    raw_rows = _load_dat_rows(str(manifest["pointwise_raw"]))
    modular_rows = _load_dat_rows(str(manifest["modular_raw"]))
    point_rows = _load_dat_rows(str(manifest["smallest_lattice_points"]))

    modular_lookup = {
        int(row["point_id"]): float(row["modular_raw"])
        for row in modular_rows
        if int(row["is_origin"]) == 0 and np.isfinite(float(row["modular_raw"]))
    }
    point_lookup = {int(row["point_id"]): row for row in point_rows}

    by_scale: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        point_id = int(row["point_id"])
        if point_id not in modular_lookup:
            continue
        by_scale[int(row["scale"])] += [row]

    scales = sorted(by_scale)
    metrics: list[dict[str, Any]] = []
    for scale in scales:
        rows = by_scale[scale]
        y = np.asarray([float(row["conn"]) for row in rows], dtype=float)
        sigma = np.asarray([float(row["conn_err"]) for row in rows], dtype=float)
        g = np.asarray([float(modular_lookup[int(row["point_id"])]) for row in rows], dtype=float)
        w = np.where(sigma > 0.0, 1.0 / np.square(sigma), 1.0)
        denom = float(np.sum(w * g * g))
        alpha = float(np.sum(w * y * g) / denom) if denom > 0.0 else float("nan")
        residual = y - alpha * g
        chi2 = float(np.sum(w * residual * residual))
        dof = max(len(y) - 1, 1)
        metrics.append(
            {
                "scale": int(scale),
                "family_size": int(rows[0]["family_size"]),
                "alpha": alpha,
                "rms": float(np.sqrt(np.mean(residual * residual))),
                "chi2_dof": float(chi2 / dof),
                "n_points": int(len(y)),
                "values": y,
                "modular_scaled": alpha * g,
                "point_ids": [int(row["point_id"]) for row in rows],
            }
        )

    largest = metrics[-1]
    coords = np.asarray(
        [
            (float(point_lookup[point_id]["nu_real"]), float(point_lookup[point_id]["nu_imag"]))
            for point_id in largest["point_ids"]
        ],
        dtype=float,
    )

    return {
        "benchmark_id": str(manifest["benchmark_id"]),
        "method": str(manifest["method"]),
        "metrics": metrics,
        "largest_coords": coords,
        "largest_residual": np.asarray(largest["values"] - largest["modular_scaled"], dtype=float),
        "largest_values": np.asarray(largest["values"], dtype=float),
        "largest_modular_scaled": np.asarray(largest["modular_scaled"], dtype=float),
    }


def _plot_benchmark(benchmark_manifest_path: str, *, output_path: str | None) -> str:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    twisted = _compute_scale_metrics(str(benchmark_manifest["methods"]["twisted"]))
    untwisted = _compute_scale_metrics(str(benchmark_manifest["methods"]["untwisted"]))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax_rms, ax_chi2, ax_scatter, ax_resid = axes.ravel()

    benchmark_id = str(benchmark_manifest["benchmark_id"])
    fig.suptitle(
        f"{benchmark_id}: finite-scale approach to the modular two-point function",
        fontsize=14,
        y=0.98,
    )

    style = {
        "twisted": {"color": "tab:blue", "marker": "o"},
        "untwisted": {"color": "tab:orange", "marker": "s"},
    }

    for data in (twisted, untwisted):
        method = str(data["method"])
        color = style[method]["color"]
        marker = style[method]["marker"]
        inv_scale = np.asarray([1.0 / float(item["scale"]) for item in data["metrics"]], dtype=float)
        rms = np.asarray([float(item["rms"]) for item in data["metrics"]], dtype=float)
        chi2 = np.asarray([float(item["chi2_dof"]) for item in data["metrics"]], dtype=float)
        sizes = [int(item["family_size"]) for item in data["metrics"]]

        label = f"{method}"
        ax_rms.plot(inv_scale, rms, marker=marker, color=color, label=label)
        ax_chi2.plot(inv_scale, chi2, marker=marker, color=color, label=label)

        largest_values = data["largest_values"]
        largest_modular = data["largest_modular_scaled"]
        ax_scatter.scatter(
            largest_modular,
            largest_values,
            s=10,
            alpha=0.45,
            color=color,
            label=f"{method} largest scale",
        )

        coords = data["largest_coords"]
        residual = data["largest_residual"]
        scatter = ax_resid.scatter(
            coords[:, 0],
            coords[:, 1],
            c=residual,
            cmap="coolwarm",
            s=12,
            alpha=0.75 if method == "twisted" else 0.45,
            vmin=-np.max(np.abs(residual)),
            vmax=np.max(np.abs(residual)),
            label=f"{method} largest-scale residual",
        )

        for x, y, size in zip(inv_scale, rms, sizes):
            ax_rms.annotate(str(size), (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7, color=color)

    lim_min = min(ax_scatter.get_xlim()[0], ax_scatter.get_ylim()[0])
    lim_max = max(ax_scatter.get_xlim()[1], ax_scatter.get_ylim()[1])
    ax_scatter.plot([lim_min, lim_max], [lim_min, lim_max], color="black", linewidth=1.0, linestyle="--")
    ax_scatter.set_xlim(lim_min, lim_max)
    ax_scatter.set_ylim(lim_min, lim_max)

    ax_rms.set_title("Amplitude-fitted RMS mismatch vs 1/scale")
    ax_rms.set_xlabel("1 / scale factor")
    ax_rms.set_ylabel("RMS(conn - alpha g_mod)")
    ax_rms.grid(True, alpha=0.3)
    ax_rms.legend()

    ax_chi2.set_title("Weighted chi2/dof vs 1/scale")
    ax_chi2.set_xlabel("1 / scale factor")
    ax_chi2.set_ylabel("chi2 / dof")
    ax_chi2.grid(True, alpha=0.3)
    ax_chi2.legend()

    ax_scatter.set_title("Largest scale: lattice vs fitted modular target")
    ax_scatter.set_xlabel("alpha_scale * modular_raw")
    ax_scatter.set_ylabel("lattice conn")
    ax_scatter.grid(True, alpha=0.3)
    ax_scatter.legend()

    ax_resid.set_title("Largest-scale residual cloud on target torus coords")
    ax_resid.set_xlabel("nu_real")
    ax_resid.set_ylabel("nu_imag")
    ax_resid.grid(True, alpha=0.3)
    ax_resid.legend(loc="upper right")

    cbar = fig.colorbar(scatter, ax=ax_resid, fraction=0.046, pad=0.04)
    cbar.set_label("conn - alpha_scale * modular_raw")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    if output_path is None:
        output_path = os.path.join(os.path.dirname(benchmark_manifest_path), f"{benchmark_id}_modular_convergence.png")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"wrote {output_path}")
    for data in (twisted, untwisted):
        method = str(data["method"])
        scales = [int(item["scale"]) for item in data["metrics"]]
        rms = [float(item["rms"]) for item in data["metrics"]]
        chi2 = [float(item["chi2_dof"]) for item in data["metrics"]]
        print(f"{method}: scales={scales} rms={rms} chi2_dof={chi2}")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot how twisted and untwisted finite-size correlators approach the modular target."
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_benchmark(
        benchmark_manifest_path=os.path.abspath(args.benchmark_manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
    )


if __name__ == "__main__":
    main()