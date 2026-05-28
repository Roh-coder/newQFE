#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import defaultdict
from typing import Any
import textwrap

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


def _compute_metrics(method_manifest_path: str, alpha_override: float | None) -> dict[str, Any]:
    manifest = _load_json(method_manifest_path)
    raw_rows = _load_dat_rows(str(manifest["pointwise_raw"]))
    modular_rows = _load_dat_rows(str(manifest["modular_raw"]))
    smallest_rows = _load_dat_rows(str(manifest["smallest_lattice_points"]))
    alignment = _load_json(str(manifest["modular_alignment"]))

    alpha = float(alpha_override) if alpha_override is not None else float(alignment["alpha"])
    modular_lookup = {
        int(row["point_id"]): float(row["modular_raw"])
        for row in modular_rows
        if int(row["is_origin"]) == 0 and np.isfinite(float(row["modular_raw"]))
    }
    coords_lookup = {
        int(row["point_id"]): (float(row["nu_real"]), float(row["nu_imag"]))
        for row in smallest_rows
    }

    payload_lookup: dict[int, dict[str, Any]] = {}
    for payload in manifest.get("payloads", []):
        scale = int(payload["scale"])
        metadata_path = payload.get("metadata_path")
        metadata: dict[str, Any] = {}
        if isinstance(metadata_path, str) and os.path.exists(metadata_path):
            metadata = _load_json(metadata_path)
        lattice = payload.get("lattice", [metadata.get("Lx"), metadata.get("Ly"), metadata.get("Tx"), metadata.get("Ty")])
        payload_lookup[scale] = {
            "lattice": tuple(int(value) for value in lattice),
            "n_traj": int(metadata.get("mc", {}).get("n_traj", 0)) if metadata else 0,
            "production_wall_seconds": float(metadata.get("production_wall_seconds", float("nan"))) if metadata else float("nan"),
        }

    by_scale: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        point_id = int(row["point_id"])
        if point_id not in modular_lookup:
            continue
        by_scale[int(row["scale"])].append(row)

    metrics: list[dict[str, Any]] = []
    largest_payload: dict[str, Any] | None = None
    for scale in sorted(by_scale):
        rows = by_scale[scale]
        values = np.asarray([float(row["conn"]) for row in rows], dtype=float)
        modular = np.asarray([float(modular_lookup[int(row["point_id"])]) for row in rows], dtype=float)
        residual = values - alpha * modular
        point_ids = [int(row["point_id"]) for row in rows]
        payload_meta = payload_lookup.get(scale, {})
        lattice = tuple(payload_meta.get("lattice", (np.nan, np.nan, np.nan, np.nan)))
        item = {
            "scale": int(scale),
            "family_size": int(rows[0]["family_size"]),
            "n_points": int(len(rows)),
            "Lx": int(lattice[0]),
            "Ly": int(lattice[1]),
            "Tx": int(lattice[2]),
            "Ty": int(lattice[3]),
            "n_traj": int(payload_meta.get("n_traj", 0)),
            "production_wall_seconds": float(payload_meta.get("production_wall_seconds", float("nan"))),
            "mean_residual": float(np.mean(residual)),
            "mean_abs_residual": float(np.mean(np.abs(residual))),
            "rms_residual": float(np.sqrt(np.mean(residual * residual))),
            "max_abs_residual": float(np.max(np.abs(residual))),
            "coords": np.asarray([coords_lookup[point_id] for point_id in point_ids], dtype=float),
            "point_ids": point_ids,
            "residual": residual,
            "values": values,
            "modular_scaled": alpha * modular,
        }
        metrics.append(item)
        largest_payload = item

    if largest_payload is None:
        raise ValueError(f"no comparable non-origin points found in {method_manifest_path}")

    return {
        "benchmark_id": str(manifest["benchmark_id"]),
        "method": str(manifest["method"]),
        "alpha": alpha,
        "metrics": metrics,
        "largest": largest_payload,
    }


def _write_table(metrics: dict[str, Any], table_path: str) -> None:
    with open(table_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scale",
                "family_size",
                "n_points",
                "Lx",
                "Ly",
                "Tx",
                "Ty",
                "n_traj",
                "production_wall_seconds",
                "mean_residual",
                "mean_abs_residual_per_site",
                "rms_residual",
                "max_abs_residual",
                "alpha_fixed",
            ]
        )
        for item in metrics["metrics"]:
            writer.writerow(
                [
                    int(item["scale"]),
                    int(item["family_size"]),
                    int(item["n_points"]),
                    int(item["Lx"]),
                    int(item["Ly"]),
                    int(item["Tx"]),
                    int(item["Ty"]),
                    int(item["n_traj"]),
                    f"{float(item['production_wall_seconds']):.6f}",
                    f"{float(item['mean_residual']):.10f}",
                    f"{float(item['mean_abs_residual']):.10f}",
                    f"{float(item['rms_residual']):.10f}",
                    f"{float(item['max_abs_residual']):.10f}",
                    f"{float(metrics['alpha']):.10f}",
                ]
            )


def _plot_metrics(metrics: dict[str, Any], output_path: str) -> None:
    benchmark_id = str(metrics["benchmark_id"])
    method = str(metrics["method"])
    alpha = float(metrics["alpha"])

    items = metrics["metrics"]
    inv_scale = np.asarray([1.0 / float(item["scale"]) for item in items], dtype=float)
    mean_abs = np.asarray([float(item["mean_abs_residual"]) for item in items], dtype=float)
    rms = np.asarray([float(item["rms_residual"]) for item in items], dtype=float)
    mean_signed = np.asarray([float(item["mean_residual"]) for item in items], dtype=float)
    sizes = [int(item["family_size"]) for item in items]

    largest = metrics["largest"]
    coords = np.asarray(largest["coords"], dtype=float)
    residual = np.asarray(largest["residual"], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))
    ax_abs, ax_rms, ax_signed, ax_cloud = axes.ravel()
    fig.suptitle(
        f"{benchmark_id} {method}: fixed-alpha modular residual collapse",
        fontsize=14,
        y=0.98,
    )

    ax_abs.plot(inv_scale, mean_abs, marker="o", color="tab:blue")
    ax_abs.set_title("Mean |conn - alpha g_mod| per site")
    ax_abs.set_xlabel("1 / scale factor")
    ax_abs.set_ylabel("mean absolute residual")
    ax_abs.grid(True, alpha=0.3)

    ax_rms.plot(inv_scale, rms, marker="o", color="tab:orange")
    ax_rms.set_title("RMS(conn - alpha g_mod) per site")
    ax_rms.set_xlabel("1 / scale factor")
    ax_rms.set_ylabel("RMS residual")
    ax_rms.grid(True, alpha=0.3)

    ax_signed.plot(inv_scale, mean_signed, marker="o", color="tab:green")
    ax_signed.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax_signed.set_title(f"Mean signed residual per site (alpha={alpha:.6f})")
    ax_signed.set_xlabel("1 / scale factor")
    ax_signed.set_ylabel("mean residual")
    ax_signed.grid(True, alpha=0.3)

    vmax = float(np.max(np.abs(residual))) if residual.size else 1.0
    scatter = ax_cloud.scatter(
        coords[:, 0],
        coords[:, 1],
        c=residual,
        cmap="coolwarm",
        s=14,
        vmin=-vmax,
        vmax=vmax,
        alpha=0.8,
    )
    ax_cloud.set_title(f"Largest scale residual cloud (size={int(largest['family_size'])})")
    ax_cloud.set_xlabel("nu_real")
    ax_cloud.set_ylabel("nu_imag")
    ax_cloud.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax_cloud, fraction=0.046, pad=0.04)
    cbar.set_label("conn - alpha g_mod")

    for axis, values in ((ax_abs, mean_abs), (ax_rms, rms), (ax_signed, mean_signed)):
        for x_value, y_value, size in zip(inv_scale, values, sizes):
            axis.annotate(str(size), (x_value, y_value), textcoords="offset points", xytext=(4, 4), fontsize=8)

    caption_entries = [
        (
            f"s{int(item['scale'])}: family={int(item['family_size'])}, "
            f"(Lx,Ly,Tx,Ty)=({int(item['Lx'])},{int(item['Ly'])},{int(item['Tx'])},{int(item['Ty'])}), "
            f"n_traj={int(item['n_traj'])}"
        )
        for item in items
    ]
    caption_text = (
        f"Fixed-alpha caption: alpha={alpha:.9f}; residuals averaged over {int(largest['n_points'])} non-origin matched points.\n"
        + "\n".join(textwrap.wrap("; ".join(caption_entries), width=120, break_long_words=False, break_on_hyphens=False))
    )
    fig.text(0.01, 0.01, caption_text, ha="left", va="bottom", fontsize=8)

    fig.tight_layout(rect=[0.0, 0.12, 1.0, 0.95])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot fixed-alpha modular residual collapse for a single method manifest across lattice scales."
        )
    )
    parser.add_argument("--method-manifest", required=True, help="Path to manifest_geometry_*_{method}.json")
    parser.add_argument("--alpha", type=float, default=None, help="Optional fixed alpha override")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--table-output", default=None, help="Optional TSV output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest_path = os.path.abspath(args.method_manifest)
    metrics = _compute_metrics(manifest_path, args.alpha)

    benchmark_id = str(metrics["benchmark_id"])
    method = str(metrics["method"])
    output_path = args.output
    table_output_path = args.table_output
    if output_path is None:
        output_path = os.path.join(os.path.dirname(manifest_path), f"{benchmark_id}_{method}_fixed_alpha_modular_collapse.png")
    if table_output_path is None:
        table_output_path = os.path.join(os.path.dirname(manifest_path), f"{benchmark_id}_{method}_fixed_alpha_modular_collapse.tsv")

    output_path = os.path.abspath(output_path)
    table_output_path = os.path.abspath(table_output_path)
    _plot_metrics(metrics, output_path)
    _write_table(metrics, table_output_path)

    print(f"wrote {output_path}")
    print(f"wrote {table_output_path}")
    print(f"alpha_fixed={metrics['alpha']:.10f}")
    for item in metrics["metrics"]:
        print(
            "scale="
            f"{int(item['scale'])} family_size={int(item['family_size'])} n_points={int(item['n_points'])} "
            f"lattice=({int(item['Lx'])},{int(item['Ly'])},{int(item['Tx'])},{int(item['Ty'])}) "
            f"n_traj={int(item['n_traj'])} "
            f"mean_abs_residual={float(item['mean_abs_residual']):.8f} "
            f"rms_residual={float(item['rms_residual']):.8f} "
            f"mean_residual={float(item['mean_residual']):.8f}"
        )


if __name__ == "__main__":
    main()