#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import textwrap
from collections import defaultdict
from math import ceil
from typing import Any

import matplotlib
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay, QhullError

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_REPO_ROOT = os.path.normpath(os.path.join(_PROJECT_ROOT, ".."))
_KFC_ROOT = os.path.join(_REPO_ROOT, "K_from_continuum")
if _KFC_ROOT not in sys.path:
    sys.path.insert(0, _KFC_ROOT)

from workflow_common import evaluate_observable_fit, fit_observable_continuum_power  # noqa: E402


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


def _format_lattice(lattice: tuple[int, int, int, int]) -> str:
    lx, ly, tx, ty = [int(value) for value in lattice]
    return f"({lx},{ly},{tx},{ty})"


def _summarize_payload_lattices(payloads: list[dict[str, Any]]) -> str:
    parts: list[str] = []
    for payload in sorted(payloads, key=lambda item: int(item["scale"])):
        lattice = tuple(int(value) for value in payload["lattice"])
        parts.append(f"s{int(payload['scale'])}={_format_lattice(lattice)}")
    return ", ".join(parts)


def _summarize_couplings(manifest: dict[str, Any]) -> str:
    couplings = manifest.get("couplings", {})
    beta_value = manifest.get("beta")
    parts: list[str] = []
    if isinstance(beta_value, (int, float)) and np.isfinite(float(beta_value)):
        parts.append(f"beta={float(beta_value):.6f}")
    if isinstance(couplings, dict):
        if all(name in couplings for name in ("k1", "k2", "k3")):
            parts.append(
                "k=("
                + ",".join(f"{float(couplings[name]):.3f}" for name in ("k1", "k2", "k3"))
                + ")"
            )
        if all(name in couplings for name in ("r1", "r2")):
            parts.append(
                "r=("
                + ",".join(f"{float(couplings[name]):.3f}" for name in ("r1", "r2"))
                + ")"
            )
    return ", ".join(parts) if parts else "couplings unavailable"


def _summarize_fss_cfg(fss_cfg: dict[str, Any]) -> str:
    return (
        f"fit={str(fss_cfg.get('fit_method', 'taylor2'))}, "
        f"c_min={float(fss_cfg.get('c_min', 0.05)):.2f}, "
        f"c_max={float(fss_cfg.get('c_max', 3.5)):.2f}, "
        f"c0={float(fss_cfg.get('c_initial', 1.0)):.2f}, "
        f"min_free_C={int(fss_cfg.get('min_sizes_for_free_C', 8))}"
    )


def _format_fit_parameters(prefix: str, fit_payload: dict[str, Any]) -> str:
    return (
        f"{prefix}: A={float(fit_payload['A']):+.4f}±{float(fit_payload['sigma_A']):.4f}, "
        f"B={float(fit_payload['B']):+.4f}, C={float(fit_payload['C']):+.4f}, "
        f"n={int(fit_payload['n_sizes_used'])}, mode={str(fit_payload['fit_mode'])}"
    )


def _format_method_base_lattice(method_label: str, method_manifest: dict[str, Any]) -> str:
    payloads = sorted(list(method_manifest.get("payloads", [])), key=lambda payload: int(payload["scale"]))
    if len(payloads) == 0:
        return f"{method_label}: base lattice unavailable"

    lattice = list(payloads[0].get("lattice", [0, 0, 0, 0]))
    while len(lattice) < 4:
        lattice.append(0)
    lx, ly, tx, ty = (int(lattice[index]) for index in range(4))

    couplings = dict(method_manifest.get("couplings", {}))
    k1 = float(couplings.get("k1", couplings.get("r1", float("nan"))))
    k2 = float(couplings.get("k2", couplings.get("r2", float("nan"))))
    k3 = float(couplings.get("k3", 1.0))
    return (
        f"{method_label}: Lx={lx} Ly={ly} Tx={tx} Ty={ty} | "
        f"k1={k1:.4f} k2={k2:.4f} k3={k3:.4f}"
    )


def _point_panel_lattice_caption(payload: dict[str, Any]) -> str:
    return f"{payload['untwisted_base_lattice']}\n{payload['twisted_base_lattice']}"


def _lattice_volume(lx: int, ly: int, tx: int, ty: int) -> int:
    return abs(int(lx) * int(ly) + int(tx) * int(ty))


def _sqrt_lattice_volume(lx: int, ly: int, tx: int, ty: int) -> float:
    return float(np.sqrt(float(_lattice_volume(lx, ly, tx, ty))))


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


def _load_fss_cfg(benchmark_manifest: dict[str, Any]) -> dict[str, Any]:
    fss_cfg = {
        "fit_method": "taylor2",
        "c_min": 0.05,
        "c_max": 3.5,
        "c_initial": 1.0,
        "min_sizes_for_free_C": 8,
    }
    config_path = benchmark_manifest.get("config_path")
    if isinstance(config_path, str) and os.path.exists(config_path):
        cfg = _load_json(config_path)
        if isinstance(cfg.get("fss"), dict):
            fss_cfg.update(cfg["fss"])
    return fss_cfg


def _build_periodic_interpolator(
    payload: dict[str, Any],
    *,
    embedding_cycles: tuple[int, int],
) -> dict[str, Any]:
    lx, ly, tx, ty = [int(value) for value in payload["lattice"]]
    rows = _load_dat_rows(str(payload["data_path"]))

    unique_points: dict[tuple[float, float], dict[str, list[float]]] = {}
    for row in rows:
        conn_value = float(row.get("corr_conn", row.get("conn", float("nan"))))
        sigma_value = float(row.get("err_conn", row.get("conn_err", float("nan"))))
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
        point = unique_points.setdefault(key, {"value": [], "sigma": []})
        point["value"].append(conn_value)
        if np.isfinite(sigma_value):
            point["sigma"].append(abs(sigma_value))

    a_tiled: list[float] = []
    b_tiled: list[float] = []
    value_tiled: list[float] = []
    sigma_tiled: list[float] = []
    a_core: list[float] = []
    b_core: list[float] = []
    value_core: list[float] = []
    sigma_core: list[float] = []

    for (a_wrap, b_wrap), values in sorted(unique_points.items()):
        value_mean = float(np.mean(values["value"]))
        sigma_mean = float(np.mean(values["sigma"])) if len(values["sigma"]) > 0 else float("nan")
        a_core.append(a_wrap)
        b_core.append(b_wrap)
        value_core.append(value_mean)
        sigma_core.append(sigma_mean)
        for delta_a in (-1.0, 0.0, 1.0):
            for delta_b in (-1.0, 0.0, 1.0):
                a_tiled.append(a_wrap + delta_a)
                b_tiled.append(b_wrap + delta_b)
                value_tiled.append(value_mean)
                sigma_tiled.append(sigma_mean)

    points_tiled = np.column_stack([np.asarray(a_tiled, dtype=float), np.asarray(b_tiled, dtype=float)])
    value_tiled_array = np.asarray(value_tiled, dtype=float)
    sigma_tiled_array = np.asarray(sigma_tiled, dtype=float)
    try:
        triangulation = Delaunay(points_tiled, qhull_options="QJ Qc Qbb Q12")
        value_interpolator = LinearNDInterpolator(triangulation, value_tiled_array, fill_value=np.nan)
        sigma_interpolator = LinearNDInterpolator(triangulation, sigma_tiled_array, fill_value=np.nan)
    except QhullError:
        value_interpolator = LinearNDInterpolator(points_tiled, value_tiled_array, fill_value=np.nan)
        sigma_interpolator = LinearNDInterpolator(points_tiled, sigma_tiled_array, fill_value=np.nan)
    nearest_value = NearestNDInterpolator(points_tiled, value_tiled_array)
    nearest_sigma = NearestNDInterpolator(points_tiled, sigma_tiled_array)

    return {
        "scale": int(payload["scale"]),
        "family_size": int(payload["family_size"]),
        "lattice": tuple(int(value) for value in payload["lattice"]),
        "value_interpolator": value_interpolator,
        "sigma_interpolator": sigma_interpolator,
        "nearest_value": nearest_value,
        "nearest_sigma": nearest_sigma,
        "a_core": np.asarray(a_core, dtype=float),
        "b_core": np.asarray(b_core, dtype=float),
        "value_core": np.asarray(value_core, dtype=float),
        "sigma_core": np.asarray(sigma_core, dtype=float),
    }


def _evaluate_periodic_interpolator(
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


def _fit_series(series_rows: list[dict[str, Any]], fss_cfg: dict[str, Any]) -> dict[str, Any]:
    sqrt_volumes = np.asarray([float(row["sqrt_volume"]) for row in series_rows], dtype=float)
    values = np.asarray([float(row["value"]) for row in series_rows], dtype=float)
    sigmas = np.asarray([float(row["sigma"]) for row in series_rows], dtype=float)
    A, sigma_A, B, C, n_sizes_used, fit_mode = fit_observable_continuum_power(
        sqrt_volumes,
        values,
        sigmas,
        fit_method=str(fss_cfg.get("fit_method", "taylor2")),
        c_min=float(fss_cfg.get("c_min", 0.05)),
        c_max=float(fss_cfg.get("c_max", 3.5)),
        c_initial=float(fss_cfg.get("c_initial", 1.0)),
        min_sizes_for_free_C=int(fss_cfg.get("min_sizes_for_free_C", 8)),
    )
    return {
        "A": float(A),
        "sigma_A": float(sigma_A),
        "B": float(B),
        "C": float(C),
        "n_sizes_used": int(n_sizes_used),
        "fit_mode": str(fit_mode),
    }


def _panel_title(display_label: str, point_row: dict[str, Any]) -> str:
    return (
        f"{display_label} | id={int(point_row['point_id'])} | "
        f"(m,n)=({int(point_row['m'])},{int(point_row['n'])})\n"
        f"(a,b)=({float(point_row['a_wrap']):.3f}, {float(point_row['b_wrap']):.3f})"
    )


def _display_title(benchmark_manifest: dict[str, Any]) -> str:
    description = str(benchmark_manifest.get("description", "")).strip()
    if description:
        return description
    return str(benchmark_manifest.get("benchmark_id", "benchmark comparison"))


def _selected_point_metrics(summary_rows: list[dict[str, Any]]) -> str:
    metrics = _point_metric_values(summary_rows)
    if metrics is None:
        return "selected-point z summary unavailable"
    return (
        f"selected-point RMS z={metrics['rms_z']:.2f}\n"
        f"mean |z|={metrics['mean_abs_z']:.2f}, max |z|={metrics['max_abs_z']:.2f}"
    )


def _point_metric_values(summary_rows: list[dict[str, Any]]) -> dict[str, float] | None:
    z_values = np.asarray([float(row["z_delta"]) for row in summary_rows if np.isfinite(float(row["z_delta"]))], dtype=float)
    if z_values.size == 0:
        return None
    return {
        "rms_z": float(np.sqrt(np.mean(np.square(z_values)))),
        "mean_abs_z": float(np.mean(np.abs(z_values))),
        "max_abs_z": float(np.max(np.abs(z_values))),
    }


def _selected_point_listing(summary_rows: list[dict[str, Any]]) -> str:
    items = [
        f"{str(row['display_label'])} (id={int(row['point_id'])}): ΔA={float(row['delta_A']):+.4f}, z={float(row['z_delta']):+.2f}"
        for row in summary_rows
    ]
    if len(items) <= 2:
        return "\n".join(items)
    split_index = int(ceil(len(items) / 2.0))
    return "\n".join([";  ".join(items[:split_index]), ";  ".join(items[split_index:])])


def _add_info_box(
    axis: Any,
    *,
    x_value: float,
    title: str,
    body: str,
    facecolor: str,
    edgecolor: str,
) -> None:
    axis.text(
        x_value,
        0.60,
        f"{title}\n{body}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9.2,
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.45,rounding_size=0.12",
            "facecolor": facecolor,
            "edgecolor": edgecolor,
            "linewidth": 1.0,
        },
    )


def _render_info_card(
    axis: Any,
    *,
    title: str,
    rows: list[tuple[str, str]],
    accent: str,
    facecolor: str,
) -> None:
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_facecolor(facecolor)
    axis.add_patch(
        plt.Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=axis.transAxes,
            facecolor=facecolor,
            edgecolor=accent,
            linewidth=1.0,
            zorder=0,
        )
    )
    axis.add_patch(
        plt.Rectangle(
            (0.0, 0.83),
            1.0,
            0.17,
            transform=axis.transAxes,
            facecolor=accent,
            edgecolor="none",
            zorder=1,
        )
    )
    axis.text(
        0.035,
        0.915,
        title,
        transform=axis.transAxes,
        ha="left",
        va="center",
        fontsize=10.0,
        fontweight="bold",
        color="white",
        zorder=2,
    )
    wrapped_rows: list[tuple[str, list[str]]] = []
    for label, value in rows:
        lines = textwrap.wrap(str(value), width=34, break_long_words=False, break_on_hyphens=False)
        if not lines:
            lines = [""]
        wrapped_rows.append((label, lines))
    total_lines = sum(len(lines) for _, lines in wrapped_rows)
    if total_lines == 0:
        return
    body_top = 0.76
    body_bottom = 0.06
    line_step = (body_top - body_bottom) / float(total_lines)
    current_y = body_top - line_step * 0.5
    for label, lines in wrapped_rows:
        axis.text(
            0.035,
            current_y,
            label,
            transform=axis.transAxes,
            ha="left",
            va="center",
            fontsize=8.4,
            color="#5C6670",
            fontweight="bold",
        )
        for line_index, line in enumerate(lines):
            axis.text(
                0.34,
                current_y - line_step * line_index,
                line,
                transform=axis.transAxes,
                ha="left",
                va="center",
                fontsize=8.4,
                color="#1A2330",
                family="monospace",
            )
        current_y -= line_step * len(lines)


def _render_point_panel(
    axis: plt.Axes,
    *,
    payload: dict[str, Any],
    color: str,
    untwisted_color: str,
    twisted_color: str,
    fit_band_alpha: float,
    title_pad: float = 24.0,
    fit_caption_y: float = 1.012,
) -> None:
    display_label = str(payload["display_label"])
    point_row = payload["point_row"]
    untwisted_series = payload["untwisted_series"]
    twisted_series = payload["twisted_series"]
    untwisted_fit = payload["untwisted_fit"]
    twisted_fit = payload["twisted_fit"]
    summary = payload["summary"]

    untwisted_sqrt_volume = np.asarray([float(row["sqrt_volume"]) for row in untwisted_series], dtype=float)
    untwisted_inv_sqrt_volume = 1.0 / untwisted_sqrt_volume
    untwisted_values = np.asarray([float(row["value"]) for row in untwisted_series], dtype=float)
    untwisted_sigmas = np.asarray([float(row["sigma"]) for row in untwisted_series], dtype=float)

    twisted_sqrt_volume = np.asarray([float(row["sqrt_volume"]) for row in twisted_series], dtype=float)
    twisted_inv_sqrt_volume = 1.0 / twisted_sqrt_volume
    twisted_values = np.asarray([float(row["value"]) for row in twisted_series], dtype=float)
    twisted_sigmas = np.asarray([float(row["sigma"]) for row in twisted_series], dtype=float)

    axis.errorbar(
        untwisted_inv_sqrt_volume,
        untwisted_values,
        yerr=untwisted_sigmas,
        fmt="o",
        capsize=3,
        color=untwisted_color,
        ecolor=untwisted_color,
        markersize=6.0,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.1,
        zorder=3,
    )
    axis.errorbar(
        twisted_inv_sqrt_volume,
        twisted_values,
        yerr=twisted_sigmas,
        fmt="s",
        capsize=3,
        color=twisted_color,
        ecolor=twisted_color,
        markersize=5.7,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.1,
        zorder=3,
    )

    max_inv_sqrt_volume = max(
        float(np.max(untwisted_inv_sqrt_volume)),
        float(np.max(twisted_inv_sqrt_volume)),
    )
    x_fit = np.linspace(0.0, max_inv_sqrt_volume * 1.05, 300)
    axis.plot(
        x_fit,
        evaluate_observable_fit(x_fit, untwisted_fit["A"], untwisted_fit["B"], untwisted_fit["C"], untwisted_fit["fit_mode"]),
        color=untwisted_color,
        linewidth=2.0,
        zorder=2,
    )
    axis.plot(
        x_fit,
        evaluate_observable_fit(x_fit, twisted_fit["A"], twisted_fit["B"], twisted_fit["C"], twisted_fit["fit_mode"]),
        color=twisted_color,
        linewidth=2.0,
        linestyle="--",
        zorder=2,
    )
    axis.axhline(untwisted_fit["A"], color=untwisted_color, linewidth=1.0, alpha=0.85)
    axis.axhspan(
        untwisted_fit["A"] - untwisted_fit["sigma_A"],
        untwisted_fit["A"] + untwisted_fit["sigma_A"],
        color=untwisted_color,
        alpha=fit_band_alpha,
    )
    axis.axhline(twisted_fit["A"], color=twisted_color, linewidth=1.0, alpha=0.85, linestyle=":")
    axis.axhspan(
        twisted_fit["A"] - twisted_fit["sigma_A"],
        twisted_fit["A"] + twisted_fit["sigma_A"],
        color=twisted_color,
        alpha=fit_band_alpha,
    )

    panel_title_main = _panel_title(display_label, point_row)
    fit_caption = (
        f"{_format_fit_parameters('untw', untwisted_fit)}\n"
        f"{_format_fit_parameters('twst', twisted_fit)}"
    )
    axis.set_title(panel_title_main, fontsize=11, pad=title_pad)
    axis.text(
        0.5,
        fit_caption_y,
        fit_caption,
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.4,
        family="monospace",
        color="#42505E",
    )
    axis.text(
        0.025,
        0.965,
        f"{display_label}   ΔA = {summary['delta_A']:+.4f}   z = {summary['z_delta']:+.2f}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10.0,
        fontweight="bold",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.35,rounding_size=0.10",
            "facecolor": color,
            "edgecolor": "none",
        },
    )
    axis.set_xlabel("1 / √(lattice volume)")
    axis.set_ylabel("Connected correlator")
    axis.grid(True, alpha=0.42)
    axis.set_xlim(0.0, max_inv_sqrt_volume * 1.08)


def _write_point_panels(
    *,
    panel_output_dir: str,
    benchmark_id: str,
    fit_payloads: list[dict[str, Any]],
    point_palette: list[str],
    untwisted_color: str,
    twisted_color: str,
    neutral_grid: str,
    fit_band_alpha: float,
) -> list[str]:
    os.makedirs(panel_output_dir, exist_ok=True)
    written_paths: list[str] = []
    style_context = {
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "axes.facecolor": "#FBFCFD",
        "figure.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": neutral_grid,
        "grid.linewidth": 0.7,
    }
    with plt.rc_context(style_context):
        for index, payload in enumerate(fit_payloads):
            point_id = int(payload["summary"]["point_id"])
            display_label = str(payload["display_label"])
            color = point_palette[index % len(point_palette)]
            fig = plt.figure(figsize=(8.8, 6.1), constrained_layout=True)
            grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.0, 0.14])
            axis = fig.add_subplot(grid[0, 0])
            axis_footer = fig.add_subplot(grid[1, 0])
            axis_footer.axis("off")
            _render_point_panel(
                axis,
                payload=payload,
                color=color,
                untwisted_color=untwisted_color,
                twisted_color=twisted_color,
                fit_band_alpha=fit_band_alpha,
                title_pad=20.0,
                fit_caption_y=1.018,
            )
            axis_footer.text(
                0.0,
                0.92,
                _point_panel_lattice_caption(payload),
                transform=axis_footer.transAxes,
                ha="left",
                va="top",
                fontsize=7.6,
                family="monospace",
                color="#42505E",
            )
            output_path = os.path.join(
                panel_output_dir,
                f"{benchmark_id}_{display_label}_point{point_id:03d}_fss.png",
            )
            fig.savefig(output_path, dpi=200)
            plt.close(fig)
            written_paths.append(output_path)
    return written_paths


def _write_summary_table(rows: list[dict[str, Any]], table_path: str) -> None:
    with open(table_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "display_label",
                "point_id",
                "m",
                "n",
                "a_wrap",
                "b_wrap",
                "untwisted_A",
                "untwisted_sigma_A",
                "twisted_interp_A",
                "twisted_interp_sigma_A",
                "delta_A",
                "combined_sigma",
                "z_delta",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    str(row["display_label"]),
                    int(row["point_id"]),
                    int(row["m"]),
                    int(row["n"]),
                    f"{float(row['a_wrap']):.8f}",
                    f"{float(row['b_wrap']):.8f}",
                    f"{float(row['untwisted_A']):.8f}",
                    f"{float(row['untwisted_sigma_A']):.8f}",
                    f"{float(row['twisted_A']):.8f}",
                    f"{float(row['twisted_sigma_A']):.8f}",
                    f"{float(row['delta_A']):.8f}",
                    f"{float(row['combined_sigma']):.8f}",
                    f"{float(row['z_delta']):.8f}",
                ]
            )


def _plot_benchmark(
    benchmark_manifest_path: str,
    *,
    output_path: str | None,
    table_output_path: str | None,
    panel_output_dir: str | None,
    point_ids: list[int] | None,
    n_panels: int,
) -> tuple[str, str, list[int]]:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    fss_cfg = _load_fss_cfg(benchmark_manifest)
    untwisted_manifest = _load_json(str(benchmark_manifest["methods"]["untwisted"]))
    twisted_manifest = _load_json(str(benchmark_manifest["methods"]["twisted"]))
    untwisted_base_lattice = _format_method_base_lattice("untw base", untwisted_manifest)
    twisted_base_lattice = _format_method_base_lattice("twst base", twisted_manifest)

    untwisted_raw_rows = _load_dat_rows(str(untwisted_manifest["pointwise_raw"]))
    untwisted_continuum_rows = _load_dat_rows(str(untwisted_manifest["pointwise_continuum"]))
    untwisted_point_rows = _load_dat_rows(str(untwisted_manifest["smallest_lattice_points"]))

    untwisted_raw_by_id = _group_raw_rows(untwisted_raw_rows)
    untwisted_point_by_id = {int(row["point_id"]): row for row in untwisted_point_rows}

    selected_point_ids = point_ids or _select_representative_points(untwisted_continuum_rows, n_panels=n_panels)
    selected_point_ids = [
        point_id for point_id in selected_point_ids if point_id in untwisted_raw_by_id and point_id in untwisted_point_by_id
    ]
    if len(selected_point_ids) == 0:
        raise ValueError("no selected points were present in the untwisted raw table")

    twisted_embedding_cycles = tuple(int(value) for value in twisted_manifest.get("embedding_cycles", [0, 1]))
    twisted_payload_interpolators = [
        _build_periodic_interpolator(payload, embedding_cycles=twisted_embedding_cycles)
        for payload in sorted(twisted_manifest["payloads"], key=lambda payload: int(payload["scale"]))
    ]

    all_point_ids = [
        point_id
        for point_id, point_row in sorted(
            untwisted_point_by_id.items(),
            key=lambda item: (float(item[1]["b_wrap"]), float(item[1]["a_wrap"]), int(item[0])),
        )
        if point_id in untwisted_raw_by_id
    ]
    if len(all_point_ids) == 0:
        raise ValueError("no base-lattice points were available for continuum comparison")

    point_payloads: dict[int, dict[str, Any]] = {}
    all_summary_rows: list[dict[str, Any]] = []
    for point_id in all_point_ids:
        point_row = untwisted_point_by_id[point_id]
        a_wrap = float(point_row["a_wrap"])
        b_wrap = float(point_row["b_wrap"])

        untwisted_series = [
            {
                "scale": int(row["scale"]),
                "family_size": int(row["family_size"]),
                "sqrt_volume": _sqrt_lattice_volume(int(row["Lx"]), int(row["Ly"]), int(row["Tx"]), int(row["Ty"])),
                "value": float(row["conn"]),
                "sigma": float(row["conn_err"]),
            }
            for row in untwisted_raw_by_id[point_id]
        ]
        twisted_series: list[dict[str, Any]] = []
        for payload_interp in twisted_payload_interpolators:
            twisted_value = _evaluate_periodic_interpolator(
                payload_interp["value_interpolator"],
                payload_interp["nearest_value"],
                a_wrap,
                b_wrap,
                a_core=payload_interp["a_core"],
                b_core=payload_interp["b_core"],
                z_core=payload_interp["value_core"],
            )
            twisted_sigma = _evaluate_periodic_interpolator(
                payload_interp["sigma_interpolator"],
                payload_interp["nearest_sigma"],
                a_wrap,
                b_wrap,
                a_core=payload_interp["a_core"],
                b_core=payload_interp["b_core"],
                z_core=payload_interp["sigma_core"],
            )
            twisted_series.append(
                {
                    "scale": int(payload_interp["scale"]),
                    "family_size": int(payload_interp["family_size"]),
                    "sqrt_volume": _sqrt_lattice_volume(*payload_interp["lattice"]),
                    "value": float(twisted_value),
                    "sigma": float(twisted_sigma),
                }
            )

        untwisted_fit = _fit_series(untwisted_series, fss_cfg)
        twisted_fit = _fit_series(twisted_series, fss_cfg)
        combined_sigma = float(np.sqrt(max(untwisted_fit["sigma_A"], 0.0) ** 2 + max(twisted_fit["sigma_A"], 0.0) ** 2))
        delta = float(untwisted_fit["A"] - twisted_fit["A"])
        z_delta = delta / combined_sigma if np.isfinite(combined_sigma) and combined_sigma > 0.0 else float("nan")

        summary_row = {
            "display_label": f"pt{int(point_id):03d}",
            "point_id": int(point_id),
            "m": int(point_row["m"]),
            "n": int(point_row["n"]),
            "a_wrap": float(point_row["a_wrap"]),
            "b_wrap": float(point_row["b_wrap"]),
            "untwisted_A": float(untwisted_fit["A"]),
            "untwisted_sigma_A": float(untwisted_fit["sigma_A"]),
            "twisted_A": float(twisted_fit["A"]),
            "twisted_sigma_A": float(twisted_fit["sigma_A"]),
            "delta_A": delta,
            "combined_sigma": combined_sigma,
            "z_delta": z_delta,
        }
        all_summary_rows.append(summary_row)
        point_payloads[int(point_id)] = {
            "display_label": str(summary_row["display_label"]),
            "point_row": point_row,
            "untwisted_series": untwisted_series,
            "twisted_series": twisted_series,
            "untwisted_fit": untwisted_fit,
            "twisted_fit": twisted_fit,
            "summary": summary_row,
            "untwisted_base_lattice": untwisted_base_lattice,
            "twisted_base_lattice": twisted_base_lattice,
        }

    fit_payloads: list[dict[str, Any]] = []
    all_fit_payloads = [point_payloads[int(point_id)] for point_id in all_point_ids]
    selected_summary_rows: list[dict[str, Any]] = []
    for display_index, point_id in enumerate(selected_point_ids, start=1):
        payload = dict(point_payloads[int(point_id)])
        display_label = f"p{display_index}"
        summary_row = dict(payload["summary"])
        summary_row["display_label"] = display_label
        payload["display_label"] = display_label
        payload["summary"] = summary_row
        fit_payloads.append(payload)
        selected_summary_rows.append(summary_row)

    tau_real = float(untwisted_manifest["target_tau"]["real"])
    tau_imag = float(untwisted_manifest["target_tau"]["imag"])
    outline_x, outline_y = _unit_cell_outline(tau_real, tau_imag)

    benchmark_id = str(benchmark_manifest["benchmark_id"])
    untwisted_lattice_summary = _summarize_payload_lattices(list(untwisted_manifest["payloads"]))
    twisted_lattice_summary = _summarize_payload_lattices(list(twisted_manifest["payloads"]))
    untwisted_coupling_summary = _summarize_couplings(untwisted_manifest)
    twisted_coupling_summary = _summarize_couplings(twisted_manifest)
    fss_cfg_summary = _summarize_fss_cfg(fss_cfg)

    point_palette = ["#1D4E89", "#D97706", "#2F855A", "#C05621", "#7B2CBF", "#B83280"]
    untwisted_color = "#1F5A91"
    twisted_color = "#C66A2B"
    neutral_grid = "#D4DBE3"
    fit_band_alpha = 0.10

    style_context = {
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "axes.facecolor": "#FBFCFD",
        "figure.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": neutral_grid,
        "grid.linewidth": 0.7,
    }

    with plt.rc_context(style_context):
        fig_width = 16.5
        fig_height = 9.2
        fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
        fig.get_layout_engine().set(w_pad=0.12, h_pad=0.18, hspace=0.05, wspace=0.10)

        grid = fig.add_gridspec(
            nrows=3,
            ncols=6,
            height_ratios=[0.55, 1.30, 2.25],
        )

        ax_title = fig.add_subplot(grid[0, :])
        ax_title.axis("off")
        ax_info_un = fig.add_subplot(grid[1, 0:2])
        ax_info_tw = fig.add_subplot(grid[1, 2:4])
        ax_info_fit = fig.add_subplot(grid[1, 4:6])
        ax_map = fig.add_subplot(grid[2, 0:2])
        ax_summary = fig.add_subplot(grid[2, 2:6])

        all_metrics = _point_metric_values(all_summary_rows)
        selected_metrics = _point_metric_values(selected_summary_rows)
        all_rms_text = f"{all_metrics['rms_z']:.2f}" if all_metrics is not None else "n/a"
        all_abs_text = (
            f"mean={all_metrics['mean_abs_z']:.2f}, max={all_metrics['max_abs_z']:.2f}"
            if all_metrics is not None
            else "n/a"
        )
        summary_text_lines: list[str] = []
        if all_metrics is not None:
            summary_text_lines.append(f"all-point RMS z={all_metrics['rms_z']:.2f}")
            summary_text_lines.append(
                f"mean |z|={all_metrics['mean_abs_z']:.2f}, max |z|={all_metrics['max_abs_z']:.2f}"
            )
        if selected_metrics is not None:
            summary_text_lines.append(f"selected-panel RMS z={selected_metrics['rms_z']:.2f}")
            summary_text_lines.append(
                f"selected mean |z|={selected_metrics['mean_abs_z']:.2f}, max |z|={selected_metrics['max_abs_z']:.2f}"
            )
        summary_text = "\n".join(summary_text_lines) if summary_text_lines else "z summary unavailable"

        raw_title = _display_title(benchmark_manifest)
        wrapped_title = "\n".join(textwrap.wrap(raw_title, width=110)) or raw_title
        subtitle_text = f"{benchmark_id}   ·   Shared-coordinate continuum comparison on the untwisted target torus"

        ax_title.text(
            0.0,
            0.95,
            wrapped_title,
            transform=ax_title.transAxes,
            ha="left",
            va="top",
            fontsize=15,
            fontweight="bold",
            color="#1A2330",
        )
        ax_title.text(
            0.0,
            0.18,
            subtitle_text,
            transform=ax_title.transAxes,
            ha="left",
            va="bottom",
            fontsize=10.0,
            color="#5C6670",
        )

        legend_handles = [
            Line2D([], [], color=untwisted_color, marker="o", linestyle="-", linewidth=2.0, markersize=6.5, label="Untwisted MC + fit"),
            Line2D([], [], color=twisted_color, marker="s", linestyle="--", linewidth=2.0, markersize=6.0, label="Twisted interp + fit"),
        ]
        ax_title.legend(
            handles=legend_handles,
            loc="lower right",
            bbox_to_anchor=(1.0, 0.0),
            frameon=True,
            facecolor="white",
            edgecolor="#D5DBE3",
            ncol=2,
            handlelength=2.4,
            columnspacing=1.6,
            fontsize=10.0,
        )

        _render_info_card(
            ax_info_un,
            title="UNTWISTED TEST",
            rows=[
                ("couplings", untwisted_coupling_summary),
                ("lattices", untwisted_lattice_summary),
            ],
            accent=untwisted_color,
            facecolor="#F2F7FC",
        )
        _render_info_card(
            ax_info_tw,
            title="TWISTED REFERENCE",
            rows=[
                ("couplings", twisted_coupling_summary),
                ("source", twisted_lattice_summary),
            ],
            accent=twisted_color,
            facecolor="#FCF4ED",
        )
        _render_info_card(
            ax_info_fit,
            title="FIT & DISCREPANCY",
            rows=[
                ("FSS", fss_cfg_summary),
                ("τ", f"{tau_real:.3f} + {tau_imag:.3f} i"),
                ("all RMS z", all_rms_text),
                ("all |z|", all_abs_text),
                ("panels", f"{len(all_fit_payloads)} exported"),
                ("highlights", _selected_point_listing(selected_summary_rows).replace("\n", "  ·  ")),
            ],
            accent="#5A6470",
            facecolor="#F4F6F9",
        )

        all_x = np.asarray([float(row["nu_real"]) for row in untwisted_point_rows], dtype=float)
        all_y = np.asarray([float(row["nu_imag"]) for row in untwisted_point_rows], dtype=float)
        ax_map.fill(outline_x, outline_y, color="#F7F2E9", alpha=0.70, zorder=0)
        ax_map.scatter(all_x, all_y, s=13, color="#C8CFD8", alpha=0.60, linewidths=0.0, zorder=1)
        ax_map.plot(outline_x, outline_y, color="#222222", linewidth=1.2, zorder=2)
        for index, payload in enumerate(fit_payloads):
            color = point_palette[index % len(point_palette)]
            display_label = str(payload["display_label"])
            point_row = payload["point_row"]
            x_value = float(point_row["nu_real"])
            y_value = float(point_row["nu_imag"])
            ax_map.scatter(
                [x_value],
                [y_value],
                s=62,
                color=color,
                edgecolors="white",
                linewidths=0.9,
                zorder=3,
            )
            ax_map.text(
                x_value + 0.015,
                y_value + 0.015,
                display_label,
                fontsize=9.5,
                fontweight="bold",
                color=color,
                va="bottom",
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
            )
        ax_map.set_aspect("equal", adjustable="box")
        ax_map.set_title("Base-lattice coordinates on the target torus")
        ax_map.set_xlabel("Re(ν)")
        ax_map.set_ylabel("Im(ν)")
        ax_map.grid(True, alpha=0.35)
        x_pad = 0.08 * (float(np.max(outline_x)) - float(np.min(outline_x)) + 1.0e-12)
        y_pad = 0.08 * (float(np.max(outline_y)) - float(np.min(outline_y)) + 1.0e-12)
        ax_map.set_xlim(float(np.min(outline_x)) - x_pad, float(np.max(outline_x)) + x_pad)
        ax_map.set_ylim(float(np.min(outline_y)) - y_pad, float(np.max(outline_y)) + y_pad)

        x_summary = np.asarray([row["twisted_A"] for row in all_summary_rows], dtype=float)
        y_summary = np.asarray([row["untwisted_A"] for row in all_summary_rows], dtype=float)
        lim_min = min(float(np.min(x_summary)), float(np.min(y_summary)))
        lim_max = max(float(np.max(x_summary)), float(np.max(y_summary)))
        padding = 0.05 * (lim_max - lim_min if lim_max > lim_min else 1.0)
        ax_summary.plot(
            [lim_min - padding, lim_max + padding],
            [lim_min - padding, lim_max + padding],
            linestyle="--",
            color="#434A54",
            linewidth=1.2,
            zorder=1,
        )
        for row in all_summary_rows:
            ax_summary.scatter(
                float(row["twisted_A"]),
                float(row["untwisted_A"]),
                s=24,
                color="#AAB7C5",
                alpha=0.58,
                linewidths=0.0,
                zorder=2,
            )
        for index, row in enumerate(selected_summary_rows):
            color = point_palette[index % len(point_palette)]
            ax_summary.errorbar(
                float(row["twisted_A"]),
                float(row["untwisted_A"]),
                xerr=float(row["twisted_sigma_A"]),
                yerr=float(row["untwisted_sigma_A"]),
                fmt="o",
                color=color,
                ecolor=color,
                elinewidth=1.2,
                capsize=3.0,
                markersize=6.5,
                markeredgecolor="white",
                markeredgewidth=0.8,
                zorder=2,
            )
            ax_summary.text(
                float(row["twisted_A"]) + 0.0015,
                float(row["untwisted_A"]) + 0.0015,
                str(row["display_label"]),
                fontsize=10,
                fontweight="bold",
                color=color,
            )
        ax_summary.text(
            0.03,
            0.97,
            summary_text,
            transform=ax_summary.transAxes,
            ha="left",
            va="top",
            fontsize=9.2,
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "#D5DBE3", "pad": 0.35},
        )
        ax_summary.set_xlim(lim_min - padding, lim_max + padding)
        ax_summary.set_ylim(lim_min - padding, lim_max + padding)
        ax_summary.set_title("Continuum intercepts for all base-lattice coordinates")
        ax_summary.set_xlabel("Twisted interpolated continuum A")
        ax_summary.set_ylabel("Untwisted continuum A")
        ax_summary.grid(True, alpha=0.40)

        if output_path is None:
            output_path = os.path.join(
                os.path.dirname(benchmark_manifest_path),
                f"{benchmark_id}_shared_coordinate_twisted_untwisted_fss.png",
            )
        if table_output_path is None:
            table_output_path = os.path.join(
                os.path.dirname(benchmark_manifest_path),
                f"{benchmark_id}_shared_coordinate_twisted_untwisted_fss.tsv",
            )

        fig.savefig(output_path, dpi=200)
        plt.close(fig)
    panel_output_paths: list[str] = []
    if panel_output_dir is not None:
        panel_output_paths = _write_point_panels(
            panel_output_dir=panel_output_dir,
            benchmark_id=benchmark_id,
            fit_payloads=all_fit_payloads,
            point_palette=point_palette,
            untwisted_color=untwisted_color,
            twisted_color=twisted_color,
            neutral_grid=neutral_grid,
            fit_band_alpha=fit_band_alpha,
        )
    _write_summary_table(all_summary_rows, table_output_path)

    print(f"wrote {output_path}")
    print(f"wrote {table_output_path}")
    if panel_output_paths:
        print(f"wrote {len(panel_output_paths)} point panels under {panel_output_dir}")
    print(f"selected_point_ids={selected_point_ids}")
    return output_path, table_output_path, selected_point_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot untwisted FSS and twisted interpolated FSS at the same untwisted coordinates for a benchmark."
        )
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--table-output", default=None, help="Optional TSV output path")
    parser.add_argument(
        "--panel-output-dir",
        default=None,
        help="Optional directory for separate per-point FSS PNG outputs",
    )
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
        panel_output_dir=os.path.abspath(args.panel_output_dir) if args.panel_output_dir else None,
        point_ids=list(args.point_ids or []),
        n_panels=max(int(args.n_panels), 1),
    )


if __name__ == "__main__":
    main()