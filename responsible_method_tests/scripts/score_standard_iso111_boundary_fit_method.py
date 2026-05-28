#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Any

import matplotlib
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_REPO_ROOT = os.path.normpath(os.path.join(_PROJECT_ROOT, ".."))
_KFC_ROOT = os.path.join(_REPO_ROOT, "K_from_continuum")
if _KFC_ROOT not in sys.path:
    sys.path.insert(0, _KFC_ROOT)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from workflow_common import evaluate_observable_fit, fit_observable_continuum_power  # noqa: E402

from plot_standard_iso111_boundary_midpoint_fss import (
    DEFAULT_FRACTION,
    DEFAULT_SWEEP_ROOT,
    STANDARD_ROOT,
    _aggregate_by_fraction,
    _boundary_point_specs,
    _boundary_ratio_series,
    _build_ratio_payloads,
    _build_series,
    _connected_ratio,
    _connected_value_with_true_jackknife,
    _discover_sweep_dirs,
    _fit_blind_power_model,
    _fraction_slug,
    _load_dat_rows,
    _parse_coupling_from_dir,
    _require_point,
    _sqrt_volume,
    _write_sweep_summary_table,
)


DEFAULT_TWISTED_REFERENCE_ROOT = os.path.join(
    STANDARD_ROOT,
    "data",
    "iso111",
    "twisted",
    "reference",
)
DEFAULT_FIT_FAMILY = "continuum"
DEFAULT_CONTINUUM_FIT_METHOD = "taylor2"
DEFAULT_C_MIN = 0.05
DEFAULT_C_MAX = 3.5
DEFAULT_C_INITIAL = 1.0
DEFAULT_MIN_SIZES_FOR_FREE_C = 8


def _default_output_path(fraction: float) -> str:
    return os.path.join(STANDARD_ROOT, f"iso111_boundary_{_fraction_slug(fraction)}_fit_method_zscores.tsv")


def _default_reference_fit_output_path(fraction: float) -> str:
    return os.path.join(STANDARD_ROOT, f"iso111_boundary_{_fraction_slug(fraction)}_fit_method_twisted_reference_fits.tsv")


def _parse_reference_lattice(dir_name: str) -> tuple[int, int, int, int]:
    match = re.fullmatch(r"Lx(-?\d+)_Ly(-?\d+)_Tx(-?\d+)_Ty(-?\d+)", dir_name)
    if match is None:
        raise ValueError(f"could not parse twisted reference lattice from {dir_name}")
    return tuple(int(value) for value in match.groups())


def _discover_twisted_reference_entries(reference_root: str) -> list[dict[str, Any]]:
    if not os.path.isdir(reference_root):
        raise FileNotFoundError(f"twisted reference root not found: {reference_root}")

    entries: list[dict[str, Any]] = []
    for name in sorted(os.listdir(reference_root)):
        reference_dir = os.path.join(reference_root, name)
        if not os.path.isdir(reference_dir):
            continue
        dat_path = os.path.join(reference_dir, "two_point_all_to_all.dat")
        if not os.path.isfile(dat_path):
            continue
        lattice = _parse_reference_lattice(name)
        entries.append(
            {
                "name": name,
                "dat_path": dat_path,
                "lattice": lattice,
                "x": 1.0 / _sqrt_volume(lattice),
            }
        )

    entries.sort(key=lambda entry: _sqrt_volume(entry["lattice"]))
    if not entries:
        raise FileNotFoundError(f"no twisted reference all-to-all files found under {reference_root}")
    return entries


def _fit_series_payload(
    *,
    x_values: np.ndarray,
    y_values: np.ndarray,
    sigma_values: np.ndarray,
    fit_family: str,
    continuum_fit_method: str,
    c_min: float,
    c_max: float,
    c_initial: float,
    min_sizes_for_free_c: int,
) -> dict[str, Any]:
    fit_family = str(fit_family).strip().lower()
    if fit_family == "blind_power":
        legacy_fit = _fit_blind_power_model(x_values, y_values, sigma_values)
        return {
            "A": float(legacy_fit["A"]),
            "sigma_A": float(legacy_fit["sigma_A"]),
            "B": float(legacy_fit["B"]),
            "C": float(legacy_fit["omega"]),
            "sigma_C": float(legacy_fit["sigma_omega"]),
            "n_used": int(legacy_fit["n_used"]),
            "fit_mode": str(legacy_fit["fit_mode"]),
            "fit_family": fit_family,
        }
    if fit_family != "continuum":
        raise ValueError(f"unknown fit family: {fit_family}")

    inv_sqrt_volume = np.asarray(x_values, dtype=float)
    sqrt_volume = np.full_like(inv_sqrt_volume, np.nan, dtype=float)
    np.reciprocal(inv_sqrt_volume, out=sqrt_volume, where=inv_sqrt_volume != 0.0)
    A_value, sigma_A, B_value, C_value, n_used, fit_mode = fit_observable_continuum_power(
        sqrt_volume,
        np.asarray(y_values, dtype=float),
        np.asarray(sigma_values, dtype=float),
        fit_method=str(continuum_fit_method),
        c_min=float(c_min),
        c_max=float(c_max),
        c_initial=float(c_initial),
        min_sizes_for_free_C=int(min_sizes_for_free_c),
    )
    return {
        "A": float(A_value),
        "sigma_A": float(sigma_A),
        "B": float(B_value),
        "C": float(C_value),
        "sigma_C": float("nan"),
        "n_used": int(n_used),
        "fit_mode": str(fit_mode),
        "fit_family": fit_family,
    }


def _apply_fit_family(
    payloads: list[dict[str, Any]],
    *,
    fit_family: str,
    continuum_fit_method: str,
    c_min: float,
    c_max: float,
    c_initial: float,
    min_sizes_for_free_c: int,
) -> None:
    for payload in payloads:
        payload["fit"] = _fit_series_payload(
            x_values=np.asarray(payload["x"], dtype=float),
            y_values=np.asarray(payload["y"], dtype=float),
            sigma_values=np.asarray(payload["sigma"], dtype=float),
            fit_family=fit_family,
            continuum_fit_method=continuum_fit_method,
            c_min=c_min,
            c_max=c_max,
            c_initial=c_initial,
            min_sizes_for_free_c=min_sizes_for_free_c,
        )


def _build_twisted_point_payloads(
    *,
    point_specs: list[dict[str, Any]],
    twisted_entries: list[dict[str, Any]],
    fit_family: str,
    continuum_fit_method: str,
    c_min: float,
    c_max: float,
    c_initial: float,
    min_sizes_for_free_c: int,
    allow_regeneration: bool = False,
) -> list[dict[str, Any]]:
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] = {}
    value_cache: dict[tuple[str, int, int], tuple[float, float]] = {}
    aggregated_by_dat: dict[str, dict[tuple[float, float], dict[str, Any]]] = {}

    for entry in twisted_entries:
        aggregated_by_dat[str(entry["dat_path"])] = _aggregate_by_fraction(
            _load_dat_rows(str(entry["dat_path"])),
            entry["lattice"],
            embedding_cycles=(0, 1),
        )

    payloads: list[dict[str, Any]] = []
    for point_spec in point_specs:
        x_values: list[float] = []
        y_values: list[float] = []
        sigma_values: list[float] = []
        lattice_labels: list[str] = []
        entry_points: list[dict[str, Any]] = []
        representative_payload: dict[str, Any] | None = None

        for entry in twisted_entries:
            dat_path = str(entry["dat_path"])
            payload = _require_point(
                aggregated_by_dat[dat_path],
                key=point_spec["key"],
                dat_path=dat_path,
            )
            if representative_payload is None:
                representative_payload = payload
            conn_value, conn_sigma = _connected_value_with_true_jackknife(
                dat_path,
                point_m=int(payload["m"]),
                point_n=int(payload["n"]),
                sample_cache=sample_cache,
                value_cache=value_cache,
                allow_regeneration=allow_regeneration,
                fallback_value=(float(payload["value"]), float(payload["sigma"])),
            )
            x_values.append(float(entry["x"]))
            y_values.append(float(conn_value))
            sigma_values.append(float(conn_sigma))
            lattice_labels.append(f"L={int(entry['lattice'][0])}")
            entry_points.append(
                {
                    "dat_path": dat_path,
                    "lattice": entry["lattice"],
                    "m": int(payload["m"]),
                    "n": int(payload["n"]),
                }
            )

        if representative_payload is None:
            raise ValueError(f"no twisted payload found for point {point_spec['key']}")
        fit_payload = _fit_series_payload(
            x_values=np.asarray(x_values, dtype=float),
            y_values=np.asarray(y_values, dtype=float),
            sigma_values=np.asarray(sigma_values, dtype=float),
            fit_family=fit_family,
            continuum_fit_method=continuum_fit_method,
            c_min=c_min,
            c_max=c_max,
            c_initial=c_initial,
            min_sizes_for_free_c=min_sizes_for_free_c,
        )
        payloads.append(
            {
                **point_spec,
                "x": np.asarray(x_values, dtype=float),
                "y": np.asarray(y_values, dtype=float),
                "sigma": np.asarray(sigma_values, dtype=float),
                "lattice_labels": lattice_labels,
                "entry_points": entry_points,
                "m": int(point_spec.get("m", representative_payload["m"])),
                "n": int(point_spec.get("n", representative_payload["n"])),
                "a_wrap": float(point_spec.get("a_wrap", representative_payload["a_wrap"])),
                "b_wrap": float(point_spec.get("b_wrap", representative_payload["b_wrap"])),
                "y_axis_label": str(point_spec.get("y_axis_label", "Connected correlator")),
                "fit": fit_payload,
            }
        )
    return payloads


def _build_twisted_ratio_payloads(
    point_payloads: list[dict[str, Any]],
    *,
    ratio_series: list[dict[str, str]],
    fit_family: str,
    continuum_fit_method: str,
    c_min: float,
    c_max: float,
    c_initial: float,
    min_sizes_for_free_c: int,
    allow_regeneration: bool = False,
) -> list[dict[str, Any]]:
    point_by_label = {str(payload["short_label"]): payload for payload in point_payloads}
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] = {}
    value_cache: dict[tuple[str, int, int], tuple[float, float]] = {}
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]] = {}
    ratio_payloads: list[dict[str, Any]] = []

    for ratio_cfg in ratio_series:
        numerator = point_by_label[str(ratio_cfg["numerator"])]
        denominator = point_by_label[str(ratio_cfg["denominator"])]
        x_values = np.asarray(numerator["x"], dtype=float)
        if x_values.shape != np.asarray(denominator["x"], dtype=float).shape or not np.allclose(x_values, denominator["x"]):
            raise ValueError(f"mismatched twisted x grids for ratio panel {ratio_cfg['short_label']}")

        y_values: list[float] = []
        sigma_values: list[float] = []
        for index, (numerator_point, denominator_point) in enumerate(zip(numerator["entry_points"], denominator["entry_points"])):
            numerator_dat_path = str(numerator_point["dat_path"])
            denominator_dat_path = str(denominator_point["dat_path"])
            if numerator_dat_path != denominator_dat_path:
                raise ValueError(
                    f"mismatched twisted data paths for ratio panel {ratio_cfg['short_label']}: {numerator_dat_path} vs {denominator_dat_path}"
                )
            ratio_value, ratio_sigma = _connected_ratio(
                dat_path=numerator_dat_path,
                point_m=int(numerator_point["m"]),
                point_n=int(numerator_point["n"]),
                anchor_m=int(denominator_point["m"]),
                anchor_n=int(denominator_point["n"]),
                point_value=float(numerator["y"][index]),
                point_sigma=float(numerator["sigma"][index]),
                anchor_value=float(denominator["y"][index]),
                anchor_sigma=float(denominator["sigma"][index]),
                sample_cache=sample_cache,
                value_cache=value_cache,
                ratio_cache=ratio_cache,
                allow_regeneration=allow_regeneration,
            )
            y_values.append(float(ratio_value))
            sigma_values.append(float(ratio_sigma))

        fit_payload = _fit_series_payload(
            x_values=x_values,
            y_values=np.asarray(y_values, dtype=float),
            sigma_values=np.asarray(sigma_values, dtype=float),
            fit_family=fit_family,
            continuum_fit_method=continuum_fit_method,
            c_min=c_min,
            c_max=c_max,
            c_initial=c_initial,
            min_sizes_for_free_c=min_sizes_for_free_c,
        )
        ratio_payloads.append(
            {
                **ratio_cfg,
                "x": x_values,
                "y": np.asarray(y_values, dtype=float),
                "sigma": np.asarray(sigma_values, dtype=float),
                "y_axis_label": "Correlator ratio",
                "fit": fit_payload,
            }
        )
    return ratio_payloads


def _fit_method_summary_row(
    *,
    untwisted_dir: str,
    untwisted_payloads: list[dict[str, Any]],
    twisted_payload_by_label: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    r1_value, r2_value = _parse_coupling_from_dir(untwisted_dir)
    row: dict[str, Any] = {
        "coupling_tag": os.path.basename(os.path.normpath(untwisted_dir)),
        "r1": float(r1_value),
        "r2": float(r2_value),
    }
    corr_chi2 = 0.0
    ratio_chi2 = 0.0
    for payload in untwisted_payloads:
        label = str(payload["short_label"])
        twisted_payload = twisted_payload_by_label[label]
        untwisted_fit = payload["fit"]
        twisted_fit = twisted_payload["fit"]
        delta = float(untwisted_fit["A"]) - float(twisted_fit["A"])
        combined_sigma = float(
            math.sqrt(max(float(untwisted_fit["sigma_A"]), 0.0) ** 2 + max(float(twisted_fit["sigma_A"]), 0.0) ** 2)
        )
        z_value = delta / combined_sigma if np.isfinite(combined_sigma) and combined_sigma > 0.0 else float("nan")
        chi2_value = float(z_value ** 2)
        row[f"{label}_z"] = float(z_value)
        row[f"{label}_chi2"] = float(chi2_value)
        if str(payload.get("y_axis_label", "Connected correlator")) == "Correlator ratio":
            ratio_chi2 += chi2_value
        else:
            corr_chi2 += chi2_value
    row["corr_chi2"] = float(corr_chi2)
    row["ratio_chi2"] = float(ratio_chi2)
    row["chi2_sum"] = float(corr_chi2 + ratio_chi2)
    row["corr_z2"] = float(corr_chi2)
    row["ratio_z2"] = float(ratio_chi2)
    row["z2_sum"] = float(corr_chi2 + ratio_chi2)
    return row


def _evaluate_fit_payload(x_values: np.ndarray, fit_payload: dict[str, Any]) -> np.ndarray:
    x_array = np.asarray(x_values, dtype=float)
    if str(fit_payload.get("fit_family", DEFAULT_FIT_FAMILY)) == "blind_power":
        return float(fit_payload["A"]) + float(fit_payload["B"]) * np.power(x_array, float(fit_payload["C"]))
    return evaluate_observable_fit(
        x_array,
        float(fit_payload["A"]),
        float(fit_payload["B"]),
        float(fit_payload["C"]),
        str(fit_payload["fit_mode"]),
    )


def _render_case_plot(
    *,
    output_path: str,
    fraction: float,
    coupling_tag: str,
    summary_row: dict[str, Any],
    untwisted_payloads: list[dict[str, Any]],
    twisted_payload_by_label: dict[str, dict[str, Any]],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(18, 9.8), squeeze=False)
    axes_flat = list(axes.ravel())
    fig.suptitle(
        f"Iso111 boundary fit_method | fraction={fraction:.3f} | {coupling_tag} | chi2_sum={float(summary_row['chi2_sum']):.3f}",
        fontsize=15,
        y=0.98,
    )
    fit_family = str(untwisted_payloads[0]["fit"].get("fit_family", DEFAULT_FIT_FAMILY)) if untwisted_payloads else DEFAULT_FIT_FAMILY
    fig.text(
        0.5,
        0.945,
        f"untwisted circles/solid fit vs twisted squares/dashed fit | fit_family={fit_family}",
        ha="center",
        va="top",
        fontsize=10,
        color="#444444",
    )

    for axis, payload in zip(axes_flat, untwisted_payloads):
        twisted_payload = twisted_payload_by_label[str(payload["short_label"])]
        color = str(payload["color"])
        untwisted_fit = payload["fit"]
        twisted_fit = twisted_payload["fit"]

        x_min = min(float(np.min(payload["x"])), float(np.min(twisted_payload["x"])))
        x_max = max(float(np.max(payload["x"])), float(np.max(twisted_payload["x"])))
        x_plot = np.geomspace(x_min * 0.95, x_max * 1.05, 300)

        axis.errorbar(
            payload["x"],
            payload["y"],
            yerr=payload["sigma"],
            fmt="o",
            color=color,
            ecolor=color,
            capsize=3,
            markersize=5,
            markeredgecolor="white",
            markeredgewidth=0.8,
            linewidth=1.1,
            label="untwisted",
            zorder=3,
        )
        axis.errorbar(
            twisted_payload["x"],
            twisted_payload["y"],
            yerr=twisted_payload["sigma"],
            fmt="s",
            color="#111827",
            ecolor="#111827",
            capsize=3,
            markersize=4.7,
            markeredgecolor="white",
            markeredgewidth=0.8,
            linewidth=1.1,
            label="twisted",
            zorder=3,
        )
        axis.plot(x_plot, _evaluate_fit_payload(x_plot, untwisted_fit), color=color, linewidth=1.8, zorder=2)
        axis.plot(x_plot, _evaluate_fit_payload(x_plot, twisted_fit), color="#111827", linewidth=1.6, linestyle="--", zorder=2)
        axis.axhspan(
            float(untwisted_fit["A"]) - float(untwisted_fit["sigma_A"]),
            float(untwisted_fit["A"]) + float(untwisted_fit["sigma_A"]),
            color=color,
            alpha=0.08,
        )
        axis.axhspan(
            float(twisted_fit["A"]) - float(twisted_fit["sigma_A"]),
            float(twisted_fit["A"]) + float(twisted_fit["sigma_A"]),
            color="#111827",
            alpha=0.05,
        )
        axis.set_xscale("log")
        axis.grid(True, alpha=0.25)
        axis.set_title(str(payload["short_label"]))
        axis.set_xlabel("1 / sqrt(lattice volume)")
        axis.set_ylabel(str(payload.get("y_axis_label", "Connected correlator")))

        z_value = float(summary_row[f"{str(payload['short_label'])}_z"])
        caption = (
            f"untw A={float(untwisted_fit['A']):+.4f}±{float(untwisted_fit['sigma_A']):.4f}\n"
            f"twst A={float(twisted_fit['A']):+.4f}±{float(twisted_fit['sigma_A']):.4f}\n"
            f"z={z_value:+.2f}"
        )
        axis.text(
            0.03,
            0.97,
            caption,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.4,
            family="monospace",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#d1d5db", "alpha": 0.92},
        )

    for axis in axes_flat[len(untwisted_payloads):]:
        axis.axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.985, 0.955), frameon=False)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _write_reference_fit_table(output_path: str, payloads: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    headers = [
        "panel",
        "y_axis_label",
        "fit_family",
        "fit_A",
        "fit_sigma_A",
        "fit_B",
        "fit_C",
        "fit_sigma_C",
        "fit_mode",
        "n_used",
    ]
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for payload in payloads:
            fit = payload["fit"]
            handle.write(
                "\t".join(
                    [
                        str(payload["short_label"]),
                        str(payload.get("y_axis_label", "Connected correlator")),
                        str(fit.get("fit_family", DEFAULT_FIT_FAMILY)),
                        f"{float(fit['A']):.10f}",
                        f"{float(fit['sigma_A']):.10f}",
                        f"{float(fit['B']):.10f}",
                        f"{float(fit['C']):.10f}",
                        f"{float(fit['sigma_C']):.10f}",
                        str(fit["fit_mode"]),
                        str(int(fit["n_used"])),
                    ]
                )
                + "\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute the standard iso111 boundary fit_method score by fitting the twisted reference size ladder and the untwisted candidate ladder separately, then comparing their continuum intercepts panel-by-panel."
        )
    )
    parser.add_argument("--fraction", type=float, default=DEFAULT_FRACTION)
    parser.add_argument("--sweep-root", default=DEFAULT_SWEEP_ROOT)
    parser.add_argument("--twisted-reference-root", default=DEFAULT_TWISTED_REFERENCE_ROOT)
    parser.add_argument("--min-twisted-sizes", type=int, default=3)
    parser.add_argument("--fit-family", choices=["continuum", "blind_power"], default=DEFAULT_FIT_FAMILY)
    parser.add_argument("--continuum-fit-method", choices=["taylor2", "power"], default=DEFAULT_CONTINUUM_FIT_METHOD)
    parser.add_argument("--c-min", type=float, default=DEFAULT_C_MIN)
    parser.add_argument("--c-max", type=float, default=DEFAULT_C_MAX)
    parser.add_argument("--c-initial", type=float, default=DEFAULT_C_INITIAL)
    parser.add_argument("--min-sizes-for-free-c", type=int, default=DEFAULT_MIN_SIZES_FOR_FREE_C)
    parser.add_argument("--output", default=None)
    parser.add_argument("--reference-fit-output", default=None)
    parser.add_argument("--plot-output-dir", default=None)
    parser.add_argument("--plot-coupling-tag", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fraction = float(args.fraction)
    sweep_root = os.path.abspath(args.sweep_root)
    twisted_reference_root = os.path.abspath(args.twisted_reference_root)
    output_path = os.path.abspath(args.output or _default_output_path(fraction))
    reference_fit_output_path = os.path.abspath(args.reference_fit_output or _default_reference_fit_output_path(fraction))

    twisted_entries = _discover_twisted_reference_entries(twisted_reference_root)
    if len(twisted_entries) < int(args.min_twisted_sizes):
        raise SystemExit(
            f"fit_method needs at least {int(args.min_twisted_sizes)} twisted sizes, found only {len(twisted_entries)} under {twisted_reference_root}"
        )

    point_specs = _boundary_point_specs(fraction)
    ratio_series = _boundary_ratio_series(fraction)
    twisted_point_payloads = _build_twisted_point_payloads(
        point_specs=point_specs,
        twisted_entries=twisted_entries,
        fit_family=str(args.fit_family),
        continuum_fit_method=str(args.continuum_fit_method),
        c_min=float(args.c_min),
        c_max=float(args.c_max),
        c_initial=float(args.c_initial),
        min_sizes_for_free_c=int(args.min_sizes_for_free_c),
        allow_regeneration=False,
    )
    twisted_ratio_payloads = _build_twisted_ratio_payloads(
        twisted_point_payloads,
        ratio_series=ratio_series,
        fit_family=str(args.fit_family),
        continuum_fit_method=str(args.continuum_fit_method),
        c_min=float(args.c_min),
        c_max=float(args.c_max),
        c_initial=float(args.c_initial),
        min_sizes_for_free_c=int(args.min_sizes_for_free_c),
        allow_regeneration=False,
    )
    twisted_payloads = [*twisted_point_payloads, *twisted_ratio_payloads]
    twisted_payload_by_label = {str(payload["short_label"]): payload for payload in twisted_payloads}
    _write_reference_fit_table(reference_fit_output_path, twisted_payloads)

    largest_entry = twisted_entries[-1]
    summary_rows: list[dict[str, Any]] = []
    panel_order: list[str] = [f"{str(payload['short_label'])}_z" for payload in twisted_payloads]
    requested_plot_tags = [str(tag) for tag in args.plot_coupling_tag]
    plot_cases: dict[str, dict[str, Any]] = {}
    for untwisted_dir in _discover_sweep_dirs(sweep_root):
        untwisted_point_payloads = _build_series(
            point_specs=point_specs,
            untwisted_dir=untwisted_dir,
            twisted_dat=str(largest_entry["dat_path"]),
            twisted_lattice=largest_entry["lattice"],
            allow_regeneration=False,
        )
        _apply_fit_family(
            untwisted_point_payloads,
            fit_family=str(args.fit_family),
            continuum_fit_method=str(args.continuum_fit_method),
            c_min=float(args.c_min),
            c_max=float(args.c_max),
            c_initial=float(args.c_initial),
            min_sizes_for_free_c=int(args.min_sizes_for_free_c),
        )
        untwisted_ratio_payloads = _build_ratio_payloads(
            untwisted_point_payloads,
            ratio_series=ratio_series,
            allow_regeneration=False,
        )
        _apply_fit_family(
            untwisted_ratio_payloads,
            fit_family=str(args.fit_family),
            continuum_fit_method=str(args.continuum_fit_method),
            c_min=float(args.c_min),
            c_max=float(args.c_max),
            c_initial=float(args.c_initial),
            min_sizes_for_free_c=int(args.min_sizes_for_free_c),
        )
        all_untwisted_payloads = [*untwisted_point_payloads, *untwisted_ratio_payloads]
        summary_row = _fit_method_summary_row(
            untwisted_dir=untwisted_dir,
            untwisted_payloads=all_untwisted_payloads,
            twisted_payload_by_label=twisted_payload_by_label,
        )
        summary_rows.append(summary_row)
        if requested_plot_tags and str(summary_row["coupling_tag"]) in requested_plot_tags:
            plot_cases[str(summary_row["coupling_tag"])] = {
                "summary_row": dict(summary_row),
                "untwisted_payloads": all_untwisted_payloads,
            }

    _write_sweep_summary_table(
        output_path=output_path,
        rows=summary_rows,
        panel_order=panel_order,
    )
    if args.plot_output_dir:
        if not requested_plot_tags:
            raise SystemExit("--plot-output-dir requires at least one --plot-coupling-tag")
        for coupling_tag in requested_plot_tags:
            if coupling_tag not in plot_cases:
                raise SystemExit(f"requested plot coupling tag not found in sweep: {coupling_tag}")
            plot_output_path = os.path.join(
                os.path.abspath(args.plot_output_dir),
                f"iso111_boundary_{_fraction_slug(fraction)}_fit_method_{coupling_tag}.png",
            )
            _render_case_plot(
                output_path=plot_output_path,
                fraction=fraction,
                coupling_tag=coupling_tag,
                summary_row=plot_cases[coupling_tag]["summary_row"],
                untwisted_payloads=plot_cases[coupling_tag]["untwisted_payloads"],
                twisted_payload_by_label=twisted_payload_by_label,
            )
            print(f"wrote {plot_output_path}")
    print(f"wrote {reference_fit_output_path}")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()