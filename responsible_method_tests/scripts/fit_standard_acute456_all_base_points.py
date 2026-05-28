#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
from typing import Any

import numpy as np
from scipy.optimize import curve_fit

from plot_standard_acute456_center_fss import (
    DEFAULT_TWISTED_DAT,
    DEFAULT_UNTWISTED_DIR,
    TWISTED_LATTICE,
    UNTWISTED_LATTICES,
    _aggregate_by_fraction,
    _build_periodic_interpolator,
    _compute_aggregate_target_score,
    _evaluate_periodic,
    _evaluate_power_model_on_x,
    _find_point_by_mn,
    _infer_family_label,
    _load_dat_rows,
    _ratio_with_uncertainty,
    _sqrt_volume,
    _to_ab,
    _untwisted_dat_path,
    _wrap_unit,
)


HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(HERE, ".."))
STANDARD_ROOT = os.path.join(PROJECT_ROOT, "standard")
DEFAULT_OUTPUT = os.path.join(STANDARD_ROOT, "acute456_all_base_points_individual_fits.tsv")


def _collect_all_base_points(
    rows: list[dict[str, Any]],
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    include_origin: bool,
) -> list[dict[str, Any]]:
    points_by_key: dict[tuple[float, float], dict[str, Any]] = {}
    for row in rows:
        a_raw, b_raw = _to_ab(
            int(row["m"]),
            int(row["n"]),
            lattice[0],
            lattice[1],
            lattice[2],
            lattice[3],
            embedding_cycles=embedding_cycles,
        )
        a_wrap = _wrap_unit(a_raw)
        b_wrap = _wrap_unit(b_raw)
        if not include_origin and abs(a_wrap) < 1.0e-12 and abs(b_wrap) < 1.0e-12:
            continue
        key = (round(a_wrap, 12), round(b_wrap, 12))
        points_by_key.setdefault(
            key,
            {
                "m": int(row["m"]),
                "n": int(row["n"]),
                "a_wrap": float(a_wrap),
                "b_wrap": float(b_wrap),
            },
        )
    return sorted(
        points_by_key.values(),
        key=lambda item: (
            float(item["b_wrap"]),
            float(item["a_wrap"]),
            int(item["m"]),
            int(item["n"]),
        ),
    )


def _canonical_orbit_key(a_wrap: float, b_wrap: float) -> tuple[float, float]:
    candidates: list[tuple[float, float]] = []
    for a_value in (float(a_wrap), (1.0 - float(a_wrap)) % 1.0):
        for b_value in (float(b_wrap), (1.0 - float(b_wrap)) % 1.0):
            candidates.append((round(a_value, 12), round(b_value, 12)))
            candidates.append((round(b_value, 12), round(a_value, 12)))
    return min(candidates)


def _compute_orbit_reduced_score(
    summary_rows: list[dict[str, Any]],
    *,
    include_origin: bool = False,
    z_key: str = "target_z",
) -> dict[str, Any]:
    filtered_rows: list[dict[str, Any]] = []
    for row in summary_rows:
        a_value = float(row["a_wrap"])
        b_value = float(row["b_wrap"])
        if not include_origin and abs(a_value) < 1.0e-12 and abs(b_value) < 1.0e-12:
            continue
        filtered_rows.append(row)

    if not filtered_rows:
        return {
            "n_points": 0,
            "n_orbits": 0,
            "rms_z": float("nan"),
            "mean_orbit_rms_z": float("nan"),
            "mean_orbit_abs_z": float("nan"),
            "max_abs_z": float("nan"),
            "orbits": [],
        }

    orbit_groups: dict[tuple[float, float], list[dict[str, Any]]] = {}
    for row in filtered_rows:
        orbit_key = _canonical_orbit_key(float(row["a_wrap"]), float(row["b_wrap"]))
        orbit_groups.setdefault(orbit_key, []).append(row)

    orbit_rows: list[dict[str, Any]] = []
    orbit_rms_values: list[float] = []
    orbit_abs_values: list[float] = []
    orbit_max_values: list[float] = []
    for orbit_key, group in sorted(orbit_groups.items()):
        z_values = np.asarray(
            [float(item[z_key]) for item in group if np.isfinite(float(item[z_key]))],
            dtype=float,
        )
        if z_values.size == 0:
            continue
        rms_z = float(np.sqrt(np.mean(np.square(z_values))))
        mean_abs_z = float(np.mean(np.abs(z_values)))
        max_abs_z = float(np.max(np.abs(z_values)))
        orbit_rms_values.append(rms_z)
        orbit_abs_values.append(mean_abs_z)
        orbit_max_values.append(max_abs_z)
        orbit_rows.append(
            {
                "a_canonical": float(orbit_key[0]),
                "b_canonical": float(orbit_key[1]),
                "multiplicity": int(len(group)),
                "rms_z": rms_z,
                "mean_abs_z": mean_abs_z,
                "max_abs_z": max_abs_z,
            }
        )

    orbit_rms_array = np.asarray(orbit_rms_values, dtype=float)
    orbit_abs_array = np.asarray(orbit_abs_values, dtype=float)
    orbit_max_array = np.asarray(orbit_max_values, dtype=float)
    return {
        "n_points": int(len(filtered_rows)),
        "n_orbits": int(len(orbit_rows)),
        "rms_z": float(np.sqrt(np.mean(np.square(orbit_rms_array)))) if orbit_rms_array.size else float("nan"),
        "mean_orbit_rms_z": float(np.mean(orbit_rms_array)) if orbit_rms_array.size else float("nan"),
        "mean_orbit_abs_z": float(np.mean(orbit_abs_array)) if orbit_abs_array.size else float("nan"),
        "max_abs_z": float(np.max(orbit_max_array)) if orbit_max_array.size else float("nan"),
        "orbits": orbit_rows,
    }


def _compute_summary_score(
    summary_rows: list[dict[str, Any]],
    *,
    z_key: str = "target_z",
    abs_delta_key: str = "target_abs_delta",
) -> dict[str, float]:
    z_values = np.asarray(
        [float(row[z_key]) for row in summary_rows if np.isfinite(float(row[z_key]))],
        dtype=float,
    )
    abs_delta_values = np.asarray(
        [float(row[abs_delta_key]) for row in summary_rows if np.isfinite(float(row[abs_delta_key]))],
        dtype=float,
    )
    if z_values.size == 0:
        return {
            "n_points": 0.0,
            "chi2": float("nan"),
            "rms_z": float("nan"),
            "mean_abs_z": float("nan"),
            "max_abs_z": float("nan"),
            "mean_abs_delta": float("nan"),
        }
    abs_z_values = np.abs(z_values)
    return {
        "n_points": float(z_values.size),
        "chi2": float(np.sum(np.square(z_values))),
        "rms_z": float(np.sqrt(np.mean(np.square(z_values)))),
        "mean_abs_z": float(np.mean(abs_z_values)),
        "max_abs_z": float(np.max(abs_z_values)),
        "mean_abs_delta": float(np.mean(abs_delta_values)) if abs_delta_values.size > 0 else float("nan"),
    }


def _compute_leave_one_out_jackknife(
    rows: list[dict[str, Any]],
    *,
    statistic_builder: Any,
    value_key: str,
) -> dict[str, float]:
    n_rows = int(len(rows))
    if n_rows <= 1:
        return {
            "n_jackknife": float(n_rows),
            "jackknife_mean": float("nan"),
            "jackknife_sigma": float("nan"),
        }

    leave_one_out_values: list[float] = []
    for index in range(n_rows):
        reduced_rows = rows[:index] + rows[index + 1 :]
        statistic_payload = statistic_builder(reduced_rows)
        value = float(statistic_payload[value_key])
        if np.isfinite(value):
            leave_one_out_values.append(value)

    if not leave_one_out_values:
        return {
            "n_jackknife": float(n_rows),
            "jackknife_mean": float("nan"),
            "jackknife_sigma": float("nan"),
        }

    jackknife_array = np.asarray(leave_one_out_values, dtype=float)
    jackknife_mean = float(np.mean(jackknife_array))
    jackknife_sigma = float(
        math.sqrt(
            max(
                0.0,
                (float(jackknife_array.size) - 1.0)
                / float(jackknife_array.size)
                * float(np.sum(np.square(jackknife_array - jackknife_mean))),
            )
        )
    )
    return {
        "n_jackknife": float(jackknife_array.size),
        "jackknife_mean": jackknife_mean,
        "jackknife_sigma": jackknife_sigma,
    }


def _compute_summary_score_jackknife(
    summary_rows: list[dict[str, Any]],
    *,
    z_key: str = "target_z",
    abs_delta_key: str = "target_abs_delta",
) -> dict[str, float]:
    score = _compute_summary_score(summary_rows, z_key=z_key, abs_delta_key=abs_delta_key)
    jackknife = _compute_leave_one_out_jackknife(
        summary_rows,
        statistic_builder=lambda rows: _compute_summary_score(rows, z_key=z_key, abs_delta_key=abs_delta_key),
        value_key="rms_z",
    )
    return {
        **score,
        "jackknife_mean": float(jackknife["jackknife_mean"]),
        "jackknife_sigma": float(jackknife["jackknife_sigma"]),
        "n_jackknife": float(jackknife["n_jackknife"]),
    }


def _compute_orbit_reduced_score_jackknife(
    summary_rows: list[dict[str, Any]],
    *,
    include_origin: bool = False,
    z_key: str = "target_z",
) -> dict[str, Any]:
    score = _compute_orbit_reduced_score(summary_rows, include_origin=include_origin, z_key=z_key)
    jackknife = _compute_leave_one_out_jackknife(
        summary_rows,
        statistic_builder=lambda rows: _compute_orbit_reduced_score(rows, include_origin=include_origin, z_key=z_key),
        value_key="rms_z",
    )
    return {
        **score,
        "jackknife_mean": float(jackknife["jackknife_mean"]),
        "jackknife_sigma": float(jackknife["jackknife_sigma"]),
        "n_jackknife": float(jackknife["n_jackknife"]),
    }


def _combine_sigmas(*sigma_values: float) -> float:
    variance = 0.0
    has_finite = False
    for sigma_value in sigma_values:
        sigma_float = float(sigma_value)
        if not np.isfinite(sigma_float) or sigma_float <= 0.0:
            continue
        variance += sigma_float ** 2
        has_finite = True
    return float(math.sqrt(variance)) if has_finite else float("nan")


def _fit_blind_power_model_with_cov(
    x_values: np.ndarray,
    y_values: np.ndarray,
    sigma_values: np.ndarray,
) -> dict[str, Any]:
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)
    sigma_array = np.asarray(sigma_values, dtype=float)
    sigma_fit = np.where(np.isfinite(sigma_array) & (sigma_array > 0.0), sigma_array, np.nan)
    finite_sigma = sigma_fit[np.isfinite(sigma_fit)]
    sigma_floor = float(np.median(finite_sigma)) if finite_sigma.size else 1.0
    sigma_fit = np.where(np.isfinite(sigma_fit) & (sigma_fit > 0.0), sigma_fit, sigma_floor)

    y_span = float(np.max(y_array) - np.min(y_array)) if y_array.size else 0.0
    A0 = float(y_array[0] - 0.25 * max(y_span, 1.0e-3))
    B0 = float(max(y_array[-1] - A0, 1.0e-4))
    omega0 = 1.0

    try:
        popt, pcov = curve_fit(
            _evaluate_power_model_on_x,
            x_array,
            y_array,
            p0=[A0, B0, omega0],
            sigma=sigma_fit,
            absolute_sigma=True,
            bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 6.0]),
            maxfev=20000,
        )
        A_value = float(popt[0])
        B_value = float(popt[1])
        omega_value = float(popt[2])
        sigma_A = (
            float(np.sqrt(max(float(pcov[0, 0]), 0.0)))
            if np.isfinite(pcov[0, 0]) else float("nan")
        )
        fit_mode = "power_fit"
    except Exception:
        A_value = float(y_array[0])
        sigma_A = float(sigma_fit[0]) if sigma_fit.size else float("nan")
        B_value = 0.0
        omega_value = 1.0
        fit_mode = "fit_failed"
        pcov = np.full((3, 3), np.nan, dtype=float)

    return {
        "A": float(A_value),
        "sigma_A": float(sigma_A),
        "B": float(B_value),
        "omega": float(omega_value),
        "n_used": int(len(x_array)),
        "fit_mode": str(fit_mode),
        "pcov": np.asarray(pcov, dtype=float),
    }


def _individual_fit_target_prediction(
    fit_payload: dict[str, Any],
    *,
    target_x: float,
    target_value: float,
    target_sigma: float,
) -> dict[str, float]:
    A_value = float(fit_payload["A"])
    B_value = float(fit_payload["B"])
    omega_value = float(fit_payload["omega"])
    pred_value = float(
        _evaluate_power_model_on_x(
            np.asarray([target_x], dtype=float),
            A_value,
            B_value,
            omega_value,
        )[0]
    )

    pred_sigma = float("nan")
    pcov = np.asarray(fit_payload.get("pcov"), dtype=float)
    if pcov.shape == (3, 3) and np.all(np.isfinite(pcov)):
        x_pow = float(target_x ** omega_value)
        grad = np.asarray(
            [
                1.0,
                x_pow,
                B_value * x_pow * math.log(target_x),
            ],
            dtype=float,
        )
        pred_var = float(grad @ pcov @ grad)
        if pred_var >= 0.0:
            pred_sigma = float(math.sqrt(pred_var))

    delta = pred_value - float(target_value)
    abs_delta = abs(delta)
    denom = float("nan")
    if np.isfinite(pred_sigma) and np.isfinite(float(target_sigma)):
        denom = float(math.sqrt(pred_sigma ** 2 + float(target_sigma) ** 2))
    elif np.isfinite(float(target_sigma)):
        denom = float(target_sigma)
    z_value = delta / denom if np.isfinite(denom) and denom > 0.0 else float("nan")

    return {
        "pred_value": float(pred_value),
        "pred_sigma": float(pred_sigma),
        "delta": float(delta),
        "abs_delta": float(abs_delta),
        "z_value": float(z_value),
    }


def _largest_size_holdout_diagnostics(
    x_values: np.ndarray,
    y_values: np.ndarray,
    sigma_values: np.ndarray,
    size_labels: list[int],
) -> dict[str, float]:
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)
    sigma_array = np.asarray(sigma_values, dtype=float)
    if x_array.size < 4:
        return {
            "holdout_size": float("nan"),
            "holdout_pred_value": float("nan"),
            "holdout_pred_sigma": float("nan"),
            "holdout_delta": float("nan"),
            "holdout_abs_delta": float("nan"),
            "holdout_z": float("nan"),
        }
    holdout_fit = _fit_blind_power_model_with_cov(x_array[1:], y_array[1:], sigma_array[1:])
    holdout_prediction = _individual_fit_target_prediction(
        holdout_fit,
        target_x=float(x_array[0]),
        target_value=float(y_array[0]),
        target_sigma=float(sigma_array[0]),
    )
    return {
        "holdout_size": float(int(size_labels[0])),
        "holdout_pred_value": float(holdout_prediction["pred_value"]),
        "holdout_pred_sigma": float(holdout_prediction["pred_sigma"]),
        "holdout_delta": float(holdout_prediction["delta"]),
        "holdout_abs_delta": float(holdout_prediction["abs_delta"]),
        "holdout_z": float(holdout_prediction["z_value"]),
    }


def _write_table(output_path: str, summary_rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "row_index",
                    "m",
                    "n",
                    "a_wrap",
                    "b_wrap",
                    "twisted_value",
                    "twisted_sigma",
                    "fit_A",
                    "fit_sigma_A",
                    "fit_B",
                    "fit_omega",
                    "fit_mode",
                    "pred_target",
                    "pred_target_sigma",
                    "target_abs_delta",
                    "target_z",
                    "target_calibrated_sigma",
                    "target_calibrated_z",
                    "holdout_size",
                    "holdout_pred_value",
                    "holdout_pred_sigma",
                    "holdout_abs_delta",
                    "holdout_z",
                ]
            )
            + "\n"
        )
        for index, row in enumerate(summary_rows, start=1):
            handle.write(
                "\t".join(
                    [
                        str(index),
                        str(int(row["m"])),
                        str(int(row["n"])),
                        f"{float(row['a_wrap']):.8f}",
                        f"{float(row['b_wrap']):.8f}",
                        f"{float(row['twisted_value']):.10f}",
                        f"{float(row['twisted_sigma']):.10f}",
                        f"{float(row['fit_A']):.10f}",
                        f"{float(row['fit_sigma_A']):.10f}",
                        f"{float(row['fit_B']):.10f}",
                        f"{float(row['fit_omega']):.10f}",
                        str(row["fit_mode"]),
                        f"{float(row['pred_target']):.10f}",
                        f"{float(row['pred_target_sigma']):.10f}",
                        f"{float(row['target_abs_delta']):.10f}",
                        f"{float(row['target_z']):.10f}",
                        f"{float(row['target_calibrated_sigma']):.10f}",
                        f"{float(row['target_calibrated_z']):.10f}",
                        str(int(row['holdout_size'])) if np.isfinite(float(row['holdout_size'])) else "-1",
                        f"{float(row['holdout_pred_value']):.10f}",
                        f"{float(row['holdout_pred_sigma']):.10f}",
                        f"{float(row['holdout_abs_delta']):.10f}",
                        f"{float(row['holdout_z']):.10f}",
                    ]
                )
                + "\n"
            )


def _build_summary_rows(
    *,
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    untwisted_embedding_cycles: tuple[int, int],
    twisted_embedding_cycles: tuple[int, int],
    include_origin: bool,
    normalization_mode: str = "raw",
    anchor_m: int = 0,
    anchor_n: int = -1,
) -> list[dict[str, Any]]:
    normalization_mode = str(normalization_mode).strip().lower()
    if normalization_mode not in {"raw", "anchor_ratio", "l8_ratio"}:
        raise ValueError(f"unknown normalization mode: {normalization_mode}")

    smallest_lattice = UNTWISTED_LATTICES[0]
    smallest_rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, smallest_lattice))
    anchor_point: dict[str, Any] | None = None
    anchor_key: tuple[float, float] | None = None
    if normalization_mode == "anchor_ratio":
        anchor_point = _find_point_by_mn(
            smallest_rows,
            smallest_lattice,
            embedding_cycles=untwisted_embedding_cycles,
            target_m=int(anchor_m),
            target_n=int(anchor_n),
        )
        anchor_key = (round(float(anchor_point["a_wrap"]), 12), round(float(anchor_point["b_wrap"]), 12))

    base_points = _collect_all_base_points(
        smallest_rows,
        smallest_lattice,
        embedding_cycles=untwisted_embedding_cycles,
        include_origin=include_origin,
    )
    if anchor_point is not None:
        base_points = [
            point
            for point in base_points
            if not (int(point["m"]) == int(anchor_point["m"]) and int(point["n"]) == int(anchor_point["n"]))
        ]
        if not base_points:
            raise ValueError("no non-anchor points remain for normalization-free scoring")

    untwisted_maps: dict[tuple[int, int, int, int], dict[tuple[float, float], dict[str, Any]]] = {}
    for lattice in UNTWISTED_LATTICES:
        rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, lattice))
        untwisted_maps[lattice] = _aggregate_by_fraction(rows, lattice, embedding_cycles=untwisted_embedding_cycles)

    twisted_rows = _load_dat_rows(twisted_dat)
    twisted_interpolator = _build_periodic_interpolator(
        twisted_rows,
        twisted_lattice,
        embedding_cycles=twisted_embedding_cycles,
    )
    target_x = 1.0 / _sqrt_volume(twisted_lattice)
    twisted_anchor_value = float("nan")
    twisted_anchor_sigma = float("nan")
    if anchor_point is not None:
        twisted_anchor_value = _evaluate_periodic(
            twisted_interpolator["value_interpolator"],
            twisted_interpolator["nearest_value"],
            float(anchor_point["a_wrap"]),
            float(anchor_point["b_wrap"]),
            a_core=twisted_interpolator["a_core"],
            b_core=twisted_interpolator["b_core"],
            z_core=twisted_interpolator["value_core"],
        )
        twisted_anchor_sigma = _evaluate_periodic(
            twisted_interpolator["sigma_interpolator"],
            twisted_interpolator["nearest_sigma"],
            float(anchor_point["a_wrap"]),
            float(anchor_point["b_wrap"]),
            a_core=twisted_interpolator["a_core"],
            b_core=twisted_interpolator["b_core"],
            z_core=twisted_interpolator["sigma_core"],
        )

    summary_rows: list[dict[str, Any]] = []
    for point in base_points:
        key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
        baseline_payload = untwisted_maps[smallest_lattice].get(key)
        if baseline_payload is None:
            raise KeyError(f"missing fractional coordinate {key} for untwisted lattice {smallest_lattice}")
        x_values: list[float] = []
        y_values: list[float] = []
        y_errors: list[float] = []
        for lattice in UNTWISTED_LATTICES:
            payload = untwisted_maps[lattice].get(key)
            if payload is None:
                raise KeyError(f"missing fractional coordinate {key} for untwisted lattice {lattice}")
            x_values.append(1.0 / _sqrt_volume(lattice))
            if normalization_mode == "raw":
                y_values.append(float(payload["value"]))
                y_errors.append(float(payload["sigma"]))
            elif normalization_mode == "anchor_ratio":
                if anchor_key is None:
                    raise RuntimeError("anchor key missing for anchor_ratio mode")
                anchor_payload = untwisted_maps[lattice].get(anchor_key)
                if anchor_payload is None:
                    raise KeyError(f"missing anchor coordinate {anchor_key} for untwisted lattice {lattice}")
                ratio, ratio_sigma = _ratio_with_uncertainty(
                    float(payload["value"]),
                    float(payload["sigma"]),
                    float(anchor_payload["value"]),
                    float(anchor_payload["sigma"]),
                )
                y_values.append(ratio)
                y_errors.append(ratio_sigma)
            elif normalization_mode == "l8_ratio":
                if lattice == smallest_lattice:
                    y_values.append(1.0)
                    y_errors.append(0.0)
                else:
                    ratio, ratio_sigma = _ratio_with_uncertainty(
                        float(payload["value"]),
                        float(payload["sigma"]),
                        float(baseline_payload["value"]),
                        float(baseline_payload["sigma"]),
                    )
                    y_values.append(ratio)
                    y_errors.append(ratio_sigma)

        x_array = np.asarray(x_values, dtype=float)
        y_array = np.asarray(y_values, dtype=float)
        yerr_array = np.asarray(y_errors, dtype=float)
        size_array = np.asarray([int(lattice[0]) for lattice in UNTWISTED_LATTICES], dtype=int)
        order = np.argsort(x_array)
        ordered_size_labels = [int(size_array[idx]) for idx in order]
        fit_payload = _fit_blind_power_model_with_cov(
            x_array[order],
            y_array[order],
            yerr_array[order],
        )
        holdout_payload = _largest_size_holdout_diagnostics(
            x_array[order],
            y_array[order],
            yerr_array[order],
            ordered_size_labels,
        )

        twisted_value = _evaluate_periodic(
            twisted_interpolator["value_interpolator"],
            twisted_interpolator["nearest_value"],
            float(point["a_wrap"]),
            float(point["b_wrap"]),
            a_core=twisted_interpolator["a_core"],
            b_core=twisted_interpolator["b_core"],
            z_core=twisted_interpolator["value_core"],
        )
        twisted_sigma = _evaluate_periodic(
            twisted_interpolator["sigma_interpolator"],
            twisted_interpolator["nearest_sigma"],
            float(point["a_wrap"]),
            float(point["b_wrap"]),
            a_core=twisted_interpolator["a_core"],
            b_core=twisted_interpolator["b_core"],
            z_core=twisted_interpolator["sigma_core"],
        )
        if normalization_mode == "anchor_ratio":
            twisted_value, twisted_sigma = _ratio_with_uncertainty(
                float(twisted_value),
                float(twisted_sigma),
                float(twisted_anchor_value),
                float(twisted_anchor_sigma),
            )
        elif normalization_mode == "l8_ratio":
            twisted_value, twisted_sigma = _ratio_with_uncertainty(
                float(twisted_value),
                float(twisted_sigma),
                float(baseline_payload["value"]),
                float(baseline_payload["sigma"]),
            )
        target_prediction = _individual_fit_target_prediction(
            fit_payload,
            target_x=target_x,
            target_value=float(twisted_value),
            target_sigma=float(twisted_sigma),
        )
        target_calibrated_sigma = _combine_sigmas(
            float(target_prediction["pred_sigma"]),
            float(twisted_sigma),
            float(holdout_payload["holdout_abs_delta"]),
        )
        target_calibrated_z = (
            float(target_prediction["delta"]) / float(target_calibrated_sigma)
            if np.isfinite(float(target_calibrated_sigma)) and float(target_calibrated_sigma) > 0.0
            else float("nan")
        )
        summary_rows.append(
            {
                "m": int(point["m"]),
                "n": int(point["n"]),
                "a_wrap": float(point["a_wrap"]),
                "b_wrap": float(point["b_wrap"]),
                "twisted_value": float(twisted_value),
                "twisted_sigma": float(twisted_sigma),
                "fit_A": float(fit_payload["A"]),
                "fit_sigma_A": float(fit_payload["sigma_A"]),
                "fit_B": float(fit_payload["B"]),
                "fit_omega": float(fit_payload["omega"]),
                "fit_mode": str(fit_payload["fit_mode"]),
                "pred_target": float(target_prediction["pred_value"]),
                "pred_target_sigma": float(target_prediction["pred_sigma"]),
                "target_abs_delta": float(target_prediction["abs_delta"]),
                "target_z": float(target_prediction["z_value"]),
                "target_calibrated_sigma": float(target_calibrated_sigma),
                "target_calibrated_z": float(target_calibrated_z),
                "holdout_size": float(holdout_payload["holdout_size"]),
                "holdout_pred_value": float(holdout_payload["holdout_pred_value"]),
                "holdout_pred_sigma": float(holdout_payload["holdout_pred_sigma"]),
                "holdout_abs_delta": float(holdout_payload["holdout_abs_delta"]),
                "holdout_z": float(holdout_payload["holdout_z"]),
            }
        )
    return summary_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit every point on the acute456 8x8 untwisted base lattice independently with "
            "A + B x^omega, evaluate the fit at the large twisted reference volume, and write "
            "a TSV table of fit parameters and target z-scores."
        )
    )
    parser.add_argument("--untwisted-dir", default=DEFAULT_UNTWISTED_DIR)
    parser.add_argument("--twisted-dat", default=DEFAULT_TWISTED_DAT)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--untwisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--twisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--twisted-lattice", nargs=4, type=int, default=list(TWISTED_LATTICE))
    parser.add_argument("--exclude-origin", action="store_true")
    parser.add_argument("--normalization-mode", choices=["raw", "anchor_ratio", "l8_ratio"], default="raw")
    parser.add_argument("--anchor-m", type=int, default=0)
    parser.add_argument("--anchor-n", type=int, default=-1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    untwisted_dir = os.path.abspath(args.untwisted_dir)
    twisted_dat = os.path.abspath(args.twisted_dat)
    output_path = os.path.abspath(args.output)
    summary_rows = _build_summary_rows(
        untwisted_dir=untwisted_dir,
        twisted_dat=twisted_dat,
        twisted_lattice=(
            int(args.twisted_lattice[0]),
            int(args.twisted_lattice[1]),
            int(args.twisted_lattice[2]),
            int(args.twisted_lattice[3]),
        ),
        untwisted_embedding_cycles=(
            int(args.untwisted_embedding_cycles[0]),
            int(args.untwisted_embedding_cycles[1]),
        ),
        twisted_embedding_cycles=(
            int(args.twisted_embedding_cycles[0]),
            int(args.twisted_embedding_cycles[1]),
        ),
        include_origin=not bool(args.exclude_origin),
        normalization_mode=str(args.normalization_mode),
        anchor_m=int(args.anchor_m),
        anchor_n=int(args.anchor_n),
    )
    _write_table(output_path, summary_rows)

    family_label = _infer_family_label(untwisted_dir)
    aggregate_score = _compute_summary_score(summary_rows, z_key="target_z", abs_delta_key="target_abs_delta")
    origin_dropped_rows = [
        row
        for row in summary_rows
        if abs(float(row["a_wrap"])) >= 1.0e-12 or abs(float(row["b_wrap"])) >= 1.0e-12
    ]
    origin_dropped_aggregate_score = _compute_summary_score(
        origin_dropped_rows,
        z_key="target_z",
        abs_delta_key="target_abs_delta",
    )
    calibrated_score = _compute_summary_score(summary_rows, z_key="target_calibrated_z", abs_delta_key="target_abs_delta")
    holdout_score = _compute_summary_score(summary_rows, z_key="holdout_z", abs_delta_key="holdout_abs_delta")
    orbit_score = _compute_orbit_reduced_score(summary_rows, include_origin=False, z_key="target_z")
    orbit_calibrated_score = _compute_orbit_reduced_score(summary_rows, include_origin=False, z_key="target_calibrated_z")
    orbit_holdout_score = _compute_orbit_reduced_score(summary_rows, include_origin=False, z_key="holdout_z")
    ranked_rows = sorted(
        summary_rows,
        key=lambda row: abs(float(row["target_z"])) if np.isfinite(float(row["target_z"])) else -1.0,
        reverse=True,
    )

    print(f"family {family_label}")
    print(f"wrote {output_path}")
    print(
        "aggregate_score "
        f"n={int(round(float(aggregate_score['n_points'])))} "
        f"rms_z={float(aggregate_score['rms_z']):.6f} "
        f"chi2={float(aggregate_score['chi2']):.6f} "
        f"mean_abs_z={float(aggregate_score['mean_abs_z']):.6f} "
        f"max_abs_z={float(aggregate_score['max_abs_z']):.6f} "
        f"mean_abs_delta={float(aggregate_score['mean_abs_delta']):.6f}"
    )
    print(
        "origin_dropped_aggregate_score "
        f"n={int(round(float(origin_dropped_aggregate_score['n_points'])))} "
        f"rms_z={float(origin_dropped_aggregate_score['rms_z']):.6f} "
        f"chi2={float(origin_dropped_aggregate_score['chi2']):.6f} "
        f"mean_abs_z={float(origin_dropped_aggregate_score['mean_abs_z']):.6f} "
        f"max_abs_z={float(origin_dropped_aggregate_score['max_abs_z']):.6f} "
        f"mean_abs_delta={float(origin_dropped_aggregate_score['mean_abs_delta']):.6f}"
    )
    print(
        "orbit_reduced_score "
        f"n_points={int(orbit_score['n_points'])} "
        f"n_orbits={int(orbit_score['n_orbits'])} "
        f"rms_z={float(orbit_score['rms_z']):.6f} "
        f"mean_orbit_rms_z={float(orbit_score['mean_orbit_rms_z']):.6f} "
        f"mean_orbit_abs_z={float(orbit_score['mean_orbit_abs_z']):.6f} "
        f"max_abs_z={float(orbit_score['max_abs_z']):.6f}"
    )
    print(
        "calibrated_target_score "
        f"n={int(round(float(calibrated_score['n_points'])))} "
        f"rms_z={float(calibrated_score['rms_z']):.6f} "
        f"mean_abs_z={float(calibrated_score['mean_abs_z']):.6f} "
        f"max_abs_z={float(calibrated_score['max_abs_z']):.6f} "
        f"orbit_rms_z={float(orbit_calibrated_score['rms_z']):.6f}"
    )
    print(
        "holdout_score "
        f"n={int(round(float(holdout_score['n_points'])))} "
        f"rms_z={float(holdout_score['rms_z']):.6f} "
        f"mean_abs_z={float(holdout_score['mean_abs_z']):.6f} "
        f"max_abs_z={float(holdout_score['max_abs_z']):.6f} "
        f"orbit_rms_z={float(orbit_holdout_score['rms_z']):.6f}"
    )
    for row in ranked_rows[:10]:
        print(
            "worst_point "
            f"m={row['m']} n={row['n']} a={row['a_wrap']:.6f} b={row['b_wrap']:.6f} "
            f"fit_A={row['fit_A']:.6f} fit_sigma_A={row['fit_sigma_A']:.6f} "
            f"fit_B={row['fit_B']:.6f} fit_omega={row['fit_omega']:.6f} "
            f"pred_target={row['pred_target']:.6f} target={row['twisted_value']:.6f} "
            f"abs_delta={row['target_abs_delta']:.6f} z={row['target_z']:.6f} "
            f"cal_z={row['target_calibrated_z']:.6f} holdout_z={row['holdout_z']:.6f}"
        )


if __name__ == "__main__":
    main()