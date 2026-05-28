#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
import math
import os
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_standard_acute456_center_fss import (
    STANDARD_ROOT,
    UNTWISTED_LATTICES,
    _aggregate_by_fraction,
    _connected_ratio_with_true_jackknife,
    _ensure_single_disp_samples,
    _evaluate_power_model_on_x,
    _fit_blind_power_model,
    _infer_family_label,
    _load_single_disp_samples,
    _load_dat_rows,
    _sqrt_volume,
    JACKKNIFE_CACHE_ROOT,
    _load_json,
    _meta_path_for_dat,
    _parse_run_id,
    _single_disp_sample_name,
)


DEFAULT_UNTWISTED_DIR = os.path.join(
    STANDARD_ROOT,
    "data",
    "iso111",
    "untwisted",
    "r1_1p000000__r2_1p000000",
)
DEFAULT_FRACTION = 0.5
DEFAULT_TWISTED_DAT = os.path.join(
    STANDARD_ROOT,
    "data",
    "iso111",
    "twisted",
    "reference",
    "Lx144_Ly144_Tx0_Ty0",
    "two_point_all_to_all.dat",
)
DEFAULT_OUTPUT = os.path.join(STANDARD_ROOT, "iso111_boundary_midpoint_fss.png")
DEFAULT_TABLE_OUTPUT = os.path.join(STANDARD_ROOT, "iso111_boundary_midpoint_fss.tsv")
DEFAULT_ALL_POINTS_OUTPUT_DIR = os.path.join(STANDARD_ROOT, "1-1-1 FSS")
DEFAULT_SWEEP_ROOT = os.path.dirname(DEFAULT_UNTWISTED_DIR)
DEFAULT_TWISTED_LATTICE = (144, 144, 0, 0)
CHECK_POINT_COLOR = "#1d4ed8"

MIDPOINTS: list[dict[str, Any]] = [
    {
        "label": "v midpoint",
        "short_label": "v/2",
        "key": (0.5, 0.0),
        "color": "#1d4ed8",
    },
    {
        "label": "u midpoint",
        "short_label": "u/2",
        "key": (0.0, 0.5),
        "color": "#047857",
    },
    {
        "label": "w midpoint",
        "short_label": "w/2",
        "key": (0.5, 0.5),
        "color": "#b45309",
    },
]

RATIO_SERIES: list[dict[str, str]] = [
    {
        "label": "v/u midpoint ratio",
        "short_label": "(v/2)/(u/2)",
        "numerator": "v/2",
        "denominator": "u/2",
        "color": "#7c2d12",
    },
    {
        "label": "w/u midpoint ratio",
        "short_label": "(w/2)/(u/2)",
        "numerator": "w/2",
        "denominator": "u/2",
        "color": "#7c3aed",
    },
]


def _fraction_short_label(direction: str, fraction: float) -> str:
    frac = Fraction(fraction).limit_denominator(32)
    if frac.numerator == 1 and math.isclose(float(frac), fraction, rel_tol=0.0, abs_tol=1.0e-12):
        return f"{direction}/{frac.denominator}"
    return f"{direction}@{fraction:.3f}"


def _fraction_point_label(fraction: float) -> str:
    if math.isclose(fraction, 0.5, rel_tol=0.0, abs_tol=1.0e-12):
        return "midpoint"
    if math.isclose(fraction, 0.25, rel_tol=0.0, abs_tol=1.0e-12):
        return "quarter-point"
    frac = Fraction(fraction).limit_denominator(32)
    if frac.numerator == 1 and math.isclose(float(frac), fraction, rel_tol=0.0, abs_tol=1.0e-12):
        return f"1/{frac.denominator} point"
    return f"{fraction:.3f}-fraction point"


def _fraction_slug(fraction: float) -> str:
    if math.isclose(fraction, 0.5, rel_tol=0.0, abs_tol=1.0e-12):
        return "midpoint"
    if math.isclose(fraction, 0.25, rel_tol=0.0, abs_tol=1.0e-12):
        return "quarter"
    frac = Fraction(fraction).limit_denominator(32)
    if frac.numerator == 1 and math.isclose(float(frac), fraction, rel_tol=0.0, abs_tol=1.0e-12):
        return f"one_over_{frac.denominator}"
    return f"frac_{fraction:.3f}".replace(".", "p")


def _default_output_path(fraction: float) -> str:
    return os.path.join(STANDARD_ROOT, f"iso111_boundary_{_fraction_slug(fraction)}_fss.png")


def _default_table_output_path(fraction: float) -> str:
    return os.path.join(STANDARD_ROOT, f"iso111_boundary_{_fraction_slug(fraction)}_fss.tsv")


def _boundary_point_specs(fraction: float) -> list[dict[str, Any]]:
    if math.isclose(fraction, DEFAULT_FRACTION, rel_tol=0.0, abs_tol=1.0e-12):
        return MIDPOINTS
    point_label = _fraction_point_label(fraction)
    return [
        {
            "label": f"v {point_label}",
            "short_label": _fraction_short_label("v", fraction),
            "key": (fraction, 0.0),
            "color": "#1d4ed8",
        },
        {
            "label": f"u {point_label}",
            "short_label": _fraction_short_label("u", fraction),
            "key": (0.0, fraction),
            "color": "#047857",
        },
        {
            "label": f"w {point_label}",
            "short_label": _fraction_short_label("w", fraction),
            "key": (fraction, fraction),
            "color": "#b45309",
        },
    ]


def _boundary_ratio_series(fraction: float) -> list[dict[str, str]]:
    if math.isclose(fraction, DEFAULT_FRACTION, rel_tol=0.0, abs_tol=1.0e-12):
        return RATIO_SERIES
    point_label = _fraction_point_label(fraction)
    v_label = _fraction_short_label("v", fraction)
    u_label = _fraction_short_label("u", fraction)
    w_label = _fraction_short_label("w", fraction)
    return [
        {
            "label": f"v/u {point_label} ratio",
            "short_label": f"({v_label})/({u_label})",
            "numerator": v_label,
            "denominator": u_label,
            "color": "#7c2d12",
        },
        {
            "label": f"w/u {point_label} ratio",
            "short_label": f"({w_label})/({u_label})",
            "numerator": w_label,
            "denominator": u_label,
            "color": "#7c3aed",
        },
    ]


def _untwisted_dat_path(root: str, lattice: tuple[int, int, int, int]) -> str:
    lx, ly, tx, ty = lattice
    return os.path.join(root, f"Lx{lx}_Ly{ly}_Tx{tx}_Ty{ty}", "two_point_all_to_all.dat")


def _require_point(
    aggregated: dict[tuple[float, float], dict[str, Any]],
    *,
    key: tuple[float, float],
    dat_path: str,
    nearest_tolerance: float | None = None,
) -> dict[str, Any]:
    payload = aggregated.get((round(float(key[0]), 12), round(float(key[1]), 12)))
    if payload is None:
        if nearest_tolerance is None:
            raise KeyError(f"missing boundary point key {key} in {dat_path}")
        best_key, best_payload = min(
            aggregated.items(),
            key=lambda item: (float(item[0][0]) - float(key[0])) ** 2 + (float(item[0][1]) - float(key[1])) ** 2,
        )
        best_distance = math.sqrt(
            (float(best_key[0]) - float(key[0])) ** 2 + (float(best_key[1]) - float(key[1])) ** 2
        )
        if best_distance > float(nearest_tolerance):
            raise KeyError(f"missing boundary point key {key} in {dat_path}")
        payload = best_payload
    return payload


def _connected_value_with_true_jackknife(
    dat_path: str,
    *,
    point_m: int,
    point_n: int,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]],
    value_cache: dict[tuple[str, int, int], tuple[float, float]],
    allow_regeneration: bool = True,
    fallback_value: tuple[float, float] | None = None,
) -> tuple[float, float]:
    cache_key = (dat_path, int(point_m), int(point_n))
    cached = value_cache.get(cache_key)
    if cached is not None:
        return cached

    sample_payload = sample_cache.get(cache_key)
    if sample_payload is None:
        sample_path = _cached_single_disp_sample_path(dat_path, point_m=int(point_m), point_n=int(point_n))
        if os.path.exists(sample_path):
            sample_payload = _load_single_disp_samples(sample_path)
        elif allow_regeneration:
            sample_payload = _load_single_disp_samples(
                _ensure_single_disp_samples(dat_path, m_value=int(point_m), n_value=int(point_n))
            )
        elif fallback_value is not None:
            value_cache[cache_key] = (float(fallback_value[0]), float(fallback_value[1]))
            return value_cache[cache_key]
        else:
            raise FileNotFoundError(
                f"no cached single-displacement samples for {dat_path} @ ({point_m},{point_n})"
            )
        sample_cache[cache_key] = sample_payload

    corr_samples = np.asarray(sample_payload["corr"], dtype=float)
    mag_samples = np.asarray(sample_payload["mag"], dtype=float)
    n_samples = int(corr_samples.size)
    if n_samples <= 1:
        raise ValueError(f"need at least two samples for midpoint jackknife: {dat_path} @ ({point_m},{point_n})")

    mean_mag = float(np.mean(mag_samples))
    conn_value = float(np.mean(corr_samples) - mean_mag * mean_mag)

    sum_corr = float(np.sum(corr_samples))
    sum_mag = float(np.sum(mag_samples))
    leave_one_out: list[float] = []
    for idx in range(n_samples):
        mean_mag_leave = (sum_mag - float(mag_samples[idx])) / float(n_samples - 1)
        conn_leave = (sum_corr - float(corr_samples[idx])) / float(n_samples - 1) - mean_mag_leave * mean_mag_leave
        leave_one_out.append(float(conn_leave))

    jack_values = np.asarray(leave_one_out, dtype=float)
    jack_var = float((float(n_samples) - 1.0) / float(n_samples) * np.sum(np.square(jack_values - conn_value)))
    conn_sigma = float(math.sqrt(max(jack_var, 0.0)))
    result = (conn_value, conn_sigma)
    value_cache[cache_key] = result
    return result


def _cached_single_disp_sample_path(
    dat_path: str,
    *,
    point_m: int,
    point_n: int,
) -> str:
    meta = _load_json(_meta_path_for_dat(dat_path))
    source_all_to_all = str(meta["source_all_to_all_file"])
    source_run_dir = os.path.dirname(source_all_to_all)
    run_info = _parse_run_id(source_run_dir)
    cache_root = os.path.join(JACKKNIFE_CACHE_ROOT, str(meta["label"]))
    return os.path.join(cache_root, str(run_info["run_id"]), _single_disp_sample_name(point_m, point_n))


def _ratio_from_summary_values(
    *,
    numerator_value: float,
    numerator_sigma: float,
    denominator_value: float,
    denominator_sigma: float,
) -> tuple[float, float]:
    if denominator_value == 0.0:
        raise ZeroDivisionError("anchor connected correlator vanished in summary ratio estimate")
    ratio_value = float(numerator_value) / float(denominator_value)
    sigma_sq = 0.0
    if np.isfinite(float(numerator_sigma)):
        sigma_sq += (float(numerator_sigma) / float(denominator_value)) ** 2
    if np.isfinite(float(denominator_sigma)):
        sigma_sq += (
            float(numerator_value) * float(denominator_sigma) / (float(denominator_value) ** 2)
        ) ** 2
    return float(ratio_value), float(math.sqrt(max(sigma_sq, 0.0)))


def _connected_ratio(
    dat_path: str,
    *,
    point_m: int,
    point_n: int,
    anchor_m: int,
    anchor_n: int,
    point_value: float,
    point_sigma: float,
    anchor_value: float,
    anchor_sigma: float,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]],
    value_cache: dict[tuple[str, int, int], tuple[float, float]],
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]],
    allow_regeneration: bool,
) -> tuple[float, float]:
    cache_key = (dat_path, int(point_m), int(point_n), int(anchor_m), int(anchor_n))
    cached = ratio_cache.get(cache_key)
    if cached is not None:
        return cached

    if allow_regeneration:
        result = _connected_ratio_with_true_jackknife(
            dat_path=dat_path,
            point_m=int(point_m),
            point_n=int(point_n),
            anchor_m=int(anchor_m),
            anchor_n=int(anchor_n),
            sample_cache=sample_cache,
            ratio_cache=ratio_cache,
        )
        ratio_cache[cache_key] = result
        return result

    point_conn, point_conn_sigma = _connected_value_with_true_jackknife(
        dat_path,
        point_m=int(point_m),
        point_n=int(point_n),
        sample_cache=sample_cache,
        value_cache=value_cache,
        allow_regeneration=False,
        fallback_value=(float(point_value), float(point_sigma)),
    )
    anchor_conn, anchor_conn_sigma = _connected_value_with_true_jackknife(
        dat_path,
        point_m=int(anchor_m),
        point_n=int(anchor_n),
        sample_cache=sample_cache,
        value_cache=value_cache,
        allow_regeneration=False,
        fallback_value=(float(anchor_value), float(anchor_sigma)),
    )
    result = _ratio_from_summary_values(
        numerator_value=float(point_conn),
        numerator_sigma=float(point_conn_sigma),
        denominator_value=float(anchor_conn),
        denominator_sigma=float(anchor_conn_sigma),
    )
    ratio_cache[cache_key] = result
    return result


def _build_series(
    *,
    point_specs: list[dict[str, Any]],
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    nearest_tolerance: float | None = None,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] | None = None,
    value_cache: dict[tuple[str, int, int], tuple[float, float]] | None = None,
    allow_regeneration: bool = True,
) -> list[dict[str, Any]]:
    return _build_point_payloads(
        point_specs=point_specs,
        untwisted_dir=untwisted_dir,
        twisted_dat=twisted_dat,
        twisted_lattice=twisted_lattice,
        nearest_tolerance=nearest_tolerance,
        sample_cache=sample_cache,
        value_cache=value_cache,
        allow_regeneration=allow_regeneration,
    )


def _build_point_payloads(
    *,
    point_specs: list[dict[str, Any]],
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    nearest_tolerance: float | None = None,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] | None = None,
    value_cache: dict[tuple[str, int, int], tuple[float, float]] | None = None,
    allow_regeneration: bool = True,
) -> list[dict[str, Any]]:
    sample_cache = {} if sample_cache is None else sample_cache
    value_cache = {} if value_cache is None else value_cache
    untwisted_maps: dict[tuple[int, int, int, int], dict[tuple[float, float], dict[str, Any]]] = {}
    for lattice in UNTWISTED_LATTICES:
        dat_path = _untwisted_dat_path(untwisted_dir, lattice)
        untwisted_maps[lattice] = _aggregate_by_fraction(
            _load_dat_rows(dat_path),
            lattice,
            embedding_cycles=(0, 1),
        )

    twisted_rows = _load_dat_rows(twisted_dat)
    twisted_map = _aggregate_by_fraction(twisted_rows, twisted_lattice, embedding_cycles=(0, 1))

    point_payloads: list[dict[str, Any]] = []
    for point_spec in point_specs:
        x_values: list[float] = []
        y_values: list[float] = []
        sigma_values: list[float] = []
        lattice_labels: list[str] = []
        lattice_points: list[dict[str, Any]] = []
        representative_payload: dict[str, Any] | None = None
        for lattice in UNTWISTED_LATTICES:
            dat_path = _untwisted_dat_path(untwisted_dir, lattice)
            payload = _require_point(untwisted_maps[lattice], key=point_spec["key"], dat_path=dat_path)
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
            x_values.append(1.0 / _sqrt_volume(lattice))
            y_values.append(float(conn_value))
            sigma_values.append(float(conn_sigma))
            lattice_labels.append(f"L={int(lattice[0])}")
            lattice_points.append(
                {
                    "dat_path": dat_path,
                    "lattice": lattice,
                    "m": int(payload["m"]),
                    "n": int(payload["n"]),
                }
            )

        twisted_payload = _require_point(
            twisted_map,
            key=point_spec["key"],
            dat_path=twisted_dat,
            nearest_tolerance=nearest_tolerance,
        )
        twisted_value, twisted_sigma = _connected_value_with_true_jackknife(
            twisted_dat,
            point_m=int(twisted_payload["m"]),
            point_n=int(twisted_payload["n"]),
            sample_cache=sample_cache,
            value_cache=value_cache,
            allow_regeneration=allow_regeneration,
            fallback_value=(float(twisted_payload["value"]), float(twisted_payload["sigma"])),
        )
        fit_payload = _fit_blind_power_model(
            np.asarray(x_values, dtype=float),
            np.asarray(y_values, dtype=float),
            np.asarray(sigma_values, dtype=float),
        )
        if representative_payload is None:
            raise ValueError(f"no representative payload found for point {point_spec['key']}")
        point_payloads.append(
            {
                **point_spec,
                "x": np.asarray(x_values, dtype=float),
                "y": np.asarray(y_values, dtype=float),
                "sigma": np.asarray(sigma_values, dtype=float),
                "lattice_labels": lattice_labels,
                "lattice_points": lattice_points,
                "m": int(point_spec.get("m", representative_payload["m"])),
                "n": int(point_spec.get("n", representative_payload["n"])),
                "a_wrap": float(point_spec.get("a_wrap", representative_payload["a_wrap"])),
                "b_wrap": float(point_spec.get("b_wrap", representative_payload["b_wrap"])),
                "twisted_x": 1.0 / _sqrt_volume(twisted_lattice),
                "twisted_y": float(twisted_value),
                "twisted_sigma": float(twisted_sigma),
                "twisted_dat_path": twisted_dat,
                "twisted_m": int(twisted_payload["m"]),
                "twisted_n": int(twisted_payload["n"]),
                "y_axis_label": str(point_spec.get("y_axis_label", "Connected correlator")),
                "fit": fit_payload,
            }
        )
    return point_payloads


def _select_check_point_specs(
    aggregated: dict[tuple[float, float], dict[str, Any]],
) -> list[dict[str, Any]]:
    point_specs: list[dict[str, Any]] = []
    for key, payload in sorted(
        aggregated.items(),
        key=lambda item: (float(item[0][1]), float(item[0][0])),
    ):
        a_wrap = float(key[0])
        b_wrap = float(key[1])
        if a_wrap > 0.5 + 1.0e-12 or b_wrap > 0.5 + 1.0e-12:
            continue
        point_specs.append(
            {
                "label": f"(a,b)=({a_wrap:.3f},{b_wrap:.3f})",
                "short_label": f"pt({a_wrap:.3f},{b_wrap:.3f})",
                "key": (a_wrap, b_wrap),
                "color": CHECK_POINT_COLOR,
                "m": int(payload["m"]),
                "n": int(payload["n"]),
                "a_wrap": a_wrap,
                "b_wrap": b_wrap,
                "y_axis_label": "Connected correlator",
            }
        )
    return point_specs


def _build_ratio_payloads(
    midpoint_payloads: list[dict[str, Any]],
    *,
    ratio_series: list[dict[str, str]],
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] | None = None,
    value_cache: dict[tuple[str, int, int], tuple[float, float]] | None = None,
    allow_regeneration: bool = True,
) -> list[dict[str, Any]]:
    midpoint_by_label = {str(payload["short_label"]): payload for payload in midpoint_payloads}
    ratio_payloads: list[dict[str, Any]] = []
    sample_cache = {} if sample_cache is None else sample_cache
    value_cache = {} if value_cache is None else value_cache
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]] = {}
    for ratio_cfg in ratio_series:
        numerator = midpoint_by_label[str(ratio_cfg["numerator"])]
        denominator = midpoint_by_label[str(ratio_cfg["denominator"])]
        x_values = np.asarray(numerator["x"], dtype=float)
        if x_values.shape != np.asarray(denominator["x"], dtype=float).shape or not np.allclose(x_values, denominator["x"]):
            raise ValueError(f"mismatched x grids for ratio panel {ratio_cfg['short_label']}")

        y_values: list[float] = []
        sigma_values: list[float] = []
        for numerator_point, denominator_point in zip(
            numerator["lattice_points"],
            denominator["lattice_points"],
        ):
            numerator_dat_path = str(numerator_point["dat_path"])
            denominator_dat_path = str(denominator_point["dat_path"])
            if numerator_dat_path != denominator_dat_path:
                raise ValueError(
                    f"mismatched data paths for ratio panel {ratio_cfg['short_label']}: {numerator_dat_path} vs {denominator_dat_path}"
                )
            ratio_value, ratio_sigma = _connected_ratio(
                dat_path=numerator_dat_path,
                point_m=int(numerator_point["m"]),
                point_n=int(numerator_point["n"]),
                anchor_m=int(denominator_point["m"]),
                anchor_n=int(denominator_point["n"]),
                point_value=float(numerator["y"][len(y_values)]),
                point_sigma=float(numerator["sigma"][len(y_values)]),
                anchor_value=float(denominator["y"][len(y_values)]),
                anchor_sigma=float(denominator["sigma"][len(y_values)]),
                sample_cache=sample_cache,
                value_cache=value_cache,
                ratio_cache=ratio_cache,
                allow_regeneration=allow_regeneration,
            )
            y_values.append(float(ratio_value))
            sigma_values.append(float(ratio_sigma))

        numerator_twisted_dat = str(numerator["twisted_dat_path"])
        denominator_twisted_dat = str(denominator["twisted_dat_path"])
        if numerator_twisted_dat != denominator_twisted_dat:
            raise ValueError(
                f"mismatched twisted data paths for ratio panel {ratio_cfg['short_label']}: {numerator_twisted_dat} vs {denominator_twisted_dat}"
            )
        twisted_ratio, twisted_ratio_sigma = _connected_ratio(
            dat_path=numerator_twisted_dat,
            point_m=int(numerator["twisted_m"]),
            point_n=int(numerator["twisted_n"]),
            anchor_m=int(denominator["twisted_m"]),
            anchor_n=int(denominator["twisted_n"]),
            point_value=float(numerator["twisted_y"]),
            point_sigma=float(numerator["twisted_sigma"]),
            anchor_value=float(denominator["twisted_y"]),
            anchor_sigma=float(denominator["twisted_sigma"]),
            sample_cache=sample_cache,
            value_cache=value_cache,
            ratio_cache=ratio_cache,
            allow_regeneration=allow_regeneration,
        )
        fit_payload = _fit_blind_power_model(
            x_values,
            np.asarray(y_values, dtype=float),
            np.asarray(sigma_values, dtype=float),
        )
        ratio_payloads.append(
            {
                **ratio_cfg,
                "x": x_values,
                "y": np.asarray(y_values, dtype=float),
                "sigma": np.asarray(sigma_values, dtype=float),
                "twisted_x": float(numerator["twisted_x"]),
                "twisted_y": float(twisted_ratio),
                "twisted_sigma": float(twisted_ratio_sigma),
                "y_axis_label": "Correlator ratio",
                "fit": fit_payload,
            }
        )
    return ratio_payloads


def _predict_at_target_x(
    fit_payload: dict[str, Any],
    *,
    target_x: float,
) -> tuple[float, float]:
    pred_value = float(
        _evaluate_power_model_on_x(
            np.asarray([target_x], dtype=float),
            float(fit_payload["A"]),
            float(fit_payload["B"]),
            float(fit_payload["omega"]),
        )[0]
    )

    pred_sigma = float("nan")
    pcov = np.asarray(fit_payload.get("pcov"), dtype=float)
    if pcov.shape == (3, 3) and np.all(np.isfinite(pcov)):
        omega_value = float(fit_payload["omega"])
        x_pow = float(target_x ** omega_value)
        grad = np.asarray(
            [
                1.0,
                x_pow,
                float(fit_payload["B"]) * x_pow * math.log(target_x),
            ],
            dtype=float,
        )
        pred_var = float(grad @ pcov @ grad)
        if pred_var >= 0.0:
            pred_sigma = float(math.sqrt(pred_var))
    return pred_value, pred_sigma


def _predict_band(
    fit_payload: dict[str, Any],
    *,
    x_values: np.ndarray,
) -> np.ndarray:
    x_array = np.asarray(x_values, dtype=float)
    band_sigma = np.full_like(x_array, np.nan, dtype=float)
    pcov = np.asarray(fit_payload.get("pcov"), dtype=float)
    if pcov.shape != (3, 3) or not np.all(np.isfinite(pcov)):
        return band_sigma

    omega_value = float(fit_payload["omega"])
    B_value = float(fit_payload["B"])
    x_pow = np.power(x_array, omega_value)
    grad = np.column_stack(
        [
            np.ones_like(x_array, dtype=float),
            x_pow,
            B_value * x_pow * np.log(x_array),
        ]
    )
    pred_var = np.einsum("ij,jk,ik->i", grad, pcov, grad)
    valid = np.isfinite(pred_var) & (pred_var >= 0.0)
    band_sigma[valid] = np.sqrt(pred_var[valid])
    return band_sigma


def _prediction_delta_summary(
    fit_payload: dict[str, Any],
    *,
    target_x: float,
    target_y: float,
    target_sigma: float,
) -> dict[str, float]:
    pred_value, pred_sigma = _predict_at_target_x(fit_payload, target_x=target_x)
    delta = float(pred_value) - float(target_y)
    if np.isfinite(pred_sigma) and pred_sigma > 0.0 and np.isfinite(float(target_sigma)):
        denom = float(math.sqrt(pred_sigma ** 2 + float(target_sigma) ** 2))
    elif np.isfinite(float(target_sigma)) and float(target_sigma) > 0.0:
        denom = float(target_sigma)
    else:
        denom = float("nan")
    z_value = delta / denom if np.isfinite(denom) and denom > 0.0 else float("nan")
    return {
        "pred_target": float(pred_value),
        "pred_target_sigma": float(pred_sigma),
        "delta_target": float(delta),
        "z_target": float(z_value),
    }


def _write_table(output_path: str, payloads: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(
            "panel\tseries\tx\ty\tsigma\tfit_A\tfit_sigma_A\tfit_omega\tpred_target\tpred_target_sigma\tdelta_target\tz_target\ttwisted_x\ttwisted_y\ttwisted_sigma\n"
        )
        for payload in payloads:
            fit = payload["fit"]
            prediction_summary = _prediction_delta_summary(
                fit,
                target_x=float(payload["twisted_x"]),
                target_y=float(payload["twisted_y"]),
                target_sigma=float(payload["twisted_sigma"]),
            )
            for x_value, y_value, sigma_value in zip(payload["x"], payload["y"], payload["sigma"]):
                handle.write(
                    f"{payload['short_label']}\tuntwisted\t{float(x_value):.10f}\t{float(y_value):.10f}\t"
                    f"{float(sigma_value):.10f}\t{float(fit['A']):.10f}\t{float(fit['sigma_A']):.10f}\t"
                    f"{float(fit['omega']):.10f}\t{float(prediction_summary['pred_target']):.10f}\t"
                    f"{float(prediction_summary['pred_target_sigma']):.10f}\t{float(prediction_summary['delta_target']):.10f}\t"
                    f"{float(prediction_summary['z_target']):.10f}\t{float(payload['twisted_x']):.10f}\t"
                    f"{float(payload['twisted_y']):.10f}\t{float(payload['twisted_sigma']):.10f}\n"
                )
            handle.write(
                f"{payload['short_label']}\ttwisted\t{float(payload['twisted_x']):.10f}\t{float(payload['twisted_y']):.10f}\t"
                f"{float(payload['twisted_sigma']):.10f}\t{float(fit['A']):.10f}\t{float(fit['sigma_A']):.10f}\t"
                f"{float(fit['omega']):.10f}\t{float(prediction_summary['pred_target']):.10f}\t"
                f"{float(prediction_summary['pred_target_sigma']):.10f}\t{float(prediction_summary['delta_target']):.10f}\t"
                f"{float(prediction_summary['z_target']):.10f}\t{float(payload['twisted_x']):.10f}\t"
                f"{float(payload['twisted_y']):.10f}\t{float(payload['twisted_sigma']):.10f}\n"
            )


def _parse_coupling_from_dir(untwisted_dir: str) -> tuple[float, float]:
    family_name = os.path.basename(os.path.normpath(untwisted_dir))
    if not (family_name.startswith("r1_") and "__r2_" in family_name):
        raise ValueError(f"could not parse coupling tag from {untwisted_dir}")
    r1_token, r2_token = family_name.split("__r2_", 1)
    return (
        float(r1_token[len("r1_"):].replace("p", ".")),
        float(r2_token.replace("p", ".")),
    )


def _prediction_summary_for_target(
    payload: dict[str, Any],
    *,
    target_prefix: str,
) -> dict[str, float]:
    return _prediction_delta_summary(
        payload["fit"],
        target_x=float(payload[f"{target_prefix}_x"]),
        target_y=float(payload[f"{target_prefix}_y"]),
        target_sigma=float(payload[f"{target_prefix}_sigma"]),
    )


def _sweep_summary_row(
    *,
    untwisted_dir: str,
    payloads: list[dict[str, Any]],
    target_prefix: str = "twisted",
) -> dict[str, Any]:
    r1_value, r2_value = _parse_coupling_from_dir(untwisted_dir)
    row: dict[str, Any] = {
        "coupling_tag": os.path.basename(os.path.normpath(untwisted_dir)),
        "r1": float(r1_value),
        "r2": float(r2_value),
    }
    corr_z2 = 0.0
    ratio_z2 = 0.0
    for payload in payloads:
        prediction_summary = _prediction_summary_for_target(payload, target_prefix=target_prefix)
        z_value = float(prediction_summary["z_target"])
        chi2_value = float(z_value ** 2)
        row[f"{str(payload['short_label'])}_z"] = z_value
        row[f"{str(payload['short_label'])}_chi2"] = chi2_value
        if str(payload.get("y_axis_label", "Connected correlator")) == "Correlator ratio":
            ratio_z2 += chi2_value
        else:
            corr_z2 += chi2_value
    row["corr_chi2"] = float(corr_z2)
    row["ratio_chi2"] = float(ratio_z2)
    row["chi2_sum"] = float(corr_z2 + ratio_z2)
    row["corr_z2"] = float(corr_z2)
    row["ratio_z2"] = float(ratio_z2)
    row["z2_sum"] = float(corr_z2 + ratio_z2)
    return row


def _write_sweep_summary_table(
    *,
    output_path: str,
    rows: list[dict[str, Any]],
    panel_order: list[str],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    chi2_panel_order = [
        f"{header[:-2]}_chi2" if header.endswith("_z") else f"{header}_chi2"
        for header in panel_order
    ]
    headers = [
        "coupling_tag",
        "r1",
        "r2",
        *panel_order,
        *chi2_panel_order,
        "corr_chi2",
        "ratio_chi2",
        "chi2_sum",
        "corr_z2",
        "ratio_z2",
        "z2_sum",
    ]
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            values: list[str] = []
            for header in headers:
                value = row[header]
                if isinstance(value, float):
                    values.append(f"{value:.10f}")
                else:
                    values.append(str(value))
            handle.write("\t".join(values) + "\n")


def _panel_title(payload: dict[str, Any]) -> str:
    if "a_wrap" in payload and "b_wrap" in payload and "m" in payload and "n" in payload:
        return (
            f"check point (a,b)=({float(payload['a_wrap']):.3f},{float(payload['b_wrap']):.3f})  "
            f"(m,n)=({int(payload['m'])},{int(payload['n'])})"
        )
    if "twisted_m" in payload and "twisted_n" in payload:
        return f"{payload['label']}  (m,n)=({int(payload['twisted_m'])},{int(payload['twisted_n'])})"
    return str(payload["label"])


def _render_panel(axis: plt.Axes, *, payload: dict[str, Any]) -> None:
    fit = payload["fit"]
    prediction_summary = _prediction_delta_summary(
        fit,
        target_x=float(payload["twisted_x"]),
        target_y=float(payload["twisted_y"]),
        target_sigma=float(payload["twisted_sigma"]),
    )
    color = str(payload["color"])
    x_min = min(float(np.min(payload["x"])), float(payload["twisted_x"]))
    x_max = max(float(np.max(payload["x"])), float(payload["twisted_x"]))
    secondary_twisted_x = payload.get("secondary_twisted_x")
    if secondary_twisted_x is not None:
        x_min = min(x_min, float(secondary_twisted_x))
        x_max = max(x_max, float(secondary_twisted_x))
    x_plot_min = x_min * 0.9
    x_plot_max = x_max * 1.08
    axis.errorbar(
        payload["x"],
        payload["y"],
        yerr=payload["sigma"],
        fmt="o",
        color=color,
        ecolor=color,
        capsize=3,
        markersize=6,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.1,
        zorder=3,
        label="untwisted",
    )
    axis.errorbar(
        [payload["twisted_x"]],
        [payload["twisted_y"]],
        yerr=[payload["twisted_sigma"]],
        fmt="s",
        color="#7c3aed",
        ecolor="#7c3aed",
        capsize=3,
        markersize=6,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.1,
        zorder=4,
        label="large twisted target",
    )
    if secondary_twisted_x is not None:
        axis.errorbar(
            [float(payload["secondary_twisted_x"])],
            [float(payload["secondary_twisted_y"])],
            yerr=[float(payload["secondary_twisted_sigma"])],
            fmt="D",
            color="#9467bd",
            ecolor="#9467bd",
            capsize=3,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.8,
            linewidth=1.1,
            zorder=4,
            label=str(payload.get("secondary_twisted_label", "smaller twisted target")),
        )
    x_fit = np.geomspace(x_plot_min, x_plot_max, 300)
    y_fit = _evaluate_power_model_on_x(x_fit, float(fit["A"]), float(fit["B"]), float(fit["omega"]))
    y_fit_sigma = _predict_band(fit, x_values=x_fit)
    finite_band = np.isfinite(y_fit_sigma)
    if np.any(finite_band):
        axis.fill_between(
            x_fit[finite_band],
            (y_fit - y_fit_sigma)[finite_band],
            (y_fit + y_fit_sigma)[finite_band],
            color=color,
            alpha=0.16,
            zorder=1,
        )
    axis.plot(x_fit, y_fit, color=color, linewidth=2.0, zorder=2)
    axis.set_title(_panel_title(payload))
    axis.text(
        0.03,
        0.97,
        (
            f"A = {float(fit['A']):.6f} +/- {float(fit['sigma_A']):.6f}\n"
            f"omega = {float(fit['omega']):.3f}\n"
            f"pred@target = {float(prediction_summary['pred_target']):.6f} +/- {float(prediction_summary['pred_target_sigma']):.6f}\n"
            f"pred - target = {float(prediction_summary['delta_target']):+.6f}  z = {float(prediction_summary['z_target']):+.2f}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
    )
    axis.set_xlabel("1 / sqrt(lattice volume)")
    axis.set_ylabel(str(payload.get("y_axis_label", "Connected correlator")))
    axis.set_xscale("log")
    axis.set_xlim(x_plot_min, x_plot_max)
    axis.grid(True, which="both", alpha=0.35)
    axis.legend(loc="lower right", fontsize=8)


def _panel_sheet_output_name(index: int, row_payloads: list[dict[str, Any]]) -> str:
    b_tag = f"{float(row_payloads[0]['b_wrap']):.3f}".replace(".", "p")
    return f"iso111_row_{index:02d}_b{b_tag}_fss.png"


def _group_point_payloads(point_payloads: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    grouped_payloads: list[list[dict[str, Any]]] = []
    current_group: list[dict[str, Any]] = []
    current_b_wrap: float | None = None
    for payload in point_payloads:
        payload_b_wrap = float(payload["b_wrap"])
        if current_b_wrap is None or abs(payload_b_wrap - current_b_wrap) <= 1.0e-12:
            current_group.append(payload)
            current_b_wrap = payload_b_wrap
            continue
        grouped_payloads.append(current_group)
        current_group = [payload]
        current_b_wrap = payload_b_wrap
    if current_group:
        grouped_payloads.append(current_group)
    return grouped_payloads


def _write_all_point_panels(
    *,
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    output_dir: str,
) -> list[str]:
    smallest_lattice = UNTWISTED_LATTICES[0]
    smallest_map = _aggregate_by_fraction(
        _load_dat_rows(_untwisted_dat_path(untwisted_dir, smallest_lattice)),
        smallest_lattice,
        embedding_cycles=(0, 1),
    )
    point_specs = _select_check_point_specs(smallest_map)
    point_payloads = _build_point_payloads(
        point_specs=point_specs,
        untwisted_dir=untwisted_dir,
        twisted_dat=twisted_dat,
        twisted_lattice=twisted_lattice,
    )
    grouped_payloads = _group_point_payloads(point_payloads)
    family_label = _infer_family_label(untwisted_dir)
    os.makedirs(output_dir, exist_ok=True)
    for name in os.listdir(output_dir):
        if name.startswith("iso111_") and name.endswith(".png"):
            os.remove(os.path.join(output_dir, name))

    written_paths: list[str] = []
    for index, row_payloads in enumerate(grouped_payloads, start=1):
        row_b_wrap = float(row_payloads[0]["b_wrap"])
        row_n_value = int(row_payloads[0]["n"])
        fig, axes = plt.subplots(2, 3, figsize=(18, 9.8), squeeze=False)
        axes_flat = list(axes.ravel())
        fig.suptitle("Iso111 check-point correlator FSS", fontsize=15, y=0.98)
        fig.text(
            0.5,
            0.945,
            (
                f"{family_label}: row {index}/5 with b={row_b_wrap:.3f}, n={row_n_value}; "
                f"untwisted L=8,16,24,32,48,64; target {twisted_lattice}"
            ),
            ha="center",
            va="top",
            fontsize=9.2,
            color="#444444",
        )
        for axis, payload in zip(axes_flat, row_payloads):
            _render_panel(axis, payload=payload)
        for axis in axes_flat[len(row_payloads):]:
            axis.axis("off")
        fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.92])
        output_path = os.path.join(output_dir, _panel_sheet_output_name(index, row_payloads))
        fig.savefig(output_path, dpi=180)
        plt.close(fig)
        written_paths.append(output_path)
    print(f"wrote {len(written_paths)} panel sheets under {output_dir}")
    return written_paths


def _discover_sweep_dirs(untwisted_root: str) -> list[str]:
    if not os.path.isdir(untwisted_root):
        raise FileNotFoundError(f"untwisted sweep root not found: {untwisted_root}")

    sweep_dirs: list[str] = []
    for name in sorted(os.listdir(untwisted_root)):
        if not name.startswith("r1_"):
            continue
        candidate_dir = os.path.join(untwisted_root, name)
        if not os.path.isdir(candidate_dir):
            continue
        if os.path.isfile(_untwisted_dat_path(candidate_dir, UNTWISTED_LATTICES[0])):
            sweep_dirs.append(candidate_dir)

    if not sweep_dirs:
        raise FileNotFoundError(
            f"no coupling-point directories with untwisted data found under {untwisted_root}"
        )
    return sweep_dirs


def _dataset_slug_from_untwisted_dir(path: str) -> str:
    normalized = os.path.normpath(os.path.abspath(path))
    base_name = os.path.basename(normalized)
    if base_name in {"untwisted", "twisted"}:
        return os.path.basename(os.path.dirname(normalized))
    return os.path.basename(os.path.dirname(os.path.dirname(normalized)))


def _sweep_output_stem(untwisted_dir: str, *, fraction: float) -> str:
    coupling_tag = os.path.basename(os.path.normpath(untwisted_dir)).replace("__", "_")
    dataset_slug = _dataset_slug_from_untwisted_dir(untwisted_dir)
    return f"{dataset_slug}_{coupling_tag}_boundary_{_fraction_slug(fraction)}_fss"


def _write_sweep_plots(
    *,
    fraction: float,
    untwisted_root: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    secondary_twisted_dat: str | None,
    secondary_twisted_lattice: tuple[int, int, int, int] | None,
    secondary_twisted_label: str | None,
    output_dir: str,
    allow_regeneration: bool = False,
) -> list[dict[str, Any]]:
    os.makedirs(output_dir, exist_ok=True)
    written_results: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    secondary_summary_rows: list[dict[str, Any]] = []
    summary_panel_order: list[str] | None = None
    for untwisted_dir in _discover_sweep_dirs(untwisted_root):
        stem = _sweep_output_stem(untwisted_dir, fraction=fraction)
        plot_result = _plot(
            fraction=fraction,
            untwisted_dir=untwisted_dir,
            twisted_dat=twisted_dat,
            twisted_lattice=twisted_lattice,
            secondary_twisted_dat=secondary_twisted_dat,
            secondary_twisted_lattice=secondary_twisted_lattice,
            secondary_twisted_label=secondary_twisted_label,
            output_path=os.path.join(output_dir, f"{stem}.png"),
            table_output_path=os.path.join(output_dir, f"{stem}.tsv"),
            allow_regeneration=allow_regeneration,
        )
        written_results.append(plot_result)
        summary_rows.append(dict(plot_result["summary_row"]))
        secondary_summary_row = plot_result.get("secondary_summary_row")
        if secondary_summary_row is not None:
            secondary_summary_rows.append(dict(secondary_summary_row))
        if summary_panel_order is None:
            summary_panel_order = list(plot_result["summary_panel_order"])

    if summary_panel_order is None:
        raise ValueError("missing sweep summary panel order")
    dataset_slug = _dataset_slug_from_untwisted_dir(untwisted_root)
    summary_output_path = os.path.join(
        output_dir,
        f"{dataset_slug}_boundary_{_fraction_slug(fraction)}_fss_zscores.tsv",
    )
    _write_sweep_summary_table(
        output_path=summary_output_path,
        rows=summary_rows,
        panel_order=summary_panel_order,
    )
    primary_target_output_path = os.path.join(
        output_dir,
        f"{dataset_slug}_boundary_{_fraction_slug(fraction)}_fss_large_target_zscores.tsv",
    )
    _write_sweep_summary_table(
        output_path=primary_target_output_path,
        rows=summary_rows,
        panel_order=summary_panel_order,
    )
    print(f"wrote {primary_target_output_path}")
    if secondary_summary_rows:
        secondary_target_output_path = os.path.join(
            output_dir,
            f"{dataset_slug}_boundary_{_fraction_slug(fraction)}_fss_small_target_zscores.tsv",
        )
        _write_sweep_summary_table(
            output_path=secondary_target_output_path,
            rows=secondary_summary_rows,
            panel_order=summary_panel_order,
        )
        print(f"wrote {secondary_target_output_path}")
    print(f"wrote {summary_output_path}")
    print(f"wrote {len(written_results)} sweep figures under {output_dir}")
    return written_results


def _plot(
    *,
    fraction: float,
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    secondary_twisted_dat: str | None,
    secondary_twisted_lattice: tuple[int, int, int, int] | None,
    secondary_twisted_label: str | None,
    output_path: str,
    table_output_path: str,
    allow_regeneration: bool = True,
) -> dict[str, Any]:
    family_label = _infer_family_label(untwisted_dir)
    point_specs = _boundary_point_specs(fraction)
    ratio_series = _boundary_ratio_series(fraction)
    point_label = _fraction_point_label(fraction)
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] = {}
    value_cache: dict[tuple[str, int, int], tuple[float, float]] = {}
    midpoint_payloads = _build_series(
        point_specs=point_specs,
        untwisted_dir=untwisted_dir,
        twisted_dat=twisted_dat,
        twisted_lattice=twisted_lattice,
        sample_cache=sample_cache,
        value_cache=value_cache,
        allow_regeneration=allow_regeneration,
    )
    ratio_payloads = _build_ratio_payloads(
        midpoint_payloads,
        ratio_series=ratio_series,
        sample_cache=sample_cache,
        value_cache=value_cache,
        allow_regeneration=allow_regeneration,
    )
    all_payloads = [*midpoint_payloads, *ratio_payloads]
    if secondary_twisted_dat is not None and secondary_twisted_lattice is not None:
        secondary_midpoint_payloads = _build_series(
            point_specs=point_specs,
            untwisted_dir=untwisted_dir,
            twisted_dat=secondary_twisted_dat,
            twisted_lattice=secondary_twisted_lattice,
            nearest_tolerance=0.02,
            sample_cache=sample_cache,
            value_cache=value_cache,
            allow_regeneration=allow_regeneration,
        )
        secondary_ratio_payloads = _build_ratio_payloads(
            secondary_midpoint_payloads,
            ratio_series=ratio_series,
            sample_cache=sample_cache,
            value_cache=value_cache,
            allow_regeneration=allow_regeneration,
        )
        secondary_by_label = {
            str(payload["short_label"]): payload
            for payload in [*secondary_midpoint_payloads, *secondary_ratio_payloads]
        }
        for payload in all_payloads:
            secondary_payload = secondary_by_label.get(str(payload["short_label"]))
            if secondary_payload is None:
                continue
            payload["secondary_twisted_x"] = float(secondary_payload["twisted_x"])
            payload["secondary_twisted_y"] = float(secondary_payload["twisted_y"])
            payload["secondary_twisted_sigma"] = float(secondary_payload["twisted_sigma"])
            payload["secondary_twisted_label"] = str(secondary_twisted_label or "smaller twisted target")
    _write_table(table_output_path, all_payloads)

    dataset_slug = _dataset_slug_from_untwisted_dir(untwisted_dir)
    fig, axes = plt.subplots(2, 3, figsize=(18, 9.8), squeeze=False)
    axes_flat = list(axes.ravel())
    fig.suptitle(f"{dataset_slug.capitalize()} {point_label} correlator and ratio FSS", fontsize=15, y=0.98)
    fig.text(
        0.5,
        0.945,
        (
            f"standard {family_label}: untwisted L=8,16,24,32,48,64 with primary twisted target at {twisted_lattice}"
            + (
                f" and secondary target at {secondary_twisted_lattice}; "
                if secondary_twisted_dat is not None and secondary_twisted_lattice is not None
                else "; "
            )
            + f"top row = {point_label} correlators, bottom row = {point_label} ratios; "
            + ", ".join(
                f"{spec['short_label']}=({float(spec['key'][0]):.3f},{float(spec['key'][1]):.3f})"
                for spec in point_specs
            )
        ),
        ha="center",
        va="top",
        fontsize=10,
        color="#444444",
    )

    for axis, payload in zip(axes_flat, all_payloads):
        _render_panel(axis, payload=payload)

    for axis in axes_flat[len(all_payloads):]:
        axis.axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.92])
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    summary_row = _sweep_summary_row(untwisted_dir=untwisted_dir, payloads=all_payloads)
    secondary_summary_row = None
    if any("secondary_twisted_x" in payload for payload in all_payloads):
        secondary_summary_row = _sweep_summary_row(
            untwisted_dir=untwisted_dir,
            payloads=all_payloads,
            target_prefix="secondary_twisted",
        )
    summary_panel_order = [f"{str(payload['short_label'])}_z" for payload in all_payloads]
    print(f"wrote {output_path}")
    print(f"wrote {table_output_path}")
    result = {
        "output_path": output_path,
        "table_output_path": table_output_path,
        "summary_row": summary_row,
        "summary_panel_order": summary_panel_order,
    }
    if secondary_summary_row is not None:
        result["secondary_summary_row"] = secondary_summary_row
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot FSS of the iso111 boundary correlators for the three boundary directions with twisted target markers."
    )
    parser.add_argument(
        "--fraction",
        type=float,
        default=DEFAULT_FRACTION,
        help="Wrapped boundary fraction to plot, e.g. 0.5 for midpoint or 0.25 for quarter-point.",
    )
    parser.add_argument("--untwisted-dir", default=DEFAULT_UNTWISTED_DIR)
    parser.add_argument(
        "--sweep-root",
        default=DEFAULT_SWEEP_ROOT,
        help="Root containing the r1_*__r2_* untwisted sweep directories.",
    )
    parser.add_argument("--twisted-dat", default=DEFAULT_TWISTED_DAT)
    parser.add_argument("--twisted-lattice", nargs=4, type=int, default=list(DEFAULT_TWISTED_LATTICE))
    parser.add_argument("--secondary-twisted-dat", default=None)
    parser.add_argument("--secondary-twisted-lattice", nargs=4, type=int, default=None)
    parser.add_argument("--secondary-twisted-label", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--table-output", default=None)
    parser.add_argument(
        "--all-points-output-dir",
        default=None,
        help="Optional directory for separate PNGs of all 25 points in the 1-1-1 check block.",
    )
    parser.add_argument(
        "--sweep-output-dir",
        default=None,
        help="Optional directory for one boundary-point/ratio FSS PNG + TSV per coupling point in the 1-1-1 parameter sweep.",
    )
    parser.add_argument(
        "--force-jackknife-sweep",
        action="store_true",
        help="When used with --sweep-output-dir, regenerate single-displacement samples as needed so sweep error bars use true jackknife values instead of summary fallbacks.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.all_points_output_dir and args.sweep_output_dir:
        raise SystemExit("--all-points-output-dir and --sweep-output-dir are mutually exclusive")
    if (args.secondary_twisted_dat is None) != (args.secondary_twisted_lattice is None):
        raise SystemExit("--secondary-twisted-dat and --secondary-twisted-lattice must be used together")

    fraction = float(args.fraction)
    untwisted_dir = os.path.abspath(args.untwisted_dir)
    sweep_root = os.path.abspath(args.sweep_root)
    twisted_dat = os.path.abspath(args.twisted_dat)
    twisted_lattice = (
        int(args.twisted_lattice[0]),
        int(args.twisted_lattice[1]),
        int(args.twisted_lattice[2]),
        int(args.twisted_lattice[3]),
    )
    secondary_twisted_dat = os.path.abspath(args.secondary_twisted_dat) if args.secondary_twisted_dat else None
    secondary_twisted_lattice = (
        (
            int(args.secondary_twisted_lattice[0]),
            int(args.secondary_twisted_lattice[1]),
            int(args.secondary_twisted_lattice[2]),
            int(args.secondary_twisted_lattice[3]),
        )
        if args.secondary_twisted_lattice
        else None
    )
    if args.sweep_output_dir:
        _write_sweep_plots(
            fraction=fraction,
            untwisted_root=sweep_root,
            twisted_dat=twisted_dat,
            twisted_lattice=twisted_lattice,
            secondary_twisted_dat=secondary_twisted_dat,
            secondary_twisted_lattice=secondary_twisted_lattice,
            secondary_twisted_label=(str(args.secondary_twisted_label) if args.secondary_twisted_label else None),
            output_dir=os.path.abspath(args.sweep_output_dir),
            allow_regeneration=bool(args.force_jackknife_sweep),
        )
        return

    if args.all_points_output_dir:
        _write_all_point_panels(
            untwisted_dir=untwisted_dir,
            twisted_dat=twisted_dat,
            twisted_lattice=twisted_lattice,
            output_dir=os.path.abspath(args.all_points_output_dir),
        )
        return

    _plot(
        fraction=fraction,
        untwisted_dir=untwisted_dir,
        twisted_dat=twisted_dat,
        twisted_lattice=twisted_lattice,
        secondary_twisted_dat=secondary_twisted_dat,
        secondary_twisted_lattice=secondary_twisted_lattice,
        secondary_twisted_label=(str(args.secondary_twisted_label) if args.secondary_twisted_label else None),
        output_path=os.path.abspath(args.output or _default_output_path(fraction)),
        table_output_path=os.path.abspath(args.table_output or _default_table_output_path(fraction)),
    )


if __name__ == "__main__":
    main()