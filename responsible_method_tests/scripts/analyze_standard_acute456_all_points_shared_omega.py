#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
from typing import Any

import numpy as np

from plot_standard_acute456_center_fss import (
    DEFAULT_UNTWISTED_DIR,
    UNTWISTED_LATTICES,
    _aggregate_by_fraction,
    _fit_blind_power_model,
    _fit_shared_blind_power_model,
    _infer_family_label,
    _load_dat_rows,
    _ratio_with_uncertainty,
    _to_ab,
    _untwisted_dat_path,
    _wrap_unit,
    _sqrt_volume,
    _evaluate_power_model_on_x,
)


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
        key=lambda item: (float(item["a_wrap"]), float(item["b_wrap"]), int(item["m"]), int(item["n"])),
    )


def _build_series_payloads(
    *,
    untwisted_dir: str,
    embedding_cycles: tuple[int, int],
    normalization_mode: str,
    include_smallest: bool,
    anchor_m: int,
    anchor_n: int,
    include_origin: bool,
) -> list[dict[str, Any]]:
    smallest_lattice = UNTWISTED_LATTICES[0]
    smallest_rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, smallest_lattice))
    smallest_map = _aggregate_by_fraction(smallest_rows, smallest_lattice, embedding_cycles=embedding_cycles)
    base_points = _collect_all_base_points(
        smallest_rows,
        smallest_lattice,
        embedding_cycles=embedding_cycles,
        include_origin=include_origin,
    )

    untwisted_maps: dict[tuple[int, int, int, int], dict[tuple[float, float], dict[str, Any]]] = {
        smallest_lattice: smallest_map,
    }
    for lattice in UNTWISTED_LATTICES[1:]:
        rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, lattice))
        untwisted_maps[lattice] = _aggregate_by_fraction(rows, lattice, embedding_cycles=embedding_cycles)

    anchor_key: tuple[float, float] | None = None
    if normalization_mode == "anchor_ratio":
        anchor_raw = None
        for row in smallest_rows:
            if int(row["m"]) == int(anchor_m) and int(row["n"]) == int(anchor_n):
                anchor_raw = row
                break
        if anchor_raw is None:
            raise KeyError(f"anchor point (m,n)=({anchor_m},{anchor_n}) not found")
        a_raw, b_raw = _to_ab(
            int(anchor_raw["m"]),
            int(anchor_raw["n"]),
            smallest_lattice[0],
            smallest_lattice[1],
            smallest_lattice[2],
            smallest_lattice[3],
            embedding_cycles=embedding_cycles,
        )
        anchor_key = (round(_wrap_unit(a_raw), 12), round(_wrap_unit(b_raw), 12))

    series_payloads: list[dict[str, Any]] = []
    lattices_to_use = UNTWISTED_LATTICES if include_smallest else UNTWISTED_LATTICES[1:]
    for point in base_points:
        point_key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
        baseline_payload = smallest_map.get(point_key)
        if baseline_payload is None:
            raise KeyError(f"missing smallest-lattice baseline for point {point_key}")

        x_values: list[float] = []
        y_values: list[float] = []
        sigma_values: list[float] = []
        size_labels: list[int] = []
        for lattice in lattices_to_use:
            payload = untwisted_maps[lattice].get(point_key)
            if payload is None:
                raise KeyError(f"missing point {point_key} on lattice {lattice}")

            x_values.append(1.0 / _sqrt_volume(lattice))
            size_labels.append(int(lattice[0]))
            if normalization_mode == "raw":
                y_values.append(float(payload["value"]))
                sigma_values.append(float(payload["sigma"]))
            elif normalization_mode == "l8_ratio":
                ratio, ratio_sigma = _ratio_with_uncertainty(
                    float(payload["value"]),
                    float(payload["sigma"]),
                    float(baseline_payload["value"]),
                    float(baseline_payload["sigma"]),
                )
                y_values.append(ratio)
                sigma_values.append(ratio_sigma)
            elif normalization_mode == "anchor_ratio":
                if anchor_key is None:
                    raise RuntimeError("anchor key missing for anchor_ratio mode")
                anchor_payload = untwisted_maps[lattice].get(anchor_key)
                if anchor_payload is None:
                    raise KeyError(f"missing anchor point {anchor_key} on lattice {lattice}")
                ratio, ratio_sigma = _ratio_with_uncertainty(
                    float(payload["value"]),
                    float(payload["sigma"]),
                    float(anchor_payload["value"]),
                    float(anchor_payload["sigma"]),
                )
                y_values.append(ratio)
                sigma_values.append(ratio_sigma)
            else:
                raise ValueError(f"unknown normalization mode: {normalization_mode}")

        order = np.argsort(np.asarray(x_values, dtype=float))
        series_payloads.append(
            {
                "point": point,
                "x": np.asarray(x_values, dtype=float)[order],
                "y": np.asarray(y_values, dtype=float)[order],
                "sigma": np.asarray(sigma_values, dtype=float)[order],
                "size_labels": [size_labels[idx] for idx in order],
            }
        )
    return series_payloads


def _compute_fit_quality(series_payloads: list[dict[str, Any]], shared_fit: dict[str, Any]) -> dict[str, Any]:
    point_summaries: list[dict[str, Any]] = []
    all_z: list[float] = []
    total_chi2 = 0.0
    total_points = 0
    for series_index, payload in enumerate(series_payloads):
        fit_payload = dict(shared_fit["series"][series_index])
        omega_value = float(shared_fit["omega"])
        y_fit = _evaluate_power_model_on_x(
            np.asarray(payload["x"], dtype=float),
            float(fit_payload["A"]),
            float(fit_payload["B"]),
            omega_value,
        )
        sigma_values = np.asarray(payload["sigma"], dtype=float)
        sigma_fit = np.where(np.isfinite(sigma_values) & (sigma_values > 0.0), sigma_values, np.nan)
        finite_sigma = sigma_fit[np.isfinite(sigma_fit)]
        sigma_floor = float(np.median(finite_sigma)) if finite_sigma.size else 1.0
        sigma_fit = np.where(np.isfinite(sigma_fit) & (sigma_fit > 0.0), sigma_fit, sigma_floor)

        residuals = np.asarray(payload["y"], dtype=float) - y_fit
        z_values = residuals / sigma_fit
        chi2 = float(np.sum(np.square(z_values)))
        total_chi2 += chi2
        total_points += int(z_values.size)
        all_z.extend([float(value) for value in z_values])

        point_summaries.append(
            {
                "m": int(payload["point"]["m"]),
                "n": int(payload["point"]["n"]),
                "a_wrap": float(payload["point"]["a_wrap"]),
                "b_wrap": float(payload["point"]["b_wrap"]),
                "chi2": chi2,
                "rms_z": float(np.sqrt(np.mean(np.square(z_values)))),
                "max_abs_z": float(np.max(np.abs(z_values))),
                "n_sizes": int(z_values.size),
            }
        )

    n_series = len(series_payloads)
    n_params = 1 + 2 * n_series
    dof = total_points - n_params
    all_z_array = np.asarray(all_z, dtype=float)
    worst_by_max = sorted(point_summaries, key=lambda item: float(item["max_abs_z"]), reverse=True)
    worst_by_rms = sorted(point_summaries, key=lambda item: float(item["rms_z"]), reverse=True)
    return {
        "n_series": int(n_series),
        "n_data": int(total_points),
        "n_params": int(n_params),
        "dof": int(dof),
        "chi2": float(total_chi2),
        "chi2_per_dof": float(total_chi2 / dof) if dof > 0 else float("nan"),
        "global_rms_z": float(np.sqrt(np.mean(np.square(all_z_array)))) if all_z_array.size else float("nan"),
        "global_mean_abs_z": float(np.mean(np.abs(all_z_array))) if all_z_array.size else float("nan"),
        "max_abs_z_any": float(np.max(np.abs(all_z_array))) if all_z_array.size else float("nan"),
        "count_max_abs_z_gt_2": int(sum(float(item["max_abs_z"]) > 2.0 for item in point_summaries)),
        "count_max_abs_z_gt_3": int(sum(float(item["max_abs_z"]) > 3.0 for item in point_summaries)),
        "worst_by_max": worst_by_max,
        "worst_by_rms": worst_by_rms,
    }


def _compute_individual_fit_quality(series_payloads: list[dict[str, Any]]) -> dict[str, Any]:
    point_summaries: list[dict[str, Any]] = []
    all_z: list[float] = []
    omegas: list[float] = []
    total_chi2 = 0.0
    total_points = 0
    for payload in series_payloads:
        fit_payload = _fit_blind_power_model(payload["x"], payload["y"], payload["sigma"])
        y_fit = _evaluate_power_model_on_x(
            np.asarray(payload["x"], dtype=float),
            float(fit_payload["A"]),
            float(fit_payload["B"]),
            float(fit_payload["omega"]),
        )
        sigma_values = np.asarray(payload["sigma"], dtype=float)
        sigma_fit = np.where(np.isfinite(sigma_values) & (sigma_values > 0.0), sigma_values, np.nan)
        finite_sigma = sigma_fit[np.isfinite(sigma_fit)]
        sigma_floor = float(np.median(finite_sigma)) if finite_sigma.size else 1.0
        sigma_fit = np.where(np.isfinite(sigma_fit) & (sigma_fit > 0.0), sigma_fit, sigma_floor)

        residuals = np.asarray(payload["y"], dtype=float) - y_fit
        z_values = residuals / sigma_fit
        chi2 = float(np.sum(np.square(z_values)))
        total_chi2 += chi2
        total_points += int(z_values.size)
        all_z.extend([float(value) for value in z_values])
        omegas.append(float(fit_payload["omega"]))

        point_summaries.append(
            {
                "m": int(payload["point"]["m"]),
                "n": int(payload["point"]["n"]),
                "a_wrap": float(payload["point"]["a_wrap"]),
                "b_wrap": float(payload["point"]["b_wrap"]),
                "chi2": chi2,
                "rms_z": float(np.sqrt(np.mean(np.square(z_values)))),
                "max_abs_z": float(np.max(np.abs(z_values))),
                "fit_omega": float(fit_payload["omega"]),
                "fit_mode": str(fit_payload["fit_mode"]),
            }
        )

    n_series = len(series_payloads)
    n_params = 3 * n_series
    dof = total_points - n_params
    all_z_array = np.asarray(all_z, dtype=float)
    omega_array = np.asarray(omegas, dtype=float)
    return {
        "n_series": int(n_series),
        "n_data": int(total_points),
        "n_params": int(n_params),
        "dof": int(dof),
        "chi2": float(total_chi2),
        "chi2_per_dof": float(total_chi2 / dof) if dof > 0 else float("nan"),
        "global_rms_z": float(np.sqrt(np.mean(np.square(all_z_array)))) if all_z_array.size else float("nan"),
        "omega_mean": float(np.mean(omega_array)) if omega_array.size else float("nan"),
        "omega_std": float(np.std(omega_array)) if omega_array.size else float("nan"),
        "omega_min": float(np.min(omega_array)) if omega_array.size else float("nan"),
        "omega_max": float(np.max(omega_array)) if omega_array.size else float("nan"),
        "count_omega_at_lower_bound": int(sum(abs(float(value) - 0.05) < 1.0e-8 for value in omega_array)),
        "point_summaries": point_summaries,
    }


def _analyze_family(
    *,
    untwisted_dir: str,
    embedding_cycles: tuple[int, int],
    normalization_mode: str,
    include_smallest: bool,
    anchor_m: int,
    anchor_n: int,
    include_origin: bool,
    top_k: int,
) -> None:
    series_payloads = _build_series_payloads(
        untwisted_dir=untwisted_dir,
        embedding_cycles=embedding_cycles,
        normalization_mode=normalization_mode,
        include_smallest=include_smallest,
        anchor_m=anchor_m,
        anchor_n=anchor_n,
        include_origin=include_origin,
    )
    shared_fit = _fit_shared_blind_power_model(series_payloads)
    fit_quality = _compute_fit_quality(series_payloads, shared_fit)
    individual_quality = _compute_individual_fit_quality(series_payloads)

    family_label = _infer_family_label(untwisted_dir)
    print(f"family {family_label}")
    print(
        "shared_omega_summary "
        f"normalization_mode={normalization_mode} "
        f"include_smallest={int(include_smallest)} "
        f"include_origin={int(include_origin)} "
        f"n_series={fit_quality['n_series']} "
        f"n_data={fit_quality['n_data']} "
        f"n_params={fit_quality['n_params']} "
        f"dof={fit_quality['dof']} "
        f"omega={float(shared_fit['omega']):.6f} "
        f"sigma_omega={float(shared_fit['sigma_omega']):.6f} "
        f"fit_mode={shared_fit['fit_mode']}"
    )
    print(
        "fit_quality "
        f"chi2={float(fit_quality['chi2']):.6f} "
        f"chi2_per_dof={float(fit_quality['chi2_per_dof']):.6f} "
        f"global_rms_z={float(fit_quality['global_rms_z']):.6f} "
        f"global_mean_abs_z={float(fit_quality['global_mean_abs_z']):.6f} "
        f"max_abs_z_any={float(fit_quality['max_abs_z_any']):.6f} "
        f"count_max_abs_z_gt_2={int(fit_quality['count_max_abs_z_gt_2'])} "
        f"count_max_abs_z_gt_3={int(fit_quality['count_max_abs_z_gt_3'])}"
    )
    print(
        "individual_omega_quality "
        f"chi2={float(individual_quality['chi2']):.6f} "
        f"chi2_per_dof={float(individual_quality['chi2_per_dof']):.6f} "
        f"global_rms_z={float(individual_quality['global_rms_z']):.6f} "
        f"omega_mean={float(individual_quality['omega_mean']):.6f} "
        f"omega_std={float(individual_quality['omega_std']):.6f} "
        f"omega_min={float(individual_quality['omega_min']):.6f} "
        f"omega_max={float(individual_quality['omega_max']):.6f} "
        f"count_omega_at_lower_bound={int(individual_quality['count_omega_at_lower_bound'])}"
    )
    print(
        "shared_vs_individual "
        f"delta_chi2={float(fit_quality['chi2'] - individual_quality['chi2']):.6f} "
        f"delta_params={int(fit_quality['n_params'] - individual_quality['n_params'])}"
    )
    for item in fit_quality["worst_by_max"][: max(1, int(top_k))]:
        print(
            "worst_point "
            f"m={item['m']} n={item['n']} a={item['a_wrap']:.6f} b={item['b_wrap']:.6f} "
            f"chi2={float(item['chi2']):.6f} rms_z={float(item['rms_z']):.6f} max_abs_z={float(item['max_abs_z']):.6f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify whether all non-origin points on the acute456 8x8 base lattice obey a shared-omega "
            "finite-size scaling ansatz across the standard untwisted size ladder."
        )
    )
    parser.add_argument("--untwisted-dir", nargs="+", default=[DEFAULT_UNTWISTED_DIR])
    parser.add_argument("--untwisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--normalization-mode", choices=["raw", "anchor_ratio", "l8_ratio"], default="l8_ratio")
    parser.add_argument("--include-smallest", action="store_true")
    parser.add_argument("--include-origin", action="store_true")
    parser.add_argument("--anchor-m", type=int, default=0)
    parser.add_argument("--anchor-n", type=int, default=-1)
    parser.add_argument("--top-k", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for idx, untwisted_dir in enumerate(args.untwisted_dir):
        if idx > 0:
            print()
        _analyze_family(
            untwisted_dir=os.path.abspath(untwisted_dir),
            embedding_cycles=(int(args.untwisted_embedding_cycles[0]), int(args.untwisted_embedding_cycles[1])),
            normalization_mode=str(args.normalization_mode),
            include_smallest=bool(args.include_smallest),
            anchor_m=int(args.anchor_m),
            anchor_n=int(args.anchor_n),
            include_origin=bool(args.include_origin),
            top_k=max(1, int(args.top_k)),
        )


if __name__ == "__main__":
    main()