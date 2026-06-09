#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from run_reweight_gradient_flow import GradientFlowRunner
from run_acute456_pow2_blind_holdout import exact_triangular_ising_beta, lattice_tag, run_simulation
from test_geometry_match_grid_interpolation import CouplingPoint, _parse_selected_bundle
from plot_standard_acute456_center_fss import (
    DEFAULT_TWISTED_DAT as ACUTE456_DEFAULT_TWISTED_DAT,
    TWISTED_LATTICE as ACUTE456_DEFAULT_TWISTED_LATTICE,
    UNTWISTED_LATTICES as ACUTE456_DEFAULT_UNTWISTED_LATTICES,
    _aggregate_by_fraction,
    _boundary_paths,
    _build_twisted_target_payload,
    _compute_aggregate_target_score,
    _find_point_by_mn,
    _fit_blind_power_model,
    _fit_shared_blind_power_model,
    _individual_fit_target_prediction,
    _load_dat_rows,
    _ratio_with_uncertainty,
    _select_base_points,
    _shared_fit_target_prediction,
    _sqrt_volume,
    _target_value_for_point,
    _to_ab,
    _untwisted_dat_path,
    _wrap_unit,
)


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent
ACUTE456_TARGET_R1 = 4.702782819756
ACUTE456_TARGET_R2 = 7.353910143333
ACUTE456_SPARSE_DAT_COLUMNS = ("d", "m", "n", "corr", "err", "corr_conn", "err_conn")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a current-stack CMA-ES search on the direct total_score objective and render "
            "how the search distribution learns over generations."
        )
    )
    parser.add_argument("--start-r1", type=float, default=1.0)
    parser.add_argument("--start-r2", type=float, default=1.0)
    parser.add_argument("--step-size", type=float, default=0.025)
    parser.add_argument("--num-steps", type=int, default=1)
    parser.add_argument("--step-mode", default="normalized")
    parser.add_argument("--step-adaptation", default="fixed")
    parser.add_argument("--adaptive-step-snr-reference", type=float, default=3.0)
    parser.add_argument("--adaptive-step-scale-min", type=float, default=0.5)
    parser.add_argument("--adaptive-step-scale-max", type=float, default=3.0)
    parser.add_argument("--min-gradient-snr", type=float, default=0.0)
    parser.add_argument("--gradient-jackknife-blocks", type=int, default=16)
    parser.add_argument("--gradient-step", type=float, default=None)
    parser.add_argument("--accept-rule", default="gradient_only")
    parser.add_argument("--balance-mode", default="none")
    parser.add_argument("--metric-key", default=None)
    parser.add_argument("--target-mode", choices=["iso111", "acute456"], default="iso111")
    parser.add_argument("--target-size", type=int, default=64)
    parser.add_argument("--target-r1", type=float, default=1.0)
    parser.add_argument("--target-r2", type=float, default=1.0)
    parser.add_argument("--fit-sizes", nargs="+", type=int, default=[32, 28, 24, 20, 16, 12, 8, 4])
    parser.add_argument("--acute456-fit-sizes", nargs="+", type=int, default=[lattice[0] for lattice in ACUTE456_DEFAULT_UNTWISTED_LATTICES])
    parser.add_argument("--acute456-target-r1", type=float, default=ACUTE456_TARGET_R1)
    parser.add_argument("--acute456-target-r2", type=float, default=ACUTE456_TARGET_R2)
    parser.add_argument(
        "--acute456-eval-workflow",
        choices=["selected_bundle", "legacy_all_to_all"],
        default="selected_bundle",
        help="How acute456 candidate families are generated before scoring.",
    )
    parser.add_argument("--acute456-fit-mode", choices=["individual", "shared"], default="individual")
    parser.add_argument("--acute456-n-panels", type=int, default=4)
    parser.add_argument("--acute456-normalization-mode", choices=["raw", "anchor_ratio", "l8_ratio"], default="raw")
    parser.add_argument("--acute456-anchor-m", type=int, default=0)
    parser.add_argument("--acute456-anchor-n", type=int, default=-1)
    parser.add_argument("--acute456-untwisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--acute456-twisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--acute456-twisted-dat", default=ACUTE456_DEFAULT_TWISTED_DAT)
    parser.add_argument("--acute456-twisted-lattice", nargs=4, type=int, default=list(ACUTE456_DEFAULT_TWISTED_LATTICE))
    parser.add_argument("--n-traj", type=int, default=100000)
    parser.add_argument("--n-therm", type=int, default=10000)
    parser.add_argument("--n-skip", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=2026060801)
    parser.add_argument("--r1-min", type=float, default=0.1)
    parser.add_argument("--r1-max", type=float, default=10.0)
    parser.add_argument("--r2-min", type=float, default=0.1)
    parser.add_argument("--r2-max", type=float, default=10.0)
    parser.add_argument("--sigma0", type=float, default=0.15)
    parser.add_argument("--popsize", type=int, default=6)
    parser.add_argument("--max-evals", type=int, default=24)
    parser.add_argument("--max-gens", type=int, default=0)
    parser.add_argument("--tolx", type=float, default=5.0e-3)
    parser.add_argument("--tolfun", type=float, default=1.0e-6)
    parser.add_argument("--cma-dsigma", type=float, default=None)
    parser.add_argument("--cma-csigma", type=float, default=None)
    parser.add_argument("--cma-seed", type=int, default=2026060802)
    parser.add_argument(
        "--plot-every-gens",
        type=int,
        default=1,
        help="Write a frame every N generations when --save-frames is enabled.",
    )
    parser.add_argument(
        "--ghost-generations",
        type=int,
        default=10,
        help="How many past Gaussian ellipses to keep visible in the trajectory panel.",
    )
    parser.add_argument("--save-frames", action="store_true")
    parser.add_argument(
        "--resume-root",
        default=None,
        help="Existing CMA-ES output root to resume in place from its saved evals/generations logs.",
    )
    parser.add_argument(
        "--additional-gens",
        type=int,
        default=0,
        help="When resuming, append this many more CMA-ES generations to the saved run.",
    )
    parser.add_argument(
        "--output-root",
        default=str(RESPONSIBLE_ROOT / "reweighting" / "reweight_cmaes_iso111"),
        help="Root directory for CMA-ES outputs and raw Monte Carlo payloads.",
    )
    return parser.parse_args()


def _format_float(value: float, fmt: str = ".10e") -> str:
    return format(float(value), fmt) if math.isfinite(float(value)) else "nan"


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _token_precise(value: float) -> str:
    return f"{float(value):.12f}".replace("-", "m").replace(".", "p")


def _unique_ints_preserve_order(values: list[int]) -> list[int]:
    ordered: list[int] = []
    seen: set[int] = set()
    for raw_value in values:
        value = int(raw_value)
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _acute456_wrap_distance(lhs: float, rhs: float) -> float:
    delta = abs(float(lhs) - float(rhs)) % 1.0
    return min(delta, 1.0 - delta)


def _acute456_point_from_mn(
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    m_value: int,
    n_value: int,
) -> dict[str, Any]:
    a_raw, b_raw = _to_ab(
        int(m_value),
        int(n_value),
        lattice[0],
        lattice[1],
        lattice[2],
        lattice[3],
        embedding_cycles=embedding_cycles,
    )
    return {
        "m": int(m_value),
        "n": int(n_value),
        "a_wrap": _wrap_unit(a_raw),
        "b_wrap": _wrap_unit(b_raw),
    }


def _acute456_point_from_fraction(
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    a_wrap: float,
    b_wrap: float,
) -> dict[str, Any]:
    paths = _boundary_paths(lattice[0], lattice[1], lattice[2], lattice[3])
    (dm_a, dn_a), (dm_b, dn_b) = [paths[idx] for idx in embedding_cycles]
    m_value = int(round(float(a_wrap) * float(dm_a) + float(b_wrap) * float(dm_b)))
    n_value = int(round(float(a_wrap) * float(dn_a) + float(b_wrap) * float(dn_b)))
    point = _acute456_point_from_mn(
        lattice,
        embedding_cycles=embedding_cycles,
        m_value=m_value,
        n_value=n_value,
    )
    if _acute456_wrap_distance(float(point["a_wrap"]), float(a_wrap)) > 1.0e-9:
        raise ValueError(f"acute456 a_wrap mapping failed for lattice {lattice}: {a_wrap} -> {point}")
    if _acute456_wrap_distance(float(point["b_wrap"]), float(b_wrap)) > 1.0e-9:
        raise ValueError(f"acute456 b_wrap mapping failed for lattice {lattice}: {b_wrap} -> {point}")
    return point


def _acute456_select_base_points_for_lattice(
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    n_panels: int,
) -> list[dict[str, Any]]:
    search_radius = max(max(abs(int(value)) for value in lattice), 1)
    candidates: dict[tuple[float, float], dict[str, Any]] = {}
    for n_value in range(-search_radius, search_radius + 1):
        point = _acute456_point_from_mn(
            lattice,
            embedding_cycles=embedding_cycles,
            m_value=0,
            n_value=n_value,
        )
        if 0.0 < float(point["b_wrap"]) <= 0.5 + 1.0e-12:
            key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
            previous = candidates.get(key)
            if previous is None or abs(int(point["n"])) < abs(int(previous["n"])):
                candidates[key] = point
    if not candidates:
        raise ValueError(f"no acute456 representative base points found on lattice {lattice}")
    ordered = sorted(candidates.values(), key=lambda item: (float(item["b_wrap"]), float(item["a_wrap"]), int(item["n"])))
    return ordered[: max(1, min(int(n_panels), len(ordered)))]


def _acute456_jackknife_sigma(leave_one_out: list[float]) -> float:
    if len(leave_one_out) < 2:
        return float("nan")
    jk = np.asarray(leave_one_out, dtype=float)
    return float(np.sqrt((jk.size - 1.0) / jk.size * np.sum((jk - np.mean(jk)) ** 2)))


def _acute456_block_slices(n_samples: int, n_blocks: int = 16) -> list[tuple[int, int]]:
    block_count = max(2, min(int(n_blocks), max(int(n_samples) // 4, 2)))
    block = int(n_samples) // block_count
    slices: list[tuple[int, int]] = []
    for index in range(block_count):
        lo = index * block
        hi = (index + 1) * block if index < block_count - 1 else int(n_samples)
        if hi > lo:
            slices.append((lo, hi))
    return slices


def _acute456_mean_with_sigma(samples: np.ndarray) -> tuple[float, float]:
    sample_array = np.asarray(samples, dtype=float)
    if sample_array.size == 0:
        return float("nan"), float("nan")
    value = float(np.mean(sample_array))
    leave_one_out: list[float] = []
    for lo, hi in _acute456_block_slices(int(sample_array.size)):
        mask = np.ones(int(sample_array.size), dtype=bool)
        mask[lo:hi] = False
        if np.count_nonzero(mask) == 0:
            continue
        leave_one_out.append(float(np.mean(sample_array[mask])))
    return value, _acute456_jackknife_sigma(leave_one_out)


def _acute456_connected_with_sigma(corr_samples: np.ndarray, mag_samples: np.ndarray) -> tuple[float, float]:
    corr_array = np.asarray(corr_samples, dtype=float)
    mag_array = np.asarray(mag_samples, dtype=float)
    if corr_array.size == 0 or mag_array.size == 0 or corr_array.size != mag_array.size:
        return float("nan"), float("nan")
    mean_corr = float(np.mean(corr_array))
    mean_mag = float(np.mean(mag_array))
    connected = float(mean_corr - mean_mag * mean_mag)
    leave_one_out: list[float] = []
    for lo, hi in _acute456_block_slices(int(corr_array.size)):
        mask = np.ones(int(corr_array.size), dtype=bool)
        mask[lo:hi] = False
        if np.count_nonzero(mask) == 0:
            continue
        corr_leave = np.asarray(corr_array[mask], dtype=float)
        mag_leave = np.asarray(mag_array[mask], dtype=float)
        mean_corr_leave = float(np.mean(corr_leave))
        mean_mag_leave = float(np.mean(mag_leave))
        leave_one_out.append(float(mean_corr_leave - mean_mag_leave * mean_mag_leave))
    return connected, _acute456_jackknife_sigma(leave_one_out)


def _acute456_sparse_rows_from_bundle(
    bundle_payload: dict[str, Any],
    lattice_specs: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    mag_samples = np.asarray(bundle_payload["mag"], dtype=float)
    rows: list[dict[str, Any]] = []
    for index, spec in enumerate(lattice_specs, start=1):
        corr_samples = np.asarray(bundle_payload["corr"][str(spec["label"])], dtype=float)
        corr_value, corr_sigma = _acute456_mean_with_sigma(corr_samples)
        corr_conn, corr_conn_sigma = _acute456_connected_with_sigma(corr_samples, mag_samples)
        rows.append(
            {
                "d": int(index),
                "m": int(spec["m"]),
                "n": int(spec["n"]),
                "corr": float(corr_value),
                "err": float(corr_sigma),
                "corr_conn": float(corr_conn),
                "err_conn": float(corr_conn_sigma),
            }
        )
    return rows


def _write_acute456_sparse_dat(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# d m n corr err corr_conn err_conn\n")
        for row in rows:
            handle.write(
                " ".join(
                    [
                        str(int(row["d"])),
                        str(int(row["m"])),
                        str(int(row["n"])),
                        _format_float(float(row["corr"])),
                        _format_float(float(row["err"])),
                        _format_float(float(row["corr_conn"])),
                        _format_float(float(row["err_conn"])),
                    ]
                )
                + "\n"
            )


class Acute456Runner:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_root = Path(args.output_root).resolve()
        self.output_root.mkdir(parents=True, exist_ok=True)
        fit_sizes = _unique_ints_preserve_order(list(args.acute456_fit_sizes))
        if not fit_sizes:
            raise ValueError("acute456_fit_sizes must not be empty")
        self.fit_lattices = tuple((int(size), int(size), 0, 0) for size in fit_sizes)
        self.fit_mode = str(args.acute456_fit_mode)
        self.normalization_mode = str(args.acute456_normalization_mode)
        self.anchor_m = int(args.acute456_anchor_m)
        self.anchor_n = int(args.acute456_anchor_n)
        self.untwisted_embedding_cycles = (
            int(args.acute456_untwisted_embedding_cycles[0]),
            int(args.acute456_untwisted_embedding_cycles[1]),
        )
        self.twisted_embedding_cycles = (
            int(args.acute456_twisted_embedding_cycles[0]),
            int(args.acute456_twisted_embedding_cycles[1]),
        )
        self.twisted_dat = Path(str(args.acute456_twisted_dat)).resolve()
        if not self.twisted_dat.is_file():
            raise FileNotFoundError(f"acute456 twisted target data not found: {self.twisted_dat}")
        self.twisted_lattice = tuple(int(value) for value in args.acute456_twisted_lattice)
        self.eval_workflow = str(args.acute456_eval_workflow)
        self.metric_key = "total_score"
        self.target_point = CouplingPoint(float(args.acute456_target_r1), float(args.acute456_target_r2))
        self.family_root = (
            self.output_root
            / f"acute456_families_traj{int(args.n_traj)}_therm{int(args.n_therm)}_skip{int(args.n_skip)}"
        )
        self.family_root.mkdir(parents=True, exist_ok=True)
        self.anchor_point = _acute456_point_from_mn(
            self.fit_lattices[0],
            embedding_cycles=self.untwisted_embedding_cycles,
            m_value=self.anchor_m,
            n_value=self.anchor_n,
        )
        selected_points = _acute456_select_base_points_for_lattice(
            self.fit_lattices[0],
            embedding_cycles=self.untwisted_embedding_cycles,
            n_panels=int(self.args.acute456_n_panels),
        )
        anchor_key = (round(float(self.anchor_point["a_wrap"]), 12), round(float(self.anchor_point["b_wrap"]), 12))
        self.selected_points: list[dict[str, Any]] = []
        self.bundle_specs: list[dict[str, Any]] = [dict(self.anchor_point, label="anchor")]
        panel_index = 1
        for point in selected_points:
            point_key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
            if point_key == anchor_key:
                self.selected_points.append(dict(self.anchor_point, label="anchor"))
                continue
            label = f"panel_{panel_index:02d}"
            spec = dict(point, label=label)
            self.selected_points.append(spec)
            self.bundle_specs.append(spec)
            panel_index += 1
        self._lattice_bundle_specs: dict[tuple[int, int, int, int], list[dict[str, Any]]] = {}
        self._target_payload: dict[str, Any] | None = None
        self.score_cache: dict[tuple[float, float], dict[str, Any]] = {}
        self.sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] = {}
        self.ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]] = {}

    def _family_key(self, point: CouplingPoint) -> tuple[float, float]:
        return (round(float(point.r1), 12), round(float(point.r2), 12))

    def _family_dir(self, point: CouplingPoint) -> Path:
        return self.family_root / f"r1_{_token_precise(point.r1)}__r2_{_token_precise(point.r2)}"

    def _flat_lattice_dir(self, family_dir: Path, lattice: tuple[int, int, int, int]) -> Path:
        lx, ly, tx, ty = lattice
        return family_dir / f"Lx{lx}_Ly{ly}_Tx{tx}_Ty{ty}"

    def _flatten_dat(self, source_dat: Path, dest_dir: Path) -> None:
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_dat, dest_dir / "two_point_all_to_all.dat")
        source_meta = source_dat.with_suffix(".meta.json")
        if source_meta.is_file():
            shutil.copy2(source_meta, dest_dir / "two_point_all_to_all.meta.json")

    def _bundle_specs_for_lattice(self, lattice: tuple[int, int, int, int]) -> list[dict[str, Any]]:
        cached = self._lattice_bundle_specs.get(lattice)
        if cached is not None:
            return [dict(spec) for spec in cached]
        lattice_specs: list[dict[str, Any]] = []
        for spec in self.bundle_specs:
            lattice_point = _acute456_point_from_fraction(
                lattice,
                embedding_cycles=self.untwisted_embedding_cycles,
                a_wrap=float(spec["a_wrap"]),
                b_wrap=float(spec["b_wrap"]),
            )
            lattice_specs.append(dict(lattice_point, label=str(spec["label"])))
        self._lattice_bundle_specs[lattice] = [dict(spec) for spec in lattice_specs]
        return [dict(spec) for spec in lattice_specs]

    def _ensure_family_data(self, point: CouplingPoint) -> Path:
        family_dir = self._family_dir(point)
        couplings = {"k1": float(point.r1), "k2": float(point.r2), "k3": 1.0}
        beta_c = float(exact_triangular_ising_beta(couplings))
        for lattice in self.fit_lattices:
            flat_dir = self._flat_lattice_dir(family_dir, lattice)
            flat_dat = flat_dir / "two_point_all_to_all.dat"
            if flat_dat.is_file():
                continue
            raw_dir = flat_dir / "_raw"
            label = f"acute456_cmaes_{lattice_tag(lattice)}_{family_dir.name}"
            if self.eval_workflow == "legacy_all_to_all":
                output_file, _ = run_simulation(
                    lattice=lattice,
                    couplings=couplings,
                    beta=beta_c,
                    n_traj=int(self.args.n_traj),
                    n_skip=int(self.args.n_skip),
                    n_therm=int(self.args.n_therm),
                    data_dir=raw_dir,
                    label=label,
                    single_disp=None,
                    force=False,
                )
                self._flatten_dat(Path(output_file), flat_dir)
                continue

            lattice_specs = self._bundle_specs_for_lattice(lattice)
            bundle_path, _ = run_simulation(
                lattice=lattice,
                couplings=couplings,
                beta=beta_c,
                n_traj=int(self.args.n_traj),
                n_skip=int(self.args.n_skip),
                n_therm=int(self.args.n_therm),
                data_dir=raw_dir,
                label=label,
                selected_disp_list=tuple((int(spec["m"]), int(spec["n"])) for spec in lattice_specs),
                selected_disp_bundle_name="selected_observables_bundle.dat",
                force=False,
            )
            bundle_payload = _parse_selected_bundle(Path(bundle_path), [str(spec["label"]) for spec in lattice_specs])
            sparse_rows = _acute456_sparse_rows_from_bundle(bundle_payload, lattice_specs)
            _write_acute456_sparse_dat(flat_dat, sparse_rows)
        return family_dir

    def _ensure_target_payload(self) -> dict[str, Any]:
        if self._target_payload is None:
            self._target_payload = _build_twisted_target_payload(
                dat_path=str(self.twisted_dat),
                lattice=self.twisted_lattice,
                embedding_cycles=self.twisted_embedding_cycles,
                anchor_point=self.anchor_point,
                label="acute456 twisted target",
            )
        return self._target_payload

    def _score_family(self, family_dir: Path) -> dict[str, Any]:
        smallest_lattice = self.fit_lattices[0]
        anchor_point = dict(self.anchor_point)
        selected_points = [dict(point) for point in self.selected_points]
        if self.normalization_mode == "anchor_ratio":
            selected_points = [
                point
                for point in selected_points
                if not (int(point["m"]) == int(anchor_point["m"]) and int(point["n"]) == int(anchor_point["n"]))
            ]
        if not selected_points:
            raise ValueError("acute456 score has no selected points after normalization filtering")

        target_payload = self._ensure_target_payload()
        untwisted_maps: dict[tuple[int, int, int, int], dict[tuple[float, float], dict[str, Any]]] = {}
        for lattice in self.fit_lattices:
            rows = _load_dat_rows(_untwisted_dat_path(str(family_dir), lattice))
            untwisted_maps[lattice] = _aggregate_by_fraction(
                rows,
                lattice,
                embedding_cycles=self.untwisted_embedding_cycles,
            )

        anchor_key = (round(float(anchor_point["a_wrap"]), 12), round(float(anchor_point["b_wrap"]), 12))
        plot_payloads: list[dict[str, Any]] = []
        for point in selected_points:
            point_key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
            baseline_payload = untwisted_maps[smallest_lattice][point_key]

            x_values: list[float] = []
            y_values: list[float] = []
            sigma_values: list[float] = []
            for lattice in self.fit_lattices:
                payload = untwisted_maps[lattice][point_key]
                x_values.append(1.0 / _sqrt_volume(lattice))
                if self.normalization_mode == "raw":
                    y_values.append(float(payload["value"]))
                    sigma_values.append(float(payload["sigma"]))
                elif self.normalization_mode == "l8_ratio":
                    ratio, ratio_sigma = _ratio_with_uncertainty(
                        float(payload["value"]),
                        float(payload["sigma"]),
                        float(baseline_payload["value"]),
                        float(baseline_payload["sigma"]),
                    )
                    y_values.append(float(ratio))
                    sigma_values.append(float(ratio_sigma))
                else:
                    anchor_payload = untwisted_maps[lattice][anchor_key]
                    ratio, ratio_sigma = _ratio_with_uncertainty(
                        float(payload["value"]),
                        float(payload["sigma"]),
                        float(anchor_payload["value"]),
                        float(anchor_payload["sigma"]),
                    )
                    y_values.append(float(ratio))
                    sigma_values.append(float(ratio_sigma))

            target_value, target_sigma = _target_value_for_point(
                target_payload=target_payload,
                point=point,
                normalization_mode=self.normalization_mode,
                baseline_value=float(baseline_payload["value"]),
                baseline_sigma=float(baseline_payload["sigma"]),
                sample_cache=self.sample_cache,
                ratio_cache=self.ratio_cache,
            )
            plot_payloads.append(
                {
                    "point": point,
                    "x": np.asarray(x_values, dtype=float),
                    "y": np.asarray(y_values, dtype=float),
                    "sigma": np.asarray(sigma_values, dtype=float),
                    "target_value": float(target_value),
                    "target_sigma": float(target_sigma),
                }
            )

        summary_rows: list[dict[str, Any]] = []
        if self.fit_mode == "shared":
            shared_fit = _fit_shared_blind_power_model(plot_payloads)
            for series_index, payload in enumerate(plot_payloads):
                target_prediction = _shared_fit_target_prediction(
                    shared_fit,
                    series_index=series_index,
                    target_x=float(target_payload["target_x"]),
                    target_value=float(payload["target_value"]),
                    target_sigma=float(payload["target_sigma"]),
                )
                summary_rows.append(
                    {
                        "target_z": float(target_prediction["z_value"]),
                        "target_abs_delta": float(target_prediction["abs_delta"]),
                    }
                )
        else:
            for payload in plot_payloads:
                fit_payload = _fit_blind_power_model(payload["x"], payload["y"], payload["sigma"])
                target_prediction = _individual_fit_target_prediction(
                    fit_payload,
                    target_x=float(target_payload["target_x"]),
                    target_value=float(payload["target_value"]),
                    target_sigma=float(payload["target_sigma"]),
                )
                summary_rows.append(
                    {
                        "target_z": float(target_prediction["z_value"]),
                        "target_abs_delta": float(target_prediction["abs_delta"]),
                    }
                )

        aggregate = _compute_aggregate_target_score(summary_rows)
        return {
            "total_score": float(aggregate["rms_z"]),
            "acute456_chi2": float(aggregate["chi2"]),
            "acute456_rms_z": float(aggregate["rms_z"]),
            "acute456_mean_abs_z": float(aggregate["mean_abs_z"]),
            "acute456_max_abs_z": float(aggregate["max_abs_z"]),
            "acute456_mean_abs_delta": float(aggregate["mean_abs_delta"]),
            "acute456_n_points": float(aggregate["n_points"]),
        }

    def score_row(self, donor_point: CouplingPoint, eval_point: CouplingPoint, *, step_index: int) -> dict[str, Any]:
        del donor_point, step_index
        key = self._family_key(eval_point)
        cached = self.score_cache.get(key)
        if cached is not None:
            return dict(cached)
        family_dir = self._ensure_family_data(eval_point)
        row = self._score_family(family_dir)
        self.score_cache[key] = dict(row)
        return dict(row)


def _clip_point(point: CouplingPoint, args: argparse.Namespace) -> CouplingPoint:
    return CouplingPoint(
        r1=min(max(float(point.r1), float(args.r1_min)), float(args.r1_max)),
        r2=min(max(float(point.r2), float(args.r2_min)), float(args.r2_max)),
    )


def _point_to_json(point: CouplingPoint) -> dict[str, float]:
    return {"r1": float(point.r1), "r2": float(point.r2)}


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _normalized_covariance_from_row(row: dict[str, Any]) -> np.ndarray:
    return np.asarray(
        [
            [float(row["norm_cov_11"]), float(row["norm_cov_12"])],
            [float(row["norm_cov_12"]), float(row["norm_cov_22"])],
        ],
        dtype=float,
    )


def _decompose_covariance(covariance: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        eigenvalues, basis = np.linalg.eigh(0.5 * (covariance + covariance.T))
    except np.linalg.LinAlgError:
        basis = np.eye(covariance.shape[0], dtype=float)
        eigenvalues = np.ones(covariance.shape[0], dtype=float)
    eigenvalues = np.clip(np.asarray(eigenvalues, dtype=float), 1.0e-30, None)
    diag = np.sqrt(eigenvalues)
    bd = basis * diag
    return basis, diag, bd


def _resume_args_from_summary(cli_args: argparse.Namespace) -> tuple[argparse.Namespace, dict[str, Any]]:
    resume_root = Path(str(cli_args.resume_root)).resolve()
    summary_path = resume_root / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Cannot resume without summary.json in {resume_root}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary_args = dict(summary.get("args", {}))
    if not summary_args:
        raise ValueError(f"summary.json in {resume_root} does not contain saved args")
    if summary_args.get("target_mode") == "acute456" and "acute456_eval_workflow" not in summary_args:
        summary_args["acute456_eval_workflow"] = "legacy_all_to_all"
    summary_args["output_root"] = str(resume_root)
    summary_args["save_frames"] = bool(summary_args.get("save_frames", False) or bool(cli_args.save_frames))
    return argparse.Namespace(**summary_args), summary


def _replay_resume_state(
    *,
    args: argparse.Namespace,
    summary: dict[str, Any],
    params: dict[str, Any],
    eval_rows: list[dict[str, str]],
    generation_rows: list[dict[str, str]],
) -> dict[str, Any]:
    n_dim = 2
    generation_count = int(summary.get("n_generations", max(len(generation_rows) - 1, 0)))
    generation_row_map = {int(row["generation"]): row for row in generation_rows}
    if 0 not in generation_row_map:
        raise ValueError("generations.tsv must contain generation 0 for resume")

    mean = np.asarray(
        [float(generation_row_map[0]["mean_r1"]), float(generation_row_map[0]["mean_r2"])],
        dtype=float,
    )
    sigma = float(generation_row_map[0]["sigma"])
    covariance = _normalized_covariance_from_row(generation_row_map[0])
    p_sigma = np.zeros(n_dim, dtype=float)
    p_c = np.zeros(n_dim, dtype=float)

    rows_by_generation: dict[int, list[dict[str, str]]] = {}
    for row in eval_rows:
        rows_by_generation.setdefault(int(row["generation"]), []).append(row)

    for generation in range(1, generation_count + 1):
        generation_eval_rows = sorted(
            rows_by_generation.get(generation, []),
            key=lambda row: int(row["candidate_index"]),
        )
        if not generation_eval_rows:
            raise ValueError(f"Missing eval rows for saved generation {generation}")

        basis, _, bd = _decompose_covariance(covariance)
        mean_prev = np.asarray(mean, dtype=float)
        sigma_prev = float(sigma)
        covariance_prev = np.asarray(covariance, dtype=float)

        completed_y: list[np.ndarray] = []
        completed_z: list[np.ndarray] = []
        scores: list[float] = []
        for row in generation_eval_rows:
            sampled = np.asarray([float(row["sampled_r1"]), float(row["sampled_r2"])], dtype=float)
            y_value = (sampled - mean_prev) / sigma_prev
            z_value = np.linalg.solve(bd, y_value)
            completed_y.append(y_value)
            completed_z.append(z_value)
            scores.append(float(row["score"]))

        fs = np.asarray(scores, dtype=float)
        completed_y_array = np.asarray(completed_y, dtype=float)
        completed_z_array = np.asarray(completed_z, dtype=float)

        mu_used = max(1, min(int(params["mu"]), len(fs) // 2 + (len(fs) % 2)))
        if mu_used != int(params["mu"]):
            weights = np.asarray(params["weights"][:mu_used], dtype=float)
            weights = weights / np.sum(weights)
            mu_eff = float(1.0 / np.sum(np.square(weights)))
        else:
            weights = np.asarray(params["weights"], dtype=float)
            mu_eff = float(params["mu_eff"])

        logged_selected = [
            local_index for local_index, row in enumerate(generation_eval_rows) if _as_bool(row.get("selected", "False"))
        ]
        if len(logged_selected) == mu_used:
            selected = np.asarray(sorted(logged_selected, key=lambda local_index: float(fs[local_index])), dtype=int)
        else:
            selected = np.argsort(fs)[:mu_used]

        y_selected = completed_y_array[selected]
        z_selected = completed_z_array[selected]
        y_weighted = weights @ y_selected
        z_weighted = weights @ z_selected
        mean = mean_prev + sigma_prev * y_weighted

        c_sigma = float(params["c_sigma"])
        d_sigma = float(params["d_sigma"])
        c_c = float(params["c_c"])
        c_1 = float(params["c_1"])
        c_mu = float(params["c_mu"])
        chi_n = float(params["chi_n"])

        p_sigma = (1.0 - c_sigma) * p_sigma + math.sqrt(c_sigma * (2.0 - c_sigma) * mu_eff) * (basis @ z_weighted)
        sigma = float(sigma_prev * math.exp((c_sigma / d_sigma) * (np.linalg.norm(p_sigma) / chi_n - 1.0)))

        hs_threshold = (1.4 + 2.0 / (n_dim + 1.0)) * chi_n
        denom = math.sqrt(max(1.0 - (1.0 - c_sigma) ** (2 * generation), 1.0e-30))
        h_sigma = 1.0 if (np.linalg.norm(p_sigma) / denom) < hs_threshold else 0.0
        p_c = (1.0 - c_c) * p_c + h_sigma * math.sqrt(c_c * (2.0 - c_c) * mu_eff) * y_weighted

        rank_mu = y_selected.T @ (weights[:, None] * y_selected)
        rank_1 = np.outer(p_c, p_c) + (1.0 - h_sigma) * c_c * (2.0 - c_c) * covariance_prev
        covariance = (1.0 - c_1 - c_mu) * covariance_prev + c_1 * rank_1 + c_mu * rank_mu
        covariance = 0.5 * (covariance + covariance.T)

        saved_row = generation_row_map.get(generation)
        if saved_row is None:
            raise ValueError(f"Missing saved Gaussian row for generation {generation}")
        saved_mean = np.asarray([float(saved_row["mean_r1"]), float(saved_row["mean_r2"])], dtype=float)
        saved_sigma = float(saved_row["sigma"])
        saved_covariance = _normalized_covariance_from_row(saved_row)
        if not np.allclose(mean, saved_mean, rtol=1.0e-9, atol=1.0e-9):
            raise ValueError(f"Resume replay mean mismatch at generation {generation}")
        if not math.isclose(sigma, saved_sigma, rel_tol=1.0e-9, abs_tol=1.0e-9):
            raise ValueError(f"Resume replay sigma mismatch at generation {generation}")
        if not np.allclose(covariance, saved_covariance, rtol=1.0e-9, atol=1.0e-9):
            raise ValueError(f"Resume replay covariance mismatch at generation {generation}")

    basis, diag, bd = _decompose_covariance(covariance)
    final_gaussian = summary.get("final_gaussian", {})
    expected_mean = np.asarray(final_gaussian.get("mean", mean.tolist()), dtype=float)
    expected_sigma = float(final_gaussian.get("sigma", sigma))
    expected_covariance = np.asarray(final_gaussian.get("normalized_covariance", covariance.tolist()), dtype=float)
    if not np.allclose(mean, expected_mean, rtol=1.0e-9, atol=1.0e-9):
        raise ValueError("Resume replay final mean does not match summary.json")
    if not math.isclose(sigma, expected_sigma, rel_tol=1.0e-9, abs_tol=1.0e-9):
        raise ValueError("Resume replay final sigma does not match summary.json")
    if not np.allclose(covariance, expected_covariance, rtol=1.0e-9, atol=1.0e-9):
        raise ValueError("Resume replay final covariance does not match summary.json")

    best_point_dict = dict(summary.get("best_point", {}))
    best_point = CouplingPoint(float(best_point_dict["r1"]), float(best_point_dict["r2"]))
    return {
        "mean": mean,
        "sigma": float(sigma),
        "covariance": covariance,
        "basis": basis,
        "diag": diag,
        "bd": bd,
        "p_sigma": p_sigma,
        "p_c": p_c,
        "generation": generation_count,
        "best_score": float(summary["best_score"]),
        "best_point": best_point,
        "best_eval_id": int(summary["best_eval_id"]),
        "eval_rows": list(eval_rows),
        "generation_rows": list(generation_rows),
    }


def _cma_default_params(
    n_dim: int,
    popsize: int,
    *,
    dsigma: float | None = None,
    csigma: float | None = None,
) -> dict[str, Any]:
    lam = int(popsize)
    mu = lam // 2
    raw_weights = np.array([math.log((lam + 1) / 2.0) - math.log(index + 1) for index in range(mu)], dtype=float)
    weights = raw_weights / np.sum(raw_weights)
    mu_eff = float(1.0 / np.sum(np.square(weights)))

    c_sigma = float((mu_eff + 2.0) / (n_dim + mu_eff + 5.0))
    if csigma is not None:
        c_sigma = float(np.clip(float(csigma), 1.0e-3, 1.0 - 1.0e-3))
    d_sigma = 1.0 + 2.0 * max(0.0, math.sqrt((mu_eff - 1.0) / (n_dim + 1.0)) - 1.0) + c_sigma
    if dsigma is not None:
        d_sigma = float(max(float(dsigma), 0.1))
    c_c = float((4.0 + mu_eff / n_dim) / (n_dim + 4.0 + 2.0 * mu_eff / n_dim))
    c_1 = float(2.0 / ((n_dim + 1.3) ** 2 + mu_eff))
    c_mu = float(
        min(1.0 - c_1, 2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((n_dim + 2.0) ** 2 + mu_eff))
    )
    chi_n = float(math.sqrt(n_dim) * (1.0 - 1.0 / (4.0 * n_dim) + 1.0 / (21.0 * n_dim * n_dim)))
    return {
        "lam": lam,
        "mu": mu,
        "weights": weights,
        "mu_eff": mu_eff,
        "c_sigma": c_sigma,
        "d_sigma": d_sigma,
        "c_c": c_c,
        "c_1": c_1,
        "c_mu": c_mu,
        "chi_n": chi_n,
    }


def _ellipse_metrics(covariance: np.ndarray) -> dict[str, float]:
    if covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        return {
            "major_radius": float("nan"),
            "minor_radius": float("nan"),
            "angle_deg": float("nan"),
            "condition": float("nan"),
            "corr": float("nan"),
        }
    try:
        eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    except np.linalg.LinAlgError:
        return {
            "major_radius": float("nan"),
            "minor_radius": float("nan"),
            "angle_deg": float("nan"),
            "condition": float("nan"),
            "corr": float("nan"),
        }
    eigenvalues = np.clip(np.asarray(eigenvalues, dtype=float), 0.0, None)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    major_radius = float(math.sqrt(eigenvalues[0])) if eigenvalues.size else float("nan")
    minor_radius = float(math.sqrt(eigenvalues[1])) if eigenvalues.size > 1 else float("nan")
    angle_deg = float(np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])))
    if minor_radius > 0.0 and math.isfinite(minor_radius):
        condition = float(major_radius / minor_radius)
    else:
        condition = float("inf") if math.isfinite(major_radius) else float("nan")
    denom = math.sqrt(max(float(covariance[0, 0]), 0.0) * max(float(covariance[1, 1]), 0.0))
    corr = float(covariance[0, 1] / denom) if denom > 0.0 else float("nan")
    return {
        "major_radius": major_radius,
        "minor_radius": minor_radius,
        "angle_deg": angle_deg,
        "condition": condition,
        "corr": corr,
    }


def _distribution_row(
    generation: int,
    *,
    mean: np.ndarray,
    sigma: float,
    normalized_covariance: np.ndarray,
    generation_best_score: float,
    global_best_score: float,
    n_completed_evals: int,
) -> dict[str, Any]:
    actual_covariance = float(sigma) * float(sigma) * np.asarray(normalized_covariance, dtype=float)
    metrics = _ellipse_metrics(actual_covariance)
    return {
        "generation": int(generation),
        "mean_r1": float(mean[0]),
        "mean_r2": float(mean[1]),
        "sigma": float(sigma),
        "norm_cov_11": float(normalized_covariance[0, 0]),
        "norm_cov_12": float(normalized_covariance[0, 1]),
        "norm_cov_22": float(normalized_covariance[1, 1]),
        "cov_11": float(actual_covariance[0, 0]),
        "cov_12": float(actual_covariance[0, 1]),
        "cov_22": float(actual_covariance[1, 1]),
        "major_radius": float(metrics["major_radius"]),
        "minor_radius": float(metrics["minor_radius"]),
        "angle_deg": float(metrics["angle_deg"]),
        "condition": float(metrics["condition"]),
        "corr": float(metrics["corr"]),
        "generation_best_score": float(generation_best_score),
        "global_best_score": float(global_best_score),
        "n_completed_evals": int(n_completed_evals),
    }


def _write_eval_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "eval_id",
        "generation",
        "candidate_index",
        "sampled_r1",
        "sampled_r2",
        "eval_r1",
        "eval_r2",
        "score",
        "selected",
        "generation_best",
        "out_of_bounds",
        "best_score_so_far",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(header, "") for header in headers])


def _write_generation_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "generation",
        "mean_r1",
        "mean_r2",
        "sigma",
        "norm_cov_11",
        "norm_cov_12",
        "norm_cov_22",
        "cov_11",
        "cov_12",
        "cov_22",
        "major_radius",
        "minor_radius",
        "angle_deg",
        "condition",
        "corr",
        "generation_best_score",
        "global_best_score",
        "n_completed_evals",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(header, "") for header in headers])


def _distribution_bounds(
    eval_rows: list[dict[str, Any]],
    generation_rows: list[dict[str, Any]],
    target_point: CouplingPoint,
) -> tuple[float, float, float, float]:
    xs = [float(target_point.r1)]
    ys = [float(target_point.r2)]
    for row in generation_rows:
        xs.append(float(row["mean_r1"]))
        ys.append(float(row["mean_r2"]))
        major = _safe_float(row.get("major_radius"))
        if math.isfinite(major):
            xs.extend([float(row["mean_r1"]) - 2.2 * major, float(row["mean_r1"]) + 2.2 * major])
            ys.extend([float(row["mean_r2"]) - 2.2 * major, float(row["mean_r2"]) + 2.2 * major])
    for row in eval_rows:
        x_value = _safe_float(row.get("sampled_r1"))
        y_value = _safe_float(row.get("sampled_r2"))
        if math.isfinite(x_value) and math.isfinite(y_value):
            xs.append(x_value)
            ys.append(y_value)
    x_min = min(xs)
    x_max = max(xs)
    y_min = min(ys)
    y_max = max(ys)
    x_pad = max(0.02, 0.08 * max(x_max - x_min, 1.0e-9))
    y_pad = max(0.02, 0.08 * max(y_max - y_min, 1.0e-9))
    return (x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad)


def _plot_state(
    path: Path,
    *,
    eval_rows: list[dict[str, Any]],
    generation_rows: list[dict[str, Any]],
    target_point: CouplingPoint,
    metric_key: str,
    ghost_generations: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0), constrained_layout=True)
    trajectory_ax = axes[0, 0]
    score_ax = axes[0, 1]
    scale_ax = axes[1, 0]
    diag_ax = axes[1, 1]

    for axis in axes.ravel():
        axis.set_facecolor("#fffdfa")
        axis.grid(True, color="#d1d5db", linewidth=0.8, alpha=0.8)

    finite_eval_rows = [row for row in eval_rows if math.isfinite(_safe_float(row.get("score")))]
    sampled_x = [_safe_float(row.get("sampled_r1")) for row in finite_eval_rows]
    sampled_y = [_safe_float(row.get("sampled_r2")) for row in finite_eval_rows]
    sampled_scores = np.asarray([_safe_float(row.get("score")) for row in finite_eval_rows], dtype=float)
    generations = np.asarray([int(row.get("generation", 0)) for row in finite_eval_rows], dtype=int)
    if finite_eval_rows:
        scatter = trajectory_ax.scatter(
            sampled_x,
            sampled_y,
            c=generations,
            cmap="viridis",
            s=34,
            alpha=0.75,
            edgecolors="none",
            zorder=1,
        )
        fig.colorbar(scatter, ax=trajectory_ax, fraction=0.046, pad=0.04, label="generation")

    ghosts = generation_rows[-max(int(ghost_generations), 1) :]
    for index, row in enumerate(ghosts):
        actual_covariance = np.asarray(
            [[float(row["cov_11"]), float(row["cov_12"])], [float(row["cov_12"]), float(row["cov_22"])]],
            dtype=float,
        )
        metrics = _ellipse_metrics(actual_covariance)
        if not math.isfinite(metrics["major_radius"]) or not math.isfinite(metrics["minor_radius"]):
            continue
        is_current = index == len(ghosts) - 1
        alpha_face = 0.20 if is_current else 0.04 + 0.10 * (index / max(len(ghosts) - 1, 1))
        alpha_edge = 0.85 if is_current else 0.22
        width_1sigma = 2.0 * float(metrics["major_radius"])
        height_1sigma = 2.0 * float(metrics["minor_radius"])
        ellipse_1sigma = mpatches.Ellipse(
            (float(row["mean_r1"]), float(row["mean_r2"])),
            width=width_1sigma,
            height=height_1sigma,
            angle=float(metrics["angle_deg"]),
            facecolor="#f59e0b",
            edgecolor="#d97706",
            alpha=alpha_face,
            linewidth=1.2 if is_current else 0.6,
            zorder=3 if is_current else 0,
        )
        trajectory_ax.add_patch(ellipse_1sigma)
        if is_current:
            ellipse_2sigma = mpatches.Ellipse(
                (float(row["mean_r1"]), float(row["mean_r2"])),
                width=2.0 * width_1sigma,
                height=2.0 * height_1sigma,
                angle=float(metrics["angle_deg"]),
                facecolor="none",
                edgecolor="#b45309",
                alpha=alpha_edge,
                linewidth=1.0,
                linestyle="--",
                zorder=4,
            )
            trajectory_ax.add_patch(ellipse_2sigma)
        trajectory_ax.plot(
            float(row["mean_r1"]),
            float(row["mean_r2"]),
            "x",
            color="#92400e" if is_current else "#d97706",
            ms=6 if is_current else 4,
            mew=1.4 if is_current else 0.8,
            zorder=5 if is_current else 2,
        )

    mean_x = [float(row["mean_r1"]) for row in generation_rows]
    mean_y = [float(row["mean_r2"]) for row in generation_rows]
    if generation_rows:
        trajectory_ax.plot(mean_x, mean_y, "-", color="#111827", linewidth=1.4, alpha=0.85, zorder=2)
        for row in generation_rows[-max(int(ghost_generations), 1) :]:
            trajectory_ax.annotate(
                str(int(row["generation"])),
                (float(row["mean_r1"]), float(row["mean_r2"])),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                color="#0f172a",
            )

    if finite_eval_rows:
        best_row = min(finite_eval_rows, key=lambda row: float(row["score"]))
        trajectory_ax.scatter(
            [float(best_row["eval_r1"])],
            [float(best_row["eval_r2"])],
            marker="*",
            s=170,
            color="#059669",
            edgecolors="#064e3b",
            linewidths=0.8,
            zorder=6,
        )
    trajectory_ax.scatter(
        [float(target_point.r1)],
        [float(target_point.r2)],
        marker="x",
        s=90,
        color="#dc2626",
        linewidths=1.8,
        zorder=7,
    )
    trajectory_ax.set_title("distribution learning in parameter space")
    trajectory_ax.set_xlabel("r1")
    trajectory_ax.set_ylabel("r2")
    x_min, x_max, y_min, y_max = _distribution_bounds(eval_rows, generation_rows, target_point)
    trajectory_ax.set_xlim(x_min, x_max)
    trajectory_ax.set_ylim(y_min, y_max)

    if finite_eval_rows:
        eval_ids = [int(row["eval_id"]) for row in finite_eval_rows]
        global_best = np.minimum.accumulate(sampled_scores)
        score_ax.plot(eval_ids, sampled_scores, "o", color="#94a3b8", markersize=4, alpha=0.8, label="all evals")
        score_ax.plot(eval_ids, global_best, "-", color="#111827", linewidth=1.8, label="best so far")
    generation_ids = [int(row["generation"]) for row in generation_rows if int(row["generation"]) > 0]
    generation_best_scores = [
        _safe_float(row.get("generation_best_score")) for row in generation_rows if int(row["generation"]) > 0
    ]
    if generation_ids:
        score_ax.plot(generation_ids, generation_best_scores, "-s", color="#ea580c", linewidth=1.4, label="generation best")
    score_ax.set_title(f"{metric_key} by evaluation")
    score_ax.set_xlabel("eval or generation index")
    score_ax.set_ylabel(metric_key)
    score_ax.legend(fontsize=8)

    scale_generation_ids = [int(row["generation"]) for row in generation_rows]
    sigmas = [_safe_float(row.get("sigma")) for row in generation_rows]
    major_radii = [_safe_float(row.get("major_radius")) for row in generation_rows]
    minor_radii = [_safe_float(row.get("minor_radius")) for row in generation_rows]
    scale_ax.plot(scale_generation_ids, sigmas, "-o", color="#1d4ed8", label="global sigma")
    scale_ax.plot(scale_generation_ids, major_radii, "-s", color="#7c3aed", label="major 1σ radius")
    scale_ax.plot(scale_generation_ids, minor_radii, "-^", color="#059669", label="minor 1σ radius")
    scale_ax.set_title("distribution scale by generation")
    scale_ax.set_xlabel("generation")
    scale_ax.set_ylabel("scale")
    scale_ax.legend(fontsize=8)

    correlations = [_safe_float(row.get("corr")) for row in generation_rows]
    angles = [_safe_float(row.get("angle_deg")) for row in generation_rows]
    diag_ax.axhline(0.0, color="#94a3b8", linewidth=0.8, alpha=0.8)
    diag_ax.plot(scale_generation_ids, correlations, "-o", color="#2563eb", label="corr")
    diag_ax.set_title("distribution orientation diagnostics")
    diag_ax.set_xlabel("generation")
    diag_ax.set_ylabel("corr")
    diag_ax.set_ylim(-1.05, 1.05)
    diag_ax_right = diag_ax.twinx()
    diag_ax_right.plot(scale_generation_ids, angles, "--s", color="#ea580c", label="angle")
    diag_ax_right.set_ylabel("angle (deg)")

    summary_lines: list[str] = []
    if finite_eval_rows:
        best_row = min(finite_eval_rows, key=lambda row: float(row["score"]))
        summary_lines.append(f"best score: {float(best_row['score']):.3e}")
        summary_lines.append(f"best point: ({float(best_row['eval_r1']):.4f}, {float(best_row['eval_r2']):.4f})")
    if generation_rows:
        current = generation_rows[-1]
        summary_lines.append(f"current sigma: {float(current['sigma']):.4f}")
        summary_lines.append(f"axes: {float(current['major_radius']):.4f}/{float(current['minor_radius']):.4f}")
        summary_lines.append(f"corr: {float(current['corr']):+.3f}")
    if summary_lines:
        trajectory_ax.text(
            0.02,
            0.98,
            "\n".join(summary_lines),
            transform=trajectory_ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="#0f172a",
            bbox={"boxstyle": "round,pad=0.35", "fc": "#fff7ed", "ec": "#fdba74", "lw": 1.0},
        )

    fig.suptitle("Current-stack CMA-ES with distribution diagnostics", fontsize=16, y=0.995)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _evaluate_direct_score(
    runner: Any,
    point: CouplingPoint,
    *,
    eval_id: int,
    generation: int,
    candidate_index: int,
) -> dict[str, Any]:
    point_row = runner.score_row(point, point, step_index=int(eval_id))
    score = float(point_row[runner.metric_key])
    return {
        "eval_id": int(eval_id),
        "generation": int(generation),
        "candidate_index": int(candidate_index),
        "sampled_r1": f"{float(point.r1):.10f}",
        "sampled_r2": f"{float(point.r2):.10f}",
        "eval_r1": f"{float(point.r1):.10f}",
        "eval_r2": f"{float(point.r2):.10f}",
        "score": _format_float(float(score)),
        "selected": "False",
        "generation_best": "False",
        "out_of_bounds": "False",
        "best_score_so_far": _format_float(float(score)),
    }


def main() -> None:
    cli_args = parse_args()
    if int(cli_args.additional_gens) < 0:
        raise ValueError("additional_gens must be non-negative")

    resume_summary: dict[str, Any] | None = None
    resume_state: dict[str, Any] | None = None
    args = cli_args
    if cli_args.resume_root:
        args, resume_summary = _resume_args_from_summary(cli_args)

    if int(args.popsize) < 2:
        raise ValueError("popsize must be at least 2")
    if int(args.max_evals) < 1:
        raise ValueError("max_evals must be at least 1")
    if not math.isfinite(float(args.sigma0)) or float(args.sigma0) <= 0.0:
        raise ValueError(f"sigma0 must be positive, got {args.sigma0}")
    if int(args.plot_every_gens) < 1:
        raise ValueError("plot_every_gens must be at least 1")

    if str(args.target_mode) == "acute456":
        runner = Acute456Runner(args)
    else:
        runner = GradientFlowRunner(args)
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    eval_tsv = output_root / "evals.tsv"
    generation_tsv = output_root / "generations.tsv"
    output_png = output_root / "trajectory.png"
    frames_dir = output_root / "frames_history"
    if bool(args.save_frames):
        frames_dir.mkdir(parents=True, exist_ok=True)

    n_dim = 2
    params = _cma_default_params(n_dim, int(args.popsize), dsigma=args.cma_dsigma, csigma=args.cma_csigma)
    rng = np.random.default_rng(int(args.cma_seed))

    if resume_summary is not None:
        eval_tsv_existing = output_root / "evals.tsv"
        generation_tsv_existing = output_root / "generations.tsv"
        if not eval_tsv_existing.exists() or not generation_tsv_existing.exists():
            raise FileNotFoundError(f"Cannot resume without evals.tsv and generations.tsv in {output_root}")
        resume_state = _replay_resume_state(
            args=args,
            summary=resume_summary,
            params=params,
            eval_rows=_read_tsv_rows(eval_tsv_existing),
            generation_rows=_read_tsv_rows(generation_tsv_existing),
        )
        for _ in range(int(resume_state["generation"])):
            rng.standard_normal((int(params["lam"]), n_dim))

    if resume_state is None:
        mean = np.asarray([float(args.start_r1), float(args.start_r2)], dtype=float)
        sigma = float(args.sigma0)
        covariance = np.eye(n_dim, dtype=float)
        p_sigma = np.zeros(n_dim, dtype=float)
        p_c = np.zeros(n_dim, dtype=float)
        basis = np.eye(n_dim, dtype=float)
        diag = np.ones(n_dim, dtype=float)
        bd = basis * diag
        eval_rows: list[dict[str, Any]] = []
        generation_rows: list[dict[str, Any]] = [
            _distribution_row(
                0,
                mean=mean,
                sigma=sigma,
                normalized_covariance=covariance,
                generation_best_score=float("nan"),
                global_best_score=float("nan"),
                n_completed_evals=0,
            )
        ]
        best_score = float("inf")
        best_point = CouplingPoint(float(args.start_r1), float(args.start_r2))
        best_eval_id = -1
        generation = 0
    else:
        mean = np.asarray(resume_state["mean"], dtype=float)
        sigma = float(resume_state["sigma"])
        covariance = np.asarray(resume_state["covariance"], dtype=float)
        p_sigma = np.asarray(resume_state["p_sigma"], dtype=float)
        p_c = np.asarray(resume_state["p_c"], dtype=float)
        basis = np.asarray(resume_state["basis"], dtype=float)
        diag = np.asarray(resume_state["diag"], dtype=float)
        bd = np.asarray(resume_state["bd"], dtype=float)
        eval_rows = list(resume_state["eval_rows"])
        generation_rows = list(resume_state["generation_rows"])
        best_score = float(resume_state["best_score"])
        best_point = resume_state["best_point"]
        best_eval_id = int(resume_state["best_eval_id"])
        generation = int(resume_state["generation"])

    previous_generation = int(generation)
    previous_eval_count = int(len(eval_rows))
    generation_limit = int(args.max_gens)
    eval_limit = int(args.max_evals)
    if resume_state is not None:
        generation_limit = int(previous_generation + int(cli_args.additional_gens))
        eval_limit = int(previous_eval_count + int(cli_args.additional_gens) * int(params["lam"]))
        if int(cli_args.additional_gens) > 0:
            summary_path = output_root / "summary.json"
            backup_path = output_root / "summary_before_resume.json"
            if summary_path.exists() and not backup_path.exists():
                backup_path.write_text(summary_path.read_text(encoding="utf-8"), encoding="utf-8")

    status = "completed"
    oob_penalty = 1.0e30

    try:
        while len(eval_rows) < eval_limit:
            generation += 1
            if generation_limit > 0 and generation > generation_limit:
                status = f"max_gens reached ({generation_limit})"
                break

            zs = rng.standard_normal((int(params["lam"]), n_dim))
            ys = zs @ bd.T
            xs = mean + sigma * ys

            batch_row_indices: list[int] = []
            batch_scores: list[float] = []
            completed_y: list[np.ndarray] = []
            completed_z: list[np.ndarray] = []
            for candidate_index in range(int(params["lam"])):
                if len(eval_rows) >= eval_limit:
                    break
                eval_id = len(eval_rows) + 1
                sampled = CouplingPoint(float(xs[candidate_index, 0]), float(xs[candidate_index, 1]))
                clipped = _clip_point(sampled, args)
                out_of_bounds = not (
                    math.isclose(float(sampled.r1), float(clipped.r1), abs_tol=1.0e-12)
                    and math.isclose(float(sampled.r2), float(clipped.r2), abs_tol=1.0e-12)
                )
                if out_of_bounds:
                    row = {
                        "eval_id": int(eval_id),
                        "generation": int(generation),
                        "candidate_index": int(candidate_index),
                        "sampled_r1": f"{float(sampled.r1):.10f}",
                        "sampled_r2": f"{float(sampled.r2):.10f}",
                        "eval_r1": "",
                        "eval_r2": "",
                        "score": _format_float(float(oob_penalty)),
                        "selected": "False",
                        "generation_best": "False",
                        "out_of_bounds": "True",
                        "best_score_so_far": _format_float(float(best_score)),
                    }
                    batch_scores.append(float(oob_penalty))
                else:
                    row = _evaluate_direct_score(
                        runner,
                        clipped,
                        eval_id=int(eval_id),
                        generation=int(generation),
                        candidate_index=int(candidate_index),
                    )
                    score = float(row["score"])
                    if score < best_score:
                        best_score = float(score)
                        best_point = clipped
                        best_eval_id = int(eval_id)
                    row["best_score_so_far"] = _format_float(float(best_score))
                    batch_scores.append(float(score))
                eval_rows.append(row)
                batch_row_indices.append(len(eval_rows) - 1)
                completed_y.append(np.asarray(ys[candidate_index], dtype=float))
                completed_z.append(np.asarray(zs[candidate_index], dtype=float))

                _write_eval_tsv(eval_tsv, eval_rows)
                _write_generation_tsv(generation_tsv, generation_rows)
                _plot_state(
                    output_png,
                    eval_rows=eval_rows,
                    generation_rows=generation_rows,
                    target_point=runner.target_point,
                    metric_key=runner.metric_key,
                    ghost_generations=int(args.ghost_generations),
                )

            if not batch_scores:
                status = "no completed evaluations"
                break

            fs = np.asarray(batch_scores, dtype=float)
            completed_y_array = np.asarray(completed_y, dtype=float)
            completed_z_array = np.asarray(completed_z, dtype=float)
            mu_used = max(1, min(int(params["mu"]), len(fs) // 2 + (len(fs) % 2)))
            if mu_used != int(params["mu"]):
                weights = np.asarray(params["weights"][:mu_used], dtype=float)
                weights = weights / np.sum(weights)
                mu_eff = float(1.0 / np.sum(np.square(weights)))
            else:
                weights = np.asarray(params["weights"], dtype=float)
                mu_eff = float(params["mu_eff"])
            order = np.argsort(fs)
            selected = order[:mu_used]
            selected_set = {int(index) for index in selected}
            generation_best_local = float(fs[int(order[0])])
            for local_index, eval_row_index in enumerate(batch_row_indices):
                eval_rows[eval_row_index]["selected"] = str(local_index in selected_set)
                eval_rows[eval_row_index]["generation_best"] = str(local_index == int(order[0]))

            y_selected = completed_y_array[selected]
            z_selected = completed_z_array[selected]
            y_weighted = weights @ y_selected
            z_weighted = weights @ z_selected
            old_mean = np.asarray(mean, dtype=float)
            mean = old_mean + sigma * y_weighted

            c_sigma = float(params["c_sigma"])
            d_sigma = float(params["d_sigma"])
            c_c = float(params["c_c"])
            c_1 = float(params["c_1"])
            c_mu = float(params["c_mu"])
            chi_n = float(params["chi_n"])

            p_sigma = (1.0 - c_sigma) * p_sigma + math.sqrt(c_sigma * (2.0 - c_sigma) * mu_eff) * (basis @ z_weighted)
            sigma = float(sigma * math.exp((c_sigma / d_sigma) * (np.linalg.norm(p_sigma) / chi_n - 1.0)))

            hs_threshold = (1.4 + 2.0 / (n_dim + 1.0)) * chi_n
            denom = math.sqrt(max(1.0 - (1.0 - c_sigma) ** (2 * generation), 1.0e-30))
            h_sigma = 1.0 if (np.linalg.norm(p_sigma) / denom) < hs_threshold else 0.0
            p_c = (1.0 - c_c) * p_c + h_sigma * math.sqrt(c_c * (2.0 - c_c) * mu_eff) * y_weighted

            rank_mu = y_selected.T @ (weights[:, None] * y_selected)
            rank_1 = np.outer(p_c, p_c) + (1.0 - h_sigma) * c_c * (2.0 - c_c) * covariance
            covariance = (1.0 - c_1 - c_mu) * covariance + c_1 * rank_1 + c_mu * rank_mu
            covariance = 0.5 * (covariance + covariance.T)

            try:
                eigenvalues, basis = np.linalg.eigh(covariance)
            except np.linalg.LinAlgError:
                covariance = np.eye(n_dim, dtype=float)
                basis = np.eye(n_dim, dtype=float)
                eigenvalues = np.ones(n_dim, dtype=float)
                p_sigma = np.zeros(n_dim, dtype=float)
                p_c = np.zeros(n_dim, dtype=float)
            eigenvalues = np.clip(np.asarray(eigenvalues, dtype=float), 1.0e-30, None)
            diag = np.sqrt(eigenvalues)
            bd = basis * diag

            generation_rows.append(
                _distribution_row(
                    generation,
                    mean=mean,
                    sigma=sigma,
                    normalized_covariance=covariance,
                    generation_best_score=generation_best_local,
                    global_best_score=best_score,
                    n_completed_evals=len(eval_rows),
                )
            )
            _write_eval_tsv(eval_tsv, eval_rows)
            _write_generation_tsv(generation_tsv, generation_rows)
            _plot_state(
                output_png,
                eval_rows=eval_rows,
                generation_rows=generation_rows,
                target_point=runner.target_point,
                metric_key=runner.metric_key,
                ghost_generations=int(args.ghost_generations),
            )
            if bool(args.save_frames) and (generation % int(args.plot_every_gens) == 0):
                frame_path = frames_dir / f"frame_{generation:04d}.png"
                _plot_state(
                    frame_path,
                    eval_rows=eval_rows,
                    generation_rows=generation_rows,
                    target_point=runner.target_point,
                    metric_key=runner.metric_key,
                    ghost_generations=int(args.ghost_generations),
                )

            sigma_step = float(sigma * np.max(diag))
            if sigma_step < float(args.tolx):
                status = f"converged_tolx ({sigma_step:.3e} < {float(args.tolx):.3e})"
                break
            if len(fs) >= 2:
                spread = float(np.max(fs) - np.min(fs))
                if spread < float(args.tolfun):
                    status = f"converged_tolfun ({spread:.3e} < {float(args.tolfun):.3e})"
                    break
    except KeyboardInterrupt:
        status = "interrupted"

    if status == "completed" and resume_state is not None and int(cli_args.additional_gens) > 0:
        status = f"resumed_completed (+{int(cli_args.additional_gens)} gens)"
    if status == "completed" and resume_state is not None and int(cli_args.additional_gens) == 0:
        status = "resume_state_validated"

    final_actual_covariance = float(sigma) * float(sigma) * np.asarray(covariance, dtype=float)
    summary = {
        "args": vars(args),
        "metric_key": runner.metric_key,
        "target_point": _point_to_json(runner.target_point),
        "start_point": {"r1": float(args.start_r1), "r2": float(args.start_r2)},
        "best_point": _point_to_json(best_point),
        "best_eval_id": int(best_eval_id),
        "best_score": float(best_score),
        "status": status,
        "n_evals": len(eval_rows),
        "n_generations": max(len(generation_rows) - 1, 0),
        "evals_tsv": str(eval_tsv),
        "generations_tsv": str(generation_tsv),
        "trajectory_png": str(output_png),
        "frames_dir": str(frames_dir) if bool(args.save_frames) else None,
        "final_gaussian": {
            "mean": [float(mean[0]), float(mean[1])],
            "sigma": float(sigma),
            "normalized_covariance": covariance.tolist(),
            "actual_covariance": final_actual_covariance.tolist(),
        },
    }
    if resume_state is not None:
        summary["resume_source_root"] = str(Path(str(cli_args.resume_root)).resolve())
        summary["resume_from_generation"] = int(previous_generation)
        summary["resume_from_evals"] = int(previous_eval_count)
        summary["continued_generations"] = int(max(generation - previous_generation, 0))
        summary["continued_evals"] = int(max(len(eval_rows) - previous_eval_count, 0))
        summary["requested_additional_gens"] = int(cli_args.additional_gens)
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()