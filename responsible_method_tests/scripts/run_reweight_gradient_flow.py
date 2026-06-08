#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent
REPO_ROOT = RESPONSIBLE_ROOT.parent
KFC_ROOT = REPO_ROOT / "K_from_continuum"
if str(KFC_ROOT) not in sys.path:
    sys.path.insert(0, str(KFC_ROOT))

from workflow_common import ensure_simulator  # noqa: E402

from plot_raw_reweighting_fss import _ratio_observables  # noqa: E402
from plot_reweighting_ratio_heatmaps import (  # noqa: E402
    BALANCE_CHOICES,
    _aggregate_plot_keys,
)
from plot_reweighting_refined_surfaces import (  # noqa: E402
    _finite_difference,
    _score_row_from_single_donor,
)
from test_geometry_match_grid_interpolation import (  # noqa: E402
    DEFAULT_EXECUTION,
    CouplingPoint,
    _block_slices,
    _parse_selected_bundle,
    _run_point,
    _selected_specs_for_size,
    _unique_ints_preserve_order,
)


STEP_MODE_CHOICES = ("normalized", "raw")
STEP_ADAPTATION_CHOICES = ("fixed", "scalar_snr", "covariance")
ACCEPT_RULE_CHOICES = ("gradient_only", "always", "direct_decrease")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a discrete gradient-flow loop driven by direct Monte Carlo at the current point and "
            "local reweight-derived score gradients around that point."
        )
    )
    parser.add_argument("--start-r1", type=float, default=1.0)
    parser.add_argument("--start-r2", type=float, default=1.0)
    parser.add_argument("--step-size", type=float, default=0.025)
    parser.add_argument("--num-steps", type=int, default=5)
    parser.add_argument("--step-mode", choices=STEP_MODE_CHOICES, default="normalized")
    parser.add_argument(
        "--step-adaptation",
        choices=STEP_ADAPTATION_CHOICES,
        default="fixed",
        help=(
            "Proposal-step adaptation mode. 'fixed' preserves the current behavior, 'scalar_snr' rescales "
            "the base step length using the local gradient SNR, and 'covariance' also rotates the descent "
            "direction using the jackknife gradient covariance."
        ),
    )
    parser.add_argument(
        "--adaptive-step-snr-reference",
        type=float,
        default=3.0,
        help=(
            "Reference SNR for adaptive step scaling. Scalar and covariance modes use "
            "clip(step_snr / adaptive_step_snr_reference, adaptive_step_scale_min, adaptive_step_scale_max)."
        ),
    )
    parser.add_argument(
        "--adaptive-step-scale-min",
        type=float,
        default=0.5,
        help="Minimum multiplicative factor applied to step_size in adaptive step modes.",
    )
    parser.add_argument(
        "--adaptive-step-scale-max",
        type=float,
        default=3.0,
        help="Maximum multiplicative factor applied to step_size in adaptive step modes.",
    )
    parser.add_argument(
        "--min-gradient-snr",
        type=float,
        default=1.0,
        help=(
            "Require the local gradient norm to exceed this jackknife signal-to-noise ratio "
            "before a proposal is evaluated with direct Monte Carlo. Use 0 to disable the gate."
        ),
    )
    parser.add_argument(
        "--gradient-jackknife-blocks",
        type=int,
        default=16,
        help="Nominal number of leave-one-block-out replicas used for the gradient SNR estimate.",
    )
    parser.add_argument(
        "--gradient-step",
        type=float,
        default=None,
        help="Centered-difference offset for the local reweight gradient. Defaults to step_size/2.",
    )
    parser.add_argument(
        "--accept-rule",
        choices=ACCEPT_RULE_CHOICES,
        default="gradient_only",
        help=(
            "Update acceptance mode. Gradient-driven updates now always advance once the local gradient passes "
            "the SNR gate; legacy values are accepted for CLI compatibility and ignored."
        ),
    )
    parser.add_argument("--balance-mode", choices=BALANCE_CHOICES, default="none")
    parser.add_argument(
        "--metric-key",
        default=None,
        help="Objective metric to flow downhill on. Defaults to the total aggregate metric for the chosen balance mode.",
    )
    parser.add_argument("--target-size", type=int, default=64)
    parser.add_argument("--target-r1", type=float, default=1.0)
    parser.add_argument("--target-r2", type=float, default=1.0)
    parser.add_argument("--fit-sizes", nargs="+", type=int, default=[32, 28, 24, 20, 16, 12, 8, 4])
    parser.add_argument("--n-traj", type=int, default=100000)
    parser.add_argument("--n-therm", type=int, default=10000)
    parser.add_argument("--n-skip", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=2026060801)
    parser.add_argument("--r1-min", type=float, default=0.1)
    parser.add_argument("--r1-max", type=float, default=10.0)
    parser.add_argument("--r2-min", type=float, default=0.1)
    parser.add_argument("--r2-max", type=float, default=10.0)
    parser.add_argument(
        "--output-root",
        default=str(RESPONSIBLE_ROOT / "results" / "reweight_gradient_flow_iso111"),
        help="Root directory for trajectory outputs and raw Monte Carlo payloads.",
    )
    return parser.parse_args()


def _resolve_gradient_step(step_size: float, requested_step: float | None) -> float:
    if requested_step is not None:
        gradient_step = float(requested_step)
    else:
        gradient_step = 0.5 * float(step_size)
    if not math.isfinite(gradient_step) or gradient_step <= 0.0:
        raise ValueError(f"gradient step must be positive, got {gradient_step}")
    return gradient_step


def _clip_point(point: CouplingPoint, args: argparse.Namespace) -> CouplingPoint:
    return CouplingPoint(
        r1=min(max(float(point.r1), float(args.r1_min)), float(args.r1_max)),
        r2=min(max(float(point.r2), float(args.r2_min)), float(args.r2_max)),
    )


def _point_to_json(point: CouplingPoint) -> dict[str, float]:
    return {"r1": float(point.r1), "r2": float(point.r2)}


def _subset_payload(payload: dict[str, Any], lo: int, hi: int) -> dict[str, Any]:
    n_samples = int(np.asarray(payload["mag"], dtype=float).size)
    mask = np.ones(n_samples, dtype=bool)
    mask[int(lo) : int(hi)] = False
    return {
        "path": str(payload["path"]),
        "corr": {key: np.asarray(values, dtype=float)[mask] for key, values in dict(payload["corr"]).items()},
        "mag": np.asarray(payload["mag"], dtype=float)[mask],
        "e1": np.asarray(payload["e1"], dtype=float)[mask],
        "e2": np.asarray(payload["e2"], dtype=float)[mask],
        "e3": np.asarray(payload["e3"], dtype=float)[mask],
        "beta": float(payload["beta"]),
        "K": tuple(float(value) for value in payload["K"]),
        "n_sites": int(payload["n_sites"]),
        "header": dict(payload["header"]),
    }


def _jackknife_sigma(values: list[float]) -> float:
    array = np.asarray([float(value) for value in values if math.isfinite(float(value))], dtype=float)
    if array.size < 2:
        return float("nan")
    mean_value = float(np.mean(array))
    return float(np.sqrt((array.size - 1.0) / array.size * np.sum(np.square(array - mean_value))))


def _jackknife_covariance(vectors: list[tuple[float, float]]) -> np.ndarray:
    array = np.asarray(vectors, dtype=float)
    if array.ndim != 2 or array.shape[0] < 2:
        return np.full((2, 2), np.nan, dtype=float)
    centered = array - np.mean(array, axis=0, keepdims=True)
    return (array.shape[0] - 1.0) / array.shape[0] * (centered.T @ centered)


def _signal_to_noise(value: float, sigma: float) -> float:
    if not math.isfinite(value) or not math.isfinite(sigma) or sigma <= 0.0:
        return float("nan")
    return float(abs(value) / sigma)


def _vector_snr(covariance: np.ndarray, grad_r1: float, grad_r2: float) -> float:
    if covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        return float("nan")
    determinant = float(np.linalg.det(covariance))
    if determinant <= 0.0:
        return float("nan")
    vector = np.asarray([float(grad_r1), float(grad_r2)], dtype=float)
    return float(math.sqrt(vector @ np.linalg.inv(covariance) @ vector))


def _adaptive_step_scale(*, grad_norm_snr: float, gradient_vector_snr: float, args: argparse.Namespace) -> tuple[float, float]:
    if str(args.step_adaptation) == "fixed":
        return 1.0, float("nan")
    step_snr = float(gradient_vector_snr) if math.isfinite(float(gradient_vector_snr)) else float(grad_norm_snr)
    if not math.isfinite(step_snr):
        step_snr = 0.0
    scale = step_snr / float(args.adaptive_step_snr_reference)
    scale = min(max(float(scale), float(args.adaptive_step_scale_min)), float(args.adaptive_step_scale_max))
    return float(scale), float(step_snr)


def _regularized_precision(covariance: np.ndarray) -> np.ndarray | None:
    if covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        return None
    try:
        eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    except np.linalg.LinAlgError:
        return None
    max_eigenvalue = float(np.max(np.abs(eigenvalues)))
    if not math.isfinite(max_eigenvalue) or max_eigenvalue <= 0.0:
        return None
    floor = max(1.0e-12, 1.0e-6 * max_eigenvalue)
    clipped = np.maximum(eigenvalues, floor)
    return eigenvectors @ np.diag(1.0 / clipped) @ eigenvectors.T


def _descent_vector(
    *,
    grad_r1: float,
    grad_r2: float,
    gradient_covariance: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    gradient = np.asarray([float(grad_r1), float(grad_r2)], dtype=float)
    if str(args.step_adaptation) != "covariance":
        return gradient
    precision = _regularized_precision(gradient_covariance)
    if precision is None:
        return gradient
    preconditioned = precision @ gradient
    if not np.all(np.isfinite(preconditioned)):
        return gradient
    norm = float(np.linalg.norm(preconditioned))
    if not math.isfinite(norm) or norm <= 0.0:
        return gradient
    return preconditioned


def _gradient_from_payloads(
    donor_point: CouplingPoint,
    *,
    fit_sizes: list[int],
    payloads: dict[tuple[int, str], dict[str, Any]],
    target_rows: dict[str, dict[str, float]],
    target_x: float,
    gradient_step: float,
    metric_key: str,
) -> dict[str, float]:
    direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]] = {}
    reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]] = {}

    def _score(eval_point: CouplingPoint) -> dict[str, Any]:
        return _score_row_from_single_donor(
            donor_point,
            eval_point,
            fit_sizes=fit_sizes,
            target_rows=target_rows,
            target_x=float(target_x),
            payloads=payloads,
            direct_family_cache=direct_family_cache,
            reweighted_family_cache=reweighted_family_cache,
        )

    center_row = _score(donor_point)
    r1_minus_row = _score(CouplingPoint(float(donor_point.r1) - float(gradient_step), float(donor_point.r2)))
    r1_plus_row = _score(CouplingPoint(float(donor_point.r1) + float(gradient_step), float(donor_point.r2)))
    r2_minus_row = _score(CouplingPoint(float(donor_point.r1), float(donor_point.r2) - float(gradient_step)))
    r2_plus_row = _score(CouplingPoint(float(donor_point.r1), float(donor_point.r2) + float(gradient_step)))

    center_score = float(center_row[metric_key])
    grad_r1 = float(
        _finite_difference(center_score, float(r1_minus_row[metric_key]), float(r1_plus_row[metric_key]), float(gradient_step))
    )
    grad_r2 = float(
        _finite_difference(center_score, float(r2_minus_row[metric_key]), float(r2_plus_row[metric_key]), float(gradient_step))
    )
    return {
        "center_score": center_score,
        "r1_minus_score": float(r1_minus_row[metric_key]),
        "r1_plus_score": float(r1_plus_row[metric_key]),
        "r2_minus_score": float(r2_minus_row[metric_key]),
        "r2_plus_score": float(r2_plus_row[metric_key]),
        "grad_r1": float(grad_r1),
        "grad_r2": float(grad_r2),
        "grad_norm": float(math.hypot(grad_r1, grad_r2)),
    }


def _format_float(value: float, fmt: str) -> str:
    return format(float(value), fmt) if math.isfinite(float(value)) else "nan"


class GradientFlowRunner:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_root = Path(args.output_root).resolve()
        self.output_root.mkdir(parents=True, exist_ok=True)
        self.raw_root = self.output_root / "raw"
        self.target_root = self.output_root / "_target"
        self.exe = ensure_simulator(DEFAULT_EXECUTION)
        self.fit_sizes = _unique_ints_preserve_order(list(args.fit_sizes))
        if not self.fit_sizes:
            raise ValueError("fit_sizes must not be empty")
        self.labels = [label for label, _ in _selected_specs_for_size(self.fit_sizes[0])]
        self.target_point = CouplingPoint(float(args.target_r1), float(args.target_r2))
        self.target_x = 1.0 / float(args.target_size)
        self.metric_keys = _aggregate_plot_keys(str(args.balance_mode))
        self.metric_key = str(args.metric_key) if args.metric_key else str(self.metric_keys[-1])
        if self.metric_key not in self.metric_keys:
            raise ValueError(f"metric_key must be one of {self.metric_keys}, got {self.metric_key}")
        self.gradient_step = _resolve_gradient_step(float(args.step_size), args.gradient_step)
        self.payload_cache: dict[tuple[int, str], dict[str, Any]] = {}
        self.direct_family_cache: dict[tuple[int, str], dict[str, dict[str, float]]] = {}
        self.reweighted_family_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]] = {}
        self.target_payload = self._load_target_payload()
        self.target_rows = _ratio_observables(self.target_payload)

    def _seed_for(self, *, step_index: int, size: int, kind_offset: int) -> int:
        size_index = self.fit_sizes.index(int(size))
        return int(self.args.seed_base) + 100000 * int(step_index) + 1000 * int(kind_offset) + int(size_index)

    def _ensure_payload(self, point: CouplingPoint, *, step_index: int) -> None:
        for size in self.fit_sizes:
            key = (int(size), point.tag)
            if key in self.payload_cache:
                continue
            bundle_path = _run_point(
                self.exe,
                self.output_root,
                size=int(size),
                point=point,
                n_traj=int(self.args.n_traj),
                n_therm=int(self.args.n_therm),
                n_skip=int(self.args.n_skip),
                seed=self._seed_for(step_index=step_index, size=int(size), kind_offset=0),
            )
            self.payload_cache[key] = _parse_selected_bundle(bundle_path, self.labels)

    def _load_target_payload(self) -> dict[str, Any]:
        bundle_path = _run_point(
            self.exe,
            self.target_root,
            size=int(self.args.target_size),
            point=self.target_point,
            n_traj=int(self.args.n_traj),
            n_therm=int(self.args.n_therm),
            n_skip=int(self.args.n_skip),
            seed=int(self.args.seed_base) + 990000,
        )
        return _parse_selected_bundle(bundle_path, self.labels)

    def score_row(self, donor_point: CouplingPoint, eval_point: CouplingPoint, *, step_index: int) -> dict[str, Any]:
        self._ensure_payload(donor_point, step_index=step_index)
        return _score_row_from_single_donor(
            donor_point,
            eval_point,
            fit_sizes=self.fit_sizes,
            target_rows=self.target_rows,
            target_x=float(self.target_x),
            payloads=self.payload_cache,
            direct_family_cache=self.direct_family_cache,
            reweighted_family_cache=self.reweighted_family_cache,
        )

    def gradient(self, point: CouplingPoint, *, step_index: int) -> dict[str, Any]:
        self._ensure_payload(point, step_index=step_index)
        full_payload = _gradient_from_payloads(
            point,
            fit_sizes=self.fit_sizes,
            payloads=self.payload_cache,
            target_rows=self.target_rows,
            target_x=float(self.target_x),
            gradient_step=float(self.gradient_step),
            metric_key=self.metric_key,
        )
        center_row = self.score_row(point, point, step_index=step_index)

        point_slices = {
            int(size): _block_slices(
                int(np.asarray(self.payload_cache[(int(size), point.tag)]["mag"], dtype=float).size),
                int(self.args.gradient_jackknife_blocks),
            )
            for size in self.fit_sizes
        }
        target_slices = _block_slices(
            int(np.asarray(self.target_payload["mag"], dtype=float).size),
            int(self.args.gradient_jackknife_blocks),
        )
        n_replicas = min(len(target_slices), *(len(slices) for slices in point_slices.values()))

        jk_vectors: list[tuple[float, float]] = []
        jk_norms: list[float] = []
        for replica_index in range(n_replicas):
            subset_payloads = {
                (int(size), point.tag): _subset_payload(
                    self.payload_cache[(int(size), point.tag)],
                    *point_slices[int(size)][replica_index],
                )
                for size in self.fit_sizes
            }
            target_subset = _subset_payload(self.target_payload, *target_slices[replica_index])
            subset_target_rows = _ratio_observables(target_subset)
            try:
                subset_gradient = _gradient_from_payloads(
                    point,
                    fit_sizes=self.fit_sizes,
                    payloads=subset_payloads,
                    target_rows=subset_target_rows,
                    target_x=float(self.target_x),
                    gradient_step=float(self.gradient_step),
                    metric_key=self.metric_key,
                )
            except Exception:
                continue
            jk_vectors.append((float(subset_gradient["grad_r1"]), float(subset_gradient["grad_r2"])))
            jk_norms.append(float(subset_gradient["grad_norm"]))

        covariance = _jackknife_covariance(jk_vectors)
        grad_r1_sigma = _jackknife_sigma([value[0] for value in jk_vectors])
        grad_r2_sigma = _jackknife_sigma([value[1] for value in jk_vectors])
        grad_norm_sigma = _jackknife_sigma(jk_norms)
        grad_r1 = float(full_payload["grad_r1"])
        grad_r2 = float(full_payload["grad_r2"])
        grad_norm = float(full_payload["grad_norm"])
        return {
            "center_row": center_row,
            "grad_r1": grad_r1,
            "grad_r2": grad_r2,
            "grad_norm": grad_norm,
            "grad_r1_sigma": float(grad_r1_sigma),
            "grad_r2_sigma": float(grad_r2_sigma),
            "grad_norm_sigma": float(grad_norm_sigma),
            "grad_r1_snr": float(_signal_to_noise(grad_r1, float(grad_r1_sigma))),
            "grad_r2_snr": float(_signal_to_noise(grad_r2, float(grad_r2_sigma))),
            "grad_norm_snr": float(_signal_to_noise(grad_norm, float(grad_norm_sigma))),
            "gradient_vector_snr": float(_vector_snr(covariance, grad_r1, grad_r2)),
            "gradient_covariance": covariance,
            "gradient_jackknife_replicas": int(len(jk_vectors)),
            "r1_minus_score": float(full_payload["r1_minus_score"]),
            "r1_plus_score": float(full_payload["r1_plus_score"]),
            "r2_minus_score": float(full_payload["r2_minus_score"]),
            "r2_plus_score": float(full_payload["r2_plus_score"]),
        }


def _propose_step(
    point: CouplingPoint,
    *,
    grad_r1: float,
    grad_r2: float,
    grad_norm: float,
    grad_norm_snr: float,
    gradient_vector_snr: float,
    gradient_covariance: np.ndarray,
    args: argparse.Namespace,
) -> dict[str, Any]:
    if grad_norm <= 0.0 or not math.isfinite(grad_norm):
        return {
            "proposal_point": CouplingPoint(float(point.r1), float(point.r2)),
            "step_r1": 0.0,
            "step_r2": 0.0,
            "step_norm": 0.0,
            "step_scale": 1.0,
            "step_scale_snr": float("nan"),
        }

    descent = _descent_vector(
        grad_r1=float(grad_r1),
        grad_r2=float(grad_r2),
        gradient_covariance=gradient_covariance,
        args=args,
    )
    descent_norm = float(np.linalg.norm(descent))
    if not math.isfinite(descent_norm) or descent_norm <= 0.0:
        return {
            "proposal_point": CouplingPoint(float(point.r1), float(point.r2)),
            "step_r1": 0.0,
            "step_r2": 0.0,
            "step_norm": 0.0,
            "step_scale": 1.0,
            "step_scale_snr": float("nan"),
        }

    if str(args.step_mode) == "normalized":
        base_step = -float(args.step_size) * descent / descent_norm
    else:
        base_step = -float(args.step_size) * descent

    step_scale, step_scale_snr = _adaptive_step_scale(
        grad_norm_snr=float(grad_norm_snr),
        gradient_vector_snr=float(gradient_vector_snr),
        args=args,
    )
    step = float(step_scale) * base_step
    unclipped_point = CouplingPoint(float(point.r1) + float(step[0]), float(point.r2) + float(step[1]))
    proposal_point = _clip_point(unclipped_point, args)
    actual_step_r1 = float(proposal_point.r1) - float(point.r1)
    actual_step_r2 = float(proposal_point.r2) - float(point.r2)
    return {
        "proposal_point": proposal_point,
        "step_r1": float(actual_step_r1),
        "step_r2": float(actual_step_r2),
        "step_norm": float(math.hypot(actual_step_r1, actual_step_r2)),
        "step_scale": float(step_scale),
        "step_scale_snr": float(step_scale_snr),
    }


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "step_index",
        "current_r1",
        "current_r2",
        "current_score",
        "grad_r1",
        "grad_r1_sigma",
        "grad_r1_snr",
        "grad_r2",
        "grad_r2_sigma",
        "grad_r2_snr",
        "grad_norm",
        "grad_norm_sigma",
        "grad_norm_snr",
        "gradient_vector_snr",
        "gradient_cov_11",
        "gradient_cov_12",
        "gradient_cov_22",
        "gradient_cov_det",
        "gradient_cov_corr",
        "gradient_jackknife_replicas",
        "step_adaptation",
        "step_scale",
        "step_scale_snr",
        "proposal_step_r1",
        "proposal_step_r2",
        "proposal_step_norm",
        "proposal_r1",
        "proposal_r2",
        "predicted_score_from_current_donor",
        "direct_score_at_proposal",
        "accepted",
        "stop_reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(header, "") for header in headers])


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _covariance_correlation(cov_11: float, cov_12: float, cov_22: float) -> float:
    if not math.isfinite(cov_11) or not math.isfinite(cov_12) or not math.isfinite(cov_22):
        return float("nan")
    if cov_11 <= 0.0 or cov_22 <= 0.0:
        return float("nan")
    denom = math.sqrt(cov_11 * cov_22)
    if denom <= 0.0:
        return float("nan")
    return float(cov_12 / denom)


def _plot_trajectory(path: Path, rows: list[dict[str, Any]], *, target_point: CouplingPoint, metric_key: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.6), constrained_layout=True)
    trajectory_ax = axes[0, 0]
    score_ax = axes[0, 1]
    covariance_ax = axes[1, 0]
    latest_cov_ax = axes[1, 1]

    accepted_points: list[tuple[float, float]] = []
    if rows:
        accepted_points.append((float(rows[0]["current_r1"]), float(rows[0]["current_r2"])))
    for row in rows:
        if str(row.get("accepted", "")).lower() == "true":
            accepted_points.append((float(row["proposal_r1"]), float(row["proposal_r2"])))

    rejected_point: tuple[float, float] | None = None
    if rows and str(rows[-1].get("accepted", "")).lower() != "true":
        rejected_point = (float(rows[-1]["proposal_r1"]), float(rows[-1]["proposal_r2"]))

    if accepted_points:
        xs = [point[0] for point in accepted_points]
        ys = [point[1] for point in accepted_points]
        trajectory_ax.plot(xs, ys, "-o", color="#2563eb", linewidth=1.8, markersize=4)
        for index, (x_value, y_value) in enumerate(accepted_points):
            trajectory_ax.text(x_value, y_value, f" {index}", fontsize=8, ha="left", va="bottom")
    if rejected_point is not None and accepted_points:
        last_point = accepted_points[-1]
        trajectory_ax.plot(
            [last_point[0], rejected_point[0]],
            [last_point[1], rejected_point[1]],
            linestyle=(0, (4, 3)),
            color="#dc2626",
            linewidth=1.5,
            alpha=0.85,
        )
        trajectory_ax.scatter(
            [rejected_point[0]],
            [rejected_point[1]],
            marker="D",
            s=54,
            facecolors="#fee2e2",
            edgecolors="#dc2626",
            linewidths=1.0,
        )
    trajectory_ax.scatter(
        [float(target_point.r1)],
        [float(target_point.r2)],
        marker="x",
        s=80,
        c="#dc2626",
        linewidths=1.8,
    )
    trajectory_ax.set_title("parameter-space trajectory")
    trajectory_ax.set_xlabel("r1")
    trajectory_ax.set_ylabel("r2")
    trajectory_ax.grid(True, color="#d1d5db", linewidth=0.8, alpha=0.8)

    step_ids = [int(row["step_index"]) for row in rows]
    current_scores = [_safe_float(row["current_score"]) for row in rows]
    predicted_scores = [_safe_float(row["predicted_score_from_current_donor"]) for row in rows]
    direct_scores = [_safe_float(row["direct_score_at_proposal"]) for row in rows]
    score_ax.plot(step_ids, current_scores, "-o", color="#111827", label="current direct score")
    score_ax.plot(step_ids, predicted_scores, "--s", color="#ea580c", label="predicted proposal score")
    score_ax.plot(step_ids, direct_scores, ":^", color="#2563eb", label="direct proposal score")
    score_ax.set_title(f"{metric_key} by step")
    score_ax.set_xlabel("step index")
    score_ax.set_ylabel(metric_key)
    score_ax.grid(True, color="#d1d5db", linewidth=0.8, alpha=0.8)
    score_ax.legend(fontsize=8)

    cov_11_values = [_safe_float(row.get("gradient_cov_11", float("nan"))) for row in rows]
    cov_12_values = [_safe_float(row.get("gradient_cov_12", float("nan"))) for row in rows]
    cov_22_values = [_safe_float(row.get("gradient_cov_22", float("nan"))) for row in rows]
    covariance_ax.axhline(0.0, color="#94a3b8", linewidth=0.9, alpha=0.8)
    covariance_ax.plot(step_ids, cov_11_values, "-o", color="#7c3aed", label="var(g_r1)")
    covariance_ax.plot(step_ids, cov_22_values, "-s", color="#059669", label="var(g_r2)")
    covariance_ax.plot(step_ids, cov_12_values, "-^", color="#dc2626", label="cov(g_r1,g_r2)")
    covariance_ax.set_title("gradient covariance by step")
    covariance_ax.set_xlabel("step index")
    covariance_ax.set_ylabel("covariance")
    covariance_ax.grid(True, color="#d1d5db", linewidth=0.8, alpha=0.8)
    covariance_ax.legend(fontsize=8)

    latest_cov_ax.set_title("latest gradient covariance")
    latest_cov_ax.set_xticks([0, 1], labels=["g_r1", "g_r2"])
    latest_cov_ax.set_yticks([0, 1], labels=["g_r1", "g_r2"])
    if rows:
        latest_cov_11 = cov_11_values[-1]
        latest_cov_12 = cov_12_values[-1]
        latest_cov_22 = cov_22_values[-1]
        latest_matrix = np.asarray(
            [[latest_cov_11, latest_cov_12], [latest_cov_12, latest_cov_22]],
            dtype=float,
        )
        scale = float(np.nanmax(np.abs(latest_matrix))) if np.isfinite(latest_matrix).any() else 1.0
        if not math.isfinite(scale) or scale <= 0.0:
            scale = 1.0
        latest_cov_ax.imshow(latest_matrix, cmap="coolwarm", vmin=-scale, vmax=scale)
        for row_index in range(2):
            for col_index in range(2):
                latest_cov_ax.text(
                    col_index,
                    row_index,
                    f"{latest_matrix[row_index, col_index]:.2e}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="#0f172a",
                    bbox={"boxstyle": "round,pad=0.15", "fc": "#ffffff", "ec": "none", "alpha": 0.80},
                )
        latest_corr = _safe_float(rows[-1].get("gradient_cov_corr", float("nan")))
        latest_cov_ax.text(
            0.50,
            -0.22,
            f"corr={latest_corr:.3f}" if math.isfinite(latest_corr) else "corr=nan",
            transform=latest_cov_ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
            color="#334155",
        )
    else:
        latest_cov_ax.text(0.5, 0.5, "no covariance yet", ha="center", va="center", fontsize=10, color="#64748b")

    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if int(args.num_steps) < 1:
        raise ValueError("num_steps must be at least 1")
    if not math.isfinite(float(args.step_size)) or float(args.step_size) <= 0.0:
        raise ValueError(f"step_size must be positive, got {args.step_size}")
    if not math.isfinite(float(args.adaptive_step_snr_reference)) or float(args.adaptive_step_snr_reference) <= 0.0:
        raise ValueError(
            f"adaptive_step_snr_reference must be positive, got {args.adaptive_step_snr_reference}"
        )
    if not math.isfinite(float(args.adaptive_step_scale_min)) or float(args.adaptive_step_scale_min) <= 0.0:
        raise ValueError(f"adaptive_step_scale_min must be positive, got {args.adaptive_step_scale_min}")
    if not math.isfinite(float(args.adaptive_step_scale_max)) or float(args.adaptive_step_scale_max) <= 0.0:
        raise ValueError(f"adaptive_step_scale_max must be positive, got {args.adaptive_step_scale_max}")
    if float(args.adaptive_step_scale_min) > float(args.adaptive_step_scale_max):
        raise ValueError(
            "adaptive_step_scale_min must be <= adaptive_step_scale_max, got "
            f"{args.adaptive_step_scale_min} > {args.adaptive_step_scale_max}"
        )
    args.accept_rule = "gradient_only"

    runner = GradientFlowRunner(args)
    current_point = _clip_point(CouplingPoint(float(args.start_r1), float(args.start_r2)), args)
    trajectory_rows: list[dict[str, Any]] = []
    final_reason = "completed_requested_steps"
    output_tsv = runner.output_root / "trajectory.tsv"
    output_png = runner.output_root / "trajectory.png"

    for step_index in range(int(args.num_steps)):
        gradient_payload = runner.gradient(current_point, step_index=step_index)
        current_row = dict(gradient_payload["center_row"])
        grad_r1 = float(gradient_payload["grad_r1"])
        grad_r2 = float(gradient_payload["grad_r2"])
        grad_norm = float(gradient_payload["grad_norm"])
        grad_r1_sigma = float(gradient_payload["grad_r1_sigma"])
        grad_r2_sigma = float(gradient_payload["grad_r2_sigma"])
        grad_norm_sigma = float(gradient_payload["grad_norm_sigma"])
        grad_norm_snr = float(gradient_payload["grad_norm_snr"])
        gradient_vector_snr = float(gradient_payload["gradient_vector_snr"])
        gradient_covariance = np.asarray(gradient_payload["gradient_covariance"], dtype=float)
        gradient_jackknife_replicas = int(gradient_payload["gradient_jackknife_replicas"])
        covariance_11 = float(gradient_covariance[0, 0])
        covariance_12 = float(gradient_covariance[0, 1])
        covariance_22 = float(gradient_covariance[1, 1])
        covariance_det = float((covariance_11 * covariance_22) - (covariance_12 * covariance_12))
        covariance_corr = _covariance_correlation(covariance_11, covariance_12, covariance_22)

        proposal_payload = _propose_step(
            current_point,
            grad_r1=grad_r1,
            grad_r2=grad_r2,
            grad_norm=grad_norm,
            grad_norm_snr=grad_norm_snr,
            gradient_vector_snr=gradient_vector_snr,
            gradient_covariance=gradient_covariance,
            args=args,
        )
        proposal_point = proposal_payload["proposal_point"]

        if float(args.min_gradient_snr) > 0.0 and (
            not math.isfinite(grad_norm_snr) or grad_norm_snr < float(args.min_gradient_snr)
        ):
            predicted_row = dict(current_row)
            if proposal_point.tag != current_point.tag:
                predicted_row = runner.score_row(current_point, proposal_point, step_index=step_index)
            direct_proposal_value = float("nan")
            accepted = False
            stop_reason = "gradient_snr_below_threshold"
        elif proposal_point.tag == current_point.tag:
            predicted_row = dict(current_row)
            direct_proposal_value = float(current_row[runner.metric_key])
            accepted = False
            stop_reason = "zero_or_clipped_step"
        else:
            predicted_row = runner.score_row(current_point, proposal_point, step_index=step_index)
            direct_proposal_row = runner.score_row(proposal_point, proposal_point, step_index=step_index + 1)
            direct_proposal_value = float(direct_proposal_row[runner.metric_key])
            accepted = True
            stop_reason = "accepted"

        row = {
            "step_index": int(step_index),
            "current_r1": f"{float(current_point.r1):.10f}",
            "current_r2": f"{float(current_point.r2):.10f}",
            "current_score": _format_float(float(current_row[runner.metric_key]), ".10e"),
            "grad_r1": _format_float(float(grad_r1), ".10e"),
            "grad_r1_sigma": _format_float(float(grad_r1_sigma), ".10e"),
            "grad_r1_snr": _format_float(float(gradient_payload["grad_r1_snr"]), ".10e"),
            "grad_r2": _format_float(float(grad_r2), ".10e"),
            "grad_r2_sigma": _format_float(float(grad_r2_sigma), ".10e"),
            "grad_r2_snr": _format_float(float(gradient_payload["grad_r2_snr"]), ".10e"),
            "grad_norm": _format_float(float(grad_norm), ".10e"),
            "grad_norm_sigma": _format_float(float(grad_norm_sigma), ".10e"),
            "grad_norm_snr": _format_float(float(grad_norm_snr), ".10e"),
            "gradient_vector_snr": _format_float(float(gradient_vector_snr), ".10e"),
            "gradient_cov_11": _format_float(float(covariance_11), ".10e"),
            "gradient_cov_12": _format_float(float(covariance_12), ".10e"),
            "gradient_cov_22": _format_float(float(covariance_22), ".10e"),
            "gradient_cov_det": _format_float(float(covariance_det), ".10e"),
            "gradient_cov_corr": _format_float(float(covariance_corr), ".10e"),
            "gradient_jackknife_replicas": int(gradient_jackknife_replicas),
            "step_adaptation": str(args.step_adaptation),
            "step_scale": _format_float(float(proposal_payload["step_scale"]), ".10e"),
            "step_scale_snr": _format_float(float(proposal_payload["step_scale_snr"]), ".10e"),
            "proposal_step_r1": _format_float(float(proposal_payload["step_r1"]), ".10e"),
            "proposal_step_r2": _format_float(float(proposal_payload["step_r2"]), ".10e"),
            "proposal_step_norm": _format_float(float(proposal_payload["step_norm"]), ".10e"),
            "proposal_r1": f"{float(proposal_point.r1):.10f}",
            "proposal_r2": f"{float(proposal_point.r2):.10f}",
            "predicted_score_from_current_donor": _format_float(float(predicted_row[runner.metric_key]), ".10e"),
            "direct_score_at_proposal": _format_float(float(direct_proposal_value), ".10e"),
            "accepted": str(bool(accepted)),
            "stop_reason": stop_reason,
        }
        trajectory_rows.append(row)
        _write_tsv(output_tsv, trajectory_rows)
        _plot_trajectory(output_png, trajectory_rows, target_point=runner.target_point, metric_key=runner.metric_key)

        if not accepted:
            final_reason = stop_reason
            break
        current_point = proposal_point

    _plot_trajectory(output_png, trajectory_rows, target_point=runner.target_point, metric_key=runner.metric_key)

    final_point = current_point
    final_row = runner.score_row(final_point, final_point, step_index=int(args.num_steps) + 1)
    summary = {
        "args": vars(args),
        "metric_key": runner.metric_key,
        "gradient_step": float(runner.gradient_step),
        "fit_sizes": [int(size) for size in runner.fit_sizes],
        "target_point": _point_to_json(runner.target_point),
        "start_point": {"r1": float(args.start_r1), "r2": float(args.start_r2)},
        "final_point": _point_to_json(final_point),
        "final_score": float(final_row[runner.metric_key]),
        "trajectory_tsv": str(output_tsv),
        "trajectory_png": str(output_png),
        "n_attempted_steps": len(trajectory_rows),
        "n_accepted_steps": sum(1 for row in trajectory_rows if str(row["accepted"]).lower() == "true"),
        "final_reason": final_reason,
    }
    (runner.output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()