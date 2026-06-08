#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from scipy.optimize import curve_fit

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent
REPO_ROOT = RESPONSIBLE_ROOT.parent
KFC_ROOT = REPO_ROOT / "K_from_continuum"
if str(KFC_ROOT) not in sys.path:
    sys.path.insert(0, str(KFC_ROOT))

from workflow_common import ensure_dir, ensure_simulator, exact_triangular_ising_beta  # noqa: E402


DEFAULT_EXECUTION = {
    "exe": None,
    "auto_build": True,
    "force_rebuild": False,
    "build_command": ["make"],
    "build_timeout_s": 600,
}

DEFAULT_GRID_VALUES = [0.99, 0.995, 1.0, 1.005, 1.01]
DEFAULT_DONOR_VALUES = [0.99, 1.0, 1.01]
DEFAULT_FIT_SIZES = [8, 12, 16, 24]


@dataclass(frozen=True)
class CouplingPoint:
    r1: float
    r2: float

    @property
    def tag(self) -> str:
        return f"r1_{_token(self.r1)}__r2_{_token(self.r2)}"


def _token(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def _run_id(lx: int, ly: int, tx: int, ty: int, k1: float, k2: float, k3: float, kt: float = 0.5) -> str:
    return (
        f"{lx}x{ly}_t{tx}x{ty}_"
        f"k{int(round(k1 * 1000))}_{int(round(k2 * 1000))}_{int(round(k3 * 1000))}_{int(round(kt * 1000))}"
    )


def _evaluate_power_model_on_x(x_values: np.ndarray, A_value: float, B_value: float, omega_value: float) -> np.ndarray:
    return A_value + B_value * np.power(x_values, omega_value)


def _fit_blind_power_model(x_values: np.ndarray, y_values: np.ndarray, sigma_values: np.ndarray) -> dict[str, Any]:
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
        sigma_A = float(np.sqrt(max(float(pcov[0, 0]), 0.0))) if np.isfinite(pcov[0, 0]) else float("nan")
        sigma_omega = float(np.sqrt(max(float(pcov[2, 2]), 0.0))) if np.isfinite(pcov[2, 2]) else float("nan")
        fit_mode = "power_fit"
    except Exception:
        A_value = float(y_array[0])
        sigma_A = float(sigma_fit[0]) if sigma_fit.size else float("nan")
        B_value = 0.0
        omega_value = 1.0
        sigma_omega = float("nan")
        fit_mode = "fit_failed"
        pcov = np.full((3, 3), np.nan, dtype=float)

    return {
        "A": float(A_value),
        "sigma_A": float(sigma_A),
        "B": float(B_value),
        "omega": float(omega_value),
        "sigma_omega": float(sigma_omega),
        "fit_mode": str(fit_mode),
        "pcov": np.asarray(pcov, dtype=float),
    }


def _predict_at_target_x(fit_payload: dict[str, Any], target_x: float) -> tuple[float, float]:
    pred_value = float(
        _evaluate_power_model_on_x(
            np.asarray([float(target_x)], dtype=float),
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
        grad = np.asarray([1.0, x_pow, float(fit_payload["B"]) * x_pow * math.log(target_x)], dtype=float)
        pred_var = float(grad @ pcov @ grad)
        if pred_var >= 0.0:
            pred_sigma = float(math.sqrt(pred_var))
    return pred_value, pred_sigma


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a small direct grid and compare the direct geometry-match proxy surface to a nearest-donor reweighted interpolation."
    )
    parser.add_argument("--grid-values", nargs="+", type=float, default=DEFAULT_GRID_VALUES)
    parser.add_argument("--donor-values", nargs="+", type=float, default=DEFAULT_DONOR_VALUES)
    parser.add_argument("--fit-sizes", nargs="+", type=int, default=DEFAULT_FIT_SIZES)
    parser.add_argument("--target-size", type=int, default=32)
    parser.add_argument("--target-r1", type=float, default=1.0)
    parser.add_argument("--target-r2", type=float, default=1.0)
    parser.add_argument("--n-traj", type=int, default=20000)
    parser.add_argument("--n-therm", type=int, default=2000)
    parser.add_argument("--n-skip", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=2026060417)
    parser.add_argument(
        "--output-root",
        default=str(
            RESPONSIBLE_ROOT
            / "results"
            / "geometry_match_reweight_interp_iso111_grid5x5_donor3x3_sizes8_12_16_24_target32_20260604"
        ),
    )
    return parser.parse_args()


def _unique_sorted(values: list[float]) -> list[float]:
    return sorted({round(float(value), 12) for value in values})


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


def _selected_specs_for_size(size: int) -> list[tuple[str, tuple[int, int]]]:
    if size % 4 != 0:
        raise ValueError(f"size must be divisible by 4 for quarter observables, got {size}")
    half = size // 2
    quarter = size // 4
    return [
        ("mid_v", (half, 0)),
        ("mid_u", (0, half)),
        ("mid_w", (half, half)),
        ("q_v", (quarter, 0)),
        ("q_u", (0, quarter)),
        ("q_w", (quarter, quarter)),
        ("anchor", (0, -1)),
    ]


def _panel_specs() -> list[dict[str, str]]:
    return [
        {"label": "mid_v_anchor", "numerator": "mid_v", "denominator": "anchor"},
        {"label": "mid_u_anchor", "numerator": "mid_u", "denominator": "anchor"},
        {"label": "mid_w_anchor", "numerator": "mid_w", "denominator": "anchor"},
        {"label": "mid_v_over_u", "numerator": "mid_v", "denominator": "mid_u"},
        {"label": "mid_w_over_u", "numerator": "mid_w", "denominator": "mid_u"},
        {"label": "q_v_anchor", "numerator": "q_v", "denominator": "anchor"},
        {"label": "q_u_anchor", "numerator": "q_u", "denominator": "anchor"},
        {"label": "q_w_anchor", "numerator": "q_w", "denominator": "anchor"},
        {"label": "q_v_over_u", "numerator": "q_v", "denominator": "q_u"},
        {"label": "q_w_over_u", "numerator": "q_w", "denominator": "q_u"},
    ]


def _point_beta(point: CouplingPoint) -> float:
    return float(exact_triangular_ising_beta({"k1": point.r1, "k2": point.r2, "k3": 1.0}))


def _selected_disp_arg(specs: list[tuple[str, tuple[int, int]]]) -> str:
    return ",".join(f"{int(coord[0])}:{int(coord[1])}" for _, coord in specs)


def _parse_selected_bundle(path: Path, labels: list[str]) -> dict[str, Any]:
    header: dict[str, Any] = {}
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split()
                if not parts:
                    continue
                key = parts[0]
                if key in {"seed", "L_x", "L_y", "T_x", "T_y", "n_therm", "n_traj", "n_skip", "n_selected"}:
                    header[key] = int(parts[1])
                elif key in {"K1", "K2", "K3", "Kt", "beta"}:
                    header[key] = float(parts[1])
                continue
            cols = line.split()
            expected_cols = 1 + len(labels) + 4
            if len(cols) != expected_cols:
                raise ValueError(f"expected {expected_cols} columns in {path}, got {len(cols)}")
            rows.append([float(value) for value in cols[1:]])
    if not rows:
        raise ValueError(f"no rows found in {path}")
    data = np.asarray(rows, dtype=float)
    corr: dict[str, np.ndarray] = {}
    for idx, label in enumerate(labels):
        corr[str(label)] = np.asarray(data[:, idx], dtype=float)
    offset = len(labels)
    return {
        "path": str(path),
        "corr": corr,
        "mag": np.asarray(data[:, offset], dtype=float),
        "e1": np.asarray(data[:, offset + 1], dtype=float),
        "e2": np.asarray(data[:, offset + 2], dtype=float),
        "e3": np.asarray(data[:, offset + 3], dtype=float),
        "beta": float(header["beta"]),
        "K": (float(header["K1"]), float(header["K2"]), float(header["K3"])),
        "n_sites": int(header["L_x"]) * int(header["L_y"]),
        "header": header,
    }


def _block_slices(n_samples: int, n_blocks: int) -> list[tuple[int, int]]:
    block_count = max(2, min(int(n_blocks), max(n_samples // 4, 2)))
    block = n_samples // block_count
    slices: list[tuple[int, int]] = []
    for index in range(block_count):
        lo = index * block
        hi = (index + 1) * block if index < block_count - 1 else n_samples
        if hi > lo:
            slices.append((lo, hi))
    return slices


def _normalize(log_w: np.ndarray) -> np.ndarray:
    shifted = log_w - np.max(log_w)
    w = np.exp(shifted)
    total = np.sum(w)
    if not np.isfinite(total) or total <= 0.0:
        return np.zeros_like(w)
    return w / total


def _connected_from_weights(w: np.ndarray, corr: np.ndarray, mag: np.ndarray) -> float:
    corr_mean = float(np.sum(w * corr))
    mag_mean = float(np.sum(w * mag))
    return corr_mean - mag_mean * mag_mean


def _n_eff(w: np.ndarray) -> float:
    total = float(np.sum(w))
    if total <= 0.0:
        return 0.0
    return total * total / max(float(np.dot(w, w)), 1.0e-300)


def _log_weights(payload: dict[str, Any], target_point: CouplingPoint) -> np.ndarray:
    beta_target = _point_beta(target_point)
    k_target = (float(target_point.r1), float(target_point.r2), 1.0)
    k_parent = tuple(float(value) for value in payload["K"])
    parent = k_parent[0] * payload["e1"] + k_parent[1] * payload["e2"] + k_parent[2] * payload["e3"]
    target = k_target[0] * payload["e1"] + k_target[1] * payload["e2"] + k_target[2] * payload["e3"]
    return -float(payload["n_sites"]) * (float(beta_target) * target - float(payload["beta"]) * parent)


def _ratio_stat(
    corr_num: np.ndarray,
    corr_den: np.ndarray,
    mag: np.ndarray,
    *,
    log_w: np.ndarray | None = None,
    n_blocks: int = 16,
) -> dict[str, float]:
    n_samples = int(mag.size)
    if n_samples < 16:
        raise ValueError("need at least 16 measurements for the interpolation test")
    base_w = np.full(n_samples, 1.0 / float(n_samples), dtype=float) if log_w is None else _normalize(log_w)
    num_conn = _connected_from_weights(base_w, corr_num, mag)
    den_conn = _connected_from_weights(base_w, corr_den, mag)
    if abs(den_conn) < 1.0e-12:
        raise ZeroDivisionError("denominator connected correlator is too small for a stable ratio")
    value = float(num_conn / den_conn)

    slices = _block_slices(n_samples, n_blocks)
    leave_one_out: list[float] = []
    for lo, hi in slices:
        mask = np.ones(n_samples, dtype=bool)
        mask[lo:hi] = False
        if np.count_nonzero(mask) < 8:
            continue
        sub_w = np.full(np.count_nonzero(mask), 1.0 / float(np.count_nonzero(mask)), dtype=float)
        if log_w is not None:
            sub_w = _normalize(log_w[mask])
        num_leave = _connected_from_weights(sub_w, corr_num[mask], mag[mask])
        den_leave = _connected_from_weights(sub_w, corr_den[mask], mag[mask])
        if abs(den_leave) < 1.0e-12:
            continue
        leave_one_out.append(float(num_leave / den_leave))
    if len(leave_one_out) < 2:
        sigma = float("nan")
    else:
        jk = np.asarray(leave_one_out, dtype=float)
        sigma = float(np.sqrt((jk.size - 1.0) / jk.size * np.sum((jk - np.mean(jk)) ** 2)))

    n_eff = float(n_samples) if log_w is None else _n_eff(base_w * n_samples)
    return {"value": value, "sigma": sigma, "n_eff_fraction": float(n_eff) / float(n_samples)}


def _panel_values(payload: dict[str, Any], *, target_point: CouplingPoint | None = None) -> dict[str, dict[str, float]]:
    log_w = _log_weights(payload, target_point) if target_point is not None else None
    corr = payload["corr"]
    rows: dict[str, dict[str, float]] = {}
    for panel in _panel_specs():
        rows[str(panel["label"])] = _ratio_stat(
            np.asarray(corr[str(panel["numerator"])], dtype=float),
            np.asarray(corr[str(panel["denominator"])], dtype=float),
            np.asarray(payload["mag"], dtype=float),
            log_w=log_w,
        )
    return rows


def _nearest_donor(point: CouplingPoint, donors: list[CouplingPoint]) -> CouplingPoint:
    return min(donors, key=lambda donor: (donor.r1 - point.r1) ** 2 + (donor.r2 - point.r2) ** 2)


def _score_family(
    family_rows: dict[int, dict[str, dict[str, float]]],
    *,
    fit_sizes: list[int],
    target_x: float,
    target_rows: dict[str, dict[str, float]],
) -> dict[str, Any]:
    panel_z: dict[str, float] = {}
    panel_neff: dict[str, float] = {}
    total = 0.0
    for panel in _panel_specs():
        label = str(panel["label"])
        x_values = np.asarray([1.0 / float(size) for size in fit_sizes], dtype=float)
        y_values = np.asarray([float(family_rows[size][label]["value"]) for size in fit_sizes], dtype=float)
        sigma_values = np.asarray([float(family_rows[size][label]["sigma"]) for size in fit_sizes], dtype=float)
        fit_payload = _fit_blind_power_model(x_values, y_values, sigma_values)
        pred_value, pred_sigma = _predict_at_target_x(fit_payload, target_x)
        target_value = float(target_rows[label]["value"])
        target_sigma = float(target_rows[label]["sigma"])
        denom = math.sqrt(
            max(pred_sigma ** 2 if np.isfinite(pred_sigma) else 0.0, 0.0) + max(target_sigma ** 2, 0.0)
        )
        z_value = (pred_value - target_value) / denom if denom > 0.0 else float("nan")
        panel_z[label] = float(z_value)
        panel_neff[label] = float(min(family_rows[size][label]["n_eff_fraction"] for size in fit_sizes))
        total += float(z_value ** 2)
    return {"z2_sum": float(total), "panel_z": panel_z, "panel_min_neff_fraction": panel_neff}


def _run_point(
    exe: str,
    output_root: Path,
    *,
    size: int,
    point: CouplingPoint,
    n_traj: int,
    n_therm: int,
    n_skip: int,
    seed: int,
) -> Path:
    size_dir = output_root / "raw" / f"L{size}" / point.tag
    ensure_dir(str(size_dir))
    specs = _selected_specs_for_size(size)
    bundle_name = "selected_observables_bundle.dat"
    run_dir = size_dir / _run_id(size, size, 0, 0, point.r1, point.r2, 1.0)
    bundle_path = run_dir / bundle_name
    if bundle_path.is_file():
        return bundle_path

    beta = _point_beta(point)
    cmd = [
        exe,
        "--L_x",
        str(size),
        "--L_y",
        str(size),
        "--T_x",
        "0",
        "--T_y",
        "0",
        "--k1",
        f"{point.r1:.12f}",
        "--k2",
        f"{point.r2:.12f}",
        "--k3",
        "1.000000000000",
        "--beta",
        f"{beta:.12f}",
        "--n_traj",
        str(int(n_traj)),
        "--n_therm",
        str(int(n_therm)),
        "--n_skip",
        str(int(n_skip)),
        "--seed",
        str(int(seed)),
        "--data_dir",
        str(size_dir),
        "--selected_disp_list",
        _selected_disp_arg(specs),
        "--selected_disp_bundle_name",
        bundle_name,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"simulator failed for size={size} point={point.tag}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
    if not bundle_path.is_file():
        raise FileNotFoundError(f"missing selected bundle after successful run: {bundle_path}")
    return bundle_path


def _write_summary(output_root: Path, rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    headers = [
        "r1",
        "r2",
        "donor_r1",
        "donor_r2",
        "direct_z2_sum",
        "reweighted_z2_sum",
        "delta_z2_sum",
    ]
    with (output_root / "grid_scores.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        f"{float(row['r1']):.6f}",
                        f"{float(row['r2']):.6f}",
                        f"{float(row['donor_r1']):.6f}",
                        f"{float(row['donor_r2']):.6f}",
                        f"{float(row['direct_z2_sum']):.16e}",
                        f"{float(row['reweighted_z2_sum']):.16e}",
                        f"{float(row['delta_z2_sum']):.16e}",
                    ]
                )
                + "\n"
            )


def _save_plot(output_root: Path, rows: list[dict[str, Any]], donors: list[CouplingPoint], target: CouplingPoint) -> None:
    r1_values = sorted({float(row["r1"]) for row in rows})
    r2_values = sorted({float(row["r2"]) for row in rows})
    index_r1 = {value: idx for idx, value in enumerate(r1_values)}
    index_r2 = {value: idx for idx, value in enumerate(r2_values)}
    direct_grid = np.full((len(r2_values), len(r1_values)), np.nan, dtype=float)
    reweighted_grid = np.full_like(direct_grid, np.nan)
    delta_grid = np.full_like(direct_grid, np.nan)
    donor_grid = np.full_like(direct_grid, np.nan)
    donor_index = {donor.tag: idx for idx, donor in enumerate(donors)}

    for row in rows:
        iy = index_r2[float(row["r2"])]
        ix = index_r1[float(row["r1"])]
        direct_grid[iy, ix] = float(row["direct_z2_sum"])
        reweighted_grid[iy, ix] = float(row["reweighted_z2_sum"])
        delta_grid[iy, ix] = float(row["delta_z2_sum"])
        donor_grid[iy, ix] = float(donor_index[f"r1_{_token(float(row['donor_r1']))}__r2_{_token(float(row['donor_r2']))}"])

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0), constrained_layout=True)
    extent = [min(r1_values), max(r1_values), min(r2_values), max(r2_values)]

    im0 = axes[0, 0].imshow(direct_grid, origin="lower", extent=extent, aspect="auto", cmap="viridis")
    axes[0, 0].set_title("Direct holdout z2_sum")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.9)

    im1 = axes[0, 1].imshow(reweighted_grid, origin="lower", extent=extent, aspect="auto", cmap="viridis")
    axes[0, 1].set_title("Nearest-donor reweighted z2_sum")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9)

    vmax = float(np.nanmax(np.abs(delta_grid)))
    im2 = axes[1, 0].imshow(delta_grid, origin="lower", extent=extent, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    axes[1, 0].set_title("Reweighted - direct")
    fig.colorbar(im2, ax=axes[1, 0], shrink=0.9)

    im3 = axes[1, 1].imshow(donor_grid, origin="lower", extent=extent, aspect="auto", cmap="tab10")
    axes[1, 1].set_title("Nearest donor patch")
    fig.colorbar(im3, ax=axes[1, 1], shrink=0.9)

    for axis in axes.ravel():
        axis.scatter([donor.r1 for donor in donors], [donor.r2 for donor in donors], marker="s", s=55, c="white", edgecolors="black")
        axis.scatter([target.r1], [target.r2], marker="*", s=120, c="#ef4444", edgecolors="black")
        axis.set_xlabel("r1")
        axis.set_ylabel("r2")

    fig.suptitle("Local geometry-match proxy grid: direct vs reweighted interpolation", fontsize=13)
    fig.savefig(output_root / "geometry_match_grid_interpolation.png", dpi=180)
    plt.close(fig)


def _panel_fit_details(
    family_rows: dict[int, dict[str, dict[str, float]]],
    *,
    fit_sizes: list[int],
    target_x: float,
    target_rows: dict[str, dict[str, float]],
) -> dict[str, dict[str, Any]]:
    details: dict[str, dict[str, Any]] = {}
    x_values = np.asarray([1.0 / math.sqrt(float(size) * float(size)) for size in fit_sizes], dtype=float)
    for panel in _panel_specs():
        label = str(panel["label"])
        y_values = np.asarray([float(family_rows[size][label]["value"]) for size in fit_sizes], dtype=float)
        sigma_values = np.asarray([float(family_rows[size][label]["sigma"]) for size in fit_sizes], dtype=float)
        fit_payload = _fit_blind_power_model(x_values, y_values, sigma_values)
        pred_value, pred_sigma = _predict_at_target_x(fit_payload, target_x)
        details[label] = {
            "x_values": [float(value) for value in x_values],
            "y_values": [float(value) for value in y_values],
            "sigma_values": [float(value) for value in sigma_values],
            "fit": {
                "A": float(fit_payload["A"]),
                "B": float(fit_payload["B"]),
                "omega": float(fit_payload["omega"]),
                "sigma_A": float(fit_payload["sigma_A"]),
                "sigma_omega": float(fit_payload["sigma_omega"]),
                "fit_mode": str(fit_payload["fit_mode"]),
            },
            "fit_payload": fit_payload,
            "pred": {"value": float(pred_value), "sigma": float(pred_sigma)},
            "target": {
                "value": float(target_rows[label]["value"]),
                "sigma": float(target_rows[label]["sigma"]),
            },
        }
    return details


def _predict_band(fit_payload: dict[str, Any], *, x_values: np.ndarray) -> np.ndarray:
    x_array = np.asarray(x_values, dtype=float)
    band_sigma = np.full_like(x_array, np.nan, dtype=float)
    pcov = np.asarray(fit_payload.get("pcov"), dtype=float)
    if pcov.shape != (3, 3) or not np.all(np.isfinite(pcov)):
        return band_sigma

    omega_value = float(fit_payload["omega"])
    b_value = float(fit_payload["B"])
    x_pow = np.power(x_array, omega_value)
    grad = np.column_stack(
        [
            np.ones_like(x_array, dtype=float),
            x_pow,
            b_value * x_pow * np.log(x_array),
        ]
    )
    pred_var = np.einsum("ij,jk,ik->i", grad, pcov, grad)
    valid = np.isfinite(pred_var) & (pred_var >= 0.0)
    band_sigma[valid] = np.sqrt(pred_var[valid])
    return band_sigma


def _prediction_delta_summary(detail: dict[str, Any]) -> dict[str, float]:
    pred_value = float(detail["pred"]["value"])
    pred_sigma = float(detail["pred"]["sigma"])
    target_value = float(detail["target"]["value"])
    target_sigma = float(detail["target"]["sigma"])
    delta = pred_value - target_value
    if np.isfinite(pred_sigma) and pred_sigma > 0.0 and np.isfinite(target_sigma):
        denom = float(math.sqrt(pred_sigma**2 + target_sigma**2))
    elif np.isfinite(target_sigma) and target_sigma > 0.0:
        denom = float(target_sigma)
    else:
        denom = float("nan")
    z_value = delta / denom if np.isfinite(denom) and denom > 0.0 else float("nan")
    return {
        "pred_target": pred_value,
        "pred_target_sigma": pred_sigma,
        "delta_target": delta,
        "z_target": float(z_value),
    }


def _acute_style_panel_groups() -> list[dict[str, Any]]:
    return [
        {
            "slug": "midpoint",
            "title": "midpoint anchored-ratio and ratio FSS",
            "panels": [
                {"label": "mid_v_anchor", "title": "mid_v_anchor", "color": "#1d4ed8"},
                {"label": "mid_u_anchor", "title": "mid_u_anchor", "color": "#047857"},
                {"label": "mid_w_anchor", "title": "mid_w_anchor", "color": "#b45309"},
                {"label": "mid_v_over_u", "title": "mid_v_over_u", "color": "#7c2d12"},
                {"label": "mid_w_over_u", "title": "mid_w_over_u", "color": "#7c3aed"},
            ],
        },
        {
            "slug": "quarter",
            "title": "quarter-point anchored-ratio and ratio FSS",
            "panels": [
                {"label": "q_v_anchor", "title": "q_v_anchor", "color": "#1d4ed8"},
                {"label": "q_u_anchor", "title": "q_u_anchor", "color": "#047857"},
                {"label": "q_w_anchor", "title": "q_w_anchor", "color": "#b45309"},
                {"label": "q_v_over_u", "title": "q_v_over_u", "color": "#7c2d12"},
                {"label": "q_w_over_u", "title": "q_w_over_u", "color": "#7c3aed"},
            ],
        },
    ]


def _render_acute_style_panel(
    axis: plt.Axes,
    *,
    panel_meta: dict[str, Any],
    direct_detail: dict[str, Any],
    reweight_detail: dict[str, Any],
    target_x: float,
    target_size: int,
    direct_z: float,
    reweighted_z: float,
    reweight_neff: float,
) -> None:
    color = str(panel_meta["color"])
    direct_fit = direct_detail["fit"]
    reweight_fit = reweight_detail["fit"]
    direct_fit_payload = dict(direct_detail.get("fit_payload", direct_fit))
    reweight_fit_payload = dict(reweight_detail.get("fit_payload", reweight_fit))
    reweight_color = "#ea580c"
    direct_summary = _prediction_delta_summary(direct_detail)
    reweight_summary = _prediction_delta_summary(reweight_detail)

    direct_x = np.asarray(direct_detail["x_values"], dtype=float)
    direct_y = np.asarray(direct_detail["y_values"], dtype=float)
    direct_sigma = np.asarray(direct_detail["sigma_values"], dtype=float)
    reweight_x = np.asarray(reweight_detail["x_values"], dtype=float)
    reweight_y = np.asarray(reweight_detail["y_values"], dtype=float)
    reweight_sigma = np.asarray(reweight_detail["sigma_values"], dtype=float)

    x_min = min(float(np.min(direct_x)), float(np.min(reweight_x)), float(target_x))
    x_max = max(float(np.max(direct_x)), float(np.max(reweight_x)), float(target_x))
    x_plot_min = max(x_min * 0.9, 1.0e-6)
    x_plot_max = x_max * 1.08
    x_fit = np.geomspace(x_plot_min, x_plot_max, 300)

    direct_fit_y = _evaluate_power_model_on_x(x_fit, float(direct_fit["A"]), float(direct_fit["B"]), float(direct_fit["omega"]))
    direct_fit_sigma = _predict_band(direct_fit_payload, x_values=x_fit)
    reweight_fit_y = _evaluate_power_model_on_x(
        x_fit,
        float(reweight_fit["A"]),
        float(reweight_fit["B"]),
        float(reweight_fit["omega"]),
    )
    reweight_fit_sigma = _predict_band(reweight_fit_payload, x_values=x_fit)

    if np.any(np.isfinite(direct_fit_sigma)):
        valid = np.isfinite(direct_fit_sigma)
        axis.fill_between(
            x_fit[valid],
            (direct_fit_y - direct_fit_sigma)[valid],
            (direct_fit_y + direct_fit_sigma)[valid],
            color=color,
            alpha=0.16,
            zorder=1,
        )
    if np.any(np.isfinite(reweight_fit_sigma)):
        valid = np.isfinite(reweight_fit_sigma)
        axis.fill_between(
            x_fit[valid],
            (reweight_fit_y - reweight_fit_sigma)[valid],
            (reweight_fit_y + reweight_fit_sigma)[valid],
            color=reweight_color,
            alpha=0.10,
            zorder=1,
        )

    axis.errorbar(
        direct_x,
        direct_y,
        yerr=direct_sigma,
        fmt="o",
        color=color,
        ecolor=color,
        capsize=3,
        markersize=5,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.0,
        label="direct",
        zorder=3,
    )
    axis.errorbar(
        reweight_x,
        reweight_y,
        yerr=reweight_sigma,
        fmt="s",
        color=reweight_color,
        ecolor=reweight_color,
        capsize=3,
        markersize=5,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.0,
        linestyle="--",
        label="reweighted",
        zorder=4,
    )
    axis.plot(x_fit, direct_fit_y, color=color, linewidth=2.0, zorder=2)
    axis.plot(x_fit, reweight_fit_y, color=reweight_color, linewidth=2.0, linestyle="--", zorder=2)

    target = direct_detail["target"]
    axis.errorbar(
        [target_x],
        [float(target["value"])],
        yerr=[float(target["sigma"])],
        fmt="*",
        color="#dc2626",
        ecolor="#dc2626",
        capsize=3,
        markersize=9,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linewidth=1.0,
        label=f"target L={target_size}",
        zorder=5,
    )
    axis.set_title(
        f"{str(panel_meta['title'])}\nzd={float(direct_z):+.2f}  zrw={float(reweighted_z):+.2f}  min N_eff/N={float(reweight_neff):.3f}",
        fontsize=9,
    )
    axis.text(
        0.03,
        0.97,
        (
            f"direct A = {float(direct_fit['A']):.6f} +/- {float(direct_fit.get('sigma_A', float('nan'))):.6f}\n"
            f"direct omega = {float(direct_fit['omega']):.3f}\n"
            f"direct pred@target = {float(direct_summary['pred_target']):.6f} +/- {float(direct_summary['pred_target_sigma']):.6f}  z = {float(direct_summary['z_target']):+.2f}\n"
            f"rw A = {float(reweight_fit['A']):.6f} +/- {float(reweight_fit.get('sigma_A', float('nan'))):.6f}\n"
            f"rw omega = {float(reweight_fit['omega']):.3f}\n"
            f"rw pred@target = {float(reweight_summary['pred_target']):.6f} +/- {float(reweight_summary['pred_target_sigma']):.6f}  z = {float(reweight_summary['z_target']):+.2f}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
    )
    axis.set_xlabel("1 / sqrt(lattice volume)")
    axis.set_ylabel("Correlator ratio")
    axis.set_xscale("log")
    axis.set_xlim(x_plot_min, x_plot_max)
    axis.grid(True, which="both", alpha=0.35)
    axis.legend(loc="lower right", fontsize=7.5)


def _save_acute_style_fss_sheet(
    output_path: Path,
    *,
    title: str,
    subtitle: str,
    group_meta: dict[str, Any],
    direct_details: dict[str, dict[str, Any]],
    reweighted_details: dict[str, dict[str, Any]],
    row: dict[str, Any],
    target_x: float,
    target_size: int,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(18, 9.8), squeeze=False)
    axes_flat = list(axes.ravel())
    fig.suptitle(title, fontsize=15, y=0.98)
    fig.text(0.5, 0.945, subtitle, ha="center", va="top", fontsize=9.5, color="#444444")

    for axis, panel_meta in zip(axes_flat, list(group_meta["panels"])):
        label = str(panel_meta["label"])
        _render_acute_style_panel(
            axis,
            panel_meta=panel_meta,
            direct_detail=direct_details[label],
            reweight_detail=reweighted_details[label],
            target_x=target_x,
            target_size=target_size,
            direct_z=float(row["direct_panel_z"][label]),
            reweighted_z=float(row["reweighted_panel_z"][label]),
            reweight_neff=float(row["reweighted_panel_min_neff_fraction"][label]),
        )
    for axis in axes_flat[len(group_meta["panels"]):]:
        axis.axis("off")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.92])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _make_row_lookup(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row["point_tag"]): row for row in rows}


def _select_fss_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(row: dict[str, Any] | None) -> None:
        if row is None:
            return
        tag = str(row["point_tag"])
        if tag in seen:
            return
        seen.add(tag)
        selected.append(row)

    row_lookup = _make_row_lookup(rows)
    add(row_lookup.get("r1_1p000000__r2_1p000000"))
    add(min(rows, key=lambda row: float(row["direct_z2_sum"])))
    add(min(rows, key=lambda row: float(row["reweighted_z2_sum"])))
    add(max(rows, key=lambda row: abs(float(row["delta_z2_sum"]))))
    return selected


def _save_fss_plots(
    output_root: Path,
    *,
    rows: list[dict[str, Any]],
    direct_family_by_tag: dict[str, dict[int, dict[str, dict[str, float]]]],
    reweighted_family_by_tag: dict[str, dict[int, dict[str, dict[str, float]]]],
    donor_by_tag: dict[str, CouplingPoint],
    fit_sizes: list[int],
    target_size: int,
    target_x: float,
    target_rows: dict[str, dict[str, float]],
) -> None:
    selected_rows = _select_fss_rows(rows)
    if not selected_rows:
        return

    fss_dir = output_root / "fss_panels"
    ensure_dir(str(fss_dir))
    x_curve = np.linspace(0.0, max(1.0 / float(size) for size in fit_sizes) * 1.05, 300)
    panel_specs = _panel_specs()
    acute_style_groups = _acute_style_panel_groups()

    for row in selected_rows:
        point_tag = str(row["point_tag"])
        direct_family = direct_family_by_tag[point_tag]
        reweighted_family = reweighted_family_by_tag[point_tag]
        direct_details = _panel_fit_details(direct_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        reweighted_details = _panel_fit_details(reweighted_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        donor = donor_by_tag[point_tag]

        fig, axes = plt.subplots(5, 2, figsize=(13.5, 17.5), constrained_layout=True)
        for axis, panel in zip(axes.ravel(), panel_specs):
            label = str(panel["label"])
            direct_detail = direct_details[label]
            reweight_detail = reweighted_details[label]

            direct_x = np.asarray(direct_detail["x_values"], dtype=float)
            direct_y = np.asarray(direct_detail["y_values"], dtype=float)
            direct_sigma = np.asarray(direct_detail["sigma_values"], dtype=float)
            rw_x = np.asarray(reweight_detail["x_values"], dtype=float)
            rw_y = np.asarray(reweight_detail["y_values"], dtype=float)
            rw_sigma = np.asarray(reweight_detail["sigma_values"], dtype=float)

            axis.errorbar(direct_x, direct_y, yerr=direct_sigma, fmt="o", color="#2563eb", label="direct", capsize=3)
            axis.errorbar(rw_x, rw_y, yerr=rw_sigma, fmt="s", color="#ea580c", label="reweighted", capsize=3)

            direct_fit = direct_detail["fit"]
            rw_fit = reweight_detail["fit"]
            axis.plot(
                x_curve,
                _evaluate_power_model_on_x(x_curve, direct_fit["A"], direct_fit["B"], direct_fit["omega"]),
                color="#2563eb",
                linewidth=1.5,
            )
            axis.plot(
                x_curve,
                _evaluate_power_model_on_x(x_curve, rw_fit["A"], rw_fit["B"], rw_fit["omega"]),
                color="#ea580c",
                linewidth=1.5,
                linestyle="--",
            )

            target = direct_detail["target"]
            axis.errorbar([target_x], [target["value"]], yerr=[target["sigma"]], fmt="*", color="#dc2626", markersize=10, label=f"target L={target_size}")

            direct_z = float(row["direct_panel_z"][label])
            reweighted_z = float(row["reweighted_panel_z"][label])
            neff = float(row["reweighted_panel_min_neff_fraction"][label])
            axis.set_title(f"{label}\nzd={direct_z:+.2f}  zrw={reweighted_z:+.2f}  min N_eff/N={neff:.3f}", fontsize=9)
            axis.set_xlabel("1 / L")
            axis.set_ylabel("ratio")
            axis.grid(alpha=0.25)

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)

        point_r1 = float(row["r1"])
        point_r2 = float(row["r2"])
        fig.suptitle(
            (
                f"FSS panels for point ({point_r1:.3f}, {point_r2:.3f}) using donor ({float(donor.r1):.3f}, {float(donor.r2):.3f})\n"
                f"direct z2={float(row['direct_z2_sum']):.4f}  reweighted z2={float(row['reweighted_z2_sum']):.4f}  delta={float(row['delta_z2_sum']):+.4f}"
            ),
            fontsize=12,
        )
        fig.savefig(fss_dir / f"{point_tag}_fss.png", dpi=180)
        plt.close(fig)

        title_base = (
            f"Point ({point_r1:.3f}, {point_r2:.3f}) using donor ({float(donor.r1):.3f}, {float(donor.r2):.3f})"
        )
        subtitle_base = (
            f"direct z2={float(row['direct_z2_sum']):.4f}  reweighted z2={float(row['reweighted_z2_sum']):.4f}  delta={float(row['delta_z2_sum']):+.4f}; "
            f"untwisted sizes={','.join(str(int(size)) for size in fit_sizes)}  target L={int(target_size)}"
        )
        for group_meta in acute_style_groups:
            _save_acute_style_fss_sheet(
                fss_dir / f"{point_tag}_{str(group_meta['slug'])}_acute_style.png",
                title=f"{title_base} {str(group_meta['title'])}",
                subtitle=subtitle_base,
                group_meta=group_meta,
                direct_details=direct_details,
                reweighted_details=reweighted_details,
                row=row,
                target_x=target_x,
                target_size=target_size,
            )


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    ensure_dir(str(output_root))
    exe = ensure_simulator(DEFAULT_EXECUTION)

    grid_values = _unique_sorted(list(args.grid_values))
    donor_values = _unique_sorted(list(args.donor_values))
    execution_sizes = _unique_ints_preserve_order(list(args.fit_sizes))
    fit_sizes = sorted(set(execution_sizes))
    target_size = int(args.target_size)
    target_point = CouplingPoint(float(args.target_r1), float(args.target_r2))
    scheduled_sizes = [*execution_sizes, *([] if target_size in execution_sizes else [target_size])]

    fine_points = [CouplingPoint(r1, r2) for r1 in grid_values for r2 in grid_values]
    donors = [CouplingPoint(r1, r2) for r1 in donor_values for r2 in donor_values]
    payloads: dict[tuple[int, str], dict[str, Any]] = {}

    labels = [label for label, _ in _selected_specs_for_size(execution_sizes[0])]
    for size_index, size in enumerate(scheduled_sizes):
        points = fine_points if size in fit_sizes else [target_point]
        for point_index, point in enumerate(points):
            bundle_path = _run_point(
                exe,
                output_root,
                size=size,
                point=point,
                n_traj=int(args.n_traj),
                n_therm=int(args.n_therm),
                n_skip=int(args.n_skip),
                seed=int(args.seed_base) + size_index * 1000 + point_index,
            )
            payloads[(size, point.tag)] = _parse_selected_bundle(bundle_path, labels)

    target_rows = _panel_values(payloads[(target_size, target_point.tag)])
    target_x = 1.0 / float(target_size)

    direct_panel_cache: dict[tuple[int, str], dict[str, dict[str, float]]] = {}
    reweighted_panel_cache: dict[tuple[int, str, str], dict[str, dict[str, float]]] = {}
    direct_family_by_tag: dict[str, dict[int, dict[str, dict[str, float]]]] = {}
    reweighted_family_by_tag: dict[str, dict[int, dict[str, dict[str, float]]]] = {}
    donor_by_tag: dict[str, CouplingPoint] = {}
    rows: list[dict[str, Any]] = []
    for point in fine_points:
        donor = _nearest_donor(point, donors)

        direct_family: dict[int, dict[str, dict[str, float]]] = {}
        reweighted_family: dict[int, dict[str, dict[str, float]]] = {}
        for size in fit_sizes:
            direct_key = (size, point.tag)
            if direct_key not in direct_panel_cache:
                direct_panel_cache[direct_key] = _panel_values(payloads[direct_key])
            direct_family[size] = direct_panel_cache[direct_key]

            reweight_key = (size, donor.tag, point.tag)
            if reweight_key not in reweighted_panel_cache:
                reweighted_panel_cache[reweight_key] = _panel_values(payloads[(size, donor.tag)], target_point=point)
            reweighted_family[size] = reweighted_panel_cache[reweight_key]

        direct_family_by_tag[point.tag] = dict(direct_family)
        reweighted_family_by_tag[point.tag] = dict(reweighted_family)
        donor_by_tag[point.tag] = donor

        direct_score = _score_family(direct_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        reweighted_score = _score_family(reweighted_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        rows.append(
            {
                "point_tag": point.tag,
                "r1": float(point.r1),
                "r2": float(point.r2),
                "donor_r1": float(donor.r1),
                "donor_r2": float(donor.r2),
                "direct_z2_sum": float(direct_score["z2_sum"]),
                "reweighted_z2_sum": float(reweighted_score["z2_sum"]),
                "delta_z2_sum": float(reweighted_score["z2_sum"] - direct_score["z2_sum"]),
                "direct_panel_z": direct_score["panel_z"],
                "reweighted_panel_z": reweighted_score["panel_z"],
                "reweighted_panel_min_neff_fraction": reweighted_score["panel_min_neff_fraction"],
            }
        )

    best_direct = min(rows, key=lambda row: float(row["direct_z2_sum"]))
    best_reweighted = min(rows, key=lambda row: float(row["reweighted_z2_sum"]))
    delta_values = np.asarray([float(row["delta_z2_sum"]) for row in rows], dtype=float)
    summary = {
        "grid_values": grid_values,
        "donor_values": donor_values,
        "fit_sizes": fit_sizes,
        "target_size": target_size,
        "target_point": {"r1": float(target_point.r1), "r2": float(target_point.r2)},
        "n_traj": int(args.n_traj),
        "n_therm": int(args.n_therm),
        "n_skip": int(args.n_skip),
        "rows": rows,
        "summary": {
            "best_direct": best_direct,
            "best_reweighted": best_reweighted,
            "delta_rmse": float(np.sqrt(np.mean(np.square(delta_values)))),
            "delta_max_abs": float(np.max(np.abs(delta_values))),
        },
    }
    _write_summary(output_root, rows, summary)
    _save_plot(output_root, rows, donors, target_point)
    _save_fss_plots(
        output_root,
        rows=rows,
        direct_family_by_tag=direct_family_by_tag,
        reweighted_family_by_tag=reweighted_family_by_tag,
        donor_by_tag=donor_by_tag,
        fit_sizes=fit_sizes,
        target_size=target_size,
        target_x=target_x,
        target_rows=target_rows,
    )
    print(json.dumps(summary["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()