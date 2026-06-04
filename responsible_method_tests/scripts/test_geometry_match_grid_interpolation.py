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


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    ensure_dir(str(output_root))
    exe = ensure_simulator(DEFAULT_EXECUTION)

    grid_values = _unique_sorted(list(args.grid_values))
    donor_values = _unique_sorted(list(args.donor_values))
    fit_sizes = sorted({int(value) for value in args.fit_sizes})
    target_size = int(args.target_size)
    target_point = CouplingPoint(float(args.target_r1), float(args.target_r2))

    fine_points = [CouplingPoint(r1, r2) for r1 in grid_values for r2 in grid_values]
    donors = [CouplingPoint(r1, r2) for r1 in donor_values for r2 in donor_values]
    payloads: dict[tuple[int, str], dict[str, Any]] = {}

    labels = [label for label, _ in _selected_specs_for_size(fit_sizes[0])]
    for size_index, size in enumerate([*fit_sizes, target_size]):
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

        direct_score = _score_family(direct_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        reweighted_score = _score_family(reweighted_family, fit_sizes=fit_sizes, target_x=target_x, target_rows=target_rows)
        rows.append(
            {
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
    print(json.dumps(summary["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()