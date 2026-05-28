#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
from collections import defaultdict
from typing import Any

import matplotlib
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay, QhullError

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(HERE, ".."))
WORKSPACE_ROOT = os.path.normpath(os.path.join(PROJECT_ROOT, ".."))
STANDARD_ROOT = os.path.join(PROJECT_ROOT, "standard")
K_FROM_CONTINUUM_ROOT = os.path.join(WORKSPACE_ROOT, "K_from_continuum")
SIMULATOR_BIN = os.path.join(K_FROM_CONTINUUM_ROOT, "bin", "ising_tri_twisted_parallelogram")
JACKKNIFE_CACHE_ROOT = os.path.join(STANDARD_ROOT, "_jackknife_samples")

DEFAULT_UNTWISTED_DIR = os.path.join(
    STANDARD_ROOT,
    "data",
    "acute456",
    "untwisted",
    "r1_4p702783__r2_7p353910",
)
DEFAULT_TWISTED_DAT = os.path.join(
    STANDARD_ROOT,
    "data",
    "acute456",
    "twisted",
    "reference",
    "Lx144_Ly144_Tx72_Ty24",
    "two_point_all_to_all.dat",
)
DEFAULT_OUTPUT = os.path.join(STANDARD_ROOT, "acute456_center_pointwise_fss.png")
DEFAULT_OUTPUT_NORMFREE = os.path.join(STANDARD_ROOT, "acute456_center_pointwise_fss_normfree.png")
DEFAULT_OUTPUT_L8_RATIO = os.path.join(STANDARD_ROOT, "acute456_center_pointwise_fss_l8_ratio.png")

UNTWISTED_LATTICES = (
    (8, 8, 0, 0),
    (16, 16, 0, 0),
    (24, 24, 0, 0),
    (32, 32, 0, 0),
    (48, 48, 0, 0),
    (64, 64, 0, 0),
)
TWISTED_LATTICE = (144, 144, 72, 24)

INT_COLUMNS = {"d", "m", "n", "sample"}

RUN_ID_PATTERN = re.compile(
    r"^(?P<lx>\d+)x(?P<ly>\d+)_t(?P<tx>-?\d+)x(?P<ty>-?\d+)_k"
    r"(?P<k1>-?\d+)_(?P<k2>-?\d+)_(?P<k3>-?\d+)_(?P<kt>-?\d+)$"
)
SEED_FILE_PATTERN = re.compile(r"_(?P<seed>[0-9A-Fa-f]{8})\.dat$")


def _infer_family_label(untwisted_dir: str) -> str:
    family_name = os.path.basename(os.path.normpath(untwisted_dir))
    if family_name.startswith("r1_") and "__r2_" in family_name:
        r1_token, r2_token = family_name.split("__r2_", 1)
        r1_value = r1_token[len("r1_"):].replace("p", ".")
        r2_value = r2_token.replace("p", ".")
        return f"r=({r1_value}, {r2_value})"
    return family_name


def _load_json(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _meta_path_for_dat(dat_path: str) -> str:
    if not dat_path.endswith(".dat"):
        raise ValueError(f"expected .dat path, got {dat_path}")
    return dat_path[:-4] + ".meta.json"


def _parse_run_id(source_run_dir: str) -> dict[str, Any]:
    run_id = os.path.basename(os.path.normpath(source_run_dir))
    match = RUN_ID_PATTERN.match(run_id)
    if match is None:
        raise ValueError(f"could not parse run id from {source_run_dir}")
    return {
        "run_id": run_id,
        "kt": float(match.group("kt")) / 1000.0,
    }


def _find_seed_from_run_dir(source_run_dir: str) -> int:
    for name in sorted(os.listdir(source_run_dir)):
        if name.startswith("two_point_") or name.startswith("traces_"):
            continue
        match = SEED_FILE_PATTERN.search(name)
        if match is not None:
            return int(match.group("seed"), 16)
    raise FileNotFoundError(f"could not find seed-tagged run file in {source_run_dir}")


def _single_disp_sample_name(m_value: int, n_value: int) -> str:
    return f"single_disp_m{int(m_value)}_n{int(n_value)}_samples.dat"


def _load_single_disp_samples(sample_path: str) -> dict[str, np.ndarray]:
    rows = _load_dat_rows(sample_path)
    corr_samples = np.asarray([float(row["corr"]) for row in rows], dtype=float)
    mag_samples = np.asarray([float(row["mag"]) for row in rows], dtype=float)
    return {
        "corr": corr_samples,
        "mag": mag_samples,
    }


def _ensure_single_disp_samples(
    dat_path: str,
    *,
    m_value: int,
    n_value: int,
) -> str:
    meta = _load_json(_meta_path_for_dat(dat_path))
    source_all_to_all = str(meta["source_all_to_all_file"])
    source_run_dir = os.path.dirname(source_all_to_all)
    run_info = _parse_run_id(source_run_dir)
    sample_name = _single_disp_sample_name(m_value, n_value)
    cache_root = os.path.join(JACKKNIFE_CACHE_ROOT, str(meta["label"]))
    sample_path = os.path.join(cache_root, str(run_info["run_id"]), sample_name)
    if os.path.exists(sample_path):
        return sample_path

    if not os.path.exists(SIMULATOR_BIN):
        raise FileNotFoundError(f"simulator binary not found: {SIMULATOR_BIN}")

    os.makedirs(cache_root, exist_ok=True)
    seed = _find_seed_from_run_dir(source_run_dir)
    mc_payload = dict(meta.get("mc", {}))
    cmd = [
        SIMULATOR_BIN,
        "-X", str(int(meta["Lx"])),
        "-Y", str(int(meta["Ly"])),
        "-P", str(int(meta["Tx"])),
        "-Q", str(int(meta["Ty"])),
        "-a", str(float(meta["k1"])),
        "-b", str(float(meta["k2"])),
        "-c", str(float(meta["k3"])),
        "-g", str(float(run_info["kt"])),
        "-B", str(float(meta["production_beta"])),
        "-S", str(int(seed)),
        "-h", str(int(mc_payload["n_therm"])),
        "-t", str(int(mc_payload["n_traj"])),
        "-s", str(int(mc_payload["n_skip"])),
        "-d", cache_root,
        "-m", str(int(m_value)),
        "-n", str(int(n_value)),
        "-J", sample_name,
    ]
    completed = subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
        cwd=WORKSPACE_ROOT,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "single-displacement sample regeneration failed for "
            f"{dat_path} @ (m,n)=({int(m_value)},{int(n_value)}):\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    if not os.path.exists(sample_path):
        raise FileNotFoundError(f"expected sample file was not written: {sample_path}")
    return sample_path


def _connected_ratio_with_true_jackknife(
    dat_path: str,
    *,
    point_m: int,
    point_n: int,
    anchor_m: int,
    anchor_n: int,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]],
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]],
) -> tuple[float, float]:
    cache_key = (dat_path, int(point_m), int(point_n), int(anchor_m), int(anchor_n))
    cached = ratio_cache.get(cache_key)
    if cached is not None:
        return cached

    point_key = (dat_path, int(point_m), int(point_n))
    anchor_key = (dat_path, int(anchor_m), int(anchor_n))
    point_payload = sample_cache.get(point_key)
    if point_payload is None:
        point_payload = _load_single_disp_samples(
            _ensure_single_disp_samples(dat_path, m_value=int(point_m), n_value=int(point_n))
        )
        sample_cache[point_key] = point_payload
    anchor_payload = sample_cache.get(anchor_key)
    if anchor_payload is None:
        anchor_payload = _load_single_disp_samples(
            _ensure_single_disp_samples(dat_path, m_value=int(anchor_m), n_value=int(anchor_n))
        )
        sample_cache[anchor_key] = anchor_payload

    point_corr = np.asarray(point_payload["corr"], dtype=float)
    point_mag = np.asarray(point_payload["mag"], dtype=float)
    anchor_corr = np.asarray(anchor_payload["corr"], dtype=float)
    anchor_mag = np.asarray(anchor_payload["mag"], dtype=float)
    if point_corr.shape != anchor_corr.shape or point_mag.shape != anchor_mag.shape:
        raise ValueError(f"sample length mismatch for true ratio jackknife from {dat_path}")
    if not np.allclose(point_mag, anchor_mag, rtol=0.0, atol=1.0e-14):
        raise ValueError(f"magnetization samples do not align for true ratio jackknife from {dat_path}")

    n_samples = int(point_corr.size)
    if n_samples <= 1:
        raise ValueError(f"need at least two samples for true ratio jackknife: {dat_path}")

    mean_mag = float(np.mean(point_mag))
    point_conn = float(np.mean(point_corr) - mean_mag * mean_mag)
    anchor_conn = float(np.mean(anchor_corr) - mean_mag * mean_mag)
    if anchor_conn == 0.0:
        raise ZeroDivisionError(f"anchor connected correlator vanished in {dat_path}")
    ratio_value = point_conn / anchor_conn

    sum_point = float(np.sum(point_corr))
    sum_anchor = float(np.sum(anchor_corr))
    sum_mag = float(np.sum(point_mag))
    leave_one_out: list[float] = []
    for idx in range(n_samples):
        mean_mag_leave = (sum_mag - float(point_mag[idx])) / float(n_samples - 1)
        point_conn_leave = (sum_point - float(point_corr[idx])) / float(n_samples - 1) - mean_mag_leave * mean_mag_leave
        anchor_conn_leave = (sum_anchor - float(anchor_corr[idx])) / float(n_samples - 1) - mean_mag_leave * mean_mag_leave
        if anchor_conn_leave == 0.0:
            raise ZeroDivisionError(f"anchor connected correlator vanished in leave-one-out for {dat_path}")
        leave_one_out.append(point_conn_leave / anchor_conn_leave)

    jack_values = np.asarray(leave_one_out, dtype=float)
    jack_var = float((float(n_samples) - 1.0) / float(n_samples) * np.sum(np.square(jack_values - ratio_value)))
    ratio_sigma = math.sqrt(max(jack_var, 0.0))
    result = (float(ratio_value), float(ratio_sigma))
    ratio_cache[cache_key] = result
    return result


def _load_dat_rows(path: str) -> list[dict[str, Any]]:
    columns: list[str] | None = None
    rows: list[dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line[1:].strip()
                if header.startswith("columns:"):
                    columns = header.split(":", 1)[1].strip().split()
                else:
                    maybe_columns = header.split()
                    if maybe_columns:
                        columns = maybe_columns
                continue
            if columns is None:
                raise ValueError(f"missing column header in {path}")
            parts = line.split()
            if len(parts) != len(columns):
                raise ValueError(
                    f"row in {path} has {len(parts)} fields, expected {len(columns)}: {line}"
                )
            row: dict[str, Any] = {}
            for key, value in zip(columns, parts):
                if key in INT_COLUMNS:
                    row[key] = int(value)
                else:
                    row[key] = float(value)
            rows.append(row)
    return rows


def _boundary_paths(lx: int, ly: int, tx: int, ty: int) -> list[tuple[int, int]]:
    return [
        (lx, ty),
        (tx, -ly),
        (-lx - tx, ly - ty),
    ]


def _to_ab(
    m_value: int,
    n_value: int,
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
    a_value = (float(dn_b) * float(m_value) - float(dm_b) * float(n_value)) / det
    b_value = (-float(dn_a) * float(m_value) + float(dm_a) * float(n_value)) / det
    return a_value, b_value


def _wrap_unit(value: float) -> float:
    return value % 1.0


def _volume(lattice: tuple[int, int, int, int]) -> int:
    lx, ly, tx, ty = lattice
    return abs(int(lx) * int(ly) + int(tx) * int(ty))


def _sqrt_volume(lattice: tuple[int, int, int, int]) -> float:
    return float(math.sqrt(float(_volume(lattice))))


def _evaluate_power_model_on_x(x_values: np.ndarray, A_value: float, B_value: float, omega_value: float) -> np.ndarray:
    return A_value + B_value * np.power(x_values, omega_value)


def _fit_blind_power_model(
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
        sigma_omega = (
            float(np.sqrt(max(float(pcov[2, 2]), 0.0)))
            if np.isfinite(pcov[2, 2]) else float("nan")
        )
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
        "n_used": int(len(x_array)),
        "fit_mode": str(fit_mode),
        "pcov": np.asarray(pcov, dtype=float),
    }


def _evaluate_shared_power_model(
    xdata: tuple[np.ndarray, np.ndarray],
    omega_value: float,
    *series_params: float,
) -> np.ndarray:
    series_index, x_values = xdata
    index_array = np.asarray(np.rint(series_index), dtype=int)
    x_array = np.asarray(x_values, dtype=float)
    y_fit = np.empty_like(x_array, dtype=float)
    n_series = len(series_params) // 2
    for idx in range(n_series):
        A_value = float(series_params[2 * idx])
        B_value = float(series_params[2 * idx + 1])
        mask = index_array == idx
        y_fit[mask] = _evaluate_power_model_on_x(x_array[mask], A_value, B_value, omega_value)
    return y_fit


def _fit_shared_blind_power_model(series_payloads: list[dict[str, Any]]) -> dict[str, Any]:
    if len(series_payloads) == 0:
        raise ValueError("shared blind power fit requires at least one series")

    initial_series_fits = [
        _fit_blind_power_model(payload["x"], payload["y"], payload["sigma"])
        for payload in series_payloads
    ]
    initial_omega = float(np.median([fit["omega"] for fit in initial_series_fits]))

    x_stack: list[np.ndarray] = []
    y_stack: list[np.ndarray] = []
    sigma_stack: list[np.ndarray] = []
    index_stack: list[np.ndarray] = []
    p0: list[float] = [initial_omega]
    lower_bounds: list[float] = [0.05]
    upper_bounds: list[float] = [6.0]
    for idx, (payload, fit_payload) in enumerate(zip(series_payloads, initial_series_fits)):
        x_values = np.asarray(payload["x"], dtype=float)
        y_values = np.asarray(payload["y"], dtype=float)
        sigma_values = np.asarray(payload["sigma"], dtype=float)
        sigma_fit = np.where(np.isfinite(sigma_values) & (sigma_values > 0.0), sigma_values, np.nan)
        finite_sigma = sigma_fit[np.isfinite(sigma_fit)]
        sigma_floor = float(np.median(finite_sigma)) if finite_sigma.size else 1.0
        sigma_fit = np.where(np.isfinite(sigma_fit) & (sigma_fit > 0.0), sigma_fit, sigma_floor)

        x_stack.append(x_values)
        y_stack.append(y_values)
        sigma_stack.append(sigma_fit)
        index_stack.append(np.full_like(x_values, float(idx), dtype=float))

        p0.extend([float(fit_payload["A"]), float(fit_payload["B"])])
        lower_bounds.extend([-np.inf, -np.inf])
        upper_bounds.extend([np.inf, np.inf])

    x_all = np.concatenate(x_stack)
    y_all = np.concatenate(y_stack)
    sigma_all = np.concatenate(sigma_stack)
    index_all = np.concatenate(index_stack)

    try:
        popt, pcov = curve_fit(
            _evaluate_shared_power_model,
            (index_all, x_all),
            y_all,
            p0=p0,
            sigma=sigma_all,
            absolute_sigma=True,
            bounds=(lower_bounds, upper_bounds),
            maxfev=50000,
        )
        omega_value = float(popt[0])
        sigma_omega = (
            float(np.sqrt(max(float(pcov[0, 0]), 0.0)))
            if np.isfinite(pcov[0, 0]) else float("nan")
        )
        fit_mode = "shared_power_fit"
        series_results: list[dict[str, Any]] = []
        for idx in range(len(series_payloads)):
            A_index = 1 + 2 * idx
            B_index = 2 + 2 * idx
            sigma_A = (
                float(np.sqrt(max(float(pcov[A_index, A_index]), 0.0)))
                if np.isfinite(pcov[A_index, A_index]) else float("nan")
            )
            series_results.append(
                {
                    "A": float(popt[A_index]),
                    "sigma_A": sigma_A,
                    "B": float(popt[B_index]),
                    "param_indices": (0, A_index, B_index),
                }
            )
    except Exception:
        omega_value = initial_omega
        sigma_omega = float("nan")
        fit_mode = "shared_fit_failed"
        pcov = np.full((1 + 2 * len(series_payloads), 1 + 2 * len(series_payloads)), np.nan, dtype=float)
        series_results = [
            {
                "A": float(fit_payload["A"]),
                "sigma_A": float(fit_payload["sigma_A"]),
                "B": float(fit_payload["B"]),
                "param_indices": (0, 1 + 2 * idx, 2 + 2 * idx),
            }
            for idx, fit_payload in enumerate(initial_series_fits)
        ]

    return {
        "omega": omega_value,
        "sigma_omega": sigma_omega,
        "fit_mode": fit_mode,
        "series": series_results,
        "pcov": np.asarray(pcov, dtype=float),
    }


def _shared_fit_target_prediction(
    shared_fit: dict[str, Any],
    *,
    series_index: int,
    target_x: float,
    target_value: float,
    target_sigma: float,
) -> dict[str, float]:
    fit_payload = dict(shared_fit["series"][series_index])
    omega_value = float(shared_fit["omega"])
    A_value = float(fit_payload["A"])
    B_value = float(fit_payload["B"])
    pred_value = float(_evaluate_power_model_on_x(np.asarray([target_x], dtype=float), A_value, B_value, omega_value)[0])

    pred_sigma = float("nan")
    pcov = np.asarray(shared_fit.get("pcov"), dtype=float)
    if pcov.ndim == 2 and pcov.size > 0:
        omega_index, A_index, B_index = [int(value) for value in fit_payload["param_indices"]]
        sub_cov = pcov[np.ix_([omega_index, A_index, B_index], [omega_index, A_index, B_index])]
        if np.all(np.isfinite(sub_cov)):
            x_pow = float(target_x ** omega_value)
            grad = np.asarray([
                B_value * x_pow * math.log(target_x),
                1.0,
                x_pow,
            ], dtype=float)
            pred_var = float(grad @ sub_cov @ grad)
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
        "pred_value": pred_value,
        "pred_sigma": pred_sigma,
        "delta": delta,
        "abs_delta": abs_delta,
        "z_value": z_value,
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


def _find_point_by_mn(
    rows: list[dict[str, Any]],
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    target_m: int,
    target_n: int,
) -> dict[str, Any]:
    for row in rows:
        if int(row["m"]) != int(target_m) or int(row["n"]) != int(target_n):
            continue
        a_raw, b_raw = _to_ab(
            int(row["m"]),
            int(row["n"]),
            lattice[0],
            lattice[1],
            lattice[2],
            lattice[3],
            embedding_cycles=embedding_cycles,
        )
        return {
            "m": int(row["m"]),
            "n": int(row["n"]),
            "a_wrap": _wrap_unit(a_raw),
            "b_wrap": _wrap_unit(b_raw),
        }
    raise KeyError(f"anchor point (m,n)=({target_m},{target_n}) not found on lattice {lattice}")


def _ratio_with_uncertainty(
    numerator: float,
    numerator_sigma: float,
    denominator: float,
    denominator_sigma: float,
) -> tuple[float, float]:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0.0:
        return float("nan"), float("nan")
    ratio = float(numerator) / float(denominator)
    if numerator == 0.0:
        rel_num = 0.0
    else:
        rel_num = float(numerator_sigma) / abs(float(numerator)) if np.isfinite(float(numerator_sigma)) else 0.0
    rel_den = float(denominator_sigma) / abs(float(denominator)) if np.isfinite(float(denominator_sigma)) else 0.0
    sigma = abs(ratio) * math.sqrt(rel_num ** 2 + rel_den ** 2)
    return float(ratio), float(sigma)


def _compute_aggregate_target_score(summary_rows: list[dict[str, Any]]) -> dict[str, float]:
    z_values = np.asarray(
        [float(row["target_z"]) for row in summary_rows if np.isfinite(float(row["target_z"]))],
        dtype=float,
    )
    abs_delta_values = np.asarray(
        [
            float(row["target_abs_delta"])
            for row in summary_rows
            if np.isfinite(float(row["target_abs_delta"]))
        ],
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
        "mean_abs_delta": (
            float(np.mean(abs_delta_values)) if abs_delta_values.size > 0 else float("nan")
        ),
    }


def _untwisted_dat_path(root: str, lattice: tuple[int, int, int, int]) -> str:
    lx, ly, tx, ty = lattice
    return os.path.join(root, f"Lx{lx}_Ly{ly}_Tx{tx}_Ty{ty}", "two_point_all_to_all.dat")


def _aggregate_by_fraction(
    rows: list[dict[str, Any]],
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
) -> dict[tuple[float, float], dict[str, Any]]:
    lx, ly, tx, ty = lattice
    grouped: dict[tuple[float, float], dict[str, Any]] = {}
    collectors: dict[tuple[float, float], dict[str, list[float]]] = defaultdict(
        lambda: {"value": [], "sigma": []}
    )
    for row in rows:
        conn_value = float(row["corr_conn"])
        sigma_value = float(row["err_conn"])
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
        grouped.setdefault(
            key,
            {
                "m": int(row["m"]),
                "n": int(row["n"]),
                "a_wrap": a_wrap,
                "b_wrap": b_wrap,
            },
        )
        collectors[key]["value"].append(conn_value)
        collectors[key]["sigma"].append(abs(sigma_value))

    for key, payload in grouped.items():
        payload["value"] = float(np.mean(collectors[key]["value"]))
        payload["sigma"] = float(np.mean(collectors[key]["sigma"]))
    return grouped


def _select_base_points(
    rows: list[dict[str, Any]],
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
    n_panels: int,
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
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
        if int(row["m"]) == 0 and 0.0 < b_wrap <= 0.5 + 1.0e-12:
            candidates.append(
                {
                    "m": int(row["m"]),
                    "n": int(row["n"]),
                    "a_wrap": a_wrap,
                    "b_wrap": b_wrap,
                }
            )
    if not candidates:
        raise ValueError("no representative non-origin base-lattice points found")
    candidates.sort(key=lambda item: (float(item["b_wrap"]), float(item["a_wrap"]), int(item["n"])))
    return candidates[: max(1, min(int(n_panels), len(candidates)))]


def _build_periodic_interpolator(
    rows: list[dict[str, Any]],
    lattice: tuple[int, int, int, int],
    *,
    embedding_cycles: tuple[int, int],
) -> dict[str, Any]:
    aggregated = _aggregate_by_fraction(rows, lattice, embedding_cycles=embedding_cycles)
    a_core: list[float] = []
    b_core: list[float] = []
    value_core: list[float] = []
    sigma_core: list[float] = []
    a_tiled: list[float] = []
    b_tiled: list[float] = []
    value_tiled: list[float] = []
    sigma_tiled: list[float] = []

    for key in sorted(aggregated):
        payload = aggregated[key]
        a_value = float(payload["a_wrap"])
        b_value = float(payload["b_wrap"])
        value = float(payload["value"])
        sigma = float(payload["sigma"])
        a_core.append(a_value)
        b_core.append(b_value)
        value_core.append(value)
        sigma_core.append(sigma)
        for delta_a in (-1.0, 0.0, 1.0):
            for delta_b in (-1.0, 0.0, 1.0):
                a_tiled.append(a_value + delta_a)
                b_tiled.append(b_value + delta_b)
                value_tiled.append(value)
                sigma_tiled.append(sigma)

    tiled_points = np.column_stack([np.asarray(a_tiled, dtype=float), np.asarray(b_tiled, dtype=float)])
    value_array = np.asarray(value_tiled, dtype=float)
    sigma_array = np.asarray(sigma_tiled, dtype=float)
    try:
        triangulation = Delaunay(tiled_points, qhull_options="QJ Qc Qbb Q12")
        value_interpolator = LinearNDInterpolator(triangulation, value_array, fill_value=np.nan)
        sigma_interpolator = LinearNDInterpolator(triangulation, sigma_array, fill_value=np.nan)
    except QhullError:
        value_interpolator = LinearNDInterpolator(tiled_points, value_array, fill_value=np.nan)
        sigma_interpolator = LinearNDInterpolator(tiled_points, sigma_array, fill_value=np.nan)

    return {
        "value_interpolator": value_interpolator,
        "sigma_interpolator": sigma_interpolator,
        "nearest_value": NearestNDInterpolator(tiled_points, value_array),
        "nearest_sigma": NearestNDInterpolator(tiled_points, sigma_array),
        "a_core": np.asarray(a_core, dtype=float),
        "b_core": np.asarray(b_core, dtype=float),
        "value_core": np.asarray(value_core, dtype=float),
        "sigma_core": np.asarray(sigma_core, dtype=float),
    }


def _evaluate_periodic(
    interpolator: Any,
    nearest: Any,
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
    nearest_value = np.asarray(nearest(np.asarray([[a_value, b_value]])), dtype=float).reshape(-1)
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


def _normalization_target_label(base_label: str, normalization_mode: str) -> str:
    if normalization_mode == "anchor_ratio":
        return f"{base_label} ratio"
    if normalization_mode == "l8_ratio":
        return f"{base_label} / L=8"
    return str(base_label)


def _build_twisted_target_payload(
    *,
    dat_path: str,
    lattice: tuple[int, int, int, int],
    embedding_cycles: tuple[int, int],
    anchor_point: dict[str, Any],
    label: str,
) -> dict[str, Any]:
    target_rows = _load_dat_rows(dat_path)
    exact_map = _aggregate_by_fraction(target_rows, lattice, embedding_cycles=embedding_cycles)
    interpolator = _build_periodic_interpolator(
        target_rows,
        lattice,
        embedding_cycles=embedding_cycles,
    )
    anchor_key = (round(float(anchor_point["a_wrap"]), 12), round(float(anchor_point["b_wrap"]), 12))
    anchor_value = _evaluate_periodic(
        interpolator["value_interpolator"],
        interpolator["nearest_value"],
        float(anchor_point["a_wrap"]),
        float(anchor_point["b_wrap"]),
        a_core=interpolator["a_core"],
        b_core=interpolator["b_core"],
        z_core=interpolator["value_core"],
    )
    anchor_sigma = _evaluate_periodic(
        interpolator["sigma_interpolator"],
        interpolator["nearest_sigma"],
        float(anchor_point["a_wrap"]),
        float(anchor_point["b_wrap"]),
        a_core=interpolator["a_core"],
        b_core=interpolator["b_core"],
        z_core=interpolator["sigma_core"],
    )
    return {
        "dat_path": str(dat_path),
        "lattice": tuple(int(v) for v in lattice),
        "label": str(label),
        "exact_map": exact_map,
        "interpolator": interpolator,
        "anchor_exact": exact_map.get(anchor_key),
        "anchor_value": float(anchor_value),
        "anchor_sigma": float(anchor_sigma),
        "target_x": 1.0 / _sqrt_volume(lattice),
    }


def _target_value_for_point(
    *,
    target_payload: dict[str, Any],
    point: dict[str, Any],
    normalization_mode: str,
    baseline_value: float,
    baseline_sigma: float,
    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]],
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]],
) -> tuple[float, float]:
    key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
    exact_payload = target_payload["exact_map"].get(key)
    interpolator = target_payload["interpolator"]
    value = _evaluate_periodic(
        interpolator["value_interpolator"],
        interpolator["nearest_value"],
        float(point["a_wrap"]),
        float(point["b_wrap"]),
        a_core=interpolator["a_core"],
        b_core=interpolator["b_core"],
        z_core=interpolator["value_core"],
    )
    sigma = _evaluate_periodic(
        interpolator["sigma_interpolator"],
        interpolator["nearest_sigma"],
        float(point["a_wrap"]),
        float(point["b_wrap"]),
        a_core=interpolator["a_core"],
        b_core=interpolator["b_core"],
        z_core=interpolator["sigma_core"],
    )
    if normalization_mode == "anchor_ratio":
        if exact_payload is not None and target_payload["anchor_exact"] is not None:
            value, sigma = _connected_ratio_with_true_jackknife(
                target_payload["dat_path"],
                point_m=int(exact_payload["m"]),
                point_n=int(exact_payload["n"]),
                anchor_m=int(target_payload["anchor_exact"]["m"]),
                anchor_n=int(target_payload["anchor_exact"]["n"]),
                sample_cache=sample_cache,
                ratio_cache=ratio_cache,
            )
        else:
            value, sigma = _ratio_with_uncertainty(
                float(value),
                float(sigma),
                float(target_payload["anchor_value"]),
                float(target_payload["anchor_sigma"]),
            )
    elif normalization_mode == "l8_ratio":
        value, sigma = _ratio_with_uncertainty(
            float(value),
            float(sigma),
            float(baseline_value),
            float(baseline_sigma),
        )
    return float(value), float(sigma)


def _plot_series(
    *,
    untwisted_dir: str,
    twisted_dat: str,
    twisted_lattice: tuple[int, int, int, int],
    secondary_twisted_dat: str | None,
    secondary_twisted_lattice: tuple[int, int, int, int] | None,
    secondary_twisted_label: str | None,
    dataset_label: str,
    fit_mode: str,
    output_path: str,
    untwisted_embedding_cycles: tuple[int, int],
    twisted_embedding_cycles: tuple[int, int],
    n_panels: int,
    normalization_mode: str,
    anchor_m: int,
    anchor_n: int,
) -> dict[str, Any]:
    normalization_mode = str(normalization_mode).strip().lower()
    if normalization_mode not in {"raw", "anchor_ratio", "l8_ratio"}:
        raise ValueError(f"unknown normalization mode: {normalization_mode}")
    fit_mode = str(fit_mode).strip().lower()
    if fit_mode not in {"individual", "shared"}:
        raise ValueError(f"unknown fit mode: {fit_mode}")
    family_label = _infer_family_label(untwisted_dir)

    smallest_lattice = UNTWISTED_LATTICES[0]
    smallest_rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, smallest_lattice))
    anchor_point = _find_point_by_mn(
        smallest_rows,
        smallest_lattice,
        embedding_cycles=untwisted_embedding_cycles,
        target_m=int(anchor_m),
        target_n=int(anchor_n),
    )
    selected_points = _select_base_points(
        smallest_rows,
        smallest_lattice,
        embedding_cycles=untwisted_embedding_cycles,
        n_panels=n_panels,
    )
    if normalization_mode == "anchor_ratio":
        selected_points = [
            point
            for point in selected_points
            if not (int(point["m"]) == int(anchor_point["m"]) and int(point["n"]) == int(anchor_point["n"]))
        ]
        if len(selected_points) == 0:
            raise ValueError("no non-anchor points remain for normalization-free plotting")

    untwisted_maps: dict[tuple[int, int, int, int], dict[tuple[float, float], dict[str, Any]]] = {}
    for lattice in UNTWISTED_LATTICES:
        rows = _load_dat_rows(_untwisted_dat_path(untwisted_dir, lattice))
        untwisted_maps[lattice] = _aggregate_by_fraction(rows, lattice, embedding_cycles=untwisted_embedding_cycles)

    sample_cache: dict[tuple[str, int, int], dict[str, np.ndarray]] = {}
    ratio_cache: dict[tuple[str, int, int, int, int], tuple[float, float]] = {}
    primary_target = _build_twisted_target_payload(
        dat_path=twisted_dat,
        lattice=twisted_lattice,
        embedding_cycles=twisted_embedding_cycles,
        anchor_point=anchor_point,
        label="large twisted target",
    )
    secondary_target: dict[str, Any] | None = None
    if secondary_twisted_dat is not None and secondary_twisted_lattice is not None:
        secondary_target = _build_twisted_target_payload(
            dat_path=secondary_twisted_dat,
            lattice=secondary_twisted_lattice,
            embedding_cycles=twisted_embedding_cycles,
            anchor_point=anchor_point,
            label=str(secondary_twisted_label or "smaller twisted target"),
        )

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    n_points = len(selected_points)
    n_cols = 2 if n_points > 1 else 1
    n_rows = int(math.ceil(float(n_points) / float(n_cols)))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.5 * n_cols, 4.8 * n_rows), squeeze=False)
    axes_flat = list(axes.ravel())

    if normalization_mode == "anchor_ratio":
        subtitle = (
            f"{dataset_label} standard {family_label}: normalization-free ratios G(p)/G(anchor), "
            f"anchor=(m,n)=({int(anchor_point['m'])},{int(anchor_point['n'])}), "
            f"sizes 8,16,24,32,48,64; twisted reference {twisted_lattice}"
        )
        fig.suptitle("Pointwise FSS of Normalization-Free Correlator Ratios", fontsize=15, y=0.98)
    elif normalization_mode == "l8_ratio":
        subtitle = (
            f"{dataset_label} standard {family_label}: normalization-free ratios G_L(p)/G_8(p) at matched points, "
            "sizes 8,16,24,32,48,64; twisted targets also normalized by untwisted L=8"
        )
        fig.suptitle("Pointwise FSS of Correlator Ratios to L=8", fontsize=15, y=0.98)
    else:
        subtitle = (
            f"{dataset_label} standard {family_label}: untwisted, "
            f"sizes 8,16,24,32,48,64; primary twisted reference {twisted_lattice}"
        )
        fig.suptitle("Pointwise FSS of Connected Correlators", fontsize=15, y=0.98)
    if secondary_target is not None:
        subtitle += f"; secondary target {secondary_target['lattice']}"
    fig.text(0.5, 0.952, subtitle, ha="center", va="top", fontsize=10, color="#444444")

    summary_rows: list[dict[str, Any]] = []
    holdout_z_values: list[float] = []
    target_x = float(primary_target["target_x"])
    anchor_key = (round(float(anchor_point["a_wrap"]), 12), round(float(anchor_point["b_wrap"]), 12))

    plot_payloads: list[dict[str, Any]] = []
    for point in selected_points:
        key = (round(float(point["a_wrap"]), 12), round(float(point["b_wrap"]), 12))
        baseline_payload = untwisted_maps[smallest_lattice].get(key)
        if baseline_payload is None:
            raise KeyError(
                f"missing L=8 baseline coordinate {key} for untwisted lattice {smallest_lattice}"
            )
        x_values: list[float] = []
        y_values: list[float] = []
        y_errors: list[float] = []
        size_labels: list[int] = []
        for lattice in UNTWISTED_LATTICES:
            payload = untwisted_maps[lattice].get(key)
            anchor_payload = untwisted_maps[lattice].get(anchor_key)
            if payload is None:
                raise KeyError(
                    f"missing fractional coordinate {key} for untwisted lattice {lattice}"
                )
            if anchor_payload is None:
                raise KeyError(
                    f"missing anchor coordinate {anchor_key} for untwisted lattice {lattice}"
                )
            x_values.append(1.0 / _sqrt_volume(lattice))
            if normalization_mode == "anchor_ratio":
                ratio, ratio_sigma = _connected_ratio_with_true_jackknife(
                    _untwisted_dat_path(untwisted_dir, lattice),
                    point_m=int(payload["m"]),
                    point_n=int(payload["n"]),
                    anchor_m=int(anchor_payload["m"]),
                    anchor_n=int(anchor_payload["n"]),
                    sample_cache=sample_cache,
                    ratio_cache=ratio_cache,
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
            else:
                y_values.append(float(payload["value"]))
                y_errors.append(float(payload["sigma"]))
            size_labels.append(int(lattice[0]))

        x_array = np.asarray(x_values, dtype=float)
        y_array = np.asarray(y_values, dtype=float)
        yerr_array = np.asarray(y_errors, dtype=float)
        order = np.argsort(x_array)
        plot_payloads.append(
            {
                "point": point,
                "x": x_array[order],
                "y": y_array[order],
                "sigma": yerr_array[order],
                "size_labels": [size_labels[idx] for idx in order],
                "baseline_value": float(baseline_payload["value"]),
                "baseline_sigma": float(baseline_payload["sigma"]),
            }
        )

    shared_fit: dict[str, Any] | None = None
    if fit_mode == "shared":
        shared_fit = _fit_shared_blind_power_model(plot_payloads)
        fit_summary_text = (
            f"blind shared-omega fit on untwisted sizes only: omega={float(shared_fit['omega']):.4f}"
            + (
                f" +/- {float(shared_fit['sigma_omega']):.4f}"
                if np.isfinite(float(shared_fit['sigma_omega'])) else ""
            )
            + f" ({str(shared_fit['fit_mode'])}); target shown only as holdout"
        )
    else:
        fit_summary_text = "blind individual fits on untwisted sizes only; target shown only as holdout"
    if normalization_mode == "anchor_ratio":
        fit_summary_text += "; ratio errors from true leave-one-out jackknife"
    if secondary_target is not None:
        fit_summary_text += "; secondary twisted target shown as comparison only"
    fig.text(
        0.5,
        0.928,
        fit_summary_text,
        ha="center",
        va="top",
        fontsize=9.5,
        color="#5a2a83",
    )

    legend_handles = None
    legend_labels = None

    for series_index, (axis, plot_payload) in enumerate(zip(axes_flat, plot_payloads)):
        point = plot_payload["point"]
        twisted_value, twisted_sigma = _target_value_for_point(
            target_payload=primary_target,
            point=point,
            normalization_mode=normalization_mode,
            baseline_value=float(plot_payload["baseline_value"]),
            baseline_sigma=float(plot_payload["baseline_sigma"]),
            sample_cache=sample_cache,
            ratio_cache=ratio_cache,
        )
        secondary_value_sigma: tuple[float, float] | None = None
        if secondary_target is not None:
            secondary_value_sigma = _target_value_for_point(
                target_payload=secondary_target,
                point=point,
                normalization_mode=normalization_mode,
                baseline_value=float(plot_payload["baseline_value"]),
                baseline_sigma=float(plot_payload["baseline_sigma"]),
                sample_cache=sample_cache,
                ratio_cache=ratio_cache,
            )

        x_array = np.asarray(plot_payload["x"], dtype=float)
        y_array = np.asarray(plot_payload["y"], dtype=float)
        yerr_array = np.asarray(plot_payload["sigma"], dtype=float)
        size_labels = list(plot_payload["size_labels"])
        holdout_prediction: dict[str, float] | None = None
        holdout_size: int | None = None
        if fit_mode == "shared":
            if shared_fit is None:
                raise RuntimeError("shared fit state was not initialized")
            fit_payload = dict(shared_fit["series"][series_index])
            fit_payload["omega"] = float(shared_fit["omega"])
            fit_payload["sigma_omega"] = float(shared_fit["sigma_omega"])
            fit_payload["fit_mode"] = str(shared_fit["fit_mode"])
            target_prediction = _shared_fit_target_prediction(
                shared_fit,
                series_index=series_index,
                target_x=target_x,
                target_value=float(twisted_value),
                target_sigma=float(twisted_sigma),
            )
        else:
            fit_payload = _fit_blind_power_model(x_array, y_array, yerr_array)
            target_prediction = _individual_fit_target_prediction(
                fit_payload,
                target_x=target_x,
                target_value=float(twisted_value),
                target_sigma=float(twisted_sigma),
            )
            if x_array.size >= 4:
                holdout_fit = _fit_blind_power_model(x_array[1:], y_array[1:], yerr_array[1:])
                holdout_prediction = _individual_fit_target_prediction(
                    holdout_fit,
                    target_x=float(x_array[0]),
                    target_value=float(y_array[0]),
                    target_sigma=float(yerr_array[0]),
                )
                holdout_size = int(size_labels[0])
                if np.isfinite(float(holdout_prediction["z_value"])):
                    holdout_z_values.append(float(holdout_prediction["z_value"]))
        x_fit = np.logspace(np.log10(target_x), np.log10(float(np.max(x_array)) * 1.05), 300)
        y_fit = _evaluate_power_model_on_x(
            x_fit,
            float(fit_payload["A"]),
            float(fit_payload["B"]),
            float(fit_payload["omega"]),
        )

        axis.errorbar(
            x_array,
            y_array,
            yerr=yerr_array,
            fmt="o-",
            color="#1f77b4",
            capsize=3,
            linewidth=1.6,
            markersize=5,
            label=(
                "untwisted"
                if normalization_mode == "raw"
                else "untwisted ratio"
                if normalization_mode == "anchor_ratio"
                else "untwisted / L=8"
            ),
        )
        axis.plot(
            x_fit,
            y_fit,
            color="#ff7f0e",
            linewidth=1.4,
            label=(
                "blind power fit"
                if normalization_mode == "raw"
                else "blind power fit (ratio)"
                if normalization_mode == "anchor_ratio"
                else "blind power fit (L=8 ratio)"
            ),
        )
        axis.axhline(twisted_value, color="#d62728", linestyle="--", linewidth=1.2, alpha=0.8)
        axis.errorbar(
            [target_x],
            [twisted_value],
            yerr=[twisted_sigma],
            fmt="*",
            color="#d62728",
            capsize=4,
            markersize=11,
            label=_normalization_target_label(str(primary_target["label"]), normalization_mode),
        )
        if secondary_target is not None and secondary_value_sigma is not None:
            secondary_value, secondary_sigma = secondary_value_sigma
            axis.axhline(secondary_value, color="#9467bd", linestyle=":", linewidth=1.2, alpha=0.9)
            axis.errorbar(
                [float(secondary_target["target_x"])],
                [secondary_value],
                yerr=[secondary_sigma],
                fmt="D",
                color="#9467bd",
                capsize=4,
                markersize=7,
                label=_normalization_target_label(str(secondary_target["label"]), normalization_mode),
            )

        for x_value, y_value, size_value in zip(x_array, y_array, size_labels):
            axis.annotate(
                str(size_value),
                (x_value, y_value),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
                color="#1f77b4",
            )

        axis.set_title(
            (
                "ratio to anchor | "
                if normalization_mode == "anchor_ratio"
                else "ratio to L=8 | "
                if normalization_mode == "l8_ratio"
                else "base 8x8 point "
            )
            +
            "" 
            f"(m,n)=({int(point['m'])},{int(point['n'])})\n"
            f"(a,b)=({float(point['a_wrap']):.3f}, {float(point['b_wrap']):.3f})",
            pad=20,
        )
        axis.set_xscale("log")
        axis.set_xlabel("1 / sqrt(V)")
        axis.set_ylabel(
            "connected correlator"
            if normalization_mode == "raw"
            else "normalized correlator ratio"
            if normalization_mode == "anchor_ratio"
            else "correlator ratio to L=8"
        )
        axis.grid(True, alpha=0.25)
        axis.text(
            0.02,
            1.01,
            f"blind fit: A={float(fit_payload['A']):.4f}±{float(fit_payload['sigma_A']):.4f}, "
            f"B={float(fit_payload['B']):+.4f}, ω={float(fit_payload['omega']):.4f}"
            + (
                f"±{float(fit_payload['sigma_omega']):.4f}" if np.isfinite(float(fit_payload['sigma_omega'])) else ""
            )
            + "\n"
            f"pred@target={float(target_prediction['pred_value']):.4f}, "
            f"|Δ|={float(target_prediction['abs_delta']):.4f}, "
            f"z={float(target_prediction['z_value']):+.2f}"
            + (
                "\n"
                + f"holdout L={holdout_size}: |Δ|={float(holdout_prediction['abs_delta']):.4f}, "
                + f"z={float(holdout_prediction['z_value']):+.2f}"
                if holdout_prediction is not None and holdout_size is not None
                else ""
            ),
            transform=axis.transAxes,
            fontsize=7.5,
            va="bottom",
            ha="left",
            color="#555555",
            clip_on=False,
        )
        if legend_handles is None or legend_labels is None:
            legend_handles, legend_labels = axis.get_legend_handles_labels()

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
                "fit_sigma_omega": float(fit_payload["sigma_omega"]),
                "fit_mode": str(fit_payload["fit_mode"]),
                "pred_target": float(target_prediction["pred_value"]),
                "pred_target_sigma": float(target_prediction["pred_sigma"]),
                "target_abs_delta": float(target_prediction["abs_delta"]),
                "target_z": float(target_prediction["z_value"]),
                "holdout_size": int(holdout_size) if holdout_size is not None else -1,
                "holdout_abs_delta": (
                    float(holdout_prediction["abs_delta"]) if holdout_prediction is not None else float("nan")
                ),
                "holdout_z": (
                    float(holdout_prediction["z_value"]) if holdout_prediction is not None else float("nan")
                ),
            }
        )

    for axis in axes_flat[len(selected_points):]:
        axis.axis("off")

    aggregate_score = _compute_aggregate_target_score(summary_rows)
    holdout_rms_z = float("nan")
    holdout_mean_abs_z = float("nan")
    if holdout_z_values:
        holdout_array = np.asarray(holdout_z_values, dtype=float)
        holdout_rms_z = float(np.sqrt(np.mean(np.square(holdout_array))))
        holdout_mean_abs_z = float(np.mean(np.abs(holdout_array)))
    fig.text(
        0.5,
        0.905,
        f"4-point target score: RMS z={float(aggregate_score['rms_z']):.3f}, "
        f"chi2={float(aggregate_score['chi2']):.3f}, "
        f"mean|z|={float(aggregate_score['mean_abs_z']):.3f}, "
        f"mean|Δ|={float(aggregate_score['mean_abs_delta']):.4f}"
        + (
            f"; largest-size holdout RMS z={holdout_rms_z:.3f}, mean|z|={holdout_mean_abs_z:.3f}"
            if np.isfinite(holdout_rms_z)
            else ""
        ),
        ha="center",
        va="top",
        fontsize=8.8,
        color="#444444",
    )

    if legend_handles is not None and legend_labels is not None:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.878),
            ncol=3,
            frameon=False,
        )

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.84], h_pad=2.6)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {
        "summary_rows": summary_rows,
        "aggregate_score": aggregate_score,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot pointwise finite-size-scaling panels for the acute456 center sweep in "
            "responsible_method_tests/standard, overlaid with the corresponding value from the "
            "large twisted reference lattice."
        )
    )
    parser.add_argument("--untwisted-dir", default=DEFAULT_UNTWISTED_DIR)
    parser.add_argument("--twisted-dat", default=DEFAULT_TWISTED_DAT)
    parser.add_argument("--secondary-twisted-dat", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--dataset-label", default="acute456")
    parser.add_argument("--fit-mode", choices=["individual", "shared"], default="individual")
    parser.add_argument("--n-panels", type=int, default=4)
    parser.add_argument("--untwisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--twisted-embedding-cycles", nargs=2, type=int, default=[0, 1])
    parser.add_argument("--twisted-lattice", nargs=4, type=int, default=list(TWISTED_LATTICE))
    parser.add_argument("--secondary-twisted-lattice", nargs=4, type=int, default=None)
    parser.add_argument("--secondary-twisted-label", default=None)
    parser.add_argument("--normalization-mode", choices=["raw", "anchor_ratio", "l8_ratio"], default="raw")
    parser.add_argument("--anchor-m", type=int, default=0)
    parser.add_argument("--anchor-n", type=int, default=-1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if (args.secondary_twisted_dat is None) != (args.secondary_twisted_lattice is None):
        raise SystemExit("--secondary-twisted-dat and --secondary-twisted-lattice must be used together")
    if args.output:
        output_path = os.path.abspath(args.output)
    elif args.normalization_mode == "raw":
        output_path = os.path.abspath(DEFAULT_OUTPUT)
    elif args.normalization_mode == "anchor_ratio":
        output_path = os.path.abspath(DEFAULT_OUTPUT_NORMFREE)
    else:
        output_path = os.path.abspath(DEFAULT_OUTPUT_L8_RATIO)
    plot_result = _plot_series(
        untwisted_dir=os.path.abspath(args.untwisted_dir),
        twisted_dat=os.path.abspath(args.twisted_dat),
        twisted_lattice=(
            int(args.twisted_lattice[0]),
            int(args.twisted_lattice[1]),
            int(args.twisted_lattice[2]),
            int(args.twisted_lattice[3]),
        ),
        secondary_twisted_dat=(os.path.abspath(args.secondary_twisted_dat) if args.secondary_twisted_dat else None),
        secondary_twisted_lattice=(
            (
                int(args.secondary_twisted_lattice[0]),
                int(args.secondary_twisted_lattice[1]),
                int(args.secondary_twisted_lattice[2]),
                int(args.secondary_twisted_lattice[3]),
            )
            if args.secondary_twisted_lattice
            else None
        ),
        secondary_twisted_label=(str(args.secondary_twisted_label) if args.secondary_twisted_label else None),
        dataset_label=str(args.dataset_label),
        fit_mode=str(args.fit_mode),
        output_path=output_path,
        untwisted_embedding_cycles=(int(args.untwisted_embedding_cycles[0]), int(args.untwisted_embedding_cycles[1])),
        twisted_embedding_cycles=(int(args.twisted_embedding_cycles[0]), int(args.twisted_embedding_cycles[1])),
        n_panels=max(1, int(args.n_panels)),
        normalization_mode=str(args.normalization_mode),
        anchor_m=int(args.anchor_m),
        anchor_n=int(args.anchor_n),
    )
    summary_rows = list(plot_result["summary_rows"])
    aggregate_score = dict(plot_result["aggregate_score"])
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
    for row in summary_rows:
        print(
            "selected_point "
            f"m={row['m']} n={row['n']} a={row['a_wrap']:.6f} b={row['b_wrap']:.6f} "
            f"twisted={row['twisted_value']:.6f} sigma={row['twisted_sigma']:.6f} "
            f"fit_A={row['fit_A']:.6f} fit_sigma_A={row['fit_sigma_A']:.6f} "
            f"fit_B={row['fit_B']:.6f} fit_omega={row['fit_omega']:.6f} mode={row['fit_mode']} "
            f"pred_target={row['pred_target']:.6f} pred_sigma={row['pred_target_sigma']:.6f} "
            f"abs_delta={row['target_abs_delta']:.6f} z={row['target_z']:.6f} "
            f"holdout_L={int(row['holdout_size'])} holdout_abs_delta={row['holdout_abs_delta']:.6f} "
            f"holdout_z={row['holdout_z']:.6f}"
        )


if __name__ == "__main__":
    main()