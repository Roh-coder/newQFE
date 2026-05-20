#!/usr/bin/env python3
"""Estimate beta_c(infinity) from pseudo-critical beta(L) data.

Reads a diag_fss_autocorr ``summary.dat`` table, keeps the original Taylor fit
in ``1/L``, and also compares a few cheap, fairly model-blind extrapolators:

* Taylor polynomial in ``1/L`` (baseline)
* free power law ``beta_c + a L^{-p}``
* shifted power law ``beta_c + a (L + L0)^{-p}``
* Bulirsch-Stoer sequence extrapolation
"""
from __future__ import annotations

import argparse
import json
import math
import os
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SUMMARY = os.path.join(
    HERE, "results", "diag_fss_autocorr_findbeta", "raw", "summary.dat"
)
EXACT_TRIANGULAR_ISING_BETA_C = math.log(3.0) / 4.0
BST_OMEGA_MIN = 0.2
BST_OMEGA_MAX = 4.0
BST_OMEGA_COUNT = 191


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit beta_pc(L) with several cheap extrapolators and report the "
            "continuum intercept."
        )
    )
    parser.add_argument(
        "--summary",
        default=DEFAULT_SUMMARY,
        help="Path to a diag_fss_autocorr raw/summary.dat table.",
    )
    parser.add_argument(
        "--degree",
        type=int,
        default=2,
        help="Polynomial degree in x = 1/L for the Taylor baseline. Default: 2.",
    )
    parser.add_argument(
        "--lmin",
        type=float,
        default=None,
        help="Minimum L value to include in the fit. Default: use all rows.",
    )
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Use 1/sigma weights when beta_sigma > 0 is available.",
    )
    parser.add_argument(
        "--window-scan",
        action="store_true",
        help="Also report Taylor-fit results for every allowed minimum-L window.",
    )
    parser.add_argument(
        "--json-out",
        default=None,
        help="Optional Taylor-fit JSON output path. Default: sibling of summary.dat.",
    )
    parser.add_argument(
        "--fit-out",
        default=None,
        help="Optional Taylor-fit table output path. Default: sibling of summary.dat.",
    )
    parser.add_argument(
        "--compare-json-out",
        default=None,
        help="Optional comparison-summary JSON path. Default: sibling of summary.dat.",
    )
    parser.add_argument(
        "--compare-plot-out",
        default=None,
        help="Optional comparison-plot PNG path. Default: sibling of summary.dat.",
    )
    return parser.parse_args()


def load_summary(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with open(path) as f:
        for line in f:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split()
            if len(parts) < 9:
                raise ValueError(
                    f"expected at least 9 columns in {path}, got: {text}"
                )
            rows.append((float(parts[0]), float(parts[7]), float(parts[8])))
    if not rows:
        raise ValueError(f"no data rows found in {path}")
    data = np.asarray(rows, dtype=float)
    order = np.argsort(data[:, 0])
    data = data[order]
    return data[:, 0], data[:, 1], data[:, 2]


def select_rows(
    L: np.ndarray, beta: np.ndarray, sigma: np.ndarray, lmin: float | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if lmin is None:
        mask = np.ones_like(L, dtype=bool)
    else:
        mask = L >= float(lmin)
    if not np.any(mask):
        raise ValueError(f"no rows satisfy L >= {lmin}")
    return L[mask], beta[mask], sigma[mask]


def format_tag_value(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def maybe_float(value: object) -> float | None:
    return None if value is None else float(value)


def covariance_sigma(covariance: np.ndarray | None, index: int) -> float | None:
    if covariance is None:
        return None
    if not np.all(np.isfinite(covariance)):
        return None
    variance = float(covariance[index, index])
    if variance < 0.0 or not np.isfinite(variance):
        return None
    return float(math.sqrt(variance))


def to_jsonable(value: object) -> object:
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): to_jsonable(val) for key, val in value.items()}
    if isinstance(value, np.ndarray):
        return [to_jsonable(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    return value


def filter_curve_fit_inputs(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    min_points: int,
    weighted: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    L_fit = L.copy()
    beta_fit = beta.copy()
    sigma_fit = sigma.copy()
    sigma_arg = None
    if weighted:
        positive = sigma_fit > 0.0
        if int(positive.sum()) < min_points:
            raise ValueError(
                "weighted fit requested, but not enough rows have beta_sigma > 0"
            )
        L_fit = L_fit[positive]
        beta_fit = beta_fit[positive]
        sigma_fit = sigma_fit[positive]
        sigma_arg = sigma_fit
    if len(L_fit) < min_points:
        raise ValueError(f"need at least {min_points} points; got {len(L_fit)}")
    return L_fit, beta_fit, sigma_fit, sigma_arg


def initial_beta_c_guess(L: np.ndarray, beta: np.ndarray) -> float:
    inv_L = 1.0 / L
    coeffs = np.polyfit(inv_L, beta, 1)
    return float(max(beta.max(), coeffs[-1]))


def fit_taylor(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    degree: int,
    weighted: bool,
) -> dict[str, object]:
    if degree < 0:
        raise ValueError("degree must be non-negative")
    if len(L) < degree + 1:
        raise ValueError(
            f"need at least {degree + 1} points for degree-{degree} fit; "
            f"got {len(L)}"
        )

    x = 1.0 / L
    y = beta.copy()
    s = sigma.copy()
    weights = None
    if weighted:
        positive = s > 0.0
        if positive.sum() < degree + 1:
            raise ValueError(
                "weighted fit requested, but not enough rows have beta_sigma > 0"
            )
        x = x[positive]
        y = y[positive]
        s = s[positive]
        L = L[positive]
        weights = 1.0 / s

    if len(L) == degree + 1:
        if weights is None:
            coeffs = np.polyfit(x, y, degree)
        else:
            coeffs = np.polyfit(x, y, degree, w=weights)
        cov = None
        intercept_sigma = None
    else:
        if weights is None:
            coeffs, cov = np.polyfit(x, y, degree, cov=True)
        else:
            coeffs, cov = np.polyfit(x, y, degree, w=weights, cov=True)
        intercept_sigma = covariance_sigma(cov, -1)
    y_fit = np.polyval(coeffs, x)
    residuals = y - y_fit

    return {
        "method": f"taylor_deg{degree}",
        "label": f"Taylor deg {degree}",
        "ansatz": "beta_pc(L) = sum_n c_n * (1/L)^n with c_0 = beta_c_continuum",
        "L": L,
        "inv_L": x,
        "beta": y,
        "sigma": s,
        "beta_fit": y_fit,
        "residual": residuals,
        "degree": degree,
        "weighted": weighted,
        "coefficients_descending_powers": coeffs,
        "covariance": cov,
        "parameters": {
            "coefficients_descending_powers": [float(value) for value in coeffs],
        },
        "beta_c_continuum": float(coeffs[-1]),
        "beta_c_continuum_sigma": intercept_sigma,
        "rms_residual": float(np.sqrt(np.mean(residuals**2))),
        "n_points": int(len(L)),
        "lmin": float(np.min(L)),
        "lmax": float(np.max(L)),
    }


def model_power_law(
    L: np.ndarray, beta_c: float, amplitude: float, exponent: float
) -> np.ndarray:
    return beta_c + amplitude * np.power(L, -exponent)


def model_shifted_power_law(
    L: np.ndarray,
    beta_c: float,
    amplitude: float,
    exponent: float,
    shift: float,
) -> np.ndarray:
    return beta_c + amplitude * np.power(L + shift, -exponent)


def evaluate_shifted_power_law_on_inv_L(
    inv_L: np.ndarray,
    beta_c: float,
    amplitude: float,
    exponent: float,
    shift: float,
) -> np.ndarray:
    values = np.full_like(inv_L, beta_c, dtype=float)
    positive = inv_L > 0.0
    if np.any(positive):
        L = 1.0 / inv_L[positive]
        values[positive] = model_shifted_power_law(
            L, beta_c, amplitude, exponent, shift
        )
    return values


def fit_power_law(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    weighted: bool,
) -> dict[str, object]:
    L_fit, y, s, sigma_arg = filter_curve_fit_inputs(
        L, beta, sigma, min_points=3, weighted=weighted
    )
    beta_c_guess = initial_beta_c_guess(L_fit, y)
    amplitude_guess = float((y[0] - beta_c_guess) * L_fit[0])
    lower = [float(y.min() - 0.10), -10.0, 0.05]
    upper = [float(y.max() + 0.10), 10.0, 6.0]
    p0 = [
        min(max(beta_c_guess, lower[0] + 1e-6), upper[0] - 1e-6),
        min(max(amplitude_guess, lower[1] + 1e-6), upper[1] - 1e-6),
        1.0,
    ]
    popt, pcov = curve_fit(
        model_power_law,
        L_fit,
        y,
        p0=p0,
        bounds=(lower, upper),
        sigma=sigma_arg,
        absolute_sigma=bool(weighted and sigma_arg is not None),
        maxfev=100000,
    )
    y_fit = model_power_law(L_fit, *popt)
    residuals = y - y_fit
    return {
        "method": "power_law",
        "label": "Free power",
        "ansatz": "beta_pc(L) = beta_c_continuum + a * L^(-p)",
        "L": L_fit,
        "inv_L": 1.0 / L_fit,
        "beta": y,
        "sigma": s,
        "beta_fit": y_fit,
        "residual": residuals,
        "covariance": pcov,
        "parameters": {
            "beta_c": float(popt[0]),
            "amplitude": float(popt[1]),
            "exponent": float(popt[2]),
        },
        "beta_c_continuum": float(popt[0]),
        "beta_c_continuum_sigma": covariance_sigma(pcov, 0),
        "rms_residual": float(np.sqrt(np.mean(residuals**2))),
        "n_points": int(len(L_fit)),
        "lmin": float(np.min(L_fit)),
        "lmax": float(np.max(L_fit)),
    }


def fit_shifted_power_law(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    weighted: bool,
    seed: dict[str, object],
) -> dict[str, object]:
    L_fit, y, s, sigma_arg = filter_curve_fit_inputs(
        L, beta, sigma, min_points=4, weighted=weighted
    )
    beta_c_guess = float(seed["beta_c_continuum"])
    amplitude_guess = float(seed["parameters"]["amplitude"])
    exponent_guess = float(seed["parameters"]["exponent"])
    shift_bound_lo = -0.75 * float(np.min(L_fit)) + 1e-6
    shift_bound_hi = 2.0 * float(np.min(L_fit))
    lower = [float(y.min() - 0.10), -10.0, 0.05, shift_bound_lo]
    upper = [float(y.max() + 0.10), 10.0, 6.0, shift_bound_hi]
    p0 = [
        min(max(beta_c_guess, lower[0] + 1e-6), upper[0] - 1e-6),
        min(max(amplitude_guess, lower[1] + 1e-6), upper[1] - 1e-6),
        min(max(exponent_guess, lower[2] + 1e-6), upper[2] - 1e-6),
        0.0,
    ]
    popt, pcov = curve_fit(
        model_shifted_power_law,
        L_fit,
        y,
        p0=p0,
        bounds=(lower, upper),
        sigma=sigma_arg,
        absolute_sigma=bool(weighted and sigma_arg is not None),
        maxfev=100000,
    )
    y_fit = model_shifted_power_law(L_fit, *popt)
    residuals = y - y_fit
    return {
        "method": "shifted_power_law",
        "label": "Shifted power",
        "ansatz": "beta_pc(L) = beta_c_continuum + a * (L + L0)^(-p)",
        "L": L_fit,
        "inv_L": 1.0 / L_fit,
        "beta": y,
        "sigma": s,
        "beta_fit": y_fit,
        "residual": residuals,
        "covariance": pcov,
        "parameters": {
            "beta_c": float(popt[0]),
            "amplitude": float(popt[1]),
            "exponent": float(popt[2]),
            "shift": float(popt[3]),
        },
        "beta_c_continuum": float(popt[0]),
        "beta_c_continuum_sigma": covariance_sigma(pcov, 0),
        "rms_residual": float(np.sqrt(np.mean(residuals**2))),
        "n_points": int(len(L_fit)),
        "lmin": float(np.min(L_fit)),
        "lmax": float(np.max(L_fit)),
    }


def bst_table(inv_L: np.ndarray, beta: np.ndarray, omega: float) -> np.ndarray:
    n = len(beta)
    table = np.full((n, n), np.nan, dtype=float)
    table[:, 0] = beta
    tiny = 1e-14
    for col in range(1, n):
        for row in range(0, n - col):
            if col == 1:
                denom = (inv_L[row] / inv_L[row + 1]) ** omega - 1.0
            else:
                delta_here = table[row + 1, col - 1] - table[row, col - 1]
                delta_prev = table[row + 1, col - 1] - table[row + 1, col - 2]
                if (
                    not np.isfinite(delta_here)
                    or not np.isfinite(delta_prev)
                    or abs(delta_prev) < tiny
                ):
                    continue
                denom = (
                    (inv_L[row] / inv_L[row + col]) ** omega
                    * (1.0 - delta_here / delta_prev)
                    - 1.0
                )
            if not np.isfinite(denom) or abs(denom) < tiny:
                continue
            table[row, col] = (
                table[row + 1, col - 1]
                + (table[row + 1, col - 1] - table[row, col - 1]) / denom
            )
    return table


def bst_estimate(inv_L: np.ndarray, beta: np.ndarray, omega: float) -> float | None:
    table = bst_table(inv_L, beta, omega)
    estimate = float(table[0, len(beta) - 1])
    return estimate if np.isfinite(estimate) else None


def fit_bst(L: np.ndarray, beta: np.ndarray) -> dict[str, object]:
    inv_L = 1.0 / L
    omega_grid = np.linspace(BST_OMEGA_MIN, BST_OMEGA_MAX, BST_OMEGA_COUNT)
    records = []
    max_start = max(0, len(L) - 4)
    for omega in omega_grid:
        window_estimates = []
        ok = True
        for start in range(0, max_start + 1):
            estimate = bst_estimate(inv_L[start:], beta[start:], omega)
            if estimate is None:
                ok = False
                break
            window_estimates.append(float(estimate))
        if not ok or not window_estimates:
            continue
        full_estimate = window_estimates[0]
        window_std = float(np.std(window_estimates))
        window_spread = float(np.max(window_estimates) - np.min(window_estimates))
        records.append(
            {
                "omega": float(omega),
                "estimate": float(full_estimate),
                "window_std": window_std,
                "window_spread": window_spread,
                "window_estimates": window_estimates,
            }
        )
    if not records:
        raise RuntimeError("BST scan failed for every omega in the search grid")
    best = min(records, key=lambda row: (row["window_std"], row["window_spread"]))
    return {
        "method": "bst",
        "label": "BST",
        "ansatz": "Bulirsch-Stoer sequence extrapolation in x = 1/L",
        "L": L,
        "inv_L": inv_L,
        "beta": beta,
        "sigma": np.zeros_like(beta),
        "beta_fit": None,
        "residual": None,
        "covariance": None,
        "parameters": {
            "best_omega": float(best["omega"]),
            "omega_grid_min": BST_OMEGA_MIN,
            "omega_grid_max": BST_OMEGA_MAX,
            "omega_grid_count": BST_OMEGA_COUNT,
        },
        "beta_c_continuum": float(best["estimate"]),
        "beta_c_continuum_sigma": None,
        "rms_residual": None,
        "window_stability": float(best["window_std"]),
        "window_spread": float(best["window_spread"]),
        "window_estimates": [float(value) for value in best["window_estimates"]],
        "omega_scan_best": float(best["omega"]),
        "n_points": int(len(L)),
        "lmin": float(np.min(L)),
        "lmax": float(np.max(L)),
    }


def default_output_paths(summary_path: str, degree: int, lmin: float) -> tuple[str, str]:
    root = os.path.dirname(os.path.abspath(summary_path))
    tag = f"deg{degree}_lmin{format_tag_value(lmin)}"
    json_out = os.path.join(root, f"beta_c_continuum_taylor_{tag}.json")
    fit_out = os.path.join(root, f"beta_c_continuum_taylor_{tag}.dat")
    return json_out, fit_out


def default_compare_output_paths(
    summary_path: str, degree: int, lmin: float
) -> tuple[str, str]:
    root = os.path.dirname(os.path.abspath(summary_path))
    tag = f"deg{degree}_lmin{format_tag_value(lmin)}"
    json_out = os.path.join(root, f"beta_c_continuum_compare_{tag}.json")
    plot_out = os.path.join(root, f"beta_c_continuum_compare_{tag}.png")
    return json_out, plot_out


def write_fit_table(path: str, fit: dict[str, object]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("# L inv_L beta_used beta_sigma beta_fit residual\n")
        rows = zip(
            fit["L"],
            fit["inv_L"],
            fit["beta"],
            fit["sigma"],
            fit["beta_fit"],
            fit["residual"],
        )
        for row in rows:
            f.write(" ".join(f"{float(value):.12g}" for value in row) + "\n")


def serialize_fit(
    fit: dict[str, object], summary_path: str, window_scan: list[dict[str, object]]
) -> dict[str, object]:
    payload = {
        "summary_file": os.path.abspath(summary_path),
        "ansatz": fit["ansatz"],
        "degree": int(fit["degree"]),
        "weighted": bool(fit["weighted"]),
        "n_points": int(fit["n_points"]),
        "lmin": float(fit["lmin"]),
        "lmax": float(fit["lmax"]),
        "beta_c_continuum": float(fit["beta_c_continuum"]),
        "beta_c_continuum_sigma": maybe_float(fit["beta_c_continuum_sigma"]),
        "coefficients_descending_powers": to_jsonable(
            fit["parameters"]["coefficients_descending_powers"]
        ),
        "rms_residual": float(fit["rms_residual"]),
        "exact_triangular_ising_beta_c": EXACT_TRIANGULAR_ISING_BETA_C,
        "delta_vs_exact": float(
            fit["beta_c_continuum"] - EXACT_TRIANGULAR_ISING_BETA_C
        ),
        "window_scan": window_scan,
    }
    return payload


def serialize_method_result(result: dict[str, object]) -> dict[str, object]:
    payload = {
        "method": result["method"],
        "label": result["label"],
        "ansatz": result["ansatz"],
        "n_points": int(result["n_points"]),
        "lmin": float(result["lmin"]),
        "lmax": float(result["lmax"]),
        "beta_c_continuum": float(result["beta_c_continuum"]),
        "beta_c_continuum_sigma": maybe_float(result["beta_c_continuum_sigma"]),
        "rms_residual": maybe_float(result["rms_residual"]),
        "delta_vs_exact": float(
            result["beta_c_continuum"] - EXACT_TRIANGULAR_ISING_BETA_C
        ),
        "parameters": to_jsonable(result.get("parameters", {})),
    }
    if result.get("window_stability") is not None:
        payload["window_stability"] = float(result["window_stability"])
    if result.get("window_spread") is not None:
        payload["window_spread"] = float(result["window_spread"])
    if result.get("window_estimates") is not None:
        payload["window_estimates"] = to_jsonable(result["window_estimates"])
    return payload


def build_window_scan(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    degree: int,
    weighted: bool,
) -> list[dict[str, object]]:
    results = []
    for lmin in sorted(set(float(value) for value in L)):
        subL, subbeta, subsigma = select_rows(L, beta, sigma, lmin)
        if len(subL) < degree + 1:
            continue
        fit = fit_taylor(subL, subbeta, subsigma, degree, weighted)
        results.append(
            {
                "lmin": float(lmin),
                "n_points": int(fit["n_points"]),
                "beta_c_continuum": float(fit["beta_c_continuum"]),
                "beta_c_continuum_sigma": maybe_float(
                    fit["beta_c_continuum_sigma"]
                ),
                "rms_residual": float(fit["rms_residual"]),
            }
        )
    return results


def build_comparison_methods(
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    weighted: bool,
    taylor_fit: dict[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    methods = [taylor_fit]
    failures: list[dict[str, str]] = []

    try:
        power_fit = fit_power_law(L, beta, sigma, weighted)
        methods.append(power_fit)
    except Exception as exc:  # pragma: no cover - runtime path
        power_fit = None
        failures.append({"method": "power_law", "error": str(exc)})

    if power_fit is not None:
        try:
            shifted_fit = fit_shifted_power_law(L, beta, sigma, weighted, power_fit)
            methods.append(shifted_fit)
        except Exception as exc:  # pragma: no cover - runtime path
            failures.append({"method": "shifted_power_law", "error": str(exc)})
    else:
        failures.append(
            {
                "method": "shifted_power_law",
                "error": "skipped because free power-law fit failed",
            }
        )

    try:
        bst_fit = fit_bst(L, beta)
        methods.append(bst_fit)
    except Exception as exc:  # pragma: no cover - runtime path
        failures.append({"method": "bst", "error": str(exc)})

    return methods, failures


def serialize_comparison(
    summary_path: str,
    methods: list[dict[str, object]],
    failures: list[dict[str, str]],
) -> dict[str, object]:
    return {
        "summary_file": os.path.abspath(summary_path),
        "exact_triangular_ising_beta_c": EXACT_TRIANGULAR_ISING_BETA_C,
        "methods": [serialize_method_result(method) for method in methods],
        "failures": failures,
    }


def print_window_scan(rows: Iterable[dict[str, object]]) -> None:
    print("window scan:")
    for row in rows:
        sigma = row["beta_c_continuum_sigma"]
        sigma_text = f"{sigma:.12f}" if sigma is not None else "n/a"
        print(
            "  "
            f"Lmin={int(row['lmin']) if float(row['lmin']).is_integer() else row['lmin']} "
            f"n={row['n_points']} "
            f"beta_c={row['beta_c_continuum']:.12f} "
            f"sigma={sigma_text} "
            f"rms={row['rms_residual']:.12e}"
        )


def print_method_summary(method: dict[str, object]) -> None:
    sigma = method["beta_c_continuum_sigma"]
    sigma_text = f"{sigma:.12f}" if sigma is not None else "n/a"
    rms = method.get("rms_residual")
    rms_text = f"{rms:.12e}" if rms is not None else "n/a"
    extras = []
    params = method.get("parameters", {})
    if method["method"] == "power_law":
        extras.append(f"p={params['exponent']:.4f}")
    elif method["method"] == "shifted_power_law":
        extras.append(f"p={params['exponent']:.4f}")
        extras.append(f"L0={params['shift']:.4f}")
    elif method["method"] == "bst":
        extras.append(f"omega*={method['omega_scan_best']:.4f}")
        extras.append(f"window_std={method['window_stability']:.3e}")
    extra_text = ("  " + "  ".join(extras)) if extras else ""
    print(
        f"{method['label']}: beta_c={method['beta_c_continuum']:.12f} "
        f"sigma={sigma_text} rms={rms_text}{extra_text}"
    )


def evaluate_method_on_inv_L(
    method: dict[str, object], inv_L: np.ndarray
) -> np.ndarray | None:
    params = method.get("parameters", {})
    if method["method"].startswith("taylor_deg"):
        coeffs = np.asarray(params["coefficients_descending_powers"], dtype=float)
        return np.polyval(coeffs, inv_L)
    if method["method"] == "power_law":
        return method["beta_c_continuum"] + params["amplitude"] * np.power(
            inv_L, params["exponent"]
        )
    if method["method"] == "shifted_power_law":
        return evaluate_shifted_power_law_on_inv_L(
            inv_L,
            method["beta_c_continuum"],
            params["amplitude"],
            params["exponent"],
            params["shift"],
        )
    return None


def make_comparison_plot(
    path: str,
    L: np.ndarray,
    beta: np.ndarray,
    sigma: np.ndarray,
    methods: list[dict[str, object]],
) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    inv_L = 1.0 / L
    x_curve = np.linspace(0.0, float(inv_L.max()) * 1.03, 500)
    colors = {
        methods[0]["method"]: "#1f77b4",
        "power_law": "#ff7f0e",
        "shifted_power_law": "#2ca02c",
        "bst": "#9467bd",
    }

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(9.0, 8.0),
        gridspec_kw={"height_ratios": [2.3, 1.2]},
    )

    if np.any(sigma > 0.0):
        ax_top.errorbar(
            inv_L,
            beta,
            yerr=sigma,
            fmt="o",
            color="#111111",
            mfc="white",
            label="Pseudo-critical data",
            zorder=5,
        )
    else:
        ax_top.scatter(
            inv_L,
            beta,
            color="#111111",
            s=55,
            label="Pseudo-critical data",
            zorder=5,
        )
    ax_top.axhline(
        EXACT_TRIANGULAR_ISING_BETA_C,
        color="#444444",
        ls=":",
        lw=1.4,
        label="Exact ln(3)/4",
    )
    for method in methods:
        color = colors.get(method["method"], None)
        curve = evaluate_method_on_inv_L(method, x_curve)
        if curve is not None:
            ax_top.plot(x_curve, curve, lw=2.0, color=color, label=method["label"])
            ax_top.scatter(
                [0.0],
                [method["beta_c_continuum"]],
                marker="*",
                s=90,
                color=color,
                zorder=6,
            )
        else:
            ax_top.hlines(
                method["beta_c_continuum"],
                -0.0015,
                0.010,
                colors=color,
                linestyles="--",
                lw=1.3,
            )
            ax_top.scatter(
                [0.0],
                [method["beta_c_continuum"]],
                marker="D",
                s=65,
                color=color,
                label=f"{method['label']} (limit)",
                zorder=6,
            )
    ax_top.set_xlim(-0.004, float(inv_L.max()) * 1.05)
    ax_top.set_xlabel("1 / L")
    ax_top.set_ylabel("beta_pc(L)")
    ax_top.set_title("Cheap continuum extrapolators for beta_pc(L)")
    ax_top.grid(alpha=0.25)
    ax_top.legend(ncol=2, fontsize=9, frameon=True)

    method_names = [method["label"] for method in methods]
    method_values = [float(method["beta_c_continuum"]) for method in methods]
    method_sigmas = [maybe_float(method["beta_c_continuum_sigma"]) for method in methods]
    xpos = np.arange(len(methods), dtype=float)
    ax_bottom.axhline(
        EXACT_TRIANGULAR_ISING_BETA_C,
        color="#444444",
        ls=":",
        lw=1.4,
        label="Exact ln(3)/4",
    )
    for index, method in enumerate(methods):
        color = colors.get(method["method"], None)
        sigma_val = method_sigmas[index]
        if sigma_val is not None:
            ax_bottom.errorbar(
                xpos[index],
                method_values[index],
                yerr=sigma_val,
                fmt="o",
                color=color,
                capsize=4,
            )
        else:
            ax_bottom.scatter(xpos[index], method_values[index], color=color, s=55)
        annotation = [f"Delta={method_values[index] - EXACT_TRIANGULAR_ISING_BETA_C:+.2e}"]
        if method["rms_residual"] is not None:
            annotation.append(f"RMS={method['rms_residual']:.1e}")
        if method["method"] == "bst":
            annotation.append(f"omega*={method['omega_scan_best']:.2f}")
        ax_bottom.text(
            xpos[index],
            method_values[index] + 0.0004,
            "\n".join(annotation),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax_bottom.set_xticks(xpos, method_names)
    ax_bottom.set_ylabel("beta_c(infinity)")
    ax_bottom.set_title("Extrapolated continuum limits")
    ax_bottom.grid(alpha=0.25)
    y_all = method_values + [EXACT_TRIANGULAR_ISING_BETA_C]
    y_margin = max(0.002, 0.18 * (max(y_all) - min(y_all) + 1e-6))
    ax_bottom.set_ylim(min(y_all) - y_margin, max(y_all) + y_margin)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    L_all, beta_all, sigma_all = load_summary(args.summary)
    L, beta, sigma = select_rows(L_all, beta_all, sigma_all, args.lmin)

    fit = fit_taylor(L, beta, sigma, args.degree, args.weighted)
    window_scan = []
    if args.window_scan:
        window_scan = build_window_scan(L, beta, sigma, args.degree, args.weighted)

    default_json, default_fit = default_output_paths(
        args.summary, args.degree, float(fit["lmin"])
    )
    json_out = args.json_out or default_json
    fit_out = args.fit_out or default_fit

    payload = serialize_fit(fit, args.summary, window_scan)
    os.makedirs(os.path.dirname(json_out), exist_ok=True)
    with open(json_out, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    write_fit_table(fit_out, fit)

    methods, failures = build_comparison_methods(L, beta, sigma, args.weighted, fit)
    default_compare_json, default_compare_plot = default_compare_output_paths(
        args.summary, args.degree, float(fit["lmin"])
    )
    compare_json_out = args.compare_json_out or default_compare_json
    compare_plot_out = args.compare_plot_out or default_compare_plot
    compare_payload = serialize_comparison(args.summary, methods, failures)
    with open(compare_json_out, "w") as f:
        json.dump(compare_payload, f, indent=2, sort_keys=True)
        f.write("\n")
    make_comparison_plot(compare_plot_out, L, beta, sigma, methods)

    print(f"summary: {os.path.abspath(args.summary)}")
    print(
        f"fit: beta_pc(L) = sum_n c_n * (1/L)^n, degree={args.degree}, "
        f"L in [{fit['lmin']:.0f}, {fit['lmax']:.0f}], n={fit['n_points']}"
    )
    sigma_text = (
        f"{fit['beta_c_continuum_sigma']:.12f}"
        if fit["beta_c_continuum_sigma"] is not None
        else "n/a"
    )
    print(
        f"beta_c_continuum = {fit['beta_c_continuum']:.12f} "
        f"+- {sigma_text}"
    )
    print(f"rms_residual = {fit['rms_residual']:.12e}")
    print(
        "delta_vs_exact = "
        f"{fit['beta_c_continuum'] - EXACT_TRIANGULAR_ISING_BETA_C:+.12f} "
        f"(exact={EXACT_TRIANGULAR_ISING_BETA_C:.12f})"
    )
    print(f"wrote {json_out}")
    print(f"wrote {fit_out}")
    print(f"wrote {compare_json_out}")
    print(f"wrote {compare_plot_out}")
    print("comparison methods:")
    for method in methods:
        print_method_summary(method)
    if failures:
        print("method failures:")
        for failure in failures:
            print(f"  {failure['method']}: {failure['error']}")
    if window_scan:
        print_window_scan(window_scan)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())