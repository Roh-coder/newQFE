#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[2]
STORY_ROOT = Path(__file__).resolve().parent
FIGURES_DIR = STORY_ROOT / "figures"
OUT_ROOT = STORY_ROOT / "gc_three_pass_demo"
K_FROM_CONT_ROOT = REPO_ROOT / "K_from_continuum"
if str(K_FROM_CONT_ROOT) not in sys.path:
    sys.path.insert(0, str(K_FROM_CONT_ROOT))
if str(K_FROM_CONT_ROOT / "lib") not in sys.path:
    sys.path.insert(0, str(K_FROM_CONT_ROOT / "lib"))

import mc_engine  # noqa: E402
from workflow_common import exact_triangular_ising_beta  # noqa: E402


SIZES_DEFAULT = [8, 16, 24, 32, 48, 64]
CASE_DEFS = {
    "iso111": {
        "title": "GC three-pass beta finder: square lattice, k=(1,1,1)",
        "k1": 1.0,
        "k2": 1.0,
        "k3": 1.0,
        "figure": FIGURES_DIR / "08_gc_three_pass_iso111.png",
        "fss_figure": FIGURES_DIR / "10_gc_three_pass_iso111_fss.png",
    },
    "r12_r22": {
        "title": "GC three-pass beta finder: square lattice, k=(2,2,1)",
        "k1": 2.0,
        "k2": 2.0,
        "k3": 1.0,
        "figure": FIGURES_DIR / "09_gc_three_pass_r12_r22.png",
        "fss_figure": FIGURES_DIR / "11_gc_three_pass_r12_r22_fss.png",
    },
}
PASS_MARKERS = {0: "o", 1: "s", 2: "D", 3: "^"}
PASS_LABELS = {0: "coarse pass", 1: "refine pass 1", 2: "refine pass 2", 3: "refine pass 3"}


def json_default(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_json(path: Path, payload: dict[str, Any]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=json_default)
        handle.write("\n")


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object in {path}")
    return value


def beta_window(exact_beta: float) -> tuple[float, float]:
    half_width = max(0.05, 0.25 * exact_beta)
    return max(0.01, exact_beta - half_width), exact_beta + half_width


def cache_path(case_name: str, size: int) -> Path:
    return OUT_ROOT / case_name / f"L{size:03d}_betac_gc_scan.json"


def clone_progress_state(state: dict[str, Any]) -> dict[str, Any]:
    gc_params = state.get("gc_params")
    return {
        "pass_num": int(state["pass_num"]),
        "all_betas": [float(value) for value in state.get("all_betas", [])],
        "all_chis": [float(value) for value in state.get("all_chis", [])],
        "all_chi_errs": [float(value) for value in state.get("all_chi_errs", [])],
        "pass_ids": [int(value) for value in state.get("pass_ids", [])],
        "gc_params": [float(value) for value in gc_params] if gc_params is not None else None,
        "beta_estimate": float(state["beta_estimate"]) if state.get("beta_estimate") is not None else None,
        "traj_done": int(state.get("traj_done", 0)),
    }


def run_scan(
    *,
    case_name: str,
    size: int,
    n_traj_coarse: int,
    n_traj_fine: int,
    n_coarse: int,
    n_refine: int,
    force: bool,
) -> dict[str, Any]:
    path = cache_path(case_name, size)
    if path.is_file() and not force:
        return load_json(path)

    case = CASE_DEFS[case_name]
    couplings = {"k1": case["k1"], "k2": case["k2"], "k3": case["k3"]}
    exact_beta = float(exact_triangular_ising_beta(couplings))
    beta_lo, beta_hi = beta_window(exact_beta)
    snapshots: list[dict[str, Any]] = []

    def progress_cb(state: dict[str, Any]) -> None:
        if state.get("scan_progress") is None:
            snapshots.append(clone_progress_state(state))

    beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis, scan_chi_errs = mc_engine.find_beta_c(
        str(K_FROM_CONT_ROOT / "bin" / "ising_tri_twisted_parallelogram"),
        size,
        size,
        0,
        0,
        couplings["k1"],
        couplings["k2"],
        couplings["k3"],
        beta_lo,
        beta_hi,
        n_coarse=n_coarse,
        n_refine=n_refine,
        n_refine2=n_refine,
        n_refine3=n_refine,
        n_traj_coarse=n_traj_coarse,
        n_traj_fine=n_traj_fine,
        max_shifts=4,
        progress_cb=progress_cb,
        jackknife=True,
        data_dir=str(OUT_ROOT / case_name / f"L{size:03d}" / "scan_data"),
    )

    payload = {
        "case_name": case_name,
        "title": case["title"],
        "L": size,
        "Lx": size,
        "Ly": size,
        "Tx": 0,
        "Ty": 0,
        "k1": couplings["k1"],
        "k2": couplings["k2"],
        "k3": couplings["k3"],
        "r1": couplings["k1"] / couplings["k3"],
        "r2": couplings["k2"] / couplings["k3"],
        "exact_beta": exact_beta,
        "beta_lo": beta_lo,
        "beta_hi": beta_hi,
        "beta_c": float(beta_c),
        "beta_c_sigma": float(beta_c_sigma),
        "chi_peak": float(chi_peak),
        "scan_betas": [float(value) for value in scan_betas],
        "scan_chis": [float(value) for value in scan_chis],
        "scan_chi_errs": [float(value) for value in scan_chi_errs],
        "pass_ids": snapshots[-1]["pass_ids"] if snapshots else [],
        "pass_snapshots": snapshots,
        "n_traj_coarse": int(n_traj_coarse),
        "n_traj_fine": int(n_traj_fine),
        "n_coarse": int(n_coarse),
        "n_refine": int(n_refine),
    }
    save_json(path, payload)
    return payload


def gc_curve(grid: np.ndarray, gc_params: list[float] | None) -> np.ndarray | None:
    if gc_params is None:
        return None
    mu, sigma, skew, kurt, amplitude = gc_params
    z = (grid - mu) / sigma
    h3 = z ** 3 - 3.0 * z
    h4 = z ** 4 - 6.0 * z ** 2 + 3.0
    return amplitude * np.exp(-0.5 * z ** 2) * (1.0 + (skew / 6.0) * h3 + (kurt / 24.0) * h4)


def fit_sigmas(y_sigmas: np.ndarray) -> np.ndarray:
    sigma = np.asarray(y_sigmas, dtype=float)
    positive = sigma > 0.0
    if np.any(positive):
        fallback_sigma = float(np.min(sigma[positive]))
        return np.where(positive, sigma, fallback_sigma)
    return np.ones_like(sigma)


def weighted_linear_fit(x_values: np.ndarray, y_values: np.ndarray, y_sigmas: np.ndarray) -> tuple[float, float, np.ndarray]:
    sigma = fit_sigmas(y_sigmas)

    design = np.column_stack((np.ones_like(x_values), x_values))
    weights = 1.0 / np.square(sigma)
    xtwx = design.T @ (weights[:, None] * design)
    xtwy = design.T @ (weights * y_values)
    coeffs = np.linalg.solve(xtwx, xtwy)
    covariance = np.linalg.inv(xtwx)
    intercept, slope = coeffs
    return float(intercept), float(slope), covariance


def linear_fit_chi_squared(
    x_values: np.ndarray,
    y_values: np.ndarray,
    y_sigmas: np.ndarray,
    intercept: float,
    slope: float,
) -> tuple[float, int, float]:
    sigma = fit_sigmas(y_sigmas)
    residuals = (y_values - (intercept + slope * x_values)) / sigma
    chi_squared = float(np.sum(np.square(residuals)))
    dof = max(int(x_values.size) - 2, 0)
    reduced_chi_squared = chi_squared / dof if dof > 0 else float("nan")
    return chi_squared, dof, float(reduced_chi_squared)


def plot_case(case_name: str, payloads: list[dict[str, Any]]) -> None:
    case = CASE_DEFS[case_name]
    exact_beta = float(payloads[0]["exact_beta"])
    payloads = sorted(payloads, key=lambda payload: int(payload["L"]))
    colors = plt.cm.viridis(np.linspace(0.10, 0.95, len(payloads)))

    fig, ax = plt.subplots(figsize=(12.2, 7.4), constrained_layout=True)
    size_handles: list[Line2D] = []
    for index, (payload, color) in enumerate(zip(payloads, colors)):
        betas = np.asarray(payload["scan_betas"], dtype=float)
        chis = np.asarray(payload["scan_chis"], dtype=float)
        chi_errs = np.asarray(payload.get("scan_chi_errs", []), dtype=float)
        pass_ids = np.asarray(payload["pass_ids"], dtype=int)
        if chi_errs.size != betas.size:
            chi_errs = np.zeros_like(chis)
        order = np.argsort(betas)
        betas = betas[order]
        chis = chis[order]
        chi_errs = chi_errs[order]
        pass_ids = pass_ids[order]
        peak = float(np.max(chis))
        normalized = chis / peak
        normalized_errs = chi_errs / peak

        for pass_num in range(4):
            mask = pass_ids == pass_num
            if np.any(mask):
                if np.any(normalized_errs[mask] > 0.0):
                    ax.errorbar(
                        betas[mask],
                        normalized[mask],
                        yerr=normalized_errs[mask],
                        fmt=PASS_MARKERS[pass_num],
                        linestyle="None",
                        color=color,
                        ecolor=color,
                        elinewidth=0.9,
                        capsize=2.0,
                        markerfacecolor=color,
                        markeredgecolor="#111827",
                        markeredgewidth=0.4,
                        markersize=4.3,
                        alpha=0.95,
                    )
                else:
                    ax.plot(
                        betas[mask],
                        normalized[mask],
                        linestyle="None",
                        marker=PASS_MARKERS[pass_num],
                        markerfacecolor=color,
                        markeredgecolor="#111827",
                        markeredgewidth=0.4,
                        markersize=4.3,
                        alpha=0.95,
                    )

        final_gc = payload["pass_snapshots"][-1]["gc_params"] if payload["pass_snapshots"] else None
        grid = None
        curve_norm = None
        if final_gc is not None:
            grid = np.linspace(float(np.min(betas)), float(np.max(betas)), 500)
            curve = gc_curve(grid, final_gc)
            if curve is not None and np.max(curve) > 0.0:
                curve_norm = np.maximum(curve / peak, 0.0)
                ax.plot(grid, curve_norm, color=color, linewidth=2.1)

        beta_c = float(payload["beta_c"])
        beta_sigma = float(payload["beta_c_sigma"])
        if grid is not None and curve_norm is not None:
            star_y = float(np.interp(beta_c, grid, curve_norm))
        else:
            star_y = float(np.interp(beta_c, betas, normalized))
        ax.vlines(beta_c, 0.0, star_y, color=color, linewidth=1.1, linestyle=":", alpha=0.8)
        ax.scatter([beta_c], [star_y], marker="*", s=100, color=color, edgecolor="#111827", linewidth=0.4, zorder=5)
        size_handles.append(
            Line2D([0], [0], color=color, linewidth=2.1, label=f"L={payload['L']}  beta_c={beta_c:.5f}")
        )

    ax.axvline(exact_beta, color="#b91c1c", linestyle="--", linewidth=2.0)
    ax.set_title(case["title"])
    ax.set_xlabel("beta")
    ax.set_ylabel("normalized susceptibility")
    ax.set_ylim(-0.02, 1.14)
    ax.grid(True, alpha=0.22)
    ax.text(
        0.015,
        0.985,
        (
            f"exact continuum coupling = {exact_beta:.8f}\n"
            "each curve is normalized by its own peak\n"
            "markers encode the GC beta-finder passes\n"
            "points show 1 sigma susceptibility errors\n"
            "solid curve = final Gram-Charlier fit"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )
    legend_handles = [
        Line2D([0], [0], color="#2563eb", linewidth=2.1, label="final GC fit"),
        Line2D([0], [0], color="#b91c1c", linestyle="--", linewidth=2.0, label="exact continuum coupling"),
        Line2D([0], [0], linestyle="None", marker="*", markerfacecolor="#94a3b8", markeredgecolor="#111827", markeredgewidth=0.4, markersize=8, label="fitted beta_c"),
    ]
    for pass_num in range(4):
        legend_handles.append(
            Line2D(
                [0],
                [0],
                linestyle="None",
                marker=PASS_MARKERS[pass_num],
                markerfacecolor="#94a3b8",
                markeredgecolor="#111827",
                markeredgewidth=0.4,
                markersize=5,
                label=PASS_LABELS[pass_num],
            )
        )
    pass_legend = ax.legend(handles=legend_handles, loc="lower left", frameon=True, fontsize=9)
    ax.add_artist(pass_legend)
    ax.legend(handles=size_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=8.5, title="lattice size")
    ensure_dir(FIGURES_DIR)
    fig.savefig(case["figure"], dpi=220)
    plt.close(fig)


def plot_beta_fss(case_name: str, payloads: list[dict[str, Any]]) -> None:
    case = CASE_DEFS[case_name]
    exact_beta = float(payloads[0]["exact_beta"])
    payloads = sorted(payloads, key=lambda payload: int(payload["L"]))
    colors = plt.cm.viridis(np.linspace(0.10, 0.95, len(payloads)))

    lengths = np.asarray([int(payload["L"]) for payload in payloads], dtype=float)
    x_values = 1.0 / lengths
    beta_values = np.asarray([float(payload["beta_c"]) for payload in payloads], dtype=float)
    beta_sigmas = np.asarray([float(payload["beta_c_sigma"]) for payload in payloads], dtype=float)

    beta_inf, slope, covariance = weighted_linear_fit(x_values, beta_values, beta_sigmas)
    beta_inf_sigma = math.sqrt(max(float(covariance[0, 0]), 0.0))
    chi_squared, dof, reduced_chi_squared = linear_fit_chi_squared(x_values, beta_values, beta_sigmas, beta_inf, slope)
    x_fit = np.linspace(0.0, float(np.max(x_values)) * 1.06, 300)
    y_fit = beta_inf + slope * x_fit

    fig, ax = plt.subplots(figsize=(9.2, 5.8), constrained_layout=True)
    ax.plot(x_fit, y_fit, color="#2563eb", linewidth=2.0, label=r"weighted fit: $\beta_c(L)=\beta_\infty + a/L$")
    ax.plot(x_values, beta_values, color="#94a3b8", linewidth=1.1, alpha=0.7)

    for payload, x_value, beta_value, beta_sigma, color in zip(payloads, x_values, beta_values, beta_sigmas, colors):
        ax.errorbar(
            [x_value],
            [beta_value],
            yerr=[beta_sigma],
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=1.1,
            capsize=3.0,
            markeredgecolor="#111827",
            markeredgewidth=0.4,
            markersize=6.0,
            zorder=4,
        )
        ax.annotate(f"L={payload['L']}", (x_value, beta_value), textcoords="offset points", xytext=(6, 4), fontsize=8, color="#111827")

    ax.errorbar(
        [0.0],
        [beta_inf],
        yerr=[beta_inf_sigma],
        fmt="*",
        color="#1d4ed8",
        ecolor="#1d4ed8",
        elinewidth=1.2,
        capsize=3.5,
        markersize=12.0,
        markeredgecolor="#111827",
        markeredgewidth=0.5,
        zorder=5,
        label=r"fit extrapolation $\beta_\infty$",
    )
    ax.axhline(exact_beta, color="#b91c1c", linestyle="--", linewidth=2.0, label="exact continuum coupling")
    ax.axvline(0.0, color="#94a3b8", linestyle=":", linewidth=1.0)

    ax.set_title(f"Beta finite-size scaling: {case['title'].split(': ', 1)[1]}")
    ax.set_xlabel(r"$1 / L$")
    ax.set_ylabel(r"$\beta_c(L)$")
    ax.set_xlim(-0.004, float(np.max(x_values)) * 1.10)
    y_all = np.concatenate((beta_values, np.asarray([beta_inf, exact_beta])))
    y_pad = max(0.0025, 0.08 * float(np.max(y_all) - np.min(y_all)))
    ax.set_ylim(float(np.min(y_all)) - y_pad, float(np.max(y_all)) + y_pad)
    ax.grid(True, alpha=0.22)
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    ax.text(
        0.98,
        0.98,
        (
            f"exact continuum coupling = {exact_beta:.8f}\n"
            f"weighted extrapolation = {beta_inf:.8f} +/- {beta_inf_sigma:.6f}\n"
            f"chi^2 = {chi_squared:.3f},  chi^2/dof = {reduced_chi_squared:.3f} ({dof} dof)\n"
            f"fit form: beta_c(L) = beta_inf + a / L\n"
            "error bars are the jackknife beta_c uncertainties"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )

    ensure_dir(FIGURES_DIR)
    fig.savefig(case["fss_figure"], dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run and plot GC three-pass beta-finder scans for selected coupling points.")
    parser.add_argument("--sizes", type=int, nargs="+", default=SIZES_DEFAULT)
    parser.add_argument("--cases", nargs="+", choices=tuple(CASE_DEFS), default=list(CASE_DEFS))
    parser.add_argument("--n-traj-coarse", type=int, default=2000)
    parser.add_argument("--n-traj-fine", type=int, default=6000)
    parser.add_argument("--n-coarse", type=int, default=11)
    parser.add_argument("--n-refine", type=int, default=5)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(OUT_ROOT)
    ensure_dir(FIGURES_DIR)

    for case_name in args.cases:
        payloads: list[dict[str, Any]] = []
        print(f"\n=== case {case_name} ===", flush=True)
        for size in args.sizes:
            payload = run_scan(
                case_name=case_name,
                size=size,
                n_traj_coarse=int(args.n_traj_coarse),
                n_traj_fine=int(args.n_traj_fine),
                n_coarse=int(args.n_coarse),
                n_refine=int(args.n_refine),
                force=args.force,
            )
            payloads.append(payload)
            print(
                f"L={size:>3d} beta_c={payload['beta_c']:.8f} +/- {payload['beta_c_sigma']:.3g} "
                f"exact={payload['exact_beta']:.8f}",
                flush=True,
            )
        plot_case(case_name, payloads)
        plot_beta_fss(case_name, payloads)
        print(f"wrote {CASE_DEFS[case_name]['figure']}", flush=True)
        print(f"wrote {CASE_DEFS[case_name]['fss_figure']}", flush=True)


if __name__ == "__main__":
    main()