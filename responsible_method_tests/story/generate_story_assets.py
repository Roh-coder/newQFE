#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
STORY_ROOT = Path(__file__).resolve().parent
FIGURES_DIR = STORY_ROOT / "figures"
ACUTE456_MATCHED_VOLUME_ROOT = REPO_ROOT / "responsible_method_tests" / "results" / "raw_manifold_fss_acute_456_matched_volume_20260522"
ACUTE456_MATCHED_VOLUME_MANIFEST = ACUTE456_MATCHED_VOLUME_ROOT / "manifest_geometry_acute_456_matched_volume.json"

EXACT_TRIANGULAR_BETA = 0.274653072167027

COPIED_FIGURES = (
    (
        REPO_ROOT / "responsible_method_tests" / "results" / "raw_manifold_fss_acute_456_matched_volume_20260522" / "geometry_acute_456_matched_volume_twisted_vs_untwisted_common_grid.png",
        FIGURES_DIR / "01_acute456_matched_volume_common_grid.png",
    ),
    (
        REPO_ROOT / "responsible_method_tests" / "results" / "raw_manifold_fss_acute_456_matched_volume_20260522" / "geometry_acute_456_matched_volume_shared_coordinate_twisted_untwisted_fss.png",
        FIGURES_DIR / "02_acute456_shared_coordinate_fss.png",
    ),
    (
        REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_hi40000_notaylor" / "acute456_literal_sizes_blind_holdout_selected_models_logx.png",
        FIGURES_DIR / "05_acute456_blind_holdout_literal_sizes.png",
    ),
    (
        REPO_ROOT / "responsible_method_tests" / "stupid_method" / "iso111_base4_size4to100_step4" / "model_fit_comparison" / "holdout_L1000_pt002_selected_models_logx.png",
        FIGURES_DIR / "06_iso111_L1000_blind_holdout.png",
    ),
)


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_dat_rows(path: Path) -> list[dict[str, str]]:
    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            if columns is None:
                raise ValueError(f"missing '# columns:' header in {path}")
            values = line.split()
            if len(values) != len(columns):
                raise ValueError(f"expected {len(columns)} columns in {path}, got {len(values)}")
            rows.append(dict(zip(columns, values)))
    return rows


def beta_scan_volume(payload: dict) -> int:
    return abs(int(payload["Lx"]) * int(payload["Ly"]) + int(payload["Tx"]) * int(payload["Ty"]))


def copy_existing_figures() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for source, destination in COPIED_FIGURES:
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, destination)


def build_holdout_scan_plot() -> None:
    holdout_json = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_betac_hi40000" / "beta_scans" / "holdout_twisted" / "acute456_holdout_smallbase024_m-108_n60_beta_scan.json"
    payload = load_json(holdout_json)

    betas = list(map(float, payload["scan_betas"]))
    chis = list(map(float, payload["scan_chis"]))
    chi_errs = list(map(float, payload.get("scan_chi_errs", [])))
    points = sorted(zip(betas, chis, chi_errs if len(chi_errs) == len(betas) else [0.0] * len(betas)), key=lambda row: row[0])

    beta_sorted = [row[0] for row in points]
    chi_sorted = [row[1] for row in points]
    chi_err_sorted = [row[2] for row in points]

    beta_c = float(payload["beta_c"])
    beta_sigma = float(payload["beta_c_sigma"])
    chi_peak = float(payload["chi_peak"])

    fig, ax = plt.subplots(figsize=(9.6, 5.4), constrained_layout=True)
    if any(error > 0.0 for error in chi_err_sorted):
        ax.errorbar(beta_sorted, chi_sorted, yerr=chi_err_sorted, fmt="o", color="#111827", ms=4.0, capsize=2.5, alpha=0.9)
    else:
        ax.plot(beta_sorted, chi_sorted, "o", color="#111827", ms=4.0, alpha=0.9)
    ax.plot(beta_sorted, chi_sorted, color="#2563eb", linewidth=1.8, alpha=0.75, label="measured susceptibility")
    ax.axvline(EXACT_TRIANGULAR_BETA, color="#7c3aed", linestyle=":", linewidth=2.0, label="exact $\\beta_c = \\ln(3)/4$")
    ax.axvline(beta_c, color="#b91c1c", linestyle="--", linewidth=2.0, label="finder $\\beta_c$")
    if beta_sigma > 0.0:
        ax.axvspan(beta_c - beta_sigma, beta_c + beta_sigma, color="#fca5a5", alpha=0.35)
    ax.scatter([beta_c], [chi_peak], marker="*", s=180, color="#b91c1c", zorder=5)

    ax.set_title("Acute 4-5-6 holdout susceptibility scan from the beta_c finder")
    ax.set_xlabel("beta")
    ax.set_ylabel("magnetic susceptibility")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", frameon=True)
    ax.text(
        0.02,
        0.98,
        (
            f"holdout lattice = [144, 144, 72, 24]\n"
            f"beta_c = {beta_c:.8f} +/- {beta_sigma:.6f}\n"
            f"exact beta = {EXACT_TRIANGULAR_BETA:.8f}\n"
            f"peak chi = {chi_peak:.5f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )
    fig.savefig(FIGURES_DIR / "03_acute456_holdout_betac_scan.png", dpi=220)
    plt.close(fig)


def build_partial_trend_plot() -> None:
    holdout_json = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_betac_hi40000" / "beta_scans" / "holdout_twisted" / "acute456_holdout_smallbase024_m-108_n60_beta_scan.json"
    twisted_dir = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_betac_hi40000" / "beta_scans" / "twisted"

    holdout_payload = load_json(holdout_json)
    training_payloads = [load_json(path) for path in sorted(twisted_dir.glob("acute456_twisted_*_beta_scan.json"))]
    training_payloads.sort(key=beta_scan_volume)

    fig, ax = plt.subplots(figsize=(9.6, 5.4), constrained_layout=True)
    ax.axhline(EXACT_TRIANGULAR_BETA, color="#7c3aed", linestyle=":", linewidth=2.0, label="exact $\\beta_c = \\ln(3)/4$")

    x_values = []
    y_values = []
    labels = []
    for payload in training_payloads:
        volume = beta_scan_volume(payload)
        x_value = 1.0 / math.sqrt(volume)
        y_value = float(payload["beta_c"])
        y_sigma = float(payload["beta_c_sigma"])
        x_values.append(x_value)
        y_values.append(y_value)
        labels.append(f"{payload['Lx']}x{payload['Ly']}")
        ax.errorbar([x_value], [y_value], yerr=[y_sigma], fmt="o", color="#111827", ms=6, capsize=3)
        ax.text(x_value + 0.0008, y_value + 0.00025, labels[-1], fontsize=8, color="#111827")

    holdout_volume = beta_scan_volume(holdout_payload)
    holdout_x = 1.0 / math.sqrt(holdout_volume)
    holdout_y = float(holdout_payload["beta_c"])
    holdout_sigma = float(holdout_payload["beta_c_sigma"])
    ax.errorbar([holdout_x], [holdout_y], yerr=[holdout_sigma], fmt="*", color="#b91c1c", ms=16, capsize=3, label="held-out twisted lattice")
    ax.plot(x_values, y_values, color="#2563eb", linewidth=1.5, alpha=0.7, label="completed twisted training ladder")

    ax.set_title("Partial acute 4-5-6 beta_c trend from the stopped production rerun")
    ax.set_xlabel("$1 / \\sqrt{V}$")
    ax.set_ylabel("$\\beta_c$")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", frameon=True)
    ax.text(
        0.02,
        0.98,
        (
            "completed training sizes: 12, 24, 36, 48, 60, 72\n"
            f"holdout: beta_c = {holdout_y:.8f} +/- {holdout_sigma:.6f}\n"
            "96x96_t48x16 was in progress when the run was stopped"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )
    fig.savefig(FIGURES_DIR / "04_acute456_partial_betac_trend.png", dpi=220)
    plt.close(fig)


def build_peak_drift_plot() -> None:
    holdout_json = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_betac_hi40000" / "beta_scans" / "holdout_twisted" / "acute456_holdout_smallbase024_m-108_n60_beta_scan.json"
    twisted_dir = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_holdout20k_betac_hi40000" / "beta_scans" / "twisted"

    payloads = [load_json(path) for path in sorted(twisted_dir.glob("acute456_twisted_*_beta_scan.json"))]
    payloads.sort(key=beta_scan_volume)
    payloads.append(load_json(holdout_json))

    fig, ax = plt.subplots(figsize=(10.8, 7.8), constrained_layout=True)
    colors = plt.cm.viridis([index / max(len(payloads) - 1, 1) for index in range(len(payloads))])

    x_min = float("inf")
    x_max = float("-inf")
    for index, (payload, color) in enumerate(zip(payloads, colors)):
        points = sorted(zip(map(float, payload["scan_betas"]), map(float, payload["scan_chis"])), key=lambda row: row[0])
        betas = [row[0] for row in points]
        chis = [row[1] for row in points]
        chi_max = max(chis) if chis else 1.0
        y_offset = 1.18 * index
        normalized = [(value / chi_max) + y_offset for value in chis]
        beta_c = float(payload["beta_c"])
        label = f"{payload['Lx']}x{payload['Ly']}"
        if int(payload["Tx"]) != 0 or int(payload["Ty"]) != 0:
            label += f"_t{payload['Tx']}x{payload['Ty']}"
        if payload["label"].startswith("acute456_holdout"):
            label = f"holdout {label}"

        ax.plot(betas, normalized, "-o", color=color, linewidth=1.6, markersize=3.3)
        ax.scatter([beta_c], [1.02 + y_offset], marker="*", s=120, color=color, edgecolor="black", linewidth=0.4, zorder=5)
        ax.text(beta_c + 0.00055, 0.18 + y_offset, f"{label}\n$\\beta_c={beta_c:.4f}$", fontsize=8.5, color=color, va="bottom")

        x_min = min(x_min, min(betas))
        x_max = max(x_max, max(betas))

    ax.axvline(EXACT_TRIANGULAR_BETA, color="#b91c1c", linestyle="--", linewidth=2.0, label="exact continuum coupling")
    ax.set_xlim(max(0.248, x_min - 0.002), min(0.2795, x_max + 0.002))
    ax.set_ylim(-0.15, 1.18 * len(payloads))
    ax.set_xlabel("beta")
    ax.set_ylabel("normalized susceptibility (stacked by lattice size)")
    ax.set_title("Acute 4-5-6 twisted susceptibility peaks drift toward the continuum coupling")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper left", frameon=True)
    ax.text(
        0.985,
        0.02,
        "each curve is divided by its own peak susceptibility\ncompleted training sizes plus the held-out twisted lattice",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )
    fig.savefig(FIGURES_DIR / "07_acute456_twisted_peak_drift.png", dpi=220)
    plt.close(fig)


def build_twisted_modular_validation_plot() -> None:
    benchmark_manifest = load_json(ACUTE456_MATCHED_VOLUME_MANIFEST)
    twisted_manifest = load_json(Path(benchmark_manifest["methods"]["twisted"]))
    alignment = load_json(Path(twisted_manifest["modular_alignment"]))
    continuum_rows = load_dat_rows(Path(twisted_manifest["pointwise_continuum"]))
    modular_rows = load_dat_rows(Path(twisted_manifest["modular_aligned"]))
    scale_values = sorted(int(payload["scale"]) for payload in twisted_manifest["payloads"])
    fit_modes = sorted({str(row["fit_mode"]) for row in continuum_rows if int(row["is_origin"]) == 0})
    n_sizes_used_values = sorted({int(row["n_sizes_used"]) for row in continuum_rows if int(row["is_origin"]) == 0})

    modular_lookup: dict[int, float] = {}
    for row in modular_rows:
        if int(row["is_origin"]) != 0:
            continue
        value = float(row["modular_aligned"])
        if math.isfinite(value):
            modular_lookup[int(row["point_id"])] = value

    point_ids: list[int] = []
    nu_real: list[float] = []
    nu_imag: list[float] = []
    continuum_values: list[float] = []
    continuum_sigmas: list[float] = []
    modular_values: list[float] = []
    for row in continuum_rows:
        if int(row["is_origin"]) != 0:
            continue
        point_id = int(row["point_id"])
        modular_value = modular_lookup.get(point_id)
        if modular_value is None:
            continue
        point_ids.append(point_id)
        nu_real.append(float(row["nu_real"]))
        nu_imag.append(float(row["nu_imag"]))
        continuum_values.append(float(row["A"]))
        continuum_sigmas.append(float(row["sigma_A"]))
        modular_values.append(modular_value)

    modular_array = np.asarray(modular_values, dtype=float)
    continuum_array = np.asarray(continuum_values, dtype=float)
    sigma_array = np.asarray(continuum_sigmas, dtype=float)
    nu_real_array = np.asarray(nu_real, dtype=float)
    nu_imag_array = np.asarray(nu_imag, dtype=float)
    residual_array = continuum_array - modular_array
    z_array = residual_array / sigma_array

    correlation = float(np.corrcoef(modular_array, continuum_array)[0, 1])
    rmse = float(np.sqrt(np.mean(np.square(residual_array))))
    relative_rmse = rmse / float(np.sqrt(np.mean(np.square(modular_array))))
    mean_abs_residual = float(np.mean(np.abs(residual_array)))
    rms_z = float(np.sqrt(np.mean(np.square(z_array))))
    chi2_per_dof = float(alignment["chi2_per_dof"])
    alpha = float(alignment["alpha"])
    fit_mode_text = fit_modes[0] if len(fit_modes) == 1 else ", ".join(fit_modes)
    if len(n_sizes_used_values) == 1:
        n_sizes_text = str(n_sizes_used_values[0])
    else:
        n_sizes_text = f"{n_sizes_used_values[0]}-{n_sizes_used_values[-1]}"

    fig, (ax_scatter, ax_residual) = plt.subplots(1, 2, figsize=(12.2, 5.8), constrained_layout=True)

    sample_stride = max(len(point_ids) // 28, 1)
    sample_idx = np.arange(0, len(point_ids), sample_stride, dtype=int)
    ax_scatter.errorbar(
        modular_array[sample_idx],
        continuum_array[sample_idx],
        yerr=sigma_array[sample_idx],
        fmt="none",
        ecolor="#94a3b8",
        elinewidth=0.9,
        capsize=2.0,
        alpha=0.75,
        zorder=1,
    )
    scatter = ax_scatter.scatter(
        modular_array,
        continuum_array,
        c=nu_imag_array,
        cmap="viridis",
        s=32,
        edgecolors="#111827",
        linewidths=0.25,
        alpha=0.92,
        zorder=2,
    )
    value_min = float(min(np.min(modular_array), np.min(continuum_array)))
    value_max = float(max(np.max(modular_array), np.max(continuum_array)))
    value_pad = 0.04 * (value_max - value_min)
    ax_scatter.plot(
        [value_min - value_pad, value_max + value_pad],
        [value_min - value_pad, value_max + value_pad],
        linestyle="--",
        linewidth=1.8,
        color="#b91c1c",
        label="perfect agreement",
    )
    ax_scatter.set_xlim(value_min - value_pad, value_max + value_pad)
    ax_scatter.set_ylim(value_min - value_pad, value_max + value_pad)
    ax_scatter.set_xlabel("aligned modular correlator")
    ax_scatter.set_ylabel("twisted continuum correlator")
    ax_scatter.set_title("Pointwise agreement with the modular solution")
    ax_scatter.grid(True, alpha=0.22)
    ax_scatter.legend(loc="upper left", frameon=True)
    colorbar = fig.colorbar(scatter, ax=ax_scatter, shrink=0.92)
    colorbar.set_label(r"$\nu_{\rm imag}$")

    z_limit = max(1.5, float(np.max(np.abs(z_array))) * 1.05)
    residual_plot = ax_residual.scatter(
        nu_real_array,
        nu_imag_array,
        c=z_array,
        cmap="coolwarm",
        vmin=-z_limit,
        vmax=z_limit,
        s=46,
        edgecolors="none",
        alpha=0.95,
    )
    ax_residual.set_title("Residual map across the twisted torus unit cell")
    ax_residual.set_xlabel("nu_real")
    ax_residual.set_ylabel("nu_imag")
    ax_residual.set_aspect("equal", adjustable="box")
    ax_residual.grid(True, alpha=0.22)
    residual_cbar = fig.colorbar(residual_plot, ax=ax_residual, shrink=0.92)
    residual_cbar.set_label(r"$(G_{\rm tw}-\alpha G_{\rm mod}) / \sigma_G$")
    ax_residual.text(
        0.98,
        0.02,
        (
            f"continuum fit = {fit_mode_text}\n"
            r"$G(s)=A + B/s + C/s^2$, plotted value $=A$" + "\n"
            f"scale range = {scale_values[0]}-{scale_values[-1]}, n_sizes_used = {n_sizes_text}\n"
            f"points = {len(point_ids)}\n"
            f"alpha = {alpha:.6f}\n"
            f"chi2/dof = {chi2_per_dof:.4f}\n"
            f"correlation = {correlation:.4f}\n"
            f"RMSE = {rmse:.5f}\n"
            f"relative RMSE = {relative_rmse:.3%}\n"
            f"mean |delta| = {mean_abs_residual:.5f}\n"
            f"RMS z = {rms_z:.3f}"
        ),
        transform=ax_residual.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.7,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )

    fig.suptitle("Acute 4-5-6 twisted continuum correlator vs aligned modular solution", fontsize=15)
    fig.savefig(FIGURES_DIR / "12_acute456_twisted_vs_modular.png", dpi=220)
    plt.close(fig)


def main() -> None:
    copy_existing_figures()
    build_holdout_scan_plot()
    build_partial_trend_plot()
    build_peak_drift_plot()
    build_twisted_modular_validation_plot()


if __name__ == "__main__":
    main()