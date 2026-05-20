#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from workflow_common import (
    check_required_sections,
    ensure_dir,
    load_json,
    load_payload_from_dat,
    log,
    metadata_path_for_data,
    resolve_path,
    sample_directional_channels,
    save_json,
    timestamp,
    write_dat,
)

_HERE = os.path.dirname(os.path.abspath(__file__))

# Edit this when running directly from an IDE without CLI args.
DEFAULT_CONFIG_PATH = "configs/score_example.json"


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


def _ratio_key(value: float) -> float:
    return round(float(value), 12)


def _power_model(L: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
    return A + B * np.power(1.0 / np.asarray(L, dtype=float), C)


def _fit_fixed_c(L: np.ndarray, y: np.ndarray, sigma: np.ndarray, C_fixed: float) -> tuple[float, float, float, float, str]:
    x = np.power(1.0 / L, C_fixed)
    w = np.ones_like(x)
    mask = np.isfinite(sigma) & (sigma > 0)
    if np.any(mask):
        w = np.where(mask, 1.0 / (sigma * sigma), 1.0)

    X = np.column_stack([np.ones_like(x), x])
    XtW = X.T * w
    try:
        cov = np.linalg.inv(XtW @ X)
        beta = cov @ (XtW @ y)
        A = float(beta[0])
        B = float(beta[1])
        sigma_A = float(np.sqrt(max(cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        A = float(np.mean(y))
        B = 0.0
        sigma_A = np.nan
    return A, sigma_A, B, float(C_fixed), "fixed_C"


def _robust_sigma(sigma: np.ndarray) -> np.ndarray:
    """Fill non-finite or non-positive sigmas with the median of the finite,
    positive ones; if none are usable, return all-ones (i.e. unweighted).

    This avoids the all-or-nothing weighting in the previous implementation,
    where a single missing sigma silently disabled weighting altogether.
    """
    sigma = np.asarray(sigma, dtype=float)
    good = np.isfinite(sigma) & (sigma > 0.0)
    if not np.any(good):
        return np.ones_like(sigma)
    med = float(np.median(sigma[good]))
    out = np.where(good, sigma, med)
    return out


def _fit_power_continuum(
    L_values: np.ndarray,
    y_values: np.ndarray,
    sigma_values: np.ndarray,
    c_min: float,
    c_max: float,
    c_initial: float,
    min_sizes_for_free_C: int = 8,
    boundary_tol: float = 1e-3,
    sigma_A_factor: float = 5.0,
    A_extreme_factor: float = 10.0,
) -> tuple[float, float, float, float, int, str]:
    """Fit y(L) = A + B*(1/L)^C with guards against degenerate fits.

    Strategy:
      * always compute a baseline weighted fit at fixed C = 1.0;
      * additionally attempt a free-C fit only when there are enough sizes
        (``n >= min_sizes_for_free_C``) — otherwise the 3-parameter model
        is unidentifiable and ``curve_fit`` parks C at a bound;
      * accept the free-C fit only if it passes sanity guards
        (C not pinned, sigma_A not exploding, |A| not absurd).
    Otherwise fall back to the fixed-C result.
    """
    mask = np.isfinite(L_values) & np.isfinite(y_values)
    L = np.asarray(L_values[mask], dtype=float)
    y = np.asarray(y_values[mask], dtype=float)
    sigma = np.asarray(sigma_values[mask], dtype=float)
    n = len(L)
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan, 0, "none"
    if n == 1:
        return float(y[0]), float(abs(sigma[0])) if np.isfinite(sigma[0]) else np.nan, 0.0, 1.0, 1, "single_point"

    sigma_fit = _robust_sigma(sigma)

    # Baseline: fixed C = 1.0 (always well-conditioned for n >= 2).
    A_fc, sA_fc, B_fc, C_fc, mode_fc = _fit_fixed_c(L, y, sigma_fit, 1.0)

    if n < max(3, min_sizes_for_free_C):
        return A_fc, sA_fc, B_fc, C_fc, n, "fixed_C_n<min"

    A0 = float(A_fc) if np.isfinite(A_fc) else float(np.mean(y[-2:]))
    B0 = float(B_fc) if np.isfinite(B_fc) else float(y[0] - A0)
    C0 = float(np.clip(c_initial, c_min, c_max))

    try:
        popt, pcov = curve_fit(
            _power_model,
            L,
            y,
            p0=[A0, B0, C0],
            sigma=sigma_fit,
            absolute_sigma=True,
            bounds=([-np.inf, -np.inf, c_min], [np.inf, np.inf, c_max]),
            maxfev=12000,
        )
        A = float(popt[0])
        B = float(popt[1])
        C = float(popt[2])
        sigma_A = (
            float(np.sqrt(max(float(pcov[0, 0]), 0.0)))
            if np.isfinite(pcov[0, 0]) else np.nan
        )
    except Exception:
        return A_fc, sA_fc, B_fc, C_fc, n, "fallback_fixed_C"

    # Sanity guards — reject obviously broken free-C fits and fall back.
    y_scale = float(np.max(np.abs(y))) if y.size else 1.0
    pinned_low = abs(C - c_min) < boundary_tol
    pinned_high = abs(C - c_max) < boundary_tol
    sigma_A_blown = np.isfinite(sigma_A) and sigma_A > sigma_A_factor * (abs(A) + 1e-6)
    A_extreme = (not np.isfinite(A)) or abs(A) > A_extreme_factor * (y_scale + 1e-12)

    if pinned_low or pinned_high or sigma_A_blown or A_extreme:
        reason = (
            "Cmin" if pinned_low else
            "Cmax" if pinned_high else
            "sigA" if sigma_A_blown else
            "Aext"
        )
        return A_fc, sA_fc, B_fc, C_fc, n, f"fallback_fixed_C_{reason}"

    return A, sigma_A, B, C, n, "power_fit"


def _draw_heatmap(ax, grid: np.ndarray, r1_values: list[float], r2_values: list[float],
                  title: str, cbar_label: str, cmap: str) -> None:
    x_lo = float(min(r1_values))
    x_hi = float(max(r1_values))
    y_lo = float(min(r2_values))
    y_hi = float(max(r2_values))
    if abs(x_hi - x_lo) < 1e-12:
        x_lo -= 0.5
        x_hi += 0.5
    if abs(y_hi - y_lo) < 1e-12:
        y_lo -= 0.5
        y_hi += 0.5
    extent = [x_lo, x_hi, y_lo, y_hi]
    ax.set_xlabel("k1/k3")
    ax.set_ylabel("k2/k3")
    ax.set_title(title)

    if np.any(np.isfinite(grid)):
        im = ax.imshow(grid, origin="lower", aspect="auto", extent=extent, cmap=cmap)
        jj, ii = np.unravel_index(np.nanargmin(grid), grid.shape)
        ax.plot(r1_values[ii], r2_values[jj], marker="*", markersize=12,
                markerfacecolor="white", markeredgecolor="k")
    else:
        im = ax.imshow(np.zeros_like(grid), origin="lower", aspect="auto", extent=extent,
                       cmap=cmap, vmin=0.0, vmax=1.0)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cbar_label)


def _save_heatmaps(r1_values: list[float], r2_values: list[float],
                   score_grid: np.ndarray, zscore_grid: np.ndarray,
                   comparison_data_dir: str, run_tag: str) -> None:
    if len(r1_values) == 0 or len(r2_values) == 0:
        return

    score_path = os.path.join(comparison_data_dir, "score_heatmap.png")
    zscore_path = os.path.join(comparison_data_dir, "zscore_heatmap.png")
    combined_path = os.path.join(comparison_data_dir, "score_zscore_heatmaps.png")

    fig1, ax1 = plt.subplots(figsize=(6.5, 5.3))
    _draw_heatmap(ax1, score_grid, r1_values, r2_values,
                  f"Score heatmap (tag={run_tag})", "score", "viridis")
    fig1.tight_layout()
    fig1.savefig(score_path, dpi=150)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(6.5, 5.3))
    _draw_heatmap(ax2, zscore_grid, r1_values, r2_values,
                  f"RMS z-score heatmap (tag={run_tag})", "z_rms", "magma")
    fig2.tight_layout()
    fig2.savefig(zscore_path, dpi=150)
    plt.close(fig2)

    fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(13.0, 5.3))
    _draw_heatmap(ax3, score_grid, r1_values, r2_values,
                  "Score", "score", "viridis")
    _draw_heatmap(ax4, zscore_grid, r1_values, r2_values,
                  "RMS z-score", "z_rms", "magma")
    fig3.suptitle(f"Score and z-score heatmaps (tag={run_tag})")
    fig3.tight_layout()
    fig3.savefig(combined_path, dpi=150)
    plt.close(fig3)


def _save_fss_plot_for_point(
    *,
    comparison_data_dir: str,
    run_tag: str,
    r1: float,
    r2: float,
    k_values: list[int],
    point_fss_channels: list[dict[str, Any]],
    dpi: int,
) -> str:
    n_rows = 3
    n_cols = max(len(k_values), 1)
    fig_w = max(4.0 * n_cols, 10.0)
    fig_h = max(3.3 * n_rows, 9.5)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    for channel in point_fss_channels:
        cycle = int(channel["cycle"])
        ik = int(channel["k_index"])
        ax = axes[cycle][ik]

        L_vals = np.asarray(channel["L_vals"], dtype=float)
        y_vals = np.asarray(channel["y_vals"], dtype=float)
        s_vals = np.asarray(channel["s_vals"], dtype=float)
        if L_vals.size > 0:
            order = np.argsort(L_vals)
            L_sorted = L_vals[order]
            y_sorted = y_vals[order]
            s_sorted = s_vals[order]
            x_sorted = 1.0 / L_sorted
            if np.any(np.isfinite(s_sorted)):
                yerr = np.where(np.isfinite(s_sorted), s_sorted, 0.0)
                ax.errorbar(x_sorted, y_sorted, yerr=yerr, fmt="o", color="#1f77b4", label="test sizes")
            else:
                ax.plot(x_sorted, y_sorted, "o", color="#1f77b4", label="test sizes")

        invL_ref = float(channel["invL_ref"])
        ref_value = float(channel["ref_value"])
        ref_sigma = float(channel["ref_sigma"])
        if np.isfinite(invL_ref) and np.isfinite(ref_value):
            if np.isfinite(ref_sigma):
                ax.errorbar([invL_ref], [ref_value], yerr=[ref_sigma], fmt="s", color="#d62728", label="reference")
            else:
                ax.plot([invL_ref], [ref_value], "s", color="#d62728", label="reference")

        A = float(channel["A"])
        sigma_A = float(channel["sigma_A"])
        B = float(channel["B"])
        C = float(channel["C"])
        if np.isfinite(A):
            if np.isfinite(sigma_A):
                ax.errorbar([0.0], [A], yerr=[sigma_A], fmt="*", markersize=12, color="#2ca02c", label="continuum A")
            else:
                ax.plot([0.0], [A], "*", markersize=12, color="#2ca02c", label="continuum A")

        if np.isfinite(A) and np.isfinite(B) and np.isfinite(C):
            x_candidates = [0.0, invL_ref]
            if L_vals.size > 0:
                x_candidates.extend((1.0 / L_vals).tolist())
            x_max = max(v for v in x_candidates if np.isfinite(v))
            x_fit = np.linspace(0.0, max(x_max, 1e-3), 200)
            y_fit = A + B * np.power(x_fit, C)
            ax.plot(x_fit, y_fit, "-", color="#2ca02c", alpha=0.9, label="fit")

        ax.set_title(f"cycle={cycle}  k={channel['k']}  t={channel['t']:.3f}")
        if cycle == n_rows - 1:
            ax.set_xlabel("1/L")
        if ik == 0:
            ax.set_ylabel("G")
        ax.grid(alpha=0.25)

    for cycle in range(n_rows):
        for ik in range(n_cols):
            ax = axes[cycle][ik]
            if not ax.has_data():
                ax.set_axis_off()

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        fig.legend(list(unique.values()), list(unique.keys()), loc="upper center", ncol=min(4, len(unique)))

    fig.suptitle(f"FSS channels with reference point  (tag={run_tag}, r1={r1:.6f}, r2={r2:.6f})")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])

    plot_dir = os.path.join(comparison_data_dir, "fss_plots")
    ensure_dir(plot_dir)
    out_png = os.path.join(plot_dir, f"fss_r1_{r1:.6f}_r2_{r2:.6f}.png")
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 3: fit each test channel to A + B*(1/L)^C and score test "
            "continuum channels against the large reference channels."
        )
    )
    parser.add_argument(
        "--config",
        type=str,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the scoring config JSON (defaults to DEFAULT_CONFIG_PATH)",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Optional override for run.tag",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = resolve_path(args.config, _HERE)
    cfg = load_json(config_path)

    check_required_sections(cfg, ["run", "paths", "input", "analysis"])
    if args.tag is not None:
        cfg["run"]["tag"] = args.tag

    _check_required_keys(cfg["run"], "run", ["tag"])
    _check_required_keys(cfg["paths"], "paths", ["results_root"])
    _check_required_keys(cfg["input"], "input", ["tests_manifest"])
    _check_required_keys(cfg["analysis"], "analysis", ["k_values", "k_denominator", "weighted", "power_fit"])
    _check_required_keys(cfg["analysis"]["power_fit"], "analysis.power_fit", ["c_min", "c_max", "c_initial"])

    reference_data_value = cfg["input"].get("reference_data_file") or cfg["input"].get("reference_payload")
    if reference_data_value is None:
        raise ValueError("input.reference_data_file is required (input.reference_payload is accepted for backward compatibility)")

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _HERE)
    run_root = os.path.join(results_root, tag)
    reference_data_dir = os.path.join(run_root, "reference_data")
    test_data_dir = os.path.join(run_root, "test_data")
    comparison_data_dir = os.path.join(run_root, "comparison_analysis_data")
    ensure_dir(reference_data_dir)
    ensure_dir(test_data_dir)
    ensure_dir(comparison_data_dir)

    reference_data_path = resolve_path(str(reference_data_value), _HERE)
    reference_metadata_path = resolve_path(
        str(cfg["input"].get("reference_metadata_file") or metadata_path_for_data(reference_data_path)),
        _HERE,
    )
    tests_manifest_path = resolve_path(str(cfg["input"]["tests_manifest"]), _HERE)

    reference_payload = load_payload_from_dat(reference_data_path, reference_metadata_path)
    tests_manifest = load_json(tests_manifest_path)
    tests_manifest_dir = os.path.dirname(tests_manifest_path)

    sizes = [int(v) for v in tests_manifest["sizes"]]
    k3 = float(tests_manifest["k3"])

    k_values = [int(v) for v in cfg["analysis"]["k_values"]]
    k_denom = int(cfg["analysis"]["k_denominator"])
    if k_denom <= 0:
        raise ValueError("analysis.k_denominator must be positive")
    fractions = np.array(k_values, dtype=float) / float(k_denom)

    weighted = bool(cfg["analysis"]["weighted"])
    c_min = float(cfg["analysis"]["power_fit"]["c_min"])
    c_max = float(cfg["analysis"]["power_fit"]["c_max"])
    c_initial = float(cfg["analysis"]["power_fit"]["c_initial"])
    min_sizes_for_free_C = int(cfg["analysis"]["power_fit"].get("min_sizes_for_free_C", 8))
    if c_max <= c_min:
        raise ValueError("analysis.power_fit requires c_max > c_min")
    if min_sizes_for_free_C < 3:
        raise ValueError("analysis.power_fit.min_sizes_for_free_C must be >= 3")

    fss_cfg = cfg["analysis"].get("fss_plots") or {}
    fss_enabled = bool(fss_cfg.get("enabled", True))
    fss_dpi = int(fss_cfg.get("dpi", 150))
    fss_max_points_raw = fss_cfg.get("max_points", None)
    fss_max_points = None if fss_max_points_raw in (None, "") else int(fss_max_points_raw)
    if fss_max_points is not None and fss_max_points < 1:
        raise ValueError("analysis.fss_plots.max_points must be >= 1 when set")

    # Build a lookup for payload file by (L, r1, r2).
    payload_map: dict[tuple[int, float, float], tuple[str, str]] = {}
    valid_jobs = [job for job in tests_manifest["jobs"] if job["status"] in {"ok", "skip"}]
    for job in valid_jobs:
        key = (int(job["L"]), _ratio_key(job["r1"]), _ratio_key(job["r2"]))
        data_path_raw = str(job.get("data_path") or job.get("payload_path") or "")
        data_path = resolve_path(data_path_raw, tests_manifest_dir) if data_path_raw else ""
        if data_path == "":
            continue
        metadata_raw = str(job.get("metadata_path") or metadata_path_for_data(data_path))
        metadata_path = resolve_path(metadata_raw, tests_manifest_dir)
        payload_map[key] = (data_path, metadata_path)

    r1_values = sorted({_ratio_key(job["r1"]) for job in valid_jobs})
    r2_values = sorted({_ratio_key(job["r2"]) for job in valid_jobs})

    G_ref, sG_ref = sample_directional_channels(reference_payload, fractions)

    ref_rows: list[list[Any]] = []
    for cycle in range(3):
        for ik, kval in enumerate(k_values):
            ref_rows.append(
                [
                    cycle,
                    int(kval),
                    float(fractions[ik]),
                    float(G_ref[cycle, ik]),
                    float(sG_ref[cycle, ik]),
                ]
            )

    raw_test_rows: list[list[Any]] = []
    continuum_rows: list[list[Any]] = []
    compare_rows: list[list[Any]] = []
    score_rows: list[list[Any]] = []
    fss_rows: list[list[Any]] = []
    fss_plot_index_rows: list[list[Any]] = []

    payload_cache: dict[tuple[str, str], dict[str, Any]] = {}
    score_grid = np.full((len(r2_values), len(r1_values)), np.nan)
    zscore_grid = np.full_like(score_grid, np.nan)
    r1_index = {_ratio_key(v): i for i, v in enumerate(r1_values)}
    r2_index = {_ratio_key(v): j for j, v in enumerate(r2_values)}
    L_ref_eff = int(reference_payload.get("L", max(abs(int(reference_payload["Lx"])), abs(int(reference_payload["Ly"])))))
    invL_ref = (1.0 / float(L_ref_eff)) if L_ref_eff > 0 else np.nan
    fss_saved_count = 0

    for r2 in r2_values:
        for r1 in r1_values:
            channels_by_size: dict[int, tuple[np.ndarray, np.ndarray, dict[str, Any]]] = {}
            missing_sizes: list[int] = []
            point_fss_channels: list[dict[str, Any]] = []
            for L in sizes:
                key = (int(L), r1, r2)
                if key not in payload_map:
                    missing_sizes.append(int(L))
                    continue
                data_path, metadata_path = payload_map[key]
                cache_key = (data_path, metadata_path)
                if cache_key not in payload_cache:
                    payload_cache[cache_key] = load_payload_from_dat(data_path, metadata_path)
                payload = payload_cache[cache_key]
                Gt, sGt = sample_directional_channels(payload, fractions)
                channels_by_size[int(L)] = (Gt, sGt, payload)

                for cycle in range(3):
                    for ik, kval in enumerate(k_values):
                        raw_test_rows.append(
                            [
                                float(r1),
                                float(r2),
                                int(L),
                                int(payload["Lx"]),
                                int(payload["Ly"]),
                                int(payload["Tx"]),
                                int(payload["Ty"]),
                                float(payload["beta_c"]),
                                cycle,
                                int(kval),
                                float(fractions[ik]),
                                float(Gt[cycle, ik]),
                                float(sGt[cycle, ik]),
                            ]
                        )

            score = 0.0
            n_score = 0
            z2_sum = 0.0
            n_z = 0

            for cycle in range(3):
                for ik, kval in enumerate(k_values):
                    L_vals: list[float] = []
                    y_vals: list[float] = []
                    s_vals: list[float] = []
                    for L in sizes:
                        if L not in channels_by_size:
                            continue
                        Gt, sGt, _ = channels_by_size[L]
                        L_vals.append(float(L))
                        y_vals.append(float(Gt[cycle, ik]))
                        s_vals.append(float(sGt[cycle, ik]))

                    A, sigma_A, B, C, n_used, fit_mode = _fit_power_continuum(
                        np.array(L_vals, dtype=float),
                        np.array(y_vals, dtype=float),
                        np.array(s_vals, dtype=float),
                        c_min=c_min,
                        c_max=c_max,
                        c_initial=c_initial,
                        min_sizes_for_free_C=min_sizes_for_free_C,
                    )

                    continuum_rows.append(
                        [
                            float(r1),
                            float(r2),
                            cycle,
                            int(kval),
                            float(fractions[ik]),
                            A,
                            sigma_A,
                            B,
                            C,
                            n_used,
                            fit_mode,
                        ]
                    )

                    ref_value = float(G_ref[cycle, ik])
                    ref_sigma = float(sG_ref[cycle, ik])

                    for L_i, y_i, s_i in zip(L_vals, y_vals, s_vals):
                        invL_i = (1.0 / float(L_i)) if float(L_i) > 0.0 else np.nan
                        fss_rows.append([
                            float(r1), float(r2), cycle, int(kval), float(fractions[ik]),
                            "test", float(L_i), invL_i, float(y_i), float(s_i), "", np.nan, np.nan,
                        ])
                    fss_rows.append([
                        float(r1), float(r2), cycle, int(kval), float(fractions[ik]),
                        "reference", float(L_ref_eff), float(invL_ref), ref_value, ref_sigma, "", np.nan, np.nan,
                    ])
                    fss_rows.append([
                        float(r1), float(r2), cycle, int(kval), float(fractions[ik]),
                        "continuum_fit", 0.0, 0.0, float(A), float(sigma_A), fit_mode, float(B), float(C),
                    ])

                    point_fss_channels.append({
                        "cycle": cycle,
                        "k_index": ik,
                        "k": int(kval),
                        "t": float(fractions[ik]),
                        "L_vals": [float(v) for v in L_vals],
                        "y_vals": [float(v) for v in y_vals],
                        "s_vals": [float(v) for v in s_vals],
                        "A": float(A),
                        "sigma_A": float(sigma_A),
                        "B": float(B),
                        "C": float(C),
                        "fit_mode": fit_mode,
                        "ref_value": ref_value,
                        "ref_sigma": ref_sigma,
                        "invL_ref": float(invL_ref),
                    })

                    diff = np.nan
                    var = np.nan
                    z = np.nan

                    if np.isfinite(A) and np.isfinite(ref_value):
                        diff = A - ref_value
                        if np.isfinite(sigma_A) and np.isfinite(ref_sigma):
                            var = sigma_A * sigma_A + ref_sigma * ref_sigma
                        weight = 1.0
                        if weighted and np.isfinite(var) and var > 0.0:
                            weight = 1.0 / var
                        # Clip per-channel residual so a single pathological
                        # extrapolation cannot dominate the cell score.
                        diff2 = diff * diff
                        if np.isfinite(var) and var > 0.0:
                            diff2 = min(diff2, 25.0 * var)
                        score += weight * diff2
                        n_score += 1

                        if np.isfinite(var) and var > 0.0:
                            z = diff / np.sqrt(var)
                            z2_sum += z * z
                            n_z += 1

                    compare_rows.append(
                        [
                            float(r1),
                            float(r2),
                            cycle,
                            int(kval),
                            float(fractions[ik]),
                            A,
                            sigma_A,
                            ref_value,
                            ref_sigma,
                            diff,
                            var,
                            z,
                            fit_mode,
                        ]
                    )

            score_value = float(score) if n_score > 0 else np.nan
            z_rms = float(np.sqrt(z2_sum / n_z)) if n_z > 0 else np.nan
            score_rows.append(
                [
                    float(r1),
                    float(r2),
                    score_value,
                    int(n_score),
                    z_rms,
                    int(n_z),
                    len(missing_sizes),
                ]
            )
            score_grid[r2_index[_ratio_key(r2)], r1_index[_ratio_key(r1)]] = score_value
            zscore_grid[r2_index[_ratio_key(r2)], r1_index[_ratio_key(r1)]] = z_rms

            if fss_enabled and point_fss_channels:
                if fss_max_points is None or fss_saved_count < fss_max_points:
                    out_png = _save_fss_plot_for_point(
                        comparison_data_dir=comparison_data_dir,
                        run_tag=tag,
                        r1=float(r1),
                        r2=float(r2),
                        k_values=k_values,
                        point_fss_channels=point_fss_channels,
                        dpi=fss_dpi,
                    )
                    fss_plot_index_rows.append([float(r1), float(r2), out_png, len(point_fss_channels)])
                    fss_saved_count += 1

            if missing_sizes:
                log(
                    f"Point r1={r1:.6f}, r2={r2:.6f}: missing sizes {missing_sizes}; "
                    "continuum fits use available sizes only"
                )

    header_common = [
        "Large-reference scoring workflow",
        f"run_tag={tag}",
        f"reference_data_file={reference_data_path}",
        f"reference_metadata_file={reference_metadata_path}",
        f"tests_manifest={tests_manifest_path}",
        f"sizes={sizes}",
        f"k3={k3}",
        f"k_values={k_values}",
        f"k_denominator={k_denom}",
        f"weighted={weighted}",
        f"power_fit=A + B*(1/L)^C with C in [{c_min}, {c_max}], min_sizes_for_free_C={min_sizes_for_free_C}",
    ]

    write_dat(
        os.path.join(reference_data_dir, "reference_channels_scored.dat"),
        header_common,
        ["cycle", "k", "t", "G_ref", "sigma_ref"],
        ref_rows,
    )
    write_dat(
        os.path.join(test_data_dir, "raw_test_channels_scored.dat"),
        header_common,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_test_rows,
    )
    write_dat(
        os.path.join(comparison_data_dir, "continuum_test_channels.dat"),
        header_common,
        ["r1", "r2", "cycle", "k", "t", "A", "sigma_A", "B", "C", "n_sizes_used", "fit_mode"],
        continuum_rows,
    )
    write_dat(
        os.path.join(comparison_data_dir, "channel_comparison.dat"),
        header_common,
        [
            "r1",
            "r2",
            "cycle",
            "k",
            "t",
            "A_test_cont",
            "sigma_A_test_cont",
            "G_ref",
            "sigma_ref",
            "diff",
            "variance",
            "z",
            "fit_mode",
        ],
        compare_rows,
    )
    write_dat(
        os.path.join(comparison_data_dir, "score_map.dat"),
        header_common,
        ["r1", "r2", "score", "n_channels_scored", "z_rms", "n_channels_with_z", "n_missing_sizes"],
        score_rows,
    )
    write_dat(
        os.path.join(comparison_data_dir, "fss_plot_data.dat"),
        header_common,
        ["r1", "r2", "cycle", "k", "t", "source", "L", "invL", "G", "sigma_G", "fit_mode", "B", "C"],
        fss_rows,
    )
    write_dat(
        os.path.join(comparison_data_dir, "fss_plot_index.dat"),
        header_common,
        ["r1", "r2", "plot_file", "n_channel_panels"],
        fss_plot_index_rows,
    )

    finite_scores = [row for row in score_rows if np.isfinite(row[2])]
    best = None
    if finite_scores:
        best = min(finite_scores, key=lambda row: row[2])
        write_dat(
            os.path.join(comparison_data_dir, "score_minimum.dat"),
            header_common,
            ["r1_min", "r2_min", "score_min", "z_rms_at_min"],
            [[best[0], best[1], best[2], best[4]]],
        )
        log(f"Best point: r1={best[0]:.6f}, r2={best[1]:.6f}, score={best[2]:.8g}")
    else:
        log("No finite scores produced. Check missing payloads or fit failures.")

    _save_heatmaps(r1_values, r2_values, score_grid, zscore_grid, comparison_data_dir, tag)

    manifest = {
        "created_at": timestamp(),
        "config_path": config_path,
        "run_tag": tag,
        "run_root": run_root,
        "reference_data_dir": reference_data_dir,
        "test_data_dir": test_data_dir,
        "comparison_analysis_data_dir": comparison_data_dir,
        "reference_data_file": reference_data_path,
        "reference_metadata_file": reference_metadata_path,
        "tests_manifest": tests_manifest_path,
        "sizes": sizes,
        "k_values": k_values,
        "k_denominator": k_denom,
        "weighted": weighted,
        "power_fit": {
            "model": "A + B*(1/L)^C",
            "c_min": c_min,
            "c_max": c_max,
            "c_initial": c_initial,
            "min_sizes_for_free_C": min_sizes_for_free_C,
        },
        "fss_plots": {
            "enabled": fss_enabled,
            "dpi": fss_dpi,
            "max_points": fss_max_points,
            "saved_points": fss_saved_count,
        },
        "outputs": {
            "reference_channels_scored": os.path.join(reference_data_dir, "reference_channels_scored.dat"),
            "raw_test_channels_scored": os.path.join(test_data_dir, "raw_test_channels_scored.dat"),
            "continuum_test_channels": os.path.join(comparison_data_dir, "continuum_test_channels.dat"),
            "channel_comparison": os.path.join(comparison_data_dir, "channel_comparison.dat"),
            "score_map": os.path.join(comparison_data_dir, "score_map.dat"),
            "score_minimum": os.path.join(comparison_data_dir, "score_minimum.dat"),
            "fss_plot_data": os.path.join(comparison_data_dir, "fss_plot_data.dat"),
            "fss_plot_index": os.path.join(comparison_data_dir, "fss_plot_index.dat"),
            "fss_plot_dir": os.path.join(comparison_data_dir, "fss_plots"),
            "score_heatmap": os.path.join(comparison_data_dir, "score_heatmap.png"),
            "zscore_heatmap": os.path.join(comparison_data_dir, "zscore_heatmap.png"),
            "score_zscore_heatmaps": os.path.join(comparison_data_dir, "score_zscore_heatmaps.png"),
        },
        "summary": {
            "n_points": len(score_rows),
            "n_finite_scores": len(finite_scores),
            "best": None
            if best is None
            else {
                "r1": float(best[0]),
                "r2": float(best[1]),
                "score": float(best[2]),
                "z_rms": float(best[4]),
            },
        },
        "config": cfg,
    }
    manifest_path = os.path.join(comparison_data_dir, "manifest_score.json")
    save_json(manifest_path, manifest)
    log(f"Scoring manifest written: {manifest_path}")


if __name__ == "__main__":
    main()
