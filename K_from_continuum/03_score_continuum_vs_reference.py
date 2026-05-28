#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from workflow_common import (
    check_required_sections,
    evaluate_observable_fit,
    ensure_dir,
    fit_observable_continuum_power,
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

DEFAULT_CONFIG_PATH = "configs/score_example.json"
_SCORE_MODES = {"raw", "additive_offset", "shape_only"}


def _parse_score_mode(value: Any) -> str:
    mode = "raw" if value is None else str(value).strip().lower()
    if mode not in _SCORE_MODES:
        allowed = ", ".join(sorted(_SCORE_MODES))
        raise ValueError(f"analysis.score_mode must be one of {{{allowed}}}")
    return mode


def _weighted_average(values: list[float], variances: list[float]) -> float:
    numerator = 0.0
    denominator = 0.0
    finite_values: list[float] = []

    for value, variance in zip(values, variances):
        if not np.isfinite(value):
            continue
        finite_values.append(float(value))
        weight = 1.0
        if np.isfinite(variance) and variance > 0.0:
            weight = 1.0 / float(variance)
        numerator += weight * float(value)
        denominator += weight

    if denominator > 0.0:
        return float(numerator / denominator)
    if finite_values:
        return float(np.mean(finite_values))
    return 0.0


def _compute_score_nuisance(
    channel_entries: list[dict[str, Any]],
    score_mode: str,
) -> tuple[float, dict[int, float]]:
    global_offset = 0.0
    cycle_offsets = {0: 0.0, 1: 0.0, 2: 0.0}
    if score_mode == "raw":
        return global_offset, cycle_offsets

    if score_mode == "additive_offset":
        global_offset = _weighted_average(
            [float(entry["diff"]) for entry in channel_entries],
            [float(entry["var"]) for entry in channel_entries],
        )
        return global_offset, cycle_offsets

    if score_mode == "shape_only":
        for cycle in range(3):
            cycle_rows = [entry for entry in channel_entries if int(entry["cycle"]) == cycle]
            cycle_offsets[cycle] = _weighted_average(
                [float(entry["diff"]) for entry in cycle_rows],
                [float(entry["var"]) for entry in cycle_rows],
            )
        return global_offset, cycle_offsets

    raise ValueError(f"unsupported score mode: {score_mode}")


def _channel_nuisance(
    *,
    cycle: int,
    score_mode: str,
    global_offset: float,
    cycle_offsets: dict[int, float],
) -> float:
    if score_mode == "raw":
        return 0.0
    if score_mode == "additive_offset":
        return float(global_offset)
    if score_mode == "shape_only":
        return float(cycle_offsets.get(int(cycle), 0.0))
    raise ValueError(f"unsupported score mode: {score_mode}")


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


def _ratio_key(value: float) -> float:
    return round(float(value), 12)


def _iter_dat_rows(path: str):
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            yield line.split()


def _load_reference_continuum_channels(path: str) -> tuple[dict[tuple[int, int], dict[str, Any]], list[list[Any]]]:
    lookup: dict[tuple[int, int], dict[str, Any]] = {}
    rows: list[list[Any]] = []
    for parts in _iter_dat_rows(path):
        cycle = int(parts[0])
        kval = int(parts[1])
        row = {
            "cycle": cycle,
            "k": kval,
            "t": float(parts[2]),
            "A": float(parts[3]),
            "sigma_A": float(parts[4]),
            "B": float(parts[5]),
            "C": float(parts[6]),
            "n_sizes_used": int(parts[7]),
            "fit_mode": str(parts[8]),
        }
        lookup[(cycle, kval)] = row
        rows.append(
            [
                cycle,
                kval,
                row["t"],
                row["A"],
                row["sigma_A"],
                row["B"],
                row["C"],
                row["n_sizes_used"],
                row["fit_mode"],
            ]
        )
    return lookup, rows


def _load_reference_raw_channels(path: str) -> dict[tuple[int, int], list[dict[str, Any]]]:
    lookup: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for parts in _iter_dat_rows(path):
        L = int(parts[0])
        cycle = int(parts[6])
        kval = int(parts[7])
        row = {
            "L": L,
            "invL": 1.0 / float(L) if L > 0 else np.nan,
            "t": float(parts[8]),
            "G": float(parts[9]),
            "sigma_G": float(parts[10]),
        }
        lookup.setdefault((cycle, kval), []).append(row)
    for rows in lookup.values():
        rows.sort(key=lambda row: row["L"])
    return lookup


def _draw_heatmap(ax, grid: np.ndarray, r1_values: list[float], r2_values: list[float], title: str, cbar_label: str, cmap: str) -> None:
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
        ax.plot(r1_values[ii], r2_values[jj], marker="*", markersize=12, markerfacecolor="white", markeredgecolor="k")
    else:
        im = ax.imshow(np.zeros_like(grid), origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=0.0, vmax=1.0)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cbar_label)


def _save_heatmaps(
    r1_values: list[float],
    r2_values: list[float],
    score_grid: np.ndarray,
    zscore_grid: np.ndarray,
    comparison_data_dir: str,
    run_tag: str,
    score_mode: str,
) -> None:
    if len(r1_values) == 0 or len(r2_values) == 0:
        return

    score_path = os.path.join(comparison_data_dir, "score_heatmap.png")
    zscore_path = os.path.join(comparison_data_dir, "zscore_heatmap.png")
    combined_path = os.path.join(comparison_data_dir, "score_zscore_heatmaps.png")

    fig1, ax1 = plt.subplots(figsize=(6.5, 5.3))
    _draw_heatmap(
        ax1,
        score_grid,
        r1_values,
        r2_values,
        f"Score heatmap ({score_mode}, tag={run_tag})",
        "score",
        "viridis",
    )
    fig1.tight_layout()
    fig1.savefig(score_path, dpi=150)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(6.5, 5.3))
    _draw_heatmap(
        ax2,
        zscore_grid,
        r1_values,
        r2_values,
        f"RMS z-score heatmap ({score_mode}, tag={run_tag})",
        "z_rms",
        "magma",
    )
    fig2.tight_layout()
    fig2.savefig(zscore_path, dpi=150)
    plt.close(fig2)

    fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(13.0, 5.3))
    _draw_heatmap(ax3, score_grid, r1_values, r2_values, "Score", "score", "viridis")
    _draw_heatmap(ax4, zscore_grid, r1_values, r2_values, "RMS z-score", "z_rms", "magma")
    fig3.suptitle(f"Score and z-score heatmaps (mode={score_mode}, tag={run_tag})")
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
    fig_h = max(3.4 * n_rows, 9.8)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    for channel in point_fss_channels:
        cycle = int(channel["cycle"])
        ik = int(channel["k_index"])
        ax = axes[cycle][ik]

        test_L = np.asarray(channel["test_L_vals"], dtype=float)
        test_y = np.asarray(channel["test_y_vals"], dtype=float)
        test_s = np.asarray(channel["test_s_vals"], dtype=float)
        if test_L.size > 0:
            order = np.argsort(test_L)
            x = 1.0 / test_L[order]
            y = test_y[order]
            s = test_s[order]
            yerr = np.where(np.isfinite(s), s, 0.0)
            ax.errorbar(x, y, yerr=yerr, fmt="o", color="#1f77b4", label="test sizes")

        ref_L = np.asarray(channel["ref_L_vals"], dtype=float)
        ref_y = np.asarray(channel["ref_y_vals"], dtype=float)
        ref_s = np.asarray(channel["ref_s_vals"], dtype=float)
        if ref_L.size > 0:
            order = np.argsort(ref_L)
            x = 1.0 / ref_L[order]
            y = ref_y[order]
            s = ref_s[order]
            yerr = np.where(np.isfinite(s), s, 0.0)
            ax.errorbar(x, y, yerr=yerr, fmt="s", color="#d62728", label="reference sizes")

        test_A = float(channel["test_A"])
        test_sigma_A = float(channel["test_sigma_A"])
        test_B = float(channel["test_B"])
        test_C = float(channel["test_C"])
        test_fit_mode = str(channel.get("test_fit_mode", "power_fit"))
        ref_A = float(channel["ref_A"])
        ref_sigma_A = float(channel["ref_sigma_A"])
        ref_B = float(channel["ref_B"])
        ref_C = float(channel["ref_C"])
        ref_fit_mode = str(channel.get("ref_fit_mode", "power_fit"))

        if np.isfinite(test_A):
            if np.isfinite(test_sigma_A):
                ax.errorbar([0.0], [test_A], yerr=[test_sigma_A], fmt="*", markersize=12, color="#2ca02c", label="test continuum")
            else:
                ax.plot([0.0], [test_A], "*", markersize=12, color="#2ca02c", label="test continuum")
        if np.isfinite(ref_A):
            if np.isfinite(ref_sigma_A):
                ax.errorbar([0.0], [ref_A], yerr=[ref_sigma_A], fmt="*", markersize=12, color="#ff7f0e", label="reference continuum")
            else:
                ax.plot([0.0], [ref_A], "*", markersize=12, color="#ff7f0e", label="reference continuum")

        x_candidates = [0.0]
        if test_L.size > 0:
            x_candidates.extend((1.0 / test_L).tolist())
        if ref_L.size > 0:
            x_candidates.extend((1.0 / ref_L).tolist())
        x_max = max(v for v in x_candidates if np.isfinite(v))
        x_fit = np.linspace(0.0, max(x_max, 1e-3), 200)

        if np.isfinite(test_A) and np.isfinite(test_B) and np.isfinite(test_C):
            y_fit = evaluate_observable_fit(x_fit, test_A, test_B, test_C, test_fit_mode)
            ax.plot(x_fit, y_fit, "-", color="#2ca02c", alpha=0.9, label="test fit")
        if np.isfinite(ref_A) and np.isfinite(ref_B) and np.isfinite(ref_C):
            y_fit = evaluate_observable_fit(x_fit, ref_A, ref_B, ref_C, ref_fit_mode)
            ax.plot(x_fit, y_fit, "--", color="#ff7f0e", alpha=0.9, label="reference fit")

        ax.set_title(f"cycle={cycle}  k={channel['k']}  t={channel['t']:.3f}")
        if cycle == n_rows - 1:
            ax.set_xlabel("1/L")
        if ik == 0:
            ax.set_ylabel("G")
        ax.grid(alpha=0.25)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        fig.legend(list(unique.values()), list(unique.keys()), loc="upper center", ncol=min(6, len(unique)))

    fig.suptitle(f"Continuum-vs-continuum FSS channels  (tag={run_tag}, r1={r1:.6f}, r2={r2:.6f})")
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
            "Step 3: fit each test channel to the continuum and score it against the "
            "continuum reference channels from Step 1."
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
    _check_required_keys(cfg["input"], "input", ["reference_manifest", "tests_manifest"])
    _check_required_keys(cfg["analysis"], "analysis", ["k_values", "k_denominator", "weighted", "power_fit"])
    _check_required_keys(cfg["analysis"]["power_fit"], "analysis.power_fit", ["c_min", "c_max", "c_initial"])

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _HERE)
    run_root = os.path.join(results_root, tag)
    reference_data_dir = os.path.join(run_root, "reference_data")
    test_data_dir = os.path.join(run_root, "test_data")
    comparison_data_dir = os.path.join(run_root, "comparison_analysis_data")
    ensure_dir(reference_data_dir)
    ensure_dir(test_data_dir)
    ensure_dir(comparison_data_dir)

    reference_manifest_path = resolve_path(str(cfg["input"]["reference_manifest"]), _HERE)
    reference_manifest = load_json(reference_manifest_path)
    reference_continuum_path = resolve_path(
        str(cfg["input"].get("reference_continuum_channels") or reference_manifest["reference_continuum_channels"]),
        _HERE,
    )
    reference_raw_path = resolve_path(
        str(cfg["input"].get("reference_raw_channels") or reference_manifest["reference_raw_channels"]),
        _HERE,
    )
    tests_manifest_path = resolve_path(str(cfg["input"]["tests_manifest"]), _HERE)

    reference_lookup, reference_rows = _load_reference_continuum_channels(reference_continuum_path)
    reference_raw_lookup = _load_reference_raw_channels(reference_raw_path)
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
    score_mode = _parse_score_mode(cfg["analysis"].get("score_mode", "raw"))
    c_min = float(cfg["analysis"]["power_fit"]["c_min"])
    c_max = float(cfg["analysis"]["power_fit"]["c_max"])
    c_initial = float(cfg["analysis"]["power_fit"]["c_initial"])
    fit_method = str(cfg["analysis"]["power_fit"].get("method", "power")).strip().lower()
    min_sizes_for_free_C = int(cfg["analysis"]["power_fit"].get("min_sizes_for_free_C", 8))
    if c_max <= c_min:
        raise ValueError("analysis.power_fit requires c_max > c_min")
    if min_sizes_for_free_C < 3:
        raise ValueError("analysis.power_fit.min_sizes_for_free_C must be >= 3")
    if fit_method not in {"power", "taylor2"}:
        raise ValueError("analysis.power_fit.method must be 'power' or 'taylor2'")

    fss_cfg = cfg["analysis"].get("fss_plots") or {}
    fss_enabled = bool(fss_cfg.get("enabled", True))
    fss_dpi = int(fss_cfg.get("dpi", 150))
    fss_max_points_raw = fss_cfg.get("max_points", None)
    fss_max_points = None if fss_max_points_raw in (None, "") else int(fss_max_points_raw)
    if fss_max_points is not None and fss_max_points < 1:
        raise ValueError("analysis.fss_plots.max_points must be >= 1 when set")

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

    raw_test_rows: list[list[Any]] = []
    continuum_rows: list[list[Any]] = []
    compare_rows: list[list[Any]] = []
    score_rows: list[list[Any]] = []
    score_diagnostic_rows: list[list[Any]] = []
    fss_rows: list[list[Any]] = []
    fss_plot_index_rows: list[list[Any]] = []

    payload_cache: dict[tuple[str, str], dict[str, Any]] = {}
    score_grid = np.full((len(r2_values), len(r1_values)), np.nan)
    zscore_grid = np.full_like(score_grid, np.nan)
    r1_index = {_ratio_key(v): i for i, v in enumerate(r1_values)}
    r2_index = {_ratio_key(v): j for j, v in enumerate(r2_values)}
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
            raw_z2_sum = 0.0
            raw_n_z = 0
            channel_score_entries: list[dict[str, Any]] = []

            for cycle in range(3):
                for ik, kval in enumerate(k_values):
                    ref_key = (cycle, int(kval))
                    if ref_key not in reference_lookup:
                        raise ValueError(f"Reference continuum channels are missing cycle={cycle}, k={kval}")
                    ref_entry = reference_lookup[ref_key]
                    ref_raw_rows = reference_raw_lookup.get(ref_key, [])

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

                    A, sigma_A, B, C, n_used, fit_mode = fit_observable_continuum_power(
                        np.array(L_vals, dtype=float),
                        np.array(y_vals, dtype=float),
                        np.array(s_vals, dtype=float),
                        fit_method=fit_method,
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

                    ref_value = float(ref_entry["A"])
                    ref_sigma = float(ref_entry["sigma_A"])
                    ref_B = float(ref_entry["B"])
                    ref_C = float(ref_entry["C"])
                    ref_fit_mode = str(ref_entry["fit_mode"])

                    for L_i, y_i, s_i in zip(L_vals, y_vals, s_vals):
                        invL_i = (1.0 / float(L_i)) if float(L_i) > 0.0 else np.nan
                        fss_rows.append(
                            [
                                float(r1),
                                float(r2),
                                cycle,
                                int(kval),
                                float(fractions[ik]),
                                "test",
                                float(L_i),
                                invL_i,
                                float(y_i),
                                float(s_i),
                                "",
                                np.nan,
                                np.nan,
                            ]
                        )
                    for ref_row in ref_raw_rows:
                        fss_rows.append(
                            [
                                float(r1),
                                float(r2),
                                cycle,
                                int(kval),
                                float(fractions[ik]),
                                "reference_size",
                                float(ref_row["L"]),
                                float(ref_row["invL"]),
                                float(ref_row["G"]),
                                float(ref_row["sigma_G"]),
                                "",
                                np.nan,
                                np.nan,
                            ]
                        )
                    fss_rows.append(
                        [
                            float(r1),
                            float(r2),
                            cycle,
                            int(kval),
                            float(fractions[ik]),
                            "test_continuum_fit",
                            0.0,
                            0.0,
                            float(A),
                            float(sigma_A),
                            fit_mode,
                            float(B),
                            float(C),
                        ]
                    )
                    fss_rows.append(
                        [
                            float(r1),
                            float(r2),
                            cycle,
                            int(kval),
                            float(fractions[ik]),
                            "reference_continuum_fit",
                            0.0,
                            0.0,
                            ref_value,
                            ref_sigma,
                            ref_fit_mode,
                            ref_B,
                            ref_C,
                        ]
                    )

                    point_fss_channels.append(
                        {
                            "cycle": cycle,
                            "k_index": ik,
                            "k": int(kval),
                            "t": float(fractions[ik]),
                            "test_L_vals": [float(v) for v in L_vals],
                            "test_y_vals": [float(v) for v in y_vals],
                            "test_s_vals": [float(v) for v in s_vals],
                            "test_A": float(A),
                            "test_sigma_A": float(sigma_A),
                            "test_B": float(B),
                            "test_C": float(C),
                            "test_fit_mode": fit_mode,
                            "ref_L_vals": [float(row["L"]) for row in ref_raw_rows],
                            "ref_y_vals": [float(row["G"]) for row in ref_raw_rows],
                            "ref_s_vals": [float(row["sigma_G"]) for row in ref_raw_rows],
                            "ref_A": ref_value,
                            "ref_sigma_A": ref_sigma,
                            "ref_B": ref_B,
                            "ref_C": ref_C,
                            "ref_fit_mode": ref_fit_mode,
                        }
                    )

                    diff = np.nan
                    var = np.nan
                    z = np.nan
                    if np.isfinite(A) and np.isfinite(ref_value):
                        diff = A - ref_value
                        if np.isfinite(sigma_A) and np.isfinite(ref_sigma):
                            var = sigma_A * sigma_A + ref_sigma * ref_sigma
                        if np.isfinite(var) and var > 0.0:
                            z = diff / np.sqrt(var)
                            raw_z2_sum += z * z
                            raw_n_z += 1

                    channel_score_entries.append(
                        {
                            "r1": float(r1),
                            "r2": float(r2),
                            "cycle": int(cycle),
                            "k": int(kval),
                            "t": float(fractions[ik]),
                            "A": A,
                            "sigma_A": sigma_A,
                            "ref_A": ref_value,
                            "ref_sigma_A": ref_sigma,
                            "diff": diff,
                            "var": var,
                            "z": z,
                            "fit_mode": fit_mode,
                            "ref_fit_mode": ref_fit_mode,
                        }
                    )

            global_offset, cycle_offsets = _compute_score_nuisance(channel_score_entries, score_mode)
            raw_diff_values: list[float] = []
            residual_values: list[float] = []

            for entry in channel_score_entries:
                diff = float(entry["diff"])
                var = float(entry["var"])
                nuisance = _channel_nuisance(
                    cycle=int(entry["cycle"]),
                    score_mode=score_mode,
                    global_offset=global_offset,
                    cycle_offsets=cycle_offsets,
                )
                score_residual = np.nan
                score_z = np.nan
                if np.isfinite(diff):
                    raw_diff_values.append(diff)
                    score_residual = diff - nuisance
                    residual_values.append(float(score_residual))
                    weight = 1.0
                    if weighted and np.isfinite(var) and var > 0.0:
                        weight = 1.0 / var
                    diff2 = score_residual * score_residual
                    if np.isfinite(var) and var > 0.0:
                        diff2 = min(diff2, 25.0 * var)
                    score += weight * diff2
                    n_score += 1
                    if np.isfinite(var) and var > 0.0:
                        score_z = score_residual / np.sqrt(var)
                        z2_sum += score_z * score_z
                        n_z += 1

                compare_rows.append(
                    [
                        float(entry["r1"]),
                        float(entry["r2"]),
                        int(entry["cycle"]),
                        int(entry["k"]),
                        float(entry["t"]),
                        entry["A"],
                        entry["sigma_A"],
                        entry["ref_A"],
                        entry["ref_sigma_A"],
                        entry["diff"],
                        entry["var"],
                        entry["z"],
                        entry["fit_mode"],
                        entry["ref_fit_mode"],
                        score_residual,
                        score_z,
                        nuisance,
                    ]
                )

            score_value = float(score) if n_score > 0 else np.nan
            z_rms = float(np.sqrt(z2_sum / n_z)) if n_z > 0 else np.nan
            raw_z_rms = float(np.sqrt(raw_z2_sum / raw_n_z)) if raw_n_z > 0 else np.nan
            mean_raw_diff = float(np.mean(raw_diff_values)) if raw_diff_values else np.nan
            mean_score_residual = float(np.mean(residual_values)) if residual_values else np.nan
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
            score_diagnostic_rows.append(
                [
                    float(r1),
                    float(r2),
                    score_mode,
                    float(global_offset),
                    float(cycle_offsets[0]),
                    float(cycle_offsets[1]),
                    float(cycle_offsets[2]),
                    mean_raw_diff,
                    mean_score_residual,
                    raw_z_rms,
                    z_rms,
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

    fit_model_desc = (
        "A + B*(1/L) + C*(1/L)^2"
        if fit_method == "taylor2"
        else f"A + B*(1/L)^C with C in [{c_min}, {c_max}], min_sizes_for_free_C={min_sizes_for_free_C}"
    )

    header_common = [
        "Continuum-to-continuum scoring workflow",
        f"run_tag={tag}",
        f"reference_manifest={reference_manifest_path}",
        f"reference_continuum_channels={reference_continuum_path}",
        f"reference_raw_channels={reference_raw_path}",
        f"tests_manifest={tests_manifest_path}",
        f"sizes={sizes}",
        f"k3={k3}",
        f"k_values={k_values}",
        f"k_denominator={k_denom}",
        f"weighted={weighted}",
        f"score_mode={score_mode}",
        f"continuum_fit={fit_model_desc}",
    ]

    reference_channels_scored_path = os.path.join(reference_data_dir, "reference_channels_scored.dat")
    raw_test_channels_path = os.path.join(test_data_dir, "raw_test_channels_scored.dat")
    continuum_test_channels_path = os.path.join(comparison_data_dir, "continuum_test_channels.dat")
    channel_comparison_path = os.path.join(comparison_data_dir, "channel_comparison.dat")
    score_map_path = os.path.join(comparison_data_dir, "score_map.dat")
    score_diagnostics_path = os.path.join(comparison_data_dir, "score_diagnostics.dat")
    fss_plot_data_path = os.path.join(comparison_data_dir, "fss_plot_data.dat")
    fss_plot_index_path = os.path.join(comparison_data_dir, "fss_plot_index.dat")
    manifest_score_path = os.path.join(comparison_data_dir, "manifest_score.json")

    write_dat(
        reference_channels_scored_path,
        header_common,
        ["cycle", "k", "t", "A_ref", "sigma_A_ref", "B_ref", "C_ref", "n_sizes_used", "fit_mode"],
        reference_rows,
    )
    write_dat(
        raw_test_channels_path,
        header_common,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_test_rows,
    )
    write_dat(
        continuum_test_channels_path,
        header_common,
        ["r1", "r2", "cycle", "k", "t", "A", "sigma_A", "B", "C", "n_sizes_used", "fit_mode"],
        continuum_rows,
    )
    write_dat(
        channel_comparison_path,
        header_common,
        [
            "r1",
            "r2",
            "cycle",
            "k",
            "t",
            "A_test_cont",
            "sigma_A_test_cont",
            "A_ref_cont",
            "sigma_A_ref_cont",
            "diff",
            "variance",
            "z",
            "test_fit_mode",
            "reference_fit_mode",
            "score_residual",
            "score_z",
            "score_nuisance",
        ],
        compare_rows,
    )
    write_dat(
        score_map_path,
        header_common,
        ["r1", "r2", "score", "n_channels_scored", "z_rms", "n_channels_with_z", "n_missing_sizes"],
        score_rows,
    )
    write_dat(
        score_diagnostics_path,
        header_common,
        [
            "r1",
            "r2",
            "score_mode",
            "global_offset",
            "cycle0_offset",
            "cycle1_offset",
            "cycle2_offset",
            "mean_raw_diff",
            "mean_score_residual",
            "raw_z_rms",
            "score_z_rms",
        ],
        score_diagnostic_rows,
    )
    write_dat(
        fss_plot_data_path,
        header_common,
        ["r1", "r2", "cycle", "k", "t", "source", "L", "invL", "G", "sigma_G", "fit_mode", "B", "C"],
        fss_rows,
    )
    write_dat(
        fss_plot_index_path,
        header_common,
        ["r1", "r2", "plot_file", "n_channel_panels"],
        fss_plot_index_rows,
    )

    finite_scores = [row for row in score_rows if np.isfinite(row[2])]
    best = None
    score_minimum_path = os.path.join(comparison_data_dir, "score_minimum.dat")
    if finite_scores:
        best = min(finite_scores, key=lambda row: row[2])
        write_dat(
            score_minimum_path,
            header_common,
            ["r1_min", "r2_min", "score_min", "z_rms_at_min"],
            [[best[0], best[1], best[2], best[4]]],
        )
        log(f"Best point ({score_mode}): r1={best[0]:.6f}, r2={best[1]:.6f}, score={best[2]:.8g}")
    else:
        log("No finite scores produced. Check missing payloads or fit failures.")

    _save_heatmaps(r1_values, r2_values, score_grid, zscore_grid, comparison_data_dir, tag, score_mode)

    manifest = {
        "created_at": timestamp(),
        "config_path": config_path,
        "run_tag": tag,
        "run_root": run_root,
        "reference_data_dir": reference_data_dir,
        "test_data_dir": test_data_dir,
        "comparison_analysis_data_dir": comparison_data_dir,
        "reference_manifest": reference_manifest_path,
        "reference_continuum_channels": reference_continuum_path,
        "reference_raw_channels": reference_raw_path,
        "tests_manifest": tests_manifest_path,
        "sizes": sizes,
        "k_values": k_values,
        "k_denominator": k_denom,
        "weighted": weighted,
        "score_mode": score_mode,
        "power_fit": {
            "method": fit_method,
            "model": fit_model_desc,
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
            "reference_channels_scored": reference_channels_scored_path,
            "raw_test_channels_scored": raw_test_channels_path,
            "continuum_test_channels": continuum_test_channels_path,
            "channel_comparison": channel_comparison_path,
            "score_map": score_map_path,
            "score_diagnostics": score_diagnostics_path,
            "score_minimum": score_minimum_path,
            "fss_plot_data": fss_plot_data_path,
            "fss_plot_index": fss_plot_index_path,
            "fss_plot_dir": os.path.join(comparison_data_dir, "fss_plots"),
            "score_heatmap": os.path.join(comparison_data_dir, "score_heatmap.png"),
            "zscore_heatmap": os.path.join(comparison_data_dir, "zscore_heatmap.png"),
            "score_zscore_heatmaps": os.path.join(comparison_data_dir, "score_zscore_heatmaps.png"),
        },
        "summary": {
            "n_points": len(score_rows),
            "n_finite_scores": len(finite_scores),
            "best": None if best is None else {
                "r1": float(best[0]),
                "r2": float(best[1]),
                "score": float(best[2]),
                "z_rms": float(best[4]),
            },
        },
        "config": cfg,
    }
    save_json(manifest_score_path, manifest)
    log(f"Scoring manifest written: {manifest_score_path}")


if __name__ == "__main__":
    main()