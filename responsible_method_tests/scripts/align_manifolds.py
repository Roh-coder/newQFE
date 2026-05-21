#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from plot_fss_interpolated_manifolds import _load_dat_rows, _load_json, _tile_periodic_points, _unit_cell_outline


COMPARISON_MODES = ("symmetric", "twisted_reference")


def _to_float_array(values: Any) -> np.ndarray:
    return np.asarray(np.ma.asarray(values).filled(np.nan), dtype=float)


def _lookup_numeric(rows: list[dict[str, Any]], key: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for row in rows:
        if int(row["is_origin"]) != 0:
            continue
        value = float(row[key])
        if np.isfinite(value):
            out[int(row["point_id"])] = value
    return out


def _build_periodic_interpolator(
    point_rows: list[dict[str, Any]],
    value_lookup: dict[int, float],
) -> mtri.LinearTriInterpolator:
    a_tiled, b_tiled, z_tiled = _tile_periodic_points(point_rows, value_lookup)
    triangulation = mtri.Triangulation(a_tiled, b_tiled)
    return mtri.LinearTriInterpolator(triangulation, z_tiled)


def _build_regular_ab_grid(grid_size: int) -> tuple[np.ndarray, np.ndarray]:
    a_values = np.linspace(0.0, 1.0, grid_size)
    b_values = np.linspace(0.0, 1.0, grid_size)
    return np.meshgrid(a_values, b_values)


def _ab_to_nu(A: np.ndarray, B: np.ndarray, *, tau_real: float, tau_imag: float) -> tuple[np.ndarray, np.ndarray]:
    return A + tau_real * B, tau_imag * B


def _summarize_pointwise_residuals(diff: np.ndarray, sigma: np.ndarray) -> dict[str, Any]:
    finite = np.isfinite(diff) & np.isfinite(sigma) & (sigma > 0.0)
    if not np.any(finite):
        return {
            "n_points": 0,
            "rms_abs": float("nan"),
            "mean_abs": float("nan"),
            "max_abs": float("nan"),
            "rms_z": float("nan"),
            "mean_abs_z": float("nan"),
            "max_abs_z": float("nan"),
            "chi2_per_point": float("nan"),
        }

    abs_diff = np.abs(diff[finite])
    z_values = diff[finite] / sigma[finite]
    abs_z = np.abs(z_values)
    return {
        "n_points": int(np.count_nonzero(finite)),
        "rms_abs": float(np.sqrt(np.mean(diff[finite] ** 2))),
        "mean_abs": float(np.mean(abs_diff)),
        "max_abs": float(np.max(abs_diff)),
        "rms_z": float(np.sqrt(np.mean(z_values ** 2))),
        "mean_abs_z": float(np.mean(abs_z)),
        "max_abs_z": float(np.max(abs_z)),
        "chi2_per_point": float(np.mean(z_values ** 2)),
    }


def _summarize_common_grid(twisted_z: np.ndarray, untwisted_z: np.ndarray) -> dict[str, Any]:
    finite = np.isfinite(twisted_z) & np.isfinite(untwisted_z)
    if not np.any(finite):
        return {
            "n_grid_points": 0,
            "rms_abs": float("nan"),
            "mean_abs": float("nan"),
            "max_abs": float("nan"),
            "relative_rms": float("nan"),
            "correlation": float("nan"),
        }

    diff = twisted_z[finite] - untwisted_z[finite]
    mean_scale = 0.5 * (np.abs(twisted_z[finite]) + np.abs(untwisted_z[finite]))
    positive_scale = mean_scale > 0.0
    if np.count_nonzero(finite) > 1:
        correlation = float(np.corrcoef(twisted_z[finite], untwisted_z[finite])[0, 1])
    else:
        correlation = float("nan")
    return {
        "n_grid_points": int(np.count_nonzero(finite)),
        "rms_abs": float(np.sqrt(np.mean(diff ** 2))),
        "mean_abs": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "relative_rms": float(
            np.sqrt(np.mean(diff ** 2)) / np.mean(mean_scale[positive_scale]) if np.any(positive_scale) else float("nan")
        ),
        "correlation": correlation,
    }


def _interpret_distinguishability(
    twisted_on_untwisted: dict[str, Any],
    untwisted_on_twisted: dict[str, Any],
    *,
    comparison_mode: str,
) -> str:
    if comparison_mode == "twisted_reference":
        threshold = float(twisted_on_untwisted["rms_z"])
        if not np.isfinite(threshold):
            return "undetermined"
        if threshold <= 1.0:
            return "not distinguished within continuum point errors"
        if threshold <= 2.0:
            return "marginally separated"
        return "distinguishable at continuum-point level"

    rms_values = [
        float(twisted_on_untwisted["rms_z"]),
        float(untwisted_on_twisted["rms_z"]),
    ]
    finite = [value for value in rms_values if np.isfinite(value)]
    if not finite:
        return "undetermined"
    threshold = max(finite)
    if threshold <= 1.0:
        return "not distinguished within continuum point errors"
    if threshold <= 2.0:
        return "marginally separated"
    return "distinguishable at continuum-point level"


def load_method_data(method_manifest_path: str) -> dict[str, Any]:
    manifest = _load_json(method_manifest_path)
    continuum_rows = [row for row in _load_dat_rows(str(manifest["pointwise_continuum"])) if int(row["is_origin"]) == 0]

    continuum_lookup = _lookup_numeric(continuum_rows, "A")
    modular_rows = _load_dat_rows(str(manifest["modular_aligned"]))
    modular_lookup = _lookup_numeric(modular_rows, "modular_aligned")
    sigma_lookup = _lookup_numeric(continuum_rows, "sigma_A")

    point_ids = np.asarray([int(row["point_id"]) for row in continuum_rows], dtype=int)
    a_wrap = np.asarray([float(row["a_wrap"]) for row in continuum_rows], dtype=float)
    b_wrap = np.asarray([float(row["b_wrap"]) for row in continuum_rows], dtype=float)
    nu_real = np.asarray([float(row["nu_real"]) for row in continuum_rows], dtype=float)
    nu_imag = np.asarray([float(row["nu_imag"]) for row in continuum_rows], dtype=float)
    continuum_values = np.asarray([float(continuum_lookup[int(point_id)]) for point_id in point_ids], dtype=float)
    sigma_values = np.asarray([float(sigma_lookup.get(int(point_id), float("nan"))) for point_id in point_ids], dtype=float)

    return {
        "method_manifest": os.path.abspath(method_manifest_path),
        "manifest": manifest,
        "benchmark_id": str(manifest["benchmark_id"]),
        "method": str(manifest["method"]),
        "tau_real": float(manifest["target_tau"]["real"]),
        "tau_imag": float(manifest["target_tau"]["imag"]),
        "alignment": _load_json(str(manifest["modular_alignment"])),
        "rows": continuum_rows,
        "point_ids": point_ids,
        "a_wrap": a_wrap,
        "b_wrap": b_wrap,
        "nu_real": nu_real,
        "nu_imag": nu_imag,
        "A": continuum_values,
        "sigma_A": sigma_values,
        "continuum_interpolator": _build_periodic_interpolator(continuum_rows, continuum_lookup),
        "modular_interpolator": _build_periodic_interpolator(continuum_rows, modular_lookup),
    }


def build_method_pair_comparison(
    *,
    twisted_method_manifest_path: str,
    untwisted_method_manifest_path: str,
    grid_size: int,
    comparison_mode: str = "symmetric",
    benchmark_id: str | None = None,
    description: str = "",
    benchmark_manifest_path: str | None = None,
    run_tag: str = "",
) -> tuple[dict[str, Any], dict[str, Any]]:
    if comparison_mode not in COMPARISON_MODES:
        raise ValueError(f"comparison_mode must be one of {COMPARISON_MODES}")

    twisted = load_method_data(twisted_method_manifest_path)
    untwisted = load_method_data(untwisted_method_manifest_path)

    tau_real = float(twisted["tau_real"])
    tau_imag = float(twisted["tau_imag"])
    if not np.isclose(tau_real, float(untwisted["tau_real"])) or not np.isclose(tau_imag, float(untwisted["tau_imag"])):
        raise ValueError("twisted and untwisted manifests do not share the same target tau")

    A_grid, B_grid = _build_regular_ab_grid(grid_size)
    X_grid, Y_grid = _ab_to_nu(A_grid, B_grid, tau_real=tau_real, tau_imag=tau_imag)
    twisted_grid = _to_float_array(twisted["continuum_interpolator"](A_grid, B_grid))
    untwisted_grid = _to_float_array(untwisted["continuum_interpolator"](A_grid, B_grid))

    twisted_on_untwisted = _to_float_array(twisted["continuum_interpolator"](untwisted["a_wrap"], untwisted["b_wrap"]))
    untwisted_on_twisted = _to_float_array(untwisted["continuum_interpolator"](twisted["a_wrap"], twisted["b_wrap"]))
    z_on_untwisted = (twisted_on_untwisted - untwisted["A"]) / untwisted["sigma_A"]
    z_on_twisted = (untwisted_on_twisted - twisted["A"]) / twisted["sigma_A"]

    comparison_id = benchmark_id or str(twisted["benchmark_id"])
    twisted_method_metrics = {
        "method_manifest": str(twisted["method_manifest"]),
        "alpha": float(twisted["alignment"]["alpha"]),
        "chi2_per_dof": float(twisted["alignment"]["chi2_per_dof"]),
        "n_fit_points": int(twisted["alignment"]["n_fit_points"]),
        "couplings": dict(twisted["manifest"].get("couplings", {})),
    }
    untwisted_method_metrics = {
        "method_manifest": str(untwisted["method_manifest"]),
        "alpha": float(untwisted["alignment"]["alpha"]),
        "chi2_per_dof": float(untwisted["alignment"]["chi2_per_dof"]),
        "n_fit_points": int(untwisted["alignment"]["n_fit_points"]),
        "couplings": dict(untwisted["manifest"].get("couplings", {})),
    }
    common_grid_metrics = _summarize_common_grid(twisted_grid, untwisted_grid)
    twisted_on_untwisted_metrics = _summarize_pointwise_residuals(twisted_on_untwisted - untwisted["A"], untwisted["sigma_A"])
    untwisted_on_twisted_metrics = _summarize_pointwise_residuals(untwisted_on_twisted - twisted["A"], twisted["sigma_A"])
    if comparison_mode == "twisted_reference":
        interpretation_basis = "twisted interpolant sampled on untwisted continuum points"
        primary_metric_key = "twisted_on_untwisted"
        primary_rms_z = float(twisted_on_untwisted_metrics["rms_z"])
    else:
        interpretation_basis = "max bidirectional pointwise RMS z across both sampling directions"
        primary_metric_key = "bidirectional_max_rms_z"
        primary_rms_z = max(
            float(twisted_on_untwisted_metrics["rms_z"]),
            float(untwisted_on_twisted_metrics["rms_z"]),
        )

    metrics = {
        "benchmark_id": comparison_id,
        "description": str(description),
        "benchmark_manifest": os.path.abspath(benchmark_manifest_path) if benchmark_manifest_path else None,
        "run_tag": str(run_tag),
        "comparison_mode": comparison_mode,
        "interpretation_basis": interpretation_basis,
        "primary_metric": {
            "key": primary_metric_key,
            "rms_z": primary_rms_z,
        },
        "target_tau": {
            "real": tau_real,
            "imag": tau_imag,
        },
        "grid_size": int(grid_size),
        "twisted": twisted_method_metrics,
        "untwisted": untwisted_method_metrics,
        "common_grid": common_grid_metrics,
        "twisted_on_untwisted": twisted_on_untwisted_metrics,
        "untwisted_on_twisted": untwisted_on_twisted_metrics,
        "interpretation": _interpret_distinguishability(
            twisted_on_untwisted_metrics,
            untwisted_on_twisted_metrics,
            comparison_mode=comparison_mode,
        ),
    }

    plot_data = {
        "tau_real": tau_real,
        "tau_imag": tau_imag,
        "X_grid": X_grid,
        "Y_grid": Y_grid,
        "twisted_grid": twisted_grid,
        "untwisted_grid": untwisted_grid,
        "twisted_sample_x": twisted["nu_real"],
        "twisted_sample_y": twisted["nu_imag"],
        "untwisted_sample_x": untwisted["nu_real"],
        "untwisted_sample_y": untwisted["nu_imag"],
        "z_on_twisted": z_on_twisted,
        "z_on_untwisted": z_on_untwisted,
    }
    return metrics, plot_data


def build_benchmark_comparison(
    benchmark_manifest_path: str,
    *,
    grid_size: int,
    comparison_mode: str = "symmetric",
) -> tuple[dict[str, Any], dict[str, Any]]:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    method_paths = benchmark_manifest.get("methods", {})
    missing = [name for name in ("twisted", "untwisted") if name not in method_paths]
    if missing:
        raise ValueError(f"benchmark manifest is missing methods: {missing}")

    return build_method_pair_comparison(
        twisted_method_manifest_path=str(method_paths["twisted"]),
        untwisted_method_manifest_path=str(method_paths["untwisted"]),
        grid_size=grid_size,
        comparison_mode=comparison_mode,
        benchmark_id=str(benchmark_manifest["benchmark_id"]),
        description=str(benchmark_manifest.get("description", "")),
        benchmark_manifest_path=benchmark_manifest_path,
        run_tag=str(benchmark_manifest.get("run_tag", "")),
    )


def _plot_benchmark_comparison(
    metrics: dict[str, Any],
    plot_data: dict[str, Any],
    *,
    output_path: str,
) -> str:
    benchmark_id = str(metrics["benchmark_id"])
    comparison_mode = str(metrics.get("comparison_mode", "symmetric"))
    X_grid = np.asarray(plot_data["X_grid"], dtype=float)
    Y_grid = np.asarray(plot_data["Y_grid"], dtype=float)
    twisted_grid = np.asarray(plot_data["twisted_grid"], dtype=float)
    untwisted_grid = np.asarray(plot_data["untwisted_grid"], dtype=float)
    diff_grid = twisted_grid - untwisted_grid

    finite_common = np.isfinite(twisted_grid) & np.isfinite(untwisted_grid)
    if not np.any(finite_common):
        raise ValueError("no finite common-grid comparison points are available")

    z_min = min(float(np.nanmin(twisted_grid[finite_common])), float(np.nanmin(untwisted_grid[finite_common])))
    z_max = max(float(np.nanmax(twisted_grid[finite_common])), float(np.nanmax(untwisted_grid[finite_common])))
    if np.isclose(z_min, z_max):
        z_min -= 1.0e-6
        z_max += 1.0e-6
    levels = np.linspace(z_min, z_max, 12)

    diff_limit = float(np.nanmax(np.abs(diff_grid[finite_common])))
    if not np.isfinite(diff_limit) or diff_limit <= 0.0:
        diff_limit = 1.0e-6

    z_residual_limit = 0.0
    for key in ("z_on_twisted", "z_on_untwisted"):
        values = np.asarray(plot_data[key], dtype=float)
        finite = np.isfinite(values)
        if np.any(finite):
            z_residual_limit = max(z_residual_limit, float(np.nanmax(np.abs(values[finite]))))
    z_residual_limit = max(z_residual_limit, 3.0)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    if comparison_mode == "twisted_reference":
        fig.suptitle(
            f"{benchmark_id}: twisted reference continuum vs untwisted continuum points",
            fontsize=15,
            y=0.98,
        )
    else:
        fig.suptitle(
            f"{benchmark_id}: twisted vs untwisted continuum manifolds on a shared target cell",
            fontsize=15,
            y=0.98,
        )
    outline_x, outline_y = _unit_cell_outline(float(plot_data["tau_real"]), float(plot_data["tau_imag"]))

    ax_overlay = axes[0, 0]
    ax_diff = axes[0, 1]
    ax_untwisted = axes[1, 0]
    ax_twisted = axes[1, 1]

    ax_overlay.contour(X_grid, Y_grid, twisted_grid, levels=levels, colors=["tab:blue"], linewidths=1.2)
    ax_overlay.contour(X_grid, Y_grid, untwisted_grid, levels=levels, colors=["tab:orange"], linewidths=1.2, linestyles="--")
    ax_overlay.scatter(plot_data["twisted_sample_x"], plot_data["twisted_sample_y"], s=8, color="tab:blue", alpha=0.25)
    ax_overlay.scatter(plot_data["untwisted_sample_x"], plot_data["untwisted_sample_y"], s=8, color="tab:orange", alpha=0.25)
    ax_overlay.plot(outline_x, outline_y, color="black", linewidth=1.2)
    if comparison_mode == "twisted_reference":
        overlay_title = (
            "Contour overlay (diagnostic context)\n"
            f"corr={metrics['common_grid']['correlation']:.4f}, rel RMS={metrics['common_grid']['relative_rms']:.4f}"
        )
    else:
        overlay_title = (
            "Contour overlay\n"
            f"corr={metrics['common_grid']['correlation']:.4f}, rel RMS={metrics['common_grid']['relative_rms']:.4f}"
        )
    ax_overlay.set_title(overlay_title, fontsize=11)
    ax_overlay.set_xlabel("nu_real")
    ax_overlay.set_ylabel("nu_imag")
    ax_overlay.set_aspect("equal", adjustable="box")
    ax_overlay.grid(True, alpha=0.25)
    ax_overlay.legend(["twisted", "untwisted"], loc="best")

    diff_handle = ax_diff.pcolormesh(
        X_grid,
        Y_grid,
        diff_grid,
        shading="auto",
        cmap="coolwarm",
        vmin=-diff_limit,
        vmax=diff_limit,
    )
    ax_diff.plot(outline_x, outline_y, color="black", linewidth=1.2)
    if comparison_mode == "twisted_reference":
        diff_title = (
            "Common-grid difference (diagnostic context): twisted - untwisted\n"
            f"RMS={metrics['common_grid']['rms_abs']:.5f}, max={metrics['common_grid']['max_abs']:.5f}"
        )
    else:
        diff_title = (
            "Common-grid difference: twisted - untwisted\n"
            f"RMS={metrics['common_grid']['rms_abs']:.5f}, max={metrics['common_grid']['max_abs']:.5f}"
        )
    ax_diff.set_title(diff_title, fontsize=11)
    ax_diff.set_xlabel("nu_real")
    ax_diff.set_ylabel("nu_imag")
    ax_diff.set_aspect("equal", adjustable="box")
    ax_diff.grid(True, alpha=0.25)
    fig.colorbar(diff_handle, ax=ax_diff, fraction=0.046, pad=0.04, label="twisted - untwisted")

    untwisted_scatter = ax_untwisted.scatter(
        plot_data["untwisted_sample_x"],
        plot_data["untwisted_sample_y"],
        c=plot_data["z_on_untwisted"],
        s=18,
        cmap="coolwarm",
        vmin=-z_residual_limit,
        vmax=z_residual_limit,
    )
    ax_untwisted.plot(outline_x, outline_y, color="black", linewidth=1.2)
    if comparison_mode == "twisted_reference":
        untwisted_title = (
            "Primary residual: twisted interpolant on untwisted continuum points\n"
            f"chi2/pt={metrics['twisted_on_untwisted']['chi2_per_point']:.3f}, RMS z={metrics['twisted_on_untwisted']['rms_z']:.3f}"
        )
    else:
        untwisted_title = (
            "Twisted interpolant on untwisted continuum points\n"
            f"chi2/pt={metrics['twisted_on_untwisted']['chi2_per_point']:.3f}, RMS z={metrics['twisted_on_untwisted']['rms_z']:.3f}"
        )
    ax_untwisted.set_title(untwisted_title, fontsize=11)
    ax_untwisted.set_xlabel("nu_real")
    ax_untwisted.set_ylabel("nu_imag")
    ax_untwisted.set_aspect("equal", adjustable="box")
    ax_untwisted.grid(True, alpha=0.25)

    twisted_scatter = ax_twisted.scatter(
        plot_data["twisted_sample_x"],
        plot_data["twisted_sample_y"],
        c=plot_data["z_on_twisted"],
        s=18,
        cmap="coolwarm",
        vmin=-z_residual_limit,
        vmax=z_residual_limit,
    )
    ax_twisted.plot(outline_x, outline_y, color="black", linewidth=1.2)
    if comparison_mode == "twisted_reference":
        twisted_title = (
            "Reverse residual (diagnostic only): untwisted interpolant on twisted points\n"
            f"chi2/pt={metrics['untwisted_on_twisted']['chi2_per_point']:.3f}, RMS z={metrics['untwisted_on_twisted']['rms_z']:.3f}"
        )
    else:
        twisted_title = (
            "Untwisted interpolant on twisted continuum points\n"
            f"chi2/pt={metrics['untwisted_on_twisted']['chi2_per_point']:.3f}, RMS z={metrics['untwisted_on_twisted']['rms_z']:.3f}"
        )
    ax_twisted.set_title(twisted_title, fontsize=11)
    ax_twisted.set_xlabel("nu_real")
    ax_twisted.set_ylabel("nu_imag")
    ax_twisted.set_aspect("equal", adjustable="box")
    ax_twisted.grid(True, alpha=0.25)

    fig.colorbar(twisted_scatter, ax=[ax_untwisted, ax_twisted], fraction=0.046, pad=0.04, label="cross-method z")
    fig.text(
        0.5,
        0.02,
        f"Interpretation basis: {metrics['interpretation_basis']}\nInterpretation: {metrics['interpretation']}",
        ha="center",
        va="center",
        fontsize=11,
    )
    fig.tight_layout(rect=[0.0, 0.05, 1.0, 0.95])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _default_output_paths(benchmark_manifest_path: str, benchmark_id: str) -> tuple[str, str]:
    base_dir = os.path.dirname(os.path.abspath(benchmark_manifest_path))
    figure_path = os.path.join(base_dir, f"{benchmark_id}_twisted_vs_untwisted_common_grid.png")
    metrics_path = os.path.join(base_dir, f"{benchmark_id}_twisted_vs_untwisted_common_grid.json")
    return figure_path, metrics_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare twisted and untwisted continuum manifolds on a shared target torus cell."
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument("--metrics-output", default=None, help="Optional output JSON path")
    parser.add_argument("--grid-size", type=int, default=180, help="Interpolation grid size in each cell direction")
    parser.add_argument(
        "--comparison-mode",
        choices=COMPARISON_MODES,
        default="symmetric",
        help="How to turn pointwise residuals into the distinguishability verdict",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    benchmark_manifest_path = os.path.abspath(args.benchmark_manifest)
    metrics, plot_data = build_benchmark_comparison(
        benchmark_manifest_path,
        grid_size=max(int(args.grid_size), 20),
        comparison_mode=str(args.comparison_mode),
    )
    default_figure_path, default_metrics_path = _default_output_paths(benchmark_manifest_path, str(metrics["benchmark_id"]))
    output_path = os.path.abspath(args.output) if args.output else default_figure_path
    metrics_output_path = os.path.abspath(args.metrics_output) if args.metrics_output else default_metrics_path

    _plot_benchmark_comparison(metrics, plot_data, output_path=output_path)
    with open(metrics_output_path, "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
        handle.write("\n")

    print(f"wrote {output_path}")
    print(f"wrote {metrics_output_path}")
    print(f"comparison_mode={metrics['comparison_mode']}")
    print(f"interpretation={metrics['interpretation']}")
    print(f"common_grid_rel_rms={metrics['common_grid']['relative_rms']:.6f}")
    print(f"twisted_on_untwisted_rms_z={metrics['twisted_on_untwisted']['rms_z']:.6f}")
    print(f"untwisted_on_twisted_rms_z={metrics['untwisted_on_twisted']['rms_z']:.6f}")


if __name__ == "__main__":
    main()