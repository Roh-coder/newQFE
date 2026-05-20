#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_modular_convergence import _load_dat_rows, _load_json


def _scale_labels(method_manifest_path: str) -> dict[int, str]:
    manifest = _load_json(method_manifest_path)
    labels: dict[int, str] = {}
    for payload in manifest["payloads"]:
        scale = int(payload["scale"])
        Lx, Ly, Tx, Ty = [int(v) for v in payload["lattice"]]
        if Tx == 0 and Ty == 0:
            labels[scale] = f"{Lx}x{Ly}"
        else:
            labels[scale] = f"{Lx}x{Ly} ({Tx},{Ty})"
    return labels


def _fit_metrics(values: np.ndarray, sigma: np.ndarray, modular_raw: np.ndarray) -> dict[str, float]:
    w = np.where(sigma > 0.0, 1.0 / np.square(sigma), 1.0)
    denom = float(np.sum(w * modular_raw * modular_raw))
    alpha = float(np.sum(w * values * modular_raw) / denom) if denom > 0.0 else float("nan")
    residual = values - alpha * modular_raw
    chi2 = float(np.sum(w * residual * residual))
    dof = max(int(values.size) - 1, 1)
    return {
        "alpha": alpha,
        "rms": float(np.sqrt(np.mean(residual * residual))),
        "chi2_dof": float(chi2 / dof),
    }


def _stable_seed(label: str, base_seed: int) -> int:
    return int(base_seed + sum((idx + 1) * ord(ch) for idx, ch in enumerate(label)) % 1_000_000)


def _compute_scale_metrics_with_uncertainty(
    method_manifest_path: str,
    *,
    n_resamples: int,
    base_seed: int,
) -> dict[str, Any]:
    manifest = _load_json(method_manifest_path)
    raw_rows = _load_dat_rows(str(manifest["pointwise_raw"]))
    modular_rows = _load_dat_rows(str(manifest["modular_raw"]))

    modular_lookup = {
        int(row["point_id"]): float(row["modular_raw"])
        for row in modular_rows
        if int(row["is_origin"]) == 0 and np.isfinite(float(row["modular_raw"]))
    }

    by_scale: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        point_id = int(row["point_id"])
        if point_id not in modular_lookup:
            continue
        by_scale[int(row["scale"])] += [row]

    rng = np.random.default_rng(_stable_seed(method_manifest_path, base_seed))
    metrics: list[dict[str, Any]] = []
    for scale in sorted(by_scale):
        rows = by_scale[scale]
        values = np.asarray([float(row["conn"]) for row in rows], dtype=float)
        sigma = np.asarray([float(row["conn_err"]) for row in rows], dtype=float)
        modular_raw = np.asarray(
            [float(modular_lookup[int(row["point_id"])]) for row in rows],
            dtype=float,
        )

        fitted = _fit_metrics(values, sigma, modular_raw)
        rms_q16 = fitted["rms"]
        rms_q84 = fitted["rms"]
        chi2_q16 = fitted["chi2_dof"]
        chi2_q84 = fitted["chi2_dof"]

        if n_resamples > 0:
            draws = rng.normal(loc=values, scale=sigma, size=(n_resamples, values.size))
            w = np.where(sigma > 0.0, 1.0 / np.square(sigma), 1.0)
            denom = float(np.sum(w * modular_raw * modular_raw))
            if denom > 0.0:
                alpha_draws = np.sum(draws * (w * modular_raw), axis=1) / denom
                residual_draws = draws - alpha_draws[:, None] * modular_raw[None, :]
                rms_draws = np.sqrt(np.mean(residual_draws * residual_draws, axis=1))
                chi2_draws = np.sum(w * residual_draws * residual_draws, axis=1) / max(values.size - 1, 1)
                rms_q16, rms_q84 = [float(x) for x in np.percentile(rms_draws, [16.0, 84.0])]
                chi2_q16, chi2_q84 = [float(x) for x in np.percentile(chi2_draws, [16.0, 84.0])]

        metrics.append(
            {
                "scale": int(scale),
                "family_size": int(rows[0]["family_size"]),
                "alpha": fitted["alpha"],
                "rms": fitted["rms"],
                "chi2_dof": fitted["chi2_dof"],
                "rms_point_estimate": fitted["rms"],
                "chi2_point_estimate": fitted["chi2_dof"],
                "rms_err_low": max(fitted["rms"] - rms_q16, 0.0),
                "rms_err_high": max(rms_q84 - fitted["rms"], 0.0),
                "chi2_err_low": max(fitted["chi2_dof"] - chi2_q16, 0.0),
                "chi2_err_high": max(chi2_q84 - fitted["chi2_dof"], 0.0),
            }
        )

    return {
        "benchmark_id": str(manifest["benchmark_id"]),
        "method": str(manifest["method"]),
        "metrics": metrics,
    }


def _plot_benchmark(
    benchmark_manifest_path: str,
    *,
    output_path: str | None,
    n_resamples: int,
    base_seed: int,
) -> str:
    benchmark_manifest = _load_json(benchmark_manifest_path)
    benchmark_id = str(benchmark_manifest["benchmark_id"])

    twisted_manifest = str(benchmark_manifest["methods"]["twisted"])
    untwisted_manifest = str(benchmark_manifest["methods"]["untwisted"])
    twisted = _compute_scale_metrics_with_uncertainty(
        twisted_manifest,
        n_resamples=n_resamples,
        base_seed=base_seed,
    )
    untwisted = _compute_scale_metrics_with_uncertainty(
        untwisted_manifest,
        n_resamples=n_resamples,
        base_seed=base_seed,
    )
    twisted_labels = _scale_labels(twisted_manifest)
    untwisted_labels = _scale_labels(untwisted_manifest)

    fig, (ax_rms, ax_chi2) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    fig.suptitle(f"{benchmark_id}: residual scaling with lattice size", fontsize=14, y=0.98)

    style = {
        "twisted": {"color": "tab:blue", "marker": "o"},
        "untwisted": {"color": "tab:orange", "marker": "s"},
    }

    for data, labels in ((twisted, twisted_labels), (untwisted, untwisted_labels)):
        method = str(data["method"])
        color = style[method]["color"]
        marker = style[method]["marker"]
        scales = np.asarray([int(item["scale"]) for item in data["metrics"]], dtype=float)
        rms = np.asarray([float(item["rms"]) for item in data["metrics"]], dtype=float)
        chi2 = np.asarray([float(item["chi2_dof"]) for item in data["metrics"]], dtype=float)
        rms_err = np.asarray(
            [
                [float(item["rms_err_low"]) for item in data["metrics"]],
                [float(item["rms_err_high"]) for item in data["metrics"]],
            ],
            dtype=float,
        )
        chi2_err = np.asarray(
            [
                [float(item["chi2_err_low"]) for item in data["metrics"]],
                [float(item["chi2_err_high"]) for item in data["metrics"]],
            ],
            dtype=float,
        )

        ax_rms.errorbar(
            scales,
            rms,
            yerr=rms_err,
            marker=marker,
            color=color,
            linewidth=1.8,
            markersize=5,
            capsize=3,
            label=method,
        )
        ax_chi2.errorbar(
            scales,
            chi2,
            yerr=chi2_err,
            marker=marker,
            color=color,
            linewidth=1.8,
            markersize=5,
            capsize=3,
            label=method,
        )

        for x, y in zip(scales, rms):
            label = labels.get(int(x), str(int(x)))
            ax_rms.annotate(label, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7, color=color)

        for x, y in zip(scales, chi2):
            label = labels.get(int(x), str(int(x)))
            ax_chi2.annotate(label, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7, color=color)

    ax_rms.set_title("Amplitude-fitted RMS residual")
    ax_rms.set_xlabel("scale factor relative to smallest lattice")
    ax_rms.set_ylabel("RMS(conn - alpha g_mod)")
    ax_rms.set_yscale("log")
    ax_rms.grid(True, alpha=0.3)
    ax_rms.legend()

    ax_chi2.set_title("Weighted chi2/dof residual")
    ax_chi2.set_xlabel("scale factor relative to smallest lattice")
    ax_chi2.set_ylabel("chi2 / dof")
    ax_chi2.set_yscale("log")
    ax_chi2.grid(True, alpha=0.3)
    ax_chi2.legend()

    fig.text(
        0.5,
        0.01,
        "Error bars: 16-84% intervals from independent Gaussian resampling of pointwise conn_err; near zero they can become upper-only because the residual metrics are positive-definite. Cross-point covariance is not included.",
        ha="center",
        va="bottom",
        fontsize=8,
    )

    fig.tight_layout(rect=[0.0, 0.05, 1.0, 0.93])

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(benchmark_manifest_path),
            f"{benchmark_id}_residual_scaling.png",
        )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"wrote {output_path}")
    for data in (twisted, untwisted):
        method = str(data["method"])
        scales = [int(item["scale"]) for item in data["metrics"]]
        rms = [float(item["rms"]) for item in data["metrics"]]
        chi2 = [float(item["chi2_dof"]) for item in data["metrics"]]
        rms_point = [float(item["rms_point_estimate"]) for item in data["metrics"]]
        chi2_point = [float(item["chi2_point_estimate"]) for item in data["metrics"]]
        rms_err = [
            (
                float(item["rms_err_low"]),
                float(item["rms_err_high"]),
            )
            for item in data["metrics"]
        ]
        chi2_err = [
            (
                float(item["chi2_err_low"]),
                float(item["chi2_err_high"]),
            )
            for item in data["metrics"]
        ]
        print(
            f"{method}: scales={scales} rms={rms} rms_point={rms_point} rms_err={rms_err} "
            f"chi2_dof={chi2} chi2_point={chi2_point} chi2_err={chi2_err}"
        )
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot how modular residuals scale with lattice size for twisted and untwisted families."
    )
    parser.add_argument("--benchmark-manifest", required=True, help="Path to manifest_geometry_*.json")
    parser.add_argument("--output", default=None, help="Optional output PNG path")
    parser.add_argument(
        "--n-resamples",
        type=int,
        default=400,
        help="Number of Gaussian resamples used to estimate error bars from pointwise conn_err.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Base RNG seed for the resampling-based error bars.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _plot_benchmark(
        benchmark_manifest_path=os.path.abspath(args.benchmark_manifest),
        output_path=os.path.abspath(args.output) if args.output else None,
        n_resamples=args.n_resamples,
        base_seed=args.seed,
    )


if __name__ == "__main__":
    main()