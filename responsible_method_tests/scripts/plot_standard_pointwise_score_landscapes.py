#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fit_standard_acute456_all_base_points import (
    _build_summary_rows,
    _compute_orbit_reduced_score_jackknife,
    _compute_orbit_reduced_score,
    _compute_summary_score_jackknife,
    _compute_summary_score,
)
from plot_standard_acute456_center_fss import STANDARD_ROOT


DEFAULT_OUTPUTS = {
    "target": os.path.join(STANDARD_ROOT, "standard_pointwise_score_landscapes.png"),
    "calibrated_target": os.path.join(STANDARD_ROOT, "standard_pointwise_calibrated_target_score_landscapes.png"),
    "holdout": os.path.join(STANDARD_ROOT, "standard_pointwise_holdout_score_landscapes.png"),
}

NORMALIZATION_MODES: dict[str, dict[str, str]] = {
    "raw": {
        "title_suffix": "",
        "colorbar_suffix": "",
        "file_suffix": "",
    },
    "anchor_ratio": {
        "title_suffix": " (anchor-ratio normalized)",
        "colorbar_suffix": " (anchor-ratio normalized)",
        "file_suffix": "_anchor_ratio",
    },
    "l8_ratio": {
        "title_suffix": " (L8-ratio normalized)",
        "colorbar_suffix": " (L8-ratio normalized)",
        "file_suffix": "_l8_ratio",
    },
}

SCORE_MODES: dict[str, dict[str, str]] = {
    "target": {
        "z_key": "target_z",
        "abs_delta_key": "target_abs_delta",
        "figure_title": "Orbit-reduced pointwise target-score landscapes",
        "colorbar_label": "orbit-reduced RMS target z",
        "table_suffix": "pointwise_score_landscape.tsv",
        "score_label": "S_target",
    },
    "calibrated_target": {
        "z_key": "target_calibrated_z",
        "abs_delta_key": "target_abs_delta",
        "figure_title": "Holdout-calibrated pointwise target-score landscapes",
        "colorbar_label": "orbit-reduced RMS calibrated target z",
        "table_suffix": "pointwise_calibrated_target_landscape.tsv",
        "score_label": "S_cal",
    },
    "holdout": {
        "z_key": "holdout_z",
        "abs_delta_key": "holdout_abs_delta",
        "figure_title": "Largest-size holdout pointwise score landscapes",
        "colorbar_label": "orbit-reduced RMS holdout z",
        "table_suffix": "pointwise_holdout_score_landscape.tsv",
        "score_label": "S_holdout",
    },
}

FAMILY_CONFIGS: dict[str, dict[str, Any]] = {
    "iso111": {
        "title": "1-1-1",
        "untwisted_root": os.path.join(STANDARD_ROOT, "data", "iso111", "untwisted"),
        "twisted_dat": os.path.join(
            STANDARD_ROOT,
            "data",
            "iso111",
            "twisted",
            "reference",
            "Lx144_Ly144_Tx0_Ty0",
            "two_point_all_to_all.dat",
        ),
        "twisted_lattice": (144, 144, 0, 0),
        "truth": (1.0, 1.0),
    },
    "acute456": {
        "title": "4-5-6",
        "untwisted_root": os.path.join(STANDARD_ROOT, "data", "acute456", "untwisted"),
        "twisted_dat": os.path.join(
            STANDARD_ROOT,
            "data",
            "acute456",
            "twisted",
            "reference",
            "Lx144_Ly144_Tx72_Ty24",
            "two_point_all_to_all.dat",
        ),
        "twisted_lattice": (144, 144, 72, 24),
        "truth": (4.702782819756, 7.353910143333),
    },
}

_CANDIDATE_RE = re.compile(r"^r1_([0-9mp]+)__r2_([0-9mp]+)$")


def _decode_ratio_token(token: str) -> float:
    return float(token.replace("m", "-").replace("p", "."))


def _candidate_ratios(candidate_dir_name: str) -> tuple[float, float]:
    match = _CANDIDATE_RE.match(candidate_dir_name)
    if match is None:
        raise ValueError(f"candidate directory does not encode (r1,r2): {candidate_dir_name}")
    return _decode_ratio_token(match.group(1)), _decode_ratio_token(match.group(2))


def _with_mode_suffix(path: str, normalization_mode: str) -> str:
    mode_cfg = NORMALIZATION_MODES[str(normalization_mode)]
    if not mode_cfg["file_suffix"]:
        return path
    root, ext = os.path.splitext(path)
    return f"{root}{mode_cfg['file_suffix']}{ext}"


def _drop_origin_rows(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        row
        for row in summary_rows
        if abs(float(row["a_wrap"])) >= 1.0e-12 or abs(float(row["b_wrap"])) >= 1.0e-12
    ]


def _write_table(output_path: str, rows: list[dict[str, Any]], *, score_mode: str, normalization_mode: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "candidate_dir",
                "r1",
                "r2",
                "score_mode",
                "normalization_mode",
                "orbit_score_rms_z",
                "orbit_score_jackknife_sigma",
                "orbit_mean_rms_z",
                "orbit_mean_abs_z",
                "orbit_max_abs_z",
                "orbit_n_orbits",
                "orbit_n_points",
                "point_rms_z",
                "point_rms_z_jackknife_sigma",
                "point_mean_abs_z",
                "point_max_abs_z",
                "point_origin_dropped_rms_z",
                "point_origin_dropped_rms_z_jackknife_sigma",
                "point_origin_inclusive_rms_z",
                "point_origin_inclusive_rms_z_jackknife_sigma",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    str(row["candidate_dir"]),
                    f"{float(row['r1']):.6f}",
                    f"{float(row['r2']):.6f}",
                    score_mode,
                    normalization_mode,
                    f"{float(row['score']):.10f}",
                    f"{float(row['score_jackknife_sigma']):.10f}",
                    f"{float(row['mean_orbit_rms_z']):.10f}",
                    f"{float(row['mean_orbit_abs_z']):.10f}",
                    f"{float(row['max_abs_z']):.10f}",
                    str(int(row["n_orbits"])),
                    str(int(row["n_points"])),
                    f"{float(row['point_rms_z']):.10f}",
                    f"{float(row['point_rms_z_jackknife_sigma']):.10f}",
                    f"{float(row['point_mean_abs_z']):.10f}",
                    f"{float(row['point_max_abs_z']):.10f}",
                    f"{float(row['point_origin_dropped_rms_z']):.10f}",
                    f"{float(row['point_origin_dropped_rms_z_jackknife_sigma']):.10f}",
                    f"{float(row['point_origin_inclusive_rms_z']):.10f}",
                    f"{float(row['point_origin_inclusive_rms_z_jackknife_sigma']):.10f}",
                ]
            )


def _evaluate_family(
    family: str,
    *,
    score_mode: str,
    normalization_mode: str,
    anchor_m: int,
    anchor_n: int,
) -> dict[str, Any]:
    if family not in FAMILY_CONFIGS:
        raise KeyError(f"unknown family: {family}")
    if score_mode not in SCORE_MODES:
        raise KeyError(f"unknown score mode: {score_mode}")
    if normalization_mode not in NORMALIZATION_MODES:
        raise KeyError(f"unknown normalization mode: {normalization_mode}")
    cfg = dict(FAMILY_CONFIGS[family])
    mode_cfg = dict(SCORE_MODES[score_mode])
    untwisted_root = os.path.abspath(str(cfg["untwisted_root"]))
    twisted_dat = os.path.abspath(str(cfg["twisted_dat"]))
    twisted_lattice = tuple(int(value) for value in cfg["twisted_lattice"])

    rows: list[dict[str, Any]] = []
    for candidate_dir_name in sorted(os.listdir(untwisted_root)):
        candidate_dir = os.path.join(untwisted_root, candidate_dir_name)
        if not os.path.isdir(candidate_dir):
            continue
        r1, r2 = _candidate_ratios(candidate_dir_name)
        summary_rows = _build_summary_rows(
            untwisted_dir=candidate_dir,
            twisted_dat=twisted_dat,
            twisted_lattice=twisted_lattice,
            untwisted_embedding_cycles=(0, 1),
            twisted_embedding_cycles=(0, 1),
            include_origin=True,
            normalization_mode=normalization_mode,
            anchor_m=anchor_m,
            anchor_n=anchor_n,
        )
        origin_dropped_rows = _drop_origin_rows(summary_rows)
        orbit_score = _compute_orbit_reduced_score_jackknife(
            origin_dropped_rows,
            include_origin=False,
            z_key=str(mode_cfg["z_key"]),
        )
        point_score = _compute_summary_score_jackknife(
            origin_dropped_rows,
            z_key=str(mode_cfg["z_key"]),
            abs_delta_key=str(mode_cfg["abs_delta_key"]),
        )
        point_inclusive_score = _compute_summary_score_jackknife(
            summary_rows,
            z_key=str(mode_cfg["z_key"]),
            abs_delta_key=str(mode_cfg["abs_delta_key"]),
        )
        rows.append(
            {
                "candidate_dir": candidate_dir_name,
                "r1": float(r1),
                "r2": float(r2),
                "score": float(orbit_score["rms_z"]),
                "score_jackknife_sigma": float(orbit_score["jackknife_sigma"]),
                "mean_orbit_rms_z": float(orbit_score["mean_orbit_rms_z"]),
                "mean_orbit_abs_z": float(orbit_score["mean_orbit_abs_z"]),
                "max_abs_z": float(orbit_score["max_abs_z"]),
                "n_orbits": int(orbit_score["n_orbits"]),
                "n_points": int(orbit_score["n_points"]),
                "point_rms_z": float(point_score["rms_z"]),
                "point_rms_z_jackknife_sigma": float(point_score["jackknife_sigma"]),
                "point_mean_abs_z": float(point_score["mean_abs_z"]),
                "point_max_abs_z": float(point_score["max_abs_z"]),
                "point_origin_dropped_rms_z": float(point_score["rms_z"]),
                "point_origin_dropped_rms_z_jackknife_sigma": float(point_score["jackknife_sigma"]),
                "point_origin_inclusive_rms_z": float(point_inclusive_score["rms_z"]),
                "point_origin_inclusive_rms_z_jackknife_sigma": float(point_inclusive_score["jackknife_sigma"]),
            }
        )

    rows.sort(key=lambda row: (float(row["score"]), float(row["r1"]), float(row["r2"])))
    if not rows:
        raise ValueError(f"no candidate rows found for family {family}")

    r1_values = sorted({float(row["r1"]) for row in rows})
    r2_values = sorted({float(row["r2"]) for row in rows})
    Z = np.full((len(r1_values), len(r2_values)), np.nan, dtype=float)
    for row in rows:
        i = r1_values.index(float(row["r1"]))
        j = r2_values.index(float(row["r2"]))
        Z[i, j] = float(row["score"])

    table_name = _with_mode_suffix(f"{family}_{mode_cfg['table_suffix']}", normalization_mode)
    table_path = os.path.join(STANDARD_ROOT, table_name)
    _write_table(table_path, rows, score_mode=score_mode, normalization_mode=normalization_mode)
    return {
        "family": family,
        "title": str(cfg["title"]),
        "score_mode": score_mode,
        "normalization_mode": normalization_mode,
        "score_label": str(mode_cfg["score_label"]),
        "truth": tuple(float(value) for value in cfg["truth"]),
        "rows": rows,
        "grid_r1": r1_values,
        "grid_r2": r2_values,
        "grid_score": Z,
        "best": rows[0],
        "table_path": table_path,
    }


def _grid_extent(values: list[float]) -> tuple[float, float]:
    if len(values) <= 1:
        return values[0] - 0.5, values[0] + 0.5
    step = float(np.median(np.diff(np.asarray(values, dtype=float))))
    return float(values[0] - 0.5 * step), float(values[-1] + 0.5 * step)


def _plot_landscapes(family_payloads: list[dict[str, Any]], *, output_path: str, annotate_values: bool) -> str:
    n_families = len(family_payloads)
    fig, axes = plt.subplots(1, n_families, figsize=(6.4 * n_families, 5.6), squeeze=False)
    axes_flat = list(axes.ravel())

    all_scores = np.concatenate(
        [payload["grid_score"][np.isfinite(payload["grid_score"])] for payload in family_payloads]
    )
    vmin = float(np.nanmin(all_scores)) if all_scores.size else 0.0
    vmax = float(np.nanmax(all_scores)) if all_scores.size else 1.0
    image = None
    for axis, payload in zip(axes_flat, family_payloads):
        r1_values = list(payload["grid_r1"])
        r2_values = list(payload["grid_r2"])
        Z = np.asarray(payload["grid_score"], dtype=float)
        x0, x1 = _grid_extent(r1_values)
        y0, y1 = _grid_extent(r2_values)
        image = axis.imshow(
            Z.T,
            origin="lower",
            extent=[x0, x1, y0, y1],
            aspect="auto",
            cmap="viridis_r",
            vmin=vmin,
            vmax=vmax,
        )
        truth_r1, truth_r2 = payload["truth"]
        best = payload["best"]
        axis.scatter([truth_r1], [truth_r2], marker="x", s=140, color="black", linewidths=2.0, label="truth")
        axis.scatter(
            [float(best["r1"])],
            [float(best["r2"])],
            marker="*",
            s=250,
            color="white",
            edgecolors="black",
            linewidths=0.9,
            label="best",
        )
        if annotate_values:
            for i, r1 in enumerate(r1_values):
                for j, r2 in enumerate(r2_values):
                    score = float(Z[i, j])
                    axis.text(
                        float(r1),
                        float(r2),
                        f"{score:.2f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if score > 0.65 * vmax else "black",
                    )
        axis.set_title(
            f"{payload['title']}\n"
            f"best=({float(best['r1']):.3f}, {float(best['r2']):.3f}), "
            f"{payload['score_label']}={float(best['score']):.3f} +/- {float(best['score_jackknife_sigma']):.3f}",
            fontsize=11,
        )
        axis.set_xlabel("r1 = k1 / k3")
        axis.set_ylabel("r2 = k2 / k3")
        axis.grid(False)
        axis.legend(loc="upper left", fontsize=8, frameon=True)

    score_mode = str(family_payloads[0]["score_mode"])
    normalization_mode = str(family_payloads[0]["normalization_mode"])
    fig.suptitle(
        str(SCORE_MODES[score_mode]["figure_title"]) + str(NORMALIZATION_MODES[normalization_mode]["title_suffix"]),
        fontsize=14,
        y=0.98,
    )
    if image is not None:
        fig.colorbar(
            image,
            ax=axes_flat,
            fraction=0.035,
            pad=0.04,
            label=(
                str(SCORE_MODES[score_mode]["colorbar_label"])
                + str(NORMALIZATION_MODES[normalization_mode]["colorbar_suffix"])
            ),
        )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render orbit-reduced pointwise score landscapes for the standard iso111 and acute456 candidate grids."
    )
    parser.add_argument("--family", action="append", choices=sorted(FAMILY_CONFIGS), help="Family to include; pass once per family")
    parser.add_argument("--score-mode", choices=sorted(SCORE_MODES), default="target")
    parser.add_argument("--normalization-mode", choices=sorted(NORMALIZATION_MODES), default="raw")
    parser.add_argument("--anchor-m", type=int, default=0)
    parser.add_argument("--anchor-n", type=int, default=-1)
    parser.add_argument("--output", default=None)
    parser.add_argument("--no-annotate-values", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    families = args.family if args.family else ["iso111", "acute456"]
    payloads = [
        _evaluate_family(
            family,
            score_mode=str(args.score_mode),
            normalization_mode=str(args.normalization_mode),
            anchor_m=int(args.anchor_m),
            anchor_n=int(args.anchor_n),
        )
        for family in families
    ]
    if args.output:
        output_path = os.path.abspath(args.output)
    else:
        output_path = os.path.abspath(_with_mode_suffix(DEFAULT_OUTPUTS[str(args.score_mode)], str(args.normalization_mode)))
    _plot_landscapes(payloads, output_path=output_path, annotate_values=not bool(args.no_annotate_values))

    print(f"wrote {output_path}")
    for payload in payloads:
        best = payload["best"]
        print(f"wrote {payload['table_path']}")
        print(
            f"family={payload['family']} "
            f"best_r1={float(best['r1']):.6f} "
            f"best_r2={float(best['r2']):.6f} "
            f"score={float(best['score']):.6f} "
            f"score_jackknife_sigma={float(best['score_jackknife_sigma']):.6f}"
        )


if __name__ == "__main__":
    main()