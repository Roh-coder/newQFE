#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MC_ENGINE_ROOT = os.path.join(REPO_ROOT, "K_from_continuum", "lib")
if MC_ENGINE_ROOT not in sys.path:
    sys.path.insert(0, MC_ENGINE_ROOT)

import mc_engine  # noqa: E402


@dataclass
class CandidateEntry:
    label: str
    manifest_path: str
    r1: float
    r2: float
    score: float
    payloads_by_scale: dict[int, dict[str, Any]]


def _load_json(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _payloads_by_scale(manifest: dict[str, Any]) -> dict[int, dict[str, Any]]:
    payloads: dict[int, dict[str, Any]] = {}
    for payload in manifest.get("payloads", []):
        payloads[int(payload["scale"])] = payload
    return payloads


def _load_all_to_all_matrix(data_path: str, field: str) -> tuple[list[int], list[int], np.ndarray]:
    data = mc_engine.load_all_to_all(data_path)
    m_values = sorted({int(m) for (m, _) in data})
    n_values = sorted({int(n) for (_, n) in data})
    matrix = np.full((len(n_values), len(m_values)), np.nan, dtype=float)
    m_index = {value: idx for idx, value in enumerate(m_values)}
    n_index = {value: idx for idx, value in enumerate(n_values)}

    for (m_value, n_value), row in data.items():
        matrix[n_index[int(n_value)], m_index[int(m_value)]] = float(row[field])
    return m_values, n_values, matrix


def _load_all_to_all_value_and_error_matrices(
    data_path: str,
    value_field: str,
    error_field: str,
) -> tuple[list[int], list[int], np.ndarray, np.ndarray]:
    data = mc_engine.load_all_to_all(data_path)
    m_values = sorted({int(m) for (m, _) in data})
    n_values = sorted({int(n) for (_, n) in data})
    value_matrix = np.full((len(n_values), len(m_values)), np.nan, dtype=float)
    error_matrix = np.full((len(n_values), len(m_values)), np.nan, dtype=float)
    m_index = {value: idx for idx, value in enumerate(m_values)}
    n_index = {value: idx for idx, value in enumerate(n_values)}

    for (m_value, n_value), row in data.items():
        row_idx = n_index[int(n_value)]
        col_idx = m_index[int(m_value)]
        value_matrix[row_idx, col_idx] = float(row[value_field])
        error_matrix[row_idx, col_idx] = float(row[error_field])
    return m_values, n_values, value_matrix, error_matrix


def _candidate_entries(scan_payload: dict[str, Any]) -> list[CandidateEntry]:
    entries: list[CandidateEntry] = []
    for row in scan_payload["rows"]:
        manifest_path = os.path.abspath(str(row["candidate_manifest"]))
        manifest = _load_json(manifest_path)
        entries.append(
            CandidateEntry(
                label=str(row["candidate_label"]),
                manifest_path=manifest_path,
                r1=float(row["r1"]),
                r2=float(row["r2"]),
                score=float(row["score"]),
                payloads_by_scale=_payloads_by_scale(manifest),
            )
        )
    return entries


def _common_scales(reference_payloads: dict[int, dict[str, Any]], entries: list[CandidateEntry]) -> list[int]:
    scales = set(reference_payloads)
    for entry in entries:
        scales &= set(entry.payloads_by_scale)
    return sorted(scales)


def _make_reference_figure(
    *,
    matrix: np.ndarray,
    m_values: list[int],
    n_values: list[int],
    scale: int,
    family_size: int,
    output_path: str,
    field_label: str,
) -> str:
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    image = ax.imshow(matrix, origin="lower", aspect="equal", cmap="viridis")
    ax.set_title(f"Reference lattice, scale={scale}, family size={family_size}")
    ax.set_xlabel("m")
    ax.set_ylabel("n")
    ax.set_xticks(range(len(m_values)))
    ax.set_yticks(range(len(n_values)))
    ax.set_xticklabels(m_values, fontsize=7)
    ax.set_yticklabels(n_values, fontsize=7)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label(field_label)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _gallery_layout(entries: list[CandidateEntry]) -> tuple[list[float], list[float], dict[tuple[float, float], CandidateEntry]]:
    r1_values = sorted({entry.r1 for entry in entries})
    r2_values = sorted({entry.r2 for entry in entries}, reverse=True)
    lookup = {(entry.r1, entry.r2): entry for entry in entries}
    return r1_values, r2_values, lookup


def _decorate_axis(
    ax: Any,
    *,
    title: str | None,
    ylabel: str | None,
    is_best: bool,
    is_isotropic: bool,
) -> None:
    if title is not None:
        ax.set_title(title, fontsize=8)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis="both", which="both", length=0, labelbottom=False, labelleft=False)
    ax.set_aspect("equal")
    if is_best:
        for spine in ax.spines.values():
            spine.set_color("red")
            spine.set_linewidth(2.0)
    elif is_isotropic:
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.8)


def _make_gallery(
    *,
    lookup: dict[tuple[float, float], CandidateEntry],
    r1_values: list[float],
    r2_values: list[float],
    scale: int,
    value_matrices: dict[tuple[float, float], np.ndarray],
    output_path: str,
    figure_title: str,
    cmap: str,
    field_label: str,
    reference_matrix: np.ndarray | None = None,
) -> str:
    n_rows = len(r2_values)
    n_cols = len(r1_values)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.4 * n_cols, 2.25 * n_rows), squeeze=False)

    stacked = np.stack([value_matrices[(r1, r2)] for r2 in r2_values for r1 in r1_values], axis=0)
    finite = stacked[np.isfinite(stacked)]
    if reference_matrix is None:
        norm = None
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
    else:
        diff_finite = finite
        diff_limit = float(np.max(np.abs(diff_finite))) if diff_finite.size else 1.0
        vmin = -diff_limit
        vmax = diff_limit
        norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)

    best_entry = min(lookup.values(), key=lambda entry: entry.score)
    image = None
    for row_idx, r2_value in enumerate(r2_values):
        for col_idx, r1_value in enumerate(r1_values):
            ax = axes[row_idx][col_idx]
            matrix = value_matrices[(r1_value, r2_value)]
            image_kwargs: dict[str, Any] = {"origin": "lower", "cmap": cmap}
            if norm is None:
                image_kwargs.update({"vmin": vmin, "vmax": vmax})
            else:
                image_kwargs.update({"norm": norm})
            image = ax.imshow(matrix, **image_kwargs)
            title = f"r1={r1_value:.3f}" if row_idx == 0 else None
            ylabel = f"r2={r2_value:.3f}" if col_idx == 0 else None
            entry = lookup[(r1_value, r2_value)]
            _decorate_axis(
                ax,
                title=title,
                ylabel=ylabel,
                is_best=(entry.r1 == best_entry.r1 and entry.r2 == best_entry.r2),
                is_isotropic=(abs(entry.r1 - 1.0) < 1e-12 and abs(entry.r2 - 1.0) < 1e-12),
            )
            ax.text(
                0.04,
                0.05,
                f"S={entry.score:.2f}",
                transform=ax.transAxes,
                fontsize=7,
                color="white" if reference_matrix is not None else "black",
                bbox={"facecolor": "black" if reference_matrix is not None else "white", "alpha": 0.55, "pad": 1.5, "edgecolor": "none"},
            )
            if row_idx == n_rows - 1:
                ax.set_xlabel(f"r1={r1_value:.3f}", fontsize=8)

    fig.suptitle(f"{figure_title} (scale={scale})", fontsize=14)
    fig.subplots_adjust(top=0.92, right=0.92, wspace=0.04, hspace=0.08)
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.015, pad=0.01)
        colorbar.set_label(field_label)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _panel_limits(matrix: np.ndarray, *, symmetric: bool) -> tuple[float, float]:
    finite = matrix[np.isfinite(matrix)]
    if finite.size == 0:
        return (-1.0, 1.0) if symmetric else (0.0, 1.0)
    if symmetric:
        limit = float(np.max(np.abs(finite)))
        if limit <= 0.0:
            limit = 1.0
        return -limit, limit
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if np.isclose(vmin, vmax):
        delta = 1.0 if np.isclose(vmin, 0.0) else abs(vmin) * 0.1
        vmin -= delta
        vmax += delta
    return vmin, vmax


def _make_best_candidate_comparison(
    *,
    reference_matrix: np.ndarray,
    reference_error_matrix: np.ndarray,
    candidate_matrix: np.ndarray,
    candidate_error_matrix: np.ndarray,
    m_values: list[int],
    n_values: list[int],
    scale: int,
    family_size: int,
    candidate_label: str,
    r1: float,
    r2: float,
    score: float,
    output_path: str,
    field_label: str,
    title_prefix: str,
) -> str:
    difference_matrix = candidate_matrix - reference_matrix
    combined_error = np.sqrt(reference_error_matrix**2 + candidate_error_matrix**2)
    zscore_matrix = np.full_like(difference_matrix, np.nan)
    valid = np.isfinite(difference_matrix) & np.isfinite(combined_error) & (combined_error > 0.0)
    zscore_matrix[valid] = difference_matrix[valid] / combined_error[valid]

    value_vmin = min(
        _panel_limits(reference_matrix, symmetric=False)[0],
        _panel_limits(candidate_matrix, symmetric=False)[0],
    )
    value_vmax = max(
        _panel_limits(reference_matrix, symmetric=False)[1],
        _panel_limits(candidate_matrix, symmetric=False)[1],
    )
    diff_vmin, diff_vmax = _panel_limits(difference_matrix, symmetric=True)
    z_vmin, z_vmax = _panel_limits(zscore_matrix, symmetric=True)

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.8), constrained_layout=True)
    panel_specs = [
        (reference_matrix, "Reference", "viridis", None, value_vmin, value_vmax, field_label),
        (candidate_matrix, f"Test: {candidate_label}", "viridis", None, value_vmin, value_vmax, field_label),
        (
            difference_matrix,
            "Test - Reference",
            "coolwarm",
            TwoSlopeNorm(vcenter=0.0, vmin=diff_vmin, vmax=diff_vmax),
            None,
            None,
            f"delta {field_label}",
        ),
        (
            zscore_matrix,
            "(Test - Reference) / sigma_combined",
            "coolwarm",
            TwoSlopeNorm(vcenter=0.0, vmin=z_vmin, vmax=z_vmax),
            None,
            None,
            "z-score",
        ),
    ]

    for idx, (matrix, title, cmap, norm, vmin, vmax, colorbar_label) in enumerate(panel_specs):
        ax = axes[idx]
        image_kwargs: dict[str, Any] = {"origin": "lower", "aspect": "equal", "cmap": cmap}
        if norm is None:
            image_kwargs.update({"vmin": vmin, "vmax": vmax})
        else:
            image_kwargs.update({"norm": norm})
        image = ax.imshow(matrix, **image_kwargs)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("m")
        if idx == 0:
            ax.set_ylabel("n")
        ax.set_xticks(range(len(m_values)))
        ax.set_yticks(range(len(n_values)))
        ax.set_xticklabels(m_values, fontsize=7)
        ax.set_yticklabels(n_values, fontsize=7)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label(colorbar_label)

    fig.suptitle(
        (
            f"{title_prefix} vs reference {field_label}, scale={scale}, family size={family_size}\n"
            f"candidate={candidate_label}, r1={r1:.3f}, r2={r2:.3f}, score={score:.4f}"
        ),
        fontsize=13,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _make_head_to_head_comparison(
    *,
    reference_matrix: np.ndarray,
    reference_error_matrix: np.ndarray,
    candidate_specs: list[dict[str, Any]],
    m_values: list[int],
    n_values: list[int],
    scale: int,
    family_size: int,
    output_path: str,
    field_label: str,
) -> str:
    if len(candidate_specs) != 2:
        raise ValueError("head-to-head comparison requires exactly two candidate specs")

    value_vmin, value_vmax = _panel_limits(reference_matrix, symmetric=False)
    diff_limit = 0.0
    z_limit = 0.0
    enriched_specs: list[dict[str, Any]] = []
    for spec in candidate_specs:
        candidate_matrix = spec["candidate_matrix"]
        candidate_error_matrix = spec["candidate_error_matrix"]
        difference_matrix = candidate_matrix - reference_matrix
        combined_error = np.sqrt(reference_error_matrix**2 + candidate_error_matrix**2)
        zscore_matrix = np.full_like(difference_matrix, np.nan)
        valid = np.isfinite(difference_matrix) & np.isfinite(combined_error) & (combined_error > 0.0)
        zscore_matrix[valid] = difference_matrix[valid] / combined_error[valid]

        cand_vmin, cand_vmax = _panel_limits(candidate_matrix, symmetric=False)
        value_vmin = min(value_vmin, cand_vmin)
        value_vmax = max(value_vmax, cand_vmax)
        diff_limit = max(diff_limit, abs(_panel_limits(difference_matrix, symmetric=True)[0]))
        z_limit = max(z_limit, abs(_panel_limits(zscore_matrix, symmetric=True)[0]))
        enriched_specs.append(
            {
                **spec,
                "difference_matrix": difference_matrix,
                "zscore_matrix": zscore_matrix,
            }
        )

    if diff_limit <= 0.0:
        diff_limit = 1.0
    if z_limit <= 0.0:
        z_limit = 1.0

    fig, axes = plt.subplots(2, 4, figsize=(18.0, 9.0), constrained_layout=True)
    diff_norm = TwoSlopeNorm(vcenter=0.0, vmin=-diff_limit, vmax=diff_limit)
    z_norm = TwoSlopeNorm(vcenter=0.0, vmin=-z_limit, vmax=z_limit)

    reference_image = axes[0][0].imshow(
        reference_matrix,
        origin="lower",
        aspect="equal",
        cmap="viridis",
        vmin=value_vmin,
        vmax=value_vmax,
    )
    axes[0][0].set_title("Reference", fontsize=10)
    axes[1][0].imshow(
        reference_matrix,
        origin="lower",
        aspect="equal",
        cmap="viridis",
        vmin=value_vmin,
        vmax=value_vmax,
    )
    axes[1][0].set_title("Reference", fontsize=10)

    candidate_image = None
    diff_image = None
    zscore_image = None
    for row_idx, spec in enumerate(enriched_specs):
        candidate_image = axes[row_idx][1].imshow(
            spec["candidate_matrix"],
            origin="lower",
            aspect="equal",
            cmap="viridis",
            vmin=value_vmin,
            vmax=value_vmax,
        )
        diff_image = axes[row_idx][2].imshow(
            spec["difference_matrix"],
            origin="lower",
            aspect="equal",
            cmap="coolwarm",
            norm=diff_norm,
        )
        zscore_image = axes[row_idx][3].imshow(
            spec["zscore_matrix"],
            origin="lower",
            aspect="equal",
            cmap="coolwarm",
            norm=z_norm,
        )
        axes[row_idx][1].set_title(
            f"Test: {spec['candidate_label']}\nr1={spec['r1']:.3f}, r2={spec['r2']:.3f}, score={spec['score']:.4f}",
            fontsize=10,
        )
        axes[row_idx][2].set_title("Test - Reference", fontsize=10)
        axes[row_idx][3].set_title("(Test - Reference) / sigma_combined", fontsize=10)

    for row_idx, spec in enumerate(enriched_specs):
        for col_idx in range(4):
            ax = axes[row_idx][col_idx]
            ax.set_xticks(range(len(m_values)))
            ax.set_yticks(range(len(n_values)))
            ax.set_xticklabels(m_values, fontsize=7)
            ax.set_yticklabels(n_values, fontsize=7)
            ax.set_xlabel("m")
            if col_idx == 0:
                ax.set_ylabel(f"n\n{spec['candidate_label']}", fontsize=9)

    value_colorbar = fig.colorbar(reference_image, ax=axes[:, :2], fraction=0.018, pad=0.02)
    value_colorbar.set_label(field_label)
    diff_colorbar = fig.colorbar(diff_image, ax=axes[:, 2], fraction=0.03, pad=0.02)
    diff_colorbar.set_label(f"delta {field_label}")
    zscore_colorbar = fig.colorbar(zscore_image, ax=axes[:, 3], fraction=0.03, pad=0.02)
    zscore_colorbar.set_label("z-score")

    fig.suptitle(
        f"Head-to-head vs reference {field_label}, scale={scale}, family size={family_size}",
        fontsize=14,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _write_summary(
    *,
    output_path: str,
    reference_manifest_path: str,
    scan_json_path: str,
    scale_outputs: list[dict[str, Any]],
    all_candidate_comparisons_dir: str | None,
) -> str:
    lines = [
        "# All-to-all correlator galleries",
        "",
        f"- scan json: {scan_json_path}",
        f"- reference manifest: {reference_manifest_path}",
        f"- all-candidate comparisons: {all_candidate_comparisons_dir or 'not generated'}",
        "",
        "## Generated outputs",
    ]
    for item in scale_outputs:
        lines.extend(
            [
                "",
                f"### scale {item['scale']} (family size {item['family_size']})",
                f"- reference: {item['reference_plot']}",
                f"- gallery: {item['gallery_plot']}",
                f"- difference gallery: {item['difference_plot']}",
                f"- best-candidate 4-panel comparison: {item['best_candidate_comparison_plot']}",
            ]
        )
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render untwisted all-to-all correlator galleries against the scan reference lattice."
    )
    parser.add_argument("--scan-json", required=True, help="Path to geometry_*_coupling_fit_scan.json")
    parser.add_argument("--output-dir", default=None, help="Optional output directory")
    parser.add_argument(
        "--scale",
        action="append",
        type=int,
        dest="scales",
        help="Optional scale to render; pass multiple times to restrict the default common-scale set",
    )
    parser.add_argument(
        "--field",
        default="conn",
        choices=("conn", "corr"),
        help="All-to-all field to render. 'conn' maps to corr_conn, 'corr' maps to corr.",
    )
    parser.add_argument(
        "--emit-all-candidate-comparisons",
        action="store_true",
        help="Also emit per-candidate 4-panel comparison PNGs for every scan point.",
    )
    parser.add_argument(
        "--all-candidate-dir",
        default=None,
        help="Optional output directory for per-candidate comparison PNGs. Defaults under the main output directory.",
    )
    parser.add_argument(
        "--head-to-head-candidate",
        action="append",
        dest="head_to_head_candidates",
        default=None,
        help="Candidate label to include in a direct head-to-head comparison; pass exactly twice.",
    )
    parser.add_argument(
        "--head-to-head-output",
        default=None,
        help="Optional output path for the direct head-to-head comparison PNG.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scan_json_path = os.path.abspath(args.scan_json)
    scan_payload = _load_json(scan_json_path)
    reference_manifest_path = os.path.abspath(str(scan_payload["target"]))
    reference_manifest = _load_json(reference_manifest_path)
    reference_payloads = _payloads_by_scale(reference_manifest)
    entries = _candidate_entries(scan_payload)
    entries_by_label = {entry.label: entry for entry in entries}
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else os.path.join(os.path.dirname(scan_json_path), "all_to_all_vs_reference")
    os.makedirs(output_dir, exist_ok=True)
    all_candidate_comparisons_dir = None
    if args.emit_all_candidate_comparisons:
        all_candidate_comparisons_dir = (
            os.path.abspath(args.all_candidate_dir)
            if args.all_candidate_dir
            else os.path.join(output_dir, f"all_candidate_comparisons_{args.field}")
        )
        os.makedirs(all_candidate_comparisons_dir, exist_ok=True)

    requested_scales = sorted(set(args.scales)) if args.scales else None
    common_scales = _common_scales(reference_payloads, entries)
    if requested_scales is not None:
        missing = [scale for scale in requested_scales if scale not in common_scales]
        if missing:
            raise ValueError(f"requested scales are not shared by reference and every candidate: {missing}")
        scales = requested_scales
    else:
        scales = common_scales

    if not scales:
        raise ValueError("no common scales found between the reference and all candidate manifests")

    head_to_head_labels = args.head_to_head_candidates
    if head_to_head_labels is not None:
        if len(head_to_head_labels) != 2:
            raise ValueError("--head-to-head-candidate must be provided exactly twice")
        missing_labels = [label for label in head_to_head_labels if label not in entries_by_label]
        if missing_labels:
            raise ValueError(f"head-to-head candidates not found in scan: {missing_labels}")

    field_key = "conn" if args.field == "conn" else "corr"
    error_field_key = "conn_err" if args.field == "conn" else "err"
    field_label = "connected correlator" if args.field == "conn" else "correlator"
    r1_values, r2_values, lookup = _gallery_layout(entries)
    best_entry = min(entries, key=lambda entry: entry.score)

    scale_outputs: list[dict[str, Any]] = []
    for scale in scales:
        reference_payload = reference_payloads[scale]
        family_size = int(reference_payload["family_size"])
        m_values, n_values, reference_matrix, reference_error_matrix = _load_all_to_all_value_and_error_matrices(
            str(reference_payload["data_path"]),
            field_key,
            error_field_key,
        )
        value_matrices: dict[tuple[float, float], np.ndarray] = {}
        diff_matrices: dict[tuple[float, float], np.ndarray] = {}
        candidate_error_matrices: dict[tuple[float, float], np.ndarray] = {}

        for entry in entries:
            payload = entry.payloads_by_scale[scale]
            cand_m_values, cand_n_values, candidate_matrix, candidate_error_matrix = _load_all_to_all_value_and_error_matrices(
                str(payload["data_path"]),
                field_key,
                error_field_key,
            )
            if cand_m_values != m_values or cand_n_values != n_values:
                raise ValueError(
                    f"candidate grid shape mismatch at scale {scale}: {entry.manifest_path}"
                )
            value_matrices[(entry.r1, entry.r2)] = candidate_matrix
            diff_matrices[(entry.r1, entry.r2)] = candidate_matrix - reference_matrix
            candidate_error_matrices[(entry.r1, entry.r2)] = candidate_error_matrix

        reference_plot = os.path.join(output_dir, f"scale{scale:02d}_reference_{args.field}.png")
        gallery_plot = os.path.join(output_dir, f"scale{scale:02d}_untwisted_gallery_{args.field}.png")
        difference_plot = os.path.join(output_dir, f"scale{scale:02d}_untwisted_minus_reference_{args.field}.png")
        best_candidate_comparison_plot = os.path.join(
            output_dir,
            f"scale{scale:02d}_best_candidate_comparison_{args.field}.png",
        )

        _make_reference_figure(
            matrix=reference_matrix,
            m_values=m_values,
            n_values=n_values,
            scale=scale,
            family_size=family_size,
            output_path=reference_plot,
            field_label=field_label,
        )
        _make_gallery(
            lookup=lookup,
            r1_values=r1_values,
            r2_values=r2_values,
            scale=scale,
            value_matrices=value_matrices,
            output_path=gallery_plot,
            figure_title=f"Untwisted {field_label} gallery vs reference family size {family_size}",
            cmap="viridis",
            field_label=field_label,
        )
        _make_gallery(
            lookup=lookup,
            r1_values=r1_values,
            r2_values=r2_values,
            scale=scale,
            value_matrices=diff_matrices,
            output_path=difference_plot,
            figure_title=f"Untwisted minus reference {field_label}",
            cmap="coolwarm",
            field_label=f"delta {field_label}",
            reference_matrix=reference_matrix,
        )
        _make_best_candidate_comparison(
            reference_matrix=reference_matrix,
            reference_error_matrix=reference_error_matrix,
            candidate_matrix=value_matrices[(best_entry.r1, best_entry.r2)],
            candidate_error_matrix=candidate_error_matrices[(best_entry.r1, best_entry.r2)],
            m_values=m_values,
            n_values=n_values,
            scale=scale,
            family_size=family_size,
            candidate_label=best_entry.label,
            r1=best_entry.r1,
            r2=best_entry.r2,
            score=best_entry.score,
            output_path=best_candidate_comparison_plot,
            field_label=field_label,
            title_prefix="Best untwisted",
        )

        if head_to_head_labels is not None:
            head_to_head_specs: list[dict[str, Any]] = []
            for label in head_to_head_labels:
                entry = entries_by_label[label]
                head_to_head_specs.append(
                    {
                        "candidate_label": entry.label,
                        "candidate_matrix": value_matrices[(entry.r1, entry.r2)],
                        "candidate_error_matrix": candidate_error_matrices[(entry.r1, entry.r2)],
                        "r1": entry.r1,
                        "r2": entry.r2,
                        "score": entry.score,
                    }
                )
            head_to_head_output = (
                os.path.abspath(args.head_to_head_output)
                if args.head_to_head_output
                else os.path.join(
                    output_dir,
                    f"scale{scale:02d}_head_to_head_{head_to_head_labels[0]}_vs_{head_to_head_labels[1]}_{args.field}.png",
                )
            )
            _make_head_to_head_comparison(
                reference_matrix=reference_matrix,
                reference_error_matrix=reference_error_matrix,
                candidate_specs=head_to_head_specs,
                m_values=m_values,
                n_values=n_values,
                scale=scale,
                family_size=family_size,
                output_path=head_to_head_output,
                field_label=field_label,
            )

        if all_candidate_comparisons_dir is not None:
            scale_dir = os.path.join(all_candidate_comparisons_dir, f"scale{scale:02d}")
            os.makedirs(scale_dir, exist_ok=True)
            for entry in entries:
                candidate_output_path = os.path.join(
                    scale_dir,
                    f"{entry.label}_comparison_{args.field}.png",
                )
                _make_best_candidate_comparison(
                    reference_matrix=reference_matrix,
                    reference_error_matrix=reference_error_matrix,
                    candidate_matrix=value_matrices[(entry.r1, entry.r2)],
                    candidate_error_matrix=candidate_error_matrices[(entry.r1, entry.r2)],
                    m_values=m_values,
                    n_values=n_values,
                    scale=scale,
                    family_size=family_size,
                    candidate_label=entry.label,
                    r1=entry.r1,
                    r2=entry.r2,
                    score=entry.score,
                    output_path=candidate_output_path,
                    field_label=field_label,
                    title_prefix="Untwisted candidate",
                )
        scale_outputs.append(
            {
                "scale": scale,
                "family_size": family_size,
                "reference_plot": reference_plot,
                "gallery_plot": gallery_plot,
                "difference_plot": difference_plot,
                "best_candidate_comparison_plot": best_candidate_comparison_plot,
            }
        )

    summary_path = os.path.join(output_dir, "README.md")
    _write_summary(
        output_path=summary_path,
        reference_manifest_path=reference_manifest_path,
        scan_json_path=scan_json_path,
        scale_outputs=scale_outputs,
        all_candidate_comparisons_dir=all_candidate_comparisons_dir,
    )

    print(f"wrote {summary_path}")
    for item in scale_outputs:
        print(f"wrote {item['reference_plot']}")
        print(f"wrote {item['gallery_plot']}")
        print(f"wrote {item['difference_plot']}")
        print(f"wrote {item['best_candidate_comparison_plot']}")
    if all_candidate_comparisons_dir is not None:
        print(f"wrote {all_candidate_comparisons_dir}")


if __name__ == "__main__":
    main()