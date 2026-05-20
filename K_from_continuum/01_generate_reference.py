#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from workflow_common import (
    build_payload_summary,
    build_test_geometry_map,
    check_required_sections,
    copy_file_atomic,
    couplings_from_ratio,
    evaluate_observable_fit,
    ensure_dir,
    ensure_simulator,
    find_beta_for_lattice,
    fit_beta_continuum_free_power,
    fit_observable_continuum_power,
    load_json,
    load_payload_from_dat,
    log,
    metadata_path_for_data,
    resolve_path,
    run_production_payload,
    sample_directional_channels,
    save_json,
    timestamp,
    write_dat,
)

_HERE = os.path.dirname(os.path.abspath(__file__))

DEFAULT_CONFIG_PATH = "configs/reference_example.json"


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


def _reference_payload_file_name(L: int, r1: float, r2: float) -> str:
    return f"reference_L{int(L)}_r1{float(r1):.6f}_r2{float(r2):.6f}.dat"


def _write_payload_outputs(summary: dict[str, Any], data_path: str) -> str:
    metadata_path = metadata_path_for_data(data_path)
    legacy_pickle_path = os.path.splitext(data_path)[0] + ".pkl"
    copy_file_atomic(summary["all_to_all_file"], data_path)
    save_json(metadata_path, summary)
    if os.path.exists(legacy_pickle_path):
        os.remove(legacy_pickle_path)
    return metadata_path


def _load_existing_metadata(data_path: str) -> dict[str, Any] | None:
    metadata_path = metadata_path_for_data(data_path)
    if not (os.path.exists(data_path) and os.path.exists(metadata_path)):
        return None
    return load_json(metadata_path)


def _existing_finite_beta_summary(
    metadata: dict[str, Any],
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta_cfg: dict[str, Any],
    label: str,
) -> dict[str, Any]:
    beta_value = metadata.get("finite_size_beta_c", metadata.get("beta_c"))
    if beta_value is None:
        raise ValueError("existing metadata is missing finite-size beta information")
    beta_sigma = metadata.get(
        "finite_size_beta_c_sigma",
        metadata.get("beta_c_sigma", 0.0),
    )
    Lx, Ly, Tx, Ty = lattice
    return {
        "label": str(metadata.get("label", label)),
        "created_at": str(metadata.get("created_at", timestamp())),
        "beta_find_wall_seconds": float(metadata.get("beta_find_wall_seconds", 0.0)),
        "L": int(metadata.get("L", max(abs(Lx), abs(Ly)))),
        "Lx": int(metadata.get("Lx", Lx)),
        "Ly": int(metadata.get("Ly", Ly)),
        "Tx": int(metadata.get("Tx", Tx)),
        "Ty": int(metadata.get("Ty", Ty)),
        "k1": float(metadata.get("k1", couplings["k1"])),
        "k2": float(metadata.get("k2", couplings["k2"])),
        "k3": float(metadata.get("k3", couplings["k3"])),
        "r1": float(metadata.get("r1", couplings["r1"])),
        "r2": float(metadata.get("r2", couplings["r2"])),
        "beta_c": float(beta_value),
        "beta_c_sigma": float(beta_sigma),
        "chi_peak": float(metadata.get("chi_peak", float("nan"))),
        "scan_betas": [float(x) for x in metadata.get("scan_betas", [])],
        "scan_chis": [float(x) for x in metadata.get("scan_chis", [])],
        "scan_chi_errs": [float(x) for x in metadata.get("scan_chi_errs", [])],
        "beta_finder": dict(metadata.get("beta_finder") or beta_cfg),
    }


def _can_reuse_continuum_payload(metadata: dict[str, Any], continuum_beta: float, tol: float = 1e-12) -> bool:
    source = str(metadata.get("production_beta_source", metadata.get("beta_source", "")))
    beta_raw = metadata.get("production_beta", metadata.get("beta_c"))
    if beta_raw is None:
        return False
    try:
        beta_val = float(beta_raw)
    except (TypeError, ValueError):
        return False
    return source == "free_power_continuum" and abs(beta_val - float(continuum_beta)) <= tol


def _save_reference_fss_plot(
    *,
    reference_dir: str,
    run_tag: str,
    k_values: list[int],
    channels: list[dict[str, Any]],
    dpi: int,
) -> str:
    n_rows = 3
    n_cols = max(len(k_values), 1)
    fig_w = max(4.0 * n_cols, 10.0)
    fig_h = max(3.3 * n_rows, 9.5)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    for channel in channels:
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
                ax.errorbar(x_sorted, y_sorted, yerr=yerr, fmt="o", color="#1f77b4", label="reference sizes")
            else:
                ax.plot(x_sorted, y_sorted, "o", color="#1f77b4", label="reference sizes")

        A = float(channel["A"])
        sigma_A = float(channel["sigma_A"])
        B = float(channel["B"])
        C = float(channel["C"])
        fit_mode = str(channel.get("fit_mode", "power_fit"))
        if np.isfinite(A):
            if np.isfinite(sigma_A):
                ax.errorbar([0.0], [A], yerr=[sigma_A], fmt="*", markersize=12, color="#2ca02c", label="continuum A")
            else:
                ax.plot([0.0], [A], "*", markersize=12, color="#2ca02c", label="continuum A")

        if np.isfinite(A) and np.isfinite(B) and np.isfinite(C):
            x_max = max([0.0] + (1.0 / L_vals).tolist()) if L_vals.size > 0 else 1.0
            x_fit = np.linspace(0.0, max(x_max, 1e-3), 200)
            y_fit = evaluate_observable_fit(x_fit, A, B, C, fit_mode)
            ax.plot(x_fit, y_fit, "-", color="#2ca02c", alpha=0.9, label="fit")

        ax.set_title(f"cycle={cycle}  k={channel['k']}  t={channel['t']:.3f}")
        if cycle == n_rows - 1:
            ax.set_xlabel("1/L")
        if ik == 0:
            ax.set_ylabel("G")
        ax.grid(alpha=0.25)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        fig.legend(list(unique.values()), list(unique.keys()), loc="upper center", ncol=min(4, len(unique)))

    fig.suptitle(f"Reference family continuum fits  (tag={run_tag})")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])
    out_png = os.path.join(reference_dir, "reference_fss.png")
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 1: generate a reference size family, fit a shared continuum beta, "
            "rerun production at that beta, and fit the reference channels to the continuum."
        )
    )
    parser.add_argument(
        "--config",
        type=str,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the reference config JSON (defaults to DEFAULT_CONFIG_PATH)",
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

    check_required_sections(
        cfg,
        ["run", "paths", "execution", "reference_family", "reference_couplings", "beta_finder", "mc", "analysis"],
    )
    if args.tag is not None:
        cfg["run"]["tag"] = args.tag

    _check_required_keys(cfg["run"], "run", ["tag"])
    _check_required_keys(cfg["paths"], "paths", ["results_root", "resume"])
    _check_required_keys(cfg["execution"], "execution", ["n_workers"])
    _check_required_keys(
        cfg["beta_finder"],
        "beta_finder",
        [
            "beta_lo",
            "beta_hi",
            "n_coarse",
            "n_refine",
            "n_refine2",
            "n_refine3",
            "n_traj_coarse",
            "n_traj_fine",
            "max_shifts",
            "jackknife",
        ],
    )
    _check_required_keys(cfg["mc"], "mc", ["n_traj", "n_skip", "n_therm"])
    _check_required_keys(cfg["analysis"], "analysis", ["k_values", "k_denominator", "power_fit"])
    _check_required_keys(cfg["analysis"]["power_fit"], "analysis.power_fit", ["c_min", "c_max", "c_initial"])

    beta_strategy = dict(cfg.get("beta_strategy") or {})
    beta_strategy_mode = str(beta_strategy.get("mode", "free_power_continuum"))
    if beta_strategy_mode != "free_power_continuum":
        raise ValueError("K_from_continuum Step 1 requires beta_strategy.mode = 'free_power_continuum'")

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _HERE)
    run_root = os.path.join(results_root, tag)
    reference_dir = os.path.join(run_root, "reference_data")
    grid_root = os.path.join(reference_dir, "grid")
    scratch_root = os.path.join(run_root, "_mc_scratch", "reference_data")
    ensure_dir(reference_dir)
    ensure_dir(grid_root)
    ensure_dir(scratch_root)

    geometry_map = build_test_geometry_map(cfg["reference_family"])
    sizes = [int(v) for v in cfg["reference_family"]["sizes"]]
    if len(sizes) < 3:
        raise ValueError("reference_family.sizes must contain at least 3 sizes")
    couplings = couplings_from_ratio(cfg["reference_couplings"], "reference_couplings")
    resume = bool(cfg["paths"]["resume"])
    workers = int(cfg["execution"]["n_workers"])
    if workers < 1:
        raise ValueError("execution.n_workers must be >= 1")
    if workers != 1:
        log("Step 1 currently runs serially for one reference point; execution.n_workers > 1 is ignored")

    analysis_cfg = cfg["analysis"]
    k_values = [int(v) for v in analysis_cfg["k_values"]]
    k_denom = int(analysis_cfg["k_denominator"])
    if k_denom <= 0:
        raise ValueError("analysis.k_denominator must be positive")
    fractions = np.array(k_values, dtype=float) / float(k_denom)
    c_min = float(analysis_cfg["power_fit"]["c_min"])
    c_max = float(analysis_cfg["power_fit"]["c_max"])
    c_initial = float(analysis_cfg["power_fit"]["c_initial"])
    fit_method = str(analysis_cfg["power_fit"].get("method", "power")).strip().lower()
    min_sizes_for_free_C = int(analysis_cfg["power_fit"].get("min_sizes_for_free_C", 8))
    if c_max <= c_min:
        raise ValueError("analysis.power_fit requires c_max > c_min")
    if min_sizes_for_free_C < 3:
        raise ValueError("analysis.power_fit.min_sizes_for_free_C must be >= 3")
    if fit_method not in {"power", "taylor2"}:
        raise ValueError("analysis.power_fit.method must be 'power' or 'taylor2'")
    fss_cfg = analysis_cfg.get("fss_plots") or {}
    fss_enabled = bool(fss_cfg.get("enabled", True))
    fss_dpi = int(fss_cfg.get("dpi", 150))

    exe = ensure_simulator(cfg["execution"])

    size_jobs: list[dict[str, Any]] = []
    finite_by_L: dict[int, dict[str, Any]] = {}
    existing_metadata: dict[int, dict[str, Any]] = {}
    for L in sizes:
        lattice = tuple(geometry_map[int(L)])
        out_dir = os.path.join(grid_root, f"L{int(L)}")
        ensure_dir(out_dir)
        data_path = os.path.join(out_dir, _reference_payload_file_name(L, couplings["r1"], couplings["r2"]))
        label = f"reference_L{int(L)}_r1{couplings['r1']:.6f}_r2{couplings['r2']:.6f}"
        size_jobs.append(
            {
                "L": int(L),
                "lattice": lattice,
                "data_path": data_path,
                "label": label,
            }
        )
        metadata = _load_existing_metadata(data_path) if resume else None
        if metadata is not None:
            existing_metadata[int(L)] = metadata
            try:
                finite_by_L[int(L)] = _existing_finite_beta_summary(
                    metadata,
                    lattice,
                    couplings,
                    cfg["beta_finder"],
                    label,
                )
                continue
            except Exception:
                pass
        finite_by_L[int(L)] = find_beta_for_lattice(
            exe=exe,
            lattice=lattice,
            couplings=couplings,
            beta_cfg=cfg["beta_finder"],
            scratch_root=scratch_root,
            label=label,
        )

    beta_fit = fit_beta_continuum_free_power(
        sizes,
        [float(finite_by_L[L]["beta_c"]) for L in sizes],
        [float(finite_by_L[L].get("beta_c_sigma", 0.0)) for L in sizes],
        weighted=bool(beta_strategy.get("weighted", False)),
        exponent_initial=float(beta_strategy.get("exponent_initial", 1.0)),
        exponent_min=float(beta_strategy.get("exponent_min", 0.05)),
        exponent_max=float(beta_strategy.get("exponent_max", 3.0)),
    )
    continuum_beta = float(beta_fit["beta_c_continuum"])
    continuum_beta_sigma = beta_fit.get("beta_c_continuum_sigma")

    job_results: list[dict[str, Any]] = []
    for size_job in size_jobs:
        L = int(size_job["L"])
        data_path = size_job["data_path"]
        metadata_path = metadata_path_for_data(data_path)
        metadata = existing_metadata.get(L)
        if resume and metadata is not None and _can_reuse_continuum_payload(metadata, continuum_beta):
            job_results.append(
                {
                    "status": "skip",
                    "L": L,
                    "r1": float(couplings["r1"]),
                    "r2": float(couplings["r2"]),
                    "data_path": data_path,
                    "metadata_path": metadata_path,
                    "beta_c": continuum_beta,
                    "finite_size_beta_c": float(finite_by_L[L]["beta_c"]),
                    "message": "reused existing payload at continuum beta",
                }
            )
            continue

        production_summary = run_production_payload(
            exe=exe,
            lattice=tuple(size_job["lattice"]),
            couplings=couplings,
            beta=continuum_beta,
            mc_cfg=cfg["mc"],
            scratch_root=scratch_root,
            label=size_job["label"],
        )
        summary = build_payload_summary(
            label=size_job["label"],
            lattice=tuple(size_job["lattice"]),
            couplings=couplings,
            finite_beta_summary=finite_by_L[L],
            production_summary=production_summary,
            beta_source="free_power_continuum",
            production_beta_sigma=(None if continuum_beta_sigma is None else float(continuum_beta_sigma)),
            beta_extrapolation=beta_fit,
        )
        metadata_path = _write_payload_outputs(summary, data_path)
        job_results.append(
            {
                "status": "ok",
                "L": L,
                "r1": float(couplings["r1"]),
                "r2": float(couplings["r2"]),
                "data_path": data_path,
                "metadata_path": metadata_path,
                "beta_c": float(summary["beta_c"]),
                "finite_size_beta_c": float(summary["finite_size_beta_c"]),
                "wall_seconds": float(summary["wall_seconds"]),
                "message": "completed with shared continuum beta",
            }
        )

    payload_cache: dict[tuple[str, str], dict[str, Any]] = {}
    channel_cache: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    log("[analysis] loading reference payloads and sampling directional channels")
    for size_job in size_jobs:
        metadata_path = metadata_path_for_data(size_job["data_path"])
        cache_key = (size_job["data_path"], metadata_path)
        if cache_key not in payload_cache:
            payload_cache[cache_key] = load_payload_from_dat(size_job["data_path"], metadata_path)
        if cache_key not in channel_cache:
            channel_cache[cache_key] = sample_directional_channels(payload_cache[cache_key], fractions)

    log("[analysis] fitting reference continuum channels")
    raw_rows: list[list[Any]] = []
    continuum_rows: list[list[Any]] = []
    fss_rows: list[list[Any]] = []
    channel_plots: list[dict[str, Any]] = []
    for cycle in range(3):
        for ik, kval in enumerate(k_values):
            L_vals: list[float] = []
            y_vals: list[float] = []
            s_vals: list[float] = []
            for size_job in size_jobs:
                L = int(size_job["L"])
                metadata_path = metadata_path_for_data(size_job["data_path"])
                cache_key = (size_job["data_path"], metadata_path)
                payload = payload_cache[cache_key]
                G_ref, sG_ref = channel_cache[cache_key]
                y_val = float(G_ref[cycle, ik])
                s_val = float(sG_ref[cycle, ik])
                L_vals.append(float(L))
                y_vals.append(y_val)
                s_vals.append(s_val)
                raw_rows.append(
                    [
                        int(L),
                        int(payload["Lx"]),
                        int(payload["Ly"]),
                        int(payload["Tx"]),
                        int(payload["Ty"]),
                        float(payload["beta_c"]),
                        cycle,
                        int(kval),
                        float(fractions[ik]),
                        y_val,
                        s_val,
                    ]
                )

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
                    cycle,
                    int(kval),
                    float(fractions[ik]),
                    float(A),
                    float(sigma_A),
                    float(B),
                    float(C),
                    int(n_used),
                    fit_mode,
                ]
            )
            for L_i, y_i, s_i in zip(L_vals, y_vals, s_vals):
                invL_i = (1.0 / float(L_i)) if float(L_i) > 0.0 else np.nan
                fss_rows.append(
                    [
                        cycle,
                        int(kval),
                        float(fractions[ik]),
                        "reference_size",
                        float(L_i),
                        invL_i,
                        float(y_i),
                        float(s_i),
                        "",
                        np.nan,
                        np.nan,
                    ]
                )
            fss_rows.append(
                [
                    cycle,
                    int(kval),
                    float(fractions[ik]),
                    "reference_continuum_fit",
                    0.0,
                    0.0,
                    float(A),
                    float(sigma_A),
                    fit_mode,
                    float(B),
                    float(C),
                ]
            )
            channel_plots.append(
                {
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
                }
            )

    fit_model_desc = (
        "A + B*(1/L) + C*(1/L)^2"
        if fit_method == "taylor2"
        else f"A + B*(1/L)^C with C in [{c_min}, {c_max}], min_sizes_for_free_C={min_sizes_for_free_C}"
    )

    raw_channels_path = os.path.join(reference_dir, "reference_raw_channels.dat")
    continuum_channels_path = os.path.join(reference_dir, "reference_continuum_channels.dat")
    fss_data_path = os.path.join(reference_dir, "reference_fss_data.dat")
    beta_path = os.path.join(reference_dir, "continuum_beta_extrapolation.json")
    manifest_path = os.path.join(reference_dir, "manifest_reference.json")
    fss_plot_path = None
    if fss_enabled:
        fss_plot_path = _save_reference_fss_plot(
            reference_dir=reference_dir,
            run_tag=tag,
            k_values=k_values,
            channels=channel_plots,
            dpi=fss_dpi,
        )

    header_common = [
        "Continuum reference workflow",
        f"run_tag={tag}",
        f"sizes={sizes}",
        f"couplings=(k1={couplings['k1']}, k2={couplings['k2']}, k3={couplings['k3']})",
        f"ratios=(k1_over_k3={couplings['r1']}, k2_over_k3={couplings['r2']})",
        f"beta_continuum={continuum_beta}",
        f"k_values={k_values}",
        f"k_denominator={k_denom}",
        f"continuum_fit={fit_model_desc}",
    ]

    write_dat(
        raw_channels_path,
        header_common,
        ["L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_rows,
    )
    write_dat(
        continuum_channels_path,
        header_common,
        ["cycle", "k", "t", "A", "sigma_A", "B", "C", "n_sizes_used", "fit_mode"],
        continuum_rows,
    )
    write_dat(
        fss_data_path,
        header_common,
        ["cycle", "k", "t", "source", "L", "invL", "G", "sigma_G", "fit_mode", "B", "C"],
        fss_rows,
    )
    save_json(
        beta_path,
        {
            "created_at": timestamp(),
            "config_path": config_path,
            "run_tag": tag,
            "mode": beta_strategy_mode,
            "point": {
                "r1": float(couplings["r1"]),
                "r2": float(couplings["r2"]),
                "beta_extrapolation": beta_fit,
                "finite_beta_rows": [
                    {
                        "L": int(L),
                        "beta_c": float(finite_by_L[L]["beta_c"]),
                        "beta_c_sigma": float(finite_by_L[L].get("beta_c_sigma", 0.0)),
                    }
                    for L in sizes
                ],
            },
        },
    )

    ok_count = sum(1 for row in job_results if row["status"] == "ok")
    skip_count = sum(1 for row in job_results if row["status"] == "skip")
    err_count = sum(1 for row in job_results if row["status"] == "err")
    manifest = {
        "created_at": timestamp(),
        "config_path": config_path,
        "run_tag": tag,
        "run_root": run_root,
        "reference_data_dir": reference_dir,
        "grid_root": grid_root,
        "reference_raw_channels": raw_channels_path,
        "reference_continuum_channels": continuum_channels_path,
        "reference_fss_data": fss_data_path,
        "reference_fss_plot": fss_plot_path,
        "continuum_beta_extrapolation": beta_path,
        "k_values": k_values,
        "k_denominator": k_denom,
        "power_fit": {
            "method": fit_method,
            "model": fit_model_desc,
            "c_min": c_min,
            "c_max": c_max,
            "c_initial": c_initial,
            "min_sizes_for_free_C": min_sizes_for_free_C,
        },
        "reference_summary": {
            "k1": float(couplings["k1"]),
            "k2": float(couplings["k2"]),
            "k3": float(couplings["k3"]),
            "r1": float(couplings["r1"]),
            "r2": float(couplings["r2"]),
            "beta_c_continuum": float(continuum_beta),
            "beta_c_continuum_sigma": None if continuum_beta_sigma is None else float(continuum_beta_sigma),
            "n_sizes": int(len(sizes)),
            "lmin": int(min(sizes)),
            "lmax": int(max(sizes)),
        },
        "summary": {
            "total_jobs": len(job_results),
            "ok": ok_count,
            "skip": skip_count,
            "err": err_count,
        },
        "jobs": sorted(job_results, key=lambda row: row["L"]),
        "config": cfg,
    }
    save_json(manifest_path, manifest)
    log(f"Reference manifest written: {manifest_path}")

    if err_count > 0:
        raise RuntimeError(f"{err_count} reference jobs failed. Check manifest: {manifest_path}")


if __name__ == "__main__":
    main()