#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from workflow_common import (
    build_payload_summary,
    build_test_geometry_map,
    check_required_sections,
    copy_file_atomic,
    ensure_dir,
    ensure_simulator,
    find_beta_for_lattice,
    fit_beta_continuum_free_power,
    fit_beta_continuum_taylor2,
    load_json,
    log,
    metadata_path_for_data,
    parse_ratio_list,
    payload_file_name,
    resolve_path,
    run_one_payload,
    run_production_payload,
    save_json,
    timestamp,
)

_HERE = os.path.dirname(os.path.abspath(__file__))

# Edit this when running directly from an IDE without CLI args.
DEFAULT_CONFIG_PATH = "configs/tests_example.json"


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


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


def _can_reuse_continuum_payload(
    metadata: dict[str, Any],
    continuum_beta: float,
    beta_source: str,
    tol: float = 1e-12,
) -> bool:
    source = str(metadata.get("production_beta_source", metadata.get("beta_source", "")))
    beta_raw = metadata.get("production_beta", metadata.get("beta_c"))
    if beta_raw is None:
        return False
    try:
        beta_val = float(beta_raw)
    except (TypeError, ValueError):
        return False
    return source == str(beta_source) and abs(beta_val - float(continuum_beta)) <= tol


def _job_worker(job: dict[str, Any]) -> dict[str, Any]:
    data_path = job["data_path"]
    metadata_path = metadata_path_for_data(data_path)
    legacy_pickle_path = os.path.splitext(data_path)[0] + ".pkl"
    if job["resume"] and os.path.exists(data_path) and os.path.exists(metadata_path):
        if os.path.exists(legacy_pickle_path):
            os.remove(legacy_pickle_path)
        return {
            "status": "skip",
            "L": int(job["L"]),
            "r1": float(job["r1"]),
            "r2": float(job["r2"]),
            "data_path": data_path,
            "metadata_path": metadata_path,
            "message": "reused existing payload",
        }

    couplings = {
        "k3": float(job["k3"]),
        "r1": float(job["r1"]),
        "r2": float(job["r2"]),
        "k1": float(job["r1"]) * float(job["k3"]),
        "k2": float(job["r2"]) * float(job["k3"]),
    }
    try:
        summary = run_one_payload(
            exe=job["exe"],
            lattice=tuple(job["lattice"]),
            couplings=couplings,
            beta_cfg=job["beta_finder"],
            mc_cfg=job["mc"],
            scratch_root=job["scratch_root"],
            label=job["label"],
        )
        metadata_path = _write_payload_outputs(summary, data_path)
        return {
            "status": "ok",
            "L": int(job["L"]),
            "r1": float(job["r1"]),
            "r2": float(job["r2"]),
            "data_path": data_path,
            "metadata_path": metadata_path,
            "beta_c": float(summary["beta_c"]),
            "finite_size_beta_c": float(summary["finite_size_beta_c"]),
            "wall_seconds": float(summary["wall_seconds"]),
            "message": "completed",
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "status": "err",
            "L": int(job["L"]),
            "r1": float(job["r1"]),
            "r2": float(job["r2"]),
            "data_path": data_path,
            "metadata_path": metadata_path,
            "message": str(exc),
        }


def _point_worker(job: dict[str, Any]) -> dict[str, Any]:
    beta_strategy_mode = str(job["beta_strategy"].get("mode", "free_power_continuum"))
    couplings = {
        "k3": float(job["k3"]),
        "r1": float(job["r1"]),
        "r2": float(job["r2"]),
        "k1": float(job["r1"]) * float(job["k3"]),
        "k2": float(job["r2"]) * float(job["k3"]),
    }
    scheduled_size_jobs = sorted(job["size_jobs"], key=lambda size_job: int(size_job["L"]), reverse=True)
    finite_by_L: dict[int, dict[str, Any]] = {}
    existing_metadata: dict[int, dict[str, Any]] = {}
    errors: list[tuple[dict[str, Any], str]] = []

    for size_job in scheduled_size_jobs:
        L = int(size_job["L"])
        metadata = _load_existing_metadata(size_job["data_path"]) if job["resume"] else None
        if metadata is not None:
            existing_metadata[L] = metadata
            try:
                finite_by_L[L] = _existing_finite_beta_summary(
                    metadata,
                    tuple(size_job["lattice"]),
                    couplings,
                    job["beta_finder"],
                    size_job["label"],
                )
                continue
            except Exception:
                pass
        try:
            finite_by_L[L] = find_beta_for_lattice(
                exe=job["exe"],
                lattice=tuple(size_job["lattice"]),
                couplings=couplings,
                beta_cfg=job["beta_finder"],
                scratch_root=job["scratch_root"],
                label=size_job["label"],
            )
        except Exception as exc:  # noqa: BLE001
            errors.append((size_job, str(exc)))

    if errors:
        error_jobs = []
        for size_job in job["size_jobs"]:
            matched = next((msg for failed_job, msg in errors if int(failed_job["L"]) == int(size_job["L"])), None)
            message = matched or "continuum beta unavailable because at least one finite-size beta fit failed"
            error_jobs.append(
                {
                    "status": "err",
                    "L": int(size_job["L"]),
                    "r1": float(job["r1"]),
                    "r2": float(job["r2"]),
                    "data_path": size_job["data_path"],
                    "metadata_path": metadata_path_for_data(size_job["data_path"]),
                    "message": message,
                }
            )
        return {
            "r1": float(job["r1"]),
            "r2": float(job["r2"]),
            "jobs": error_jobs,
            "beta_extrapolation": None,
            "finite_beta_rows": [],
        }

    ordered_sizes = [int(size_job["L"]) for size_job in job["size_jobs"]]
    if beta_strategy_mode == "free_power_continuum":
        beta_fit = fit_beta_continuum_free_power(
            ordered_sizes,
            [float(finite_by_L[L]["beta_c"]) for L in ordered_sizes],
            [float(finite_by_L[L].get("beta_c_sigma", 0.0)) for L in ordered_sizes],
            weighted=bool(job["beta_strategy"].get("weighted", False)),
            exponent_initial=float(job["beta_strategy"].get("exponent_initial", 1.0)),
            exponent_min=float(job["beta_strategy"].get("exponent_min", 0.05)),
            exponent_max=float(job["beta_strategy"].get("exponent_max", 6.0)),
        )
    elif beta_strategy_mode == "taylor2_continuum":
        beta_fit = fit_beta_continuum_taylor2(
            ordered_sizes,
            [float(finite_by_L[L]["beta_c"]) for L in ordered_sizes],
            [float(finite_by_L[L].get("beta_c_sigma", 0.0)) for L in ordered_sizes],
            weighted=bool(job["beta_strategy"].get("weighted", False)),
        )
    else:
        raise ValueError(f"unsupported beta_strategy.mode: {beta_strategy_mode}")

    beta_source = str(beta_fit["method"])
    continuum_beta = float(beta_fit["beta_c_continuum"])
    continuum_beta_sigma = beta_fit.get("beta_c_continuum_sigma")

    job_results: list[dict[str, Any]] = []
    for size_job in scheduled_size_jobs:
        L = int(size_job["L"])
        data_path = size_job["data_path"]
        metadata_path = metadata_path_for_data(data_path)
        legacy_pickle_path = os.path.splitext(data_path)[0] + ".pkl"
        metadata = existing_metadata.get(L)
        if (
            job["resume"]
            and metadata is not None
            and _can_reuse_continuum_payload(metadata, continuum_beta, beta_source)
        ):
            if os.path.exists(legacy_pickle_path):
                os.remove(legacy_pickle_path)
            job_results.append(
                {
                    "status": "skip",
                    "L": L,
                    "r1": float(job["r1"]),
                    "r2": float(job["r2"]),
                    "data_path": data_path,
                    "metadata_path": metadata_path,
                    "beta_c": continuum_beta,
                    "finite_size_beta_c": float(finite_by_L[L]["beta_c"]),
                    "message": "reused existing payload at continuum beta",
                }
            )
            continue
        try:
            production_summary = run_production_payload(
                exe=job["exe"],
                lattice=tuple(size_job["lattice"]),
                couplings=couplings,
                beta=continuum_beta,
                mc_cfg=job["mc"],
                scratch_root=job["scratch_root"],
                label=size_job["label"],
            )
            summary = build_payload_summary(
                label=size_job["label"],
                lattice=tuple(size_job["lattice"]),
                couplings=couplings,
                finite_beta_summary=finite_by_L[L],
                production_summary=production_summary,
                beta_source=beta_source,
                production_beta_sigma=(
                    None if continuum_beta_sigma is None else float(continuum_beta_sigma)
                ),
                beta_extrapolation=beta_fit,
            )
            metadata_path = _write_payload_outputs(summary, data_path)
            job_results.append(
                {
                    "status": "ok",
                    "L": L,
                    "r1": float(job["r1"]),
                    "r2": float(job["r2"]),
                    "data_path": data_path,
                    "metadata_path": metadata_path,
                    "beta_c": float(summary["beta_c"]),
                    "finite_size_beta_c": float(summary["finite_size_beta_c"]),
                    "wall_seconds": float(summary["wall_seconds"]),
                    "message": "completed with shared continuum beta",
                }
            )
        except Exception as exc:  # noqa: BLE001
            job_results.append(
                {
                    "status": "err",
                    "L": L,
                    "r1": float(job["r1"]),
                    "r2": float(job["r2"]),
                    "data_path": data_path,
                    "metadata_path": metadata_path,
                    "message": str(exc),
                }
            )

    finite_beta_rows = [
        {
            "L": int(L),
            "beta_c": float(finite_by_L[L]["beta_c"]),
            "beta_c_sigma": float(finite_by_L[L].get("beta_c_sigma", 0.0)),
        }
        for L in ordered_sizes
    ]
    return {
        "r1": float(job["r1"]),
        "r2": float(job["r2"]),
        "continuum_beta": float(continuum_beta),
        "continuum_beta_sigma": (
            None if continuum_beta_sigma is None else float(continuum_beta_sigma)
        ),
        "jobs": job_results,
        "production_jobs": job_results,
        "beta_extrapolation": beta_fit,
        "finite_beta_rows": finite_beta_rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 2: generate test payloads for all requested L values and "
            "(k1/k3, k2/k3) grid points."
        )
    )
    parser.add_argument(
        "--config",
        type=str,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the test-generation config JSON (defaults to DEFAULT_CONFIG_PATH)",
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
        ["run", "paths", "execution", "test_family", "couplings", "beta_finder", "mc"],
    )
    if args.tag is not None:
        cfg["run"]["tag"] = args.tag

    _check_required_keys(cfg["run"], "run", ["tag"])
    _check_required_keys(cfg["paths"], "paths", ["results_root", "resume"])
    _check_required_keys(cfg["execution"], "execution", ["n_workers"])
    _check_required_keys(cfg["couplings"], "couplings", ["k3", "k1_over_k3_values", "k2_over_k3_values"])
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

    beta_strategy = dict(cfg.get("beta_strategy") or {})
    beta_strategy_mode = str(beta_strategy.get("mode", "per_size_beta"))
    if beta_strategy_mode not in {"per_size_beta", "free_power_continuum", "taylor2_continuum"}:
        raise ValueError(
            "beta_strategy.mode must be one of {'per_size_beta', 'free_power_continuum', 'taylor2_continuum'}"
        )

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _HERE)
    run_root = os.path.join(results_root, tag)
    test_data_dir = os.path.join(run_root, "test_data")
    grid_root = os.path.join(test_data_dir, "grid")
    scratch_root = os.path.join(run_root, "_mc_scratch", "test_data")
    ensure_dir(test_data_dir)
    ensure_dir(grid_root)
    ensure_dir(scratch_root)

    manifest_path = os.path.join(test_data_dir, "manifest_tests.json")
    extrapolation_manifest_path = os.path.join(test_data_dir, "continuum_beta_extrapolations.json")
    resume = bool(cfg["paths"]["resume"])

    geometry_map = build_test_geometry_map(cfg["test_family"])
    sizes = [int(v) for v in cfg["test_family"]["sizes"]]
    r1_values = parse_ratio_list(cfg["couplings"]["k1_over_k3_values"], "couplings.k1_over_k3_values")
    r2_values = parse_ratio_list(cfg["couplings"]["k2_over_k3_values"], "couplings.k2_over_k3_values")
    k3 = float(cfg["couplings"]["k3"])
    if k3 == 0.0:
        raise ValueError("couplings.k3 must be non-zero")
    if beta_strategy_mode in {"free_power_continuum", "taylor2_continuum"} and len(sizes) < 3:
        raise ValueError(
            "beta_strategy.mode in {'free_power_continuum', 'taylor2_continuum'} requires at least 3 test sizes"
        )

    exe = ensure_simulator(cfg["execution"])
    workers = int(cfg["execution"]["n_workers"])
    if workers < 1:
        raise ValueError("execution.n_workers must be >= 1")

    results: list[dict[str, Any]] = []
    beta_extrapolations: list[dict[str, Any]] = []

    if beta_strategy_mode == "per_size_beta":
        jobs: list[dict[str, Any]] = []
        for L in sizes:
            Lx, Ly, Tx, Ty = geometry_map[int(L)]
            out_dir = os.path.join(grid_root, f"L{int(L)}")
            ensure_dir(out_dir)
            for r2 in r2_values:
                for r1 in r1_values:
                    filename = payload_file_name(L, r1, r2)
                    data_path = os.path.join(out_dir, filename)
                    jobs.append(
                        {
                            "exe": exe,
                            "L": int(L),
                            "lattice": [int(Lx), int(Ly), int(Tx), int(Ty)],
                            "r1": float(r1),
                            "r2": float(r2),
                            "k3": float(k3),
                            "beta_finder": cfg["beta_finder"],
                            "mc": cfg["mc"],
                            "resume": resume,
                            "scratch_root": scratch_root,
                            "data_path": data_path,
                            "label": f"test_L{int(L)}_r1{float(r1):.6f}_r2{float(r2):.6f}",
                        }
                    )

        log(f"Running {len(jobs)} test jobs with workers={workers}")
        if workers == 1:
            for idx, job in enumerate(jobs, start=1):
                log(
                    f"[start {idx}/{len(jobs)}] L={job['L']} r1={job['r1']:.6f} "
                    f"r2={job['r2']:.6f}"
                )
                result = _job_worker(job)
                results.append(result)
                log(
                    f"[{idx}/{len(jobs)}] L={result['L']} r1={result['r1']:.6f} "
                    f"r2={result['r2']:.6f} status={result['status']}"
                )
        else:
            done = 0
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_job_worker, job) for job in jobs]
                log(f"[dispatch] queued {len(jobs)}/{len(jobs)} test jobs")
                for future in as_completed(futures):
                    result = future.result()
                    done += 1
                    results.append(result)
                    log(
                        f"[{done}/{len(jobs)}] L={result['L']} r1={result['r1']:.6f} "
                        f"r2={result['r2']:.6f} status={result['status']}"
                    )
    else:
        point_jobs: list[dict[str, Any]] = []
        for r2 in r2_values:
            for r1 in r1_values:
                size_jobs = []
                for L in sizes:
                    Lx, Ly, Tx, Ty = geometry_map[int(L)]
                    out_dir = os.path.join(grid_root, f"L{int(L)}")
                    ensure_dir(out_dir)
                    filename = payload_file_name(L, r1, r2)
                    size_jobs.append(
                        {
                            "L": int(L),
                            "lattice": [int(Lx), int(Ly), int(Tx), int(Ty)],
                            "data_path": os.path.join(out_dir, filename),
                            "label": f"test_L{int(L)}_r1{float(r1):.6f}_r2{float(r2):.6f}",
                        }
                    )
                point_jobs.append(
                    {
                        "exe": exe,
                        "r1": float(r1),
                        "r2": float(r2),
                        "k3": float(k3),
                        "beta_finder": cfg["beta_finder"],
                        "beta_strategy": beta_strategy,
                        "mc": cfg["mc"],
                        "n_workers": workers,
                        "resume": resume,
                        "scratch_root": scratch_root,
                        "size_jobs": size_jobs,
                    }
                )

        log(
            f"Running {len(point_jobs)} coupling-point jobs serially from large to small "
            f"with workers_per_lattice={workers} using beta_strategy.mode={beta_strategy_mode}"
        )
        for idx, point_job in enumerate(point_jobs, start=1):
            log(
                f"[start {idx}/{len(point_jobs)}] r1={point_job['r1']:.6f} "
                f"r2={point_job['r2']:.6f}"
            )
            point_result = _point_worker(point_job)
            results.extend(point_result["jobs"])
            if point_result["beta_extrapolation"] is not None:
                beta_extrapolations.append(point_result)
            point_status = "ok" if all(job["status"] in {"ok", "skip"} for job in point_result["jobs"]) else "err"
            log(
                f"[{idx}/{len(point_jobs)}] r1={point_result['r1']:.6f} "
                f"r2={point_result['r2']:.6f} status={point_status}"
            )

        extrap_payload = {
            "created_at": timestamp(),
            "config_path": config_path,
            "run_tag": tag,
            "mode": beta_strategy_mode,
            "points": sorted(beta_extrapolations, key=lambda row: (row["r2"], row["r1"])),
        }
        save_json(extrapolation_manifest_path, extrap_payload)
        log(f"Continuum beta extrapolation summary written: {extrapolation_manifest_path}")

    ok_count = sum(1 for r in results if r["status"] == "ok")
    skip_count = sum(1 for r in results if r["status"] == "skip")
    err_count = sum(1 for r in results if r["status"] == "err")

    manifest = {
        "created_at": timestamp(),
        "config_path": config_path,
        "run_tag": tag,
        "run_root": run_root,
        "test_data_dir": test_data_dir,
        "grid_root": grid_root,
        "exe": exe,
        "sizes": sizes,
        "geometry_map": {str(k): list(v) for k, v in geometry_map.items()},
        "k3": k3,
        "k1_over_k3_values": r1_values,
        "k2_over_k3_values": r2_values,
        "beta_strategy": {
            "mode": beta_strategy_mode,
            **beta_strategy,
        },
        "summary": {
            "total_jobs": len(results),
            "ok": ok_count,
            "skip": skip_count,
            "err": err_count,
        },
        "jobs": sorted(results, key=lambda r: (r["L"], r["r2"], r["r1"])),
        "config": cfg,
    }
    if beta_strategy_mode != "per_size_beta":
        manifest["continuum_beta_extrapolations"] = extrapolation_manifest_path
    save_json(manifest_path, manifest)
    log(f"Test manifest written: {manifest_path}")

    if err_count > 0:
        raise RuntimeError(f"{err_count} test jobs failed. Check manifest: {manifest_path}")


if __name__ == "__main__":
    main()
