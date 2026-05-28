#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent
REPO_ROOT = RESPONSIBLE_ROOT.parent
KFC_ROOT = REPO_ROOT / "K_from_continuum"
if str(KFC_ROOT) not in sys.path:
    sys.path.insert(0, str(KFC_ROOT))

from workflow_common import (  # noqa: E402
    copy_file_atomic,
    ensure_dir,
    ensure_simulator,
    exact_triangular_ising_beta,
    save_json,
    timestamp,
)


MULTIPLIERS = (0.8, 0.9, 1.0, 1.1, 1.2)
SIZES = (8, 16, 24, 32, 48, 64)
DEFAULT_MC = {
    "n_traj": 20000,
    "n_skip": 10,
    "n_therm": 2000,
}
DEFAULT_EXECUTION = {
    "exe": None,
    "auto_build": True,
    "force_rebuild": False,
    "build_command": ["make"],
    "build_timeout_s": 600,
}
GEOMETRIES = {
    "iso111": {
        "description": "Equilateral 1-1-1 control. Twisted reference is a large untwisted torus with isotropic couplings.",
        "twisted_lattice": (144, 144, 0, 0),
        "twisted_label": "large isotropic control",
        "center_r1": 1.0,
        "center_r2": 1.0,
        "center_source": "exact equilateral solution",
    },
    "acute456": {
        "description": "Story acute 4-5-6 benchmark. Twisted reference uses the large [144,144,72,24] torus from the holdout figures.",
        "twisted_lattice": (144, 144, 72, 24),
        "twisted_label": "story acute 4-5-6 holdout lattice",
        "center_r1": 4.702782819756,
        "center_r2": 7.353910143333,
        "center_source": "story acute 4-5-6 realized untwisted match",
    },
}


def _volume(lattice: tuple[int, int, int, int]) -> int:
    lx, ly, tx, ty = lattice
    return abs(int(lx) * int(ly) + int(tx) * int(ty))


def _ratio_token(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def _multiplier_token(value: float) -> str:
    return f"{value:.3f}".replace("-", "m").replace(".", "p")


def _couplings(r1: float, r2: float) -> dict[str, float]:
    return {
        "k1": float(r1),
        "k2": float(r2),
        "k3": 1.0,
        "r1": float(r1),
        "r2": float(r2),
    }


def _build_job_plan(root: Path) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    scratch_root = root / "_mc_scratch"
    data_root = root / "data"

    for geometry_id, geometry in GEOMETRIES.items():
        twisted_lattice = tuple(int(v) for v in geometry["twisted_lattice"])
        twisted_couplings = _couplings(1.0, 1.0)
        twisted_beta = exact_triangular_ising_beta(twisted_couplings)
        twisted_dir = data_root / geometry_id / "twisted" / "reference" / (
            f"Lx{twisted_lattice[0]}_Ly{twisted_lattice[1]}_"
            f"Tx{twisted_lattice[2]}_Ty{twisted_lattice[3]}"
        )
        jobs.append(
            {
                "job_id": f"{geometry_id}__twisted__reference",
                "geometry_id": geometry_id,
                "geometry_description": geometry["description"],
                "method": "twisted",
                "candidate_label": "reference",
                "candidate_source": geometry["twisted_label"],
                "center_r1": float(geometry["center_r1"]),
                "center_r2": float(geometry["center_r2"]),
                "multiplier_r1": 1.0,
                "multiplier_r2": 1.0,
                "lattice": list(twisted_lattice),
                "volume": _volume(twisted_lattice),
                "couplings": twisted_couplings,
                "beta": float(twisted_beta),
                "critical_source": "exact triangular sinh-rule critical point",
                "data_path": str(twisted_dir / "two_point_all_to_all.dat"),
                "meta_path": str(twisted_dir / "two_point_all_to_all.meta.json"),
                "scratch_root": str(scratch_root / geometry_id / "twisted"),
            }
        )

        center_r1 = float(geometry["center_r1"])
        center_r2 = float(geometry["center_r2"])
        for mul_r1 in MULTIPLIERS:
            for mul_r2 in MULTIPLIERS:
                r1 = center_r1 * mul_r1
                r2 = center_r2 * mul_r2
                couplings = _couplings(r1, r2)
                beta = exact_triangular_ising_beta(couplings)
                candidate_dir = data_root / geometry_id / "untwisted" / (
                    f"r1_{_ratio_token(r1)}__r2_{_ratio_token(r2)}"
                )
                for size in SIZES:
                    lattice = (int(size), int(size), 0, 0)
                    lattice_dir = candidate_dir / f"Lx{size}_Ly{size}_Tx0_Ty0"
                    jobs.append(
                        {
                            "job_id": (
                                f"{geometry_id}__untwisted__"
                                f"r1x{_multiplier_token(mul_r1)}__r2x{_multiplier_token(mul_r2)}__L{size:03d}"
                            ),
                            "geometry_id": geometry_id,
                            "geometry_description": geometry["description"],
                            "method": "untwisted",
                            "candidate_label": (
                                f"r1x{mul_r1:.3f}_r2x{mul_r2:.3f}"
                            ),
                            "candidate_source": geometry["center_source"],
                            "center_r1": center_r1,
                            "center_r2": center_r2,
                            "multiplier_r1": float(mul_r1),
                            "multiplier_r2": float(mul_r2),
                            "lattice": list(lattice),
                            "volume": _volume(lattice),
                            "couplings": couplings,
                            "beta": float(beta),
                            "critical_source": "exact triangular sinh-rule critical point",
                            "data_path": str(lattice_dir / "two_point_all_to_all.dat"),
                            "meta_path": str(lattice_dir / "two_point_all_to_all.meta.json"),
                            "scratch_root": str(scratch_root / geometry_id / "untwisted"),
                        }
                    )
    return jobs


def _existing_result(job: dict[str, Any]) -> bool:
    return Path(job["data_path"]).is_file() and Path(job["meta_path"]).is_file()


def _collect_named_files(root_dir: str, target_name: str) -> list[str]:
    root_path = Path(root_dir)
    if not root_path.is_dir():
        return []
    return [str(path) for path in root_path.rglob(target_name)]


def _tail_lines(text: str, n_lines: int = 8) -> list[str]:
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) <= n_lines:
        return lines
    return lines[-n_lines:]


def _read_text(path: Path) -> str:
    if not path.is_file():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def _manual_wall_seconds(source_dir: Path, all_to_all_file: Path) -> float:
    start_path = source_dir / "start_epoch.txt"
    if not start_path.is_file():
        return 0.0
    raw_value = _read_text(start_path).strip()
    if not raw_value:
        return 0.0
    try:
        start_epoch = float(raw_value)
    except ValueError:
        return 0.0
    return max(0.0, float(all_to_all_file.stat().st_mtime) - start_epoch)


def _run_production_payload_direct(
    *,
    exe: str,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta: float,
    mc_cfg: dict[str, Any],
    scratch_root: str,
    label: str,
) -> dict[str, Any]:
    lx, ly, tx, ty = lattice
    k1 = float(couplings["k1"])
    k2 = float(couplings["k2"])
    k3 = float(couplings["k3"])
    prod_data_dir = Path(scratch_root) / label / "prod"
    ensure_dir(str(prod_data_dir))
    before_files = set(_collect_named_files(str(prod_data_dir), "two_point_all_to_all.dat"))
    seed = int(time.time() * 1000) & 0xFFFFFFFF
    cmd = [
        exe,
        "--L_x", str(lx),
        "--L_y", str(ly),
        "--T_x", str(tx),
        "--T_y", str(ty),
        "--k1", f"{k1:.12f}",
        "--k2", f"{k2:.12f}",
        "--k3", f"{k3:.12f}",
        "--beta", f"{float(beta):.12f}",
        "--n_traj", str(int(mc_cfg["n_traj"])),
        "--n_skip", str(int(mc_cfg["n_skip"])),
        "--n_therm", str(int(mc_cfg["n_therm"])),
        "--seed", str(seed),
        "--data_dir", str(prod_data_dir),
    ]
    start_t = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "simulator failed with non-zero exit code\n"
            f"stdout tail: {_tail_lines(result.stdout)}\n"
            f"stderr tail: {_tail_lines(result.stderr)}"
        )

    all_files = _collect_named_files(str(prod_data_dir), "two_point_all_to_all.dat")
    new_files = [path for path in all_files if path not in before_files]
    candidate_files = new_files if new_files else all_files
    if not candidate_files:
        raise RuntimeError(
            "simulator completed but no two_point_all_to_all.dat was found under "
            f"{prod_data_dir}"
        )
    all_to_all_file = max(candidate_files, key=lambda path: os.path.getmtime(path))
    return {
        "production_beta": float(beta),
        "production_wall_seconds": float(time.time() - start_t),
        "all_to_all_file": str(all_to_all_file),
        "stdout_tail": _tail_lines(result.stdout),
        "stderr_tail": _tail_lines(result.stderr),
        "mc": {
            "n_traj": int(mc_cfg["n_traj"]),
            "n_skip": int(mc_cfg["n_skip"]),
            "n_therm": int(mc_cfg["n_therm"]),
        },
    }


def _latest_logged_status(root: Path) -> dict[str, str]:
    log_path = root / "job_log.jsonl"
    latest: dict[str, str] = {}
    if not log_path.is_file():
        return latest
    with log_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            job_id = payload.get("job_id")
            status = payload.get("status")
            if isinstance(job_id, str) and isinstance(status, str):
                latest[job_id] = status
    return latest


def _write_plan(root: Path, all_jobs: list[dict[str, Any]], selected_jobs: list[dict[str, Any]], workers: int) -> None:
    payload = {
        "campaign": "responsible_method_tests/standard",
        "created_at": timestamp(),
        "root": str(root),
        "mc": DEFAULT_MC,
        "workers_requested": int(workers),
        "geometries": GEOMETRIES,
        "sizes": list(SIZES),
        "multipliers": list(MULTIPLIERS),
        "job_count": len(all_jobs),
        "jobs": all_jobs,
    }
    save_json(str(root / "campaign_plan.json"), payload)
    selection_payload = {
        "updated_at": timestamp(),
        "selected_job_count": len(selected_jobs),
        "selected_job_ids": [str(job["job_id"]) for job in selected_jobs],
    }
    save_json(str(root / "current_selection.json"), selection_payload)


def _append_log(root: Path, payload: dict[str, Any]) -> None:
    log_path = root / "job_log.jsonl"
    ensure_dir(str(log_path.parent))
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True))
        handle.write("\n")


def _write_summary(root: Path, jobs: list[dict[str, Any]]) -> None:
    latest_status = _latest_logged_status(root)
    completed = sum(1 for job in jobs if _existing_result(job))
    skipped = sum(1 for job in jobs if latest_status.get(str(job["job_id"])) == "skipped")
    failed = sum(
        1
        for job in jobs
        if latest_status.get(str(job["job_id"])) == "failed" and not _existing_result(job)
    )
    payload = {
        "updated_at": timestamp(),
        "job_count": len(jobs),
        "completed": int(completed),
        "skipped": int(skipped),
        "failed": int(failed),
        "pending": int(len(jobs) - completed - failed),
    }
    save_json(str(root / "status_summary.json"), payload)


def _import_manual_job(root: Path, jobs: list[dict[str, Any]], job_id: str, source_dir: Path) -> dict[str, Any]:
    job = next((candidate for candidate in jobs if str(candidate["job_id"]) == job_id), None)
    if job is None:
        raise ValueError(f"unknown job id: {job_id}")

    all_files = _collect_named_files(str(source_dir), "two_point_all_to_all.dat")
    if not all_files:
        raise FileNotFoundError(
            f"no two_point_all_to_all.dat found under {source_dir}"
        )
    all_to_all_file = Path(max(all_files, key=lambda path: os.path.getmtime(path)))

    data_path = Path(job["data_path"])
    meta_path = Path(job["meta_path"])
    ensure_dir(str(data_path.parent))
    copy_file_atomic(str(all_to_all_file), str(data_path))

    lattice = tuple(int(v) for v in job["lattice"])
    couplings = dict(job["couplings"])
    stdout_tail = _tail_lines(_read_text(source_dir / "stdout.log"))
    stderr_tail = _tail_lines(_read_text(source_dir / "stderr.log"))
    metadata = {
        "label": str(job["job_id"]),
        "created_at": timestamp(),
        "campaign": "responsible_method_tests/standard",
        "geometry_id": str(job["geometry_id"]),
        "geometry_description": str(job["geometry_description"]),
        "method": str(job["method"]),
        "candidate_label": str(job["candidate_label"]),
        "candidate_source": str(job["candidate_source"]),
        "Lx": int(lattice[0]),
        "Ly": int(lattice[1]),
        "Tx": int(lattice[2]),
        "Ty": int(lattice[3]),
        "volume": int(job["volume"]),
        "k1": float(couplings["k1"]),
        "k2": float(couplings["k2"]),
        "k3": float(couplings["k3"]),
        "r1": float(couplings["r1"]),
        "r2": float(couplings["r2"]),
        "center_r1": float(job["center_r1"]),
        "center_r2": float(job["center_r2"]),
        "multiplier_r1": float(job["multiplier_r1"]),
        "multiplier_r2": float(job["multiplier_r2"]),
        "beta_c": float(job["beta"]),
        "beta_c_sigma": 0.0,
        "production_beta": float(job["beta"]),
        "production_beta_source": "exact_triangular_sinh_rule",
        "critical_source": str(job["critical_source"]),
        "production_wall_seconds": _manual_wall_seconds(source_dir, all_to_all_file),
        "stdout_tail": stdout_tail,
        "stderr_tail": stderr_tail,
        "mc": dict(DEFAULT_MC),
        "all_to_all_file": str(data_path),
        "source_all_to_all_file": str(all_to_all_file),
    }
    save_json(str(meta_path), metadata)

    result = {
        "job_id": str(job["job_id"]),
        "status": "completed",
        "data_path": str(data_path),
        "meta_path": str(meta_path),
        "finished_at": timestamp(),
        "production_wall_seconds": float(metadata["production_wall_seconds"]),
    }
    _append_log(root, result)
    _write_summary(root, jobs)
    return result


def _run_one_job(job: dict[str, Any], exe: str) -> dict[str, Any]:
    try:
        data_path = Path(job["data_path"])
        meta_path = Path(job["meta_path"])
        ensure_dir(str(data_path.parent))

        if data_path.is_file() and meta_path.is_file():
            return {
                "job_id": job["job_id"],
                "status": "skipped",
                "data_path": str(data_path),
                "meta_path": str(meta_path),
                "finished_at": timestamp(),
            }

        lattice = tuple(int(v) for v in job["lattice"])
        couplings = dict(job["couplings"])
        summary = _run_production_payload_direct(
            exe=exe,
            lattice=lattice,
            couplings=couplings,
            beta=float(job["beta"]),
            mc_cfg=DEFAULT_MC,
            scratch_root=str(job["scratch_root"]),
            label=str(job["job_id"]),
        )
        copy_file_atomic(summary["all_to_all_file"], str(data_path))
        metadata = {
            "label": str(job["job_id"]),
            "created_at": timestamp(),
            "campaign": "responsible_method_tests/standard",
            "geometry_id": str(job["geometry_id"]),
            "geometry_description": str(job["geometry_description"]),
            "method": str(job["method"]),
            "candidate_label": str(job["candidate_label"]),
            "candidate_source": str(job["candidate_source"]),
            "Lx": int(lattice[0]),
            "Ly": int(lattice[1]),
            "Tx": int(lattice[2]),
            "Ty": int(lattice[3]),
            "volume": int(job["volume"]),
            "k1": float(couplings["k1"]),
            "k2": float(couplings["k2"]),
            "k3": float(couplings["k3"]),
            "r1": float(couplings["r1"]),
            "r2": float(couplings["r2"]),
            "center_r1": float(job["center_r1"]),
            "center_r2": float(job["center_r2"]),
            "multiplier_r1": float(job["multiplier_r1"]),
            "multiplier_r2": float(job["multiplier_r2"]),
            "beta_c": float(job["beta"]),
            "beta_c_sigma": 0.0,
            "production_beta": float(job["beta"]),
            "production_beta_source": "exact_triangular_sinh_rule",
            "critical_source": str(job["critical_source"]),
            "production_wall_seconds": float(summary.get("production_wall_seconds", 0.0)),
            "stdout_tail": list(summary.get("stdout_tail", [])),
            "stderr_tail": list(summary.get("stderr_tail", [])),
            "mc": dict(summary["mc"]),
            "all_to_all_file": str(data_path),
            "source_all_to_all_file": str(summary["all_to_all_file"]),
        }
        save_json(str(meta_path), metadata)
        return {
            "job_id": job["job_id"],
            "status": "completed",
            "data_path": str(data_path),
            "meta_path": str(meta_path),
            "finished_at": timestamp(),
            "production_wall_seconds": float(summary.get("production_wall_seconds", 0.0)),
        }
    except Exception as exc:
        return {
            "job_id": job["job_id"],
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "finished_at": timestamp(),
        }


def _filtered_jobs(
    jobs: list[dict[str, Any]],
    geometries: set[str] | None,
    methods: set[str] | None,
    max_jobs: int | None,
) -> list[dict[str, Any]]:
    selected = []
    for job in jobs:
        if geometries is not None and job["geometry_id"] not in geometries:
            continue
        if methods is not None and job["method"] not in methods:
            continue
        selected.append(job)
    if max_jobs is not None:
        return selected[: max(int(max_jobs), 0)]
    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the responsible_method_tests/standard raw two-point campaign directly from the source-built simulator."
    )
    parser.add_argument("--workers", type=int, default=2, help="Parallel workers for independent jobs")
    parser.add_argument(
        "--geometry",
        action="append",
        choices=sorted(GEOMETRIES.keys()),
        help="Restrict to one or more geometry ids",
    )
    parser.add_argument(
        "--method",
        action="append",
        choices=["twisted", "untwisted"],
        help="Restrict to one or more method families",
    )
    parser.add_argument("--max-jobs", type=int, default=None, help="Optional cap for debugging or chunked restarts")
    parser.add_argument("--dry-run", action="store_true", help="Write the plan and print counts without launching jobs")
    parser.add_argument(
        "--import-manual-job",
        type=str,
        default=None,
        help="Import a completed direct simulator output directory for the given job id",
    )
    parser.add_argument(
        "--manual-source-dir",
        type=str,
        default=None,
        help="Directory containing a completed manual simulator run and its logs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = HERE
    ensure_dir(str(root))
    all_jobs = _build_job_plan(root)

    if args.import_manual_job is not None:
        if not args.manual_source_dir:
            raise SystemExit("--manual-source-dir is required with --import-manual-job")
        result = _import_manual_job(
            root,
            all_jobs,
            str(args.import_manual_job),
            Path(args.manual_source_dir).resolve(),
        )
        print(json.dumps(result, sort_keys=True))
        return

    jobs = _filtered_jobs(
        all_jobs,
        set(args.geometry) if args.geometry else None,
        set(args.method) if args.method else None,
        args.max_jobs,
    )
    _write_plan(root, all_jobs, jobs, max(int(args.workers), 1))

    existing = sum(1 for job in jobs if _existing_result(job))
    print(f"planned_jobs={len(jobs)}")
    print(f"already_present={existing}")
    print(f"new_jobs={len(jobs) - existing}")
    if args.dry_run:
        return

    exe = ensure_simulator(DEFAULT_EXECUTION)
    print(f"simulator={exe}")

    workers = max(int(args.workers), 1)
    _write_summary(root, all_jobs)

    if workers == 1:
        for index, job in enumerate(jobs, start=1):
            result = _run_one_job(job, exe)
            _append_log(root, result)
            _write_summary(root, all_jobs)
            print(f"[{index}/{len(jobs)}] {result['job_id']} -> {result['status']}", flush=True)
        return

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(_run_one_job, job, exe): job for job in jobs}
        for index, future in enumerate(as_completed(future_map), start=1):
            result = future.result()
            _append_log(root, result)
            _write_summary(root, all_jobs)
            print(f"[{index}/{len(jobs)}] {result['job_id']} -> {result['status']}", flush=True)


if __name__ == "__main__":
    main()