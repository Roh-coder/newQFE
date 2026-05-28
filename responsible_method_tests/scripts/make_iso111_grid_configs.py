#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import sys
from datetime import date
from typing import Any

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_REPO_ROOT = os.path.normpath(os.path.join(_PROJECT_ROOT, ".."))
_KFC_ROOT = os.path.join(_REPO_ROOT, "K_from_continuum")
if _KFC_ROOT not in sys.path:
    sys.path.insert(0, _KFC_ROOT)

from workflow_common import couplings_from_ratio, exact_triangular_ising_beta, save_json  # noqa: E402


def parse_args() -> argparse.Namespace:
    today = date.today().strftime("%Y%m%d")
    parser = argparse.ArgumentParser(
        description="Generate production configs for the 1-1-1 high-stat twisted reference and untwisted local grid scan."
    )
    parser.add_argument("--output-dir", default=os.path.join(_PROJECT_ROOT, "configs"), help="Directory for generated config JSON files")
    parser.add_argument("--results-root", default="results", help="results_root value to write into the generated configs")
    parser.add_argument("--reference-tag", default=None, help="Optional run.tag for the twisted reference campaign")
    parser.add_argument("--grid-tag", default=None, help="Optional run.tag for the untwisted grid campaign")
    parser.add_argument("--reference-config", default=None, help="Optional explicit path for the twisted reference config JSON")
    parser.add_argument("--grid-config", default=None, help="Optional explicit path for the untwisted grid config JSON")
    parser.add_argument("--r1-min", type=float, default=0.5, help="Minimum r1 = k1 / k3")
    parser.add_argument("--r1-max", type=float, default=1.5, help="Maximum r1 = k1 / k3")
    parser.add_argument("--r2-min", type=float, default=0.5, help="Minimum r2 = k2 / k3")
    parser.add_argument("--r2-max", type=float, default=1.5, help="Maximum r2 = k2 / k3")
    parser.add_argument("--step", type=float, default=0.1, help="Grid step for both r1 and r2")
    parser.add_argument("--sizes", type=int, nargs="+", default=[12, 24, 36, 48, 60], help="Family sizes for the twisted reference and untwisted candidates")
    parser.add_argument("--n-traj", type=int, default=10000, help="Production trajectories per lattice")
    parser.add_argument("--n-skip", type=int, default=10, help="Measurement thinning interval")
    parser.add_argument("--n-therm", type=int, default=2000, help="Thermalization trajectories")
    parser.add_argument("--n-workers", type=int, default=2, help="Parallel workers for generate_pointwise_manifold_dataset.py")
    parser.add_argument("--benchmark-id", default="geometry_111", help="Benchmark id for the twisted reference family")
    parser.add_argument("--candidate-prefix", default="iso111_scan", help="Prefix for untwisted grid benchmark ids")
    parser.add_argument("--k3", type=float, default=1.0, help="Base coupling k3 for all generated points")
    parser.add_argument("--fit-method", default="taylor2", help="FSS fit method")
    parser.add_argument("--c-min", type=float, default=0.05, help="Minimum continuum-fit exponent/parameter")
    parser.add_argument("--c-max", type=float, default=3.5, help="Maximum continuum-fit exponent/parameter")
    parser.add_argument("--c-initial", type=float, default=1.0, help="Initial continuum-fit exponent/parameter")
    parser.add_argument("--min-sizes-for-free-c", type=int, default=3, help="Minimum sizes for the free-C fit path")
    parser.add_argument("--twisted-description", default="Fresh high-stat twisted reference for the equilateral 1-1-1 control.", help="Description written into the twisted reference benchmark")
    parser.add_argument("--candidate-description-prefix", default="Equilateral responsible wide scan candidate", help="Description prefix for each untwisted grid point")
    parser.add_argument("--critical-source", default="exact anisotropic triangular critical point", help="critical_source string for untwisted benchmarks")
    args = parser.parse_args()

    if args.reference_tag is None:
        args.reference_tag = f"raw_manifold_fss_iso111_twisted_target_hi{args.n_traj}_{today}"
    if args.grid_tag is None:
        step_token = _format_step_token(args.step)
        args.grid_tag = (
            f"raw_manifold_fss_iso111_grid_r{_compact_bound(args.r1_min)}to{_compact_bound(args.r1_max)}"
            f"_r{_compact_bound(args.r2_min)}to{_compact_bound(args.r2_max)}_step{step_token}_hi{args.n_traj}_{today}"
        )
    if args.reference_config is None:
        args.reference_config = os.path.join(args.output_dir, f"{args.reference_tag}.json")
    if args.grid_config is None:
        args.grid_config = os.path.join(args.output_dir, f"{args.grid_tag}.json")
    return args


def _format_decimal_token(value: float) -> str:
    return f"{value:.3f}".replace("-", "m").replace(".", "p")


def _format_step_token(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def _compact_bound(value: float) -> str:
    return f"{value:.1f}".replace(".", "p")


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    if step <= 0.0:
        raise ValueError("step must be positive")
    if stop < start:
        raise ValueError("grid maximum must be >= minimum")
    count_float = (stop - start) / step
    count = int(round(count_float))
    if not math.isclose(count_float, float(count), rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError("grid range must be an integer number of steps")
    values = [round(start + idx * step, 10) for idx in range(count + 1)]
    values[-1] = round(stop, 10)
    return values


def _base_config(*, tag: str, results_root: str, n_workers: int, n_traj: int, n_skip: int, n_therm: int, fit_method: str, c_min: float, c_max: float, c_initial: float, min_sizes_for_free_c: int) -> dict[str, Any]:
    return {
        "run": {"tag": tag},
        "paths": {
            "results_root": results_root,
            "resume": True,
        },
        "execution": {
            "n_workers": n_workers,
            "exe": None,
            "auto_build": True,
            "force_rebuild": False,
            "build_command": ["make"],
            "build_timeout_s": 600,
        },
        "mc": {
            "n_traj": n_traj,
            "n_skip": n_skip,
            "n_therm": n_therm,
        },
        "fss": {
            "fit_method": fit_method,
            "c_min": c_min,
            "c_max": c_max,
            "c_initial": c_initial,
            "min_sizes_for_free_C": min_sizes_for_free_c,
        },
        "benchmarks": [],
    }


def _reference_benchmark(*, benchmark_id: str, description: str, sizes: list[int], k3: float) -> dict[str, Any]:
    couplings_cfg = {
        "k3": float(k3),
        "k1_over_k3": 1.0,
        "k2_over_k3": 1.0,
    }
    couplings = couplings_from_ratio(couplings_cfg, f"{benchmark_id}.twisted.couplings")
    beta = exact_triangular_ising_beta(couplings)
    return {
        "id": benchmark_id,
        "description": description,
        "modular_target_geometry": [12, 12, 0, 0],
        "twisted": {
            "base_geometry": [12, 12, 0, 0],
            "scales": [int(size // 12) for size in sizes],
            "couplings": couplings_cfg,
            "beta": float(beta),
            "critical_source": "exact isotropic triangular critical point",
        },
    }


def _candidate_benchmark(*, candidate_prefix: str, description_prefix: str, r1: float, r2: float, sizes: list[int], k3: float, critical_source: str) -> dict[str, Any]:
    couplings_cfg = {
        "k3": float(k3),
        "k1_over_k3": float(r1),
        "k2_over_k3": float(r2),
    }
    benchmark_id = f"{candidate_prefix}_r1{_format_decimal_token(r1)}_r2{_format_decimal_token(r2)}"
    couplings = couplings_from_ratio(couplings_cfg, f"{benchmark_id}.untwisted.couplings")
    beta = exact_triangular_ising_beta(couplings)
    return {
        "id": benchmark_id,
        "description": f"{description_prefix} r1={r1:.3f}, r2={r2:.3f}.",
        "modular_target_geometry": [12, 12, 0, 0],
        "untwisted": {
            "sizes": [int(size) for size in sizes],
            "geometry_defaults": {
                "Lx_mult": 1.0,
                "Ly_mult": 1.0,
                "Tx_frac": 0.0,
                "Ty_frac": 0.0,
            },
            "couplings": couplings_cfg,
            "beta": float(beta),
            "critical_source": critical_source,
        },
    }


def main() -> None:
    args = parse_args()
    if any(size <= 0 for size in args.sizes):
        raise ValueError("sizes must all be positive")
    if any(size % 12 != 0 for size in args.sizes):
        raise ValueError("all sizes must be multiples of 12 so the twisted scales remain integral")

    os.makedirs(args.output_dir, exist_ok=True)

    reference_config = _base_config(
        tag=args.reference_tag,
        results_root=args.results_root,
        n_workers=args.n_workers,
        n_traj=args.n_traj,
        n_skip=args.n_skip,
        n_therm=args.n_therm,
        fit_method=args.fit_method,
        c_min=args.c_min,
        c_max=args.c_max,
        c_initial=args.c_initial,
        min_sizes_for_free_c=args.min_sizes_for_free_c,
    )
    reference_config["benchmarks"] = [
        _reference_benchmark(
            benchmark_id=args.benchmark_id,
            description=args.twisted_description,
            sizes=list(args.sizes),
            k3=args.k3,
        )
    ]

    r1_values = _grid_values(args.r1_min, args.r1_max, args.step)
    r2_values = _grid_values(args.r2_min, args.r2_max, args.step)

    grid_config = _base_config(
        tag=args.grid_tag,
        results_root=args.results_root,
        n_workers=args.n_workers,
        n_traj=args.n_traj,
        n_skip=args.n_skip,
        n_therm=args.n_therm,
        fit_method=args.fit_method,
        c_min=args.c_min,
        c_max=args.c_max,
        c_initial=args.c_initial,
        min_sizes_for_free_c=args.min_sizes_for_free_c,
    )
    grid_config["benchmarks"] = [
        _candidate_benchmark(
            candidate_prefix=args.candidate_prefix,
            description_prefix=args.candidate_description_prefix,
            r1=r1,
            r2=r2,
            sizes=list(args.sizes),
            k3=args.k3,
            critical_source=args.critical_source,
        )
        for r1 in r1_values
        for r2 in r2_values
    ]

    save_json(os.path.abspath(args.reference_config), reference_config)
    save_json(os.path.abspath(args.grid_config), grid_config)

    print(f"wrote {os.path.abspath(args.reference_config)}")
    print(f"wrote {os.path.abspath(args.grid_config)}")
    print(f"reference_tag={args.reference_tag}")
    print(f"grid_tag={args.grid_tag}")
    print(f"grid_points={len(r1_values)}x{len(r2_values)}={len(grid_config['benchmarks'])}")


if __name__ == "__main__":
    main()