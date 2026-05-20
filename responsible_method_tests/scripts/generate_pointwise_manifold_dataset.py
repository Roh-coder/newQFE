#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from typing import Any

import numpy as np
from mpmath import mp

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_REPO_ROOT = os.path.normpath(os.path.join(_PROJECT_ROOT, ".."))
_KFC_ROOT = os.path.join(_REPO_ROOT, "K_from_continuum")
if _KFC_ROOT not in sys.path:
    sys.path.insert(0, _KFC_ROOT)

from workflow_common import (  # noqa: E402
    build_test_geometry_map,
    check_required_sections,
    copy_file_atomic,
    couplings_from_ratio,
    ensure_dir,
    ensure_simulator,
    fit_observable_continuum_power,
    load_json,
    load_payload_from_dat,
    log,
    metadata_path_for_data,
    resolve_path,
    run_production_payload,
    save_json,
    timestamp,
    write_dat,
)

DEFAULT_CONFIG_PATH = "../configs/raw_manifold_fss_campaign.json"
METHOD_NAMES = ("twisted", "untwisted")


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


def _boundary_paths(lx: int, ly: int, tx: int, ty: int) -> list[tuple[int, int]]:
    return [
        (lx, ty),
        (tx, -ly),
        (-lx - tx, ly - ty),
    ]


def _normalize_embedding_cycles(value: Any) -> tuple[int, int]:
    if value is None:
        return (0, 1)
    cycles = tuple(int(v) for v in value)
    if len(cycles) != 2:
        raise ValueError("embedding_cycles must contain exactly two cycle indices")
    if len(set(cycles)) != 2 or any(cycle not in (0, 1, 2) for cycle in cycles):
        raise ValueError("embedding_cycles must be two distinct cycle indices chosen from {0,1,2}")
    return cycles


def _to_ab(
    m: int,
    n: int,
    lx: int,
    ly: int,
    tx: int,
    ty: int,
    *,
    embedding_cycles: tuple[int, int] = (0, 1),
) -> tuple[float, float]:
    paths = _boundary_paths(lx, ly, tx, ty)
    (dm_a, dn_a), (dm_b, dn_b) = [paths[idx] for idx in embedding_cycles]
    det = float(dm_a * dn_b - dn_a * dm_b)
    if det == 0.0:
        raise ValueError(
            f"embedding cycle basis {embedding_cycles} is degenerate for lattice {(lx, ly, tx, ty)}"
        )
    a = (float(dn_b) * float(m) - float(dm_b) * float(n)) / det
    b = (-float(dn_a) * float(m) + float(dm_a) * float(n)) / det
    return a, b


def _wrap_unit(value: float) -> float:
    return value % 1.0


def _modular_tau(lx: int, ly: int, tx: int, ty: int) -> tuple[complex, int]:
    omega = 0.5 + 0.5j * np.sqrt(3.0)
    v = complex(lx + ty * omega)
    u = complex(tx - ly * omega)
    tau = u / v
    b_sign = 1
    if tau.imag < 0:
        tau = -tau
        b_sign = -1
    return tau, b_sign


def _ising_torus_shape(nu: complex, tau: complex, theta1p0: mp.mpf, q: mp.mpc) -> float:
    z = mp.pi * mp.mpc(nu.real, nu.imag)
    z2 = 0.5 * z
    th1 = mp.jtheta(1, z, q)
    if abs(th1) == 0:
        return float("nan")

    pref = abs(theta1p0 / th1) ** mp.mpf("0.25")
    num = (
        abs(mp.jtheta(1, z2, q))
        + abs(mp.jtheta(2, z2, q))
        + abs(mp.jtheta(3, z2, q))
        + abs(mp.jtheta(4, z2, q))
    )
    den = (
        abs(mp.jtheta(2, mp.mpf("0.0"), q))
        + abs(mp.jtheta(3, mp.mpf("0.0"), q))
        + abs(mp.jtheta(4, mp.mpf("0.0"), q))
    )
    if den == 0:
        return float("nan")
    return float(pref * (num / den))


def _integer_scale(lattice: tuple[int, int, int, int], smallest: tuple[int, int, int, int]) -> int:
    ratios: list[int] = []
    for current, base in zip(lattice, smallest):
        if base == 0:
            if current != 0:
                raise ValueError(f"lattice {lattice} is not an integer scaling of {smallest}")
            continue
        if current % base != 0:
            raise ValueError(f"lattice {lattice} is not an integer scaling of {smallest}")
        ratios.append(current // base)

    if not ratios:
        return 1
    scale = ratios[0]
    if scale <= 0 or any(r != scale for r in ratios[1:]):
        raise ValueError(f"lattice {lattice} is not a positive uniform scaling of {smallest}")
    return int(scale)


def _normalize_positive_int_sequence(values: Any, field_name: str) -> list[int]:
    items = [int(v) for v in values]
    if len(items) == 0:
        raise ValueError(f"{field_name} must not be empty")
    if any(item <= 0 for item in items):
        raise ValueError(f"{field_name} must contain only positive integers")
    if len(set(items)) != len(items):
        raise ValueError(f"{field_name} must not contain duplicates")
    return sorted(items)


def _build_family_entries(method_cfg: dict[str, Any]) -> list[dict[str, Any]]:
    if "base_geometry" in method_cfg:
        _check_required_keys(method_cfg, "method", ["base_geometry", "scales"])
        base = tuple(int(v) for v in method_cfg["base_geometry"])
        if len(base) != 4:
            raise ValueError("base_geometry must contain 4 integers")
        scales = _normalize_positive_int_sequence(method_cfg["scales"], "scales")
        entries = []
        for scale in scales:
            lattice = tuple(scale * v for v in base)
            entries.append(
                {
                    "scale": int(scale),
                    "family_size": int(scale),
                    "lattice": lattice,
                }
            )
        return entries

    if "sizes" not in method_cfg:
        raise ValueError("method must provide either base_geometry/scales or sizes/geometry_defaults")
    geometry_map = build_test_geometry_map(method_cfg)
    sizes = _normalize_positive_int_sequence(method_cfg["sizes"], "sizes")
    smallest = geometry_map[sizes[0]]
    entries = []
    for size in sizes:
        lattice = geometry_map[size]
        entries.append(
            {
                "scale": _integer_scale(lattice, smallest),
                "family_size": int(size),
                "lattice": lattice,
            }
        )
    return entries


def _point_sort_key(item: tuple[tuple[int, int], dict[str, Any]]) -> tuple[int, int, int]:
    (m, n), row = item
    return int(row["d"]), int(m), int(n)


def _write_payload_fixed_beta(
    *,
    benchmark_id: str,
    method_name: str,
    entry: dict[str, Any],
    couplings: dict[str, float],
    beta: float,
    critical_source: str,
    exe: str,
    mc_cfg: dict[str, Any],
    scratch_root: str,
    grid_root: str,
    resume: bool,
) -> dict[str, Any]:
    scale = int(entry["scale"])
    family_size = int(entry["family_size"])
    lattice = tuple(int(v) for v in entry["lattice"])
    Lx, Ly, Tx, Ty = lattice
    stem = (
        f"{benchmark_id}_{method_name}_scale{scale:02d}_"
        f"Lx{Lx}_Ly{Ly}_Tx{Tx}_Ty{Ty}_"
        f"r1{couplings['r1']:.6f}_r2{couplings['r2']:.6f}"
    )
    data_path = os.path.join(grid_root, stem + ".dat")
    meta_path = metadata_path_for_data(data_path)
    if resume and os.path.exists(data_path) and os.path.exists(meta_path):
        metadata = load_json(meta_path)
        metadata["all_to_all_file"] = data_path
        return metadata

    label = f"{benchmark_id}_{method_name}_scale{scale:02d}"
    production_summary = run_production_payload(
        exe=exe,
        lattice=lattice,
        couplings=couplings,
        beta=float(beta),
        mc_cfg=mc_cfg,
        scratch_root=scratch_root,
        label=label,
    )
    copy_file_atomic(production_summary["all_to_all_file"], data_path)
    metadata = {
        "label": label,
        "created_at": timestamp(),
        "benchmark_id": benchmark_id,
        "method": method_name,
        "scale": scale,
        "family_size": family_size,
        "L": family_size,
        "Lx": Lx,
        "Ly": Ly,
        "Tx": Tx,
        "Ty": Ty,
        "k1": float(couplings["k1"]),
        "k2": float(couplings["k2"]),
        "k3": float(couplings["k3"]),
        "r1": float(couplings["r1"]),
        "r2": float(couplings["r2"]),
        "beta_c": float(beta),
        "beta_c_sigma": 0.0,
        "production_beta": float(beta),
        "production_beta_source": "exact_critical",
        "critical_source": str(critical_source),
        "production_wall_seconds": float(production_summary.get("production_wall_seconds", 0.0)),
        "mc": dict(production_summary["mc"]),
        "all_to_all_file": data_path,
        "source_all_to_all_file": production_summary["all_to_all_file"],
    }
    save_json(meta_path, metadata)
    return metadata


def _fit_modular_alignment(
    continuum_rows: list[list[Any]],
    modular_rows: list[list[Any]],
) -> dict[str, Any]:
    modular_by_point = {int(row[0]): float(row[8]) for row in modular_rows}
    numerator = 0.0
    denominator = 0.0
    chi2 = 0.0
    n_fit = 0
    for row in continuum_rows:
        point_id = int(row[0])
        A = float(row[10])
        sigma_A = float(row[11])
        g = modular_by_point.get(point_id, float("nan"))
        if not np.isfinite(A) or not np.isfinite(sigma_A) or sigma_A <= 0.0 or not np.isfinite(g):
            continue
        weight = 1.0 / (sigma_A * sigma_A)
        numerator += weight * A * g
        denominator += weight * g * g
        n_fit += 1

    alpha = numerator / denominator if denominator > 0.0 else float("nan")
    if np.isfinite(alpha):
        for row in continuum_rows:
            point_id = int(row[0])
            A = float(row[10])
            sigma_A = float(row[11])
            g = modular_by_point.get(point_id, float("nan"))
            if not np.isfinite(A) or not np.isfinite(sigma_A) or sigma_A <= 0.0 or not np.isfinite(g):
                continue
            chi2 += ((A - alpha * g) / sigma_A) ** 2

    dof = max(n_fit - 1, 0)
    return {
        "alpha": alpha,
        "n_fit_points": int(n_fit),
        "chi2": float(chi2) if np.isfinite(alpha) else float("nan"),
        "chi2_per_dof": float(chi2 / dof) if dof > 0 and np.isfinite(alpha) else float("nan"),
    }


def _build_family_dataset(
    *,
    benchmark_id: str,
    benchmark_description: str,
    target_geometry: tuple[int, int, int, int],
    method_name: str,
    method_cfg: dict[str, Any],
    exe: str,
    mc_cfg: dict[str, Any],
    fss_cfg: dict[str, Any],
    run_root: str,
    resume: bool,
    limit_sizes: int | None,
) -> str:
    method_root = os.path.join(run_root, benchmark_id, method_name)
    grid_root = os.path.join(method_root, "grid")
    scratch_root = os.path.join(run_root, "_mc_scratch", benchmark_id, method_name)
    ensure_dir(method_root)
    ensure_dir(grid_root)
    ensure_dir(scratch_root)

    couplings = couplings_from_ratio(method_cfg["couplings"], f"{benchmark_id}.{method_name}.couplings")
    beta = float(method_cfg["beta"])
    critical_source = str(method_cfg.get("critical_source", "exact_critical"))
    entries = _build_family_entries(method_cfg)
    if limit_sizes is not None:
        entries = entries[:max(int(limit_sizes), 1)]
    if len(entries) == 0:
        raise ValueError(f"{benchmark_id}.{method_name} has no sizes after filtering")

    payloads: list[dict[str, Any]] = []
    for entry in entries:
        metadata = _write_payload_fixed_beta(
            benchmark_id=benchmark_id,
            method_name=method_name,
            entry=entry,
            couplings=couplings,
            beta=beta,
            critical_source=critical_source,
            exe=exe,
            mc_cfg=mc_cfg,
            scratch_root=scratch_root,
            grid_root=grid_root,
            resume=resume,
        )
        payload = load_payload_from_dat(metadata["all_to_all_file"], metadata_path_for_data(metadata["all_to_all_file"]))
        payload["scale"] = int(metadata["scale"])
        payload["family_size"] = int(metadata["family_size"])
        payloads.append(payload)

    payloads.sort(key=lambda row: int(row["scale"]))
    smallest = payloads[0]
    smallest_lattice = (
        int(smallest["Lx"]),
        int(smallest["Ly"]),
        int(smallest["Tx"]),
        int(smallest["Ty"]),
    )
    embedding_cycles = _normalize_embedding_cycles(method_cfg.get("embedding_cycles"))
    target_tau, target_b_sign = _modular_tau(*target_geometry)

    mp.dps = 50
    q = mp.e ** (mp.pi * 1j * mp.mpc(target_tau.real, target_tau.imag))
    theta1p0 = mp.diff(lambda zz: mp.jtheta(1, zz, q), mp.mpf("0.0"))

    smallest_rows = sorted(smallest["data"].items(), key=_point_sort_key)
    smallest_point_rows: list[list[Any]] = []
    raw_rows: list[list[Any]] = []
    continuum_rows: list[list[Any]] = []
    modular_rows: list[list[Any]] = []

    scale_values = np.asarray([int(payload["scale"]) for payload in payloads], dtype=float)

    for point_id, ((m0, n0), row0) in enumerate(smallest_rows, start=1):
        a_raw, b_raw = _to_ab(
            m0,
            n0,
            *smallest_lattice,
            embedding_cycles=embedding_cycles,
        )
        a_wrap = _wrap_unit(a_raw)
        b_wrap = _wrap_unit(target_b_sign * b_raw)
        nu_real = float(a_wrap + b_wrap * target_tau.real)
        nu_imag = float(b_wrap * target_tau.imag)
        is_origin = int(abs(nu_real) < 1.0e-14 and abs(nu_imag) < 1.0e-14)
        modular_value = float("nan") if is_origin else _ising_torus_shape(complex(nu_real, nu_imag), target_tau, theta1p0, q)

        smallest_point_rows.append(
            [
                point_id,
                int(row0["d"]),
                int(m0),
                int(n0),
                a_raw,
                b_raw,
                a_wrap,
                b_wrap,
                nu_real,
                nu_imag,
                float(row0["corr"]),
                float(row0["err"]),
                float(row0["conn"]),
                float(row0["conn_err"]),
                is_origin,
            ]
        )
        modular_rows.append(
            [
                point_id,
                int(row0["d"]),
                int(m0),
                int(n0),
                a_wrap,
                b_wrap,
                nu_real,
                nu_imag,
                modular_value,
                is_origin,
            ]
        )

        y_values: list[float] = []
        sigma_values: list[float] = []
        for payload in payloads:
            scale = int(payload["scale"])
            m = int(scale * m0)
            n = int(scale * n0)
            point = payload["data"].get((m, n))
            if point is None:
                raise KeyError(
                    f"missing scaled point {(m, n)} for point {(m0, n0)} at scale={scale} "
                    f"in {benchmark_id}.{method_name}"
                )
            raw_rows.append(
                [
                    point_id,
                    scale,
                    int(payload["family_size"]),
                    int(payload["Lx"]),
                    int(payload["Ly"]),
                    int(payload["Tx"]),
                    int(payload["Ty"]),
                    m,
                    n,
                    int(point["d"]),
                    float(point["corr"]),
                    float(point["err"]),
                    float(point["conn"]),
                    float(point["conn_err"]),
                ]
            )
            y_values.append(float(point["conn"]))
            sigma_values.append(float(point["conn_err"]))

        A, sigma_A, B, C, n_sizes_used, fit_mode = fit_observable_continuum_power(
            scale_values,
            np.asarray(y_values, dtype=float),
            np.asarray(sigma_values, dtype=float),
            fit_method=str(fss_cfg.get("fit_method", "taylor2")),
            c_min=float(fss_cfg["c_min"]),
            c_max=float(fss_cfg["c_max"]),
            c_initial=float(fss_cfg["c_initial"]),
            min_sizes_for_free_C=int(fss_cfg.get("min_sizes_for_free_C", 8)),
        )
        continuum_rows.append(
            [
                point_id,
                int(row0["d"]),
                int(m0),
                int(n0),
                a_raw,
                b_raw,
                a_wrap,
                b_wrap,
                nu_real,
                nu_imag,
                A,
                sigma_A,
                B,
                C,
                int(n_sizes_used),
                str(fit_mode),
                is_origin,
            ]
        )

    alignment = _fit_modular_alignment(continuum_rows, modular_rows)
    modular_aligned_rows: list[list[Any]] = []
    alpha = float(alignment.get("alpha", float("nan")))
    for row in modular_rows:
        raw_value = float(row[8])
        aligned_value = alpha * raw_value if np.isfinite(alpha) and np.isfinite(raw_value) else float("nan")
        modular_aligned_rows.append(list(row) + [aligned_value])

    smallest_points_path = os.path.join(method_root, "smallest_lattice_points.dat")
    raw_points_path = os.path.join(method_root, "pointwise_raw.dat")
    continuum_path = os.path.join(method_root, "pointwise_continuum.dat")
    modular_raw_path = os.path.join(method_root, "modular_raw.dat")
    modular_aligned_path = os.path.join(method_root, "modular_aligned.dat")
    alignment_path = os.path.join(method_root, "modular_alignment.json")
    manifest_path = os.path.join(method_root, f"manifest_{benchmark_id}_{method_name}.json")

    write_dat(
        smallest_points_path,
        [
            f"benchmark={benchmark_id}",
            f"method={method_name}",
            f"description={benchmark_description}",
            "smallest-lattice master point set",
        ],
        [
            "point_id",
            "d",
            "m",
            "n",
            "a_raw",
            "b_raw",
            "a_wrap",
            "b_wrap",
            "nu_real",
            "nu_imag",
            "corr",
            "err",
            "conn",
            "conn_err",
            "is_origin",
        ],
        smallest_point_rows,
    )
    write_dat(
        raw_points_path,
        [
            f"benchmark={benchmark_id}",
            f"method={method_name}",
            "raw connected correlator values for every scaled copy of the smallest-lattice points",
        ],
        [
            "point_id",
            "scale",
            "family_size",
            "Lx",
            "Ly",
            "Tx",
            "Ty",
            "m",
            "n",
            "d",
            "corr",
            "err",
            "conn",
            "conn_err",
        ],
        raw_rows,
    )
    write_dat(
        continuum_path,
        [
            f"benchmark={benchmark_id}",
            f"method={method_name}",
            "pointwise continuum extrapolation of the connected correlator",
        ],
        [
            "point_id",
            "d",
            "m",
            "n",
            "a_raw",
            "b_raw",
            "a_wrap",
            "b_wrap",
            "nu_real",
            "nu_imag",
            "A",
            "sigma_A",
            "B",
            "C",
            "n_sizes_used",
            "fit_mode",
            "is_origin",
        ],
        continuum_rows,
    )
    write_dat(
        modular_raw_path,
        [
            f"benchmark={benchmark_id}",
            f"method={method_name}",
            "raw modular-theta shape sampled on the smallest-lattice normalized coordinates",
        ],
        [
            "point_id",
            "d",
            "m",
            "n",
            "a_wrap",
            "b_wrap",
            "nu_real",
            "nu_imag",
            "modular_raw",
            "is_origin",
        ],
        modular_rows,
    )
    write_dat(
        modular_aligned_path,
        [
            f"benchmark={benchmark_id}",
            f"method={method_name}",
            "modular-theta shape with a single weighted amplitude fit to the continuum lattice data",
        ],
        [
            "point_id",
            "d",
            "m",
            "n",
            "a_wrap",
            "b_wrap",
            "nu_real",
            "nu_imag",
            "modular_raw",
            "is_origin",
            "modular_aligned",
        ],
        modular_aligned_rows,
    )
    save_json(alignment_path, alignment)

    method_manifest = {
        "created_at": timestamp(),
        "benchmark_id": benchmark_id,
        "benchmark_description": benchmark_description,
        "method": method_name,
        "target_geometry": list(target_geometry),
        "target_tau": {
            "real": float(target_tau.real),
            "imag": float(target_tau.imag),
        },
        "critical_source": critical_source,
        "beta": float(beta),
        "couplings": couplings,
        "embedding_cycles": list(embedding_cycles),
        "resume": bool(resume),
        "payloads": [
            {
                "scale": int(payload["scale"]),
                "family_size": int(payload["family_size"]),
                "lattice": [int(payload["Lx"]), int(payload["Ly"]), int(payload["Tx"]), int(payload["Ty"])],
                "data_path": str(payload["all_to_all_file"]),
                "metadata_path": metadata_path_for_data(str(payload["all_to_all_file"])),
            }
            for payload in payloads
        ],
        "smallest_point_count": int(len(smallest_point_rows)),
        "comparable_point_count": int(sum(1 for row in modular_rows if int(row[9]) == 0 and np.isfinite(float(row[8])))),
        "smallest_lattice_points": smallest_points_path,
        "pointwise_raw": raw_points_path,
        "pointwise_continuum": continuum_path,
        "modular_raw": modular_raw_path,
        "modular_aligned": modular_aligned_path,
        "modular_alignment": alignment_path,
    }
    save_json(manifest_path, method_manifest)
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate raw all-to-all datasets plus pointwise continuum extrapolations "
            "for the 1-1-1 and 4-5-6 manifold benchmarks."
        )
    )
    parser.add_argument(
        "--config",
        type=str,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the raw-manifold campaign config JSON",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Optional override for run.tag",
    )
    parser.add_argument(
        "--benchmark",
        action="append",
        dest="benchmarks",
        help="Restrict to one or more benchmark ids",
    )
    parser.add_argument(
        "--method",
        action="append",
        dest="methods",
        choices=list(METHOD_NAMES),
        help="Restrict to one or more method names",
    )
    parser.add_argument(
        "--limit-sizes",
        type=int,
        default=None,
        help="Use only the first N sizes in each selected family",
    )
    return parser.parse_args()


def _resolve_input_path(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return os.path.normpath(path)
    if os.path.exists(path):
        return os.path.normpath(os.path.abspath(path))
    return resolve_path(path, base_dir)


def main() -> None:
    args = parse_args()
    config_path = _resolve_input_path(args.config, _HERE)
    cfg = load_json(config_path)

    check_required_sections(cfg, ["run", "paths", "execution", "mc", "fss", "benchmarks"])
    if args.tag is not None:
        cfg["run"]["tag"] = args.tag

    _check_required_keys(cfg["run"], "run", ["tag"])
    _check_required_keys(cfg["paths"], "paths", ["results_root", "resume"])
    _check_required_keys(cfg["execution"], "execution", ["n_workers"])
    _check_required_keys(cfg["mc"], "mc", ["n_traj", "n_skip", "n_therm"])
    _check_required_keys(cfg["fss"], "fss", ["fit_method", "c_min", "c_max", "c_initial"])

    benchmarks_cfg = list(cfg["benchmarks"])
    if len(benchmarks_cfg) == 0:
        raise ValueError("benchmarks must not be empty")

    selected_benchmarks = set(args.benchmarks or [str(item["id"]) for item in benchmarks_cfg])
    selected_methods = tuple(args.methods or METHOD_NAMES)

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _PROJECT_ROOT)
    run_root = os.path.join(results_root, tag)
    ensure_dir(run_root)
    resume = bool(cfg["paths"]["resume"])

    exe = ensure_simulator(cfg["execution"])
    benchmark_manifest_paths: list[str] = []

    for benchmark_cfg in benchmarks_cfg:
        benchmark_id = str(benchmark_cfg["id"])
        if benchmark_id not in selected_benchmarks:
            continue
        benchmark_description = str(benchmark_cfg.get("description", ""))
        _check_required_keys(benchmark_cfg, benchmark_id, ["modular_target_geometry"])
        target_geometry = tuple(int(v) for v in benchmark_cfg["modular_target_geometry"])
        if len(target_geometry) != 4:
            raise ValueError(f"{benchmark_id}.modular_target_geometry must have 4 integers")

        method_manifests: dict[str, str] = {}
        for method_name in selected_methods:
            if method_name not in benchmark_cfg:
                raise ValueError(f"{benchmark_id} is missing method section: {method_name}")
            log(f"[{benchmark_id}] building {method_name} dataset")
            method_manifests[method_name] = _build_family_dataset(
                benchmark_id=benchmark_id,
                benchmark_description=benchmark_description,
                target_geometry=target_geometry,
                method_name=method_name,
                method_cfg=dict(benchmark_cfg[method_name]),
                exe=exe,
                mc_cfg=dict(cfg["mc"]),
                fss_cfg=dict(cfg["fss"]),
                run_root=run_root,
                resume=resume,
                limit_sizes=args.limit_sizes,
            )

        benchmark_manifest_path = os.path.join(run_root, f"manifest_{benchmark_id}.json")
        existing_method_manifests: dict[str, str] = {}
        if os.path.exists(benchmark_manifest_path):
            existing_benchmark_manifest = load_json(benchmark_manifest_path)
            existing_methods = existing_benchmark_manifest.get("methods")
            if isinstance(existing_methods, dict):
                existing_method_manifests = {str(key): str(value) for key, value in existing_methods.items()}

        benchmark_manifest = {
            "created_at": timestamp(),
            "config_path": config_path,
            "run_tag": tag,
            "benchmark_id": benchmark_id,
            "description": benchmark_description,
            "target_geometry": list(target_geometry),
            "methods": {**existing_method_manifests, **method_manifests},
        }
        save_json(benchmark_manifest_path, benchmark_manifest)
        benchmark_manifest_paths.append(benchmark_manifest_path)

    campaign_manifest_path = os.path.join(run_root, "campaign_manifest.json")
    save_json(
        campaign_manifest_path,
        {
            "created_at": timestamp(),
            "config_path": config_path,
            "run_tag": tag,
            "run_root": run_root,
            "benchmark_manifests": benchmark_manifest_paths,
        },
    )
    log(f"Campaign manifest written: {campaign_manifest_path}")


if __name__ == "__main__":
    main()