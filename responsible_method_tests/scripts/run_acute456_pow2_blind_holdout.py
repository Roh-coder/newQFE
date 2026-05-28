#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


REPO_ROOT = Path(__file__).resolve().parents[2]
K_FROM_CONTINUUM_ROOT = REPO_ROOT / "K_from_continuum"
if str(K_FROM_CONTINUUM_ROOT) not in sys.path:
    sys.path.insert(0, str(K_FROM_CONTINUUM_ROOT))

from workflow_common import exact_triangular_ising_beta, find_beta_for_lattice

SIM_EXE = REPO_ROOT / "K_from_continuum" / "bin" / "ising_tri_twisted_parallelogram"
DEFAULT_OUT_ROOT = REPO_ROOT / "responsible_method_tests" / "stupid_method" / "acute456_literal_sizes_blind_holdout"

SMALL_BASE_TWISTED = (6, 6, 3, 1)
DEFAULT_TWISTED_SMALL_BASE_SCALES = (2, 4, 6, 8, 10, 12, 16)
DEFAULT_UNTWISTED_SIZES = (8, 16, 24, 32, 48, 64, 96)
DEFAULT_HOLDOUT_MIN_VOLUME = 1_000_000
DEFAULT_EMBEDDING_CYCLES = (0, 1)

TWISTED_CFG = {
    "name": "twisted",
    "couplings": {"k1": 1.0, "k2": 1.0, "k3": 1.0},
    "beta": 0.274653072167027,
}
UNTWISTED_CFG = {
    "name": "untwisted",
    "couplings": {"k1": 4.702782819756, "k2": 7.353910143333, "k3": 1.0},
    "beta": 0.066421804006,
}
METHODS = (TWISTED_CFG, UNTWISTED_CFG)

MODEL_ORDER = ("power_free", "pade11", "pade21")
MODEL_LABELS = {
    "power_free": "A + B / V^(omega/2)",
    "pade11": "Pade [1/1]",
    "pade21": "Pade [2/1]",
}
MODEL_STYLES = {
    "power_free": ("#b91c1c", "--"),
    "pade11": ("#c05621", "-"),
    "pade21": ("#65a30d", "-."),
}


def parse_positive_int_list(value: str, *, field_name: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError(f"{field_name} must be a comma-separated list of positive integers")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a shared-point blind-holdout study for the acute 4-5-6 geometry with "
            "literal untwisted square sizes and a twisted family built from a very small "
            "4-5-6 base lattice, then compare both fits against one large twisted blind point."
        )
    )
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--train-traj", type=int, default=3000)
    parser.add_argument("--train-therm", type=int, default=1000)
    parser.add_argument("--train-skip", type=int, default=10)
    parser.add_argument(
        "--untwisted-sizes",
        type=lambda value: parse_positive_int_list(value, field_name="untwisted sizes"),
        default=DEFAULT_UNTWISTED_SIZES,
        help="Comma-separated untwisted square sizes, e.g. 8,16,24,32,48,64,96",
    )
    parser.add_argument(
        "--twisted-small-base-scales",
        type=lambda value: parse_positive_int_list(value, field_name="twisted small-base scales"),
        default=DEFAULT_TWISTED_SMALL_BASE_SCALES,
        help=(
            "Comma-separated multiples of the small twisted base [6,6,3,1]. "
            "Use even scales when the shared point is a boundary midpoint."
        ),
    )
    parser.add_argument("--holdout-traj", type=int, default=500)
    parser.add_argument("--holdout-therm", type=int, default=200)
    parser.add_argument("--holdout-skip", type=int, default=10)
    parser.add_argument(
        "--beta-method",
        choices=("fixed", "beta_finder"),
        default="fixed",
        help="Use the existing fixed beta values or estimate beta_c per lattice from susceptibility scans.",
    )
    parser.add_argument(
        "--beta-window-frac",
        type=float,
        default=0.25,
        help="Half-width of the beta_c finder bracket as a fraction of the nominal beta.",
    )
    parser.add_argument(
        "--beta-window-min",
        type=float,
        default=0.015,
        help="Minimum absolute half-width used for the beta_c finder bracket.",
    )
    parser.add_argument("--beta-n-coarse", type=int, default=11)
    parser.add_argument("--beta-n-refine", type=int, default=5)
    parser.add_argument("--beta-n-refine2", type=int, default=5)
    parser.add_argument("--beta-n-refine3", type=int, default=5)
    parser.add_argument("--beta-n-traj-coarse", type=int, default=80)
    parser.add_argument("--beta-n-traj-fine", type=int, default=160)
    parser.add_argument("--beta-max-shifts", type=int, default=4)
    parser.add_argument(
        "--beta-jackknife",
        action="store_true",
        help="Estimate beta_c uncertainty by leave-one-out jackknife over the scan points.",
    )
    parser.add_argument(
        "--holdout-min-volume",
        type=int,
        default=DEFAULT_HOLDOUT_MIN_VOLUME,
        help="Minimum twisted holdout volume; the first on-lattice boundary-midpoint multiple above this is used.",
    )
    parser.add_argument("--probe-traj", type=int, default=1)
    parser.add_argument("--probe-therm", type=int, default=0)
    parser.add_argument("--probe-skip", type=int, default=1)
    parser.add_argument(
        "--boundary-index",
        type=int,
        choices=(0, 1, 2),
        default=None,
        help="Physical boundary direction whose midpoint is used for the shared point; defaults to the longest side.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def base_volume(lattice: tuple[int, int, int, int]) -> int:
    lx, ly, tx, ty = lattice
    return abs(lx * ly + tx * ty)


def lattice_for_scale(base_lattice: tuple[int, int, int, int], scale: int) -> tuple[int, int, int, int]:
    return tuple(scale * value for value in base_lattice)


def square_lattice(size: int) -> tuple[int, int, int, int]:
    return (size, size, 0, 0)


def lattice_tag(lattice: tuple[int, int, int, int]) -> str:
    lx, ly, tx, ty = lattice
    return f"{lx}x{ly}_t{tx}x{ty}"


def build_twisted_training_lattices(scales: tuple[int, ...]) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(lattice_for_scale(SMALL_BASE_TWISTED, scale) for scale in scales)


def build_untwisted_training_lattices(sizes: tuple[int, ...]) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(square_lattice(size) for size in sizes)


def wrapped_run_id(lattice: tuple[int, int, int, int], couplings: dict[str, float]) -> str:
    lx, ly, tx, ty = lattice
    kt = 0.5 * couplings["k3"]
    return (
        f"{lx}x{ly}_t{tx}x{ty}_k"
        f"{int(round(couplings['k1'] * 1000))}_"
        f"{int(round(couplings['k2'] * 1000))}_"
        f"{int(round(couplings['k3'] * 1000))}_"
        f"{int(round(kt * 1000))}"
    )


def stable_seed(label: str) -> int:
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big") or 1


def latest_point_file(data_dir: Path) -> Path | None:
    candidates = sorted(data_dir.rglob("two_point_all_to_all.dat"), key=lambda path: path.stat().st_mtime)
    return candidates[-1] if candidates else None


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected object JSON in {path}")
    return value


def save_json(path: Path, payload: dict[str, Any]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def couplings_with_ratios(method_cfg: dict[str, Any]) -> dict[str, float]:
    couplings = dict(method_cfg["couplings"])
    k3 = float(couplings["k3"])
    if k3 == 0.0:
        raise ValueError("k3 must be non-zero")
    couplings["r1"] = float(couplings["k1"]) / k3
    couplings["r2"] = float(couplings["k2"]) / k3
    return couplings


def beta_finder_cfg(nominal_beta: float, args: argparse.Namespace) -> dict[str, Any]:
    half_width = max(abs(float(nominal_beta)) * float(args.beta_window_frac), float(args.beta_window_min))
    beta_lo = max(1.0e-6, float(nominal_beta) - half_width)
    beta_hi = float(nominal_beta) + half_width
    return {
        "mode": "scan",
        "beta_lo": beta_lo,
        "beta_hi": beta_hi,
        "n_coarse": int(args.beta_n_coarse),
        "n_refine": int(args.beta_n_refine),
        "n_refine2": int(args.beta_n_refine2),
        "n_refine3": int(args.beta_n_refine3),
        "n_traj_coarse": int(args.beta_n_traj_coarse),
        "n_traj_fine": int(args.beta_n_traj_fine),
        "max_shifts": int(args.beta_max_shifts),
        "jackknife": bool(args.beta_jackknife),
    }


def resolve_beta_summary(
    *,
    method_cfg: dict[str, Any],
    lattice: tuple[int, int, int, int],
    args: argparse.Namespace,
    summary_dir: Path,
    label: str,
    force: bool,
) -> tuple[dict[str, Any], Path | None]:
    couplings = couplings_with_ratios(method_cfg)
    exact_beta = float(exact_triangular_ising_beta(couplings))
    nominal_beta = float(method_cfg["beta"])

    if args.beta_method == "fixed":
        summary = {
            "label": label,
            "method_name": method_cfg["name"],
            "beta_method": "fixed",
            "baseline_beta": nominal_beta,
            "exact_beta": exact_beta,
            "beta_c": nominal_beta,
            "beta_c_sigma": 0.0,
            "chi_peak": float("nan"),
            "scan_betas": [],
            "scan_chis": [],
            "scan_chi_errs": [],
            "beta_finder": {"mode": "fixed"},
        }
        return summary, None

    summary_path = summary_dir / f"{label}_beta_scan.json"
    if summary_path.is_file() and not force:
        cached = load_json(summary_path)
        return cached, summary_path

    beta_cfg = beta_finder_cfg(nominal_beta, args)
    scratch_root = summary_dir / "scratch"
    summary = find_beta_for_lattice(
        str(SIM_EXE),
        lattice,
        couplings,
        beta_cfg,
        str(scratch_root),
        label,
    )
    summary["method_name"] = method_cfg["name"]
    summary["beta_method"] = "beta_finder"
    summary["baseline_beta"] = nominal_beta
    summary["exact_beta"] = exact_beta
    summary["scan_dir"] = str(scratch_root / label / "scan")
    save_json(summary_path, summary)
    return summary, summary_path


def run_simulation(
    *,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta: float,
    n_traj: int,
    n_skip: int,
    n_therm: int,
    data_dir: Path,
    label: str,
    single_disp: tuple[int, int] | None = None,
    force: bool = False,
) -> tuple[Path, str]:
    ensure_dir(data_dir)
    existing = None if force else latest_point_file(data_dir)
    if existing is not None:
        print(f"reusing {label}: {existing}", flush=True)
        return existing, ""

    if not SIM_EXE.is_file():
        raise FileNotFoundError(f"simulator not found: {SIM_EXE}")

    lx, ly, tx, ty = lattice
    cmd = [
        str(SIM_EXE),
        "--L_x", str(lx),
        "--L_y", str(ly),
        "--T_x", str(tx),
        "--T_y", str(ty),
        "--k1", f"{couplings['k1']:.12f}",
        "--k2", f"{couplings['k2']:.12f}",
        "--k3", f"{couplings['k3']:.12f}",
        "--beta", f"{beta:.15f}",
        "--seed", str(stable_seed(label)),
        "--n_therm", str(n_therm),
        "--n_traj", str(n_traj),
        "--n_skip", str(n_skip),
        "--data_dir", str(data_dir),
    ]
    if single_disp is not None:
        cmd.extend(["--single_disp_m", str(single_disp[0]), "--single_disp_n", str(single_disp[1])])

    print(
        (
            f"running {label}: lattice={list(lattice)} "
            f"disp={list(single_disp) if single_disp is not None else 'all'} "
            f"mc=(n_therm={n_therm}, n_traj={n_traj}, n_skip={n_skip})"
        ),
        flush=True,
    )

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"simulation failed for {label}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    output_file = latest_point_file(data_dir)
    if output_file is None:
        raise RuntimeError(f"no two_point_all_to_all.dat found under {data_dir} for {label}")
    print(f"finished {label}: {output_file}", flush=True)
    return output_file, result.stdout


def load_all_points(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            rows.append(
                {
                    "d": int(parts[0]),
                    "m": int(parts[1]),
                    "n": int(parts[2]),
                    "corr": float(parts[3]),
                    "err": float(parts[4]),
                    "conn": float(parts[5]),
                    "conn_err": float(parts[6]),
                }
            )
    return rows


def load_single_point(path: Path) -> dict[str, float | int]:
    rows = load_all_points(path)
    if len(rows) != 1:
        raise ValueError(f"expected exactly one selected point in {path}, found {len(rows)}")
    return rows[0]


def boundary_paths(lx: int, ly: int, tx: int, ty: int) -> list[tuple[int, int]]:
    return [(lx, ty), (tx, -ly), (-lx - tx, ly - ty)]


def triangular_length_sq(m: int, n: int) -> int:
    return m * m + n * n + m * n


def longest_boundary_index(lattice: tuple[int, int, int, int]) -> int:
    paths = boundary_paths(*lattice)
    return max(range(3), key=lambda idx: triangular_length_sq(paths[idx][0], paths[idx][1]))


def midpoint_disp(lattice: tuple[int, int, int, int], boundary_index: int) -> tuple[int, int]:
    dm, dn = boundary_paths(*lattice)[boundary_index]
    if dm % 2 != 0 or dn % 2 != 0:
        raise ValueError(
            f"boundary midpoint is not a lattice point for lattice={lattice} boundary_index={boundary_index}"
        )
    return (dm // 2, dn // 2)


def midpoint_coord(boundary_index: int) -> tuple[float, float]:
    if boundary_index == 0:
        return (0.5, 0.0)
    if boundary_index == 1:
        return (0.0, 0.5)
    return (0.5, 0.5)


def to_ab(m: int, n: int, lx: int, ly: int, tx: int, ty: int, cycles: tuple[int, int] = (0, 1)) -> tuple[float, float]:
    paths = boundary_paths(lx, ly, tx, ty)
    (dm_a, dn_a), (dm_b, dn_b) = [paths[index] for index in cycles]
    det = float(dm_a * dn_b - dn_a * dm_b)
    a = (float(dn_b) * float(m) - float(dm_b) * float(n)) / det
    b = (-float(dn_a) * float(m) + float(dm_a) * float(n)) / det
    return a, b


def wrap_unit(value: float) -> float:
    return value % 1.0


def coordinate_map(rows: list[dict[str, float | int]], lattice: tuple[int, int, int, int]) -> dict[tuple[float, float], dict[str, float | int]]:
    lx, ly, tx, ty = lattice
    mapped: dict[tuple[float, float], dict[str, float | int]] = {}
    for row in rows:
        a_raw, b_raw = to_ab(int(row["m"]), int(row["n"]), lx, ly, tx, ty)
        key = (round(wrap_unit(a_raw), 12), round(wrap_unit(b_raw), 12))
        mapped[key] = {
            **row,
            "a_wrap": wrap_unit(a_raw),
            "b_wrap": wrap_unit(b_raw),
        }
    return mapped


def choose_shared_point(
    boundary_index: int,
    *,
    twisted_reference_lattice: tuple[int, int, int, int],
    untwisted_reference_lattice: tuple[int, int, int, int],
) -> dict[str, Any]:
    key = midpoint_coord(boundary_index)
    return {
        "coord": {"a_wrap": key[0], "b_wrap": key[1]},
        "twisted": {
            "lattice": list(twisted_reference_lattice),
            "m": int(midpoint_disp(twisted_reference_lattice, boundary_index)[0]),
            "n": int(midpoint_disp(twisted_reference_lattice, boundary_index)[1]),
        },
        "untwisted": {
            "lattice": list(untwisted_reference_lattice),
            "m": int(midpoint_disp(untwisted_reference_lattice, boundary_index)[0]),
            "n": int(midpoint_disp(untwisted_reference_lattice, boundary_index)[1]),
        },
        "boundary_index": int(boundary_index),
    }


def validate_boundary_midpoint_family(
    lattices: tuple[tuple[int, int, int, int], ...],
    *,
    boundary_index: int,
    family_name: str,
) -> None:
    for lattice in lattices:
        midpoint_disp(lattice, boundary_index)


def robust_sigma(values: np.ndarray) -> np.ndarray:
    good = np.isfinite(values) & (values > 0)
    median = float(np.median(values[good])) if np.any(good) else 1.0
    return np.where(good, values, median)


def fit_basis_powers(x: np.ndarray, y: np.ndarray, sigma: np.ndarray, powers: tuple[int, ...]) -> dict[str, Any]:
    sigma = robust_sigma(sigma)
    design = np.column_stack([x ** power for power in powers])
    design_w = design / sigma[:, None]
    y_w = y / sigma
    beta, *_ = np.linalg.lstsq(design_w, y_w, rcond=None)
    y_hat = design @ beta
    chi2 = float(np.sum(((y - y_hat) / sigma) ** 2))
    dof = max(len(y) - len(powers), 1)
    return {
        "coef": [float(value) for value in beta],
        "chi2dof_train": chi2 / dof,
        "predict": lambda x_new: float(sum(beta[idx] * float(x_new) ** power for idx, power in enumerate(powers))),
        "curve": lambda x_grid: sum(beta[idx] * x_grid ** power for idx, power in enumerate(powers)),
    }


def model_power_free(x: np.ndarray, A: float, B: float, omega: float) -> np.ndarray:
    return A + B * np.power(x, omega)


def fit_power_free(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict[str, Any]:
    sigma = robust_sigma(sigma)
    A0 = float(np.mean(y[-2:]))
    B0 = float(y[0] - A0)
    popt, _ = curve_fit(
        model_power_free,
        x,
        y,
        p0=[A0, B0, 1.0],
        sigma=sigma,
        absolute_sigma=True,
        bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 4.0]),
        maxfev=50000,
    )
    y_hat = model_power_free(x, *popt)
    chi2 = float(np.sum(((y - y_hat) / sigma) ** 2))
    dof = max(len(y) - 3, 1)
    return {
        "coef": [float(value) for value in popt],
        "chi2dof_train": chi2 / dof,
        "predict": lambda x_new: float(model_power_free(np.asarray([x_new]), *popt)[0]),
        "curve": lambda x_grid: model_power_free(x_grid, *popt),
    }


def model_pade11(x: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
    return A + (B * x) / (1.0 + C * x)


def fit_pade11(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict[str, Any]:
    sigma = robust_sigma(sigma)
    A0 = float(np.mean(y[-2:]))
    B0 = float((y[0] - A0) / max(x[0], 1.0e-12))
    c_low = -0.95 / float(np.max(x))
    popt, _ = curve_fit(
        model_pade11,
        x,
        y,
        p0=[A0, B0, 0.0],
        sigma=sigma,
        absolute_sigma=True,
        bounds=([-np.inf, -np.inf, c_low], [np.inf, np.inf, 500.0]),
        maxfev=50000,
    )
    y_hat = model_pade11(x, *popt)
    chi2 = float(np.sum(((y - y_hat) / sigma) ** 2))
    dof = max(len(y) - 3, 1)
    return {
        "coef": [float(value) for value in popt],
        "chi2dof_train": chi2 / dof,
        "predict": lambda x_new: float(model_pade11(np.asarray([x_new]), *popt)[0]),
        "curve": lambda x_grid: model_pade11(x_grid, *popt),
    }


def model_pade21(x: np.ndarray, A: float, B: float, C: float, D: float) -> np.ndarray:
    return A + (B * x + C * x * x) / (1.0 + D * x)


def fit_pade21(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict[str, Any]:
    sigma = robust_sigma(sigma)
    A0 = float(np.mean(y[-2:]))
    B0 = float((y[0] - A0) / max(x[0], 1.0e-12))
    popt, _ = curve_fit(
        model_pade21,
        x,
        y,
        p0=[A0, B0, 0.0, 0.0],
        sigma=sigma,
        absolute_sigma=True,
        bounds=([-np.inf, -np.inf, -np.inf, -0.95 / float(np.max(x))], [np.inf, np.inf, np.inf, 500.0]),
        maxfev=100000,
    )
    y_hat = model_pade21(x, *popt)
    chi2 = float(np.sum(((y - y_hat) / sigma) ** 2))
    dof = max(len(y) - 4, 1)
    return {
        "coef": [float(value) for value in popt],
        "chi2dof_train": chi2 / dof,
        "predict": lambda x_new: float(model_pade21(np.asarray([x_new]), *popt)[0]),
        "curve": lambda x_grid: model_pade21(x_grid, *popt),
    }


FITTERS = {
    "taylor2": lambda x, y, sigma: fit_basis_powers(x, y, sigma, (0, 1, 2)),
    "power_free": fit_power_free,
    "pade11": fit_pade11,
    "pade21": fit_pade21,
}


def format_params(model_name: str, coef: list[float]) -> str:
    if model_name == "taylor2":
        return f"A={coef[0]:.4f}, B={coef[1]:.4f}, C={coef[2]:.4f}"
    if model_name == "power_free":
        return f"A={coef[0]:.4f}, B={coef[1]:.4f}, omega={coef[2]:.4f}"
    if model_name == "pade11":
        return f"A={coef[0]:.4f}, B={coef[1]:.4f}, C={coef[2]:.4f}"
    return f"A={coef[0]:.4f}, B={coef[1]:.4f}, C={coef[2]:.4f}, D={coef[3]:.4f}"


def collect_method_series(
    *,
    method_cfg: dict[str, Any],
    lattices: tuple[tuple[int, int, int, int], ...],
    boundary_index: int,
    mc_cfg: dict[str, int],
    out_root: Path,
    args: argparse.Namespace,
    force: bool,
) -> list[dict[str, Any]]:
    series: list[dict[str, Any]] = []
    beta_summary_dir = out_root / "beta_scans" / method_cfg["name"]
    for family_index, lattice in enumerate(lattices, start=1):
        disp = midpoint_disp(lattice, boundary_index)
        run_dir = out_root / "mc_runs" / method_cfg["name"] / lattice_tag(lattice)
        label = f"acute456_{method_cfg['name']}_{lattice_tag(lattice)}_m{disp[0]}_n{disp[1]}"
        beta_summary, beta_summary_path = resolve_beta_summary(
            method_cfg=method_cfg,
            lattice=lattice,
            args=args,
            summary_dir=beta_summary_dir,
            label=label,
            force=force,
        )
        point_path, stdout = run_simulation(
            lattice=lattice,
            couplings=couplings_with_ratios(method_cfg),
            beta=float(beta_summary["beta_c"]),
            data_dir=run_dir,
            label=label,
            single_disp=disp,
            force=force,
            **mc_cfg,
        )
        point = load_single_point(point_path)
        series.append(
            {
                "family_index": family_index,
                "label": lattice_tag(lattice),
                "lattice": list(lattice),
                "disp": [disp[0], disp[1]],
                "volume": float(base_volume(lattice)),
                "sqrt_volume": math.sqrt(float(base_volume(lattice))),
                "conn": float(point["conn"]),
                "conn_err": float(point["conn_err"]),
                "beta_c": float(beta_summary["beta_c"]),
                "beta_c_sigma": float(beta_summary.get("beta_c_sigma", 0.0)),
                "beta_summary": beta_summary,
                "beta_summary_path": str(beta_summary_path) if beta_summary_path is not None else None,
                "path": str(point_path),
                "stdout_tail": stdout.strip().splitlines()[-8:] if stdout else [],
            }
        )
    return series


def run_holdout(
    *,
    boundary_index: int,
    holdout_min_volume: int,
    mc_cfg: dict[str, int],
    out_root: Path,
    args: argparse.Namespace,
    force: bool,
) -> dict[str, Any]:
    scale = 1
    while True:
        lattice = lattice_for_scale(SMALL_BASE_TWISTED, scale)
        if base_volume(lattice) <= holdout_min_volume:
            scale += 1
            continue
        try:
            disp = midpoint_disp(lattice, boundary_index)
        except ValueError:
            scale += 1
            continue
        break

    run_dir = out_root / "mc_runs" / "holdout_twisted"
    label = f"acute456_holdout_smallbase{scale:03d}_m{disp[0]}_n{disp[1]}"
    beta_summary, beta_summary_path = resolve_beta_summary(
        method_cfg=TWISTED_CFG,
        lattice=lattice,
        args=args,
        summary_dir=out_root / "beta_scans" / "holdout_twisted",
        label=label,
        force=force,
    )
    point_path, stdout = run_simulation(
        lattice=lattice,
        couplings=couplings_with_ratios(TWISTED_CFG),
        beta=float(beta_summary["beta_c"]),
        data_dir=run_dir,
        label=label,
        single_disp=disp,
        force=force,
        **mc_cfg,
    )
    point = load_single_point(point_path)
    return {
        "small_base_scale": scale,
        "min_volume_target": int(holdout_min_volume),
        "label": lattice_tag(lattice),
        "lattice": list(lattice),
        "disp": [disp[0], disp[1]],
        "volume": float(base_volume(lattice)),
        "sqrt_volume": math.sqrt(float(base_volume(lattice))),
        "conn": float(point["conn"]),
        "conn_err": float(point["conn_err"]),
        "beta_c": float(beta_summary["beta_c"]),
        "beta_c_sigma": float(beta_summary.get("beta_c_sigma", 0.0)),
        "beta_summary": beta_summary,
        "beta_summary_path": str(beta_summary_path) if beta_summary_path is not None else None,
        "path": str(point_path),
        "stdout_tail": stdout.strip().splitlines()[-8:] if stdout else [],
    }


def fit_method_series(series: list[dict[str, Any]], holdout: dict[str, Any]) -> dict[str, Any]:
    volumes = np.asarray([item["volume"] for item in series], dtype=float)
    x = 1.0 / np.sqrt(volumes)
    y = np.asarray([item["conn"] for item in series], dtype=float)
    sigma = np.asarray([item["conn_err"] for item in series], dtype=float)
    x_holdout = 1.0 / math.sqrt(float(holdout["volume"]))
    result: dict[str, Any] = {"training_series": series, "fits": {}}
    for model_name in MODEL_ORDER:
        fit = FITTERS[model_name](x, y, sigma)
        pred = float(fit["predict"](x_holdout))
        delta = pred - float(holdout["conn"])
        z_holdout = delta / float(holdout["conn_err"])
        result["fits"][model_name] = {
            "prediction_holdout": pred,
            "delta_vs_holdout": delta,
            "z_holdout": z_holdout,
            "chi2dof_train": float(fit["chi2dof_train"]),
            "coef": list(fit["coef"]),
        }
        fit["prediction_holdout"] = pred
        fit["z_holdout"] = z_holdout
        result["fits"][model_name]["_fit_obj"] = fit
    ranking = sorted(result["fits"].items(), key=lambda item: abs(item[1]["z_holdout"]))
    result["ranking"] = [name for name, _ in ranking]
    result["best_by_abs_z"] = ranking[0][0]
    return result


def strip_fit_objects(summary: dict[str, Any]) -> dict[str, Any]:
    encoded = json.loads(json.dumps(summary, default=str))
    for method_name in ("twisted", "untwisted"):
        for model_name in MODEL_ORDER:
            encoded["methods"][method_name]["fits"][model_name].pop("_fit_obj", None)
    return encoded


def collapse_duplicate_scan_points(
    betas: np.ndarray,
    chis: np.ndarray,
    chi_errs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if betas.size == 0:
        return betas, chis, chi_errs

    rounded = np.round(betas, 12)
    unique_betas = np.unique(rounded)
    merged_betas: list[float] = []
    merged_chis: list[float] = []
    merged_errs: list[float] = []
    has_errs = chi_errs.size == betas.size

    for beta_key in unique_betas:
        mask = rounded == beta_key
        merged_betas.append(float(np.mean(betas[mask])))

        if has_errs:
            errs = chi_errs[mask]
            valid = np.isfinite(errs) & (errs > 0.0)
            if np.any(valid):
                weights = 1.0 / np.square(errs[valid])
                merged_chis.append(float(np.average(chis[mask][valid], weights=weights)))
                merged_errs.append(float(math.sqrt(1.0 / np.sum(weights))))
                continue

        merged_chis.append(float(np.mean(chis[mask])))
        merged_errs.append(float("nan"))

    return (
        np.asarray(merged_betas, dtype=float),
        np.asarray(merged_chis, dtype=float),
        np.asarray(merged_errs, dtype=float),
    )


def build_susceptibility_grid(entries: list[dict[str, Any]], *, title: str, plot_path: Path) -> None:
    if not entries:
        return

    ncols = 3
    nrows = int(math.ceil(len(entries) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.3 * ncols, 3.4 * nrows), constrained_layout=True)
    flat_axes = np.atleast_1d(axes).ravel()
    for ax, entry in zip(flat_axes, entries):
        beta_summary = entry["beta_summary"]
        betas = np.asarray(beta_summary.get("scan_betas", []), dtype=float)
        chis = np.asarray(beta_summary.get("scan_chis", []), dtype=float)
        chi_errs = np.asarray(beta_summary.get("scan_chi_errs", []), dtype=float)
        exact_beta = float(beta_summary.get("exact_beta", np.nan))
        beta_c = float(beta_summary["beta_c"])
        beta_sigma = float(beta_summary.get("beta_c_sigma", 0.0))

        if betas.size:
            order = np.argsort(betas)
            betas = betas[order]
            chis = chis[order]
            if chi_errs.size == order.size:
                chi_errs = chi_errs[order]
            betas, chis, chi_errs = collapse_duplicate_scan_points(betas, chis, chi_errs)
            if chi_errs.size == betas.size and np.any(np.isfinite(chi_errs) & (chi_errs > 0.0)):
                ax.errorbar(betas, chis, yerr=chi_errs, fmt="o-", color="#1f2937", ms=3.5, capsize=2.0)
            else:
                ax.plot(betas, chis, "o-", color="#1f2937", ms=3.5)
        else:
            ax.text(0.5, 0.5, "fixed beta\n(no scan)", transform=ax.transAxes, ha="center", va="center")

        if np.isfinite(exact_beta):
            ax.axvline(exact_beta, color="#2563eb", linestyle=":", linewidth=1.5, label="exact beta")
        ax.axvline(beta_c, color="#b91c1c", linestyle="--", linewidth=1.5, label="beta_c")
        ax.set_title(entry["title"], fontsize=10)
        ax.set_xlabel("beta")
        ax.set_ylabel("m_susc")
        ax.grid(True, alpha=0.25)
        ax.text(
            0.02,
            0.98,
            f"beta_c={beta_c:.8f}\nσ={beta_sigma:.3g}\nV={entry['volume']:.0f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cbd5e1"},
        )

    for ax in flat_axes[len(entries):]:
        ax.axis("off")

    handles, labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", frameon=True)
    fig.suptitle(title, fontsize=14)
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)


def build_beta_pc_continuum_plot(summary: dict[str, Any], plot_path: Path) -> dict[str, Any]:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), constrained_layout=True)
    fit_payload: dict[str, Any] = {}

    for ax, method_name in zip(axes, ("untwisted", "twisted")):
        series = summary["methods"][method_name]["training_series"]
        x = np.asarray([1.0 / math.sqrt(float(item["volume"])) for item in series], dtype=float)
        beta_pc = np.asarray([float(item["beta_c"]) for item in series], dtype=float)
        beta_sigma = np.asarray([float(item.get("beta_c_sigma", 0.0)) for item in series], dtype=float)
        exact_beta = float(series[0]["beta_summary"].get("exact_beta", np.nan)) if series else float("nan")

        fit = fit_basis_powers(x, beta_pc, beta_sigma, (0, 1, 2))
        x_grid = np.linspace(0.0, float(np.max(x)) * 1.05, 500)
        y_grid = np.asarray(fit["curve"](x_grid), dtype=float)
        ax.errorbar(
            x,
            beta_pc,
            yerr=np.where(beta_sigma > 0.0, beta_sigma, np.nan),
            fmt="o",
            color="#111827",
            ms=5.0,
            capsize=2.0,
            label="training beta_pc",
        )
        ax.plot(x_grid, y_grid, color="#b91c1c", linewidth=2.0, label="taylor2 fit")
        if np.isfinite(exact_beta):
            ax.axhline(exact_beta, color="#2563eb", linestyle=":", linewidth=1.6, label="exact beta")

        holdout_row: dict[str, Any] | None = None
        if method_name == "twisted":
            holdout_row = summary["holdout"]
            x_holdout = 1.0 / math.sqrt(float(holdout_row["volume"]))
            ax.errorbar(
                [x_holdout],
                [float(holdout_row["beta_c"])],
                yerr=[float(holdout_row.get("beta_c_sigma", 0.0))],
                fmt="*",
                color="#7c3aed",
                ms=14,
                capsize=3.0,
                label="holdout beta_pc",
            )
        ax.set_xlabel("1 / sqrt(volume)")
        ax.set_ylabel("beta_pc")
        ax.set_title(f"Acute 4-5-6 {method_name} beta_pc continuum fit")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8.0, loc="best", frameon=True)
        ax.text(
            0.02,
            0.98,
            (
                f"beta_inf={fit['coef'][0]:.8f}\n"
                f"chi2/dof={fit['chi2dof_train']:.3f}\n"
                f"B={fit['coef'][1]:.5f}, C={fit['coef'][2]:.5f}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cbd5e1"},
        )

        fit_payload[method_name] = {
            "exact_beta": exact_beta,
            "beta_inf": float(fit["coef"][0]),
            "coef": [float(value) for value in fit["coef"]],
            "chi2dof_train": float(fit["chi2dof_train"]),
        }
        if holdout_row is not None:
            x_holdout = 1.0 / math.sqrt(float(holdout_row["volume"]))
            holdout_pred = float(fit["predict"](x_holdout))
            holdout_sigma = float(holdout_row.get("beta_c_sigma", 0.0))
            fit_payload[method_name]["holdout_beta_pc"] = float(holdout_row["beta_c"])
            fit_payload[method_name]["holdout_prediction"] = holdout_pred
            fit_payload[method_name]["holdout_delta"] = holdout_pred - float(holdout_row["beta_c"])
            fit_payload[method_name]["holdout_z"] = (
                (holdout_pred - float(holdout_row["beta_c"])) / holdout_sigma
                if holdout_sigma > 0.0 else float("nan")
            )

    fig.suptitle("Acute 4-5-6 beta_pc fits to the continuum", fontsize=15)
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    return fit_payload


def build_plot(summary: dict[str, Any], plot_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.4), constrained_layout=True)
    holdout = summary["holdout"]
    x_holdout = 1.0 / math.sqrt(float(holdout["volume"]))
    for ax, method_name in zip(axes, ("untwisted", "twisted")):
        method = summary["methods"][method_name]
        series = method["training_series"]
        x = np.asarray([1.0 / math.sqrt(float(item["volume"])) for item in series], dtype=float)
        y = np.asarray([item["conn"] for item in series], dtype=float)
        sigma = np.asarray([item["conn_err"] for item in series], dtype=float)
        ax.errorbar(x, y, yerr=sigma, fmt="o", color="#111827", ms=5.0, capsize=2.0, label="training data")
        ax.errorbar([x_holdout], [holdout["conn"]], yerr=[holdout["conn_err"]], fmt="*", color="#b91c1c", ms=15, capsize=3, label="twisted blind holdout")

        x_min = min(float(np.min(x)), x_holdout)
        x_max = max(float(np.max(x)), x_holdout)
        x_grid = np.geomspace(x_min * 0.9, x_max * 1.05, 600)

        param_lines = []
        for model_name in MODEL_ORDER:
            fit = method["fits"][model_name]["_fit_obj"]
            color, linestyle = MODEL_STYLES[model_name]
            ax.plot(x_grid, fit["curve"](x_grid), color=color, linestyle=linestyle, linewidth=2.0, label=f"{MODEL_LABELS[model_name]} (z={fit['z_holdout']:+.2f})")
            param_lines.append(
                f"{MODEL_LABELS[model_name]}: chi2/dof={fit['chi2dof_train']:.3f}, z={fit['z_holdout']:+.2f}\n  {format_params(model_name, fit['coef'])}"
            )

        ax.set_xscale("log")
        ax.set_xlabel("1 / sqrt(volume)")
        ax.set_ylabel("connected correlator")
        ax.set_title(f"Acute 4-5-6 {method_name} fit vs twisted blind holdout")
        ax.grid(True, alpha=0.3, which="both")
        ax.legend(fontsize=8.0, loc="best", frameon=True)
        ax.text(
            0.02,
            0.98,
            "\n".join(param_lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.0,
            family="monospace",
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cbd5e1"},
        )
        ax.text(
            0.98,
            0.02,
            (
                f"holdout small-base scale={holdout['small_base_scale']}, V={holdout['volume']:.0f}, conn={holdout['conn']:.4f} +/- {holdout['conn_err']:.4f}\n"
                f"min target volume={holdout['min_volume_target']}\n"
                f"twisted blind lattice={holdout['lattice']} disp={holdout['disp']}"
            ),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.0,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cbd5e1"},
        )

    fig.suptitle("Acute 4-5-6 blind holdout with literal untwisted sizes and a small-base twisted family", fontsize=15)
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dir(args.out_root)

    untwisted_sizes = tuple(int(size) for size in args.untwisted_sizes)
    twisted_small_base_scales = tuple(int(scale) for scale in args.twisted_small_base_scales)
    train_cfg = {"n_traj": args.train_traj, "n_skip": args.train_skip, "n_therm": args.train_therm}
    holdout_cfg = {"n_traj": args.holdout_traj, "n_skip": args.holdout_skip, "n_therm": args.holdout_therm}
    boundary_index = longest_boundary_index(lattice_for_scale(SMALL_BASE_TWISTED, 2)) if args.boundary_index is None else int(args.boundary_index)

    untwisted_lattices = build_untwisted_training_lattices(untwisted_sizes)
    twisted_lattices = build_twisted_training_lattices(twisted_small_base_scales)
    validate_boundary_midpoint_family(untwisted_lattices, boundary_index=boundary_index, family_name="untwisted")
    validate_boundary_midpoint_family(twisted_lattices, boundary_index=boundary_index, family_name="twisted")

    shared_point = choose_shared_point(
        boundary_index,
        twisted_reference_lattice=twisted_lattices[0],
        untwisted_reference_lattice=untwisted_lattices[0],
    )
    print(
        (
            "shared boundary midpoint selected: "
            f"boundary_index={shared_point['boundary_index']} "
            f"twisted_ref={shared_point['twisted']['lattice']} "
            f"twisted=({shared_point['twisted']['m']},{shared_point['twisted']['n']}) "
            f"untwisted_ref={shared_point['untwisted']['lattice']} "
            f"untwisted=({shared_point['untwisted']['m']},{shared_point['untwisted']['n']})"
        ),
        flush=True,
    )
    print(f"twisted training lattices: {[list(lattice) for lattice in twisted_lattices]}", flush=True)
    print(f"untwisted training lattices: {[list(lattice) for lattice in untwisted_lattices]}", flush=True)
    holdout = run_holdout(
        boundary_index=boundary_index,
        holdout_min_volume=int(args.holdout_min_volume),
        mc_cfg=holdout_cfg,
        out_root=args.out_root,
        args=args,
        force=args.force,
    )
    print(
        (
            f"blind holdout finished: small_base_scale={holdout['small_base_scale']} lattice={holdout['lattice']} "
            f"volume={holdout['volume']:.0f}"
        ),
        flush=True,
    )

    methods: dict[str, Any] = {}
    for method_cfg in METHODS:
        method_lattices = twisted_lattices if method_cfg["name"] == "twisted" else untwisted_lattices
        series = collect_method_series(
            method_cfg=method_cfg,
            lattices=method_lattices,
            boundary_index=boundary_index,
            mc_cfg=train_cfg,
            out_root=args.out_root,
            args=args,
            force=args.force,
        )
        methods[method_cfg["name"]] = fit_method_series(series, holdout)

    summary = {
        "benchmark": "geometry_acute_456_literal_sizes_blind_holdout",
        "description": (
            "Acute 4-5-6 shared-boundary-midpoint blind holdout study. The untwisted training family uses literal square sizes, "
            "the twisted training family uses explicit multiples of the small 4-5-6 base [6,6,3,1], and the blind datum is the smallest "
            "twisted multiple with boundary midpoint on-lattice and volume above the requested minimum."
        ),
        "shared_point": shared_point,
        "untwisted_sizes": list(untwisted_sizes),
        "twisted_small_base": list(SMALL_BASE_TWISTED),
        "twisted_small_base_scales": list(twisted_small_base_scales),
        "holdout_min_volume": int(args.holdout_min_volume),
        "beta_method": args.beta_method,
        "training_lattices": {
            "twisted": [list(lattice) for lattice in twisted_lattices],
            "untwisted": [list(lattice) for lattice in untwisted_lattices],
        },
        "train_mc": train_cfg,
        "holdout_mc": holdout_cfg,
        "beta_finder_cfg": beta_finder_cfg(float(TWISTED_CFG["beta"]), args) if args.beta_method == "beta_finder" else {"mode": "fixed"},
        "holdout": holdout,
        "methods": methods,
    }

    plot_path = args.out_root / "acute456_literal_sizes_blind_holdout_selected_models_logx.png"
    twisted_scan_plot_path = args.out_root / "acute456_literal_sizes_blind_holdout_twisted_susceptibilities.png"
    untwisted_scan_plot_path = args.out_root / "acute456_literal_sizes_blind_holdout_untwisted_susceptibilities.png"
    beta_pc_plot_path = args.out_root / "acute456_literal_sizes_blind_holdout_beta_pc_continuum.png"
    summary_path = args.out_root / "acute456_literal_sizes_blind_holdout_selected_models_logx.json"

    twisted_scan_entries = [
        {
            "title": item["label"],
            "volume": item["volume"],
            "beta_summary": item["beta_summary"],
        }
        for item in methods["twisted"]["training_series"]
    ]
    twisted_scan_entries.append(
        {
            "title": f"holdout {holdout['label']}",
            "volume": holdout["volume"],
            "beta_summary": holdout["beta_summary"],
        }
    )
    untwisted_scan_entries = [
        {
            "title": item["label"],
            "volume": item["volume"],
            "beta_summary": item["beta_summary"],
        }
        for item in methods["untwisted"]["training_series"]
    ]

    build_susceptibility_grid(
        twisted_scan_entries,
        title="Acute 4-5-6 twisted susceptibility scans",
        plot_path=twisted_scan_plot_path,
    )
    build_susceptibility_grid(
        untwisted_scan_entries,
        title="Acute 4-5-6 untwisted susceptibility scans",
        plot_path=untwisted_scan_plot_path,
    )
    summary["beta_pc_fits"] = build_beta_pc_continuum_plot(summary, beta_pc_plot_path)
    summary["plot_paths"] = {
        "holdout_fit": str(plot_path),
        "twisted_susceptibilities": str(twisted_scan_plot_path),
        "untwisted_susceptibilities": str(untwisted_scan_plot_path),
        "beta_pc_continuum": str(beta_pc_plot_path),
    }

    build_plot(summary, plot_path)
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(strip_fit_objects(summary), handle, indent=2)

    print(plot_path)
    print(twisted_scan_plot_path)
    print(untwisted_scan_plot_path)
    print(beta_pc_plot_path)
    print(summary_path)
    print(
        f"shared_point twisted=({shared_point['twisted']['m']},{shared_point['twisted']['n']}) "
        f"untwisted=({shared_point['untwisted']['m']},{shared_point['untwisted']['n']})"
    )
    print(
        f"holdout small_base_scale={holdout['small_base_scale']} lattice={holdout['lattice']} volume={holdout['volume']:.0f} "
        f"conn={holdout['conn']:.9f} err={holdout['conn_err']:.9f}"
    )
    for method_name in ("untwisted", "twisted"):
        print(f"\n{method_name}")
        beta_fit = summary["beta_pc_fits"][method_name]
        print(
            f"  beta_pc continuum: beta_inf={beta_fit['beta_inf']:.9f} "
            f"exact={beta_fit['exact_beta']:.9f} chi2/dof={beta_fit['chi2dof_train']:.3f}"
        )
        for model_name in methods[method_name]["ranking"]:
            row = methods[method_name]["fits"][model_name]
            print(
                f"  {model_name:10s} pred={row['prediction_holdout']:.9f} "
                f"z={row['z_holdout']:+.3f} chi2/dof={row['chi2dof_train']:.3f} "
                f"{format_params(model_name, row['coef'])}"
            )


if __name__ == "__main__":
    main()