#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
RESPONSIBLE_ROOT = HERE.parent
REPO_ROOT = RESPONSIBLE_ROOT.parent
KFC_ROOT = REPO_ROOT / "K_from_continuum"
if str(KFC_ROOT) not in sys.path:
    sys.path.insert(0, str(KFC_ROOT))

from workflow_common import ensure_dir, ensure_simulator, exact_triangular_ising_beta  # noqa: E402


DEFAULT_EXECUTION = {
    "exe": None,
    "auto_build": True,
    "force_rebuild": False,
    "build_command": ["make"],
    "build_timeout_s": 600,
}


@dataclass(frozen=True)
class CouplingPoint:
    label: str
    r1: float
    r2: float


def _resolve_sizes(values: list[int] | None) -> list[int]:
    if not values:
        return []
    resolved: list[int] = []
    seen: set[int] = set()
    for raw in values:
        size = int(raw)
        if size <= 0:
            raise ValueError(f"lattice sizes must be positive, got {size}")
        if size in seen:
            continue
        seen.add(size)
        resolved.append(size)
    return resolved


def _resolve_displacement(lattice: tuple[int, int, int, int], args: argparse.Namespace) -> tuple[int, int]:
    if args.disp_fraction is None:
        return (int(args.disp_m), int(args.disp_n))
    raw_m = float(lattice[0]) * float(args.disp_fraction)
    disp_m = int(round(raw_m))
    if not math.isclose(raw_m, float(disp_m), rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(
            f"disp-fraction={args.disp_fraction} is not exactly representable for Lx={lattice[0]}"
        )
    return (disp_m, int(args.disp_n))


def _run_id(lx: int, ly: int, tx: int, ty: int, k1: float, k2: float, k3: float, kt: float = 0.5) -> str:
    return (
        f"{lx}x{ly}_t{tx}x{ty}_"
        f"k{int(round(k1 * 1000))}_{int(round(k2 * 1000))}_{int(round(k3 * 1000))}_{int(round(kt * 1000))}"
    )


def _point_beta(point: CouplingPoint) -> float:
    return float(exact_triangular_ising_beta({"k1": point.r1, "k2": point.r2, "k3": 1.0}))


def _parse_sample_file(path: Path) -> dict[str, Any]:
    header: dict[str, Any] = {}
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split()
                if len(parts) >= 2:
                    key = parts[0]
                    if key in {"seed", "L_x", "L_y", "T_x", "T_y", "n_therm", "n_traj", "n_skip", "d", "m", "n"}:
                        header[key] = int(parts[1])
                    elif key in {"K1", "K2", "K3", "Kt", "beta"}:
                        header[key] = float(parts[1])
                continue
            cols = line.split()
            if len(cols) < 6:
                raise ValueError(f"expected 6 sample columns in {path}, got: {line}")
            rows.append([float(cols[1]), float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5])])
    if not rows:
        raise ValueError(f"no samples found in {path}")
    arr = np.asarray(rows, dtype=float)
    return {
        "path": str(path),
        "corr": arr[:, 0],
        "mag": arr[:, 1],
        "e1": arr[:, 2],
        "e2": arr[:, 3],
        "e3": arr[:, 4],
        "beta": float(header["beta"]),
        "K": (float(header["K1"]), float(header["K2"]), float(header["K3"])),
        "n_sites": int(header["L_x"]) * int(header["L_y"]),
        "header": header,
    }


class CorrelatorReweighter:
    def __init__(self, payload: dict[str, Any], *, n_blocks: int = 16, n_eff_floor: float = 0.1):
        self._corr = np.asarray(payload["corr"], dtype=float)
        self._mag = np.asarray(payload["mag"], dtype=float)
        self._e1 = np.asarray(payload["e1"], dtype=float)
        self._e2 = np.asarray(payload["e2"], dtype=float)
        self._e3 = np.asarray(payload["e3"], dtype=float)
        self._beta_parent = float(payload["beta"])
        self._k_parent = tuple(float(value) for value in payload["K"])
        self._n_sites = int(payload["n_sites"])
        self._n = int(self._corr.size)
        if self._n < 16:
            raise ValueError("need at least 16 measurements for block jackknife")
        self._n_blocks = max(2, min(int(n_blocks), self._n // 4))
        self._n_eff_floor = float(n_eff_floor)

    def _log_weights(self, beta_target: float, k_target: tuple[float, float, float]) -> np.ndarray:
        parent = self._k_parent[0] * self._e1 + self._k_parent[1] * self._e2 + self._k_parent[2] * self._e3
        target = k_target[0] * self._e1 + k_target[1] * self._e2 + k_target[2] * self._e3
        return -self._n_sites * (float(beta_target) * target - self._beta_parent * parent)

    @staticmethod
    def _normalize(log_w: np.ndarray) -> np.ndarray:
        shifted = log_w - np.max(log_w)
        w = np.exp(shifted)
        total = np.sum(w)
        if not np.isfinite(total) or total <= 0.0:
            return np.zeros_like(w)
        return w / total

    @staticmethod
    def _connected_from_weights(w: np.ndarray, corr: np.ndarray, mag: np.ndarray) -> float:
        corr_mean = float(np.sum(w * corr))
        mag_mean = float(np.sum(w * mag))
        return corr_mean - mag_mean * mag_mean

    @staticmethod
    def _n_eff(w: np.ndarray) -> float:
        total = float(np.sum(w))
        if total <= 0.0:
            return 0.0
        return total * total / max(float(np.dot(w, w)), 1.0e-300)

    def connected(self, beta_target: float, k_target: tuple[float, float, float]) -> dict[str, float]:
        log_w = self._log_weights(beta_target, k_target)
        w = self._normalize(log_w)
        n_eff = self._n_eff(w * self._n)
        connected = self._connected_from_weights(w, self._corr, self._mag)
        if n_eff < self._n_eff_floor * self._n:
            return {"connected": connected, "sigma": float("nan"), "n_eff": n_eff}

        block = self._n // self._n_blocks
        leave_one_out: list[float] = []
        for index in range(self._n_blocks):
            lo = index * block
            hi = (index + 1) * block if index < self._n_blocks - 1 else self._n
            mask = np.ones(self._n, dtype=bool)
            mask[lo:hi] = False
            if np.count_nonzero(mask) < 8:
                continue
            w_sub = self._normalize(log_w[mask])
            leave_one_out.append(self._connected_from_weights(w_sub, self._corr[mask], self._mag[mask]))
        if len(leave_one_out) < 2:
            sigma = float("nan")
        else:
            jk = np.asarray(leave_one_out, dtype=float)
            sigma = float(np.sqrt((jk.size - 1.0) / jk.size * np.sum((jk - np.mean(jk)) ** 2)))
        return {"connected": connected, "sigma": sigma, "n_eff": n_eff}


def _direct_connected(payload: dict[str, Any]) -> dict[str, float]:
    corr = np.asarray(payload["corr"], dtype=float)
    mag = np.asarray(payload["mag"], dtype=float)
    connected = float(np.mean(corr) - np.mean(mag) ** 2)
    n = int(corr.size)
    if n <= 1:
        return {"connected": connected, "sigma": 0.0}
    jk = np.empty(n, dtype=float)
    sum_corr = float(np.sum(corr))
    sum_mag = float(np.sum(mag))
    for index in range(n):
        corr_leave = (sum_corr - float(corr[index])) / (n - 1)
        mag_leave = (sum_mag - float(mag[index])) / (n - 1)
        jk[index] = corr_leave - mag_leave * mag_leave
    sigma = float(np.sqrt((n - 1.0) / n * np.sum((jk - np.mean(jk)) ** 2)))
    return {"connected": connected, "sigma": sigma}


def _run_point(
    exe: str,
    output_root: Path,
    lattice: tuple[int, int, int, int],
    disp: tuple[int, int],
    point: CouplingPoint,
    n_traj: int,
    n_therm: int,
    n_skip: int,
    seed: int,
) -> dict[str, Any]:
    lx, ly, tx, ty = lattice
    beta = _point_beta(point)
    point_dir = output_root / point.label
    ensure_dir(str(point_dir))
    sample_name = "single_disp_reweight_samples.dat"
    cmd = [
        exe,
        "--L_x",
        str(lx),
        "--L_y",
        str(ly),
        "--T_x",
        str(tx),
        "--T_y",
        str(ty),
        "--k1",
        f"{point.r1:.12f}",
        "--k2",
        f"{point.r2:.12f}",
        "--k3",
        "1.000000000000",
        "--beta",
        f"{beta:.12f}",
        "--n_traj",
        str(int(n_traj)),
        "--n_therm",
        str(int(n_therm)),
        "--n_skip",
        str(int(n_skip)),
        "--seed",
        str(int(seed)),
        "--data_dir",
        str(point_dir),
        "--single_disp_m",
        str(int(disp[0])),
        "--single_disp_n",
        str(int(disp[1])),
        "--single_disp_samples_name",
        sample_name,
    ]
    start_t = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"simulator failed for {point.label}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
    run_id = _run_id(lx, ly, tx, ty, point.r1, point.r2, 1.0)
    sample_path = point_dir / run_id / sample_name
    if not sample_path.is_file():
        raise FileNotFoundError(f"missing sample file for {point.label}: {sample_path}")
    payload = _parse_sample_file(sample_path)
    payload["wall_seconds"] = float(time.time() - start_t)
    payload["label"] = point.label
    payload["r1"] = float(point.r1)
    payload["r2"] = float(point.r2)
    payload["beta"] = float(beta)
    return payload


def _build_direct_rows(payloads: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for payload in payloads:
        direct = _direct_connected(payload)
        rows.append(
            {
                "label": payload["label"],
                "r1": payload["r1"],
                "r2": payload["r2"],
                "beta": payload["beta"],
                "wall_seconds": payload["wall_seconds"],
                "connected": direct["connected"],
                "sigma": direct["sigma"],
                "n_samples": int(len(payload["corr"])),
            }
        )
    return rows


def _build_reweight_rows(payloads: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    direct_lookup = {str(payload["label"]): _direct_connected(payload) for payload in payloads}
    for source in payloads:
        rw = CorrelatorReweighter(source)
        for target in payloads:
            predicted = rw.connected(float(target["beta"]), (float(target["r1"]), float(target["r2"]), 1.0))
            target_direct = direct_lookup[str(target["label"])]
            sigma_combined = math.sqrt(
                max(float(predicted["sigma"]) ** 2 if math.isfinite(float(predicted["sigma"])) else 0.0, 0.0)
                + float(target_direct["sigma"]) ** 2
            )
            delta = float(predicted["connected"]) - float(target_direct["connected"])
            z = delta / sigma_combined if sigma_combined > 0.0 else float("nan")
            rows.append(
                {
                    "source": source["label"],
                    "target": target["label"],
                    "predicted_connected": float(predicted["connected"]),
                    "predicted_sigma": float(predicted["sigma"]),
                    "target_connected": float(target_direct["connected"]),
                    "target_sigma": float(target_direct["sigma"]),
                    "delta": delta,
                    "z": z,
                    "n_eff": float(predicted["n_eff"]),
                    "n_eff_fraction": float(predicted["n_eff"]) / max(len(source["corr"]), 1),
                }
            )
    return rows


def _write_single_summary(
    output_root: Path,
    lattice: tuple[int, int, int, int],
    disp: tuple[int, int],
    args: argparse.Namespace,
    direct_rows: list[dict[str, Any]],
    reweight_rows: list[dict[str, Any]],
) -> None:
    summary = {
        "lattice": {"Lx": lattice[0], "Ly": lattice[1], "Tx": lattice[2], "Ty": lattice[3]},
        "displacement": {"m": disp[0], "n": disp[1]},
        "n_traj": int(args.n_traj),
        "n_therm": int(args.n_therm),
        "n_skip": int(args.n_skip),
        "direct": direct_rows,
        "reweight": reweight_rows,
    }
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    tsv_lines = [
        "source\ttarget\tpredicted_connected\tpredicted_sigma\ttarget_connected\ttarget_sigma\tdelta\tz\tn_eff\tn_eff_fraction"
    ]
    for row in reweight_rows:
        tsv_lines.append(
            "\t".join(
                [
                    str(row["source"]),
                    str(row["target"]),
                    f"{row['predicted_connected']:.16e}",
                    f"{row['predicted_sigma']:.16e}",
                    f"{row['target_connected']:.16e}",
                    f"{row['target_sigma']:.16e}",
                    f"{row['delta']:.16e}",
                    f"{row['z']:.16e}",
                    f"{row['n_eff']:.16e}",
                    f"{row['n_eff_fraction']:.16e}",
                ]
            )
        )
    (output_root / "reweight_vs_direct.tsv").write_text("\n".join(tsv_lines) + "\n", encoding="utf-8")


def _write_multi_size_summary(
    output_root: Path,
    args: argparse.Namespace,
    points: list[CouplingPoint],
    runs: list[dict[str, Any]],
) -> None:
    summary = {
        "sizes": [int(run["size"]) for run in runs],
        "disp_fraction": (float(args.disp_fraction) if args.disp_fraction is not None else None),
        "n_traj": int(args.n_traj),
        "n_therm": int(args.n_therm),
        "n_skip": int(args.n_skip),
        "points": [{"label": point.label, "r1": point.r1, "r2": point.r2} for point in points],
        "runs": runs,
    }
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    direct_lines = [
        "size\tLx\tLy\tTx\tTy\tdisp_m\tdisp_n\tlabel\tr1\tr2\tbeta\tconnected\tsigma\tn_samples\twall_seconds"
    ]
    reweight_lines = [
        "size\tLx\tLy\tTx\tTy\tdisp_m\tdisp_n\tsource\ttarget\tpredicted_connected\tpredicted_sigma\ttarget_connected\ttarget_sigma\tdelta\tz\tn_eff\tn_eff_fraction"
    ]
    for run in runs:
        lattice = run["lattice"]
        disp = run["displacement"]
        for row in run["direct"]:
            direct_lines.append(
                "\t".join(
                    [
                        str(run["size"]),
                        str(lattice["Lx"]),
                        str(lattice["Ly"]),
                        str(lattice["Tx"]),
                        str(lattice["Ty"]),
                        str(disp["m"]),
                        str(disp["n"]),
                        str(row["label"]),
                        f"{row['r1']:.6f}",
                        f"{row['r2']:.6f}",
                        f"{row['beta']:.16e}",
                        f"{row['connected']:.16e}",
                        f"{row['sigma']:.16e}",
                        str(row["n_samples"]),
                        f"{row['wall_seconds']:.6f}",
                    ]
                )
            )
        for row in run["reweight"]:
            reweight_lines.append(
                "\t".join(
                    [
                        str(run["size"]),
                        str(lattice["Lx"]),
                        str(lattice["Ly"]),
                        str(lattice["Tx"]),
                        str(lattice["Ty"]),
                        str(disp["m"]),
                        str(disp["n"]),
                        str(row["source"]),
                        str(row["target"]),
                        f"{row['predicted_connected']:.16e}",
                        f"{row['predicted_sigma']:.16e}",
                        f"{row['target_connected']:.16e}",
                        f"{row['target_sigma']:.16e}",
                        f"{row['delta']:.16e}",
                        f"{row['z']:.16e}",
                        f"{row['n_eff']:.16e}",
                        f"{row['n_eff_fraction']:.16e}",
                    ]
                )
            )
    (output_root / "direct.tsv").write_text("\n".join(direct_lines) + "\n", encoding="utf-8")
    (output_root / "reweight_vs_direct.tsv").write_text("\n".join(reweight_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run three nearby coupling points and test single-displacement correlator reweighting between them."
    )
    parser.add_argument("--Lx", type=int, default=32)
    parser.add_argument("--Ly", type=int, default=32)
    parser.add_argument("--Tx", type=int, default=0)
    parser.add_argument("--Ty", type=int, default=0)
    parser.add_argument("--sizes", nargs="*", type=int, default=None)
    parser.add_argument("--disp-m", type=int, default=8)
    parser.add_argument("--disp-n", type=int, default=0)
    parser.add_argument(
        "--disp-fraction",
        type=float,
        default=None,
        help="If provided, use m = disp_fraction * Lx for each size; requires exact integer representability.",
    )
    parser.add_argument("--n-traj", type=int, default=100000)
    parser.add_argument("--n-therm", type=int, default=10000)
    parser.add_argument("--n-skip", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=2026060401)
    parser.add_argument(
        "--output-root",
        default=str(RESPONSIBLE_ROOT / "results" / "correlator_reweight_nearby_iso111_L32_hi100k_20260604"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    ensure_dir(str(output_root))
    exe = ensure_simulator(DEFAULT_EXECUTION)
    size_list = _resolve_sizes(args.sizes)
    points = [
        CouplingPoint("center", 1.000000, 1.000000),
        CouplingPoint("r1_plus_0p01", 1.010000, 1.000000),
        CouplingPoint("r2_plus_0p01", 1.000000, 1.010000),
    ]

    lattices = (
        [(int(size), int(size), int(args.Tx), int(args.Ty)) for size in size_list]
        if size_list
        else [(int(args.Lx), int(args.Ly), int(args.Tx), int(args.Ty))]
    )

    runs = []
    for lattice_index, lattice in enumerate(lattices):
        disp = _resolve_displacement(lattice, args)
        print(
            f"[run] lattice=({lattice[0]},{lattice[1]},{lattice[2]},{lattice[3]}) "
            f"disp=({disp[0]},{disp[1]}) n_traj={args.n_traj}"
        )
        payloads = []
        for point_index, point in enumerate(points):
            print(f"  [point] {point.label} r1={point.r1:.6f} r2={point.r2:.6f}")
            payloads.append(
                _run_point(
                    exe,
                    output_root,
                    lattice,
                    disp,
                    point,
                    int(args.n_traj),
                    int(args.n_therm),
                    int(args.n_skip),
                    int(args.seed_base) + lattice_index * 100 + point_index,
                )
            )

        direct_rows = _build_direct_rows(payloads)
        reweight_rows = _build_reweight_rows(payloads)
        runs.append(
            {
                "size": int(lattice[0]),
                "lattice": {"Lx": lattice[0], "Ly": lattice[1], "Tx": lattice[2], "Ty": lattice[3]},
                "displacement": {"m": disp[0], "n": disp[1]},
                "direct": direct_rows,
                "reweight": reweight_rows,
            }
        )

    if len(runs) == 1:
        run = runs[0]
        _write_single_summary(
            output_root,
            (run["lattice"]["Lx"], run["lattice"]["Ly"], run["lattice"]["Tx"], run["lattice"]["Ty"]),
            (run["displacement"]["m"], run["displacement"]["n"]),
            args,
            list(run["direct"]),
            list(run["reweight"]),
        )
        print(
            json.dumps(
                {
                    "lattice": run["lattice"],
                    "displacement": run["displacement"],
                    "n_traj": int(args.n_traj),
                    "n_therm": int(args.n_therm),
                    "n_skip": int(args.n_skip),
                    "direct": run["direct"],
                    "reweight": run["reweight"],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return

    _write_multi_size_summary(output_root, args, points, runs)
    print(
        json.dumps(
            {
                "sizes": [int(run["size"]) for run in runs],
                "disp_fraction": (float(args.disp_fraction) if args.disp_fraction is not None else None),
                "n_traj": int(args.n_traj),
                "n_therm": int(args.n_therm),
                "n_skip": int(args.n_skip),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
