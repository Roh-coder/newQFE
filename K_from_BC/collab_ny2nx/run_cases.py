#!/usr/bin/env python3
"""Run the curated Ny=2Nx comparison cases for collaborators.

This script compiles the K_from_BC executable and runs three high-stat cases:
1) Equilateral sinh-rule reference
2) 45-45-90 square-limit (k3=0)
3) Right triangle with height/base=1/2 (k3=0)

Outputs are stored under K_from_BC/collab_ny2nx/output/.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List


THETA_EQ = 1.0471975512


@dataclass
class Case:
    name: str
    title: str
    k1: float
    k2: float
    k3: float


CASES: List[Case] = [
    Case(
        name="equilateral",
        title="Equilateral sinh-rule",
        k1=0.27465307216702745,
        k2=0.27465307216702745,
        k3=0.27465307216702745,
    ),
    Case(
        name="right_45_45_90_k3zero",
        title="45-45-90 square-limit (k3=0)",
        k1=0.4406867935097715,
        k2=0.4406867935097715,
        k3=0.0,
    ),
    Case(
        name="right_h_over_b_half_k3zero",
        title="Right triangle h/b=1/2 (k3=0)",
        k1=0.24060591252980174,
        k2=0.7218177375894052,
        k3=0.0,
    ),
]


def sh(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=str(cwd), check=True)


def sh_capture(cmd: list[str], cwd: Path) -> str:
    cp = subprocess.run(cmd, cwd=str(cwd), check=True, text=True, capture_output=True)
    return cp.stdout


def format_dat_name(nx: int, ny: int, k1: float, k2: float, k3: float) -> str:
    return (
        f"{nx}_{ny}_t1_{THETA_EQ:.3f}_t2_{THETA_EQ:.3f}_"
        f"k_{k1:.3f}_{k2:.3f}_{k3:.3f}.dat"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=32)
    ap.add_argument("--ny", type=int, default=64)
    ap.add_argument("--n-therm", type=int, default=5000)
    ap.add_argument("--n-traj", type=int, default=120000)
    ap.add_argument("--n-skip", type=int, default=20)
    ap.add_argument("--n-wolff", type=int, default=3)
    ap.add_argument("--n-metropolis", type=int, default=5)
    ap.add_argument(
        "--source",
        default="K_from_BC/collab_ny2nx/ising_flat_crit_k_from_bc.cc",
        help="C++ source file to compile",
    )
    ap.add_argument("--binary", default="K_from_BC/collab_ny2nx/bin/ising_kbc")
    ap.add_argument("--skip-compile", action="store_true")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    out_dir = repo / "K_from_BC" / "collab_ny2nx" / "output"
    results_dir = out_dir / "results"
    log_dir = out_dir / "logs"
    data_dir = results_dir / "runs"
    bin_path = repo / args.binary
    src = repo / args.source

    log_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "plots").mkdir(parents=True, exist_ok=True)
    bin_path.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_compile:
        print(f"Compiling: {bin_path}")
        sh(
            [
                "g++",
                "-std=c++14",
                "-O3",
                "-Wall",
                "-Wno-sign-compare",
                "-I",
                "include",
                str(src),
                "-o",
                str(bin_path),
            ],
            cwd=repo,
        )

    manifest = {
        "nx": args.nx,
        "ny": args.ny,
        "theta1": THETA_EQ,
        "theta2": THETA_EQ,
        "n_therm": args.n_therm,
        "n_traj": args.n_traj,
        "n_skip": args.n_skip,
        "n_wolff": args.n_wolff,
        "n_metropolis": args.n_metropolis,
        "cases": [],
    }

    root_dat_dir = repo / "ising_flat_crit_k_from_bc"
    root_dat_dir.mkdir(exist_ok=True)

    for case in CASES:
        log_path = log_dir / f"{case.name}.log"
        print(f"Running case: {case.name}")
        run_cmd = [
            str(bin_path),
            "--n_x",
            str(args.nx),
            "--n_y",
            str(args.ny),
            "--theta1",
            str(THETA_EQ),
            "--theta2",
            str(THETA_EQ),
            "--k1",
            str(case.k1),
            "--k2",
            str(case.k2),
            "--k3",
            str(case.k3),
            "--n_therm",
            str(args.n_therm),
            "--n_traj",
            str(args.n_traj),
            "--n_skip",
            str(args.n_skip),
            "--n_wolff",
            str(args.n_wolff),
            "--n_metropolis",
            str(args.n_metropolis),
        ]
        stdout = sh_capture(run_cmd, cwd=repo)
        log_path.write_text(stdout)

        dat_name = format_dat_name(args.nx, args.ny, case.k1, case.k2, case.k3)
        src_dat = root_dat_dir / dat_name
        if not src_dat.exists():
            raise FileNotFoundError(f"Expected output not found: {src_dat}")

        dst_dat = data_dir / f"{case.name}.dat"
        shutil.copy2(src_dat, dst_dat)

        manifest["cases"].append(
            {
                **asdict(case),
                "source_dat": str(src_dat.relative_to(repo)),
                "copied_dat": str(dst_dat.relative_to(repo)),
                "log": str(log_path.relative_to(repo)),
            }
        )

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote manifest: {manifest_path.relative_to(repo)}")


if __name__ == "__main__":
    main()
