#!/usr/bin/env python3
"""Spyder-friendly launcher for the local L128 continuum workflow bundle.

Open this file in Spyder, adjust the settings below, and press Run.
The script uses the *current* Python interpreter (the one Spyder is attached to).

Modes:
  - "reference" : generate the reference family only
  - "tests"     : generate the 7x7 test grid only
  - "score"     : fit test vs reference and make the comparison plots
  - "plots"     : same as "score"; convenient alias for plot regeneration
  - "all"       : run reference, tests, then score in sequence
"""

from __future__ import annotations

import json
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

MODE = "all"
AUTO_INSTALL_PYTHON_DEPS = True
N_WORKERS_OVERRIDE = None
SHOW_COMPARISON_PLOTS = False
MONITOR_INTERVAL_SECONDS = 300

ROOT = Path(__file__).resolve().parent
WORKDIR = ROOT / "K_from_continuum"
CONFIG_DIR = WORKDIR / "configs"
RESULTS_DIR = WORKDIR / "results"
PYTHON = sys.executable
SRC_PATH = WORKDIR / "src" / "ising_tri_twisted_parallelogram.cc"
INCLUDE_DIR = WORKDIR / "include"
BIN_DIR = WORKDIR / "bin"
EXE_PATH = BIN_DIR / ("ising_tri_twisted_parallelogram.exe" if os.name == "nt" else "ising_tri_twisted_parallelogram")

REFERENCE_CONFIG = CONFIG_DIR / "local_L128_reference.json"
TESTS_CONFIG = CONFIG_DIR / "local_L128_tests.json"
SCORE_CONFIG = CONFIG_DIR / "local_L128_score.json"

SCORE_OUTPUT_DIR = (
    RESULTS_DIR
    / "local_iso111_ref4x_qtr_score_L128_grid7x7"
    / "comparison_analysis_data"
)

VALID_MODES = {"reference", "tests", "score", "plots", "all"}
REQUIRED_MODULES = ("numpy", "scipy", "matplotlib")

RUN_SIM_START_RE = re.compile(
    r"run_simulator start "
    r"L=\((?P<Lx>-?\d+),(?P<Ly>-?\d+)\) T=\((?P<Tx>-?\d+),(?P<Ty>-?\d+)\) "
    r"k=\((?P<k1>[^,]+),(?P<k2>[^,]+),(?P<k3>[^)]+)\) beta=(?P<beta>\S+) "
    r"n_traj=(?P<n_traj>\d+) n_skip=(?P<n_skip>\d+) n_therm=(?P<n_therm>\d+)"
)
FIND_BETA_START_RE = re.compile(
    r"find_beta_c start "
    r"L=\((?P<Lx>-?\d+),(?P<Ly>-?\d+)\) T=\((?P<Tx>-?\d+),(?P<Ty>-?\d+)\) "
    r"k=\((?P<k1>[^,]+),(?P<k2>[^,]+),(?P<k3>[^)]+)\) "
    r"beta=\[(?P<beta_lo>[^,]+),(?P<beta_hi>[^\]]+)\]"
)
FIND_BETA_DONE_RE = re.compile(
    r"find_beta_c done .*beta_c=(?P<beta_c>\S+).*traj_total=(?P<traj_total>\d+)"
)
SCAN_PROGRESS_RE = re.compile(
    r"beta scan pass=(?P<pass_num>\d+) pt=(?P<pt_done>\d+)/(?P<pt_total>\d+) "
    r"beta=(?P<beta>\S+) chi=(?P<chi>\S+) \+/- (?P<chi_err>\S+) "
    r"traj_done=(?P<traj_done>\d+)"
)
POINT_PROGRESS_RE = re.compile(
    r"Point r1=(?P<r1>\S+), r2=(?P<r2>\S+):"
)


def format_seconds(seconds: float) -> str:
    total = max(0, int(round(seconds)))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def timestamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{timestamp()}] {message}", flush=True)


def _format_status_snapshot(status: dict[str, object]) -> str:
    parts: list[str] = []
    workflow_stage = status.get("workflow_stage")
    if workflow_stage:
        parts.append(f"workflow={workflow_stage}")
    step_name = status.get("step_name")
    if step_name:
        parts.append(f"step={step_name}")
    mc_phase = status.get("mc_phase")
    if mc_phase:
        parts.append(f"phase={mc_phase}")
    lattice = status.get("lattice")
    if lattice:
        parts.append(f"lattice={lattice}")
    couplings = status.get("couplings")
    if couplings:
        parts.append(f"k={couplings}")
    beta = status.get("beta")
    if beta:
        parts.append(f"beta={beta}")
    beta_window = status.get("beta_window")
    if beta_window:
        parts.append(f"beta_window={beta_window}")
    scan_progress = status.get("scan_progress")
    if scan_progress:
        parts.append(f"scan={scan_progress}")
    traj_done = status.get("traj_done")
    if traj_done is not None:
        parts.append(f"traj_done={traj_done}")
    current_n_traj = status.get("current_n_traj")
    if current_n_traj is not None:
        parts.append(f"run_n_traj={current_n_traj}")
    score_progress = status.get("score_progress")
    if score_progress:
        parts.append(f"score={score_progress}")
    if not parts:
        return "no parsed workflow details yet"
    return " | ".join(parts)


def _update_status_from_line(status: dict[str, object], line: str, label: str) -> str | None:
    stripped = line.strip()
    if not stripped:
        return None

    if label == "01_generate_reference.py":
        status["workflow_stage"] = "reference generation"
    elif label == "02_generate_tests.py":
        status["workflow_stage"] = "test generation"
    elif label == "03_score_continuum_vs_reference.py":
        status["workflow_stage"] = "scoring"

    if label == "03_score_continuum_vs_reference.py":
        status["mc_phase"] = "scoring"
        point_match = POINT_PROGRESS_RE.search(stripped)
        if point_match:
            status["score_progress"] = (
                f"processing point r1={point_match.group('r1')}, r2={point_match.group('r2')}"
            )
            return _format_status_snapshot(status)
        if "comparison_analysis_data" in stripped or "score_heatmap" in stripped or "zscore_heatmap" in stripped:
            status["score_progress"] = "writing score outputs"
            return _format_status_snapshot(status)
        return None

    match = FIND_BETA_START_RE.search(stripped)
    if match:
        status["step_name"] = label
        status["mc_phase"] = "beta MC"
        status["lattice"] = f"({match.group('Lx')},{match.group('Ly')}) T=({match.group('Tx')},{match.group('Ty')})"
        status["couplings"] = f"({match.group('k1')},{match.group('k2')},{match.group('k3')})"
        status["beta_window"] = f"[{match.group('beta_lo')},{match.group('beta_hi')}]"
        status["beta"] = None
        status["scan_progress"] = None
        status["traj_done"] = 0
        status["current_n_traj"] = None
        return _format_status_snapshot(status)

    match = SCAN_PROGRESS_RE.search(stripped)
    if match:
        status["mc_phase"] = "beta MC"
        status["beta"] = match.group("beta")
        status["scan_progress"] = f"pass {match.group('pass_num')} point {match.group('pt_done')}/{match.group('pt_total')}"
        status["traj_done"] = int(match.group("traj_done"))
        return _format_status_snapshot(status)

    match = FIND_BETA_DONE_RE.search(stripped)
    if match:
        status["mc_phase"] = "beta MC complete"
        status["beta"] = match.group("beta_c")
        status["traj_done"] = int(match.group("traj_total"))
        status["scan_progress"] = "done"
        return _format_status_snapshot(status)

    match = RUN_SIM_START_RE.search(stripped)
    if match:
        current_n_traj = int(match.group("n_traj"))
        status["step_name"] = label
        status["lattice"] = f"({match.group('Lx')},{match.group('Ly')}) T=({match.group('Tx')},{match.group('Ty')})"
        status["couplings"] = f"({match.group('k1')},{match.group('k2')},{match.group('k3')})"
        status["beta"] = match.group("beta")
        status["current_n_traj"] = current_n_traj
        if status.get("mc_phase") not in {"beta MC", "beta MC complete"}:
            status["mc_phase"] = "production MC"
        elif status.get("scan_progress") == "done":
            status["mc_phase"] = "production MC"
        if status.get("mc_phase") == "production MC":
            status["scan_progress"] = None
            status["traj_done"] = None
        return _format_status_snapshot(status)

    if stripped.startswith("[scan]"):
        scan_progress = stripped.replace("[scan]", "").strip()
        status["scan_progress"] = scan_progress
        if status.get("mc_phase") in {"beta MC", "beta MC complete"}:
            status["mc_phase"] = "beta MC"
        return _format_status_snapshot(status)

    if "production run" in stripped.lower() or "prod_retry" in stripped.lower():
        status["mc_phase"] = "production MC"
        return _format_status_snapshot(status)

    return None


def run_monitored_subprocess(cmd: list[str], *, cwd: Path, label: str) -> None:
    log(f"Starting {label}")
    log(f"Command: {' '.join(cmd)}")
    start = time.perf_counter()
    status: dict[str, object] = {
        "step_name": label,
        "workflow_stage": (
            "reference generation" if label == "01_generate_reference.py"
            else "test generation" if label == "02_generate_tests.py"
            else "scoring" if label == "03_score_continuum_vs_reference.py"
            else None
        ),
        "mc_phase": "scoring" if label == "03_score_continuum_vs_reference.py" else None,
        "lattice": None,
        "couplings": None,
        "beta": None,
        "beta_window": None,
        "scan_progress": None,
        "traj_done": None,
        "current_n_traj": None,
        "score_progress": "loading scoring inputs" if label == "03_score_continuum_vs_reference.py" else None,
    }
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    line_queue: queue.Queue[str] = queue.Queue()

    def _reader() -> None:
        assert proc.stdout is not None
        for raw_line in proc.stdout:
            line_queue.put(raw_line)
        proc.stdout.close()

    reader = threading.Thread(target=_reader, daemon=True)
    reader.start()

    last_status_line: str | None = None
    try:
        while True:
            try:
                raw_line = line_queue.get(timeout=MONITOR_INTERVAL_SECONDS)
                print(raw_line, end="", flush=True)
                status_line = _update_status_from_line(status, raw_line, label)
                if status_line and status_line != last_status_line:
                    log(f"Status | {status_line}")
                    last_status_line = status_line
            except queue.Empty:
                if proc.poll() is not None:
                    break
                elapsed = time.perf_counter() - start
                log(
                    f"{label} still running | elapsed {format_seconds(elapsed)} | "
                    f"{_format_status_snapshot(status)}"
                )
    except KeyboardInterrupt:
        proc.terminate()
        raise

    while not line_queue.empty():
        raw_line = line_queue.get_nowait()
        print(raw_line, end="", flush=True)
        status_line = _update_status_from_line(status, raw_line, label)
        if status_line and status_line != last_status_line:
            log(f"Status | {status_line}")
            last_status_line = status_line

    reader.join(timeout=1.0)
    returncode = proc.wait()

    elapsed = time.perf_counter() - start
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)
    log(f"Finished {label} | elapsed {format_seconds(elapsed)}")


def require_command(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Missing required command: {name}")


def find_cpp_compiler() -> tuple[str, str] | None:
    for name in ("g++", "c++", "cl"):
        path = shutil.which(name)
        if path is not None:
            return name, path
    return None


def ensure_local_tools() -> None:
    compiler = find_cpp_compiler()
    if compiler is None:
        raise RuntimeError(
            "Missing C++ build tool. Install one of: g++, c++, or cl "
            "(Visual Studio Build Tools). 'make' is optional."
        )
    if shutil.which("make") is None:
        print(f"[setup] 'make' not found; will compile directly with {compiler[0]}")
    else:
        print(f"[setup] Found make and compiler {compiler[0]}")


def missing_python_modules() -> list[str]:
    missing = []
    for name in REQUIRED_MODULES:
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    return missing


def ensure_python_deps() -> None:
    missing = missing_python_modules()
    if not missing:
        log("Python dependencies already available")
        return
    if not AUTO_INSTALL_PYTHON_DEPS:
        raise RuntimeError(
            "Missing Python packages in the current interpreter: "
            f"{', '.join(missing)}. Install them manually or set "
            "AUTO_INSTALL_PYTHON_DEPS = True."
        )
    requirements = WORKDIR / "requirements.txt"
    log(f"Installing Python dependencies from {requirements}")
    try:
        run_monitored_subprocess(
            [PYTHON, "-m", "pip", "install", "-r", str(requirements)],
            cwd=ROOT,
            label="pip install",
        )
    except subprocess.CalledProcessError as exc:
        still_missing = missing_python_modules()
        if not still_missing:
            log("pip install failed, but required Python modules are available; continuing")
            return
        raise RuntimeError(
            "Automatic pip install failed. Install these packages into the "
            f"Spyder interpreter manually: {', '.join(still_missing)}"
        ) from exc


def _input_paths() -> list[Path]:
    return [SRC_PATH, *sorted(INCLUDE_DIR.glob("*.h"))]


def _needs_rebuild(output_path: Path) -> bool:
    if not output_path.exists():
        return True
    out_mtime = output_path.stat().st_mtime
    return any(path.stat().st_mtime > out_mtime for path in _input_paths())


def _build_with_make() -> None:
    run_monitored_subprocess(["make"], cwd=WORKDIR, label="simulator build (make)")


def _build_with_compiler() -> None:
    compiler = find_cpp_compiler()
    if compiler is None:
        raise RuntimeError("No supported C++ compiler found")

    name, path = compiler
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    if not _needs_rebuild(EXE_PATH):
        log(f"Simulator already up to date: {EXE_PATH}")
        return

    if name in {"g++", "c++"}:
        cmd = [
            path,
            "-std=c++14",
            "-O3",
            "-Wall",
            "-Wno-sign-compare",
            "-Wno-deprecated-declarations",
            "-I",
            str(INCLUDE_DIR),
            str(SRC_PATH),
            "-o",
            str(EXE_PATH),
        ]
    elif name == "cl":
        cmd = [
            path,
            "/nologo",
            "/O2",
            "/EHsc",
            "/std:c++14",
            f"/I{INCLUDE_DIR}",
            str(SRC_PATH),
            f"/Fe:{EXE_PATH}",
        ]
    else:
        raise RuntimeError(f"Unsupported compiler: {name}")

    run_monitored_subprocess(cmd, cwd=WORKDIR, label=f"simulator build ({name})")


def build_simulator() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if shutil.which("make") is not None:
        _build_with_make()
    else:
        _build_with_compiler()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_temp_config(base_path: Path, overrides: dict) -> Path:
    data = load_json(base_path)
    for key, value in overrides.items():
        cursor = data
        parts = key.split(".")
        for part in parts[:-1]:
            cursor = cursor.setdefault(part, {})
        cursor[parts[-1]] = value

    temp_dir = Path(tempfile.mkdtemp(prefix="local_l128_cfg_"))
    temp_path = temp_dir / base_path.name
    with temp_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
    return temp_path


def run_step(script_name: str, config_path: Path, *, overrides: dict | None = None) -> None:
    actual_config = config_path
    if overrides:
        actual_config = write_temp_config(config_path, overrides)
        log(f"Using temporary config override: {actual_config}")

    cmd = [PYTHON, "-u", script_name, "--config", str(actual_config)]
    run_monitored_subprocess(cmd, cwd=WORKDIR, label=script_name)


def maybe_show_comparison_plots() -> None:
    if not SHOW_COMPARISON_PLOTS:
        return
    main_paths = [
        SCORE_OUTPUT_DIR / "score_heatmap.png",
        SCORE_OUTPUT_DIR / "zscore_heatmap.png",
        SCORE_OUTPUT_DIR / "score_zscore_heatmaps.png",
    ]
    if not all(path.exists() for path in main_paths):
        print("[plots] Comparison heatmaps not found yet; skipping display")
        return

    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    for path in main_paths:
        image = mpimg.imread(path)
        plt.figure(figsize=(8, 6))
        plt.imshow(image)
        plt.axis("off")
        plt.title(path.name)
    plt.show()


def print_plot_paths() -> None:
    log("Main comparison outputs:")
    print(f"  {SCORE_OUTPUT_DIR / 'score_heatmap.png'}")
    print(f"  {SCORE_OUTPUT_DIR / 'zscore_heatmap.png'}")
    print(f"  {SCORE_OUTPUT_DIR / 'score_zscore_heatmaps.png'}")
    print(f"  {SCORE_OUTPUT_DIR / 'fss_plots'}")


def main() -> None:
    mode = MODE.strip().lower()
    if mode not in VALID_MODES:
        raise ValueError(f"MODE must be one of {sorted(VALID_MODES)}, got {MODE!r}")

    total_start = time.perf_counter()
    log(f"Bundle root: {ROOT}")
    log(f"Python interpreter: {PYTHON}")
    log(f"Selected mode: {mode}")
    ensure_local_tools()
    ensure_python_deps()
    build_simulator()

    worker_overrides = None
    if N_WORKERS_OVERRIDE is not None:
        worker_overrides = {"execution.n_workers": int(N_WORKERS_OVERRIDE)}

    if mode in {"reference", "all"}:
        run_step("01_generate_reference.py", REFERENCE_CONFIG, overrides=worker_overrides)

    if mode in {"tests", "all"}:
        run_step("02_generate_tests.py", TESTS_CONFIG, overrides=worker_overrides)

    if mode in {"score", "plots", "all"}:
        run_step("03_score_continuum_vs_reference.py", SCORE_CONFIG)
        print_plot_paths()
        maybe_show_comparison_plots()

    log(f"Workflow finished | total elapsed {format_seconds(time.perf_counter() - total_start)}")


if __name__ == "__main__":
    main()
