#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import hashlib
from datetime import datetime
from typing import Any

import numpy as np
from scipy.optimize import curve_fit

_WORKFLOW_ROOT = os.path.dirname(os.path.abspath(__file__))
_RUNTIME_ROOT = _WORKFLOW_ROOT
_LIB_ROOT = os.path.join(_RUNTIME_ROOT, "lib")
if _LIB_ROOT not in sys.path:
    sys.path.insert(0, _LIB_ROOT)

import mc_engine  # noqa: E402
from cost import _SQRT3_2, _tile_interp, boundary_paths  # noqa: E402


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{timestamp()}] {message}", flush=True)


def json_default(obj: Any):
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def resolve_path(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return os.path.normpath(path)
    return os.path.normpath(os.path.join(base_dir, path))


def _strip_json_comments(text: str) -> str:
    out: list[str] = []
    i = 0
    n = len(text)
    in_string = False
    escaped = False

    while i < n:
        ch = text[i]
        if in_string:
            out.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue

        if ch == '"':
            in_string = True
            out.append(ch)
            i += 1
            continue

        if ch == "/" and i + 1 < n:
            nxt = text[i + 1]
            if nxt == "/":
                i += 2
                while i < n and text[i] not in "\r\n":
                    i += 1
                continue
            if nxt == "*":
                i += 2
                while i < n:
                    if i + 1 < n and text[i] == "*" and text[i + 1] == "/":
                        i += 2
                        break
                    if text[i] == "\n":
                        out.append("\n")
                    i += 1
                continue

        out.append(ch)
        i += 1

    return "".join(out)


def _remove_trailing_json_commas(text: str) -> str:
    out: list[str] = []
    i = 0
    n = len(text)
    in_string = False
    escaped = False

    while i < n:
        ch = text[i]
        if in_string:
            out.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue

        if ch == '"':
            in_string = True
            out.append(ch)
            i += 1
            continue

        if ch == ",":
            j = i + 1
            while j < n and text[j] in " \t\r\n":
                j += 1
            if j < n and text[j] in "]}":
                i += 1
                continue

        out.append(ch)
        i += 1

    return "".join(out)


def _decode_json_relaxed(text: str) -> dict[str, Any]:
    without_comments = _strip_json_comments(text)
    without_trailing = _remove_trailing_json_commas(without_comments)
    value = json.loads(without_trailing)
    if not isinstance(value, dict):
        raise ValueError("Top-level JSON value must be an object")
    return value


def load_json(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        raw = handle.read()

    # Fast path: strict JSON.
    try:
        value = json.loads(raw)
        if not isinstance(value, dict):
            raise ValueError(f"Top-level JSON value must be an object: {path}")
        return value
    except json.JSONDecodeError as strict_err:
        # Compatibility path for collaborators who add comments or trailing commas.
        try:
            value = _decode_json_relaxed(raw)
            log(f"WARNING: Loaded non-strict JSON (comments/trailing commas) from {path}")
            return value
        except Exception as relaxed_err:
            lines = raw.splitlines()
            bad_line = ""
            if 1 <= strict_err.lineno <= len(lines):
                bad_line = lines[strict_err.lineno - 1].strip()
            raise ValueError(
                f"Invalid JSON in {path} at line {strict_err.lineno}, "
                f"column {strict_err.colno}: {strict_err.msg}. "
                "Use double-quoted property names and strings; remove malformed commas/comments. "
                f"Offending line: {bad_line}"
            ) from relaxed_err


def save_json(path: str, payload: dict[str, Any]) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=json_default)
        handle.write("\n")


def copy_file_atomic(src_path: str, dst_path: str) -> None:
    ensure_dir(os.path.dirname(dst_path))
    tmp = dst_path + ".tmp"
    shutil.copy2(src_path, tmp)
    os.replace(tmp, dst_path)


def metadata_path_for_data(data_path: str) -> str:
    if data_path.lower().endswith(".dat"):
        return data_path[:-4] + ".meta.json"
    return data_path + ".meta.json"


def load_payload_from_dat(data_path: str, metadata_path: str | None = None) -> dict[str, Any]:
    meta_path = metadata_path or metadata_path_for_data(data_path)
    metadata = load_json(meta_path)
    payload = dict(metadata)
    payload["all_to_all_file"] = data_path
    payload["data"] = mc_engine.load_all_to_all(data_path)
    return payload


def default_exe_path() -> str:
    exe = os.path.join(_RUNTIME_ROOT, "bin", "ising_tri_twisted_parallelogram")
    if os.name == "nt":
        exe += ".exe"
    return exe


def _collect_named_files(root_dir: str, target_name: str) -> list[str]:
    files: list[str] = []
    if not os.path.isdir(root_dir):
        return files
    for root, _, names in os.walk(root_dir):
        if target_name in names:
            files.append(os.path.join(root, target_name))
    return files


def _latest_named_file(root_dir: str, target_name: str) -> str | None:
    candidates = _collect_named_files(root_dir, target_name)
    if not candidates:
        return None
    return max(candidates, key=lambda p: os.path.getmtime(p))


def _safe_mc_scratch_dir(scratch_root: str, label: str) -> str:
    scratch = os.path.join(scratch_root, label)
    if os.name != "nt":
        return scratch
    if len(scratch) <= 180:
        return scratch
    digest = hashlib.sha1(scratch.encode("utf-8")).hexdigest()[:16]
    short_scratch = os.path.join(tempfile.gettempdir(), "large_reference_workflow_mc", digest)
    log(
        "MC scratch path is long for Windows; using short temp scratch: "
        f"{short_scratch}"
    )
    return short_scratch


def _newest_source_mtime() -> float:
    newest = os.path.getmtime(os.path.join(_RUNTIME_ROOT, "Makefile"))
    for folder in ("src", "include"):
        root_dir = os.path.join(_RUNTIME_ROOT, folder)
        for root, _, files in os.walk(root_dir):
            for name in files:
                newest = max(newest, os.path.getmtime(os.path.join(root, name)))
    return newest


def _direct_compile_command(exe: str, compiler: str) -> list[str]:
    src = os.path.join(_RUNTIME_ROOT, "src", "ising_tri_twisted_parallelogram.cc")
    inc = os.path.join(_RUNTIME_ROOT, "include")
    return [
        compiler,
        "-std=c++14",
        "-O3",
        "-Wall",
        "-Wno-sign-compare",
        "-Wno-unknown-warning-option",
        "-Wno-deprecated-declarations",
        f"-I{inc}",
        src,
        "-o",
        exe,
    ]


def _default_build_commands(execution_cfg: dict[str, Any], exe: str) -> list[list[str] | str]:
    commands: list[list[str] | str] = []
    seen: set[str] = set()

    def _append_command(cmd: list[str] | str) -> None:
        key = cmd if isinstance(cmd, str) else "\x1f".join(cmd)
        if key in seen:
            return
        seen.add(key)
        commands.append(cmd)

    configured = execution_cfg.get("build_command")
    if configured is not None:
        _append_command(configured)

    # Prefer direct compiler invocation so Windows runs do not depend on make.
    for compiler in ("g++", "clang++", "c++"):
        if shutil.which(compiler) is not None:
            _append_command(_direct_compile_command(exe, compiler))

    if os.name == "nt":
        _append_command(["mingw32-make"])
    _append_command(["make"])
    return commands


def ensure_simulator(execution_cfg: dict[str, Any]) -> str:
    exe = execution_cfg.get("exe") or default_exe_path()
    if not os.path.isabs(exe):
        exe = os.path.join(_RUNTIME_ROOT, exe)
    exe = os.path.normpath(exe)

    auto_build = bool(execution_cfg.get("auto_build", True))
    force_rebuild = bool(execution_cfg.get("force_rebuild", False))
    needs_build = force_rebuild or (not os.path.exists(exe))
    if auto_build and os.path.exists(exe):
        needs_build = needs_build or (_newest_source_mtime() > os.path.getmtime(exe))

    if auto_build and needs_build:
        timeout_s = int(execution_cfg.get("build_timeout_s", 600))
        ensure_dir(os.path.dirname(exe))
        build_errors: list[str] = []
        build_commands = _default_build_commands(execution_cfg, exe)
        build_ok = False

        for command in build_commands:
            log(f"Building simulator with command: {command}")
            try:
                result = subprocess.run(
                    command,
                    cwd=_RUNTIME_ROOT,
                    capture_output=True,
                    text=True,
                    timeout=timeout_s,
                    shell=isinstance(command, str),
                )
            except FileNotFoundError as exc:
                build_errors.append(
                    f"Command not found: {command}\n"
                    f"error: {exc}"
                )
                continue
            except subprocess.TimeoutExpired as exc:
                build_errors.append(
                    f"Build command timed out after {timeout_s}s: {command}\n"
                    f"error: {exc}"
                )
                continue

            if result.returncode == 0:
                build_ok = True
                break

            build_errors.append(
                f"Build command failed: {command}\n"
                f"stdout:\n{result.stdout}\n\n"
                f"stderr:\n{result.stderr}"
            )

        if not build_ok:
            details = "\n\n".join(build_errors) if build_errors else "No build commands were attempted."
            raise RuntimeError(
                "Simulator build failed.\n"
                f"{details}\n\n"
                "On Windows, install a C++ toolchain (for example MinGW-w64) and ensure "
                "`mingw32-make` or `g++` is available on PATH."
            )

    if not os.path.exists(exe):
        raise FileNotFoundError(
            f"Simulator executable not found at {exe}. "
            "Either set execution.auto_build=true or provide execution.exe."
        )
    return exe


def require_lattice(cfg: dict[str, Any], section_name: str) -> tuple[int, int, int, int]:
    missing = [key for key in ("Lx", "Ly", "Tx", "Ty") if key not in cfg]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")
    return (int(cfg["Lx"]), int(cfg["Ly"]), int(cfg["Tx"]), int(cfg["Ty"]))


def couplings_from_ratio(couplings_cfg: dict[str, Any], section_name: str) -> dict[str, float]:
    missing = [key for key in ("k3", "k1_over_k3", "k2_over_k3") if key not in couplings_cfg]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")

    k3 = float(couplings_cfg["k3"])
    if k3 == 0.0:
        raise ValueError(f"{section_name}.k3 must be non-zero")
    r1 = float(couplings_cfg["k1_over_k3"])
    r2 = float(couplings_cfg["k2_over_k3"])
    return {
        "k3": k3,
        "r1": r1,
        "r2": r2,
        "k1": r1 * k3,
        "k2": r2 * k3,
    }


def build_test_geometry_map(test_family_cfg: dict[str, Any]) -> dict[int, tuple[int, int, int, int]]:
    if "sizes" not in test_family_cfg:
        raise ValueError("test_family.sizes is required")
    sizes = [int(v) for v in test_family_cfg["sizes"]]
    if len(sizes) == 0:
        raise ValueError("test_family.sizes must not be empty")

    geometry_map = test_family_cfg.get("geometry_map")
    if geometry_map is not None:
        out: dict[int, tuple[int, int, int, int]] = {}
        for key, value in geometry_map.items():
            if len(value) != 4:
                raise ValueError(f"test_family.geometry_map[{key}] must have 4 ints")
            out[int(key)] = tuple(int(x) for x in value)
        missing = [L for L in sizes if L not in out]
        if missing:
            raise ValueError(f"test_family.geometry_map is missing sizes: {missing}")
        return {L: out[L] for L in sizes}

    defaults = test_family_cfg.get("geometry_defaults") or {}
    lx_mult = float(defaults.get("Lx_mult", 1.0))
    ly_mult = float(defaults.get("Ly_mult", 1.0))
    tx_frac = float(defaults.get("Tx_frac", 0.0))
    ty_frac = float(defaults.get("Ty_frac", 0.0))

    out = {}
    for L in sizes:
        out[L] = (
            int(round(lx_mult * L)),
            int(round(ly_mult * L)),
            int(round(tx_frac * L)),
            int(round(ty_frac * L)),
        )
    return out


def parse_ratio_list(values: Any, field_name: str) -> list[float]:
    if values is None:
        raise ValueError(f"{field_name} is required")
    out = [float(v) for v in values]
    if len(out) == 0:
        raise ValueError(f"{field_name} must not be empty")
    return out


def _beta_continuum_power_model(
    L: np.ndarray, beta_c: float, amplitude: float, exponent: float
) -> np.ndarray:
    return beta_c + amplitude * np.power(np.asarray(L, dtype=float), -exponent)


def fit_beta_continuum_free_power(
    L_values: list[int] | np.ndarray,
    beta_values: list[float] | np.ndarray,
    sigma_values: list[float] | np.ndarray | None = None,
    *,
    weighted: bool = False,
    exponent_initial: float = 1.0,
    exponent_min: float = 0.05,
    exponent_max: float = 6.0,
) -> dict[str, Any]:
    L = np.asarray(L_values, dtype=float)
    beta = np.asarray(beta_values, dtype=float)
    sigma = np.zeros_like(beta) if sigma_values is None else np.asarray(sigma_values, dtype=float)
    mask = np.isfinite(L) & np.isfinite(beta)
    if sigma.shape == beta.shape:
        mask &= np.isfinite(sigma)
    L = L[mask]
    beta = beta[mask]
    sigma = sigma[mask]

    if len(L) < 3:
        raise ValueError("free-power continuum beta fit requires at least 3 sizes")
    if exponent_max <= exponent_min:
        raise ValueError("free-power continuum beta fit requires exponent_max > exponent_min")

    sigma_fit = None
    if weighted:
        good = sigma > 0.0
        if int(np.count_nonzero(good)) < 3:
            raise ValueError(
                "weighted free-power continuum beta fit requires at least 3 positive sigma values"
            )
        L = L[good]
        beta = beta[good]
        sigma = sigma[good]
        sigma_fit = sigma

    inv_L = 1.0 / L
    linear_coeffs = np.polyfit(inv_L, beta, 1)
    beta_c_guess = float(max(float(beta.max()), float(linear_coeffs[-1])))
    amplitude_guess = float((beta[0] - beta_c_guess) * L[0])
    lower = [float(beta.min() - 0.10), -10.0, float(exponent_min)]
    upper = [float(beta.max() + 0.10), 10.0, float(exponent_max)]
    p0 = [
        min(max(beta_c_guess, lower[0] + 1e-6), upper[0] - 1e-6),
        min(max(amplitude_guess, lower[1] + 1e-6), upper[1] - 1e-6),
        min(max(float(exponent_initial), lower[2] + 1e-6), upper[2] - 1e-6),
    ]
    popt, pcov = curve_fit(
        _beta_continuum_power_model,
        L,
        beta,
        p0=p0,
        bounds=(lower, upper),
        sigma=sigma_fit,
        absolute_sigma=bool(weighted and sigma_fit is not None),
        maxfev=100000,
    )
    beta_fit = _beta_continuum_power_model(L, *popt)
    residual = beta - beta_fit
    beta_c_sigma = None
    if np.all(np.isfinite(pcov)):
        variance = float(pcov[0, 0])
        if variance >= 0.0:
            beta_c_sigma = float(np.sqrt(variance))
    return {
        "method": "free_power_continuum",
        "ansatz": "beta_c(L) = beta_c(infinity) + a * L^(-p)",
        "weighted": bool(weighted),
        "parameters": {
            "beta_c": float(popt[0]),
            "amplitude": float(popt[1]),
            "exponent": float(popt[2]),
            "exponent_initial": float(exponent_initial),
            "exponent_min": float(exponent_min),
            "exponent_max": float(exponent_max),
        },
        "L_values": [int(x) for x in L],
        "beta_values": [float(x) for x in beta],
        "sigma_values": [float(x) for x in sigma],
        "beta_fit_values": [float(x) for x in beta_fit],
        "residual_values": [float(x) for x in residual],
        "beta_c_continuum": float(popt[0]),
        "beta_c_continuum_sigma": beta_c_sigma,
        "rms_residual": float(np.sqrt(np.mean(residual * residual))),
        "n_points": int(len(L)),
        "lmin": int(np.min(L)),
        "lmax": int(np.max(L)),
    }


def find_beta_for_lattice(
    exe: str,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta_cfg: dict[str, Any],
    scratch_root: str,
    label: str,
) -> dict[str, Any]:
    Lx, Ly, Tx, Ty = lattice
    k1 = couplings["k1"]
    k2 = couplings["k2"]
    k3 = couplings["k3"]

    scratch = _safe_mc_scratch_dir(scratch_root, label)
    ensure_dir(scratch)
    t0 = time.time()
    beta_c, beta_c_sigma, chi_peak, scan_betas, scan_chis, scan_chi_errs = mc_engine.find_beta_c(
        exe,
        Lx,
        Ly,
        Tx,
        Ty,
        k1,
        k2,
        k3,
        float(beta_cfg["beta_lo"]),
        float(beta_cfg["beta_hi"]),
        n_coarse=int(beta_cfg["n_coarse"]),
        n_refine=int(beta_cfg["n_refine"]),
        n_refine2=int(beta_cfg["n_refine2"]),
        n_refine3=int(beta_cfg["n_refine3"]),
        n_traj_coarse=int(beta_cfg["n_traj_coarse"]),
        n_traj_fine=int(beta_cfg["n_traj_fine"]),
        data_dir=os.path.join(scratch, "scan"),
        max_shifts=int(beta_cfg.get("max_shifts", 4)),
        jackknife=bool(beta_cfg.get("jackknife", False)),
    )
    return {
        "label": label,
        "created_at": timestamp(),
        "beta_find_wall_seconds": float(time.time() - t0),
        "L": int(max(abs(Lx), abs(Ly))),
        "Lx": int(Lx),
        "Ly": int(Ly),
        "Tx": int(Tx),
        "Ty": int(Ty),
        "k1": float(k1),
        "k2": float(k2),
        "k3": float(k3),
        "r1": float(couplings["r1"]),
        "r2": float(couplings["r2"]),
        "beta_c": float(beta_c),
        "beta_c_sigma": float(beta_c_sigma),
        "chi_peak": float(chi_peak),
        "scan_betas": [float(x) for x in scan_betas],
        "scan_chis": [float(x) for x in scan_chis],
        "scan_chi_errs": [float(x) for x in scan_chi_errs],
        "beta_finder": {
            "beta_lo": float(beta_cfg["beta_lo"]),
            "beta_hi": float(beta_cfg["beta_hi"]),
            "n_coarse": int(beta_cfg["n_coarse"]),
            "n_refine": int(beta_cfg["n_refine"]),
            "n_refine2": int(beta_cfg["n_refine2"]),
            "n_refine3": int(beta_cfg["n_refine3"]),
            "n_traj_coarse": int(beta_cfg["n_traj_coarse"]),
            "n_traj_fine": int(beta_cfg["n_traj_fine"]),
            "max_shifts": int(beta_cfg.get("max_shifts", 4)),
            "jackknife": bool(beta_cfg.get("jackknife", False)),
        },
    }


def run_production_payload(
    exe: str,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta: float,
    mc_cfg: dict[str, Any],
    scratch_root: str,
    label: str,
) -> dict[str, Any]:
    Lx, Ly, Tx, Ty = lattice
    k1 = couplings["k1"]
    k2 = couplings["k2"]
    k3 = couplings["k3"]
    scratch = _safe_mc_scratch_dir(scratch_root, label)
    ensure_dir(scratch)
    prod_data_dir = os.path.join(scratch, "prod")
    t0 = time.time()
    _, subdir = mc_engine.run_simulator(
        exe,
        Lx,
        Ly,
        Tx,
        Ty,
        k1,
        k2,
        k3,
        float(beta),
        n_traj=int(mc_cfg["n_traj"]),
        n_skip=int(mc_cfg["n_skip"]),
        n_therm=int(mc_cfg["n_therm"]),
        data_dir=prod_data_dir,
    )
    if subdir is None:
        recovered = _latest_named_file(prod_data_dir, "two_point_all_to_all.dat")
        if recovered is not None:
            subdir = os.path.dirname(recovered)

    if subdir is None and os.name == "nt":
        retry_data_dir = os.path.join(
            tempfile.gettempdir(),
            "large_reference_workflow_prod_retry",
            f"{int(time.time())}_{os.getpid()}",
        )
        log(
            "Simulator output directory was not detected; retrying production run "
            f"with short temp data_dir: {retry_data_dir}"
        )
        _, subdir = mc_engine.run_simulator(
            exe,
            Lx,
            Ly,
            Tx,
            Ty,
            k1,
            k2,
            k3,
            float(beta),
            n_traj=int(mc_cfg["n_traj"]),
            n_skip=int(mc_cfg["n_skip"]),
            n_therm=int(mc_cfg["n_therm"]),
            data_dir=retry_data_dir,
        )
        if subdir is None:
            recovered = _latest_named_file(retry_data_dir, "two_point_all_to_all.dat")
            if recovered is not None:
                subdir = os.path.dirname(recovered)

    if subdir is None:
        raise RuntimeError(
            "Simulator did not report output directory and no two_point_all_to_all.dat "
            f"was found under {prod_data_dir}."
        )

    all_to_all_path = os.path.join(subdir, "two_point_all_to_all.dat")
    if not os.path.isfile(all_to_all_path):
        recovered = _latest_named_file(subdir, "two_point_all_to_all.dat")
        if recovered is None:
            raise RuntimeError(
                "Simulator output directory resolved but missing two_point_all_to_all.dat: "
                f"{subdir}"
            )
        all_to_all_path = recovered
    return {
        "production_beta": float(beta),
        "production_wall_seconds": float(time.time() - t0),
        "all_to_all_file": all_to_all_path,
        "mc": {
            "n_traj": int(mc_cfg["n_traj"]),
            "n_skip": int(mc_cfg["n_skip"]),
            "n_therm": int(mc_cfg["n_therm"]),
        },
    }


def build_payload_summary(
    *,
    label: str,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    finite_beta_summary: dict[str, Any],
    production_summary: dict[str, Any],
    beta_source: str,
    production_beta_sigma: float | None = None,
    beta_extrapolation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    Lx, Ly, Tx, Ty = lattice
    production_beta = float(production_summary["production_beta"])
    summary = {
        "label": label,
        "created_at": timestamp(),
        "wall_seconds": float(
            float(finite_beta_summary.get("beta_find_wall_seconds", 0.0))
            + float(production_summary.get("production_wall_seconds", 0.0))
        ),
        "L": int(max(abs(Lx), abs(Ly))),
        "Lx": int(Lx),
        "Ly": int(Ly),
        "Tx": int(Tx),
        "Ty": int(Ty),
        "k1": float(couplings["k1"]),
        "k2": float(couplings["k2"]),
        "k3": float(couplings["k3"]),
        "r1": float(couplings["r1"]),
        "r2": float(couplings["r2"]),
        "beta_c": production_beta,
        "beta_c_sigma": None if production_beta_sigma is None else float(production_beta_sigma),
        "production_beta": production_beta,
        "production_beta_sigma": None if production_beta_sigma is None else float(production_beta_sigma),
        "production_beta_source": str(beta_source),
        "finite_size_beta_c": float(finite_beta_summary["beta_c"]),
        "finite_size_beta_c_sigma": float(finite_beta_summary.get("beta_c_sigma", 0.0)),
        "chi_peak": float(finite_beta_summary.get("chi_peak", float("nan"))),
        "scan_betas": [float(x) for x in finite_beta_summary.get("scan_betas", [])],
        "scan_chis": [float(x) for x in finite_beta_summary.get("scan_chis", [])],
        "scan_chi_errs": [float(x) for x in finite_beta_summary.get("scan_chi_errs", [])],
        "beta_find_wall_seconds": float(finite_beta_summary.get("beta_find_wall_seconds", 0.0)),
        "production_wall_seconds": float(production_summary.get("production_wall_seconds", 0.0)),
        "mc": production_summary["mc"],
        "beta_finder": finite_beta_summary["beta_finder"],
        "all_to_all_file": production_summary["all_to_all_file"],
    }
    if beta_extrapolation is not None:
        summary["beta_extrapolation"] = beta_extrapolation
    return summary


def run_one_payload(
    exe: str,
    lattice: tuple[int, int, int, int],
    couplings: dict[str, float],
    beta_cfg: dict[str, Any],
    mc_cfg: dict[str, Any],
    scratch_root: str,
    label: str,
) -> dict[str, Any]:
    finite_beta_summary = find_beta_for_lattice(
        exe=exe,
        lattice=lattice,
        couplings=couplings,
        beta_cfg=beta_cfg,
        scratch_root=scratch_root,
        label=label,
    )
    production_summary = run_production_payload(
        exe=exe,
        lattice=lattice,
        couplings=couplings,
        beta=float(finite_beta_summary["beta_c"]),
        mc_cfg=mc_cfg,
        scratch_root=scratch_root,
        label=label,
    )
    return build_payload_summary(
        label=label,
        lattice=lattice,
        couplings=couplings,
        finite_beta_summary=finite_beta_summary,
        production_summary=production_summary,
        beta_source="per_size_beta",
        production_beta_sigma=maybe_float(finite_beta_summary.get("beta_c_sigma")),
    )


def sample_directional_channels(
    payload: dict[str, Any],
    fractions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    fractions = np.asarray(fractions, dtype=float)
    if np.any((fractions < 0.0) | (fractions > 1.0)):
        raise ValueError("fractions must lie in [0, 1]")

    Lx = int(payload["Lx"])
    Ly = int(payload["Ly"])
    Tx = int(payload["Tx"])
    Ty = int(payload["Ty"])
    data = payload["data"]

    interp_g = _tile_interp(data, Lx, Ly, Tx, Ty, "conn", copies=2)
    interp_s = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)

    G = np.full((3, len(fractions)), np.nan, dtype=float)
    sG = np.full_like(G, np.nan)

    for cycle, (dm, dn) in enumerate(boundary_paths(Lx, Ly, Tx, Ty)):
        ex = dm + 0.5 * dn
        ey = _SQRT3_2 * dn
        pts = np.column_stack([fractions * ex, fractions * ey])
        G[cycle] = np.asarray(interp_g(pts), dtype=float).ravel()
        sG[cycle] = np.abs(np.asarray(interp_s(pts), dtype=float).ravel())

    return G, sG


def write_dat(path: str, header_lines: list[str], columns: list[str], rows: list[list[Any]]) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as handle:
        for line in header_lines:
            handle.write(f"# {line}\n")
        handle.write("# columns: " + " ".join(columns) + "\n")
        for row in rows:
            chunks = []
            for value in row:
                if isinstance(value, (int, np.integer)):
                    chunks.append(str(int(value)))
                elif isinstance(value, (float, np.floating)):
                    chunks.append(f"{float(value):.10g}")
                else:
                    chunks.append(str(value))
            handle.write(" ".join(chunks) + "\n")


def payload_file_name(L: int, r1: float, r2: float) -> str:
    return f"test_L{int(L)}_r1{float(r1):.6f}_r2{float(r2):.6f}.dat"


def check_required_sections(config: dict[str, Any], section_names: list[str]) -> None:
    missing = [name for name in section_names if name not in config]
    if missing:
        raise ValueError(f"Config is missing top-level sections: {missing}")
