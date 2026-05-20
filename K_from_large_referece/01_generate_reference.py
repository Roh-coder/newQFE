#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from typing import Any

import numpy as np

from workflow_common import (
    check_required_sections,
    copy_file_atomic,
    couplings_from_ratio,
    ensure_dir,
    ensure_simulator,
    load_json,
    load_payload_from_dat,
    log,
    metadata_path_for_data,
    require_lattice,
    resolve_path,
    run_one_payload,
    sample_directional_channels,
    save_json,
    timestamp,
    write_dat,
)

_HERE = os.path.dirname(os.path.abspath(__file__))

# Edit this when running directly from an IDE without CLI args.
DEFAULT_CONFIG_PATH = "configs/reference_example.json"


def _check_required_keys(section: dict[str, Any], section_name: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in section]
    if missing:
        raise ValueError(f"{section_name} is missing keys: {missing}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 1: generate one large reference all-to-all two-point payload "
            "with configured lattice, beta-finder, and MC settings."
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
        ["run", "paths", "execution", "reference_lattice", "reference_couplings", "beta_finder", "mc"],
    )

    if args.tag is not None:
        cfg["run"]["tag"] = args.tag

    _check_required_keys(cfg["run"], "run", ["tag"])
    _check_required_keys(cfg["paths"], "paths", ["results_root", "resume"])
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

    tag = str(cfg["run"]["tag"])
    results_root = resolve_path(str(cfg["paths"]["results_root"]), _HERE)
    run_root = os.path.join(results_root, tag)
    reference_dir = os.path.join(run_root, "reference_data")
    scratch_root = os.path.join(run_root, "_mc_scratch", "reference_data")
    ensure_dir(reference_dir)
    ensure_dir(scratch_root)

    reference_data_path = os.path.join(reference_dir, "reference_all_to_all.dat")
    reference_meta_path = metadata_path_for_data(reference_data_path)
    legacy_pickle_path = os.path.join(reference_dir, "reference_payload.pkl")
    channels_path = os.path.join(reference_dir, "reference_channels.dat")
    manifest_path = os.path.join(reference_dir, "manifest_reference.json")

    resume = bool(cfg["paths"]["resume"])
    if resume and os.path.exists(reference_data_path) and os.path.exists(reference_meta_path):
        log(f"Reusing existing reference data: {reference_data_path}")
        if os.path.exists(legacy_pickle_path):
            os.remove(legacy_pickle_path)
            log(f"Removed legacy pickle artifact: {legacy_pickle_path}")
        summary = None
    else:
        exe = ensure_simulator(cfg["execution"])
        lattice = require_lattice(cfg["reference_lattice"], "reference_lattice")
        couplings = couplings_from_ratio(cfg["reference_couplings"], "reference_couplings")

        label = (
            f"reference_Lx{lattice[0]}_Ly{lattice[1]}_Tx{lattice[2]}_Ty{lattice[3]}"
            f"_r1{couplings['r1']:.6f}_r2{couplings['r2']:.6f}"
        )
        log("Running reference beta finder + production simulation")
        summary = run_one_payload(
            exe=exe,
            lattice=lattice,
            couplings=couplings,
            beta_cfg=cfg["beta_finder"],
            mc_cfg=cfg["mc"],
            scratch_root=scratch_root,
            label=label,
        )
        copy_file_atomic(summary["all_to_all_file"], reference_data_path)
        save_json(reference_meta_path, summary)
        if os.path.exists(legacy_pickle_path):
            os.remove(legacy_pickle_path)
            log(f"Removed legacy pickle artifact: {legacy_pickle_path}")
        log(f"Reference data written: {reference_data_path}")

    payload = load_payload_from_dat(reference_data_path, reference_meta_path)

    analysis_points = cfg.get("analysis_points") or {}
    k_values = [int(v) for v in analysis_points.get("k_values", [1, 2, 3, 4, 5, 6, 7])]
    k_denom = int(analysis_points.get("k_denominator", 8))
    if k_denom <= 0:
        raise ValueError("analysis_points.k_denominator must be positive")
    fractions = np.array(k_values, dtype=float) / float(k_denom)

    G_ref, sG_ref = sample_directional_channels(payload, fractions)
    channel_rows: list[list[Any]] = []
    for cycle in range(3):
        for ik, kval in enumerate(k_values):
            channel_rows.append(
                [
                    cycle,
                    int(kval),
                    float(fractions[ik]),
                    float(G_ref[cycle, ik]),
                    float(sG_ref[cycle, ik]),
                ]
            )

    header = [
        "Large reference channels sampled at directional fractions",
        f"run_tag={tag}",
        f"lattice=[{payload['Lx']}, {payload['Ly']}, {payload['Tx']}, {payload['Ty']}]",
        f"couplings=(k1={payload['k1']}, k2={payload['k2']}, k3={payload['k3']})",
        f"ratios=(k1_over_k3={payload['r1']}, k2_over_k3={payload['r2']})",
        f"beta_c={payload['beta_c']}",
        f"k_values={k_values}",
        f"k_denominator={k_denom}",
    ]
    write_dat(
        channels_path,
        header,
        ["cycle", "k", "t", "G_ref", "sigma_ref"],
        channel_rows,
    )

    manifest = {
        "created_at": timestamp(),
        "config_path": config_path,
        "run_tag": tag,
        "run_root": run_root,
        "reference_data_dir": reference_dir,
        "reference_data_file": reference_data_path,
        "reference_metadata_file": reference_meta_path,
        "reference_channels": channels_path,
        "k_values": k_values,
        "k_denominator": k_denom,
        "reference_summary": {
            "Lx": int(payload["Lx"]),
            "Ly": int(payload["Ly"]),
            "Tx": int(payload["Tx"]),
            "Ty": int(payload["Ty"]),
            "k1": float(payload["k1"]),
            "k2": float(payload["k2"]),
            "k3": float(payload["k3"]),
            "r1": float(payload["r1"]),
            "r2": float(payload["r2"]),
            "beta_c": float(payload["beta_c"]),
            "beta_c_sigma": float(payload["beta_c_sigma"]),
            "chi_peak": float(payload["chi_peak"]),
            "wall_seconds": float(payload["wall_seconds"]),
        },
        "config": cfg,
    }
    save_json(manifest_path, manifest)
    log(f"Reference manifest written: {manifest_path}")


if __name__ == "__main__":
    main()
