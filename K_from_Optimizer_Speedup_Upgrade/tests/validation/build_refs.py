"""build_refs.py — Build reference correlator caches for every
geometry in REFERENCE_GEOMETRIES.

Phase V.2 of the validation panel.  This script loops over the 5
reference geometries (R-equi, R-equiL, R-twist, R-rect, R-shear) and
dispatches a `build_reference()` call for each one with the
panel-budget MC stats (ref_n_traj=40k, scan stats from `_base_cfg`
in run_panel.py).

Usage
-----
    python tests/validation/build_refs.py                 # all 5 geoms
    python tests/validation/build_refs.py --geoms R-equi  # subset
    python tests/validation/build_refs.py --dry-run

Outputs
-------
    results/_reference_Lx<..>_Ly<..>_Tx<..>_Ty<..>/
        two_point_all_to_all.dat
        reference_meta.json
        reference_betac_scan.json

Cached entries are reused on a re-run (build_reference checks the
geometry hash), so this script is idempotent.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKG = HERE.parent.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

from tests.validation.geometries import (  # noqa: E402
    REFERENCE_GEOMETRIES, list_geometries,
)


def _ref_cfg(geom: dict) -> dict:
    """Build a minimal CONFIG dict for build_reference()."""
    return {
        "ref_Lx": geom["ref_Lx"], "ref_Ly": geom["ref_Ly"],
        "ref_Tx": geom["ref_Tx"], "ref_Ty": geom["ref_Ty"],
        "ref_n_traj":            40_000,
        "ref_scan_n_traj_coarse": 4_000,
        "ref_scan_n_traj_fine":  10_000,
        "ref_beta_seed":         (0.20, 0.32),
        "ref_scan_n_coarse":     11,
        "ref_scan_n_refine":     5,
        "ref_scan_n_refine2":    5,
        "ref_scan_n_refine3":    5,
        "ref_scan_max_shifts":   4,
        "scan_jackknife":        False,
        "exe":                   "bin/ising_tri_twisted_parallelogram",
    }


def build_one(tag: str, dry_run: bool = False) -> dict:
    geom = REFERENCE_GEOMETRIES[tag]
    cfg = _ref_cfg(geom)
    ref_tag = (f"Lx{cfg['ref_Lx']}_Ly{cfg['ref_Ly']}"
               f"_Tx{cfg['ref_Tx']}_Ty{cfg['ref_Ty']}")
    ref_dir = PKG / "results" / f"_reference_{ref_tag}"
    print(f"[refs] {tag:<10} -> {ref_dir}")
    if dry_run:
        return {"tag": tag, "status": "DRY-RUN", "ref_dir": str(ref_dir)}
    ref_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(PKG)
    import run as _run
    _run.ensure_binary(cfg["exe"])
    t0 = time.time()
    _, meta = _run.build_reference(cfg, str(ref_dir))
    wall = time.time() - t0
    return {"tag": tag, "status": "OK", "ref_dir": str(ref_dir),
            "beta_c": meta.get("beta_c"), "wall_s": wall}


def _main():
    p = argparse.ArgumentParser()
    p.add_argument("--geoms", default=None,
                   help="Comma-separated subset of geometry tags")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    tags = (args.geoms.split(",") if args.geoms
            else list_geometries())
    results = []
    for tag in tags:
        if tag not in REFERENCE_GEOMETRIES:
            print(f"[refs] unknown tag {tag!r}, skipping")
            continue
        results.append(build_one(tag, dry_run=args.dry_run))
    print("\n[refs] summary:")
    for r in results:
        bc = r.get("beta_c")
        wall = r.get("wall_s", 0)
        bc_s = f"β_c={bc:.6f}" if bc is not None else "—"
        print(f"  {r['tag']:<10} {r['status']:<8} {bc_s}  "
              f"wall={wall:.0f}s")


if __name__ == "__main__":
    _main()
