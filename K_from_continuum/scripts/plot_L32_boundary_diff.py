#!/usr/bin/env python3
"""Plot test vs reference *boundary* connected two-point correlators at L=32.

For each test all-to-all output at L=32 in
``production_iso111_ref4x_qtr_tests_20260515``, extract the three boundary
slices (v, u, w) of the parallelogram via the existing FEM linear
interpolation in ``mc_engine.extract_boundary_slices`` and produce a 3x3
figure:

  rows : boundary direction (v, u, w)
  cols : (test vs reference overlay)  (diff = test - ref)  (z-score)

Reference is the L=32, (r1, r2) = (1, 1) production run.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "K_from_continuum"))
from lib.mc_engine import (  # noqa: E402
    extract_boundary_slices,
    load_all_to_all,
)

RESULTS = REPO_ROOT / "K_from_continuum" / "results"
TEST_RUN = RESULTS / "production_iso111_ref4x_qtr_tests_20260515"
REF_RUN = RESULTS / "production_iso111_ref4x_qtr_reference_20260515"
OUT_DIR = RESULTS / "L32_boundary_diff_plots"

L = 32
N_SAMPLES = 400

META_RE = re.compile(r"_L(\d+)_r1([0-9.]+)_r2([0-9.]+)\.meta\.json$")


def find_test_metas() -> list[Path]:
    grid_dir = TEST_RUN / "test_data" / "grid" / f"L{L}"
    metas = []
    for p in sorted(grid_dir.glob("test_L*_r1*_r2*.meta.json")):
        m = META_RE.search(p.name)
        if m and int(m.group(1)) == L:
            metas.append(p)
    return metas


def load_meta(meta_path: Path) -> tuple[dict, dict]:
    meta = json.loads(meta_path.read_text())
    data = load_all_to_all(meta["all_to_all_file"])
    return data, meta


def plot_one(test_meta: Path, ref_data: dict) -> None:
    test_data, meta = load_meta(test_meta)
    r1, r2 = meta["r1"], meta["r2"]
    Lx = int(meta.get("Lx", L)); Ly = int(meta.get("Ly", L))
    Tx = int(meta.get("Tx", 0)); Ty = int(meta.get("Ty", 0))

    slices = extract_boundary_slices(
        ref_data, test_data,
        Lx, Ly, Tx, Ty,
        L, L, 0, 0,
        copies=2, n_samples=N_SAMPLES,
    )

    fig, axes = plt.subplots(3, 3, figsize=(13, 10), constrained_layout=True,
                             sharex=True)
    fig.suptitle(
        f"L=32 boundary connected two-point: test (r1={r1:.4f}, r2={r2:.4f})"
        f" vs reference (1,1) — FEM linear interpolation along boundary paths",
        fontsize=12,
    )

    max_abs_z = 0.0
    max_abs_diff = 0.0
    for row, sl in enumerate(slices):
        t = sl["t"]
        gr = sl["g_ref"]; gt = sl["g_test"]
        diff = sl["diff"]; derr = sl["diff_err"]
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(derr > 0, diff / derr, 0.0)
        max_abs_diff = max(max_abs_diff, float(np.max(np.abs(diff))) if diff.size else 0.0)
        max_abs_z = max(max_abs_z, float(np.max(np.abs(z))) if z.size else 0.0)

        ax = axes[row, 0]
        ax.plot(t, gr, "-", color="C0", label="reference (1,1)")
        ax.plot(t, gt, "-", color="C3", label=f"test ({r1:.3f},{r2:.3f})")
        ax.set_ylabel(f"direction {sl['label']}\nG_conn")
        ax.grid(True, alpha=0.3)
        if row == 0:
            ax.legend(loc="best", fontsize=8)
            ax.set_title("test vs reference")

        ax = axes[row, 1]
        ax.axhline(0.0, color="k", lw=0.5)
        ax.fill_between(t, diff - derr, diff + derr, color="C2", alpha=0.3,
                        label="±1σ")
        ax.plot(t, diff, "-", color="C2")
        ax.set_ylabel("test − ref")
        ax.grid(True, alpha=0.3)
        if row == 0:
            ax.set_title("difference")
            ax.legend(loc="best", fontsize=8)

        ax = axes[row, 2]
        ax.axhline(0.0, color="k", lw=0.5)
        for zl in (-3, -1, 1, 3):
            ax.axhline(zl, color="gray", lw=0.4, ls="--", alpha=0.6)
        ax.plot(t, z, "-", color="C4")
        ax.set_ylabel("z-score")
        ax.grid(True, alpha=0.3)
        if row == 0:
            ax.set_title("z-score = diff / σ")

    for ax in axes[-1, :]:
        ax.set_xlabel("path parameter t ∈ [0, 1]")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"boundary_diff_L32_r1{r1:.4f}_r2{r2:.4f}.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"wrote {out_path.relative_to(REPO_ROOT)}  "
          f"|max diff|={max_abs_diff:.3e}  |max z|={max_abs_z:.2f}")


def main() -> int:
    ref_meta_path = (REF_RUN / "reference_data" / "grid" / f"L{L}"
                     / f"reference_L{L}_r11.000000_r21.000000.meta.json")
    ref_data, _ = load_meta(ref_meta_path)

    metas = find_test_metas()
    if not metas:
        print(f"no test metas found under {TEST_RUN}", file=sys.stderr)
        return 1
    print(f"reference: {ref_meta_path.relative_to(REPO_ROOT)}")
    print(f"{len(metas)} test points at L={L}")
    for mp in metas:
        plot_one(mp, ref_data)
    print(f"\nall plots written to {OUT_DIR.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
