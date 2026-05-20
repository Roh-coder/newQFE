#!/usr/bin/env python3
"""Plot test vs reference connected two-point functions at L=32.

For every test all-to-all output at L=32 in
``production_iso111_ref4x_qtr_tests_20260515``, build a four-panel figure:

  1. test G(m,n)
  2. reference G(m,n) (1,1 reference at L=32)
  3. difference (test - reference)
  4. z-score (test - reference) / sqrt(err_test^2 + err_ref^2)

Each panel is shown as a 2x2 periodic tiling of the torus, with the underlying
32x32 lattice values upsampled by bilinear / Q1 FEM linear interpolation onto a
finer grid for smooth heatmaps.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "K_from_continuum"))
from lib.mc_engine import load_all_to_all  # noqa: E402

RESULTS = REPO_ROOT / "K_from_continuum" / "results"
TEST_RUN = RESULTS / "production_iso111_ref4x_qtr_tests_20260515"
REF_RUN = RESULTS / "production_iso111_ref4x_qtr_reference_20260515"
OUT_DIR = RESULTS / "L32_diff_plots"

L = 32
UPSAMPLE = 4  # 32 -> 128 fine grid per tile; 2x2 tiling => 256x256 image

META_RE = re.compile(r"_L(\d+)_r1([0-9.]+)_r2([0-9.]+)\.meta\.json$")


def to_field(data: dict, L: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (G, sigma) as L x L arrays indexed by (m, n)."""
    G = np.zeros((L, L), dtype=float)
    sig = np.zeros((L, L), dtype=float)
    seen = np.zeros((L, L), dtype=bool)
    for (m, n), rec in data.items():
        mi = int(m) % L
        ni = int(n) % L
        G[mi, ni] = rec["conn"]
        sig[mi, ni] = rec["conn_err"]
        seen[mi, ni] = True
    if not seen.all():
        missing = int((~seen).sum())
        raise RuntimeError(f"missing {missing} of {L*L} (m,n) entries")
    return G, sig


def fem_q1_upsample(field: np.ndarray, up: int) -> np.ndarray:
    """Periodic bilinear (Q1 FEM) upsample of an LxL field by factor ``up``."""
    L = field.shape[0]
    ext = np.empty((L + 1, L + 1), dtype=field.dtype)
    ext[:L, :L] = field
    ext[L, :L] = field[0, :]
    ext[:L, L] = field[:, 0]
    ext[L, L] = field[0, 0]

    nfine = L * up
    t = np.arange(nfine) / up  # 0, 1/up, ..., L - 1/up
    i = np.floor(t).astype(int)
    f = t - i
    # rows
    row_lo = ext[i, :]
    row_hi = ext[i + 1, :]
    by_row = (1.0 - f)[:, None] * row_lo + f[:, None] * row_hi  # (nfine, L+1)
    # cols
    col_lo = by_row[:, i]
    col_hi = by_row[:, i + 1]
    fine = (1.0 - f)[None, :] * col_lo + f[None, :] * col_hi
    return fine


def tile_2x2(field: np.ndarray) -> np.ndarray:
    return np.tile(field, (2, 2))


def find_test_metas() -> list[Path]:
    grid_dir = TEST_RUN / "test_data" / "grid" / f"L{L}"
    metas = []
    for p in sorted(grid_dir.glob("test_L*_r1*_r2*.meta.json")):
        m = META_RE.search(p.name)
        if m and int(m.group(1)) == L:
            metas.append(p)
    return metas


def load_meta_field(meta_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    import json
    meta = json.loads(meta_path.read_text())
    a2a = meta["all_to_all_file"]
    data = load_all_to_all(a2a)
    G, sig = to_field(data, L)
    return G, sig, meta


def plot_one(test_meta: Path, ref_G: np.ndarray, ref_sig: np.ndarray,
             ref_label: str) -> None:
    G_t, sig_t, meta = load_meta_field(test_meta)
    r1, r2 = meta["r1"], meta["r2"]

    diff = G_t - ref_G
    denom = np.sqrt(sig_t * sig_t + ref_sig * ref_sig)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(denom > 0, diff / denom, 0.0)

    fields = [G_t, ref_G, diff, z]
    titles = [
        f"test G  (r1={r1:.4f}, r2={r2:.4f})",
        f"reference G  ({ref_label})",
        "diff = test - ref",
        "z-score = diff / sqrt(err_t^2 + err_r^2)",
    ]
    cmaps = ["viridis", "viridis", "RdBu_r", "RdBu_r"]

    # Symmetric color limits for diff and z; max-abs from interpolated tile is
    # essentially the same as raw, so use raw arrays for limits.
    diff_lim = float(np.max(np.abs(diff))) if np.any(diff) else 1.0
    z_lim = float(np.max(np.abs(z))) if np.any(z) else 1.0

    fig, axes = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
    for ax, field, title, cmap in zip(axes.flat, fields, titles, cmaps):
        fine = fem_q1_upsample(field, UPSAMPLE)
        tiled = tile_2x2(fine)
        if cmap == "RdBu_r":
            lim = diff_lim if "diff" in title else z_lim
            im = ax.imshow(
                tiled, origin="lower", cmap=cmap,
                extent=(0, 2 * L, 0, 2 * L),
                vmin=-lim, vmax=lim, interpolation="nearest",
            )
        else:
            im = ax.imshow(
                tiled, origin="lower", cmap=cmap,
                extent=(0, 2 * L, 0, 2 * L), interpolation="nearest",
            )
        fig.colorbar(im, ax=ax, shrink=0.85)
        for k in (L,):
            ax.axhline(k, color="white", lw=0.5, alpha=0.4)
            ax.axvline(k, color="white", lw=0.5, alpha=0.4)
        ax.set_title(title)
        ax.set_xlabel("m (2x2 tile)")
        ax.set_ylabel("n (2x2 tile)")

    fig.suptitle(
        f"L=32 connected two-point: test (r1={r1:.4f}, r2={r2:.4f})"
        f" vs reference (1,1) — FEM Q1 upsample x{UPSAMPLE}, periodic 2x2 tiling",
        fontsize=12,
    )

    out_name = f"diff_L32_r1{r1:.4f}_r2{r2:.4f}.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / out_name
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"wrote {out_path.relative_to(REPO_ROOT)}  "
          f"|max diff|={diff_lim:.3e}  |max z|={z_lim:.2f}")


def main() -> int:
    ref_meta = REF_RUN / "reference_data" / "grid" / f"L{L}" \
        / f"reference_L{L}_r11.000000_r21.000000.meta.json"
    ref_G, ref_sig, _ = load_meta_field(ref_meta)
    ref_label = f"L={L}, r1=1, r2=1"

    metas = find_test_metas()
    if not metas:
        print(f"no test metas found under {TEST_RUN}", file=sys.stderr)
        return 1
    print(f"reference: {ref_meta.relative_to(REPO_ROOT)}")
    print(f"{len(metas)} test points at L={L}")
    for mp in metas:
        plot_one(mp, ref_G, ref_sig, ref_label)
    print(f"\nall plots written to {OUT_DIR.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
