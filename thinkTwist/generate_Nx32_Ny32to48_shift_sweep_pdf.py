#!/usr/bin/env python3
"""Generate 7-panel plots for Nx=32, Ny in range, shifts 0..Nx, and bundle into PDF.

Each page corresponds to one (Ny, shift) pair.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def canonicalize(x: int, y: int, nx: int, ny: int, shift: int) -> tuple[int, int]:
    y0 = y % ny
    b = (y - y0) // ny
    x0 = (x - b * shift) % nx
    return x0, y0


def orbit_period(step: tuple[int, int], nx: int, ny: int, shift: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, nx, ny, shift)
    for n in range(1, nx * ny + 1):
        x, y = canonicalize(x + sx, y + sy, nx, ny, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError


def run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def copy_set(prefix: Path, suffix: str) -> None:
    for d in ["base", "left", "diag"]:
        (Path(f"{prefix}_{d}_{suffix}.png")).write_bytes(Path(f"{prefix}_{d}.png").read_bytes())
        (Path(f"{prefix}_{d}_{suffix}.svg")).write_bytes(Path(f"{prefix}_{d}.svg").read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Nx", type=int, default=32)
    parser.add_argument("--Ny-min", type=int, default=32)
    parser.add_argument("--Ny-max", type=int, default=48)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=32)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument("--decay-scale", type=float, default=10.0)
    parser.add_argument("--out-pdf", default="thinkTwist/two_point_panels_Nx32_Ny32to48_shifts_0_to_32.pdf")
    parser.add_argument("--summary-csv", default="thinkTwist/two_point_panels_Nx32_Ny32to48_shifts_0_to_32_summary.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    think = root / "thinkTwist"

    rows = []
    panel_paths: list[Path] = []

    for ny in range(args.Ny_min, args.Ny_max + 1):
        for shift in range(args.shift_min, args.shift_max + 1):
            px = orbit_period((1, 0), args.Nx, ny, shift)
            pl = orbit_period((0, 1), args.Nx, ny, shift)
            pd = orbit_period((-1, 1), args.Nx, ny, shift)

            left_len = args.x_length * px / pl
            diag_len = args.x_length * px / pd

            # base
            run(
                [
                    str(root / ".venv/bin/python"),
                    "thinkTwist/plot_twisted_directional_two_point_schematic.py",
                    "--Nx",
                    str(args.Nx),
                    "--Ny",
                    str(ny),
                    "--shift",
                    str(shift),
                    "--fractional-one-period",
                    "--x-length",
                    str(args.x_length),
                    "--left-length",
                    str(args.x_length),
                    "--diag-length",
                    str(args.x_length),
                    "--decay-scale",
                    str(args.decay_scale),
                ],
                root,
            )

            if ny == args.Nx:
                prefix = think / f"twisted_base_L{args.Nx}_shift{shift}_one_period_fractional_two_point"
            else:
                prefix = think / f"twisted_base_Nx{args.Nx}_Ny{ny}_shift{shift}_one_period_fractional_two_point"
            copy_set(prefix, "BASE_equal")

            # compensated
            run(
                [
                    str(root / ".venv/bin/python"),
                    "thinkTwist/plot_twisted_directional_two_point_schematic.py",
                    "--Nx",
                    str(args.Nx),
                    "--Ny",
                    str(ny),
                    "--shift",
                    str(shift),
                    "--fractional-one-period",
                    "--x-length",
                    str(args.x_length),
                    "--left-length",
                    str(left_len),
                    "--diag-length",
                    str(diag_len),
                    "--decay-scale",
                    str(args.decay_scale),
                ],
                root,
            )
            copy_set(prefix, "COMPENSATED")

            panel_png = think / f"two_point_base_vs_comp_panel_Nx{args.Nx}_Ny{ny}_shift{shift}.png"
            run(
                [
                    str(root / ".venv/bin/python"),
                    "thinkTwist/make_base_vs_comp_panel.py",
                    "--Nx",
                    str(args.Nx),
                    "--Ny",
                    str(ny),
                    "--shift",
                    str(shift),
                    "--x-length",
                    str(args.x_length),
                    "--out",
                    str(panel_png),
                ],
                root,
            )

            panel_paths.append(panel_png)
            rows.append(
                {
                    "Nx": args.Nx,
                    "Ny": ny,
                    "shift": shift,
                    "p_x": px,
                    "p_left": pl,
                    "p_diag": pd,
                    "x_length": args.x_length,
                    "left_length": left_len,
                    "diag_length": diag_len,
                    "decay_scale": args.decay_scale,
                    "panel_png": str(panel_png.relative_to(root)),
                }
            )

    summary = root / args.summary_csv
    with summary.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "Nx",
                "Ny",
                "shift",
                "p_x",
                "p_left",
                "p_diag",
                "x_length",
                "left_length",
                "diag_length",
                "decay_scale",
                "panel_png",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    pdf_path = root / args.out_pdf
    with PdfPages(pdf_path) as pdf:
        for panel in panel_paths:
            img = mpimg.imread(panel)
            h, w = img.shape[0], img.shape[1]
            fig = plt.figure(figsize=(w / 150.0, h / 150.0))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0)
            plt.close(fig)

    print(f"Wrote {summary}")
    print(f"Wrote {pdf_path}")
    print(f"Pages: {len(panel_paths)}")


if __name__ == "__main__":
    main()
