#!/usr/bin/env python3
"""Generate 7-panel plots for L=32, shifts 0..32, and bundle into one PDF.

Each page in the output PDF is one shift value.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def canonicalize(x: int, y: int, L: int, shift: int) -> tuple[int, int]:
    y0 = y % L
    b = (y - y0) // L
    x0 = (x - b * shift) % L
    return x0, y0


def orbit_period(step: tuple[int, int], L: int, shift: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, L, shift)
    for n in range(1, L * L + 1):
        x, y = canonicalize(x + sx, y + sy, L, shift)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period found for step={step}, L={L}, shift={shift}")


def run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def copy_set(src_prefix: Path, dst_suffix: str) -> None:
    for direction in ["base", "left", "diag"]:
        src_png = Path(f"{src_prefix}_{direction}.png")
        src_svg = Path(f"{src_prefix}_{direction}.svg")
        dst_png = Path(f"{src_prefix}_{direction}_{dst_suffix}.png")
        dst_svg = Path(f"{src_prefix}_{direction}_{dst_suffix}.svg")
        dst_png.write_bytes(src_png.read_bytes())
        dst_svg.write_bytes(src_svg.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=32)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=32)
    parser.add_argument("--x-length", type=float, default=1.0)
    parser.add_argument("--decay-scale", type=float, default=10.0)
    parser.add_argument("--out-pdf", default="thinkTwist/two_point_panels_L32_shifts_0_to_32.pdf")
    parser.add_argument("--summary-csv", default="thinkTwist/two_point_panels_L32_shifts_0_to_32_summary.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    think = root / "thinkTwist"
    panel_paths: list[Path] = []
    rows: list[dict[str, float | int]] = []

    for shift in range(args.shift_min, args.shift_max + 1):
        px = orbit_period((1, 0), args.L, shift)
        pl = orbit_period((0, 1), args.L, shift)
        pd = orbit_period((-1, 1), args.L, shift)

        left_len = args.x_length * px / pl
        diag_len = args.x_length * px / pd

        # Base (equal lengths)
        run(
            [
                str(root / ".venv/bin/python"),
                "thinkTwist/plot_twisted_directional_two_point_schematic.py",
                "--L",
                str(args.L),
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
        prefix = think / f"twisted_base_L{args.L}_shift{shift}_one_period_fractional_two_point"
        copy_set(prefix, "BASE_equal")

        # Compensated
        run(
            [
                str(root / ".venv/bin/python"),
                "thinkTwist/plot_twisted_directional_two_point_schematic.py",
                "--L",
                str(args.L),
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

        # 7-panel per shift
        panel_png = think / f"two_point_base_vs_comp_panel_L{args.L}_shift{shift}.png"
        run(
            [
                str(root / ".venv/bin/python"),
                "thinkTwist/make_base_vs_comp_panel.py",
                "--L",
                str(args.L),
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
                "L": args.L,
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

    # Write summary CSV.
    summary_csv = root / args.summary_csv
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "L",
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

    # Build multi-page PDF.
    out_pdf = root / args.out_pdf
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_pdf) as pdf:
        for panel in panel_paths:
            img = mpimg.imread(panel)
            h, w = img.shape[0], img.shape[1]
            fig = plt.figure(figsize=(w / 140.0, h / 140.0))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0)
            plt.close(fig)

    print(f"Wrote {summary_csv}")
    print(f"Wrote {out_pdf}")
    print(f"Pages: {len(panel_paths)}")


if __name__ == "__main__":
    main()
