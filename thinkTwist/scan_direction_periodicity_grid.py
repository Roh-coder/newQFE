#!/usr/bin/env python3
"""Scan directional orbit periodicities over (Nx, Ny, Ns) parameter space.

Boundary conditions:
  (i + Nx, j) ~ (i, j)
  (i, j + Ny) ~ (i + Ns, j)

Directions scanned:
  x      = (1, 0)
  y      = (0, 1)
  y_minus_x = (-1, 1)
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def canonicalize(x: int, y: int, nx: int, ny: int, ns: int) -> tuple[int, int]:
    y0 = y % ny
    b = (y - y0) // ny
    x0 = (x - b * ns) % nx
    return x0, y0


def period(step: tuple[int, int], nx: int, ny: int, ns: int) -> int:
    sx, sy = step
    x, y = 0, 0
    x, y = canonicalize(x, y, nx, ny, ns)

    # Finite quotient has nx*ny states, so period is <= nx*ny.
    for n in range(1, nx * ny + 1):
        x, y = canonicalize(x + sx, y + sy, nx, ny, ns)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period found for step={step}, nx={nx}, ny={ny}, ns={ns}")


def parse_int_list(text: str) -> list[int]:
    out: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        out.append(int(item))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nx-list",
        default="4,8,12,16",
        help="Comma-separated Nx values",
    )
    parser.add_argument(
        "--ny-list",
        default="4,8,12,16",
        help="Comma-separated Ny values",
    )
    parser.add_argument(
        "--ns-mode",
        choices=["full", "simple"],
        default="full",
        help="full: Ns=0..Nx-1, simple: Ns in {0,1,2,Nx//2} filtered to valid unique",
    )
    parser.add_argument(
        "--ns-list",
        default="",
        help="Optional comma-separated Ns values; when provided, overrides --ns-mode",
    )
    parser.add_argument(
        "--out-prefix",
        default="thinkTwist/periodicity_scan_NxNyNs",
        help="Output prefix (without extension)",
    )
    args = parser.parse_args()

    nx_list = parse_int_list(args.nx_list)
    ny_list = parse_int_list(args.ny_list)
    ns_list = parse_int_list(args.ns_list) if args.ns_list.strip() else []

    rows: list[dict[str, int]] = []

    for nx in nx_list:
        if nx <= 0:
            raise ValueError(f"Invalid Nx={nx}")
        for ny in ny_list:
            if ny <= 0:
                raise ValueError(f"Invalid Ny={ny}")

            if ns_list:
                ns_values = ns_list
            elif args.ns_mode == "full":
                ns_values = list(range(nx))
            else:
                candidates = [0, 1, 2, nx // 2]
                ns_values = sorted({ns for ns in candidates if 0 <= ns < nx})

            for ns in ns_values:
                row = {
                    "Nx": nx,
                    "Ny": ny,
                    "Ns": ns,
                    "p_x": period((1, 0), nx, ny, ns),
                    "p_y": period((0, 1), nx, ny, ns),
                    "p_y_minus_x": period((-1, 1), nx, ny, ns),
                }
                rows.append(row)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    txt_path = out_prefix.with_suffix(".txt")

    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Nx", "Ny", "Ns", "p_x", "p_y", "p_y_minus_x"],
        )
        w.writeheader()
        w.writerows(rows)

    lines: list[str] = []
    lines.append("Directional periodicity scan over (Nx, Ny, Ns)")
    lines.append("============================================")
    lines.append("")
    lines.append(f"nx_list = {nx_list}")
    lines.append(f"ny_list = {ny_list}")
    if ns_list:
        lines.append(f"ns_list = {ns_list}")
    else:
        lines.append(f"ns_mode = {args.ns_mode}")
    lines.append(f"rows = {len(rows)}")
    lines.append("")
    lines.append("Columns: Nx Ny Ns p_x p_y p_y_minus_x")
    lines.append("")
    lines.append(f"{'Nx':>4} {'Ny':>4} {'Ns':>4} {'p_x':>6} {'p_y':>6} {'p_y-x':>8}")
    lines.append("-" * 42)

    for r in rows:
        lines.append(
            f"{r['Nx']:>4d} {r['Ny']:>4d} {r['Ns']:>4d} {r['p_x']:>6d} {r['p_y']:>6d} {r['p_y_minus_x']:>8d}"
        )

    txt_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {txt_path}")
    print(f"Total rows: {len(rows)}")


if __name__ == "__main__":
    main()
