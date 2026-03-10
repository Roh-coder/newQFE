#!/usr/bin/env python3
"""Verify closed-form periodicity rules for three directions.

Boundary conditions:
  (i + Nx, j) ~ (i, j)
  (i, j + Ny) ~ (i + Ns, j)

Directions:
  x         = (1, 0)
  y         = (0, 1)
  y_minus_x = (-1, 1)
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def canonicalize(x: int, y: int, nx: int, ny: int, ns: int) -> tuple[int, int]:
    y0 = y % ny
    b = (y - y0) // ny
    x0 = (x - b * ns) % nx
    return x0, y0


def period_exact(step: tuple[int, int], nx: int, ny: int, ns: int) -> int:
    sx, sy = step
    x, y = canonicalize(0, 0, nx, ny, ns)
    for n in range(1, nx * ny + 1):
        x, y = canonicalize(x + sx, y + sy, nx, ny, ns)
        if x == 0 and y == 0:
            return n
    raise RuntimeError(f"No period found for step={step}, nx={nx}, ny={ny}, ns={ns}")


def period_formula_x(nx: int, ny: int, ns: int) -> int:
    _ = ny, ns
    return nx


def period_formula_y(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ns)


def period_formula_y_minus_x(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ny + ns)


def parse_int_list(text: str) -> list[int]:
    vals: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if item:
            vals.append(int(item))
    return vals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx-list", default="1,2,3,4,5,6,7,8")
    parser.add_argument("--ny-list", default="1,2,3,4,5,6,7,8")
    parser.add_argument("--ns-list", default="1,2,3,4,5,6,7,8")
    parser.add_argument(
        "--out-prefix",
        default="thinkTwist/periodicity_rules_verification_Nx1to8_Ny1to8_Ns1to8",
    )
    args = parser.parse_args()

    nx_list = parse_int_list(args.nx_list)
    ny_list = parse_int_list(args.ny_list)
    ns_list = parse_int_list(args.ns_list)

    total = 0
    mismatches = 0
    mismatch_rows: list[dict[str, int]] = []

    for nx in nx_list:
        if nx <= 0:
            raise ValueError(f"Invalid Nx={nx}")
        for ny in ny_list:
            if ny <= 0:
                raise ValueError(f"Invalid Ny={ny}")
            for ns in ns_list:
                total += 1

                exact_x = period_exact((1, 0), nx, ny, ns)
                exact_y = period_exact((0, 1), nx, ny, ns)
                exact_d = period_exact((-1, 1), nx, ny, ns)

                form_x = period_formula_x(nx, ny, ns)
                form_y = period_formula_y(nx, ny, ns)
                form_d = period_formula_y_minus_x(nx, ny, ns)

                if (exact_x, exact_y, exact_d) != (form_x, form_y, form_d):
                    mismatches += 1
                    mismatch_rows.append(
                        {
                            "Nx": nx,
                            "Ny": ny,
                            "Ns": ns,
                            "exact_x": exact_x,
                            "formula_x": form_x,
                            "exact_y": exact_y,
                            "formula_y": form_y,
                            "exact_y_minus_x": exact_d,
                            "formula_y_minus_x": form_d,
                        }
                    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    txt_path = out_prefix.with_suffix(".txt")
    csv_path = out_prefix.with_suffix(".csv")

    lines: list[str] = []
    lines.append("Closed-form periodicity rules")
    lines.append("===========================")
    lines.append("")
    lines.append("Boundary conditions:")
    lines.append("  (i + Nx, j) ~ (i, j)")
    lines.append("  (i, j + Ny) ~ (i + Ns, j)")
    lines.append("")
    lines.append("Directions and rules:")
    lines.append("  p_x       = Nx")
    lines.append("  p_y       = Nx*Ny / gcd(Nx, Ns)")
    lines.append("  p_y_minus_x = Nx*Ny / gcd(Nx, Ny + Ns)")
    lines.append("")
    lines.append(f"nx_list = {nx_list}")
    lines.append(f"ny_list = {ny_list}")
    lines.append(f"ns_list = {ns_list}")
    lines.append(f"total_cases = {total}")
    lines.append(f"mismatches = {mismatches}")
    lines.append("status = PASS" if mismatches == 0 else "status = FAIL")

    if mismatches:
        lines.append("")
        lines.append("First mismatches:")
        lines.append("Nx Ny Ns exact_x formula_x exact_y formula_y exact_y-x formula_y-x")
        for r in mismatch_rows[:20]:
            lines.append(
                f"{r['Nx']:>2} {r['Ny']:>2} {r['Ns']:>2} "
                f"{r['exact_x']:>7} {r['formula_x']:>9} "
                f"{r['exact_y']:>7} {r['formula_y']:>9} "
                f"{r['exact_y_minus_x']:>9} {r['formula_y_minus_x']:>11}"
            )

    txt_path.write_text("\n".join(lines) + "\n")

    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "Nx",
                "Ny",
                "Ns",
                "exact_x",
                "formula_x",
                "exact_y",
                "formula_y",
                "exact_y_minus_x",
                "formula_y_minus_x",
            ],
        )
        w.writeheader()
        w.writerows(mismatch_rows)

    print(f"Wrote {txt_path}")
    print(f"Wrote {csv_path}")
    print(f"total_cases={total} mismatches={mismatches}")


if __name__ == "__main__":
    main()
