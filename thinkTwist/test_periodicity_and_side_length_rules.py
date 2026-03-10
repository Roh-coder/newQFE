#!/usr/bin/env python3
"""Test periodicity and side-length rules on a finite parameter grid.

Rules tested:
- p_x = Nx
- p_y = Nx*Ny / gcd(Nx, Ns)
- p_y_minus_x = Nx*Ny / gcd(Nx, Ny + Ns)

Square-target side lengths (up to scale s=1):
- l_x = 1
- l_y = p_x / p_y
- l_r = sqrt(2) * p_x / p_y_minus_x

Feasibility condition:
- |l_x - l_y| <= l_r <= l_x + l_y
"""

from __future__ import annotations

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


def p_formula_x(nx: int, ny: int, ns: int) -> int:
    _ = ny, ns
    return nx


def p_formula_y(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ns)


def p_formula_diag(nx: int, ny: int, ns: int) -> int:
    return nx * ny // math.gcd(nx, ny + ns)


def main() -> None:
    nx_vals = range(1, 9)
    ny_vals = range(1, 9)
    ns_vals = range(1, 9)

    periodicity_mismatches: list[dict[str, int]] = []
    half_orbit_mismatches: list[dict[str, float]] = []
    infeasible_examples: list[dict[str, float]] = []

    total = 0
    feasible_count = 0
    infeasible_below = 0
    infeasible_above = 0

    for nx in nx_vals:
        for ny in ny_vals:
            for ns in ns_vals:
                total += 1

                pxe = period_exact((1, 0), nx, ny, ns)
                pye = period_exact((0, 1), nx, ny, ns)
                pde = period_exact((-1, 1), nx, ny, ns)

                pxf = p_formula_x(nx, ny, ns)
                pyf = p_formula_y(nx, ny, ns)
                pdf = p_formula_diag(nx, ny, ns)

                if (pxe, pye, pde) != (pxf, pyf, pdf):
                    periodicity_mismatches.append(
                        {
                            "Nx": nx,
                            "Ny": ny,
                            "Ns": ns,
                            "exact_x": pxe,
                            "formula_x": pxf,
                            "exact_y": pye,
                            "formula_y": pyf,
                            "exact_y_minus_x": pde,
                            "formula_y_minus_x": pdf,
                        }
                    )

                # Side-length rule for square-target ratios with scale s=1.
                lx = 1.0
                ly = pxf / pyf
                lr = math.sqrt(2.0) * pxf / pdf

                # Check that half-orbit physical lengths match targets.
                # Targets are Dx=1, Dy=1, Dr=sqrt(2).
                dx = 0.5 * pxf * lx
                dy = 0.5 * pyf * ly
                dr = 0.5 * pdf * lr

                if not (
                    abs(dx - 0.5 * pxf) < 1e-12
                    and abs(dy - 0.5 * pxf) < 1e-12
                    and abs(dr - 0.5 * math.sqrt(2.0) * pxf) < 1e-12
                ):
                    half_orbit_mismatches.append(
                        {
                            "Nx": nx,
                            "Ny": ny,
                            "Ns": ns,
                            "dx": dx,
                            "dy": dy,
                            "dr": dr,
                        }
                    )

                feasible = (abs(lx - ly) <= lr <= lx + ly)
                if feasible:
                    feasible_count += 1
                else:
                    lower = abs(lx - ly)
                    upper = lx + ly
                    if lr < lower:
                        infeasible_below += 1
                    elif lr > upper:
                        infeasible_above += 1
                    if len(infeasible_examples) < 10:
                        infeasible_examples.append(
                            {
                                "Nx": nx,
                                "Ny": ny,
                                "Ns": ns,
                                "p_x": pxf,
                                "p_y": pyf,
                                "p_y_minus_x": pdf,
                                "l_y": ly,
                                "l_r": lr,
                                "lower": lower,
                                "upper": upper,
                            }
                        )

    out_dir = Path("thinkTwist")
    out_dir.mkdir(parents=True, exist_ok=True)

    txt_path = out_dir / "rules_test_report.txt"
    per_csv = out_dir / "rules_test_periodicity_mismatches.csv"
    hor_csv = out_dir / "rules_test_half_orbit_mismatches.csv"

    with per_csv.open("w", newline="") as f:
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
        w.writerows(periodicity_mismatches)

    with hor_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Nx", "Ny", "Ns", "dx", "dy", "dr"])
        w.writeheader()
        w.writerows(half_orbit_mismatches)

    lines = []
    lines.append("Rule Test Report")
    lines.append("================")
    lines.append("")
    lines.append("Grid: Nx,Ny,Ns in {1..8}")
    lines.append(f"total_cases = {total}")
    lines.append("")
    lines.append("Periodicity formula test:")
    lines.append(f"mismatches = {len(periodicity_mismatches)}")
    lines.append("status = PASS" if not periodicity_mismatches else "status = FAIL")
    lines.append("")
    lines.append("Side-length half-orbit target test (square-target ratio):")
    lines.append(f"mismatches = {len(half_orbit_mismatches)}")
    lines.append("status = PASS" if not half_orbit_mismatches else "status = FAIL")
    lines.append("")
    lines.append("Triangle feasibility count for square-target side lengths:")
    lines.append(f"feasible = {feasible_count}")
    lines.append(f"infeasible = {total - feasible_count}")
    lines.append(f"infeasible_below_lower_bound = {infeasible_below}")
    lines.append(f"infeasible_above_upper_bound = {infeasible_above}")
    lines.append("")
    lines.append("Infeasibility interpretation:")
    lines.append("For square-target side lengths with scale s=1,")
    lines.append("  l_x = 1")
    lines.append("  l_y = p_x/p_y = gcd(Nx, Ns)/Ny")
    lines.append("  l_r = sqrt(2)*p_x/p_y_minus_x = sqrt(2)*gcd(Nx, Ny+Ns)/Ny")
    lines.append("Feasibility requires |1-l_y| <= l_r <= 1+l_y.")
    lines.append("Most failures occur when l_r is too small to bridge |1-l_y|.")
    lines.append("")
    lines.append("Sample infeasible points:")
    lines.append("Nx Ny Ns p_x p_y p_y-x l_y l_r lower upper")
    for ex in infeasible_examples:
        lines.append(
            f"{int(ex['Nx']):>2} {int(ex['Ny']):>2} {int(ex['Ns']):>2} "
            f"{int(ex['p_x']):>3} {int(ex['p_y']):>3} {int(ex['p_y_minus_x']):>5} "
            f"{ex['l_y']:.5f} {ex['l_r']:.5f} {ex['lower']:.5f} {ex['upper']:.5f}"
        )

    txt_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {txt_path}")
    print(f"Wrote {per_csv}")
    print(f"Wrote {hor_csv}")


if __name__ == "__main__":
    main()
