#!/usr/bin/env python3
"""Search twisted-parallelogram integer parameters for target side ratios.

Conventions follow thinkDoubleTwist/ising_tri_twisted_parallelogram.cc:
  v = (Lx, Ty)
  u = (Tx, -Ly)
  w = -(v + u) = (-Lx - Tx, Ly - Ty)

On a triangular lattice basis e1=(1,0), e2=(1/2,sqrt(3)/2),
the squared length of integer direction (m,n) is:
  |(m,n)|^2 = m^2 + n^2 + m*n

This script finds integer tuples (Lx, Ly, Tx, Ty) such that
|v|:|u|:|w| approximates a user-provided target ratio.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Candidate:
    lx: int
    ly: int
    tx: int
    ty: int
    ncell: int
    lv2: int
    lu2: int
    lw2: int
    score: float
    mapping: tuple[str, str, str]
    achieved_ratio: tuple[float, float, float]


def tri_len_sq(m: int, n: int) -> int:
    return m * m + n * n + m * n


def compute_lengths_sq(lx: int, ly: int, tx: int, ty: int) -> tuple[int, int, int]:
    # v = (Lx, Ty), u = (Tx, -Ly), w = (-Lx-Tx, Ly-Ty)
    lv2 = tri_len_sq(lx, ty)
    lu2 = tri_len_sq(tx, -ly)
    lw2 = tri_len_sq(-lx - tx, ly - ty)
    return lv2, lu2, lw2


def ratio_score(lengths_sq: tuple[int, int, int], target: tuple[int, int, int]) -> float:
    # Compare lv2:lu2:lw2 with target^2 ratio via normalized spread of scales.
    scales = [
        lengths_sq[0] / float(target[0] * target[0]),
        lengths_sq[1] / float(target[1] * target[1]),
        lengths_sq[2] / float(target[2] * target[2]),
    ]
    mean = sum(scales) / 3.0
    if mean == 0.0:
        return float("inf")
    spread = max(scales) - min(scales)
    return spread / mean


def achieved_ratio(lengths_sq: tuple[int, int, int]) -> tuple[float, float, float]:
    lv = math.sqrt(lengths_sq[0])
    lu = math.sqrt(lengths_sq[1])
    lw = math.sqrt(lengths_sq[2])
    base = lv if lv > 0.0 else 1.0
    return (1.0, lu / base, lw / base)


def find_candidates(
    target: tuple[int, int, int],
    max_l: int,
    max_t: int,
    top_k: int,
    max_score: float,
    max_ncell: int | None,
    fixed_vuw_order: bool,
) -> list[Candidate]:
    candidates: list[Candidate] = []
    seen: set[tuple[int, int, int, int]] = set()

    side_names = ("v", "u", "w")
    target_perms = [target] if fixed_vuw_order else list(itertools.permutations(target))

    for lx in range(1, max_l + 1):
        for ly in range(1, max_l + 1):
            for tx in range(-max_t, max_t + 1):
                for ty in range(-max_t, max_t + 1):
                    ncell = lx * ly + tx * ty
                    if ncell <= 0:
                        continue
                    if max_ncell is not None and ncell > max_ncell:
                        continue

                    lv2, lu2, lw2 = compute_lengths_sq(lx, ly, tx, ty)
                    if lv2 == 0 or lu2 == 0 or lw2 == 0:
                        continue

                    best_score = float("inf")
                    best_perm: tuple[int, int, int] | None = None
                    best_map: tuple[str, str, str] | None = None

                    # Allow any assignment of target sides to (v,u,w).
                    for perm in target_perms:
                        s = ratio_score((lv2, lu2, lw2), perm)
                        if s < best_score:
                            best_score = s
                            best_perm = perm
                            best_map = (
                                f"v->{perm[0]}",
                                f"u->{perm[1]}",
                                f"w->{perm[2]}",
                            )

                    assert best_perm is not None and best_map is not None
                    if best_score > max_score:
                        continue

                    key = (lx, ly, tx, ty)
                    if key in seen:
                        continue
                    seen.add(key)

                    candidates.append(
                        Candidate(
                            lx=lx,
                            ly=ly,
                            tx=tx,
                            ty=ty,
                            ncell=ncell,
                            lv2=lv2,
                            lu2=lu2,
                            lw2=lw2,
                            score=best_score,
                            mapping=best_map,
                            achieved_ratio=achieved_ratio((lv2, lu2, lw2)),
                        )
                    )

    candidates.sort(key=lambda c: (c.score, c.ncell, abs(c.tx) + abs(c.ty), c.lx + c.ly))
    return candidates[:top_k]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Find (Lx, Ly, Tx, Ty) for target triangle side ratios.")
    ap.add_argument(
        "--ratio",
        type=int,
        nargs=3,
        default=[4, 5, 6],
        metavar=("A", "B", "C"),
        help="Target side ratio A B C (default: 4 5 6).",
    )
    ap.add_argument("--max-L", type=int, default=40, help="Max absolute value for Lx and Ly (positive scan).")
    ap.add_argument("--max-T", type=int, default=40, help="Max absolute value for Tx and Ty (signed scan).")
    ap.add_argument("--top-k", type=int, default=25, help="Number of best candidates to print.")
    ap.add_argument(
        "--max-score",
        type=float,
        default=0.02,
        help="Keep candidates with score <= max-score (0 means exact ratio in squared-length sense).",
    )
    ap.add_argument(
        "--max-ncell",
        type=int,
        default=0,
        help="Optional max Ncell = Lx*Ly + Tx*Ty; 0 disables this cut.",
    )
    ap.add_argument(
        "--fixed-vuw-order",
        action="store_true",
        help="Require |v|:|u|:|w| to match the ratio order exactly, without permutation.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    target = tuple(args.ratio)
    if min(target) <= 0:
        raise SystemExit("All ratio entries must be positive.")

    max_ncell = None if args.max_ncell <= 0 else args.max_ncell
    out = find_candidates(
        target=target,
        max_l=args.max_L,
        max_t=args.max_T,
        top_k=args.top_k,
        max_score=args.max_score,
        max_ncell=max_ncell,
        fixed_vuw_order=args.fixed_vuw_order,
    )

    print("Target side ratio:", f"{target[0]}:{target[1]}:{target[2]}")
    print("Conventions: v=(Lx,Ty), u=(Tx,-Ly), w=-(v+u)")
    print("Score: relative spread of scale factors in squared-length matching (smaller is better).")
    print()

    if not out:
        print("No candidates found in the requested search window.")
        print("Try increasing --max-L/--max-T or relaxing --max-score.")
        return

    header = (
        f"{'rank':>4}  {'Lx':>3} {'Ly':>3} {'Tx':>4} {'Ty':>4} {'Ncell':>6} "
        f"{'|v|^2':>8} {'|u|^2':>8} {'|w|^2':>8} {'score':>10}  {'map':<20} {'ratio(1:..:..)':<20}"
    )
    print(header)
    print("-" * len(header))

    for i, c in enumerate(out, start=1):
        r1, r2, r3 = c.achieved_ratio
        mapping = ",".join(c.mapping)
        print(
            f"{i:4d}  {c.lx:3d} {c.ly:3d} {c.tx:4d} {c.ty:4d} {c.ncell:6d} "
            f"{c.lv2:8d} {c.lu2:8d} {c.lw2:8d} {c.score:10.6f}  "
            f"{mapping:<20} {r1:4.2f}:{r2:4.2f}:{r3:4.2f}"
        )


if __name__ == "__main__":
    main()
