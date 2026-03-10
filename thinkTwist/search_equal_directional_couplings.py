#!/usr/bin/env python3
import argparse
import json
import math
import os
import subprocess
from pathlib import Path


def run_case(binary, nx, ny, shift, k1, k2, k3, therm, traj, skip, wolff, metro, seed, out_path):
    cmd = [
        binary,
        "--n_x", str(nx),
        "--n_y", str(ny),
        "--shift", str(shift),
        "--k1", f"{k1:.8f}",
        "--k2", f"{k2:.8f}",
        "--k3", f"{k3:.8f}",
        "--n_therm", str(therm),
        "--n_traj", str(traj),
        "--n_skip", str(skip),
        "--n_wolff", str(wolff),
        "--n_metropolis", str(metro),
        "--seed", str(seed),
        "--out", str(out_path),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_dat(path):
    meta = {"period_x": None, "period_left": None, "period_right": None}
    data = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# periods"):
                parts = line.split()
                meta["period_x"] = int(parts[3])
                meta["period_left"] = int(parts[5])
                meta["period_right"] = int(parts[7])
                continue
            if line.startswith("#"):
                continue
            p = line.split()
            d = int(p[0])
            r = int(p[1])
            m = float(p[2])
            if d in data:
                data[d][r] = m
    return meta, data


def directional_mismatch(meta, data, use_connected=True):
    dx = 3 if use_connected else 0
    dl = 4 if use_connected else 1
    dr = 5 if use_connected else 2

    px = meta["period_x"]
    pl = meta["period_left"]
    pr = meta["period_right"]
    if px is None or pl is None or pr is None:
        raise RuntimeError("Missing period metadata")

    step_l = pl // px
    step_r = pr // px

    # Match the three directions at equal fractional position over one x-period.
    rs = list(range(1, px // 2 + 1))
    sq = []
    lr_sq = []
    for rx in rs:
        rl = (rx * step_l) % pl
        rr = (rx * step_r) % pr
        gx = data[dx][rx]
        gl = data[dl][rl]
        gr = data[dr][rr]
        sq.append((gx - gl) ** 2)
        sq.append((gx - gr) ** 2)
        sq.append((gl - gr) ** 2)
        lr_sq.append((gl - gr) ** 2)

    rms3 = math.sqrt(sum(sq) / len(sq))
    rms_lr = math.sqrt(sum(lr_sq) / len(lr_sq))
    return rms3, rms_lr


def scan_grid(args, nx, ny, shift, k1_fixed, k2_vals, k3_vals, tag):
    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    run_id = 0

    for k2 in k2_vals:
        for k3 in k3_vals:
            run_id += 1
            out_path = tmp_dir / f"scan_{tag}_k2_{k2:.6f}_k3_{k3:.6f}.dat"
            run_case(
                args.binary,
                nx,
                ny,
                shift,
                k1_fixed,
                k2,
                k3,
                args.n_therm,
                args.n_traj,
                args.n_skip,
                args.n_wolff,
                args.n_metropolis,
                args.seed + run_id,
                out_path,
            )
            meta, data = parse_dat(out_path)
            rms3, rms_lr = directional_mismatch(meta, data, use_connected=True)
            rows.append(
                {
                    "Nx": nx,
                    "Ny": ny,
                    "shift": shift,
                    "k1": k1_fixed,
                    "k2": k2,
                    "k3": k3,
                    "rms3": rms3,
                    "rms_lr": rms_lr,
                    "period_x": meta["period_x"],
                    "period_left": meta["period_left"],
                    "period_right": meta["period_right"],
                }
            )
    rows.sort(key=lambda r: r["rms3"])
    return rows


def coarse_vals(center, half_width, step):
    vals = []
    x = center - half_width
    while x <= center + half_width + 1e-12:
        vals.append(round(x, 6))
        x += step
    return vals


def write_summary(path, best, top_rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps({"best": best, "top": top_rows}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="/workspaces/newQFE/bin/ising_twist_two_point_fast")
    parser.add_argument("--out-dir", default="/workspaces/newQFE/K_from_BC/results")
    parser.add_argument("--tmp-dir", default="/workspaces/newQFE/K_from_BC/results/_scan_tmp")
    parser.add_argument("--n-therm", type=int, default=1500)
    parser.add_argument("--n-traj", type=int, default=6000)
    parser.add_argument("--n-skip", type=int, default=20)
    parser.add_argument("--n-wolff", type=int, default=2)
    parser.add_argument("--n-metropolis", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260310)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        (32, 64, 0, "Nx32_Ny64_shift0"),
        (32, 32, 16, "Nx32_Ny32_shift16"),
    ]

    all_best = []

    for nx, ny, shift, tag in configs:
        # Coarse 2D scan around isotropic point with k1 fixed by scale convention.
        coarse_k2 = coarse_vals(1.0, 0.6, 0.2)
        coarse_k3 = coarse_vals(1.0, 0.6, 0.2)
        coarse_rows = scan_grid(args, nx, ny, shift, 1.0, coarse_k2, coarse_k3, tag + "_coarse")
        best_coarse = coarse_rows[0]

        # Fine scan around coarse minimum.
        c2 = best_coarse["k2"]
        c3 = best_coarse["k3"]
        fine_k2 = coarse_vals(c2, 0.12, 0.04)
        fine_k3 = coarse_vals(c3, 0.12, 0.04)
        fine_rows = scan_grid(args, nx, ny, shift, 1.0, fine_k2, fine_k3, tag + "_fine")
        best_fine = fine_rows[0]

        # Validation run with higher statistics for best candidate.
        val_path = out_dir / f"ising_twist_equalized_{tag}_best.dat"
        run_case(
            args.binary,
            nx,
            ny,
            shift,
            best_fine["k1"],
            best_fine["k2"],
            best_fine["k3"],
            therm=5000,
            traj=50000,
            skip=20,
            wolff=3,
            metro=5,
            seed=args.seed + 999,
            out_path=val_path,
        )
        meta, data = parse_dat(val_path)
        val_rms3, val_rms_lr = directional_mismatch(meta, data, use_connected=True)

        summary = {
            "tag": tag,
            "coarse_best": best_coarse,
            "fine_best": best_fine,
            "validation": {
                "path": str(val_path),
                "rms3_connected": val_rms3,
                "rms_lr_connected": val_rms_lr,
                "period_x": meta["period_x"],
                "period_left": meta["period_left"],
                "period_right": meta["period_right"],
            },
        }
        all_best.append(summary)

        write_summary(
            out_dir / f"coupling_scan_summary_{tag}.json",
            summary,
            fine_rows[:10],
        )

    print(json.dumps(all_best, indent=2))


if __name__ == "__main__":
    main()
