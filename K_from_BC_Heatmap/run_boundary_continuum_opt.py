#!/usr/bin/env python3
"""
run_boundary_continuum_opt.py

Continuum-calibrated boundary-correlator optimization over (r1, r2).

Workflow:
1) Precompute MC payloads for test geometries on the (r1, r2) grid, and for
    reference geometries as fixed isotropic baselines.
2) Sample boundary correlators at t_k = k/8 (k=1..7), cycles c=0,1,2.
3) Extrapolate each channel to L->infinity with quadratic fit in 1/L.
4) Build score S(r1,r2) = sum_{c,k} w(c,k) * [G_test_inf - G_ref_inf]^2.
5) Write reusable .dat files for raw channels, continuum channels, and score map.

Notes:
- Internally this reuses cached pickle payloads from precompute_grid._run_one.
- All analysis outputs are written to well-formatted .dat files.
- Reference methodology: fixed lattice geometry with isotropic couplings
    (k1=k2=k3), computed once per reference level and reused across (r1, r2).
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from precompute_grid import _point_path, _run_one  # noqa: E402
from cost import boundary_paths, _tile_interp, _SQRT3_2  # noqa: E402


def _write_dat(path, header_lines, columns, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for h in header_lines:
            f.write(f"# {h}\n")
        f.write("# columns: " + " ".join(columns) + "\n")
        for r in rows:
            out = []
            for v in r:
                if isinstance(v, (int, np.integer)):
                    out.append(str(int(v)))
                elif isinstance(v, (float, np.floating)):
                    out.append(f"{float(v):.10g}")
                else:
                    out.append(str(v))
            f.write(" ".join(out) + "\n")


def _parse_ref_geom_entries(entries):
    """Parse repeated --ref-geom entries of form 'L:Lx:Ly:Tx:Ty'."""
    out = {}
    for s in entries:
        parts = s.split(":")
        if len(parts) != 5:
            raise ValueError(f"bad --ref-geom '{s}', expected L:Lx:Ly:Tx:Ty")
        L, Lx, Ly, Tx, Ty = map(int, parts)
        out[L] = (Lx, Ly, Tx, Ty)
    return out


def _sample_boundary_from_payload(payload, t_fracs):
    Lx, Ly = payload["Lx"], payload["Ly"]
    Tx, Ty = payload["Tx"], payload["Ty"]
    data = payload["data"]

    iG = _tile_interp(data, Lx, Ly, Tx, Ty, "conn", copies=2)
    iE = _tile_interp(data, Lx, Ly, Tx, Ty, "conn_err", copies=2)

    paths = boundary_paths(Lx, Ly, Tx, Ty)
    G = np.full((3, len(t_fracs)), np.nan)
    sG = np.full_like(G, np.nan)

    for c, (dm, dn) in enumerate(paths):
        ex, ey = dm + 0.5 * dn, _SQRT3_2 * dn
        pts = np.column_stack([t_fracs * ex, t_fracs * ey])
        G[c] = np.asarray(iG(pts), dtype=float).ravel()
        sG[c] = np.abs(np.asarray(iE(pts), dtype=float).ravel())
    return G, sG


def _fit_quad_in_invL(Larr, y, sigma):
    """Weighted fit y = a + b*(1/L) + c*(1/L)^2. Returns (a, sa, b, c, n_used)."""
    x = 1.0 / np.asarray(Larr, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.asarray(sigma, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    if np.count_nonzero(m) < 3:
        return np.nan, np.nan, np.nan, np.nan, int(np.count_nonzero(m))

    xv, yv, sv = x[m], y[m], s[m]
    X = np.column_stack([np.ones_like(xv), xv, xv * xv])
    w = 1.0 / (sv * sv)
    XtW = X.T * w
    try:
        cov = np.linalg.inv(XtW @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan, int(np.count_nonzero(m))

    beta = cov @ (XtW @ yv)
    a, b, c = float(beta[0]), float(beta[1]), float(beta[2])
    sa = float(np.sqrt(max(cov[0, 0], 0.0)))
    return a, sa, b, c, int(np.count_nonzero(m))


def _geom_test(L):
    return (int(L), int(L), 0, 0)


def _resolve_exe(path):
    if path:
        return path
    base = os.path.join(_HERE, "bin", "ising_tri_twisted_parallelogram")
    cands = (base + ".exe", base) if os.name == "nt" else (base, base + ".exe")
    for c in cands:
        if os.path.exists(c):
            return c
    raise FileNotFoundError("simulator binary not found in bin/")


def _build_jobs(args, out_root, exe, rs, ref_geom_map):
    jobs = []
    points = [(float(r1), float(r2)) for r1 in rs for r2 in rs]
    scratch = os.path.join(out_root, "_mc_scratch")
    os.makedirs(scratch, exist_ok=True)

    for L in args.test_sizes:
        Lx, Ly, Tx, Ty = _geom_test(L)
        grid_dir = os.path.join(out_root, "grid", f"L{L}", "test")
        os.makedirs(grid_dir, exist_ok=True)
        for r1, r2 in points:
            out_pkl = _point_path(grid_dir, r1, r2)
            label = f"test_L{L}_r1{r1:.3f}_r2{r2:.3f}"
            jobs.append((exe, label, Lx, Ly, Tx, Ty, r1, r2,
                         args.n_traj, args.n_traj_coarse, args.n_traj_fine,
                         args.beta_seed[0], args.beta_seed[1],
                         scratch, out_pkl))

    ref_r1 = float(args.ref_fixed_r1)
    ref_r2 = float(args.ref_fixed_r2)

    if args.reference_mode == "continuum":
        for L in args.ref_sizes:
            if L not in ref_geom_map:
                raise ValueError(f"Missing ref geometry for L={L}. Use --ref-geom.")
            Lx, Ly, Tx, Ty = ref_geom_map[L]
            grid_dir = os.path.join(out_root, "grid", f"L{L}", "ref")
            os.makedirs(grid_dir, exist_ok=True)
            out_pkl = _point_path(grid_dir, ref_r1, ref_r2)
            label = f"ref_L{L}_fixed_r1{ref_r1:.3f}_r2{ref_r2:.3f}"
            jobs.append((exe, label, Lx, Ly, Tx, Ty, ref_r1, ref_r2,
                         args.n_traj, args.n_traj_coarse, args.n_traj_fine,
                         args.beta_seed[0], args.beta_seed[1],
                         scratch, out_pkl))
    else:
        Lx, Ly, Tx, Ty = args.ref_large
        grid_dir = os.path.join(out_root, "grid", "Lref_large", "ref")
        os.makedirs(grid_dir, exist_ok=True)
        out_pkl = _point_path(grid_dir, ref_r1, ref_r2)
        label = f"ref_large_fixed_r1{ref_r1:.3f}_r2{ref_r2:.3f}"
        jobs.append((exe, label, Lx, Ly, Tx, Ty, ref_r1, ref_r2,
                     args.n_traj, args.n_traj_coarse, args.n_traj_fine,
                     args.beta_seed[0], args.beta_seed[1],
                     scratch, out_pkl))

    return jobs, points


def _run_jobs(jobs, workers):
    n_total = len(jobs)
    n_cached = sum(1 for j in jobs if os.path.exists(j[-1]))
    print(f"[precompute] total={n_total} cached={n_cached} todo={n_total-n_cached} workers={workers}")

    ok = skip = err = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_run_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label, status, info = fut.result()
            if status == "ok":
                ok += 1
                print(f"[{i:4d}/{n_total}] ok   {label}  beta_c={info['beta_c']:.5f} wall={info['wall']:.1f}s")
            elif status == "skip":
                skip += 1
            else:
                err += 1
                print(f"[{i:4d}/{n_total}] ERR  {label}  {info.get('msg','?')}")
    print(f"[precompute] done in {time.time()-t0:.1f}s  ok={ok} skip={skip} err={err}")
    return err


def _score_heatmap_plot(rs, score_grid, out_png, title):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(score_grid, origin="lower", aspect="auto",
                   extent=[rs.min(), rs.max(), rs.min(), rs.max()])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("score S(r1,r2)")
    j, i = np.unravel_index(np.nanargmin(score_grid), score_grid.shape)
    r1_min = rs[i]
    r2_min = rs[j]
    ax.plot(r1_min, r2_min, marker="*", markersize=14,
            markerfacecolor="white", markeredgecolor="k")
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _score_and_zscore_heatmap_plot(rs, score_grid, zscore_grid, out_png, title):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.8), sharex=True, sharey=True)

    im1 = ax1.imshow(score_grid, origin="lower", aspect="auto",
                     extent=[rs.min(), rs.max(), rs.min(), rs.max()])
    cb1 = fig.colorbar(im1, ax=ax1)
    cb1.set_label("score S(r1,r2)")
    if np.any(np.isfinite(score_grid)):
        j, i = np.unravel_index(np.nanargmin(score_grid), score_grid.shape)
        ax1.plot(rs[i], rs[j], marker="*", markersize=12,
                 markerfacecolor="white", markeredgecolor="k")
    ax1.set_title("Score Heatmap")
    ax1.set_xlabel("r1")
    ax1.set_ylabel("r2")

    im2 = ax2.imshow(zscore_grid, origin="lower", aspect="auto",
                     extent=[rs.min(), rs.max(), rs.min(), rs.max()])
    cb2 = fig.colorbar(im2, ax=ax2)
    cb2.set_label("RMS z-score")
    if np.any(np.isfinite(zscore_grid)):
        jz, iz = np.unravel_index(np.nanargmin(zscore_grid), zscore_grid.shape)
        ax2.plot(rs[iz], rs[jz], marker="*", markersize=12,
                 markerfacecolor="white", markeredgecolor="k")
    ax2.set_title("Test-vs-Reference RMS z-score")
    ax2.set_xlabel("r1")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", type=str, default="boundary_opt")
    ap.add_argument("--exe", type=str, default=None)

    ap.add_argument("--r-min", type=float, default=0.5)
    ap.add_argument("--r-max", type=float, default=1.5)
    ap.add_argument("--r-step", type=float, default=0.25)

    ap.add_argument("--test-sizes", type=int, nargs="+",
                    default=[8, 16, 24, 32, 40])
    ap.add_argument("--reference-mode", choices=["continuum", "large"], default="continuum")
    ap.add_argument("--ref-sizes", type=int, nargs="+", default=[8, 16, 24, 32, 40])
    ap.add_argument("--ref-geom", action="append", default=[],
                    help="Reference geometry map entry: L:Lx:Ly:Tx:Ty (repeatable)")
    ap.add_argument("--ref-large", type=int, nargs=4, default=[96, 96, 24, 24],
                    metavar=("Lx", "Ly", "Tx", "Ty"))

    ap.add_argument("--n-traj", type=int, default=20000)
    ap.add_argument("--n-traj-coarse", type=int, default=2000)
    ap.add_argument("--n-traj-fine", type=int, default=4000)
    ap.add_argument("--beta-seed", type=float, nargs=2, default=[0.05, 0.40])
    ap.add_argument("--n-workers", type=int, default=4)

    ap.add_argument("--ref-fixed-r1", type=float, default=1.0,
                    help="Fixed isotropic r1 used for reference payload generation")
    ap.add_argument("--ref-fixed-r2", type=float, default=1.0,
                    help="Fixed isotropic r2 used for reference payload generation")

    ap.add_argument("--weighted", action="store_true",
                    help="Use inverse-variance weighting in score sum")

    args = ap.parse_args()
    args.ref_large = tuple(int(v) for v in args.ref_large)

    rs = np.arange(args.r_min, args.r_max + 1e-12, args.r_step)
    ref_geom_map = _parse_ref_geom_entries(args.ref_geom)

    # Default reference geometry mapping if none provided.
    if args.reference_mode == "continuum" and not ref_geom_map:
        for L in args.ref_sizes:
            T = int(round(0.25 * L))
            ref_geom_map[L] = (int(L), int(L), T, T)

    exe = _resolve_exe(args.exe)
    out_root = os.path.join(_HERE, "results", args.tag)
    dat_dir = os.path.join(out_root, "dat")
    os.makedirs(dat_dir, exist_ok=True)

    jobs, points = _build_jobs(args, out_root, exe, rs, ref_geom_map)
    err = _run_jobs(jobs, args.n_workers)
    if err > 0:
        print("[fatal] precompute errors encountered; aborting analysis")
        return 2

    k_values = np.arange(1, 8, dtype=int)
    t_fracs = k_values.astype(float) / 8.0

    raw_test_rows = []
    raw_ref_rows = []
    cont_test_rows = []
    cont_ref_rows = []
    score_rows = []
    zscore_rows = []

    score_grid = np.full((len(rs), len(rs)), np.nan)
    zscore_grid = np.full((len(rs), len(rs)), np.nan)

    # Build quick lookup from r-value to index in the mesh.
    r_index = {round(float(v), 12): i for i, v in enumerate(rs)}

    for r1, r2 in points:
        # -------- collect test channels --------
        test_by_L = {}
        for L in args.test_sizes:
            pkl = _point_path(os.path.join(out_root, "grid", f"L{L}", "test"), r1, r2)
            if not os.path.exists(pkl):
                continue
            with open(pkl, "rb") as f:
                payload = pickle.load(f)
            G, sG = _sample_boundary_from_payload(payload, t_fracs)
            test_by_L[L] = (payload, G, sG)
            for c in range(3):
                for ik, k in enumerate(k_values):
                    t = float(t_fracs[ik])
                    raw_test_rows.append([
                        r1, r2, L, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                        payload["beta_c"], c, int(k), t, G[c, ik], sG[c, ik]
                    ])

        # -------- collect reference channels --------
        ref_by_L = {}
        if args.reference_mode == "continuum":
            for L in args.ref_sizes:
                pkl = _point_path(
                    os.path.join(out_root, "grid", f"L{L}", "ref"),
                    float(args.ref_fixed_r1),
                    float(args.ref_fixed_r2),
                )
                if not os.path.exists(pkl):
                    continue
                with open(pkl, "rb") as f:
                    payload = pickle.load(f)
                G, sG = _sample_boundary_from_payload(payload, t_fracs)
                ref_by_L[L] = (payload, G, sG)
                for c in range(3):
                    for ik, k in enumerate(k_values):
                        t = float(t_fracs[ik])
                        raw_ref_rows.append([
                            r1, r2, L, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                            payload["beta_c"], c, int(k), t, G[c, ik], sG[c, ik]
                        ])
        else:
            pkl = _point_path(
                os.path.join(out_root, "grid", "Lref_large", "ref"),
                float(args.ref_fixed_r1),
                float(args.ref_fixed_r2),
            )
            if os.path.exists(pkl):
                with open(pkl, "rb") as f:
                    payload = pickle.load(f)
                G, sG = _sample_boundary_from_payload(payload, t_fracs)
                ref_by_L["large"] = (payload, G, sG)
                for c in range(3):
                    for ik, k in enumerate(k_values):
                        t = float(t_fracs[ik])
                        raw_ref_rows.append([
                            r1, r2, -1, payload["Lx"], payload["Ly"], payload["Tx"], payload["Ty"],
                            payload["beta_c"], c, int(k), t, G[c, ik], sG[c, ik]
                        ])

        # -------- continuum channels --------
        Gt_inf = np.full((3, len(t_fracs)), np.nan)
        sGt_inf = np.full_like(Gt_inf, np.nan)
        Gr_inf = np.full_like(Gt_inf, np.nan)
        sGr_inf = np.full_like(Gt_inf, np.nan)

        for c in range(3):
            for ik, k in enumerate(k_values):
                t = float(t_fracs[ik])
                # test fit
                Lt = sorted(test_by_L)
                yt = np.array([test_by_L[L][1][c, ik] for L in Lt], dtype=float)
                st = np.array([test_by_L[L][2][c, ik] for L in Lt], dtype=float)
                at, sat, bt, ct, nt = _fit_quad_in_invL(np.array(Lt, dtype=float), yt, st)
                Gt_inf[c, ik], sGt_inf[c, ik] = at, sat
                cont_test_rows.append([r1, r2, c, int(k), t, at, sat, bt, ct, nt])

                # reference fit or large-L pass-through
                if args.reference_mode == "continuum":
                    Lr = sorted(ref_by_L)
                    yr = np.array([ref_by_L[L][1][c, ik] for L in Lr], dtype=float)
                    sr = np.array([ref_by_L[L][2][c, ik] for L in Lr], dtype=float)
                    ar, sar, br, cr, nr = _fit_quad_in_invL(np.array(Lr, dtype=float), yr, sr)
                    Gr_inf[c, ik], sGr_inf[c, ik] = ar, sar
                    cont_ref_rows.append([r1, r2, c, int(k), t, ar, sar, br, cr, nr])
                else:
                    if "large" in ref_by_L:
                        ar = float(ref_by_L["large"][1][c, ik])
                        sar = float(ref_by_L["large"][2][c, ik])
                    else:
                        ar, sar = np.nan, np.nan
                    Gr_inf[c, ik], sGr_inf[c, ik] = ar, sar
                    cont_ref_rows.append([r1, r2, c, int(k), t, ar, sar, 0.0, 0.0, 1])

        # -------- score --------
        score = 0.0
        n_used = 0
        z2_sum = 0.0
        n_z = 0
        for c in range(3):
            for ik in range(len(t_fracs)):
                gt, st = Gt_inf[c, ik], sGt_inf[c, ik]
                gr, sr = Gr_inf[c, ik], sGr_inf[c, ik]
                if not (np.isfinite(gt) and np.isfinite(gr)):
                    continue
                diff = gt - gr
                w = 1.0
                if args.weighted and np.isfinite(st) and np.isfinite(sr) and (st > 0 or sr > 0):
                    var = st * st + sr * sr
                    if var > 0:
                        w = 1.0 / var
                score += w * diff * diff
                n_used += 1

                if np.isfinite(st) and np.isfinite(sr):
                    var = st * st + sr * sr
                    if var > 0:
                        z = diff / np.sqrt(var)
                        z2_sum += z * z
                        n_z += 1

        if n_used == 0:
            score = np.nan
        score_rows.append([r1, r2, score, n_used])

        z_rms = np.sqrt(z2_sum / n_z) if n_z > 0 else np.nan
        zscore_rows.append([r1, r2, z_rms, n_z])

        i = r_index[round(float(r1), 12)]
        j = r_index[round(float(r2), 12)]
        score_grid[j, i] = score
        zscore_grid[j, i] = z_rms

    # -------- write dat outputs --------
    hdr = [
        "Continuum boundary optimization outputs",
        f"tag={args.tag}",
        f"reference_mode={args.reference_mode}",
        f"test_sizes={args.test_sizes}",
        f"ref_sizes={args.ref_sizes if args.reference_mode=='continuum' else 'large_only'}",
        f"r_grid=[{args.r_min},{args.r_max}] step={args.r_step}",
        "channels: cycle c=0,1,2 and t_k=k/8 with k=1..7",
    ]

    _write_dat(
        os.path.join(dat_dir, "raw_test_channels.dat"),
        hdr,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_test_rows,
    )

    _write_dat(
        os.path.join(dat_dir, "raw_ref_channels.dat"),
        hdr,
        ["r1", "r2", "L", "Lx", "Ly", "Tx", "Ty", "beta_c", "cycle", "k", "t", "G", "sigma_G"],
        raw_ref_rows,
    )

    _write_dat(
        os.path.join(dat_dir, "continuum_test.dat"),
        hdr,
        ["r1", "r2", "cycle", "k", "t", "G_inf", "sigma_inf", "b_1overL", "c_1overL2", "n_used"],
        cont_test_rows,
    )

    _write_dat(
        os.path.join(dat_dir, "continuum_ref.dat"),
        hdr,
        ["r1", "r2", "cycle", "k", "t", "G_inf", "sigma_inf", "b_1overL", "c_1overL2", "n_used"],
        cont_ref_rows,
    )

    _write_dat(
        os.path.join(dat_dir, "score_map.dat"),
        hdr + [f"weighted={args.weighted}"],
        ["r1", "r2", "score", "n_channels_used"],
        score_rows,
    )

    _write_dat(
        os.path.join(dat_dir, "zscore_map.dat"),
        hdr,
        ["r1", "r2", "zscore_rms", "n_channels_used"],
        zscore_rows,
    )

    if np.any(np.isfinite(score_grid)):
        jj, ii = np.unravel_index(np.nanargmin(score_grid), score_grid.shape)
        r1_min, r2_min = rs[ii], rs[jj]
        s_min = score_grid[jj, ii]
        _write_dat(
            os.path.join(dat_dir, "score_minimum.dat"),
            hdr,
            ["r1_min", "r2_min", "score_min"],
            [[r1_min, r2_min, s_min]],
        )
        _score_heatmap_plot(
            rs,
            score_grid,
            os.path.join(out_root, "score_heatmap.png"),
            f"Boundary continuum score map  (tag={args.tag})",
        )
        _score_and_zscore_heatmap_plot(
            rs,
            score_grid,
            zscore_grid,
            os.path.join(out_root, "score_and_zscore_heatmaps.png"),
            f"Boundary continuum diagnostics  (tag={args.tag})",
        )
        print(f"[done] min score at r1={r1_min:.4f}, r2={r2_min:.4f}, S={s_min:.6g}")
    else:
        print("[warn] score map has no finite entries")

    manifest = {
        "tag": args.tag,
        "reference_mode": args.reference_mode,
        "test_sizes": args.test_sizes,
        "ref_sizes": args.ref_sizes if args.reference_mode == "continuum" else [],
        "ref_large": list(args.ref_large),
        "ref_geom_map": {str(k): list(v) for k, v in ref_geom_map.items()},
        "r_min": args.r_min,
        "r_max": args.r_max,
        "r_step": args.r_step,
        "n_traj": args.n_traj,
        "n_traj_coarse": args.n_traj_coarse,
        "n_traj_fine": args.n_traj_fine,
        "weighted": args.weighted,
        "ref_fixed_r1": float(args.ref_fixed_r1),
        "ref_fixed_r2": float(args.ref_fixed_r2),
        "dat_dir": dat_dir,
    }
    with open(os.path.join(out_root, "manifest_boundary_opt.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"[done] dat outputs in {dat_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
