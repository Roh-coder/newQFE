#!/usr/bin/env python3
"""
run.py — single-entry orchestrator for the K_from_BC_Heatmap workflow.

Wraps three stages (precompute → cost-replay → plot) into one script.
Every knob is exposed as both a CONFIG dict (edit at the top of this
file) and a CLI flag.  CLI values override CONFIG.

Stages
------
  precompute : run MC for test and ref geometries at each (L, r1, r2),
               find β_c, and dump pickles.
  compute    : replay cost(r1,r2;L) and fit cost(L) vs 1/L → cost_inf.
  plot       : heatmaps per L + continuum + sigma.

Use --stages to select a subset:
  --stages precompute              # long overnight MC run
  --stages compute plot            # replay / re-plot with no new MC

Geometry parameterisation
-------------------------
Test and reference lattice shapes are independently specified as
multipliers (Lx, Ly) and fractions (Tx, Ty) of the FSS scale L:

  Lx = round(Lx_mult * L)          Ly = round(Ly_mult * L)
  Tx = round(Tx_frac * L)          Ty = round(Ty_frac * L)

Defaults produce the standard untwisted-vs-twisted comparison:
  test  : Lx_mult=1, Ly_mult=1, Tx_frac=0,    Ty_frac=0
  ref   : Lx_mult=1, Ly_mult=1, Tx_frac=0.25, Ty_frac=0.25

For a sanity check (cost → 0) set all ref fracs = test fracs.
For a non-square ref (e.g. 4:3 aspect), set ref_Ly_mult=0.75.

Progress output
---------------
Timestamped + elapsed per job.  Periodic heartbeat every
--status-interval seconds: done/total, ok/err, mean wall, ETA.

    nohup python3 run.py --tag prod > run_prod.log 2>&1 &
    tail -f run_prod.log
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine                                    # noqa: E402
import cost as cost_module                          # noqa: E402
from precompute_grid import _point_path, _run_one   # noqa: E402
from compute_landscape import _fit_continuum        # noqa: E402


# ===========================================================================
# DEFAULT CONFIG — edit here, or override via CLI flags
# ===========================================================================
CONFIG = {
    # --- output / bookkeeping -----------------------------------------------
    "tag":    "prod",   # results/<tag>/
    "exe":    None,     # None → bin/ising_tri_twisted_parallelogram

    # which stages to run
    "stages": ["precompute", "compute", "plot"],

    # --- coupling grid (r1, r2) ---------------------------------------------
    "r_min":  0.5,
    "r_max":  3.0,
    "r_step": 0.25,

    # --- FSS size levels (≥2 for a continuum fit) ---------------------------
    "sizes":  [12, 16, 24],

    # --- test-lattice geometry: (Lx,Ly,Tx,Ty) = (round(mult*L), ...) -------
    "test_Lx_mult": 1.0,   # Lx = round(mult * L)
    "test_Ly_mult": 1.0,   # Ly = round(mult * L)
    "test_Tx_frac": 0.0,   # Tx = round(frac * L)  [0 = untwisted]
    "test_Ty_frac": 0.0,   # Ty = round(frac * L)

    # --- reference-lattice geometry -----------------------------------------
    "ref_Lx_mult": 1.0,
    "ref_Ly_mult": 1.0,
    "ref_Tx_frac": 0.25,   # default moderate twist
    "ref_Ty_frac": 0.25,

    # --- Monte Carlo --------------------------------------------------------
    "n_traj":        50000,      # production trajectories per grid point
    "n_traj_coarse": 2000,       # beta_c finder coarse-pass trajectories
    "n_traj_fine":   4000,       # beta_c finder fine-pass trajectories
    "n_skip":        20,         # measurement interval (sweeps)
    "n_therm":       2000,       # thermalisation sweeps before measurements
    "beta_seed":     (0.05, 0.40),  # initial beta bracket [lo, hi]
    "n_workers":     4,          # parallel MC processes

    # --- cost replay --------------------------------------------------------
    "cost_modes":  ["test_native"],
    # available: test_native | huber_log | affine_fit | spectral | l4mean_both_interp
    "cost_power":  2,            # residual exponent for test_native (2 or 4)
    "fit":         "linear",     # continuum extrapolation: linear | quadratic

    # --- plotting -----------------------------------------------------------
    "plot_log":  False,    # render log10(cost)
    "plot_vmax": None,     # cap colour scale (None → 99th percentile)

    # --- monitoring ---------------------------------------------------------
    "status_interval": 60.0,   # seconds between heartbeat lines
}


# ===========================================================================
# Geometry helpers
# ===========================================================================

def _make_geom(L, lx_mult, ly_mult, tx_frac, ty_frac):
    """Return (Lx, Ly, Tx, Ty) for scale L and geometry parameters."""
    return (int(round(lx_mult * L)),
            int(round(ly_mult * L)),
            int(round(tx_frac * L)),
            int(round(ty_frac * L)))


def _geom_label(Lx, Ly, Tx, Ty):
    return f"({Lx},{Ly},{Tx:+d},{Ty:+d})"


# ===========================================================================
# Pretty-print helpers
# ===========================================================================
_T0_GLOBAL = time.time()


def _ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _elapsed(t0=None):
    s = time.time() - (t0 if t0 is not None else _T0_GLOBAL)
    return str(timedelta(seconds=int(s)))


def log(msg):
    print(f"[{_ts()}  +{_elapsed()}]  {msg}", flush=True)


def banner(title):
    bar = "=" * 72
    print(f"\n{bar}\n[{_ts()}  +{_elapsed()}]  {title}\n{bar}", flush=True)


def fmt_eta(t_per_job, remaining):
    if t_per_job <= 0 or remaining <= 0:
        return "—"
    return str(timedelta(seconds=int(t_per_job * remaining)))


# ===========================================================================
# Stage 1 — precompute (MC)
# ===========================================================================

def stage_precompute(cfg):
    banner(f"STAGE 1/3 — PRECOMPUTE  tag={cfg['tag']}")
    exe = cfg["exe"] or os.path.join(_HERE, "bin",
                                     "ising_tri_twisted_parallelogram")
    if not os.path.exists(exe):
        log(f"ERROR: simulator not found: {exe}")
        sys.exit(1)
    log(f"simulator: {exe}")

    out_root     = os.path.join(_HERE, "results", cfg["tag"])
    scratch_root = os.path.join(out_root, "_mc_scratch")
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(scratch_root, exist_ok=True)

    rs     = np.arange(cfg["r_min"], cfg["r_max"] + 1e-9, cfg["r_step"])
    points = [(float(r1), float(r2)) for r1 in rs for r2 in rs]

    # Build per-L geometry table and log it.
    geom_table = {}
    for L in cfg["sizes"]:
        tg = _make_geom(L, cfg["test_Lx_mult"], cfg["test_Ly_mult"],
                        cfg["test_Tx_frac"], cfg["test_Ty_frac"])
        rg = _make_geom(L, cfg["ref_Lx_mult"],  cfg["ref_Ly_mult"],
                        cfg["ref_Tx_frac"],  cfg["ref_Ty_frac"])
        geom_table[L] = {"test": tg, "ref": rg}
        log(f"  L={L:3d}  test={_geom_label(*tg)}  ref={_geom_label(*rg)}")

    jobs = []
    for L in cfg["sizes"]:
        for kind in ("test", "ref"):
            Lx, Ly, Tx, Ty = geom_table[L][kind]
            grid_dir = os.path.join(out_root, "grid", f"L{L}", kind)
            os.makedirs(grid_dir, exist_ok=True)
            for r1, r2 in points:
                out_pkl = _point_path(grid_dir, r1, r2)
                label   = f"L{L}_{kind}_r1{r1:.3f}_r2{r2:.3f}"
                jobs.append((exe, label, Lx, Ly, Tx, Ty, r1, r2,
                              cfg["n_traj"],
                              cfg["n_traj_coarse"],
                              cfg["n_traj_fine"],
                              cfg["beta_seed"][0],
                              cfg["beta_seed"][1],
                              scratch_root, out_pkl))

    n_total    = len(jobs)
    n_existing = sum(1 for j in jobs if os.path.exists(j[-1]))
    n_to_do    = n_total - n_existing
    log(f"plan: {len(points)} (r1,r2) × {len(cfg['sizes'])} L × 2 geoms "
        f"= {n_total} jobs ({n_existing} cached, {n_to_do} to run)")
    log(f"plan: n_traj_prod={cfg['n_traj']}  "
        f"scan=({cfg['n_traj_coarse']},{cfg['n_traj_fine']})  "
        f"workers={cfg['n_workers']}")
    log(f"plan: beta_seed=[{cfg['beta_seed'][0]},{cfg['beta_seed'][1]}]  "
        f"n_skip={cfg['n_skip']}  n_therm={cfg['n_therm']}")
    log(f"plan: out → {out_root}")

    if n_to_do == 0:
        log("nothing to do — all points already cached.")
        _write_manifest(out_root, cfg, geom_table,
                        n_total, 0, n_total, 0, 0.0)
        return {"out_root": out_root, "n_ok": 0,
                "n_skip": n_existing, "n_err": 0}

    t_start       = time.time()
    t_last_status = t_start
    n_ok = n_skip = n_err = 0
    completed = 0
    walls = []

    with ProcessPoolExecutor(max_workers=cfg["n_workers"]) as ex:
        futs = {ex.submit(_run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            label, status, info = fut.result()
            completed += 1
            if status == "ok":
                n_ok += 1
                walls.append(info["wall"])
                log(f"[{completed:4d}/{n_total}] ok   {label}  "
                    f"beta_c={info['beta_c']:.5f}  wall={info['wall']:.1f}s")
            elif status == "skip":
                n_skip += 1
            else:
                n_err += 1
                log(f"[{completed:4d}/{n_total}] ERR  {label}  "
                    f"{info.get('msg', '?')}")

            now = time.time()
            if now - t_last_status >= cfg["status_interval"]:
                done_now   = n_ok + n_err
                remaining  = n_to_do - done_now
                mean_wall  = sum(walls) / len(walls) if walls else 0.0
                throughput = mean_wall / max(cfg["n_workers"], 1)
                pct        = 100.0 * done_now / max(n_to_do, 1)
                log(f"--- status  done={done_now}/{n_to_do} ({pct:.1f}%)  "
                    f"ok={n_ok}  err={n_err}  "
                    f"mean_wall={mean_wall:.1f}s  "
                    f"throughput~{throughput:.1f}s/job  "
                    f"ETA~{fmt_eta(throughput, remaining)} ---")
                t_last_status = now

    wall = time.time() - t_start
    _write_manifest(out_root, cfg, geom_table,
                    n_total, n_ok, n_skip, n_err, wall)
    shutil.rmtree(scratch_root, ignore_errors=True)
    log(f"DONE precompute: wall={wall:.1f}s  "
        f"ok={n_ok}  skip={n_skip}  err={n_err}")
    log(f"manifest -> {os.path.join(out_root, 'manifest.json')}")
    return {"out_root": out_root,
            "n_ok": n_ok, "n_skip": n_skip, "n_err": n_err}


def _write_manifest(out_root, cfg, geom_table,
                    n_total, n_ok, n_skip, n_err, wall):
    m = {
        "tag":   cfg["tag"],
        "sizes": cfg["sizes"],
        # all 8 geometry parameters
        "test_Lx_mult": cfg["test_Lx_mult"],
        "test_Ly_mult": cfg["test_Ly_mult"],
        "test_Tx_frac": cfg["test_Tx_frac"],
        "test_Ty_frac": cfg["test_Ty_frac"],
        "ref_Lx_mult":  cfg["ref_Lx_mult"],
        "ref_Ly_mult":  cfg["ref_Ly_mult"],
        "ref_Tx_frac":  cfg["ref_Tx_frac"],
        "ref_Ty_frac":  cfg["ref_Ty_frac"],
        # expanded per-L geometries for quick inspection
        "test_geoms": {str(L): list(geom_table[L]["test"])
                       for L in cfg["sizes"]},
        "ref_geoms":  {str(L): list(geom_table[L]["ref"])
                       for L in cfg["sizes"]},
        # r-grid
        "r_min":  cfg["r_min"],
        "r_max":  cfg["r_max"],
        "r_step": cfg["r_step"],
        # MC parameters
        "n_traj_prod":         cfg["n_traj"],
        "n_traj_scan_coarse":  cfg["n_traj_coarse"],
        "n_traj_scan_fine":    cfg["n_traj_fine"],
        "n_skip":   cfg["n_skip"],
        "n_therm":  cfg["n_therm"],
        "beta_seed": list(cfg["beta_seed"]),
        # run stats
        "wall_s":          round(wall, 1),
        "n_total":         n_total,
        "n_ok_this_run":   n_ok,
        "n_skip_this_run": n_skip,
        "n_err_this_run":  n_err,
    }
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(m, f, indent=2)


# ===========================================================================
# Stage 2 — compute cost & continuum extrapolation
# ===========================================================================

def stage_compute(cfg, mode):
    banner(f"STAGE 2/3 — COMPUTE  tag={cfg['tag']}  "
           f"mode={mode}  fit={cfg['fit']}")
    root  = os.path.join(_HERE, "results", cfg["tag"])
    mpath = os.path.join(root, "manifest.json")
    if not os.path.exists(mpath):
        log(f"ERROR: manifest missing: {mpath}")
        sys.exit(1)
    with open(mpath) as f:
        manifest = json.load(f)
    sizes = list(manifest["sizes"])

    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n  = len(rs)
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")

    COST_L  = {L: np.full((n, n), np.nan) for L in sizes}
    SIGMA_L = {L: np.full((n, n), np.nan) for L in sizes}
    BETAC_T = {L: np.full((n, n), np.nan) for L in sizes}
    BETAC_R = {L: np.full((n, n), np.nan) for L in sizes}
    COST_INF  = np.full((n, n), np.nan)
    SIGMA_INF = np.full((n, n), np.nan)
    SLOPE_INF = np.full((n, n), np.nan)
    OK_INF    = np.zeros((n, n), dtype=bool)

    fit_records = []
    n_full = n_partial = n_missing = 0
    t0          = time.time()
    last_status = t0
    total_pts   = n * n
    done_pts    = 0

    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            done_pts += 1
            costs_per_L = []
            sigmas_per_L = []
            Ls_used = []
            for L in sizes:
                test_pkl = _point_path(
                    os.path.join(root, "grid", f"L{L}", "test"), r1, r2)
                ref_pkl  = _point_path(
                    os.path.join(root, "grid", f"L{L}", "ref"),  r1, r2)
                if not (os.path.exists(test_pkl) and
                        os.path.exists(ref_pkl)):
                    continue
                with open(test_pkl, "rb") as f:
                    pt = pickle.load(f)
                with open(ref_pkl, "rb") as f:
                    pr = pickle.load(f)
                BETAC_T[L][i, j] = pt["beta_c"]
                BETAC_R[L][i, j] = pr["beta_c"]
                try:
                    c, s, _pd, _pds = cost_module.l2_cost(
                        pr["data"], pt["data"],
                        pt["Lx"], pt["Ly"], pt["Tx"], pt["Ty"],
                        pr["Lx"], pr["Ly"], pr["Tx"], pr["Ty"],
                        cost_mode=mode,
                        cost_power=cfg["cost_power"],
                    )
                except Exception as e:  # noqa: BLE001
                    log(f"  [skip cost] L={L} r=({r1:.2f},{r2:.2f}): {e}")
                    continue
                COST_L[L][i, j]  = c
                SIGMA_L[L][i, j] = s
                costs_per_L.append(c)
                sigmas_per_L.append(s)
                Ls_used.append(L)

            if not costs_per_L:
                n_missing += 1
            elif len(costs_per_L) < len(sizes):
                n_partial += 1
            else:
                n_full += 1

            if costs_per_L:
                cinf, sinf, slope, ok = _fit_continuum(
                    Ls_used, costs_per_L, sigmas_per_L, mode=cfg["fit"])
                COST_INF[i, j]  = cinf
                SIGMA_INF[i, j] = sinf
                SLOPE_INF[i, j] = slope
                OK_INF[i, j]    = ok
                fit_records.append({
                    "r1": float(r1), "r2": float(r2),
                    "Ls":     [int(L)   for L in Ls_used],
                    "costs":  [float(c) for c in costs_per_L],
                    "sigmas": [float(s) for s in sigmas_per_L],
                    "cost_inf":  float(cinf),
                    "sigma_inf": float(sinf),
                    "slope":     float(slope),
                    "ok":        bool(ok),
                })

            now = time.time()
            if now - last_status >= cfg["status_interval"]:
                pct = 100.0 * done_pts / total_pts
                log(f"--- compute status  "
                    f"{done_pts}/{total_pts} pts ({pct:.1f}%)  "
                    f"full={n_full}  partial={n_partial}  "
                    f"missing={n_missing} ---")
                last_status = now

    wall = time.time() - t0
    log(f"DONE compute: full={n_full}  partial={n_partial}  "
        f"missing={n_missing}  wall={wall:.1f}s")

    out_name = f"cost_{mode}"
    npz = os.path.join(root, f"{out_name}.npz")
    np.savez(
        npz,
        rs=rs, R1=R1, R2=R2,
        sizes=np.array(sizes, dtype=int),
        # geometry params so plot stage is self-contained
        test_Lx_mult=cfg["test_Lx_mult"], test_Ly_mult=cfg["test_Ly_mult"],
        test_Tx_frac=cfg["test_Tx_frac"], test_Ty_frac=cfg["test_Ty_frac"],
        ref_Lx_mult=cfg["ref_Lx_mult"],   ref_Ly_mult=cfg["ref_Ly_mult"],
        ref_Tx_frac=cfg["ref_Tx_frac"],   ref_Ty_frac=cfg["ref_Ty_frac"],
        cost_inf=COST_INF, sigma_inf=SIGMA_INF,
        slope_inf=SLOPE_INF, ok_inf=OK_INF,
        **{f"cost_L{L}":       COST_L[L]  for L in sizes},
        **{f"sigma_L{L}":      SIGMA_L[L] for L in sizes},
        **{f"betac_test_L{L}": BETAC_T[L] for L in sizes},
        **{f"betac_ref_L{L}":  BETAC_R[L] for L in sizes},
    )
    json_path = os.path.join(root, f"{out_name}_fits.json")
    with open(json_path, "w") as f:
        json.dump({
            "tag": cfg["tag"], "cost_mode": mode,
            "cost_power": cfg["cost_power"], "fit": cfg["fit"],
            "sizes": sizes,
            "test_geoms": manifest.get("test_geoms", {}),
            "ref_geoms":  manifest.get("ref_geoms",  {}),
            "points": fit_records,
        }, f, indent=2)
    log(f"wrote {npz}")
    log(f"wrote {json_path}")

    finite = COST_INF[np.isfinite(COST_INF)]
    if finite.size:
        idx = np.unravel_index(np.nanargmin(COST_INF), COST_INF.shape)
        log(f"continuum cost: min={float(np.nanmin(COST_INF)):.4e}  "
            f"max={float(np.nanmax(COST_INF)):.4e}  "
            f"argmin=(r1={rs[idx[0]]:.3f}, r2={rs[idx[1]]:.3f})")
    return {"npz": npz, "fits": json_path}


# ===========================================================================
# Stage 3 — plot
# ===========================================================================

def stage_plot(cfg, mode):
    banner(f"STAGE 3/3 — PLOT  tag={cfg['tag']}  mode={mode}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root     = os.path.join(_HERE, "results", cfg["tag"])
    npz_path = os.path.join(root, f"cost_{mode}.npz")
    if not os.path.exists(npz_path):
        log(f"ERROR: {npz_path} missing — run compute stage first.")
        sys.exit(1)
    z     = np.load(npz_path)
    rs    = z["rs"]
    sizes = list(z["sizes"])
    cost_inf  = z["cost_inf"]
    sigma_inf = z["sigma_inf"]

    def _gv(key, fallback):
        return float(z[key]) if key in z else fallback
    test_desc = (f"Lx={_gv('test_Lx_mult',1)}*L "
                 f"Ly={_gv('test_Ly_mult',1)}*L "
                 f"Tx={_gv('test_Tx_frac',0)}*L "
                 f"Ty={_gv('test_Ty_frac',0)}*L")
    ref_desc  = (f"Lx={_gv('ref_Lx_mult',1)}*L "
                 f"Ly={_gv('ref_Ly_mult',1)}*L "
                 f"Tx={_gv('ref_Tx_frac',0.25)}*L "
                 f"Ty={_gv('ref_Ty_frac',0.25)}*L")

    plot_dir = os.path.join(root, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    def _heatmap(ax, Z, title, log_scale=False, vmax=None, cmap="viridis"):
        Zp = Z.copy()
        if log_scale:
            Zp = np.log10(np.where(Zp > 0, Zp, np.nan))
        finite = Zp[np.isfinite(Zp)]
        if vmax is None and finite.size > 0:
            vmax = float(np.nanpercentile(finite, 99))
        vmin = float(np.nanmin(finite)) if finite.size > 0 else 0.0
        extent = [rs[0], rs[-1], rs[0], rs[-1]]
        im = ax.imshow(Zp.T, origin="lower", extent=extent, aspect="auto",
                       cmap=cmap, vmin=vmin, vmax=vmax)
        if finite.size > 0:
            idx = np.unravel_index(np.nanargmin(Z), Z.shape)
            ax.plot(rs[idx[0]], rs[idx[1]], "w*", ms=10,
                    markeredgecolor="k",
                    label=f"min=({rs[idx[0]]:.2f},{rs[idx[1]]:.2f})")
            ax.legend(loc="upper right", fontsize=7)
        ax.set_xlabel("r1")
        ax.set_ylabel("r2")
        ax.set_title(title, fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    n_panels = len(sizes) + 2
    n_cols   = min(3, n_panels)
    n_rows   = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.6 * n_cols, 4.2 * n_rows),
                             squeeze=False)
    axes_flat = axes.flatten()
    for k, L in enumerate(sizes):
        _heatmap(axes_flat[k], z[f"cost_L{L}"],
                 f"L={L}  cost={mode}",
                 log_scale=cfg["plot_log"], vmax=cfg["plot_vmax"])
    _heatmap(axes_flat[len(sizes)], cost_inf,
             f"continuum (L->inf)  mode={mode}",
             log_scale=cfg["plot_log"], vmax=cfg["plot_vmax"])
    _heatmap(axes_flat[len(sizes) + 1], sigma_inf,
             "sigma (continuum extrapolation)", cmap="magma")
    for ax in axes_flat[n_panels:]:
        ax.axis("off")
    fig.suptitle(
        f"tag={cfg['tag']}  mode={mode}  sizes={sizes}\n"
        f"test: {test_desc}\nref:  {ref_desc}",
        fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    suffix = "_log" if cfg["plot_log"] else ""
    panel_out = os.path.join(plot_dir, f"landscape_{mode}{suffix}.png")
    fig.savefig(panel_out, dpi=140)
    plt.close(fig)
    log(f"wrote {panel_out}")

    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    _heatmap(ax, cost_inf,
             f"Continuum cost  mode={mode}\ntest: {test_desc}\nref: {ref_desc}",
             log_scale=cfg["plot_log"], vmax=cfg["plot_vmax"])
    fig.tight_layout()
    cont_out = os.path.join(plot_dir, f"continuum_{mode}{suffix}.png")
    fig.savefig(cont_out, dpi=160)
    plt.close(fig)
    log(f"wrote {cont_out}")
    return panel_out


# ===========================================================================
# CLI
# ===========================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # bookkeeping
    p.add_argument("--tag",    default=CONFIG["tag"])
    p.add_argument("--exe",    default=CONFIG["exe"],
                   help="path to simulator (default: bin/ in this folder)")
    p.add_argument("--stages", nargs="+", default=CONFIG["stages"],
                   choices=["precompute", "compute", "plot"])

    # coupling grid
    g = p.add_argument_group("coupling grid")
    g.add_argument("--r-min",  type=float, default=CONFIG["r_min"])
    g.add_argument("--r-max",  type=float, default=CONFIG["r_max"])
    g.add_argument("--r-step", type=float, default=CONFIG["r_step"])

    # FSS sizes
    g = p.add_argument_group("FSS lattice sizes")
    g.add_argument("--sizes", type=int, nargs="+", default=CONFIG["sizes"],
                   metavar="L",
                   help="scale levels for continuum extrapolation (>=2 values)")

    # test-lattice geometry
    g = p.add_argument_group("test-lattice geometry  [Lx=round(mult*L), Tx=round(frac*L)]")
    g.add_argument("--test-Lx-mult", type=float, default=CONFIG["test_Lx_mult"],
                   metavar="M", help="Lx = round(M*L)  [default 1.0]")
    g.add_argument("--test-Ly-mult", type=float, default=CONFIG["test_Ly_mult"],
                   metavar="M", help="Ly = round(M*L)  [default 1.0]")
    g.add_argument("--test-Tx-frac", type=float, default=CONFIG["test_Tx_frac"],
                   metavar="F", help="Tx = round(F*L)  [default 0.0 = untwisted]")
    g.add_argument("--test-Ty-frac", type=float, default=CONFIG["test_Ty_frac"],
                   metavar="F", help="Ty = round(F*L)  [default 0.0]")

    # reference-lattice geometry
    g = p.add_argument_group("reference-lattice geometry")
    g.add_argument("--ref-Lx-mult",  type=float, default=CONFIG["ref_Lx_mult"],
                   metavar="M", help="Lx = round(M*L)  [default 1.0]")
    g.add_argument("--ref-Ly-mult",  type=float, default=CONFIG["ref_Ly_mult"],
                   metavar="M", help="Ly = round(M*L)  [default 1.0]")
    g.add_argument("--ref-Tx-frac",  type=float, default=CONFIG["ref_Tx_frac"],
                   metavar="F", help="Tx = round(F*L)  [default 0.25]")
    g.add_argument("--ref-Ty-frac",  type=float, default=CONFIG["ref_Ty_frac"],
                   metavar="F", help="Ty = round(F*L)  [default 0.25]")

    # Monte Carlo
    g = p.add_argument_group("Monte Carlo")
    g.add_argument("--n-traj",        type=int,   default=CONFIG["n_traj"],
                   help="production trajectories per grid point")
    g.add_argument("--n-traj-coarse", type=int,   default=CONFIG["n_traj_coarse"],
                   help="beta_c finder coarse-pass trajectories")
    g.add_argument("--n-traj-fine",   type=int,   default=CONFIG["n_traj_fine"],
                   help="beta_c finder fine-pass trajectories")
    g.add_argument("--n-skip",        type=int,   default=CONFIG["n_skip"],
                   help="measurement interval in sweeps")
    g.add_argument("--n-therm",       type=int,   default=CONFIG["n_therm"],
                   help="thermalisation sweeps")
    g.add_argument("--beta-seed", type=float, nargs=2,
                   default=list(CONFIG["beta_seed"]),
                   metavar=("LO", "HI"),
                   help="initial beta bracket for beta_c finder")
    g.add_argument("--n-workers", type=int, default=CONFIG["n_workers"],
                   help="parallel MC processes")

    # cost replay
    g = p.add_argument_group("cost replay")
    g.add_argument("--cost-modes", nargs="+", default=CONFIG["cost_modes"],
                   choices=["l4mean_both_interp", "test_native",
                            "spectral", "affine_fit", "huber_log"],
                   help="one or more cost modes (each gets its own .npz + plots)")
    g.add_argument("--cost-power", type=int, default=CONFIG["cost_power"],
                   help="residual exponent for test_native: 2 or 4")
    g.add_argument("--fit", default=CONFIG["fit"],
                   choices=["linear", "quadratic"],
                   help="continuum extrapolation model")

    # plotting
    g = p.add_argument_group("plotting")
    g.add_argument("--plot-log",  action="store_true", default=CONFIG["plot_log"],
                   help="render log10(cost)")
    g.add_argument("--plot-vmax", type=float, default=CONFIG["plot_vmax"],
                   help="cap colour scale (default: 99th percentile)")

    # monitoring
    g = p.add_argument_group("monitoring")
    g.add_argument("--status-interval", type=float,
                   default=CONFIG["status_interval"],
                   help="seconds between heartbeat status lines (default 60)")
    return p


def cfg_from_args(args):
    cfg = dict(CONFIG)
    cfg.update({
        "tag":    args.tag,
        "exe":    args.exe,
        "stages": args.stages,
        "r_min":  args.r_min,
        "r_max":  args.r_max,
        "r_step": args.r_step,
        "sizes":  args.sizes,
        # test geometry
        "test_Lx_mult": args.test_Lx_mult,
        "test_Ly_mult": args.test_Ly_mult,
        "test_Tx_frac": args.test_Tx_frac,
        "test_Ty_frac": args.test_Ty_frac,
        # ref geometry
        "ref_Lx_mult": args.ref_Lx_mult,
        "ref_Ly_mult": args.ref_Ly_mult,
        "ref_Tx_frac": args.ref_Tx_frac,
        "ref_Ty_frac": args.ref_Ty_frac,
        # MC
        "n_traj":        args.n_traj,
        "n_traj_coarse": args.n_traj_coarse,
        "n_traj_fine":   args.n_traj_fine,
        "n_skip":        args.n_skip,
        "n_therm":       args.n_therm,
        "beta_seed":     tuple(args.beta_seed),
        "n_workers":     args.n_workers,
        # cost
        "cost_modes":  args.cost_modes,
        "cost_power":  args.cost_power,
        "fit":         args.fit,
        # plot
        "plot_log":  args.plot_log,
        "plot_vmax": args.plot_vmax,
        # monitoring
        "status_interval": args.status_interval,
    })
    return cfg


# ===========================================================================
# Entry point
# ===========================================================================

def main(argv=None):
    args = build_parser().parse_args(argv)
    cfg  = cfg_from_args(args)

    banner(f"K_from_BC_Heatmap run.py  tag={cfg['tag']}  "
           f"stages={cfg['stages']}")
    log(f"config:\n{json.dumps({k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()}, indent=2)}")

    run_dir = os.path.join(_HERE, "results", cfg["tag"])
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "run_config.json"), "w") as f:
        json.dump({k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in cfg.items()}, f, indent=2)

    if "precompute" in cfg["stages"]:
        stage_precompute(cfg)
    if "compute" in cfg["stages"]:
        for mode in cfg["cost_modes"]:
            stage_compute(cfg, mode)
    if "plot" in cfg["stages"]:
        for mode in cfg["cost_modes"]:
            stage_plot(cfg, mode)

    banner(f"ALL DONE  total wall={_elapsed()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
