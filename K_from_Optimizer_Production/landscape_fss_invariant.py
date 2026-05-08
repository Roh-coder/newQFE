"""
landscape_fss_invariant.py — FSS-invariant scale-aware cost variants.

Variants (all replayed on the existing precompute grid):

A) huber_log baseline                         (broken — bowl to corner)
B) huber_log + λ·(β_c·L_test − β_c_ref·L_ref)²  scale-aware via FSS
C) huber_log against ANALYTIC CFT target       (use cft_torus.ising_torus_F
   at the test's bare (16,16,0,0) modular τ — no MC reference needed)
D) (B) + (C)                                    (combined scale-aware)

The key physical idea: β_c·L is the standard FSS-invariant combination
for ν=1 (2d Ising). For two geometries to represent the same continuum
CFT, their β_c·L_eff values must agree. L_eff = √area of the torus.
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
import landscape_dozen as ld
from cft_torus import (torus_periods, ising_torus_F, modular_tau,
                       lattice_sites, _SQRT3_2)
from cost import (boundary_paths, _direction_lattice_steps)

_REF_GEOMS = {"_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3)}
_TRUTH_456 = (5.0652, 7.7429)
_BETA_C_REF = 0.25605    # measured ref β_c (from reference_meta.json)


def _torus_area(w1, w2):
    return abs(w1.real * w2.imag - w1.imag * w2.real)


def _analytic_cft_samples(test_Lx, test_Ly, test_Tx, test_Ty,
                          ref_pack):
    """For each direction, evaluate the analytic Ising CFT correlator
    at the integer lattice sites visited by the test along its boundary
    path, using the test's BARE modular τ.
    """
    w1, w2 = torus_periods(test_Lx, test_Ly, test_Tx, test_Ty)
    tau = modular_tau(test_Lx, test_Ly, test_Tx, test_Ty)
    test_paths = boundary_paths(test_Lx, test_Ly, test_Tx, test_Ty)

    out = []
    for d, (tdm, tdn) in zip(ref_pack, test_paths):
        if d is None: out.append(None); continue
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        Gs = np.empty(N)
        for i, (m, n) in enumerate(zip(ms, ns)):
            x = m + 0.5 * n
            y = _SQRT3_2 * n
            zeta = complex(x, y) / w1
            Gs[i] = ising_torus_F(zeta, tau)
        # Subtract value at the largest |zeta| as a "G_far" baseline so
        # connected-correlator conventions match the measured data.
        out.append(Gs)
    return out


def k_huber_log(dirs, delta=0.5, eps=1e-12):
    pd = []
    for d in dirs:
        if d is None: pd.append(None); continue
        m = (np.isfinite(d["G_test"]) & np.isfinite(d["G_ref"])
             & (d["G_test"] > eps) & (d["G_ref"] > eps))
        if m.sum() < 2: pd.append(None); continue
        r = np.log(d["G_test"][m]) - np.log(d["G_ref"][m])
        a = np.abs(r)
        rho = np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))
        pd.append(float(np.mean(rho)))
    arr = [c for c in pd if c is not None and np.isfinite(c)]
    return float(np.mean(arr)) if arr else np.nan


def k_huber_log_vs_analytic(dirs, analytic, delta=0.5, eps=1e-12):
    """huber_log but against the analytic CFT correlator (per-dir)."""
    pd = []
    # Normalisation: scale measured Gt to match analytic at first sample.
    for d, A in zip(dirs, analytic):
        if d is None or A is None: pd.append(None); continue
        Gt = d["G_test"]
        m = np.isfinite(Gt) & np.isfinite(A) & (Gt > eps) & (A > eps)
        if m.sum() < 2: pd.append(None); continue
        Gt_, A_ = Gt[m], A[m]
        # Match overall scale (the analytic F has arbitrary normalisation).
        scale = float(np.median(Gt_ / A_))
        A_scaled = A_ * scale
        r = np.log(Gt_) - np.log(A_scaled)
        a = np.abs(r)
        rho = np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))
        pd.append(float(np.mean(rho)))
    arr = [c for c in pd if c is not None and np.isfinite(c)]
    return float(np.mean(arr)) if arr else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag", default="_reference_Lx13_Ly16_Tx-3_Ty3")
    ap.add_argument("--lambda-fss", type=float, default=200.0,
                    help="weight on (β_c·L_test - β_c_ref·L_ref)² penalty")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    manifest = json.load(open(os.path.join(root, "manifest.json")))
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    rL, rL2, rT, rT2 = _REF_GEOMS[args.ref_tag]
    ref_data = mc_engine.load_all_to_all(
        os.path.join(_HERE, "results", args.ref_tag,
                     "two_point_all_to_all.dat"))
    ref_pack = ld._ref_samples(ref_data, rL, rL2, rT, rT2, Lx, Ly, Tx, Ty)

    # FSS scale handle.
    w1_test, w2_test = torus_periods(Lx, Ly, Tx, Ty)
    L_test = float(np.sqrt(_torus_area(w1_test, w2_test)))
    w1_ref,  w2_ref  = torus_periods(rL, rL2, rT, rT2)
    L_ref  = float(np.sqrt(_torus_area(w1_ref, w2_ref)))
    fss_target = _BETA_C_REF * L_ref
    print(f"[fss] L_test={L_test:.3f}  L_ref={L_ref:.3f}  "
          f"β_c_ref·L_ref={fss_target:.4f}  truth β_c·L_test≈{0.0628*L_test:.4f}")

    # Analytic CFT samples on test bare τ (constant per grid point — depends
    # only on test geometry, not on r1,r2).
    analytic = _analytic_cft_samples(Lx, Ly, Tx, Ty, ref_pack)
    print(f"[analytic] τ_test = {modular_tau(Lx,Ly,Tx,Ty)}")

    print("\n[load + eval]")
    t0 = time.time()
    Z_huber = np.full((n, n), np.nan)
    Z_anal  = np.full((n, n), np.nan)
    BC = np.full((n, n), np.nan)
    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            pt = pickle.load(open(pkl, "rb"))
            dirs = ld._test_samples(pt["test_data"], ref_pack, Lx, Ly, Tx, Ty)
            Z_huber[i, j] = k_huber_log(dirs)
            Z_anal[i, j]  = k_huber_log_vs_analytic(dirs, analytic)
            BC[i, j]      = pt["beta_c"]
    print(f"  {time.time()-t0:.1f}s for {n*n} pts")

    # Build composite costs.
    Z_fss_pen = (BC * L_test - fss_target) ** 2
    Z_combo_fss = Z_huber + args.lambda_fss * Z_fss_pen
    Z_combo_all = Z_anal  + args.lambda_fss * Z_fss_pen

    panels = [
        ("A) huber_log (baseline)",                           Z_huber),
        (f"B) huber_log + {args.lambda_fss:.0f}·(β_c·L_test − β_c_ref·L_ref)²", Z_combo_fss),
        ("C) huber_log vs analytic CFT (no MC ref)",          Z_anal),
        (f"D) (C) + {args.lambda_fss:.0f}·FSS penalty",       Z_combo_all),
    ]

    print("\n[summary]")
    for title, Z in panels:
        if np.all(np.isnan(Z)):
            print(f"  {title:60s} all-NaN"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        d_truth = float(np.hypot(rs[ij[0]] - _TRUTH_456[0],
                                 rs[ij[1]] - _TRUTH_456[1]))
        ti = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        tj = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        print(f"  {title:60s} argmin=({rs[ij[0]]:5.2f},{rs[ij[1]]:5.2f})  "
              f"d_truth={d_truth:.2f}  cost(min)={Z[ij]:.3e}  "
              f"cost(truth)={Z[ti,tj]:.3e}")

    # Render
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (title, Z) in zip(axes, panels):
        Zp = np.log10(np.maximum(Z, 1e-15))
        vmin, vmax = np.nanpercentile(Zp, [2, 98])
        im = ax.pcolormesh(R1, R2, Zp, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log10 cost")
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=14, mec="k", mew=1.2,
                label="truth")
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        ax.plot(rs[ij[0]], rs[ij[1]], "w*", ms=14, mec="k", mew=1.2,
                label="argmin")
        ax.set_xlabel("r1"); ax.set_ylabel("r2")
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"FSS-invariant cost diagnostics  test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)  "
                 f"L_test={L_test:.2f}  L_ref={L_ref:.2f}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(root, "fss_invariant.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\n[done] figure → {out}")
    np.savez(os.path.join(root, "fss_invariant.npz"),
             R1=R1, R2=R2, BC=BC,
             Z_huber=Z_huber, Z_anal=Z_anal,
             Z_combo_fss=Z_combo_fss, Z_combo_all=Z_combo_all,
             L_test=L_test, L_ref=L_ref, fss_target=fss_target,
             lambda_fss=args.lambda_fss)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
