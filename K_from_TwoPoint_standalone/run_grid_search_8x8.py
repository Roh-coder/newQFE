#!/usr/bin/env python3
"""
run_grid_search_8x8.py — Fast 8×8 equilateral convergence test.
Uses the hybrid GC3 + Z²/Z⁴ strategy (same engine as run_grid_search.py).
Runs in parallel with the 32×32 test; uses a separate output directory.
"""

# ===========================================================================
# USER CONFIGURATION
# ===========================================================================
Lx, Ly = 8, 8
Tx, Ty = 0, 0

REF_Lx, REF_Ly = 8, 8
REF_Tx, REF_Ty = 0, 0

R1_INIT    = 1.5
R2_INIT    = 1.5
HALF_SPAN  = 0.30

BETA_INIT = 0.27

N_TRAJ_PROD         = 50000
N_TRAJ_SCAN_COARSE  = 50000
N_TRAJ_SCAN_FINE    = 100000

MAX_ITER = 6
N_GRID   = 5
TOP_N    = 3           # re-score top N with Z⁴ per level

OUTPUT_DIR = "results/equilateral_8x8_from_1.5"

REFERENCE  = None
REF_N_TRAJ = 100000   # 8×8 is cheap; 100k is fast and sufficient

# Exact beta_c for the reference.  Set to skip susceptibility scan.
# ln(3)/4 ≈ 0.27465307216702745 for equilateral triangular Ising.
REF_BETA_C = None

import sys as _sys
EXE = ("bin/ising_tri_twisted_parallelogram.exe"
       if _sys.platform == "win32"
       else "bin/ising_tri_twisted_parallelogram")
# ===========================================================================
# END OF USER CONFIGURATION
# ===========================================================================

import os, sys, json, time, warnings
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
os.chdir(_HERE)

import optimise_couplings as oc
from scipy.optimize import curve_fit, minimize_scalar

if not os.path.isfile(EXE):
    print(f"Simulator not found at '{EXE}' — build it first (make).")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Reference generation / loading
# ---------------------------------------------------------------------------
if REFERENCE is None:
    _ref_tag = f"ref_{REF_Lx}x{REF_Ly}" + (f"_t{REF_Tx}x{REF_Ty}" if REF_Tx or REF_Ty else "")
    REFERENCE = os.path.join("ref_data", _ref_tag, "ref_metadata.json")

if not os.path.isfile(REFERENCE):
    print(f"Reference not found at '{REFERENCE}' — generating it now.")
    print(f"  Lattice: {REF_Lx}x{REF_Ly}  Tx={REF_Tx} Ty={REF_Ty}")
    print(f"  Trajectories: {REF_N_TRAJ}")
    print()
    ref_output_dir = os.path.dirname(REFERENCE) or "."
    sys.argv = [
        "optimise_couplings.py",
        "--exe", EXE,
        "--Lx", str(REF_Lx), "--Ly", str(REF_Ly),
        "--Tx", str(REF_Tx), "--Ty", str(REF_Ty),
        "--beta_init", str(BETA_INIT),
        "--ref_n_traj", str(REF_N_TRAJ),
        "--output_dir", ref_output_dir,
        "--gen_ref",
    ]
    if REF_BETA_C is not None:
        sys.argv += ["--ref_beta_c", str(REF_BETA_C)]
    import importlib
    importlib.reload(oc)
    oc.main()
    if not os.path.isfile(REFERENCE):
        print(f"ERROR: Reference generation did not produce '{REFERENCE}'")
        sys.exit(1)
    print(f"Reference generated.\n")

with open(REFERENCE) as f:
    ref_meta = json.load(f)
ref_data = oc.load_all_to_all(ref_meta["a2a_path"])
print(f"Reference loaded: {REFERENCE}")
print(f"  β_c = {ref_meta['beta_c']:.10f}, Lx={ref_meta['Lx']}, Ly={ref_meta['Ly']}")

L_eff = min(Lx, Ly)
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# GC3 β_c finder (3-pass Gram-Charlier)
# ---------------------------------------------------------------------------
def _gc_refit(sb, sc, beta_init):
    """Fit Gram-Charlier to scan data and return mode."""
    def _gram_charlier(b, mu, sigma, skew, kurt, A):
        z = (b - mu) / sigma
        H3 = z**3 - 3.0 * z
        H4 = z**4 - 6.0 * z**2 + 3.0
        return A * np.exp(-0.5 * z**2) * (1.0 + (skew / 6.0) * H3 + (kurt / 24.0) * H4)

    sb_arr, sc_arr = np.array(sb), np.array(sc)
    b_rng = sb_arr.max() - sb_arr.min()
    p0 = [beta_init, b_rng / 2, 0.0, 0.0, float(sc_arr.max())]
    bounds = ([sb_arr.min(), 1e-8, -5.0, -3.0, 0.0],
              [sb_arr.max(), b_rng * 2, 5.0, 10.0, np.inf])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        popt, _ = curve_fit(_gram_charlier, sb_arr, sc_arr,
                            p0=p0, bounds=bounds, maxfev=8000)
    mu_fit, sigma_fit = popt[0], popt[1]
    mode_lo = max(sb_arr.min(), mu_fit - 3.0 * sigma_fit)
    mode_hi = min(sb_arr.max(), mu_fit + 3.0 * sigma_fit)
    res = minimize_scalar(lambda b: -_gram_charlier(b, *popt),
                          bounds=(mode_lo, mode_hi), method='bounded')
    return float(np.clip(res.x, sb_arr.min(), sb_arr.max())), tuple(popt)


def find_beta_gc3(Lx, Ly, Tx, Ty, k1, k2, k3, beta_guess, scan_dir):
    """3-pass Gram-Charlier β_c finder."""
    margin = max(0.20 * beta_guess, 0.04)
    beta_lo = max(0.01, beta_guess - margin)
    beta_hi = beta_guess + margin

    beta_c_2, chi_peak, sb, sc, se, fit_params_12 = oc.find_beta_c(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3,
        beta_lo, beta_hi,
        n_coarse=11, n_refine=5, n_refine2=5,
        n_traj_coarse=N_TRAJ_SCAN_COARSE, n_traj_fine=N_TRAJ_SCAN_FINE,
        data_dir=scan_dir,
    )
    b_range = max(sb) - min(sb)
    half3 = b_range / (2 * 11) / 4
    fine3_betas = np.linspace(beta_c_2 - half3, beta_c_2 + half3, 7)
    for b in fine3_betas:
        scan_dir3 = scan_dir + "_p3"
        stdout_f, _ = oc.run_simulator(
            EXE, Lx, Ly, Tx, Ty, k1, k2, k3, float(b),
            n_traj=N_TRAJ_SCAN_FINE, data_dir=scan_dir3)
        chi_f, chi_f_err = oc.parse_stdout_value_with_err(stdout_f, "m_susc:")
        sb.append(float(b))
        sc.append(chi_f)
        se.append(chi_f_err)

    fit_params_3 = None
    try:
        beta_c_3, fit_params_3 = _gc_refit(sb, sc, beta_c_2)
    except Exception:
        beta_c_3 = beta_c_2
    # Use pass-2 fit from find_beta_c as "coarse" and pass-3 as "fine"
    fp_coarse = fit_params_12[1] if fit_params_12 else None
    fp_fine = fit_params_3
    return beta_c_3, chi_peak, sb, sc, se, (fp_coarse, fp_fine)


# ---------------------------------------------------------------------------
# Evaluate one grid point: GC3 → production → Z² cost
# ---------------------------------------------------------------------------
def evaluate_point_gc3(r1, r2, beta_guess, level, i, j, data_dir):
    """Return (beta_c, z2_cost, z2_ndof, test_data_path)."""
    k3 = 1.0
    k1, k2 = r1 * k3, r2 * k3
    label = f"L{level}_i{i}_j{j}"

    level_dir = os.path.join(data_dir, f"level{level:02d}")
    os.makedirs(level_dir, exist_ok=True)
    scan_dir = os.path.join(level_dir, f"scan_{label}")

    plotter.set_current_point(r1, r2)
    print(f"  [{label}] r1={r1:.6f} r2={r2:.6f}: GC3 scan ...", end="", flush=True)
    t0 = time.time()
    beta_c, chi_peak, scan_b, scan_c, scan_e, gc_fit_params = find_beta_gc3(
        Lx, Ly, Tx, Ty, k1, k2, k3, beta_guess, scan_dir)
    t_scan = time.time() - t0
    print(f" β_c={beta_c:.8f} ({t_scan:.0f}s)", end="", flush=True)

    prod_dir = os.path.join(level_dir, f"prod_{label}")
    stdout, subdir = oc.run_simulator(
        EXE, Lx, Ly, Tx, Ty, k1, k2, k3, beta_c,
        n_traj=N_TRAJ_PROD, n_therm=3000, data_dir=prod_dir)
    a2a_path = os.path.join(subdir, "two_point_all_to_all.dat")
    test_data = oc.load_all_to_all(a2a_path)

    chi2, ndof, chi2_err = oc.compute_chi2(
        ref_data, test_data,
        r_min=0.0, r_max_frac=0.33, L_eff=L_eff, Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
        ref_Lx=REF_Lx, ref_Ly=REF_Ly, ref_Tx=REF_Tx, ref_Ty=REF_Ty,
        cost="boundary_slices_normed",
    )
    z2_cost = chi2 / max(ndof, 1)
    t_total = time.time() - t0
    print(f"  Z²={z2_cost:.4f} ({t_total:.0f}s)")
    sys.stdout.flush()

    slices = None
    try:
        slices = oc.extract_boundary_slices(
            ref_data, test_data, Lx, Ly, Tx, Ty,
            ref_Lx=REF_Lx, ref_Ly=REF_Ly,
            ref_Tx=REF_Tx, ref_Ty=REF_Ty)
    except Exception:
        pass

    plotter.record_point(r1, r2, beta_c, chi2, ndof,
                         scan_b, scan_c, scan_e,
                         parabola=gc_fit_params,
                         chi2_err=chi2_err, slices=slices)

    return beta_c, z2_cost, ndof, a2a_path


def rescore_z4(a2a_path):
    """Re-score an already-run production with Z⁴ (no extra MC needed)."""
    test_data = oc.load_all_to_all(a2a_path)
    chi2, ndof, chi2_err = oc.compute_chi2(
        ref_data, test_data,
        r_min=0.0, r_max_frac=0.33, L_eff=L_eff, Lx=Lx, Ly=Ly, Tx=Tx, Ty=Ty,
        ref_Lx=REF_Lx, ref_Ly=REF_Ly, ref_Tx=REF_Tx, ref_Ty=REF_Ty,
        cost="boundary_slices_normed_quartic",
    )
    return chi2 / max(ndof, 1)


# ---------------------------------------------------------------------------
# Hybrid grid search: GC3 + Z² grid / Z⁴ re-score
# ---------------------------------------------------------------------------
print(f"\n{'='*70}")
print(f"Hybrid Grid Search: GC3 + Z² grid / Z⁴ re-score")
print(f"  Lattice: {Lx}×{Ly}  Grid: {N_GRID}×{N_GRID}  Levels: {MAX_ITER}")
print(f"  Start: r1={R1_INIT}, r2={R2_INIT}, half_span={HALF_SPAN}")
print(f"  MC: prod={N_TRAJ_PROD}, coarse={N_TRAJ_SCAN_COARSE}, fine={N_TRAJ_SCAN_FINE}")
print(f"  Top-{TOP_N} Z⁴ re-scoring per level")
print(f"{'='*70}\n")

r1_center = R1_INIT
r2_center = R2_INIT
hs_r1 = HALF_SPAN
hs_r2 = HALF_SPAN
best_beta_c = BETA_INIT

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plotter = oc.LivePlotter(OUTPUT_DIR, cost="boundary_slices_normed")

all_results = []
run_dir = os.path.join(OUTPUT_DIR, "runs")

for level in range(MAX_ITER):
    r1_lo = max(r1_center - hs_r1, 0.1)
    r1_hi = r1_center + hs_r1
    r2_lo = max(r2_center - hs_r2, 0.1)
    r2_hi = r2_center + hs_r2

    r1_vals = np.linspace(r1_lo, r1_hi, N_GRID)
    r2_vals = np.linspace(r2_lo, r2_hi, N_GRID)

    print(f"\n{'='*70}")
    print(f"Level {level}: centre=({r1_center:.4f}, {r2_center:.4f})"
          f"  hs=({hs_r1:.4f}, {hs_r2:.4f})")
    print(f"  r1 ∈ [{r1_lo:.4f}, {r1_hi:.4f}]  r2 ∈ [{r2_lo:.4f}, {r2_hi:.4f}]")
    print(f"{'='*70}")
    sys.stdout.flush()

    plotter.start_level(level, r1_vals, r2_vals)

    grid = {}
    for i, r1 in enumerate(r1_vals):
        for j, r2 in enumerate(r2_vals):
            beta_c, z2, ndof, a2a_path = evaluate_point_gc3(
                r1, r2, best_beta_c, level, i, j, run_dir)
            grid[(i, j)] = {
                "r1": r1, "r2": r2, "beta_c": beta_c,
                "z2": z2, "a2a_path": a2a_path,
            }
            all_results.append({
                "level": level, "i": i, "j": j,
                "r1": r1, "r2": r2, "beta_c": beta_c,
                "z2_cost": z2,
            })

    ranked = sorted(grid.keys(), key=lambda k: grid[k]["z2"])
    top_keys = ranked[:TOP_N]

    print(f"\n  Re-scoring top {TOP_N} with Z⁴:")
    for rank, key in enumerate(top_keys):
        pt = grid[key]
        z4 = rescore_z4(pt["a2a_path"])
        pt["z4"] = z4
        print(f"    #{rank+1}: ({key[0]},{key[1]}) r1={pt['r1']:.4f} r2={pt['r2']:.4f}"
              f"  Z²={pt['z2']:.4f}  Z⁴={z4:.4f}")

    best_key = min(top_keys, key=lambda k: grid[k].get("z4", grid[k]["z2"]))
    bi, bj = best_key
    best = grid[best_key]
    r1_best, r2_best = best["r1"], best["r2"]
    best_beta_c = best["beta_c"]
    best_z2 = best["z2"]
    best_z4 = best.get("z4", best_z2)

    print(f"\n  Level {level} best: ({bi},{bj})  r1={r1_best:.6f} r2={r2_best:.6f}"
          f"  Z²={best_z2:.4f}  Z⁴={best_z4:.4f}  β_c={best_beta_c:.8f}")

    r1_border = (bi == 0) or (bi == N_GRID - 1)
    r2_border = (bj == 0) or (bj == N_GRID - 1)

    level_summary = {
        "level": level,
        "r1_center": r1_center, "r2_center": r2_center,
        "hs_r1": hs_r1, "hs_r2": hs_r2,
        "r1_vals": r1_vals.tolist(), "r2_vals": r2_vals.tolist(),
        "best_ij": [bi, bj],
        "best_r1": r1_best, "best_r2": r2_best,
        "best_z2": best_z2, "best_z4": best_z4,
        "best_beta_c": best_beta_c,
        "r1_border": r1_border, "r2_border": r2_border,
        "top_n_z4": [
            {"ij": list(k), "r1": grid[k]["r1"], "r2": grid[k]["r2"],
             "z2": grid[k]["z2"], "z4": grid[k].get("z4")}
            for k in top_keys
        ],
    }
    with open(os.path.join(OUTPUT_DIR, f"grid_level_{level:02d}.json"), "w") as f:
        json.dump(level_summary, f, indent=2)

    # Centre on best; halve span only for interior axes.
    r1_center = r1_best
    if r1_border:
        print(f"  r1: border → translate to {r1_center:.4f} (span unchanged)")
    else:
        hs_r1 /= 2.0
        print(f"  r1: interior → refine, hs_r1={hs_r1:.4f}")

    r2_center = r2_best
    if r2_border:
        print(f"  r2: border → translate to {r2_center:.4f} (span unchanged)")
    else:
        hs_r2 /= 2.0
        print(f"  r2: interior → refine, hs_r2={hs_r2:.4f}")

    sys.stdout.flush()


# ---------------------------------------------------------------------------
# Final results
# ---------------------------------------------------------------------------
print(f"\n{'='*70}")
print("OPTIMISATION COMPLETE")
print(f"{'='*70}")

best_overall = min(all_results, key=lambda r: r["z2_cost"])
print(f"Best point: r1={best_overall['r1']:.6f} r2={best_overall['r2']:.6f}"
      f"  Z²={best_overall['z2_cost']:.6f}  β_c={best_overall['beta_c']:.8f}")

fit_result = {
    "best_r1": best_overall["r1"],
    "best_r2": best_overall["r2"],
    "best_k3": 1.0,
    "best_beta_c": best_overall["beta_c"],
    "best_z2_cost": best_overall["z2_cost"],
    "method": "gc3 + Z²/Z⁴ hybrid",
    "levels_completed": min(level + 1, MAX_ITER),
    "all_results": all_results,
}

with open(os.path.join(OUTPUT_DIR, "fit_result.json"), "w") as f:
    json.dump(fit_result, f, indent=2)

plotter.save_final()
print(f"\nResults saved: {OUTPUT_DIR}/fit_result.json")
print("Done.")
