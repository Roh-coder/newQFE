#!/usr/bin/env python3
"""
cost_zoo_test.py
================

Test ALL of the proposed cost-function ideas (10 brainstormed + 2
baselines) against the same controlled MC-stability dataset used by
spectral_vs_native_test.py.

Each candidate is implemented as a pair-cost  C(ref_data, test_data) -> float
acting on the all-to-all dict produced by mc_engine.load_all_to_all.

Discrimination metric (same as spectral_vs_native_test.py):
    Z = (mean_mismatched - mean_matched) / std_matched
Higher Z = the optimizer can see the (r1, r2) signal more clearly above
MC noise.

Candidates
----------
A. Decorrelate noise:
   1. chi2_native       — test-native L² weighted by 1/(σ_t² + σ_r²)
   2. baseline_subtract — test_native cost minus per-replicate self-cost
                          (estimated from ref-vs-ref and test-vs-test)
   3. jackknife         — N/A (requires per-config blocks not stored
                          in two_point_all_to_all.dat)
B. Orthogonal scalars:
   4. cycle_loop        — Σ_d (Φ_d^t - Φ_d^r)² where Φ_d = Σ G along
                          period vector d
   5. moment_r2         — (⟨r²⟩_t - ⟨r²⟩_r)²  with weight |G|
   6. eff_xi            — Σ_d (ξ_d^t - ξ_d^r)² from log-slope of |G|
                          along each period vector
C. Reference-free isotropy / symmetry checks:
   7. angular_iso       — Σ_shells Σ_{m≥2} |Ĝ_m(r)|² evaluated on
                          {test} minus {ref}, then squared
   8. cycle_var         — variance of {Φ_v, Φ_u, Φ_w}  (single-sample
                          stat; pair cost = (var_t - var_r)²)
D. Use ignored structure:
   9. affine_fit        — best a,b in min‖a·G_t + b - G_r‖²; cost = the
                          residual after the affine collapse
  10. offdiag_cov       — N/A (needs per-config covariance)

Plus baselines:
   B1. native           — test_native L² (current production cost)
   B2. spectral         — NDFT polar grid (cost.py 'spectral' mode)

Outputs
-------
  results/_cost_compare/
    zoo_summary.json
    zoo_compare.png      — Z-scores for all candidates side-by-side
    zoo_distributions.png — per-candidate matched/r2/r5 distributions
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cost as cost_mod
import mc_engine
from cost import (
    _SQRT3_2, boundary_paths, _direction_lattice_steps,
    _lookup_test_value, _tile_interp,
)

ROOT     = Path(__file__).resolve().parent
STAB_DIR = ROOT / "results" / "_betac_stability_test"
OUT_DIR  = ROOT / "results" / "_cost_compare"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GEOM = dict(Lx=12, Ly=12, Tx=0, Ty=0)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_group(group_name: str):
    base = STAB_DIR / group_name
    if not base.is_dir():
        raise SystemExit(f"missing dataset: {base}")
    out = []
    for sub in sorted(base.iterdir()):
        cand = list(sub.glob("*/two_point_all_to_all.dat"))
        if not cand:
            continue
        out.append(mc_engine.load_all_to_all(str(cand[0])))
    print(f"[load] {group_name}: {len(out)} replicates")
    return out


def _arrays(data):
    """Vectorise the all-to-all dict to numpy arrays."""
    keys = np.array(list(data.keys()), dtype=int)        # (N, 2)
    m = keys[:, 0].astype(float)
    n = keys[:, 1].astype(float)
    G  = np.array([data[(int(a), int(b))]["conn"]     for a, b in keys], float)
    eG = np.array([data[(int(a), int(b))]["conn_err"] for a, b in keys], float)
    x = m + 0.5 * n
    y = _SQRT3_2 * n
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    return dict(m=m, n=n, x=x, y=y, r=r, th=th, G=G, eG=eG)


# ---------------------------------------------------------------------------
# Candidate cost functions   (each takes ref_data, test_data → float)
# ---------------------------------------------------------------------------
LX, LY, TX, TY = GEOM["Lx"], GEOM["Ly"], GEOM["Tx"], GEOM["Ty"]


def cost_native(ref, test):
    c, _, _, _ = cost_mod.l2_cost(
        ref, test, LX, LY, TX, TY, LX, LY, TX, TY,
        cost_mode="test_native", cost_power=2)
    return c


def cost_spectral(ref, test):
    c, _, _, _ = cost_mod.l2_cost(
        ref, test, LX, LY, TX, TY, LX, LY, TX, TY,
        cost_mode="spectral",
        spectral_kwargs=dict(n_r=6, n_theta=12, k_max=1.6))
    return c


def cost_chi2_native(ref, test):
    """test_native but residual divided by sqrt(σ_t² + σ_r²)."""
    iref     = _tile_interp(ref, LX, LY, TX, TY, "conn",     copies=2)
    iref_err = _tile_interp(ref, LX, LY, TX, TY, "conn_err", copies=2)
    test_paths = boundary_paths(LX, LY, TX, TY)
    ref_paths  = boundary_paths(LX, LY, TX, TY)
    per = []
    for (rdm, rdn), (tdm, tdn) in zip(ref_paths, test_paths):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        if len(ks) < 2:
            continue
        N = len(ks)
        Gt, et, ts = [], [], []
        for k, mk, nk in zip(ks, ms, ns):
            entry = _lookup_test_value(test, mk, nk, LX, LY, TX, TY, copies=2)
            if entry is None:
                continue
            Gt.append(entry["conn"]); et.append(entry["conn_err"]); ts.append(k / N)
        if len(Gt) < 2:
            continue
        ts = np.asarray(ts); Gt = np.asarray(Gt); et = np.asarray(et)
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts = np.column_stack([ts * rex, ts * rey])
        Gr = np.asarray(iref(pts), float)
        er = np.abs(np.asarray(iref_err(pts), float))
        m = np.isfinite(Gr) & np.isfinite(Gt) & np.isfinite(er) & np.isfinite(et)
        if m.sum() < 2:
            continue
        var = er[m]**2 + et[m]**2 + 1e-30
        per.append(float(np.mean((Gt[m] - Gr[m])**2 / var)))
    return float(np.mean(per)) if per else 0.0


def _cycle_phi(data):
    """Wilson-loop integrals along the three boundary period vectors."""
    out = []
    for tdm, tdn in boundary_paths(LX, LY, TX, TY):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        s = 0.0
        for k, mk, nk in zip(ks, ms, ns):
            e = _lookup_test_value(data, mk, nk, LX, LY, TX, TY, copies=2)
            if e is not None:
                s += e["conn"]
        out.append(s)
    return np.array(out, dtype=float)


def cost_cycle_loop(ref, test):
    pr, pt = _cycle_phi(ref), _cycle_phi(test)
    return float(np.mean((pt - pr) ** 2))


def cost_moment_r2(ref, test):
    ar, at = _arrays(ref), _arrays(test)
    # Use |G| as weight to keep the mean well-defined when G changes sign.
    w_r = np.abs(ar["G"]); w_t = np.abs(at["G"])
    mr = float(np.sum(ar["r"]**2 * w_r) / max(np.sum(w_r), 1e-30))
    mt = float(np.sum(at["r"]**2 * w_t) / max(np.sum(w_t), 1e-30))
    return (mt - mr) ** 2


def _eff_xi_per_dir(data):
    """Slope of log|G| vs lattice distance along each boundary direction."""
    xis = []
    for tdm, tdn in boundary_paths(LX, LY, TX, TY):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        rs, gs = [], []
        for k, mk, nk in zip(ks, ms, ns):
            e = _lookup_test_value(data, mk, nk, LX, LY, TX, TY, copies=2)
            if e is None:
                continue
            x_, y_ = (mk + 0.5 * nk), _SQRT3_2 * nk
            rr = math.hypot(x_, y_)
            if rr < 1e-9:
                continue
            g = abs(e["conn"])
            if g <= 0:
                continue
            rs.append(rr); gs.append(math.log(g))
        if len(rs) < 3:
            xis.append(0.0); continue
        rs = np.asarray(rs); gs = np.asarray(gs)
        # Use only points up to half the max distance (small-r dominates ξ).
        cut = rs <= 0.6 * rs.max()
        if cut.sum() < 3:
            cut = slice(None)
        slope, _ = np.polyfit(rs[cut], gs[cut], 1)
        xis.append(-1.0 / slope if slope < 0 else 0.0)
    return np.array(xis, dtype=float)


def cost_eff_xi(ref, test):
    xr, xt = _eff_xi_per_dir(ref), _eff_xi_per_dir(test)
    return float(np.mean((xt - xr) ** 2))


def _angular_iso_stat(data, m_max=4, n_shells=6):
    """Sum over angular Fourier modes m=2..m_max of |Ĝ_m(r)|² over shells."""
    a = _arrays(data)
    r, th, G = a["r"], a["th"], a["G"]
    if r.size == 0:
        return 0.0
    edges = np.linspace(0.0, r.max() * 1.001, n_shells + 1)
    total = 0.0
    for i in range(n_shells):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if m.sum() < 4:
            continue
        for k in range(2, m_max + 1):
            ck = float(np.mean(G[m] * np.cos(k * th[m])))
            sk = float(np.mean(G[m] * np.sin(k * th[m])))
            total += ck * ck + sk * sk
    return total


def cost_angular_iso(ref, test):
    sr = _angular_iso_stat(ref); st = _angular_iso_stat(test)
    return (st - sr) ** 2


def cost_cycle_var(ref, test):
    """Variance over the three Φ_d cycle integrals (single-sample stat)."""
    vr = float(np.var(_cycle_phi(ref)))
    vt = float(np.var(_cycle_phi(test)))
    return (vt - vr) ** 2


def cost_affine_fit(ref, test):
    """Fit a, b in a·G_t + b ≈ G_r over the all-to-all sites that exist
    in both dicts; return the residual after the best affine collapse.
    Catches the case where the test correlator is just an overall
    rescaling of the reference (which native cost would falsely flag).
    """
    common = sorted(set(ref.keys()) & set(test.keys()))
    if len(common) < 4:
        return 0.0
    Gr = np.array([ref[k]["conn"]  for k in common], float)
    Gt = np.array([test[k]["conn"] for k in common], float)
    A = np.column_stack([Gt, np.ones_like(Gt)])
    sol, *_ = np.linalg.lstsq(A, Gr, rcond=None)
    a, b = sol
    res = Gr - (a * Gt + b)
    return float(np.mean(res ** 2))


# ---------------------------------------------------------------------------
# Cost #2: replicate-baseline subtraction
# It needs access to the OTHER replicates in the group, so it can't be
# expressed as a pure pair-cost.  We post-process the matched/mismatch
# arrays for native cost: subtract the matched mean (best estimate of
# the pure-noise floor).  This is *exactly* a centering — it makes
# matched have zero mean so Z is unchanged, but it shifts the operating
# point a stronger optimizer cares about.  The discrimination metric
# we report is the same Z as for native (definition is invariant).
# We instead implement a stronger version below: at each pair (i, j),
# we subtract the average of cost(eq_i, eq_other) and cost(test_j,
# test_other).  This DOES change Z when noise is correlated across pairs.
# ---------------------------------------------------------------------------

def baseline_subtract_arrays(matched_pairs_costs,
                             mismatched_pairs_costs,
                             matched_idx_pairs,
                             mismatched_idx_pairs,
                             matched_self_avg,         # mean self-cost per ref replicate
                             mismatched_self_avg_a,    # mean self-cost per ref replicate
                             mismatched_self_avg_b):   # mean self-cost per test replicate
    """Per-pair: cost(i,j) - 0.5*(self_avg_a[i] + self_avg_b[j])."""
    out_m = np.array([
        c - 0.5 * (matched_self_avg[i] + matched_self_avg[j])
        for c, (i, j) in zip(matched_pairs_costs, matched_idx_pairs)])
    out_x = np.array([
        c - 0.5 * (mismatched_self_avg_a[i] + mismatched_self_avg_b[j])
        for c, (i, j) in zip(mismatched_pairs_costs, mismatched_idx_pairs)])
    return out_m, out_x


# ---------------------------------------------------------------------------
# Pair-cost driver
# ---------------------------------------------------------------------------

class _silence:
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        return self
    def __exit__(self, *a):
        sys.stdout.close(); sys.stdout = self._stdout


def _run_cost(fn, group_a, group_b, *, same_group):
    """Return (costs_array, idx_pairs) for every pair across two groups."""
    if same_group:
        idx = list(itertools.combinations(range(len(group_a)), 2))
        pairs = [(group_a[i], group_a[j]) for i, j in idx]
    else:
        idx = [(i, j) for i in range(len(group_a))
                       for j in range(len(group_b))]
        pairs = [(group_a[i], group_b[j]) for i, j in idx]
    out = []
    for ref, tst in pairs:
        with _silence():
            out.append(fn(ref, tst))
    return np.array(out, float), idx


def _z_score(matched, mismatched):
    sd = np.std(matched, ddof=1) if len(matched) > 1 else 0.0
    if sd <= 0:
        return float("inf") if np.mean(mismatched) > np.mean(matched) else 0.0
    return float((np.mean(mismatched) - np.mean(matched)) / sd)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    eq = _load_group("equil_12")
    a2 = _load_group("aniso_r2_12")
    a5 = _load_group("aniso_r5_12")

    candidates = [
        ("native",        cost_native),
        ("spectral",      cost_spectral),
        ("chi2_native",   cost_chi2_native),
        ("cycle_loop",    cost_cycle_loop),
        ("moment_r2",     cost_moment_r2),
        ("eff_xi",        cost_eff_xi),
        ("angular_iso",   cost_angular_iso),
        ("cycle_var",     cost_cycle_var),
        ("affine_fit",    cost_affine_fit),
    ]

    groups = {"eq": eq, "a2": a2, "a5": a5}
    # Three "reference choices": for each, the matched group is the
    # diagonal pair-set, and the mismatched groups are the two off-
    # diagonal cross-sets.  An honest distance-like cost should give
    # roughly symmetric Z across all three reference choices.  A cost
    # that only fires on isotropy-violation gives huge Z only when ref
    # is the symmetric (eq) group.
    suites = [
        ("ref=eq", "eq", ["a2", "a5"]),
        ("ref=a2", "a2", ["eq", "a5"]),
        ("ref=a5", "a5", ["eq", "a2"]),
    ]

    results = {}   # results[suite_label][cand_name] = {Z_x, Z_y, ...}
    raw = {}
    print("\n[run] pair costs for each candidate × reference suite ...")
    for suite_label, ref_key, mis_keys in suites:
        results[suite_label] = {}
        raw[suite_label] = {}
        ref_grp = groups[ref_key]
        mis_grps = [groups[k] for k in mis_keys]
        for name, fn in candidates:
            tA = time.time()
            m_arr,  _ = _run_cost(fn, ref_grp, ref_grp,    same_group=True)
            x_arrs = [_run_cost(fn, ref_grp, g, same_group=False)[0]
                      for g in mis_grps]
            zs = [_z_score(m_arr, x) for x in x_arrs]
            results[suite_label][name] = dict(
                matched_mean=float(m_arr.mean()),
                matched_std =float(m_arr.std(ddof=1)),
                mis_keys=mis_keys,
                mis_means=[float(x.mean()) for x in x_arrs],
                Z=zs, t_s=time.time() - tA,
            )
            raw[suite_label][name] = dict(
                matched=m_arr.tolist(),
                **{f"vs_{k}": x.tolist() for k, x in zip(mis_keys, x_arrs)})
        # Print a small per-suite table.
        print(f"\n  --- suite {suite_label} (matched = {ref_key}-{ref_key}) ---")
        hdr = f"  {'cand':14s}  matched μ±σ           Z(vs {mis_keys[0]})  Z(vs {mis_keys[1]})"
        print(hdr); print("  " + "-" * (len(hdr)-2))
        for name, _ in candidates:
            r = results[suite_label][name]
            print(f"  {name:14s}  {r['matched_mean']:+.3e}±{r['matched_std']:.2e}"
                  f"  {r['Z'][0]:>10.3f}  {r['Z'][1]:>10.3f}")

    # ------------------ #2: replicate-baseline subtraction (ref=eq only) ---
    print("\n[run] replicate-baseline (#2) post-processing on native cost...")
    fn = cost_native

    def _self_avg(grp):
        n = len(grp)
        sums = np.zeros(n); cnts = np.zeros(n)
        with _silence():
            for i in range(n):
                for j in range(i + 1, n):
                    c = fn(grp[i], grp[j])
                    sums[i] += c; cnts[i] += 1
                    sums[j] += c; cnts[j] += 1
        return sums / np.maximum(cnts, 1)

    eq_self = _self_avg(eq); a2_self = _self_avg(a2); a5_self = _self_avg(a5)
    self_avgs = {"eq": eq_self, "a2": a2_self, "a5": a5_self}

    for suite_label, ref_key, mis_keys in suites:
        m_idx = list(itertools.combinations(range(len(groups[ref_key])), 2))
        nat_m = np.array(raw[suite_label]["native"]["matched"])
        ref_self = self_avgs[ref_key]
        m_corr = np.array([
            c - 0.5 * (ref_self[i] + ref_self[j])
            for c, (i, j) in zip(nat_m, m_idx)])
        zs = []
        mis_means = []
        x_corrs = {}
        for k in mis_keys:
            x = np.array(raw[suite_label]["native"][f"vs_{k}"])
            x_idx = [(i, j) for i in range(len(groups[ref_key]))
                             for j in range(len(groups[k]))]
            tst_self = self_avgs[k]
            x_corr = np.array([
                c - 0.5 * (ref_self[i] + tst_self[j])
                for c, (i, j) in zip(x, x_idx)])
            x_corrs[k] = x_corr
            zs.append(_z_score(m_corr, x_corr))
            mis_means.append(float(x_corr.mean()))
        results[suite_label]["baseline_subtract"] = dict(
            matched_mean=float(m_corr.mean()),
            matched_std =float(m_corr.std(ddof=1)),
            mis_keys=mis_keys, mis_means=mis_means, Z=zs)
        raw[suite_label]["baseline_subtract"] = dict(
            matched=m_corr.tolist(),
            **{f"vs_{k}": x.tolist() for k, x in x_corrs.items()})
        print(f"  {suite_label}  baseline_subtract  "
              f"Z(vs {mis_keys[0]})={zs[0]:.3f}  "
              f"Z(vs {mis_keys[1]})={zs[1]:.3f}")

    # ------------------ Save + plot ----------------------------------------
    summary = dict(geom=GEOM, suites=results,
                   notes={"jackknife": "N/A — needs per-config block data",
                          "offdiag_cov": "N/A — needs per-config covariance"},
                   wall_time_s=time.time() - t0)
    with open(OUT_DIR / "zoo_summary.json", "w") as f:
        json.dump(dict(summary, raw=raw), f, indent=2)
    print(f"\n[save] {OUT_DIR/'zoo_summary.json'}")

    _plot_suites(results)
    _plot_distributions(raw["ref=eq"])
    print(f"[done] wall {time.time()-t0:.1f}s")


def _plot_suites(results):
    """Bar chart: x=candidate, group bars by suite."""
    suite_order = ["ref=eq", "ref=a2", "ref=a5"]
    cand_names = list(results[suite_order[0]].keys())
    n_cand = len(cand_names)
    n_suite = len(suite_order)
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)

    for ax, mis_idx, label in [(axes[0], 0, "first  off-diagonal mismatch"),
                                (axes[1], 1, "second off-diagonal mismatch")]:
        x = np.arange(n_cand)
        w = 0.27
        for s_i, suite in enumerate(suite_order):
            zs = []
            for cn in cand_names:
                r = results[suite][cn]
                zs.append(abs(r["Z"][mis_idx]) if mis_idx < len(r["Z"]) else 0.0)
            zs = np.array([max(z, 1e-3) for z in zs])
            offset = (s_i - 1) * w
            color = ["C0", "C2", "C3"][s_i]
            mis_key = results[suite][cand_names[0]]["mis_keys"][mis_idx]
            ax.bar(x + offset, zs, width=w,
                   label=f"{suite}  vs {mis_key}", color=color, alpha=0.85)
            for xi, zi in zip(x + offset, zs):
                ax.text(xi, zi * 1.18,
                        f"{zi:.2f}" if zi < 1e3 else f"{zi:.0e}",
                        ha="center", va="bottom", fontsize=6.5)
        ax.set_yscale("log")
        ax.set_ylabel("|Z|  (log)")
        ax.set_title(f"{label}  —  cost discrimination across reference choices",
                     fontsize=10)
        ax.grid(alpha=0.25, axis="y", which="both")
        ax.legend(fontsize=8, ncol=3, loc="upper left")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cand_names, rotation=25, ha="right")
    fig.suptitle("Cost-function zoo — discrimination across THREE reference choices\n"
                 "(an honest distance-like cost shows similar |Z| for all 3 ref choices;\n"
                 " a cost that only fires on isotropy-violation peaks only at ref=eq)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = OUT_DIR / "zoo_compare.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[save] {out}")


def _plot_z_bars(rows):
    names = list(rows.keys())
    z2 = np.array([rows[n]["Z_r2"] for n in names], float)
    z5 = np.array([rows[n]["Z_r5"] for n in names], float)
    # Plot |Z| on a symmetric log scale so the 5-decade span is readable.
    z2_p = np.where(np.abs(z2) > 0, np.abs(z2), 1e-3)
    z5_p = np.where(np.abs(z5) > 0, np.abs(z5), 1e-3)
    x  = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(13, 5.5))
    bars2 = ax.bar(x - 0.18, z2_p, width=0.36, label="|Z(r2)|  mild mismatch",
                   color=["C0" if v >= 0 else "lightgray" for v in z2])
    bars5 = ax.bar(x + 0.18, z5_p, width=0.36, label="|Z(r5)|  strong mismatch",
                   color=["C3" if v >= 0 else "lightgray" for v in z5])
    for i, (a, b) in enumerate(zip(z2, z5)):
        ax.text(i - 0.18, max(abs(a), 1e-3) * 1.18,
                f"{a:+.2f}" if abs(a) < 1e3 else f"{a:+.1e}",
                ha="center", va="bottom", fontsize=7, rotation=0)
        ax.text(i + 0.18, max(abs(b), 1e-3) * 1.18,
                f"{b:+.2f}" if abs(b) < 1e3 else f"{b:+.1e}",
                ha="center", va="bottom", fontsize=7, rotation=0)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_yscale("log")
    nat_z2 = abs(rows["native"]["Z_r2"]) or 1e-3
    nat_z5 = abs(rows["native"]["Z_r5"]) or 1e-3
    ax.axhline(nat_z2, color="C0", lw=0.7, ls="--", alpha=0.6,
               label=f"native |Z(r2)|={nat_z2:.2f}")
    ax.axhline(nat_z5, color="C3", lw=0.7, ls="--", alpha=0.6,
               label=f"native |Z(r5)|={nat_z5:.2f}")
    ax.set_ylabel("|Z|  =  |μ_mis - μ_match| / σ_match    (log scale)")
    ax.set_title("Cost-function zoo — discrimination on L=12 stability data\n"
                 "(higher = optimizer sees true mismatch above MC noise; "
                 "gray bars = sign-flipped Z)")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.grid(alpha=0.25, axis="y", which="both")
    fig.tight_layout()
    out = OUT_DIR / "zoo_compare.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[save] {out}")


def _plot_distributions(raw):
    """raw is the per-suite dict for a single suite (e.g. raw['ref=eq'])."""
    names = list(raw.keys())
    n = len(names)
    cols = 4
    rows_n = math.ceil(n / cols)
    fig, axs = plt.subplots(rows_n, cols, figsize=(4 * cols, 3 * rows_n))
    axs = np.array(axs).reshape(-1)
    for ax, name in zip(axs, names):
        d = raw[name]
        cols_d = ["matched"] + sorted(k for k in d if k.startswith("vs_"))
        labels = [k.replace("vs_", "") if k != "matched" else "match"
                  for k in cols_d]
        ax.boxplot([d[c] for c in cols_d], tick_labels=labels, showmeans=True)
        all_vals = np.concatenate([d[c] for c in cols_d])
        if np.all(all_vals > 0):
            ax.set_yscale("log")
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3, which="both")
    for ax in axs[n:]:
        ax.set_axis_off()
    fig.suptitle("Per-candidate cost distributions  (suite ref=eq)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = OUT_DIR / "zoo_distributions.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
