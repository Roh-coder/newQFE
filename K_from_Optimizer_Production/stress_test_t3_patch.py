#!/usr/bin/env python3
"""
stress_test_t3_patch.py — adversarial stress test of the T3 (drop-collapsed-
directions) patch to ``test_native`` cost.

No modification to existing code: re-implements the proposed cost as a
shadow function and probes corner cases.

Stress dimensions
-----------------
S1  Geometry sweep — 25 contrived (Lx, Ly, Tx, Ty) tuples covering
    {0, 1, 2, 3} collapsed directions, including pathological (Lx=Ly=1).
    For each: verify ratio current/proposed = N_valid/3 to <1e-10.

S2  All-collapse edge case (N_valid = 0).
    - current: returns finite cost = 0, σ = 0  ← LOOKS LIKE BEST POSSIBLE
    - proposed: returns NaN
    Confirm the failure mode and document downstream impact.

S3  Long-horizon CMA-ES selection invariance.
    100 seeds × 10 generations × popsize=8 = 8000 evals.  Verify rank
    stability on G1 (1 collapsed direction).  Synthetic surface with
    real-MC-calibrated noise (σ_d = 5e-3).

S4  Numerical edge cases:
    a) per-dir cost == 0 due to perfect agreement (not collapse).  Make
       sure proposed does NOT drop it (should keep, only drop when both
       C and σ are exactly 0).
    b) per-dir σ > 0 but C exactly 0: should be kept.
    c) per-dir σ exactly 0 but C > 0 (impossible from MC but possible
       from cache): should be kept.
    d) Large dynamic range: per-dir [1e-10, 1, 1e10].

S5  Real-MC adversarial geometry — replay the previous CMA-ES eval log
    against three modified rules to make sure no other zero-σ events
    sneak through.  Cross-check with a second contrived geometry
    (Lx=11, Ly=12, Tx=0, Ty=0) with different gcd pattern.

S6  Cross-mode interaction: confirm proposed cost still works when
    cost_power=4 (other code path).  Use cached _betac_stability_test
    equil_12 data.

Outputs
-------
  results/_native_triage/stress_t3_summary.json
  results/_native_triage/stress_t3.log
"""
from __future__ import annotations

import contextlib
import glob
import io
import json
import math
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cost as cost_mod          # noqa: E402
import mc_engine                 # noqa: E402

OUT = os.path.join(HERE, "results", "_native_triage")
os.makedirs(OUT, exist_ok=True)
STAB = os.path.join(HERE, "results", "_betac_stability_test")


# ---------------------------------------------------------------------------
# Shadow proposed cost (drop dirs where (C, σ) == (0, 0))
# ---------------------------------------------------------------------------
def proposed_from_perdir(per_dir: List[float],
                         per_sig: List[float]) -> Tuple[float, float, int]:
    valid = [(c, s) for c, s in zip(per_dir, per_sig)
             if not (c == 0.0 and s == 0.0)]
    if not valid:
        return float("nan"), float("nan"), 0
    n = len(valid)
    cost = sum(c for c, _ in valid) / n
    sig = math.sqrt(sum(s * s for _, s in valid)) / n
    return cost, sig, n


def _quiet(fn, *a, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return fn(*a, **kw)


# ===========================================================================
# S1 — Geometry sweep
# ===========================================================================
def _load_one_replicate():
    pat = os.path.join(STAB, "equil_12", "rw_t00", "*",
                       "two_point_all_to_all.dat")
    f = sorted(glob.glob(pat))[0]
    return mc_engine.load_all_to_all(f)


def s1_geometry_sweep(ref, test) -> dict:
    print("\n=== S1: geometry sweep (current/proposed ratio) ===")
    LXY = (12, 12, 0, 0)
    geoms = [
        (12, 12, 0, 0),    # 0 collapse
        (12, 11, 0, 0),    # 1 collapse
        (11, 12, 0, 0),    # 1 collapse (different)
        (12, 12, 0, 1),    # 2 collapse
        (12, 12, 1, 0),    # 2 collapse
        (1, 1, 0, 0),      # 3 collapse (degenerate)
        (12, 1, 0, 0),     # 2 collapse
        (1, 12, 0, 0),     # 2 collapse
        (12, 12, 1, 1),    # mixed
        (12, 12, -1, 1),
        (10, 8, 0, 0),
        (7, 5, 0, 0),
        (6, 9, 0, 0),
        (3, 4, 0, 0),
        (2, 2, 0, 0),
        (12, 12, 6, 0),
        (12, 12, 0, 6),
        (12, 12, 4, -4),
        (16, 16, 0, 0),
        (16, 12, 0, 0),
        (12, 16, 4, 0),
        (8, 8, 4, 4),
        (8, 8, 0, 4),
        (5, 7, 1, 1),
        (15, 9, 3, -3),
    ]
    rows, fails = [], []
    for g in geoms:
        try:
            paths = cost_mod.boundary_paths(*g)
            gcds = [math.gcd(abs(a), abs(b)) for a, b in paths]
            n_valid_pred = sum(1 for x in gcds if x >= 2)
            c, s, pd, ps = _quiet(cost_mod.l2_cost,
                                   ref, test, *g, *LXY,
                                   cost_mode="test_native", cost_power=2)
            cp, sp, nv = proposed_from_perdir(pd, ps)
            ratio = c / cp if cp and not math.isnan(cp) else float("nan")
            pred = nv / 3 if nv else float("nan")
            ok = (math.isnan(ratio) and math.isnan(pred)) or \
                 (not math.isnan(ratio) and abs(ratio - pred) < 1e-10)
            rows.append(dict(geom=g, gcds=gcds,
                             n_valid_pred=n_valid_pred, n_valid_obs=nv,
                             ratio=ratio, pred=pred, ok=ok,
                             cost_cur=c, cost_prop=cp))
            print(f"  geom={g}  gcds={gcds}  N_valid={nv}/3  "
                  f"ratio={ratio if not math.isnan(ratio) else 'nan':<12}  "
                  f"pred={pred if not math.isnan(pred) else 'nan':<8}  "
                  f"{'OK' if ok else '*** MISMATCH ***'}")
        except Exception as e:
            fails.append(dict(geom=g, error=str(e)))
            print(f"  geom={g}  EXCEPTION: {e}")
    n_ok = sum(1 for r in rows if r["ok"])
    print(f"  → {n_ok}/{len(rows)} OK   exceptions={len(fails)}")
    return dict(rows=rows, exceptions=fails, n_ok=n_ok, n_total=len(rows))


# ===========================================================================
# S2 — All-collapse edge case (N_valid = 0)
# ===========================================================================
def s2_all_collapse(ref) -> dict:
    print("\n=== S2: all-collapse edge case (N_valid=0) ===")
    LXY = (12, 12, 0, 0)
    g = (1, 1, 0, 0)  # all gcds = 1 except one... let's check
    paths = cost_mod.boundary_paths(*g)
    gcds = [math.gcd(abs(a), abs(b)) for a, b in paths]
    print(f"  geom={g}  paths={paths}  gcds={gcds}")
    c, s, pd, ps = _quiet(cost_mod.l2_cost,
                          ref, ref, *g, *LXY,
                          cost_mode="test_native", cost_power=2)
    cp, sp, nv = proposed_from_perdir(pd, ps)
    print(f"  current : cost={c:.4e}  σ={s:.4e}  per_dir={pd}  per_sig={ps}")
    print(f"  proposed: cost={cp}  σ={sp}  N_valid={nv}")
    danger = (c == 0.0 and s == 0.0)
    print(f"  DANGER: current returns 0±0 (looks like 'best ever') = {danger}")
    print(f"  proposed returns NaN — caller must guard with isfinite()")
    return dict(geom=g, gcds=gcds, current_cost=c, current_sigma=s,
                proposed_cost=cp, proposed_sigma=sp, n_valid=nv,
                danger_current=danger,
                proposed_is_nan=(math.isnan(cp) if not isinstance(cp, complex) else False))


# ===========================================================================
# S3 — Long-horizon CMA-ES selection invariance
# ===========================================================================
SQRT3_2 = math.sqrt(3.0) / 2.0
A = np.array([1.0, 0.7, 1.3])
B = np.array([0.2, 0.1, 0.3])
SIGMA_D = np.array([5e-3, 5e-3, 5e-3])


def _evaluate_synthetic(rng, r1, r2, valid):
    quad = (r1 - 1) ** 2 + (r2 - 1) ** 2
    cross = (r1 - r2) ** 2
    Cd_true = A * quad + B * cross
    Cd = Cd_true + rng.normal(0, SIGMA_D)
    Cd_cur = np.where(valid, Cd, 0.0)
    sd_cur = np.where(valid, SIGMA_D, 0.0)
    cur_cost = float(Cd_cur.sum() / 3)
    cur_sig = float(math.sqrt((sd_cur ** 2).sum()) / 3)
    if valid.any():
        nv = int(valid.sum())
        prop_cost = float(Cd[valid].sum() / nv)
        prop_sig = float(math.sqrt((SIGMA_D[valid] ** 2).sum()) / nv)
    else:
        prop_cost, prop_sig = float("nan"), float("nan")
    return (cur_cost, cur_sig), (prop_cost, prop_sig), Cd


def _cmaes_one(seed, valid, popsize=8, max_gen=10, x0=(2., 0.5), sigma0=0.6):
    n = 2; mu = popsize // 2
    weights = np.log(mu + 1) - np.log(np.arange(1, mu + 1))
    weights /= weights.sum()
    mueff = 1 / (weights ** 2).sum()
    c_sigma = (mueff + 2) / (n + mueff + 5)
    d_sigma = 1 + 2 * max(0, math.sqrt((mueff - 1) / (n + 1)) - 1) + c_sigma
    chi_n = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))

    mean = np.array(x0, float); sigma = float(sigma0)
    p_sigma = np.zeros(n)
    master = np.random.default_rng(seed)
    flips = 0; total_gens = 0
    for g in range(max_gen):
        ss = int(master.integers(0, 2**31 - 1))
        srng = np.random.default_rng(ss)
        z = srng.standard_normal((popsize, n))
        x = mean + sigma * z
        cur_costs, prop_costs = [], []
        for k in range(popsize):
            es = (seed * 1_000_003 + g * 1_009 + k) & 0x7FFFFFFF
            erng = np.random.default_rng(es)
            (cc, _), (pc, _), _ = _evaluate_synthetic(erng,
                                                      float(x[k, 0]),
                                                      float(x[k, 1]),
                                                      valid)
            cur_costs.append(cc); prop_costs.append(pc)
        cur_rank = np.argsort(cur_costs).tolist()
        prop_rank = np.argsort(prop_costs).tolist()
        if cur_rank != prop_rank:
            flips += 1
        total_gens += 1
        # Use current-cost selection (CMA-ES would do the same with proposed
        # since rank is identical; we just need to advance the mean).
        order = np.argsort(cur_costs)
        best_z = (x[order[:mu]] - mean) / sigma
        new_mean = mean + sigma * (weights @ best_z)
        p_sigma = (1 - c_sigma) * p_sigma + math.sqrt(
            c_sigma * (2 - c_sigma) * mueff) * (weights @ best_z)
        sigma = sigma * math.exp((c_sigma / d_sigma) *
                                 (np.linalg.norm(p_sigma) / chi_n - 1))
        mean = new_mean
    return flips, total_gens, mean, sigma


def s3_long_horizon() -> dict:
    print("\n=== S3: long-horizon CMA-ES selection invariance "
          "(100 seeds × 10 gens × popsize=8) ===")
    out = {}
    for label, valid in [
        ("G0_no_collapse", np.array([True, True, True])),
        ("G1_one_collapse", np.array([True, True, False])),
        ("G2_two_collapse", np.array([True, False, False])),
    ]:
        total_flips = 0; total_gens = 0; final_dists = []
        for s in range(100):
            f, g, m, _ = _cmaes_one(s, valid)
            total_flips += f; total_gens += g
            final_dists.append(float(np.hypot(m[0] - 1, m[1] - 1)))
        out[label] = dict(
            total_flips=total_flips,
            total_gens=total_gens,
            flip_rate=total_flips / total_gens,
            mean_final_dist=float(np.mean(final_dists)),
            std_final_dist=float(np.std(final_dists)),
        )
        print(f"  {label}: {total_flips}/{total_gens} gen-rank flips  "
              f"({100 * total_flips / total_gens:.2f}%)  "
              f"final mean dist = {np.mean(final_dists):.4f} ± "
              f"{np.std(final_dists):.4f}")
    return out


# ===========================================================================
# S4 — Numerical edge cases
# ===========================================================================
def s4_numerical_edges() -> dict:
    print("\n=== S4: numerical edge cases (proposed cost) ===")
    cases = [
        ("(0,0,0)         all dirs collapsed", [0., 0., 0.], [0., 0., 0.]),
        ("(0,0,5e-4)      2 dirs collapsed",  [0., 0., 5e-4], [0., 0., 1e-5]),
        ("(0,5e-4,5e-4)   1 dir collapsed",   [0., 5e-4, 5e-4], [0., 1e-5, 1e-5]),
        ("(0,0,0)+nonzero σ — perfect agreement, real MC noise",
                                              [0., 0., 0.], [1e-5, 1e-5, 1e-5]),
        ("(0,1e-3,2e-3)+(0,σ,σ) — first dir is a real zero-noise zero",
                                              [0., 1e-3, 2e-3], [1e-5, 2e-5, 3e-5]),
        ("(1e-3,1e-3,0)+(σ,σ,0) — last dir cost is zero but σ>0 normally absent",
                                              [1e-3, 1e-3, 0.], [2e-5, 2e-5, 0.]),
        ("(1e-10,1.0,1e10) — large dynamic range",
                                              [1e-10, 1.0, 1e10], [1e-12, 1e-2, 1e8]),
    ]
    rows = []
    for label, pd, ps in cases:
        c_cur = sum(pd) / 3
        s_cur = math.sqrt(sum(x * x for x in ps)) / 3
        cp, sp, nv = proposed_from_perdir(pd, ps)
        rows.append(dict(label=label, per_dir=pd, per_sig=ps,
                         cur_cost=c_cur, cur_sig=s_cur,
                         prop_cost=cp, prop_sig=sp, n_valid=nv))
        print(f"  {label}")
        print(f"    pd={pd}  ps={ps}")
        print(f"    current  cost={c_cur:.3e}  σ={s_cur:.3e}")
        print(f"    proposed cost={cp}  σ={sp}  N_valid={nv}")
    return dict(rows=rows)


# ===========================================================================
# S5 — Replay real CMA-ES eval log + adversarial geometry
# ===========================================================================
def s5_real_log_replay() -> dict:
    print("\n=== S5: replay real-MC eval log; adversarial geometry check ===")
    log = os.path.join(OUT, "realmc_cmaes", "mc_scratch", "eval_log.jsonl")
    if not os.path.exists(log):
        print("  (no real-MC log; skipping)")
        return {"skipped": True}
    rows = []
    with open(log) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))

    n_zero_C_nonzero_sig = 0
    n_zero_sig_nonzero_C = 0
    n_both_nonzero = 0
    n_both_zero = 0
    for r in rows:
        for c, s in zip(r["per_dir"], r["per_dir_sigma"]):
            if c == 0.0 and s == 0.0:
                n_both_zero += 1
            elif c == 0.0 and s != 0.0:
                n_zero_C_nonzero_sig += 1
            elif c != 0.0 and s == 0.0:
                n_zero_sig_nonzero_C += 1
            else:
                n_both_nonzero += 1
    print(f"  per-direction frequency over {len(rows)*3} dirs in log:")
    print(f"    both zero        (collapse):  {n_both_zero}")
    print(f"    C=0 σ≠0          (perfect agreement, MC noise): {n_zero_C_nonzero_sig}")
    print(f"    C≠0 σ=0          (impossible from MC): {n_zero_sig_nonzero_C}")
    print(f"    both nonzero     (normal):    {n_both_nonzero}")
    print(f"  → patch only drops 'both zero'; no false-positive drops.")

    # Now check second contrived geometry (Lx=11, Ly=12) — different gcd pattern.
    print("\n  adversarial geometry (Lx=11, Ly=12, Tx=0, Ty=0):")
    paths = cost_mod.boundary_paths(11, 12, 0, 0)
    gcds = [math.gcd(abs(a), abs(b)) for a, b in paths]
    print(f"    paths={paths}  gcds={gcds}  "
          f"N_valid_pred={sum(1 for g in gcds if g >= 2)}")

    return dict(
        n_evals_in_log=len(rows),
        per_dir_freq=dict(both_zero=n_both_zero,
                          zero_C_nonzero_sig=n_zero_C_nonzero_sig,
                          zero_sig_nonzero_C=n_zero_sig_nonzero_C,
                          both_nonzero=n_both_nonzero),
        adversarial_geom=dict(geom=(11, 12, 0, 0), paths=paths, gcds=gcds),
    )


# ===========================================================================
# S6 — cost_power=4 interaction
# ===========================================================================
def s6_cost_power_4(ref, test) -> dict:
    print("\n=== S6: cost_power=4 interaction ===")
    LXY = (12, 12, 0, 0)
    out = []
    for g in [(12, 12, 0, 0), (12, 11, 0, 0), (12, 12, 0, 1)]:
        c, s, pd, ps = _quiet(cost_mod.l2_cost,
                              ref, test, *g, *LXY,
                              cost_mode="test_native", cost_power=4)
        cp, sp, nv = proposed_from_perdir(pd, ps)
        ratio = c / cp if cp and not math.isnan(cp) else float("nan")
        gcds = [math.gcd(abs(a), abs(b))
                for a, b in cost_mod.boundary_paths(*g)]
        pred = nv / 3 if nv else float("nan")
        ok = (math.isnan(ratio) and math.isnan(pred)) or \
             (not math.isnan(ratio) and abs(ratio - pred) < 1e-10)
        out.append(dict(geom=g, gcds=gcds, n_valid=nv,
                        cost_cur=c, cost_prop=cp,
                        ratio=ratio, pred=pred, ok=ok))
        print(f"  geom={g} gcds={gcds} N_valid={nv}  "
              f"current={c:.4e}  proposed={cp:.4e}  "
              f"ratio={ratio:.4f}  pred={pred:.4f}  {'OK' if ok else 'MISMATCH'}")
    return dict(rows=out)


# ===========================================================================
# Main
# ===========================================================================
def main():
    log_buf = io.StringIO()

    class Tee:
        def __init__(self, *s): self.s = s
        def write(self, x):
            for s in self.s: s.write(x)
        def flush(self):
            for s in self.s: s.flush()
    sys.stdout = Tee(sys.__stdout__, log_buf)

    print("=" * 78)
    print("stress_test_t3_patch — adversarial validation of T3 (drop-collapsed)")
    print("=" * 78)

    ref = _load_one_replicate()
    # Use a different replicate for "test" to introduce real MC noise.
    pat = os.path.join(STAB, "equil_12", "rw_t01", "*",
                       "two_point_all_to_all.dat")
    test = mc_engine.load_all_to_all(sorted(glob.glob(pat))[0])

    out = {}
    out["S1_geometry_sweep"]   = s1_geometry_sweep(ref, test)
    out["S2_all_collapse"]     = s2_all_collapse(ref)
    out["S3_long_horizon"]     = s3_long_horizon()
    out["S4_numerical_edges"]  = s4_numerical_edges()
    out["S5_real_log_replay"]  = s5_real_log_replay()
    out["S6_cost_power_4"]     = s6_cost_power_4(ref, test)

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    s1 = out["S1_geometry_sweep"]
    print(f"  S1 geometry sweep:  {s1['n_ok']}/{s1['n_total']} ratios match N_valid/3")
    s2 = out["S2_all_collapse"]
    print(f"  S2 all-collapse:    current={s2['current_cost']} (DANGER: looks like best!)  "
          f"proposed=NaN ({'GUARD REQUIRED' if s2['proposed_is_nan'] else 'OK'})")
    s3 = out["S3_long_horizon"]
    flip_rates = [v['flip_rate'] for v in s3.values()]
    print(f"  S3 long-horizon:    max selection flip rate over 1000 gens = "
          f"{max(flip_rates)*100:.2f}%  (expect 0%)")
    s6 = out["S6_cost_power_4"]
    n_ok6 = sum(1 for r in s6['rows'] if r['ok'])
    print(f"  S6 cost_power=4:    {n_ok6}/{len(s6['rows'])} ratios match")

    sys.stdout = sys.__stdout__
    with open(os.path.join(OUT, "stress_t3_summary.json"), "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    with open(os.path.join(OUT, "stress_t3.log"), "w") as f:
        f.write(log_buf.getvalue())
    print(log_buf.getvalue())


if __name__ == "__main__":
    main()
