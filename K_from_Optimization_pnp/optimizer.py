"""
optimizer.py — Hand-coded Nelder-Mead simplex minimiser for the
2-D coupling-ratio search (r1, r2).

Standard NM coefficients
------------------------
    alpha    = 1.0    reflection
    beta_exp = 2.0    expansion
    gamma    = 0.5    contraction
    shrink   = 0.75   shrink (larger than scipy's 0.5 for noise tolerance)

Stopping rules
--------------
    1. Hard cap:  ``max_evals``
    2. Geometric: simplex diameter < ``xatol`` AND cost spread < ``fatol``
    3. SNR floor: running best SNR < ``snr_floor`` after ≥ 5 evals
    4. Indistinguishable-from-reference: best SNR < ``indist_stop_snr``
       (only triggered when adaptive stats have been maxed out)

Adaptive MC statistics
----------------------
If ``snr_target > 0``, after every evaluation with SNR < snr_target the
evaluator's ``n_traj_prod`` is multiplied by 1.5 (capped at
``snr_max_traj_factor × starting value``).  This automatically tightens
the cost measurement as the simplex nears the minimum.
"""

from __future__ import annotations

import numpy as np


class _ResultProxy:
    """Lightweight stand-in for evaluator.EvalResult built from a dict.

    Used only by the parallel CMA-ES path: workers return ``asdict()``-ed
    EvalResults, which the main process wraps here so the existing
    ``history`` machinery (sorted by ``r.cost``, indexed by attribute) keeps
    working without importing/copying the dataclass into this module.
    """
    __slots__ = ("__dict__",)

    def __init__(self, d: dict):
        self.__dict__.update(d)


def _summary(history: list, status: str) -> dict:
    if not history:
        return {"status": status, "n_evals": 0}
    best = min(history, key=lambda r: r.cost)
    return {
        "status": status,
        "n_evals": len(history),
        "best_eval_id": best.eval_id,
        "best_r1": best.r1,
        "best_r2": best.r2,
        "best_cost": best.cost,
        "best_sigma": best.sigma_cost,
        "best_snr": best.snr,
        "best_beta_c": best.beta_c,
        "total_traj": sum(r.n_traj_prod + r.n_traj_scan_total
                          for r in history),
        "total_wall_s": sum(r.wall_time_s for r in history),
    }


def run_nelder_mead(evaluator, x0=(1.0, 1.0), *,
                    max_evals=40, sigma0=0.1,
                    xatol=0.005, fatol=1e-6,
                    alpha=1.0, beta_exp=2.0, gamma=0.5, shrink=0.75,
                    snr_floor=0.0,
                    indist_stop_snr=0.0,
                    snr_target=0.0, snr_max_traj_factor=4):
    history: list = []
    _n_traj_prod_base = int(evaluator.n_traj_prod)
    _n_traj_max = _n_traj_prod_base * max(1, int(snr_max_traj_factor))

    def _eval(pt, simplex, slot=-1):
        if len(history) >= max_evals:
            raise StopIteration("max_evals reached")
        if hasattr(evaluator, "current_simplex"):
            display_simplex = [np.asarray(v).tolist() for v in simplex]
            if slot is not None:
                display_simplex[slot] = [float(pt[0]), float(pt[1])]
            evaluator.current_simplex = display_simplex
        res = evaluator(float(pt[0]), float(pt[1]))
        history.append(res)

        # Adaptive stats
        if snr_target > 0 and res.snr < snr_target:
            new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
            if new_n > evaluator.n_traj_prod:
                print(f"[nm] SNR={res.snr:.2f} < {snr_target}; "
                      f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n

        best = min(history, key=lambda r: r.cost)

        if snr_floor > 0 and best.snr < snr_floor and len(history) >= 5:
            print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}; "
                  "need more statistics.")
            raise StopIteration("snr floor reached")

        if indist_stop_snr > 0 and len(history) >= 5 and best.snr < indist_stop_snr:
            if snr_target > 0 and evaluator.n_traj_prod < _n_traj_max:
                new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
                print(f"[nm] best SNR={best.snr:.2f} < {indist_stop_snr} "
                      f"but stats not maxed; boosting "
                      f"n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n
            else:
                print(f"[stop] running best SNR={best.snr:.2f} < {indist_stop_snr}; "
                      f"test indistinguishable from reference at "
                      f"r1={best.r1:.4f} r2={best.r2:.4f}.")
                raise StopIteration("indistinguishable from reference")
        return res.cost

    s = float(sigma0)
    verts = np.array([
        [x0[0],     x0[1]    ],
        [x0[0] + s, x0[1]    ],
        [x0[0],     x0[1] + s],
    ], dtype=float)
    costs = np.full(3, np.inf)

    status = "completed"
    try:
        for i in range(len(verts)):
            costs[i] = _eval(verts[i], verts, slot=None)

        while len(history) < max_evals:
            order = np.argsort(costs)
            verts, costs = verts[order], costs[order]

            diam = max(np.linalg.norm(v - verts[0]) for v in verts[1:])
            if diam < xatol and (costs[-1] - costs[0]) < fatol:
                status = "converged"
                break

            xbar = verts[:-1].mean(axis=0)

            # Reflection
            xr = xbar + alpha * (xbar - verts[-1])
            fr = _eval(xr, verts)

            if costs[0] <= fr < costs[-2]:
                verts[-1], costs[-1] = xr, fr
                continue

            if fr < costs[0]:
                # Expansion
                xe = xbar + beta_exp * (xr - xbar)
                fe = _eval(xe, verts)
                if fe < fr:
                    verts[-1], costs[-1] = xe, fe
                else:
                    verts[-1], costs[-1] = xr, fr
                continue

            # Contraction
            if fr < costs[-1]:
                xc = xbar + gamma * (xr - xbar)
                fc = _eval(xc, verts)
                if fc <= fr:
                    verts[-1], costs[-1] = xc, fc
                    continue
            else:
                xc = xbar + gamma * (verts[-1] - xbar)
                fc = _eval(xc, verts)
                if fc < costs[-1]:
                    verts[-1], costs[-1] = xc, fc
                    continue

            # Shrink
            new_verts = np.empty_like(verts)
            new_costs = np.empty_like(costs)
            new_verts[0], new_costs[0] = verts[0], costs[0]
            for i in range(1, len(verts)):
                vs = verts[0] + shrink * (verts[i] - verts[0])
                new_costs[i] = _eval(vs, verts, slot=i)
                new_verts[i] = vs
            verts, costs = new_verts, new_costs

    except StopIteration as e:
        status = str(e)

    if hasattr(evaluator, "current_simplex"):
        evaluator.current_simplex = None

    return _summary(history, status)


# ---------------------------------------------------------------------------
# CMA-ES backend (hand-coded, pure NumPy — no `cma` package needed)
# ---------------------------------------------------------------------------
#
# Reference implementation of (μ/μ_w, λ)-CMA-ES from
#     Hansen, "The CMA Evolution Strategy: A Tutorial" (2016), arXiv:1604.00772
# with the standard rank-μ + rank-1 covariance update and cumulative step-size
# adaptation (CSA).  No active update, no restarts (TPA / IPOP), no boundary
# handling — those are easy to add later but unnecessary for a 2-D unbounded
# problem with ≤ a few hundred evals.
#
# Notation matches the tutorial:
#     m         : mean of the search distribution                (n,)
#     sigma     : global step size                               scalar
#     C         : covariance matrix                              (n, n)
#     p_sigma   : evolution path for sigma (CSA)                 (n,)
#     p_c       : evolution path for C (rank-1 update)           (n,)
#     B, D      : eigendecomposition C = B diag(D)^2 B.T         (n,n), (n,)


def _cma_default_params(n: int, popsize: int):
    """Standard CMA-ES strategy parameters (Hansen 2016, Table 1)."""
    import math as _math
    lam = int(popsize)
    mu = lam // 2

    # Recombination weights w_i ∝ log((lam+1)/2) − log(i),  i = 1..mu
    raw_w = np.array([_math.log((lam + 1) / 2.0) - _math.log(i + 1)
                      for i in range(mu)])
    weights = raw_w / raw_w.sum()
    mu_eff  = 1.0 / np.sum(weights ** 2)        # variance-effective mu

    # Adaptation rates
    c_sigma = (mu_eff + 2.0) / (n + mu_eff + 5.0)
    d_sigma = 1.0 + 2.0 * max(0.0, _math.sqrt((mu_eff - 1.0) / (n + 1.0)) - 1.0) + c_sigma
    c_c     = (4.0 + mu_eff / n) / (n + 4.0 + 2.0 * mu_eff / n)
    c_1     = 2.0 / ((n + 1.3) ** 2 + mu_eff)
    c_mu    = min(1.0 - c_1,
                  2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((n + 2.0) ** 2 + mu_eff))

    # E[||N(0,I)||] (Hansen eq. 31)
    chiN = _math.sqrt(n) * (1.0 - 1.0 / (4.0 * n) + 1.0 / (21.0 * n * n))

    return {
        "lam": lam, "mu": mu, "weights": weights, "mu_eff": mu_eff,
        "c_sigma": c_sigma, "d_sigma": d_sigma,
        "c_c": c_c, "c_1": c_1, "c_mu": c_mu, "chiN": chiN,
    }


def run_cmaes(evaluator, x0=(1.0, 1.0), *,
              max_evals=40, max_gens=0, sigma0=0.1,
              popsize=6, tolx=0.005, tolfun=1e-6,
              snr_floor=0.0,
              indist_stop_snr=0.0,
              snr_target=0.0, snr_max_traj_factor=4,
              seed=None, pool=None):
    """Covariance Matrix Adaptation Evolution Strategy (hand-coded).

    Drop-in replacement for ``run_nelder_mead`` with the same evaluator
    interface and summary dict.  Pure-NumPy implementation of the
    (μ/μ_w, λ)-CMA-ES with rank-μ + rank-1 covariance updates and
    cumulative step-size adaptation; no external `cma` dependency.

    Run-budget parameters
    ---------------------
    ``popsize`` (λ)
        Number of evaluations per generation.  Default 6 (≈ 4 + 3·ln N for
        N=2).  Larger λ gives more accurate gradients per generation but
        costs more MC evals.
    ``max_gens``
        Hard cap on the number of generations.  ``0`` (default) means no
        per-generation cap — only ``max_evals`` limits the run.  When set,
        the effective evaluation budget is ``min(max_evals, max_gens *
        popsize)``.
    ``max_evals``
        Hard cap on the total number of evaluator calls.  Always honoured.

    Why CMA-ES on noisy / flat cost surfaces:
      • adapts a full n×n covariance + global step size σ each generation,
        so the search direction follows the descent ridge,
      • uses a population of λ evals per generation and weighted
        recombination of the μ best, so isolated noisy evaluations don't
        derail the mean,
      • CSA naturally shrinks σ near the optimum and grows it on plateaus.

    Stopping rules (any one triggers):
      • ``max_evals`` reached
      • ``max_gens`` generations completed (if > 0)
      • σ * sqrt(max eigenvalue of C) < ``tolx``         (search step tiny)
      • spread of last λ costs < ``tolfun``              (cost flat)
      • SNR-based rules from the evaluator (same as NM backend)
    """
    history: list = []
    _n_traj_prod_base = int(evaluator.n_traj_prod)
    _n_traj_max = _n_traj_prod_base * max(1, int(snr_max_traj_factor))

    rng = np.random.default_rng(seed)

    def _eval(x):
        if len(history) >= max_evals:
            raise StopIteration("max_evals reached")
        # CMA-ES has no simplex; clear any leftover from a previous NM run
        # so the visualizer doesn't draw a stale triangle.
        if hasattr(evaluator, "current_simplex"):
            evaluator.current_simplex = None
        # current_gaussian is set by the outer loop before each generation;
        # the evaluator forwards it to the optimizer plotter.
        res = evaluator(float(x[0]), float(x[1]))
        history.append(res)

        if snr_target > 0 and res.snr < snr_target:
            new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
            if new_n > evaluator.n_traj_prod:
                print(f"[cma] SNR={res.snr:.2f} < {snr_target}; "
                      f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n

        best = min(history, key=lambda r: r.cost)

        if snr_floor > 0 and best.snr < snr_floor and len(history) >= 5:
            print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}; "
                  "need more statistics.")
            raise StopIteration("snr floor reached")

        if indist_stop_snr > 0 and len(history) >= 5 and best.snr < indist_stop_snr:
            if snr_target > 0 and evaluator.n_traj_prod < _n_traj_max:
                new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
                print(f"[cma] best SNR={best.snr:.2f} < {indist_stop_snr} "
                      f"but stats not maxed; boosting "
                      f"n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n
            else:
                print(f"[stop] running best SNR={best.snr:.2f} < {indist_stop_snr}; "
                      f"test indistinguishable from reference at "
                      f"r1={best.r1:.4f} r2={best.r2:.4f}.")
                raise StopIteration("indistinguishable from reference")
        return res.cost

    # ----- initialise state -----
    n     = len(x0)
    m     = np.array(x0, dtype=float)
    sigma = float(sigma0)
    C     = np.eye(n)
    p_sigma = np.zeros(n)
    p_c     = np.zeros(n)

    p = _cma_default_params(n, popsize)
    lam     = p["lam"];  mu      = p["mu"]
    weights = p["weights"]
    mu_eff  = p["mu_eff"]
    c_sigma = p["c_sigma"];  d_sigma = p["d_sigma"]
    c_c     = p["c_c"];      c_1     = p["c_1"];  c_mu = p["c_mu"]
    chiN    = p["chiN"]

    # Initial eigendecomposition (C = I → trivial).
    B = np.eye(n)
    D = np.ones(n)
    BD = B * D                                    # n x n: cols are B[:,i]*D[i]

    status = "completed"
    gen    = 0
    last_pop_costs = None

    try:
        while len(history) < max_evals:
            gen += 1

            # Publish current search distribution so the optimizer plotter
            # can draw the Gaussian (1σ / 2σ ellipses) for this generation.
            if hasattr(evaluator, "current_gaussian"):
                evaluator.current_gaussian = {
                    "mean":  [float(m[0]), float(m[1])],
                    "cov":   [[float(C[0, 0]), float(C[0, 1])],
                              [float(C[1, 0]), float(C[1, 1])]],
                    "sigma": float(sigma),
                    "gen":   int(gen),
                }

            # ----- sample λ offspring  x_k = m + σ B D z_k,  z_k ~ N(0, I) -----
            zs = rng.standard_normal((lam, n))    # (lam, n)
            ys = zs @ BD.T                        # (lam, n) = z @ (BD).T
            xs = m + sigma * ys                   # (lam, n)

            # Evaluate (may raise StopIteration mid-generation).
            # When `pool` is provided (Speedup 3 / parallel CMA-ES) the
            # whole generation is dispatched in one batch and we replay
            # the per-eval housekeeping (history append, SNR boost,
            # stop-rule checks) on the returned EvalResult-shaped dicts
            # in submission order.  When `pool` is None we keep the
            # serial loop intact so legacy behaviour is bit-equivalent.
            fs = []
            if pool is None:
                for k in range(lam):
                    if len(history) >= max_evals:
                        break
                    fs.append(_eval(xs[k]))
            else:
                # Cap the batch at the remaining eval budget.
                budget = max(0, max_evals - len(history))
                batch = [(float(xs[k][0]), float(xs[k][1]))
                         for k in range(min(lam, budget))]
                if not batch:
                    break
                if hasattr(evaluator, "current_simplex"):
                    evaluator.current_simplex = None
                # Pool returns plain dicts (asdict(EvalResult)).  Replay
                # them through the same housekeeping that the serial
                # _eval does — adaptive stats, SNR boost, stop checks.
                eid_base = len(history) + 1
                results = pool.map_generation(batch, eid_base)
                for res_dict in results:
                    res = _ResultProxy(res_dict)
                    history.append(res)
                    fs.append(res.cost)
                    # Adaptive stats
                    if snr_target > 0 and res.snr < snr_target:
                        new_n = min(int(evaluator.n_traj_prod * 1.5),
                                    _n_traj_max)
                        if new_n > evaluator.n_traj_prod:
                            print(f"[cma|par] SNR={res.snr:.2f} < {snr_target}; "
                                  f"boosting n_traj_prod "
                                  f"{evaluator.n_traj_prod}→{new_n}")
                            evaluator.n_traj_prod = new_n
                    best = min(history, key=lambda r: r.cost)
                    if (snr_floor > 0 and best.snr < snr_floor
                            and len(history) >= 5):
                        print(f"[stop] running best SNR={best.snr:.2f} "
                              f"< {snr_floor}; need more statistics.")
                        raise StopIteration("snr floor reached")
                    if (indist_stop_snr > 0 and len(history) >= 5
                            and best.snr < indist_stop_snr):
                        if (snr_target > 0
                                and evaluator.n_traj_prod < _n_traj_max):
                            new_n = min(int(evaluator.n_traj_prod * 1.5),
                                        _n_traj_max)
                            evaluator.n_traj_prod = new_n
                        else:
                            print(f"[stop] running best SNR={best.snr:.2f} "
                                  f"< {indist_stop_snr}; "
                                  "test indistinguishable from reference at "
                                  f"r1={best.r1:.4f} r2={best.r2:.4f}.")
                            raise StopIteration(
                                "indistinguishable from reference")
            fs = np.array(fs, dtype=float)

            if len(fs) < mu:
                # Not enough completed evals to do a recombination — bail out.
                if len(fs) == 0:
                    break
                # Fall through with whatever we have; mu_used = len(fs)//2 + 1
                # so we still update something meaningful.
                pass

            # ----- selection: take the μ best of the actually-evaluated -----
            n_eval = len(fs)
            mu_used = max(1, min(mu, n_eval // 2 + (n_eval % 2)))
            # If we got a partial population, renormalise the weights.
            if mu_used != mu:
                w = weights[:mu_used]
                w = w / w.sum()
                mueff_used = 1.0 / np.sum(w ** 2)
            else:
                w = weights
                mueff_used = mu_eff

            order = np.argsort(fs)
            sel   = order[:mu_used]
            ys_sel = ys[sel]                      # (mu_used, n)
            zs_sel = zs[sel]                      # (mu_used, n)

            last_pop_costs = fs.copy()

            # ----- recombination -----
            y_w = w @ ys_sel                      # (n,)  weighted mean of y
            z_w = w @ zs_sel
            m_old = m
            m = m_old + sigma * y_w               # update mean

            # ----- step-size path & σ update (CSA) -----
            # C^{-1/2} (m - m_old) / σ  =  B z_w
            p_sigma = (1.0 - c_sigma) * p_sigma \
                      + np.sqrt(c_sigma * (2.0 - c_sigma) * mueff_used) * (B @ z_w)
            sigma = sigma * np.exp((c_sigma / d_sigma) *
                                   (np.linalg.norm(p_sigma) / chiN - 1.0))

            # ----- covariance path & C update -----
            # h_sigma: indicator that ||p_sigma|| is "not too large" — prevents
            # spurious axis growth in early generations.  (Hansen eq. 33.)
            hs_thresh = (1.4 + 2.0 / (n + 1.0)) * chiN
            denom = np.sqrt(1.0 - (1.0 - c_sigma) ** (2 * gen))
            h_sigma = 1.0 if (np.linalg.norm(p_sigma) / max(denom, 1e-30)) < hs_thresh else 0.0

            p_c = (1.0 - c_c) * p_c \
                  + h_sigma * np.sqrt(c_c * (2.0 - c_c) * mueff_used) * y_w

            # rank-μ update term: Σ w_i y_i y_i^T
            artmp = ys_sel                        # (mu_used, n)
            rank_mu = artmp.T @ (w[:, None] * artmp)
            # rank-1 update term: p_c p_c^T  (with h_sigma correction)
            rank_1  = np.outer(p_c, p_c) \
                      + (1.0 - h_sigma) * c_c * (2.0 - c_c) * C   # variance loss correction

            C = (1.0 - c_1 - c_mu) * C \
                + c_1 * rank_1 \
                + c_mu * rank_mu

            # Keep C symmetric (numerical hygiene) and re-decompose.
            C = 0.5 * (C + C.T)
            try:
                eigvals, B = np.linalg.eigh(C)
            except np.linalg.LinAlgError:
                # Reset on numerical failure.
                C = np.eye(n)
                eigvals, B = np.linalg.eigh(C)
                p_sigma = np.zeros(n); p_c = np.zeros(n)
            eigvals = np.maximum(eigvals, 1e-30)  # guard against tiny negatives
            D = np.sqrt(eigvals)
            BD = B * D

            # ----- termination tests -----
            if max_gens and gen >= int(max_gens):
                status = f"max_gens reached ({gen})"
                break
            sigma_step = sigma * float(D.max())
            if sigma_step < tolx:
                status = f"converged (tolx: σ·max(D)={sigma_step:.2e} < {tolx})"
                break
            if last_pop_costs is not None and len(last_pop_costs) >= 2:
                spread = float(last_pop_costs.max() - last_pop_costs.min())
                if spread < tolfun:
                    status = f"converged (tolfun: pop spread={spread:.2e} < {tolfun})"
                    break

    except StopIteration as e:
        status = str(e)

    if hasattr(evaluator, "current_simplex"):
        evaluator.current_simplex = None
    if hasattr(evaluator, "current_gaussian"):
        evaluator.current_gaussian = None

    return _summary(history, status)
