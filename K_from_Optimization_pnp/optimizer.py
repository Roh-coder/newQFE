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
