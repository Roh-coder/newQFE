"""
optimizers.py — five classical optimizer drivers for the (r1, r2) search.

All drivers share the same signature:

    run_<method>(evaluator, x0, *, max_evals, sigma0=...) -> dict

where *evaluator* is an Evaluator instance (callable returning EvalResult).
Each driver calls evaluator(r1, r2) repeatedly and returns a summary dict
with the converged best point, total evaluations, and the method name.

Stopping rules (applied uniformly):
  - Hard cap on the number of evaluations (`max_evals`).
  - Soft "statistics adequacy" check: if the running best has SNR < 1
    (i.e. the cost is already consistent with zero), stop early — the
    optimizer cannot make further progress without more MC statistics.
"""

from __future__ import annotations

from typing import Callable
import numpy as np

# Optional dependencies; imported lazily inside their drivers.
# scikit-optimize for GP/Bayesian, cma for CMA-ES.


# ---------------------------------------------------------------------------
#  Shared helpers
# ---------------------------------------------------------------------------

def _wrap(evaluator, history: list, snr_floor: float, max_evals: int):
    """Return a scalar cost function f(x) that records EvalResults.

    Raises ``StopIteration`` to break the optimizer loop when the budget or
    the SNR floor is reached.
    """
    def f(x):
        if len(history) >= max_evals:
            raise StopIteration("max_evals reached")
        r1, r2 = float(x[0]), float(x[1])
        res = evaluator(r1, r2)
        history.append(res)
        # SNR floor check on the running best.
        best = min(history, key=lambda r: r.cost)
        if best.snr < snr_floor and len(history) >= 5:
            print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}; "
                  f"need more statistics.")
            raise StopIteration("snr floor reached")
        return res.cost
    return f


def _summary(method: str, history: list, status: str) -> dict:
    if not history:
        return {"method": method, "status": status, "n_evals": 0}
    best = min(history, key=lambda r: r.cost)
    return {
        "method": method,
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


# ---------------------------------------------------------------------------
#  1. Nelder-Mead simplex  (hand-coded; no scipy dependency)
# ---------------------------------------------------------------------------

def run_nelder_mead(evaluator, x0=(1.0, 1.0), *,
                    max_evals=40, sigma0=0.1, snr_floor=1.0,
                    xatol=0.005, fatol=1e-6,
                    alpha=1.0, beta_exp=2.0, gamma=0.5, shrink=0.75,
                    noise_stop_snr=0.0, noise_stop_n=5,
                    indist_stop_snr=0.0,
                    snr_target=0.0, snr_max_traj_factor=4):
    """
    Hand-coded Nelder-Mead downhill-simplex minimiser (2-D).

    Standard NM coefficients
    ------------------------
    alpha    = 1.0    reflection   — mirror worst vertex through centroid
    beta_exp = 2.0    expansion    — stretch further along reflection direction
    gamma    = 0.5    contraction  — pull toward centroid (inside or outside)
    shrink   = 0.75   shrink       — collapse all non-best vertices toward best

    scipy defaults to shrink=0.5, which halves the simplex on every failed
    contraction step.  Under MC noise, failed contractions are common (a
    reflected/contracted point can look better than it really is due to a
    lucky fluctuation), so shrink=0.5 collapses the triangle well before it
    has walked to the minimum.  shrink=0.75 keeps the simplex larger longer,
    giving the method more room to explore before it freezes.

    Noise-dominated early stop
    ---------------------------
    noise_stop_snr  = 0.0   threshold SNR; 0 disables the check entirely.
    noise_stop_n    = 5     stop if the last *noise_stop_n* consecutive evals
                            all have SNR < noise_stop_snr.  This signals that
                            the simplex has descended into a region where MC
                            noise completely swamps the cost signal and further
                            NM steps are meaningless.
    indist_stop_snr = 0.0   stop if the running best has SNR < this value
                            (i.e. cost < indist_stop_snr × sigma_cost).
                            SNR < 1 means the test correlator is statistically
                            indistinguishable from the reference — the
                            optimizer has found the target coupling.

    Adaptive MC statistics
    ----------------------
    snr_target          = 0.0   when > 0: after each eval, if SNR < snr_target
                                multiply n_traj_prod by 1.5 (capped at
                                snr_max_traj_factor × starting value).  This
                                automatically tightens the cost measurement
                                as the simplex nears the minimum where cost
                                values get small.
    snr_max_traj_factor = 4     maximum multiple of starting n_traj_prod.

    Simplex exposure for the visualizer
    ------------------------------------
    Before every evaluator call the current simplex vertices are stored on
    ``evaluator.current_simplex`` (list of three [r1, r2] lists) so that
    OptimizerPlotter can draw it as a filled triangle in the trajectory panel.
    """
    history: list = []
    _n_traj_prod_base = int(evaluator.n_traj_prod)
    _n_traj_max = _n_traj_prod_base * max(1, int(snr_max_traj_factor))

    def _eval(pt, simplex, slot=-1):
        """Evaluate pt=[r1,r2]; push current simplex to plotter; check limits.

        The simplex pushed to the visualizer has *pt* substituted at *slot*
        (default = worst-vertex index -1), so the displayed triangle reflects
        what is actually being tried right now — not the stale pre-iteration
        triangle.  Pass slot=None to display the simplex unchanged (used for
        the initial vertex evaluations where pt already IS a simplex vertex).
        """
        if len(history) >= max_evals:
            raise StopIteration("max_evals reached")
        # Store simplex on evaluator so optimizer_plot.update() picks it up.
        if hasattr(evaluator, 'current_simplex'):
            display_simplex = [np.asarray(v).tolist() for v in simplex]
            if slot is not None:
                display_simplex[slot] = [float(pt[0]), float(pt[1])]
            evaluator.current_simplex = display_simplex
        res = evaluator(float(pt[0]), float(pt[1]))
        history.append(res)
        # Adaptive stats: boost n_traj_prod when SNR is too low.
        if snr_target > 0 and res.snr < snr_target:
            new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
            if new_n > evaluator.n_traj_prod:
                print(f"[nm] SNR={res.snr:.2f} < {snr_target}; "
                      f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n
        best = min(history, key=lambda r: r.cost)
        if best.snr < snr_floor and len(history) >= 5:
            print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}; "
                  f"need more statistics.")
            raise StopIteration("snr floor reached")
        # Noise-dominated stop: last noise_stop_n evals all below threshold.
        if noise_stop_snr > 0 and len(history) >= noise_stop_n:
            recent = history[-noise_stop_n:]
            if all(r.snr < noise_stop_snr for r in recent):
                avg_snr = sum(r.snr for r in recent) / noise_stop_n
                print(f"[stop] last {noise_stop_n} evals all noise-dominated "
                      f"(avg SNR={avg_snr:.2f} < {noise_stop_snr}); "
                      f"curve collapse swamped by MC noise.")
                raise StopIteration("noise dominated")
        # Indistinguishable-from-reference stop: running best cost < noise.
        # Only fire when stats are already at maximum — if there is still room
        # to boost, boost instead so the best point can be re-confirmed later.
        if indist_stop_snr > 0 and len(history) >= 5:
            if best.snr < indist_stop_snr:
                if snr_target > 0 and evaluator.n_traj_prod < _n_traj_max:
                    new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
                    print(f"[nm] best SNR={best.snr:.2f} < {indist_stop_snr} but stats not maxed; "
                          f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                    evaluator.n_traj_prod = new_n
                else:
                    print(f"[stop] running best SNR={best.snr:.2f} < {indist_stop_snr}; "
                          f"test lattice indistinguishable from reference at "
                          f"r1={best.r1:.4f} r2={best.r2:.4f}.")
                    raise StopIteration("indistinguishable from reference")
        return res.cost

    # Initial right-triangle simplex centred on x0.
    s = float(sigma0)
    verts = np.array([
        [x0[0],       x0[1]    ],
        [x0[0] + s,   x0[1]    ],
        [x0[0],       x0[1] + s],
    ], dtype=float)
    costs = np.full(3, np.inf)

    status = "completed"
    try:
        # --- Evaluate the three initial vertices -------------------------
        for i in range(len(verts)):
            costs[i] = _eval(verts[i], verts, slot=None)

        while len(history) < max_evals:
            # Sort: verts[0] = best (lowest cost), verts[-1] = worst.
            order = np.argsort(costs)
            verts, costs = verts[order], costs[order]

            # Convergence: simplex diameter and cost spread both small.
            diam = max(np.linalg.norm(v - verts[0]) for v in verts[1:])
            if diam < xatol and (costs[-1] - costs[0]) < fatol:
                status = "converged"
                break

            # Centroid of all but worst vertex.
            xbar = verts[:-1].mean(axis=0)

            # ---- Reflection ----
            xr = xbar + alpha * (xbar - verts[-1])
            fr = _eval(xr, verts)

            if costs[0] <= fr < costs[-2]:
                # Better than second-worst but not best: accept reflection.
                verts[-1], costs[-1] = xr, fr
                continue

            if fr < costs[0]:
                # Better than best: try expansion before accepting.
                xe = xbar + beta_exp * (xr - xbar)
                fe = _eval(xe, verts)
                if fe < fr:
                    verts[-1], costs[-1] = xe, fe
                else:
                    verts[-1], costs[-1] = xr, fr
                continue

            # ---- Contraction (reflection was not useful) ----
            if fr < costs[-1]:
                # Outside contraction: reflected point was better than worst.
                xc = xbar + gamma * (xr - xbar)
                fc = _eval(xc, verts)
                if fc <= fr:
                    verts[-1], costs[-1] = xc, fc
                    continue
            else:
                # Inside contraction: reflected point was worse than worst.
                xc = xbar + gamma * (verts[-1] - xbar)
                fc = _eval(xc, verts)
                if fc < costs[-1]:
                    verts[-1], costs[-1] = xc, fc
                    continue

            # ---- Shrink: pull all non-best vertices toward best ---------
            # Show old (pre-shrink) simplex in the visualizer while each
            # shrunk vertex is being evaluated.
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

    # Clear the simplex reference so plotter doesn't keep drawing it.
    if hasattr(evaluator, 'current_simplex'):
        evaluator.current_simplex = None

    return _summary("nelder_mead", history, status)


# ---------------------------------------------------------------------------
#  2. Powell's method
# ---------------------------------------------------------------------------

def run_powell(evaluator, x0=(1.0, 1.0), *,
               max_evals=40, sigma0=0.1, snr_floor=1.0):
    from scipy.optimize import minimize
    history: list = []
    f = _wrap(evaluator, history, snr_floor, max_evals)
    direc = np.array([[sigma0, 0.0], [0.0, sigma0]])
    status = "completed"
    try:
        minimize(f, np.asarray(x0, dtype=float),
                 method="Powell",
                 options={"direc": direc,
                          "xtol": 0.005, "ftol": 1e-6,
                          "maxiter": max_evals, "maxfev": max_evals,
                          "disp": False})
    except StopIteration as e:
        status = str(e)
    return _summary("powell", history, status)


# ---------------------------------------------------------------------------
#  3. BFGS with finite-difference gradients (3-point)
# ---------------------------------------------------------------------------

def run_bfgs_fd(evaluator, x0=(1.0, 1.0), *,
                max_evals=40, sigma0=0.05, snr_floor=1.0):
    from scipy.optimize import minimize
    history: list = []
    f = _wrap(evaluator, history, snr_floor, max_evals)
    status = "completed"
    try:
        minimize(f, np.asarray(x0, dtype=float),
                 method="L-BFGS-B",
                 jac="3-point",
                 options={"eps": sigma0,
                          "maxiter": max_evals, "maxfun": max_evals,
                          "ftol": 1e-7, "gtol": 1e-5,
                          "disp": False})
    except StopIteration as e:
        status = str(e)
    return _summary("bfgs_fd", history, status)


# ---------------------------------------------------------------------------
#  4. Gaussian-process Bayesian optimization (scikit-optimize)
# ---------------------------------------------------------------------------

def run_gp(evaluator, x0=(1.0, 1.0), *,
           max_evals=30, bounds=((0.5, 1.5), (0.5, 1.5)),
           snr_floor=1.0, n_initial=8):
    try:
        from skopt import Optimizer as SkOptimizer
    except ImportError:
        return {"method": "gp", "status": "skipped: scikit-optimize not installed",
                "n_evals": 0}
    history: list = []
    opt = SkOptimizer(dimensions=list(bounds),
                      base_estimator="GP",
                      acq_func="EI",
                      n_initial_points=n_initial,
                      random_state=42)
    # Seed with x0 first.
    seed_pts = [list(x0)]
    status = "completed"
    try:
        for x_seed in seed_pts:
            res = evaluator(*x_seed)
            history.append(res)
            opt.tell(x_seed, res.cost)
        while len(history) < max_evals:
            x = opt.ask()
            res = evaluator(float(x[0]), float(x[1]))
            history.append(res)
            opt.tell(x, res.cost)
            best = min(history, key=lambda r: r.cost)
            if best.snr < snr_floor and len(history) >= 5:
                print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}.")
                status = "snr floor reached"
                break
    except StopIteration as e:
        status = str(e)
    return _summary("gp", history, status)


# ---------------------------------------------------------------------------
#  5. CMA-ES (covariance matrix adaptation)
# ---------------------------------------------------------------------------

def run_cma(evaluator, x0=(1.0, 1.0), *,
            max_evals=40, sigma0=0.1, snr_floor=1.0,
            bounds=((0.5, 1.5), (0.5, 1.5))):
    try:
        import cma
    except ImportError:
        return {"method": "cma", "status": "skipped: cma not installed",
                "n_evals": 0}
    history: list = []
    es = cma.CMAEvolutionStrategy(
        list(x0), sigma0,
        {"bounds": [[bounds[0][0], bounds[1][0]],
                    [bounds[0][1], bounds[1][1]]],
         "popsize": 6, "maxfevals": max_evals,
         "verbose": -9, "tolx": 0.005, "tolfun": 1e-6, "seed": 42})
    status = "completed"
    try:
        while not es.stop() and len(history) < max_evals:
            xs = es.ask()
            costs = []
            for x in xs:
                if len(history) >= max_evals:
                    break
                res = evaluator(float(x[0]), float(x[1]))
                history.append(res)
                costs.append(res.cost)
                best = min(history, key=lambda r: r.cost)
                if best.snr < snr_floor and len(history) >= 5:
                    print(f"[stop] running best SNR={best.snr:.2f}.")
                    status = "snr floor reached"
                    raise StopIteration
            if len(costs) == len(xs):
                es.tell(xs, costs)
    except StopIteration:
        pass
    return _summary("cma", history, status)


# ---------------------------------------------------------------------------
#  Registry
# ---------------------------------------------------------------------------

def run_bfgs_hand(evaluator, x0=(1.0, 1.0), *,
                  max_evals=40, sigma0=0.1, snr_floor=1.0,
                  fd_eps=0.05, max_step=0.4,
                  noise_stop_snr=0.0, noise_stop_n=5,
                  indist_stop_snr=0.0,
                  snr_target=3.0, snr_max_traj_factor=4):
    """
    Hand-coded BFGS with central-difference FD gradient (2-D).

    Algorithm
    ---------
    Maintains a rank-2 BFGS approximation H of the inverse Hessian.
    At each outer iteration:

      1. Estimate gradient g via central FD (4 evals: ±fd_eps along each axis).
      2. Compute search direction  p = −H g   (pure Newton-like descent).
      3. Armijo backtracking line search along p, starting from the
         step that moves ‖step‖ = min(1, sigma0/‖p‖) in parameter space.
         At most ``max_backtracks`` halvings; accepts the last candidate
         regardless (MC noise can block the Armijo condition permanently).
      4. BFGS curvature update:
           s = x_new − x,   y = g_new − g
           If  s·y > 0  (positive curvature), perform the standard update:
               ρ = 1/(y·s)
               H ← (I − ρ s yᵀ) H (I − ρ y sᵀ) + ρ s sᵀ
           Otherwise reset H ← I (noise destroyed curvature info).
      5. Move to x_new; reuse g_new as gradient for next iteration.

    Noise robustness
    ----------------
    fd_eps              : central-FD step.  Should be ≳ MC noise scale on the
                          cost (~σ_cost).  0.10 keeps gradient signal well
                          above the ~1e-4 cost noise at far-start points.
    max_step            : hard cap on ‖step‖ to prevent runaway when H is
                          poorly conditioned early in the run.
    snr_target          : adaptive stats threshold.  After each eval, if
                          SNR < snr_target the evaluator's n_traj_prod is
                          multiplied by 1.5 (capped at snr_max_traj_factor ×
                          starting value).  This automatically tightens the
                          cost measurement when the gradient signal is weak.
    snr_max_traj_factor : maximum multiple of starting n_traj_prod allowed
                          by the adaptive boost (default 4).
    noise_stop_snr / noise_stop_n / indist_stop_snr : same semantics as
                          run_nelder_mead.

    Eval budget breakdown (approx)
    --------------------------------
    Initial eval   :  1
    Initial grad   :  4
    Per BFGS step  :  1–4 (line search)  +  4 (new gradient)  ≈  8
    → budget 25  gives ~2–3 full BFGS steps after initial setup.
    → budget 60  gives ~6 full BFGS steps.
    """
    history: list = []
    _n_traj_prod_base = int(evaluator.n_traj_prod)
    _n_traj_max = _n_traj_prod_base * max(1, int(snr_max_traj_factor))

    def _eval(pt):
        if len(history) >= max_evals:
            raise StopIteration("max_evals reached")
        res = evaluator(float(pt[0]), float(pt[1]))
        # Adaptive stats: if SNR too low, boost n_traj_prod for next evals
        # (up to snr_max_traj_factor × base) so the cost signal stays useful.
        if snr_target > 0 and res.snr < snr_target:
            new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
            if new_n > evaluator.n_traj_prod:
                print(f"[bfgs] SNR={res.snr:.2f} < {snr_target}; "
                      f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                evaluator.n_traj_prod = new_n
        history.append(res)
        best = min(history, key=lambda r: r.cost)
        if best.snr < snr_floor and len(history) >= 5:
            print(f"[stop] running best SNR={best.snr:.2f} < {snr_floor}; "
                  "need more statistics.")
            raise StopIteration("snr floor reached")
        if noise_stop_snr > 0 and len(history) >= noise_stop_n:
            recent = history[-noise_stop_n:]
            if all(r.snr < noise_stop_snr for r in recent):
                avg = sum(r.snr for r in recent) / noise_stop_n
                print(f"[stop] last {noise_stop_n} evals noise-dominated "
                      f"(avg SNR={avg:.2f} < {noise_stop_snr}).")
                raise StopIteration("noise dominated")
        if indist_stop_snr > 0 and len(history) >= 5:
            if best.snr < indist_stop_snr:
                if snr_target > 0 and evaluator.n_traj_prod < _n_traj_max:
                    new_n = min(int(evaluator.n_traj_prod * 1.5), _n_traj_max)
                    print(f"[bfgs] best SNR={best.snr:.2f} < {indist_stop_snr} but stats not maxed; "
                          f"boosting n_traj_prod {evaluator.n_traj_prod}→{new_n}")
                    evaluator.n_traj_prod = new_n
                else:
                    print(f"[stop] running best SNR={best.snr:.2f} < {indist_stop_snr}; "
                          f"indistinguishable from reference at "
                          f"r1={best.r1:.4f} r2={best.r2:.4f}.")
                    raise StopIteration("indistinguishable from reference")
        return res.cost

    def _fd_grad(pt):
        """Central-difference gradient; 4 evaluations."""
        g = np.zeros(2)
        for i in range(2):
            ei = np.zeros(2); ei[i] = fd_eps
            fp = _eval(pt + ei)
            fm = _eval(pt - ei)
            g[i] = (fp - fm) / (2.0 * fd_eps)
        return g

    I2 = np.eye(2)
    H = I2.copy()           # approximate inverse Hessian
    x = np.asarray(x0, dtype=float).copy()
    status = "completed"

    try:
        fx = _eval(x)
        g = _fd_grad(x)     # 4 evals

        while len(history) < max_evals:
            # --- Search direction ------------------------------------------
            p = -H @ g
            gp = float(g @ p)
            if gp >= 0.0:
                # Not a descent direction (can happen when H is bad); reset.
                H = I2.copy()
                p = -g
                gp = float(g @ p)

            pnorm = np.linalg.norm(p)
            if pnorm < 1e-12:
                status = "zero_gradient"
                break

            # Scale initial step so ‖α p‖ = min(1, sigma0) in param space.
            alpha = min(1.0, sigma0) / pnorm

            # --- Armijo backtracking line search ---------------------------
            c1 = 1e-4
            max_backtracks = 4
            x_new = np.clip(x + alpha * p, 0.3, 2.5)
            fx_new = _eval(x_new)
            for _ in range(max_backtracks):
                if fx_new <= fx + c1 * alpha * gp:
                    break
                alpha *= 0.5
                x_new = np.clip(x + alpha * p, 0.3, 2.5)
                fx_new = _eval(x_new)
            # Accept regardless; MC noise can permanently block Armijo.

            # Hard cap on actual step length.
            step = x_new - x
            if np.linalg.norm(step) > max_step:
                step = step * (max_step / np.linalg.norm(step))
                x_new = x + step
                fx_new = _eval(x_new)

            # --- Gradient at new point ------------------------------------
            g_new = _fd_grad(x_new)   # 4 evals

            # --- BFGS inverse Hessian update ------------------------------
            s = x_new - x
            y = g_new - g
            sy = float(s @ y)
            sy_thresh = 1e-10 * np.linalg.norm(s) * np.linalg.norm(y)
            if sy > sy_thresh:
                rho = 1.0 / sy
                A = I2 - rho * np.outer(s, y)
                H = A @ H @ A.T + rho * np.outer(s, s)
            else:
                # Curvature condition violated (noise); reset.
                print(f"[bfgs] curvature reset at eval {len(history)} "
                      f"(s·y={sy:.2e})")
                H = I2.copy()

            x, fx, g = x_new, fx_new, g_new

            # Noise-limited convergence: gradient norm below FD noise floor.
            gnorm = np.linalg.norm(g)
            sigma_g = (history[-1].sigma_cost if history else 0.0) / fd_eps
            if gnorm < max(1e-5, sigma_g):
                status = "converged"
                break

    except StopIteration as e:
        status = str(e)

    return _summary("bfgs_hand", history, status)


# ---------------------------------------------------------------------------
#  Registry
# ---------------------------------------------------------------------------

ALL_METHODS: dict = {
    "nelder_mead": run_nelder_mead,
    "powell":      run_powell,
    "bfgs_fd":     run_bfgs_fd,
    "bfgs_hand":   run_bfgs_hand,
    "gp":          run_gp,
    "cma":         run_cma,
}
