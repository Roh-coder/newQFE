"""
surrogate.py — Cost-surface surrogate (Speedup 5).

A pure-NumPy heteroscedastic Gaussian-process regressor for the noisy
2-D cost surface ``c(r1, r2)`` produced by the Monte Carlo pipeline.
The surrogate is the engine behind the ``"bo"`` optimiser backend: at
each step we

    1. fit / refit the GP on the existing ``(r1, r2, cost, sigma_cost)``
       evaluations,
    2. propose the next ``(r1, r2)`` by minimising an acquisition
       function (default: lower-confidence bound (LCB)) over a coarse
       grid + a local Nelder-Mead polish,
    3. dispatch the proposed point to the real evaluator,
    4. append the result and loop.

The GP uses a Matérn-5/2 kernel — the standard choice for noisy
optimisation.  Its smoothness assumption (twice mean-square
differentiable) is a reasonable fit for the curve-collapse cost
surface, which is locally quadratic near the optimum but has noisier
flanks.

Hyper-parameters (signal variance σ_f², length scales ℓ_d) are fit by
maximising the log marginal likelihood with `scipy.optimize.minimize`
using L-BFGS-B and 5 random restarts.  The per-point measurement
noise σ_i is taken from ``sigma_cost`` directly (heteroscedastic),
floored at ``noise_floor`` to keep the kernel matrix well-conditioned.

This module is independent of the evaluator: it accepts arrays and
returns arrays.  The ``"bo"`` backend in :mod:`optimizer` glues it to
the production loop.

Public API
----------
* ``MaternGP`` — heteroscedastic Matérn-5/2 GP regressor with `fit`,
  `predict`, `log_marginal_likelihood`, `sample_posterior`.
* ``lcb`` — lower-confidence-bound acquisition (κ-weighted).
* ``ei``  — expected-improvement acquisition (rarely needed; LCB is
  the workhorse for noisy opt).
* ``propose_next`` — search the acquisition function over a domain
  rectangle + a local polish; returns ``(r1, r2)``.

Calibration
-----------
``MaternGP.calibration_residuals(X_held, y_held, sigma_held)``
returns standardised residuals ``(y - μ) / sqrt(σ_pred² + σ_obs²)``
on a held-out set; the BO panel acceptance check requires the 95th
percentile to sit within 2.0 (i.e. error bars are honest).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Tuple

import numpy as np
from scipy.optimize import minimize


# ---------------------------------------------------------------------------
# Matérn-5/2 kernel
# ---------------------------------------------------------------------------

def _scaled_dist(X1: np.ndarray, X2: np.ndarray, length_scales: np.ndarray) -> np.ndarray:
    """Pairwise scaled Euclidean distances ‖(x_i − y_j)/ℓ‖₂."""
    diff = (X1[:, None, :] - X2[None, :, :]) / length_scales
    return np.sqrt(np.sum(diff * diff, axis=-1) + 1e-300)


def matern52(X1: np.ndarray, X2: np.ndarray, sigma_f: float,
             length_scales: np.ndarray) -> np.ndarray:
    """k(x, y) = σ_f² · (1 + √5 r + 5/3 r²) · exp(−√5 r), r = ‖(x−y)/ℓ‖."""
    r = _scaled_dist(X1, X2, length_scales)
    sqrt5_r = np.sqrt(5.0) * r
    return (sigma_f * sigma_f) * (1.0 + sqrt5_r + (5.0 / 3.0) * r * r) * np.exp(-sqrt5_r)


# ---------------------------------------------------------------------------
# Heteroscedastic GP regressor
# ---------------------------------------------------------------------------

@dataclass
class _GPState:
    sigma_f: float
    length_scales: np.ndarray
    X: np.ndarray
    y: np.ndarray
    sigma_obs: np.ndarray
    L: np.ndarray         # Cholesky of K + diag(σ²) + jitter
    alpha: np.ndarray     # K⁻¹ (y − ȳ)
    y_mean: float


class MaternGP:
    """Heteroscedastic Matérn-5/2 Gaussian-process regressor (n_dim=2 default).

    Parameters
    ----------
    n_dim : int
        Input dimensionality (use 2 for ``(r1, r2)``).
    noise_floor : float
        Lower bound applied to the per-point ``sigma_obs`` before it goes
        on the diagonal.  Prevents singular kernels when an MC eval
        accidentally reports ``σ ≈ 0``.
    jitter : float
        Diagonal jitter added on top of σ² for numerical stability.
    """

    def __init__(self, n_dim: int = 2, noise_floor: float = 1e-6,
                 jitter: float = 1e-8):
        self.n_dim = int(n_dim)
        self.noise_floor = float(noise_floor)
        self.jitter = float(jitter)
        self._state: Optional[_GPState] = None

    # ------------------------------------------------------------------
    # Internal: build Cholesky + α for given hyper-parameters
    # ------------------------------------------------------------------

    def _build(self, X, y, sigma_obs, sigma_f, length_scales):
        K = matern52(X, X, sigma_f, length_scales)
        diag = np.maximum(sigma_obs, self.noise_floor) ** 2 + self.jitter
        K = K + np.diag(diag)
        try:
            L = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            # Bump jitter until decomposition succeeds.
            for boost in (1e-6, 1e-4, 1e-2):
                try:
                    L = np.linalg.cholesky(K + boost * np.eye(len(X)))
                    break
                except np.linalg.LinAlgError:
                    continue
            else:
                raise
        y_mean = float(np.mean(y))
        z = y - y_mean
        alpha = np.linalg.solve(L.T, np.linalg.solve(L, z))
        return L, alpha, y_mean

    def _nll(self, theta, X, y, sigma_obs):
        sigma_f = float(np.exp(theta[0]))
        ls = np.exp(theta[1:1 + self.n_dim])
        try:
            L, alpha, y_mean = self._build(X, y, sigma_obs, sigma_f, ls)
        except np.linalg.LinAlgError:
            return 1e12
        z = y - y_mean
        n = len(y)
        nll = 0.5 * z @ alpha + np.sum(np.log(np.diag(L))) + 0.5 * n * np.log(2.0 * np.pi)
        return float(nll)

    # ------------------------------------------------------------------
    # Public: fit / predict
    # ------------------------------------------------------------------

    def fit(self, X: np.ndarray, y: np.ndarray, sigma_obs: np.ndarray,
            *, n_restarts: int = 5, seed: int = 0) -> "MaternGP":
        """Fit hyper-parameters by L-BFGS-B with random restarts.

        Returns ``self`` so calls can be chained.
        """
        X = np.asarray(X, dtype=float).reshape(-1, self.n_dim)
        y = np.asarray(y, dtype=float).reshape(-1)
        sigma_obs = np.asarray(sigma_obs, dtype=float).reshape(-1)
        if len(X) < 2:
            raise ValueError("MaternGP.fit needs ≥ 2 points")
        rng = np.random.default_rng(seed)
        # Reasonable scale priors: log σ_f around log std(y); log ℓ around
        # log of input range / 4.
        y_scale = max(float(np.std(y)), 1e-6)
        x_ranges = np.maximum(X.max(axis=0) - X.min(axis=0), 1e-3)
        log_sigma_f0 = np.log(y_scale)
        log_ls0 = np.log(x_ranges / 4.0)
        bounds = [(np.log(1e-4 * y_scale), np.log(1e2 * y_scale))]
        bounds += [(np.log(1e-3 * r), np.log(10.0 * r)) for r in x_ranges]
        best = (np.inf, None)
        for k in range(max(1, int(n_restarts))):
            if k == 0:
                theta0 = np.r_[log_sigma_f0, log_ls0]
            else:
                theta0 = np.r_[
                    log_sigma_f0 + 0.5 * rng.standard_normal(),
                    log_ls0 + 0.5 * rng.standard_normal(self.n_dim),
                ]
            theta0 = np.clip(theta0, [b[0] for b in bounds],
                             [b[1] for b in bounds])
            try:
                res = minimize(self._nll, theta0, args=(X, y, sigma_obs),
                               method="L-BFGS-B", bounds=bounds,
                               options={"maxiter": 100, "ftol": 1e-6})
            except Exception:
                continue
            if res.fun < best[0]:
                best = (float(res.fun), res.x)
        if best[1] is None:
            raise RuntimeError("MaternGP.fit: all restarts failed")
        theta = best[1]
        sigma_f = float(np.exp(theta[0]))
        length_scales = np.exp(theta[1:1 + self.n_dim])
        L, alpha, y_mean = self._build(X, y, sigma_obs, sigma_f, length_scales)
        self._state = _GPState(
            sigma_f=sigma_f,
            length_scales=length_scales,
            X=X.copy(), y=y.copy(), sigma_obs=sigma_obs.copy(),
            L=L, alpha=alpha, y_mean=y_mean,
        )
        return self

    def predict(self, X_star: np.ndarray, *,
                return_std: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """Posterior mean (and std if ``return_std``) at ``X_star``."""
        if self._state is None:
            raise RuntimeError("call fit() before predict()")
        s = self._state
        Xs = np.asarray(X_star, dtype=float).reshape(-1, self.n_dim)
        K_s = matern52(s.X, Xs, s.sigma_f, s.length_scales)  # (N, M)
        mu = s.y_mean + K_s.T @ s.alpha
        if not return_std:
            return mu, None
        v = np.linalg.solve(s.L, K_s)                         # (N, M)
        K_ss_diag = (s.sigma_f * s.sigma_f) * np.ones(Xs.shape[0])
        var = K_ss_diag - np.sum(v * v, axis=0)
        var = np.maximum(var, 0.0)
        return mu, np.sqrt(var)

    def log_marginal_likelihood(self) -> float:
        if self._state is None:
            raise RuntimeError("call fit() first")
        s = self._state
        z = s.y - s.y_mean
        n = len(s.y)
        return float(-0.5 * z @ s.alpha
                     - np.sum(np.log(np.diag(s.L)))
                     - 0.5 * n * np.log(2.0 * np.pi))

    def calibration_residuals(self, X_held: np.ndarray, y_held: np.ndarray,
                              sigma_held: np.ndarray) -> np.ndarray:
        """Standardised residuals on held-out points.

        ``r_i = (y_i − μ_pred(x_i)) / sqrt(σ_pred(x_i)² + σ_obs(x_i)²)``.
        Honest error bars give ``r ~ N(0, 1)``.
        """
        mu, sd = self.predict(X_held, return_std=True)
        sigma_obs = np.maximum(np.asarray(sigma_held, dtype=float),
                               self.noise_floor)
        denom = np.sqrt(sd * sd + sigma_obs * sigma_obs)
        return (np.asarray(y_held, dtype=float) - mu) / np.maximum(denom, 1e-12)

    # ------------------------------------------------------------------
    # Convenience: hyper-parameter accessors (for diagnostics + tests)
    # ------------------------------------------------------------------

    @property
    def sigma_f(self) -> float:
        return self._state.sigma_f if self._state else float("nan")

    @property
    def length_scales(self) -> np.ndarray:
        return (self._state.length_scales.copy()
                if self._state else np.array([float("nan")] * self.n_dim))

    @property
    def n_train(self) -> int:
        return len(self._state.X) if self._state else 0


# ---------------------------------------------------------------------------
# Acquisition functions
# ---------------------------------------------------------------------------

def lcb(gp: MaternGP, X: np.ndarray, kappa: float = 2.0) -> np.ndarray:
    """Lower confidence bound: μ − κ·σ.  Lower is more attractive."""
    mu, sd = gp.predict(X, return_std=True)
    return mu - kappa * sd


def ei(gp: MaternGP, X: np.ndarray, y_best: float, xi: float = 0.0) -> np.ndarray:
    """Expected improvement (negated so LOWER is more attractive).

    Returns ``-E[max(0, y_best − ξ − f(X))]`` so callers can use the same
    "minimise me" convention as :func:`lcb`.
    """
    from math import sqrt
    from scipy.special import erf  # local: avoid hard dep at module load

    mu, sd = gp.predict(X, return_std=True)
    sd = np.maximum(sd, 1e-12)
    z = (y_best - xi - mu) / sd
    Phi = 0.5 * (1.0 + erf(z / sqrt(2.0)))
    phi = np.exp(-0.5 * z * z) / sqrt(2.0 * np.pi)
    improvement = (y_best - xi - mu) * Phi + sd * phi
    return -improvement


# ---------------------------------------------------------------------------
# Acquisition optimiser
# ---------------------------------------------------------------------------

def propose_next(acq_fn: Callable[[np.ndarray], np.ndarray],
                 bounds: Tuple[Tuple[float, float], Tuple[float, float]],
                 *, n_grid: int = 41, n_polish: int = 8,
                 seed: int = 0) -> Tuple[float, float]:
    """Coarse grid + local polish of an acquisition function on a 2-D box.

    ``acq_fn`` must accept an ``(N, 2)`` array and return an ``(N,)``
    array; lower is better.  ``bounds = ((r1_lo, r1_hi), (r2_lo, r2_hi))``.
    Returns the chosen ``(r1, r2)``.
    """
    (r1_lo, r1_hi), (r2_lo, r2_hi) = bounds
    g1 = np.linspace(r1_lo, r1_hi, n_grid)
    g2 = np.linspace(r2_lo, r2_hi, n_grid)
    G1, G2 = np.meshgrid(g1, g2)
    grid = np.column_stack([G1.ravel(), G2.ravel()])
    vals = acq_fn(grid)
    i = int(np.argmin(vals))
    x0 = grid[i]
    # Local polish via Nelder-Mead — tiny budget, no derivatives needed.
    rng = np.random.default_rng(seed)
    best_x, best_v = x0, float(vals[i])
    for k in range(max(0, int(n_polish))):
        try:
            res = minimize(lambda x: float(acq_fn(np.atleast_2d(x))[0]),
                           x0=best_x + 0.01 * rng.standard_normal(2),
                           method="Nelder-Mead",
                           options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 50})
        except Exception:
            continue
        # Clip to bounds.
        cand = np.clip(res.x, [r1_lo, r2_lo], [r1_hi, r2_hi])
        cand_v = float(acq_fn(np.atleast_2d(cand))[0])
        if cand_v < best_v:
            best_x, best_v = cand, cand_v
    return float(best_x[0]), float(best_x[1])
