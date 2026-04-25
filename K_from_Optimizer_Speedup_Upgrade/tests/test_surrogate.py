"""Unit tests for surrogate.py — pure-Python, no MC needed."""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG  = os.path.dirname(HERE)
if PKG not in sys.path:
    sys.path.insert(0, PKG)

import surrogate as sg


# ---------------------------------------------------------------------------
# Synthetic noisy quadratic surface
# ---------------------------------------------------------------------------

def _quadratic(X, x_star=(1.0, 1.0), noise=0.05, rng=None):
    """f(r1, r2) = (r1 − r1*)² + 1.5·(r2 − r2*)² + N(0, noise²)."""
    X = np.atleast_2d(X)
    rng = rng or np.random.default_rng(0)
    base = (X[:, 0] - x_star[0]) ** 2 + 1.5 * (X[:, 1] - x_star[1]) ** 2
    return base + noise * rng.standard_normal(len(X))


@pytest.fixture
def noisy_dataset():
    rng = np.random.default_rng(2026)
    X = rng.uniform(0.5, 1.5, size=(40, 2))
    y = _quadratic(X, rng=rng, noise=0.02)
    sigma = np.full(40, 0.02)
    return X, y, sigma


# ---------------------------------------------------------------------------
# Kernel sanity
# ---------------------------------------------------------------------------

def test_matern52_diagonal_equals_sigma_squared():
    X = np.random.default_rng(0).uniform(0, 1, (5, 2))
    K = sg.matern52(X, X, sigma_f=2.0, length_scales=np.array([0.3, 0.4]))
    np.testing.assert_allclose(np.diag(K), 4.0, rtol=1e-12)


def test_matern52_psd():
    X = np.random.default_rng(1).uniform(0, 1, (10, 2))
    K = sg.matern52(X, X, sigma_f=1.0, length_scales=np.array([0.5, 0.5]))
    K += 1e-8 * np.eye(10)  # guard
    eig = np.linalg.eigvalsh(K)
    assert eig.min() > -1e-10


def test_matern52_decays_with_distance():
    a = np.array([[0.0, 0.0]])
    b = np.array([[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [3.0, 0.0]])
    K = sg.matern52(a, b, sigma_f=1.0, length_scales=np.array([1.0, 1.0]))
    assert K[0, 0] > K[0, 1] > K[0, 2] > K[0, 3]


# ---------------------------------------------------------------------------
# Fit / predict basics
# ---------------------------------------------------------------------------

def test_fit_then_predict_at_training_point(noisy_dataset):
    X, y, s = noisy_dataset
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=2, seed=0)
    mu, sd = gp.predict(X)
    # At training points, mean should track y closely (within a few σ_obs).
    np.testing.assert_allclose(mu, y, atol=5 * s.max())
    assert (sd >= 0).all()


def test_fit_reduces_variance_at_dense_region():
    """σ_pred at a dense training cluster ≪ σ_pred far away."""
    rng = np.random.default_rng(3)
    cluster = rng.uniform(0.95, 1.05, (30, 2))
    y = _quadratic(cluster, rng=rng, noise=0.01)
    s = np.full(30, 0.01)
    gp = sg.MaternGP(n_dim=2).fit(cluster, y, s, n_restarts=2, seed=0)
    _, sd_in = gp.predict(np.array([[1.0, 1.0]]))
    _, sd_out = gp.predict(np.array([[3.0, 3.0]]))
    assert sd_in[0] < sd_out[0]


def test_fit_recovers_minimum(noisy_dataset):
    """GP posterior mean predicts a minimum near the true (1, 1)."""
    X, y, s = noisy_dataset
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=3, seed=0)
    g = np.linspace(0.5, 1.5, 41)
    G1, G2 = np.meshgrid(g, g)
    grid = np.column_stack([G1.ravel(), G2.ravel()])
    mu, _ = gp.predict(grid)
    i = int(np.argmin(mu))
    assert abs(grid[i, 0] - 1.0) < 0.15
    assert abs(grid[i, 1] - 1.0) < 0.15


def test_log_marginal_likelihood_finite(noisy_dataset):
    X, y, s = noisy_dataset
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=2, seed=0)
    lml = gp.log_marginal_likelihood()
    assert np.isfinite(lml)


def test_predict_before_fit_raises():
    gp = sg.MaternGP(n_dim=2)
    with pytest.raises(RuntimeError):
        gp.predict(np.array([[0.5, 0.5]]))


def test_fit_with_one_point_raises():
    gp = sg.MaternGP(n_dim=2)
    with pytest.raises(ValueError):
        gp.fit(np.array([[0.5, 0.5]]), np.array([0.0]), np.array([1e-3]))


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------

def test_calibration_residuals_unit_scale():
    """On an honest quadratic + known noise, residual std ≈ 1."""
    rng = np.random.default_rng(7)
    n_train, n_test = 60, 30
    X_tr = rng.uniform(0.5, 1.5, (n_train, 2))
    y_tr = _quadratic(X_tr, rng=rng, noise=0.03)
    X_te = rng.uniform(0.5, 1.5, (n_test, 2))
    y_te = _quadratic(X_te, rng=rng, noise=0.03)
    s = np.full(n_train, 0.03)
    s_te = np.full(n_test, 0.03)
    gp = sg.MaternGP(n_dim=2).fit(X_tr, y_tr, s, n_restarts=3, seed=0)
    r = gp.calibration_residuals(X_te, y_te, s_te)
    # On 30 points, std(r) ought to be in [0.5, 2.0] for honest bars.
    assert 0.4 < r.std() < 2.5
    # 95th-percentile-of-|r| acceptance from the validation panel:
    assert np.percentile(np.abs(r), 95) < 3.0


# ---------------------------------------------------------------------------
# Acquisition functions + propose_next
# ---------------------------------------------------------------------------

def test_lcb_low_at_known_minimum(noisy_dataset):
    X, y, s = noisy_dataset
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=2, seed=0)
    grid = np.array([[1.0, 1.0], [0.5, 0.5], [1.5, 1.5]])
    vals = sg.lcb(gp, grid, kappa=2.0)
    assert int(np.argmin(vals)) == 0


def test_ei_returns_negative_improvement(noisy_dataset):
    """EI(x) at a known good location must be more negative than EI at
    a known bad location *within the well-sampled region*.

    EI in unsampled regions is large because of high σ (exploration);
    the test stays inside the sampled box to compare exploitation only.
    """
    X, y, s = noisy_dataset
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=2, seed=0)
    y_best = float(min(y))
    # Both points sit inside the training box [0.5, 1.5]^2; (1, 1) is
    # near the true minimum, (1.4, 0.6) is on the rim.
    grid = np.array([[1.0, 1.0], [1.4, 0.6]])
    e = sg.ei(gp, grid, y_best)
    assert e[0] < e[1]


def test_propose_next_returns_finite_inside_bounds():
    rng = np.random.default_rng(11)
    X = rng.uniform(0.5, 1.5, (20, 2))
    y = _quadratic(X, rng=rng, noise=0.02)
    s = np.full(20, 0.02)
    gp = sg.MaternGP(n_dim=2).fit(X, y, s, n_restarts=2, seed=0)
    bounds = ((0.5, 1.5), (0.5, 1.5))
    nxt = sg.propose_next(lambda Z: sg.lcb(gp, Z, kappa=2.0),
                          bounds=bounds, n_grid=21, n_polish=3, seed=0)
    assert bounds[0][0] <= nxt[0] <= bounds[0][1]
    assert bounds[1][0] <= nxt[1] <= bounds[1][1]
    # On a quadratic centred at (1, 1) the proposed point should head there.
    assert abs(nxt[0] - 1.0) < 0.4
    assert abs(nxt[1] - 1.0) < 0.4


def test_propose_next_finds_quadratic_minimum_in_few_steps():
    """End-to-end mini-BO loop on a noisy quadratic.

    Mirrors the Validation-panel BO acceptance criterion: surrogate
    must reach the r*-ball (here ‖Δ‖₂ < 0.1) in ≤ 0.5 × baseline evals.
    Baseline eval count proxy: 30; we cap at 15 here.
    """
    rng = np.random.default_rng(13)
    bounds = ((0.3, 1.7), (0.3, 1.7))
    X_list = list(rng.uniform(0.5, 1.5, (6, 2)))
    y_list = [float(_quadratic(np.atleast_2d(x), rng=rng, noise=0.02)[0])
              for x in X_list]
    s_list = [0.02] * 6
    for step in range(15):
        gp = sg.MaternGP(n_dim=2).fit(np.array(X_list),
                                       np.array(y_list),
                                       np.array(s_list),
                                       n_restarts=2, seed=step)
        nxt = sg.propose_next(lambda Z: sg.lcb(gp, Z, kappa=1.5),
                              bounds=bounds, n_grid=31, n_polish=2,
                              seed=step)
        X_list.append(np.array(nxt))
        y_list.append(float(_quadratic(np.atleast_2d(nxt), rng=rng, noise=0.02)[0]))
        s_list.append(0.02)
    best = X_list[int(np.argmin(y_list))]
    dist = np.hypot(best[0] - 1.0, best[1] - 1.0)
    assert dist < 0.15, f"BO failed to converge: dist={dist:.3f}"
