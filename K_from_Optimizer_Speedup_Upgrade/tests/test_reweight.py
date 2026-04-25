"""Unit tests for reweight.py — pure-Python, no C++ simulator needed.

Strategy: synthesise trace data from a Gaussian-like Ising-onset
distribution where we know χ(β) analytically.  Specifically, we draw
samples (E_i, m_i) from a joint distribution where the marginal of E
is normal and m is correlated with E, then:

  * reweighting from β₀ to β₀ + Δβ should recover the Boltzmann shift
    of ⟨m²⟩ to leading order in Δβ;
  * χ(β₀) from the reweighter must match the direct estimator
    V·(⟨m²⟩ − ⟨|m|⟩²) on the same samples (Δβ = 0 sanity);
  * N_eff at Δβ = 0 must equal N (unit weights);
  * N_eff falls off and the validity gate triggers as |Δβ|·σ_E grows;
  * jackknife σ on χ is positive and bounded.
"""
from __future__ import annotations

import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG  = os.path.dirname(HERE)
if PKG not in sys.path:
    sys.path.insert(0, PKG)

import reweight as rw


# ---------------------------------------------------------------------------
# Synthetic trace generator
# ---------------------------------------------------------------------------

def _write_traces(path, e_per_site, abs_m, m2, n_sites, beta,
                  K=(1.0, 1.0, 1.0, 0.0)):
    with open(path, "w") as f:
        f.write("# E_per_site abs_m m^2  (one line per measured trajectory)\n")
        f.write(
            f"# n_sites {n_sites} beta {beta} K1 {K[0]} K2 {K[1]} K3 {K[2]} Kt {K[3]}\n"
        )
        for e, am, m_sq in zip(e_per_site, abs_m, m2):
            f.write(f"{e:.10g} {am:.10g} {m_sq:.10g}\n")


@pytest.fixture
def synthetic_traces(tmp_path):
    """Return (path, dict-of-arrays) for a single synthetic trace file."""
    rng = np.random.default_rng(2026)
    n = 4000
    n_sites = 144
    beta0 = 0.27
    # E_per_site ~ N(-1.2, 0.05); E_total = E_per_site * n_sites.
    e_ps = rng.normal(loc=-1.2, scale=0.05, size=n)
    # |m| weakly anti-correlated with E (lower energy → more order).
    am = np.clip(0.6 - 2.0 * (e_ps + 1.2) + rng.normal(0, 0.03, n), 0.0, 1.0)
    m2 = am * am + rng.normal(0, 0.005, n)
    m2 = np.clip(m2, 0.0, 1.0)
    path = tmp_path / "traces_test.dat"
    _write_traces(str(path), e_ps, am, m2, n_sites, beta0)
    return str(path), dict(E_per_site=e_ps, abs_m=am, m2=m2,
                           n_sites=n_sites, beta_parent=beta0)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_parse_traces_roundtrip(synthetic_traces):
    path, ref = synthetic_traces
    parsed = rw.parse_traces(path)
    assert parsed["n_sites"] == ref["n_sites"]
    assert parsed["beta_parent"] == pytest.approx(ref["beta_parent"])
    assert len(parsed["E_per_site"]) == len(ref["E_per_site"])
    np.testing.assert_allclose(parsed["E_per_site"], ref["E_per_site"], rtol=1e-6)
    np.testing.assert_allclose(parsed["abs_m"], ref["abs_m"], rtol=1e-6)
    np.testing.assert_allclose(parsed["m2"], ref["m2"], rtol=1e-6)


def test_parse_traces_rejects_empty(tmp_path):
    p = tmp_path / "empty.dat"
    p.write_text("# header only\n# n_sites 4 beta 0.1 K1 1 K2 1 K3 1 Kt 0\n")
    with pytest.raises(ValueError):
        rw.parse_traces(str(p))


def test_parse_traces_rejects_missing_header(tmp_path):
    p = tmp_path / "noheader.dat"
    p.write_text("# stray comment\n-1.0 0.5 0.25\n-1.1 0.6 0.36\n" * 5)
    with pytest.raises(ValueError):
        rw.parse_traces(str(p))


def test_chi_zero_dbeta_matches_direct(synthetic_traces):
    """At Δβ = 0 the reweighter must reproduce the direct χ estimator
    (per-site convention, no V factor — matches simulator m_susc)."""
    path, ref = synthetic_traces
    rw_obj = rw.Reweighter(path)
    est = rw_obj.chi(rw_obj.beta_parent)
    direct = np.mean(ref["m2"]) - np.mean(ref["abs_m"]) ** 2
    assert est.chi == pytest.approx(direct, rel=1e-10)
    # All weights are unity → N_eff equals N.
    assert est.n_eff == pytest.approx(rw_obj.n_samples, rel=1e-10)
    assert est.sigma > 0  # jackknife produces a finite bar


def test_chi_grid_matches_scalar(synthetic_traces):
    """Vectorised chi_grid must agree with point-by-point chi() calls."""
    path, _ = synthetic_traces
    rw_obj = rw.Reweighter(path)
    grid = np.linspace(rw_obj.beta_parent - 0.005,
                       rw_obj.beta_parent + 0.005, 7)
    chis_vec, _, neff_vec = rw_obj.chi_grid(grid)
    for i, b in enumerate(grid):
        est = rw_obj.chi(float(b))
        if np.isfinite(chis_vec[i]):
            assert chis_vec[i] == pytest.approx(est.chi, rel=1e-10)
            assert neff_vec[i] == pytest.approx(est.n_eff, rel=1e-10)


def test_neff_drops_with_dbeta(synthetic_traces):
    """N_eff must monotonically decrease as |Δβ| grows."""
    path, _ = synthetic_traces
    rw_obj = rw.Reweighter(path)
    # Probe out to several σ_E/N: with σ_E_total = 0.05*n_sites = 7.2,
    # Δβ = 0.5 gives Δβ·σ_E ≈ 3.6 — plenty of weight spread.
    dbs = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5])
    neffs = []
    for db in dbs:
        w = rw_obj._weights(rw_obj.beta_parent + db)
        neffs.append(rw_obj._n_eff(w * rw_obj.n_samples))
    neffs = np.asarray(neffs)
    assert neffs[0] == pytest.approx(rw_obj.n_samples, rel=1e-10)
    assert all(neffs[i] >= neffs[i + 1] - 1e-6 for i in range(len(neffs) - 1))
    assert neffs[-1] < 0.25 * neffs[0]


def test_validity_gate_returns_nan(synthetic_traces):
    """Far enough from β_parent the gate must mark χ as NaN."""
    path, _ = synthetic_traces
    rw_obj = rw.Reweighter(path, n_eff_floor=0.5)  # aggressive floor
    # Push far enough that N_eff certainly fails the 0.5 fraction.
    est = rw_obj.chi(rw_obj.beta_parent + 1.0)
    assert np.isnan(est.chi)
    assert np.isnan(est.sigma)
    assert est.n_eff < 0.5 * rw_obj.n_samples


def test_valid_window_contains_parent(synthetic_traces):
    path, _ = synthetic_traces
    rw_obj = rw.Reweighter(path)
    lo, hi = rw_obj.valid_window()
    assert lo <= rw_obj.beta_parent <= hi
    assert hi > lo


def test_jackknife_sigma_positive(synthetic_traces):
    path, _ = synthetic_traces
    rw_obj = rw.Reweighter(path, n_blocks=8)
    sigma = rw_obj._jackknife_sigma(rw_obj.beta_parent)
    assert sigma > 0.0
    # Trivially: σ on a 4000-sample histogram is bounded by χ itself.
    chi = rw_obj.chi(rw_obj.beta_parent).chi
    assert sigma < abs(chi)


def test_combined_reweighter_picks_best_parent(tmp_path):
    """CombinedReweighter must pick the parent with largest N_eff."""
    rng = np.random.default_rng(7)
    n_sites = 100
    paths = []
    for i, beta_p in enumerate([0.20, 0.27, 0.34]):
        e_ps = rng.normal(-1.2, 0.05, 2000)
        am = np.clip(0.5 + rng.normal(0, 0.05, 2000), 0, 1)
        m2 = am * am
        p = tmp_path / f"traces_{i}.dat"
        _write_traces(str(p), e_ps, am, m2, n_sites, beta_p)
        paths.append(str(p))
    rws = [rw.Reweighter(p) for p in paths]
    combo = rw.CombinedReweighter(rws)
    # At each parent's own β, that parent must win (full N_eff).
    for parent in rws:
        est = combo.chi(parent.beta_parent)
        assert est.n_eff == pytest.approx(parent.n_samples, rel=1e-6)


def test_combined_reweighter_falls_back_when_all_fail(tmp_path):
    rng = np.random.default_rng(11)
    p = tmp_path / "traces.dat"
    # Wide histogram so even a far Δβ collapses N_eff well below the floor.
    e_ps = rng.normal(-1.2, 0.5, 200)
    _write_traces(str(p), e_ps, np.full(200, 0.5), np.full(200, 0.25), 100, 0.27)
    combo = rw.CombinedReweighter(
        [rw.Reweighter(str(p), n_eff_floor=0.5)],
        n_eff_floor=0.5,
    )
    est = combo.chi(10.0)  # absurdly far → no parent passes the 0.5 floor
    assert np.isnan(est.chi)
    assert est.n_eff == 0.0


def test_too_few_samples_rejected(tmp_path):
    p = tmp_path / "tiny.dat"
    _write_traces(str(p), [-1.0] * 4, [0.5] * 4, [0.25] * 4, 100, 0.27)
    with pytest.raises(ValueError):
        rw.Reweighter(str(p))


def test_find_traces_file(tmp_path):
    sub = tmp_path / "run"
    sub.mkdir()
    (sub / "traces_DEADBEEF.dat").write_text("# stub\n")
    (sub / "two_point_all_to_all.dat").write_text("# unrelated\n")
    found = rw.find_traces_file(str(sub))
    assert found is not None and found.endswith("traces_DEADBEEF.dat")
    assert rw.find_traces_file(str(tmp_path / "nope")) is None
