#!/usr/bin/env python3
"""
modular_data.py — three-direction torus correlator profiles for the
modular testbed.

For a fixed reference period |ω₁| and a target τ = Re τ + i Im τ, build
the analytic Ising spin-spin two-point function F(z, τ) along the three
shortest torus geodesics:

    "re"   direction — along ω₁ axis  (length |ω₁|)
                       = the **Re τ** cycle
    "im"   direction — perpendicular axis (i)  (length |Im τ|·|ω₁|)
                       = the **Im τ** cycle
    "diag" direction — along (ω₂ − ω₁)
                       = the **Im τ − Re τ** cycle (after one twist)

Each profile is a 1-D array G[s] for s = 1, 2, … K_min lattice spacings,
where K_min is bounded by the shortest half-period.  This sidesteps the
boundary-path / gcd machinery entirely: every τ point on the modular
grid stores three short vectors of analytic correlator values, all
referred to the SAME physical units (lattice spacing 1).

Costs:
    * l2          : mean( (G_t - G_r)² )                        averaged over 3 dirs
    * relative    : mean( ((G_t - G_r) / G_r)² )
    * log         : mean( (log G_t - log G_r)² )
    * effmass     : mean( (m_eff_t - m_eff_r)² ) with m_eff_k = -log(G_{k+1}/G_k)
    * pmean4      : ( mean_d (mean_k diff_k²)² )^(1/2)          p=4 over directions

Noise:
    add_noise(profile, sigma, rng)  — multiplicative Gaussian on every G entry
    A deterministic seed per (i, j) grid cell keeps the cost reproducible.
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np

from cft_torus import ising_torus_F

DIRS = ("re", "im", "diag")


# ---------------------------------------------------------------------------
# Build one three-direction profile at a given τ
# ---------------------------------------------------------------------------
def make_threedir_profile(tau: complex, *, w1_norm: float = 1.0,
                          K_max: int = 8, eps: float = 1e-12) -> dict:
    """Return {profiles: {dir: G[K]}, dists: {dir: s[K]}, tau, omegas}.

    The three direction unit vectors are:
        u_re   = ω₁ / |ω₁|
        u_im   = i   (perpendicular to a real ω₁)
        u_diag = (ω₂ - ω₁) / |ω₂ - ω₁|

    Sample at integer physical distances s = 1, 2, …, K_eff with
        K_eff = min(K_max, floor( min(|ω₁|, |Im τ·ω₁|, |ω₂-ω₁|) / 2 ))
    """
    w1 = complex(w1_norm, 0.0)
    w2 = complex(tau.real, tau.imag) * w1

    L_re   = abs(w1)
    L_im   = abs(complex(0.0, tau.imag) * w1)        # |Im τ| · |ω₁|
    L_diag = abs(w2 - w1)
    half   = min(L_re, L_im, L_diag) / 2.0
    K_eff  = int(min(K_max, max(2, math.floor(half))))

    if L_re   < eps: u_re   = complex(1.0, 0.0)
    else:            u_re   = w1 / L_re
    u_im = complex(0.0, 1.0)
    if L_diag < eps: u_diag = complex(1.0, 0.0)
    else:            u_diag = (w2 - w1) / L_diag

    s = np.arange(1, K_eff + 1, dtype=float)
    profiles, dists = {}, {}
    for label, u in (("re", u_re), ("im", u_im), ("diag", u_diag)):
        zs = s * u
        Gs = np.empty(K_eff, dtype=float)
        for k, z in enumerate(zs):
            zeta = z / w1
            Gs[k] = ising_torus_F(zeta, tau)
        profiles[label] = Gs
        dists[label] = s.copy()
    return dict(profiles=profiles, dists=dists,
                tau=complex(tau.real, tau.imag),
                omega1=w1, omega2=w2,
                K=K_eff)


# ---------------------------------------------------------------------------
# Synthetic noise injection
# ---------------------------------------------------------------------------
def add_noise(profile: dict, sigma_rel: float,
              rng: np.random.Generator) -> dict:
    """Return a NEW profile with G_k → G_k · (1 + σ·N(0,1)) per entry,
    plus an ``err`` array storing the imposed σ_rel · |G_k|.  Pass
    sigma_rel = 0 to get a noise-free copy.
    """
    out_prof, out_err = {}, {}
    for label in DIRS:
        G = profile["profiles"][label]
        if sigma_rel <= 0:
            out_prof[label] = G.copy()
            out_err[label]  = np.zeros_like(G)
        else:
            noise = rng.standard_normal(G.shape) * sigma_rel
            out_prof[label] = G * (1.0 + noise)
            out_err[label]  = sigma_rel * np.abs(G)
    new = dict(profile)
    new["profiles"] = out_prof
    new["errors"]   = out_err
    new["sigma_rel"] = float(sigma_rel)
    return new


# ---------------------------------------------------------------------------
# Costs operating on three-direction profiles
# ---------------------------------------------------------------------------
def _per_dir_l2(Gt, Gr, et=None, er=None, p=2):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    diff = Gt[:K] - Gr[:K]
    if p == 2:
        return float(np.mean(diff ** 2))
    return float(np.mean(diff ** p))


def _per_dir_relative(Gt, Gr, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    denom = np.abs(Gr[:K]) + 1e-30
    return float(np.mean(((Gt[:K] - Gr[:K]) / denom) ** 2))


def _per_dir_log(Gt, Gr, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 2: return float("nan")
    d = np.log(Gt_[m]) - np.log(Gr_[m])
    return float(np.mean(d ** 2))


def _per_dir_effmass(Gt, Gr, **_):
    K = min(len(Gt), len(Gr))
    if K < 3: return float("nan")
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 3: return float("nan")
    Gt_ = Gt_[m]; Gr_ = Gr_[m]
    mt = -np.log(Gt_[1:] / Gt_[:-1])
    mr = -np.log(Gr_[1:] / Gr_[:-1])
    return float(np.mean((mt - mr) ** 2))


# --- 10 additional cost kernels for stress testing -----------------------
def _per_dir_l1(Gt, Gr, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    return float(np.mean(np.abs(Gt[:K] - Gr[:K])))


def _per_dir_l1_log(Gt, Gr, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 2: return float("nan")
    return float(np.mean(np.abs(np.log(Gt_[m]) - np.log(Gr_[m]))))


def _per_dir_huber(Gt, Gr, delta=None, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    r = Gt[:K] - Gr[:K]
    if delta is None: delta = max(1e-6, np.median(np.abs(r)) * 1.4826)
    a = np.abs(r)
    quad = np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))
    return float(np.mean(quad))


def _per_dir_huber_log(Gt, Gr, delta=0.5, **_):
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 2: return float("nan")
    r = np.log(Gt_[m]) - np.log(Gr_[m])
    a = np.abs(r)
    quad = np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))
    return float(np.mean(quad))


def _per_dir_chi2(Gt, Gr, **_):
    """Poisson-like χ²: (G_t - G_r)² / |G_r| — emphasises tail bins where
    G is small (hardest to measure)."""
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    denom = np.abs(Gr[:K]) + 1e-30
    return float(np.mean((Gt[:K] - Gr[:K]) ** 2 / denom))


def _per_dir_cosine(Gt, Gr, **_):
    """Scale-invariant shape match: 1 - <G_t,G_r>/(|G_t||G_r|)."""
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    a = Gt[:K]; b = Gr[:K]
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na < 1e-30 or nb < 1e-30: return float("nan")
    return float(1.0 - np.dot(a, b) / (na * nb))


def _per_dir_slope(Gt, Gr, **_):
    """Squared diff of fitted slope of log G vs s — global m_eff."""
    K = min(len(Gt), len(Gr))
    if K < 3: return float("nan")
    s = np.arange(1, K + 1, dtype=float)
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 3: return float("nan")
    sm = s[m]; lt = np.log(Gt_[m]); lr = np.log(Gr_[m])
    # least-squares slope
    st = np.polyfit(sm, lt, 1)[0]
    sr = np.polyfit(sm, lr, 1)[0]
    return float((st - sr) ** 2)


def _per_dir_trim_rel(Gt, Gr, frac=0.2, **_):
    """Relative L2 with top `frac` residuals trimmed."""
    K = min(len(Gt), len(Gr))
    if K < 3: return float("nan")
    denom = np.abs(Gr[:K]) + 1e-30
    r2 = ((Gt[:K] - Gr[:K]) / denom) ** 2
    n_keep = max(2, int(np.ceil(K * (1 - frac))))
    kept = np.sort(r2)[:n_keep]
    return float(np.mean(kept))


def _per_dir_median_log(Gt, Gr, **_):
    """Median squared log residual — very robust to outliers."""
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    Gt_, Gr_ = Gt[:K], Gr[:K]
    m = (Gt_ > 1e-12) & (Gr_ > 1e-12)
    if m.sum() < 2: return float("nan")
    d = np.log(Gt_[m]) - np.log(Gr_[m])
    return float(np.median(d ** 2))


def _per_dir_rel_l1(Gt, Gr, **_):
    """Relative L1: mean |Δ/G_r|."""
    K = min(len(Gt), len(Gr))
    if K < 2: return float("nan")
    denom = np.abs(Gr[:K]) + 1e-30
    return float(np.mean(np.abs((Gt[:K] - Gr[:K]) / denom)))


_PER_DIR = {
    "l2":         _per_dir_l2,
    "relative":   _per_dir_relative,
    "log":        _per_dir_log,
    "effmass":    _per_dir_effmass,
    "l1":         _per_dir_l1,
    "l1_log":     _per_dir_l1_log,
    "huber":      _per_dir_huber,
    "huber_log":  _per_dir_huber_log,
    "chi2":       _per_dir_chi2,
    "cosine":     _per_dir_cosine,
    "slope":      _per_dir_slope,
    "trim_rel":   _per_dir_trim_rel,
    "median_log": _per_dir_median_log,
    "rel_l1":     _per_dir_rel_l1,
}


def threedir_cost(ref: dict, test: dict, *, mode: str = "l2",
                  pmean: int = 1) -> Tuple[float, Dict[str, float]]:
    """Compare two three-direction profiles.  Returns (cost, per_dir_dict).

    ``mode`` ∈ {l2, relative, log, effmass}.
    ``pmean`` aggregates per-direction values via ( mean(C_d^p) )^(1/p).
    pmean=1 → arithmetic mean; pmean=4 → anisotropy-penalising p4-mean.
    """
    fn = _PER_DIR.get(mode)
    if fn is None:
        raise KeyError(f"unknown threedir mode {mode!r}")
    per = {}
    vals = []
    for label in DIRS:
        Gr = ref["profiles"].get(label)
        Gt = test["profiles"].get(label)
        if Gr is None or Gt is None:
            per[label] = float("nan"); continue
        c = fn(Gt, Gr)
        per[label] = c
        if math.isfinite(c) and c >= 0:
            vals.append(c)
    if not vals:
        return float("nan"), per
    if pmean == 1:
        cost = sum(vals) / len(vals)
    else:
        cost = (sum(v ** pmean for v in vals) / len(vals)) ** (1.0 / pmean)
    return float(cost), per


# Convenience registry mirroring cost_zoo.ZOO
THREEDIR_COSTS = {
    "l2":       lambda r, t: threedir_cost(r, t, mode="l2"),
    "relative": lambda r, t: threedir_cost(r, t, mode="relative"),
    "log":      lambda r, t: threedir_cost(r, t, mode="log"),
    "effmass":  lambda r, t: threedir_cost(r, t, mode="effmass"),
    "pmean4":   lambda r, t: threedir_cost(r, t, mode="l2", pmean=4),
}


# ---------------------------------------------------------------------------
# Modular-aware distance: account for T : τ → τ ± 1 (and S : τ → -1/τ)
# equivalences in the upper half plane.  CMA-ES converging to a modular
# image of truth is a SUCCESS for the cost-surface evaluation, even though
# the Euclidean distance is large.
# ---------------------------------------------------------------------------
def modular_distance(tau1: complex, tau2: complex, *,
                      n_T: int = 2, include_S: bool = True) -> float:
    """Minimum Euclidean distance between τ1 and any image of τ2 under the
    modular group generators T^k (k ∈ [-n_T, n_T]) and S∘T^k.
    """
    z1 = complex(tau1.real, tau1.imag)
    z2 = complex(tau2.real, tau2.imag)
    candidates = []
    for k in range(-n_T, n_T + 1):
        candidates.append(z2 + complex(k, 0))
    if include_S:
        for k in range(-n_T, n_T + 1):
            try:
                Sz = -1.0 / (z2 + complex(k, 0))
                if Sz.imag > 0: candidates.append(Sz)
            except ZeroDivisionError:
                pass
    return float(min(abs(z1 - c) for c in candidates))


def add_noise_frozen(profile: dict, sigma_rel: float, cell_seed: int) -> dict:
    """Like add_noise() but with a deterministic seed tied ONLY to the
    cell index, so the same (i, j) cell always returns the same noisy
    profile across CMA-ES calls.  Models 'static landscape roughness'
    rather than re-sampled stochastic eval noise.
    """
    rng = np.random.default_rng(cell_seed)
    return add_noise(profile, sigma_rel, rng)
