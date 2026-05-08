#!/usr/bin/env python3
"""
cft_torus.py — analytic torus 2-point function for testing cost variants
without MC noise.

Uses a massive free scalar on a parallelogram torus, computed by direct
momentum-space sum:

    G_m(r; ω₁, ω₂) = (1/A) Σ_{n₁,n₂ ∈ Z}  exp(i k·r) / (k² + m²)

where A = |Im(conj(ω₁)·ω₂)| is the torus area and the momenta live on
the reciprocal lattice  k = 2π · (n₁ b₁ + n₂ b₂),  with  b_i · ω_j = δ_ij.

The mass m > 0 acts as a soft IR regulator so the n=0 mode is finite.
For m → 0 this reduces to the GFF (up to the missing zero-mode contact);
the modular structure (dependence on τ = ω₂/ω₁) is identical to the
Ising CFT spin-spin correlator's leading term for the cost-discrimination
test we care about.

The output dict has the same shape as ``mc_engine.load_all_to_all`` so the
existing cost functions (production + zoo) can consume it directly.
"""
from __future__ import annotations

import math
import cmath
from typing import Dict, Tuple

import numpy as np

_SQRT3_2 = math.sqrt(3.0) / 2.0


def torus_periods(Lx: int, Ly: int, Tx: int, Ty: int) -> Tuple[complex, complex]:
    """Return torus periods (ω₁, ω₂) as complex numbers in physical xy.

    Same convention as ``cost.boundary_paths``: triangular-lattice
    fundamental periods are
        v = (Lx + Ty/2,  Ty · √3/2)
        u = (Tx − Ly/2, −Ly · √3/2)
    treated as z = x + i y.
    """
    w1 = complex(Lx + 0.5 * Ty,  _SQRT3_2 * Ty)
    w2 = complex(Tx - 0.5 * Ly, -_SQRT3_2 * Ly)
    return w1, w2


def _reciprocal(w1: complex, w2: complex) -> Tuple[complex, complex]:
    """Reciprocal-lattice basis (b₁, b₂) with b_i · ω_j = δ_ij in ℝ²."""
    # Treat 2-D real vectors as complex; the real dot product is
    # Re(conj(b)·ω). Solve the linear system.
    A = np.array([[w1.real, w2.real],
                  [w1.imag, w2.imag]], dtype=float)
    inv = np.linalg.inv(A)        # rows of inv are (b1, b2)
    b1 = complex(inv[0, 0], inv[0, 1])
    b2 = complex(inv[1, 0], inv[1, 1])
    return b1, b2


def _torus_area(w1: complex, w2: complex) -> float:
    return abs(w1.real * w2.imag - w1.imag * w2.real)


def torus_green(z: complex, w1: complex, w2: complex,
                m_sq: float = 0.01, n_cut: int = 20) -> float:
    """Massive free-field Green function at separation z = x + i y.

    G_m(z) = (1/A) Σ_{n₁,n₂}  cos(k·r) / (k² + m²)
    (use cos because the sum is symmetric n → −n; G is real)
    """
    b1, b2 = _reciprocal(w1, w2)
    A = _torus_area(w1, w2)
    total = 0.0
    rx, ry = z.real, z.imag
    for n1 in range(-n_cut, n_cut + 1):
        kx = 2.0 * math.pi * (n1 * b1.real)
        ky = 2.0 * math.pi * (n1 * b1.imag)
        for n2 in range(-n_cut, n_cut + 1):
            kxn = kx + 2.0 * math.pi * (n2 * b2.real)
            kyn = ky + 2.0 * math.pi * (n2 * b2.imag)
            k2 = kxn * kxn + kyn * kyn
            phase = kxn * rx + kyn * ry
            total += math.cos(phase) / (k2 + m_sq)
    return total / A


def torus_green_array(zs_x: np.ndarray, zs_y: np.ndarray,
                      w1: complex, w2: complex,
                      m_sq: float = 0.01, n_cut: int = 20) -> np.ndarray:
    """Vectorised G_m(z) over an array of separations."""
    b1, b2 = _reciprocal(w1, w2)
    A = _torus_area(w1, w2)
    n1s = np.arange(-n_cut, n_cut + 1)
    n2s = np.arange(-n_cut, n_cut + 1)
    N1, N2 = np.meshgrid(n1s, n2s, indexing="ij")
    kx = 2.0 * math.pi * (N1 * b1.real + N2 * b2.real)
    ky = 2.0 * math.pi * (N1 * b1.imag + N2 * b2.imag)
    k2 = kx * kx + ky * ky
    denom = k2 + m_sq                            # shape (2N+1, 2N+1)
    # Outer over (sites, modes)
    rx = zs_x[:, None, None]
    ry = zs_y[:, None, None]
    phase = kx[None, :, :] * rx + ky[None, :, :] * ry
    return np.sum(np.cos(phase) / denom[None, :, :], axis=(1, 2)) / A


def lattice_sites(Lx: int, Ly: int, Tx: int, Ty: int) -> Dict[Tuple[int, int],
                                                              Tuple[float, float]]:
    """Integer (m, n) sites in the test fundamental domain.

    Uses the same enumeration as the simulator's two_point_all_to_all.dat:
    a square block of (m, n) ∈ [-something, something]² that covers the
    fundamental cell after torus reduction.  For the zoo cost test we
    only need sites along the three boundary directions (the cost
    functions index test_data via _lookup_test_value), so we generate a
    generous (m, n) box around the origin.
    """
    w1, w2 = torus_periods(Lx, Ly, Tx, Ty)
    # Generate (m, n) in a box; positions = m + n/2 + i n √3/2
    M = 2 * max(Lx, Ly, abs(Tx), abs(Ty))
    sites = {}
    for m in range(-M, M + 1):
        for n in range(-M, M + 1):
            sites[(m, n)] = (m + 0.5 * n, _SQRT3_2 * n)
    return sites


def make_data_dict(Lx: int, Ly: int, Tx: int, Ty: int, *,
                   omega1_override: complex | None = None,
                   omega2_override: complex | None = None,
                   m_sq: float = 0.01, n_cut: int = 18,
                   err: float = 1e-6) -> Dict[Tuple[int, int], dict]:
    """Build a {(m, n): {'conn': G, 'conn_err': err}} dict with G computed
    at each integer site (m, n) under the (possibly overridden) torus
    periods.

    To scan a 2-D family of effective τ for the cost-zoo test, override
    ω₁ and/or ω₂ with the desired effective periods while keeping the
    site enumeration on the literal (Lx, Ly, Tx, Ty) integer lattice.
    """
    w1, w2 = torus_periods(Lx, Ly, Tx, Ty)
    if omega1_override is not None: w1 = omega1_override
    if omega2_override is not None: w2 = omega2_override

    sites = lattice_sites(Lx, Ly, Tx, Ty)
    keys  = list(sites.keys())
    xs = np.array([sites[k][0] for k in keys], dtype=float)
    ys = np.array([sites[k][1] for k in keys], dtype=float)
    Gs = torus_green_array(xs, ys, w1, w2, m_sq=m_sq, n_cut=n_cut)
    return {k: {"conn": float(g), "conn_err": err} for k, g in zip(keys, Gs)}


def modular_tau(Lx: int, Ly: int, Tx: int, Ty: int) -> complex:
    """τ = ω₂/ω₁ of the torus, conjugated as needed so Im(τ) > 0."""
    w1, w2 = torus_periods(Lx, Ly, Tx, Ty)
    tau = w2 / w1
    if tau.imag < 0:
        tau = complex(tau.real, -tau.imag)
    return tau


def modular_reduce(tau: complex, max_iter: int = 100) -> complex:
    """Reduce τ to the standard fundamental domain |τ|≥1, |Re(τ)|≤1/2,
    Im(τ)>0 via PSL(2, ℤ) generators T (τ→τ+1) and S (τ→−1/τ)."""
    if tau.imag <= 0:
        tau = complex(tau.real, abs(tau.imag) + 1e-12)
    for _ in range(max_iter):
        # Translate Re(τ) into [-1/2, 1/2]
        t = tau.real - round(tau.real)
        tau = complex(t, tau.imag)
        if abs(tau) < 1.0:
            tau = -1.0 / tau
        else:
            return tau
    return tau


# ===========================================================================
# Ising CFT spin-spin correlator on the torus (Jacobi theta functions)
# ===========================================================================
#
# Ref: Di Francesco, Mathieu, Sénéchal §12.6, Eq. (12.106-117).  The
# Ising spin operator σ has weights h = h̄ = 1/16.  On a torus with
# modular parameter τ = ω₂/ω₁ and complex separation z (in units of ω₁,
# i.e. ζ = z/ω₁), the spin-spin correlator summed over the three even
# spin structures (a, b) ∈ {(0,0), (0,1/2), (1/2,0)} is
#
#   ⟨σ(z) σ(0)⟩_τ  =  (1/Z_Ising) · |2π η(τ)|^(-1/4) ·
#                     [|θ₁(ζ|τ)|^(-1/4) ·
#                      ( |θ₂(ζ|τ)|^(1/2) · |θ₂(0|τ)|^(1/2)
#                      + |θ₃(ζ|τ)|^(1/2) · |θ₃(0|τ)|^(1/2)
#                      + |θ₄(ζ|τ)|^(1/2) · |θ₄(0|τ)|^(1/2) ) ]
#                     / (2 |θ₂(0|τ) θ₃(0|τ) θ₄(0|τ)|^(1/4))
#
# where η(τ) is the Dedekind eta function and θ_a(ζ|τ) are the Jacobi
# theta functions in the q = exp(2πiτ) convention.  The overall
# normalisation Z_Ising and the η factor cancel in the *cost*
# (which is a difference of correlators across geometries) so we drop
# them and return the unnormalised
#
#   F(ζ; τ) = |θ₁(ζ|τ)|^(-1/4) ·
#             ( |θ₂(ζ|τ)|^(1/2) · |θ₂(0|τ)|^(1/2)
#             + |θ₃(ζ|τ)|^(1/2) · |θ₃(0|τ)|^(1/2)
#             + |θ₄(ζ|τ)|^(1/2) · |θ₄(0|τ)|^(1/2) )
#
# This has the EXACT modular structure of the Ising spin correlator —
# what the cost-collapse hypothesis is supposed to detect — without the
# overall constants that the cost is invariant to.

# ---------------------------------------------------------------------------
# Jacobi theta functions θ_a(z|τ),  q = exp(iπτ),  z complex.
# ---------------------------------------------------------------------------
# Series form (DMS Eq. 10.158-161), truncated at |n| ≤ N_THETA.
# At Im(τ) ≳ 0.3 the series converges to ~16-digit precision by N=20.
N_THETA = 24


def _jtheta_all(z: complex, tau: complex):
    """Return (θ₁, θ₂, θ₃, θ₄) at (z|τ) using nome q = exp(iπτ).

    θ₁(z|τ) =  2 q^(1/4) Σ_{n=0..∞} (-1)^n q^(n(n+1))    sin((2n+1) z)
    θ₂(z|τ) =  2 q^(1/4) Σ_{n=0..∞}            q^(n(n+1)) cos((2n+1) z)
    θ₃(z|τ) =  1 + 2 Σ_{n=1..∞} q^(n²)                 cos(2 n z)
    θ₄(z|τ) =  1 + 2 Σ_{n=1..∞} (-1)^n q^(n²)          cos(2 n z)
    """
    q = cmath.exp(1j * cmath.pi * tau)        # nome
    q14 = cmath.exp(1j * cmath.pi * tau / 4)  # q^(1/4)
    th1 = 0 + 0j
    th2 = 0 + 0j
    th3 = 1 + 0j
    th4 = 1 + 0j
    for n in range(N_THETA + 1):
        # θ₁, θ₂ use q^{n(n+1)} cos/sin((2n+1) z)
        qpow_a = q ** (n * (n + 1))
        sign1  = (-1) ** n
        ang_a  = (2 * n + 1) * z
        th1   += sign1 * qpow_a * cmath.sin(ang_a)
        th2   +=         qpow_a * cmath.cos(ang_a)
        if n >= 1:
            qpow_b = q ** (n * n)
            ang_b  = 2 * n * z
            th3   += 2 * qpow_b * cmath.cos(ang_b)
            th4   += 2 * ((-1) ** n) * qpow_b * cmath.cos(ang_b)
    th1 = 2 * q14 * th1
    th2 = 2 * q14 * th2
    return th1, th2, th3, th4


def ising_torus_F(zeta: complex, tau: complex) -> float:
    """Unnormalised Ising spin-spin two-point function on the torus.

    zeta = z/ω₁ is the complex separation in units of ω₁ (real part along
    the "spatial" cycle, imag part along the "time" cycle).  τ = ω₂/ω₁
    must have Im(τ) > 0.

    Returns a positive real number.  Cost functions only see relative
    differences across geometries, so overall normalisation is irrelevant.
    """
    # Convert zeta from "fraction of ω₁" to the theta-function convention
    # in which the period is π (DMS uses θ_a(z|τ) periodic in z → z + π
    # for θ₃, θ₄; θ₁, θ₂ pick up minus signs).
    z = cmath.pi * zeta
    th1, th2, th3, th4 = _jtheta_all(z, tau)
    # |θ₁(0|τ)| = 0, so the |θ₁(ζ|τ)|^(-1/4) factor diverges as ζ→0.
    # Cap it at a small floor to avoid numerical inf at the contact point.
    a1 = max(abs(th1), 1e-30)
    th2_0, th3_0, th4_0 = abs(th2.real if abs(th2.imag) < 1e-12 else th2), \
                          abs(th3.real if abs(th3.imag) < 1e-12 else th3), \
                          abs(th4.real if abs(th4.imag) < 1e-12 else th4)
    # θ_a(0|τ) for a=2,3,4
    _, th2_at_0, th3_at_0, th4_at_0 = _jtheta_all(0 + 0j, tau)
    s = (math.sqrt(abs(th2)) * math.sqrt(abs(th2_at_0))
       + math.sqrt(abs(th3)) * math.sqrt(abs(th3_at_0))
       + math.sqrt(abs(th4)) * math.sqrt(abs(th4_at_0)))
    return float(a1 ** (-0.25) * s)


def ising_torus_array(zs_x: np.ndarray, zs_y: np.ndarray,
                      w1: complex, w2: complex) -> np.ndarray:
    """Vectorised Ising spin correlator at lattice sites z = x + iy."""
    tau = w2 / w1
    if tau.imag < 0:
        tau = complex(tau.real, -tau.imag)
        w1, w2 = complex(w1.real, -w1.imag), complex(w2.real, -w2.imag)
    out = np.empty(zs_x.shape[0], dtype=float)
    for i, (x, y) in enumerate(zip(zs_x, zs_y)):
        zeta = complex(x, y) / w1
        out[i] = ising_torus_F(zeta, tau)
    return out


def make_ising_data_dict(Lx: int, Ly: int, Tx: int, Ty: int, *,
                         omega1_override: complex | None = None,
                         omega2_override: complex | None = None,
                         err: float = 1e-6) -> Dict[Tuple[int, int], dict]:
    """Same shape as ``make_data_dict`` but uses the Ising CFT
    spin-spin correlator (Jacobi theta) instead of the massive scalar.

    Subtracts the value at the largest |z| in the box so the result has
    the same "connected" sign convention used by the production cost
    functions (G_conn → 0 at large separation).
    """
    w1, w2 = torus_periods(Lx, Ly, Tx, Ty)
    if omega1_override is not None: w1 = omega1_override
    if omega2_override is not None: w2 = omega2_override

    sites = lattice_sites(Lx, Ly, Tx, Ty)
    keys  = list(sites.keys())
    xs = np.array([sites[k][0] for k in keys], dtype=float)
    ys = np.array([sites[k][1] for k in keys], dtype=float)
    Gs = ising_torus_array(xs, ys, w1, w2)

    # Identify and remove the contact divergence at z = 0.  The Ising
    # spin correlator F(ζ) ~ |ζ|^{-1/4} blows up at the origin and that
    # one site otherwise dominates every cost.  Replace it by the value
    # at its *nearest* lattice neighbour (smallest |z|>0).
    rs = np.hypot(xs, ys)
    near0 = rs < 1e-9
    if near0.any():
        finite_mask = (~near0)
        if finite_mask.any():
            i_near = int(np.argmin(rs[finite_mask]))
            Gs[near0] = float(Gs[finite_mask][i_near])

    # Cap any remaining numerical outliers at the 99-th percentile so a
    # single near-singular site cannot dominate cost.
    cap = float(np.percentile(np.abs(Gs), 99.0))
    Gs = np.clip(Gs, -cap * 5.0, cap * 5.0)

    # The Ising correlator is positive; we keep it so (do NOT subtract a
    # G_far) so that log/effective-mass cost variants work natively.
    # Cost functions only see relative differences, so the absolute
    # offset has no effect on shape-comparison metrics.
    return {k: {"conn": float(g), "conn_err": err} for k, g in zip(keys, Gs)}

