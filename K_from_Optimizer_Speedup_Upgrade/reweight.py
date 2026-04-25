"""
reweight.py — Single-histogram Ferrenberg–Swendsen reweighting (Speedup 4).

Given one MC ensemble at parent inverse temperature β₀ with per-trajectory
samples of (E_total, |m|, m²), construct an unbiased estimate of
χ(β) = V * (⟨m²⟩_β − ⟨|m|⟩_β²) on a dense β grid via canonical
reweighting:

    w_i(β) = exp( −(β − β₀) · E_total,i ) / Z(β)

with Z(β) = Σ_i exp(−(β−β₀) · E_i).  Effective sample size
N_eff(β) = (Σ wᵢ)² / Σ wᵢ² is reported alongside every estimate so the
caller can detect when β has drifted out of the validity window.

Public API
----------
    Reweighter(traces_path, n_sites)
        .chi(beta)              → (chi, sigma, n_eff)
        .chi_grid(beta_grid)    → (chis, sigmas, n_effs)
        .energy_per_site_mean   property
        .beta_parent            property
        .n_samples              property
        .valid_window(threshold=0.1)  → (β_lo, β_hi)

A ``CombinedReweighter`` builds on multiple parent histograms (multi-
histogram / WHAM extension); it picks the parent with the largest
N_eff at each query β and falls back to the next-best when below the
N_eff floor.  Suitable for stitching reweighted estimates across the
3-pass scan when one parent histogram doesn't span the whole bracket.

Notes
-----
* The error bar on χ(β) is computed from a *block jackknife* on the
  parent histogram (default 16 blocks).  Per-β χ values returned at
  different β share the same parent and are therefore highly
  correlated; downstream fits (e.g. the GC β_c finder) MUST treat
  them as such — see ``find_beta_c_reweight`` in mc_engine.py.
* Pass ``per_site_energy=True`` when the trace stores energy already
  divided by V (this is what ``traces.dat`` writes); the reweighter
  multiplies internally by ``n_sites`` to recover the extensive
  weight-argument E_total = E_per_site · V.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Trace file I/O
# ---------------------------------------------------------------------------

def parse_traces(path: str) -> dict:
    """Read a traces.dat file written by --dump_traces 1.

    Returns a dict::

        {"E_per_site": np.ndarray,  # (N,)
         "abs_m":      np.ndarray,  # (N,)
         "m2":         np.ndarray,  # (N,)
         "n_sites":    int,
         "beta_parent":float,
         "K":          (K1, K2, K3, Kt),
         "path":       str}
    """
    n_sites = None
    beta = None
    K = (None, None, None, None)
    rows: list = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Look for the metadata header; tolerate other comment lines.
                if "n_sites" in line and "beta" in line:
                    parts = line.split()
                    # tokens like ["#", "n_sites", "36", "beta", "0.27", "K1", "1.0", ...]
                    try:
                        n_sites = int(parts[parts.index("n_sites") + 1])
                        beta = float(parts[parts.index("beta") + 1])
                        K = tuple(
                            float(parts[parts.index(name) + 1])
                            for name in ("K1", "K2", "K3", "Kt")
                        )
                    except (ValueError, IndexError):
                        pass
                continue
            cols = line.split()
            if len(cols) < 3:
                continue
            try:
                rows.append([float(cols[0]), float(cols[1]), float(cols[2])])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"no trajectory rows in {path}")
    if n_sites is None or beta is None:
        raise ValueError(
            f"traces file {path} missing required header (n_sites/beta); "
            "rebuild simulator binary if older than the --dump_traces patch."
        )
    arr = np.asarray(rows, dtype=float)
    return {
        "E_per_site": arr[:, 0],
        "abs_m":      arr[:, 1],
        "m2":         arr[:, 2],
        "n_sites":    int(n_sites),
        "beta_parent": float(beta),
        "K":          K,
        "path":       str(path),
    }


# ---------------------------------------------------------------------------
# Single-histogram reweighter
# ---------------------------------------------------------------------------

@dataclass
class _ChiEstimate:
    chi: float
    sigma: float
    n_eff: float


class Reweighter:
    """Single-histogram Ferrenberg–Swendsen reweighting.

    Parameters
    ----------
    traces_path : str | dict
        Path to a traces.dat file (or a dict already returned by
        :func:`parse_traces`).
    n_blocks : int, default 16
        Number of contiguous blocks used for the jackknife error bar
        on χ(β).  Set higher for smoother error curves at the cost of
        slightly larger runtime per query (negligible vs MC).
    n_eff_floor : float, default 0.1
        Reject χ(β) when ``n_eff(β) / n_samples < n_eff_floor``.
        Returns NaN for both ``chi`` and ``sigma``.  Use
        ``valid_window()`` to discover the usable Δβ from this floor.
    """

    def __init__(self, traces_path, *, n_blocks: int = 16,
                 n_eff_floor: float = 0.1):
        if isinstance(traces_path, dict):
            t = traces_path
        else:
            t = parse_traces(str(traces_path))
        self._E_total = t["E_per_site"] * t["n_sites"]   # extensive
        self._abs_m = t["abs_m"]
        self._m2 = t["m2"]
        self._n_sites = int(t["n_sites"])
        self._beta = float(t["beta_parent"])
        self._K = tuple(t["K"])
        self._n = len(self._E_total)
        if self._n < 8:
            raise ValueError(
                f"reweighter needs ≥ 8 samples, got {self._n} from {t['path']}"
            )
        self._n_blocks = max(2, min(int(n_blocks), self._n // 4))
        self._n_eff_floor = float(n_eff_floor)
        # Pre-shift the energy so weights stay numerically tame even at
        # large |Δβ| · E.  The shift cancels in all expectation values.
        self._E_shifted = self._E_total - float(np.mean(self._E_total))

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def beta_parent(self) -> float:
        return self._beta

    @property
    def n_sites(self) -> int:
        return self._n_sites

    @property
    def n_samples(self) -> int:
        return self._n

    @property
    def couplings(self) -> Tuple[float, float, float, float]:
        return self._K

    @property
    def energy_per_site_mean(self) -> float:
        return float(np.mean(self._E_total) / self._n_sites)

    # ------------------------------------------------------------------
    # Core reweighting kernel
    # ------------------------------------------------------------------

    def _weights(self, beta: float) -> np.ndarray:
        """Normalised weights w_i(β); shape (N,), sum = 1."""
        dbeta = float(beta) - self._beta
        # log w_i ∝ −Δβ · E_i  (extensive); subtract max for stability.
        log_w = -dbeta * self._E_shifted
        log_w -= log_w.max()
        w = np.exp(log_w)
        Z = w.sum()
        if Z <= 0.0 or not np.isfinite(Z):
            return np.zeros_like(w)
        return w / Z

    def _chi_from_weights(self, w: np.ndarray) -> float:
        """χ = ⟨m²⟩ − ⟨|m|⟩² at given weights.

        Matches the C++ simulator's ``m_susc`` convention (per-site,
        i.e. NO V factor).  Multiply by ``self.n_sites`` if you want
        the extensive susceptibility.
        """
        if w.sum() <= 0.0:
            return float("nan")
        m2_mean = float((w * self._m2).sum())
        am_mean = float((w * self._abs_m).sum())
        return m2_mean - am_mean * am_mean

    @staticmethod
    def _n_eff(w: np.ndarray) -> float:
        s = w.sum()
        if s <= 0.0:
            return 0.0
        return float(s * s / max(np.dot(w, w), 1e-300))

    # ------------------------------------------------------------------
    # Public estimators
    # ------------------------------------------------------------------

    def chi(self, beta: float) -> _ChiEstimate:
        """Return a ``_ChiEstimate(chi, sigma, n_eff)`` at ``beta``."""
        w = self._weights(float(beta))
        n_eff = self._n_eff(w * self._n)   # un-normalised N_eff
        if n_eff < self._n_eff_floor * self._n:
            return _ChiEstimate(float("nan"), float("nan"), n_eff)
        chi_hat = self._chi_from_weights(w)
        sigma = self._jackknife_sigma(float(beta))
        return _ChiEstimate(chi_hat, sigma, n_eff)

    def chi_grid(self, betas: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Vectorised χ(β) on a 1-D grid.

        Returns ``(chis, sigmas, n_effs)`` of shape ``(len(betas),)``.
        Uses one combined matrix exp() instead of looping over betas
        for the central estimate; jackknife errors still loop over
        blocks.
        """
        betas = np.asarray(betas, dtype=float)
        dbetas = betas - self._beta
        # log_w[k, i] = -dbetas[k] * E_shifted[i]
        log_w = -np.outer(dbetas, self._E_shifted)
        log_w -= log_w.max(axis=1, keepdims=True)
        w = np.exp(log_w)
        Z = w.sum(axis=1, keepdims=True)
        # Normalised weights wn[k, i]
        wn = w / np.where(Z > 0, Z, 1.0)
        m2_mean = wn @ self._m2
        am_mean = wn @ self._abs_m
        chis = m2_mean - am_mean * am_mean
        n_effs = (w.sum(axis=1) ** 2) / np.maximum((w * w).sum(axis=1), 1e-300)
        # NaN-out points below the N_eff floor.
        mask = n_effs < self._n_eff_floor * self._n
        chis = chis.copy()
        chis[mask] = np.nan
        # Per-point sigma via the jackknife loop (cheap: <ms per beta).
        sigmas = np.array([self._jackknife_sigma(float(b)) for b in betas])
        sigmas[mask] = np.nan
        return chis, sigmas, n_effs

    # ------------------------------------------------------------------
    # Jackknife error
    # ------------------------------------------------------------------

    def _jackknife_sigma(self, beta: float) -> float:
        """Block-jackknife σ on χ(β) over the parent histogram.

        The parent samples are partitioned into ``n_blocks`` contiguous
        blocks; for each block we re-reweight on the remaining N − N_b
        samples and record χ_b.  The standard formula
            σ² = (B−1)/B · Σ (χ_b − ⟨χ⟩)²
        gives the error bar.
        """
        B = self._n_blocks
        if B < 2 or self._n < 2 * B:
            return 0.0
        block_size = self._n // B
        chis = np.empty(B, dtype=float)
        for b in range(B):
            lo = b * block_size
            hi = (b + 1) * block_size if b < B - 1 else self._n
            mask = np.ones(self._n, dtype=bool)
            mask[lo:hi] = False
            E_sub = self._E_shifted[mask]
            log_w = -(float(beta) - self._beta) * E_sub
            log_w -= log_w.max()
            w = np.exp(log_w)
            Z = w.sum()
            if Z <= 0.0 or not np.isfinite(Z):
                chis[b] = np.nan
                continue
            wn = w / Z
            m2_mean = float((wn * self._m2[mask]).sum())
            am_mean = float((wn * self._abs_m[mask]).sum())
            chis[b] = m2_mean - am_mean * am_mean
        chis = chis[np.isfinite(chis)]
        if len(chis) < 2:
            return 0.0
        var = (len(chis) - 1) / len(chis) * np.sum(
            (chis - chis.mean()) ** 2
        )
        return float(np.sqrt(var))

    def valid_window(self, threshold: Optional[float] = None,
                     n_probe: int = 41,
                     max_dbeta: Optional[float] = None) -> Tuple[float, float]:
        """Estimate the (β_lo, β_hi) window where N_eff stays above the floor.

        Probes ``n_probe`` βs spaced symmetrically around β_parent and
        returns the largest contiguous interval containing β_parent on
        which the N_eff fraction stays above ``threshold`` (defaults to
        ``self._n_eff_floor``).  ``max_dbeta`` defaults to a heuristic
        based on σ_E / √N — the natural scale of the Boltzmann weight
        spread.
        """
        thr = float(threshold) if threshold is not None else self._n_eff_floor
        if max_dbeta is None:
            sigma_E = float(np.std(self._E_total))
            max_dbeta = max(0.005, 5.0 * sigma_E / max(np.sqrt(self._n), 1.0))
        bs = np.linspace(self._beta - max_dbeta, self._beta + max_dbeta, n_probe)
        _, _, n_effs = self.chi_grid(bs)
        ok = (n_effs / self._n) >= thr
        if not ok.any():
            return (self._beta, self._beta)
        # Find the contiguous "True" run that contains the centre.
        centre = n_probe // 2
        if not ok[centre]:
            return (self._beta, self._beta)
        lo = centre
        while lo > 0 and ok[lo - 1]:
            lo -= 1
        hi = centre
        while hi < n_probe - 1 and ok[hi + 1]:
            hi += 1
        return (float(bs[lo]), float(bs[hi]))


# ---------------------------------------------------------------------------
# Multi-histogram (WHAM-lite) combiner
# ---------------------------------------------------------------------------

class CombinedReweighter:
    """Pick the most informative parent at each β; fallback to next-best.

    Not full WHAM (no iterative free-energy solve) — for our use case
    (3-pass scan around the susceptibility peak) the parent histograms
    overlap heavily and the simple "max-N_eff parent wins" rule gives
    almost the same answer as full WHAM at a fraction of the bookkeeping.
    """

    def __init__(self, reweighters: Iterable[Reweighter], *,
                 n_eff_floor: float = 0.1):
        self._rws = list(reweighters)
        if not self._rws:
            raise ValueError("CombinedReweighter needs ≥ 1 parent reweighter")
        self._n_eff_floor = float(n_eff_floor)

    def chi(self, beta: float) -> _ChiEstimate:
        # For each parent: compute weights and N_eff at this β, pick the
        # one with the largest N_eff that also passes the floor.
        best = None
        for rw in self._rws:
            w = rw._weights(float(beta))
            n_eff = Reweighter._n_eff(w * rw._n)
            if n_eff < self._n_eff_floor * rw._n:
                continue
            if best is None or n_eff > best[1]:
                best = (rw, n_eff, w)
        if best is None:
            return _ChiEstimate(float("nan"), float("nan"), 0.0)
        rw, n_eff, w = best
        chi_hat = rw._chi_from_weights(w)
        sigma = rw._jackknife_sigma(float(beta))
        return _ChiEstimate(chi_hat, sigma, n_eff)

    def chi_grid(self, betas: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        out_chi = np.full(len(betas), np.nan)
        out_sig = np.full(len(betas), np.nan)
        out_neff = np.zeros(len(betas))
        for i, b in enumerate(betas):
            est = self.chi(float(b))
            out_chi[i] = est.chi
            out_sig[i] = est.sigma
            out_neff[i] = est.n_eff
        return out_chi, out_sig, out_neff

    def add(self, rw: Reweighter) -> None:
        self._rws.append(rw)

    def __len__(self) -> int:
        return len(self._rws)


# ---------------------------------------------------------------------------
# Convenience: locate the traces file produced by run_simulator
# ---------------------------------------------------------------------------

def find_traces_file(subdir: str) -> Optional[str]:
    """Return the path to the lone traces_*.dat in *subdir*, or None.

    The simulator writes exactly one such file per run.
    """
    if not os.path.isdir(subdir):
        return None
    for name in sorted(os.listdir(subdir)):
        if name.startswith("traces_") and name.endswith(".dat"):
            return os.path.join(subdir, name)
    return None
