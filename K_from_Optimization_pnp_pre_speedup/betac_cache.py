"""
betac_cache.py — Persistent β_c(r1, r2) interpolation cache (Speedup 2).

Each test geometry gets its own cache directory:

    results/_betac_cache_Lx{Lx}_Ly{Ly}_Tx{Tx}_Ty{Ty}/
        samples.jsonl    append-only log of every (r1, r2, β_c, σ) ever computed
        surface.npz      cached Delaunay + LinearND interpolant (rebuilt lazily)

Hit rule (a query (r1, r2) returns from the cache iff ALL hold):

    1. cache has at least `min_neighbours` samples
    2. point lies inside the convex hull of cached samples
    3. enclosing simplex is NOT a hull-boundary triangle (would extrapolate)
    4. local mesh size at the simplex vertices is ≤ tol_r
    5. precomputed leave-one-out residual at the nearest neighbour < tol_beta

Concurrent writers from multiple `run.py` processes are coordinated by a
POSIX `fcntl.flock` (Windows: `msvcrt.locking`) on a sentinel file in the
cache directory.  Reads stay lock-free because samples.jsonl is
append-only — partial last lines are tolerated and skipped.

Public API
----------
    BetacCache(geom, root, *, tol_r, tol_beta, min_neighbours, rebuild_every)
        .lookup(r1, r2) -> (beta_c, sigma) | None
        .add(r1, r2, beta_c, sigma, source_run="")
        .stats() -> dict
"""
from __future__ import annotations

import json
import math
import os
import time
from contextlib import contextmanager
from typing import Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Cross-platform append-write lock
# ---------------------------------------------------------------------------

@contextmanager
def _file_lock(path: str):
    """Exclusive lock around the cache directory; auto-releases on exit."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        if os.name == "nt":
            import msvcrt
            # Block until the lock is acquired (1 byte is enough).
            msvcrt.locking(fd, msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                try:
                    msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
                except OSError:
                    pass
        else:
            import fcntl
            fcntl.flock(fd, fcntl.LOCK_EX)
            try:
                yield
            finally:
                try:
                    fcntl.flock(fd, fcntl.LOCK_UN)
                except OSError:
                    pass
    finally:
        os.close(fd)


# ---------------------------------------------------------------------------
# Cache class
# ---------------------------------------------------------------------------

class BetacCache:
    """Persistent (r1, r2) → β_c interpolation cache for one test geometry."""

    def __init__(self, geom, root="results", *,
                 tol_r: float = 0.05,
                 tol_beta: float = 2e-3,
                 min_neighbours: int = 4,
                 rebuild_every: int = 8,
                 verbose: bool = True):
        Lx, Ly, Tx, Ty = (int(g) for g in geom)
        self.geom = (Lx, Ly, Tx, Ty)
        self.tol_r = float(tol_r)
        self.tol_beta = float(tol_beta)
        self.min_neighbours = int(min_neighbours)
        self.rebuild_every = max(1, int(rebuild_every))
        self.verbose = bool(verbose)

        tag = f"Lx{Lx}_Ly{Ly}_Tx{Tx}_Ty{Ty}"
        self.dir = os.path.join(root, f"_betac_cache_{tag}")
        os.makedirs(self.dir, exist_ok=True)
        self.samples_path = os.path.join(self.dir, "samples.jsonl")
        self.surface_path = os.path.join(self.dir, "surface.npz")
        self.lock_path    = os.path.join(self.dir, ".lock")

        # In-memory state.  Reload from disk; refresh on every lookup so
        # concurrent writers in other processes are picked up.
        self._samples: list = []
        self._samples_mtime: float = 0.0
        self._tri = None              # scipy.spatial.Delaunay
        self._interp = None           # scipy.interpolate.LinearNDInterpolator
        self._hull_simplex_mask = None  # bool array, True = touches hull
        self._loo_residuals = None    # per-sample LOO residual on β_c
        self._surface_n = 0           # how many samples were in the last build
        self._stats = {"lookups": 0, "hits": 0, "misses": 0, "adds": 0,
                       "rebuilds": 0}
        self._refresh_samples()

    # ----- disk I/O -----

    def _refresh_samples(self) -> None:
        """Reload samples.jsonl if it changed on disk (cheap stat call)."""
        if not os.path.exists(self.samples_path):
            self._samples = []
            self._samples_mtime = 0.0
            return
        try:
            mtime = os.path.getmtime(self.samples_path)
        except OSError:
            return
        if mtime == self._samples_mtime and self._samples:
            return
        out = []
        with open(self.samples_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    # Tolerate a truncated tail (concurrent writer mid-line).
                    continue
                if "r1" in rec and "r2" in rec and "beta_c" in rec:
                    out.append(rec)
        self._samples = out
        self._samples_mtime = mtime

    def _maybe_rebuild(self) -> None:
        """Rebuild the interpolant when ≥ rebuild_every new samples added."""
        n = len(self._samples)
        if n < max(self.min_neighbours, 4):
            return
        if self._interp is not None and (n - self._surface_n) < self.rebuild_every:
            return
        self._build_surface()

    def _build_surface(self) -> None:
        from scipy.spatial import Delaunay
        from scipy.interpolate import LinearNDInterpolator

        pts   = np.array([[s["r1"], s["r2"]] for s in self._samples], float)
        betas = np.array([s["beta_c"]        for s in self._samples], float)

        # Drop any duplicate (r1, r2) — Delaunay would fail.  Average β_c
        # at duplicates so we keep all the information.
        keys = [(round(p[0], 9), round(p[1], 9)) for p in pts]
        seen: dict = {}
        for k, b in zip(keys, betas):
            seen.setdefault(k, []).append(float(b))
        u_pts   = np.array([list(k) for k in seen.keys()], float)
        u_betas = np.array([float(np.mean(v)) for v in seen.values()], float)

        if len(u_pts) < max(self.min_neighbours, 4):
            self._tri = None
            self._interp = None
            return

        try:
            tri = Delaunay(u_pts)
        except Exception:
            # Co-linear samples or other degeneracy — bail; try again next time.
            self._tri = None
            self._interp = None
            return

        interp = LinearNDInterpolator(tri, u_betas)

        # Mark each simplex that has any vertex on the hull.  Hits served
        # from those simplices would have one anchor on the boundary of
        # the cached region and so can extrapolate badly.
        hull_vertices = np.unique(tri.convex_hull.ravel())
        on_hull = np.zeros(len(u_pts), dtype=bool)
        on_hull[hull_vertices] = True
        simplex_on_hull = on_hull[tri.simplices].any(axis=1)

        # Leave-one-out residuals — for each sample i, refit the interpolant
        # on the others and record |β_c(i) − β_interp(i)|.  Cheap on the
        # ~tens-of-points scale we run at; gates the tol_beta rule.
        loo = np.full(len(u_pts), np.inf, dtype=float)
        if len(u_pts) >= max(self.min_neighbours + 1, 5):
            for i in range(len(u_pts)):
                mask = np.ones(len(u_pts), dtype=bool)
                mask[i] = False
                try:
                    sub_tri = Delaunay(u_pts[mask])
                    sub_interp = LinearNDInterpolator(sub_tri, u_betas[mask])
                    pred = float(sub_interp(u_pts[i:i + 1])[0])
                    if math.isfinite(pred):
                        loo[i] = abs(pred - u_betas[i])
                except Exception:
                    pass

        self._tri = tri
        self._interp = interp
        self._hull_simplex_mask = simplex_on_hull
        self._loo_residuals = loo
        self._u_pts = u_pts
        self._u_betas = u_betas
        self._surface_n = len(self._samples)
        self._stats["rebuilds"] += 1
        if self.verbose:
            n_interior_simplices = int((~simplex_on_hull).sum())
            finite_loo = loo[np.isfinite(loo)]
            med = float(np.median(finite_loo)) if finite_loo.size else float("nan")
            print(f"[betac_cache] rebuilt surface "
                  f"(n={len(u_pts)} unique pts, {n_interior_simplices} "
                  f"interior simplices, median LOO={med:.2e})")

    # ----- public API -----

    def lookup(self, r1: float, r2: float) -> Optional[Tuple[float, float]]:
        """Return (β_c, σ_interp) if the query qualifies for a cache hit."""
        self._stats["lookups"] += 1
        self._refresh_samples()
        if len(self._samples) < self.min_neighbours:
            self._stats["misses"] += 1
            return None
        self._maybe_rebuild()
        if self._tri is None or self._interp is None:
            self._stats["misses"] += 1
            return None

        q = np.array([[float(r1), float(r2)]])

        # 2. inside convex hull?
        simplex = int(self._tri.find_simplex(q)[0])
        if simplex < 0:
            self._stats["misses"] += 1
            return None

        # 3. simplex must not touch the hull
        if bool(self._hull_simplex_mask[simplex]):
            self._stats["misses"] += 1
            return None

        # 4. all vertices of the enclosing simplex within tol_r in some
        #    metric — use the simplex circumradius bound (max edge length).
        verts_idx = self._tri.simplices[simplex]
        verts = self._u_pts[verts_idx]
        edges = [
            float(np.linalg.norm(verts[a] - verts[b]))
            for a in range(len(verts)) for b in range(a + 1, len(verts))
        ]
        if max(edges) > self.tol_r:
            self._stats["misses"] += 1
            return None

        # 5. nearest-neighbour LOO residual must be small
        d = np.linalg.norm(self._u_pts[verts_idx] - q, axis=1)
        nn = int(verts_idx[int(np.argmin(d))])
        loo = float(self._loo_residuals[nn])
        if not math.isfinite(loo) or loo > self.tol_beta:
            self._stats["misses"] += 1
            return None

        beta = float(self._interp(q)[0])
        if not math.isfinite(beta):
            self._stats["misses"] += 1
            return None

        # σ for the interpolated value: take the LOO residual as the
        # local error scale.  Always conservative because the live
        # neighbour was *included* in the actual interpolant.
        self._stats["hits"] += 1
        return (beta, loo)

    def add(self, r1: float, r2: float, beta_c: float, sigma: float = 0.0,
            source_run: str = "", n_traj_total: int = 0) -> None:
        """Append a new (r1, r2, β_c) sample to the cache."""
        rec = {
            "r1":          float(r1),
            "r2":          float(r2),
            "beta_c":      float(beta_c),
            "sigma":       float(sigma),
            "source_run":  str(source_run),
            "n_traj_total": int(n_traj_total),
            "ts":          time.time(),
        }
        # Keep the lock scope tiny: only the append is contended.
        with _file_lock(self.lock_path):
            with open(self.samples_path, "a") as f:
                f.write(json.dumps(rec) + "\n")
        self._stats["adds"] += 1
        # Force a refresh so subsequent lookups in this process see the
        # new sample without waiting for the next mtime tick.
        self._samples_mtime = 0.0

    def stats(self) -> dict:
        s = dict(self._stats)
        s["n_samples"]    = len(self._samples)
        s["surface_n"]    = self._surface_n
        s["hit_rate"]     = (s["hits"] / s["lookups"]) if s["lookups"] else 0.0
        return s

    def __len__(self) -> int:
        self._refresh_samples()
        return len(self._samples)
