"""Tests for betac_cache.BetacCache (pure Python, no MC)."""
from __future__ import annotations

import json
import os
import threading

import numpy as np
import pytest

from betac_cache import BetacCache


# Smooth analytic β_c surface used to seed synthetic cache content.
def _truth(r1, r2):
    return 0.27 + 0.03 * r1 - 0.02 * r2 + 0.005 * (r1 - 1.0) * (r2 - 1.0)


@pytest.fixture
def empty_cache(tmp_path):
    return BetacCache((24, 24, 0, 0), root=str(tmp_path),
                      verbose=False)


@pytest.fixture
def grid_cache(tmp_path):
    cache = BetacCache((24, 24, 0, 0), root=str(tmp_path),
                       tol_r=0.30, tol_beta=1e-2,
                       min_neighbours=4, verbose=False)
    # 7×7 grid in [0.5, 1.5]^2
    for r1 in np.linspace(0.5, 1.5, 7):
        for r2 in np.linspace(0.5, 1.5, 7):
            cache.add(float(r1), float(r2), _truth(r1, r2), sigma=0.0)
    return cache


def test_empty_cache_misses(empty_cache):
    assert empty_cache.lookup(1.0, 1.0) is None
    assert len(empty_cache) == 0


def test_lookup_below_min_neighbours(tmp_path):
    cache = BetacCache((8, 8, 0, 0), root=str(tmp_path),
                       min_neighbours=10, verbose=False)
    for r1 in np.linspace(0.5, 1.5, 5):
        cache.add(float(r1), 1.0, _truth(r1, 1.0))
    assert cache.lookup(1.0, 1.0) is None  # < min_neighbours


def test_interior_hit_accuracy(grid_cache):
    # Several interior queries should all return values close to truth.
    queries = [(1.0, 1.0), (0.83, 1.17), (1.21, 0.79), (0.95, 1.05)]
    n_hits = 0
    for r1, r2 in queries:
        hit = grid_cache.lookup(r1, r2)
        if hit is not None:
            beta_pred, _sigma = hit
            assert abs(beta_pred - _truth(r1, r2)) < 5e-3
            n_hits += 1
    assert n_hits >= 3, f"expected ≥3 hits among {queries}, got {n_hits}"


def test_lookup_outside_hull_misses(grid_cache):
    assert grid_cache.lookup(2.5, 1.0) is None
    assert grid_cache.lookup(0.0, 0.0) is None


def test_hull_facet_rejected(grid_cache):
    # Points just inside the hull boundary should fall on a hull-touching
    # simplex and be rejected (no extrapolation risk allowed).
    rejected = 0
    for r1 in np.linspace(0.51, 0.55, 4):
        for r2 in np.linspace(0.51, 0.55, 4):
            if grid_cache.lookup(float(r1), float(r2)) is None:
                rejected += 1
    assert rejected >= 8  # majority of corner-region queries should be rejected


def test_persistence_across_instances(tmp_path):
    geom = (16, 16, 0, 0)
    a = BetacCache(geom, root=str(tmp_path), verbose=False)
    a.add(1.0, 1.0, 0.27)
    a.add(0.9, 1.0, 0.30)
    b = BetacCache(geom, root=str(tmp_path), verbose=False)
    assert len(b) == 2


def test_cross_geometry_isolation(tmp_path):
    a = BetacCache((24, 24, 0, 0), root=str(tmp_path), verbose=False)
    b = BetacCache((16, 16, 0, 0), root=str(tmp_path), verbose=False)
    a.add(1.0, 1.0, 0.27)
    assert len(a) == 1
    assert len(b) == 0
    assert a.dir != b.dir


def test_concurrent_add_lock(tmp_path):
    cache = BetacCache((8, 8, 0, 0), root=str(tmp_path), verbose=False)
    n_per = 25

    def worker(offset):
        for i in range(n_per):
            cache.add(1.0 + offset + 0.001 * i, 1.0,
                      0.27 + 0.001 * (offset + i))

    threads = [threading.Thread(target=worker, args=(0.1 * t,)) for t in range(4)]
    for t in threads: t.start()
    for t in threads: t.join()

    # Re-read fresh — total lines should be exactly 4 * n_per with no
    # truncated/garbled entries.
    fresh = BetacCache((8, 8, 0, 0), root=str(tmp_path), verbose=False)
    assert len(fresh) == 4 * n_per
    with open(cache.samples_path) as f:
        for line in f:
            json.loads(line.strip())  # raises if any line is malformed


def test_rebuild_cadence(tmp_path):
    cache = BetacCache((8, 8, 0, 0), root=str(tmp_path),
                       min_neighbours=4, rebuild_every=3, verbose=False)
    for r1 in np.linspace(0.5, 1.5, 6):
        for r2 in np.linspace(0.5, 1.5, 6):
            cache.add(float(r1), float(r2), _truth(r1, r2))
    # Force one rebuild then check stat counter increments.
    cache.lookup(1.0, 1.0)
    n0 = cache.stats()["rebuilds"]
    # No new samples → no new rebuild.
    cache.lookup(1.0, 1.0)
    assert cache.stats()["rebuilds"] == n0
    # Add < rebuild_every new samples → still no rebuild.
    cache.add(1.0, 1.01, _truth(1.0, 1.01))
    cache.lookup(1.0, 1.0)
    assert cache.stats()["rebuilds"] == n0
    # Add enough new samples to trigger rebuild.
    for k in range(3):
        cache.add(1.0 + 0.001 * k, 1.02, _truth(1.0, 1.02))
    cache.lookup(1.0, 1.0)
    assert cache.stats()["rebuilds"] == n0 + 1


def test_stats_hit_rate(grid_cache):
    grid_cache.lookup(1.0, 1.0)        # likely hit
    grid_cache.lookup(2.5, 1.0)        # miss
    s = grid_cache.stats()
    assert s["lookups"] == 2
    assert s["hits"] + s["misses"] == 2
    assert 0.0 <= s["hit_rate"] <= 1.0
