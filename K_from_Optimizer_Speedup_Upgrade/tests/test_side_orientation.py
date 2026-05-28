from __future__ import annotations

import os
import sys

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
if PKG not in sys.path:
    sys.path.insert(0, PKG)

import cost


def test_l2_cost_pairs_cycles_by_physical_side(monkeypatch):
    ref_geom = (13, 16, 3, -3)
    test_geom = (16, 16, 0, 0)

    role_level = {"short": 1.0, "mid": 2.0, "long": 3.0}

    def fake_tile_interp(_data, Lx, Ly, Tx, Ty, field, copies=2):
        geom = (Lx, Ly, Tx, Ty)
        paths = cost.boundary_paths(*geom)
        roles = cost.cycle_roles(*geom)
        xy_dirs = [np.array(cost._triangular_xy(dm, dn), dtype=float)
                   for dm, dn in paths]

        def interp(pts):
            pts = np.asarray(pts, dtype=float)
            if field == "conn_err":
                return np.full(len(pts), 1.0e-3)
            direction = pts[-1]
            norms = [abs(vec[0] * direction[1] - vec[1] * direction[0]) /
                     max(np.dot(vec, vec), 1.0e-12)
                     for vec in xy_dirs]
            cycle_idx = int(np.argmin(norms))
            vec = xy_dirs[cycle_idx]
            scale = max(np.dot(vec, vec), 1.0e-12)
            t = np.clip((pts @ vec) / scale, 0.0, 1.0)
            base = role_level[roles[cycle_idx]]
            return base + 0.25 * t

        return interp

    monkeypatch.setattr(cost, "_tile_interp", fake_tile_interp)

    matched_cost, _, _, _ = cost.l2_cost(
        {}, {},
        *test_geom,
        *ref_geom,
        n_samples=32,
        copies=1,
    )
    raw_index_cost, _, _, _ = cost.l2_cost(
        {}, {},
        *test_geom,
        *ref_geom,
        n_samples=32,
        copies=1,
        ref_cycle_roles=(0, 1, 2),
        test_cycle_roles=(0, 1, 2),
    )

    assert matched_cost < 1.0e-12
    assert raw_index_cost > 0.5