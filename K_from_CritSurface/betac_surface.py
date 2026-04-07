#!/usr/bin/env python3
"""
betac_surface.py — FEM interpolant for the pre-computed β_c(r₁, r₂) surface
============================================================================

Loads a β_c surface JSON file produced by build_betac_surface.py and provides
a smooth bilinear FEM interpolant over the (r₁, r₂) domain.

The surface data consists of β_c values measured by MC susceptibility scans
at discrete (r₁, r₂) grid points.  This module triangulates the grid and
builds a piecewise-linear (P1 finite element) interpolant, so that any
query point (r₁, r₂) inside the convex hull of the data returns a smooth,
deterministic β_c estimate with no MC jitter.

Usage:
    from betac_surface import BetacSurface

    surf = BetacSurface("betac_surfaces/surface_32x32_T0x0.json")
    beta_c = surf(1.05, 0.97)       # interpolate at a point
    beta_c = surf.interpolate(r1, r2)  # same thing
"""

import json
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


class BetacSurface:
    """
    Piecewise-linear (P1 FEM) interpolant of β_c(r₁, r₂).

    Built from a JSON surface file containing scattered (r₁, r₂, β_c) data
    produced by build_betac_surface.py.

    Parameters
    ----------
    path : str
        Path to the surface JSON file.

    Attributes
    ----------
    meta : dict
        Lattice metadata (Lx, Ly, Tx, Ty) from the surface file.
    r1_vals, r2_vals : np.ndarray
        The 1D grid coordinates used to build the surface.
    points : np.ndarray, shape (N, 2)
        All (r₁, r₂) data points.
    beta_c_vals : np.ndarray, shape (N,)
        Corresponding β_c values.
    """

    def __init__(self, path):
        with open(path) as f:
            data = json.load(f)

        self.meta = {
            "Lx": data["Lx"], "Ly": data["Ly"],
            "Tx": data["Tx"], "Ty": data["Ty"],
        }
        self.r1_vals = np.array(data["r1_vals"])
        self.r2_vals = np.array(data["r2_vals"])

        pts = data["points"]
        self.points = np.array([[p["r1"], p["r2"]] for p in pts])
        self.beta_c_vals = np.array([p["beta_c"] for p in pts])

        # Build the Delaunay triangulation and P1 interpolant
        self._tri = Delaunay(self.points)
        self._interp = LinearNDInterpolator(self._tri, self.beta_c_vals)

        # Domain bounds (for extrapolation warnings)
        self._r1_lo = self.r1_vals.min()
        self._r1_hi = self.r1_vals.max()
        self._r2_lo = self.r2_vals.min()
        self._r2_hi = self.r2_vals.max()

    def interpolate(self, r1, r2):
        """
        Interpolate β_c at (r₁, r₂).

        Returns the piecewise-linear interpolated value.  If (r₁, r₂) is
        outside the convex hull of the surface data, falls back to
        nearest-neighbor extrapolation with a warning.

        Parameters
        ----------
        r1, r2 : float
            Coupling ratios.

        Returns
        -------
        float
            Interpolated β_c.
        """
        val = self._interp(r1, r2)
        if np.isnan(val):
            # Outside convex hull — nearest-neighbor fallback
            dists = np.sum((self.points - [r1, r2]) ** 2, axis=1)
            idx = np.argmin(dists)
            import warnings
            warnings.warn(
                f"BetacSurface: ({r1:.4f}, {r2:.4f}) is outside the surface domain "
                f"[{self._r1_lo:.3f}–{self._r1_hi:.3f}] × [{self._r2_lo:.3f}–{self._r2_hi:.3f}]; "
                f"using nearest neighbor β_c={self.beta_c_vals[idx]:.8f}"
            )
            return float(self.beta_c_vals[idx])
        return float(val)

    def __call__(self, r1, r2):
        """Shorthand for interpolate(r1, r2)."""
        return self.interpolate(r1, r2)

    def domain_contains(self, r1, r2):
        """Check whether (r₁, r₂) is inside the surface convex hull."""
        return self._tri.find_simplex([r1, r2]) >= 0

    def summary(self):
        """Print a summary of the surface."""
        print(f"β_c surface: {self.meta['Lx']}×{self.meta['Ly']} "
              f"Tx={self.meta['Tx']} Ty={self.meta['Ty']}")
        print(f"  r₁ ∈ [{self._r1_lo:.4f}, {self._r1_hi:.4f}]  "
              f"({len(self.r1_vals)} pts)")
        print(f"  r₂ ∈ [{self._r2_lo:.4f}, {self._r2_hi:.4f}]  "
              f"({len(self.r2_vals)} pts)")
        print(f"  {len(self.beta_c_vals)} total data points, "
              f"{len(self._tri.simplices)} triangles")
        print(f"  β_c range: [{self.beta_c_vals.min():.8f}, "
              f"{self.beta_c_vals.max():.8f}]")


# ---------------------------------------------------------------------------
# CLI: inspect a surface file
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inspect a β_c surface file")
    parser.add_argument("surface", help="Path to surface JSON file")
    parser.add_argument("--r1", type=float, default=None,
                        help="Query r₁ (interpolate at this point)")
    parser.add_argument("--r2", type=float, default=None,
                        help="Query r₂")
    parser.add_argument("--plot", action="store_true",
                        help="Plot the surface as a heatmap")
    args = parser.parse_args()

    surf = BetacSurface(args.surface)
    surf.summary()

    if args.r1 is not None and args.r2 is not None:
        bc = surf(args.r1, args.r2)
        inside = surf.domain_contains(args.r1, args.r2)
        print(f"\n  β_c({args.r1:.4f}, {args.r2:.4f}) = {bc:.10f}"
              f"  {'(interpolated)' if inside else '(extrapolated/nearest)'}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.tri import Triangulation

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: scatter colored by β_c
        sc = axes[0].scatter(surf.points[:, 0], surf.points[:, 1],
                             c=surf.beta_c_vals, cmap="viridis", s=40)
        axes[0].set_xlabel("r₁ = k₁/k₃")
        axes[0].set_ylabel("r₂ = k₂/k₃")
        axes[0].set_title("β_c scan data")
        fig.colorbar(sc, ax=axes[0], label="β_c")

        # Right: interpolated surface on a fine grid
        r1_fine = np.linspace(surf._r1_lo, surf._r1_hi, 100)
        r2_fine = np.linspace(surf._r2_lo, surf._r2_hi, 100)
        R1, R2 = np.meshgrid(r1_fine, r2_fine)
        Bc = surf._interp(R1.ravel(), R2.ravel()).reshape(R1.shape)
        im = axes[1].pcolormesh(R1, R2, Bc, cmap="viridis", shading="auto")
        axes[1].set_xlabel("r₁ = k₁/k₃")
        axes[1].set_ylabel("r₂ = k₂/k₃")
        axes[1].set_title("β_c surface (P1 FEM interpolant)")
        fig.colorbar(im, ax=axes[1], label="β_c")

        plt.tight_layout()
        out = args.surface.replace(".json", "_surface.png")
        plt.savefig(out, dpi=150)
        print(f"\n  Plot saved: {out}")
