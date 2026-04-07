#!/usr/bin/env python3
"""
betac_surface.py — Piecewise-linear FEM interpolant over a pre-computed β_c grid.

Usage as a module:
    from betac_surface import BetacSurface
    surf = BetacSurface("betac_surfaces/surface_32x32_T0x0.json")
    beta_c = surf(r1, r2)

Usage from the command line:
    python betac_surface.py surface.json --r1 1.0 --r2 1.0
    python betac_surface.py surface.json --plot
"""

import json
import sys
import warnings

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


class BetacSurface:
    """Piecewise-linear interpolant of β_c(r₁, r₂) from a JSON surface file."""

    def __init__(self, path):
        with open(path) as f:
            raw = json.load(f)

        self.meta = {
            "Lx": raw["Lx"], "Ly": raw["Ly"],
            "Tx": raw["Tx"], "Ty": raw["Ty"],
        }
        pts = raw["points"]
        self.r1_arr = np.array([p["r1"] for p in pts])
        self.r2_arr = np.array([p["r2"] for p in pts])
        self.beta_c_arr = np.array([p["beta_c"] for p in pts])
        self.points = np.column_stack([self.r1_arr, self.r2_arr])

        self._tri = Delaunay(self.points)
        self._interp = LinearNDInterpolator(self._tri, self.beta_c_arr)
        self._nn_vals = self.beta_c_arr  # for fallback

    def __call__(self, r1, r2):
        return self.interpolate(r1, r2)

    def interpolate(self, r1, r2):
        """Return interpolated β_c.  Falls back to nearest-neighbour outside hull."""
        val = float(self._interp(r1, r2))
        if np.isnan(val):
            dists = np.sqrt((self.r1_arr - r1)**2 + (self.r2_arr - r2)**2)
            idx = np.argmin(dists)
            val = float(self.beta_c_arr[idx])
            warnings.warn(
                f"({r1:.4f}, {r2:.4f}) is outside the surface domain; "
                f"using nearest-neighbour β_c={val:.8f}")
        return val

    def domain_contains(self, r1, r2):
        return self._tri.find_simplex(np.array([[r1, r2]]))[0] >= 0

    def summary(self):
        print(f"BetacSurface: {len(self.beta_c_arr)} points")
        print(f"  Lattice: {self.meta['Lx']}×{self.meta['Ly']}  "
              f"Tx={self.meta['Tx']} Ty={self.meta['Ty']}")
        print(f"  r1 ∈ [{self.r1_arr.min():.4f}, {self.r1_arr.max():.4f}]")
        print(f"  r2 ∈ [{self.r2_arr.min():.4f}, {self.r2_arr.max():.4f}]")
        print(f"  β_c ∈ [{self.beta_c_arr.min():.8f}, {self.beta_c_arr.max():.8f}]")

    def plot(self, output_path=None):
        """Generate a 2-panel visualisation of the surface."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.tri import Triangulation

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Left: scatter of data points
        sc = ax1.scatter(self.r1_arr, self.r2_arr, c=self.beta_c_arr,
                         cmap="viridis", s=40, edgecolors="k", linewidths=0.5)
        ax1.set_xlabel("r₁ = k₁/k₃")
        ax1.set_ylabel("r₂ = k₂/k₃")
        ax1.set_title("Data points")
        fig.colorbar(sc, ax=ax1, label="β_c")

        # Right: interpolated surface on fine grid
        r1_fine = np.linspace(self.r1_arr.min(), self.r1_arr.max(), 100)
        r2_fine = np.linspace(self.r2_arr.min(), self.r2_arr.max(), 100)
        R1, R2 = np.meshgrid(r1_fine, r2_fine)
        Bc = self._interp(R1.ravel(), R2.ravel()).reshape(R1.shape)
        im = ax2.pcolormesh(R1, R2, Bc, cmap="viridis", shading="auto")
        ax2.set_xlabel("r₁ = k₁/k₃")
        ax2.set_ylabel("r₂ = k₂/k₃")
        ax2.set_title("Interpolated surface")
        fig.colorbar(im, ax=ax2, label="β_c")

        fig.tight_layout()
        out = output_path or "betac_surface.png"
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
        plt.close(fig)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Inspect a β_c surface file")
    parser.add_argument("surface", help="Path to surface JSON file")
    parser.add_argument("--r1", type=float, default=None)
    parser.add_argument("--r2", type=float, default=None)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    surf = BetacSurface(args.surface)
    surf.summary()

    if args.r1 is not None and args.r2 is not None:
        bc = surf(args.r1, args.r2)
        inside = surf.domain_contains(args.r1, args.r2)
        print(f"\n  β_c({args.r1:.4f}, {args.r2:.4f}) = {bc:.10f}"
              f"  {'(inside hull)' if inside else '(OUTSIDE hull — nearest-neighbour)'}")

    if args.plot:
        surf.plot()


if __name__ == "__main__":
    main()
