#!/usr/bin/env python3
"""
Test of the boundary-slice concept for the K_from_TwoPoint optimizer.

For a twisted parallelogram with basis vectors
  v = (Lx, Ty),  u = (Tx, -Ly)
the three *boundary* paths on the torus are along v, u, and w = -(v+u).
At the critical point these three 1-D correlators should be identical
(isotropic) for an equilateral geometry.

Approach
--------
1. Load the all-to-all dataset.
2. Tile 3x3 periodic copies (shifting by multiples of v and u in (m,n)).
3. Convert (m,n) -> Cartesian (x,y) using triangular-lattice embedding,
   then build a 2D scipy linear interpolant.
4. Sample 500 points along each of the three boundary paths.
5. Plot all three on the same axes (should collapse to one curve for
   an isotropic lattice at critical coupling).

Usage:
  python plot_boundary_slices_test.py        # uses default 4-5-6 scale-3 data
  python plot_boundary_slices_test.py --dat PATH --Lx N --Ly N --Tx N --Ty N
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def mn_to_xy(m, n):
    """Triangular lattice: e1=(1,0), e2=(1/2, sqrt(3)/2)."""
    m = np.asarray(m, dtype=float)
    n = np.asarray(n, dtype=float)
    x = m + 0.5 * n
    y = (np.sqrt(3.0) / 2.0) * n
    return x, y


def load_a2a(path: Path):
    """Load two_point_all_to_all.dat -> arrays (m, n, corr, err, conn, conn_err)."""
    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 7:
                continue
            rows.append((int(parts[1]), int(parts[2]),
                         float(parts[3]), float(parts[4]),
                         float(parts[5]), float(parts[6])))
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4], arr[:, 5]


def tile_periodic(m, n, cc, ec, Lx, Ly, Tx, Ty, copies=1):
    """
    Tile data by [-copies..+copies] shifts along v=(Lx,Ty) and u=(Tx,-Ly).
    Removes exact-duplicate (x,y) positions after conversion to Cartesian.
    Returns (x_tiled, y_tiled, cc_tiled, ec_tiled).
    """
    m_list, n_list, cc_list, ec_list = [], [], [], []
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            dm = a * Lx + b * Tx
            dn = a * Ty - b * Ly
            m_list.append(m + dm)
            n_list.append(n + dn)
            cc_list.append(cc)
            ec_list.append(ec)

    m_all = np.concatenate(m_list)
    n_all = np.concatenate(n_list)
    x_all, y_all = mn_to_xy(m_all, n_all)
    cc_all = np.concatenate(cc_list)
    ec_all = np.concatenate(ec_list)

    # Remove duplicates (can arise from degenerate lattice geometry)
    pts = np.column_stack([x_all, y_all])
    _, unique_idx = np.unique(np.round(pts, 8), axis=0, return_index=True)
    return x_all[unique_idx], y_all[unique_idx], cc_all[unique_idx], ec_all[unique_idx]


def make_interp(x, y, z):
    """Build a 2D linear interpolant from scattered (x,y) -> z."""
    pts = np.column_stack([x, y])
    return LinearNDInterpolator(pts, z)


def sample_path(interp_cc, interp_ec, dm_end, dn_end, n_samples=500):
    """
    Sample the interpolant from (0,0) to (dm_end, dn_end) in (m,n) space.
    Returns (t, cc_vals, ec_vals) where t is in [0,1].
    """
    t = np.linspace(0.0, 1.0, n_samples)
    m_path = t * dm_end
    n_path = t * dn_end
    x_path, y_path = mn_to_xy(m_path, n_path)
    pts = np.column_stack([x_path, y_path])
    cc_vals = interp_cc(pts)
    ec_vals = interp_ec(pts)
    mask = np.isfinite(cc_vals) & np.isfinite(ec_vals)
    return t, cc_vals, ec_vals, mask


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", default=None,
                    help="Path to two_point_all_to_all.dat")
    ap.add_argument("--Lx", type=int, default=39)
    ap.add_argument("--Ly", type=int, default=48)
    ap.add_argument("--Tx", type=int, default=9)
    ap.add_argument("--Ty", type=int, default=-9)
    ap.add_argument("--samples", type=int, default=500)
    ap.add_argument("--copies", type=int, default=2,
                    help="Number of periodic tile copies in each direction (default 2)")
    ap.add_argument("--output", default="/tmp/boundary_slices_456.png")
    args = ap.parse_args()

    if args.dat is None:
        base = Path(__file__).parent
        args.dat = str(base / "out_456_highprec_scale3" /
                       "Lx39_Ly48_Tx9_Ty-9_k0.275_0.275_0.275" /
                       "two_point_all_to_all.dat")

    Lx, Ly, Tx, Ty = args.Lx, args.Ly, args.Tx, args.Ty
    Ncell = Lx * Ly + Tx * Ty

    print(f"Geometry: Lx={Lx}, Ly={Ly}, Tx={Tx}, Ty={Ty}  (Ncell={Ncell})")
    print(f"Loading: {args.dat}")

    m, n, corr, err, conn, conn_err = load_a2a(Path(args.dat))
    print(f"  {len(m)} displacement vectors loaded")

    # Three boundary directions
    # v = (Lx, Ty), u = (Tx, -Ly), w = -(v+u)
    v = (Lx, Ty)
    u = (Tx, -Ly)
    w = (-(Lx + Tx), Ly - Ty)

    print(f"\nBoundary vectors:")
    print(f"  v = {v}   (Lx steps in e1 + Ty steps in e2)")
    print(f"  u = {u}   (Tx steps in e1 + (-Ly) steps in e2)")
    print(f"  w = {w}   = -(v+u), closes the triangle")

    # Build interpolant from tiled data
    print(f"\nTiling data ({2*args.copies+1}x{2*args.copies+1} copies)...")
    x_t, y_t, cc_t, ec_t = tile_periodic(m, n, conn, conn_err, Lx, Ly, Tx, Ty,
                                          copies=args.copies)
    print(f"  {len(x_t)} unique points after tiling and de-duplication")

    print("Building 2D interpolant...")
    interp_cc = make_interp(x_t, y_t, cc_t)
    interp_ec = make_interp(x_t, y_t, ec_t)

    # Sample along each boundary path
    N = args.samples
    paths = {"v": v, "u": u, "w": w}
    colors = {"v": "tab:red", "u": "tab:blue", "w": "tab:green"}

    # Physical lengths (for x-axis label: fraction vs actual length)
    def path_length(dm, dn):
        dx, dy = mn_to_xy(dm, dn)
        return np.sqrt(dx**2 + dy**2)

    lengths = {k: path_length(*paths[k]) for k in paths}
    print(f"\nPath physical lengths:")
    for k, (dm, dn) in paths.items():
        print(f"  {k} = ({dm},{dn}):  |L| = {lengths[k]:.4f}")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), dpi=150)
    fig.suptitle(
        f"Boundary-direction correlators: 4-5-6 triangle\n"
        f"Lx={Lx}, Ly={Ly}, Tx={Tx}, Ty={Ty}  "
        f"k1=k2=k3=ln(3)/4 (isotropic critical)",
        fontsize=11
    )

    nan_warned = set()
    for name, (dm, dn) in paths.items():
        t, cc_vals, ec_vals, mask = sample_path(interp_cc, interp_ec, dm, dn, N)
        frac_valid = mask.sum() / len(mask)
        if frac_valid < 0.9 and name not in nan_warned:
            print(f"  WARNING: {name} direction only {100*frac_valid:.0f}% points "
                  f"inside convex hull — try --copies 3")
            nan_warned.add(name)
        else:
            print(f"  {name}: {100*frac_valid:.1f}% of samples inside interpolation hull")

        # Top panel: connected correlator vs fraction of path
        axes[0].plot(t[mask], cc_vals[mask], color=colors[name], lw=2.0, label=f"{name} = {(dm,dn)}")
        axes[0].fill_between(t[mask],
                             cc_vals[mask] - ec_vals[mask],
                             cc_vals[mask] + ec_vals[mask],
                             color=colors[name], alpha=0.2)

        # Bottom panel: log(connected correlator) to see shape distortions
        pos = mask & (cc_vals > 0)
        axes[1].plot(t[pos], np.log(cc_vals[pos]), color=colors[name], lw=2.0, label=f"{name}")
        # error bars on log:  d(log G) ≈ dG/G
        log_err = np.where(cc_vals[pos] > 0, ec_vals[pos] / cc_vals[pos], np.nan)
        axes[1].fill_between(t[pos],
                             np.log(cc_vals[pos]) - log_err,
                             np.log(cc_vals[pos]) + log_err,
                             color=colors[name], alpha=0.2)

    axes[0].set_ylabel("G_conn(t)")
    axes[0].set_title("Connected correlator along each boundary path")
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(alpha=0.25)
    axes[0].set_xlim(0, 1)

    axes[1].set_ylabel("log G_conn(t)")
    axes[1].set_xlabel("fraction of boundary path  (0 = origin, 1 = full period)")
    axes[1].set_title("Log correlator — shape differences indicate anisotropy")
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(alpha=0.25)
    axes[1].set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"\nSaved: {args.output}")

    # Print endpoint values for sanity
    print("\nEndpoint check (t=0 and t~1 should match by periodicity):")
    for name, (dm, dn) in paths.items():
        t, cc_vals, ec_vals, mask = sample_path(interp_cc, interp_ec, dm, dn, N)
        valid = np.where(mask)[0]
        if len(valid) > 2:
            i0 = valid[0]; i1 = valid[-1]
            print(f"  {name}: G(0) = {cc_vals[i0]:.5f},  G(1) = {cc_vals[i1]:.5f}  "
                  f"(diff = {abs(cc_vals[i0]-cc_vals[i1]):.2e})")


if __name__ == "__main__":
    main()
