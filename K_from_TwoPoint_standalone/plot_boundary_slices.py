#!/usr/bin/env python3
"""
plot_boundary_slices.py — Visualise and compare boundary-direction correlators.

For a twisted parallelogram on torus basis vectors
  v = (Lx, Ty),   u = (Tx, -Ly),   w = -(v+u)
the three boundary paths each have a distinct physical period length
  L_v, L_u, L_w  (from the triangular metric: |a,b|^2 = a^2+ab+b^2)

Reference and test must share the SAME geometry (same Lx, Ly, Tx, Ty), so
they share the same physical boundary periods L_v, L_u, L_w.  The isotropic
reference at criticality (k1=k2=k3=Kc) gives three non-collapsed slices
because L_v != L_u != L_w.  This non-collapse pattern is the target signal:
the optimizer finds (k1,k2,k3) such that the test slices reproduce it exactly.

Modes
-----
Single dataset (--dat only): plot the three boundary slices for one run.
Comparison (--ref + --test):  overlay reference and test per direction; show
  residuals and integrated cost.

Usage examples
--------------
# Single: just the isotropic reference slices
python plot_boundary_slices.py \\
    --dat  results/test_456_iso/.../two_point_all_to_all.dat \\
    --Lx 13 --Ly 16 --Tx 3 --Ty -3 \\
    --title "4-5-6 isotropic reference  k1=k2=k3" \\
    --output results/boundary_456_ref.png

# Comparison: reference vs anisotropic test
python plot_boundary_slices.py \\
    --ref  results/test_456_iso/.../two_point_all_to_all.dat \\
    --test results/test_456_aniso/.../two_point_all_to_all.dat \\
    --Lx 13 --Ly 16 --Tx 3 --Ty -3 \\
    --title "4-5-6  ref(iso) vs test(k1=1.2,k2=0.8,k3=1)" \\
    --output results/boundary_456_compare.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_a2a(path: Path):
    """Load two_point_all_to_all.dat → arrays (m,n,corr,err,conn,conn_err)."""
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
    if not rows:
        raise SystemExit(f"No data in {path}")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4], arr[:, 5]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def mn_to_xy(m, n):
    """Triangular lattice: e1=(1,0), e2=(1/2,sqrt(3)/2)."""
    m = np.asarray(m, dtype=float)
    n = np.asarray(n, dtype=float)
    return m + 0.5 * n, (np.sqrt(3.0) / 2.0) * n


def path_length(dm, dn):
    """Physical length of displacement (dm,dn) in triangular metric."""
    dx, dy = mn_to_xy(float(dm), float(dn))
    return float(np.sqrt(dx**2 + dy**2))


def tile_and_interp(m, n, conn, conn_err, Lx, Ly, Tx, Ty, copies=2):
    """Tile ±copies along v=(Lx,Ty), u=(Tx,-Ly); return interpolants."""
    m_list, n_list, cc_list, ec_list = [], [], [], []
    for a in range(-copies, copies + 1):
        for b in range(-copies, copies + 1):
            dm = a * Lx + b * Tx
            dn = a * Ty - b * Ly
            m_list.append(m + dm)
            n_list.append(n + dn)
            cc_list.append(conn)
            ec_list.append(conn_err)
    m_all = np.concatenate(m_list)
    n_all = np.concatenate(n_list)
    x_all, y_all = mn_to_xy(m_all, n_all)
    cc_all = np.concatenate(cc_list)
    ec_all = np.concatenate(ec_list)
    pts = np.column_stack([x_all, y_all])
    _, uid = np.unique(np.round(pts, 8), axis=0, return_index=True)
    pts = pts[uid]; cc_all = cc_all[uid]; ec_all = ec_all[uid]
    return (LinearNDInterpolator(pts, cc_all),
            LinearNDInterpolator(pts, ec_all))


def sample_path(icc, iec, dm_end, dn_end, n_samples=400):
    """Sample from (0,0) to (dm_end,dn_end); return (t, cc, ec, mask)."""
    t = np.linspace(0.0, 1.0, n_samples)
    x_p, y_p = mn_to_xy(t * dm_end, t * dn_end)
    pts = np.column_stack([x_p, y_p])
    cc = np.array(icc(pts), dtype=float)
    ec = np.array(iec(pts), dtype=float)
    mask = np.isfinite(cc) & np.isfinite(ec) & (cc > 0)
    return t, cc, ec, mask


def get_slices(m, n, conn, conn_err, Lx, Ly, Tx, Ty, paths, copies=2, n_samples=400):
    """Build interpolants and extract slices along each boundary path."""
    icc, iec = tile_and_interp(m, n, conn, conn_err, Lx, Ly, Tx, Ty, copies)
    slices = {}
    for name, (dm, dn) in paths.items():
        t, cc, ec, mask = sample_path(icc, iec, dm, dn, n_samples)
        frac = mask.sum() / len(mask)
        if frac < 0.9:
            print(f"  WARNING {name}: only {100*frac:.0f}% inside hull — "
                  f"try --copies {copies+1}")
        slices[name] = (t, cc, ec, mask)
    return slices


def slice_cost(t, cc_test, mask_test, cc_ref, mask_ref):
    """
    Integrated squared log-ratio cost for one direction.
    Uses the interior only (t in [0.05, 0.95]) to avoid endpoint artefacts.
    Returns scalar cost.
    """
    interior = (t >= 0.05) & (t <= 0.95)
    both = interior & mask_test & mask_ref & (cc_test > 0) & (cc_ref > 0)
    if both.sum() < 4:
        return float("nan")
    log_ratio = np.log(cc_test[both]) - np.log(cc_ref[both])
    # Subtract mean (amplitude-free comparison)
    log_ratio -= log_ratio.mean()
    dt = t[both][1:] - t[both][:-1]
    integrand = 0.5 * (log_ratio[:-1]**2 + log_ratio[1:]**2)
    return float(np.sum(integrand * dt))


# ---------------------------------------------------------------------------
# Single-dataset plot
# ---------------------------------------------------------------------------

def plot_single(slices, path_defs, title, output):
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), dpi=150, sharex=True)
    fig.suptitle(f"Boundary-direction correlator slices\n{title}",
                 fontsize=11, y=0.98)

    tv, ccv, ecv, maskv = slices["v"]
    for name, (dm, dn), color, label in path_defs:
        t, cc, ec, mask = slices[name]

        axes[0].plot(t[mask], cc[mask], color=color, lw=2.0, label=label)
        axes[0].fill_between(t[mask], cc[mask]-ec[mask], cc[mask]+ec[mask],
                             color=color, alpha=0.20)

        pos = mask
        axes[1].plot(t[pos], np.log(cc[pos]), color=color, lw=2.0, label=label)
        log_err = ec[pos] / np.where(cc[pos] > 0, cc[pos], np.nan)
        axes[1].fill_between(t[pos],
                             np.log(cc[pos]) - log_err,
                             np.log(cc[pos]) + log_err,
                             color=color, alpha=0.20)

        both = mask & maskv & (ccv > 0)
        ratio = cc[both] / ccv[both]
        rerr  = ratio * np.sqrt((ec[both]/cc[both])**2 + (ecv[both]/ccv[both])**2)
        axes[2].plot(t[both], ratio, color=color, lw=2.0, label=label)
        axes[2].fill_between(t[both], ratio-rerr, ratio+rerr, color=color, alpha=0.20)

    axes[2].axhline(1.0, color="gray", lw=0.9, ls="--")
    axes[0].set_ylabel("G_conn(t)"); axes[0].set_title("(a) Connected correlator")
    axes[1].set_ylabel("log G_conn(t)"); axes[1].set_title("(b) Log correlator [shape, amplitude-free]")
    axes[2].set_ylabel("G / G_v"); axes[2].set_title("(c) Ratio to v-direction")
    axes[2].set_xlabel("fraction of boundary path  (0 = origin → 1 = full period)")
    for ax in axes:
        ax.legend(frameon=False, fontsize=8); ax.grid(alpha=0.22)
    axes[2].set_xlim(0, 1)
    plt.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=150)
    print(f"Saved → {output}")


# ---------------------------------------------------------------------------
# Comparison plot: reference vs test
# ---------------------------------------------------------------------------

def plot_comparison(slices_ref, slices_test, path_defs, title, output):
    """
    3 rows (one per direction v, u, w), 2 columns:
      Left:  log G_ref and log G_test overlaid — the curves to match
      Right: residual log(G_test/G_ref), mean-subtracted — anisotropy signal
    Also prints the per-direction and total integrated cost.
    """
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), dpi=150)
    fig.suptitle(f"Boundary-slice comparison  (ref = dashed, test = solid)\n{title}",
                 fontsize=11, y=0.99)

    dir_colors = {"v": "tab:red", "u": "tab:blue", "w": "tab:green"}
    total_cost = 0.0

    for row, (name, (dm, dn), color, label) in enumerate(path_defs):
        ax_L = axes[row, 0]   # log-correlator overlay
        ax_R = axes[row, 1]   # residual

        t_r, cc_r, ec_r, mask_r = slices_ref[name]
        t_t, cc_t, ec_t, mask_t = slices_test[name]

        # --- Left panel: log G overlay ---
        for (t, cc, ec, mask, ls, lbl) in [
            (t_r, cc_r, ec_r, mask_r, "--", "ref (iso)"),
            (t_t, cc_t, ec_t, mask_t, "-",  "test"),
        ]:
            pos = mask
            ax_L.plot(t[pos], np.log(cc[pos]), color=color, lw=1.8, ls=ls, label=lbl)
            log_err = ec[pos] / np.where(cc[pos] > 0, cc[pos], np.nan)
            ax_L.fill_between(t[pos],
                              np.log(cc[pos]) - log_err,
                              np.log(cc[pos]) + log_err,
                              color=color, alpha=0.18)

        Lphys = path_length(dm, dn)
        ax_L.set_title(f"Direction {name} = ({dm},{dn})  |L|={Lphys:.3f}\nlog G")
        ax_L.legend(frameon=False, fontsize=8)
        ax_L.grid(alpha=0.22)
        ax_L.set_xlim(0, 1)

        # --- Right panel: residual log(G_test/G_ref), mean-subtracted ---
        interior = (t_t >= 0.05) & (t_t <= 0.95)
        both = interior & mask_t & mask_r & (cc_t > 0) & (cc_r > 0)
        if both.sum() > 4:
            log_ratio = np.log(cc_t[both]) - np.log(cc_r[both])
            log_ratio_err = np.sqrt((ec_t[both]/cc_t[both])**2 +
                                    (ec_r[both]/cc_r[both])**2)
            mean_lr = log_ratio.mean()
            resid = log_ratio - mean_lr          # amplitude-free residual
            ax_R.axhline(0.0, color="gray", lw=0.8, ls="--")
            ax_R.plot(t_t[both], resid, color=color, lw=1.8)
            ax_R.fill_between(t_t[both],
                              resid - log_ratio_err,
                              resid + log_ratio_err,
                              color=color, alpha=0.25)
            # integrated cost
            dt = t_t[both][1:] - t_t[both][:-1]
            integrand = 0.5 * (resid[:-1]**2 + resid[1:]**2)
            cost_d = float(np.sum(integrand * dt))
            total_cost += cost_d
            ax_R.set_title(f"Residual log(G_test/G_ref) − mean\ncost = {cost_d:.4e}")
        else:
            ax_R.set_title("Residual  (insufficient valid points)")

        ax_R.grid(alpha=0.22)
        ax_R.set_xlim(0, 1)

    for ax in axes[2, :]:
        ax.set_xlabel("fraction of boundary path  (0 → 1 full period)")

    plt.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=150)
    print(f"Saved → {output}")
    print(f"Total integrated boundary-slice cost = {total_cost:.4e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    # Single-dataset mode
    ap.add_argument("--dat",  default=None, help="Single two_point_all_to_all.dat")
    # Comparison mode
    ap.add_argument("--ref",  default=None, help="Reference (isotropic) dat file")
    ap.add_argument("--test", default=None, help="Test (anisotropic) dat file")
    # Geometry (shared by ref and test — they must have the same Lx,Ly,Tx,Ty)
    ap.add_argument("--Lx",  type=int, required=True)
    ap.add_argument("--Ly",  type=int, required=True)
    ap.add_argument("--Tx",  type=int, required=True)
    ap.add_argument("--Ty",  type=int, required=True)
    ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--copies",  type=int, default=2)
    ap.add_argument("--title",   default="")
    ap.add_argument("--output",  required=True)
    args = ap.parse_args()

    Lx, Ly, Tx, Ty = args.Lx, args.Ly, args.Tx, args.Ty
    Ncell = Lx * Ly + Tx * Ty

    # Boundary vectors and their physical lengths
    v = (Lx,         Ty)
    u = (Tx,        -Ly)
    w = (-(Lx + Tx), Ly - Ty)

    paths = {"v": v, "u": u, "w": w}
    path_defs = [
        ("v", v, "tab:red",   f"v=({v[0]},{v[1]})  |L|={path_length(*v):.3f}"),
        ("u", u, "tab:blue",  f"u=({u[0]},{u[1]})  |L|={path_length(*u):.3f}"),
        ("w", w, "tab:green", f"w=({w[0]},{w[1]})  |L|={path_length(*w):.3f}"),
    ]

    print(f"Geometry  Lx={Lx} Ly={Ly} Tx={Tx} Ty={Ty}  Ncell={Ncell}")
    for name, (dm, dn), _, _ in path_defs:
        print(f"  {name} = ({dm},{dn})  physical period L = {path_length(dm,dn):.4f}")

    N = args.samples
    C = args.copies

    if args.dat and not args.ref and not args.test:
        # ---- single-dataset mode ----
        print(f"Loading {args.dat}")
        m, n, _, _, conn, conn_err = load_a2a(Path(args.dat))
        print(f"  {len(m)} vectors  →  tiling & interpolating ...")
        slices = get_slices(m, n, conn, conn_err, Lx, Ly, Tx, Ty, paths, C, N)
        title = args.title or f"Lx={Lx} Ly={Ly} Tx={Tx} Ty={Ty}"
        plot_single(slices, path_defs, title, args.output)

    elif args.ref and args.test:
        # ---- comparison mode ----
        print(f"Loading reference: {args.ref}")
        mr, nr, _, _, ccr, ecr = load_a2a(Path(args.ref))
        print(f"  {len(mr)} vectors  →  tiling & interpolating ...")
        slices_ref = get_slices(mr, nr, ccr, ecr, Lx, Ly, Tx, Ty, paths, C, N)

        print(f"Loading test:      {args.test}")
        mt, nt, _, _, cct, ect = load_a2a(Path(args.test))
        print(f"  {len(mt)} vectors  →  tiling & interpolating ...")
        slices_test = get_slices(mt, nt, cct, ect, Lx, Ly, Tx, Ty, paths, C, N)

        title = args.title or f"Lx={Lx} Ly={Ly} Tx={Tx} Ty={Ty}  ref vs test"
        plot_comparison(slices_ref, slices_test, path_defs, title, args.output)

    else:
        ap.error("Provide either --dat (single) or both --ref and --test (comparison)")


if __name__ == "__main__":
    main()
