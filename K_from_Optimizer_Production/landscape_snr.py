"""
landscape_snr.py — heatmaps of signal-to-noise across the (r1,r2) grid.

For each precompute pkl on the test landscape, walk the integer-lattice
sites along the u/v/w cycles and compute |conn|/conn_err. Then plot:

  median SNR (all sites)
  min SNR (worst site)
  median SNR with t≈0 dropped (the t=0 spike has SNR~1000+)
  median SNR for t∈[0.25, 0.75]   (the actually-informative midrange)
  median |conn|                   (signal magnitude)
  median conn_err                 (noise magnitude)
"""
from __future__ import annotations
import argparse, json, os, pickle, sys, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from cost import boundary_paths, _direction_lattice_steps, _lookup_test_value

_TRUTH_456 = (5.0652, 7.7429)


def per_pt(td, Lx, Ly, Tx, Ty):
    conns = []; errs = []; ts = []
    for (dm, dn) in boundary_paths(Lx, Ly, Tx, Ty):
        ks, ms, ns = _direction_lattice_steps(dm, dn)
        N = len(ks)
        if N < 3: continue
        for k, m, n in zip(ks, ms, ns):
            e = _lookup_test_value(td, int(m), int(n), Lx, Ly, Tx, Ty, copies=2)
            if e is None: continue
            conns.append(abs(e["conn"]))
            errs.append(e["conn_err"])
            ts.append(k / N)
    return (np.asarray(conns), np.asarray(errs), np.asarray(ts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    print(f"[load+eval] {n*n} points on {args.landscape}")
    t0 = time.time()
    SNR_med    = np.full((n, n), np.nan)
    SNR_min    = np.full((n, n), np.nan)
    SNR_no_t0  = np.full((n, n), np.nan)
    SNR_mid    = np.full((n, n), np.nan)
    SIG_med    = np.full((n, n), np.nan)
    ERR_med    = np.full((n, n), np.nan)

    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            c, e, t = per_pt(pt["test_data"], Lx, Ly, Tx, Ty)
            m = (e > 0) & np.isfinite(c) & np.isfinite(e)
            if not m.any(): continue
            snr_all = c[m] / e[m]
            SNR_med[i, j] = np.median(snr_all)
            SNR_min[i, j] = np.min(snr_all)
            mask_no_t0 = m & (t > 0.02)  # drop t≈0
            if mask_no_t0.any():
                snr0 = c[mask_no_t0] / e[mask_no_t0]
                SNR_no_t0[i, j] = np.median(snr0)
            mask_mid = m & (t >= 0.25) & (t <= 0.75)
            if mask_mid.any():
                snr_mid = c[mask_mid] / e[mask_mid]
                SNR_mid[i, j] = np.median(snr_mid)
            SIG_med[i, j] = np.median(c[m])
            ERR_med[i, j] = np.median(e[m])
    print(f"[load+eval] done in {time.time()-t0:.1f}s")

    panels = [
        ("median SNR (all sites)",  SNR_med, "log"),
        ("min SNR (worst site)",    SNR_min, "log"),
        ("median SNR (t > 0)",      SNR_no_t0, "log"),
        ("median SNR (0.25 ≤ t ≤ 0.75)", SNR_mid, "log"),
        ("median |conn|",           SIG_med, "linear"),
        ("median conn_err",         ERR_med, "log"),
    ]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    print("\n[summary]")
    for ax, (name, Z, scale) in zip(axes.flatten(), panels):
        if np.all(np.isnan(Z)):
            ax.axis("off"); continue
        finZ = Z[np.isfinite(Z)]
        vmin = np.nanpercentile(finZ, 2)
        vmax = np.nanpercentile(finZ, 98)
        if scale == "log":
            vmin = max(vmin, 1e-12)
            norm = LogNorm(vmin=vmin, vmax=vmax)
            im = ax.pcolormesh(R1, R2, Z, shading="auto", cmap="viridis",
                               norm=norm)
        else:
            im = ax.pcolormesh(R1, R2, Z, shading="auto", cmap="viridis",
                               vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=11, mec="k", mew=1)
        ax.set_title(name)
        ax.set_xlabel("r1"); ax.set_ylabel("r2"); ax.set_aspect("equal")
        # report value at truth (nearest grid)
        i_t = int(np.argmin(np.abs(rs - _TRUTH_456[0])))
        j_t = int(np.argmin(np.abs(rs - _TRUTH_456[1])))
        v_truth = Z[i_t, j_t]
        v_med   = np.nanmedian(Z)
        print(f"  {name:35s}  truth(grid)≈{v_truth:9.3g}   "
              f"grid-median={v_med:9.3g}")

    fig.suptitle(f"Signal-to-noise heatmaps   geom=({Lx},{Ly},{Tx},{Ty})  "
                 f"n_traj_prod={manifest['n_traj_prod']}   "
                 f"truth=(5.07,7.74)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(root, "snr_heatmaps.png")
    fig.savefig(out, dpi=130)
    print("→", out)


if __name__ == "__main__":
    main()
