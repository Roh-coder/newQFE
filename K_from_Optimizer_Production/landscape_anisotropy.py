"""
landscape_anisotropy.py — heatmaps of the Tier-1 anisotropy observables
from the development plan.

For each (r1, r2) pkl on the precompute grid, walk integer lattice sites
along the three boundary cycles u, v, w. For both the test correlator
and the reference (interpolated at the same physical points), compute
three direction triplets:

  A. half-decay length t_half_i      (proxy for the cycle decay rate)
  B. integrated susceptibility χ_i = Σ G(t·p_i)
  C. multiplicative amplitude c_i = exp〈log G_test - log G_ref〉

Each triplet gives two scale-free ratios (·_u/·_w, ·_v/·_w) — match
between test and ref defines the cost. Also computes:

  F. heteroscedastic χ² in log-G  (Tier-2 wrapper, drop-in)

Plots a 3×3 grid of heatmaps.
"""
from __future__ import annotations

import argparse, json, os, pickle, sys, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mc_engine
from cost import (boundary_paths, _direction_lattice_steps,
                  _lookup_test_value, _tile_interp, _SQRT3_2)

_REF_GEOMS = {
    "_reference_Lx13_Ly16_Tx-3_Ty3": (13, 16, -3, 3),
    "_reference_Lx39_Ly48_Tx-9_Ty9": (39, 48, -9, 9),
    "_reference_Lx16_Ly16_Tx0_Ty0":  (16, 16, 0, 0),
}
_TRUTH_456 = (5.0652, 7.7429)


# --------------------------- ref pack ---------------------------------

def _ref_pack(ref_data, rL, rL2, rT, rT2, tL, tL2, tT, tT2, copies=2):
    """Pre-compute, for each test cycle, the t-grid and ref values
    at those physical-space points."""
    iref = _tile_interp(ref_data, rL, rL2, rT, rT2, "conn", copies)
    rp = boundary_paths(rL, rL2, rT, rT2)
    tp = boundary_paths(tL, tL2, tT, tT2)
    out = []
    for (rdm, rdn), (tdm, tdn) in zip(rp, tp):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 4: out.append(None); continue
        t = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5*rdn, _SQRT3_2*rdn
        pts = np.column_stack([t*rex, t*rey])
        Gr = np.asarray(iref(pts), dtype=float)
        out.append(dict(t=t, Gr=Gr, ms=ms, ns=ns, N=N))
    return out


# ---------------------- per-direction observables ----------------------

def _half_decay(t, G, drop=2):
    """t value where G crosses 0.5*G[drop] (the early-t reference).
    Linear interp; folded (use min of forward/backward). Returns nan
    if the curve never crosses."""
    if len(G) < 4:
        return np.nan
    # use min along the cycle as the floor, half-life = where curve
    # drops to (G[drop] + G_min)/2
    G = np.asarray(G, float); t = np.asarray(t, float)
    g0 = G[drop]
    g_min = G.min()
    target = 0.5 * (g0 + g_min)
    # only consider t>=t[drop]
    for i in range(drop, len(G)-1):
        if (G[i] - target) * (G[i+1] - target) <= 0:
            # linear interp
            f = (G[i] - target) / max(G[i] - G[i+1], 1e-12)
            return float(t[i] + f * (t[i+1] - t[i]))
    return np.nan


def _per_pt_observables(test_data, ref_pack, tL, tL2, tT, tT2,
                        copies=2, drop_first=1, drop_last=1, eps=1e-12):
    """Return a dict of test/ref triplets and stats for one (r1,r2) pt."""
    chi_t, chi_r = [], []
    th_t,  th_r  = [], []
    c_arr        = []
    chi2_het     = 0.0
    n_chi2       = 0
    sh_resid     = []   # for an aggregate "shape residual after c"
    for d in ref_pack:
        if d is None:
            chi_t.append(np.nan); chi_r.append(np.nan)
            th_t.append(np.nan);  th_r.append(np.nan)
            c_arr.append(np.nan); continue
        Gr = d["Gr"]; ms = d["ms"]; ns = d["ns"]
        N  = d["N"]; t = d["t"]
        Gt = np.full(N, np.nan); Et = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       tL, tL2, tT, tT2, copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
                Et[i] = entry["conn_err"]
        m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps)
        # half-decay (use full curve, drop just the t=0 sample)
        th_t.append(_half_decay(t[m], Gt[m], drop=1) if m.sum()>=4 else np.nan)
        th_r.append(_half_decay(t[m], Gr[m], drop=1) if m.sum()>=4 else np.nan)

        # susceptibility (sum over the cycle, full)
        chi_t.append(float(np.sum(Gt[m])) if m.sum() else np.nan)
        chi_r.append(float(np.sum(Gr[m])) if m.sum() else np.nan)

        # amplitude offset c — drop both endpoints
        idx = np.where(m)[0]
        if len(idx) > drop_first + drop_last + 1:
            idx_c = idx[drop_first:-drop_last] if drop_last>0 else idx[drop_first:]
            log_diff = np.log(Gt[idx_c]) - np.log(Gr[idx_c])
            c = float(np.exp(np.mean(log_diff)))
            c_arr.append(c)
            sh_resid.append(float(np.sqrt(np.mean(
                (log_diff - np.log(c))**2))))
        else:
            c_arr.append(np.nan); sh_resid.append(np.nan)

        # heteroscedastic chi2 in log-G
        if m.sum() > 2:
            sig_t2 = (Et[m]/Gt[m])**2
            # we don't have ref errors per point post-interp -> use constant
            sig_r2 = 0.005**2
            res = (np.log(Gt[m]) - np.log(Gr[m]))**2
            chi2_het += float(np.sum(res / (sig_t2 + sig_r2 + 1e-12)))
            n_chi2 += int(m.sum())

    chi_t = np.asarray(chi_t); chi_r = np.asarray(chi_r)
    th_t  = np.asarray(th_t);  th_r  = np.asarray(th_r)
    c_arr = np.asarray(c_arr)

    out = dict(chi_t=chi_t, chi_r=chi_r, th_t=th_t, th_r=th_r,
               c=c_arr, chi2_het=chi2_het / max(n_chi2, 1),
               shape_resid=float(np.nanmean(sh_resid)) if len(sh_resid) else np.nan)

    # Triplet ratios (u/w, v/w).  We pick whichever direction is
    # last (w) as the reference index = 2.
    def _ratios(triple):
        if np.any(~np.isfinite(triple)) or abs(triple[2]) < 1e-12:
            return (np.nan, np.nan)
        return (float(triple[0]/triple[2]), float(triple[1]/triple[2]))
    Rchi_t = _ratios(chi_t); Rchi_r = _ratios(chi_r)
    Rth_t  = _ratios(th_t);  Rth_r  = _ratios(th_r)
    Rc_t   = _ratios(c_arr); Rc_r   = (1.0, 1.0)   # ref-vs-ref is unity

    # cost contributions: sum of (R_test - R_ref)^2 / R_ref^2
    def _cost(Rt, Rr):
        if not np.all(np.isfinite(Rt + Rr)) or any(abs(r) < 1e-9 for r in Rr):
            return np.nan
        return float(((Rt[0]-Rr[0])/Rr[0])**2 + ((Rt[1]-Rr[1])/Rr[1])**2)

    out.update(dict(
        cost_chi   = _cost(Rchi_t, Rchi_r),     # B
        cost_thalf = _cost(Rth_t,  Rth_r),      # A
        cost_c     = _cost(Rc_t,   Rc_r),       # C
        Rchi_t=Rchi_t, Rchi_r=Rchi_r,
        Rth_t=Rth_t,   Rth_r=Rth_r,
        Rc_t=Rc_t,
    ))
    return out


# ------------------------------ main ----------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--landscape", default="Lx16_Ly16_Tx0_Ty0")
    ap.add_argument("--ref-tag",   default="_reference_Lx39_Ly48_Tx-9_Ty9")
    ap.add_argument("--out",       default=None)
    args = ap.parse_args()

    root = os.path.join(_HERE, "results", "_landscape", args.landscape)
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    Lx, Ly, Tx, Ty = manifest["geom"]
    rs = np.arange(manifest["r_min"], manifest["r_max"] + 1e-9,
                   manifest["r_step"])
    n = len(rs)

    rg = _REF_GEOMS[args.ref_tag]
    ref_data = mc_engine.load_all_to_all(
        os.path.join(_HERE, "results", args.ref_tag, "two_point_all_to_all.dat"))
    ref_pack = _ref_pack(ref_data, *rg, Lx, Ly, Tx, Ty)

    print(f"[load+eval] {n*n} pts  test=({Lx},{Ly},{Tx},{Ty})  ref={args.ref_tag}")
    t0 = time.time()
    COST_A = np.full((n, n), np.nan)
    COST_B = np.full((n, n), np.nan)
    COST_C = np.full((n, n), np.nan)
    COST_F = np.full((n, n), np.nan)
    SHAPE  = np.full((n, n), np.nan)
    COMBO  = np.full((n, n), np.nan)

    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            o = _per_pt_observables(pt["test_data"], ref_pack,
                                    Lx, Ly, Tx, Ty)
            COST_A[i, j] = o["cost_thalf"]
            COST_B[i, j] = o["cost_chi"]
            COST_C[i, j] = o["cost_c"]
            COST_F[i, j] = o["chi2_het"]
            SHAPE[i, j]  = o["shape_resid"]
    # combined: simple sum of the three Tier-1 costs (each scale-free)
    finite = (np.isfinite(COST_A) & np.isfinite(COST_B) & np.isfinite(COST_C))
    COMBO = np.where(finite, COST_A + COST_B + COST_C, np.nan)
    print(f"[load+eval] done in {time.time()-t0:.1f}s")

    panels = [
        ("A: half-decay ratio cost",     COST_A),
        ("B: susceptibility ratio cost", COST_B),
        ("C: amplitude ratio cost",      COST_C),
        ("F: heteroscedastic χ²(log G)", COST_F),
        ("shape residual (post-c)",      SHAPE),
        ("Combined  A+B+C",              COMBO),
    ]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R1, R2 = np.meshgrid(rs, rs, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
    print("\n[summary]")
    for ax, (name, Z) in zip(axes.flatten(), panels):
        if np.all(np.isnan(Z)):
            ax.axis("off"); continue
        ij = np.unravel_index(np.nanargmin(Z), Z.shape)
        am = (rs[ij[0]], rs[ij[1]])
        d_t = float(np.hypot(am[0]-_TRUTH_456[0], am[1]-_TRUTH_456[1]))
        # use log color for the cost panels (positive, wide range)
        Zp = np.log10(np.maximum(Z, 1e-9))
        vmin, vmax = np.nanpercentile(Zp, [2, 98])
        im = ax.pcolormesh(R1, R2, Zp, shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.plot(_TRUTH_456[0], _TRUTH_456[1], "rX", ms=11, mec="k", mew=1)
        ax.plot(am[0], am[1], "w*", ms=12, mec="k", mew=1)
        ax.set_title(f"{name}\nargmin=({am[0]:.1f},{am[1]:.1f})  "
                     f"d_truth={d_t:.2f}", fontsize=10)
        ax.set_xlabel("r1"); ax.set_ylabel("r2"); ax.set_aspect("equal")
        print(f"  {name:38s} argmin=({am[0]:5.2f},{am[1]:5.2f})  "
              f"d_truth={d_t:.2f}  log10[min]={np.nanmin(Zp):.2f}")

    fig.suptitle(f"Tier-1 anisotropy costs   test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = args.out or os.path.join(root, "anisotropy_costs.png")
    fig.savefig(out, dpi=130)
    print("→", out)


if __name__ == "__main__":
    main()
