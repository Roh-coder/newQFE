"""
landscape_isotropy_dev.py — Section 7.4 step 2 of plan.md.

Re-cast Tier-1 anisotropy probes B (susceptibility) and C (amplitude
offset) as **deviation-from-isotropy** observables instead of raw
triplet ratios. The motivation is the failure mode documented in the
plan (Section 7.2): at small (r1, r2) the test correlator decays so
fast that all per-cycle observables collapse toward each other, so
their ratios become trivially close to the reference's, producing a
spurious minimum at the small-r corner.

For each cycle i ∈ {u, v, w} we form a per-direction observable
X_i ∈ {χ_i, c_i, λ_i (≡ −log(t_half))}. We then form the scale-free
"isotropy-violation" scalar

    Δ_X = (1/X̄²) · Σ_i (X_i − X̄)²

separately on the test and the reference. The cost contribution is

    J_X = (Δ_X^test − Δ_X^ref)² / Δ_X^ref²

so a perfect anisotropy match scores zero.  Both terms vanish only at
genuine isotropy (which the 4-5-6 reference is NOT — its Δ is large)
and the small-r corner now scores very far from the reference's Δ,
not very close.

Plotted panels:
    A : λ-deviation cost            (decay-rate isotropy violation)
    B : χ-deviation cost            (susceptibility isotropy violation)
    C : c-deviation cost            (amplitude isotropy violation)
    Δ_χ^test, Δ_c^test              (raw test-side scalars)
    A + B + C                       (naive equal-weight sum)
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


# ---------------------------------------------------------------------
# ref pack: one entry per cycle with the t-grid, ms,ns indices, and
# the pre-interpolated reference correlator at the matching physical
# points on the test torus.
# ---------------------------------------------------------------------
def _ref_pack(ref_data, rL, rL2, rT, rT2, tL, tL2, tT, tT2, copies=2):
    iref = _tile_interp(ref_data, rL, rL2, rT, rT2, "conn", copies)
    rp = boundary_paths(rL, rL2, rT, rT2)
    tp = boundary_paths(tL, tL2, tT, tT2)
    out = []
    for (rdm, rdn), (tdm, tdn) in zip(rp, tp):
        ks, ms, ns = _direction_lattice_steps(tdm, tdn)
        N = len(ks)
        if N < 4:
            out.append(None); continue
        t = np.asarray([k / N for k in ks], dtype=float)
        rex, rey = rdm + 0.5 * rdn, _SQRT3_2 * rdn
        pts = np.column_stack([t * rex, t * rey])
        Gr = np.asarray(iref(pts), dtype=float)
        out.append(dict(t=t, Gr=Gr, ms=ms, ns=ns, N=N))
    return out


def _half_decay(t, G, drop=1):
    if len(G) < 4:
        return np.nan
    G = np.asarray(G, float); t = np.asarray(t, float)
    g0 = G[drop]; g_min = G.min()
    target = 0.5 * (g0 + g_min)
    for i in range(drop, len(G) - 1):
        if (G[i] - target) * (G[i + 1] - target) <= 0:
            f = (G[i] - target) / max(G[i] - G[i + 1], 1e-12)
            return float(t[i] + f * (t[i + 1] - t[i]))
    return np.nan


def _per_pt_triplets(test_data, ref_pack, tL, tL2, tT, tT2,
                     copies=2, drop_first=1, drop_last=1, eps=1e-12):
    """Return (chi_t, chi_r, lam_t, lam_r, c_arr) — each length-3 arrays."""
    chi_t, chi_r = [], []
    lam_t, lam_r = [], []
    c_arr        = []
    for d in ref_pack:
        if d is None:
            for L in (chi_t, chi_r, lam_t, lam_r, c_arr): L.append(np.nan)
            continue
        Gr = d["Gr"]; ms = d["ms"]; ns = d["ns"]
        N  = d["N"]; t = d["t"]
        Gt = np.full(N, np.nan)
        for i, (mk, nk) in enumerate(zip(ms, ns)):
            entry = _lookup_test_value(test_data, int(mk), int(nk),
                                       tL, tL2, tT, tT2, copies=copies)
            if entry is not None:
                Gt[i] = entry["conn"]
        m = np.isfinite(Gr) & np.isfinite(Gt) & (Gr > eps) & (Gt > eps)

        # susceptibility (sum over cycle, no endpoint drop — cycle is closed)
        chi_t.append(float(np.sum(Gt[m])) if m.any() else np.nan)
        chi_r.append(float(np.sum(Gr[m])) if m.any() else np.nan)

        # half-decay → λ ≡ −ln(t_half)  (so larger λ = faster decay)
        th_t = _half_decay(t[m], Gt[m], drop=1) if m.sum() >= 4 else np.nan
        th_r = _half_decay(t[m], Gr[m], drop=1) if m.sum() >= 4 else np.nan
        lam_t.append(-np.log(th_t) if (np.isfinite(th_t) and th_t > 0) else np.nan)
        lam_r.append(-np.log(th_r) if (np.isfinite(th_r) and th_r > 0) else np.nan)

        # multiplicative offset c  (drop both singular endpoints)
        idx = np.where(m)[0]
        if len(idx) > drop_first + drop_last + 1:
            idx_c = idx[drop_first:-drop_last] if drop_last > 0 else idx[drop_first:]
            log_diff = np.log(Gt[idx_c]) - np.log(Gr[idx_c])
            c_arr.append(float(np.exp(np.mean(log_diff))))
        else:
            c_arr.append(np.nan)

    return (np.asarray(chi_t), np.asarray(chi_r),
            np.asarray(lam_t), np.asarray(lam_r),
            np.asarray(c_arr))


def _delta(triple, eps=1e-12):
    """Scale-free isotropy-violation:  Σ_i (X_i − X̄)² / X̄²."""
    if not np.all(np.isfinite(triple)):
        return np.nan
    bar = float(np.mean(triple))
    if abs(bar) < eps:
        return np.nan
    return float(np.sum((triple - bar) ** 2) / bar ** 2)


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
        os.path.join(_HERE, "results", args.ref_tag,
                     "two_point_all_to_all.dat"))
    ref_pack = _ref_pack(ref_data, *rg, Lx, Ly, Tx, Ty)

    print(f"[load+eval] {n*n} pts  test=({Lx},{Ly},{Tx},{Ty})  ref={args.ref_tag}")
    t0 = time.time()
    COST_A = np.full((n, n), np.nan)
    COST_B = np.full((n, n), np.nan)
    COST_C = np.full((n, n), np.nan)
    DLT_X  = np.full((n, n), np.nan)   # Δ_χ^test  (diagnostic)
    DLT_C  = np.full((n, n), np.nan)   # Δ_c^test  (diagnostic)

    # Reference Δ (single set of three numbers per observable, computed
    # by walking the ref correlator on its OWN cycles is overkill —
    # ref_pack already exposes Gr at the test-tied physical sample
    # points, and the (chi_r, lam_r) triplets it returns are the same
    # for every test (r1, r2). We compute them once.
    chi_r_ref = lam_r_ref = None
    DR_X = DR_L = None

    for i, r1 in enumerate(rs):
        for j, r2 in enumerate(rs):
            pkl = os.path.join(root, "grid",
                               f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            if not os.path.exists(pkl): continue
            with open(pkl, "rb") as f:
                pt = pickle.load(f)
            chi_t, chi_r, lam_t, lam_r, c_arr = _per_pt_triplets(
                pt["test_data"], ref_pack, Lx, Ly, Tx, Ty)
            if chi_r_ref is None:
                # Cache the ref-side triplets / Δ once.
                chi_r_ref = chi_r; lam_r_ref = lam_r
                DR_X = _delta(chi_r); DR_L = _delta(lam_r)
                # Reference c is by construction (1,1,1) on itself
                # → Δ_c^ref = 0. Replace with floor for the relative cost.
                print(f"[ref] Δ_χ^ref={DR_X:.4g}  Δ_λ^ref={DR_L:.4g}  "
                      f"(Δ_c^ref ≡ 0 by construction)")
            DT_X = _delta(chi_t); DT_L = _delta(lam_t); DT_C = _delta(c_arr)
            DLT_X[i, j] = DT_X; DLT_C[i, j] = DT_C

            # Costs: relative squared deviation between Δ_test and Δ_ref.
            # For c, the ref Δ is 0 by construction; we replace the
            # denominator with a Δ_c^* = max(Δ_c^test_grid_median, 1e-3)
            # AFTER the loop, so for now store the raw difference.
            if np.isfinite(DT_L) and np.isfinite(DR_L) and DR_L > 1e-12:
                COST_A[i, j] = (DT_L - DR_L) ** 2 / DR_L ** 2
            if np.isfinite(DT_X) and np.isfinite(DR_X) and DR_X > 1e-12:
                COST_B[i, j] = (DT_X - DR_X) ** 2 / DR_X ** 2
            if np.isfinite(DT_C):
                COST_C[i, j] = DT_C  # store raw, normalize after

    # Δ_c^ref ≡ 0 → cost_C = (Δ_c^test)^2 / Δ_c^scale^2 where the
    # scale is the median Δ_c^test over the grid (target: 0).
    scale_C = float(np.nanmedian(COST_C))
    if not np.isfinite(scale_C) or scale_C < 1e-12: scale_C = 1.0
    COST_C = (COST_C ** 2) / scale_C ** 2

    finite = (np.isfinite(COST_A) & np.isfinite(COST_B) & np.isfinite(COST_C))
    COMBO = np.where(finite, COST_A + COST_B + COST_C, np.nan)
    print(f"[load+eval] done in {time.time()-t0:.1f}s")

    panels = [
        ("A: Δλ-deviation cost",         COST_A),
        ("B: Δχ-deviation cost",         COST_B),
        ("C: Δc-deviation cost",         COST_C),
        ("Δχ^test (raw)",                DLT_X),
        ("Δc^test (raw)",                DLT_C),
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
        d_t = float(np.hypot(am[0] - _TRUTH_456[0],
                             am[1] - _TRUTH_456[1]))
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
        print(f"  {name:32s} argmin=({am[0]:5.2f},{am[1]:5.2f})  "
              f"d_truth={d_t:.2f}  log10[min]={np.nanmin(Zp):.2f}")

    fig.suptitle(f"Δ-from-isotropy costs   test=({Lx},{Ly},{Tx},{Ty})  "
                 f"ref={args.ref_tag}  truth=(5.07,7.74)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = args.out or os.path.join(root, "isotropy_dev_costs.png")
    fig.savefig(out, dpi=130)
    print("→", out)


if __name__ == "__main__":
    main()
