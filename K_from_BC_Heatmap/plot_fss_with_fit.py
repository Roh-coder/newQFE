#!/usr/bin/env python3
"""Generate FSS channel plots with the same weighted quadratic fit used in analysis.

Uses k=1..7 only (k=0 excluded).
"""

import glob
import os
import pickle
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_boundary_continuum_opt import _fit_quad_in_invL, _sample_boundary_from_payload


def main():
    base = os.path.join("K_from_BC_Heatmap", "results", "boundary_opt_prod_phase1")
    grid = os.path.join(base, "grid")
    out_dir = os.path.join(base, "plots", "fss_with_fit_all_points")
    os.makedirs(out_dir, exist_ok=True)

    k_values = np.arange(1, 8, dtype=int)
    t_fracs = k_values.astype(float) / 8.0
    pat = re.compile(r"r1_([0-9]+\.[0-9]+)_r2_([0-9]+\.[0-9]+)\.pkl$")

    sizes = sorted(int(d[1:]) for d in os.listdir(grid) if d.startswith("L") and d[1:].isdigit())
    inv_l = np.array([1.0 / L for L in sizes], dtype=float)
    x_fit = np.linspace(inv_l.min(), inv_l.max(), 200)
    x_cont = 0.0

    first_l = sizes[0]
    points = []
    for p in glob.glob(os.path.join(grid, f"L{first_l}", "test", "*.pkl")):
        m = pat.search(os.path.basename(p))
        if m:
            points.append((float(m.group(1)), float(m.group(2))))
    points = sorted(set(points))

    ref_files = sorted(glob.glob(os.path.join(grid, f"L{max(sizes)}", "ref", "*.pkl")))
    if not ref_files:
        raise RuntimeError("No reference pkl files found")
    ref_name = "r1_1.000_r2_1.000.pkl"
    if not any(os.path.basename(p) == ref_name for p in ref_files):
        ref_name = os.path.basename(ref_files[0])

    rm = pat.search(ref_name)
    ref_r1 = float(rm.group(1)) if rm else float("nan")
    ref_r2 = float(rm.group(2)) if rm else float("nan")

    cmap = plt.get_cmap("tab10", 7)

    def _collect_point_channels(r1, r2, mode):
        """Collect channel arrays over sizes once to avoid repeated file I/O."""
        G_by_size = {}
        E_by_size = {}
        for L in sizes:
            if mode == "test":
                pkl = os.path.join(grid, f"L{L}", "test", f"r1_{r1:.3f}_r2_{r2:.3f}.pkl")
            else:
                pkl = os.path.join(grid, f"L{L}", "ref", ref_name)
            if not os.path.exists(pkl):
                continue

            with open(pkl, "rb") as f:
                payload = pickle.load(f)
            G, sG = _sample_boundary_from_payload(payload, t_fracs)
            G_by_size[L] = G
            E_by_size[L] = np.abs(sG)
        return G_by_size, E_by_size

    def plot_one(dataset_label, r1, r2, mode="test"):
        fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.8), sharex=True)

        G_by_size, E_by_size = _collect_point_channels(r1, r2, mode)

        for cyc in range(3):
            ax = axs[cyc]
            for ik, k in enumerate(k_values):
                x = []
                y = []
                s = []

                for L in sizes:
                    if L not in G_by_size:
                        continue
                    x.append(1.0 / float(L))
                    y.append(float(G_by_size[L][cyc, ik]))
                    s.append(float(E_by_size[L][cyc, ik]))

                if not x:
                    continue

                x = np.array(x, dtype=float)
                y = np.array(y, dtype=float)
                s = np.array(s, dtype=float)
                order = np.argsort(x)
                x = x[order]
                y = y[order]
                s = s[order]

                color = cmap(ik)
                ax.errorbar(
                    x,
                    y,
                    yerr=s,
                    fmt="o",
                    markersize=3.3,
                    elinewidth=0.9,
                    capsize=2,
                    color=color,
                    alpha=0.95,
                    label=f"k={k}" if cyc == 0 else None,
                )

                a, sa, b, c2, n = _fit_quad_in_invL(1.0 / x, y, s)
                if np.isfinite(a):
                    y_fit = a + b * x_fit + c2 * (x_fit * x_fit)
                    ax.plot(x_fit, y_fit, "-", lw=1.1, color=color, alpha=0.9)
                    if np.isfinite(sa):
                        # Continuum prediction from fit intercept with propagated error.
                        ax.errorbar(
                            [x_cont],
                            [a],
                            yerr=[sa],
                            fmt="D",
                            markersize=3.8,
                            elinewidth=1.0,
                            capsize=2,
                            color=color,
                            alpha=0.95,
                        )

            # Leave room at x=0 for continuum-point markers.
            ax.set_xlim(-0.004, float(inv_l.max()) * 1.06)

            ax.set_title(f"cycle c={cyc}")
            ax.set_xlabel("1/L")
            ax.grid(alpha=0.25)

        axs[0].set_ylabel("G_conn(t_k)  [jackknife sigma]")

        handles, labels = axs[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, ncol=7, loc="upper center", bbox_to_anchor=(0.5, 1.04), fontsize=8)

        fit_txt = (
            r"fit: $G(L)=a+b/L+c/L^2$ (weighted by $1/\sigma^2$); "
            r"circles: jackknife $\sigma_G$; diamonds at $1/L=0$: $a\pm\sigma_a$"
        )
        fig.suptitle(f"{dataset_label}: r1={r1:.1f}, r2={r2:.1f}   |   {fit_txt}", y=1.10)
        fig.tight_layout()

        if mode == "test":
            out = os.path.join(out_dir, f"fss_fit_test_r1_{r1:.3f}_r2_{r2:.3f}.png")
        else:
            out = os.path.join(out_dir, f"fss_fit_reference_r1_{r1:.3f}_r2_{r2:.3f}.png")
        fig.savefig(out, dpi=170, bbox_inches="tight")
        plt.close(fig)

    for i, (r1, r2) in enumerate(points, 1):
        print(f"[plot] {i}/{len(points)}  test r1={r1:.3f} r2={r2:.3f}")
        plot_one("FSS boundary channels (test)", r1, r2, mode="test")
    print("[plot] reference")
    plot_one("FSS boundary channels (reference)", ref_r1, ref_r2, mode="ref")

    fig, ax = plt.subplots(figsize=(8, 7.5))
    rvals = sorted(set(r for r, _ in points))
    ax.set_xlim(-0.5, len(rvals) - 0.5)
    ax.set_ylim(-0.5, len(rvals) - 0.5)
    ax.set_xticks(range(len(rvals)))
    ax.set_yticks(range(len(rvals)))
    ax.set_xticklabels([f"{v:.1f}" for v in rvals])
    ax.set_yticklabels([f"{v:.1f}" for v in rvals])
    ax.set_xlabel("r1")
    ax.set_ylabel("r2")
    ax.set_title("FSS with fit index")
    for r1, r2 in points:
        ix = rvals.index(r1)
        iy = rvals.index(r2)
        ax.text(ix, iy, f"{r1:.1f},{r2:.1f}", ha="center", va="center", fontsize=7)
    ax.grid(alpha=0.35)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fss_fit_test_index.png"), dpi=170)
    plt.close(fig)

    print(f"Generated {len(points)} test FSS+fit plots and 1 reference plot in {out_dir}")
    print("Reference key used:", ref_name)


if __name__ == "__main__":
    main()
