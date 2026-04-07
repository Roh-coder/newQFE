#!/usr/bin/env python3
"""
Plot β_c finder diagnostic results.
Reads results/betac_diagnostic/diagnostic_results.json and produces
a multi-panel figure comparing the 5 β_c methods.
"""

import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

DATA = "results/betac_diagnostic/diagnostic_results.json"
OUTDIR = "results/betac_diagnostic"
EXACT_BC_INF = np.log(3) / 4  # 0.27465... (infinite volume)

with open(DATA) as f:
    data = json.load(f)

methods_order = ["exact", "gc3", "gc2", "gc2_boot", "gc2_histat"]
labels = {
    "exact": "Exact (ref β_c)",
    "gc3": "GC 3-pass",
    "gc2": "GC 2-pass",
    "gc2_boot": "GC2 + bootstrap",
    "gc2_histat": "GC2 hi-stat (2×)",
}
colors = {
    "exact": "#2ecc71",
    "gc3": "#e74c3c",
    "gc2": "#3498db",
    "gc2_boot": "#9b59b6",
    "gc2_histat": "#e67e22",
}

# Gather per-method arrays
md = {}
for m in methods_order:
    entries = [d for d in data if d["method"] == m]
    md[m] = {
        "beta_c": np.array([d["beta_c"] for d in entries]),
        "raw": np.array([d["boundary_slices"]["chi2_ndof"] for d in entries]),
        "z2": np.array([d["boundary_slices_normed"]["chi2_ndof"] for d in entries]),
        "z4": np.array([d["boundary_slices_normed_quartic"]["chi2_ndof"] for d in entries]),
        "scan_t": np.array([d["scan_time_s"] for d in entries]),
        "total_t": np.array([d["total_time_s"] for d in entries]),
    }

exact_bc = md["exact"]["beta_c"][0]
n_trials = len(md["exact"]["beta_c"])

# =========================================================================
# Figure: 2×3 grid
# =========================================================================
fig = plt.figure(figsize=(18, 11))
fig.suptitle(
    f"β_c Finder Diagnostic — {n_trials} trials × 5 methods  "
    f"(8×8 lattice, r₁=r₂=1)",
    fontsize=14, fontweight="bold",
)
gs = GridSpec(2, 3, hspace=0.35, wspace=0.30,
              left=0.06, right=0.97, top=0.92, bottom=0.06)

# ---------------------------------------------------------------------------
# Panel 1: β_c distributions (box + strip)
# ---------------------------------------------------------------------------
ax1 = fig.add_subplot(gs[0, 0])
positions = np.arange(len(methods_order))
bp_data = [md[m]["beta_c"] for m in methods_order]
bp = ax1.boxplot(bp_data, positions=positions, widths=0.5, patch_artist=True,
                 showfliers=False, medianprops=dict(color="k", linewidth=1.5))
for patch, m in zip(bp["boxes"], methods_order):
    patch.set_facecolor(colors[m])
    patch.set_alpha(0.5)
for i, m in enumerate(methods_order):
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(md[m]["beta_c"]))
    ax1.scatter(positions[i] + jitter, md[m]["beta_c"], s=18, alpha=0.7,
                color=colors[m], edgecolor="none", zorder=3)
ax1.axhline(exact_bc, ls="--", color="#2ecc71", alpha=0.7, lw=1,
            label=f"ref β_c = {exact_bc:.6f}")
ax1.set_xticks(positions)
ax1.set_xticklabels([labels[m] for m in methods_order], rotation=25, ha="right", fontsize=8)
ax1.set_ylabel("β_c found")
ax1.set_title("β_c distributions")
ax1.legend(fontsize=7, loc="upper left")

# ---------------------------------------------------------------------------
# Panel 2: β_c bias (mean ± std)
# ---------------------------------------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
scan_methods = [m for m in methods_order if m != "exact"]
biases = [md[m]["beta_c"].mean() - exact_bc for m in scan_methods]
stds = [md[m]["beta_c"].std() for m in scan_methods]
x = np.arange(len(scan_methods))
bars = ax2.bar(x, biases, yerr=stds, capsize=5, width=0.6,
               color=[colors[m] for m in scan_methods], alpha=0.7,
               edgecolor="k", linewidth=0.5)
ax2.axhline(0, ls="-", color="k", lw=0.5)
ax2.set_xticks(x)
ax2.set_xticklabels([labels[m] for m in scan_methods], rotation=25, ha="right", fontsize=8)
ax2.set_ylabel("β_c bias  (found − ref)")
ax2.set_title("β_c bias ± std  (lower = better)")
for i, (b, s) in enumerate(zip(biases, stds)):
    ax2.text(i, b + s + 0.0001, f"+{b:.5f}", ha="center", va="bottom", fontsize=7)

# ---------------------------------------------------------------------------
# Panel 3: Z⁴ cost distributions (log scale, violin-ish)
# ---------------------------------------------------------------------------
ax3 = fig.add_subplot(gs[0, 2])
for i, m in enumerate(methods_order):
    vals = md[m]["z4"]
    jitter = np.random.default_rng(7 + i).uniform(-0.15, 0.15, len(vals))
    ax3.scatter(i + jitter, vals, s=20, alpha=0.6, color=colors[m],
                edgecolor="none", zorder=3)
    ax3.plot([i - 0.25, i + 0.25], [np.median(vals)] * 2,
             color="k", lw=2, zorder=4)
ax3.set_yscale("log")
ax3.set_xticks(positions)
ax3.set_xticklabels([labels[m] for m in methods_order], rotation=25, ha="right", fontsize=8)
ax3.set_ylabel("Normed Z⁴ cost (χ²/ndof)")
ax3.set_title("Z⁴ cost by method  (lower = better)")

# ---------------------------------------------------------------------------
# Panel 4: |β_c bias| vs Z⁴ cost (scatter, all methods)
# ---------------------------------------------------------------------------
ax4 = fig.add_subplot(gs[1, 0])
for m in methods_order:
    bias = np.abs(md[m]["beta_c"] - exact_bc)
    z4 = md[m]["z4"]
    ax4.scatter(bias, z4, s=25, alpha=0.6, color=colors[m],
                label=labels[m], edgecolor="none")
ax4.set_xlabel("|β_c − ref β_c|")
ax4.set_ylabel("Normed Z⁴ cost")
ax4.set_yscale("log")
ax4.set_title("|β_c bias| vs Z⁴ cost")
ax4.legend(fontsize=7, loc="upper left")

# Pearson on log-log (excluding exact which has zero bias)
all_bias = np.concatenate([np.abs(md[m]["beta_c"] - exact_bc) for m in scan_methods])
all_z4 = np.concatenate([md[m]["z4"] for m in scan_methods])
mask = all_bias > 0
from scipy.stats import pearsonr
r, p = pearsonr(np.log10(all_bias[mask]), np.log10(all_z4[mask]))
ax4.text(0.95, 0.05, f"r(log|bias|, logZ⁴) = {r:.2f}\np = {p:.1e}",
         transform=ax4.transAxes, ha="right", va="bottom", fontsize=8,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

# ---------------------------------------------------------------------------
# Panel 5: Head-to-head win rates vs gc2
# ---------------------------------------------------------------------------
ax5 = fig.add_subplot(gs[1, 1])
trials = sorted(set(d["trial"] for d in data))
competitors = ["gc3", "gc2_boot", "gc2_histat"]
cost_names = ["z2", "z4"]
cost_labels = ["Normed Z²", "Normed Z⁴"]
bar_width = 0.25
xpos = np.arange(len(competitors))
for ci, (cn, cl) in enumerate(zip(cost_names, cost_labels)):
    win_rates = []
    for comp in competitors:
        wins = 0
        for t in trials:
            gc2_val = [d[f"boundary_slices{'_normed' if cn != 'raw' else ''}{'_quartic' if cn == 'z4' else ''}"]["chi2_ndof"]
                       for d in data if d["trial"] == t and d["method"] == "gc2"]
            comp_val = [d[f"boundary_slices{'_normed' if cn != 'raw' else ''}{'_quartic' if cn == 'z4' else ''}"]["chi2_ndof"]
                        for d in data if d["trial"] == t and d["method"] == comp]
            if gc2_val and comp_val and comp_val[0] < gc2_val[0]:
                wins += 1
        win_rates.append(wins / len(trials) * 100)
    ax5.bar(xpos + ci * bar_width, win_rates, width=bar_width,
            color=["#2ecc71", "#e74c3c"][ci], alpha=0.7, label=cl,
            edgecolor="k", linewidth=0.5)
    for i, wr in enumerate(win_rates):
        ax5.text(i + ci * bar_width, wr + 1, f"{wr:.0f}%", ha="center",
                 va="bottom", fontsize=8)

ax5.axhline(50, ls="--", color="gray", lw=1, alpha=0.5)
ax5.set_xticks(xpos + bar_width / 2)
ax5.set_xticklabels([labels[m] for m in competitors], fontsize=8)
ax5.set_ylabel("Win rate vs GC2 (%)")
ax5.set_title(f"Head-to-head vs GC2 ({len(trials)} paired trials)")
ax5.set_ylim(0, 105)
ax5.legend(fontsize=8)

# ---------------------------------------------------------------------------
# Panel 6: Median cost × time (efficiency)
# ---------------------------------------------------------------------------
ax6 = fig.add_subplot(gs[1, 2])
for m in methods_order:
    med_z4 = np.median(md[m]["z4"])
    mean_t = md[m]["total_t"].mean()
    ax6.scatter(mean_t, med_z4, s=120, color=colors[m], edgecolor="k",
                linewidth=1, zorder=3, label=labels[m])
    ax6.annotate(labels[m], (mean_t, med_z4), fontsize=7,
                 xytext=(8, 4), textcoords="offset points")
ax6.set_xlabel("Mean wall time per trial (s)")
ax6.set_ylabel("Median Z⁴ cost")
ax6.set_yscale("log")
ax6.set_title("Efficiency: cost vs time  (lower-left = best)")

out = os.path.join(OUTDIR, "diagnostic_plots.png")
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.close()
