#!/usr/bin/env python3
"""
Plot optimizer score comparison results.
Reads results/stability_normed_compare/compare_results.json and produces
a multi-panel figure comparing the 3 cost functions evaluated on the same MC data.
"""

import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, spearmanr

DATA = "results/stability_normed_compare/compare_results.json"
OUTDIR = "results/stability_normed_compare"
EXACT_BC = 0.23888464199556  # 8×8 ref β_c

with open(DATA) as f:
    data = json.load(f)

n = len(data)
betas = np.array([d["beta_c"] for d in data])
raw = np.array([d["boundary_slices"]["chi2_ndof"] for d in data])
z2 = np.array([d["boundary_slices_normed"]["chi2_ndof"] for d in data])
z4 = np.array([d["boundary_slices_normed_quartic"]["chi2_ndof"] for d in data])
bias = np.abs(betas - EXACT_BC)

cost_data = {"Raw quartic (diff⁴)": raw, "Normed Z²": z2, "Normed Z⁴": z4}
cost_colors = {"Raw quartic (diff⁴)": "#3498db", "Normed Z²": "#2ecc71", "Normed Z⁴": "#e74c3c"}

# =========================================================================
fig = plt.figure(figsize=(18, 12))
fig.suptitle(
    f"Optimizer Score Comparison — {n} paired trials  "
    f"(8×8, r₁=r₂=1, gc2 finder, same MC data per trial)",
    fontsize=14, fontweight="bold",
)
gs = GridSpec(2, 3, hspace=0.35, wspace=0.35,
              left=0.06, right=0.97, top=0.92, bottom=0.06)

# ---------------------------------------------------------------------------
# Panel 1: Cost distributions (histograms in log space)
# ---------------------------------------------------------------------------
ax1 = fig.add_subplot(gs[0, 0])
for name, vals in cost_data.items():
    lv = np.log10(vals)
    bins = np.linspace(lv.min() - 0.2, lv.max() + 0.2, 25)
    ax1.hist(lv, bins=bins, alpha=0.5, color=cost_colors[name], label=name,
             edgecolor="k", linewidth=0.3)
ax1.set_xlabel("log₁₀(χ²/ndof)")
ax1.set_ylabel("Count")
ax1.set_title("Cost distributions (log scale)")
ax1.legend(fontsize=8)

# ---------------------------------------------------------------------------
# Panel 2: Summary stats (CoV bar chart)
# ---------------------------------------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
names = list(cost_data.keys())
covs = [cost_data[n].std() / cost_data[n].mean() for n in names]
bars = ax2.bar(range(len(names)), covs, color=[cost_colors[n] for n in names],
               alpha=0.7, edgecolor="k", linewidth=0.5, width=0.6)
for i, (c, n_) in enumerate(zip(covs, names)):
    ax2.text(i, c + 0.03, f"{c:.2f}", ha="center", va="bottom", fontsize=10,
             fontweight="bold")
ax2.set_xticks(range(len(names)))
ax2.set_xticklabels(names, fontsize=9)
ax2.set_ylabel("Coefficient of Variation (σ/μ)")
ax2.set_title("Stability: CoV  (lower = more reliable)")
ax2.set_ylim(0, max(covs) * 1.3)

# Add median and IQR as text
for i, n_ in enumerate(names):
    v = cost_data[n_]
    q25, med, q75 = np.percentile(v, [25, 50, 75])
    ax2.text(i, -0.02, f"med={med:.2e}\nIQR={q75-q25:.2e}", ha="center",
             va="top", fontsize=6, transform=ax2.get_xaxis_transform())

# ---------------------------------------------------------------------------
# Panel 3: log(Z²) vs log(Z⁴) scatter (paired correlation)
# ---------------------------------------------------------------------------
ax3 = fig.add_subplot(gs[0, 2])
ax3.scatter(np.log10(z2), np.log10(z4), s=20, alpha=0.6, c=betas,
            cmap="coolwarm", edgecolor="none")
cb = plt.colorbar(ax3.collections[0], ax=ax3, pad=0.02)
cb.set_label("β_c", fontsize=8)
# Fit line
lz2, lz4 = np.log10(z2), np.log10(z4)
m_fit, b_fit = np.polyfit(lz2, lz4, 1)
xx = np.linspace(lz2.min(), lz2.max(), 50)
ax3.plot(xx, m_fit * xx + b_fit, "k--", lw=1, alpha=0.7)
rp, pp = pearsonr(lz2, lz4)
rs, ps = spearmanr(lz2, lz4)
ax3.set_xlabel("log₁₀(Normed Z²)")
ax3.set_ylabel("log₁₀(Normed Z⁴)")
ax3.set_title("Z² vs Z⁴ (paired, same MC data)")
ax3.text(0.05, 0.95, f"Pearson r = {rp:.3f}\nSpearman ρ = {rs:.3f}\nslope = {m_fit:.2f}",
         transform=ax3.transAxes, va="top", fontsize=9,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

# ---------------------------------------------------------------------------
# Panel 4: |β_c bias| vs each cost
# ---------------------------------------------------------------------------
ax4 = fig.add_subplot(gs[1, 0])
for name, vals in cost_data.items():
    ax4.scatter(bias, vals, s=15, alpha=0.5, color=cost_colors[name],
                label=name, edgecolor="none")
ax4.set_xlabel("|β_c − ref β_c|")
ax4.set_ylabel("Cost (χ²/ndof)")
ax4.set_yscale("log")
ax4.set_title("|β_c bias| drives ALL cost functions")
ax4.legend(fontsize=7)
# Annotate correlations
for name, vals in cost_data.items():
    r, p = pearsonr(bias, np.log10(vals))
    ax4.text(0.98, {"Raw quartic (diff⁴)": 0.35, "Normed Z²": 0.22, "Normed Z⁴": 0.09}[name],
             f"{name}: r = {r:.2f}", transform=ax4.transAxes, ha="right",
             fontsize=8, color=cost_colors[name], fontweight="bold")

# ---------------------------------------------------------------------------
# Panel 5: β_c distribution
# ---------------------------------------------------------------------------
ax5 = fig.add_subplot(gs[1, 1])
ax5.hist(betas, bins=20, color="#3498db", alpha=0.6, edgecolor="k", linewidth=0.5)
ax5.axvline(EXACT_BC, ls="--", color="#2ecc71", lw=2, label=f"ref β_c = {EXACT_BC:.6f}")
ax5.axvline(betas.mean(), ls="-", color="#e74c3c", lw=1.5,
            label=f"mean = {betas.mean():.6f}")
ax5.set_xlabel("β_c")
ax5.set_ylabel("Count")
ax5.set_title(f"β_c distribution ({n} trials, gc2 finder)")
ax5.legend(fontsize=8)
ax5.text(0.95, 0.85, f"std = {betas.std():.6f}\nbias = +{betas.mean()-EXACT_BC:.6f}",
         transform=ax5.transAxes, ha="right", fontsize=9,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

# ---------------------------------------------------------------------------
# Panel 6: Z⁴/Z² ratio stability
# ---------------------------------------------------------------------------
ax6 = fig.add_subplot(gs[1, 2])
ratio = z4 / z2
ax6.hist(ratio, bins=20, color="#9b59b6", alpha=0.6, edgecolor="k", linewidth=0.5)
ax6.axvline(ratio.mean(), ls="-", color="k", lw=1.5, label=f"mean = {ratio.mean():.2f}")
ax6.axvline(np.median(ratio), ls="--", color="k", lw=1, label=f"median = {np.median(ratio):.2f}")
ax6.axvline(3.0, ls=":", color="gray", lw=1, label="Null expectation (E[Z⁴]/E[Z²] = 3)")
ax6.set_xlabel("Z⁴ / Z²")
ax6.set_ylabel("Count")
ax6.set_title("Z⁴/Z² ratio  (constant ⇒ redundant)")
ax6.legend(fontsize=7)
ax6.text(0.95, 0.70, f"CoV = {ratio.std()/ratio.mean():.2f}",
         transform=ax6.transAxes, ha="right", fontsize=10,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

out = os.path.join(OUTDIR, "score_comparison_plots.png")
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.close()
