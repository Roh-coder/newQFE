#!/usr/bin/env python3
"""Re-render thermalization plot using running means; write summary md."""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "results", "diag_fss_autocorr")
PLOTS = os.path.join(ROOT, "plots")
RAW = os.path.join(ROOT, "raw")

# Reload traces from scratch dirs
SCRATCH = os.path.join(ROOT, "_scratch")
sizes = sorted(int(d[1:]) for d in os.listdir(SCRATCH) if d.startswith("L"))

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax_run, ax_inst = axes
for L in sizes:
    tfiles = sorted(glob.glob(os.path.join(SCRATCH, f"L{L}", "*", "traces_*.dat")),
                    key=os.path.getmtime)
    if not tfiles:
        continue
    rows = []
    with open(tfiles[-1]) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            rows.append(float(p[1]))  # |m|
    s = np.asarray(rows)
    n = len(s)
    # running mean
    cum = np.cumsum(s) / np.arange(1, n + 1)
    ax_run.plot(np.arange(n), cum, lw=1.0, label=f"L={L}")
    # block-mean every 100 sweeps for visibility
    nb = n // 100
    bm = s[: nb * 100].reshape(nb, 100).mean(axis=1)
    ax_inst.plot(np.arange(nb) * 100, bm, lw=0.8, alpha=0.8, label=f"L={L}")

ax_run.set_ylabel(r"running mean $\langle |m|\rangle$")
ax_run.set_title("Thermalization diagnostics (post-thermalization production chain)")
ax_run.legend(ncol=4, fontsize=8)
ax_run.grid(alpha=0.3)
ax_inst.set_xlabel("measured trajectory index")
ax_inst.set_ylabel(r"block-averaged $|m|$ (block=100)")
ax_inst.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(PLOTS, "thermalization_running.png")
fig.savefig(out, dpi=120)
print(f"wrote {out}")
