#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_typed_dat(path: Path):
    data = np.loadtxt(path)
    out = {}
    for t in sorted(set(data[:, 0].astype(int))):
        rows = data[data[:, 0].astype(int) == t]
        out[t] = {
            "r": rows[:, 1],
            "c": rows[:, 2],
            "e": rows[:, 3],
        }
    return out


def main():
    typed = Path("K_from_BC/results/runs/main/aniso_h3over4_geometry_rule_connected_jk_highstat_L48x64_typed.dat")
    out_png = Path("K_from_BC/results/plots/main/aniso_h3over4_geometry_rule_connected_jk_highstat_L48x64_xy_oblique.png")
    out_pdf = Path("K_from_BC/results/plots/main/aniso_h3over4_geometry_rule_connected_jk_highstat_L48x64_xy_oblique.pdf")

    if not typed.exists():
        raise FileNotFoundError(f"Missing input file: {typed}")

    d = load_typed_dat(typed)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)

    plotted = []
    if 4 in d:
        ax.errorbar(
            d[4]["r"],
            d[4]["c"],
            yerr=d[4]["e"],
            fmt="o-",
            lw=1.5,
            ms=3,
            capsize=2,
            label="x direction (connected jk)",
        )
        plotted.append("x")

    if 5 in d:
        ax.errorbar(
            d[5]["r"],
            d[5]["c"],
            yerr=d[5]["e"],
            fmt="o-",
            lw=1.5,
            ms=3,
            capsize=2,
            label="y direction (connected jk)",
        )
        plotted.append("y")

    if 3 in d:
        ax.errorbar(
            d[3]["r"],
            d[3]["c"],
            yerr=d[3]["e"],
            fmt="o-",
            lw=1.5,
            ms=3,
            capsize=2,
            label="oblique (x-y)",
        )
        plotted.append("oblique")
    else:
        ax.text(
            0.02,
            0.03,
            "Oblique (x-y) channel is not present in\nL48x64 typed data (only types 4,5).",
            transform=ax.transAxes,
            fontsize=10,
            color="crimson",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="crimson", alpha=0.9),
        )

    ax.set_title("Largest h=3/4 Run: Directional Correlators (L=48x64)")
    ax.set_xlabel("lattice separation bin r")
    ax.set_ylabel("connected correlation C(r)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    fig.savefig(out_pdf)

    print(f"Plotted channels: {', '.join(plotted) if plotted else 'none'}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
