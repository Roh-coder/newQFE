#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_critical_points(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()
            if len(tok) < 6:
                continue
            rows.append(
                {
                    "L": int(tok[0]),
                    "beta_peak": float(tok[1]),
                    "beta_err": float(tok[2]),
                    "chi_peak": float(tok[3]),
                    "chi_err": float(tok[4]),
                    "idx": int(tok[5]),
                }
            )
    return rows


def main():
    out_dir = Path("Potts/output")
    crit_path = out_dir / "critical_points.dat"

    if not crit_path.exists():
        raise FileNotFoundError(f"Missing {crit_path}")

    peaks = load_critical_points(crit_path)
    if not peaks:
        raise RuntimeError("No peak rows found in critical_points.dat")

    fig, ax = plt.subplots(figsize=(9, 6), dpi=140)

    for p in sorted(peaks, key=lambda x: x["L"]):
        L = p["L"]
        data_path = out_dir / f"susceptibility_L{L}.dat"
        if not data_path.exists():
            print(f"Skipping L={L}: missing {data_path}")
            continue

        data = np.loadtxt(data_path)
        beta = data[:, 0]
        chi = data[:, 3]
        chi_err = data[:, 4]

        eb = ax.errorbar(
            beta,
            chi,
            yerr=chi_err,
            marker="o",
            ms=4,
            lw=1.7,
            elinewidth=1.0,
            capsize=2,
            label=f"L={L}",
        )
        line = eb.lines[0]
        color = line.get_color()

        ax.errorbar(
            [p["beta_peak"]],
            [p["chi_peak"]],
            xerr=[p["beta_err"]],
            yerr=[p["chi_err"]],
            fmt="*",
            ms=14,
            color=color,
            ecolor=color,
            elinewidth=1.2,
            capsize=3,
            zorder=6,
        )

        ax.annotate(
            f"L={L}\n$\\beta_c(L)$={p['beta_peak']:.4f}",
            xy=(p["beta_peak"], p["chi_peak"]),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=9,
            color=color,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=color, alpha=0.9),
        )

    ax.axvline(np.log(3.0), color="black", ls="--", lw=1.2, label=r"$\beta_c(\infty)=\ln 3$")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"Susceptibility $\chi$")
    ax.set_title("q=4 Potts Susceptibility Curves with Labeled Peaks")
    ax.grid(alpha=0.25)
    ax.legend(frameon=True)
    fig.tight_layout()

    png_path = out_dir / "susceptibility_peaks.png"
    pdf_path = out_dir / "susceptibility_peaks.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
