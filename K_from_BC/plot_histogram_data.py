#!/usr/bin/env python3
"""
Visualize histogram reweighting data from triangular Ising simulations.
Creates histograms of magnetization, m^2, and Binder cumulant by lattice size.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_data(filepath):
    """Load histogram data file and return data organized by (N, beta) pairs"""
    data_by_N_beta = {}
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split()
            N = int(parts[0])
            beta = float(parts[1])
            m = float(parts[2])
            m2 = float(parts[3])
            action = float(parts[4])
            
            key = (N, beta)
            if key not in data_by_N_beta:
                data_by_N_beta[key] = {"m": [], "m2": [], "action": []}
            
            data_by_N_beta[key]["m"].append(m)
            data_by_N_beta[key]["m2"].append(m2)
            data_by_N_beta[key]["action"].append(action)
    
    return data_by_N_beta


def compute_binder(m2_array, m_array):
    """Compute Binder cumulant U4 for array of configurations"""
    m4 = m2_array ** 2
    # U4 = 1 - m^4 / (3 * m2^2)
    # Avoid division by zero for disordered phase
    u4 = np.ones_like(m2_array)
    nonzero = m2_array > 1e-10
    u4[nonzero] = 1.0 - m4[nonzero] / (3.0 * m2_array[nonzero] ** 2)
    return u4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="K_from_BC/histogram_data.dat",
        help="Input data file"
    )
    parser.add_argument(
        "--output",
        default="K_from_BC/histogram_plots.png",
        help="Output plot file"
    )
    args = parser.parse_args()
    
    data_by_N_beta = load_data(args.input)
    
    # Organize by size and beta
    sizes = sorted(set(N for N, beta in data_by_N_beta.keys()))
    betas = sorted(set(beta for N, beta in data_by_N_beta.keys()))
    
    print(f"Lattice sizes: {sizes}")
    print(f"Beta values: {len(betas)} points from {min(betas):.4f} to {max(betas):.4f}")
    print(f"Total (N,beta) pairs: {len(data_by_N_beta)}")
    
    # Create massive figure with all histograms
    # 3 rows (m, m2, U4) x len(sizes) x len(betas) plots
    fig, axes = plt.subplots(
        3 * len(sizes),
        len(betas),
        figsize=(3 * len(betas), 4 * len(sizes))
    )
    
    for size_idx, N in enumerate(sizes):
        for beta_idx, beta in enumerate(betas):
            key = (N, beta)
            if key not in data_by_N_beta:
                continue
            
            m = np.array(data_by_N_beta[key]["m"])
            m2 = np.array(data_by_N_beta[key]["m2"])
            u4 = compute_binder(m2, m)
            
            # Row for m histograms
            ax = axes[3 * size_idx, beta_idx]
            n_bins = max(5, min(30, len(np.unique(m))))
            ax.hist(m, bins=n_bins, alpha=0.7, color="blue", edgecolor="black")
            if beta_idx == 0:
                ax.set_ylabel(f"L={N}\n|m|", fontsize=8)
            if size_idx == 0:
                ax.set_title(f"β={beta:.3f}", fontsize=8)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=6)
            
            # Row for m2 histograms
            ax = axes[3 * size_idx + 1, beta_idx]
            n_bins = max(5, min(30, len(np.unique(m2))))
            ax.hist(m2, bins=n_bins, alpha=0.7, color="green", edgecolor="black")
            if beta_idx == 0:
                ax.set_ylabel(f"L={N}\nm²", fontsize=8)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=6)
            
            # Row for U4 histograms
            ax = axes[3 * size_idx + 2, beta_idx]
            u4_ordered = u4[m2 > 0.01]
            if len(u4_ordered) > 1:
                unique_vals = len(np.unique(np.round(u4_ordered, 8)))
                n_bins = max(1, min(10, unique_vals))
                try:
                    ax.hist(u4_ordered, bins=n_bins, alpha=0.7, color="red", edgecolor="black")
                    ax.set_xlim([0, 1])
                except:
                    if len(u4_ordered) > 0:
                        ax.scatter(range(len(u4_ordered)), u4_ordered, alpha=0.5, color="red", s=5)
            else:
                n_disord = len(u4) - len(u4_ordered)
                ax.text(0.5, 0.5, f"{n_disord}", ha="center", va="center", fontsize=8)
            if beta_idx == 0:
                ax.set_ylabel(f"L={N}\nU₄", fontsize=8)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=6)
            if size_idx == len(sizes) - 1:
                ax.set_xlabel(f"", fontsize=6)
            
            # Print stats
            if beta_idx == 0 and size_idx == 0:
                print(f"\n{'L':>3} {'Beta':>8} {'<m>':>8} {'<m²>':>8} {'<U₄>':>8} {'n_ord':>6}")
            print(f"{N:3d} {beta:8.4f} {np.mean(m):8.4f} {np.mean(m2):8.4f} {np.mean(u4):8.4f} "
                  f"{len(u4_ordered):6d}")
    
    plt.tight_layout()
    print(f"\nSaved comprehensive plot to: {args.output}")
    plt.savefig(args.output, dpi=100, bbox_inches="tight")
    plt.show()
    print(f"Plot dimensions: {3*len(sizes)} rows x {len(betas)} columns")


if __name__ == "__main__":
    main()
