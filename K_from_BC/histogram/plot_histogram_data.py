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


def compute_binder_block_samples(m_array, m2_array, n_blocks=20):
    """Compute block-level Binder estimates U4 = 1 - <m^4>/(3< m^2 >^2)."""
    m_array = np.asarray(m_array)
    m2_array = np.asarray(m2_array)
    n = len(m_array)
    if n < 4:
        return np.array([], dtype=float)

    # Use contiguous blocks to show a Binder estimate distribution from one run.
    n_blocks = max(2, min(n_blocks, n // 10 if n >= 20 else 2))
    block_size = n // n_blocks
    if block_size < 2:
        return np.array([], dtype=float)

    u4_blocks = []
    m4_array = m_array ** 4
    for b in range(n_blocks):
        start = b * block_size
        stop = (b + 1) * block_size if b < n_blocks - 1 else n
        if stop - start < 2:
            continue
        mean_m2 = float(np.mean(m2_array[start:stop]))
        if mean_m2 <= 1e-14:
            continue
        mean_m4 = float(np.mean(m4_array[start:stop]))
        u4 = 1.0 - mean_m4 / (3.0 * mean_m2 * mean_m2)
        u4_blocks.append(u4)

    return np.array(u4_blocks, dtype=float)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="K_from_BC/histogram/histogram_data.dat",
        help="Input data file"
    )
    parser.add_argument(
        "--output",
        default="K_from_BC/histogram/histogram_plots.png",
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
            u4_blocks = compute_binder_block_samples(m, m2)
            
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
            
            # Row for U4 histograms (block-level Binder estimates)
            ax = axes[3 * size_idx + 2, beta_idx]
            u4_finite = u4_blocks[np.isfinite(u4_blocks)]
            if len(u4_finite) > 1:
                n_bins = max(4, min(20, len(u4_finite) // 2))
                try:
                    ax.hist(u4_finite, bins=n_bins, alpha=0.7, color="red", edgecolor="black")
                    ax.set_xlim([0, 1])
                except:
                    if len(u4_finite) > 0:
                        ax.scatter(range(len(u4_finite)), u4_finite, alpha=0.5, color="red", s=5)
            else:
                ax.text(0.5, 0.5, "insufficient\nblocks", ha="center", va="center", fontsize=8)
            if beta_idx == 0:
                ax.set_ylabel(f"L={N}\nU₄", fontsize=8)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=6)
            if size_idx == len(sizes) - 1:
                ax.set_xlabel(f"", fontsize=6)
            
            # Print stats
            if beta_idx == 0 and size_idx == 0:
                print(f"\n{'L':>3} {'Beta':>8} {'<m>':>8} {'<m²>':>8} {'<U₄>':>8} {'n_blk':>6}")
            mean_u4 = np.mean(u4_finite) if len(u4_finite) > 0 else float("nan")
            print(f"{N:3d} {beta:8.4f} {np.mean(m):8.4f} {np.mean(m2):8.4f} {mean_u4:8.4f} "
                  f"{len(u4_finite):6d}")
    
    plt.tight_layout()
    print(f"\nSaved comprehensive plot to: {args.output}")
    plt.savefig(args.output, dpi=100, bbox_inches="tight")
    plt.show()
    print(f"Plot dimensions: {3*len(sizes)} rows x {len(betas)} columns")


if __name__ == "__main__":
    main()
