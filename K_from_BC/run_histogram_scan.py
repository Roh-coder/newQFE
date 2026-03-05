#!/usr/bin/env python3
"""
Histogram reweighting data collection for equilateral triangular Ising model.

Compiles and runs tri_histogram_reweight.cc across multiple lattice sizes
and beta ranges, collecting configuration-level magnetization data for later
histogram reweighting analysis to determine the critical point.

Output format (one per configuration):
  N beta m m2 action
where m is |<s>| and m2 is m^2, allowing posterior computation of binder
cumulant U4 = 1 - m^4/(3*m2^2).
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def compile_program(workspace_root):
    """Compile tri_histogram_reweight.cc"""
    exe = "/tmp/tri_histogram_reweight"
    source = workspace_root / "K_from_BC" / "tri_histogram_reweight.cc"
    
    print(f"Compiling {source}...")
    result = subprocess.run(
        [
            "g++",
            "-std=c++14",
            "-O3",
            "-Wall",
            "-Wno-sign-compare",
            f"-I{workspace_root / 'include'}",
            str(source),
            "-o",
            exe,
        ],
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        print(f"Compilation failed:\n{result.stderr}")
        sys.exit(1)
    
    print(f"Compiled to {exe}")
    return exe


def run_simulation(exe, N, beta_min, beta_max, n_beta, n_therm, n_traj, n_skip, n_wolff, n_metropolis):
    """Run a single simulation and return raw output"""
    cmd = [
        exe,
        "--N", str(N),
        "--beta_min", str(beta_min),
        "--beta_max", str(beta_max),
        "--n_beta", str(n_beta),
        "--n_therm", str(n_therm),
        "--n_traj", str(n_traj),
        "--n_skip", str(n_skip),
        "--n_wolff", str(n_wolff),
        "--n_metropolis", str(n_metropolis),
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Simulation failed:\n{result.stderr}")
        sys.exit(1)
    
    return result.stdout


def parse_dat_files(output_dir):
    """Parse all .dat files from a simulation run and return array of (N, beta, m, m2, action)"""
    data = []
    dat_dir = Path(output_dir) / "tri_histogram_reweight"
    
    if not dat_dir.exists():
        print(f"Warning: {dat_dir} not found")
        return data
    
    for dat_file in sorted(dat_dir.glob("*.dat")):
        # Extract N and beta from filename: N16_beta_0.2500.dat
        stem = dat_file.stem
        parts = stem.split("_")
        N = int(parts[0][1:])  # Remove 'N' prefix
        beta = float(parts[2])  # beta value is after "beta_"
        
        with open(dat_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                try:
                    m, m2, action = map(float, line.split())
                    data.append((N, beta, m, m2, action))
                except ValueError:
                    continue
    
    return data


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--workspace",
        default="/workspaces/newQFE",
        help="Workspace root directory"
    )
    parser.add_argument(
        "--output",
        default="K_from_BC/histogram_data.dat",
        help="Output data file (relative to workspace)"
    )
    parser.add_argument(
        "--sizes",
        default="16,24,32",
        help="Comma-separated lattice sizes to scan"
    )
    parser.add_argument(
        "--beta_min",
        type=float,
        default=0.246,
        help="Minimum beta value"
    )
    parser.add_argument(
        "--beta_max",
        type=float,
        default=0.274,
        help="Maximum beta value"
    )
    parser.add_argument(
        "--n_beta",
        type=int,
        default=8,
        help="Number of beta steps"
    )
    parser.add_argument(
        "--n_therm",
        type=int,
        default=1000,
        help="Thermalization steps"
    )
    parser.add_argument(
        "--n_traj",
        type=int,
        default=5000,
        help="Trajectory steps"
    )
    parser.add_argument(
        "--n_skip",
        type=int,
        default=10,
        help="Measurement skip interval"
    )
    parser.add_argument(
        "--n_wolff",
        type=int,
        default=3,
        help="Wolff cluster updates per step"
    )
    parser.add_argument(
        "--n_metropolis",
        type=int,
        default=5,
        help="Metropolis updates per step"
    )
    
    args = parser.parse_args()
    
    workspace = Path(args.workspace)
    output_path = workspace / args.output
    sizes = [int(s.strip()) for s in args.sizes.split(",")]
    
    print("=" * 70)
    print("Equilateral Triangular Ising - Histogram Reweighting Data Collection")
    print("=" * 70)
    print(f"Workspace: {workspace}")
    print(f"Output: {output_path}")
    print(f"Lattice sizes: {sizes}")
    print(f"Beta range: [{args.beta_min}, {args.beta_max}] with {args.n_beta} points")
    print(f"Thermalization: {args.n_therm}, Trajectories: {args.n_traj}, Skip: {args.n_skip}")
    print()
    
    # Compile once
    os.chdir(workspace)
    exe = compile_program(workspace)
    print()
    
    # Collect all data
    all_data = []
    
    for N in sizes:
        print(f"\nRunning lattice size N={N}...")
        try:
            stdout = run_simulation(
                exe,
                N,
                args.beta_min,
                args.beta_max,
                args.n_beta,
                args.n_therm,
                args.n_traj,
                args.n_skip,
                args.n_wolff,
                args.n_metropolis,
            )
            print(stdout)
            
            # Parse generated data files
            data = parse_dat_files(".")
            all_data.extend(data)
            
        except Exception as e:
            print(f"Error running N={N}: {e}")
            continue
    
    # Write master data file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w") as f:
        f.write("# Histogram reweighting data for equilateral triangular Ising\n")
        f.write("# Format: N beta m m2 action\n")
        f.write("# where m = |<s>|, m2 = m^2\n")
        f.write("# Binder cumulant (posterior): U4 = 1 - m^4 / (3*m2^2)\n")
        f.write("# Magnetic susceptibility (posterior): chi = (m2 - m^2) * N^2\n")
        f.write("#\n")
        f.write(f"# Data collection parameters:\n")
        f.write(f"# beta_min={args.beta_min}, beta_max={args.beta_max}, n_beta={args.n_beta}\n")
        f.write(f"# n_therm={args.n_therm}, n_traj={args.n_traj}, n_skip={args.n_skip}\n")
        f.write(f"# n_wolff={args.n_wolff}, n_metropolis={args.n_metropolis}\n")
        f.write("#\n")
        f.write("# Data columns:\n")
        f.write("#   N      : lattice size\n")
        f.write("#   beta   : inverse temperature\n")
        f.write("#   m      : absolute magnetization per spin\n")
        f.write("#   m2     : m squared\n")
        f.write("#   action : action per spin\n")
        f.write("#\n")
        
        for N, beta, m, m2, action in sorted(all_data):
            f.write(f"{N:4d} {beta:.12e} {m:.12e} {m2:.12e} {action:.12e}\n")
    
    print()
    print("=" * 70)
    print(f"Data collection complete!")
    print(f"Total configurations: {len(all_data)}")
    print(f"Output written to: {output_path}")
    print()
    print("File contents (first 20 lines):")
    with open(output_path) as f:
        for i, line in enumerate(f):
            if i >= 20:
                break
            print(line.rstrip())
    
    if len(all_data) > 20:
        print(f"... ({len(all_data) - 20} more configurations)")
    
    print()
    print("Ready for histogram reweighting analysis.")
    print("Usage example:")
    print("  python analyze_histogram.py --input", output_path)


if __name__ == "__main__":
    main()
