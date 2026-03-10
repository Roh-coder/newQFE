# Ny=2Nx Collaboration Workflow (Clean Reproducible Package)

This folder is a collaborator-facing, minimal workflow for the directional-correlation comparison study.

It runs **three curated cases** on the same lattice shape (`Nx=32, Ny=64`) and produces a single clean figure of connected correlators vs fractional separation `r/L`.

## What this package runs

Included C++ source entry point in this folder:

- `K_from_BC/collab_ny2nx/ising_flat_crit_k_from_bc.cc`
- `K_from_BC/collab_ny2nx/k_from_bc_collab_ny2nx.cc` (wrapper)

The full source is now copied into this folder and is the default file
compiled by `run_cases.sh`.

All runs use:

- triangular-lattice code: `K_from_BC/ising_flat_crit_k_from_bc.cc`
- geometry angles: `theta1=theta2=pi/3`
- default stats: `n_therm=5000`, `n_traj=120000`, `n_skip=20`, `n_wolff=3`, `n_metropolis=5`

Cases:

1. `equilateral`
- `k1=k2=k3=0.27465307216702745`

2. `right_45_45_90_k3zero`
- `k1=k2=0.4406867935097715`, `k3=0`

3. `right_h_over_b_half_k3zero`
- right triangle with `height/base = 1/2` under sinh-rule mapping
- `k1=0.24060591252980174`, `k2=0.7218177375894052`, `k3=0`

## Quick start

From repository root:

```bash
bash K_from_BC/collab_ny2nx/run_all.sh
```

Outputs appear in:

- `K_from_BC/collab_ny2nx/output/manifest.json`
- `K_from_BC/collab_ny2nx/output/logs/*.log`
- `K_from_BC/collab_ny2nx/output/results/runs/*.dat`
- `K_from_BC/collab_ny2nx/output/results/plots/correlators_fractional.png`

## Run scripts separately

```bash
bash K_from_BC/collab_ny2nx/run_cases.sh
python3 K_from_BC/collab_ny2nx/plot_cases.py
```

Useful options:

```bash
bash K_from_BC/collab_ny2nx/run_cases.sh --nx 32 --ny 64 --n-traj 120000
bash K_from_BC/collab_ny2nx/run_cases.sh --skip-compile
python3 K_from_BC/collab_ny2nx/plot_cases.py --output K_from_BC/collab_ny2nx/output/results/plots/custom.png
```

## Direct compile/run (single case data generation)

If collaborators want raw `.dat` output directly from the C++ file:

```bash
g++ -std=c++14 -O3 -Wall -Wno-sign-compare -I include \
  K_from_BC/collab_ny2nx/ising_flat_crit_k_from_bc.cc \
  -o K_from_BC/collab_ny2nx/bin/ising_kbc_direct

K_from_BC/collab_ny2nx/bin/ising_kbc_direct \
  --n_x 32 --n_y 64 \
  --theta1 1.0471975512 --theta2 1.0471975512 \
  --k1 0.27465307216702745 --k2 0.27465307216702745 --k3 0.27465307216702745 \
  --n_therm 5000 --n_traj 120000 --n_skip 20 --n_wolff 3 --n_metropolis 5
```

This writes a typed data file in:

- `ising_flat_crit_k_from_bc/`

For the command above, the filename is:

- `ising_flat_crit_k_from_bc/32_64_t1_1.047_t2_1.047_k_0.275_0.275_0.275.dat`

## Plot conventions

- plotted channels are connected correlators:
  - type `4`: `x` (base)
  - type `5`: `left/y`
  - type `7`: `right` (`y-x`)
- x-axis is fractional separation `r/L`
- normalizations used:
  - `Lx = Nx`
  - `Ly = Ny`
  - `Lright = lcm(Nx, Ny)`
- for `Nx=32, Ny=64`, this means `Lright=64` (full periodic right-direction orbit)

## Reproducibility notes

- `run_cases.sh` compiles `K_from_BC/collab_ny2nx/ising_flat_crit_k_from_bc.cc`
  into a local binary at `K_from_BC/collab_ny2nx/bin/ising_kbc`.
- Raw run output data still originates from the code's standard output folder:
  - `ising_flat_crit_k_from_bc/*.dat`
- A copy of each case data file is saved into this package under:
  - `K_from_BC/collab_ny2nx/output/results/runs/`

This makes collaborator handoff straightforward: point collaborators to this folder and `run_all.sh`.
