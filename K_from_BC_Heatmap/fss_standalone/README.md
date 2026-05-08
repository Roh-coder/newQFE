# fss_standalone — boundary-correlator FSS driver (zip-and-go)

Self-contained directory to reproduce the test/ref boundary two-point
continuum-limit plots on any laptop. No dependency on the surrounding
`newQFE` tree.

## Contents
```
Makefile                     # builds bin/ising_tri_twisted_parallelogram
src/, include/               # C++ source for the simulator
bin/                         # prebuilt Linux x86_64 binary (rebuild on macOS/Win)
lib/mc_engine.py             # β_c finder + production MC wrapper
lib/cost.py                  # boundary path geometry + tile interpolant
run_fss.py                   # the driver: MC + cache + two plots
config.example.json          # commented example with every knob exposed
```

## Quick start
```bash
make                         # rebuild the simulator (needed on macOS / Windows)
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib

cp config.example.json my_run.json
# edit my_run.json (lattice list, n_traj, β_c finder grid, …)
python3 run_fss.py --config my_run.json
```
Plots land under `results/<tag>/plots/`. MC pickles are cached under
`results/<tag>/grid/L<L>/{test,ref}/r1_*_r2_*.pkl`; rerunning the same
config is a no-op.

## All settings
Every knob is in the JSON:
- `couplings.r1, r2, k3` — bond couplings (k1=r1, k2=r2, k3=k3).
- `lattices[]` — list of `{L, test:[Lx,Ly,Tx,Ty], ref:[Lx,Ly,Tx,Ty]}`.
- `mc.{n_traj, n_skip, n_therm}` — production MC settings.
- `beta_c_finder.*` — β_c scan (Gram-Charlier 3-pass): bracket, points
  per pass, traj per pass, max bracket-translation shifts, jackknife.
- `execution.{n_workers, exe}` — parallel pool size, optional override
  path to a prebuilt simulator.
- `plots.{eighths, midpoints}` — toggle each plot.

CLI overrides: `--n-workers N`, `--n-traj N`, `--tag NAME`.

## Notes
- The Linux binary in `bin/` is shipped for convenience; `make` rebuilds
  cleanly on any POSIX platform with a C++14 compiler.
- The β_c finder runs a 3-pass GC-filtered susceptibility scan per
  lattice; each lattice is sampled at its own pseudo-β_c.
- For ref geometries with `Tx = Ty = round(L/4)`, integer aliasing
  causes outliers at L=6 (`round(1.5)=2`); use multiples of 4 if you
  want clean T/L = 1/4 across the whole sweep.
