# 4-State Potts (Swendsen-Wang) Prototype

This directory contains a prototype scan for the 2D square-lattice 4-state Potts model using Swendsen-Wang cluster updates.

## Files

- `potts_sw4_scan.cc`: runs scans in inverse temperature `beta` for multiple lattice sizes `L`.
- `Makefile`: builds `../bin/potts_sw4_scan`.
- `output/`: default location for generated scan data.

## Build

```bash
cd Potts
make
```

## Run

Example moderate scan:

```bash
cd ..
./bin/potts_sw4_scan \
  --L_list 16,24,32,48,64 \
  --beta_min 1.02 \
  --beta_max 1.15 \
  --n_beta 26 \
  --n_therm 3000 \
  --n_meas 12000 \
  --n_skip 1 \
  --seed 12345
```

The exact infinite-volume value on the square lattice is `beta_c = ln(3) ~= 1.098612`.

## Outputs

- `Potts/output/susceptibility_L< L >.dat`
  - Columns: `beta m_mean m_err susceptibility chi_err`
- `Potts/output/critical_points.dat`
  - Columns: `L beta_peak beta_err chi_peak chi_err peak_grid_index`

The `beta_peak` values are pseudocritical points `beta_c(L)` from the susceptibility peak for each finite `L`.

## Quick Plot

From repo root:

```bash
python3 - <<'PY'
import glob
import numpy as np
import matplotlib.pyplot as plt

for fn in sorted(glob.glob('Potts/output/susceptibility_L*.dat')):
    data = np.loadtxt(fn)
    L = fn.split('L')[-1].split('.dat')[0]
    plt.plot(data[:,0], data[:,3], marker='o', ms=3, label=f'L={L}')

plt.xlabel(r'$\beta$')
plt.ylabel(r'$\chi$')
plt.title('q=4 Potts susceptibility scan')
plt.legend()
plt.tight_layout()
plt.show()
PY
```
