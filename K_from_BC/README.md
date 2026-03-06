# Couplings from boundary conditions (K_from_BC)

## Can `examples/ising_flat_crit.cc` be adapted for this idea?

Short answer: **yes, as a strong starting point**.

`ising_flat_crit.cc` already has the core ingredients you need:

- independent directional couplings `k1`, `k2`, `k3` exposed as run-time options
- automatic critical-point estimate via `find_crit(k1,k2,k3)`
- directional two-point correlator measurements (`corr_x`, `corr_y`, `corr_z`, etc.)
- finite-size and thermodynamic observables (`U4`, susceptibility) for consistency checks

So the code is already close to a "coupling calibration" workflow. What is still missing is
explicit boundary-shape control and a fitting loop that maps target geometry
(lengths/aspect ratio) to couplings.

## Recommended adaptation path

1. **Introduce explicit skew/shape boundary parameters** (torus period vectors)
   in lattice initialization, so boundary geometry can be scanned directly.
2. **Measure correlation lengths / decay rates by direction** from the existing
   correlators and build an affine-rescaling diagnostic.
3. **Fit (`k1`,`k2`,`k3`) to target shape data** with multi-observable matching
   (directional correlators + Binder + susceptibility + finite-size scaling).
4. **Validate with collapse tests** after affine normalization at multiple `N`.

## Important caveat

The affine transformation gives continuum anisotropy information, but generally
not a unique microscopic `k_i = f(l_1,l_2,l_3)` by itself. Treat the method as a
calibration procedure, not a one-shot algebraic inversion.

## Relevant hooks in `ising_flat_crit.cc`

- Coupling inputs: `--k1`, `--k2`, `--k3`
- Critical beta from couplings: `find_crit(k1,k2,k3)`
- Directional correlators for anisotropy diagnostics:
  - `corr_x`, `corr_y`, `corr_z`
  - cross/derived channels `corr_w`, `corr_xz`, `corr_yz`
  - zero-momentum projected channels `corr_zero_x`, `corr_zero_y`, `corr_zero_z`

These make `ising_flat_crit.cc` a practical base for a new `K_from_BC`
analysis pipeline.

## Prototype driver in this folder

- `K_from_BC/ising_flat_crit_k_from_bc.cc` is a copy/adaptation of
  `examples/ising_flat_crit.cc`.
- New shape inputs: `--theta1`, `--theta2`, `--n_x`, `--n_y`.
- The code includes a local helper block that builds an oblique torus basis from
  `(theta1, theta2)` and uses that basis for minimum-image distance binning
  (`corr_skew_radial`).
- No `lattice.h` changes are required; the helper logic is embedded directly in
  the `.cc` file as requested.

- Two-point direction plotting helper: `K_from_BC/plot_two_point_directions.py`
  (plots x, y, and oblique `(x-y)` channels from the output `.dat`).

## Aspect=1 scan (fixed `k1`, vary `k2=k3`)

- Use `K_from_BC/sweep_aspect1_k_ratio.py` for the focused scan requested here.
- It fixes the target isosceles geometry to `height/base = 1` via
  `theta1 = theta2 = atan(2)`.
- It holds an overall scale `K = k1` fixed and scans a ratio list where
  `k2/k1 = k3/k1`.
- For each ratio it runs `ising_flat_crit_k_from_bc.cc`, writes the raw `.dat`,
  extracts a compact x/y `.dat` file for local plotting, and generates **PDF**
  plots (optionally PNG) showing `corr_x` and `corr_y`.
- It also prints midpoint-bin matching diagnostics (`corr_y/corr_x`) to support
  the coupling-tuning criterion.

Example:

```bash
r### g++ -std=c++14 -O3 -Wall -Wno-sign-compare -I include \
  K_from_BC/ising_flat_crit_k_from_bc.cc -o /tmp/ising_kbc_aspect1
python K_from_BC/sweep_aspect1_k_ratio.py --exe /tmp/ising_kbc_aspect1
# optional: also emit PNG alongside PDF
python K_from_BC/sweep_aspect1_k_ratio.py --exe /tmp/ising_kbc_aspect1 --formats pdf,png
```

# print the extracted x/y correlator files created by the sweep
for f in K_from_BC/aspect1_scan_xy_dat/*.dat; do
  echo "==== $f ===="
  cat "$f"
done
```

The extracted files use columns:

`r corr_x err_x corr_y err_y`

## Triangle-side twisted-BC symmetry test

New executable:

- `K_from_BC/ising_flat_crit_k_from_bc_triangle_sides.cc`

Goal:

- Test only the two side directions of the imposed right-triangle geometry.
- Use twisted boundary conditions so these two real-geometry directions are
  periodic and directly comparable.

Twisted BC convention:

- `(x + Nx, y) ~ (x, y + twist_shift)`

Measured channels in output typed `.dat`:

- `0`: raw `x` direction
- `1`: raw oblique side direction (the imposed real-`y` direction)
- `4`: connected `x`
- `5`: connected oblique side

Build and run example:

```bash
g++ -std=c++14 -O3 -Wall -Wno-sign-compare -I include \
  K_from_BC/ising_flat_crit_k_from_bc_triangle_sides.cc \
  -o /tmp/ising_kbc_triangle_sides

/tmp/ising_kbc_triangle_sides \
  --n_x 48 --n_y 64 \
  --lx 1.0 --ly 0.75 --use_sinh_rule 1 \
  --twist_shift 48 --oblique_sign 1 \
  --n_therm 2000 --n_traj 20000 --n_skip 20
```

This mode computes couplings from a right-triangle `sinh(2k)=l*/l` mapping via
`--lx --ly` and reports `symmetry_check_rms_z_connected` in stdout.

Plot helper:

- `K_from_BC/plot_triangle_side_directions.py`

Example:

```bash
python K_from_BC/plot_triangle_side_directions.py \
  --typed ising_flat_crit_k_from_bc/triangle_sides_twist_Nx48_Ny64_s48_ob1_lx1.000_ly0.750_k0.347_0.281_0.575.dat \
  --output K_from_BC/triangle_side_twist_connected.png
```


## Downloadable results bundle

To package outputs for download/share (PDF/PNG/log/dat + manifest), run:

```bash
bash K_from_BC/package_results.sh
```

This writes `k_from_bc_results.tar.gz` at the repository root.
You can choose a custom filename:

```bash
bash K_from_BC/package_results.sh my_results.tar.gz
```