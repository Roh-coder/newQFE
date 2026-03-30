# K_from_BC twisted-parallelogram 3d extension

This subdirectory is reserved for the `K_from_BC` version of the twisted
parallelogram project with an appended third direction (`N_t`).

Intended scope:

- drivers adapted from the existing `K_from_BC` executables for `N_t > 1`
- scripts that analyze or tune couplings for the extruded geometry
- run notes, plots, and small utilities specific to the 3d extension

Current reference point:

- simulation and plotting work for this geometry currently lives under
  `thinkDoubleTwist/`

Suggested next files:

- a 3d-capable driver source file
- a local Makefile or run script
- a plotting script for the equal-time or full 3d two-point function