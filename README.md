# Quantum Finite Elements

Code for lattice simulations using the quantum finite element method.

We currently have Monte Carlo routines for ising and phi4 field theory,
including Metropolis and Wolff cluster updates, and overrelaxation for phi4.
The lattice can be defined on the following manifolds:

- Flat triangular and rectangular lattices with periodic boundary conditions.
- AdS2 lattice with open boundary conditions.
- S2 lattice from Platonic solids triangulated and projected onto a unit sphere.
- S3 lattice from 5-cell, 16-cell, and 600-cell projected onto a unit sphere.
- Any of the above lattices with an orthogonal extra dimension.

## Notes

- Coupling-from-boundary-condition adaptation notes: `K_from_BC/README.md`.
