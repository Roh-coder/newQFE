# Symmetry-Constrained Triangle Coupling Fit

This note describes a local polynomial ansatz for triangle edge couplings near the equilateral triangle. The goal is to fit a shared coupling function from measured or optimized edge couplings while enforcing the discrete symmetries of the triangle.

The setup is:

- one Euclidean triangle
- vertex angles `theta_i = pi/3 + delta_i`
- angular constraint `delta_1 + delta_2 + delta_3 = 0`
- one unoriented edge coupling `K_ij = K_ji` per edge
- cyclic covariance under `1 -> 2 -> 3 -> 1`

The three edge couplings are:

- `K12`
- `K23`
- `K31`

The main object to fit is a single shared function

```text
F(s, a^2)