# Picking 10 Geometry Benchmarks

## Goal

Pick 10 target triangle shapes with simple integer side ratios, approximate each one by a small twisted lattice `(Lx, Ly, Tx, Ty)`, and record the exact realized side ratios from that twisted cell. Those realized ratios can then be converted into untwisted anisotropic couplings using the exact 2D triangular Ising critical-surface rule.

## Important Constraints

- I excluded `1-1-2` because it is degenerate: `1 + 1 = 2`, so it does not define a genuine triangle.
- I restricted the target list to acute triangles. With the conventions already used in this repo, acute triangles give positive ferromagnetic untwisted couplings. Right or obtuse targets would drive one of the dimensionless critical couplings to zero or negative values.
- The side-orientation bug from the `4-5-6` case still matters here. The untwisted branch has to use the `embedding_cycles` that match the same physical sides used by the twisted branch.

## Geometry Conventions

For a twisted cell `(Lx, Ly, Tx, Ty)`, define the three cycle vectors in the triangular basis as

```text
c0 = Lx + Ty * omega
c1 = Tx - Ly * omega
c2 = -(Lx + Tx) + (Ly - Ty) * omega

omega = 1/2 + i * sqrt(3) / 2
```

The three realized side lengths are then

```text
s0 = |c0|,  s1 = |c1|,  s2 = |c2|
```

For each target ratio, I searched over small twisted cells with

```text
Lx, Ly in [4, 20]
Tx in [-(Lx - 1), ..., Lx - 1]
Ty in [0, ..., Ly - 1]
```

and minimized the mismatch between the normalized sorted side lengths and the target integer ratio.

## Coupling Conversion

Let the sorted realized side lengths be

```text
A <= B <= C
```

with opposite angles `alpha`, `beta`, `gamma`. In the untwisted anisotropic conventions already used here,

- `k2` is paired with the smallest side `A`
- `k1` is paired with the middle side `B`
- `k3` is paired with the largest side `C`

The exact critical products are then defined by the sinh rule

```text
sinh(2 * beta_c * k2) = cot(alpha)
sinh(2 * beta_c * k1) = cot(beta)
sinh(2 * beta_c * k3) = cot(gamma)
```

Equivalently,

```text
beta_c * k_i = -1/2 * log(tan(theta_i / 2))
```

where `theta_i` is the matching angle.

In the tables below I normalize the untwisted couplings by setting

```text
k3 = 1
r1 = k1 / k3
r2 = k2 / k3
beta = beta_c
```

## Orientation Handling

For the untwisted `(L, L, 0, 0)` branch in this workflow, the physical-side ordering is fixed as

```text
untwisted cycle 0 -> middle side
untwisted cycle 1 -> shortest side
untwisted cycle 2 -> longest side
```

So after sorting the twisted side lengths from shortest to longest, the untwisted `embedding_cycles` must be chosen to match the same two physical sides that the twisted default basis `(cycle 0, cycle 1)` is using.

When two sides are equal, the listed `embedding_cycles` is only a representative choice; swapping equal sides is harmless.

## Recommended 10-Geometry Set

| target integer ratio | twisted cell `(Lx,Ly,Tx,Ty)` | realized sorted lengths | realized ratio | max ratio error | twisted cycle order short->long | untwisted `embedding_cycles` | `beta_c` | `r1` | `r2` |
| --- | --- | --- | --- | ---: | --- | --- | ---: | ---: | ---: |
| `1-1-1` | `(4, 4, 0, 0)` | `(4.0000, 4.0000, 4.0000)` | `(1.0000, 1.0000, 1.0000)` | `0.000%` | `[0, 1, 2]` | `[1, 0]` | `0.274653072167` | `1.000000000000` | `1.000000000000` |
| `1-2-2` | `(8, 5, 2, 1)` | `(4.3589, 8.5440, 8.7178)` | `(1.0000, 1.9601, 2.0000)` | `1.994%` | `[1, 0, 2]` | `[0, 1]` | `0.109103599502` | `1.365636502427` | `6.159868968409` |
| `2-3-3` | `(4, 4, -1, 1)` | `(3.0000, 4.5826, 4.5826)` | `(1.0000, 1.5275, 1.5275)` | `1.835%` | `[2, 0, 1]` | `[0, 2]` | `0.169915682051` | `1.000000000000` | `3.119581887046` |
| `3-3-4` | `(4, 6, 0, 1)` | `(4.5826, 4.5826, 6.0000)` | `(1.0000, 1.0000, 1.3093)` | `1.802%` | `[0, 2, 1]` | `[1, 2]` | `0.071920518113` | `5.446287367229` | `5.446287367229` |
| `3-4-4` | `(5, 7, 1, 0)` | `(5.0000, 6.5574, 6.5574)` | `(1.0000, 1.3115, 1.3115)` | `1.638%` | `[0, 1, 2]` | `[1, 0]` | `0.200758621832` | `1.000000000000` | `2.206078057497` |
| `4-4-5` | `(5, 7, 0, 1)` | `(5.5678, 5.5678, 7.0000)` | `(1.0000, 1.0000, 1.2572)` | `0.579%` | `[0, 2, 1]` | `[1, 2]` | `0.106416953856` | `3.472796597450` | `3.472796597450` |
| `4-5-5` | `(6, 8, 1, 0)` | `(6.0000, 7.5498, 7.5498)` | `(1.0000, 1.2583, 1.2583)` | `0.664%` | `[0, 1, 2]` | `[1, 0]` | `0.210254830503` | `1.000000000000` | `1.990413763108` |
| `5-5-6` | `(7, 9, 0, 1)` | `(7.5498, 7.5498, 9.0000)` | `(1.0000, 1.0000, 1.1921)` | `0.660%` | `[2, 0, 1]` | `[0, 2]` | `0.148995858027` | `2.305372163697` | `2.305372163697` |
| `4-5-6` | `(6, 6, 3, 1)` | `(5.1962, 6.5574, 7.8102)` | `(1.0000, 1.2620, 1.5031)` | `0.958%` | `[1, 0, 2]` | `[0, 1]` | `0.066421804006` | `4.702782819756` | `7.353910143333` |
| `5-6-7` | `(4, 5, 2, 3)` | `(4.3589, 5.2915, 6.0828)` | `(1.0000, 1.2140, 1.3955)` | `1.163%` | `[1, 2, 0]` | `[2, 1]` | `0.110136103641` | `2.666681064033` | `4.069832921022` |

## Notes On This Set

- This list is intentionally biased toward acute triangles with moderate anisotropy, because those are the most straightforward next-step benchmarks for the twisted-vs-untwisted continuum check.
- The realized ratios, not the original target integers, are what should be fed into the untwisted branch. The twisted cell only approximates the target integer triangle; the untwisted couplings should match the exact realized twisted geometry.
- The strongest-shape case in this set is still `4-5-6`, which remains useful as a stress benchmark because it forces a clear side-ordering test.
- If we want an even broader scan later, the next natural extension is to add skinnier acute isosceles targets such as `1-3-3` or `1-4-4`, but those push `r2` much larger and are a worse first production pass.

## Immediate Next Step

Build a 10-benchmark campaign config directly from this table:

- use the listed twisted cells as `modular_target_geometry`
- use isotropic twisted couplings `(r1, r2) = (1, 1)` at `beta = ln(3) / 4`
- use the listed `r1`, `r2`, `beta_c`, and `embedding_cycles` for the untwisted branch

That gives a clean first pass for the "same continuum?" check across 10 distinct geometries while explicitly handling the side-orientation issue.

## Twisted-Reference 4x Rerun

For the coupling-to-geometry workflow, the more relevant comparison is not a symmetric common-grid score. The intended production-style check is:

- build a denser twisted reference family
- fit its continuum manifold first
- interpolate that twisted continuum manifold onto the untwisted continuum points
- use only those twisted-on-untwisted residuals to decide whether the untwisted candidate matches

For that rerun, I keep the untwisted size ladder at the base multiples `{1,2,3,4,5}` of the `12x12` untwisted cell, so the untwisted family remains `(12,24,36,48,60)`.

The twisted family should also use base multiples `{1,2,3,4,5}`, but with a base cell chosen as the nearest integer multiple of the original pickGeo twisted cell whose volume is approximately `4 x 12 x 12 = 576`. Because the allowed twisted cells are restricted to integer multiples of the base geometry, that `4x` target is approximate rather than exact for some shapes.

| benchmark | original twisted cell | chosen scale to define twisted reference base | twisted reference base `(Lx,Ly,Tx,Ty)` | base volume ratio vs untwisted `12x12` |
| --- | --- | ---: | --- | ---: |
| `1-1-1` | `(4,4,0,0)` | `6` | `(24,24,0,0)` | `4.000` |
| `1-2-2` | `(8,5,2,1)` | `4` | `(32,20,8,4)` | `4.667` |
| `2-3-3` | `(4,4,-1,1)` | `6` | `(24,24,-6,6)` | `3.750` |
| `3-3-4` | `(4,6,0,1)` | `5` | `(20,30,0,5)` | `4.167` |
| `3-4-4` | `(5,7,1,0)` | `4` | `(20,28,4,0)` | `3.889` |
| `4-4-5` | `(5,7,0,1)` | `4` | `(20,28,0,4)` | `3.889` |
| `4-5-5` | `(6,8,1,0)` | `3` | `(18,24,3,0)` | `3.000` |
| `5-5-6` | `(7,9,0,1)` | `3` | `(21,27,0,3)` | `3.938` |
| `4-5-6` | `(6,6,3,1)` | `4` | `(24,24,12,4)` | `4.333` |
| `5-6-7` | `(4,5,2,3)` | `5` | `(20,25,10,15)` | `4.514` |

The corresponding production config is `responsible_method_tests/configs/raw_manifold_fss_pickgeo10_twisted_reference4x_20260521.json`, and the matching analysis mode is `align_manifolds.py --comparison-mode twisted_reference`.