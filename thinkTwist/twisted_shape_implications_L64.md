# Twisted-BC Shape Compensation: What The Sweep Means

## Scope
This note summarizes the full shift sweep for L = 64 from:
- thinkTwist/twisted_shape_sweep_L64.csv
- thinkTwist/twisted_shape_sweep_L64.png

Model setup:
- base direction = x
- left direction = y
- right direction = y - x
- twisted boundary condition: (i + L, j) ~ (i, j), (i, j + L) ~ (i + shift, j)
- compensate local triangle shape so folded half-orbit physical scales match square targets

Targets:
- x : y : (y-x) physical half-orbit scales = 1 : 1 : sqrt(2)

## Main Numerical Findings
For L = 64:
- p_x is always 64
- p_y and p_(y-x) are always equal to each other
- p_y takes values 64, 128, 256, 512, 1024, 2048, or 4096 depending on shift

Representative cases:
- shift = 0: right-isosceles local shape (theta = 90 deg)
- shift = 32 (half-twist): |y| = 0.5, |y-x| = 0.7071, theta ~= 41.41 deg
- shift = 1 (or most odd shifts): |y| = 1/64, |y-x| ~= 0.0221, theta = 0 deg (degenerate limit)

## Why Theta Was Often Zero (Important Correction)
In the first pass, theta was computed after clamping cos(theta) into [-1, 1].
That made many impossible cases appear as theta = 0 deg.

In the reanalysis, we now mark these as infeasible instead of assigning a fake angle.

For L = 64 sweep (shift = 0..63):
- feasible exact Euclidean triangle cases: 2 / 64
- feasible shifts: 0 and 32

So most shifts do not admit any exact triangle shape that satisfies all three matching constraints simultaneously.

Reanalysis files:
- thinkTwist/twisted_shape_sweep_L64_reanalysis.txt
- thinkTwist/twisted_shape_sweep_L64_reanalysis.csv
- thinkTwist/twisted_shape_sweep_L64_reanalysis.png

## Geometric Interpretation
The compensation rule forces:
- |y| = p_x / p_y
- |y-x| = sqrt(2) * p_x / p_(y-x)

Because p_x is fixed at 64 while p_y can be much larger under twist, the solver must shrink y strongly.
For many shifts, the requested |y-x| is then smaller than the triangle-inequality lower bound |x-y|,
so no Euclidean triangle exists.

That is the true reason many earlier entries looked like theta = 0 deg: they were boundary-projected artifacts,
not valid geometric solutions.

In plain terms:
- strong twists make the global topology dominate
- local metric deformation cannot fully restore a square-like continuum geometry without extreme distortion

## What This Says About The Idea
The idea is still useful, but with an important limit:

1. Local shape tuning can compensate some global period imbalance.
2. For many shifts, exact matching requires unphysical or nearly singular triangles.
3. The half-twist case is less extreme than generic odd shifts, but still far from right-isosceles.
4. Therefore, shape compensation alone is not a robust universal fix for twisted BC artifacts.

## Practical Guidance
If the goal is continuum-like square behavior in x, y, and y-x simultaneously:

1. Prefer shifts that keep p_y close to p_x (small period inflation).
2. Treat large-shift or odd-shift cases as topology-dominated regimes.
3. Use compensation as a diagnostic, not only as a prescription.
4. Expect to combine shape tuning with either:
- adjusted observables (direction-specific normalization), or
- different BC design if exact isotropy is required.

## Files Produced In This Pass
- thinkTwist/solve_compensating_triangle_cases.py
- thinkTwist/twisted_shape_sweep_L64.txt
- thinkTwist/twisted_shape_sweep_L64.csv
- thinkTwist/twisted_shape_sweep_L64.png
- thinkTwist/twisted_shape_sweep_L64.svg
