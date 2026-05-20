# Results

Store validation artifacts here, grouped by benchmark.

## Suggested Layout

- `B0_equilateral_control/`
- `B1_quarter_twist/`
- `B2_456_stress/`

Within each benchmark directory, keep the following together:

- copied or linked manifests used for the comparison
- aligned-manifold plots
- line-cut plots
- score tables
- one short summary note with the pass/fail decision

The purpose of this folder is reproducibility. A collaborator should be able to
open a benchmark subdirectory and see exactly which inputs produced the stated
comparison result.