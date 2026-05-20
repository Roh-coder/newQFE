# First Validation Campaign

This is the first-pass runbook for collaborator briefing and internal execution.

## Objective

Demonstrate, with the smallest convincing benchmark set, that the modular,
twisted, and untwisted-anisotropic routes can be made continuum-equivalent at
the level of the full all-to-all correlator manifold, and that nearby wrong
points are score-separable.

## Execution Order

1. Run B0 equilateral control.
2. Freeze the alignment and score definitions after B0 passes.
3. Run B1 quarter-twist benchmark with a tuned untwisted anisotropic target.
4. Run B2 4-5-6 stress benchmark only after B1 is stable.

## Data Products Needed

- Modular manifold samples on a common normalized cell.
- Twisted-lattice size-family outputs with continuum extrapolation metadata.
- Untwisted anisotropic size-family outputs with continuum extrapolation
  metadata.
- One manifest per benchmark collecting the three matched routes and the local
  deformation stencil.

## Minimum Figure Set

1. One figure for B0 showing matched surfaces and the difference surface.
2. One figure for B1 showing the same comparison for a nontrivial torus shape.
3. One figure for B2 showing whether the strong-shape case still collapses.
4. One score-separation figure per benchmark showing matched vs wrong points.

## Decision Rules

- Do not promote a score from exploratory to default use until it passes B0 and
  B1.
- Do not interpret B2 failure as physics failure until B2 is rechecked with the
  same alignment and cycle-pairing logic that passed B1.
- If a score is zero-like at the truth but does not separate local wrong points,
  it is not sufficient for downstream optimization.

## Immediate Next Artifacts

- Fill `configs/campaign_template.json` with real manifest paths for B0.
- Create the first comparison script that ingests one modular source and two
  lattice manifests.
- Save all first-run outputs under `results/B0_equilateral_control/`.