# Collaborator Story Outline

This folder is a presentation-ready subset of the current `responsible_method_tests` evidence. The goal is to support two claims:

1. the `beta_c` finder is a real, usable critical-tuning tool rather than a noisy convenience wrapper
2. blind holdouts plus matched-geometry comparisons make a credible case that twisted-vs-untwisted correlator matching can pin couplings to continuum geometry

## Core Thesis

The project is building a geometry-inference pipeline for lattice field theory. The central claim is that a twisted isotropic lattice and an untwisted anisotropic lattice can realize the same continuum torus, and that this equivalence is visible in the full continuum-extrapolated two-point correlator manifold rather than just in a few line cuts. If that claim is right, then correlator matching becomes a practical route from couplings to geometry.

## Recommended Figure Set

1. [Matched 4-5-6 common-grid comparison](figures/01_acute456_matched_volume_common_grid.png)
Claim supported: when geometry is matched, twisted and untwisted continuum manifolds can agree within current errors.
Key numbers: `rms_z = 0.750`, `relative RMS = 0.0288`, `correlation = 0.9949`, interpretation `not distinguished within continuum point errors`.

2. [Matched 4-5-6 shared-coordinate FSS comparison](figures/02_acute456_shared_coordinate_fss.png)
Claim supported: the agreement is not just a heatmap artifact; pointwise shared-coordinate continuum intercepts also line up.
Key numbers: selected-point Z values are about `+0.31`, `-0.41`, `-0.35`, `-0.37`.

3. [Acute 4-5-6 holdout beta_c scan](figures/03_acute456_holdout_betac_scan.png)
Claim supported: the `beta_c` finder sees a real interior susceptibility peak on the actual stressed geometry.
Key numbers: held-out twisted lattice `[144,144,72,24]`, `beta_c = 0.27350859 +/- 0.00099267`, exact isotropic triangular value `0.27465307`.

4. [Partial acute 4-5-6 beta_c trend from the stopped rerun](figures/04_acute456_partial_betac_trend.png)
Claim supported: the beta finder behaves physically across the twisted size ladder and trends toward the exact large-volume value.
Key numbers: completed twisted ladder moves from `0.25699835` at `12x12_t6x2` to `0.27003786` at `72x72_t36x12`; the held-out larger lattice is `0.27350859`.

7. [Acute 4-5-6 twisted susceptibility peak drift](figures/07_acute456_twisted_peak_drift.png)
Claim supported: the raw susceptibility peak locations themselves drift monotonically toward the continuum coupling as lattice size grows.
Use this when you want a single slide that shows the beta finder on the actual scan data rather than on fitted summaries.

5. [Acute 4-5-6 literal-size blind holdout](figures/05_acute456_blind_holdout_literal_sizes.png)
Claim supported: blind extrapolation from both twisted and untwisted training families can predict an unseen large twisted datum.
Key numbers: twisted best fit `power_free`, `z = -0.43`; untwisted best fit `pade21`, `z = -0.33`; both branches stay within `|z| < 1`.

6. [1000x1000 blind holdout](figures/06_iso111_L1000_blind_holdout.png)
Claim supported: blind holdout success is not limited to the 4-5-6 acute benchmark; large out-of-sample points can still be predicted.
Key numbers: untwisted `power_free` gives `z = -0.43`; twisted `power_free` gives `z = +0.90`; the other fit families miss badly.

## Suggested 7-Slide Talk Arc

## Slide 1: The Question

Title: `Can correlators determine continuum geometry?`

Message:
The project is testing whether different microscopic lattice realizations of the same continuum torus collapse onto the same continuum two-point manifold. If they do, then geometry can be inferred by matching correlators rather than by guessing couplings directly.

Say out loud:
“This is not just a tuning problem. It is a validation problem first. We need to know whether twisted and untwisted lattices really land on the same continuum object before we trust any optimization built on top of them.”

## Slide 2: The Validation Standard

Title: `We compare full manifolds, not just line cuts`

Figure: [Matched 4-5-6 common-grid comparison](figures/01_acute456_matched_volume_common_grid.png)

Message:
The object of interest is the full continuum correlator manifold over the target cell. In the matched acute `4-5-6` benchmark, twisted and untwisted manifolds are highly correlated and not statistically distinguished.

Numbers to say:
`correlation = 0.9949`, `relative RMS = 0.0288`, `bidirectional RMS z = 0.75`.

## Slide 3: The Agreement Survives Pointwise Scrutiny

Title: `Matched geometry also agrees at shared coordinates`

Figure: [Matched 4-5-6 shared-coordinate FSS comparison](figures/02_acute456_shared_coordinate_fss.png)

Message:
The agreement is visible not only on a common-grid heatmap but also in continuum intercepts built from shared coordinate points. This is the concrete evidence that twisted and untwisted realizations can describe the same continuum geometry.

Say out loud:
“This is the slide that lets us say the comparison is not just qualitatively close. At the representative shared points, the continuum intercept differences sit comfortably inside present uncertainties.”

## Slide 4: Why The beta_c Finder Matters

Title: `Critical tuning has to be stable before geometry matching means anything`

Figures:

- [Acute 4-5-6 holdout beta_c scan](figures/03_acute456_holdout_betac_scan.png)
- [Partial acute 4-5-6 beta_c trend](figures/04_acute456_partial_betac_trend.png)
- [Acute 4-5-6 twisted susceptibility peak drift](figures/07_acute456_twisted_peak_drift.png)

Message:
The `beta_c` finder is doing the right physics job on the stressed twisted geometry: it sees a clear interior susceptibility peak on the held-out lattice and the completed twisted size ladder trends toward the exact isotropic triangular critical value.

Numbers to say:

- holdout `beta_c = 0.27350859 +/- 0.00099267`
- exact `ln(3)/4 = 0.27465307`
- completed twisted sizes rise from about `0.2570` to `0.2700`, with the larger held-out lattice even closer to exact

Important caveat:
This rerun was stopped before the untwisted branch and final blind-fit plots finished, so use these slides as evidence that the beta finder works, not as the final acute-456 beta-finder science result.

## Slide 5: Blind Holdout On The Acute 4-5-6 Benchmark

Title: `Twisted and untwisted training families both predict an unseen twisted datum`

Figure: [Acute 4-5-6 literal-size blind holdout](figures/05_acute456_blind_holdout_literal_sizes.png)

Message:
This is the best current direct evidence that correlator comparison is a viable couplings-to-geometry strategy. The blind point is not fit during training. Both branches extrapolate to it successfully.

Numbers to say:

- twisted best model: `power_free`, `z = -0.43`
- untwisted best model: `pade21`, `z = -0.33`
- the bad model is visibly bad: untwisted `pade11` misses at `z = +3.24`

Takeaway line:
“When we ask both constructions to predict a large unseen twisted datum, they do not just look qualitatively close; they land within one sigma.”

## Slide 6: Blind Holdout Still Works Far Outside The Training Window

Title: `The idea survives a 1000x1000 out-of-sample point`

Figure: [1000x1000 blind holdout](figures/06_iso111_L1000_blind_holdout.png)

Message:
The same logic survives on a much larger blind holdout. Both branches again favor the same asymptotic model family and predict the held-out point within acceptable Z.

Numbers to say:

- untwisted `power_free`: `z = -0.43`
- twisted `power_free`: `z = +0.90`
- Taylor and Padé alternatives miss badly, so the holdout really does discriminate models

## Slide 7: What We Can Claim Right Now

Title: `Current claim, honest scope, next step`

Message:

- Matched twisted and untwisted lattices can already agree within continuum-level uncertainties on a nontrivial `4-5-6` benchmark.
- Blind holdouts show that both branches carry predictive information about unseen large-lattice data.
- The `beta_c` finder is operational and physically well-behaved on the stressed twisted branch.
- The remaining task is to finish the full acute `4-5-6` beta-finder rerun, including the untwisted branch, and then rerun the blind comparison with per-lattice pseudo-critical tuning on both sides.

Good closing sentence:
“The project is now past the stage of asking whether there is any signal at all. The real question is how strong and how geometry-specific that signal becomes once the critical tuning is fully controlled.”

## Recommended Opening Paragraph

This project is trying to make continuum geometry measurable from lattice correlators. We simulate the same torus in different microscopic ways, twisted isotropic lattices, untwisted anisotropic lattices, and analytic modular data, and ask whether their continuum-extrapolated two-point manifolds collapse onto the same object. The evidence so far says yes in matched controls, and the blind-holdout tests say the comparison has real predictive power. That is why we think correlator matching is a viable route from couplings to geometry rather than just a diagnostic curiosity.

## Build Note

If any source figure is regenerated later, rerun:

```bash
source /workspaces/newQFE/.venv/bin/activate
python /workspaces/newQFE/responsible_method_tests/story/generate_story_assets.py
```