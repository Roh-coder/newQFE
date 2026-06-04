# Modular anchored-ratio r-grid heatmaps

Target geometry: (12, 12, 0, 0)

Target tau: 0.500000000000 + 0.866025403784 i

Target couplings: r1 = 1.000000000000, r2 = 1.000000000000

r1 range: [0.600000, 1.400000]
r2 range: [0.600000, 1.400000]

Anchor choice: nu_anchor = tau / 8.000

Point fraction: 3/8
Aggregate balance mode: mean
Observable set:
- 3v/8 anchored to nu_anchor
- 3u/8 anchored to nu_anchor
- 3w/8 anchored to nu_anchor
- (3v/8)/(3u/8)
- (3w/8)/(3u/8)

Balanced positive panels divide each observable chi2 map by its own grid scale before aggregation.
Per-observable chi2 scales:
- v4_chi2: 0.00021303864
- u4_chi2: 1.079935e-05
- w4_chi2: 9.8168909e-05
- v4_over_u4_chi2: 0.00013123402
- w4_over_u4_chi2: 4.4912121e-05

Sinh-rule mapping used:

- sinh(2 beta_c k2) = cot(alpha)
- sinh(2 beta_c k1) = cot(beta)
- sinh(2 beta_c k3) = cot(gamma)
- tau_eff = (sin(alpha) / sin(beta)) * exp(i gamma)
