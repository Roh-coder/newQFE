# Modular anchored-ratio r-grid heatmaps

Target geometry: (12, 12, 0, 0)

Target tau: 0.500000000000 + 0.866025403784 i

Target couplings: r1 = 1.000000000000, r2 = 1.000000000000

r1 range: [0.600000, 1.400000]
r2 range: [0.600000, 1.400000]

Anchor choice: nu_anchor = tau / 8.000

Point fraction: 3/8
Aggregate balance mode: none
Balance rule: raw observable chi2 maps are summed without rescaling
Observable set:
- 3v/8 anchored to nu_anchor
- 3u/8 anchored to nu_anchor
- 3w/8 anchored to nu_anchor
- (3v/8)/(3u/8)
- (3w/8)/(3u/8)

Additional derived metrics:
- v+w score = chi2[3v/8 anchor-ratio] + chi2[3w/8 anchor-ratio]
- v-w score = (log residual[3v/8 anchor-ratio] - log residual[3w/8 anchor-ratio])^2

Sinh-rule mapping used:

- sinh(2 beta_c k2) = cot(alpha)
- sinh(2 beta_c k1) = cot(beta)
- sinh(2 beta_c k3) = cot(gamma)
- tau_eff = (sin(alpha) / sin(beta)) * exp(i gamma)
