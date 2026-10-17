# MACS0417 h002 time-delay re-check: plan (before any new number)

Owner instruction 2026-10-17 (after #308). Budget 1 h. Only the delays of #308 item 3; tools and the other items
unchanged. A reviewer's claim (fixed-β0 Fermat potential) is checked, not assumed.

1. Units and conventions: compare ∇ψ (finite differences of the glafic ψ map) with the deflection maps (sign, axes,
   units); confirm D_LS/D_S scaling applies to α and ψ alike; z_l 0.443 cited from the arXiv:1810.13439 abstract.
2. Fixed source position: β0 = θ0 − α(θ0) at h002; φ(θ; β0) = ½|θ − β0|² − f ψ(θ). Image positions refined from the
   grid candidates by root-finding β(θ) = β0 with bilinear-interpolated α (scipy.optimize.root); check convergence
   (residual < 1e-4″, and stability against interpolation order 1 vs 3 and the map grid) and that h002's delay
   relative to itself is 0.
3. Old vs new delays stored side by side for the best fit and the same 10 realisations; dates and the "not testable
   with current epochs" statement are updated only if the new delays change them.
4. The realisation range is a sample range of the first 10 realisations, not a confidence interval.
5. If not settled within the budget: delays and dates are marked "to be re-checked"; the real-change record is not
   affected.
