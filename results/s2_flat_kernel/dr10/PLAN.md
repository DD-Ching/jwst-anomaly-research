# S2 with a deeper galaxy column: LS DR10, dereddened z < 22 (frozen plan, 2026-10-18)

Plan written before any DR10 galaxy was read. D-070 and its DES-SN5YR addendum: the flat-kernel term is null, but the
lensing positive control is only 1.7σ (γ_L = −0.018 ± 0.011) with LS DR9 z < 21 galaxies, so the null is not yet
calibrated by a detected ordinary effect.

**Hypothesis.** A deeper foreground column (LS DR10 Tractor + DR10 photo-z, dereddened z < 22, same chain and
selection otherwise) raises the shot-noise dilution λ_L and the lensing S/N, so the ordinary weak-lensing
magnification of DES-SN5YR SNe is detected (negative γ_L: more foreground matter, brighter SN, smaller residual).

**Ordinary explanations / systematics to watch.** Photo-z failures at 21–22 mag (scatter, catastrophic outliers move
galaxies across the z_s − 0.05 host gap); depth variations inside the DES SN fields; dust in the SN hosts correlated
with foreground counts (colour term); the same SNe as the DR9 run (not independent data: exploratory re-analysis).

**Gates (ASSUMPTION: thresholds).**
- **Q (photo-z, before the fit):** galaxies in the fetched discs with a DR10 `z_spec`, split at z = 21: in 21–22,
  σ_NMAD ≤ 0.10 (1 + z) and outlier fraction (|Δz| > 0.15 (1 + z)) ≤ 0.20. If Q fails, the fit is run at z < 21.5 only
  and labelled post hoc.
- **C (positive control):** γ_L < 0 at ≥ 3σ (|γ_L| / err ≥ 3) **and** scramble p_L < 0.01.
  - C passes → the α = 2 flat-kernel limit from this run is the calibrated S2 limit (replaces 0.019 if smaller).
  - C fails → S2 is paused: LS-depth photo-z columns cannot calibrate it on DES-SN5YR; reopen only with DES Y3 Gold
    or a deeper/spec-z galaxy column.
- The flat-kernel result is reported either way, with the chain recovery ratio and λ_F, as in D-070.

Everything else (2′ disc, Δz 0.05 shells, α = 2, P(Ia) ≥ 0.5, coverage cut, 1000 scrambles, 20 chain injections,
seed 69) is unchanged from `fit_des_alpha2.json`.
