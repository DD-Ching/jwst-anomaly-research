# Two-epoch check of SMACS `f356w_1243` (2736 o001 → VENUS 6882 o057, F444W)

Frozen before the comparison is computed (claim commit). Open question in
`docs/candidates/jw02736-o001_t001_nircam_f356w_1243.md`: "Is there a second epoch of NIRCam imaging that
could test the transient hypothesis?"

**Disclosed before freezing:** the coverage check read the nearest-row F444W `aper50_abmag` in both CATs
(26.114 in 2736, 26.070 in 6882; separations 0.006″ / 0.009″). No frame tie, peer scatter or centroid
residual was computed before this file was written.

## Coverage (observed)
- 2736 o001 (2022-06): F150W, F356W, F444W CATs; the target has F356W and F444W rows only.
- 6882 o057 (2026-06-05): F150W, F182M, F210M, F300M, F410M, F444W; **F444W is the only shared band with a
  target row**. 12774 o002 does not cover the target (nearest row 18.7″).

## Hypotheses
- H-gal: dusty compact galaxy at z_phot ≈ 5.6–5.75 (record verdict): constant flux, no motion.
- H-AGN: red compact AGN: flux change possible (rest-frame ~0.6 yr baseline at z 5.7), no motion.
- H-SN: a transient in 2022 (2023ApJ...947L...1Y's question for point-like candidates): faded by 2026
  for ordinary SNe (4.0 yr observed ≈ 0.6 yr rest frame); long-lived superluminous / pair-instability
  events need not have faded (model_prediction, order of magnitude).
- H-BD: a cool brown dwarf (known red point-source contaminant): proper motion ~ v_t / d, e.g. 100 km/s at
  2 kpc ≈ 10 mas/yr ≈ 40 mas over 4 yr (model_prediction, order of magnitude).

## Ordinary explanations of a measured change
Frame-tie error, centroid noise at faint S/N, deblending differences between pipeline versions, zero-point
differences (cancel in the frame-tie median, D-017), and the PSF wings of `f200w_1874` 1.7″ away.

## Method (D-017 tool, with the peer null of D-017 addendum 1)
`scripts/epoch_compare.py` `compare_epochs` on the two F444W CATs (defaults: references S/N ≥ 30 within 60″,
match 0.3″). Peers: mutual matches within 60″ of the target, |aper50 mag − target| ≤ 0.5, size within 30 %,
target rows dropped (`smacs_epochs_1816_384.empirical_null` selection). Peer Δm and positional residuals use the
same frame tie.

## Frozen outcome rules (ASSUMPTIONs; F444W only; computed by `decide` in the script)
- Target rows must be mutual nearest neighbours; else "not measured". Fewer than 10 peers → "undetermined".
- Variability, with σ_peer = 1.4826 × MAD of peer (Δm − frame median): |Δm_target − frame median| ≥ 5 σ_peer
  → "variable"; ≤ 2 σ_peer → "constant" (H-SN for ordinary SNe disfavoured; fade bound 3 σ_peer stated);
  else "undetermined".
- Motion: target residual > peer p99 → "moves" (H-BD supported); ≤ peer p95 → "no motion detected" (bound
  p95 / 4.0 yr); else "undetermined".

## What would change the plan
"moves" → brown dwarf; record closes as a star. "variable" → `/vet-candidate` (AGN or transient). "constant" +
"no motion detected" → H-gal / H-AGN stay, ordinary SNe disfavoured; the record stays "known population".
