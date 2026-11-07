# Image-level two-epoch centroids of SMACS `f200w_1816` / `f200w_384` (2736 o001 → VENUS 6882 o057)

Frozen in the claim commit, before any centroid is measured. Follows `../epochs_1816_384/RESULT.md` "Next":
the catalog route was undetermined (no F444W CAT rows; D-017 error model ~10× too small for faint compact
sources, D-017 addendum).

## Hypotheses
- H-star: faint Galactic stars, ≳ 1–10 mas/yr (≳ 4–40 mas over 4.0 yr; model_prediction, order of magnitude).
- H-gal: compact background galaxies (no motion).

## Ordinary explanations
Frame-tie error, centroid noise, resampling/drizzle differences between programs, PSF differences.

## Method
- Level-3 `_i2d.fits` read from `s3://stpubdata` by byte range (`jwst_anomaly.cutouts.make_cutouts`), F150W and
  F444W, both epochs; 1.5″ cutouts at the epoch-1 F150W CAT position of each target and of every peer.
- Centroid: `photutils.centroids.centroid_2dg` on a 7 × 7 px box around the brightest pixel within 0.15″ of
  the cutout centre, after subtracting the cutout's sigma-clipped median; pixel → sky with each cutout's WCS.
- Peers (empirical null, chosen before measuring): epoch-1 F150W CAT sources within 60″ of the target,
  aper50 AB within ±1.0 mag of the target's, `semimajor_sigma` within ±50 % of the target's, isolated
  (no other CAT source within 1″), target rows excluded, mutual F150W matches within 0.3″ across epochs.
- Frame tie per band and target: median peer offset (epoch 2 − epoch 1). Residual = offset − tie;
  each peer's residual uses the leave-one-out median of the other peers.
- Validation first (synthetic): shift each peer's epoch-2 cutout by a known (+20, 0) mas and (0, +20) mas
  (Fourier shift) and re-measure; the method passes if the median recovered shift is within 5 mas of the
  injected one in both bands. If it fails, no target conclusion is drawn.

## Frozen outcome rules (ASSUMPTIONs)
- Needs ≥ 10 peers with valid centroids in a band; otherwise that band is "not measured".
- "moves: star": in **both** bands the target residual exceeds the 95th percentile of peer residuals and the
  two directions agree within 45°.
- "no motion detected": in both bands the target residual is ≤ the peer median residual; the bound is the
  peer 95th percentile (mas, and mas/yr over 4.0 yr).
- Anything else: "undetermined".

## What would change the plan
"moves" closes the record as a star. "no motion" closes the star route at the stated bound (record stays
"compact, unresolved; needs spectroscopy"). "undetermined" retires the record with that label; no further
astrometric route is planned for these two sources.
