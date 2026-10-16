# Per-pair vetting of the 11 escalated Abell 2744 orphans: plan written before measuring (2026-10-16)

Owner request 2026-10-16: finish the per-pair check that #269's frozen escalation rule calls for. One unit of ≤ 2 h.
- **Fixed:** no threshold of #269 changes, there is no new null, no parameter scan and no large download.
- **Pairs** (`../top_orphans.ecsv`): #1, #2, #3, #5, #6, #7, #9, #10, #11, #13 and #14. The main effort goes to #2, #3, #9,
  #11 and #13. Where #269 already holds sufficient evidence, it is reused.

## Data (existing, small)
- **Images:** UNCOVER 2561 o001 t003 level-3 mosaics in F115W, F150W, F200W, F277W, F356W, F410M and F444W. They
  are read as 6″ boxes by S3 byte range (`arc_diff.read_box`).
- **Catalogue:** the #269 matched-pair table (`matched_pairs.ecsv`, scratchpad), with the DJA photo-z.
- **Public spectroscopy:** VizieR cone searches of 1.5″ around each member in published Abell 2744 redshift
  catalogues (MUSE, UNCOVER). There is no bulk download.
- **Not used:** DJA segmentation maps (a large download). Blending is judged from the catalogue `deblend` flag and
  the stamps instead.

## Measurements per pair (all derived; every threshold an ASSUMPTION fixed now)
- **E1 bridge (one system).** Take a 0.1″-wide strip along the segment joining the members, excluding 0.15″ around
  each member. Measure its mean surface brightness in the summed LW (F277W+F356W+F444W) image, against the robust σ of
  the same strip shape placed at 40 random positions in the box.
  - A bridge of ≥ 3σ counts as connected emission: knots of one galaxy, or an interacting pair.
- **E2 surface-brightness test (lensing).** Lensing conserves surface brightness, so two images of one source have
  equal surface brightness and equal colours.
  - Size: from second moments of the background-subtracted LW stamp within 0.3″ of each member.
  - Surface brightness: the 0.15″-aperture flux divided by the PSF-deconvolved area.
  - **Inconsistent with one lensed source:** the surface-brightness ratio differs from 1 by > 3σ.
  - **Undetermined:** the members are unresolved (deconvolved size below 0.5 × the PSF FWHM). Unresolved members
    carry no surface-brightness information.
- **E3 colours.** Per-band 0.15″ flux ratios, a/b, in 7 bands.
  - **Inconsistent:** the χ² of a constant ratio has p < 0.01.
  - This re-tests #269's SED match on our own photometry. It is not a new null.
- **E4 deflector light.** Measure the flux in a 0.2″ aperture at the brightest point of the joining segment (stacked
  LW) and give a 3σ upper limit on a deflector magnitude there.
  - Compare it with the SIS σ that #269 lists for the pair, using the D-048 scaling (σ* = 180 km/s at m*(F160W) = 19,
    hypothesis-level).
  - The limit only states how dark a deflector would have to be. On its own it is not evidence for one.
- **E5 spectroscopy.** If both members have public spectroscopic redshifts that differ by more than 0.01 (1 + z),
  the pair is not one source.
- **E6 artifacts and neighbours.**
  - Detection: each member must be detected at ≥ 3σ in at least 3 bands.
  - Artifacts: a spike, edge or stripe within 0.5″, or a pixel-level defect on the stamp.
  - Neighbours: count catalogue sources within 1.5″ of the pair midpoint with overlapping photo-z (from
    `matched_pairs` and the top-orphan table). With ≥ 2 such third sources, the pair is part of a compact same-z group.

## Classes (decision rules)
- **Ordinary explanation supported.** At least one of:
  - an E1 bridge;
  - an E2 surface-brightness mismatch;
  - an E3 colour mismatch;
  - an E5 redshift mismatch;
  - an E6 artifact;
  - an E6 compact group of ≥ 2 third sources **together with** E1 or E2 evidence.
  A group on its own is not enough.
- **Insufficient data.** None of the above, and E2 is undetermined (unresolved members) or a member fails the E6
  detection rule. The record names the minimal new data, e.g. a spectrum of both members, or deeper or
  higher-resolution imaging.
- **Anomaly worth independent checking.** All of these together:
  - no ordinary evidence;
  - E2 equal surface brightness and E3 equal colours, both measured (not undetermined);
  - an E4 deflector limit fainter than the SIS-scaled magnitude by ≥ 3 mag.
  Even then it is only a candidate for `/vet-candidate` and independent data.
