# MACS J1931 counter-image route: proposal (not frozen)

Status: **proposal, not frozen**. Nothing here has been run; no predicted position, map value or pixel was read for it.
Follows `../macs1931_gatek/RESULT.md` (gate K passes; spec-z systems 1, 5, 6, 7, 100, 101 have more predicted than
catalogued images) and the footprint step (`RESULT.md` here).

## Why a new reference is needed
The extra-image test (`../bullet_extra/PLAN.md` route) needs an arc-aware reference catalogue and passes gate C2 per
field first. MACS J1931 has none public: no DJA field, no CANUCS field (D-034 addendum aa). Pipeline isophotal
segments were rejected as a reference in cluster cores (D-034 addendum 2026-10-22: MACS0417 gate C2 10 / 27).
That addendum names the two alternatives: a deblended arc catalogue, or an arc-aperture sum on member-subtracted
images.

## Proposed route (to be frozen in its own commit before any run)
1. **Own deblended catalogue.** On an S3 byte-range cutout of the GO 6882 F200W level-3 mosaic covering the
   catalogued images plus a margin, subtract a smooth background that carries the cluster light: photutils
   `Background2D` with a box small enough to follow the bCG / member halos (ASSUMPTION, to be fixed: of order 1–2″).
   Then detect sources with `detect_sources` and split them with `deblend_sources`; measure Kron totals and shapes
   with `SourceCatalog`. Map them into the `m0417_arcref` catalogue form (`ra`, `dec`, `a`, `b`, `pa`, `flux`,
   `snr`, `ext`, `ell`, `scale`). photutils is the same library the JWST pipeline's source catalogue step uses; the
   change from the rejected pipeline segments is the background, which keeps cluster light out of the arc segments.
   Adding photutils as an explicit dependency needs `/reuse-check` and a D-NNN entry first.
2. **Validation on a known case (stage V), before MACS J1931.** Run the same builder on the MACS0417 GO pipeline
   F200W mosaic and pass it to the unchanged frozen MACS0417 route. Gate C2 there is 37 / 40 with CANUCS DR1
   (passes) and 10 / 27 with pipeline segments (fails). The builder is accepted for MACS J1931 only if gate C2
   (shared-pair rule) passes on MACS0417 with it; a failure ends the route without a MACS J1931 reading.
3. **MACS J1931 stage**, only after V passes: gate K-WLS repeated (`../macs1931_gatek`, spec-z systems only),
   frame tables with the WLS source, local noise on the same F200W mosaic (constants of `../bullet_extra/PLAN.md`
   rule 5), gate C2 (shared pairs), flags = predicted targets classed `absent` while C2 passes. Photo-z and mixed-z
   systems (8, 10, 11, 13, 15, 16, 17) stay out until their z is settled (`../macs1931_gatek/RESULT.md` reading).

## Hypothesis, ordinary explanations, what would change the plan
- Question (science, stage 3): do the VENUS maps predict a magnified, unmatched counter-image of a spec-z system
  where the F200W mosaic shows no source at the predicted flux?
- Ordinary explanations for any flag: model μ / position error near critical curves, low surface brightness or
  deblending, cluster-member light left by the background filter, a reference flux that includes unrelated light,
  wrong system z. Each goes to `/vet-candidate` before any wording beyond "flag". Exotic readings are not considered.
- Stage V fails → try the other alternative (arc-aperture sums on `member_sub` residuals) as a new plan, or wait
  for a public VENUS catalogue (the VENUS HLSP on archive.stsci.edu could not be checked: TLS fails through the
  cloud proxy; a local session can check it).
