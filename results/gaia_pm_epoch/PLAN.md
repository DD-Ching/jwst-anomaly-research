# Plan: does Gaia proper motion cost SMACS 0723 star classifications? (D-012, TASKS "Next (M2)")

Written before any query (2026-11-06).

- **Question.** `classify` matches JWST sources to Gaia DR3 positions at epoch J2016.0. SMACS 0723 (program 2736)
  was observed 2022-06 (~6.45 yr later). A star moving faster than ~46 mas/yr shifts > 0.3″ (`gaia_position`
  radius) and > ~78 mas/yr shifts > 0.5″ (query radius), so it could be missed and ranked as a galaxy.
- **Hypothesis H0 (ordinary expectation).** At b ≈ −18° the 0.1° field holds few Gaia stars with
  6.45 yr × μ > 0.3″; none of them lacks its star label.
- **Measurement.** One Gaia DR3 TAP cone around the field (the `field` backend's circle). For every source with
  pmra/pmdec: displacement d = 6.45 yr × μ. For those with d > 0.15″: nearest jwst 3.0.0 F200W catalogue
  source to the J2016 and to the propagated position, and its star label in run `20261106T171253Z-28af4751`.
- **Decision rule (ASSUMPTION).** If ≥ 1 star with d > 0.3″ has its JWST counterpart labelled `other`, or sits in
  the galaxy top 50, build propagation into `crossmatch` (next PR). Otherwise record the bound and close the
  TASKS item as "not needed for SMACS"; revisit for other fields when a field is at low latitude or > 10 yr from J2016.
