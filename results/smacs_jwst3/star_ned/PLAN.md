# SMACS 0723 jwst 3.0.0: star-stratum sources with NED galaxy matches — plan (before the run)

TASKS "Then" item 4, re-posed on the jwst 3.0.0 catalogs. The two 2026-10-07 examples (`f090w_1719`,
`f200w_2054`, run `20261007T162133Z-2511977c`) came from the jwst 2.0.1 catalogs, which MAST no longer serves.

- **Question:** which quality-passing star-stratum sources of the SIMBAD-mirror re-run (`../simbad_mirror/`,
  `configs/smacs_jwst3_simbad.yaml`; 36 SIMBAD/Gaia astrometry, 97 stellar locus, 6 Gaia position only) have a NED
  object of a galaxy-like type (G, G_Lens, QSO, ...) within 0.5″ (ASSUMPTION, the classify radius)?
- **Hypothesis:** astrometric stars have none; any NED galaxy match sits on stellar-locus stars (photometric label
  only) and is a compact galaxy (cluster member or lensed image) that the locus mislabels.
- **Ordinary explanations:** a NED galaxy offset from the star but within the radius (blend or neighbour); NED
  position catalogues built on the same JWST images (photometric-redshift catalogues label every source "G");
  a real compact galaxy inside the locus box.
- **Test:** NED TAP cone per star; per match: type, separation, reference designation, redshift; DJA `r50` and
  photo-z of the star; contact sheet of the matched stars in F200W.
- **What changes the plan:** any astrometric star with a NED galaxy match (then the Gaia/SIMBAD star rule is
  re-examined), or ≥ 10 locus stars with a NED galaxy match within 0.3″ (then the locus rule is re-examined, and the
  affected sources are re-ranked in the galaxy stratum as a paired check).
