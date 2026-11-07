# SMACS f200w_3670 is a foreground star on the end of lensed image 2.3; no lens anomaly

Plan: `PLAN.md` (pushed before any image was read). Data: jwst 3.0.0 program 2736 o001 pipeline catalogs
(`*_cat.ecsv`, streamed, not stored) and S3 byte-range cutouts of the i2d mosaics (`scripts/smacs_f3670_sheet.py`,
sheet `sheet.png`). System 2 positions: Mahler+2023 ICLv2 `arcs.dat` (sha256 pinned in `lensmodel.py`). Numbers in
`summary.json`.

## Result (observed / derived)
- **Separate object.** In every band the pipeline gives f200w_3670 its own segment; the 2.3 arc is the neighbouring
  segment (F200W label 3671) 0.81″ away. The NED 0.369″ separation is to the Mahler+2023 NED position of 2.3, not to the arc
  light.
- **Colours** (aper_total, AB): f200w_3670 F090W−F444W −0.95, F150W−F277W −0.97, brightest in F150W (22.23); arc 2.3
  +2.28 / +0.91; counter-image 2.1 +2.41 / +1.04. The source differs from 2.3 by 3.2 / 1.9 mag (rule threshold 0.3,
  ASSUMPTION); 2.1 and 2.3 agree within 0.13 mag, as lensing requires.
- **Looked at** (`sheet.png`): f200w_3670 shows the six-spike NIRCam diffraction pattern in F150W and F200W; the
  2.3 arc is a separate elongated object extending from it (mosaic pixel frame, orientation not checked), and 2.1 (compact) and 2.2 (long tangential arc) look
  like normal images of the same source.
- **Outcome A** (plan rule 2): a Galactic star projected onto one end of arc image 2.3, not a lensed knot.
  The stellar-locus label is correct. No Gaia or SIMBAD object within 0.5″ (`../star_ned/followup.json`); at F200W
  22.5 it is likely below Gaia's limit (not checked).

## Consequence
- No counterpart check is needed (rule 5 applies only to outcome B). TASKS Then 4 is closed.
- The star's light can bias 2.3 photometry near the star; the pipeline arc segment's colour still matches 2.1, so no
  system-2 colour inconsistency follows from it.

## Limitations
- Colours from pipeline `aper_total` apertures; the star and the arc tip partly overlap, so either can contaminate
  the other at the 0.1-mag level (not quantified).
- Star class from morphology and colour only; no spectrum or proper motion.
- **Deviation from PLAN rule 1:** the cutouts are 12″ and 4″ in F150W / F200W / F444W instead of 6″ in F200W / F444W,
  and there is no colour composite; colours come from the catalog instead. The rule-2 inputs (the spike pattern in
  F200W and the two colours) are unchanged.
