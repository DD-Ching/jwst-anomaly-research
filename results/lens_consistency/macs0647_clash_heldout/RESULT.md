# MACS0647: CLASH HST-era maps on the JWST-first systems (2026-10-23)

Frozen plan `PLAN.md` (pushed in 031af83 before any model read). Script `scripts/macs0647_clash_heldout.py`,
rule functions tested in `tests/test_macs0647_clash_heldout.py`. Numbers from `summary.json`; per-image
residuals in `imageplane_<model>.ecsv` (provenance `model_prediction`; positions observed, Meena+2023).
No imaging was read: this is a position test of published maps against a published image list.

## Frozen result
| | CLASH LTM-Gauss | CLASH PIEMD + eNFW |
|---|---|---|
| Gate K0, pre-JWST families 1–12 (38 images, fitted z) | **passes**: rms 1.66″, 37 / 38 within 3″ | **fails**: rms 2.34″ (> 2.0″), 31 / 38 |
| Frame check F (median observed − predicted) | passes: +0.13″ RA·cos δ, +0.10″ Dec | passes: −0.03″, +0.08″ |
| Held-out families 13–23 (31 images, fitted z) | rms 1.41″, 31 / 31 within 3″ | rms 2.67″, 21 / 31 |
| Usable for the reading (K0 and F) | yes | no |

Rule 4 reading: **no image flagged** (`flags` is empty). With its family redshifts fitted, the LTM-Gauss model
reproduces every certain image of the eleven systems first identified in JWST data within 3″, with a held-out
rms (1.41″) no larger than the pre-JWST one (1.66″). Families with more predicted than catalogued images are
listed in `summary.json` (`extra_predicted`, 17 families for LTM-Gauss; stage-2 inputs, not read).

## Post hoc: family z at the photometric z (not in the plan)
The fitted z carries much of the fit: LTM-Gauss puts MACS0647-JD (6-1 / 6-2 / 6-3) at z 5.0 against a photo-z
median of 9.2–10.6, and family 16 at the grid edge (12) against 6.69 (`z_fit`, `z_phot_median`). Re-run with
each family at the median of its image photo-z (`post_hoc_zphot`; 6-2 has none and is left out):

| | LTM-Gauss | PIEMD + eNFW |
|---|---|---|
| pre-JWST (35 images) | rms 4.17″, 20 / 35 within 3″ | rms 2.71″, 27 / 35 |
| held-out (31 images) | rms 2.88″, 25 / 31 | rms 2.91″, 22 / 31 |

At the photo-z the models do worse on the pre-JWST systems than on the held-out ones, so the per-image photo-z
(e.g. 0.52 and 9.19 within family 6-3, 0.39 / 3.38 / 11.4 within family 15; observed, Meena+2023) are not a
usable redshift constraint here; the unmatched images at the photo-z (LTM-Gauss: 13.2, 18.1, 19.2, 21.2, 22.2,
22.3) are not readings.

## Reading
- The frozen test finds no lensing-residual input: the HST-era LTM-Gauss maps predict the positions of the
  JWST-first systems to the 3″ match level once each family's z is free.
- **Limit of the test:** a free z per family absorbs one degree of freedom of a two- or three-image family, so
  this test can only find images that miss by > 3″ at every z on the grid. It cannot test the redshift
  structure, which needs spectroscopic redshifts (none in the list).
- The PIEMD + eNFW map fails its own validation gate (rms 2.34″ > 2.0″, ASSUMPTION) and is not used.
- No candidate. A JWST-era model with public maps (D-034 addendum v "Revisit if") or spectroscopic redshifts
  for the new systems would make a stronger test.
