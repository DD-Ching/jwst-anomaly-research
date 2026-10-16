# PEARLS 1176 registration check on Abell 2744: plan written before measuring (committed 2026-10-16 04:55 UTC)

Owner request 2026-10-16. One unit of ≤ 2 h. Context: #263 saw a PEARLS 1176 cutout offset of about 0.3″ on control
c63; #264 and #265 left this offset unmeasured. Nothing below has been measured when this file is committed.

## Question
Are the PEARLS 1176 o221 (t010) level-3 mosaics offset against UNCOVER 2561 o001 (t003)? If so, is the offset
one global shift per mosaic, or does it depend on position?

## Data and tool (reuse first)
- **Data:** official MAST level-3 source catalogues (`*_cat.ecsv`, `sky_centroid`), produced by the pipeline from each
  mosaic's own WCS, for both programs in every band they share: F115W, F150W, F200W, F277W, F356W, F410M and F444W.
  They are read from S3. Nothing is stored beyond the scratchpad.
- **Baseline pair:** UNCOVER 2561 o001 against ALT 3516 o002 in F356W. Those programs have no reported offset, so
  this pair sets the match noise of catalogue centroids.
- **Image check:** if the catalogue route fails, centroids are measured on `read_box` cutouts.

## Control sources (ASSUMPTIONs, fixed now)
A source is a control if it meets all of these:
- `is_extended == False` in both catalogues;
- AB magnitude 21.0–25.5 in that band (bright enough to centroid, below the saturation regime);
- no other catalogue source within 1.0″ in either catalogue;
- a nearest-neighbour match within 0.6″ (the expected offset is ~0.3″);
- outside 3″ of any mosaic edge, i.e. inside both footprints.

Stars are kept. A proper-motion outlier is removed by the 3σ clip below.

**Split.** Each control goes to calibration set A or validation set B by the parity of crc32 of its rounded UNCOVER
position. B is not used for fitting.

## Models (fit on A, robust: iterative 3 × 1.4826 MAD clip)
- **M0:** one shift (ΔRA cosδ, ΔDec) per mosaic.
- **M1:** affine in (ΔRA, ΔDec) against position (6 parameters), for position-dependent errors.

## Decision rules (fixed now)
1. **Position dependence.** It is declared only if M1 lowers the robust RMS of B's residuals by ≥ 20 % against M0
   and its slope terms are ≥ 3σ. Otherwise M0 is used.
2. **Reliable correction (pass).** In each band, the corrected B residuals must have:
   - |median| ≤ 0.010″ on each axis;
   - robust RMS ≤ 1.5 × the robust RMS of the baseline pair (2561 against 3516, F356W, same control rules).
3. **Band consistency.** The bands of 1176 agree in shift within 0.02″, or the per-band shifts are used.
4. **If rule 2 fails in a band,** that band is not corrected. The diagnostics are saved and the work does not
   expand.

## After a pass (step 3)
- The cross-epoch photometry that uses 1176 is re-measured at corrected positions: the h182 arc-control run of #263,
  in `vet/ctrl182arc`.
- The verdicts of #261 and #263 are then re-read.
- **Exploratory label:** h182, and the controls c63, c18 and c12 that were looked at in #263, are exploratory; they
  are not a new blind test.
- No threshold of #261–#265 is changed.
