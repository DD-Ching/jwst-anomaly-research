# MIRI single-band detections confirmed by a NIRCam counterpart (TASKS Then 7, D-014 "Revisit if")

Status: plan, written before measurement (2026-11-07).

- **Question.** D-014 drops every single-band MIRI detection of SMACS 0723 (program 2736), because MIRI has
  no independent confirmation; 103 of 530 pass. How many of the dropped single-band MIRI sources have a NIRCam
  counterpart in the DJA v7.4 SMACS catalog (NIRCam detection image), and would confirmed ones change the
  MIRI ranking?
- **Hypothesis.** Most single-band MIRI drops that sit inside the NIRCam footprint are real galaxies seen only
  in F770W (the deepest MIRI band); the outskirts drops are noise.
- **Ordinary explanations for a match.** Chance coincidence with the dense NIRCam catalog (estimated from
  shifted-position matches), MIRI PSF (0.27″ FWHM at F770W) blending several NIRCam sources.
- **Rule (ASSUMPTION).** Confirmed = a DJA source within 0.3″ (about one F770W FWHM), and the chance rate
  from 10 shifted copies (offsets 5–10″) below 20 % of the real match rate. If the chance rate is higher,
  the confirmation is not adopted.
- **What changes the plan.** Chance rate ≥ 20 % of the real rate → no MIRI confirmation; record as failed.
  Otherwise wire it as a position-only confirmation join for MIRI.
