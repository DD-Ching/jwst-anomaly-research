# Batch-2 follow-up: result

Owner instruction 2026-10-17 (after #304); budget 13:47Z → 15:47Z; plan `PLAN.md` (rules written before each new check).
Exploratory: items 2 and 3 were seen before. Every number is `derived` unless labelled; lens-model quantities are
`model_prediction`. Session-clock times are not stored.

## 1. AS1063 1 h002 (batch 2; 342.1798224, −44.5339836): artifact (edge); no significant change
- The pending batch-2 entry, not batch-1 PLCK G165 h002 or batch-1 AS1063 h005 (`as1063_h002/result.json`).
- Screen re-run with the batch settings reproduced the committed hit table exactly; sheet and stamps looked at.
- The hit is 12.54″ from the box centre in Dec (box half-width 12″): coded class `artifact (edge)`.
- Independently of the edge: a faint compact source present in 1840, 3293 and VENUS; search-epoch z F150W −1.68,
  F444W −0.64; all-mosaic |z| < 3 in every band. No significant change.
- The box is undeterminable; it is not a reliable zero-candidate search.

## 2. MS0451 h011 (73.548773, −3.020842): artifact (detector-fixed signal)
- 45 level-2 `_cal` LW exposures of 5058 (5 visits × 3 filters × 3 dithers) cover the position, all on NRCALONG;
  each exposure's header WCS was corrected by the median offset of 12 compact mosaic sources (≤ 0.26 px).
- **At the sky position:** in visits 1–2, S/N 9–25 only in dither 1 of each filter (detector pixel ≈ (586, 954));
  |S/N| < 1.6 in dithers 2–3, where the same sky lands ≈ 6″ away on the detector. Visits 3–5: |S/N| < 2.2.
- **At that detector pixel:** S/N 9–28 in all 18 exposures of visits 1–2, including dithers 2–3 that see different
  sky; 68–103 nJy without strong decay (`ms0451_h011/per_exposure_detector_pixel.json`).
- Pre-fixed rule outcome: consistent with a detector effect. Class: artifact. Mechanism (persistence or a temporarily
  hot pixel) not distinguished. "LW only" was not used as evidence.
- arXiv:2601.10097: its abstract (two swallowtail-lensed arcs) gives no transient-search area, dates or depth; not used.

## 3. MACS0417 h002 (64.396407, −11.907040): counter-image test feasible in principle, not testable with current epochs
- Inputs found: public RELICS glafic v3 (deflection and potential maps, plus realisations) and lenstool v2 (deflection
  only) on MAST. Assumptions: z_l 0.443, z_s 0.8718, flat ΛCDM (H0 70, Ωm 0.3).
- Both models put h002 on image 1.1 (μ ≈ 5.2) and predict counter-images at 7.3–7.6″ (0.29–0.41″ from catalogued
  1c.2), 8.9–9.3″ (near 1b.4) and 25.4–25.5″ (0.2″ from 1c.3) (`macs0417_h002/predicted_images.json`).
  The transient's exact counterpart position within the knot is not set by a galaxy-scale model.
- Delays vs h002 (glafic; 10 realisations): +649 to +817 d at 7.5–7.9″, +1,068 to +1,183 d at 9.3″, −4,505 to
  −5,259 d at 25″. Model-to-model systematics unknown (lenstool has no potential map).
- Coverage: 1208 (2022-10) and VENUS (2025-10). No new point source at the 7.3″ and 25.4″ positions in VENUS (stamps
  looked at). Because the 7.5″ images lag by ≈ 2 yr, in 2025-10 they show image 1.1's state of ≈ 2023-07 to 2023-12,
  for which no data exist: **the non-detection is uninformative.** The 25″ image leads by ≈ 12–14 yr (pre-JWST).
- Expected counter-image signal when it arrives: ≈ 200–330 nJy in F150W (|μ| 2–3), far above the 5–12 nJy errors.
- A counterpart, if found later, would support an ordinary lensed transient; a non-detection would not by itself
  support anything exotic.

## Deviations
- Item 2: the per-exposure WCS offset correction and the detector-pixel measurement were added after the first
  uncorrected pass returned ≈ 0 everywhere (a method check, not a threshold change).
- Item 3: realisation maps use a 0.4″ grid; image search there used a 3-pixel window and 1.0″ source-plane tolerance
  (best-fit maps: 7 pixels, 0.5″). The first 10 realisations in file order were used.

## Not done
- Item 3: lenstool delays; archival HST epochs of the 25″ image (≈ 2011–2013) not checked.
- Item 2: the mechanism (persistence vs hot pixel) and the preceding exposure history.
