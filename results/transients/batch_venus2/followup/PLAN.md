# Batch-2 follow-up: plan (written before any new check is read)

Owner instruction 2026-10-17 (after #304). Budget 13:47Z → **15:47Z**. Exploratory: the sheets and photometry of
items 2 and 3 were already seen in #304. No new solver, no new sky, no re-run of undeterminable boxes, no threshold
changes. An item that lacks data is saved with what is missing, and the next item starts.

## 1. AS1063 1 h002 (batch 2; 342.1798224, −44.5339836; `../boxes/pending.json`)
- Not batch-1 PLCK G165 h002 or batch-1 AS1063 h005. Box `AS1063_1` of `../boxes.json` (3293 o001 → VENUS o059), an
  undeterminable box (6/19 control false hits): its result is never counted as a searched zero.
- Method (existing, unchanged): one `arc_diff.py` screen of that box with the batch settings (the session maps were
  deleted after #304), its sheet looked at, `hit_lightcurve.py` forced photometry, `batch_review.classify`, and the
  stored #296 post-review checks and full-rule function (`../boxes/review/final_classes.json` label).
- Outcomes: artifact / no significant change / real change / undetermined (incl. data insufficient).

## 2. MS0451 h011 persistence (73.548773, −3.020842; LW-only in 5058 visits at MJD 60699 and 60703)
- Persistence is detector-fixed; a sky source is sky-fixed. Within one visit the dithers put the sky position on
  different detector pixels.
- Check (minimal): list the level-2 `_cal` LW exposures of the 5058 visits covering the position (MAST), compute the
  position's detector pixel in each (header WCS), and measure a 0.15″-radius forced flux at the sky position in each
  exposure from small byte-range cutouts (no full-file download beyond what range reads need).
- Decision rule (fixed now):
  - **sky-fixed signal** if, in a visit where the mosaic detects it, ≥ 2 exposures on detector pixels ≥ 1″ apart each
    give ≥ 3σ positive flux at the sky position, with fluxes consistent within 3σ;
  - **consistent with detector persistence** if the flux is ≥ 3σ in only one exposure (or only exposures on the same
    detector pixel ± 0.5″) and absent (< 2σ) in dithers elsewhere;
  - otherwise, or if the exposures cannot be read: undetermined, with the missing input recorded.
  - "LW only" is not used as evidence either way.
- arXiv:2601.10097: read the abstract/text only to check whether h011's position and dates fall inside its search
  area and sensitivity; its null result is not used beyond that.

## 3. MACS0417 h002 counter-images (64.396407, −11.907040; on image 1.1 of system 1, z = 0.8718)
- Feasibility first: a public lens model (RELICS HLSP on MAST or published tables) with deflection maps for z = 0.8718
  (or scalable), the predicted positions of the transient's counter-images with an uncertainty, time delays and
  magnification ratios, and the epochs/depths of NIRCam coverage at those positions.
- If deflection maps exist: map h002 to the source plane with the model's α, find the image-plane positions that map
  back within the model's quoted image-position RMS, and run `hit_lightcurve.py` once at those positions (stamps looked
  at). Without time delays the result is reported as "no detection at epoch X", not as a contradiction.
- If any of these inputs is missing: record which, and stop. No lens model is fitted or rebuilt.
- Outcomes: counterpart seen at a predicted position (supports ordinary lensed transient), not seen (uninformative
  unless delays and depths say it should be seen), or not testable.

## Stop
- 15:30Z: no new checks; records, review and PR only. 15:47Z: stop.
