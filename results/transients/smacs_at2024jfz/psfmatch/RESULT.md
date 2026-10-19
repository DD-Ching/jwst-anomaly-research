# AT 2024jfz z_s 1.55–1.70: PSF-matched NIRISS 4368 vs ERO F150W (result, 2026-10-19)

Plan: `PLAN.md` (frozen in `claude/at2024jfz-psfmatch` before any pixel read; rule `outcome_psfmatch` unit-tested
before the run). Every number is in `result.json`. Positions `model_prediction` (Mahler+22 ICLv2, best model);
fluxes `derived`; thresholds ASSUMPTIONs. `stamps_psfmatch.png` was looked at before this was written.

## Frozen outcome: **untestable** (all four rows; subtraction gate)
- Kernel K (9 × 9 E pixels, sum 0.9654, fit χ² 0.779 per pixel) maps T (resampled onto E's 0.0656″ grid) to E;
  σ_D = 6.48 nJy (aperture r 0.15″).
- Positive control **passes**: D at the TNS position 228.6 nJy (plain E − T in #360: 223.4 nJy).
- Injection gate passes in every row (D_inj 179.3–216.8 nJy against thresholds 96.9–99.7 nJy).
- Subtraction gate **fails**: D at the static source centroid is −131.9 nJy, |D| above every row threshold. The
  disc maxima (97.2–167.2 nJy) lie 0.203–0.258″ from the source (computed from `delta_max_radec` and `static_source.radec`; distances not stored), on its
  core residual.

## Post hoc (not in PLAN.md; does not change the outcome)
- Stamps: the static source shows diffraction spikes in T, so it is star-like (by eye; no catalogue check). D
  shows a compact positive / negative pattern within about 0.25″ of its core and is flat at the four predictions.
- Disc points more than 0.35″ from the source (69–75 of 81 per row), subtraction gate waived: max D 9.9, 31.5, 45.8
  and 41.8 nJy for z_s 1.55, 1.60, 1.65, 1.70, all below the thresholds (96.9–99.7 nJy); the same rule gives
  "absent" in every row (`posthoc.outcome_outside_core`). D at the predicted positions themselves: −16.7 to
  17.4 nJy.
- Reading (post hoc): under a single-peaked light curve and this model, a counter-image at the predicted
  z_s 1.55–1.70 positions (0.57–0.70″ from the star) is not seen; one within ≈ 0.35″ of the star remains untested.

## Limitations
- Best model only; 0.5″ position radius assumed. Rows step by 0.05 in z_s.
- K is one kernel for the whole box; a bright star's core needs more than a 9 × 9 kernel or a sub-pixel registration
  (the frozen plan allowed no tuning).
- The T variance after resampling uses an independent-pixel area scaling (fit weights only).
- No candidate can come from this step.
