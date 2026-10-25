# SN Eos 101.2 / 101.1 flux ratio vs the VENUS model: result (2026-10-25)

Plan `PLAN.md` (frozen in a3665b33 / 8450e58f before any pixel was read). Tool `scripts/sneos_ratio.py`, tests
`tests/test_sneos_ratio.py`; functions and run driver committed before the run. Numbers from `summary.json`;
difference stamps in `contact_sheet.png`. Fluxes are observed (difference images); R_model is a model_prediction.

## Result (frozen rules)
- **Decision: "agrees".** Combined R_obs = 1.382, σ_ln 0.083, from 5 of 6 band × epoch pairs; R_model 1.162.
  |ln R_obs − ln R_model| = 0.173 (not stored; from the two stored values), below max(3σ_ln = 0.249, floor 0.30).
- Registration: every later stamp within 0.004–0.008″ of 6882 o038 (rule ≤ 0.03″). Gate D: F277W 6882 − 9493 fails
  (101.1 S/N 4.2); the other five pass. Gate I: all six pass (36 / 36 trials within 10 %).
- Per pair (6882 o038 minus): F277W 12299 1.46 ± 0.12; F356W 9493 1.37 ± 0.21, 12299 1.56 ± 0.09; F444W 9493 1.93 ±
  0.41, 12299 1.16 ± 0.07. The SN was brighter in 6882 o038 in every pair (no sign flip).
- PSF star: Gaia DR3 6766145724110284544 (G 19.5), 7.4″ from the pair.

## Reading
- Under the frozen rule the observed ratio is consistent with this model's 1.16. The rule's floor (0.30 in ln R,
  an ASSUMPTION for microlensing and macro-μ error at μ ≈ 25) decides here; on statistics alone the offset is about
  2.1σ (not stored), above the model.
- The per-pair ratios scatter more than their errors (1.16 to 1.93); the scatter enters σ_ln through rule 5. Its cause
  (chromatic microlensing, host residuals, PSF differences) is not tested.
- Contact sheet viewed: both images are point sources in every difference stamp, with no other comparable residual
  in the 6″ stamps; the direct-image panels are over-stretched (display only).

## Post hoc and deviations
- **Deviation:** the first run crashed when an injection PSF ran off the stamp edge; the clipping fix
  (commit "clip injected PSF ...") was made before any measurement was printed or written.
- **Post hoc:** gate I subtracts the un-injected aperture sums, so it tests aperture linearity and PSF placement, not
  contamination at the SN positions; it is a weak gate. The two pairs per band share the 6882 o038 stamp, so they are
  not independent and σ_ln is somewhat optimistic.
- Rule 4's "isolated unsaturated star" was implemented as the nearest Gaia DR3 source with a > 3σ parallax or proper
  motion, no other Gaia source within 3″ and a finite 5 × 5 core.
- The plan's literature check (whether Allingham et al. 2026 or the SN Eos discovery paper report this ratio) was not
  done in this run.

**Review fixes (/code-review, before merge):** a stamp that runs off a mosaic now raises; the contact sheet circles use each band's own positions. A re-run reproduces `summary.json` exactly.

## Limits
One model, no macro-μ uncertainty, no microlensing model; aperture photometry only. Data: stamps read by byte range
from 9 `_i2d` mosaics on `s3://stpubdata`; nothing stored.
