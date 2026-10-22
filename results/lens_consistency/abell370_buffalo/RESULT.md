# Abell 370: BUFFALO Niemiec+2023 model A maps, gate K (result)

Plan: `PLAN.md` (frozen; addendum K2 frozen after K1). Script: `scripts/a370_buffalo.py`. Numbers: `summary.json`;
per-image residuals: `imageplane_abell370_*.ecsv`. All residuals are `model_prediction` minus `observed` positions.
No JWST pixel, catalogue flux or counter-image prediction was read.

## Frozen gates
- **Gate L passes** (loading): κ from the deflection-map Hessian vs the published `main` κ map, median relative
  difference 3.4e-5 (threshold 0.05, ASSUMPTION).
- **Gate K1 fails** (Lagattuta+2019 gold z_spec list, 102 images): rms 7.70″ > 1.0″; 92.2 % within 3″. As frozen,
  this stops the plan's original route.
- **Gate K2 passes** (deviation: addendum K2, the model's own released gold list, 104 images, frozen after K1):
  rms 0.910″ ≤ 1.0″; 99.0 % within 3″. It agrees with the published 0.90″ (Niemiec+2023 Table A1, model A, 98 images).
  One image is outside 3″: 28.1 at 3.36″ (`readings_abell370_k2.ecsv`, reading `lensing_residual_input`; not vetted,
  and in K1 too at 3.41″).

## Post hoc (not in the plan; not used for a decision)
- K1's rms is carried by 21.1 / 21.2 (54.1″ and 52.5″: one predicted image for the family) and by three images whose
  family has a single selected image, so the solver predicts none (10.2, 29.1, 31.1; NaN). The selection rule (class
  per image) created those singletons: a plan defect, recorded as a failed approach.
- Without the 8 images `run_field` drops (outside 3″ or NaN), K1 gives rms 0.91″ (94 images, 97.9 % within 3″).
- Lagattuta images with a BUFFALO gold image within 0.5″: 85 of 102 (`post_hoc_k1_overlap`; the split depends on
  this arbitrary radius). K1 rms is 5.87″ on those 85 (21.2 is among them) and 14.07″ on the other 17 (15 finite).
  Nearest-gold separations (`nearest_gold_arcsec`): 21.1 0.585″, 21.2 0.123″, 21.3 1.442″, so family 21's positions
  differ between the two lists. The cause of the 53″ residuals is not decided.

## Reading
- observed vs model_prediction: the BUFFALO model A maps reproduce their own constraints at the published accuracy
  through our map route (gate L and K2). Abell 370 now has an image-plane-validated map model for the extra-image
  test, validated on its own constraints only.
- Not established: that the model reproduces an independent list (K1 fails as frozen).

## Next
Pin the CANUCS DR1 Abell 370 catalogue and mosaic by sha256 and run the extra-image stages of
`../macs0416_canucsref/PLAN.md` on the BUFFALO gold systems (frozen plan first; `extra_predicted` in `summary.json`
lists 19 families with more predicted than listed images).
