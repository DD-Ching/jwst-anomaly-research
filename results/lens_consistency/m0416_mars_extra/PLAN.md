# MACS0416: are the CANUCS-model extra images also predicted by the MARS maps? (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any MARS prediction at these positions was computed.
Follows D-034 addendum ap (`../mars_hff_heldout/RESULT.md`: MARS reproduces 77 held-out MACS0416 images, rms
0.514″) and the CANUCS-reference extra-image test (`../macs0416_canucsref/`: 49 extra predictions, one flag,
D115b, vetted `inconclusive` in `docs/candidates/macs0416-canucsref-D115b.md`). Script
`scripts/m0416_mars_extra.py`; test `tests/test_m0416_mars_extra.py` (before the run).

## Inputs
- model_prediction: the 49 extra (counter-) image predictions of the CANUCS Lenstool model
  (`../macs0416_canucsref/predicted_macs0416.ecsv`, as stored; not recomputed).
- model_prediction: MARS HFF free-form maps of MACS0416 (Cha & Jee 2023, the zip pinned in
  `scripts/mars_hff_inputs.py`; 1.515″ map pixels; gate K 0.513″, gate H 0.514″).
- observed: the CANUCS multiple-image list (`../macs0416_canucsref/canucs_m0416_images.ecsv`), the anchors.

## Question, ordinary explanations, what changes the plan
- Question: does an independent model (free-form, different constraints) also predict each CANUCS extra image?
- Hypothesis (ordinary): extra images near galaxy-scale critical curves (|μ| ≳ 10, D115a / D115b) are specific to
  one model's cluster-member parametrisation (vetting record: member 8515). Then MARS does not predict them.
- If MARS does predict D115b's flagged image (rule 2 "both"), the flag's ordinary model-side reading weakens and
  the record's next step becomes a deeper photometric / deblending check; still not a candidate in this stage.
- If MARS does not ("CANUCS only"), the flag gets one more ordinary explanation (model-specific prediction);
  the vetting verdict is updated with that evidence, not closed by it (the MARS map is coarse near galaxy-scale
  critical curves).

## Frozen rules (functions with a unit test before the run)
1. Frame: the CANUCS→MARS frame offset of addendum ap (`mars_hff_heldout.frame_offset`, recomputed from the same
   inputs) is applied to the anchors and to the CANUCS extra positions.
2. Per family with an extra prediction: source = mean back-traced position of all its CANUCS images in the MARS
   model (as addendum ap's type N anchors); `lensmodel.find_images` on a 0.25″ grid at the family z.
   `classify`: per extra prediction, "both" when the nearest MARS image not matched to an observed image of the
   family is within 3.0″ (ASSUMPTION; gate H's radius), "CANUCS only" otherwise, "not testable" when an anchor is
   outside the map. A MARS image is "matched to an observed image" when it is the nearest MARS image to an
   observed image of the family within 3.0″ (one-to-one, greedy by distance).
3. Reported, not a gate: counts per class overall and for |μ_CANUCS| ≥ 10; the MARS |μ| at each "both" match;
   the class of the D115b flag; MARS extra images (unmatched to observed or CANUCS extra) per family.

## Stop rules
- No rule change after the first MARS prediction; departures labelled deviation or post hoc.
- No flux or pixel reading in this stage; the flag's vetting record gets one dated line with the outcome.
- Data: the MARS zip (2.7 MB) is deleted after the run; only `summary.json` and `extra_classes.ecsv` are stored.
