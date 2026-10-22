# MACS0417 family 10: redshift scan in the CANUCS v1 maps (plan, frozen 2026-10-22 02:20Z before any scan)

**Question.** `../macs0417_canucsref` flagged a family-10 predicted image (|μ| 2.85, `candidate`, at
64.384714, −11.915405) as `absent` at the map-fitted z 2.15. Post hoc, CANUCS DR1 1115938 lies 1.49″ away with
1.07× the predicted F200W flux and photo-z 1.81 (z160–z840 1.72–1.88), matching 10.1 / 10.2. Does the CANUCS model,
at a family-10 redshift consistent with that photo-z, put the third image on 1115938?

**Hypothesis.** The flag is a source-redshift artefact: the third image moves onto 1115938 at z ≈ 1.8.
**Ordinary explanations.** (1) wrong fitted z (family 10 has no z_spec); (2) model position error (CANUCS v1 rms);
(3) 1115938 is an unrelated z ≈ 1.8 galaxy (field density at that flux).

**Method (model_prediction; image positions observed).** Same model and frame as `../macs0417_canucsref`
(`macs0417-canucs-maps`, `apply_frame_offset`). For z = 1.00–4.00 in steps of 0.025: the source is the
magnification-weighted least-squares position of the back-traced 10.1 and 10.2 (`a2744_s700.solve`), images by
`find_images` on the scan grid; recorded per z: number of images, separations of 10.1 / 10.2 to their nearest
predicted image, and the predicted image nearest to 1115938 (separation, |μ|, parity).
Known-case check first: the same procedure must reproduce 10.1 and 10.2 (each within 1.0″) at z 2.15.

**Frozen rules (ASSUMPTION thresholds).**
- *Lands at z*: 10.1 and 10.2 each within 1.0″ of a distinct predicted image, and a third distinct predicted image
  within 1.0″ of 1115938 (the D-037 pair criterion, as in the frozen 1.0″ forced-search radius).
- Outcome A "z artefact, counter-image identified": the lands-z set intersects the photo-z interval 1.72–1.88. The
  flag is explained by redshift; 1115938 is a plausible third image (model_prediction + observed; not a z_spec).
- Outcome A′ "lands only outside the photo-z interval": record the z range; flag stays explained only at a z the
  photometry disfavours: undetermined, not a candidate.
- Outcome B "never lands": record the minimum separation and its z; the flag stays a model-position question
  (CANUCS v1 rms), still not a candidate without a second model.
- If the known-case check fails, no outcome is read.

**What would change the plan.** Outcome B with a minimum separation > 3″ would make family 10's third image a
model residual worth a second-model test (RELICS Lenstool, `../macs0417_images`).
