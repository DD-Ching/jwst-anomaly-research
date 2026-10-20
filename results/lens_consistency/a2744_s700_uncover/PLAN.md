# Abell 2744 system 700.1: z scan in the UNCOVER v2.0 maps (plan, frozen 2026-10-20 21:12Z before any scan; commit time is the record)

**Question.** Bergamini+23 reproduces 700.1a/b/c at no z in 1–4 (`../a2744_s700_zscan/`, outcome B): a and b
share one predicted image. D-036 says to revisit if UNCOVER v2.0 is scanned the same way. Does an independent
model (UNCOVER v2.0 maps, D-041 second model) reproduce the three images at some z?

**Ordinary explanations.** (1) model resolution in the northern infall region, 160″ out (neither model has a
constraint there); (2) a/b are the two ends of one image catalogued as two (sheet of #412); (3) wrong z;
(4) solver grid (0.25″, refined next to critical curves).

**Method (model_prediction; positions observed, Bergamini+23 `arcs.dat`).** `scripts/a2744_s700.py --model
uncover2`: the same per-z procedure as #412 (WLS source of the three back-traced images, `find_images`, one-to-one
minimum-separation assignment), z = 1.000–4.000, step 0.025, maps scaled by D_LS/D_S at each z.

**Tool gate (known cases; ASSUMPTION thresholds).** The same rule on every UNCOVER own-list system with a
catalogued z > 0 and ≥ 2 images, at its own z, at two match radii: 1.0″ (the #412 radius) and 1.5″ (the D-041
UNCOVER match radius, 1.8 × its 0.84″ own-list rms). The scan outcome uses the smallest of the two radii at which
≥ 80 % of the known systems are reproduced (`gated_match`). If neither radius reaches 80 %, the outcome is
"tool not validated" and no reading is made.

**Frozen outcome (function `verdict` at `gated_match`).**
- A "reproduced in UNCOVER": at least one z gives all three images distinct predicted images within the radius.
  Then 700.1's failure is specific to Bergamini+23 (model-dependent); not a candidate.
- B "not reproduced at any z": record the best rms and the failure pattern (shared image, missing image).
  Together with #412's B this makes 700.1 a two-model residual, and the next step is `/vet-candidate` on it
  (one source or not: colours, photo-z, morphology), as #412's plan states.
- Neither outcome is evidence of exotic lensing.
