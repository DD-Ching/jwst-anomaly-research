# Plan: Cosmic Gems arc (SPT-CL J0615−5746) vs the HST-era RELICS Lenstool v1 and GLAFIC v2 maps

Frozen before any map is evaluated at the arc. Answers D-034 addendum bl's "Revisit if" (a JWST-era SPT0615 image
list). Tool: `scripts/spt0615_gems.py`; rules are functions with unit tests (`tests/test_spt0615_gems.py`).

## Inputs (observed)
- Adamo et al. 2024, Nature (arXiv:2401.03224v2), Extended Data Table 3: clump positions A.1–E.1 (image 1, SE half)
  and A.2, B.2, C.2+D.2 (image 2, NW half) of the Cosmic Gems arc; E.2 (2σ) excluded as in the paper.
- Bradley et al. 2025, ApJ 991, 32 (arXiv:2404.10770v2), text: counter-image candidate (93.95000245, −57.77021846),
  photo-z 10.8 (+0.6 / −1.4); JWST-era predictions Lenstool-A (93.9490607, −57.7701814), glafic
  (93.9504865, −57.7696559). Arc photo-z 10.2 ± 0.2 (95 %), no spectroscopic redshift.
- Maps: RELICS HLSP `spt0615-relics-lenstool` (v1) and `spt0615-relics-glafic` (v2), the same products and gates
  (P, frame F) as #537 (`results/lens_consistency/spt0615_images/`). Neither map used the arc's mirror symmetry
  (HST-era constraints: Paterno-Mahler et al. 2018 families 1–3), so this is a held-out test. ASSUMPTION: the
  RELICS GLAFIC v2 constraint list is not documented in the HLSP; if it included the arc as a constraint, the test is
  a reproduction for that map (stated, not tested).

## Hypothesis and ordinary explanations
- H0 (ordinary, expected): each HST-era map puts the z ≈ 10.2 tangential critical curve through the arc between the
  mirror clumps and predicts a third image within a few arcsec of the JWST candidate.
- A miss has ordinary explanations, tested or stated in this order: photo-z (rerun at z 10.0 and 10.4); frame offset
  (the #537 frame offsets are applied); HST-era under-constraint near the arc (no HST constraint within ~20″;
  stated, not tested); single plane (the z ≈ 0.4 foreground is in neither map). A miss is a model-reliability
  result, not a lensing anomaly: the JWST-era models reproduce the arc.

## Rules (per map, per redshift z ∈ {10.0, 10.2, 10.4}; z 10.2 is the headline)
1. **Crossing** (`classify_crossing`): sample 401 points on the straight segment A.1 → A.2 (model frame). Tangential
   eigenvalue λ_t = 1 − κ − γ. Count sign changes of λ_t.
   - `through`: an odd number of sign changes on the inner sub-segment E.1 → C.2+D.2 (projected onto the
     A.1 → A.2 line);
   - `between`: no odd count on the inner sub-segment, but an odd number on A.1 → A.2;
   - `miss`: an even number (incl. 0) on A.1 → A.2.
   Report the crossing position as the offset (arcsec) along the segment from the symmetry point, defined as the
   midpoint of E.1 and C.2+D.2 projected on the segment.
2. **Parity** (`parity_split`): signed magnification at A.1 and A.2 have opposite signs. Reported; must agree with
   rule 1's odd/even count on A.1 → A.2 (a disagreement is a tool error and stops the run).
3. **Counter-image** (`classify_counter`): backtrace A.1 and A.2 to the source plane at z; find all images of each
   source point (`lensmodel.find_images`, grid step 0.25″). Counter-image = predicted images more than 5″ from the
   arc midpoint (A.1, A.2 mean). Distance d to the candidate from the nearest one, per source point; the class
   uses the smaller of the two:
   - `consistent`: d ≤ 2.2″ (ASSUMPTION: the largest JWST-era model offset from the candidate, Bradley et al.);
   - `offset`: 2.2″ < d ≤ 5″;
   - `none`: d > 5″ or no counter-image.
   Also report the distances to the Lenstool-A and glafic JWST predictions.
4. **Flag** (`joint_flag`): a lens-model residual is flagged only if every usable map gives `miss` at all three z,
   or every usable map gives `none` at all three z. Gates as #537: a map with failed gate P or frame F is unusable.

## What would change the plan
- `through` / `consistent` for both maps: the HST-era models predicted the JWST configuration; SPT0615 maps stay
  usable near the arc. No follow-up.
- A flag: record as an HST-era model limitation (ordinary), list it in the held-out synthesis; no exotic screen.
- Rule 2 / rule 1 disagreement: fix the tool, re-run, record the deviation.

## Provenance
Clump and candidate positions observed (transcribed); κ, γ, μ, crossings, predicted images model_prediction;
classes derived; thresholds ASSUMPTION.
