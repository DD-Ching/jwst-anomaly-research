# Abell 2744 system 700.1: fine z scan in Bergamini+23 (result)

Plan: `PLAN.md` (frozen before the scan). Script: `scripts/a2744_s700.py --sheet`. Numbers: `summary.json`
(`known_cases`, `outcome`, `best_assignment`, per-z `rows`). All computed positions are `model_prediction`; the
catalogued positions (Bergamini+23 `obs_arcs.cat`) are `observed`.

**Known cases (tool check).** The same rule at each system's own catalogued z, on the 40 other systems with a fixed z
and ≥ 2 images: 39 reproduced. The one failure is 34.1, the fold pair already known to merge in this model
(D-030, D-045). The assignment rule therefore separates a known merged pair from reproduced systems.

**Frozen outcome: B, "not reproduced at any z"** (z 1.000–4.000, step 0.025, 121 values; 0 reproduced).
- Best rms: 8.28″ at z 1.10 (`best_assignment`): 700.1a 14.3″, 700.1b 0.79″, 700.1c 0.55″. 700.1b and 700.1c each
  get their own predicted image, but 700.1a's own image is 14.3″ away.
- At z 1.075–1.275 (9 values, which include the fitted 1.217), every catalogued image has *a* predicted image
  within 1″ (`nearest_arcsec`), but 700.1a and 700.1b share one image with μ ≈ −3. This is the D-030 `shared_match`
  over the whole fitted-z range, not only at the best fit.
- Above z ≈ 1.3 the back-traced sources spread out (`source_spread_arcsec` 1.0″ at z 1.3, 2.6″ at z 2.3, 3.7″ at
  z 4.0), and from z 1.9 on the nearest-image separations are 2.1–5.1″. The scan does not reach z < 1.0; the source spread is
  smallest at z 1.05 (0.90″).

**Sheet (`sheet.png`, looked at).** F150W / F200W / F444W stamps, 8″ wide, on sky offsets (the UNCOVER mosaic is
rotated; a first north-up rendering misplaced the markers and was replaced). 700.1a and 700.1b sit at the two ends
of one thin, straight feature about 1.7″ long, visible in all three bands; the model's single image at z 1.10
(μ −2.8) lies on it. 700.1c is a compact source. The sheet cannot tell "two images across a critical curve" from
"one image (or one edge-on galaxy) catalogued as two"; the model's low |μ| there (≈ 3) fits the second reading
better, but that is not tested.

**Reading.** The redshift is not what makes 700.1 fail: no z in 1–4 gives two separate images for a and b. Within
this model the failure is a missing critical curve between 700.1a and 700.1b, about 160″ from the cluster core in
the northern infall region. Two ordinary explanations stay untested: model resolution there (few constraints,
no local member potential), and a/b being the two ends of one image. The system stays "model-limited". This is not a
candidate.

**Deviations (labelled).** The known-case check and the sheet were added after the first scan (post hoc), and the sheet was
re-rendered once on sky offsets; the scan
numbers are the same with and without them. PLAN.md wording was corrected after the first scan (decision reference
D-036, claim time, grid description); the rule and thresholds are unchanged.

**Limitations.** One model. The UNCOVER v2.0 maps were not re-run here (D-041: their nearest 700.1 images are 3.56″
and 12.7″ away at their own redshift). The solver grid is 0.25″, refined to 0.02″ next to critical curves (D-040).
z < 1 was not scanned.

**Next (only if this region matters):** solve 700.1 in the UNCOVER v2.0 maps with this scan (≈ 250 MB of maps while
the run lasts), or find a model with a constraint near (3.579, −30.358).
