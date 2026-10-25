# Abell 2744: Bergamini+2023 and MARS vs UNCOVER v2.0-only images: both fail gate H (D-034 addendum ax)

Plan frozen in `PLAN.md` (commit f075f0e6) before any residual; script and tests committed before the run.
Numbers from `summary.json` (per-image rows in `heldout_<model>.ecsv`) unless marked post hoc (`posthoc.json`).
Labels: UNCOVER v2.0 positions and spectroscopic z are observed (as listed); residuals are model_prediction;
thresholds are ASSUMPTIONs (`../mars_hff_heldout/PLAN.md`). No image pixels were read.

## Result
172 UNCOVER v2.0 images with spectroscopic z (`../abell2744_uncoverref/uncover2_images_spec.ecsv`).

| Model | Frame offset (E, N) | Seen / held-out | Types E / N / M | Outside map | Gate H rms | Within 3″ | H |
|---|---|---|---|---|---|---|---|
| Bergamini+23 (Lenstool, JWST GLASS-era; 149 own images) | −0.033″, +0.071″ (135 pairs) | 129 / 43 | 2 / 41 / 0 | 0 | 1.78″ | 88 % (of 43) | **fails** (rms) |
| MARS (HST-era; reported only, gate K failed) | 0.000″, 0.000″ (94 pairs) | 94 / 78 | 1 / 77 / 0 | 67 | 23.97″ | 64 % (of 11) | fails |

- **Bergamini+23:** the 5 images > 3″ are all `b` images of one source, family 80 (z 3.67; 80.2b–80.6b, 3.5–4.4″).
  18 of the 43 held-out images are its six knots. No held-out image fell outside the model.
- **MARS:** 67 of 78 held-out images lie in families with an image outside the MARS map, so only 11 are tested;
  83 (z 9.76) gets one predicted image (83a / 83b at 52″ / 60″) and 86a / 86b miss by 4.9 / 5.1″. Reported only.
- No held-out image is > 3″ in both models (`over_3arcsec_in_both` is empty; MARS tests few images).

## Post hoc (not in the plan; `posthoc.json`)
| Bergamini subset | n | rms | Within 3″ |
|---|---|---|---|
| Northern group (Dec > −30.385°), families with z < 3 (67, 68, 70, 74, 78) | 14 | 0.27″ | 100 % |
| Northern group, families with z ≥ 3 (72, 80, 81) | 24 | 2.34″ | 79 % |
| Core (34, 84, 86) | 5 | 0.92″ | 100 % |
| All without family 80 | 25 | 1.14″ | 100 % |

- In the northern group the three z ≥ 3 families (72, 80, 81; z 3.48–3.75) have larger residuals (2.34″) than
  the z < 3 families (0.27″), with the nearest pair of images 1.08″ apart (`min_sep_…`). The z ≥ 3 group is
  mostly one source (18 of its 24 images are family 80's six knots) and the Dec cut was chosen after the run, so
  this is a description, not a measured z dependence. UNCOVER v2.0, which was fitted to these images, gives family
  80 0.33–1.20″ under the WLS source (`../abell2744_uncoverref/imageplane_wls_abell2744_abell2744-uncover2-maps.ecsv`).
  Ordinary explanations, not separated here: the NW-clump mass slope that Bergamini's lower-z constraints leave
  free (the distance-ratio lever arm), the UNCOVER z of 72 / 80 / 81, and the mean-source estimator on a six-knot
  source.

## MARS vs Bergamini image by image (PLAN.md reading, from the stored tables; `posthoc.json`)
Only family 86 is tested by both: 86a 1.08″ / 4.86″, 86b 0.45″ / 5.08″, 86c 0.40″ / 0.68″ (Bergamini / MARS).
The other eight MARS-tested images (61, 62, 83, 53c) have no Bergamini row (seen by Bergamini), so the comparison
has almost no overlap. Deviation (code, after the run): the comparison is written to `posthoc.json`, and the
stored `summary.json` key `canucs_images` of the MARS block was renamed `reference_images` (the list is UNCOVER
v2.0); no number changed.

## Reading (narrow)
- Bergamini+23 does not pass gate H, so Abell 2744 has no held-out-validated second model; the early two-model
  comparison (`../abell2744-uncover2/crossmodel.json`) stays unvalidated. Held-out validation now passes on one
  field only (MARS, MACS0416).
- No flag: the plan reads no flag. The larger residuals of the z ≥ 3 NW families have ordinary candidate causes
  (above), and a model that saw the images reproduces them.

## Not done / limits
- 41 of 43 Bergamini held-out images are type N (the family's own images anchor its source), so gate H tests the
  model's lens equation at the listed z more than a free prediction; only 34.1c and 84b are type E.
- Knot ids sit < 0.5″ apart; a held-out knot next to a seen one may be called seen (PLAN.md).
