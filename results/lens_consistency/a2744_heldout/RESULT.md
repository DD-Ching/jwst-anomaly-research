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

- The misses depend on redshift at the same place. Families 72, 80 and 81 (z 3.48–3.75) lie within about 10″ of
  67, 68, 70 and 74 (z 2.32–2.58), which Bergamini reproduces to ≤ 0.6″. UNCOVER v2.0, which was fitted to these
  images, gives family 80 0.33–1.20″ under the WLS source (`../abell2744_uncoverref/imageplane_wls_abell2744_abell2744-uncover2-maps.ecsv`).
  Ordinary explanations, not separated here: the NW-clump mass slope that Bergamini's lower-z constraints leave
  free (the distance-ratio lever arm), the UNCOVER z of 72 / 80 / 81, and the mean-source estimator on a six-knot
  source.

## Reading (narrow)
- Bergamini+23 does not pass gate H, so Abell 2744 has no held-out-validated second model; the early two-model
  comparison (`../abell2744-uncover2/crossmodel.json`) stays unvalidated. Held-out validation now passes on one
  field only (MARS, MACS0416).
- No flag: the plan reads no flag. The z-dependent miss in the NW clump is a model-constraint question with
  ordinary candidate causes (above); a model that saw the images reproduces them.

## Not done / limits
- 41 of 43 Bergamini held-out images are type N (the family's own images anchor its source), so gate H tests the
  model's lens equation at the listed z more than a free prediction; only 34.1c and 84b are type E.
- Knot ids sit < 0.5″ apart; a held-out knot next to a seen one may be called seen (PLAN.md).
