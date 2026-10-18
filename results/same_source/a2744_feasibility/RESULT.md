# Same-source pairing feasibility on Abell 2744: result

Plan `PLAN.md` (committed before any score); code `scripts/same_source_pairs.py` with its pass rule tested first
(`tests/test_batch_rules.py::test_same_source_verdict`). Every number is `derived` from the files named; identities
(Bergamini+23 images) were used only for scoring.

## Verdict (frozen rule): **fails** on both halves (`summary.json`)
| | development | holdout |
|---|---|---|
| queries (in-sample images with a counterpart at 2–60″) | 18 | 21 |
| P1 recall@10 (≥ 0.5) | 0.83 ✓ | 0.86 ✓ |
| recall@1 / @5 | 0.56 / 0.78 | 0.71 / 0.81 |
| median candidates per query | 585.5 | 613 |
| median false `match` per query (≤ 10) | 0 ✓ | 0 ✓ |
| positive `match` rate (≥ 0.6) | 0.17 (6 pairs with a decoy) ✗ | 0.32 (19) ✗ |
| decoy `match` rate (≤ 0.2) | 0.00 ✓ | 0.00 ✓ |

- P2 fails because true image pairs mostly fail the absolute SED-match rule: of 33 positive pairs, 20 fail the χ²
  probability and 13 the photo-z overlap; 10 pass (`positive_fail_reasons.json`, post hoc diagnosis).
- Sample: of 149 published images, 65 have a catalogue row within 0.3″ and 53 are eligible and unshared; 33 positive
  pairs at 2–60″ (1 closer). The test therefore covers bright, catalogued, S/N-selected images only.
- Morphology add-on (development, plan rule): affine fit of either handedness scores 6 positives against 6 decoys
  with AUC 0.53 (`summary_morph_dev.json`), below 0.75, so it was not run on the holdout and is recorded as not
  discriminating here. Cost 14.5 s.

## What the numbers mean (post hoc reading, labelled)
- **Ranking works, an absolute rule does not.** Relative SED similarity puts the true counterpart first among ~600
  candidates in 10 of 18 / 15 of 21 queries, but a rule that decides "same source" on its own rejects most true
  pairs, because real image pairs differ by more than achromatic magnification plus a 3 % floor (blending, ICL,
  differential magnification across the aperture, photo-z scatter).
- **Field-wide burden** (`field_burden_posthoc.json`, computed after the verdict): among 9,385 eligible sources,
  2,674,270 pairs at 2–60″; 26,145 pass the SED rule, of which 11 have both members within 0.3″ of a published image
  (10 same-system pairs pass, `positive_fail_reasons.json`); 6,439 sources sit in at least one passing pair. A field
  search with this rule would propose ~2,600 other pairs per passing true pair. That these are mostly same-redshift
  neighbours (as the orphan-pair nulls D-048, D-055 showed at 0.3–3″) is suggested, not tested here.
- The decoy match rate of 0 (25 decoys matched in brightness, colour and size, not in redshift) suggests that those
  properties alone rarely produce an SED match; the cause of the field-wide burden is not tested.

## Deviations (labelled)
- Ranking adds a χ²/ν tie-break after `match` and p_sed (not in the plan; affects only exact ties of p_sed).
- The positive and decoy match rates in P2 are over pairs that have a decoy (6 dev, 19 holdout). Over all 33 positive
  pairs the rate is 10/33 = 0.30, which also fails.
- Morphology used `os.cpu_count()` processes (4 here, as planned).
- Post hoc numbers come from the stored snippets `posthoc_script.py.txt`.
- Tool-evaluation details (install time, import results) are session observations, not stored.

## Scope
- One cluster, bright catalogued images, one catalogue, one band for morphology. Not a false-pair rate for any other
  field or the sky. No new pair was proposed or searched (branch B of the owner instruction).
