# MACS0416: the CANUCS-model extra images under the MARS maps (D-034 addendum ar)

Plan frozen in 50c7a749; `classify` / `match_observed` and tests in 0e239c56; both before the run. Numbers from
`summary.json` (per-prediction rows in `extra_classes.ecsv`). Labels: CANUCS image positions are observed; all
predicted positions and |μ| are model_prediction; the 3.0″ radius and the |μ| ≥ 10 split are ASSUMPTIONs (PLAN.md).

## Result
| Set | Extra predictions | both | CANUCS only | not testable |
|---|---|---|---|---|
| All CANUCS extra predictions | 49 | 38 | 11 | 0 |
| \|μ_CANUCS\| ≥ 10 | 5 | 0 | 5 | 0 |

- **The D115b flag is "CANUCS only":** the nearest MARS image not matched to an observed D115b image is 7.29″ away
  (|μ_MARS| 6.3). The other two D115b extra predictions: 10.18″ (CANUCS only) and 0.45″ (both).
- The five |μ_CANUCS| ≥ 10 predictions (D28 at 18.2, D115a at 30.0 and 20.4, D115b at 26.6 and 21.8) are all
  CANUCS only: the galaxy-scale critical-curve images are not reproduced by the free-form map.
- CANUCS only with no free MARS image in the family at all (`d_nearest` = inf): D8a, D20, D52a, D52b, D52c.
- MARS images matched to neither an observed image nor a CANUCS extra prediction, over the 40 families tested: one,
  in D115a (|μ_MARS| 29.2; its position was not stored; the three D115a CANUCS extras lie 17.1–18.7″ from it).
- Frame offset CANUCS→MARS (−0.218″, +0.082″), 242 pairs, as addendum ap.

## Reading (narrow)
- For the flag, this adds one ordinary explanation: the extra image is a prediction of the CANUCS Lenstool model
  that an independent model validated on 77 held-out MACS0416 images does not make. With the vetting record's
  finding (member 8515's σ makes it), the flag is best read as a model-specific prediction; the record stays
  `inconclusive` under its own frozen rule (the Lenstool re-fit has not run).
- The MARS map pixel is 1.515″ (D-034 addendum ao), coarse next to galaxy-scale critical curves, so "MARS does not
  predict it" is not a test of whether such an image exists. Neither model alone settles the D115 region.
- 38 of 49 extra predictions have a MARS counterpart within 3″; this stage reads no flux, so it says nothing
  about whether those images are seen.

## Not done / limits
- No forced photometry at the one MARS-only D115a image (position not stored); MARS-only extra images of families
  without a CANUCS extra prediction were not searched.
- Many-to-one matching is allowed (one MARS image can serve several CANUCS predictions; e.g. D25, D77).
