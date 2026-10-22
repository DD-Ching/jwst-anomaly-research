# MACS0417 family 10: catalogue colour test of CANUCS 1115938 (result)

Plan, script and tested rules: `PLAN.md`, `scripts/m0417_f10_colour.py`, `tests/test_m0417_f10_colour.py`, all
pushed before any measurement (f7541d1). Numbers: `summary.json`; post hoc specificity `posthoc.json`. Fluxes are
`observed` (CANUCS DR1 0.3″ colour fluxes, 8 JWST bands); χ² values are `derived`.

## Frozen outcome: "H1 not excluded"
| Step | Result | Rule |
|---|---|---|
| Known-case gate (first two images of other Jauzac+19 families, F200W S/N ≥ 10) | passes: 8 / 9 consistent (17.1–17.2 differs, χ² 24.7 / 7) | ≥ 8 pairs, ≥ 80 % |
| 10.1 vs 1115938 | consistent: χ² 2.80 / 7, p 0.90, flux ratio 4.69 | p < 0.01 = differ |
| 10.2 vs 1115938 | consistent: χ² 3.98 / 7, p 0.78, flux ratio 2.78 | p < 0.01 = differ |
| 10.1 vs 10.2 (reference) | consistent: χ² 5.00 / 7, p 0.66 | — |

Reading (frozen): 1115938's colours match both catalogued family-10 images, as a third image requires. Equal colours
are necessary, not sufficient, so this does not prove the identification.

## Post hoc (not frozen rules)
- **Specificity** (`posthoc.json`): of 2,033 CANUCS sources with F200W S/N ≥ 10 and 0.5–2× 1115938's F200W flux, 153
  (7.5 %) are consistent with both 10.1 and 10.2. A random galaxy of that flux passes this test with about that rate.
- **Flux ratios against the model** (from `../macs0417_f10_zscan/summary.json`, z 2.025 row, μ +13.10 / −8.38 /
  +3.01): predicted |μ| ratios 4.35 (10.1 / third) and 2.78 (10.2 / third) against the fitted colour-flux ratios 4.69
  and 2.78. Not stored as a separate number; computed from the two committed files.
- Together with the z scan (third image within 1.0″ of 1115938 at z 1.975–2.100): position, colours and flux ratios are
  consistent with 1115938 being family 10's third image at z ≈ 2.0, with its photo-z (1.72–1.88) lower by about 0.1–0.2.
  The #453 `candidate` flag has an ordinary reading: source redshift. No candidate; `/vet-candidate` not run.

## Limitations
- Deviation (after /code-review, before merge): a known pair whose two images match one CANUCS source is skipped (none
  did; `summary.json` unchanged), and the specificity sample excludes 1115938 itself.
- Catalogue photometry only; 0.3″ colour apertures on lensed images can mix components (17.1–17.2 differs).
- One model for the position and magnifications (CANUCS v1); no z_spec for family 10.
- No new contact sheet; the 1115938 neighbourhood was looked at in `../macs0417_canucsref` (post hoc there).

## Next
- MACS0417 family 10 is closed as ordinary unless a z_spec contradicts z ≈ 2.0. Next lensing step: the CANUCS DR1
  reference on MACS1423 and the other CANUCS fields, with a local-noise `undetectable` rule frozen first.
