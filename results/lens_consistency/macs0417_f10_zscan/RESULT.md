# MACS0417 family 10: redshift scan in the CANUCS v1 maps (result)

Plan: `PLAN.md` (frozen and pushed before the run, a60bdd2; rules as tested functions in
`scripts/m0417_f10_zscan.py`, `tests/test_m0417_f10_zscan.py`, pushed before the run). Numbers: `summary.json`.
Everything is `model_prediction` from the CANUCS v1 maps; 10.1 / 10.2 (Jauzac+19, shifted to the JWST frame by
`../macs0417_frame`) and CANUCS DR1 1115938 positions are `observed`.

## Frozen outcome: A′, "lands only outside the photo-z interval"
| Step | Result | Rule |
|---|---|---|
| Known-case check at z 2.15 | passes: 10.1 / 10.2 reproduced at 0.159″ / 0.285″ | each ≤ 1.0″, distinct images |
| z where 10.1, 10.2 and 1115938 all land (6 of 121 z) | z 1.975–2.100 | each ≤ 1.0″, distinct images |
| Closest approach to 1115938 | 0.464″ at z 2.025 (μ +3.01, positive parity) | — |
| Within the photo-z interval 1.72–1.88 (1115938 z160–z840) | third image 1.96″ (z 1.875) to 3.95″ (z 1.725) from 1115938 | — |

- The model's own fit of the pair is best near the same z: 10.1 / 10.2 separations 0.118″ / 0.128″ at z 2.05; the
  pair passes the 1.0″ check only for z 1.775–2.425 (27 of 121 z).
- Consistency (post hoc, computed from `summary.json` `known_case.predicted`): at z 2.15 the scan's third image lies
  0.098″ from the `../macs0417_canucsref` flagged position, so the scan reproduces that run's prediction.
- Reading (frozen A′): the flag of `../macs0417_canucsref` disappears at a family-10 z of 1.975–2.100, where the third
  image falls on 1115938, a source with 1.07× the predicted F200W flux (`../macs0417_canucsref/posthoc.json`). That z
  lies 0.095–0.22 above the upper end of 1115938's z160–z840 interval and below the map-fitted z 2.15. The flag
  is explained only at a z the photometry disfavours at the 84th-percentile level: **undetermined, not a candidate**.
  `/vet-candidate` not run.

## Limitations
- One model (CANUCS v1 maps). Point-source images; no model position error is propagated (CANUCS v1 image-plane
  rms on this field 1.12″, `../macs0417_images`), so a 1.96″ offset at the photo-z edge is under 2× that rms.
- Photo-z intervals are the catalogue's z160–z840; 10.1 / 10.2 have z_ML 1.80 / 1.92 (`../macs0417_canucsref`), so
  the family's own photometry spans the landing z. No spectroscopic z for family 10.
- No new contact sheet: the 1115938 neighbourhood was looked at in `../macs0417_canucsref` (post hoc there).
- 1115938 being an unrelated z ≈ 1.8 galaxy (ordinary explanation 3) is not tested (no colour or morphology test).

## Next
- A colour test of 1115938 against 10.1 / 10.2 (same-source colours, as `../a2744_s700_vet`), frozen first, would
  separate explanations 1 and 3; a z_spec for family 10 would decide it.
- The CANUCS DR1 reference on MACS1423 and the other CANUCS fields, with a local-noise `undetectable` rule frozen first.
