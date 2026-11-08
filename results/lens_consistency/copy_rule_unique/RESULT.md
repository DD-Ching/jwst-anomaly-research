# D-042 second copy rule as a function: gate U fails as frozen on MACS0416 (system 122)

Plan: `PLAN.md` (frozen before the run). Code: `scripts/lens_consistency.py` `model_copies_unique` (unit test
`tests/test_lens_consistency.py::test_model_copies_unique_one_to_one_and_parity`), gate script
`scripts/copy_rule_unique_gate.py`. Numbers: `summary.json` (derived), `post_hoc_122.json` (derived, post hoc).

## Route
- `lens_consistency.py --model macs0416-cats images --step 0.25 --forced-image <F277W i2d, S3 byte range>` on the
  F200W pipeline catalogue `jw01208-o004_t002` (sha256 `00cb6e18…`, as the manifest), CATS frame offset applied.
  Output: `images_macs0416.json`, `images_predicted_macs0416.ecsv.gz`.
- Route gate passes: 41 systems, 143 predicted; forced recovered 3, confused 13, absent 2, as recorded in
  docs/fields/macs0416.md. CANUCS photo-z were not passed (`--photoz`); the forced counts match the record without them.
- Unpredicted catalogued images: 1.2, 15.1, 16.2, 45.2, 55.2, 60.3, 122.2, 122.3. Radius 6.5 × 0.76″ = 4.94″.

## Result
- 15 `absent` + `confused` rows; the rule labels 2: system 1 (copy of 1.2, 1.58″) and system 15 (15.1, 2.57″), both
  in the hand list.
- System 45 (μ -8.99) pairs with 45.2 (1.89″), whose model μ is +2.28: opposite parity, not labelled, as D-042 reads it.
- System 122 ×2: not labelled; neither row has an unpredicted image within 4.94″.
- **Gate U: fails as frozen** (hand-only: 122 ×2). As the plan says, the rule is not adopted; the by-eye reading stays,
  and the post-U diagnostic on MACS0717, MACS1149 and Abell 2744 was not run.

## Post hoc (not part of the gate)
- Post hoc: the two 122 rows lie 6.44–7.00″ from 122.2 / 122.3 (`post_hoc_122.json`), outside the 1.6–4.7″ range the
  docs give for this verdict; 4.60–4.71″ is their distance to 122.1, which is matched to an observed row. So the
  recorded copy reading of 122 is not supported by the re-run separations. Both rows are `confused`
  (a brighter source in the search circle), not `absent`, so no unmet prediction appears: the row class is
  unchanged and no flag or candidate results. The system-122 verdict in docs/fields/macs0416.md is now marked as
  not reproduced.

## Limitations
One field, the rule's own defining case; the gate tests agreement with a by-eye reading, not with truth. Parity
uses the CATS map μ at the catalogued position at the system's redshift, a model prediction.
