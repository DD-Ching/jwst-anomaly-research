# D-041 model-copy rule as a function: gate G fails as frozen; not adopted as a class

Plan: `PLAN.md` (frozen before the run). Code: `scripts/lens_consistency.py` `model_copies` (unit test
`tests/test_lens_consistency.py::test_model_copies_positional_and_detection`), gate script
`scripts/model_copy_gate.py`. Numbers: `summary.json` (derived).

## Inputs
- Stored MACS0717 `macs0717-cats` forced-check table and its `unpredicted_catalogued_images`
  (`../macs0717-cats/`), CATS v4.1 `arcs.txt` (132 images, sha256 pinned in `MODELS`).
- rms 3.2117″: recomputed with `lens_consistency.py --model macs0717-cats validate` this run (not stored by
  `validate`; passed to the gate script and recorded as `rms_arcsec` in `summary.json`). Limit 1.75 × rms = 5.62″.
- Detection condition of D-041 not applied (`detection_condition_applied: false`; no per-image S/N is stored).

## Result
- 44 `absent` + `confused` rows; the positional rule flags 33 (`n_copy`), the same count as the docs' first pass.
- `confused`: 28 flagged, 28 in the hand list, per-system identical (no row only in one list).
- `absent`: 5 flagged vs 3 hand copies (5, 6, 34 all flagged). Rule-only: the two system-7 rows (both nearest to
  7.2).
- **Gate G: fails as frozen** (per-system counts differ in `absent`). As the plan says, the rule is not adopted as
  a class; the per-row reading stays.

## Post hoc (not part of the gate)
- The hand verdict for the system-7 pair (docs/fields/macs0717.md, "Flags other than `confused`") is "model
  position error and near-critical μ": CATS's 7.1/7.2 pair lies 2.3–4.6″ off and the μ there is untestable
  (D-036). So both lists blame model position error for these rows; they differ in bookkeeping (the hand list
  files the pair under μ-untestable), not in the ordinary explanation.
- Both system-7 rows map to the same image 7.2 because the rule takes the nearest image per row; a one-to-one
  assignment (each catalogued image explains at most one prediction) would be an amendment, untested here.
- No hand copy was missed (0 hand-only rows); the 3 copies of undetected images (45.3, 79.2, 80.2) are flagged
  because the detection condition was not applied, as the plan expected.

## Limitations
One field, one model; the gate tests agreement with a by-eye reading, not with truth. The detection condition
needs per-image S/N (fluxratio output) to be tested.
