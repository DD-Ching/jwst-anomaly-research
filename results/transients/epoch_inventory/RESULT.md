# Two-epoch NIRCam inventory on 19 cluster fields: result

Plan `PLAN.md` (committed before any query). Tool `scripts/epoch_inventory.py`, rule tested in
`tests/test_epoch_inventory.py`. Outputs: `inventory.json` (all epochs and new pairs), `upcoming.json` (not yet public),
`macs1931_epochs.json` (post hoc). Every row is MAST metadata (`observed`); nothing was downloaded or read as science.

## Result
- 19 fields, 117 NIRCam level-3 epochs within 3′; 96 (proposal, observation) epochs are named in `results/transients`
  ("used"); 15 public epochs are unused; 82 pairs pass the frozen rule (a)–(d).
- **Fields with new pairs:** Abell 2744 78 (unused 1324 o001 / o006 / o007, 2561 o002 / o006, 3073 o005 / o008 / o009),
  MACS J1931 3, MACS0647 1 (1433 o010 → o020, F200W only, 107 d). MACS0717 and RXC J0018 each have one unused epoch
  but no partner under the rule. The other 14 fields have 0 new pairs.
- **MACS J1931 (the useful one, post hoc reading):** besides VENUS 6882 o038 (2025-09-01), two LW follow-up epochs
  target `SN_Eos_Image_101.1`: 9493 o001 (2025-10-08) and 12299 o001 (2026-09-28, public the same day), both centred
  ≈ 8″ from the VENUS centre (s_ra / s_dec in `macs1931_epochs.json`; not stored as a number) and sharing F277W / F356W / F444W with it (`macs1931_epochs.json`). SN Eos is lensed system 101 of
  the AstroLensPy v2.1 model we pinned and gate-K-tested (D-034 addendum ab); the authors predict five SN images
  (Allingham et al. 2026, arXiv:2602.14074; spectroscopic z 5.13, arXiv:2601.04156). The 12299 epoch is 392.6 d after
  VENUS.
- **Not yet public (Abell 370):** 7345 o101–o105, 9371 o001 and 5324 o015, releases 2026-12-29 to 2027-08-11
  (`upcoming.json`).

## Limits and deviations
- "Used" means named anywhere under `results/transients`; a named epoch need not have been differenced, and the
  3′ cone does not check footprint overlap (GLASS 1324 NIRCam parallels may be off the cluster). The Abell 2744 count
  is therefore an upper bound on useful pairs.
- **Deviation:** the query first ran twice with the tool uncommitted (lint failures stopped the commit; the
  first run failed at the JSON write and wrote nothing). The rule function was unchanged from the tested version;
  only docstrings, `strict=True` and a float cast changed. The `macs1931_epochs.json` query is post hoc.

## Next
SN Eos lens-model test (new plan first): from the pinned AstroLensPy v2.1 maps, the predicted positions and time
delays of SN Eos images 3–5 relative to 101.1 / 101.2; then forced photometry at those positions in 6882 o038, 9493
o001 and 12299 o001 (F277W / F356W / F444W, S3 byte-range cutouts), with injection at each position. A predicted image
that should be bright in an epoch but is absent, or a bright point source where none is predicted, is a lens-model
residual to vet; a non-detection where the model predicts a later or fainter image is not.
