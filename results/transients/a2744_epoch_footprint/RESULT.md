# Abell 2744 unused NIRCam epochs: footprint check (result, 2026-10-26)

Plan `PLAN.md` (pushed before the query). Tool `scripts/a2744_epoch_footprint.py`, rules tested in
`tests/test_a2744_epoch_footprint.py`, both committed before the run. Numbers: `footprints.json`. Footprints are
MAST `s_region` metadata (`observed`); areas are `derived`. No pixel was read.

## Result (frozen rules)
- **Outcome O (overlap only).** Of the 78 inventory pairs, 0 are `core`, 28 `overlap`, 50 `none`.
  Every pair had footprints in at least one shared filter; rule 5 dropped no whole pair.
- No unused epoch reaches the 1′ core disc (π arcmin²) by more than 0.315 arcmin² (3073 o005 / o008); 1324 o001 /
  o006 / o007, 2561 o002 / o006 and 1176 o221 have 0 core area in every filter.
- Largest overlaps (best filter, all outside the core):
  - 1324 o001 → o007: 6.327 arcmin² in F115W, 373.3 d;
  - 2561 o002 → o006: 5.817 arcmin² in F277W, 365.6 d;
  - 1324 o001 → o006: 4.436 arcmin² in F356W, 134.7 d;
  - 1324 o006 → o007: 4.297 arcmin² in F444W, 238.7 d.
- So the inventory's "Abell 2744 78" is, under this rule, 28 off-core pairs and no pair through the strong-lensing
  core. The core is imaged only by epochs the transient runs already used.

## Post hoc (not in the plan)
- **Positive control:** the used core epochs 2561 o001 (UNCOVER), 2883 o002, 3516 o001 / o002, 3538 o007 / o008
  and 4111 o001 each cover 3.04–3.14 arcmin² of the core disc in their best filter (`epoch_filter_area`), so the
  projection and the core test see the core where it is known to be imaged.
- 1324 (GLASS) and 2561 o002 / o006 are flanking pointings, as the plan's ordinary explanation expected.

## Reading and next
- A D-027 grid on an unused Abell 2744 pair is a flanking-field transient search: low lensing magnification, not
  a test of the cluster lens model. Its value is a general two-epoch rate, not the mission's lensing question.
- If one is run, the frozen best candidates are 1324 o001 → o007 and 2561 o002 → o006 (largest overlap, ~1 yr
  apart, several shared filters).
- The Abell 2744 core needs a new epoch; the inventory's `upcoming.json` lists none for Abell 2744.

## Limits
`s_region` is the mosaic outline: chip gaps and masked regions inside it count as covered, so overlaps are upper
bounds. One centre, one core radius (ASSUMPTION 1′). The inventory's 3′ cone still bounds which epochs were seen.
