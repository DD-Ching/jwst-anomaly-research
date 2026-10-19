# W3 dimming: NEXUS `o014` footprint/depth test and forced photometry on the 12 `on_sky` targets

Run 2026-10-19 17:12–17:21 UTC under the frozen `PLAN.md`. Screen rebuilt from the S3 catalogues (same 3,793 flags
and 81,136 sources as the committed NEXUS screen; not stored), Gaia mask from the ESA TAP (#376/#377 route),
`forced` unchanged with 300 controls (`make_targets.py` leaves only the 12 targets unvetoed). Numbers are in
`summary.json` (rule outcome per target, control depth) and `targets.ecsv` (forced fluxes and flags).

## Calibration (derived)
- `calibrated` true; noise scale F200W 1.918 / F444W 2.335 (the committed `nexus_calib` run: 1.911 / 2.568).
- **Frozen requirement not met:** 8 / 300 controls flagged against the plan's ≤ 1 %. The plan said "as in the earlier
  fields", which was wrong: the committed NEXUS `forced` run flagged 10 / 300 controls
  (`nexus_calib/forced_summary.json`; Abell 2744 and MACS0416 had 0 / 300). By the frozen rule this run **does not
  count** as a calibrated null. The footprint reading below does not depend on the noise calibration (it uses image
  coverage only); the depth reading does.

## Targets (derived; rule `summarise.explanation`, unit-tested before the run)
| outcome | n | uids |
|---|---|---|
| survivor | 0 | |
| `footprint_or_image_quality` | 5 | d10423 d18183 d29073 d30717 (no other F200W epoch measurable); d29638 (vanish reappears in forced photometry; `cut_low_weight`) |
| `catalogue_depth_or_deblend` | 7 | d13073 d15524 d20291 d22021 d23045 d24470 d46213 |

- The 7 `catalogue_depth_or_deblend` targets have forced F200W S/N 4.1–28.2 in the other epochs that cover them
  (1–4 epochs each); their `o014` forced S/N is 6.0–27.0. In the sheet (`targets_sheet.png`, looked at) d13073,
  d15524, d20291, d22021, d23045 and d24470 are visibly present in every covering epoch; d46213 sits in the halo of a
  bright galaxy.
- 11 of 12 targets have at least one cutout-quality flag (`cut_edge` and / or `cut_low_weight`): in the sheet they lie
  at or beyond the edge of the other epochs' F200W mosaics. `o014` covers them; most other epochs do not.
- **Depth (derived):** median control F200W 1σ forced error 0.0115 µJy in `o014` (e2) and 0.0116 in e0, against
  0.0164–0.0171 µJy in e1 and e3–e7, so `o014` is about 1.5× deeper than six of the seven other epochs.
- d29638 (by eye, not a rule outcome): faint at the centre of the e4 stamp, noise-dominated low-weight edge in e7; its
  forced S/N is 3.9 (e4) and −0.7 (e7) against 15.1 in `o014`. It fails `cut_low_weight`, so it is not a survivor.

## Reading (narrowest statement)
- None of the 12 NEXUS `on_sky` single-epoch vanish flags survives the frozen rule. 11 of 12 lie where the other
  epochs' mosaics are absent, at their edge or at low weight; the 7 that other epochs cover are detected there in
  forced photometry. This is consistent with H0 (catalogue footprint and depth), not with a vanish.
- The screen's coverage and `edge_proxy` tests did not remove these edge positions; the persistence test cannot,
  because they are real sources.

## Next
1. Finish the persistence test on MACS0416 epochs 2–3 and Abell 2744 epochs 1–2 (#377 Next 2), which tests their 20
   targets; for each `on_sky` target, run this same forced check.
2. Before the D-039 swap, check whether the injection chain counts mosaic-edge positions as covered (if it does,
   the veto-off what-if in `persistence_scope` over-counts recoveries at edges).
3. Decide the control false-flag requirement for NEXUS from the committed runs (8 and 10 / 300) before quoting a
   NEXUS forced null.
