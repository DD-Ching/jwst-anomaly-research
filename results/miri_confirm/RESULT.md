# MIRI single-band detections vs DJA NIRCam counterparts: result (2026-11-07)

Plan: `PLAN.md` (rule frozen before measurement). Script: `scripts/miri_confirm.py`; numbers: `summary.json`.
Contact sheet: `contact_sheet.png` (`scripts/miri_confirm_sheet.py`, visual check only).

**Data (observed).** Pipeline run `20261107T001054Z-9ae92a90` (not stored; scratch outputs; catalog checksums in `data/manifests/smacs_jwst3_miri.ecsv`), sample
`smacs0723_miri` (program 2736, jwst 3.0.0 MIRI catalogs, config `smacs_jwst3_simbad.yaml` stages with the MIRI
sample of `reference_sample.yaml`): 908 sources, 113 pass the quality gate. DJA v7.4 SMACS `fix_phot.fits`
(sha256 361a4eb9…, the D-013 catalog): 12,363 rows with S/N ≥ 5 in at least one NIRCam band (0.5″ aperture).

**Frozen rule (ASSUMPTIONs).** Match within 0.3″; chance rate = mean match rate at 10 RA shifts of 5–10″; adopt
NIRCam confirmation if chance < 0.2 × real. Applied to the plan's set: sources whose only gate failure is
`single_band` (163; confirmation in `quality.assess_sources` lifts only `single_band`, not `low_snr`).

**Result (derived).**
- Single-band-only: 61 / 163 matched (0.374), chance 0.051 ± 0.013 (std over shifts); 0.051 < 0.2 × 0.374,
  so the rule says **adopt**. Expected chance matches about 8 of the 61 (0.051 × 163, not stored).
- By the band the source was detected in: F770W 45 / 66, F1000W 10 / 20, F1500W 3 / 42, F1800W 3 / 35. The
  F1500W and F1800W rates (0.071, 0.086) are close to the group's chance rate (per-band chance not computed).
- Positive control (post hoc): sources that pass the gate (mostly multi-band) match at 97 / 113 (0.858), chance
  0.048.
- Post hoc, all 756 single-band sources (any other failure too): 135 matched (0.179), chance 0.049; this group is
  not the plan's and its `verdict` field in `summary.json` is not used.

**Visual check (contact sheet, 19 seeded sources).** Matched F770W (1–8): compact galaxies in F200W, several
with a visible F770W blob. Unmatched F770W (9–12): all four lie on diffraction spikes or the halo of bright stars.
Unmatched F1800W (13–16): noise or low-weight mosaic edge. Matched F1500W (17–19): one galaxy, one inside the
halo of a large galaxy, one off the NIRCam image (no data in F200W; matched in the DJA catalogue, which includes
HST bands in its detection footprint — hypothesis, not checked).

**Reading.** NIRCam counterparts confirm real F770W/F1000W single-band MIRI sources at a rate far above chance;
single-band F1500W/F1800W detections are mostly not confirmed and their few matches are near the chance level.
The unmatched F770W singles in the sheet are spike or halo detections, which the D-014 gate correctly removes.

**Not done (deviation from the plan's last line).** The confirmation is not yet wired into the pipeline as a
position-only join for MIRI: with about 8 expected chance confirmations, mostly in F1500W/F1800W, a band
restriction (e.g. F770W/F1000W only) looks better, but that is post hoc and needs its own frozen rule. Ranking
therefore is unchanged.
