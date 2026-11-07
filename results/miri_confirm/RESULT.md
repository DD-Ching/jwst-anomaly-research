# MIRI single-band detections vs DJA NIRCam counterparts: result (2026-11-07)

Plan: `PLAN.md` (rule frozen before measurement). Script: `scripts/miri_confirm.py`; numbers: `summary.json`.
Contact sheet: `contact_sheet.png` (`scripts/miri_confirm_sheet.py <run_dir> <out_dir> summary.json`, visual check only).

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

## Per-band gate (2026-11-07, `PLAN_bands.md`)

Re-run `20261107T011004Z-0662e572` (same catalogs, checksums as the manifest; classify, cutouts and cross-match
off): 908 sources, 113 pass, and every group's matched source list equals the first run's. Numbers:
`summary_bands.json` (`scripts/miri_confirm.py`, now with per-band chance rates and `band_verdict`).

| Detection band (single-band-only) | n | matched | real rate | chance rate (± std over 10 shifts) | frozen verdict |
|---|---|---|---|---|---|
| F770W | 66 | 45 | 0.682 | 0.046 ± 0.019 | adopt |
| F1000W | 20 | 10 | 0.500 | 0.020 ± 0.033 | adopt |
| F1500W | 42 | 3 | 0.071 | 0.069 ± 0.027 | not adopted |
| F1800W | 35 | 3 | 0.086 | 0.057 ± 0.034 | not adopted |

**Reading (derived).** The plan's hypothesis holds: F770W and F1000W single-band detections are confirmed at
10–25× the chance rate (55 confirmed, about 3.4 expected by chance), while F1500W and F1800W matches are at the
chance level (6 matched vs about 4.9 expected). A pipeline confirmation join for MIRI is therefore restricted to
F770W and F1000W detections. F1000W rests on 20 sources; its chance rate is 0.4 expected matches. Not yet wired:
the existing confirmation path (`confirm_column`) comes from the matched-aperture join, which would also replace
the MIRI colours, so wiring needs a position-only join plus a `confirm_bands` restriction in the quality gate.

## Re-rank with band-restricted confirmation in the pipeline (2026-11-07, `PLAN_rerank.md`)

**Runs.** `20261107T021815Z-a98aae0e` (`configs/smacs_jwst3_miri_confirm.yaml`, position confirmation on) vs
`20261107T022631Z-036f249e` (the same config without the block); scratch outputs, not stored. Comparison:
`scripts/miri_rerank_compare.py` → `summary_rerank.json`. Contact sheet of the new top 20: `contact_sheet_rerank.png`.

**Checks (derived).**
1. Gate count 168 of 908 (baseline 113): the frozen expectation is met, and the 55 admitted sources equal the
   F770W/F1000W matched list of `summary.json`. Position join: 264 / 908 sources have a NIRCam detection within 0.3″.
2. Top-20 Jaccard 0.60 (15 shared). New: `f770w_221` (#2), `_239` (#4), `_176` (#13), `_118` (#16), `_240` (#18),
   all single-band F770W; dropped: `f1000w_192`, `f770w_117`, `_163`, `_177`, `_46`.
3. Visual (contact sheet viewed): #2 is a knot on the diffuse elongated light of a bright neighbour; #4 a compact
   source in the wing of a large bright galaxy; #18 a faint knot between bright sources; #16 an elongated faint
   feature next to a compact source; #13 a faint compact source. Four of five lead with CI_50_30/CI_70_50 at
   +3–5σ (high concentration from the neighbour's gradient, hypothesis), with 3–5 of the features missing.

**Reading.** 5 / 20 new entries are newly admitted singles, below the plan's "more than half" trigger, so no
separate single-band stratum. The new entries look like deblend fragments or sources on bright neighbours' light
(ordinary processing explanation, untested), not anomalies; none is vetted. Expected chance confirmations among
the 55: about 3.4 (plan).
