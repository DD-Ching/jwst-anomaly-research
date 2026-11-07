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

## Deblend-fragment check of the new top-20 singles (2026-11-07, `PLAN_fragments.md`)

**Data.** MAST L3 `jw02736-o002_t001_miri_f770w_cat.ecsv` + `_segm.fits` (S3 `stpubdata`, CAL_VER 3.0.0, as the
run); 264 catalog rows, 240 with finite CI_50_30 and isophotal S/N ≥ 5. Script `scripts/miri_fragments.py` →
`summary_fragments.json`. Raw files were streamed to scratch and deleted.

**Frozen rule outcome (derived).**
1. Bright-adjacency (segment 4-adjacent to one ≥ 5 × brighter, ASSUMPTION): 14 / 240 sources. Mechanism test:
   overall one-sided Mann–Whitney p = 1.1e-5 (median CI_50_30 2.13 vs 1.89), but every flux tertile holds only
   4–5 bright-adjacent sources (< 10), so the "all tested tertiles" clause has nothing to test and the frozen
   function returns **not supported**.
2. Base rates: 0 / 97 NIRCam-matched baseline gate passers, 1 / 45 admitted F770W singles, 1 / 5 new top-20
   entries (`f770w_240`, next to `_239` at 5.5 × its flux). Plan-change trigger not met.
3. Verdict per the frozen rule: all five (`f770w_221`, `_239`, `_176`, `_118`, `_240`) go to `/vet-candidate`.

**Deviation.** The baseline group is the 97 NIRCam-matched gate passers of `summary.json`, not all 113 (the full
list was not stored).

**Post hoc (not in the plan; `post_hoc_annulus` in the JSON).** The F770W CI_50_30 flat-surface-brightness ceiling
is the aperture-area ratio 2.409 (EE30 / EE50 radii 1.353 / 2.100 px, catalog meta); 9 / 240 sources reach it, including
`_221`, `_239`, `_118` and `_240` (2.52–2.53). `_221` and `_239` have isophotal footprints with equivalent radii
16.9 and 18.2 px, beyond the 8.9 px inner radius of the local-background annulus; such footprints reach the
ceiling in 3 / 18 cases vs 6 / 222 for smaller ones (small numbers; hypothesis: the annulus samples the source's
own light, so the subtracted aperture fluxes and the CIs are not profile measurements for these). Contact sheet
`contact_sheet_rerank.png` viewed again: #2 (`_221`) is marked on the diffuse band of an extended galaxy beside a
brighter knot; #4 (`_239`) is the core of a large galaxy; #18 (`_240`) a faint knot between brighter sources; #16
(`_118`) an elongated faint feature; #13 (`_176`) a faint compact source (sharpness-led, CI 1.79, below the ceiling).

**Reading.** The addendum 3 "deblend fragment" hypothesis is not supported as a group description: only `_240` is
a child next to a much brighter segment. Two of the high-CI entries are themselves large extended segments. None
is vetted; none is an anomaly claim.

## Size-aware band association (2026-11-07)

**Plan** `PLAN_sizemerge.md` (frozen before measuring). **Data:** MAST L3 `jw02736-o002_t001_miri_{f770w,f1000w,
f1500w,f1800w}_cat.ecsv` (S3 `stpubdata`, jwst 3.0.0, pixel scale 0.1109″; sha256 in the JSON), streamed to scratch
and deleted. Script `scripts/miri_sizemerge.py` → `summary_sizemerge.json`; option `catalog.merge_bands(...,
size_factor=α)` (pair accepted at `sep ≤ max(0.1″, α · r_eq)`, `r_eq` from the smaller isophotal area of the pair).

**Frozen rule outcome (derived).**

| α | pairs | added over α = 0 | chance (mean of 4 shifts) | chance fraction | passes (≤ 0.10) |
|---|---|---|---|---|---|
| 0 | 238 | 0 | 0 | — | — |
| 0.25 | 269 | 31 | 1.25 | 0.040 | yes |
| 0.5 | 312 | 74 | 8.5 | 0.115 | no |
| 0.75 | 353 | 115 | 19.5 | 0.170 | no |
| 1.0 | 399 | 161 | 31.5 | 0.196 | no |

Adopted α = 0.25. Vetted check: 1 / 4 expected merges at α = 0.25 (`_239` gains F1000W); `_176` stays single as
expected. That is below the plan's 3 / 4, so **a size-scaled radius that passes the chance rule does not explain
the vetted "single band" rows**. At α = 0.5 all four merge, but its chance fraction is 0.115.

**Re-rank** (`configs/smacs_jwst3_miri_sizemerge.yaml`, run `20261107T051243Z-b065bb6d` vs a same-session re-run of
the fixed-radius config `20261107T051306Z-bfc1921b`, which reproduces #754's gate 168 and the five new entries' ranks
2/4/13/16/18 of `summary_rerank.json`; `scripts/
miri_rerank_compare.py` → `summary_sizemerge_rerank.json`): rows 908 → 877, gate 168 → 168, top-20 Jaccard 0.739.
New: `f1500w_173` (#2), `f770w_220` (#8), `f1500w_207` (#10); dropped: `f1500w_117`, `f1500w_119`, `f770w_156`.
Contact sheet viewed (not tracked; F770W cutouts): #10 has an F770W galaxy about 0.5″ from the F1500W marker (by
eye: another split row), #2 shows no F770W source at the marker; both new F1500W rows are large F1500W/F1800W
segments (F1500W 691 and 1474 px, `post_hoc_new_f1500w_rows`) without F770W/F1000W detections, ranked by `blue_dropout`. Unvetted.

**Deviation (post hoc fix after a first, invalid run).** The plan shifted every non-reference band by the same
vector, which keeps those bands aligned with each other: 93.25 "chance" pairs at α = 0 (not stored). Each band now
gets the base shift rotated by 120° per band (same 11.3″ length). The invalid run's table (seen) had chance
fractions 0.81–0.93 and no α passing (not stored). The table above uses the fixed test.

**Post hoc (diagnostic, `post_hoc_fixed_radius`).** Fixed radii 0.1/0.15/0.2/0.3/0.4/0.6″ give 238/263/277/316/344/404
real vs 0.5/1.5/4/9.5/19.5/39 chance pairs; the chance share of the pairs added over 0.1″ is 0.090 at 0.2″ (3.5 / 39) and
0.115 at 0.3″ ((9.5 − 0.5) / (316 − 238), computed from the stored counts). This is not a rule.

**Reading.** The 0.1″ radius splits some real MIRI multi-band sources (31 added pairs at 4 % chance), but the
vetted singles need radii where chance pairs reach 12 % or more. The config with α = 0.25 brings new split-row
entries (`blue_dropout`) into the top 20, so it is not adopted as the MIRI default. Segment-overlap association
(a pair shares footprint pixels across bands) is the next test, as the plan says.

## Segment-overlap band association (2026-11-07)

**Plan** `PLAN_segoverlap.md` (frozen before measuring). **Data:** MAST L3 `jw02736-o002_t001_miri_{f770w,f1000w,
f1500w,f1800w}_segm.fits` and `_cat.ecsv` (S3 `stpubdata`, jwst 3.0.0; sha256 in the JSON), streamed to scratch and
deleted. Script `scripts/miri_segoverlap.py` → `summary_segoverlap.json`. Each other-band map is resampled onto the
F770W grid (nearest neighbour through both WCSs); an F770W ↔ band-b pair is accepted at `n_ij / min(A_i, A_j) ≥ f`,
one-to-one by decreasing overlap; "extra" pairs have catalog centroids > 0.1″ apart.

**Validation (post hoc, not in the plan; derived, `validation_r0_pairs`; every pair within 0.1″, not one-to-one, no frozen criterion).** All 232 F770W pairs within 0.1″ (F1000W 145, F1500W 55, F1800W 32)
overlap ≥ 0.75; per-band 10th percentiles 0.979 / 0.951 / 0.949. The resampling and the overlap measure reproduce the
known associations by this check.

**Frozen rule outcome (derived).**

| f | accepted within 0.1″ | extra (real) | extra (chance, mean of 4 shifts) | chance fraction | passes (≤ 0.10) |
|---|---|---|---|---|---|
| 0.1 | 231 | 124 | 49.0 | 0.395 | no |
| 0.25 | 231 | 121 | 41.75 | 0.345 | no |
| 0.5 | 231 | 114 | 34.5 | 0.303 | no |
| 0.75 | 231 | 102 | 22.5 | 0.221 | no |

No f is adopted. The four vetted rows (`_221`, `_239`, `_118`, `_240`) gain ≥ 1 band at every f (`vetted_bands_gained`),
`_176` none. The segmentation view of the five (post hoc, scratch, not tracked) shows the F770W footprints of `_221`, `_118`
overlapping their F1000W/F1500W/F1800W segments, and `_239`, `_240` partly overlapping F1000W segments of bright
neighbours. So footprint overlap does reach the vetted rows, but its extra pairs are 22–40 % chance. The chance
measure here is reference-anchored (F770W pairs only), unlike the chained merge of the size-aware test, so the two
fractions are indicative, not a strict comparison.

**Post hoc corrected counting (`post_hoc_corrected_grid`, after code review).** The frozen counting breaks overlap
ties by label, lets chance pairs use labels that the real run spends on 0.1″ pairs, measures band-b areas after
truncation at the F770W grid edge, and counts labels without a catalog row as extra; all four can inflate the chance
count. With all four fixed (native areas, ties by shared pixels, catalogued labels only, real 0.1″ labels removed
from the chance run) the chance fractions are 0.280 / 0.246 / 0.224 / 0.174 at f 0.1 / 0.25 / 0.5 / 0.75: still
none ≤ 0.10, same vetted merges. The frozen outcome stands.

**Reading.** Neither a size-scaled radius (addendum 5) nor footprint overlap separates real split MIRI rows from
chance neighbours under the 0.10 ceiling (ASSUMPTION). Per the plan's stop branch, the MIRI association stays the
fixed 0.1″ and this line is closed. Stated limitation: a MIRI row's "single band" status is unreliable for extended or
blended sources; vetting checks other-band footprints and centroids within ~0.7″ (as the five-source vetting did).

## MIRI `G_Lens` object and MIRI flux ratios of ICLv2 multiple images (2026-11-07, TASKS Then 5)

Frozen plan: `PLAN_glens.md` (rules committed 07:22Z, before any flux was read). Run `20261107T071516Z-7e9ac450`
(`configs/smacs_jwst3_miri_confirm.yaml`, jwst 3.0.0; CDS SIMBAD TAP and XMatch unreachable, NED worked).
Script `scripts/miri_fluxratio.py` → `summary_fluxratio.json`.

- **`G_Lens` object (observed):** the only top-20 row with NED `G_Lens` is rank 16 `f770w_118`
  (`[NDA2023] 01149`, 0.06″), already vetted as Mahler+2022 system 3 (`docs/candidates/smacs_miri_f770w_top20_singles.md`).
  The Oct-7 object cannot be re-identified on the pre-3.0.0 run (outputs not kept); no other `G_Lens`/`Q_Lens`
  row ranks now. Then 5 is closed: nothing new to vet.
- **Flux-ratio rule outcome (derived; ICLv2 μ at the catalogued positions, D-032 rules, 0.5″ match):** 13 / 60
  images match a MIRI row.

| Band | usable | tested (systems) | rms / max \|LOO residual\| | flagged (> 0.75 mag) |
|---|---|---|---|---|
| F770W | 9 | 6 (1, 4) | 0.36 / 0.69 mag | none |
| F1000W | 8 | 8 (1, 3, 4) | 0.64 / 1.25 mag | 4.2 (−0.97), 4.3 (+1.25) |

- **4.2 / 4.3 are chromatic, so a measurement failure under the frozen rule:** their F770W residuals are −0.37 /
  +0.69 mag (not flagged). Post hoc: the F1000W isophotes shrink to 67 px (4.2, an arc; 238 px in F770W) and 5 px
  (4.3, S/N 23.5), so isophotal F770W−F1000W is −1.57 / −1.14 / −1.91 for 4.1 / 4.2 / 4.3, while the
  aperture-corrected `aper_total` colours agree (−0.42 / −0.33 / −0.36). Lensing is achromatic; a shallower band's
  isophote truncating faint or extended images is the ordinary explanation. Cutouts (F770W and F1000W, 4″, S3 byte
  ranges; scratch, not tracked) viewed: 4.3 is a faint knot in F1000W, 4.2 an arc whose outer part drops below the
  F1000W isophote; no extra or missing image near any matched position.
- **System 3 in MIRI:** 3.1, 3.3 have no F770W/F1000W isophotal row at their nearest match (`f1500w_28`,
  `f1500w_52`, F1500W-anchored rows); `f770w_118` sits 0.40″ from ICLv2 3.3, beyond `f1500w_52` (0.16″), the same
  split-row effect as D-014 addendum 6. Only 3.2/3.4 are tested (F1000W, ±0.44 mag).
- **Result: no MIRI flux-ratio anomaly** in systems 1, 3 and 4; consistent with D-032 (NIRCam). Limits: 3 systems,
  best-fit μ only, isophotal fluxes (D-032's caveat that aperture flux follows surface brightness for arcs), 47 / 60
  images without a MIRI row within 0.5″.
