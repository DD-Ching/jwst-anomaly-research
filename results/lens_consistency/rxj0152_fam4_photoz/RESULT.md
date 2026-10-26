# RX J0152 family 4: RELICS photo-z of the three chance-level matches (result, 2026-10-26)

Plan `PLAN.md` (pushed in 11bd41f9 before any catalogue row was read; tool and rule tests pushed before the run).
Numbers: `summary.json` (frozen run) and `posthoc.json` (post hoc). Catalogue: RELICS
`hlsp_relics_hst_acs-wfc3ir_cl0152m13_multi_v1_cat.txt` (13,052 rows, sha256 in `summary.json`; read in the session,
not stored). Catalogue positions `observed`, BPZ photo-z `derived` (RELICS), family z `model_prediction` (SLICE;
family 1 spectroscopic), classes `derived`. No imaging was read in this run (bg's contact sheet covers the region).

## Frozen result
- **Gate P fails:** of the 18 SLICE images of families 1–6, 14 match a RELICS galaxy row within 0.6″ and only 5 / 14
  have their family z inside BPZ's 95 % interval (≥ 60 % needed). **Outcome U: no reading.**
- The spectroscopic family fails too: 1.1 (z_spec 3.93) matches z_b 3.785, interval 3.697–3.868, odds 1.0; 1.3's
  interval 0.273–3.693 (odds 0.094) also excludes 3.93; 1.2 and 1.4 have no match.
- Family 4's own images disagree with each other and with SLICE z 1.69: 4.1 z_b 3.266 (3.163–3.345, odds 0.999),
  4.2 z_b 0.209 (0.164–0.503, odds 0.56).
- Segment classes (reported, not read because P failed): 112 `consistent` (z_b 4.636, 0.742–4.919, odds 0.617),
  127 `excluded` (z_b 0.624, 0.534–0.661, odds 0.922), 131 `excluded` (z_b 3.131, 2.938–3.325, odds 0.907).

## Post hoc (not in the plan)
- **Offset:** the matched controls' median RELICS − SLICE offset is (−0.238″, +0.277″) in (RA cos δ, Dec)
  (`summary.json`). With it removed (`posthoc.json`, `scripts/rxj0152_fam4_photoz_posthoc.py`), gate P still fails,
  6 / 14; the three segments match the same RELICS rows within 0.04″ and keep their classes.
- 131's interval (2.938–3.325) overlaps 4.1's own (3.163–3.345), while 4.2 sits at z_b 0.209. Because the two
  images of one source get incompatible BPZ redshifts, this overlap is not evidence for H_img.

## Reading and next
RELICS HST BPZ does not reproduce the SLICE family redshifts on this field at the frozen tolerance, the
spectroscopic family included, so it cannot decide H_img vs H_field for segments 112 / 127 / 131. The line stays
paused (D-034 addendum bg "Revisit if"): it now needs JWST-resolution photo-z (e.g. a DJA or VENUS catalogue of RX
J0152, none public in the DJA v7 index on 2026-10-26) or spectroscopy, or a SLICE model of RX J0152.

## Limits
One catalogue, one match radius (ASSUMPTION 0.6″), the 95 % interval as the consistency window (ASSUMPTION). SLICE z
of families 2–6 are model / photometric values, so gate P mixes BPZ errors with errors of the reference z; the one
spectroscopic family fails on its own (1.1 at odds 1.0).

## Deviations
- Post hoc offset-corrected re-match (above); frozen outcome unchanged.
- The first run read the catalogue under a scratch file name; it was re-run once under the archive file name
  (same sha256) so `summary.json` records the archive name. Results identical (not stored separately).
