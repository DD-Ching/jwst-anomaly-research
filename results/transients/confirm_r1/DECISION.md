# Decision summary (confirmation round 1, 2026-10-18)

Evidence: `a370_controls.json`, `macs0417_feasibility.json`, `jwst_gtvt_macs0417.txt`; plan `PLAN.md`. Three questions
are kept apart: did the brightness really change, what kind of object is it, and anything exotic (no check here bears
on that, and nothing points to it).

## What is now settled
- **B3-A370-box2-h003: consistent with the matched-control scatter (rule, tail 0.61); no evidence of a source-specific
  change.** Interpretation added after the check, from the control distribution: 54 unrelated compact sources in the same
  two mosaics lose a median 15 % of their F444W 0.15″-aperture flux between 1208 and VENUS (ratio 0.852, 16–84 %:
  0.824–0.882) and h003's 0.887 lies inside that band, which points to a common flux offset between the two mosaics. VENUS → 18 days later: within the control scatter (tail 0.094).
- **B3-A370-box2-h010: existing data do not support a real change.** With 80 matched arc pixels its change is within
  their scatter (tail 0.100 for 1208 → VENUS F444W; 0.113 for the within-5324 F356W drop). This does not prove an
  instrumental cause. To say more needs a change larger than the arc-pixel scatter (calibrated |z| here 1.4–1.5; a
  clear case would need ≳ 2–3), i.e. deeper data or a PSF-matched model of the arc light; no catalogued microlensing
  event lies within 0.5″.
- **MACS0417 h002 is the lensed supernova "SN Athena"**, found in October 2025 by another team (JWST DD 12781 abstract;
  GO 4568 target 0.014″ from our position). Our detection is independent but not first. Its type will come from their
  NIRSpec PRISM spectrum.

## What is still missing, per result
| Result | Brightness change | Identity / type | Missing |
|---|---|---|---|
| Batch-1 PLCK G165 h002 | real (multi-band, multi-programme, absent in 2016) | unknown; SN-like hypothesis | a spectrum or classification; TNS query |
| Batch-2 MACS1423 h000 | real | unknown; SN-like hypothesis | spectrum; host redshift confirmation |
| Batch-2 MACS0417 h002 | real | **SN Athena** (named; type pending) | the other team's spectrum and reappearance |
| Batch-2 Bullet h019 | real at S/N 5–7 | unknown | a second epoch or spectrum |
| Batch-3 A370 h003 | consistent with control scatter; no source-specific evidence | compact galaxy z 0.36 | — |
| Batch-3 A370 h010 | not supported | on the Dragon arc | deeper or PSF-modelled data |

## Smallest new observation for MACS0417
- None needed from us. The reappearance window predicted by the other team (mid-Aug 2026 – early Mar 2027) is covered
  by HST 18470 (visits 2026-08 to 2026-12; four already public) and JWST DD 12781 (NIRCam 2027-01/02 in six filters +
  NIRSpec PRISM), at our predicted counter-image position (0.29″). Depth is not limiting (S/N ~16–73 per visit at
  1208/VENUS depth; estimate).
- Our single-model window (2027-09-14 to 2027-11-22) is later and only partly JWST-visible (to 2027-10-21); if the
  supernova has not reappeared by 2027-03, the scheduled data will say so, and only then would a later epoch matter.

## Recommendation: **not worth requesting now**
- A JWST or HST request would duplicate an approved, executing programme (duplication checks; DD urgency), and the
  object is already named and followed by its discoverers.
- What we can still contribute without new observations: a recorded, blind comparison of our model-based
  reappearance window and position with the outcome, when later HST 18470 / JWST 12781 data become public.
- **Done this round (existing data, rule fixed first):** in the four public HST 18470 WFC3/IR F105W visits
  (2026-08-11 to 10-07) no new point source is seen at the 12781 "Athena-2" position or at our 7.6″ position: changes
  vs visit 1 are +1.4σ at both (167.6 → 226.7 nJy and 10.4 → 66.9 nJy; errors ~28–32 nJy) and the stamps show only
  the lensed galaxy ~1″ away (`macs0417_hst18470/`). The expected ~200–330 nJy (F150W) would be ~5–8σ against the visit-1 difference error (~41–42 nJy), but the
  counterpart position is uncertain at the ~1″ level and the F105W brightness of a supernova may differ, so this is a
  statement about existing data to 2026-10-07, not a test failure of either prediction.
- No science-goal draft is written, because no proposal is recommended.
