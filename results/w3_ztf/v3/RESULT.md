# W3 ZTF v3: counting fix and what the catalogue non-detections of three survivors mean — result

Plan `PLAN.md` (committed before any v3 number; deviations at its end). Development half only (already seen, not
blind); the holdout half stays unread. Thresholds unchanged. Every number is `derived` from the files named; the
image work is exploratory (3 sources, 23 exposures).

## A. Counting fix (`scripts/w3_ztf.py` v3)
- Exposures come from IRSA's metadata (`ztf.ztf_current_meta_sci`: 38,026 exposures of the 26 quadrants,
  `exposures_summary.json`); light-curve times are matched to each exposure's heliocentric mid-time within 0.001 d
  (25,289 of 21,334,085 detections unmatched, 0.12 %; `v2_vs_v3.json`). Each comparison star counts once per
  exposure; the target's own group is excluded; only comparison stars with detections before and after the exposure
  are in the denominator; a run epoch with none in play fails the common-mode test. Regression tests in
  `tests/test_w3_ztf.py` cover duplicate counting, target exclusion, out-of-span comparison and target epochs,
  infobits, detection-state priority and the injection noise.
- **Changes v2 → v3** (`v2_vs_v3.json`):
  - non-detections per development r-band star: median 0.04 % of good exposures (90th percentile 4.3 %, 99th 48 %),
    so v1's "3.4 %" was a binning artefact;
  - real survivors 9 → 8 (826210200000603 dropped; the other 8 unchanged, in the windows MJD 58657–58661 and
    59106–59109); flagged runs 1,619 → 1,520;
  - W3 recovery, all cells (blend 1 / blend 0.7):
    - ρ 0.01: t_E 3 d 0.46 / 0.355; 10 d 0.605 / 0.595; 30 d 0.69 / 0.685; 100 d 0.655 / 0.575;
    - ρ 0.1: t_E 3 d 0.41 / 0.40; 10 d 0.585 / 0.45; 30 d 0.525 / 0.475; 100 d **0.30** / 0.21;
    - some cells fell from v2 (e.g. blend 1, ρ 0.1, t_E 30 d 0.62 → 0.525; blend 0.7, ρ 0.01, t_E 100 d 0.68 → 0.575);
  - controls passing vetting unchanged: eclipses 0, fades 3 / 200, PSPL 0, common mode 0 (by construction).
  - The frozen pass conditions still fail (ρ 0.1, t_E 100 d: 0.30 < 0.5; 8 survivors > 5).

## B. The catalogue non-detections of three survivors (`image_selection.json`, `forced_photometry.json`, `stamps.png`)
- Sources fixed in the plan: 848201300006405, 826215400000661 (window MJD 58655–58661) and 826215100008147
  (59106–59109); 23 exposures chosen by the plan rule before any flux (`image_selection_script.py.txt`); IRSA IBE
  positional search (`ibe_meta.json`) and 120″ science-image cutouts for the 19 exposures that contain the
  positions: 1.34 MB in all, 55–81 kB each (IBE answered 503 or reset several times; retries succeeded). Forced
  photometry at the fixed sky position (2″ aperture, 4–6″ annulus, MAGZP, local random-aperture error); 1–5 field
  stars in the same cutouts change by 0.82–1.18 × their pre-window level where measured.
- **The 11 catalogue non-detections:**
  - **class 1, no coverage or no valid observation: 6.** For 848201300006405 the 4 window exposures (field 826) do
    not contain the position in any image (IBE position search; chip gap or pointing, not determined). One exposure
    of 826215400000661 has infobits 67108864 (its flux there is 1.07 × the pre-window median, but it is invalid by
    the rule). One of 826215100008147 has non-finite pixels in the aperture (edge or mask).
  - **class 3, star present in the image: 5**, at 0.56, 0.98 and 1.15 × (826215400000661) and 0.98, 1.02 ×
    (826215100008147) of the pre-window flux; the 0.56 case has 3.9″ seeing (a 2″ aperture undercounts; borderline,
    just under the 4″ "poor" rule). In all 5 the star was imaged in a **neighbouring quadrant** of CCD 15 (q1 vs q4),
    and the release has no objectid for the star there: no member of either group lies in that quadrant
    (`class3_other_objectids.json`), so no light curve can hold those detections.
  - **class 4, a real drop: 0.**
- For these three survivors, the shared non-detections are coverage gaps (no image, an invalid exposure, an edge)
  and quadrant changes without a light curve in the other quadrant — not dimming. The other 5 survivors were not
  imaged; whether the same holds for them, and for the sample's non-detection tail, is not established.

## Answers to the two questions
1. **What the light-curve tables are still usable for after the fix:** partial dims of detected stars that stay
   above the limit (W3 with a blend, or other dips ≥ 1 mag), at the stage efficiencies above (t_E 10–100 d:
   0.575–0.69 at ρ = 0.01; t_E 3 d: 0.355–0.46; slow-edged large ρ t_E: 0.21–0.30). A catalogue non-detection by
   itself is not established as a disappearance: in the 11 cases checked it was a coverage gap, an invalid exposure
   or a quadrant change.
2. **Can a few images fill the gaps, at acceptable cost?** For the 23 exposures checked, the frozen rules put every
   exposure in class 1 or 3 (one borderline case), at 1.34 MB and a few minutes; no real drop was present to test
   the other side. Part of the gap can be closed from metadata alone: `ztf_current_meta_sci` carries each quadrant's
   corner coordinates (ra1–ra4, dec1–dec4), so "position inside the exposure" and "imaged in a quadrant without an
   objectid" can be tested for every star without images.

## If a new independent validation is run (not started)
- Scale: 1,000 new stars in a different HATS partition, with three additions frozen first: per-star footprint
  containment from the corner coordinates (exposures where the star is off-chip or in an objectid-less quadrant are
  not counted as non-detections), merging of a star's objectids across quadrants, and IBE forced photometry of every
  non-detection exposure of each vetted survivor (here 8 survivors × ≤ 8 exposures ≈ 64 cutouts, ~5 MB).
- Costs (session observations unless named): light-curve partition 57–61 s (`fetch_info.json`: 61.0 s); exposure
  metadata ~95 s for four fields; development search ~1 min; 4,000 injections ~3–4 min on one core; cutouts
  ~71 kB on average (`forced_photometry.json`: 1,344,960 bytes for 19 cutouts).
- Gaps this would not close: full disappearance is calibrated only through images (the catalogue has no upper
  limits); t_E ≈ 3 d events are recovered at 0.36–0.46; slow-edged W3 (large ρ t_E) at ≤ 0.30. The reason for the
  slow-edge loss is a hypothesis carried over from the v1 diagnostic (the sharp-edge veto), not re-measured in v3.
  W3 remains a restricted negative-mass hypothesis; nothing here is a wormhole candidate.
