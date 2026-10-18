# W3 (negative-mass point lens) feasibility trial on ZTF DR24 — result: **failed the frozen validation (v1 and v2); no holdout search**

Plan `PLAN.md` (committed before any light curve was read; labelled deviation: the per-epoch limit; labelled
version 2: bug fixes after review, thresholds unchanged). Code `scripts/w3_ztf.py`, rules unit-tested
(`tests/test_w3_ztf.py`). Every number is `derived` from the files named (`simulated` for injections); session-clock
costs are labelled "not stored".

## Data actually used (`data_versions.json`, `fetch_info.json`, `sample.json`)
- ZTF DR24 light curves, IRSA HATS collection (README version date 2026-05-13), read anonymously from
  s3://ipac-irsa-ztf/ztf/enhanced/dr24/: one objects partition (13 metadata columns) and one light-curve partition
  (952.75 MB; 38,165 rows kept for the sample and comparison stars, 61 s). Nothing else downloaded.
- Field near RA 272.65°, Dec +64.95° (ZTF fields 826, 848, 1856, 1857; PLAN deviation note). Sample: 1,000 isolated
  star groups, r reference mag 14.5–17.5, ≥ 300 good epochs, drawn with seed 7 from 4,601; 514 development, 486
  holdout; 5,736 objectids, disjoint groups (all filters and fields of a star in one group).

## What the data can and cannot show
- Light curves hold **detections only** (no NaN or upper-limit rows); a non-detection exists only by our coverage
  rule (≥ 50 % of the quadrant's comparison stars detected with catflags 0 in that exposure). No flux is known for a
  non-detection, so a vanished star is a run of missing epochs, and only partial dims above the limit have
  magnitudes. All efficiencies are for this stage on these light curves, **not** for ZTF's catalogue matching.
- Cadence: median spacing between epochs 1.0–4.0 d (10th–90th percentile of comparison stars; PLAN note); median
  1,414 good r-band epochs per development star (v1 count, `dev_diagnostics_posthoc.json`); seasonal gaps of months.

## v1 (frozen; `search_dev.json`, `injections_dev.json`): failed
- W3 recovery ≥ 0.5 at t_E 30/100 d for both ρ (blend 1): **fails** (ρ = 0.1, t_E = 100 d: 0.215). Real survivors ≤ 5:
  **fails** (6). Controls ≤ 10 %: passes (eclipses 0, fades 0.05, PSPL 0, common mode 0).
- Review then found code bugs (PLAN "Version 2"): the most serious, fixed 0.0005-d time bins, split single exposures
  and created fake non-detections — v1's "median 3.4 % of good epochs undetected" is mostly this artefact and is
  withdrawn as a data property.

## v2 (bug fixes only; same development half, so not blind; `search_dev_v2.json`, `injections_dev_v2.json`): failed
- **W3 recovery** (flagged, overlapping the umbra, not already flagged in the carrier, and vetted; 200 per cell):

| ρ, blend | t_E 3 d | 10 d | 30 d | 100 d |
|---|---|---|---|---|
| 0.01, 1 | 0.455 | 0.560 | **0.710** | **0.615** |
| 0.1, 1 | 0.335 | 0.565 | **0.620** | **0.240** |
| 0.01, 0.7 | 0.375 | 0.610 | 0.720 | 0.680 |
| 0.1, 0.7 | 0.345 | 0.505 | 0.505 | 0.240 |

- **Ordinary controls injected** (200 each; one event model in both filters): eclipses 0 pass vetting (41 flagged);
  dust/R CrB-like fades 3 (0.015; 179 flagged); ordinary positive-mass microlensing 0 (0 flagged); common-mode 0 (147
  flagged; rejected by construction, the comparison dim fraction is set to 0.5).
- **Real development light curves:** 1,023 series, 1,619 flagged runs in 95 of 514 groups, **9 groups pass**.
- Pass conditions: ρ = 0.1, t_E = 100 d recovery 0.240 < 0.5 — **fails**; 9 real survivors > 5 — **fails**; controls
  pass. The holdout half was **not searched**.

## Why (post hoc, labelled)
- **Survivors:** all 9 fall in two windows (MJD 58655–58661: 6 groups; 59106–59109: 3), and in those windows 46 of
  509 and 51 of 499 development r-band series are mostly undetected, in fields 826 and 848
  (`survivor_windows_v2_posthoc.json`; light curves in `survivors_v2_examples.png`). So they coincide with
  sample-wide non-detection episodes, not single-star events (individual images not inspected). The common-mode test
  missed them; one reason found: the comparison-star dim fraction is negative in 1,886 of 30,838 quadrant epochs
  (a star detected twice in one exposure cluster is counted twice; v2 bug, not fixed).
- **W3 losses** (v1 diagnostic, `injection_failures_posthoc.json`, three cells, seed 99): at ρ = 0.1, t_E = 100 d the
  finite-source edges last ~ρ t_E ≈ 10 d and leave intermediate detections, so the sharp-edge test (built against dust
  fades) removed 86 of the flagged injections; in the three cells analysed about a third were removed by the one-off
  test (in v1, partly from the binning artefact). Example light curves: `injection_failure_examples.png`.

## What this shows
- At this stage, with the first two pipeline versions (both failed), W3 umbrae with t_E 10–100 d and ρ = 0.01 are
  flagged and vetted in 56–71 % of injections, and ordinary fades pass vetting at ≤ 5 %; slow-edged (large ρ t_E)
  W3 events are mostly lost to the fade veto. No real candidate, no limit, nothing about wormholes in general.
- Cost (not stored unless named): light curves 61 s per 953-MB partition (~38 k objectids; `fetch_info.json`);
  sample stage 23.8 s (`sample.json`); development search ~1 min; 4,000 injections a few minutes on one core.
