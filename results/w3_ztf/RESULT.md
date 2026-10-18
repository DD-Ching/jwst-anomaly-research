# W3 (negative-mass point lens) feasibility trial on ZTF DR24 — result: **failed the frozen validation; no holdout search**

Plan `PLAN.md` (committed before any light curve was read; one labelled deviation: the per-epoch limit). Code
`scripts/w3_ztf.py`, rules unit-tested (`tests/test_w3_ztf.py`) before the run. Every number is `derived` from the
files named (`simulated` for injections).

## Data actually used
- ZTF DR24 light curves, IRSA HATS collection (README version date 2026-05-13; 5,401,134,101 rows), read anonymously
  from s3://ipac-irsa-ztf/ztf/enhanced/dr24/: the objects partition Norder 4 / Npix 959 (13 metadata columns, 3.49 M
  rows, 70 s) and the light-curve partition Norder 5 / Npix 3839 (953 MB; 38,165 rows kept for the sample and the
  comparison stars, 61 s; `fetch_info.json`). Nothing else downloaded.
- Field near RA 272.65°, Dec +64.95° (ZTF fields 826, 848, 1856, 1857). Sample (`sample.json`): 1,000 isolated
  star groups with r reference mag 14.5–17.5, ≥ 300 good epochs, drawn with seed 7 from 4,601; 514 development, 486
  holdout; 5,736 objectids (all filters and fields of each star in one group).

## What the data can and cannot show
- Light curves hold **detections only** (no NaN or upper-limit rows); a missing epoch is a non-detection only by our
  rule (≥ 50 % of the quadrant's comparison stars detected with catflags 0 at that epoch). No flux is known for a
  non-detection, so a vanished star is a run of missing epochs and only partial dims above the limit have
  magnitudes. Efficiencies below are for this stage on these light curves, **not** for ZTF's catalogue matching.
- Cadence: median spacing between epochs 1.0–4.0 d (10th–90th percentile of comparison stars); median 1,414 good
  r-band epochs per development star (`dev_diagnostics_posthoc.json`); seasonal gaps of months.
- The detections-only format produces non-detections on good epochs without any astrophysical change: median 3.4 %
  of good epochs per star (90th percentile 8.8 %, 99th 51 %); 117 of 514 development r-band series exceed 5 %.

## Validation on the development half (`search_dev.json`, `injections_dev.json`)
- **W3 recovery** (flagged and vetted; 200 injections per cell; blend 1 = unblended):

| ρ | t_E 3 d | 10 d | 30 d | 100 d |
|---|---|---|---|---|
| 0.01, blend 1 | 0.355 | 0.560 | **0.600** | **0.595** |
| 0.1, blend 1 | 0.290 | 0.385 | **0.505** | **0.215** |
| 0.01, blend 0.7 | 0.315 | 0.495 | 0.535 | 0.575 |
| 0.1, blend 0.7 | 0.325 | 0.475 | 0.445 | 0.205 |

- **Ordinary controls injected** (200 each): periodic eclipses 0 pass vetting (123 flagged); dust/R CrB-like fades
  10 pass (0.05; 184 flagged); ordinary positive-mass microlensing 0 (42 flagged); common-mode dims 0 (148 flagged).
- **Real development light curves:** 1,023 series searched, 1,753 flagged runs in 129 of 514 groups, **6 groups pass
  vetting**.
- **Frozen pass conditions:** W3 ≥ 0.5 at t_E 30 and 100 d for both ρ (blend 1) — **fails** (ρ = 0.1, t_E = 100 d:
  0.215); each control ≤ 10 % — passes; real survivors ≤ 5 — **fails** (6). So no holdout search was run.

## Why it failed (post hoc diagnostics, labelled)
- **The real survivors are shared-epoch artefacts:** the 6 passing groups fall on two epoch windows (MJD 58657–58663:
  2 groups; 59106–59109: 3 groups; plus 58457–58510: 1). Half of all 1,753 flagged runs start on epochs shared by ≥ 5
  groups (`dev_diagnostics_posthoc.json`: 859 runs on 91 such epochs). The quadrant-level common-mode test (< 20 % of
  comparison stars dim) does not catch them; a sample-wide shared-epoch test (as learned in MOA) is missing.
- **W3 losses** (`injection_failures_posthoc.json`, separate seed): for ρ = 0.1, t_E = 100 d the finite-source edges
  last ~ρ t_E ≈ 10 d and leave intermediate detections, so the sharp-edge test (built against dust fades) removes
  86 of the flagged injections; in every cell about a third are removed by the one-off test because the carrier's own
  detections-only non-detection runs count as "another dip". Examples: `injection_failure_examples.png`.
- So the two failures are (1) a known missing test (shared epochs across the sample) and (2) a real overlap between
  the model's slow edges at large ρ t_E and gradual ordinary fades, plus the non-detection noise of the data format.

## What this does and does not show
- W3 events with t_E ≈ 10–100 d and small ρ, if present in such light curves, would be flagged and pass vetting in
  about half the cases at this stage; the cadence covers the umbra (days to a year) but often not the spikes.
- Not shown: any limit on W3 events (no holdout search, small sample, stage-only efficiency); not a statement about
  wormholes in general (the model is a negative-mass point lens).
- Cost (measured): metadata 24–70 s, light curves 61 s for one 953-MB partition (~38 k objectids), development
  search 52 s, 4,000 injections ~15 min on one core (session observation).
