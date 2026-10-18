# W3 ZTF v3: counting fix and what a catalogue non-detection means (plan; fixed before any v3 number)

Owner instruction 2026-10-18 (after #336). Budget 11:47Z → **13:47Z** (session clock), ≥ 20 min for records.
Development half only (already seen in v1/v2, so not blind); the holdout half stays unread. Thresholds of
`../PLAN.md` unchanged; no injection class dropped.

## A. Fixes (v3)
- **Exposures from official metadata:** IRSA TAP `ztf.ztf_current_meta_sci` (field, ccdid, qid, filtercode, obsjd,
  expid, infobits, maglimit) for the quadrants of the sample. Each light-curve hmjd is matched to the quadrant's
  exposure whose obsjd + 15 s, converted to heliocentric MJD at the quadrant's position (astropy), is nearest within
  0.001 d; unmatched detections are counted and reported. This replaces time binning/clustering.
- **Coverage denominator:** for an exposure, the comparison stars "in play" are the quadrant's comparison objectids
  whose own light curve has detections both before and after the exposure (others are not counted as missing).
- **Common-mode count:** each comparison objectid counts at most once per exposure (duplicates removed); the target's
  own group members are excluded from its comparison set. dimfrac = (in play − detected and not dim) / in play,
  which lies in [0, 1].
- **Good exposure:** unchanged rule (≥ 50 % of in-play comparison stars detected with catflags 0) and infobits = 0.
- **Limit for injections:** the exposure's official `maglimit` (replaces the v1/v2 magerr extrapolation; labelled).
- Regression tests that fail on the v2 code: duplicate detections in one exposure, the target as its own comparison,
  and a comparison star outside its light-curve span counted as missing.
- Re-run on the development half: search and injections (same cells, seeds 1100 + c / 1200 + k). Report v2 → v3
  changes for the 9 v2 survivors and every recovery cell.

## B. Images for the non-detections (exploratory)
- Sources (fixed now): up to three v2 survivors covering both windows — 826215400000661 (window MJD 58657–58661),
  848201300006405 (same window, other field) and 826215100008147 (window 59106–59109).
- Exposures (fixed now, before any flux): per source, from its matched exposure list in the same filter: the 2
  nearest good-epoch detections before the window, every exposure of the window (all of its v2 ND epochs, up to 4)
  and the 2 nearest detections after; at most 8 per source, 24 in all.
- Data: IRSA IBE science-image cutouts (`sciimg.fits?center=..&size=40arcsec`), cumulative download ≤ 500 MB, 60 s
  timeouts, cached in the session scratch. Forced photometry at the fixed sky position with
  `hit_lightcurve.aperture_photometry` (radius 2.0″, annulus 4–6″), zero point MAGZP from the header, local error
  from 200 random apertures in the cutout, and the same measurement on up to 5 field stars in the cutout as a
  common-change check.
- Each exposure classified: (1) no coverage or no valid observation (missing cutout, infobits ≠ 0, NaN at the
  position); (2) observed but poor (seeing > 4″ or maglimit < 19.5); (3) source present in the image (flux ≥ 5σ and
  ≥ 50 % of the pre-window median) although the catalogue has no detection; (4) valid observation with flux below
  50 % of the pre-window median or below 5σ (a real dimming or below the local limit, stated with the flux and error).

## Deviations (labelled; recorded after the run)
- Cutouts were 120″, not 40″: the 40″ cutouts held no field star for the common-change check (all three sources are
  isolated); the classification rule is unchanged.
- The image windows were taken from the v2 survivors' non-detection runs (MJD 58655–58661 and 59106–59109); the v3
  run's flagged runs start at 58657.41 and 59106.19.
- After review (same day): `star_state` returns the magnitude it chose; a run epoch with no comparison star in play
  now fails the common-mode test (was 0.0); an `exposures` stage writes the TAP query and a summary. Re-running
  search and injections gave identical numbers.
