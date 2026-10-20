# W3 dimming: NaN-tolerant aperture in the D-039 persistence test: G1 fails at 2 pixels, k stays 0 (2026-10-20)

Plan frozen first (`PLAN.md`, amended once before G1 or any score). Code: `persistence_check.fill_nan_aperture`
and `measure_pixel(max_nan=k)` (default 0 = D-039 as frozen); `persistence_run/run.py ... max_nan=2`. Data: NEXUS
`_cal` exposures read from S3 by range (0 MB stored). Every number is derived; thresholds are ASSUMPTIONs.

## Gates
- **G0 (code identity): pass.** In e1, e3, e4, e5, e6 and e7, all 266 exposures with finite k = 0 flux have
  identical k = 2 flux, no row is missing, and the 89 pairs with no NaN exposure keep their verdict (`score.json`
  `G0`). k = 2 makes 81 NaN exposures finite; 7 stay NaN with 1–2 NaN pixels in the stamp's aperture (centre
  pixel, an unfillable reflection or the array edge; the cause is not split).
- **G1 (fill accuracy on 42 known on-sky apertures, e5; `g1.json`): fails.** 1 masked pixel: median
  |f_fill / f_true − 1| 0.037 (p90 0.126), S/N ≥ 5 call unchanged 42 / 42: pass. 2 masked pixels: median 0.079
  (p90 0.170) > 0.05, call unchanged 42 / 42: fail. No fill was refused (42 / 42 finite). The error is scatter, not bias (signed median +0.000 and
  +0.022).

**Decision: D2.** k stays 0 and the swap keeps ε = 0.412 for ≤ 4 exposures. No limit changes.

## Post hoc (G1 failed; descriptive only, not adopted)
Held-out (e5–e7) and train (e1, e3, e4; e0 not re-run) NEXUS monitored pairs, mosaic S/N ≥ 10, ≤ 4 exposures:

| aperture | train `on_sky` / n (Wilson 95 %) | held out `on_sky` / n (Wilson 95 %) |
|---|---|---|
| k = 0 (committed) | 40 / 86 (0.363–0.570) | 33 / 89 (0.278–0.475) |
| k = 1 (approx.: k = 2 tables minus 2-pixel exposures; earlier-exposure fills not undone) | 57 / 86 (0.558–0.754) | 58 / 89 (0.548–0.743) |
| k = 2 | 62 / 86 (0.618–0.805) | 62 / 89 (0.595–0.782) |

- 0 `persistence` verdicts under k = 2 (no illuminator invented by the fill); held-out `undetected` 11 → 2.
- Post hoc alternatives for 2 masked pixels (`g1.json`): bilinear reflection through the true centre 0.069,
  area correction 0.053; neither meets 0.05 either.
- Reading: the NaN-aperture loss is most of the ≤ 4 deficit, and k = 1 passed G1 as written. The 5 % flux gate is
  stricter than the D-039 call needs (S/N ≥ 5, unchanged 42 / 42; the 20× illuminator ratio was not tested here), but it
  was frozen and is kept.

## Failed approach
- k = 2 point-reflection fill: per-pixel scatter of an undersampled NIRCam PSF about the rounded centre gives
  8 % median flux scatter with two filled pixels; the unfillable centre pixel needed a plan amendment before G1.

## Next
1. New frozen plan for k = 1 (the case that passed G1). The held-out e5–e7 numbers above are already seen, so its
   efficiency check must use unseen data: NEXUS e0 (train in `expcount_eff`, 1076 exposures, not re-run here)
   and e2 controls; then the 37 targets in all three fields with k = 1 and a re-run `inject`.
2. Split the remaining held-out `inconclusive` pairs by cause (unchanged from `undetected_stamps`).
