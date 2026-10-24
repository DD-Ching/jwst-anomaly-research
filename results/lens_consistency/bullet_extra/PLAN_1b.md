# Bullet Cluster extra-image test, deviation 1b: CANUCS F200W mosaic WCS (frozen)

Frozen 2026-10-24, after the frozen run of `PLAN.md` stopped at rule 5's input check (`RESULT.md`) and before any
CANUCS pixel or any CANUCS flux at a predicted or control position was read. Only the mosaic **header** and the
catalogue's position columns (`X`, `Y`, `RA`, `DEC`) have been looked at. This is a post-failure deviation; it changes
only how the mosaic WCS is built. Rules 4–8 of `PLAN.md` and all their constants are unchanged.

## What the header shows (observed, header only)
`NAXIS` 10240 × 8240, `CRPIX` (5120.5, 4120.5) = the array centre, `CRVAL` (104.61608616357, −55.954676475915),
`CTYPE` RA---TAN / DEC--TAN, `PC` = diag(−1, 1) with no `CDELT`, `CD` = diag(−2.22e-5°, 2.22e-5°) and `PIXSCALE`
0.08. astropy prefers `PC` (unit `CDELT`, 3600″ per pixel). A 0.08″ `CD` would make the array 819″ wide, against a
catalogue that maps `X`, `Y` at 0.040″ per pixel over the same array extent (`RESULT.md`).
**Hypothesis:** the `CD` / `PIXSCALE` keywords were carried over from an 80 mas header; the true grid is 0.04″ with the
header's `CRPIX`, `CRVAL` and orientation.

## Frozen rules
- **W0 (candidate WCS):** `W_hdr` = the header's `CRPIX`, `CRVAL`, `CTYPE` and `PC`, with the `CD` keywords removed and
  `CDELT1 = CDELT2 = 0.04/3600` (degrees per pixel) with `PC` = diag(−1, 1). No distortion terms.
- **W1 (catalogue consistency):** map every catalogue row with finite `X`, `Y`, `RA`, `DEC` through `W_hdr` with `X`, `Y`
  read once as 0-based and once as 1-based pixel coordinates. For each reading, the sky separation to (`RA`, `DEC`):
  passes iff rms ≤ 0.01″ (0.25 pixel) **and** 99th percentile ≤ 0.04″ (1 pixel) (ASSUMPTIONs). W1 passes iff at least
  one reading passes; the reading with the smaller rms is reported as the catalogue's convention. Fails → stop.
- **W2 (pixel check of the absolute registration):** because W1 cannot tell a catalogue convention from a one-pixel
  `CRPIX` error, the image decides: the 100 catalogue rows with `FLAG_POINTSRC` true, the highest `SNR_DETIMG`, a
  `W_hdr` pixel position ≥ 10 pixels inside the array, no other catalogue row within 1.0″, and a peak pixel below
  0.5 × the array's 99.999th percentile (saturation guard, ASSUMPTION). Per source: intensity-weighted centroid in a
  9 × 9 box centred on the nearest pixel to `W_hdr`(RA, DEC) (0-based, `world_to_pixel`), after subtracting the median
  of the box's border pixels. W2 passes iff ≥ 50 sources are usable (finite, positive box sum) **and** |median Δx| and
  |median Δy| (centroid − `W_hdr` position) ≤ 0.25 pixel (ASSUMPTION). Fails → stop (no C2, no flag).
- **Use:** with W1 and W2 passing, `W_hdr` replaces `WCS(header)` in `PLAN.md` rule 5 (its 0.04″ check then holds by
  construction and is still run). `src_x`, `src_y` (the containing-source and exclusion positions) are
  `W_hdr.world_to_pixel(RA, DEC)`, as in `PLAN.md`.
- **Rules 1–3 are not re-run for their outputs.** Rule 1 is recomputed (catalogue and gold list only) so its numbers are
  stored this time; a result different from "passes" → stop. Rule 2's output is stored and sha-pinned
  (`summary_wls.json`); rule 3's frame tables are read from `frame/` as committed in 87ee789. Their inputs are
  sha-pinned, so a re-run is not expected to change them.
- Rules 4–8: `PLAN.md`, unchanged.

## What would change the plan
W1 or W2 failing → stop and record; the next step would be a TAN fit to catalogue `X`, `Y` (another deviation), not a
reading. C2 failing → as `PLAN.md` (route not validated on the Bullet Cluster; no flag is read).
