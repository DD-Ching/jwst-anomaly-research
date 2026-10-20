# W3 dimming: NEXUS-Center-Tile-6 (o015) cannot serve as the unseen k = 1 sample: not tested (R3) (2026-10-20)

Plan frozen first (`PLAN.md`). Every number below is derived and stored in `footprint.json` or
`posthoc_k0/counts.json`; thresholds are ASSUMPTIONs.

## New input (observed)
- MAST (2026-10-20 08:15 UTC): program 5105 has exactly one public NIRCam level-3 observation outside
  `configs/dimming_screen.yaml`: `jw05105-o015_t005` (NEXUS-Center-Tile-6, F200W t_exptime 1,245 s, MJD 60878.76;
  80 `_nrc` `_cal` files, 8 per detector). No NEXUS-Center epoch after o012 (MJD 61127.37) is public (not stored:
  query output in the session only).

## Frozen run
- The frozen sample (NEXUS monitored sources detected in e2 with e2 S/N ≥ 50) is **empty**: the e2 mosaic S/N of
  the 151 monitored sources peaks at 33.2 (`posthoc_k0/counts.json`), so `run_o015.py` stopped before reading any
  pixel. Decision **R3 (not tested)**; the swap keeps ε(≤ 4) = 0.412 and k = 0. The S/N ≥ 50 cut was set without
  looking at the e2 S/N distribution (a plan error, recorded as a failed approach).

## Post hoc (deviation: S/N cut removed)
- With no S/N cut (`snr_min=0`), **0 of the 151 monitored sources** fall on any of the 80 o015 exposures, at k = 0
  and at k = 1 (`posthoc_k0/counts.json`). The tool itself works on o015 (a position inside the `nrca1` polygon gets
  12 covering exposures; not stored).
- The monitored sources span RA 268.222–268.702, while 17,778 screen light-curve rows lie at RA > 268.66
  (`footprint.json`): the mosaics overlap Tile-6, but the monitored selection holds no source in that overlap.

## Failed approach
- o015 as the unseen ≤ 4 sample for the k = 1 fill: no monitored pair in its footprint. Using non-monitored e2
  detections there would change the sample definition after the fact, so it could not adopt k = 1 (exploratory only).

## Next
1. k = 1 adoption still waits for a new NEXUS-Center wide epoch (none public after MJD 61127.37) or JADES.
2. Split the held-out `inconclusive` pairs by cause (unchanged).
