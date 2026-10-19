# W3 dimming screen, Abell 2744: forced stage calibrated (D-052 addendum 2026-10-19)

Plan frozen first (`PLAN.md`). Code and thresholds unchanged (`scripts/dimming_screen.py`, `Params` defaults);
`fetch`, `screen`, `forced --n-controls 300` (one sequential S3 byte-range job, 1,347 cutout files, 666 s wall), `inject`.
Outputs here: `forced_summary.json`, `inject_summary.json`, `limits.ecsv` (derived; injections simulated).

## Result
- **Screen (derived):** reproduces D-052 exactly: 530 catalogue flags (525 / 49 / 6), 32 after the catalogue
  vetoes, 63 monitored compact sources (62 unflagged).
- **Forced stage: calibrated** (the D-052 gate: ≥ 50 controls in some epoch pair, zero points in every epoch).
  - Controls with finite forced flux per epoch (o001 / o002 / o006): F200W 39 / 295 / 252, F444W 36 / 290 / 259. Pairs:
    o002–o006 247 / 256 (pass); o001–o002 39 / 29 and o001–o006 2 / 1 (below 50: the o001 imaging footprint
    barely overlaps the two parallels). So the noise scale rests on the o002–o006 pair alone.
  - Noise scale 4.05 (F200W) / 7.21 (F444W), against 2.15 / 1.99 in MACS0416. Forced / catalogue fractional
    scatter 1.00 / 3.60. Observed, not explained here (possible causes: the 0.15″ aperture on the
    parallels' dither and depth). It makes this field's forced confirmation less sensitive, which the injection
    carries through the F444W inflation.
  - 32 flags measured: **0 forced-confirmed**, 0 controls flagged. **0 surviving events.** No contact sheet
    (nothing confirmed); no SIMBAD/NED or TNS check needed.
- **Injection (simulated):** with the measured inflation, W3 efficiency over all magnitudes is 0.000–0.013
  (ρ = 0.01 / 0.1; e.g. t_E = 0.3 yr: 0.011 / 0.013). The 2026-10-08 uncalibrated values (0.00–0.06) were
  about 4–7× too high at t_E ≥ 0.3 yr, because they used no forced inflation.
- **Headline limits, MACS0416 + Abell 2744** (ρ = 0.1; exposures add; MACS0416 exposures as rounded in
  docs/exotic_limits.md, so the sums are good to ±0.05 source yr):

  | t_E (yr) | Abell 2744 exposure (source yr) | MACS0416 | sum | per source per yr (was) | τ | limit tighter by |
  |---|---|---|---|---|---|---|
  | 0.01 | 0.09 | 1.2 | 1.29 | 2.3 (2.5) | 0.073 | 7 % |
  | 0.03 | 0.56 | 4.7 | 5.26 | 0.57 (0.64) | 0.054 | 11 % |
  | 0.1 | 1.18 | 9.9 | 11.08 | 0.27 (0.30) | 0.085 | 11 % |
  | 0.3 | 2.35 | 18.7 | 21.05 | 0.14 (0.16) | 0.13 | 11 % |
  | 1 | 2.87 | 53.7 | 56.57 | 0.053 (0.056) | 0.17 | 5 % |
  | 3 | 6.18 | 62.6 | 68.78 | 0.044 (0.048) | 0.41 | 9 % |

## Limits of this result
- Same caveats as D-052 (isolated lens, no macro-magnification; catalogue-level injection; the `single_epoch`
  veto still blinds full vanishes of two-epoch sources).
- The noise calibration uses one epoch pair; o001 enters the light curves with zero points from ≥ 10 controls
  but no pair-level noise check of its own.
- The per-deg² headline (which needs the area exposure of both fields) was not recomputed from a joint
  `combine`; it needs MACS0416's `limits.ecsv`, which a fresh clone must regenerate (`forced`, ≈ 32 min).

## Next
- NEXUS (139 sources, uncalibrated efficiency 0.07–0.24): run `forced` with ≥ 300 controls in a run with ≈ 40 min
  free; it is the field that could change the headline by a large factor, if its epoch pairs reach 50 controls.
- Then `combine` over all three fields' regenerated `limits.ecsv` for the per-deg² headline.
