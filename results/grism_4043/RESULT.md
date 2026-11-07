# Program 4043 F444W grism: covers one open SMACS record; pipeline c1d not usable for a line test

Plan: `PLAN.md` (frozen in 54100b2 before any grism flux was read). Scripts: `scripts/grism_4043_coverage.py`,
`scripts/grism_4043_lines.py` (frozen rule), `scripts/grism_4043_fluxcheck.py` (post hoc). Tables:
`coverage.json`, `lines_2915.json`, `fluxcheck.json`.

## Data (observed)
- MAST CAOM ObsCore (`obs_id LIKE 'jw04043%'`, 2026-11-07): level-3 `f444w-grismr` and `-grismc`, target
  `SMACSPROTO`, 8 footprint polygons each. S3 `stpubdata/jwst/public/jw04043/L3/t/o001/` holds the pipeline
  `_x1d` and `_c1d` files (byte sizes in SOURCES.md), last modified 2026-10-21, `CAL_VER` 3.0.0, extracted
  from the 4043 F444W imaging catalog (`SCATFILE`). The c1d files list 212 (R) and 221 (C) sources (not stored; read from the FITS headers in the session).
- Only the two c1d files and the 4043 F444W catalog were read (streamed to the scratchpad, deleted after use).

## Coverage (derived, `coverage.json`)
- 46.3 % of the 4080 2736 F200W catalog sources lie inside each grism footprint (necessary, not sufficient; ASSUMPTION).
- Of the nine candidate records and the two point-like top-20 sources: inside are `f200w_1096`, `f200w_2559`,
  `f200w_2915`, `f200w_2925`, `f277w_829`; outside are `f200w_2804` (verdict "needs NIRSpec"), `f200w_1032`,
  `f200w_438`, `f356w_1243`, `f200w_1816`, `f200w_384`.
- In the c1d source lists only `f200w_2915` has a spectrum within 0.15″ (0.03″, source 1238, both dispersions);
  `f200w_2925` is at 0.31″ (not stored; computed in the session), beyond the radius. So the grism route cannot
  address `f200w_2804`, `f200w_1816` or `f200w_384`.

## Line test for `f200w_2915` (frozen rule, `lines_2915.json`)
- Validation gate **failed**: injected S/N-7 lines are recovered at 0.71 (He I), 0.655 (Paγ) and 0.0 (Paβ) < 0.90.
  Per the plan: **no verdict**; the quiescent-galaxy verdict stands as recorded.
- **Deviations:** (1) bug fix after `/code-review`: an injected centre off a spectrum made `inject` return NaN; it
  now leaves that spectrum unchanged, and recovery is also reported over draws on both spectra (the gate numbers did
  not change). (2) The plan's zero-amplitude "false-positive rate" was dropped: without injection every draw
  repeats the same observed decision, so it is not a rate. (3) The rule functions were unit-tested before the run
  but committed with the result (2ad1461), not with the plan, so git cannot show that they were frozen first.
- Post hoc reasons (not tested as causes): Paβ lies on both spectra in only 9 / 200 draws (grismc ends at
  4.7579 µm), and 0 / 9 are recovered; the window peaks sit at band edges (grismr Paβ S/N 13.38 at 5.0815 µm, where
  the flux diverges), and the argmax over the window takes no edge mask, so an edge feature can outrank an
  injected line; both spectra have one narrow negative spike (spectrum plot looked at, not stored).

## Flux check (post hoc, `fluxcheck.json`)
- The c1d continuum (median, 4.2–4.6 µm) over the 4043 F444W `aper_total_flux`: median 23.1 (16–84 %: 7.5–86.4) for
  163 grismr sources and 26.4 (10.7–110.4) for 178 grismc sources; `f200w_2915`: 17.5 / 17.6.
- So the pipeline c1d levels are not on the imaging flux scale, and the spread is wide; the cause (calibration,
  extraction aperture, contamination by overlapping traces) is not tested.

## Limitations
- Coverage uses footprints and the pipeline's own source list; a decontaminated extraction (e.g. grizli from the
  rate files) could extract more sources, but none of the three uncovered records lies inside the footprint.
- No line limit for `f200w_2915` (gate failed). Closing the grism route for the open records is **post hoc**: it
  rests on the coverage (step 1) and the flux check, not on the plan's outcome rule.
