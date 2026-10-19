# NEXUS W3 dimming survivors vetted at exposure level (D-052 addendum 2026-10-19c)

Plan frozen first (`PLAN.md`). Tool: `scripts/dimming_exposure_check.py` (fixed-position 0.15″ aperture on every
level-2 `_cal` exposure of the flagged epoch that covers the position, JUMP_DET / DO_NOT_USE pixel counts, stamp
sheet). Numbers: `single_points/` and `d73390/` (`exposure_fluxes.ecsv`, `exposure_summary.json`,
`exposure_sheet.png`) and `vetting.json` (built by `make_vetting.py`). All values derived.

## Result: 0 of 5 survivors is a W3 dimming; each has a specific ordinary reading
| survivor | test | outcome |
|---|---|---|
| d53823 (e0) | dithers | F200W 93.6 µJy in exposure 1 with 60 / 72 aperture pixels JUMP-flagged; 0.0003 µJy in exposure 2 at the same sky position. One-exposure event: cosmic ray / detector residual. |
| d60156 (e3) | dithers; colour (post hoc) | one exposure covers it (no dither test); 46 / 75 pixels JUMP-flagged, compact core with a halo in the stamp (snowball-like). Forced F200W 10.44 µJy against F444W 0.038 ± 0.015 µJy in the same exposure: chromatic veto. |
| d79928 (e6) | dithers; colour (post hoc) | one exposure covers it; 10 / 74 pixels JUMP-flagged. Forced F200W 24.37 µJy against F444W 0.017 ± 0.017 µJy: chromatic veto. |
| d73390 | aperture / exposure level | on the `_cal` frames the F200W epoch medians are 5.20–6.01 µJy (min / max 0.865; e3 not measurable, DO_NOT_USE in the aperture), F444W 3.44–3.90 µJy (0.883). The i2d forced F200W range (max / min 2.03) is not reproduced at exposure level; no epoch fades in one band only. Why the i2d forced amplitude is larger is not established. |
| d79956 | catalogue; light curve | a brightening, not a dimming: F200W 0.019–0.042 µJy in e0–e4 (consistent within errors), then 0.096 / 0.160 / 0.118 in e5–e7, with F444W rising too (0.056 in e6); present at the same flux in the e6 `_cal` exposure with no JUMP pixels (a real sky source). SIMBAD: no object within 5″. NED and TNS not checked (timeout; 403 from TNS). An ordinary transient (supernova or AGN) is the likely reading; not tested further. |

- **Chromatic veto (post hoc; `chromatic_veto`, unit-tested):** F200W and F444W are exposed simultaneously, so any
  sky source appears in both, and a lensing spike is achromatic. A point is vetoed when
  `f_F444W + 3 σ < f_F200W / 20` (ASSUMPTION: the Rayleigh–Jeans limit is 5; ×4 margin for emission lines). It
  cannot remove an injected W3 event (achromatic by construction), so the injection efficiencies stand.
- **Controls (the same artefact class at random positions):** 4 of the 10 flagged controls are single-epoch F200W
  points. 3 are covered at exposure level, and each is in one exposure only with 31–48 JUMP pixels (ctl0290 181 µJy,
  ctl0294 1.13 µJy, ctl0298 0.55 µJy); 2 of the 4 meet the chromatic veto (ctl0290, ctl0298). The other 6 controls
  are faint (F200W ≤ 0.30 µJy) and were not examined at exposure level (deviation). In 5 of them, e0 reads below
  zero.

## Consequence
The condition in `../nexus_calib/headline.json` holds: the three-field limits stand as written there (0.020 per
source per yr at t_E = 1 yr, 0.029 at 0.3 yr, ρ = 0.1). The per-deg² headline still needs the joint `combine`.

## Deviations
- The chromatic veto is not in the frozen plan (post hoc). It is needed where only one exposure covers the
  position (d60156, d79928). For these two, the reading rests on it and on the JUMP flags.
- d73390: the nearby-galaxy comparison in the plan was not run. The exposure-level aperture test was enough to show
  that no single-band fade exists.
- d79956: TNS (403) and NED (timeout) were not checked. A transient reading does not change the headline, because
  the source brightens and does not vanish.
