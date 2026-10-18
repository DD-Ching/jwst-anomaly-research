# MACS1423 batch-2 box-2 h000: HST archive light curve — result

Plan `PLAN.md` (fixed before any HST flux; two labelled deviations: registration to the JWST frame, and an offset
quality gate). Every number is `derived` from `hst_forced.json` / `epochs.json` (registered run) unless labelled;
`unregistered/` keeps the first run. Target: **B2-MACSJ1423+2404-box2-h000 (215.949099, +24.086948)**.

## What was there (`hst_obs.ecsv`, `stamps_hst.png`)
- Every listed visit has finite data at the target (stamps looked at). Used: 43 association products in F435W,
  F606W, F775W, F814W, F850LP, F110W and F160W from 2006–2024. **Useful new epochs (visits within 7 days merged):**
  - before the JWST detection (2023-02-26): ACS F435W + F606W, 2022-02 (−380 d, 16 visits) and 2022-08 (−181 d,
    5 visits); CLASH 2013 as a decade-earlier baseline (F606W, F110W);
  - after it: ACS F435W 2023-08 (+180 d, 1 visit); 17504 WFC3/IR F110W + UVIS F606W 2023-12 (+290 d, 4 visits);
    17504 F110W + UVIS F814W 2024-04 (+409 d, 2 visits).
  - So **5 new epochs** (2 before, 3 after), 3 of them same-instrument comparisons; F606W 2023-12 and F814W 2024-04
    are cross-instrument (ACS → UVIS) and only limited.
- No point source is visible in the aperture at any HST epoch (stamps).

## Measurements (0.2″ aperture incl. host; `epochs.json` `comparisons`)
| Reference → epoch (same filter) | Days from detection | Δ (nJy) | σ (nJy) | Result |
|---|---|---|---|---|
| ACS F435W 2022-02 → 2022-08 | −183 | −13.2 | 7.2 | no change (3σ < 22) |
| ACS F606W 2022-02 → 2022-08 | −181 | −14.9 | 6.5 | no change (3σ < 19) |
| ACS F606W 2013 → 2022-08 (added comparison) | −181 | −2.2 | 13.7 | no change (3σ < 41) |
| ACS F435W 2022-02 → 2023-08 | +180 | +14.6 | 9.3 | no change (3σ < 28) |
| WFC3/IR F110W 2013 → 2023-12 | +290 | −17.0 | 19.6 | no change (3σ < 59) |
| WFC3/IR F110W 2013 → 2024-04 | +409 | −39.5 | 19.4 | no change (3σ < 58) |
| ACS F606W 2022-08 → UVIS F606W 2023-12 (cross-instrument, ZP 0.95 from 8 field sources) | +290 | +17.7 | 9.9 | no change (limited) |
| ACS F814W 2006 → UVIS F814W 2024-04 (cross-instrument, no ZP) | +409 | +69 | 72 | no change (limited) |

- Errors: per-visit local random apertures; within an epoch the larger of the local error and the visit scatter is
  used (choice made after seeing the scatter; labelled). Registration (23 of 43 visits shifted, offsets ≤ 0.27″;
  `offsets_info.json`) removed most of the IR visit scatter (F110W 2024-04: 245 vs 492 nJy unregistered → 277 ± 0.6
  scatter). Field-source zero points were applied only where ≥ 8 sources matched (ACS F606W pair, cross-instrument
  F606W); the IR and F435W pairs have no zero-point correction (labelled).
- Light curve: `lightcurve.png` (HST per filter, host included; JWST F150W/F444W beside it, not on the same scale).

## What this adds
- **Before:** no source at the position in F435W or F606W 181 and 380 days before the JWST detection, at the limits
  above (a constant source would be absorbed into the host, but the 2022-08 F606W level also equals the 2013 level).
  Within these bands, the source was not present for ≥ 6 months before the detection.
- **After:** nothing in F435W at +180 d, F110W at +290 d and +409 d (3σ < 59 nJy). For comparison (different band and
  aperture, so only indicative): the JWST F150W transient part was ≈ 277 nJy (552 − 275, 2023 − 2026, host assumed
  constant), so by +290 d it had faded by more than a factor ~4–5 in the near-infrared (if its F110W/F150W colour is
  near 1; not measured).
- **Not obtained:** no second detection, so no rise or decline curve and **no colour evolution**; the transient's
  F606W/F435W flux at detection is unknown (no HST epoch near 2023-02), so the pre-detection limits cannot be turned
  into a rise-time or a luminosity statement without assuming a spectrum. No ordinary-transient model comparison was
  run (the plan's condition, enough detections, is not met).

## Reading (narrow)
- Consistent with a **single episode that started after 2022-08-30 and was gone in F110W by 2023-12**: an observed
  duration of less than ~470 days (< ~310 rest-frame days at z 0.525, if the host redshift applies). A supernova is
  consistent with this; so is any other one-off transient of that duration.
- **Weakened:** a long-lived source present a year before and still bright after (e.g. a steadily bright AGN or star
  cluster) — in the blue bands before and the near-infrared after, at the limits above.
- **Not separated:** supernova types; supernova vs other one-off transients (e.g. a tidal disruption event); an
  intrinsic transient vs a magnification event of a faint persistent source (a single episode does not separate
  them; no colours over time). Nothing here is evidence for or against any exotic lens.
- Still no catalogued transient or spectrum (`../identity_r2/identity_check.json`).

## Archive status and next step
- The HST archive at this position is now used for this question (all 2006–2024 imaging epochs listed in MAST at
  2026-10-18 in the four comparison bands). Remaining public data with any power: the 13459 GLASS / 13063 grism
  observations are 2013–2014 (pre-transient) and do not help. No public epoch between 2022-08 and 2023-02 or near the
  peak exists in MAST, so the rise and peak cannot be recovered from archives. **No low-cost next step remains** for
  this object; its identity stays open.

## Cost
- MAST queries ~1 min; three HST byte-range runs of 43 visits ~25–40 s each; registration ~10 s (session
  observations, not stored).
