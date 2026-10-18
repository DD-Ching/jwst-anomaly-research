# MACS1423 batch-2 box-2 h000: HST archive light curve — result

Plan `PLAN.md` (committed before any HST flux; times in it are session clock, a few minutes off the git times — the
commit order is the record). **Deviations and post hoc choices (all labelled):** (1) registration of every visit to
the JWST frame; (2) an offset quality gate; (3) epoch error = larger of local error and visit scatter, single-visit
epochs scaled by the same instrument's median scatter/local ratio (WFC3/IR 1.50; ACS and UVIS 1.0); (4) three
comparisons added after the first measurement (2013 → 2022-08 F606W and F435W; UVIS F606W 2014 → 2023-12); (5) the
field-source selection is per visit (peaks above 20 × the smoothed-box noise, merged across visits at 0.15″ and kept
if seen in ≥ 2 visits), not the plan's "reference epoch" wording. Every number is `derived` from `hst_forced.json` /
`epochs.json` (registered run); `unregistered/` keeps the first run.
Target: **B2-MACSJ1423+2404-box2-h000 (215.949099, +24.086948)**.

## What was there
- 45 pipeline association products (`_drc`/`_drz`) in F435W, F606W, F775W, F814W, F850LP, F110W, F160W from
  2006–2024 were measured; all have finite data at the target (`stamps_hst.png`, looked at). `hst_obs.ecsv` is the
  MAST listing used to find them (HAP single-visit rows; the measured products are the matching pipeline
  associations). No point source is visible in the aperture at any HST epoch.
- Epochs (visits within 7 days merged): before the JWST detection (2023-02-26), ACS F435W + F606W in 2022-02 (−380 d,
  16 visits; used as the reference) and 2022-08 (−181 d, 5 visits); after it, ACS F435W 2023-08 (+180 d, 1 visit),
  17504 WFC3/IR F110W + UVIS F606W 2023-12 (+290 d, 4 visits), F110W + UVIS F814W 2024-04 (+409 d, 2 visits);
  baselines from 2006–2014 (CLASH 12790, 13386, 10493). **Four new epochs are tested** (2022-08, 2023-08, 2023-12,
  2024-04), each against a same-instrument baseline; only F814W remains cross-instrument.
- Registration: 24 of 45 visits shifted (offsets up to 0.27″; `offsets_info.json`, `offsets_used.json`); it removed
  most of the IR visit scatter (F110W 2024-04: 245 vs 492 nJy unregistered → 277 ± 17).

## Measurements (0.2″ aperture incl. host; `epochs.json` `comparisons`)
| Reference → epoch (same instrument and filter) | Days | Δ (nJy) | σ (nJy) | Result |
|---|---|---|---|---|
| ACS F435W 2022-02 → 2022-08 | −183 | −13.2 | 7.2 | no change (3σ < 22) |
| ACS F606W 2022-02 → 2022-08 | −181 | −15.6 | 6.4 | no change (3σ < 19) |
| ACS F606W 2013 → 2022-08 (added) | −181 | −2.2 | 13.7 | no change (3σ < 41) |
| ACS F435W 2013 → 2022-08 (added) | −183 | −7.9 | 21.0 | no change (3σ < 63) |
| ACS F435W 2022-02 → 2023-08 | +180 | +14.6 | 9.3 | no change (3σ < 28) |
| WFC3/IR F110W 2013 → 2023-12 | +290 | −17.0 | 22.6 | no change (3σ < 68) |
| WFC3/IR F110W 2013 → 2024-04 (ZP 1.033, 10 field sources) | +409 | −50.1 | 32.2 | no change (3σ < 97) |
| WFC3/UVIS F606W 2014 → 2023-12 (added, same instrument) | +290 | +0.7 | 13.8 | no change (3σ < 41) |
| ACS F606W 2022-08 → UVIS F606W 2023-12 (cross-instrument) | +290 | +14.6 | 9.6 | no change (limited) |
| ACS F814W 2006 → UVIS F814W 2024-04 (cross-instrument) | +409 | +69 | 72 | no change (limited) |
- Field-source zero points: only the F110W 2024-04 pair has ≥ 8 matched sources (10); every other pair is
  uncorrected (labelled). Single-visit references may still understate errors where the scaling is 1.0.
- Light curve: `lightcurve.png` (HST per filter with 2013/14 baselines; JWST F150W/F444W beside it, not on one scale).

## What this adds
- **Before:** no change at the position between −380 and −181 d (F435W/F606W, 3σ < 19–22 nJy), and the −181 d level
  matches 2013 (F606W 3σ < 41, F435W < 63 nJy). A source constant since 2013 would be absorbed into the "host" and is
  not tested; so the statement is "no new or varying source at these limits", not "nothing there".
- **After:** no change in F435W at +180 d, F606W at +290 d (same-instrument vs 2014) and F110W at +290 / +409 d against
  2013 (3σ < 68 / 97 nJy). Indicative only (different band and aperture): the JWST F150W transient part was ≈ 277 nJy
  (552 − 275, host assumed constant), so by +290 d it had faded by more than a factor ~4 in the near-infrared if its
  F110W/F150W colour is near 1 (not measured).
- **Not obtained:** no second detection, so no rise or decline curve and **no colour evolution**; the transient's
  blue flux at detection is unknown (no HST epoch near 2023-02). No ordinary-transient model comparison was run (the
  plan's condition, enough detections, is not met).

## Reading (narrow)
- Consistent with a **single episode that began after 2022-08-30 and had faded in F110W by 2023-12** (observed span
  < ~470 d, < ~310 rest-frame days at z 0.525 if the host redshift applies), for a source not much redder than the
  blue pre-detection limits allow. A supernova fits this; so does any other one-off transient of that duration.
- **Weakened:** a source that varied or newly appeared in the year before (blue bands) or stayed bright after
  (F110W, F606W), at the limits above.
- **Not separated:** supernova types; supernova vs other one-off transients (e.g. a tidal disruption event); an
  intrinsic transient vs magnification of a faint persistent source (one episode, no colours over time). Nothing here
  bears on any exotic lens.

## Archive status and next step
- Measured: every HST imaging epoch at this position from 2006–2024 in the bands compared here, including the 2013
  F110W pair and the 13386 UVIS F606W baseline; listed but not used: CLASH F105W/F125W/F140W/F160W and other single
  2012–2014 bands (pre-transient baselines in bands without a post-detection epoch), and 2013–2014 grism data.
- No public HST or JWST epoch between 2022-08-30 and 2023-02-26 or within ±180 d of the detection exists in MAST for
  this position (MAST query 2026-10-18), so the rise and peak cannot be recovered from archives. **No low-cost archive
  step remains** for this question; the identity stays open.

## Cost
- MAST queries ~1–2 min; four HST byte-range runs of 43–45 visits, ~25–40 s each; registration ~10 s (session
  observations, not stored).
