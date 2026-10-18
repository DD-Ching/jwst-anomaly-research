# Confirmation round 1: plan (fixed before any new check)

Owner instruction 2026-10-18 (after #311). Budget 00:51Z → **02:51Z**. No new batch. Three questions are kept apart: is
the brightness change real, what kind of object is it, and (not addressed by any check here) anything exotic.

## A. A370 box 2 (batch 3) — matched-control checks (existing `hit_lightcurve_controls.py`, unchanged)
Positions and mosaics from `results/transients/batch_venus3/boxes/review/lightcurves.ecsv`.
- **B3-A370-box2-h010 (39.970466, −1.585253), arc pixels:** `pick-arc` on `jw01208-o002_t001_nircam_clear-f444w`,
  `--half-arcsec 16 --match 1.5 --step-arcsec 0.2 --min-sep-arcsec 0.6 --n 80 --seed 1` (the #300 settings).
- **B3-A370-box2-h003 (39.966687, −1.583019), compact source:** `pick` on the same mosaic, `--half-arcsec 16
  --flux-range 128 510` (0.5–2 × its 255 nJy F444W flux) `--n 80 --seed 1`.
- Controls and the hit are measured together by `hit_lightcurve.py` (0.15″), then `compare`.
- **Decisive pair (both):** F444W `jw01208-o002` → `jw06882-o053` (the search pair).
- **Pre-declared secondary pairs:** h010 F356W `jw05324-o002` → `jw05324-o003` (the within-programme drop);
  h003 F444W `jw06882-o053` → `jw05324-o014` (VENUS → 18 d later).
- **Rule:** with ≥ 20 controls, control tail fraction < 0.05 ⇒ "the change exceeds the scatter of matched controls"
  (supports a real change); ≥ 0.05 ⇒ "consistent with the matched-control scatter" (this does **not** prove an
  instrumental cause). < 20 controls ⇒ not decided; the missing measurement is named.
- Brightness returning to an earlier level is not taken as evidence against variability.

## B. MACS0417 h002 — is a new observation worth requesting? (no lens-model rebuild; #309 numbers)
- Executed, approved and publicly listed future observations covering the counter-image positions (position, bands,
  dates), not only public images.
- The most discriminating question, the observation needed (dates, area, bands, depth, number of visits), JWST
  visibility windows, current GO/DD rules, and whether HST, other telescopes or existing schedules can answer it.
- Depths are scaled from the measured VENUS/1208 forced-photometry errors and labelled as estimates; an ETC is used only
  if available, with inputs and outputs saved.
- The #309 708–778 d is a single-model sample range; position offsets (0.8–1.1″ from 1c.2) and model-to-model
  uncertainty are kept explicit.

## C. A one-page decision summary (`DECISION.md`)
- Choose one of: start a formal proposal / first one specific check / not worth requesting now, with reasons.
- If a proposal is worthwhile, a short science-goal draft based on the evidence, without exotic premises.
- Nothing is submitted, and no one is contacted.

## Stop
- 02:30Z: no new checks; 02:51Z stop.

## Deviation (00:54Z, recorded before any compare result)
- The compact-source pick for h003 returned 16 controls (< 20). Fallback fixed before any comparison: the same `pick` with `--half-arcsec 30` (other settings unchanged). The 16-control result is reported as not decided.

## D. Optional existing-data check (added 01:08Z, before any HST 18470 data are read; timebox 30 min)
- Public HST 18470 WFC3/IR drizzled images, visits 1–4, one filter; forced 0.2″-aperture flux (annulus 0.6–1.0″,
  random-aperture errors, as `h002_hst2016.json` in #300) at two positions: the 12781 target "Athena-2"
  (64.394375, −11.90758) and our glafic 7.6″ image (64.394334, −11.907649).
- Rule: "a new point source is seen" only if a later visit exceeds visit 1 by ≥ 5σ at a position and the stamp shows a
  point source; otherwise "not seen in visits 1–4 (to 2026-10-07)". Either outcome is a statement about existing public
  data, not a discovery claim; the other team monitors this.
