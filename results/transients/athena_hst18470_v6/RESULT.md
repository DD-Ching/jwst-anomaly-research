# SN Athena reappearance in HST 18470 visit 6: result (2026-11-24)

Plan: `PLAN.md` (committed before the run; method and rule frozen in `../athena_hst18470_v5/PLAN.md`). Script
`scripts/athena_hst18470.py` (rule `decide` / `combine`, unit-tested in `tests/test_athena_hst18470.py`). Numbers:
`forced.json` (all six F105W visits, observed), `summary.json` (rule applied to visit 6 vs visit 1, derived).
Contact sheet `stamps.png` (looked at).

## Outcome: `not seen`
| Position | F105W visit 1 (nJy) | visit 6 (nJy) | Δ / σ_d | Point source in stamp | Decision |
|---|---|---|---|---|---|
| DD 12781 "Athena-2" | 167.6 ± 31.2 | 140.5 ± 27.3 | −27.1 / 41.5 = −0.65 | no | not seen |
| RELICS glafic 7.6″ image (#309) | 10.4 ± 29.1 | 20.2 ± 35.6 | +9.8 / 46.0 = +0.21 | no | not seen |

- Visit 6 is `ifu906010`, F105W, MJD 61368.23 (2026-11-24 05:30 UTC); public on 2026-11-24 (MAST CAOM TAP
  `dbo.ObsPointing`, proposal 18470, 31 rows at 19:10Z, 30 at 18:10Z).
- Visits 1–5 reproduce `../athena_hst18470_v5/forced.json` exactly (same method, same random seed).
- Sensitivity (derived, `summary.json`): `appeared` needs Δ ≥ 207.3 nJy at Athena-2 (229.9 nJy at the glafic
  position); `not seen` holds below 124.4 nJy (137.9 nJy). The #323 estimate for the reappearing image was
  ~200–330 nJy in F150W; its F105W brightness is not modelled, so a fainter (e.g. redder or early-rising) image is
  not excluded.
- Stamps are cut in each visit's detector frame; visit 6 has another roll angle, so the lensed galaxy ~1″ away
  sits below the marker. No compact source appears at either marker.

## Meaning
- Both predictions stay open: the discovery team's window (mid-Aug 2026 to early Mar 2027) and our pre-JWST RELICS
  glafic window (2027-09-14 to 2027-11-22, #309). No lens-model residual is implied; this covers existing data to
  2026-11-24 only.

## Deviations
- None from the frozen rule. Visits 2–5 are re-tabulated in `forced.json` (not needed by the rule).

## Next
Re-run at the next public HST 18470 visit (visit spacing 18.7–26.7 d in visits 1–6, from the `forced.json` MJDs) and the JWST 12781 epochs
(2027-01/02): `python scripts/athena_hst18470.py <data_dir> <out_dir> ifu901010 <new visit>`.
