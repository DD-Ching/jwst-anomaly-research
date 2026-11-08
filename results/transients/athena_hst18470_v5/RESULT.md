# SN Athena reappearance in HST 18470 visit 5: result (2026-11-08)

Plan and frozen rule: `PLAN.md` (committed before the run). Script `scripts/athena_hst18470.py` (rule `decide`
unit-tested in `tests/test_athena_hst18470.py`). Numbers: `forced.json` (all five F105W visits, observed),
`summary.json` (rule applied to visit 5 vs visit 1, derived). Contact sheet `stamps.png` (looked at).

## Outcome: `not seen`
| Position | F105W visit 1 (nJy) | visit 5 (nJy) | Δ / σ_d | Point source in stamp | Decision |
|---|---|---|---|---|---|
| DD 12781 "Athena-2" | 167.6 ± 31.2 | 197.4 ± 23.3 | +29.8 / 38.9 = +0.77 | no | not seen |
| RELICS glafic 7.6″ image (#309) | 10.4 ± 29.1 | 24.9 ± 24.0 | +14.5 / 37.7 = +0.38 | no | not seen |

- Visit 5 is `ifu905010`, MJD 61341.56 (2026-10-28 13:26 UTC); it was the only public HST 18470 visit not checked
  in #323 (MAST query 2026-11-08; visits 1–5 public).
- Visits 1–4 reproduce #323's `hst18470_forced.json` to the stated precision (same method, same random seed).
- Sensitivity (derived, from `summary.json`): `appeared` needs Δ ≥ 5 σ_d, i.e. Δ ≥ 194.7 nJy at Athena-2; `not seen`
  holds below 3 σ_d (116.8 nJy). The #323 estimate for the reappearing image was ~200–330 nJy in F150W; its F105W
  brightness is not modelled here, so a fainter (e.g. redder or early-rising) image is not excluded.
- Stamps are cut in each visit's detector frame; the roll angle differs between visits (visit 5 is rotated), so the
  lensed galaxy ~1″ away moves around the marker. No compact source appears at either marker.

## Meaning
- Both predictions stay open: the discovery team's window (mid-Aug 2026 to early Mar 2027, three JWST-era models;
  `results/transients/confirm_r1/DECISION.md`) and our pre-JWST RELICS glafic window (2027-09-14 to 2027-11-22,
  `results/transients/batch_venus2/followup/macs0417_h002/recheck/RESULT.md`). No lens-model residual is implied.
- This is a statement about existing data to 2026-10-28, not a test failure of either prediction.

## Deviations
- None from the frozen rule. Post hoc: all five visits are tabulated in `forced.json` (the plan needed only 1 and 5).

## Next
Re-run with each newly public HST 18470 visit (cadence ~19–21 d in visits 1–5) and the JWST 12781 epochs
(2027-01/02): `python scripts/athena_hst18470.py <data_dir> <out_dir> ifu901010 <new visit>`.
