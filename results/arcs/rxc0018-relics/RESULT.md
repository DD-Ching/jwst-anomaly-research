# `arcs` on RXC J0018.5+1626 with the RELICS Lenstool v1 maps: result (2026-10-20)

Plan: `PLAN.md` (committed before the run). Labels: shapes `observed`; photo-z `derived` (DJA eazy); predicted axes
`model_prediction`; classes, counts and P `derived`; thresholds are ASSUMPTIONs (D-024 defaults).

## Result (`summary.json`, `rxc0018-relics-lenstool/arcs.json`)
| Model | `validate` | Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) | P(≥ anti) |
|---|---|---|---|---|---|---|
| Lenstool v1 (only model) | passes (κ median relative difference 4.4 × 10⁻⁵) | 4 | 3 (75 %) | 0 | 0.095 | 1.0 |

**Frozen outcome: "does not replicate", but the test had no power.** Post hoc (`posthoc.json` "reachability"): with
4 strong-shear background rows a passing control (≥ 70 % aligned) allows at most 1 anti row, and P(≥ 1) under the
0.095 expectation is 0.091 > 0.05, so "replicates" could not have come out of any data on this field. The field
therefore counts as **no power** in the S1063 tally, not as a non-replication. The plan did not check before the
run that "replicates" was reachable (deviation from good practice, recorded here; future plans should).

## Post hoc (not in the plan; does not change the outcome)
- `posthoc.json`: 21 of the 30 strong-shear rows have no DJA photo-z row within 0.3″ (nearest 0.39–2.78″; the DJA
  table has 5,002 rows). Those rows are `unknown` and drop out of the test, which is why the background sample is
  so small. One `unknown` row matches at 0.02″ (population left `unknown` by the z-interval rule).
- `sheet_strong_shear.png` (F200W / F444W, program 6882 o001 i2d; blue observed axis, red predicted z = 2 axis) was
  looked at. The three aligned background rows (192, 1082, 1683) are elongated along the red axis; 731 (`mixed`)
  is faint and diffuse. The two member `anti` rows (2167, 2196; not part of the test) are compact faint sources
  with no clear elongation by eye.

## Inputs
- `jw06882-o001_t001_nircam_clear-f200w_cat.ecsv` sha256 `9554af81fe3b73dbd03b6463b793fd310485de9719d6259b7388f613fafe92d0`.
- `rj0018p1626-grizli-v7.5-fix.eazypy.zout.fits` sha256 `2a7f9f8eb852dea00f0ff165acc493667d161887f5fee3d255efeca4ac8d604a`.
- Maps: `lensmodel.RXC0018_RELICS_LENSTOOL` (four files of 144,005,760 B each, Content-Length on 2026-10-20, not
  stored; cached for the run and deleted after it).

## Reading
This field cannot test the S1063 anti excess under the frozen rule. A useful
re-test needs a deeper photo-z catalogue covering the faint F200W rows (or a JWST-era model with an image list).
