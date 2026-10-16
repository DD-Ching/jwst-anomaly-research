# `arcs` on MACS0647: result (2026-10-16)

Plan: `PLAN.md` (frozen before the run; outcome rule `outcome()` in `scripts/arcs_macs0647.py`, unit-tested
first). Labels: shapes `observed`; predicted directions `model_prediction`; classes, counts and P values
`derived`. All thresholds are the D-024 defaults (ASSUMPTIONs). Inputs and checksums: SOURCES.md "MACS0647".

## Model check (`epoch1/<model>/validate.json`)
Both CLASH Zitrin v2 map models reproduce their published products: κ from the deflections vs the published κ,
median relative difference 2.2 × 10⁻⁵ (PIEMD + eNFW) and 2.7 × 10⁻⁵ (LTM-Gauss); z = 2 μ 5.7 × 10⁻⁶ and
9.2 × 10⁻⁶.

## Result (`summary.json`; frozen outcome: **no result**)
| Run (epoch 1, F200W) | Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) | P(≥ anti) |
|---|---|---|---|---|---|
| primary `macs0647-clash-nfw` | 15 | 10 (67 %) | 2 | 0.646 | 0.137 |
| second `macs0647-clash-ltm` | 14 | 7 (50 %) | 0 | 0.23 | 1.0 |

Both controls fail, so under the frozen rule MACS0647 gives **no result** on the S1063 question: neither a
replication nor a non-replication. No anti count above is interpreted.

## Post-hoc and deviation numbers (`posthoc.json`, `epoch2/`; none changes the outcome)
- **Frame offset (post hoc, after review):** no frame offset between the CLASH maps and the JWST catalogue was
  measured (no image list). Shifting the model frame by 0.3″ along ±RA and ±Dec moves the primary control between
  10 / 15 (67 %) and 10 / 14 (71 %): the primary failure is marginal (one row) and passes under one of four shifts.
  Its background anti count stays 2 in all four (P not computed for the shifted runs). LTM-Gauss stays at
  43–50 %.
- **Epoch 2 (deviation: the plan pre-registers only matching the deciding run's anti rows to epoch 2; full
  epoch-2 counts were added):** primary 13 / 18 aligned (72 %), anti 2 vs 0.425 (P 0.068); LTM-Gauss 11 / 18
  (61 %), anti 1 vs 0.315 (P 0.27). Epoch 2 does not decide (the plan names epoch 1).
- **By g_min (post hoc):** epoch 1 primary, g_min < 1: 10 / 13 aligned, 1 anti; g_min ≥ 1: 0 / 2 aligned, 1 anti.
- **The two epoch-1 primary anti rows (post hoc):** 753 (g_min 0.44, F200W 26.3) has a mutual epoch-2 row 0.05″
  away that is also anti (872); 705 (g_min 6.87, 25.4) has no selected epoch-2 row within 0.3″ (nearest 2.2″). No
  conclusion is drawn about either row: neither was looked at in cutouts or vetted (the plan sends cutouts only from
  a deciding run).
- **Photo-z coverage (post hoc, not stored):** 31 of 62 selected epoch-1 rows have no DJA match within 0.3″ and fall
  into population `unknown`; the median DJA-minus-pipeline offset of the rows matched within 1″ is (+0.01″, +0.03″),
  so the missing matches are not a frame offset between the two catalogues. Small samples (14–18 background rows)
  follow from the 2,104 s exposure, r ≤ 50″ and this coverage.

## Reading
- With HST-era CLASH models (image-plane rms 1.4″ for LTM-Gauss per its readme; not stored) the predicted stretch
  does not reach the 70 % alignment control on JWST F200W shapes in this field (67 % and 50 %; the 67 % is within
  one row of the threshold and offset-sensitive, see above). All three HST-era map models tried in the
  two replication fields have now failed the control (MACS0417 RELICS v1 64 %, #280; CLASH PIEMD + eNFW and
  LTM-Gauss here); the one JWST-era model tried passed (MACS0417 CANUCS, 76 %). The HST-era CATS v4.1 model of
  S1063 passed there (20 / 24, #273), so the control failure is a per-model, per-field outcome, not a rule.
- The S1063 excess (#273, #279) stays a single-field exploratory result: one independent field does not replicate
  it (MACS0417, CANUCS), and this one cannot test it.

## Data handling
Twelve CLASH FITS maps (70.6 MB each, uncompressed; not stored: deflections, κ, γ1, γ2 and z = 2 μ of both
models) and two readmes were cached in the gitignored data root and deleted after the run; the deflection, κ and
μ maps were fetched again for the post-review re-run (identical `arcs.ecsv` tables) and deleted. The photo-z tarball
was streamed and not stored.
