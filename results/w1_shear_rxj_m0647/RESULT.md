# W1 shear screen (D-053) on RXJ2129 and MACS0647: result (2026-10-17)

Plan: `PLAN.md` (committed before the run). Labels: shapes `observed`; ε, S maps, R `derived`; cluster shear
`model_prediction` (Lenstool maps of Zenodo 21460864 for RXJ2129, CLASH Zitrin v2 PIEMD + eNFW for MACS0647);
injections `simulated`; limits `derived` from `simulated`; every threshold is an ASSUMPTION (D-053 defaults,
unchanged).

## Inputs (not stored; streamed to the gitignored cache and deleted after the run)
F200W pipeline catalogues `jw02767-o002_t001` (RXJ2129, sha256 `2ed1a2c7…`) and `jw01433-o010_t001` (MACS0647
epoch 1, `473e32d2…`), both matching SOURCES.md; DJA v7.0 zouts `rxj2129` (`651c0b24…`, matches SOURCES.md) and
`macs0647` (`ceb199c6…`, first recorded here). Gaia DR3 G < 17 within 4′ of the model reference point: 25 / 28
stars. Fields enter `inject_shear.py` as shear-only entries (`SHEAR_FIELD_OVERRIDES`). Commands:
`exotic_screens.py --model <model> --out <dir> shear --catalog … --photoz … --spike-stars … --max-radius 100`
(`screen_shear_*.json`) and `inject_shear.py --fields <field>` (seed 50, 200 per mass; `inject_summary_*.json`).

## Screen
| Field | rows used (weak shear) | R | S_max (p_rot) | max S_× / max −S_× (p_rot) | D-053 rule |
|---|---|---|---|---|---|
| RXJ2129 | 289 | 0.551 ± 0.076 | 2.91 (0.375) | 3.00 (0.245) / 2.90 (0.405) | not passed |
| MACS0647 | 620 | 0.231 ± 0.070 | 2.97 (0.325) | 3.03 (0.21) / 2.96 (0.33) | not passed |

- Known-case check: R is positive in both. RXJ2129's 0.55 lies inside the 0.41–0.60 of the ten earlier fields.
  MACS0647's 0.23 lies below it (2.6 σ below 0.41 using its own error; not stored); the plan named no action for
  this case. Its cause was not measured (the CLASH model is HST-era; no frame offset to the JWST catalogue was
  measured).
- Sheet looked at (`sheet_rxj2129_macs0647_smax.png`, the two highest-S grid centres per field, F200W 30″): the
  RXJ2129 centres show ordinary galaxies, no streak. Both MACS0647 centres lie on two crossing linear streaks
  that span the cutout (by eye, not measured; origin not determined). Neither field passes the rule, so no
  vetting test was run. With Abell S1063 and MACS1423 this is the third field where the top W1 peak sits on a
  linear streak.

## Injections (`inject_summary_*.json`; recovered / 200)
| \|M\| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| RXJ2129 | 0 | 0 | 8 | 53 | 92 |
| MACS0647 | 0 | 0 | 0 | 1 | 1 |

MACS0647 is close to blind (1 / 200 at 2 × 10¹³ M☉), with the lowest R of all fields; no causal link was tested.

## Limits (`combined_11fields.json`, `combined_10fields_no_macs0647.json`)
- Check: `--combine-only` on the nine committed headline summaries reproduces
  `results/w1_shear_more/combined_9fields.json` exactly (equality checked in the session; not stored).
- **Headline, eleven fields (the nine plus RXJ2129 and MACS0647, as planned), 64.0 arcmin²:** Σ(W1) < 4.5 × 10³ /
  7.6 × 10² / 3.8 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (border-corrected 5.0 × 10³ / 8.4 × 10² /
  4.3 × 10²); was 4.9 × 10³ / 8.2 × 10² / 4.1 × 10² from nine fields.
- Without MACS0647 (ten fields, 59.3 arcmin²) the limits differ by < 0.2 % (`combined_10fields_no_macs0647.json`),
  because its efficiency is near zero; the low-R question does not move the headline.
- MACS1423 stays out (`results/w1_shear_more/RESULT.md`). Still blind at ≤ 2 × 10¹¹ M☉ (one recovery in 2,200 at
  2 × 10¹¹ M☉, from Abell 370).

## Limits of this result
- Same caveats as D-053 (catalogue isophotal moments, no PSF-anisotropy model, one B map per field).
- Both catalogues are shallow (about 2 ks; 289 and 620 weak-shear rows), so the two fields add little efficiency.
- No frame offset between either model and the JWST catalogue was measured.
