# MACS0416: MARS-only images in every CANUCS family (D-034 addendum av)

Plan frozen in 8cb26331; rule functions (`split_families`, `gate_r`, `sheet_rows`) and tests in 874e50ff, both before
the run. Numbers from `summary.json` (frozen run), `controls.ecsv` / `targets.ecsv` and `posthoc.json` (post hoc).
Labels: CANUCS image positions and catalogue fluxes are observed; MARS and CANUCS-model positions, |μ| and critical
curves are model_prediction; F_pred is derived; the 3.0″ pairing radius, 1.5″ search radius and classify / local-noise
thresholds are ASSUMPTIONs (PLAN.md).

## Frozen result
- **Gate R passes:** the 40 addendum-as families give one target, D115a at (64.044742, −24.060901), |μ| 29.2, as
  addendum as.
- 71 new families (no CANUCS extra prediction) solved; **6 MARS-only targets** in 4 families (D15 × 2, D31b, D31c,
  D34 × 2).
- **Gate C2 passes** over all 111 families: 195 / 205 testable controls recovered or confused (0.951), arc-referenced
  60 / 63 (0.952). New families alone (reported): 156 / 164, 48 / 50.
- Targets, `class_local`: 2 absent, 1 confused, 1 undetectable, 2 undetectable_local. Readings: **2
  `lensing_residual_input`** (D15 at |μ| 11.2, F_pred / σ_loc 69.7; D34 at |μ| 66.8, F_pred / σ_loc 27.5), 4
  `no_reading`.
- Cutout sheet `targets.png` (both readings with their families' observed images) looked at: both target circles are
  empty sky at F200W; no edge, spike or masked pixels in them.

## Post hoc vetting (`posthoc.json`; not frozen rules)
| | D15 target | D34 target |
|---|---|---|
| CANUCS extra predictions in the family | 0 | 0 |
| Nearest observed image | D15.2, 12.33″ (all three paired to MARS) | **D34.2, 3.12″, not paired to MARS** |
| MARS images per source (each observed image's own β) | 5, 5, 5; one 0.05–0.13″ from the target (\|μ\| 10.3–11.9) and one 3.84–5.00″ away (\|μ\| 43.8–64.0) | 5, 4, 4; one 0.05–0.22″ from the target |
| CANUCS-model images per source (0.25″ grid) | 3, 3, 3; nearest to the target 10.89–12.33″ | 3, 3, 3; nearest 0.74″ (β of D34.1), 3.12″ (D34.2), 1.39″ (D34.3) |
| Other MARS-only image of the family | \|μ\| 51.3, 4.43″ away (`confused`) | \|μ\| 12.7, 9.86″ (`undetectable`) |
| MARS critical curve (sign change of det A) | 1.46″ | 1.14″ |
| CANUCS-model critical curve | 0.67″ | 0.64″ |

- **D34:** the target is 3.12″ from catalogued D34.2, the one D34 image MARS does not pair (radius 3.0″). The CANUCS
  model's own image in this region moves 0.74–3.12″ from the target with the choice of source (table). Reading
  (hypothesis, ordinary): it is MARS's image of D34.2, displaced just beyond the pairing radius, as the D115a target of
  addendum as was for D115a.2 (3.24″). Not a missing image.
- **D15:** MARS predicts two extra images (|μ| 51.3 and 11.2, 4.43″ apart) for each of the three source choices
  (table), the target 1.46″ from the MARS critical curve, about one map pixel (1.515″). The CANUCS model also has a
  critical curve near the spot (0.67″; whether it is the same curve is not tested) but gives 3 images for each source,
  none within 10.89″. F200W shows no source within 1.5″ at F_pred / σ 69.7. Reading (hypothesis, ordinary): the
  extra pair is a MARS-specific prediction near its critical curve; the data agree with the CANUCS model there.
  Vetting records: `docs/candidates/macs0416-mars-D15-extra.md`, `docs/candidates/macs0416-mars-D34-extra.md`.

## Reading (narrow)
Over 111 MACS0416 families, the two MARS-only predictions read `absent` are one pairing-radius case and one MARS-only
extra pair that the CANUCS model does not make. The other four targets (D15 |μ| 51.3 `confused`; D31b, D31c
catalogue-`recovered` but `undetectable_local`; D34 |μ| 12.7 `undetectable`) have no reading under rule 5 and were
not vetted. No candidate. The two-model MARS / CANUCS route on MACS0416 is complete.

## Limits
- Single band (F200W); MARS |μ| near critical curves is coarse (1.515″ pixels), so F_pred there carries an
  uncertainty this stage does not quantify.
- Rule 1 drops families whose images do not all back-trace; none was dropped here (`posthoc.json` `families`: 111
  total, 111 solved).
- The D15 test shows only that the data agree with the CANUCS prediction there; it does not test the CANUCS model
  elsewhere.
