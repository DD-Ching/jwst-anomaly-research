# MACS0717: D-041 cross-model readings with the pinned Keeton v4 and Sharon v4cor — result

Run 2026-10-21 (cloud routine) under the frozen `PLAN.md`. Script `scripts/m0717_crossmodel.py`; every number below is
in `summary.json` (this directory) or in `results/lens_consistency/macs0717-cats/images_predicted.ecsv`. Positions and
μ are `model_prediction`; catalogued positions `observed`; forced classes `derived`. Thresholds are ASSUMPTIONs
(PLAN.md).

## Frozen outcome: tool not validated (both models fail the known-case gate at the CATS redshifts)

| Model | catalogued images reproduced / tested (CATS `z_sys`) | gate ≥ 0.8 | secondary: systems with an `arcs.txt` z |
|---|---|---|---|
| `macs0717-keeton` | 72 / 132 (0.55) | fails | 29 / 31, passes |
| `macs0717-sharon` | 63 / 106 (0.59) | fails | 24 / 25, passes |

So, as frozen, no pinned D-041 reading is used: all 159 non-`observed` rows are `untestable` under the frozen rule, and
the Step 6 verdicts in `docs/fields/macs0717.md` keep their scratch-code basis. The gap between the two gates is
consistent with the redshift source: 41 of 51 systems take z from CATS `params.txt` `z_m_limit` (D-035), a CATS-fitted
value that need not suit another model; with catalogue redshifts the models reproduce 0.94 (Keeton) and 0.96 (Sharon) of the images. This
reading is untested beyond these two counts.

## CATS rows (reproduction check)

- Forced classes equal Step 4 exactly: recovered 25, confused 38, absent 6, undetectable 13, no_reference 54,
  inconsistent_reference 5, off_image 3 (ambiguous 0).
- **Deviation:** image classes differ by one row: observed 41 (Step 4: 40), 200 rows (199); the other classes are equal.
  The 2026-10-08 table was not stored, so the extra row cannot be identified. The new table is used.
- Run time 179 s (not stored), F200W catalogue and forced F277W `_i2d` by S3 byte range.

## Post hoc (deviation): the frozen rules on the 10 systems with an `arcs.txt` redshift

Both models pass the secondary gate there, so they are used (`summary.json` `post_hoc_arcs_z`): 34 rows, of which 14
are flagged.
- `explained_model_dependent` 12, `explained_step6` 2, `no_flag` 20; **`to_vet_candidate` 0.**
- The two second-model-supported flags: 6 (`absent`; Keeton has an extra image 2.12″ away; Sharon does not cover it) is
  Step 6's copy of 6.3; one confused row of system 1 (both models within 3″) is a Step 6 copy-system row.
- **Comparison with the scratch readings:** the four Step 6 "CATS-only extra image" rows in these systems (3, the 19
  extra, 68 ×2) read `model-dependent (CATS only)` again (`summary.json` rows). The other two (45, 71) are in systems
  without an `arcs.txt` z and are not read post hoc.
- **One disagreement:** Step 6 says Sharon and Keeton put a high-|μ| image 3.3″ / 1.3″ from the near-critical 68
  (−251.9) row (scratch, each model's own back-traced source). The pinned reading finds their nearest extra image
  93.98″ / 95.04″ away and reads it `model-dependent (CATS only)`. A near-critical image pair can merge or be assigned
  to a catalogued image under a different source position; not diagnosed. The Step 6 verdict (μ untestable, D-036)
  does not depend on it.

## Conclusion

No lensing residual and no candidate. The frozen test could not be validated on MACS0717 at CATS-fitted redshifts; on
the systems with catalogue redshifts (post hoc) the pinned models reproduce the Step 6 cross-model verdicts except
the 68 near-critical reading above.

## Limitations

- 41 of 51 systems lack a redshift in `arcs.txt`; testing them needs each second model's own source redshift, which
  this run does not read.
- Sharon v4cor's gate tests 106 of 132 catalogued images: systems with an image off its maps are skipped.
