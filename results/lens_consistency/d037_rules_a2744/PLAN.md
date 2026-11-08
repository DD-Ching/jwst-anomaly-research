# D-037 rules as functions: second-field gate on Abell 2744 (plan, frozen before the run)

**Question.** Do `lens_consistency.repeated_pair_topology` and `mu_model_dependent` (D-037 addendum 2026-11-08; gate R
passed on MACS1149, the defining field) reproduce hand readings on a field they were not built on? A pass is the
condition D-037's addendum sets before the labels become a reported column.

**Field and rows.** Abell 2744, the stored 2026-10-14 D-041 cross-model run (`results/lens_consistency/abell2744-uncover2/`):
two forced tables on Bergamini+23's image list, one per model (Bergamini+23 Lenstool, UNCOVER v2.0 maps). Each
model's `absent` / `confused` rows get both functions; the other model is the second model (its nearest *extra*
image, `crossmodel.ecsv` `other_sep_arcsec` / `other_magnification`, as on MACS1149).

**Inputs.**
- Forced tables `images_predicted_{bergamini23,uncover2}.ecsv`; `crossmodel.ecsv` (stored).
- `shared_match`: UNCOVER model from the stored `imageplane_bergamini23.ecsv` (UNCOVER on Bergamini's list). Bergamini
  model: the recorded `validate` list (docs/fields/abell2744.md: 3.2a/3.2b, 34.1a/34.1b, 700.1a/700.1b), because
  www.fe.infn.it (Bergamini `best.par` / `obs_arcs.cat`) resets every connection from this cloud session
  (2026-11-08 16:12Z), so `validate` cannot be re-run. Deviation from the MACS1149 route; stated in the result.
- Catalogued positions: Bergamini+23 Table 2 from the VizieR mirror vizier.cfa.harvard.edu (J/ApJ/952/84/table2,
  149 images, same ids as `obs_arcs.cat`; sha256 of the TSV recorded in the summary).

**Route check (before the gate).** Every `imageplane_bergamini23.ecsv` id is in Table 2; Table 2 has 149 rows
(= Lenstool's 149 images); the docs-table positions of 3.2a/b, 34.1a/b and 700.1a/b agree with Table 2 within 0.1″.
Otherwise stop: the gate is not run.

**Gate A (hand readings, docs/fields/abell2744.md 2026-10-08 table and D-041 section).** Exact labels wanted:

| Model | Row | Hand reading | Want |
|---|---|---|---|
| UNCOVER | 3.1 `absent` | fold pair 3.1a/b merged by UNCOVER (misplaced critical curve) | `model_topology` |
| UNCOVER | 34.1 `absent` | fold pair 34.1a/b merged by UNCOVER | `model_topology` |
| Bergamini | 700.1 `absent` ×2 | wrong fitted z / model-dependent position (D-041); rows ~15″ from 700.1a/b | none (fitted-z rule is not coded) |
| Bergamini | 18.1 `confused` | member light; galaxy-scale caustic, not robust to the source position | none |
| Bergamini | 22.1 `confused` | known 1.6″ position residual of 22.1a | none |

Passes when all six rows carry exactly the wanted labels. UNCOVER's other `confused` rows (3.3, C200.3, 301.1) have no
recorded hand reading: their labels are reported and read post hoc, outside the gate.

**Disclosure.** Before freezing this plan the stored cross-model rows (other-model separations and μ) were printed,
so the μ inputs were seen; the functions were not run on this field.

**Ordinary explanations targeted.** Model topology (fold pairs straddling a misplaced critical curve), model-dependent μ.
**What changes the plan.** Pass: both labels become a `d037_label` column for future forced runs (still diagnostics,
not verdicts). Fail: the functions stay diagnostics only; the mismatches are listed and the item closes.
