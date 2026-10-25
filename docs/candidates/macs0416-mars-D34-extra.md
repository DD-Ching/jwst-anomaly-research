# MACS0416 D34: MARS-only predicted image not seen (MARS / CANUCS two-model stage)

## 2026-10-25: vetting record

- **What was flagged.** `results/lens_consistency/m0416_mars_only_all/` (plan frozen in 8cb26331, run at 874e50ff):
  MARS predicts a D34 image (z 3.2355) at 64.047690, −24.061432, |μ| 66.8, F_pred / σ_loc 27.5 in CANUCS DR1 F200W;
  rule 6 reads it `absent` → `lensing_residual_input`.
- **source_uid / run.** `macs0416-mars:D34-extra` (no catalogue source); the frozen run above, `summary.json`.
- **Score.** None: a lens-model consistency flag (`model_prediction`).

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge, masks, spikes | A | does not explain | `targets.png` (looked at): empty sky in the 1.5″ circle; σ_loc 4.63 nJy finite | observed |
| Pairing radius (a displaced image of a catalogued one) | B | **explains** | `posthoc.json`: catalogued D34.2 is 3.12″ from the target and is the one D34 image not paired to MARS (radius 3.0″, ASSUMPTION) | derived |
| Source-position choice | C | image persists | MARS gives 5, 4, 4 images for the three observed β; one 0.05–0.22″ from the target each time | model_prediction |
| Critical curve | C | near | MARS 1.14″, CANUCS 0.64″ from the target | model_prediction |
| Second model | C | consistent | CANUCS: 0 extra predictions in D34 | model_prediction |

**Verdict:** `catalog effect: pairing radius` — the target is MARS's counterpart of catalogued D34.2, displaced 3.12″
near the critical curve (`hypothesis`, supported by the distance and the unpaired D34.2), as the D115a target of
addendum as (3.24″). Not a missing image; not a candidate.

### Regenerate
As `macs0416-mars-D15-extra.md`.
