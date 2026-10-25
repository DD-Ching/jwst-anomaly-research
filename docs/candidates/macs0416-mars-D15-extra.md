# MACS0416 D15: MARS-only predicted image not seen (MARS / CANUCS two-model stage)

## 2026-10-25: vetting record

- **What was flagged.** `results/lens_consistency/m0416_mars_only_all/` (plan frozen in 8cb26331, run at 874e50ff):
  the MARS free-form map (Cha & Jee 2023) predicts a fifth D15 image (z 2.1067) at 64.045460, −24.063220, |μ| 11.2,
  F_pred / σ_loc 69.7 in CANUCS DR1 F200W; rule 6 reads it `absent` → `lensing_residual_input` (gate C2 passes,
  195 / 205).
- **source_uid / run.** `macs0416-mars:D15-extra` (no catalogue source); the frozen run above, `summary.json`.
- **Score.** None: a lens-model consistency flag. Positions / |μ| are `model_prediction`; photometry `observed` / `derived`.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge, masks, spikes | A | does not explain | `targets.png` (looked at): empty sky in the 1.5″ circle; σ_loc 8.54 nJy is finite (inside the mosaic) | observed |
| Bright neighbour, deblending | A/B | does not apply | `targets.ecsv`: 3 compared sources, max F / F_pred 0.059 | derived |
| Pairing radius (a displaced image of a catalogued one) | B | does not explain | `posthoc.json`: nearest observed image D15.2 at 12.33″; all three observed images paired to MARS | derived |
| Source-position choice | C | prediction robust | MARS gives 5 images for each observed image's own β, one 0.05–0.13″ from the target | model_prediction |
| Map resolution / critical curve | C | consistent with model geometry | MARS critical curve 1.46″ away (map pixel 1.515″); partner MARS-only image \|μ\| 51.3, 4.43″ away: a fold pair | model_prediction |
| Second model | C | **explains** (as a model difference) | CANUCS Lenstool model: 0 extra predictions in D15 (`predicted_macs0416.ecsv`); its critical curve lies 0.67″ from the target, but the source is outside its caustic | model_prediction |

**Verdict:** `artifact: model-specific MARS fold pair` (the "artifact" is of the lens model, not of the data). Two
validated models share a critical curve here; only MARS puts the D15 source inside the caustic, and the F200W data
(no flux at F_pred / σ 69.7) agree with the CANUCS model. Ordinary reading (`hypothesis`): MARS caustic placement at
1.515″ resolution. Not a candidate.

### Regenerate
`SSL_CERT_FILE=<proxy bundle + certifi> python scripts/m0416_mars_only_all.py --zip DIR --cat <CANUCS cat> --image
<F200W bgsub>` (pinned files: `scripts/m0416_canucsref.py`), then `--zip DIR --posthoc`.

### Open questions
- None for this flag. A finer independent MACS0416 model (D-034 addendum ar "Revisit if") would test the caustic
  placement directly.
