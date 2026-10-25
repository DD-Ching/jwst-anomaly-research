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
| Source-position choice | C | prediction robust | `posthoc.json`: MARS gives 5 images for each observed image's own β, one 0.05–0.13″ from the target (\|μ\| 10.3–11.9) and one 3.84–5.00″ away (\|μ\| 43.8–64.0) | model_prediction |
| Map resolution / critical curve | C | consistent with model geometry | MARS critical curve 1.46″ away (map pixel 1.263″) | model_prediction |
| Second model | C | **explains** (as a model difference) | CANUCS Lenstool model: 0 extra predictions in D15 (`predicted_macs0416.ecsv`); 3 images for each source (0.25″ grid), none within 10.89″ of the target; a CANUCS critical curve lies 0.67″ from the target (same curve as MARS's not tested) | model_prediction |

**Verdict:** `artifact: model-specific MARS extra image pair` (an artifact of the lens model, not of the data). Of
two models validated on held-out MACS0416 images, only MARS predicts images here, and the F200W data (no flux at
F_pred / σ 69.7) agree with the CANUCS model. Ordinary reading (`hypothesis`): MARS's critical-curve / caustic
geometry at 1.263″ resolution. Not a candidate.

### Regenerate
`SSL_CERT_FILE=<proxy bundle + certifi> python scripts/m0416_mars_only_all.py --zip DIR --cat <CANUCS cat> --image
<F200W bgsub>` (pinned files: `scripts/m0416_canucsref.py`), then `--zip DIR --posthoc`.

### Open questions
- None for this flag. A finer independent MACS0416 model (D-034 addendum ar "Revisit if") would test the caustic
  placement directly.
