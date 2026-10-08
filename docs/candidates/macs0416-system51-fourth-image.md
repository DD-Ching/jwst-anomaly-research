# MACS0416 system 51: predicted fourth image not seen

## 2026-10-08: vetting record

- **What was flagged.** Both lens models predict a fourth image of MACS0416 system 51 (z 4.102; CANUCS K52, z_spec
  4.1032), and it is not seen. The three catalogued images (51.1–51.3) are matched by both models.
  - CATS v4.1 (`macs0416-cats`, D-035/D-040): μ 3.78 at 64.026803, −24.069756.
  - CANUCS Lenstool best fit (Rihtarsic et al. 2025): μ 5.49 at 64.026708, −24.069919, 0.2″ from the CATS position.
- **Run.** `lens_consistency.py --model macs0416-cats images --forced-image <F277W i2d> --photoz <CANUCS DR1 zout>`
  on CANUCS program 1208 `jw01208-o004_t002`. Repository main at the D-040 merge (c23dfe5); the field run is D-042.
- **Score.** None (this is a lens-model consistency flag, not an anomaly-score ranking). The predictions are
  `model_prediction`; the photometry is `observed` / `derived`.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge, low weight, detector artefacts | A | does not explain | Inside the full-depth area of all 8 bands; no spike, saturation or edge at the position (cutouts) | observed |
| Bright neighbour light | A | does not explain | A z_spec 0.268 galaxy (CANUCS 3101008) lies 0.88″ away. After an isophote model of it (photutils, 45–52 isophotes per band) is subtracted, the 0.2″ aperture gives S/N −1.8 / −0.9 / 0.2 (F150W/F277W/F444W, CANUCS position) and −0.8 / −1.1 / −0.5 (CATS position) | derived |
| Expected brightness | — | — | From 51.1–51.3 (0.2″ apertures, flux / \|μ\| per model): 9–10σ in F150W, 24–25σ in F277W, 11–13σ in F444W (local noise 0.6–2″ from the position) | derived |
| Dust in the neighbour | C | unlikely, not excluded | The positions lie roughly along the neighbour's minor axis, about 3.6 kpc from its centre at z 0.268, off the disk plane. No colour gradient is detectable at the position (all bands ≈ 0σ) | derived / hypothesis |
| Catalogue counterpart | B | none | No CANUCS DR1 source within 0.8″. The nearest are cluster members at 0.85–1.1″ (z 0.32–0.35) and the z 0.268 galaxy | observed |
| Independent model | C | predicts it | Both CATS and CANUCS predict the image (μ 3.8 / 5.5), so it is not a CATS-only feature (unlike the system-27 extras, D-042) | model_prediction |
| Galaxy-scale sensitivity | C | **explains the prediction's fragility** | In CANUCS the image survives ±30 % changes in σ of the nearest potential (8757, σ 102 km/s, 0.02″ from the z 0.268 galaxy): it moves 0.4–0.5″, μ 5.1–5.7. With 8757 removed it vanishes: 3 images remain, none within 0.8″. The image exists because of this galaxy's potential | model_prediction |
| Galaxy redshift | C | **untested ordinary explanation** | 8757 is modelled as a cluster member (z_lens 0.396, scaling-relation mass). The galaxy is at z_spec 0.268, in the foreground. A correct treatment is multi-plane: different distances, a mass rescaled to its true luminosity, and the cluster deflection applied behind it. Neither model does this | assumption |

**Verdict:** `inconclusive: needs a multi-plane lens model with the z 0.268 galaxy as a foreground deflector`.

The image's existence in both models depends on one galaxy that both treat as a cluster member, although it lies in
the foreground. A multi-plane model, or at least a single-plane model with that galaxy's mass and distance rescaled,
is needed before the absence means anything. Nothing here points beyond ordinary lensing. The strongest statement
allowed is "unexplained under the single-plane models tested", which is a statement about the models.

### Open questions

- Does a multi-plane model (for example Lenstool's multi-plane mode, or a two-plane composition of the CANUCS model
  plus a foreground dPIE at z 0.268) still predict a fourth image near this position, and with what μ?
- Is the neighbour's mass from the member scaling relation too high for its true luminosity? It would be fainter at
  z 0.268: L lower by about 2.4×, σ by about 20–25 %. The −30 % test kept the image, so the multi-plane geometry
  matters more.
- Does any spectroscopic data cover the predicted position (MUSE, NIRSpec)? A faint, extincted image would show
  Lyα or [O II] at z 4.10.

### Regeneration

Scripts are scratch only:
- `s51.py`: solving in both models;
- `s51pert.py`: ±30 % changes and removal of potential 8757;
- `s51iso.py`: isophote subtraction and photometry;
- `s51b.py`: high-pass sheet.

Inputs:
- the CANUCS bestparam file (SOURCES, D-042);
- the CANUCS DR1 photometry catalogue;
- the CATS v4.1 maps;
- S3 byte ranges of `jw01208-o004_t002_nircam_clear-{f150w,f277w,f444w}_i2d.fits`.

Queries and downloads are dated 2026-10-08. The contact sheets (`s51_sheet.png`, `s51_iso.png`) were inspected by eye.
