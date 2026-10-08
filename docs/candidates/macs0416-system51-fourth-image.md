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
| Galaxy redshift and mass (two-plane model) | C | **explains** | The galaxy is put on its own plane at z 0.268 and the rest of the CANUCS model is kept at z 0.396. Two-plane lens equation; source re-fitted from 51.1–51.3; image search on a 0.05″ grid. **σ 102 km/s (as fitted):** 5 images; the fourth persists (μ 6.1, 0.12″ from the single-plane position). **σ 81 km/s (rescaled to D_L(0.268), σ ∝ L^1/4):** 4 images. **σ ≤ 70 km/s:** only 3 images. The former fourth image becomes the model's counterpart of 51.3, 1.3–1.6″ from it (μ 5.8 → 2.9 as σ goes 70 → 1). No unseen image remains. 70 km/s is 14 % below the rescaled value, within the scatter of the member scaling relation | model_prediction / assumption |

**Verdict (updated the same day):** `explained: lens-model systematics, not an anomaly`.

The fourth image exists only while the z 0.268 foreground galaxy is massive enough (σ ≳ 75–80 km/s in a two-plane
model) to split 51.3. Both published models give it a cluster-member scaling-relation mass at the cluster redshift.
A two-plane model with σ ≤ 70 km/s, within the scaling relation's scatter, predicts exactly the three observed
images. The observed absence of a fourth image therefore constrains this galaxy's mass (a `hypothesis`-level
statement about the galaxy). It is not evidence for anything unusual.

The first verdict, earlier the same day, was `inconclusive: needs a multi-plane lens model`; the two-plane test
resolved it.

### Open questions

- A full multi-plane re-fit of the cluster, with the galaxy's σ free, would test whether σ ≈ 70 km/s is consistent
  with the other constraints. The two-plane test here keeps the rest of the CANUCS model fixed.
- Rule for later fields (D-042): before treating a missing extra image as a flag, check every potential that
  produces it against its spectroscopic redshift. Foreground or background galaxies modelled as members are the
  first suspect.

### Regeneration

Scripts are scratch only:
- `s51.py`: solving in both models;
- `s51pert.py`: ±30 % changes and removal of potential 8757;
- `s51iso.py`: isophote subtraction and photometry;
- `s51b.py`: high-pass sheet;
- `twoplane.py` and `twoplane_scan.py`: the two-plane lens equation with angular-diameter distances
  (β = θ − D₁ₛ/Dₛ·α₁(θ) − D₂ₛ/Dₛ·α₂(θ − D₁₂/D₂·α₁(θ))), the source traced back from 51.1–51.3, and seeds from a 0.05″
  triangle grid (under 2M cells).

Inputs:
- the CANUCS bestparam file (SOURCES, D-042);
- the CANUCS DR1 photometry catalogue;
- the CATS v4.1 maps;
- S3 byte ranges of `jw01208-o004_t002_nircam_clear-{f150w,f277w,f444w}_i2d.fits`.

Queries and downloads are dated 2026-10-08. The contact sheets (`s51_sheet.png`, `s51_iso.png`) were inspected by eye.
