# Halo-subtracted flux ratios on Abell 370 (CANUCS model) — plan (frozen before running)

- **Hypothesis:** a fifth cluster extends the D-032 flux-ratio line. No compact catalogued image deviates by
  > 0.75 mag (both routes, `OUTLIER_MAG`) from the `abell370-canucs` μ ratio, and `fluxratios --halo-phot`
  finds no 3σ `flux_tension`.
- **Ordinary explanations for any flag:** halo/ICL residual (blank-aperture checks), neighbour light, image
  misidentification in the catalogue, microlensing / substructure (known populations), μ near a critical curve.
- **Inputs:** `abell370-canucs` model; CANUCS 1208 o002 t001 F150W + F444W `_i2d` by S3 byte range (no download);
  unchanged `scripts/halo_phot.py` gates (q ≥ 0.80 star calibration, |μ| ≤ 20), ±1 mag injections.
- **Decision rule:** a `flux_outlier` or `flux_tension` goes to `/vet-candidate`; none → null for a fifth cluster
  with the injection reach stated. Exploratory: the gates were set on Abell 2744, El Gordo and SMACS, not here.
