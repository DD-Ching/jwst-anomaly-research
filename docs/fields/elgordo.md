# El Gordo (ACT-CL J0102−4915)

Field note for the `batch-clusters` expansion (D-023). El Gordo is a merging cluster at z = 0.87. The run
sections below are pipeline output (`derived` / `model_prediction`) and are unvetted; no physical
interpretation is made.

## Data (accessed 2026-10-08)

- **JWST:** PEARLS GTO program 1176 (PI Windhorst), NIRCam observation 241, association `t012`, all 8 PUBLIC
  level-3 bands: `jw01176-o241_t012_nircam_clear-{f090w,f115w,f150w,f200w,f277w,f356w,f410m,f444w}`.
  - The cluster's SIMBAD position (ACT-CL J0102-4915) lies in NIRCam module B. Module A images a flanking
    field about 3′ north. Footprints come from MAST `s_region`.
  - Every band has a `_cat.ecsv` (2.0–3.4 MB as served) and an `_i2d.fits` (0.31–1.29 GB). The `_i2d` files
    are never downloaded; cutouts use S3 byte ranges.
  - **Pipeline version:** jwst 3.0.0, CRDS `jwst_1584.pmap`. The catalogs and `_i2d` headers were reprocessed
    on 2026-10-01.
    - MAST's product listing still reports `prvversion` 2.0.1 and the old sizes, so `acquire` accepts the
      served size only after an HTTP HEAD check (PR #29).
    - For comparison, the SMACS 2736 catalogs are jwst 2.0.1 and CEERS is jwst 3.0.0 (D-010, D-022).
  - Manifests: `data/manifests/elgordo.ecsv` (8 catalogs with sha256) and `data/manifests/elgordo_products.ecsv`
    (all 16 level-3 products with S3 URIs).
- **Matched photometry (D-013):** DJA
  [`elgordo-grizli-v7.0-fix_phot.fits`](https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix_phot.fits).
  - 17,170,560 bytes, sha256 `ac1cf7064edc3f1a907c1b95b7eacd864e19810ff559e7dc39b917081ed66e31`, 6,775 rows.
  - Aperture 1 is 0.5″ in diameter.
  - v7.0 is the only El Gordo catalog on the DJA v7 S3 path (v7.1–v7.5 are absent). It is therefore older than
    the v7.4 catalogs used for SMACS and CEERS.

## Lens models (recorded, not downloaded)

Every title was checked on its arXiv abstract page and every URL was opened on 2026-10-08.

| Model | Method | Data | Public download |
|---|---|---|---|
| Caminha et al. 2023, "A MUSE view of the massive merging galaxy cluster ACT-CL J0102-4915 (El Gordo) at z = 0.87: robust strong lensing model and data release", A&A 678, A3, [arXiv:2209.02718](https://arxiv.org/abs/2209.02718) | Lenstool (parametric) | HST + MUSE | [CDS J/A+A/678/A3](https://cdsarc.cds.unistra.fr/ftp/J/A+A/678/A3/): best-fit magnification maps for z = 1–20 (9 FITS files); `files/` holds `best_fit.par`, `bayes.dat`, the cluster members and the multiple-image list `obs_arcs_v1_new_IDs.dat`; MUSE redshift catalogue (402 rows). No κ map. |
| Diego et al. 2023, "JWST's PEARLS: a new lens model for ACT-CL J0102−4915, 'El Gordo', and the first red supergiant star at cosmological distances discovered by JWST", A&A 672, A3, [arXiv:2210.06514](https://arxiv.org/abs/2210.06514) | WSLAP+ (free-form) | JWST PEARLS 1176 | none found |
| Frye et al. 2023, "The JWST PEARLS View of the El Gordo Galaxy Cluster and of the Structure It Magnifies", ApJ 952, 81, [arXiv:2303.03556](https://arxiv.org/abs/2303.03556) | light-traces-mass | JWST PEARLS 1176 | none found (imaging DOI 10.17909/x49n-d207 only) |
| Galan et al. 2024, "El Gordo needs El Anzuelo: Probing the structure of cluster members with multi-band extended arcs in JWST data", A&A 689, A304, [arXiv:2402.18636](https://arxiv.org/abs/2402.18636) | Herculens + NIFTy, galaxy-scale, built on the Caminha+2023 cluster model | JWST PEARLS 1176 | [GitHub `aymgal/ElAnzuelo_modeling_public`](https://github.com/aymgal/ElAnzuelo_modeling_public) (COOLEST format) |
| Cerny et al. 2018, "RELICS: Strong Lens Models for Five Galaxy Clusters From the Reionization Lensing Cluster Survey", ApJ 859, 159, [arXiv:1710.09329](https://arxiv.org/abs/1710.09329) | Lenstool (pre-JWST; the paper calls this model under-constrained) | HST RELICS | [MAST HLSP RELICS models](https://archive.stsci.edu/hlsps/relics/act0102m49/models/) (DOI 10.17909/T9SP45): `lenstool/` v1 and `glafic/` v2 and v3, each with κ, γ, ψ, deflection and z = 1, 2, 4, 9 magnification maps |

**For the D-023 lens-model test:**
- **Caminha+2023** is the public model with MUSE spectroscopic constraints. Its multiple-image list gives the
  positions to compare against. κ/γ would have to be regenerated from `best_fit.par` with Lenstool (untested).
- **RELICS** offers the only ready-made κ/γ/deflection grids.
- **Diego+2023 and Frye+2023**, the JWST-era free-form models, have no public files.

## Lens-model consistency: Caminha+2023 (2026-10-08, D-030)

Reproduced in the repository with `python scripts/lens_consistency.py --model elgordo-caminha23 validate` (D-030).
That gives χ² 82.53 against Lenstool's 80.22, rms 0.754″, and the magnification maps below. The counter-image
numbers further down still come from the scratch run of issue #41. Every number below is `model_prediction` or
`derived`, and every threshold is an ASSUMPTION.

- **Files:** CDS `files/best_fit.par` (sha256 `7b0153ae…`) and `obs_arcs_v1_new_IDs.dat` (sha256 `d6317439…`).
  - The model has 265 dPIE potentials, single plane at z = 0.8703. Lenstool optimised it in the image plane:
    Chi2pos 80.22, dof 52.
  - The image list has 56 images in 23 systems, all with spectroscopic z, and errors of 0.621″.
- **Magnification maps:** the port reproduces CDS `magnification_best_fit_z{2,8}.fits`. On about 20k pixels with
  |μ| < 10, the median |Δμ|/μ is 3.8e-5 / 6.7e-5, and the 99th percentile is 2.6e-3 / 5.1e-3.
- **Image positions:**
  - The exact image-plane solve gives χ² = 82.5 against Lenstool's 80.22, and rms 0.754″ against the paper's 0.75″.
  - The source-plane back-trace approximation gives 121.6 and is not valid for image-plane-optimised models.
  - The largest residuals are 23c 1.54″, 3a 1.49″ and 17a 1.44″. **No image-position anomaly.**
- **Counter-images:**
  - The model predicts 72 images: 56 match the catalogued images and 16 are not catalogued. 14 are third images
    30–60″ away with |μ| 2.4–4.6; the other two are 4th images on bright members (sys 6 and 12).
  - Most observed images of these MUSE Lyα systems are at or below the F277W noise in a 0.5″ aperture. The 5″
    search circles hold 9–32 DJA sources each, so a flux and colour match cannot discriminate.
  - **Result: inconclusive; no candidate.**
- **Frame:** DJA v7.0 sits at dRA +0.224″, dDec −0.016″ (median) relative to the RELICS/HST frame of the image list.

## Run

- Run `20261007T210354Z-9dccdf35`, config [`configs/elgordo.yaml`](../../configs/elgordo.yaml) (`elgordo_v1`,
  sha256 prefix `d647e092f21b872c`). The config's `stages:` block is identical to `configs/reference_sample.yaml`.
- Code: `b787e29`, which includes the PR #29 fix. The earlier run `20261007T205438Z-b58ac10d` differed only in a
  config comment, and its sample results are identical.
- **Sources:** 6,921 merged sources. DJA matches one-to-one for 3,081 of them (D-013).
- **Gate (D-011/D-014):** 3,397 of 6,921 sources passed.
  - Flags: edge 94, low_snr 2,132, low_weight 1,347, no_coverage 31, sharper_than_psf 8, single_band 2,968.
  - Of the passing sources, 11 are catalogued stars (D-012) and 3,386 galaxies are ranked.
- **Stellar locus (D-015/016):** not applied. Only 4 catalogued stars have 20 < mag < 22.5, against the 10
  needed. Spike screening therefore uses only the host test (D-020).
- **Screening:** 4 sources were removed from the galaxy top 20:
  - #5 `f200w_2530` and #7 `f200w_2969`: spikes, no host light (D-020);
  - #6 `f200w_2724` and #22 `f277w_47`: low cutout weight (D-021).
- **Top-k composition, galaxies:** n = 20: known object 6 (30%), lens-related 3 (15%), catalogued star 0 (0%),
  cutout-flagged 0 (0%), spikes 0 (0%); 4 screened out before selection (D-019 to D-021).
- **Top-k composition, stars:** n = 10: known object 10 (100%), lens-related 0 (0%), catalogued star 6 (60%),
  cutout-flagged 10 (100%), spikes 10 (100%).

**Contact sheets (checked visually).**
- Star stratum: all 10 are isolated point sources with hexagonal spikes. #3 `f200w_3049` sits on a no-data
  edge strip.
- Galaxy stratum:
  - No spikes remain.
  - 12 of the top 20 are faint compact sources in noisy, speckled backgrounds: #1, #2, #3, #8, #9, #11,
    #12, #14, #18, #19, #21, #23.
  - #3 `f200w_1601` is barely visible in F200W despite its −11.4σ colour, so a DJA mismatch or catalog effect
    is possible.
  - #4 `f200w_2889` is a knot on the edge of a large galaxy, likely a deblending fragment.
  - #23 `f200w_889` has a checkerboard pixel cluster about 1″ away.
  - #24 `f200w_3070` is a bright extended galaxy whose top feature is `n_gaps` = +6. This is probably the
    0.1″ merge radius splitting a large galaxy across bands (D-003).

### Top 10 (galaxy stratum, after screening; original ranks)

UIDs are short for `jw01176-o241_t012_nircam_<uid>`. Every cutout flag is `ok`.

| Rank | uid | Top features (robust z) | Cross-match best match (≤1″) |
|---|---|---|---|
| 1 | `f200w_1779` | F277W−F356W −7.7; sharpness +4.2; n_gaps +4.0 | NED ACT-CL J0102-49:[CGR2023] 04b (G_Lens, 0.24″) |
| 2 | `f200w_333` | sharpness +6.6; F200W−F277W +4.9; F090W−F115W +2.6 | none |
| 3 | `f200w_1601` | F200W−F277W −11.4; F115W−F150W −3.8; log CI50/30 +3.0 | none |
| 4 | `f200w_2889` | sharpness +11.2; roundness −2.9; red_dropout +2.7 | none |
| 8 | `f200w_2695` | log CI50/30 −4.3; n_gaps +4.0; log CI70/50 −3.6 | none |
| 9 | `f200w_793` | F277W−F356W −6.9; log CI50/30 −3.1; log CI70/50 −3.0 | none |
| 10 | `f200w_1065` | F200W−F277W +5.9; F277W−F356W +5.8; F115W−F150W +5.6 | none |
| 11 | `f200w_545` | F200W−F277W −6.6; log CI50/30 −3.7; log CI70/50 −3.5 | none |
| 12 | `f200w_1709` | F150W−F200W +5.7; F277W−F356W −5.6; F200W−F277W −5.2 | none |
| 13 | `f444w_1836` | F200W−F277W +9.3; F277W−F356W +7.0; F150W−F200W +6.0 | NED ACT-CL J0102-4915:[ALB2018] 01 (G, 0.08″; 5 matches incl. NED G_Lens, SIMBAD LeG `[ALB2018] ACTJ0102-1` = `[FKO2024] AC0102-C224`) |

The other lens-related source in the top 20 is #18 `f277w_1143`: NED [CGR2023] 19c (G_Lens, 0.10″). Its top
features are F277W−F356W −7.3 and F356W−F410M +7.2.

### Notable for follow-up (unvetted)

- **#13 `f444w_1836`:** a catalogued lensed galaxy. It appears as a thin diagonal arc in F200W, with three colour
  outliers above 6σ. It is a direct target for the Caminha+2023 model comparison (D-023).
- **#1 `f200w_1779`:** 0.24″ from Caminha+2023 multiple image 04b. It has a blue F277W−F356W outlier and detection gaps in 4 of 8 bands. The next step is to check its counter-images against the
  model's image list.
- **#18 `f277w_1143`:** Caminha+2023 image 19c. Its outliers are opposite-sign colours on either side of F356W
  (F277W−F356W −7.3, F356W−F410M +7.2).
- **#17 `f200w_2828`:** a horizontal linear feature crossing the cutout, 10–14″ west of #13 and #8 in module B.
  It has no counterpart within 1″ and an F200W−F277W outlier of −8.3. It could be an arc or an artifact (stripe/spike);
  this is untested.
- **#10 `f200w_1065`:** a clumpy, multi-knot source with three colour outliers of +5.6 to +5.9 and no
  counterpart within 1″.
