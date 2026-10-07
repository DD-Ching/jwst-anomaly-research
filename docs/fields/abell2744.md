# Abell 2744 (UNCOVER)

Field record for `configs/abell2744.yaml`. Every ranking below is a `model_prediction`, and no candidate
has been vetted. Exotic-spacetime readings are hypotheses only, to be considered after the ordinary tests
in docs/methodology.md.

## Data used (observed)

- **Program 2561 (UNCOVER), observation 1:** association `jw02561-o001_t003_nircam_clear-<filter>`, MAST
  target `ABELL2744-PREIMG`. MAST was queried on 2026-10-08 for PUBLIC level-3 NIRCam imaging within 0.5′ of
  the core.
  - obs_ids: `jw02561-o001_t003_nircam_clear-` + `f115w`, `f150w`, `f200w`, `f277w`, `f356w`,
    `f410m`, `f444w`.
  - Every obs_id has a level-3 `_cat.ecsv` (4.8–9.6 MB, 52.5 MB in total). jwst 3.0.0 and photutils 3.0.0
    produced all seven, according to MAST `prvversion` and the file headers.
  - The `_i2d.fits` files are 0.83–7.7 GB and were not downloaded. Cutouts use S3 byte ranges
    (`s3://stpubdata/jwst/public/jw02561/L3/t/o001/...`, resolved with MAST path_lookup).
  - Records: `data/manifests/abell2744.ecsv` (downloads, sha256) and `data/manifests/abell2744_products.ecsv`
    (all level-3 products with S3 URIs).
- **Why this association:** its `s_region` contains the cluster core (00:14:21.2 −30:23:50). The others were
  rejected:
  - `jw02561-o003_t006` has 7 bands, but its footprint does not contain the core.
  - `jw02756-o002/o003_t001` have 6 bands but are processed with jwst 2.0.1 and are shallower.
  - 4111 (medium bands) and 3516 (F070W/F090W/F356W) lack the required wide bands.
- **DJA matched photometry:** https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix_phot.fits
  - Size: 233,527,680 bytes. sha256 `5a9ed7bc3f5f8d41486f2475f3bedb3baa2a0068b1eab6db62fdcea216633084`.
    66,172 rows. Accessed 2026-10-08.
  - The apertures are the same as in SMACS/CEERS v7.4: ASEC_1 = 0.5″ on a 0.04″/px detection image.
  - This is grizli **v7.2**, while SMACS and CEERS use v7.4.
  - It is over the 200 MB limit; the stated reason is in the config comment (`max_bytes: 240000000`).

## Lens models (recorded, not downloaded)

| Model | Reference (title as on arXiv) | arXiv / journal / DOI | Public products |
|---|---|---|---|
| UNCOVER v1.0 (2022-12-14) / v1.1 (2023-07-25), Furtak et al. | "UNCOVERing the extended strong lensing structures of Abell 2744 with the deepest JWST imaging" | 2212.04381; MNRAS 523, 4568 (2023); 10.1093/mnras/stad1627 | Dropbox folder linked from https://jwst-uncover.github.io/DR1.html#LensingMaps and https://jwst-uncover.github.io/DR2.html#LensingMaps |
| UNCOVER v2.0 (2024-08-07), same model with 13 new spectroscopic redshifts and 132 cluster members | The release asks users to cite Furtak et al. and also Price et al., "The UNCOVER Survey: First Release of Ultradeep JWST/NIRSpec PRISM spectra for ~700 galaxies from z~0.3-13 in Abell 2744" | 2408.03920; ApJ 982, 51 (2025); 10.3847/1538-4357/adaec1 | FITS maps (WCS) of deflection, κ, γ, potential and magnification at z = 1–20, on 0.04″ and 0.1″ grids, best-fit model: https://jwst-uncover.github.io/DR4.html#LensingMaps (links a Google Drive folder; not DOI-archived) |
| Bergamini et al. 2023b (GLASS-JWST) | "The GLASS-JWST Early Release Science Program. III. Strong lensing model of Abell 2744 and its infalling regions" | 2303.10210; ApJ 952, 84 (2023); 10.3847/1538-4357/acd643 | Lenstool files and the MCMC chain, but no FITS maps: https://www.fe.infn.it/astro/lensing/A2744_Bergamini23/ . The paper points to the SLOT web tool for maps; SLOT timed out on 2026-10-08 |
| Bergamini et al. 2023a (pre-JWST; marked "obsolete" on the authors' model page) | "New high-precision strong lensing modeling of Abell 2744. Preparing for JWST observations" | 2207.09416; A&A 670, A60 (2023); 10.1051/0004-6361/202244575 | Lenstool files: https://www.fe.infn.it/astro/lensing/A2744_Bergamini22/ |
| Cha et al. 2024 (MARS, free-form) | "Precision MARS Mass Reconstruction of Abell 2744: Synergizing the Largest Strong Lensing and Densest Weak Lensing Datasets from JWST" | 2308.14805; ApJ 961, 186 (2024); 10.3847/1538-4357/ad0cbf | none found. The abstract says the maps are public, but no location was found |

## Run

Run `20261007T204844Z-200b7e32`, config `abell2744_v1`. The code was `main` at f27816c plus this config, with
no `src/` changes. It took about 4 minutes.

**Sources and photometry (derived)**
- 15,648 merged sources.
- 8,236 of them match DJA one-to-one within 0.2″ (D-013).

**Gate (D-011/D-014)**
- 8,849 of 15,648 pass.
- Flags: low_snr 3,256, low_weight 2,624, no_coverage 107, sharper_than_psf 11, single_band 4,960.
- 48 of the passing sources are stars, so 8,801 galaxies are ranked.

**Stellar locus (D-015/D-016): applied**
- r50_psf is 2.13 px, from 13 catalogued stars. SMACS had 2.86 px from 22.
- It added 16 stars, 14 of them past the gate, with r50 1.95–2.48 px.
- Catalogued stars before the gate: 39 by SIMBAD or Gaia astrometry, 8 by Gaia position only.
- The calibration colour ranges (5–95%) are F150W−F444W [−1.94, +0.06] and F200W−F356W [−1.14, +0.01].
  SMACS had [−1.94, −0.52] and [−1.67, −0.23].
- Extended "stars" in the calibration set widen the red end: 4 of the 13 have r50 4.7–13.5 px. The
  main one is `f444w_3999` (r50 13.2 px, F150W−F444W +1.03).

**Screening (D-019 to D-021)**
- D-020 (spikes, no host light) removed 2 sources from the top 20:
  - #1 `f200w_5904`: s6 3.2, host ratio 0.0009;
  - #16 `f200w_3207`: s6 4.8, host ratio 0.0022.
- There were no D-019 or D-021 removals.

**Top-k composition (derived)**
- Galaxies: n = 20: known object 19 (95%), lens-related 1 (5%), catalogued star 0 (0%), cutout-flagged 0 (0%),
  spikes 0 (0%); 2 screened out before selection.
- Stars: n = 10: known object 10 (100%), lens-related 1 (10%), catalogued star 4 (40%), cutout-flagged 2 (20%),
  spikes 2 (20%).
- "Known" mostly means an UNCOVER DR1 entry in NED. It measures how heavily the field is catalogued, not how
  ordinary the sources are.

### Contact sheets (checked by eye)

**Galaxy top 20**
- No spikes, stripes, edge or no-data artifacts are visible.
- Most entries are faint compact sources: DJA mag_auto 27.4–31.1 for 13 of the 20, and 2 more have no DJA
  mag_auto. Several sit in the core's intracluster light.
- #3, #4, #8, #9 and #15 are barely visible in F200W (mag_auto 29.3–31.1; #9 has none, F200W 28.9).
  Colours of near-noise sources probably drive their scores.
- #18 `f200w_7603` has no DJA match and lies in a bright galaxy's halo.
- #12's cutout weight, 0.50, sits right at the D-021 threshold.

**Star top 10: five are extended galaxies**
- #3 `765`, #4 `7606`, #5 `7776` (an edge-on disk), #8 `8025` and #9 `8137` are extended. Their DJA r50 is
  11.9–20.2 px, against r50_psf 2.13 px.
  - #4, #5, #8 and #9 entered by `gaia_position` (a Gaia DR3 source within 0.3″, D-012). Gaia DR3 lists
    the cores of bright cluster galaxies.
- #1 `f444w_3999` and #2 `f200w_7179` are one object split into two merged sources 0.25″ apart. Both match
  UNCOVER DR1 24980. The cutouts show a bright point source about 0.4″ from the target position, next to a
  clumpy galaxy.
- #6 `7694` (no DJA match) and #10 `7893` are spike-flagged; #10 is a point source on a spiral galaxy.
- #7 `7386` shows a PSF ring pattern.

### Galaxy top 10 (after screening)

Top features are robust z per feature.

| Rank | uid | Top features | Cutout | Best cross-match (observed) |
|---|---|---|---|---|
| 2 | `f200w_8209` | F150W−F200W +9.9; ref_mag −4.0; F115W−F150W +3.6 | ok | NED UNCOVER DR1 34426 (G, 0.04″) |
| 3 | `f200w_2397` | F356W−F410M +6.6; log CI50/30 −4.8; F277W−F356W −4.5 | ok | NED UNCOVER DR1 07090 (G, 0.02″) |
| 4 | `f200w_2390` | sharpness +8.9; roundness +3.7; ref_mag +2.3 | ok | NED UNCOVER DR1 06685 (G, 0.09″) |
| 5 | `f200w_7987` | F115W−F150W +10.7; F150W−F200W +8.9; log CI50/30 −5.4 | ok | NED [PFF2025] UNCOVER 09447 (*, 0.01″) |
| 6 | `f200w_7298` | F356W−F410M −6.8; F200W−F277W +4.6; F410M−F444W +3.3 | ok | NED ABELL 2744:[FZW2023] c65.13CI (G, 0.05″) |
| 7 | `f200w_3365` | F150W−F200W +9.4; F410M−F444W +7.3; log CI50/30 −3.2 | ok | NED UNCOVER DR1 05409 (G, 0.04″) |
| 8 | `f200w_2684` | F150W−F200W +5.6; F115W−F150W −5.3; log CI50/30 −4.7 | ok | NED UNCOVER DR1 22273 (G, 0.05″) |
| 9 | `f200w_2985` | F150W−F200W +6.7; sharpness −4.5; F115W−F150W −3.7 | ok | NED UNCOVER DR1 08389 (G, 0.03″) |
| 10 | `f200w_285` | F410M−F444W −7.0; F200W−F277W −4.3; log CI50/30 −3.0 | ok | NED UNCOVER DR1 00948 (G, 0.05″) |
| 11 | `f200w_3120` | F356W−F410M +5.8; F150W−F200W −4.9; log CI70/50 −4.5 | ok | SIMBAD [MAC2016] A2744cl 1214 (G, 0.14″); lens-related (SIMBAD LeG) |

uid prefix: `jw02561-o001_t003_nircam_`.

## Notable for follow-up (unvetted)

- **`f200w_5904` (screened, originally #1).** It is point-like (DJA r50 2.08 px) with F150W−F444W +2.64,
  outside the stellar colour box. D-020 removed it because it shows spikes and no host light.
  - It matches SIMBAD UNCOVER DR2 46803 (type EmG, 0.06″) and NED [AJO2024] MSAID 45924 (0.17″).
  - So the screening rule removed a catalogued compact source with non-stellar colours from rank 1. This
    is D-020's "Revisit if" condition.
- **`f200w_7987` (#5).** It is point-like (r50 2.03 px) with F444W 22.5, F115W−F150W +1.92 and F150W−F444W
  +2.07. NED gives it type `*` ([PFF2025] UNCOVER 09447). Its colours keep it out of the star stratum.
- **`f200w_4731` (#12).** A compact (r50 1.23 px), faint source (F200W 27.25, DJA) with no SIMBAD, NED or Gaia
  counterpart within 1″. Its cutout weight is at the D-021 threshold.
- **`f277w_264` (#13).** A clumpy, multi-knot source anchored on its F277W pipeline detection. F150W−F444W is +6.0 (F150W
  29.4, F444W 23.4), and F150W−F200W is the most extreme feature in the top 20 (robust z +20). It matches NED
  ABELL 2744:[BC2023] 02.
- **`f200w_7298` (#6).** A faint, arc-like source with NED designation ABELL 2744:[FZW2023] c65.13CI (type G).
  The pipeline does not flag it as lens-related (only SIMBAD types do that). Whether it is a multiple image in a
  published lens model was not checked.
