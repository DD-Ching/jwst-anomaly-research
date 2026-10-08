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

## Lens-model consistency: Bergamini+2023b (2026-10-08, issue #41, D-030)

**Reproduced in the repository (D-030):** `python scripts/lens_consistency.py --model <name> validate` gives the
image-plane χ² and rms below. The counter-image numbers still come from issue #41's scratch code. Every number below
is `model_prediction` or `derived`, and every threshold is an ASSUMPTION.

- **Files:** `best.par` (sha256 `7245368f…`) and `obs_arcs.cat` (sha256 `d02c231f…`) from the authors' page.
  - The model has 180 dPIE potentials at z = 0.3072 and was optimised in the image plane: Chi2pos 146.60, dof 148.
  - The image list has 149 images; 28 take their z from `z_m_limit`.
- **Loading:** fixed in D-030 (letter IDs, multi-id `z_m_limit`, `_kpc` rounding).
- **Image positions:** the exact image-plane solve gives χ² = 146.64 against Lenstool's 146.60, and rms 0.427″. The
  largest residual is 22.1a, 1.62″ at σ 0.57. **No image-position anomaly.**
  Multiplicity residual: 3.2a/b, 34.1a/b and 700.1a/b each match one predicted image (pairs on the same side of
  the critical curve; D-030). `bayes.dat` samples (D-045): all three stay merged in 13 of 13 models; CATS v4.1 splits
  34.1a/b; 3.2a/b sit on the caustic in both models. Model resolution at folds, not an anomaly.
- **Counter-images:**
  - 30 predicted images are not catalogued. 15 are near-critical-curve pairs within 3″ of an observed image (|μ|
    mostly > 9). 15 are far third images, including the z = 7.39 system A200/B200/C200 (μ ≈ 10).
  - The 2561 pipeline catalogs match only 20 of 149 images within 1″, because of deblending in the core. The test
    therefore needs the 233 MB DJA catalogue, which needs a DECISIONS entry.
  - **Result: not yet testable.**

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

## Lens-model and exotic screens (2026-10-08)

Model `abell2744-bergamini23` (Bergamini+2023b). Code: `main` at 64172c4, no `src/` or `scripts/` changes.
Every number is `derived` from `observed` images and `model_prediction` positions and magnifications. Every
threshold is an ASSUMPTION (the script defaults, recorded in each JSON output).

**Inputs**
- F200W pipeline catalogue (`jw02561-o001_t003_nircam_clear-f200w_cat.ecsv`, sha256 `9953c869…`).
- DJA eazy photo-z: `https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix.eazypy.zout.fits`
  (63,570,240 B, sha256 `436626861c7dbf370419a8787dd79140fad99c2ecbd268257941afe6d493991c`, Last-Modified
  2023-12-12, accessed 2026-10-08). The `…v7.2-fix.photoz.tar.gz` tarball named in earlier notes returns S3 403
  (no such key); the zout table is published on its own, so no tarball was needed.
- Forced photometry and cutouts: S3 byte ranges of the F150W, F200W, F277W and F444W `_i2d` (no full download).

**Commands** (`JWST_ANOMALY_DATA` = shared data root; `$CAT`, `$ZOUT` as above; `$S3` =
`s3://stpubdata/jwst/public/jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear`)

```
python scripts/lens_consistency.py --model abell2744-bergamini23 validate                       # 33 s
python scripts/lens_consistency.py --model abell2744-bergamini23 images --catalog $CAT --photoz $ZOUT \
    --half-width 190 --step 0.25 --forced-image $S3-f277w_i2d.fits                             # 60 s
python scripts/exotic_screens.py --model abell2744-bergamini23 fluxratio \
    --band1 $S3-f150w_i2d.fits --band2 $S3-f444w_i2d.fits                                       # 3 min 10 s
python scripts/exotic_screens.py --model abell2744-bergamini23 radial --catalog $CAT --photoz $ZOUT \
    --max-radius 120                                                                            # 66 s
```

**Validate.** Image-plane χ² 146.64 against Lenstool's 146.60 (149/149 images solved, rms 0.427″, max 1.62″,
none over 3σ). `shared_match_images`: 3.2a/3.2b, 34.1a/34.1b and 700.1a/700.1b (each pair matches one predicted
image).

**Counter-images (`images`).** 50 systems, 176 predicted images: observed 145, candidate 8, other_source 7,
missing 1, faint 2, no_flux_ref 9, outside 4, demagnified 0. 5σ depth proxy 27.8 mag. Catalogued images with no
predicted image within 1.5″: 22.1a (the 1.62″ residual) and 64.1a (1.50″, μ ≈ −535: on the critical curve).
Forced photometry in F277W (r = 0.2″, search 1.0″, flux-ratio window [1/3, 3], ASSUMPTIONs) on 27 predicted
images: recovered 7, confused 9, absent 4, undetectable 4, ambiguous 0, no_reference 3, off_image 0.

**Flux ratios (`fluxratio`, F150W/F444W).** 149 images: underluminous 0, overluminous 0, chromatic 3, consistent
30, resolved 99, untestable 17. The 3 chromatic ones are all of system 602.1 (z from `z_m_limit`): 602.1d
(μ 24.9) has F150W/F444W ratios 4.28/0.79 and is blended with a red compact neighbour 0.5″ south (cutout), which
also shifts its siblings' reference medians (602.1a/b 0.39/1.15 and 0.37/1.09). Blend; not a flag.

**Radial arcs (`radial`, 120″).** 943 elongated sources, 110 dropped as foreground or cluster by photo-z, 149
anti-tangential, 148 not predicted radial (μ_r < 3). 35 convergence peaks (25 with no source within 1″), max
5 lines. 200 random-angle draws give 50.4 peaks on average (95th percentile 65) and a maximum of ≥5 lines in
115 of 200 draws: p_random ≥ 0.575 for every peak. **No radial peak has p_random < 0.05.**
These numbers predate the D-034 spike veto. With it (current code, re-run for D-049 with the same catalogue,
photo-z and 120″ radius), 42 spike segments are dropped, 134 `anti` arcs remain, the maximum is 5 lines and
p_random 0.505: still null (docs/exotic_limits.md, base-screen table).

**Checks applied to each flag.** 3-band cutouts (F150W/F200W/F444W) of every flagged position, viewed by eye.
Also high-pass versions (0.6″ median filter subtracted) to remove cluster-member light. Further checks:
- 0.2″ aperture S/N at the exact predicted position on the high-pass F150W and F277W stamps;
- DJA neighbours within 2.5″ whose photo-z 95% interval contains z_sys;
- each prediction recomputed from each observed image's own back-traced source position, and for 700.1 at
  z = 1, 1.5, 2, 3, 5.

| Flag (system, μ) | RA, Dec (deg) | Screen | Forced S/N pred → best (at pred) | Flux ratio | Ordinary explanation tested | Verdict |
|---|---|---|---|---|---|---|
| 3.2c (3.43) | 3.57667, −30.40165 | confused | 40.7 → 131.4 (1.4) | 4.26 | Knot mismatch: catalogued 3.1c (DJA 17477, z_spec 3.987) is 0.57″ away; the c image is compact, so all knots fall in one aperture | Ordinary (image seen as 3.1c) |
| 4.1, 4.2 central (−0.88, −0.74) | 3.57950/3.57955, −30.40932 | confused ×2 | 8.3/6.4 → ~2190 | 1471/1807 | Cluster-member light: 0.7″ from the BCG core (6983, DJA z_spec 0.303) | Ordinary (buried; untestable) |
| **4.2c (8.71)** | **3.57972, −30.40835** | **absent** | **251.7 → 30.7 (−0.3); high-pass 0.2** | **0.14** | Catalogue: no DJA or pipeline source with 4.2's colour within 2″ (6985 at 1.14″ is diffuse halo light, absent after high-pass). Knot mismatch: the forced reference is the galaxy peak; knot-only expectation is still high-pass S/N ≈ 60 (4.2a/b give 63/68). Cluster light: removed by high-pass, nothing. Source position: robust to 0.65″. z: spectroscopic (not z_m_limit). Not near a critical curve (μ 8.7). Achromatic absence (F150W–F444W all ≈ 0σ). The c image of the companion knot (4.1c) is seen 2.5″ away. **Not tested:** model error from the galaxy-scale potential of BCG 6983, 2.9″ away; the MCMC (`bayes.dat`) spread | **Survives the cheap tests; most likely a model-prediction residual near a bright member. Candidate for `/vet-candidate`, not an anomaly** |
| 8.1c (3.33) | 3.57650, −30.40231 | confused | 11.3 → 31.6 (−0.9) | 3.28 | Search radius: the 1″ search found a bright z = 0.25 galaxy (2122). Pipeline sources 2121 and 2118, 0.80″ and 1.03″ away, have F150W−F444W −1.46, as do 8.1a/b (−1.41/−1.03) | Ordinary (counterpart at 0.8–1.0″) |
| 18.1 extra (−3.87) | 3.58861, −30.39626 | confused | 122.0 → 1941 (226 raw, 5 high-pass) | 103.9 | Cluster-member light: 0.8″ from a member core (DJA z_spec 0.302). Robustness: the image exists only for the mean and 18.1a source positions; with 18.1b or 18.1c's position the model predicts 3 images | Ordinary (galaxy-scale caustic, not robust) |
| 22.1a (5.41) | 3.58739, −30.41163 | confused; unpredicted | 30.7 → 113.1 (0.6) | 4.50 | Position residual: observed 22.1a (DJA 13659, z_spec 5.283) is 1.62–1.66″ away, the known largest residual. The 1″ search found a z = 1.1 neighbour | Ordinary |
| 33.1 extra ×2 (−32.2, −19.1) | 3.58498, −30.40352; 3.58423, −30.40280 | confused ×2 | 52.5/22.1 → 1327/1191 | 110/155 | Cluster-member light: on or next to members (DJA z_spec 0.304/0.316). Near-critical magnification. The reference 33.1a has only high-pass S/N 6–12 | Ordinary (untestable) |
| 33.1 third (2.94) | 3.60056, −30.39528 | confused | 5.4 → 152 (6.6 raw, 0.2 high-pass) | 43.7 | Predicted S/N 5.4 is marginal; crowded field next to a 20.1 mag member | Ordinary (below sensitivity) |
| 34.1 (+44.98) and pair 34.1a/b | 3.59269, −30.41106; a 3.59341, −30.41081; b 3.59380, −30.41069 | absent; shared_match | 387 → 15.5 (12.1; a compact z ≈ 0.34 dwarf) | 0.04 | Near critical curve: the observed fold pair a/b (1.3″ apart, high-pass S/N 10.6/8.1) straddles a critical curve that the model places so its +parity image falls 2.4″ away. The μ = 45 scaling of 34.1c's flux is invalid there. Robust to 0.3–1.0″ | Ordinary (the observed 34.1b is the predicted image; misplaced critical curve) |
| 700.1 ×2 (−91.95, +10.30) and pair 700.1a/b | 3.57519, −30.35561; 3.57883, −30.35360; a 3.57970, −30.35772; b 3.57917, −30.35783 | absent ×2; shared_match | 3634/392 → 2.9/8.4 | 0.00/0.02 | Wrong z (z_m_limit 1.217): at the sampled z = 1, 1.5, 2, 3, 5 the model never gives 700.1a and b separate images, so the configuration is not reproduced. The 4th image comes and goes with the source position (2 images from 700.1a's position). The −92 image is the fold partner of 700.1c, 1.9″ away. At 160″ from the core this is the model's infall region | Ordinary (model-limited system) |
| 3.2a/3.2b pair | 3.58921, −30.39382; 3.58896, −30.39380 | shared_match | — | — | Knot identification along the giant merging arc (μ ≈ ±50, cutout): knots straddle the critical curve | Ordinary |
| 64.1a | 3.58119, −30.39871 | unpredicted | — | — | On the critical curve (μ ≈ −535); nearest predicted image 1.50″ away | Ordinary |
| 28d (−4.23) | 3.58721, −30.40146 | catalog `missing` | 33.6 → 21.5 | 0.68 | Forced photometry recovers a faint compact source 0.85″ away (cutout) | Not a flag. Possible uncatalogued 28d (unconfirmed) |

**Summary.**
- Predicted images screened: 176 (27 with forced photometry). Images in the flux-ratio test: 149. Radial
  peaks: 35.
- Flags raised: 16 = 13 forced (absent 4, confused 9) + 3 shared-match pairs. Also inspected: 2 unpredicted
  images, 3 chromatic images and 1 catalogue-missing image.
- Surviving the cheap ordinary tests: one, 4.2c, which is most likely a model residual near BCG 6983. No
  under- or overluminous image, and no significant radial convergence.
- Nothing here is evidence of non-standard lensing.
- **Next:**
  - ~~`bayes.dat` spreads for 4.2c, 34.1 and 700.1~~ — done (D-045): 4.2c moves 0.2–0.4″ (μ 8.7–10.0); the pairs stay
    merged.
  - ~~Run `/vet-candidate` on 4.2c~~ — done: explained, see docs/candidates/abell2744-family4-c.md.
  - Widen the forced search to the model's positional rms × 3 (≈1.3″) together with a colour match (the
    8.1c lesson).

### Re-run under the D-034 rules and the 4.2c verdict (coordinator, 2026-10-08)

- `images --forced-image` was re-run with the code after #49 (D-034: compact, at-position, consistent references;
  residual-scaled search radius). Classes: recovered 3, confused 3, absent 2, undetectable 3, no_reference 16.
  - The only `absent` images are both of system 700.1 (z = 1.217 fitted). The worker had already shown that the
    model cannot reproduce this system at any z from 1 to 5.
  - **4.2c is now `no_reference`.** System 4.2's catalogued images (4.2a, 4.2b) fail the reference rules: they are
    resolved or neighbour-contaminated. The first pass's 252σ prediction therefore rested on aperture fluxes that
    do not scale with |μ|.
- **4.2c verdict: ordinary, not a candidate.**
  - The reference knots 4.1a and 4.2a are small, faint clumps of one thin arc (high-pass cutouts, F150W/F277W/
    F444W).
  - The predicted position lies on the brightest cluster galaxy's halo, where the high-pass residuals are strong.
  - The family's other knot image, 4.1c, is observed 2.5″ away. That is consistent with knot-level model offsets
    near the cluster core.
- 34.1 (μ 45) and 28 are recovered under the new rules (flux ratios 0.80 and 0.68).
- **Family 4 c images, vetted** (docs/candidates/abell2744-family4-c.md): 4.1c is underluminous 4–8× after BCG
  subtraction, and 4.2c is undetected. Both are explained by the model's μ and position systematics next to member
  34423: plausible changes move μ(4.1c) from 3.9 to 28.7, and CATS v4.1 gives 7.3.
- **Result: the Abell 2744 screens are a null result.** Flags raised: 16. Surviving vetting: 0.
