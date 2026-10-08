# WHL0137−08 ("Sunrise arc", Earendel)

Lensing cluster WHL J013719.8−082841 (z = 0.566). It lenses the z ≈ 6 Sunrise arc and the extremely
magnified source Earendel. Config `configs/sunrise.yaml` (sample `sunrise_nircam`), whose archive, cloud and
`stages:` blocks are the same as `configs/reference_sample.yaml`. All rankings below are unvetted
`model_prediction`s.

## Data used (observed)

- **Ranked:** program 2282 (PI Coe), observation 10, t001, NIRCam, 2022-07-30, PUBLIC. It has 8 bands:
  F090W, F115W, F150W, F200W, F277W, F356W, F410M and F444W (`jw02282-o010_t001_nircam_clear-<band>`).
  - Level-3 `_cat.ecsv` exists for every band (1.6–2.8 MB, jwst 2.0.1, the same as SMACS 2736). They are pinned in
    `data/manifests/sunrise.ecsv`, and all level-3 products with S3 URIs are in `sunrise_products.ecsv`.
  - The `_i2d.fits` files (SW 1.29 GB, LW 0.31 GB) are read only by S3 byte range.
- **All public level-3 NIRCam imaging within 6′** (MAST, queried 2026-10-07 UTC, all data rights):

  | Epoch | obs prefix | Date (t_min) | Bands | jwst | `_cat.ecsv` |
  |---|---|---|---|---|---|
  | 1 | `jw02282-o010_t001_nircam_clear-` | 2022-07-30 | F090W F115W F150W F200W F277W F356W F410M F444W | 2.0.1 | yes |
  | 2 | `jw02282-o120_t001_nircam_clear-` | 2023-01-10 | F090W F115W F277W F356W | 2.0.1 | yes |
  | 3 | `jw06882-o052_t052_nircam_clear-` (VENUS) | 2025-07-08 | F150W F210M F300M F444W | 3.0.0 | yes |

  Filters shared with epoch 1, for `scripts/epoch_compare.py` (D-017):
  - epoch 2 has F090W, F115W, F277W and F356W (a 5-month baseline, same pipeline version);
  - epoch 3 has F150W and F444W (a 3-year baseline), with a jwst 3.0.0 vs 2.0.1 catalog difference (D-010).
  - No exclusive-access NIRCam imaging was listed.
- **DJA matched photometry (D-013):** https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix_phot.fits
  - Version v7.5 (SMACS and CEERS use v7.4). 28,779,840 bytes, sha256
    `8bd178e94156d6a9a126e487dec8a850b6df9f6d76ceb9e36931d9960558ebf0`, downloaded 2026-10-07 UTC.
  - It has 10,938 rows and 19 HST+JWST bands, with ASEC_1 = 0.5″.
  - It has F210M and F300M columns, filters that MAST lists here only for 6882. So the DJA stack probably includes
    the 2025 epoch (inferred).

## Lens models (model_prediction; not downloaded)

- **RELICS HST models**, v1 (HLSP DOI [10.17909/T9SP45](https://doi.org/10.17909/T9SP45)):
  `glafic`, `lenstool`, `wslap` and `zitrin-ltm-gauss`.
  - Each has κ, γ, deflection and z = 1/2/4 magnification maps. Glafic, WSLAP and LTM also have a z = 6.2
    magnification map (`z06p2-magnif-rad`).
  - Download: https://archive.stsci.edu/hlsps/relics/whl0137m08/models/
  - These are the four modelling codes used in Welch et al. 2022 (Nature).
- **Scofield, Jee, Cha & Park 2025**, "Is Earendel a Star?: Investigating the Sunrise Arc Using JWST Strong and Weak
  Gravitational Lensing Analyses", arXiv:[2504.08879](https://arxiv.org/abs/2504.08879).
  - A JWST-era MrMARTIAN strong+weak model built on 2282 epochs 1–2.
  - Public κ, deflection and magnification maps scaled to z = 6.2: Zenodo
    [10.5281/zenodo.15110933](https://doi.org/10.5281/zenodo.15110933) (CC-BY-4.0).
- **Welch et al. 2022**, "JWST Imaging of Earendel, the Extremely Magnified Star at Redshift z = 6.2", ApJL 940, L1,
  arXiv:[2208.09007](https://arxiv.org/abs/2208.09007).
  - It uses 5 earlier models (LTM, Glafic c = 1 and c = 7, WSLAP+, Lenstool).
  - The paper points to https://cosmic-spring.github.io for products. The page fetched on 2026-10-07 (UTC) linked no
    lens-model files.

## Earendel (reference case for time-domain lensing)

- **Position:**
  - Scofield+2025 Table 3, image 1.1e: RA 24.346822, Dec −8.464511 (JWST frame). This is the reference used here.
  - Welch+2022 ApJL: 01:37:23.25 −08:27:52.27, 0.19″ from the reference.
  - Welch+2022 Nature, arXiv:[2209.14866](https://arxiv.org/abs/2209.14866), HST frame: 01:37:23.232 −08:27:52.20,
    0.09″ from the reference.
  - SIMBAD's `[WCD2022] WHL0137-LS` position is corrupt (quality D, Dec +23.2), so do not use it.
- **Not ranked:** there is no pipeline merged source within 0.84″ of any of the three positions.
  - The nearest sources are `f277w_1631` at 1.0″ (gated out: low_snr, single_band), `f200w_2083` at 1.1″ (rank 98)
    and `f200w_2054` at 1.2″ (rank 2578).
  - DJA has object 4244 at 0.14″, at 27.2 AB in F200W in its 0.5″ aperture. The per-band pipeline catalogs do not
    list it separately: it is part of an arc segment or below their threshold (an F200W/F150W cutout was checked).
  - Its nature (single star, binary or star cluster) is debated: Pascale et al. 2025, "Is Earendel a Star Cluster?",
    arXiv:[2507.05483](https://arxiv.org/abs/2507.05483).

## Run `20261007T204826Z-7f9b687b`

The code is at `f27816c`; the only uncommitted file was `configs/sunrise.yaml`.

- **Merged sources:** 5,943. DJA matched 2,993 of them one-to-one (0.2″).
- **Gate:** 2,923 of 5,943 passed. Flags: low_snr 2173, single_band 2425, low_weight 1195, edge 36, no_coverage 17,
  sharper_than_psf 9. That leaves 2,918 galaxies ranked and 5 stars.
- **Stellar locus:** not applied. Only 2 catalogued stars have 20 < mag < 22.5, and 10 are needed. So D-019
  stellar-colour screening is off, and only the D-020 host test screens.
- **Screened from the galaxy top 20:** 3, all spikes with no host light (D-020): `1569` (#5), `1691` (#7) and
  `2092` (#19). Their host ratios are below 0.004, and the D-021 low-weight screen removed nothing.
- **Top-k composition (galaxies):** n = 20: known object 2 (10%), lens-related 1 (5%), catalogued star 0 (0%),
  cutout-flagged 1 (5%), spikes 0 (0%); 3 screened. The one lens-related match is probably spurious (see `1743`
  below), so the true lens-related count is likely 0.
- **Stars stratum (n = 5):** all 5 are Gaia DR3 stars with spikes (`1977` has a NaN core). The contact sheet
  shows clean PSFs.

### Galaxy top 10 (unvetted)

uids are `jw02282-o010_t001_nircam_<short>`. Cutout flags are F200W.

| Rank | Short uid | Top features (robust z) | Cutout | Cross-match best (≤ 1″) |
|---|---|---|---|---|
| 1 | `f200w_549` | ref_sharpness −13.7; F200W−F277W −6.6; nn_dist −3.2 | ok | none |
| 2 | `f200w_1578` | CI_50_30 +6.5; F200W−F277W +6.3; CI_70_50 +6.1 | ok | none |
| 3 | `f200w_2205` | F410M−F444W +7.3; F200W−F277W −4.8; CI_70_50 −4.3 | ok | none |
| 4 | `f200w_1869` | sharpness +6.9; n_gaps +6.3; blue_dropout +2.9 | ok | none |
| 6 | `f200w_1743` | F356W−F410M +8.6; F410M−F444W −6.8; CI_70_50 −3.0 | ok | SIMBAD `NAME Sunrise Arc` (G, 0.31″) |
| 8 | `f200w_1432` | F410M−F444W −6.9; F200W−F277W −5.0; CI_50_30 −4.7 | ok | none |
| 9 | `f200w_697` | F410M−F444W −8.0; CI_70_50 −4.1; F200W−F277W −3.9 | ok | none |
| 10 | `f200w_1419` | CI_50_30 −5.4; F200W−F277W +4.5; CI_70_50 −2.9 | ok | none |
| 11 | `f200w_1032` | n_gaps +6.3; F150W−F200W −4.8; F115W−F150W +3.8 | edge, nan | none |
| 12 | `f200w_1907` | F150W−F200W +6.6; F200W−F277W +5.6; F115W−F150W +4.9 | ok | none |

Ranks 5 and 7 were screened (D-020).

### Contact sheets, checked by eye

- **Contaminants in the galaxy top 20:**
  - `1869` (#4) is a 29.3 mag F200W+F444W detection with no DJA match. It sits on the diffraction-spike streaks of
    Gaia star `1817`, 2.6″ away. This is a multi-band spike detection, the D-014 "Revisit if" trigger.
  - `1032` (#11) lies on the mosaic edge.
- **Near-noise sources:** `549` (#1) and `677` (#16) are F200W-only pipeline detections, at aper50 ≈ 29.3. A DJA
  match (dja05 F200W 27.8 and 27.3) confirms them under D-014. The cutouts show little at the target.
- **Half-depth plateau:** 1,909 of the 5,943 sources have a coarse-map rel_weight of 0.45–0.55 (842 below 0.5,
  1,067 above). 8 of the 40 pool cutouts have core weights of 0.503–0.518. The 0.5 threshold of D-011 and D-021
  therefore splits a half-depth plateau of this mosaic on noise.
- **Contaminants absent:** no stars or star+galaxy blends remain after screening.

### Notable for follow-up (unvetted)

- **`f200w_1756` (#13):** compact and bright (F200W aper50 24.7), detected in F115W, F150W, F200W and F444W, with
  no DJA object within 0.2″. A bright straight stripe passes about 1″ from it in the cutout. It is probably a spike
  of Gaia star `1817`, 14.9″ away. Compare it with epoch 2 (F115W) and the DJA segmentation.
- **`f200w_1907` (#12):** an extended two-knot, arc-like morphology with extreme red F115W→F277W colours (+4.9 to
  +6.6σ). It has no catalogue counterpart.
- **`f200w_1578` (#2):** a bright extended galaxy crossed by a dark lane-like band, with extreme concentration and a
  red F200W−F277W colour. It is uncatalogued within 1″.
- **`f200w_1743` (#6):** compact and red (dja05 F200W 26.2, F444W 25.2). Its SIMBAD `NAME Sunrise Arc` match
  comes from coordinates in 2020ApJ...889..189S. Those coordinates are 37.7″ from the arc's images (`[WCD2022]
  1.1a–e`) and from Earendel, so the "lens-related" label is probably a SIMBAD position error.
- **`f200w_1082` (#23):** this is `[BCB2023] WHL0137-08004`, a published z ≈ 9–10 candidate (Bradley et al. 2023,
  ApJ 955, 13, arXiv:[2210.01777](https://arxiv.org/abs/2210.01777)). The ranking recovers a known rare
  population here.

## Two-epoch search, o010 (2022-07-30) against o120 (2023-01-10) (D-027)

The code is at this PR's head; outputs are under `outputs/transients_sunrise/` (not committed). Every number is
`derived`, and every threshold is an ASSUMPTION.

- **Catalog stage** (`transient_search.py`, jwst 2.0.1 in both epochs):
  - 660–873 matched pairs per band (the epoch-2 footprint overlaps only part of epoch 1);
  - frame shifts of ≤ 0.018″;
  - candidates: F090W 71, F115W 65, F277W 129, F356W 107.
- **Combination** (`transient_combine.py`): 58 positions have the same kind in ≥ 2 bands (13 variable, 27
  appeared, 18 disappeared). 6 lie near Gaia DR3 sources (59 in a 4′ cone) and are excluded, leaving 52 candidates. Earendel is measured alongside as a
  reference position and is not counted among them.
- **Forced photometry,** 0.15″ aperture, recentred with a 0.15″ box: `appeared` candidates are centroided in
  epoch 2, all others in epoch 1. The median Δm of the candidates is ≤ 0.005 in every band, and the 5–95% range is −0.39…+0.16 mag.
  **0 of 52 candidates pass in all four bands.** Passing in one or two bands:
  - `c0049` (disappeared in F277W and F356W; it passes in F090W and F115W, −0.54/−0.66 mag): a linear streak crosses
    the position in epoch 2 in all four bands;
  - `c0018` (appeared in F277W and F356W; it passes in F356W only, +0.32 mag): galaxy outskirts, 0.3″ from a bright
    core, at the epoch-2 LW footprint edge.
- **Result: no credible transient** at |Δm| ≥ 0.3 mag and ≥ 5σ in all four bands, for sources the catalogs
  detect at S/N ≥ 10 in one epoch, in the overlap of the two epochs.
- **Noise calibration:** forced photometry of 150 ordinary sources (F356W 25.5–28 mag) gives a robust std of the
  flux significance of 1.18–1.49. So ERR-based significances are 1.2–1.5× too large.
- **Earendel:**
  - At the Scofield+2025 reference position, without recentring, it seemed to brighten by −0.62 (F277W) and −0.76
    mag (F356W, 5.7σ). This is an artefact: the centroid lies 0.15–0.16″ away, and the PSF rotates by about 180°
    between the epochs.
  - Recentred (centroid 0.15–0.16″ from the reference): Δm = +0.11 ± 0.18 (F090W), +0.18 ± 0.12 (F115W),
    −0.05 ± 0.08 (F277W) and −0.08 ± 0.07 (F356W), all within 1.5σ. The errors are ERR-based plus a 0.05 mag floor,
    before the 1.2–1.5× calibration. This is consistent with no change over 164 days.

## Third epoch: o010 (2022-07-30) against VENUS o052 (2025-07-08), F150W and F444W (D-027)

Outputs under `outputs/t_sunrise_e13/` (not committed). Every number is `derived`; thresholds are ASSUMPTIONs.

- **Catalog stage** (jwst 2.0.1 against 3.0.0): F150W 1,305 matched, 91 variable, 194 appeared, 39 disappeared;
  F444W 984 matched, 51/131/49. Frame shifts ≤ 0.019″.
- **Combination:** 58 positions with the same kind in both bands; 1 near a Gaia source (cone 4′ around
  24.3537, −8.4573); 57 candidates plus Earendel.
- **Forced photometry,** recentred (0.15″ box), with `--controls` (200 random epoch-1 F150W sources at
  25.5–28 mag, ≥ 1″ from any candidate; 184/175 measurable): noise scale **1.30 (F150W), 1.18 (F444W)**. The
  candidates' Δm spans −0.08…+0.12 mag (5–95%); **0 of 57 pass**, raw or calibrated. They are deblending
  differences between the pipeline versions.
- **Earendel** (recentred): Δm = +0.05 (F150W, 0.3σ calibrated) and +0.12 mag (F444W, 0.9σ) over 2.9 years.
- **Controls that change** (|flux σ| ≥ 5 after calibration). Checked on cutouts in every epoch-1 band and in
  epoch 2 (o120: F090W, F115W) at 0.1″ apertures:
  - `n0022` (24.364244, −8.433747): **detector persistence, not a transient** (D-039). In o010 it appears only in
    dithers 3 and 4 of each SW filter, at NRCB4 pixels (21, 298) and (29, 109). A bright galaxy at (24.364688,
    −8.435356) lit pixel (21, 298) in dither 2 and pixel (29, 109) in dither 1 (F150W aperture flux 256 and 268
    against 5.0 and 2.1 afterwards: 2.0 % and 0.8 %; the same ratios in F090W–F200W). The dither geometry puts both
    afterimages on one sky position, so the mosaic shows a "source" that later epochs lack.
  - `n0150` (24.339250, −8.442280): **detector persistence** (D-039). In every epoch it appears only in the
    exposure right after a saturated star at (24.340874, −8.441846) sat on the same pixel: o010 dither 4 (NRCB3
    (100, 929), all SW filters and F277W), o120 dither 2 and o052 dither 4; the preceding exposure has 19–31 saturated pixels
    in the aperture.
    The afterimage holds 0.04–0.07 % of the star's aperture flux. The "different position each epoch" is that
    epoch's dither vector.
  - `n0153` (24.333856, −8.426836): brighter by 0.63 mag in F150W (5.9σ) and 0.52 mag in F444W (3.6σ) in 2025;
    compact and unchanged in shape. Per `_cal` exposure (0.1″ apertures) F150W is 1.5–2.0 in all four o010 dithers
    and 3.0–3.8 in all four o052 dithers (×1.9; jwst 2.0.1 against 3.0.0 calibration). Not persistence. A variable
    (e.g. AGN) candidate; not yet vetted.
  - `scripts/persistence_check.py` on o010, o120 and o052 (404 `_cal` files, S3 byte ranges, 2 min): `n0022` 8 of 8
    detections suspect, `n0150` 5/5, 4/4 and 2/2 at its three positions; `n0153` 0 of 39 and Earendel 0 of 4 suspect
    (`on_sky`).
  - **Why the search missed them:** the combination step needs the same kind in two bands, and in this pair only
    F150W and F444W overlap. `n0022` and `n0150` are blue and faint in F444W. That `n0022`/`n0150` did not appear in the o010/o120
    search (F090W and F115W) is not yet understood (footprint or depth proxy, or the S/N ≥ 10 cut).
