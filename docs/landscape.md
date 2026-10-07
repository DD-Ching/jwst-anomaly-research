# Reuse landscape for future milestones (M2–M4)

This is a decision-oriented survey of the tools, models and data products we can reuse in the
milestones after the M1 catalog slice (ROADMAP.md). The rule (agent charter; DECISIONS.md D-010)
is to start from the "use first" entry. Build something new only after recording why that entry
is not enough.

- **Checked 2026-10-07.** Versions, release dates and licenses come from the GitHub releases API
  and the PyPI JSON API on that date, and match the linked pages. Each arXiv ID was checked
  against its abstract page (the title matches). Links are checked with `scripts/check_links.py`.
- **Verdicts.** *Reuse now* means use it in M0/M1. *Reuse at M2/M3/M4* means adopt it when that
  milestone starts. *Reject* comes with a reason, so nobody re-investigates without new evidence.
- **When to re-check.** Treat this survey as stale after about 6 months, or as soon as a D-010
  "Revisit if" condition holds. Re-verify the version before adopting a tool. Don't re-survey an
  area whose verdict still holds.
- **Not a result.** Nothing here is a scientific result. A new model still has to beat the M1
  baseline on injection-recovery (docs/methodology.md).

## Summary

| Need | When | Use first | Then / instead of |
|---|---|---|---|
| Consistent multi-band photometry + photo-z | M2 | DJA grizli v7.4 catalogs (cover SMACS 0723 and CEERS) | CEERS DR1.0 to validate colours; photutils forced photometry as fallback (§4, §6) |
| Image embeddings | M2 | Zoobot encoders via `timm`; DINOv2 as the generic-vision control | AstroCLIP, AstroPT and AION-1 rejected: they expect other surveys' bands and pixel scales (§1) |
| Unlabelled JWST training data | M2 | Multimodal Universe JWST subset (built from DJA) | (§1) |
| More detectors / active learning | M2 | PyOD (detector zoo), coniferest (active learning) | anomalib rejected; Astronomaly as method/UI reference (§2) |
| Lens finding | M3 | AnomalyMatch (already used on NIRCam); fine-tuned Zoobot | No public lens finder works on NIRCam without fine-tuning (§3a) |
| Lens modelling | M3 | lenstronomy (galaxy scale); PyAutoLens as cross-check | Don't refit the cluster; use published SMACS 0723 models (§3b) |
| Existing SMACS 0723 lens models | M3 | Mahler+2022 (GitHub, CC0), RELICS HLSP models, Caminha+2022 (VizieR), Noirot+2023 redshifts | (§3c) |
| Calibration | now | MAST products; re-run `jwst` only under the §5a rule | Community artifact tools only when re-running (§5c) |
| Artifact vetting | now | JDox known-issues pages + i2d `WHT`/`CON` coverage (the catalog has no flag columns) | (§5b) |
| Run tracking / data versioning | now → M2 | git + manifests; MLflow (local backend) once runs multiply | DVC and DataLad rejected for now (§7) |

**Two findings affect M1 directly** (details in §4 and §5a):

1. **CEERS position is wrong.** For CEERS `jw01345-o001_t021`, MAST `s_ra/s_dec` is the position of
   the MIRI1 prime target, and it lies *outside* the NIRCam footprint. Coverage tests and cone
   searches must use `s_region`.
2. **The two catalogs are not comparable yet.** The program 2736 catalog (jwst 2.0.1, build 12.3)
   and the CEERS catalog (jwst 3.0.0, build 13.0) were made with different `source_catalog`
   deblending defaults. Their shape and neighbour columns can't be compared directly until MAST
   reprocesses program 2736.

## 1. Pretrained astronomy and vision representations (M2)

**Bottom line.** None of the models below was trained on JWST imaging, and we found no public
JWST-finetuned checkpoint.

- **Usable without retraining.** Zoobot and DINOv2 produce embeddings for NIRCam cutouts as they
  are. They need either a 3-band arcsinh RGB composite or, for Zoobot, its greyscale encoders.
- **Not usable as is.** AstroCLIP, AstroPT and AION-1 are survey-specific foundation models. They
  expect Legacy Survey or HSC bands and pixel scales.
- **Whether embeddings help is an empirical question.** An embedding must beat the M1 baseline on
  injection-recovery (docs/methodology.md). The band-to-RGB mapping and the resampling are
  variables to test, not fixed choices.

| Option | What · maintainer | License | Latest release | Training data · JWST? | Verdict |
|---|---|---|---|---|---|
| [Zoobot](https://github.com/mwalmsley/zoobot) | Galaxy-morphology encoders (ConvNeXt, EfficientNet, ...) trained on Galaxy Zoo votes · M. Walmsley | code GPL-3.0; Hugging Face encoder weights Apache-2.0 ([example model card](https://huggingface.co/mwalmsley/zoobot-encoder-convnext_nano)) | v2.9.0 (2025-07-25) | "GZ Evo" set: GZ2, UKIDSS, Hubble, CANDELS, DECaLS/DESI, HSC ([docs](https://zoobot.readthedocs.io/en/latest/pretrained_models.html)); a Euclid encoder exists; no JWST checkpoint | **Reuse at M2** as the first backbone. Load the weights with `timm` (`hf_hub:` id), which avoids importing the GPL package |
| [DINOv2](https://github.com/facebookresearch/dinov2) | General self-supervised ViT · Meta | Apache-2.0 (code and standard weights) | no tagged releases (Torch Hub) | 142M natural images; no astronomy | **Reuse at M2** as the generic-vision control |
| [DINOv3](https://github.com/facebookresearch/dinov3) | Successor of DINOv2 · Meta | custom [DINOv3 License](https://github.com/facebookresearch/dinov3/blob/main/LICENSE.md) (gated download, attribution, use restrictions) | no tagged releases | web images (plus satellite variants); no astronomy | **Reject for now.** The gated download and custom license hurt reproducibility. Revisit if DINOv2 shows clear value |
| [AstroCLIP](https://github.com/PolymathicAI/AstroCLIP) | Image–spectrum contrastive model · Polymathic AI | MIT | tag `mnras` (2024-06-18); not on PyPI | Legacy Survey g,r,z, 144 px; no JWST | **Reject.** Bands and scale don't match ours, and the spectrum alignment adds nothing here |
| [AstroPT](https://github.com/Smith42/astroPT) | Autoregressive transformer on galaxy stamps · M. J. Smith | code AGPL-3.0; [weights](https://huggingface.co/Smith42/astroPT_v2.0) CC-BY-SA-4.0 | v2.0.9 (2025-09-09) | Legacy Survey DR8 grz (plus a Euclid variant); no JWST | **Reject for now.** Same mismatch, plus copyleft code and share-alike weights |
| [AION-1](https://github.com/PolymathicAI/AION) | Multimodal model with per-survey tokenisers · Polymathic AI | MIT | v0.0.2 (2025-06-27) | Legacy Survey, HSC, SDSS/DESI spectra, Gaia; no JWST modality | **Reject.** Its image tokenisers are tied to Legacy Survey/HSC bands |
| [Multimodal Universe](https://github.com/MultimodalUniverse/MultimodalUniverse) (dataset) | ML-ready astronomy data · MMU collaboration. Includes a [JWST subset](https://huggingface.co/datasets/MultimodalUniverse/jwst) from DJA fields (CEERS, PRIMER, GOODS-S/N, NGDEEP): F090W–F444W, 96 px cutouts | repo MIT; the JWST data asks for DJA acknowledgement (terms not fully settled) | v1.5 (2026-06-17) | is itself JWST data | **Reuse at M2** as unlabelled JWST "normal" data for finetuning or self-supervision. **Exclude CEERS** (our evaluation field) and any other field being scored from the training data, or the injection-recovery comparison is biased |
| [timm](https://github.com/huggingface/pytorch-image-models) | Model zoo and loader · Hugging Face | Apache-2.0 | v1.0.30 (2026-09-22) | n/a | **Reuse at M2** as the single loading path for Zoobot and DINO weights |

**Read before M2:**
- Galaxy Zoo JWST: CEERS volunteer morphologies,
  "Galaxy Zoo JWST: Up to 75% of discs are featureless at 3<z<7" ([arXiv:2503.21869](https://arxiv.org/abs/2503.21869)).
- "The Platonic Universe" ([arXiv:2509.19453](https://arxiv.org/abs/2509.19453)), which compares
  foundation-model embeddings on cross-matched imagery that includes JWST.

## 2. Anomaly-detection frameworks (vs the sklearn baseline)

The M1 baseline (unit 3: robust z-scores, IsolationForest and LOF in scikit-learn) stays the
reference.

| Option | What · maintainer | License | Latest release | Fit | Verdict |
|---|---|---|---|---|---|
| [PyOD](https://github.com/yzhao062/pyod) ([arXiv:1901.01588](https://arxiv.org/abs/1901.01588)) | About 60 outlier detectors behind one `fit`/`decision_function` API · Y. Zhao | BSD-2 | v3.6.6 (2026-09-17) | Superset of IF/LOF; takes feature matrices, including embeddings | **Reuse at M2** once the detector set grows beyond IF/LOF. Not needed for M1 |
| [coniferest](https://github.com/snad-space/coniferest) ([arXiv:2410.17142](https://arxiv.org/abs/2410.17142)) | Isolation forest plus active anomaly detection (AAD, PineForest) · SNAD | MIT | v0.2.1 (2026-07-18) | Re-ranks with a human in the loop, on features or embeddings; pip-installable | **Reuse at M2** for active learning on vetted candidates |
| [Astronomaly](https://github.com/MichelleLochner/astronomaly) ([arXiv:2010.11202](https://arxiv.org/abs/2010.11202); "Protege" [arXiv:2411.04188](https://arxiv.org/abs/2411.04188)) | Active anomaly-detection framework with a web labelling UI · M. Lochner | BSD-3 | v2.0 (2024-12-13); not on PyPI | Same detectors as ours plus feedback-driven ranking. Protege extends it to self-supervised features | **Reuse at M2 as a method/UI reference.** Prefer coniferest as the library dependency |
| [AnomalyMatch](https://github.com/esa/AnomalyMatch) ([arXiv:2505.03509](https://arxiv.org/abs/2505.03509)) | Image classifier combining semi-supervised learning (FixMatch) with active learning · ESA | ESA-PL Permissive v2.4 | v1.3.2 (2026-07-08); install from GitHub | Works on images (FITS input, 1–N channels). Already used to find lenses in JWST NIRCam ([arXiv:2605.03442](https://arxiv.org/abs/2605.03442)) | **Reuse at M3** for the lens/arc search; optional image-level competitor at M2 |
| [anomalib](https://github.com/open-edge-platform/anomalib) ([arXiv:2202.08341](https://arxiv.org/abs/2202.08341)) | Visual anomaly detection for industrial inspection · Open Edge Platform | Apache-2.0 | v2.7.0 (2026-10-06) | Trains on normal-only RGB images to produce pixel-level defect maps. Our samples are unlabelled and may already contain the anomalies | **Reject** |

## 3. Gravitational-lens finding and lens modelling (M3)

### 3a. Lens finders

No published lens finder runs on NIRCam cutouts as is. Every one we found was trained on
ground-based or Euclid imaging, so each needs fine-tuning on JWST data.

| Option | Weights public? | License | Fit | Verdict |
|---|---|---|---|---|
| AnomalyMatch (see §2) | Code only; it trains from a handful of labels | ESA-PL Permissive v2.4 | Already used on NIRCam lenses ([arXiv:2605.03442](https://arxiv.org/abs/2605.03442)) | **Reuse at M3** (first choice) |
| Zoobot Euclid encoder ([HF](https://huggingface.co/mwalmsley/zoobot-encoder-euclid)) + Euclid Q1 expert lens labels ([HF dataset](https://huggingface.co/datasets/mwalmsley/euclid_strong_lens_expert_judges)) | Encoder yes. The fine-tuned Q1 lens finder ([arXiv:2503.15326](https://arxiv.org/abs/2503.15326)) was not found | Apache-2.0 (encoder) | Backbone to fine-tune on JWST lens examples | **Reuse at M3** as the fine-tuning backbone |
| [ssl-legacysurvey](https://github.com/georgestein/ssl-legacysurvey) (Stein et al., [arXiv:2110.00023](https://arxiv.org/abs/2110.00023)) | Yes (via Globus) | MIT | Legacy Survey grz only | **Reject** (method reference only) |
| [CMU DeepLens](https://github.com/McWilliamsCenter/CMUDeepLens) ([arXiv:1703.02642](https://arxiv.org/abs/1703.02642)) | No | MIT | Theano/Lasagne; last push 2018 | **Reject** |
| HOLISMOKES (HSC), DESI (Huang/Storfer) finders | Not found (unverified) | n/a | n/a | Revisit if weights are released |

### 3b. Lens-modelling codes

| Option | Scope | License | Latest release | Verdict |
|---|---|---|---|---|
| [lenstronomy](https://github.com/lenstronomy/lenstronomy) ([docs](https://lenstronomy.readthedocs.io/en/latest/)) | Galaxy-scale modelling from imaging; [JAXtronomy](https://github.com/lenstronomy/JAXtronomy) is a JAX/GPU port of a subset | BSD-3 | v1.14.2 (2026-07-09) | **Reuse at M3** for galaxy-scale lenses and arcs in the field |
| [PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens) ([docs](https://pyautolens.readthedocs.io/en/latest/)) | Galaxy (and group) scale; automated pipelines; JAX/GPU | MIT | 2026.10.4.1 (2026-10-04) | **Reuse at M3** as an alternative or cross-check to lenstronomy |
| [herculens](https://github.com/Herculens/herculens) | Differentiable JAX modelling, including pixelated sources | MIT | v0.3.0 (2026-04-02) | **Reuse at M3 only** if subhalo or pixelated-potential tests are needed |
| Cluster-scale codes (details below) | Cluster mass models | mixed | n/a | **Don't refit at M3.** Use the published cluster models (§3c). If refitting is ever needed, glafic is the easiest to install |

The cluster-scale codes:
- [glafic](https://github.com/oguri/glafic2): GPL-3.0, v2.1.15 (2026-06-24).
- [Lenstool](https://projets.lam.fr/projects/lenstool/wiki): license not stated.
- Grale2.
- WSLAP+: not public.
- Zitrin LTM: no public code.

None of lenstronomy, PyAutoLens or herculens advertises cluster-scale modelling.

### 3c. Public lens models and catalogs for SMACS J0723.3-7327

Before calling a candidate "unusual", M3 can use these to test whether its geometry (arc
orientation, number of images, magnification) is consistent with existing cluster models.

| Product | Method · era | Where | Terms | Verdict |
|---|---|---|---|---|
| RELICS HLSP lens models `glafic/v2`, `lenstool/v1`, `zitrin-ltm-gauss/v2`: kappa, gamma, deflection, and magnification at z = 1, 2, 4, 9, as FITS | HST era | [HLSP page](https://archive.stsci.edu/hlsp/relics), [files](https://archive.stsci.edu/hlsps/relics/smacs0723m73/models/) | CC BY 4.0, DOI 10.17909/T9SP45 | **Reuse at M3** as the HST-era reference set |
| Mahler et al. 2022 ([arXiv:2207.07101](https://arxiv.org/abs/2207.07101)) | Lenstool, JWST | [GitHub](https://github.com/guillaumemahler/SMACS0723-mahler2022) (release v2, 2022-12-27): mass, kappa, gamma, deflection and magnification FITS, plus multiple-image catalogs | CC0-1.0 | **Reuse at M3** as the primary JWST-era model |
| Caminha et al. 2022 ([arXiv:2207.07567](https://arxiv.org/abs/2207.07567); A&A 666, L9) | Parametric; JWST + MUSE | CDS catalogue J/A+A/666/L9 ([page](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/666/L9)); file contents not verified | CDS terms | **Reuse at M3** as an independent model, with MUSE redshifts |
| Pascale et al. 2022 ([arXiv:2207.07102](https://arxiv.org/abs/2207.07102)); Golubchik et al. 2022 ([arXiv:2207.05007](https://arxiv.org/abs/2207.05007)) | Parametric; LTM (HST only) | Dropbox links in the papers (not verified to still work) | not stated | Optional cross-checks |
| Diego et al. 2023 ([arXiv:2301.03629](https://arxiv.org/abs/2301.03629)) | Free-form | No availability statement | n/a | Skip unless obtained from the authors |
| Noirot et al. 2023 spectroscopic redshifts ([arXiv:2212.07366](https://arxiv.org/abs/2212.07366)) | NIRISS; 190 secure redshifts | [niriss.github.io/smacs0723](https://niriss.github.io/smacs0723) | not stated | **Reuse at M3** for multiple-image redshifts |

**Access pattern.** RELICS and the Hubble Frontier Fields ([HFF lens models](https://archive.stsci.edu/prepds/frontier/lensmodels/))
share one directory layout, `.../{cluster}/models/{team}/{version}/hlsp_..._{product}.fits`, so
one loader can read SMACS 0723 and the HFF clusters. The RELICS deflection-map normalisation is
unverified (the README PDF was not parsed); check the FITS headers before use.

MAST DOIs (10.17909/...) appear as plain text in this file. On 2026-10-07, doi.org redirected
them to an `http://archive.stsci.edu/doi/resolve/...` URL that returned 404; the same URL over
`https://` worked. Link the HLSP pages instead.

## 4. JWST HLSP catalogs with consistent multi-band photometry and photo-z (M2)

**Bottom line.**
- **DJA covers both fields.** The [DAWN JWST Archive (DJA)](https://dawn-cph.github.io/dja/) is
  the only source we verified that covers **both** reference fields with one pipeline: grizli
  mosaics, one detection image, matched circular apertures and eazy photo-z.
- **DJA photometry is not PSF-homogenised.** It is matched-aperture only
  ([method note](https://dawn-cph.github.io/dja/blog/2023/07/14/photometric-catalog-demo/)). For
  CEERS, the team's PSF-matched DR1.0 catalog is the reference for validating DJA colours.
- **SMACS 0723 has no other JWST catalog.** Apart from DJA, the only team HLSP covering it is
  RELICS, which is HST only.

**Footprint caveat (verified 2026-10-07).**
- **CEERS t021:** for `jw01345-o001_t021_nircam_clear-f200w`, MAST gives `s_ra/s_dec` =
  (215.1620, 53.0513), which is the position of the prime MIRI1 target. The NIRCam parallel's
  `s_region` spans RA 214.913–215.053 and Dec 52.934–53.022, so that position lies outside the
  NIRCam data. Coverage tests and cone searches must use `s_region`, not `s_ra/s_dec`.
- **Program 2736:** the `s_ra/s_dec` position does lie inside the footprint.

| Survey | Field(s) | Latest release (as stated) | Forced/matched phot · photo-z | Access · terms | SMACS 0723 | CEERS t021 |
|---|---|---|---|---|---|---|
| **DJA** (grizli v7) | many public fields, including `smacs0723` and `ceers-full` | v7.4 catalogs: smacs0723 files dated 2025-03-05, ceers-full 2024-09-27 | matched apertures (not PSF-matched) · eazy | [file index](https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html); FITS; no license stated; cite [arXiv:2302.10936](https://arxiv.org/abs/2302.10936) | **yes** | **yes** (per-exposure WCS) |
| [CEERS DR1.0](https://ceers.github.io/dr1.html) | EGS, NIRCam pointings 1–10 | October 2025 | PSF-matched, 7 NIRCam + 7 HST bands · LePHARE | FITS, >80k galaxies; CC BY 4.0 per the catalog paper ([arXiv:2510.08743](https://arxiv.org/abs/2510.08743)); on MAST "soon" | no | yes, inferred (Epoch-1 NIRCam parallels to MIRI; not checked against coordinates) |
| [CEERS HLSP (MAST)](https://archive.stsci.edu/hlsp/ceers) | EGS | updated 2023-08-30 | mosaics only, no catalog | CC BY 4.0, DOI 10.17909/z7p0-8481 | no | imaging only |
| [JADES](https://archive.stsci.edu/hlsp/jades) | GOODS-S/N | MAST DR3 (2024); DR5 catalog paper [arXiv:2601.15956](https://arxiv.org/abs/2601.15956) | forced Kron/circular apertures on PSF-matched mosaics · yes | FITS; CC BY 4.0, DOI 10.17909/8tdj-8n28 | no | no |
| [COSMOS2025](https://cosmos2025.iap.fr/) (COSMOS-Web) | COSMOS, ~0.54 deg² | catalog v1.1, 2026-05-04 | PSF-homogenised apertures + SourceXtractor++ models · LePHARE | one FITS file, 784,016 sources; cite [arXiv:2506.03243](https://arxiv.org/abs/2506.03243) | no | no |
| PRIMER | COSMOS, UDS | no team catalog found (MAST `hlsp/primer` returned 404 on 2026-10-07) | via DJA only | n/a | no | no |
| [UNCOVER](https://jwst-uncover.github.io/) | Abell 2744 | DR4 (2024-08-07; v4.1 2024-12-05) | PSF-matched to F444W · EAZY, Prospector | team-site downloads; cite [arXiv:2301.02671](https://arxiv.org/abs/2301.02671) | no | no |
| [CANUCS](https://archive.stsci.edu/hlsp/canucs) | A370, MACS0416, MACS0417, MACS1149, MACS1423 | DR1, 2025-06-25 | PSF-matched to F444W · EAzY | FITS; CC BY 4.0, DOI 10.17909/18nv-np70; [arXiv:2506.21685](https://arxiv.org/abs/2506.21685) | no | no |
| [GLASS-JWST](https://archive.stsci.edu/hlsp/glass-jwst) | Abell 2744 | 2023-07-07 | forced on PSF-matched images · no photo-z | FITS, 24,389 sources; CC BY 4.0, DOI 10.17909/kw3c-n857; [arXiv:2301.02179](https://arxiv.org/abs/2301.02179) | no | no |
| [RELICS](https://archive.stsci.edu/hlsp/relics) (HST) | HST lensing clusters, including SMACS J0723.3-7327 | updated 2022-09-30 | SExtractor catalogs · BPZ-style photo-z | CC BY 4.0, DOI 10.17909/T9SP45 | yes (HST bands) | no |

DJA v7.4 files we verified (HTTP HEAD, 2026-10-07), all under
`https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/`:

| File | Size |
|---|---|
| `smacs0723-grizli-v7.4-fix_phot.fits` | 38.6 MB |
| `smacs0723-grizli-v7.4-fix.photoz.tar.gz` | 59.0 MB |
| `ceers-full-grizli-v7.4-fix_phot.fits` | 250.5 MB |
| `ceers-full-grizli-v7.4-fix.photoz.tar.gz` | 396.8 MB |

Both `ceers-full` files are above the 200 MB single-download guideline in CLAUDE.md. Before
downloading either one, state the reason in the PR, or use CEERS DR1.0 instead.

**Verdict.**
- **Reuse at M2:** the DJA v7.4 catalogs for both fields (one schema), CEERS DR1.0 to validate
  the CEERS colours, and RELICS for the HST bands in SMACS 0723.
- **Fallback:** photutils forced photometry (§6), only where the colours from these catalogs
  prove inadequate.
- **M4:** JADES, COSMOS2025, UNCOVER, CANUCS and GLASS are the natural next fields; all have
  PSF-matched catalogs and photo-z (GLASS has no photo-z).

## 5. JWST calibration and artifact handling (vetting, M1–M2)

### 5a. Archive products vs re-running the `jwst` pipeline

- **Reprocessing.** Since build 11.1, each operations build is tied to one fixed CRDS "build
  context". STScI reprocesses current and prior MAST data with each new build, "typically taking
  approximately 4–6 weeks after build installation"
  ([build information](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information);
  [quarterly CRDS](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/crds-migration-to-quarterly-calibration-updates)).
- **Pedigree and local runs.** The FITS headers `CAL_VER` (pipeline version) and `CRDS_CTX` record
  how a product was made
  ([MAST: updates to JWST data](https://outerspace.stsci.edu/display/MASTDOCS/Updates+to+JWST+Data)).
  Local runs need `CRDS_PATH` and `CRDS_SERVER_URL=https://jwst-crds.stsci.edu`, and optionally a
  pinned `CRDS_CONTEXT`
  ([pipeline docs](https://jwst-pipeline.readthedocs.io/en/latest/jwst/user_documentation/reference_files_crds.html)).
- **Our two reference catalogs come from different builds** (checked 2026-10-07):
  - Program 2736 F200W: jwst 2.0.1 = build 12.3 (archived), context `jwst_1535.pmap`.
  - CEERS t021: jwst 3.0.0 = build 13.0 (in operations since 2026-09-08), context `jwst_1584.pmap`.
  - Build 13.0 "changed the default behavior of the imaging source_catalog step to apply
    deblending instead of merging all contiguous pixels". Until 2736 is reprocessed (expected
    within ~6 weeks of 2026-09-08), the shape, neighbour and concentration columns can't be
    compared directly between the two programs.
- **Some corrections are off by default in operations**, so only a re-run applies them:
  - `clean_flicker_noise` for 1/f noise ([docs](https://jwst-pipeline.readthedocs.io/en/latest/jwst/clean_flicker_noise/main.html)).
  - Persistence flagging, which is inactive unless `persistence_time` is set
    ([docs](https://jwst-pipeline.readthedocs.io/en/latest/jwst/persistence/description.html)).
  - NIRCam wisp subtraction, which uses offline templates and is not a pipeline step.

**Rule (reuse now).** Use MAST level-3 products. Re-run the pipeline from `_rate`/`_uncal` (with
the current build's `jwst` version unless there is a reason not to) only when:
1. a product's `CAL_VER`/`CRDS_CTX` lags the operational build and MAST has not caught up;
2. a top candidate could plausibly be 1/f striping, a wisp, persistence or a saturated core, and
   the fix is off by default; or
3. a comparison needs identical catalog settings across programs.

Record `CAL_VER`/`CRDS_CTX` in the manifests so this can be checked.

### 5b. Artifacts to rule out before interpreting a candidate

| Artifact | Instr. | Appearance | Pipeline handling | Documentation |
|---|---|---|---|---|
| Snowballs / showers | NIRCam / MIRI | Round residuals (NIRCam) or elongated, diffuse ones (MIRI); they can survive when there are few dithers | `jump` large-event flagging (MIRI shower flagging is not applied in F1800W–F2550W); `outlier_detection` with ≥4 dithers | [JDox](https://jwst-docs.stsci.edu/known-issues/shower-and-snowball-artifacts), [jump arguments](https://jwst-pipeline.readthedocs.io/en/latest/jwst/jump/arguments.html) |
| Wisps, claws, dragon's breath | NIRCam (mostly SW) | Diffuse patterns at fixed detector positions; claw-shaped stray light from bright stars outside the field | None in the pipeline (wisp templates are applied offline) | [JDox scattered light](https://jwst-docs.stsci.edu/known-issues/nircam-known-issues/nircam-scattered-light-artifacts) |
| 1/f noise | NIRCam (rows), MIRI (faint columns) | Stripes across the detector | `refpix` corrects part of it; `clean_flicker_noise` is off by default | [JDox 1/f](https://jwst-docs.stsci.edu/known-issues/nircam-known-issues/nircam-1-f-noise-removal-methods) |
| Persistence | NIRCam, MIRI | Latent copies of earlier bright sources | Flagging only, and only if opted in | [JDox NIRCam persistence](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-persistence) |
| Diffraction spikes, saturation | both | Spikes aligned with the telescope V3 axis; NaN or flagged cores | `saturation` flags | [JDox NIRCam PSFs](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions), [NIRCam known issues](https://jwst-docs.stsci.edu/known-issues/nircam-known-issues) |
| Cruciform, tree rings, striping, row/column pull-up/down, edge glow | MIRI | A cross around bright sources (≤10 µm), stripes, bright edges | None for imaging | [MIRI imaging known issues](https://jwst-docs.stsci.edu/known-issues/miri-known-issues/miri-imaging-known-issues) |
| Low-coverage edges, outlier overflagging | both | Few contributing exposures, holes, doubled sources | `outlier_detection` (DQ OUTLIER) | [outlier_detection](https://jwst-pipeline.readthedocs.io/en/latest/jwst/outlier_detection/main.html) |

**The catalog alone is not enough for vetting.** The `_cat.ecsv` has **no** saturation, edge or
DQ flag columns ([column list](https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/main.html)).
Vetting must join catalog positions with the i2d `WHT` and `CON` arrays, which give the coverage
and the contributing exposures.

Useful catalog columns are `is_extended`, `sharpness`, `roundness`, `CI_*`, `nn_dist` and
`ellipticity`. Heuristics:
- Very high sharpness or a low concentration index suggests a cosmic ray or hot pixel.
- High ellipticity or roundness suggests a fragment of a spike or streak.
- A small `nn_dist` to a bright star suggests spikes, ghosts or scattered light.

### 5c. Community artifact tools

All of these work on level-1/2 products, not on i2d mosaics.

| Tool | Maintainer · license | Latest | Purpose | Verdict |
|---|---|---|---|---|
| [snowblind](https://github.com/mpi-astronomy/snowblind) | MPIA · BSD-3 | 0.2.1 (on PyPI since 2023-12; repo pushed 2026-03) | Snowball/shower masks; persistence and open-pixel flags | **Reuse at M2** if re-running stages 1–3 (compatibility with jwst 3.0 not documented) |
| [grizli](https://github.com/gbrammer/grizli) | G. Brammer · MIT | 1.14.2 (2026-07-24) | Full HST/JWST preprocessing, including wisp, 1/f, snowball and persistence handling (the DJA mosaics are made with it) | **Reuse at M2** only if we must make our own mosaics; prefer the DJA products |
| [chriswillott/jwst](https://github.com/chriswillott/jwst) | C. Willott · GPL-3.0 | no releases (pushed 2025-12) | `image1overf.py` for 1/f removal; snowball and persistence flagging scripts | **Reference only.** Use it as an external cross-check; don't vendor it (GPL) |
| [ceers-nircam](https://github.com/ceers/ceers-nircam) | CEERS team · no license | last push 2023-05 (targets jwst 1.7.2) | Wisp, 1/f and snowball scripts | **Reject as code.** The methods are described in [arXiv:2211.02495](https://arxiv.org/abs/2211.02495) |

## 6. Source extraction and forced photometry

The pipeline's `source_catalog` step already runs photutils, so our M1 catalogs are photutils
output. The open question for M2 is consistent multi-band (forced) photometry in fields that no
HLSP catalog covers (§4).

| Tool | Maintainer · license | Latest release | Fit | Verdict |
|---|---|---|---|---|
| [photutils](https://github.com/astropy/photutils) | Astropy project · BSD-3 | 3.0.0 (2026-04-17) | `SourceCatalog(..., detection_catalog=...)` reuses one detection image's segmentation, centroids and apertures in every band, i.e. forced photometry ([docs](https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceCatalog.html)). The `psf_matching` module builds PSF-matching kernels ([PSF matching](https://photutils.readthedocs.io/en/stable/user_guide/psf_matching.html)). It is the code the pipeline already uses: pure Python, installs with pip | **Reuse at M2** as the default forced-photometry tool |
| [SEP](https://github.com/sep-developers/sep) | sep-developers (originally K. Barbary) · LGPLv3 as a whole (the Python wrapper is MIT) | 1.4.1 (2025-02-18) | Source Extractor algorithms on NumPy arrays. Faster than photutils, but duplicates what photutils already gives us | **Reject for now.** Revisit only if photutils is a measured bottleneck at M4 |
| [SourceXtractor++](https://github.com/astrorama/SourceXtractorPlusPlus) | astrorama (Euclid-funded) · LGPL-3.0 | 1.1.0 (2026-07-07) | Model-fitting photometry across several images with different PSFs and pixel grids ([docs](https://astrorama.github.io/SourceXtractorPlusPlus/)). Installs only with conda (`-c astrorama -c conda-forge`); heavy | **Reuse at M3** only if blended arcs need simultaneous model fitting for deblending |

## 7. Experiment tracking and data versioning

Current practice (D-001):
- git holds the code and configs.
- Tracked manifests (URI + sha256 + size + pipeline version) record the archive data.
- The data itself stays in an immutable public archive (MAST/S3).

| Tool | Maintainer · license | Latest release | Fit | Verdict |
|---|---|---|---|---|
| git + manifests + candidate-store run records | this repo | n/a | Enough while runs are few and the inputs are public archive files | **Reuse now** |
| [MLflow](https://github.com/mlflow/mlflow) | MLflow project · Apache-2.0 | 3.16.1 (2026-09-17) | Tracks params, metrics and artifacts per run. Works locally (an `mlruns` directory or SQLite) with no server ([docs](https://mlflow.org/docs/latest/ml/tracking/)). Useful once M2 compares many models and seeds | **Reuse at M2** once the run count outgrows a results table; local backend only |
| [DVC](https://github.com/treeverse/dvc) | Treeverse (formerly iterative) · Apache-2.0 | 3.67.1 (2026-03-31) | Git metafiles plus a cache and remotes such as S3, HTTP or SSH ([docs](https://doc.dvc.org/user-guide/data-management/remote-storage)). Our raw inputs already have a permanent public "remote" (MAST), so DVC would store copies of files that our manifests only reference | **Reject for now.** Revisit at M4 if large derived sets (embeddings, cutout banks) must be shared |
| [DataLad](https://github.com/datalad/datalad) | DataLad team · MIT | 1.7.1 (2026-10-06) | Data version control on git + git-annex; `datalad run` records command provenance ([handbook](https://handbook.datalad.org/en/latest/)). git-annex is an extra system dependency on Windows | **Reject for now**, with the same revisit trigger as DVC |

## Open and unverified items

- **Lens-finder weights.** We found no public weights for HOLISMOKES, the DESI (Huang/Storfer)
  finders or the fine-tuned Euclid Q1 lens finder. This could change.
- **SMACS 0723 model details.** Not verified: the RELICS deflection-map normalisation, the
  contents of the Caminha+2022 VizieR files, and whether the Pascale and Golubchik Dropbox links
  still work.
- **Data licenses.** DJA, COSMOS2025 and UNCOVER state no explicit data license; they ask for
  citations. The Multimodal Universe JWST subset asks for a DJA acknowledgement.
- **CEERS DR1.0 coverage of t021.** This is inferred (Epoch-1 NIRCam parallels), not checked
  against catalog coordinates. Check it with the `s_region` polygon.
- **Artifact tools.** Unknown whether snowblind works with jwst 3.0; the license and date of the
  NIRCam wisp templates are unknown.
- **License notes are not legal advice.** They are factual summaries. GPL/AGPL tools (Zoobot
  code, AstroPT, glafic, chriswillott/jwst) stay optional dependencies and are never vendored into
  this BSD-3-Clause repository.
