# Sources

Enough information to recover each external source later: URL, version or DOI where relevant, and
the date it was checked. Sections marked with a unit number are owned by that bootstrap unit;
edit only your own section in parallel work.

## Data (checked 2026-10-07)

- **MAST** (Mikulski Archive for Space Telescopes), CAOM API: https://mast.stsci.edu/api/v0/_services.html ;
  Python examples: https://mast.stsci.edu/api/v0/pyex.html
- **JWST program 2736**, ERO SMACS J0723.3-7327: https://www.stsci.edu/jwst/science-execution/program-information?id=2736
  Level-3 imaging: NIRCam F090W/F150W/F200W/F277W/F356W/F444W (obs `o001`), MIRI F770W/F1000W/F1500W/F1800W
  (obs `o002`), all PUBLIC (release MJD 59773.625). NIRCam F200W: `_cat.ecsv` 3,444,394 B (3,145 rows,
  59 columns, made with jwst 2.0.1 / photutils 2.3.0), `_i2d.fits` 1,762,162,560 B. MIRI F770W `_i2d.fits` 36,083,520 B.
- **JWST program 1345** (CEERS): https://www.stsci.edu/jwst/science-execution/program-information?id=1345
  Control pointing `jw01345-o001_t021_nircam_clear-f200w`: `_cat.ecsv` 4,017,583 B (3,695 rows, jwst 3.0.0),
  `_i2d.fits` 1,248,422,400 B.
- **STScI public S3 bucket** `stpubdata` (anonymous): https://stpubdata.s3.amazonaws.com/ — level-3 key pattern
  `jwst/public/jw{PPPPP}/L3/t/o{OOO}/{obs_id}_{suffix}` (verified for 2736 and 1345).

## Core software

- **Astropy**: https://www.astropy.org — cloud FITS byte-range access:
  https://docs.astropy.org/en/stable/io/fits/usage/cloud.html
- **astroquery.mast**: https://astroquery.readthedocs.io/en/latest/mast/mast_obsquery.html
- **JWST pipeline source catalog step** (defines `_cat.ecsv` columns):
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/main.html

## Agent tooling (unit 8 extends)

- Claude Code skills: https://code.claude.com/docs/en/skills ; routines (`/schedule`):
  https://code.claude.com/docs/en/routines
- `anthropics/claude-code-action`: https://github.com/anthropics/claude-code-action

## Archive access (unit 1)

_Pending._

## Catalogs and cross-band matching (unit 2)

_Pending._

## Features and anomaly-detection methods (unit 3)

_Pending._

## Imaging, cutouts and visualization (unit 4)

_Pending._

## External catalog services (unit 5)

_Pending._

## Candidate store and run provenance (unit 6)

_Pending._

## Open-source project tooling (unit 7)

_Pending._

## Future milestones: models, frameworks, lensing, HLSPs, tracking (unit 9)

Everything here was accessed on 2026-10-07. Versions are the latest GitHub release (or PyPI
release) on that date. Verdicts and context are in [docs/landscape.md](docs/landscape.md).

**Software and models**
- Zoobot v2.9.0: https://github.com/mwalmsley/zoobot (encoder weights: https://huggingface.co/mwalmsley/zoobot-encoder-convnext_nano)
- DINOv2, no tagged release; `main` at 7764ea0f912e (2026-06-03): https://github.com/facebookresearch/dinov2
- DINOv3, no tagged release; `main` at 6876159a11b4 (2026-07-15); custom license: https://github.com/facebookresearch/dinov3
- AstroCLIP, tag `mnras`: https://github.com/PolymathicAI/AstroCLIP
- AstroPT v2.0.9: https://github.com/Smith42/astroPT ; weights (CC-BY-SA-4.0): https://huggingface.co/Smith42/astroPT_v2.0
- AION-1 v0.0.2: https://github.com/PolymathicAI/AION
- timm v1.0.30: https://github.com/huggingface/pytorch-image-models
- PyOD v3.6.6: https://github.com/yzhao062/pyod
- coniferest v0.2.1: https://github.com/snad-space/coniferest
- Astronomaly v2.0: https://github.com/MichelleLochner/astronomaly
- AnomalyMatch v1.3.2: https://github.com/esa/AnomalyMatch
- anomalib v2.7.0: https://github.com/open-edge-platform/anomalib
- lenstronomy v1.14.2: https://github.com/lenstronomy/lenstronomy
- JAXtronomy, no tagged release: https://github.com/lenstronomy/JAXtronomy
- PyAutoLens 2026.10.4.1: https://github.com/PyAutoLabs/PyAutoLens
- herculens v0.3.0: https://github.com/Herculens/herculens
- glafic v2.1.15: https://github.com/oguri/glafic2
- Lens finders (rejected; kept for the record):
  - ssl-legacysurvey: https://github.com/georgestein/ssl-legacysurvey
  - CMU DeepLens: https://github.com/McWilliamsCenter/CMUDeepLens
- jwst 3.0.0: https://github.com/spacetelescope/jwst
- CRDS 14.0.2: https://github.com/spacetelescope/crds
- snowblind 0.2.1: https://github.com/mpi-astronomy/snowblind
- grizli 1.14.2: https://github.com/gbrammer/grizli
- chriswillott/jwst (no release; pushed 2025-12-18): https://github.com/chriswillott/jwst
- ceers-nircam (no license; last push 2023-05-03): https://github.com/ceers/ceers-nircam
- photutils 3.0.0: https://github.com/astropy/photutils
- SEP 1.4.1: https://github.com/sep-developers/sep
- SourceXtractor++ 1.1.0: https://github.com/astrorama/SourceXtractorPlusPlus
- MLflow 3.16.1: https://github.com/mlflow/mlflow
- DVC 3.67.1: https://github.com/treeverse/dvc
- DataLad 1.7.1: https://github.com/datalad/datalad
- Link-checker alternatives (rejected for `scripts/check_links.py`, see D-010):
  - lychee v0.24.2: https://github.com/lycheeverse/lychee
  - lychee-action v2.9.0: https://github.com/lycheeverse/lychee-action
  - markdown-link-check v3.15.0: https://github.com/tcort/markdown-link-check
  - LinkChecker v10.6.0: https://github.com/linkchecker/linkchecker

**Data products**
- DJA grizli v7.4 catalogs (`smacs0723-*`, `ceers-full-*`): https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html
  ; overview: https://dawn-cph.github.io/dja/
- Multimodal Universe JWST subset, release v1.5: https://huggingface.co/datasets/MultimodalUniverse/jwst
- CEERS DR1.0 catalog (October 2025): https://ceers.github.io/dr1.html
- MAST HLSPs, each with its DOI:
  - RELICS, DOI 10.17909/T9SP45, includes the SMACS 0723 lens models: https://archive.stsci.edu/hlsp/relics
  - CEERS, DOI 10.17909/z7p0-8481: https://archive.stsci.edu/hlsp/ceers
  - JADES, DOI 10.17909/8tdj-8n28: https://archive.stsci.edu/hlsp/jades
  - CANUCS, DOI 10.17909/18nv-np70: https://archive.stsci.edu/hlsp/canucs
  - GLASS-JWST, DOI 10.17909/kw3c-n857: https://archive.stsci.edu/hlsp/glass-jwst
- COSMOS2025 catalog v1.1: https://cosmos2025.iap.fr/
- UNCOVER DR4: https://jwst-uncover.github.io/
- Mahler+2022 SMACS 0723 lens model, release v2 (2022-12-27): https://github.com/guillaumemahler/SMACS0723-mahler2022
- Noirot+2023 SMACS 0723 redshifts: https://niriss.github.io/smacs0723
- Caminha+2022 SMACS 0723 model and MUSE redshifts, CDS catalogue J/A+A/666/L9 (file contents not verified): https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/666/L9
- RELICS SMACS 0723 lens-model files (`glafic/v2`, `lenstool/v1`, `zitrin-ltm-gauss/v2`): https://archive.stsci.edu/hlsps/relics/smacs0723m73/models/
- Hubble Frontier Fields lens models (same layout as RELICS): https://archive.stsci.edu/prepds/frontier/lensmodels/
- Zoobot Euclid encoder (Apache-2.0): https://huggingface.co/mwalmsley/zoobot-encoder-euclid
- Euclid Q1 strong-lens expert labels: https://huggingface.co/datasets/mwalmsley/euclid_strong_lens_expert_judges

**Papers** (title as shown on arXiv)
- Lens finding and lens models:
  - arXiv:2207.07101: Mahler et al., "Precision modeling of JWST's first cluster lens SMACSJ0723.3-7327"
  - arXiv:2207.07567: Caminha et al., "First JWST observations of a gravitational lens: Mass model from new multiple images with near-infrared observations of SMACS J0723.3-7327"
  - arXiv:2207.07102: Pascale et al., "Unscrambling the lensed galaxies in JWST images behind SMACS0723"
  - arXiv:2207.05007: Golubchik et al., "HST strong-lensing model for the first JWST galaxy cluster SMACS J0723.3-7327"
  - arXiv:2301.03629: Diego et al., "On the correlation between dark matter, intracluster light and globular cluster distribution in SMACS0723"
  - arXiv:2212.07366: Noirot et al., "The first large catalogue of spectroscopic redshifts in Webb's First Deep Field, SMACS J0723.3-7327"
  - arXiv:2605.03442: Dima et al., "High-Redshift Gravitational Lens Discoveries in JWST NIRCam Using AnomalyMatch"
  - arXiv:2503.15326: "Euclid Quick Data Release (Q1). The Strong Lensing Discovery Engine C: Finding lenses with machine learning"
  - arXiv:2110.00023: Stein et al., "Mining for Strong Gravitational Lenses with Self-supervised Learning"
  - arXiv:1703.02642: Lanusse et al., "CMU DeepLens: Deep Learning For Automatic Image-based Galaxy-Galaxy Strong Lens Finding"
- Representations and anomaly detection:
  - arXiv:2404.02973: Zoobot, "Scaling Laws for Galaxy Images"
  - arXiv:2310.03024: AstroCLIP
  - arXiv:2405.14930: AstroPT
  - arXiv:2510.17960: AION-1
  - arXiv:2412.02527: Multimodal Universe
  - arXiv:2304.07193: DINOv2
  - arXiv:2508.10104: DINOv3
  - arXiv:2503.21869: Galaxy Zoo JWST
  - arXiv:2509.19453: "The Platonic Universe"
  - arXiv:2010.11202: Astronomaly
  - arXiv:2411.04188: Astronomaly Protege
  - arXiv:1901.01588: PyOD
  - arXiv:2202.08341: anomalib
  - arXiv:2410.17142: coniferest
  - arXiv:2505.03509: AnomalyMatch
- Catalogs:
  - arXiv:2302.10936: DJA mosaics, Valentino et al.
  - arXiv:2510.08743: CEERS catalog
  - arXiv:2211.02495: CEERS Epoch 1 reduction
  - arXiv:2601.15956: JADES DR5
  - arXiv:2506.03243: COSMOS2025
  - arXiv:2301.02671: UNCOVER
  - arXiv:2506.21685: CANUCS DR1
  - arXiv:2301.02179: GLASS-JWST

**Documentation**
- JWST operations build table (build 13.0 = jwst 3.0.0 = `jwst_1584.pmap`, installed 2026-09-08):
  https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information
- JDox known issues for NIRCam and MIRI (pages linked in docs/landscape.md §5b): https://jwst-docs.stsci.edu/known-issues/nircam-known-issues
- MAST note on JWST reprocessing and the `CAL_VER`/`CRDS_CTX` headers: https://outerspace.stsci.edu/display/MASTDOCS/Updates+to+JWST+Data
