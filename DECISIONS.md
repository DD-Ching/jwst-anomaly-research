# Decisions

Format: **Decision** / **Alternatives rejected** (and why) / **Evidence** / **Revisit if**.
Read the headings before researching a tool; if a decision exists, reuse it unless its
"Revisit if" condition holds. Sections D-002 to D-010 are owned by bootstrap units 1–9.

## D-001 Repository conventions and catalog-first first slice (2026-10-07)

**Decision.**
- Python package `jwst_anomaly` (src layout, hatchling, uv environments), one module per pipeline stage.
  Stages exchange `astropy.table.Table` with column contracts in `schema.py` and a mandatory
  provenance label in `meta`.
- Raw and derived data stay out of git. Tracked manifests (data URI, sha256, size, pipeline version)
  plus download scripts are the reproducibility record.
- The first vertical slice ranks sources from JWST pipeline level-3 catalogs (`_cat.ecsv`) of
  program 2736 (SMACS 0723, lensing cluster) and uses CEERS (program 1345) as a non-cluster control.
  Image cutouts are fetched only for top-ranked sources, via S3 byte-range reads.
- GitHub flow: public repo, CI on every PR, owner merges all PRs (see CLAUDE.md).
- The agent harness is built from Claude Code features (CLAUDE.md, project skills, `/loop`,
  `/schedule` routines, `/batch`, `claude-code-action`) rather than custom orchestration code.

**Alternatives rejected.**
- pandas DataFrames at stage boundaries: they lose units and metadata. Converting at the ML edge is cheap.
- FITS in git or Git LFS: size and cost. Manifests make it unnecessary.
- Starting from level-2 `_cal` or raw ramps: the calibration work is unneeded for a first ranking.
- Running our own source extraction first: the pipeline already runs photutils. Revisit it for forced photometry.
- A custom scheduler or daemon for long runs: that would rebuild `/loop` and `/schedule`.

**Evidence.** MAST product sizes measured 2026-10-07 (catalog ~3 MB vs i2d ~1.8 GB; SOURCES.md "Data").

**Revisit if.** Cross-band colors from independently detected catalogs prove too noisy (then use a
forced-photometry HLSP catalog or photutils on i2d mosaics), or the sample grows beyond memory
(then use Parquet partitions or a database).

## D-002 Archive access and acquisition (unit 1)

_Open._

## D-003 Catalog ingestion and cross-band matching (unit 2)

_Open._

## D-004 Baseline features and anomaly ranking (unit 3)

_Open._

## D-005 Image cutouts and visualization (unit 4)

_Open._

## D-006 External catalog cross-checking (unit 5)

_Open._

## D-007 Candidate store and run provenance (unit 6)

_Open._

## D-008 Open-source project tooling (unit 7)

_Open._

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9, 2026-10-07)

The full survey, with versions, licenses and per-option verdicts, is in
[docs/landscape.md](docs/landscape.md).

**Decision.**
- **M2 photometry.** Start from the DJA grizli v7.4 catalogs, which cover both SMACS 0723 and the
  CEERS t021 data. They use matched apertures (not PSF-homogenised) and come with eazy photo-z.
  Validate CEERS colours against the PSF-matched CEERS DR1.0 catalog. Use photutils forced
  photometry (`SourceCatalog(detection_catalog=...)`) only where those colours prove inadequate.
  Don't make our own mosaics. The `ceers-full` DJA files are 250–400 MB, so state the reason
  before downloading them (CLAUDE.md).
- **M2 representations.** Start with Zoobot encoders loaded through `timm` (Apache-2.0 weights).
  Use DINOv2 as the generic-vision control. Use the Multimodal Universe JWST subset as unlabelled
  JWST data, but exclude CEERS and any other field being scored from it. Any embedding must beat
  the sklearn baseline on injection-recovery.
- **M2 detectors.** Keep the M1 sklearn baseline. Add PyOD when more detectors are needed, and
  coniferest for active learning.
- **M3 lensing.**
  - Test candidates against published SMACS 0723 models (Mahler+2022, CC0; RELICS HLSP; Caminha+2022)
    and the Noirot+2023 redshifts. Don't refit the cluster.
  - Lens finding: AnomalyMatch first, then a fine-tuned Zoobot.
  - Galaxy-scale modelling: lenstronomy, with PyAutoLens as a cross-check.
- **Calibration.** Use MAST products. Re-run `jwst` only when a product's `CAL_VER`/`CRDS_CTX` lags
  the operational build, when a top candidate needs a correction that is off by default (1/f,
  wisps, persistence, saturated cores), or when a comparison needs identical catalog settings.
  Record `CAL_VER`/`CRDS_CTX` in the manifests.
- **Tracking.** git + manifests now. Add MLflow (local backend only) at M2 if the run count
  outgrows a results table.

**Alternatives rejected.**
- **AstroCLIP, AstroPT, AION-1:** they are tied to Legacy Survey/HSC bands and pixel scales.
  AstroPT also has AGPL code and CC-BY-SA weights.
- **DINOv3:** gated download and a custom license.
- **anomalib:** built for industrial, normal-only training with pixel-level defect maps.
- **SEP:** duplicates photutils.
- **SourceXtractor++ at M2:** conda-only and heavy. Reconsider at M3 for blended arcs.
- **DVC, DataLad:** our inputs already live in a permanent public archive.
- **Refitting cluster lens models:** published models exist.
- **ssl-legacysurvey, CMU DeepLens:** ground-based training, or no weights.
- **ceers-nircam scripts:** no license, and stale.
- **lychee, markdown-link-check, LinkChecker** for `scripts/check_links.py`:
  - lychee is Rust, Apache-2.0, v0.24.2.
  - markdown-link-check is Node, ISC, v3.15.0.
  - LinkChecker is Python, GPL-2.0, v10.6.0.

  The checker had to be stdlib-only and run before any environment exists, without adding a
  Node/Rust toolchain or a GPL dependency. If CI link checking is wanted later, lychee-action is
  the upgrade path.

**Evidence.** docs/landscape.md (checked 2026-10-07; versions from the GitHub/PyPI APIs; arXiv IDs
checked against their abstract pages), DJA file listing
https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html, JWST build table
https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information,
and the MAST `s_region` of `jw01345-o001_t021_nircam_clear-f200w` (queried with astroquery).

**Revisit if.**
- A JWST-trained encoder or a lens finder with public weights appears.
- DJA colours disagree with CEERS DR1.0 beyond the photometric errors.
- MAST reprocesses program 2736 with build 13.0 (then redo the cross-program comparisons).
- A study needs more than ~50 tracked runs (MLflow).
- Derived datasets of several GB must be shared (DVC/DataLad).
- docs/landscape.md is more than 6 months old when a milestone starts.
- The link checker needs more than small fixes (then switch to lychee).
