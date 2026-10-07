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

**Decision.** (2026-10-07)
- Cutouts: `astropy.io.fits` opened on an `fsspec` file (read-ahead cache, 64 KiB block) plus
  `astropy.nddata.Cutout2D` on a lazy `Section`. A ~10-line `Section` subclass reads each cutout's
  rows as one contiguous full-width byte range (one request per extension per cutout) instead of
  one read per row. No new dependency (`fsspec`/`s3fs` are the `cloud` extra). `mast:` URIs go to
  the MAST download endpoint, which serves HTTP range requests.
- Output: SCI + WHT (ERR optional) with the cutout WCS, as small FITS files under
  `outputs/cutouts/<image stem>/`. Pixel size rounded to an odd count (target on the central pixel).
- Quality flags (i2d has no DQ). No-data = non-finite SCI or WHT <= 0. `edge`: no-data on the box
  border or box beyond the array. `nan_center`: no-data within 0.2". `nan`: no-data fraction > 0.1.
  `low_weight`: median core WHT < 0.5 x the image's median positive WHT (32 bin-centre rows; the
  whole WHT if <= 8 MiB). `outside`: no data in the box (row kept, no file). Thresholds are untuned
  heuristics. The table is labeled observed; `meta["column_provenance"]` marks derived columns.
- Display: `astropy.visualization` asinh (a = 0.01) over each panel's 1-99.9 percentile range;
  matplotlib `Figure` (no pyplot). Chosen by eye on 2736 cutouts over PercentileInterval(99.5) with
  a = 0.1 (cores saturated, sky crushed) and a = 0.002 (noise-dominated).

**Alternatives rejected.** Cost of one 96x96 SCI cutout of the 2736 F200W i2d on S3 (1.76 GB):
- `fits.open(uri, use_fsspec=True)` with s3fs defaults (50 MiB read-ahead): 104.9 MB, 25.5 s.
- `astrocut` 1.3.0 `FITSCutout` (STScI; reads S3 FITS): 419.9 MB, 9 requests, 33.1 s by default,
  because it reads every HDU header through 50 MiB blocks; 5.8 MB, 22 requests, 9.6 s with
  `default_block_size=256 KiB`. It also takes one position per instance (file reopened per target)
  and adds asdf, gwcs, spherical-geometry and s3path.
- Plain `Cutout2D(hdu.section)`, 256 KiB read-ahead: 3.9 MB but 15 requests, 7.3 s (one read per
  row); without a cache, 33 s. Chosen row-strip reads: 4.5 MB, 2 requests, 2.4 s.
- Concurrent per-row ranges (`fsspec` `cat_ranges`) would cut bytes ~100x but need our own FITS
  offset/scaling code; not worth it at top-k scale. Full NIRCam i2d download: disk budget
  (CLAUDE.md).

**Evidence.** Measured 2026-10-07 on the owner's machine (throwaway scripts; output in the
unit-4 PR).
- Module run on 10 bright F200W catalog sources: SCI+WHT cutouts plus the weight sample cost
  84.7 MB in 56 requests, 34 s (4.8% of the file). MIRI F770W via MAST HTTPS range reads: 6.2 MB in
  9 requests, 5.5 s, pixel-identical to cutouts from the full 36 MB download.
- Weight reference on the F200W WHT: 64 rows give 6634; 8 rows give 4690-6695 depending on which
  rows, 16 rows 6151-6668, 32 rows 6630-6636. Hence 32 rows (~3.4 MB, ~12 s once per image).
- The i2d SCI headers carry a TAN FITS WCS (0.0312"/px F200W, 0.1109"/px F770W). Pixel positions
  from it for the catalog `sky_centroid`s equal the catalog `xcentroid/ycentroid` to < 1e-4 px. The
  light-weighted centroid is within 1.1 px of the cutout centre for 11 of 12 cutouts (2.2 px for
  the extended source, label 1845).
- MIRI F770W i2d (CAL_VER 2.0.1): SCI is NaN exactly where WHT == 0 (25.6% of pixels); positive WHT
  is bimodal (1st percentile 542, median 5225), so a relative weight flag is informative. Pipeline
  defaults `fillval` NAN, `weight_type` ivm:
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
- https://docs.astropy.org/en/stable/io/fits/usage/cloud.html ;
  https://astrocut.readthedocs.io/en/latest/astrocut/index.html

**Revisit if.** Runs move next to the data (us-east-1), where larger blocks may be cheaper than
latency; top-k grows to thousands (share strips between nearby targets, or `cat_ranges`);
astrocut gains multi-position FITS cutouts with lighter dependencies; injection tests show the flag
thresholds are badly calibrated; multi-band comparison needs a common north-up grid (`reproject`).

## D-006 External catalog cross-checking (unit 5)

_Open._

## D-007 Candidate store and run provenance (unit 6)

_Open._

## D-008 Open-source project tooling (unit 7)

_Open._

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9)

_Open._
