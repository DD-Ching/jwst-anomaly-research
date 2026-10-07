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

Checked 2026-10-07; versions are those installed and used (see D-005).

- **astropy** 8.0.1 — cloud FITS subsets (`use_fsspec`, `.section`, `Cutout2D`):
  https://docs.astropy.org/en/stable/io/fits/usage/cloud.html ; `Cutout2D`:
  https://docs.astropy.org/en/stable/api/astropy.nddata.Cutout2D.html ; asinh stretch and intervals:
  https://docs.astropy.org/en/stable/visualization/normalization.html
- **fsspec** 2026.9.0 (read-ahead cache, block size): https://filesystem-spec.readthedocs.io/en/latest/api.html ;
  **s3fs** 2026.9.0 (default 50 MiB blocks): https://s3fs.readthedocs.io/en/latest/
- **astrocut** 1.3.0 (evaluated, not used): https://astrocut.readthedocs.io/en/latest/astrocut/index.html ;
  https://pypi.org/project/astrocut/1.3.0/ ; https://github.com/spacetelescope/astrocut
- **matplotlib** 3.11.2, `Figure` without pyplot:
  https://matplotlib.org/stable/gallery/user_interfaces/web_application_server_sgskip.html
- **JWST pipeline resample** (i2d `fillval`, WHT `weight_type`):
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html ; i2d layout (SCI, ERR,
  CON, WHT, VAR_*): https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/science_products.html
- **NIRCam filters** (pupil-wheel filters F162M, F164N, F323N, F405N, F466N, F470N pair with a
  filter-wheel filter, so `PUPIL` names the band):
  https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters
- **MAST download endpoint** `https://mast.stsci.edu/api/v0.1/Download/file?uri=<mast URI>` answers
  HTTP range requests (206 Partial Content on the MIRI F770W i2d).
- E2E inputs (program 2736): `jw02736-o001_t001_nircam_clear-f200w_cat.ecsv` sha256
  `7e7b760110fe778df6b4f143e68846113df20b7c681cbbed9dac3e9064652ba3`;
  `jw02736-o002_t001_miri_f770w_i2d.fits` sha256
  `ff7d4d687941311752f5dbce8ebccd5065e6eafc7eb275a76bcedcc241e5155a`; F200W i2d read in place on S3.

## External catalog services (unit 5)

_Pending._

## Candidate store and run provenance (unit 6)

_Pending._

## Open-source project tooling (unit 7)

_Pending._

## Future milestones: models, frameworks, lensing, HLSPs, tracking (unit 9)

_Pending._
