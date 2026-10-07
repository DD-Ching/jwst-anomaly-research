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

All links checked 2026-10-07.

- **astroquery** 0.4.11 (used): https://pypi.org/project/astroquery/0.4.11/ ; cloud access section:
  https://astroquery.readthedocs.io/en/latest/mast/mast_obsquery.html#cloud-data-access
- **MAST field definitions**: observations https://mast.stsci.edu/api/v0/_c_a_o_mfields.html ;
  products (`productSubGroupDescription`, `parent_obsid`, `prvversion`, ...) https://mast.stsci.edu/api/v0/_productsfields.html
- **MAST download endpoint**: `https://mast.stsci.edu/api/v0.1/Download/file?uri=<dataURI>`
  (anonymous for PUBLIC data, `Accept-Ranges: bytes`).
- **MAST path lookup** (dataURI → S3 key, JSON): `https://mast.stsci.edu/api/v0.1/path_lookup/?uri=<dataURI>`
- **MAST API tokens** (only for exclusive-access data; read from `$MAST_API_TOKEN`): https://auth.mast.stsci.edu/token
- **MAST public data on AWS**: https://outerspace.stsci.edu/display/MASTDOCS/Public+AWS+Data ;
  registry entry https://registry.opendata.aws/mast-jwst/
- **JWST file naming** (level-3 `<obs_id>_<suffix>` names): https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
- **ECSV** (manifest format): https://docs.astropy.org/en/stable/io/ascii/ecsv.html
- **botocore pin conflict** (D-002): aiobotocore 3.9.2 https://pypi.org/project/aiobotocore/3.9.2/
  (`botocore<1.43.107`) vs boto3 1.43.108 https://pypi.org/project/boto3/1.43.108/ (`botocore>=1.43.108`).
- **pooch** (rejected alternative): https://www.fatiando.org/pooch/latest/

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

_Pending._
