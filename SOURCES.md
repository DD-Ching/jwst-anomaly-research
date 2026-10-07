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

Checked 2026-10-07. Used with astropy 8.0.1.

- **astropy ECSV reader**: https://docs.astropy.org/en/stable/io/ascii/ecsv.html
- **astropy catalog matching**: https://docs.astropy.org/en/stable/coordinates/matchsep.html ;
  `search_around_sky`: https://docs.astropy.org/en/stable/api/astropy.coordinates.search_around_sky.html
- **photutils `SourceCatalog`**, which the pipeline uses to measure the catalogs:
  https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.SourceCatalog.html (docs showed 3.0.0).
  Column meanings are in the JWST `source_catalog` link under "Core software".
- **Program 2736 NIRCam catalogs** `jw02736-o001_t001_nircam_clear-<band>_cat.ecsv`, retrieved 2026-10-07 from
  `https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/<file>`. Made with jwst 2.0.1,
  photutils 2.3.0 and astropy 7.2.0; catalog `date` 2026-08-02. Rows and sha256 prefixes:

  | Band | Rows | sha256 prefix |
  |---|---|---|
  | F090W | 2652 | `c83881377453` |
  | F150W | 2962 | `f5f7c752ace8` |
  | F200W | 3145 | `7e7b760110fe` |
  | F277W | 1860 | `9f7d00b2f1a5` |
  | F356W | 1983 | `d284126a950c` |
  | F444W | 1877 | `bf128e859be8` |

- **Test fixtures** `tests/data/catalog_jw02736-o001_t001_nircam_clear-{f200w,f444w}_cat.ecsv`: the full
  header of the F200W and F444W files above plus only the rows within 5″ of RA 110.685719,
  Dec −73.470661 (10 and 8 rows).
- **CEERS** `jw01345-o001_t021_nircam_clear-f200w_cat.ecsv`: jwst 3.0.0, photutils 3.0.0, 3,695 rows.
  Used for the loader compatibility check.
- **Alternatives considered (D-003)**:
  - STILTS `tmatchn`: https://www.star.bris.ac.uk/~mbt/stilts/sun256/tmatchn-usage.html
  - NWAY: https://github.com/JohannesBuchner/nway
  - CDS XMatch: https://cdsxmatch.u-strasbg.fr/

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
