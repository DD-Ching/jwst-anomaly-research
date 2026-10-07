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

All links checked 2026-10-07. Exact catalog identifiers are as queried by `crossmatch.py` (D-006).

- **SIMBAD** (CDS), main type + all types per object. Wenger et al. 2000, doi:10.1051/aas:2000332.
  - TAP: https://simbad.cds.unistra.fr/simbad/sim-tap ; tables `basic`, `alltypes`, `otypedef`.
  - Object types: https://simbad.cds.unistra.fr/guide/otypes.htx . The `otypedef.path` hierarchy
    was read on 2026-10-07: 226 types; star branch `*`; lensing branch `grv > gLS`.
- **CDS XMatch**: https://cdsxmatch.cds.unistra.fr/ ; docs: https://cdsxmatch.cds.unistra.fr/xmatch/doc/ .
  Targets `cat2='simbad'` and `cat2='vizier:I/355/gaiadr3'`.
- **Gaia DR3**. Gaia Collaboration, Vallenari et al. 2023, doi:10.1051/0004-6361/202243940.
  Positions at epoch J2016.0.
  - VizieR `I/355/gaiadr3`: https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=I/355
  - ESA archive `gaiadr3.gaia_source`: https://gea.esac.esa.int/archive/ ; TAP endpoint
    https://gea.esac.esa.int/tap-server/tap ; data model:
    https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
- **NED** (NASA/IPAC Extragalactic Database).
  - TAP: https://ned.ipac.caltech.edu/tap ; table `NEDTAP.objdir`, described by the service as the
    "NED object directory for N36.1". No upload method; output limit 1,000,000 rows
    (https://ned.ipac.caltech.edu/tap/capabilities).
  - Object types (`G_Lens`, `Q_Lens`, `*` = "Star or Point Source"):
    https://ned.ipac.caltech.edu/help/ui/nearposn-list_objecttypes
- **astroquery** 0.4.11. Ginsburg et al. 2019, doi:10.3847/1538-3881/aafc33. Modules used:
  - XMatch: https://astroquery.readthedocs.io/en/latest/xmatch/xmatch.html
  - SIMBAD: https://astroquery.readthedocs.io/en/latest/simbad/simbad.html
  - NED: https://astroquery.readthedocs.io/en/latest/ipac/ned/ned.html
  - Gaia: https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
- **pyvo** 1.9.1 (TAP client for NED): https://pyvo.readthedocs.io/en/latest/dal/index.html
- **Test fixtures**: `tests/data/xmatch_*.ecsv` are raw responses from the six backends above,
  recorded 2026-10-07 and trimmed to the columns read. Each header records provenance; the
  positions' origins are in `xmatch_targets.ecsv`.

## Candidate store and run provenance (unit 6)

_Pending._

## Open-source project tooling (unit 7)

_Pending._

## Future milestones: models, frameworks, lensing, HLSPs, tracking (unit 9)

_Pending._
