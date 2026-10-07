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

Checked 2026-10-07. Library versions are the ones resolved for the project environment that day.

- **scikit-learn** 1.9.1. `IsolationForest`:
  https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html ;
  `LocalOutlierFactor`: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.LocalOutlierFactor.html ;
  overview: https://scikit-learn.org/stable/modules/outlier_detection.html
- **SciPy** 1.18.1. `median_abs_deviation`:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html ;
  `rankdata`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rankdata.html
- Liu, Ting & Zhou 2008, "Isolation Forest", ICDM, 413–422, doi:10.1109/ICDM.2008.17
- Breunig, Kriegel, Ng & Sander 2000, "LOF", SIGMOD, 93–104, doi:10.1145/342009.335388
- JWST `source_catalog` step: aperture EE defaults of 30/50/70 and the CI star thresholds,
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/arguments.html ; APCORR reference
  file (NIRCam rows: filter, pupil, eefraction, radius),
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/reference_files.html
- **PyOD** 3.6.6 (PyPI, released 2026-09-17; evaluated and not adopted, see D-004):
  https://github.com/yzhao062/pyod . ECOD: Li et al. 2022, arXiv:2201.00382,
  doi:10.1109/TKDE.2022.3159580
- **ADBench**, Han et al. 2022, arXiv:2206.09426 (30 algorithms on 57 datasets)
- **Astronomaly**, Lochner & Bassett 2021, arXiv:2010.11202, doi:10.1016/j.ascom.2021.100481 ;
  https://github.com/MichelleLochner/astronomaly (last push 2026-05-12)

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

Checked 2026-10-07.

- **SQLite** 3.49.1 (bundled with CPython 3.12.10): https://www.sqlite.org/whentouse.html ;
  `PRAGMA user_version`: https://www.sqlite.org/pragma.html#pragma_user_version ;
  JSON functions: https://www.sqlite.org/json1.html
- **Python standard library**: `sqlite3` https://docs.python.org/3/library/sqlite3.html ;
  `importlib.metadata` https://docs.python.org/3/library/importlib.metadata.html ;
  `argparse` https://docs.python.org/3/library/argparse.html
- **Astropy** 8.0.1 table I/O for intermediate tables: ECSV https://docs.astropy.org/en/stable/io/ascii/ecsv.html ;
  Parquet https://docs.astropy.org/en/stable/io/unified_table_parquet.html
- Considered, not used (D-007): DuckDB SQLite extension https://duckdb.org/docs/current/core_extensions/sqlite.html ;
  Datasette https://datasette.io/ ; sqlite-utils https://sqlite-utils.datasette.io/ ;
  Alembic https://alembic.sqlalchemy.org/ ; GitPython https://gitpython.readthedocs.io/ ;
  Click https://click.palletsprojects.com/ ; Typer https://typer.tiangolo.com/ ;
  Snakemake https://snakemake.readthedocs.io/

## Open-source project tooling (unit 7)

All links checked 2026-10-07.

- **Contributor Covenant 2.1**: https://www.contributor-covenant.org/version/2/1/code_of_conduct/ ;
  the markdown copied into CODE_OF_CONDUCT.md is
  https://www.contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md
  (sha256 of the downloaded file: `977d781349351fd7c1f076e4c7dc7de2a05b40e12c773542c3815dd4ce7f37ba`.
  The repository copy differs only in the filled-in contact line and the trimmed leading and trailing
  blank lines). The contact line links GitHub "Report abuse":
  https://docs.github.com/en/communities/maintaining-your-safety-on-github/reporting-abuse-or-spam . Source repository:
  https://github.com/EthicalSource/contributor_covenant (release `2.1`). Newer version 3.0, not adopted:
  https://www.contributor-covenant.org/version/3/0/code_of_conduct/
- **Citation File Format 1.2.0**: https://citation-file-format.github.io/ ; schema guide:
  https://github.com/citation-file-format/citation-file-format/blob/1.2.0/schema-guide.md ;
  validator **cffconvert 2.0.0**: https://github.com/citation-file-format/cffconvert . GitHub
  citation support: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files ;
  DOIs through Zenodo: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
- **GitHub issue forms and templates**: syntax:
  https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests/syntax-for-issue-forms ;
  form schema: https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests/syntax-for-githubs-form-schema ;
  `config.yml`: https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests/configuring-issue-templates-for-your-repository ;
  PR template: https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests/creating-a-pull-request-template-for-your-repository
- **GitHub private vulnerability reporting**: enabling it:
  https://docs.github.com/en/code-security/how-tos/report-and-fix-vulnerabilities/configure-vulnerability-reporting/configure-for-a-repository ;
  how reporters use it: https://docs.github.com/en/code-security/how-tos/report-and-fix-vulnerabilities/report-privately
- **Dependabot options reference** (ecosystems `github-actions`, `pre-commit`, `pip`; `groups`, `cooldown`,
  `versioning-strategy`, `labels`): https://docs.github.com/en/code-security/reference/supply-chain-security/dependabot-options-reference ;
  grouping of version vs security updates (`applies-to`):
  https://docs.github.com/en/code-security/tutorials/secure-your-dependencies/optimizing-pr-creation-version-updates
- **GitHub Actions `schedule` event** (delays at the top of the hour; disabled after 60 days without
  activity in public repos): https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule
- **pre-commit 4.6.2**: https://pre-commit.com/ ; hooks pinned in `.pre-commit-config.yaml`:
  **astral-sh/ruff-pre-commit v0.16.10** (https://github.com/astral-sh/ruff-pre-commit/releases/tag/v0.16.10)
  and **pre-commit/pre-commit-hooks v6.0.0** (https://github.com/pre-commit/pre-commit-hooks/releases/tag/v6.0.0).
- Validation tools (not project dependencies): **actionlint 1.7.12** (https://github.com/rhysd/actionlint ,
  PyPI wrapper `actionlint-py` 1.7.12.25) and **check-jsonschema 0.38.2**, which bundles SchemaStore
  schemas for workflows, Dependabot, issue forms and CFF (https://github.com/python-jsonschema/check-jsonschema).
- Considered but not adopted: Renovate (https://docs.renovatebot.com/) and pre-commit.ci (https://pre-commit.ci/).
  D-008 gives the reasons.

## Future milestones: models, frameworks, lensing, HLSPs, tracking (unit 9)

_Pending._
