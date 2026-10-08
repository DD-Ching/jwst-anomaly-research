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
- **Claude Code docs**, checked 2026-10-07 against Claude Code v2.1.292:
  - `/loop` and scheduled tasks: https://code.claude.com/docs/en/scheduled-tasks
  - cloud environments (network levels, default allowlist, GitHub proxy, setup scripts):
    https://code.claude.com/docs/en/cloud-environments
  - Desktop scheduled tasks: https://code.claude.com/docs/en/desktop-scheduled-tasks
  - GitHub Actions: https://code.claude.com/docs/en/github-actions
  - permissions (rule syntax, deny before ask before allow): https://code.claude.com/docs/en/permissions
  - settings: https://code.claude.com/docs/en/settings , with schema https://json.schemastore.org/claude-code-settings.json
  - memory: https://code.claude.com/docs/en/memory
  - commands (`/batch`): https://code.claude.com/docs/en/commands
  - subagents: https://code.claude.com/docs/en/sub-agents
  - headless (`-p`): https://code.claude.com/docs/en/headless
  - `claude setup-token`: https://code.claude.com/docs/en/authentication
- **claude-code-action v1**: tag `v1`, released 2025-08-26, https://github.com/anthropics/claude-code-action/releases/tag/v1 .
  Inputs: https://github.com/anthropics/claude-code-action/blob/main/docs/usage.md . Security model (write-access
  check, PR creation by link, `include_comments_by_actor`): https://github.com/anthropics/claude-code-action/blob/main/docs/security.md
  (checked 2026-10-07). The action's git handling at tag `v1`: it embeds its token in the remote URL and fetches fork
  PRs via `refs/pull` (https://github.com/anthropics/claude-code-action/tree/v1/src/github/operations).
- Workflow actions: `astral-sh/setup-uv@v10.2.0` (https://github.com/astral-sh/setup-uv); `actions/checkout@v7`.
  Validated with actionlint 1.7.12 (https://github.com/rhysd/actionlint), checked 2026-10-07.
- GitHub docs: rulesets, https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets ;
  Actions billing (free for public repositories on standard runners),
  https://docs.github.com/en/billing/concepts/product-billing/github-actions (checked 2026-10-07).
- **Hosts used by the science services**, observed on 2026-10-07 with astroquery 0.4.11 by logging DNS lookups:
  `mast.stsci.edu` (query, product list, download), `stpubdata.s3.amazonaws.com` (S3 byte-range FITS),
  `simbad.cds.unistra.fr`, `vizier.cds.unistra.fr`, `cdsxmatch.u-strasbg.fr`, `ned.ipac.caltech.edu` and
  `gea.esac.esa.int`. These feed the routine allowlist in docs/operations.md.

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

## Vetting (cycle 3, checked 2026-10-07)

- **Mahler+2022 multiple-image catalog** (model entry above; CC0-1.0): `ICLv0/arcs.dat`, 62 images, absolute
  coordinates, pinned at commit `f36a41c365865df252a6da08cd79955f9cd68f11`:
  https://raw.githubusercontent.com/guillaumemahler/SMACS0723-mahler2022/f36a41c365865df252a6da08cd79955f9cd68f11/ICLv0/arcs.dat
  (sha256 `35d568c43423974940275bddffb163c0acd4a95c907c108a278caf6b866c9154`). Model reference centre RA 110.826750, Dec −73.454628 (`ICLv0/best.par`).
- **NED** entries `SMACS J0723:[NDA2023] 01100` (z 1.9807) and `[NDA2023] 00908` (z 1.3618), retrieved through
  NED TAP by `crossmatch.query_matches` on 2026-10-07. The original reference for the `NDA2023` designation and
  the redshift type were not verified.

## Matched-aperture photometry (cycle 4, checked 2026-10-07)

- **DJA v7.4 SMACS 0723 catalog** (grizli; no license stated; cite arXiv:2302.10936, see docs/landscape.md):
  - file https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix_phot.fits, last-modified
    2025-03-05, 38,626,560 bytes, sha256 `361a4eb944fefa3afbbbfe08192d2d5c5764ae7137e41d77c2ce63c5b4a47166`;
  - 13,554 sources;
  - apertures `ASEC_0..2` = 0.36, 0.5 and 0.7″ diameter;
  - fluxes in µJy.
- **SEP flag bits** used to filter DJA aperture measurements (`APER_TRUNC` 0x10, `APER_HASMASKED` 0x20,
  `APER_ALLMASKED` 0x40, `APER_NONPOSITIVE` 0x80): https://raw.githubusercontent.com/kbarbary/sep/main/src/sep.h
- **SIMBAD/NED designations** recovered in the top 20 of runs `20261007T120528Z-378af5a1` and `20261007T122054Z-c52935ec`, as returned by
  `crossmatch` on 2026-10-07 (papers of the first three verified in "Vetting (cycle 7)"):
  - `SMACS J0723-73:[MJR2023] 028.2` (NED G_Lens);
  - `[YML2023] F150DB-C-4` (SIMBAD);
  - `[RBI2023] 18`, `SMACS J0723-7327:[CSM2022] 78`, `[MS2023] WDF-P-6576` (SIMBAD).

## Vetting (cycle 7, checked 2026-10-07 UTC)

- **DJA v7.4 SMACS 0723 eazy photo-z** (same DJA terms as the catalog above):
  - tarball https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix.photoz.tar.gz, last-modified
    2025-03-05, 59,032,200 bytes, sha256 `13b5b6f40048b1f474827d93ca52dc94f0200261187c59ca397d70382c8c2821`;
  - only the member `smacs0723-grizli-v7.4-fix.eazypy.zout.fits` is used (13,052,160 bytes, sha256
    `82db8910962f3050ddb9f8209aafa487fb30f3cf46551d57ef3824d7d2c355ea`, 13,554 rows). Extract it with
    `tarfile` (`filter="data"`) to `$JWST_ANOMALY_DATA/cache/dja/photoz-v7.4/`;
  - header: eazy-py `VERSION` 0.8.3, templates `templates/sfhz/sorted_agn_blue_sfhz_13.param`, `PRIOR` False.
- **Papers behind catalogue designations**, as SIMBAD TAP (`ident`→`has_ref`→`ref`) and NED TAP (`objdir`) return
  them. Titles are verbatim from SIMBAD.
  - 2023ApJ...942L...9Y: "First Batch of z ≈ 11-20 Candidate Objects Revealed by the James Webb Space Telescope
    Early Release Observations on SMACS 0723-73." (`[YML2023]`)
  - 2023ApJ...945...49M: "Precision Modeling of JWST's First Cluster Lens SMACS J0723.3-7327." It is NED's only
    reference for `[MJR2023] 028.2` and the published version of arXiv:2207.07101.
  - 2023MNRAS.518L..19R: "JWST unveils heavily obscured (active and passive) sources up to z ∼ 13." (`[RBI2023]`)
  - 2023ApJS..265....5H: "A Comprehensive Study of Galaxies at z ∼ 9-16 Found in the Early JWST Data: Ultraviolet
    Luminosity Functions and Cosmic Star Formation History at the Pre-reionization Epoch."
  - 2023ApJ...947L...1Y: "Pointlike Sources among z > 11 Galaxy Candidates: Contaminants due to Supernovae at High
    Redshifts?"
  - 2023MNRAS.518.4755A: "Discovery and properties of ultra-high redshift galaxies (9 < z < 12) in the JWST ERO
    SMACS 0723 Field."
  - 2023MNRAS.525.2087B: "High-z galaxies with JWST and local analogues - it is not only star formation." It is the
    NED redshift reference (z 2.7412, flag `SLS`) for `[MS2023] WDF-C-2420`. The `[MS2023]` designation paper
    itself is not verified.
- **NED redshift flags** (first letter: S = spectroscopic, P = photometric; `SLS` = spectroscopic, several lines,
  secure): https://ned.ipac.caltech.edu/Documents/Guides/Database

## Second epochs (cycle 9, checked 2026-10-07 UTC)

- **JWST program 6882**, "Vast Exploration for Nascent, Unexplored Sources (VENUS)", PI Fujimoto (MAST metadata).
  Observation `o057_t057` covers part of SMACS 0723 on 2026-06-05 (t_min 21:55:17 UTC). Program 2736's F444W
  epoch is 2022-06-07 (t_min 05:18:42 UTC). Both are public, level 3. The VENUS catalogs are jwst 3.0.0 and
  photutils 3.0.0; sha256 values follow.
  - `jw06882-o057_t057_nircam_clear-f150w_cat.ecsv` `d71e54f235b5acb0ee7fc21acd7527bf8fdefe91eae5c771c53afa9a42f5a4a8`
  - `jw06882-o057_t057_nircam_clear-f182m_cat.ecsv` `8f76b50430d8ec90b8cb7afd268f72be6ee5367743f71e5446190cbfd761226b`
  - `jw06882-o057_t057_nircam_clear-f210m_cat.ecsv` `92b5a1c661d3e361497721071a0f6a09ebd6894f09073e96785261b04ee3d7a3`
  - `jw06882-o057_t057_nircam_clear-f300m_cat.ecsv` `f560d0f3ac43b014e7e34d47c94ba96f393e8cd8507ac1447574a9f91100fb0c`
  - `jw06882-o057_t057_nircam_clear-f410m_cat.ecsv` `876f469e51a51ae34d324a31d6b8ad1b695364a804685c1e420f2acff8782deb`
  - `jw06882-o057_t057_nircam_clear-f444w_cat.ecsv` `e431813bdaf255fbdc0cb4003e96051f7c8babfddc8223199c27582ad43efca6`
- **Little-red-dot colour criteria:** Kokorev et al., "A Census of Photometrically Selected Little Red Dots at
  4 < z < 9 in JWST Blank Fields", arXiv:2401.09981 (§3.1: `red1`/`red2` colour cuts, F444W compactness
  f(0.4″)/f(0.2″) < 1.7, F444W > 14σ and < 27.7 mag, brown-dwarf removal F115W−F200W > −0.5, z16 > 4). Read
  2026-10-07.
- **Brown-dwarf NIRCam selection:** Hainline et al. 2024, ApJ 975, 31, "Brown Dwarf Candidates in the JADES and CEERS
  Extragalactic Surveys". It selects on blue 1–2.5 µm and red 3–4.5 µm colours (abstract).
  https://experts.arizona.edu/en/publications/brown-dwarf-candidates-in-the-jades-and-ceers-extragalactic-surve/
  The colour thresholds themselves were not read.
- Langeroodi & Hjorth, arXiv:2308.10900, "Little Red Dots or Brown Dwarfs? NIRSpec Discovery of Three Distant
  Brown Dwarfs Masquerading as NIRCam-Selected Highly-Reddened AGNs": the two populations overlap in NIRCam colours.
  Only the title is read.

## Diffraction-spike detection survey (cycle 10, checked 2026-10-08)

All of these were rejected for D-018; they are listed so later sessions do not survey them again.
- STPSF 2.2.0 (2025-12-23, BSD-3-Clause; requires `numpy<2.4.0`): https://github.com/spacetelescope/stpsf. Data
  files (88 MB): https://stpsf.readthedocs.io/en/latest/installation.html
- MaxiMask 1.4.1 (2024-03-22, MIT): https://github.com/mpaillassa/MaxiMask. Paper: Paillassa, Bertin & Bouy,
  arXiv:1907.08298.
- grizli 1.14.2 `mask_IR_psf_spikes` (WFC3/IR only):
  https://grizli.readthedocs.io/en/latest/api/grizli.pipeline.auto_script.mask_IR_psf_spikes.html
- LSST pipe_tasks `DiffractionSpikeMaskTask`: https://github.com/lsst/pipe_tasks/pull/1178/files
- spike-psf: https://spike-psf.readthedocs.io/en/latest/psf.html
- JWST1PASS mask example (J. Anderson, STScI; no license stated):
  https://www.stsci.edu/~jayander/JWST1PASS/CODE/MASK_EXAMPLE/PLEASE_README.txt
- JWST spike geometry (six spikes plus two fainter ones from the secondary-mirror support):
  https://en.wikipedia.org/wiki/Diffraction_spike
- **Quiescent-galaxy atlas:** 2023ApJ...947...20V, "An Atlas of Color-selected Quiescent Galaxies at z > 3 in Public
  JWST Fields" (SIMBAD designation `[VBG2023]`; title from SIMBAD TAP, 2026-10-07).
- **JWST program 4043**, "Unveiling the build-up of large scale structure in the early Universe", PI Witten (MAST
  metadata). It has NIRCam F444W grism spectra (`jw04043-o001_t001_nircam_f444w-grismr` and `-grismc`, public,
  level 3, t_min 2024-05-11) and F090W/F115W/F444W imaging over part of SMACS 0723, including `2915`.

## Control-field matched photometry (cycle 16, checked 2026-10-07 UTC)

- **DJA v7.4 CEERS catalog** (grizli; same terms as the SMACS catalog above):
  - https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/ceers-full-grizli-v7.4-fix_phot.fits
  - last-modified 2024-09-27, 250,545,600 bytes, sha256 `741ec72f761b19bba57a6c1ef0f10934e03933084330faade88bfacb1d761c28`
  - 81,671 rows; apertures `ASEC_0..2` = 0.36, 0.5 and 0.7″ (`APER_1` = 12.5 px, 0.04″/px)
  - 23 bands, including F115W, F150W, F200W, F277W, F356W, F410M and F444W.

## Cluster fields (cycle 17, checked 2026-10-07/08 UTC; details in docs/fields/*.md)

- **Abell 2744.**
  - Data: MAST 2561 `jw02561-o001_t003_nircam_clear-*` (jwst 3.0.0).
  - DJA v7.2: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix_phot.fits
    (233,527,680 B; sha256 in `configs/abell2744.yaml`).
  - DJA v7.2 eazy photo-z (no tarball for v7.2): https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix.eazypy.zout.fits
    (63,570,240 B, sha256 `436626861c7dbf370419a8787dd79140fad99c2ecbd268257941afe6d493991c`, Last-Modified 2023-12-12,
    accessed 2026-10-08).
  - Lens models:
    - UNCOVER v2.0 maps, https://jwst-uncover.github.io/DR4.html#LensingMaps. Cite Furtak et al. arXiv:2212.04381,
      "UNCOVERing the extended strong lensing structures of Abell 2744 with the deepest JWST imaging", and Price
      et al. arXiv:2408.03920.
    - Bergamini et al. arXiv:2303.10210: Lenstool files at https://www.fe.infn.it/astro/lensing/A2744_Bergamini23/.
      Accessed 2026-10-08: `best.par` sha256 `7245368f96ad9c7159eb9c8d0045030eda0804554312ee86d84f2e052025b9fb`,
      `obs_arcs.cat` sha256 `d02c231f4ee8c81f47335a99182a9f64a4c553e14818b1c9bd07314c2f4f5e1c`.
- **El Gordo.**
  - Data: MAST 1176 `jw01176-o241_t012` (jwst 3.0.0, reprocessed 2026-10-01).
  - DJA v7.0: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix_phot.fits (17,170,560 B,
    sha256 `ac1cf7064edc3f1a907c1b95b7eacd864e19810ff559e7dc39b917081ed66e31`).
  - DJA v7.0 eazy photo-z (D-032): https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix.photoz.tar.gz
    (38,066,605 B, last-modified 2023-08-07, sha256 `9a4a39cee4b154c34040e05591a151597e85d05d0f580cbabdd15519488baade`,
    accessed 2026-10-08). Only the member `elgordo-grizli-v7.0-fix.eazypy.zout.fits` is used (6,546,240 B, sha256
    `efd051e4e688ed46f7279e08cab19ec4179e114705d38f879ffb8feff230edd0`, rows match the `fix_phot` ids). The tarball
    also holds `elgordo-grizli-v7.0-fix_phot_apcorr.fits` (not inspected).
  - Lens models:
    - Caminha et al. 2023 (A&A 678, A3, arXiv:2209.02718): magnification maps, `best_fit.par` and the
      multiple-image list at https://cdsarc.cds.unistra.fr/ftp/J/A+A/678/A3/.
      Accessed 2026-10-08: `best_fit.par` sha256 `7b0153ae0ee02f057f6aaa6f46b1b698502e6fc427266ac9a09d241ddc63a472`,
      `obs_arcs_v1_new_IDs.dat` sha256 `d631743921266c34689a1d509f08e53dc3c90bc88064393d7b8fd524a3d5c700`.
    - RELICS models: https://archive.stsci.edu/hlsps/relics/act0102m49/models/ (DOI 10.17909/T9SP45).
- **Sunrise (WHL0137−08).**
  - Data: MAST 2282 `jw02282-o010_t001` (jwst 2.0.1); later epochs 2282 o120 and 6882 o052.
  - DJA v7.5: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix_phot.fits (28,779,840 B,
    sha256 `8bd178e94156d6a9a126e487dec8a850b6df9f6d76ceb9e36931d9960558ebf0`).
  - Lens models:
    - RELICS: https://archive.stsci.edu/hlsps/relics/whl0137m08/models/. Lenstool v1 maps used (accessed
      2026-10-08, sha256 in `lensmodel.WHL0137_RELICS_LENSTOOL`): x/y-arcsec-deflect, kappa, z06p2-magnif.
    - Scofield et al. arXiv:2504.08879: maps at Zenodo 10.5281/zenodo.15110933.
  - DJA v7.5 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix.photoz.tar.gz
    (51,231,576 B, sha256 `a99b604f7d7418271d94f4665e2e8ab2839f67edf53f28d93e0b9a02cfffc6f1`). Member
    `sunrise-grizli-v7.5-fix.eazypy.zout.fits`, sha256
    `672d4bdf31926209633c1a7feeb436bab3e5b08e847a079521820f3409f548e0`.
  - Earendel: Welch et al. arXiv:2209.14866, "A highly magnified star at redshift 6.2".
- **MAST stale product sizes** after the 2026-10-01 reprocessing (PR #29): `acquire` checks the Download
  service's Content-Length.
## HFF CATS lens models (accessed 2026-10-08; D-035)

- Hubble Frontier Fields lens models, CATS team (Lenstool), HLSP https://archive.stsci.edu/prepds/frontier/lensmodels/.
  Versions: MACS0416 v4.1, MACS1149 v4.1, Abell 370 v4, MACS0717 v4.1, Abell S1063 v4.1, Abell 2744 v4.1. The files
  (x/y-arcsec-deflect, kappa, z02-magnif, arcs.txt, params.txt where present) are pinned by sha256 in
  `lensmodel.HFF_CATS`. The release readmes quote rms 0.72″, 0.63″, n/a, 2.41″, 0.48″ and n/a.

### MACS1149 field run (accessed 2026-10-08; D-037)

- JWST CANUCS program 1208, NIRCam level 3 `jw01208-o008_t004_nircam_clear-*`, 8 bands (F090W–F444W), jwst 3.0.0,
  S3 `s3://stpubdata/jwst/public/jw01208/L3/t/o008/`. Manifests: `data/manifests/macs1149*.ecsv`.
- HFF Sharon v4cor MACS1149 deflection maps (second model for the comparison only; not in `MODELS`):
  `https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/sharon/v4cor/hlsp_frontier_model_macs1149_sharon_v4cor_{x,y}-arcsec-deflect.fits`,
  64,008,000 B each. sha256: x `44f1d21a78cead051c08295d13d7c9189638971cf1e9445bb2caa641a129324f`,
  y `54508deb502d9bcafca9733c3b28f3b4d794b64f11167589b56c1fc37d66a15c`.
- DJA v7 imaging index (dawn-cph GitHub Pages `/dja/imaging/v7/`): no MACS1149 mosaic.
- CANUCS DR1 (DOI 10.17909/18nv-np70) MACS1149 cluster-field photometry and EAzY photo-z:
  https://archive.stsci.edu/hlsps/canucs/dr1/macs1149/clu/hlsp_canucs_jwst-hst_multi_macs1149-clu_multi_v1_photometry-cat.fits.gz,
  29,779,449 B, sha256 `08ab67347f2c3dfe4f743cc1b2a9d4b77a1e40eb66ddeb611648bb296adc0739`. The same directory tree has
  CANUCS lens models (`model/`: deflection, κ, γ, best and samples), not yet used.
### MACS0416 field run (accessed 2026-10-08; D-038)

- JWST CANUCS program 1208, NIRCam level 3 `jw01208-o004_t002`, 8 bands, jwst 3.0.0. Manifests:
  `data/manifests/macs0416*.ecsv`. Also available: PEARLS 1176 o211/o212/o213 `t009`.
- DJA v7 index https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html: no MACS0416 mosaic.
- CANUCS DR1 (DOI 10.17909/18nv-np70) MACS0416 cluster-field photometry and EAzY photo-z:
  https://archive.stsci.edu/hlsps/canucs/dr1/macs0416/clu/hlsp_canucs_jwst-hst_multi_macs0416-clu_multi_v1_photometry-cat.fits.gz,
  31,693,111 B, sha256 `339107a5c5041621d7bed4ecc8b4a51b5148a6913d4e513c3a5e783363f11bff`.

## Exotic-lens literature (checked 2026-10-08; summary in docs/exotic_lensing.md)

- arXiv:1009.6084 Abe, "Gravitational Microlensing by the Ellis Wormhole".
- arXiv:1211.0379 Kitamura, Nakajima & Asada, "Demagnifying gravitational lenses toward hunting a clue of exotic
  matter and energy".
- arXiv:1305.5037 Izumi et al., "Gravitational lensing shear by an exotic lens object with negative convergence or
  negative mass".
- arXiv:astro-ph/9409051 Cramer et al., "Natural Wormholes as Gravitational Lenses".
- arXiv:gr-qc/0105070 Safonova, Torres & Romero, "Microlensing by natural wormholes: theory and simulations".
- arXiv:1711.01730 Asada, "Gravitational lensing by exotic objects".
- arXiv:gr-qc/9907019 Clark, Hiscock & Larson, "Null geodesics in the Alcubierre warp drive spacetime: the view
  from the bridge".
- arXiv:2406.02466 Clough, Dietrich & Khan, gravitational waves from warp-drive collapse.
- arXiv:1303.1301 Takahashi & Asada (ApJL 768, L16): limits on negative-mass objects and Ellis wormholes from the
  SDSS quasar lens search.
- arXiv:astro-ph/9904399 Anchordoqui et al., "In search for natural wormholes".
- Degeneracies:
  - arXiv:1706.10279 Kelly et al., "Extreme magnification of a star at redshift 1.5 by a galaxy-cluster lens";
  - arXiv:astro-ph/9707187 Mao & Schneider (flux-ratio anomalies from substructure).

## Lens-model stage (D-024, checked 2026-10-08)

- **Mahler+2022 SMACS 0723 ICLv2** (model entry above; CC0-1.0), pinned at commit
  `f36a41c365865df252a6da08cd79955f9cd68f11`. URLs and sha256 are in `lensmodel.SMACS0723_MAHLER22_ICLV2`.
  - `ICLv2/best.par`: 149 dPIE potentials; reference RA 110.826989, Dec −73.454723; `Chi2pos` 30.913.
  - `ICLv2/arcs.dat`: 60 images.
  - `ICLv2/input.par`: `sigposArcsec` 0.4427.
  - `ICLv2/tmp_k/0000_k.fits.tar.xz`: 17.1 MB; a 3000 × 3000 px κ map, 0.01334″/px, centred on the reference.
- **PyAutoGalaxy** `autogalaxy==2026.10.7.1` (MIT, https://github.com/PyAutoLabs/PyAutoGalaxy;
  https://pypi.org/project/autogalaxy/), `autogalaxy/profiles/mass/total/dual_pseudo_isothermal_mass.py`: the source
  of the dPIE port.
- Profile references named by that code:
  - Kassiola & Kovner 1993, ApJ 417, 450, doi:10.1086/173325, "Elliptic Mass Distributions versus Elliptic
    Potentials in Gravitational Lenses" (Crossref, checked 2026-10-08);
  - arXiv:0710.5636, Elíasdóttir et al., "Where is the matter in the Merging Cluster Abell 2218?"
