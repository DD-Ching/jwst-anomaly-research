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
      MCMC chain `bayes.dat` 70,069,007 B, sha256 `bf6ae670…` (accessed 2026-10-08; D-045).
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
      `obs_arcs_v1_new_IDs.dat` sha256 `d631743921266c34689a1d509f08e53dc3c90bc88064393d7b8fd524a3d5c700`,
      `bayes.dat` 1,666,720 B sha256 `2d3f7362…` (MCMC chain, D-045), `to_sample.par` (read for the potfile `mag0`).
      `fits/magnification_best_fit_z2.fits` and `_z4.fits` (4,003,200 B each; sha256 `2cfe1b62…`, `6eab61e2…`; signed μ,
      0.4″/px, model frame), pinned in `lensmodel.ELGORDO_CAMINHA23_MAG_MAPS` for `validate`.
    - Lenstool source, https://git-cral.univ-lyon1.fr/lenstool/lenstool.git (v8.15.6, commit 09cf4cc4, accessed
      2026-10-08): `docs/sphinx/source/section_parfile/image.rst` ("forme"), `src/o_chi.c` (`chi2_img`, σ² = a·b for
      `forme -10`), `src/bayesapp.c` (bayes.dat ln(Lhood) and Chi2). Read for the D-045 amendment.
    - RELICS models: https://archive.stsci.edu/hlsps/relics/act0102m49/models/ (DOI 10.17909/T9SP45).
- **Sunrise (WHL0137−08).**
  - Data: MAST 2282 `jw02282-o010_t001` (jwst 2.0.1); later epochs 2282 o120 and 6882 o052.
  - Level-2 `_cal` exposures (D-039): `s3://stpubdata/jwst/public/jw02282/jw02282010001/`, `.../jw02282120001/`,
    `s3://stpubdata/jwst/public/jw06882/jw06882052001/` (404 imaging files; headers and byte-range cutouts only;
    accessed 2026-10-08).
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

### Abell 370 and Abell S1063 field runs (accessed 2026-10-08; D-043)

- JWST programs 1208 (CANUCS, PI Willott; Abell 370 `o002_t001`) and 3293 (GLIMPSE, PI Atek; S1063 `o001_t001`).
  Manifests: `data/manifests/abell370*.ecsv`, `data/manifests/abells1063*.ecsv`.
- CANUCS DR1 Abell 370 photometry and photo-z:
  https://archive.stsci.edu/hlsps/canucs/dr1/a370/clu/hlsp_canucs_jwst-hst_multi_a370-clu_multi_v1_photometry-cat.fits.gz,
  32,153,881 B, sha256 `f5622f2867aa3094df6861981ee67f7f1879b8381b348224748c7381436ae17b` (13,567 sources, 509 Z_SPEC).
- DJA v7.5 Abell S1063 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abells1063-grizli-v7.5-fix.photoz.tar.gz,
  74,481,914 B (last-modified 2025-01-09), sha256 `088a1954e357d96e9e3b0ab4af1f474eb9d35676774c1ad990f40140e81b1319`.
  The member `abells1063-grizli-v7.5-fix.eazypy.zout.fits` has sha256
  `03bc4c289982dc50a4540305b8858b279aca7c7cdb64c809d69a855e631b6428` (12,597 rows).
- Gaia DR3 (VizieR I/355/gaiadr3): positions and G magnitudes of the spike-seeding stars.

### MACS0717 field run (accessed 2026-10-08; D-041)

- JWST program 6882 (VENUS, PI Fujimoto), observation 29 `jw06882-o029_t063`, 10 NIRCam bands, jwst 3.0.0, released
  2025-11-03. Manifests: `data/manifests/macs0717*.ecsv`.
- HFF MACS0717 second models, https://archive.stsci.edu/pub/hlsp/frontier/macs0717/models/<team>/<version>/. They
  were read by byte range (κ/γ) or solved in scratch code and then deleted (Sharon and Keeton deflection maps). Full
  sha256 values were streamed on 2026-10-08. The deflection-map prefixes match those recorded during the run.
  - `hlsp_frontier_model_macs0717_sharon_v4cor_kappa.fits` 64,008,000 B `4dcb81cb70a5b429687566b850490ee9f727cb909b6c769adbcf03927930511c`
  - `hlsp_frontier_model_macs0717_sharon_v4cor_gamma.fits` 64,008,000 B `077c548ea9b19fe717b9efd4b6f25a274f7c799b73b20dee0acd3ab55b8ccbca`
  - `hlsp_frontier_model_macs0717_sharon_v4cor_x-arcsec-deflect.fits` 64,008,000 B `ded45fa9802561f764fa3c8fd427b8eebb620155ef43760b608bf7a8e719f758`
  - `hlsp_frontier_model_macs0717_sharon_v4cor_y-arcsec-deflect.fits` 64,008,000 B `913acb9272a63776b70c71eb6b3424f94c94f2392572c1e844ef733c3c8c68f8`
  - `hlsp_frontier_model_macs0717_keeton_v4_kappa.fits` 64,039,680 B `2bd093da64ea551bf49ba5026a24690ee3884d544a2b476e8ae0dce6c5b49684`
  - `hlsp_frontier_model_macs0717_keeton_v4_gamma.fits` 64,039,680 B `50f203bd22748da418810dbb8e5eaa1548297082e8010652823ede58503b075e`
  - `hlsp_frontier_model_macs0717_keeton_v4_x-arcsec-deflect.fits` 64,039,680 B `794590824dc56f12da9a261774c9488ab79a57cb792dcfebff586bd6cb078d6c`
  - `hlsp_frontier_model_macs0717_keeton_v4_y-arcsec-deflect.fits` 64,039,680 B `9820902025c448a0f760facd51f107fffc2aa61ea7e6efb519d4a6c240e207c8`
  - `hlsp_frontier_model_macs0717_glafic_v3_kappa.fits` 277,318,080 B `3c0f8b0ed0b9c7a04454c4b59640550479845b88ba47832b4cd539a19ff8319d`
  - `hlsp_frontier_model_macs0717_glafic_v3_gamma.fits` 277,318,080 B `618bd396d6c3ade6ad67f2e593c475b6a99371d320483426f757cdc2142a0bfe`
  - `hlsp_frontier_model_macs0717_williams_v4.1_kappa.fits` 3,090,240 B `cc807a9a391246e911f8e0a54a52e35139360c41949bd678c0983e7c3a5a7a40`
  - `hlsp_frontier_model_macs0717_williams_v4.1_gamma.fits` 3,090,240 B `75ed37026ecedf04b250684669ab7263a3fb20e66754a017fc272dfd94cd87d3`
  - `hlsp_frontier_model_macs0717_diego_v4.1_kappa.fits` 1,056,960 B `e3304e045a04639df725351bdfd12aec2dc530bc2b863b52e3631f927d10e69f`
  - `hlsp_frontier_model_macs0717_diego_v4.1_gamma.fits` 1,056,960 B `316ab0e8953021b4ad0cef8d085aa0bb8c5df1acb5f74bcdee66fcc95bdfe9da`

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
- CANUCS DR1 MACS0416 Lenstool model (Rihtarsic et al. 2025, A&A, doi:10.1051/0004-6361/202451117), files under
  https://archive.stsci.edu/hlsps/canucs/dr1/macs0416/model/:
  - `hlsp_canucs_jwst-hst_multi_macs0416-lenstool-bestparam_multi_v1_model.txt` 70,753 B
    `f3d9a8415044ff8d5ab4774573fcf3670c9dee9169abfca46086bdbf3dd45bed`;
  - `…-lenstool-multim_multi_v1_model.txt` 19,280 B `ce00444dcc9239f0fd72d5fb37e35cbbea281647803ef4c30eb7b187262f5507`;
  - `…-lenstool-readme_multi_v1_model.txt` 2,275 B `0d3f6796f79c9b9710223aa2ef5d3535c27e2da06c7cd89dc47f1ae3f0915cbf`.
  Deflection, κ and γ maps (best fit and 100 MCMC samples) are in the same directory and were not used.
  - `…-lenstool-param_multi_v1_model.txt` (input parameters; sigposArcsec 0.49) 4,468 B
    `0f1fb7d6947d467b28d8b74485799321b09a1cc63fd7e7df82f01620ddb9b337`.
- CANUCS DR1 Abell 370 Lenstool model (Gledhill et al. 2025, ApJ, doi:10.3847/1538-4357/ad684a), files under
  https://archive.stsci.edu/hlsps/canucs/dr1/a370/model/ (accessed 2026-10-08; D-044):
  - `hlsp_canucs_jwst-hst_multi_a370-lenstool-bestparam_multi_v1_model.txt` 90,403 B
    `3c1eea91755ea532e424b1b143e39b35a9d89beae2958d3b042876af50443580`;
  - `…-lenstool-multim_multi_v1_model.txt` 7,881 B `d72c3d98e675e5bc00cfbdd9b84d1b8528b22e36311924fea072293af23ef9e2`;
  - `…-lenstool-param_multi_v1_model.txt` (sigposArcsec 0.3) 6,254 B
    `aacc2dadd442645d2222c23ee2c3f9f6a76fddaa73996a26a6e88c0691a6bf3d`;
  - `…-lenstool-readme_multi_v1_model.txt` 1,988 B `a73bfe3b7d23c7605e16717781d94736d003e7b5fc5d08c0ef80a0358367edc0`.

## Orphan-pair search inputs (accessed 2026-10-08; D-048)

Pinned in `scripts/orphan_pairs.py` `FIELDS` and fetched by `photometry.fetch_catalog` (cache
`data/cache/external/<sha256[:12]>_<name>`):
- CANUCS DR1 photometry catalogues of MACS0416, MACS1149 and Abell 370 (URLs and sha256 above and in
  docs/fields/abell370.md).
- The Abell 370 image list is the pinned `abell370-canucs` `arcs.dat` (D-044).
- MACS0416 all-multiple-image catalogue (CANUCS `allmultim-cat`) 18,806 B
  `c8978003d8dd617cb980ed7ba5acde1485cd742db43846c25ed251f110dc2417`:
  <https://archive.stsci.edu/hlsps/canucs/dr1/macs0416/model/hlsp_canucs_jwst-hst_multi_macs0416-allmultim-cat_multi_v1_model.txt>.
- MACS1149 Lenstool readme `hlsp_canucs_jwst-hst_multi_macs1149-lenstool-readme_multi_v1_model.txt` 2,198 B
  `aef1cfffaf1865a0f8f2927f4df6fe2b3d44ed5e211aa7dc3f57a97d4651cee1` (no image list until v2).
- Photometry readme <https://archive.stsci.edu/hlsps/canucs/dr1/webpage/hlsp_canucs_jwst-hst_multi_v1_photometry-cat_readme.txt>
  (read, not cached).

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

## Exotic-lensing predictions (D-047)

All accessed 2026-10-08. Each ID was fetched from the arXiv API (export.arxiv.org); journal DOIs were checked on
Crossref. The papers marked "full text" were read in their arXiv source for the formulas used in
`src/jwst_anomaly/exotic_sim.py`; the others were checked from the abstract.
- Wormholes and negative masses:
  - arXiv:1211.0379, Kitamura, Nakajima & Asada, PRD 87, 027501 (2013), doi:10.1103/PhysRevD.87.027501 (full
    text): α = ε̄/bⁿ; 2/(n+1) is their leading-order, large-n estimate of the demagnification onset; n = 10
    onset at β = 0.187 (numerical); Fig. 2c (n = 3) minimum A ≈ 0.865, read from the figure pixels.
  - arXiv:1305.5037, Izumi et al., PRD 88, 024049 (2013), doi:10.1103/PhysRevD.88.024049 (full text): λ_±, κ, γ;
    radial images for ε < 0; voids as negative convergence.
  - arXiv:1307.6637, Kitamura et al., "Microlensed image centroid motions by an exotic lens object with negative
    convergence or negative mass", PRD 89, 084020 (2014), doi:10.1103/PhysRevD.89.084020 (full text): θ_E for
    any n and sign; centroid shifts.
  - arXiv:1009.6084, Abe, ApJ 725, 787 (2010), doi:10.1088/0004-637X/725/1/787 (full text): Ellis α = πa²/(4b²),
    R_E, gutters of about 4 %, inner-image values at β = 2 and 3.
  - arXiv:1107.5374, Toki, Kitamura, Asada & Abe, "Astrometric Image Centroid Displacements due to Gravitational
    Microlensing by the Ellis Wormhole", ApJ 740, 121 (2011), doi:10.1088/0004-637X/740/2/121: µas centroid
    shifts.
  - arXiv:gr-qc/0105070, Safonova, Torres & Romero, PRD 65, 023001 (2002), doi:10.1103/PhysRevD.65.023001 (full
    text): caustic at 2θ_E, umbra, A = (u² − 2)/(u√(u² − 4)).
  - arXiv:astro-ph/9409051, Cramer et al., PRD 51, 3117 (1995), doi:10.1103/PhysRevD.51.3117: negative-mass light
    curves differ qualitatively from MACHO ones.
  - arXiv:astro-ph/9802106, Torres, Romero & Anchordoqui, "Might some gamma ray bursts be an observable signature
    of natural wormholes?", PRD 58, 123001 (1998), doi:10.1103/PhysRevD.58.123001: negative-mass density
    ≲ O(10⁻³⁶) g cm⁻³.
  - arXiv:gr-qc/9805075, Torres, Romero & Anchordoqui, "Wormholes, Gamma Ray Bursts and the Amount of Negative Mass
    in the Universe", MPLA 13, 1575 (1998), doi:10.1142/S0217732398001650: an essay version of the same bound.
  - arXiv:1303.1301, Takahashi & Asada, ApJL 768, L16 (2013), doi:10.1088/2041-8205/768/1/L16: SDSS quasar-lens
    limits (see "Existing observational limits").
  - arXiv:1302.7170, Yoo, Harada & Tsukamoto, "Wave Effect in Gravitational Lensing by the Ellis Wormhole", PRD 87,
    084045 (2013), doi:10.1103/PhysRevD.87.084045: n ≲ 10⁻⁹ AU⁻³ for a ~ 1 cm from femtolensing.
  - arXiv:1711.04560, Tsukamoto & Gong, "Extended source effect on microlensing light curves by an Ellis
    wormhole", PRD 97, 084051 (2018), doi:10.1103/PhysRevD.97.084051: an extended source makes the gutter
    shallower.
  - arXiv:gr-qc/0104076, Eiroa, Romero & Torres, "Chromaticity effects in microlensing by wormholes", MPLA 16, 973
    (2001), doi:10.1142/S021773230100398X: finite-source colour signatures (not used numerically).
- Warp drives:
  - arXiv:gr-qc/0009013, Alcubierre, "The warp drive: hyper-fast travel within general relativity", CQG 11, L73
    (1994), doi:10.1088/0264-9381/11/5/001.
  - arXiv:gr-qc/9907019, Clark, Hiscock & Larson, CQG 16, 3965 (1999), doi:10.1088/0264-9381/16/12/313: view from
    inside the bubble.
  - arXiv:1107.5650, Müller & Weiskopf, "Detailed study of null and time-like geodesics in the Alcubierre Warp
    spacetime", GRG 44, 509 (2011), doi:10.1007/s10714-011-1289-0.
  - arXiv:gr-qc/0110086, Natário, "Warp Drive With Zero Expansion", CQG 19, 1157 (2002),
    doi:10.1088/0264-9381/19/6/308.
  - arXiv:2006.07125, Lentz, "Breaking the Warp Barrier: Hyper-Fast Solitons in Einstein-Maxwell-Plasma Theory",
    CQG 38, 075015 (2021), doi:10.1088/1361-6382/abe692.
  - arXiv:2102.06824, Bobrick & Martire, "Introducing Physical Warp Drives", CQG 38, 105009 (2021),
    doi:10.1088/1361-6382/abdf6e.
  - arXiv:2104.06488, Fell & Heisenberg, "Positive Energy Warp Drive from Hidden Geometric Structures", CQG 38,
    155020 (2021), doi:10.1088/1361-6382/ac0e47.
  - arXiv:2406.02466, Clough, Dietrich & Khan, "What no one has seen before: gravitational waveforms from warp
    drive collapse", OJAp 7 (2024), doi:10.33232/001c.121868 (full text): flat exterior; f ~ 300 kHz and
    h ~ 10⁻²¹ at 1 Mpc for a 1 km bubble.
  - arXiv:2205.15950, Schuster, Santiago & Visser, "ADM mass in warp drive spacetimes", GRG 55, 14 (2023),
    doi:10.1007/s10714-022-03061-9.
  - arXiv:1202.5708, McMonigal, Lewis & O'Byrne, "The Alcubierre Warp Drive: On the Matter of Matter", PRD 85,
    064024 (2012), doi:10.1103/PhysRevD.85.064024.
  - arXiv:2608.10800, Fell & Loeb, "Radiative Signatures from Warp Drives Traveling Through the Earth's
    Atmosphere" (2026, preprint).
  - Literature check 2026-10-08 (arXiv API, 2023–2026; INSPIRE citations of 2406.02466, none searches data):
    arXiv:2310.16067, Kuwahara & Cannon, "Development and Application of a Detection System for a Novel Class of
    Gravitational-Wave Transients" (2023): LIGO/Virgo/KAGRA O3 search for GW bursts from superluminal curvature
    sources, null. arXiv:2212.02065, Sellers, Bobrick, Martire et al. (2022): GWs from accelerating massive
    spacecraft (not a warp metric). arXiv:2405.19381, Lentz & Felton, "Motivating Emissions from Positive Energy Warp
    Bubbles" (2024): order-of-magnitude EM fluxes for a bubble 100 lyr away (Eqs. 12–13, Figs. 6–8; no template).
    arXiv:2311.12069, Pieri (2023): no quantitative prediction for a distant observer. The Clough et al. waveform is not public (no data statement or
    Zenodo record found).

## Injection-recovery limits (D-049)

- arXiv:1807.06209, Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters", A&A 641, A6 (2020),
  doi:10.1051/0004-6361/201833910 (checked on the arXiv API 2026-10-08). Used as `astropy.cosmology.Planck18`
  for the θ_E → |M| conversion in `scripts/inject_radial.py`, as in D-047.
- Inputs that were already recorded: the field catalogues, photo-z and models of docs/fields/*.md (SOURCES
  "Cluster fields", "HFF CATS lens models"); the Gaia DR3 stars from `scripts/gaia_stars.py` (D-043); and
  Takahashi & Asada (2013) for the comparison limit ("Exotic-lens literature").

## W1 shear screen (D-050, checked 2026-10-08)
- **TreeCorr** 5.1.4 (2026-09-01, BSD-3-Clause): https://rmjarvis.github.io/TreeCorr/_build/html/ng.html. Evaluated,
  not used.
- **lenspack** 1.0.0 (2020-09-04, MIT): https://github.com/CosmoStat/lenspack. Evaluated, not used.
- **Schneider 1996**, "Detection of (dark) matter concentrations via weak gravitational lensing", MNRAS 283, 837,
  doi:10.1093/mnras/283.3.837, https://arxiv.org/abs/astro-ph/9601039. Catalogue aperture-mass estimator.
- **Schirmer et al. 2007**, "GaBoDS IX. A sample of 158 shear-selected mass concentration candidates", A&A 462, 875,
  doi:10.1051/0004-6361:20065955, https://arxiv.org/abs/astro-ph/0607022. Shear-shaped filter Q_TANH, eqs. 15–16,
  read from the arXiv full text on 2026-10-08 (D-053).
- **Seitz & Schneider**, "Steps towards nonlinear cluster inversion through gravitational distortions III. Including
  a redshift distribution of the sources", https://arxiv.org/abs/astro-ph/9601079 (A&A, 1997), accessed 2026-10-08.
  Reduced-shear inversion ε_s = (ε − g)/(1 − g* ε) used by `exotic_screens.remove_cluster_shear` (D-053).

## Deep-field orphan-pair screen (accessed 2026-10-08; D-051)

- A. L. Read, "Presentation of search results: the CLs technique", J. Phys. G 28, 2693 (2002),
  doi:10.1088/0954-3899/28/10/313 (checked on Crossref 2026-10-08): the background-aware limit in
  `scripts/inject_pairs.py` (`poisson_signal_ul`).

Pinned by URL and sha256 in `scripts/orphan_pairs.py` (`DEEP_FIELDS`); downloaded with `photometry.fetch_catalog`.
- **CANUCS DR1 NIRCam flanking-field (NCF) photometry + EAzY photo-z**, same format and readme as the cluster
  catalogues (https://archive.stsci.edu/hlsps/canucs/dr1/webpage/hlsp_canucs_jwst-hst_multi_v1_photometry-cat_readme.txt,
  sha256 `978413b1052cfc6f951599e468f4083aa11efb4b19528c39a9ebcb46c486511b`). All Last-Modified 2025-04-29.
  URL pattern `https://archive.stsci.edu/hlsps/canucs/dr1/<c>/ncf/hlsp_canucs_jwst-hst_multi_<c>-ncf_multi_v1_photometry-cat.fits.gz`:
  - `macs0416`: 38,095,398 B, sha256 `72a2c609015f830af4ca05bcd49c7754723c0a55434176f780508bfe62a5afe3`, 10,804 rows;
  - `macs1149`: 37,393,732 B, sha256 `f83510e39332657bf149a08d8060b8dd63a9bb9b839ca5d157714c88ea273ceb`, 11,657 rows;
  - `a370`: 37,054,880 B, sha256 `b35c6eed37ae964d7f996a563dca5794ccfa8ccb0e8b364800fc04b54aa3039f`, 10,267 rows;
  - `macs0417`: 23,683,659 B, sha256 `f82940179ebf63ddf2ec872e7a0297f2fce1785faa0de83ee43b401e96acc49d`, 11,027 rows
    (no F356W);
  - `macs1423`: 23,191,315 B, sha256 `b4881e0ea7fc4dd11fbad76afa73a3ca56f0848d0f85384f59646070c9b7b63b`, 10,412 rows
    (no F356W).
  - Cite Sarrouh, Asada et al. 2025 (named in the readme as the catalogue paper; arXiv id not checked here).
- **DJA v7.3 GOODS-North** (grizli; same DJA terms and citation as the other DJA catalogues, arXiv:2302.10936):
  - catalogue https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/gdn-grizli-v7.3-fix_phot.fits, Last-Modified
    2024-02-19, 223,813,440 B, sha256 `9b18b41731c3a86085cb9c4fdb7a4c9f15c5477431f4eea904d410a1f173b6c1`, 70,421 rows;
    apertures `ASEC_0..2` = 0.36, 0.5, 0.7″ on 0.04″ pixels; per-band images not PSF-matched (SW bands on 0.02″
    pixels, header `<BAND>_aper_0` = 18 px);
  - photo-z tarball https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/gdn-grizli-v7.3-fix.photoz.tar.gz,
    Last-Modified 2024-02-19, 371,072,778 B, sha256 `1d89eef3f613eeb7d592ecf03d94869dfd76110c8208eb728291592ddca8db82`;
    only the member `gdn-grizli-v7.3-fix.eazypy.zout.fits` is kept (67,645,440 B, sha256
    `363176053431708410d4957a77e56826f5ced3d6ed072722aeccc9a2edde48fc`, 70,421 rows aligned with the catalogue ids;
    eazy-py `0.6.8.dev0+gf3eee26.d20230928`, templates `templates/sfhz/agn_blue_sfhz_13.param`, `PRIOR` False);
  - **do not use** the standalone https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/gdn-grizli-v7.3-fix.eazypy.zout.fits
    (Last-Modified 2024-02-08, 60,586,560 B, sha256 `3cc7a92465ed76503b1f29c621074162c446c1eb46d8d1c58583ed75a6727c56`):
    63,069 rows from an earlier catalogue; its ids point at other positions (median offset 137″).
  - Other DJA v7 deep-field catalogues checked (HEAD, 2026-10-08), all > 200 MB with ≥ 350 MB photo-z tarballs:
    `gds-grizli-v7.2` 248.3 MB (+352.8 MB), `ceers-full-grizli-v7.4` 250.5 MB (+396.8 MB), `primer-uds-north-grizli-v7.2`
    270.1 MB (+393.0 MB), `primer-cosmos-east-grizli-v7.4` 269.0 MB (+407.9 MB). The gds, ceers and primer-uds-north tarballs
    (first 300 kB read) and the GOODS-N one start with a 112–138 MB `eazypy.h5`; in the GOODS-N tarball the zout is the fifth member.
- **Cutouts**: MAST level-3 `_i2d` F150W/F277W/F444W read by S3 byte range. CANUCS NCF observations of program 1208:
  o023_t002 (MACS0416), o029_t004 (MACS1149), o020_t001 (Abell 370), o026_t003 (MACS0417), o032_t005 (MACS1423);
  GOODS-N: JADES program 1181 (the observation whose MAST footprint contains the pair).

## Multi-epoch NIRCam fields (accessed 2026-10-08; D-052)

Survey: one MAST `Observations.query_criteria` per field (`instrument_name=NIRCAM/IMAGE`, `calib_level=3`,
PUBLIC; cones of 3′ around MACS0416, El Gordo and Abell 2744; `proposal_id` 5105 for NEXUS and 1180/1210/1286/3215
for JADES). Programs, PIs and titles are MAST metadata. The epochs used, with MAST `t_min`, are in
`configs/dimming_screen.yaml`; catalogue sha256 and sizes in `data/manifests/dimming_<field>.ecsv`, level-3 products
with S3 URIs (never downloaded `_i2d.fits`) in `data/manifests/dimming_<field>_products.ecsv`. All catalogues are
jwst 3.0.0 / photutils 3.0.0 (file headers).
- **JWST 5105**, "NEXUS: the North ecliptic pole EXtragalactic Unified Survey", PI Shen: `jw05105-o001/o002/o014/
  o004/o006/o008/o010/o012_t001` (F200W, F444W), 2024-09-12 to 2026-03-28.
- **JWST 1176** (PEARLS, PI Windhorst) `jw01176-o211/o212/o213_t009` and **JWST 1208** (CANUCS, PI Willott)
  `jw01208-o004_t002`, MACS0416, 2022-10-07 to 2023-02-10; **JWST 6882** (VENUS, PI Fujimoto) `jw06882-o054_t054`
  F444W, 2026-01-10.
- **JWST 2561** (UNCOVER, PI Labbé) `jw02561-o001_t003`, `o002_t001`, `o006_t007`, Abell 2744, 2022-11-02 to
  2024-07-31.
- Not used: El Gordo (only 1176 o241 in 2022 and 6882 o051 F444W in 2026: one shared band); JADES (many programs
  and pointings; left for later, TASKS follow-up).
- Gaia DR3 (VizieR I/355/gaiadr3), one 6′ cone per field (`transient_combine.fetch_gaia`), for the star mask.

## Signature layer and time-domain archives (D-054, checked 2026-10-08)

- MulensModel 3.12.0 (2026-09-25, MIT): https://github.com/rpoleski/MulensModel ; docs https://rpoleski.github.io/MulensModel/
- VBMicrolensing 5.6.1 (2026-10-08, LGPL-3.0): https://github.com/valboz/VBMicrolensing
- pyLIMA 1.9.8 (2025-04-04, GPL-3.0-or-later; rejected): https://pylima.readthedocs.io/en/latest/
- OGLE-IV EWS 2011–2026 (publication terms on the page; cite Udalski et al. 2015, Acta Astron. 65, 1):
  https://ogle.astrouw.edu.pl/ogle4/ews/ews.html ; files https://www.astrouw.edu.pl/ogle/ogle4/ews/
- Mróz et al. 2019, OGLE-IV bulge (arXiv:1906.02210, ApJS 244, 29): https://www.astrouw.edu.pl/ogle/ogle4/microlensing_maps/
  (`table3.dat` 5,836 rows; `phot.tar.gz` 48 MB)
- Mróz et al. 2020, OGLE-IV Galactic plane (arXiv:2004.07289, ApJS 249, 16), 630 events:
  https://www.astrouw.edu.pl/ogle/ogle4/galactic_disk_microlensing/
- Gaia DR3 `vari_microlensing`, 363 events (Wyrzykowski et al., arXiv:2206.06121, A&A 674, A23):
  https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_variability_tables/ssec_dm_vari_microlensing.html
- astroquery.gaia epoch photometry (`Gaia.load_data`, astroquery 0.4.11, BSD): https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
- KMTNet event lists 2015–2025 (cite Kim et al. 2016, JKAS 49, 37; proprietary until 1 July of the following year):
  https://kmtnet.kasi.re.kr/~ulens/
- ZTF light-curve API (not used now): https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html
- lenscat `catalog.csv`, 32,838 entries (MIT; Vujeva et al., arXiv:2406.04398): https://github.com/lenscat/lenscat
- Euclid Q1 Strong Lensing Discovery Engine catalogue v0.0.3 (CC-BY-4.0; arXiv:2503.15324): https://doi.org/10.5281/zenodo.15003116
- SuGOHI public lens list (no licence stated): https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/
- Prior-art check (no survey light-curve search for negative-mass or Ellis lenses found): arXiv:1711.04560,
  arXiv:0807.2774, arXiv:1302.7170, arXiv:gr-qc/9805075.

## Published lens catalogues and deep-imaging checks (D-056, accessed 2026-10-08)

- lenscat 1.1.3 `catalog.csv` (MIT; Vujeva et al., arXiv:2406.04398), 32,838 rows, 3,824,227 bytes, sha256
  `7de5111afb6486c119198cb2df868a2e6a15f79869e6cdc8109ec7d59b06a7cc` (identical at HEAD 053719a):
  https://raw.githubusercontent.com/lenscat/lenscat/f531b8a8f3fa4bbd1ee8a58b2ddcd54e6936ea14/lenscat/data/catalog.csv
- Euclid Q1 Strong Lensing Discovery Engine v0.0.3 (CC-BY-4.0; Zenodo record 15025832, concept DOI
  10.5281/zenodo.15003116; Walmsley et al. arXiv:2503.15324, lens models arXiv:2503.15325–15328), files at
  `https://zenodo.org/api/records/15025832/files/<name>/content`:
  - `q1_discovery_engine_lens_catalog.csv`, 2,584 rows, 422,423 bytes, sha256
    `ee5e60cd507413eabf3da8ffa37a3c527212da00feb29083cf516a86d6f26877`;
  - `modeling_lens_mass.csv` (PyAutoLens SIE + shear, 336 lenses), 322,960 bytes, sha256
    `f2a52616a1ac65137b34abe0892d251e18e2c60e9a64c459914102f44ec49db5`;
  - `modeling_sersic_magnitude.csv` (lens VIS/Y/J/H magnitudes), 288,390 bytes, sha256
    `7c7506af33f27d65cc618e998cbbb0868c36a106703f303d6b0b0a74240d94f2`.
  `lens.zip` (3.05 GB), `group.zip`, `unsuccess.zip` were not downloaded.
- SuGOHI public candidate list (HSC-SSP; no licence stated; cite SuGOHI I–X and HOLISMOKES VI, VIII, XIII, XVI as
  listed on the page), 3,961 rows, 423,170 bytes, sha256
  `72fb96dc8d13851c20304b25b8889087d92df27cf1549e6522c8657befda4a55`, served by a PHP script (may change in place):
  https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/download_list.php?file=list_ra_asc_public.csv
- Legacy Surveys DR10 Tractor catalogue `ls_dr10.tractor` via NOIRLab Astro Data Lab TAP (synchronous ADQL; table
  upload and q3c functions are rejected by the ADQL front end on 2026-10-08, so box ORs are batched 300 per query):
  https://datalab.noirlab.edu/tap ; cutouts https://www.legacysurvey.org/viewer/cutout.jpg (layer `ls-dr10`).
  Cite Dey et al. 2019 (AJ 157, 168) and the DR10 acknowledgement at https://www.legacysurvey.org/acknowledgment/
- Cluster catalogues via CDS XMatch: redMaPPer SDSS DR8 v6.3 (Rykoff et al. 2016, VizieR J/ApJS/224/1) and
  Wen, Han & Liu 2012 (VizieR J/ApJS/199/34); SIMBAD (CDS XMatch `simbad`).
- Lemon et al. lensed-quasar database (https://research.ast.cam.ac.uk/lensedquasars/): HTTP 500 on 2026-10-08;
  its entries enter through lenscat.
- Prior art on dark lenses (arXiv API, 2026-10-08): Jackson, Helbig & Browne 1998, "Lensing galaxies: light or
  dark?" (astro-ph/9804136; lens galaxies found in 12 of 12 JVAS/CLASS lenses); Koopmans et al. 2000, CLASS
  B0827+525 "dark lens or binary radio-loud quasar" (astro-ph/0007286); Frey, Paragi & Campbell 2010,
  J1218+2953 (arXiv:1002.1714).
- Spingola et al. 2019, mJIVE-20 VLBI lens search (arXiv:1811.09152; MJV16999 rejected, sect. 4.1.12); SMILE
  milli-lens searches (Casadio et al. 2021, arXiv:2107.06896; Pötzl et al. 2024, arXiv:2409.15229). Cluster-survey
  references in `lenscats.CLUSTER_SURVEY_REFS` were checked by title on the arXiv API (2026-10-08); the two DOI-only
  ones (Lopes et al. 2004 NoSOCS, Gioia et al. 1990 EMSS) by their lenscat names (NSCS, MS cluster designations).
- Legacy Surveys DR10 brick summary `ls_dr10.bricks_s` via Data Lab TAP (query and sha256 of the 2026-10-08 download,
  332,581 bricks / 23.5 MB, recorded in the run's summary.json; not pinned, the service output may change): footprint
  (nexp_r, nexp_z) and per-brick 5σ galaxy depth for D-056. He et al. 2025, lensed-quasar confirmations
  (arXiv:2509.03858; HSC J2212−0103 lens-light fit).

## OGLE-IV microlensing samples (accessed 2026-10-08; D-057)

Cite Mróz et al. as their pages ask; data are pinned in `src/jwst_anomaly/ogle.py` (`FILES`) and
`data/manifests/ogle_mroz.ecsv`, fetched with `photometry.fetch_catalog`. Not OGLE EWS seasons (D-054).

- Mróz, Udalski, Szymański et al. 2019, ApJS 244, 29, arXiv:1906.02210 (v2 e-print read for the selection,
  Table 2, and the efficiency definition): https://www.astrouw.edu.pl/ogle/ogle4/microlensing_maps/
  5,790 events in 112 low-cadence bulge fields (D-054 said 5,836; the nine high-cadence fields come from
  Mróz et al. 2017, Nature 548, 183, and have no light curves here).

  | File | Bytes | sha256 |
  |---|---|---|
  | README | 2,688 | 9e3e62038163881521d5a895b27ba0980edd502c20dcdb65f772ff4c0d9b2136 |
  | table3.dat | 1,542,323 | ca47555840c808967e9a257dba9071acddd499dc078efb5f5c361b3ec8fef67f |
  | table6.dat | 7,355 | a3e27a63b597e49563a061480c0ea37ae475e569900549d04dd8a8afa3ad50da |
  | table7.dat | 10,537 | 56455aecadfea463c9ba623d7fe459a552a1017a2532601d3c610057fa24aa2d |
  | eff.tar.gz | 15,845 | dd5ffa37e4860dfb137691f90e93e11278c16bef705c8950d5f4a5bacc133ac1 |
  | phot.tar.gz | 50,731,904 | 5dafa6835b8456b00eb379d1803a6f6dae46bb5eac54fb8c2f809fe4685bfe6f |

- Mróz, Udalski, Szymański et al. 2020, ApJS 249, 16, arXiv:2004.07289 (v2 e-print read):
  https://www.astrouw.edu.pl/ogle/ogle4/galactic_disk_microlensing/ — 460 events that pass the selection
  (Table B1); the "630" of the abstract adds 170 possible events (Table B2, `data_c/`), not used.

  | File | Bytes | sha256 |
  |---|---|---|
  | README | 2,571 | a456c9f2f0041b81aa1b3dbc21e0be0b7957e29dae7051cc1d88534270b85f33 |
  | table_A1.txt | 139,773 | 60990a6d207171333c8771b9dcc1240217da80ad7f272d20c5faf3432271a954 |
  | table_B1.txt | 131,669 | 065fea96c0f0343f268c675b3132dd36c2c0243ba6c6a1690b784549fbeec410 |
  | eff21.tar.gz | 213,175 | aacd45898288269573cb73ecd936c96fb0992563f65f6752b63701129d11f563 |
  | data.tar.gz | 497,649 | 74868d53863167ebca08cbbbc2df6214433bfc7382f8ee82705529ca206c93c7 |

- Skowron et al. 2016, Acta Astron. 66, 1 (reference list of arXiv:1906.02210): error-bar correction
  already applied to the published photometry.
- Vetting catalogues via CDS XMatch (astroquery 0.4.11): AAVSO VSX `B/vsx/vsx`; Gaia DR3 variability
  classification `I/358/vclassre`. Literature: arXiv API `all:"<event name>"` (export.arxiv.org); the
  per-name queries are rate-limited and returned an error for 6 of 212 bulge names on 2026-10-08
  (recorded as −1 in the vetting record, so they can be re-run).
- Mróz et al. 2017, Nature 548, 183 (reference list of arXiv:1906.02210): the nine high-cadence bulge fields,
  whose events are in the 2019 optical-depth analysis but whose light curves are not in `phot.tar.gz`.
- Hubble Source Catalog v3 (Whitmore et al. 2016, AJ 151, 134), summary `magaper2` via the MAST catalogs API
  `https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/summary/magaper2.csv` (cone search; accessed 2026-10-08; D-060).

## Gaia DR3 microlensing candidates (accessed 2026-10-08; D-061)

- Wyrzykowski, Kruszyńska, Rybicki et al. 2023, "Gaia Data Release 3: Microlensing events from all over the sky",
  A&A 674, A23, doi:10.1051/0004-6361/202243756, arXiv:2206.06121. The v2 e-print source
  (https://arxiv.org/src/2206.06121v2, 2,461,772 bytes, sha256 40c60eee14f1a5e691cf7efb5b1233878fcabd37496eb4a0865fa5e7afa07295)
  was read for the Sample A cuts (Appendix C, Table C.1), the error rescaling (Eq. 9–10) and Table D.1 (Method
  A / B / A+B, parsed by `gaia_mulens.parse_method_table`).
- Gaia archive TAP `SELECT * FROM gaiadr3.vari_microlensing` (363 rows) and `gaiadr3.gaia_source` positions
  (https://gea.esac.esa.int/tap-server/tap/sync); DR3 epoch photometry from the DataLink service
  (https://gea.esac.esa.int/data-server/data, `RETRIEVAL_TYPE=EPOCH_PHOTOMETRY`, INDIVIDUAL CSV). sha256 of every
  cached file: `data/manifests/gaia_dr3_mulens.ecsv`.

## MOA-II 9-year bulge release (accessed 2026-10-08; D-062)

Pinned in `src/jwst_anomaly/moa.py` (`FILES`) and `data/manifests/moa_ii.ecsv`; fetched with
`photometry.fetch_catalog` into `$JWST_ANOMALY_DATA/raw/moa/`. Licence not stated; the archive asks for this
acknowledgement: "This paper makes use of data obtained by the MOA collaboration with the 1.8-metre MOA-II
telescope at the University of Canterbury Mount John Observatory, Lake Tekapo, New Zealand. The MOA
collaboration is supported by JSPS KAKENHI grant and the Royal Society of New Zealand Marsden Fund. These data
are made available using services at the NASA Exoplanet Archive, which is operated by the California Institute
of Technology, under contract with the National Aeronautics and Space Administration under the Exoplanet
Exploration Program."

- NASA Exoplanet Archive, MOA mission page (2006–2014, ~2.4 M light curves, 22 fields, selection text):
  https://exoplanetarchive.ipac.caltech.edu/docs/MOAMission.html ; columns and flux zero point (20 mag =
  691.8 counts on chip 2, 1,445 on other chips): https://exoplanetarchive.ipac.caltech.edu/docs/API_moa_columns.html

  | File | Bytes | Last-Modified | sha256 |
  |---|---|---|---|
  | bulk/metadata.ipac.tar.gz | 97,332,954 | 2023-10-13 | b339ec176933f4bfcad09d6e17f08166c7a18ffea2ddc567b6cca83fa25fe55c |
  | bulk/gb22.tar | 3,510,138,880 | 2023-10-10 | 1cb0173e676915dcf602647e5a2dc315f0ab51c8f2a612ffcdcab4e44edad2d0 |

  `metadata.ipac` has 2,409,061 rows (gb22: 18,599), equal to the Cut-0 count of Nunota et al. 2024. The
  per-object path `data/Contributed/MOA/gb{F}/R/{C}/gb{F}-R-{C}-{S}-{ID}.ipac` (used by the archive viewer;
  undocumented) also answers; not used here.
- Koshimoto, Sumi, Bennett et al. 2023, "Terrestrial and Neptune mass free-floating planet candidates from the
  MOA-II 9-year Galactic Bulge survey", arXiv:2303.08279 (e-print read for Cut-0, Table 2: S/N of SIM > 2.7,
  N_continue,8 ≥ 3, σ_x,y ≤ 1/0.8 px, positive and negative PSF profiles; the archive page still quotes the
  Sumi et al. 2011 cuts S/N > 5, N_detect,continue > 2).
- Nunota, Sumi, Koshimoto et al. 2024, "The Microlensing Event Rate and Optical Depth from MOA-II 9 year
  Survey toward the Galactic Bulge", arXiv:2410.23553 (e-print read): 2,409,061 Cut-0 objects; span
  HJD 2453824–2456970; Table 1 N_s (10 ≤ I_s ≤ 21.4) for 20 fields; gb6 and gb22 excluded (no clear RCG).
- Sumi et al. 2011, Nature 473, 349, doi:10.1038/nature10092, arXiv:1105.3544: the selection the archive page cites.
- Gaia DR3 `gaiadr3.gaia_source` via the ESA TAP service (https://gea.esac.esa.int/tap-server/tap), queried
  2026-10-08: RP histogram within 0.3° of (279.148°, −23.767°) for the luminosity-function slope (D-062).
- Rejected readers (D-062): merida 0.3.2 (https://pypi.org/project/merida/), qusi 1.5.6
  (https://pypi.org/project/qusi/).

## W5 count-deficit screen (accessed 2026-10-08; D-063)

- Legacy Surveys DR10 Tractor `ls_dr10.tractor` via NOIRLab Astro Data Lab TAP (https://datalab.noirlab.edu/tap),
  aggregated server-side per `nest4096` pixel in 2° × 2° chunks (three queries per chunk:
  `countmap.chunk_queries`); regions RA 20–40° and 50–70°, Dec −30° to −20° (DES-wide area). The service output is
  not pinned (it may change); chunk FITS files are cached under `$JWST_ANOMALY_DATA/cache/w5_counts/`. The ADQL
  front end rejected sub-selects, CASE, SIGN and GROUP BY on expressions on 2026-10-08, and returned HTTP 502 for
  ~30 min the same evening. Cite Dey et al. 2019 (AJ 157, 168) and https://www.legacysurvey.org/acknowledgment/
- Galaxy number counts N(< r) of the same selection: `results/w5_counts/numcounts.ecsv` (observed; boxes in its
  meta).
- Legacy Surveys DR10 random catalogues (Myers et al. 2023, arXiv:2208.08518, doi:10.3847/1538-3881/aca5f9;
  ~19.9 GB per file, not used): https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/randoms/
- astropy-healpix 2.0.1 (BSD-3-Clause; https://pypi.org/project/astropy-healpix/) for nest4096 geometry. Rejected:
  healpy 1.20.1 (GPL-2.0, no Windows wheels; https://pypi.org/project/healpy/), healsparse 1.15.0 and hpgeom
  (GPL-3.0-or-later; https://pypi.org/project/healsparse/), Pylians 0.12 (https://pypi.org/project/pylians/).
- Gaia DR3 `gaiadr3.gaia_source` (G < 9) via the ESA TAP service https://gea.esac.esa.int/tap-server/tap, queried
  2026-10-08, for the bright-star vetting test.
- Wen & Han 2024, galaxy clusters from the DESI Legacy Surveys and WISE, ApJS 272, 39 (VizieR J/ApJS/272/39,
  table2; M500 ≥ 3 × 10¹⁴ M☉ and z ≤ 0.6 used), queried 2026-10-08 through astroquery.vizier, for the
  cluster-depletion test.
- HyperLEDA PGC (Paturel et al. 2003, A&A 412, 45; VizieR VII/237/pgc), galaxies with D25 ≥ 1′ (logD25 ≥ 1.0 in
  log 0.1′), queried 2026-10-08, for the large-galaxy (sky over-subtraction) test.
- Safonova, Torres & Romero 2001, "Macrolensing signatures of large-scale violations of the weak energy
  condition", MPLA 16, 153, arXiv:astro-ph/0104075, doi:10.1142/S0217732301003188: the W5 "central void" in the
  background galaxy field.
- Safonova & Torres 2002, "Degeneracy in exotic gravitational lensing", MPLA 17, 1685, arXiv:gr-qc/0208039,
  doi:10.1142/S0217732302008083.
- Broadhurst, Taylor & Peacock 1995, ApJ 438, 49, arXiv:astro-ph/9406052, doi:10.1086/175053, and Umetsu &
  Broadhurst 2008, ApJ 684, 177, arXiv:0712.3441, doi:10.1086/589683: count depletion behind clusters by
  magnification bias (the main ordinary mimic).
- Amendola, Frieman & Waga 1999, MNRAS 309, 465, arXiv:astro-ph/9811458, doi:10.1046/j.1365-8711.1999.02841.x:
  lensing by voids (mimic context).
- Void finders considered and not used: VIDE (Sutter et al. 2015, arXiv:1406.1191); REVOLVER (Nadathur et al. 2019,
  arXiv:1904.01030); DES SV photometric voids (Sánchez et al. 2017, arXiv:1605.03982). DES Y6 Gold (Bechtol et al.
  2025, arXiv:2501.05739) and HSC-SSP PDR3 (Aihara et al. 2022, arXiv:2108.13045) masks not used.

- Rejected lensed-quasar candidates and controls (D-064, `data/manifests/w12_niq_inputs.json` has URLs, bytes, sha256;
  accessed 2026-10-08), VizieR ASU-TSV `https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=<ID>&-out.max=5000&-out.all`:
  Lemon et al. 2023, MNRAS 520, 3305 (arXiv:2206.07714), `J/MNRAS/520/3305/table1`; Inada et al. 2008, AJ 135, 496,
  `J/AJ/135/496/table2,table3`; Inada et al. 2010, AJ 140, 403, `J/AJ/140/403/table2,table3`; Inada et al. 2012,
  AJ 143, 119, `J/AJ/143/119/table3,table4`; Hennawi et al. 2006, AJ 131, 1, `J/AJ/131/1/binqso`.
- MAST HAP cutouts (`astroquery.mast.Hapcut`, https://mast.stsci.edu/hapcut/api/v0.1/astrocut ; accessed 2026-10-09):
  HST ACS/WFC F814W skycell cutouts from program 17308 (J0130+0725, J0728+2607, J2308+3201) and WFC3/UVIS F814W
  from program 17199 (SDSS J1515+1511), used by D-064's archival HST addendum. Not stored (outputs/, gitignored).

## Euclid Q1 MER catalogue (accessed 2026-10-09; D-065)

- Euclid Quick Data Release Q1 MER catalogue, table `euclid_q1_mer_catalogue` on the IRSA TAP service
  https://irsa.ipac.caltech.edu/TAP (sync endpoint `/TAP/sync`), queried 2026-10-09 for counts only
  (`scripts/w5_euclid_feasibility.py`). Use `CONTAINS(POINT, CIRCLE)` for spatial cuts; plain RA/Dec ranges are
  not indexed. The service output is not pinned.
- Same table, rows (shapes) in 0.2–0.3° discs, 2026-10-09 (`scripts/w5_euclid_shear.py`, D-066). The column
  `position_angle` must be quoted in ADQL; its description ("CCW/x") is wrong: it is PA east of north.
- Euclid Q1 MER VIS mosaics (IRSA SIA collection `euclid_DpdMerBksMosaic`; IBE cutouts
  `?center=RA,Dec&size=8arcsec`, gzip-compressed), tile 102022477, accessed 2026-10-09 (`pacheck`).
- Planck PSZ2 (Planck Collaboration 2016, A&A 594, A27; VizieR J/A+A/594/A27) and ACT DR5 clusters (Hilton et al.
  2021, ApJS 253, 3; VizieR J/ApJS/253/3), queried 2026-10-09 for clusters inside Q1 (`CLUSTERS` in the script).
- Q1 MER tile list: IRSA ObsCore (`ivoa.obscore`, `obs_collection = 'euclid_DpdMerBksMosaic'`, VIS), 352 tiles,
  2026-10-09; whole-tile rows fetched by the indexed `tileid` column (`w5_euclid_shear.py fetch`, D-067).
- NFW lensing (R calibration, D-067): Wright & Brainerd 2000, ApJ 534, 34 (arXiv:astro-ph/9908213); c200(M200):
  Duffy et al. 2008, MNRAS 390, L64 (arXiv:0804.2486); astropy `Planck18` cosmology.
- Pantheon+SH0ES distances (Scolnic et al. 2022, ApJ 938, 113, arXiv:2112.03863; Brout et al. 2022, ApJ 938, 110,
  arXiv:2202.04077): `Pantheon+SH0ES.dat` from https://github.com/PantheonPlusSH0ES/DataRelease (branch `main`,
  `Pantheon+_Data/4_DISTANCES_AND_COVAR/`), accessed 2026-10-09 (`scripts/s2_flat_kernel.py`, D-070).
- DES-SN5YR data release (DES Collaboration 2024, ApJL 973, L14, arXiv:2401.02929), Dovekie re-analysis files:
  https://github.com/des-science/DES-SN5YR (branch `main`, HEAD c9a4fca of 2026-01-28),
  `4_DISTANCES_COVMAT/DES-Dovekie_HD.csv`, `DES-Dovekie_Metadata.csv`, `0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz`
  (sha256 in `results/s2_flat_kernel/fit_des_alpha2.json`), accessed 2026-10-09 (D-070 addendum). The Dovekie
  paper reference is not recorded here (not verified this cycle).
- Legacy Surveys DR9 Tractor + DR9 photometric redshifts (Data Lab TAP tables `ls_dr9.tractor`, `ls_dr9.photo_z`,
  joined on `ls_id`; https://www.legacysurvey.org/dr9/), batched box queries, accessed 2026-10-09 (D-070).

## S1 burst twins: Fermi GBM (accessed 2026-10-09; D-071)

- Fermi GBM burst catalogue, HEASARC table `fermigbrst` via TAP https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync
  (ADQL, VOTable), 4,390 rows, 2026-10-09; columns kept and the response sha256 are in
  `results/s1_twins/catalogue.ecsv.gz` (meta). Catalogue papers: von Kienlin et al. 2020
  (arXiv:2002.11460); Poolakkil et al. 2021 (arXiv:2103.13528, doi:10.3847/1538-4357/abf24d).
- GBM burst-catalogue "bcat" files `glg_bcat_all_bn<id>_v<NN>.fit`, HEASARC FTP
  https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts/<YYYY>/bn<id>/current/ , newest version per burst,
  4,389 of 4,390 present (bn number, file name, size and sha256 per row in `results/s1_twins/lc_<YYYY>.ecsv.gz`);
  streamed into memory, never stored. HDU 2 `PHTFLUX`/`PHTFLUXB` used; HDU 1 `PHTCNTS` rejected (D-071).
- GBM localisation systematic: Connaughton et al. 2015 (arXiv:1411.2685), 3.7° (68 %) core plus a
  ~10 % tail to ~14°.
- Prior lensed-GRB search (method and gap): Ahlgren & Larsson 2020 (arXiv:2006.07095).
- Fermi GBM Data Tools (GDT), https://astro-gdt.readthedocs.io/projects/astro-gdt-fermi/en/latest/ (docs 2.2.x,
  opened 2026-10-09): TTE, PHAII, RSP, trigdat, poshist, scat, tcat and catalogue finders; no bcat reader listed.

## E1 causal event network: GBM × ICECAT-1 × GWTC × CHIME/FRB Cat 2 (accessed 2026-10-09; D-TBD)

Inputs are pinned by sha256 in `data/manifests/e1_events.ecsv` (4.7 MB in total; event tables only).

- Fermi GBM burst catalogue, HEASARC `fermigbrst` via TAP https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync
  (columns trigger_name, ra, dec, error_radius, trigger_time, t90, fluence, last_modified), 4,390 rows.
  Catalogue papers as in "S1 burst twins".
- IceCube ICECAT-1 v4, Harvard Dataverse doi:10.7910/DVN/SCRUCD, file `IceCube_Gold_Bronze_Tracks.tab`
  (datafile 7502710, original CSV; the Dataverse md5 cfb7a988cfd2591ba71ab95cd365a3fa matched), 348 tracks
  (340 after dropping `CR_VETO`). Paper: Abbasi et al. 2023, ApJS 269, 25 (arXiv:2304.01174).
- GWOSC cumulative GWTC confident event list,
  https://gwosc.org/api/v2/catalogs/GWTC/events?include-default-parameters=true&format=csv (391 events,
  GWTC-1 to GWTC-5.0). The CSV has no sky positions; sky maps are in the Zenodo PE releases.
- CHIME/FRB Catalog 2, CANFAR doi:10.11570/25.0066, `data/table/chimefrbcat2.csv` (4,057,396 bytes). It
  downloaded at the first attempt on 2026-10-09 (an earlier session's download had been reset by the proxy).
  Paper: arXiv:2601.09399.
- CHIME/FRB Catalog 2 exposure, same DOI, `data/exposure/chimefrbcat2_exposure.h5` (216,024,090 bytes, sha256
  cd8411f92d0ac31bd05dff47f62797638c354444de27a5c056113ca00470d514; `data/manifests/e1_chime_exposure.ecsv`). It
  holds two HEALPix nside-4096 RING maps of time-integrated exposure (s), upper and lower transit, 2018-09-04 to
  2023-09-15, with no time axis. It was reduced to `results/e1_events/chime_exposure_dec_profile.ecsv` and deleted.
  The CANFAR release has no time-resolved uptime file (directories listed 2026-10-09).
- Positive control: Abbott et al. 2017, ApJL 848, L13 (arXiv:1710.05834): GRB 170817A began 1.74 ± 0.05 s
  after the GW170817 merger. SSS17a position: Coulter et al. 2017, Science, doi:10.1126/science.aap9811
  (arXiv:1710.05452).
- GBM instrument: Meegan et al. 2009, ApJ 702, 791 (arXiv:0908.0450). The ~95.6 min orbital period used for the
  orbit-phase null is an ASSUMPTION (a ~96 min low-Earth orbit), not taken from that abstract.
- Prior coincidence searches. Both are same-direction only, and no public code was found on 2026-10-09:
  - Curtin et al. 2023, ApJ 954, 154 (arXiv:2208.00803): CHIME/FRB × GBM/BAT GRBs, 3σ position overlap,
    ≤ 1 week;
  - Masaoka et al. 2026 (arXiv:2603.24983): CHIME Cat 2 × ICECAT-1, best post-trial p = 0.076.
- Time-dependent two-point correlation of a burst catalogue (a same-direction repeater test): Brainerd et al. 1995,
  ApJL (arXiv:astro-ph/9501010, doi:10.1086/187784).
- The antipodal (~176°) BATSE correlation peak is explained by a position-determination bias: Maoz 1994, MNRAS
  269, L1 (arXiv:astro-ph/9308040).

## COSMOGRAIL XIX light curves (accessed 2026-10-09; D-072)
- Millon et al. 2020, A&A 640, A105, arXiv:2002.05736: R-band light curves of 23 lensed quasars, CDS
  J/A+A/640/A105 (https://cdsarc.cds.unistra.fr/ftp/J/A+A/640/A105/, `lcab/*.dat`); delays and redshifts from the
  paper's Tables 1 and 4 (arXiv source `tabdelay.tex`, `tabdata.tex`). `scripts/s3_hybrid.py`.

## D1 distance self-consistency (accessed 2026-10-09; D-073)

- H0LiCOW public distance posteriors, https://github.com/shsuyu/H0LiCOW-public (commit 57cf973, 2025-05-14):
  `h0licow_distance_chains/*`, `MontePython_cosmo_sampling/data/timedelay_6lenses/B1608_Dd_Ddt_params.dat` and the
  lens redshifts in `MontePython_cosmo_sampling/likelihoods/timedelay_6lenses/__init__.py`. Papers: Wong et al.
  2020, MNRAS 498, 1420 (arXiv:1907.04869); Suyu et al. 2010 (B1608 D_dt); Jee et al. 2019, Science 365, 1134
  (B1608 D_d); Chen et al. 2019, MNRAS 490, 1743 (HE0435, RXJ1131, PG1115); Birrer et al. 2019, MNRAS 484, 4726
  (J1206); Rusu et al. 2020 (WFI2033, arXiv:1905.09338).
- TDCOSMO 2025 public release, https://github.com/TDCOSMO/TDCOSMO2025_public (commit d7f38db, 2026-01-21):
  `TDCOSMO_sample/TDCOSMO_data/SDSS1206+4332/final_D_d.npy`, `final_D_dt.npy`, `TDCOSMO_sample/tdcosmo_sample.yaml`.
  Paper: TDCOSMO Collaboration 2025, A&A 704, A63 (arXiv:2506.03023).
- FRBs/FRB repository, https://github.com/FRBs/FRB (commit 996fcda, 2026-05-06, BSD-3):
  `frb/data/Galaxies/public_hosts.csv` and `frb/data/FRBs/FRB*.json` (DM, DMISM). `DMISM` is NE2001 (Cordes &
  Lazio 2002, arXiv:astro-ph/0207156) from `frb/mw.py` `ismDM` (python `ne2001` package, `ElectronDensity().DM(l, b,
  100.)`), set by `frb/builds/build_frbs.py`; used as stored, not recomputed. Macquart et al. 2020, Nature 581,
  391 (arXiv:2005.13161) for ⟨DM_cosmic⟩, p(Δ) and the host log-normal; James et al. 2022, MNRAS 516, 4862
  (arXiv:2208.00819) for F ≈ 0.32.
