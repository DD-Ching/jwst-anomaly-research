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
- PyAutoLens 2026.10.17.1 (installed in an isolated venv 2026-10-18, `docs/same_source.md`; not a dependency)
- Lenstruction 0.1.0 (MIT; last commit 2022-11-18, pinned to lenstronomy 1.10.3): https://github.com/ylilan/lenstruction
  (cloned 2026-10-18; Yang et al. 2020, MNRAS 498, 3241)
- glafic v2.1.15: https://github.com/oguri/glafic2
- Lens finders (rejected; kept for the record):
  - ssl-legacysurvey: https://github.com/georgestein/ssl-legacysurvey
  - CMU DeepLens: https://github.com/McWilliamsCenter/CMUDeepLens
- jwst 3.0.0: https://github.com/spacetelescope/jwst
- CRDS 14.0.2: https://github.com/spacetelescope/crds
- snowblind 0.2.1: https://github.com/mpi-astronomy/snowblind
- grizli 1.14.2: https://github.com/gbrammer/grizli
- grizli 1.14.3 source (sdist from PyPI, accessed 2026-10-23): `pipeline/photoz.py` (total flux = `{band}_flux_aper_k × flux_auto / flux_aper_k × {band}_tot_corr`), `prep.get_kron_tot_corr`; https://pypi.org/project/grizli/1.14.3/
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
- UNCOVER DR3 photometry (Weaver et al. 2024, ApJS 270, 7): https://jwst-uncover.github.io/DR3.html; LW SUPER catalogue
  `UNCOVER_DR3_LW_SUPER_catalog.fits` (Drive id 19Jd3Xg5BTrHX5icb81txSn_P0LqxmLsz, 83,684,160 bytes), segmentation map
  (1s7DnoX8G73sOLiMiDttlsyZLjwSB2WlX), bCG-subtracted F200W mosaic v5.2.0 (1CPjv4mJUj-4Qn8nwWS0FEPEdBfhDxzV_,
  1.84 GB gzip); headers read 2026-10-23. Pinned 2026-10-23 (stage 2a): SUPER catalogue sha256
  b8f265311d041b815e1511d74eadec2d88bf8c80ef2fee2394d7f78dbe4adb4c; `UNCOVER_DR3_SEGMAP.fits.gz` sha256
  5a64d01fcd75ecbea38efcc4ec673d2cb5fb3fb19b76ced21417fed6bc00087a; `README_DR3.txt` (Drive id
  1Xbfkl2ZxIm-XnGj1ResBSKWCDjkgULqc) sha256 93cfa19494b077a60541a34fc17db6125c5667f3a17a9dcd2c10d2915828c397. The
  mosaic is not pinned.
- Mahler+2022 SMACS 0723 lens model, release v2 (2022-12-27): https://github.com/guillaumemahler/SMACS0723-mahler2022
- Noirot+2023 SMACS 0723 redshifts: https://niriss.github.io/smacs0723
- Caminha+2022 SMACS 0723 model and MUSE redshifts, CDS catalogue J/A+A/666/L9: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/666/L9.
  `muse.dat` (table A1, 78 rows, 4,914 B, sha256 `81868f29…`, accessed 2026-10-14) is pinned in
  `scripts/offplane_potentials.py` `MUSE` for `smacs0723-iclv2`; the model files were not read.
- RELICS SMACS 0723 lens-model files (`glafic/v2`, `lenstool/v1`, `zitrin-ltm-gauss/v2`): https://archive.stsci.edu/hlsps/relics/smacs0723m73/models/
  - Lenstool v1 best maps sit directly under `lenstool/` (not `lenstool/v1/`); x / y arcsec deflection sha256 pinned in `scripts/smacs_relics_heldout.py` (accessed 2026-10-25). zitrin-ltm-gauss v2 readme: Golubchik+22, "JWST-blind".
  - Golubchik et al. 2022, arXiv:2207.05007, Table 1 (17 HST-era image rows of five families; PDF read 2026-10-25): `results/lens_consistency/smacs_relics_heldout/golubchik22_images.ecsv`.
  - Used 2026-10-14 (`scripts/smacs_s17_vet.py`): `glafic/v2` and `zitrin-ltm-gauss/v2` `{x,y}-arcsec-deflect.fits`, sha256 in the
    script (`lenstool/v1` lists no files); z_lens 0.39 and H0 70 / Ωm 0.3 assumed for the D_LS/D_S = 1 maps.
  - arXiv:2207.07102 (Pascale+22) Table 1: system 17 model z 1.72, photo-z 2.34 (accessed 2026-10-14).
- RELICS lens-model readme v3 (`hlsp_relics_hst_multi_lens-models_multi_v3_readme.pdf`, sha256 318f4c03…652d074; every
  deflection map "scaled to D_LS/D_S = 1"): https://archive.stsci.edu/hlsps/relics/docs/ (accessed 2026-10-25).
- RELICS GLAFIC v2 maps of MACS J0553.4−3342 and RXC J2211.7−0349 (Oguri; best fit, D_LS/D_S = 1):
  https://archive.stsci.edu/hlsps/relics/macs0553m33/models/glafic/v2/ and `.../rxc2211m03/models/glafic/v2/`
  (accessed 2026-10-25; sha256 in `lensmodel.MACS0553_RELICS_GLAFIC`, `lensmodel.RXC2211_RELICS_GLAFIC`).
- RELICS MACS J0417 lens models (`glafic/v3` incl. `range/` realisations, `lenstool/v2`; maps scaled to D_LS/D_S = 1;
  readme `hlsp_relics_hst_multi_lens-models_multi_v3_readme.pdf`): https://archive.stsci.edu/hlsps/relics/macs0417m11/models/
  (accessed 2026-10-17; used by `results/transients/batch_venus2/followup/macs0417_h002`).
- RELICS SPT-CL J0615−5746 lens models (`glafic/v2`, `lenstool/v1`; maps at D_LS/D_S = 1) and BPZ catalogue
  `hlsp_relics_hst_acs-wfc3ir_spt0615m57_multi_v1_cat.txt` (Coe, 2017-07-06): https://archive.stsci.edu/hlsps/relics/spt0615m57/
  (accessed 2026-10-18; sha256 in `lensmodel.SPT0615_RELICS_*` and `results/arcs/spt0615-relics/summary.json`).
- Paterno-Mahler et al. 2018, RELICS SPT-CLJ0615−5746 model, z_lens 0.972, ApJ 863, 154, arXiv:1805.09834 (abstract read 2026-10-18; Table 2 image list from the e-print, sha256 83c2bc6c…6cdb26,
  accessed 2026-10-25, `scripts/spt0615_images.py`).
- Adamo et al. 2024, Nature, "Bound star clusters observed in a lensed galaxy 460 Myr after the Big Bang",
  arXiv:2401.03224v2 (HTML read 2026-10-26; Extended Data Table 3 Cosmic Gems clump positions, transcribed in
  `scripts/spt0615_gems.py`).
- Bradley et al. 2025, ApJ 991, 32, doi:10.3847/1538-4357/adf638, arXiv:2404.10770v2 (HTML read 2026-10-26;
  counter-image candidate and JWST-era predicted positions from the text, `scripts/spt0615_gems.py`).
- SN H0pe / PLCK G165.7+67.0 (`results/lens_consistency/h0pe_mu/`; e-prints accessed 2026-10-26, sha256 of the
  arXiv e-print tarballs):
  - Pascale et al. 2025, arXiv:2403.18902v2 (sha256 cdd4d3f2…ca189144): Table 1, model μ and delays (unblinded).
  - Pierel et al. 2024, doi:10.3847/1538-4357/ad3c43, arXiv:2403.18954v2 (sha256 6bb34189…f8d4f6ab): Tables 3–4.
  - Agrawal et al. 2026, ApJ 1002, 187, doi:10.3847/1538-4357/ae5f65, arXiv:2510.07637v2 (sha256 38ada1be…b02bfb19):
    Table 1 and Sect. 4 Δμ offsets.
  - Galan et al. 2026, A&A 708, A359, arXiv:2510.20561v2, and Grayling et al., arXiv:2510.11719v2: abstracts read
    (context: A_V ≈ 0.9 at the SN; BayeSN-TD μ).
- SN Encore / MACS J0138.0−2155 (`results/lens_consistency/encore_mu/`; e-prints accessed 2026-10-26, sha256 of the
  arXiv e-print tarballs):
  - Suyu et al., A&A, doi:10.1051/0004-6361/202557235, arXiv:2509.12319 (e-print of 2026-03-03; sha256
    3cfab090…22266e33ba): Table "tab:snencore", blind model μ and Δt.
  - Pierel et al., arXiv:2509.12301 (sha256 3085083c…26c9e922ff): Tables "tab:final_delay_mag", "tab:systematics",
    "tab:fitting", "tab:h0" and Sect. "Magnifications".
  - Suyu et al. (same e-print) Table "tab:snrequiem": blind model delays of SN Requiem
    (`results/lens_consistency/requiem_dt/`).
  - Rodney et al. 2021, Nature Astronomy 5, 1118, arXiv:2106.08935v2 (sha256 eb10e253…d2ec01c1): Table 1 delays and
    positions.
- SN Refsdal SX check (`results/lens_consistency/refsdal_sx/`, accessed 2026-10-26):
  - Treu et al. 2016, ApJ 817, 60, arXiv:1510.05750 e-print source (sha256 4c2fd6b9…acb094494): Table
    "tab:predictions" (blind Δt_X1, μX/μ1).
  - Kelly et al. 2023, Science 380, eabh1322, arXiv:2305.06367 e-print source (sha256 1e4005f5…a0664b87eef2): Tables
    "tab:measurements", "tab:mod_time_delays", "tab:mod_mag_ratios"; Sharon delay revision in the supplementary text.
- Salmon et al. 2018, SPT0615-JD z ~ 10 arc, arXiv:1801.03103 (abstract read 2026-10-26: RELICS-model μ ~ 4–7,
  predicted counter-images).
- RELICS MACS J0308.9+2645 lens models (`glafic/v2`, `zitrin-ltm-gauss` v1; maps at D_LS/D_S = 1) and BPZ catalogue
  `hlsp_relics_hst_acs-wfc3ir_macs0308p26_multi_v1_cat.txt`: https://archive.stsci.edu/hlsps/relics/macs0308p26/
  (accessed 2026-10-18; sha256 in `lensmodel.MACS0308_RELICS_*` and `results/arcs/macs0308-relics/summary.json`).
- Acebron et al. 2018, RELICS strong-lensing analysis of MACS J0308.9+2645 and PLCK G171.9−40.7, ApJ 858, 42,
  arXiv:1803.00560 (Table 1 z_lens 0.356, read 2026-10-18).
- JWST program 5293 o005 NIRCam F200W level-3 catalogue `jw05293-o005_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC,
  accessed 2026-10-18; sha256 in `results/arcs/macs0308-relics/summary.json`).
- LAM Strong Lensing Cluster Atlas (Lenstool model tarballs and maps; 17 models on 2026-10-18): https://data.lam.fr/sl-cluster-atlas
  (ANIS API `https://data.lam.fr/server/search/sl-cluster-atlas/model?a=1;2;3`; needs cookie `Yogsototh_opens_the_door=1`).
- Mahler et al. 2019 (RELICS MACS J0417 lenstool model, system 1 images): arXiv:1810.13439 (abstract read 2026-10-17: "MACS J0417.5-1154 at z=0.443"; journal reference not checked).
- Meena et al. 2026, "Caught in Swallowtails" (MS 0451.6−0305), arXiv:2601.10097 (abstract read 2026-10-17).
- Palencia et al. 2026, catalogue of microlensed stars in the A370 Dragon arc, A&A 714, A22, arXiv:2604.22702
  (via a literature search 2026-10-17; used for B3-A370-box2-h010).
- Fudamoto et al., Dragon-arc transients, arXiv:2404.08045 (via a literature search 2026-10-17).
- JWST DD 12781 (PI Pascale), "Catching the reappearance of strongly lensed supernova with a 1+ year time-delay"
  (SN Athena, MACS J0417): https://www.stsci.edu/jwst-program-info/program/?program=12781 (proposal PDF read 2026-10-18);
  JWST GO 4568 (target "athena"); HST 18470 and 17474 (MAST programme pages, 2026-10-18).
- JWST GO 9811 (PI Kamieneski), "SN-ECHOES", Cycle 5, includes a PLCK G165 epoch:
  https://www.stsci.edu/jwst-program-info/download/jwst/pdf/9811/ (read 2026-10-18); HST 17504 (PI Kelly), "SNAP Survey for
  Strongly Lensed Supernovae and Magnified Stars" (proposal PDF read 2026-10-18); Legacy Survey DR8 north, VizieR VII/292
  (queried 2026-10-18). Identity checks: `results/transients/identity_r2/identity_check.json`.
- JWST General Target Visibility Tool `jwst_gtvt` 1.2 (PyPI), run 2026-10-18 for MACS0417.
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
  - tarball https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix.photoz.tar.gz (re-checked 2026-10-17, sha256 unchanged), last-modified
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
  files (88 MB): https://stpsf.readthedocs.io/en/latest/installation.html (used 2026-10-17 for `results/w1_stpsf`: stpsf 2.2.0
  from PyPI, data 2.2.0 via the "LATEST" link, https://stsci.box.com/shared/static/kqfolg2bfzqc4mjkgmujo06d3iaymahv.gz)
- jwst `tweakreg/tweakreg_catalog.py` (main branch, read 2026-10-17): source morphology is measured on data convolved
  with a Gaussian of FWHM `kernel_fwhm` (step default 2.0 px, `source_catalog_step.py`):
  https://github.com/spacetelescope/jwst/blob/main/jwst/tweakreg/tweakreg_catalog.py (used for `results/w1_psf_smooth`).
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
    (233,527,680 B; sha256 in `configs/abell2744.yaml`; re-fetched 2026-10-13, sha256 matched, used for D-032's
    Abell 2744 addendum).
  - DJA v7.2 eazy photo-z (no tarball for v7.2): https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix.eazypy.zout.fits
    (63,570,240 B, sha256 `436626861c7dbf370419a8787dd79140fad99c2ecbd268257941afe6d493991c`, Last-Modified 2023-12-12,
    accessed 2026-10-08).
  - Lens models:
    - UNCOVER v2.0 maps, https://jwst-uncover.github.io/DR4.html#LensingMaps. Cite Furtak et al. arXiv:2212.04381,
      "UNCOVERing the extended strong lensing structures of Abell 2744 with the deepest JWST imaging", and Price
      et al. arXiv:2408.03920.
    - Furtak et al. 2025, A&A 698, A227, doi:10.1051/0004-6361/202554110, arXiv:2502.07875v3 (accessed 2026-10-13):
      A2744-QSO1 (z = 7.045), Table 1 image positions = Bergamini+23 system 53; UNCOVER v2 μ A/B/C 5.4/7.5/2.7,
      delays B 855 d, C 7043 d after A; 20 % μ floor, systematics up to ~40 %.
    - Bergamini et al. arXiv:2303.10210: Lenstool files at https://www.fe.infn.it/astro/lensing/A2744_Bergamini23/.
      MCMC chain `bayes.dat` 70,069,007 B, sha256 `bf6ae670…` (accessed 2026-10-08; D-045).
      Accessed 2026-10-08: `best.par` sha256 `7245368f96ad9c7159eb9c8d0045030eda0804554312ee86d84f2e052025b9fb`,
      `obs_arcs.cat` sha256 `d02c231f4ee8c81f47335a99182a9f64a4c553e14818b1c9bd07314c2f4f5e1c`.
- **El Gordo.**
  - Data: MAST 1176 `jw01176-o241_t012` (jwst 3.0.0, reprocessed 2026-10-01).
  - Diego et al. 2023, A&A 672, A3, arXiv:2210.06514 (JWST PEARLS El Gordo lens model): Table `tab_arcs`
    (`TableArcs_NEW3.tex` in the arXiv e-print source, 218 image rows, accessed 2026-10-25), parsed to
    `results/lens_consistency/elgordo_heldout/diego23_images.ecsv`.
  - DJA v7.0: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix_phot.fits (17,170,560 B,
    sha256 `ac1cf7064edc3f1a907c1b95b7eacd864e19810ff559e7dc39b917081ed66e31`).
  - DJA v7.0 eazy photo-z (D-032): https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/elgordo-grizli-v7.0-fix.photoz.tar.gz
    (38,066,605 B, last-modified 2023-08-07, sha256 `9a4a39cee4b154c34040e05591a151597e85d05d0f580cbabdd15519488baade`,
    accessed 2026-10-08; re-checked 2026-10-17, sha256 unchanged). Only the member `elgordo-grizli-v7.0-fix.eazypy.zout.fits` is used (6,546,240 B, sha256
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
      `catalog.dat` (MUSE redshift catalogue, table 1; 26,532 B, sha256 `c4ddb6e2…`, accessed 2026-10-14), pinned in
      `scripts/offplane_potentials.py` `MUSE`.
    - Richard et al. 2021 (A&A 646, A83, arXiv:2009.09784), MUSE lensing-cluster redshift catalogues,
      https://cdsarc.cds.unistra.fr/ftp/J/A+A/646/A83/: `a2744.dat` (284,372 B, sha256 `89d4615f…`, 506 rows, accessed
      2026-10-14), pinned in `scripts/offplane_potentials.py` `MUSE`; `a370.dat` (321,464 B, sha256 `05f524e5…`, 572
      rows, accessed 2026-10-14), pinned there for `abell370-canucs`.
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
      GLAFIC v1 maps used for `arcs` (accessed 2026-10-17, sha256 in `lensmodel.WHL0137_RELICS_GLAFIC`):
      x/y-arcsec-deflect, kappa, z02-magnif.
    - DJA v7.5 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix.photoz.tar.gz
      (51,231,576 B, last-modified 2025-07-11), member `sunrise-grizli-v7.5-fix.eazypy.zout.fits` sha256
      `672d4bdf31926209633c1a7feeb436bab3e5b08e847a079521820f3409f548e0` (accessed 2026-10-17).
    - Scofield et al. arXiv:2504.08879: maps at Zenodo 10.5281/zenodo.15110933 (M5-G1, files dated 250325,
      1,112,414,400 B each; md5 in `scripts/sunrise_scofield.py`; accessed 2026-10-13 by byte range). Image list:
      Table 3 of the arXiv source (v1 e-print, accessed 2026-10-13).
  - DJA v7.5 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix.photoz.tar.gz
    (51,231,576 B, sha256 `a99b604f7d7418271d94f4665e2e8ab2839f67edf53f28d93e0b9a02cfffc6f1`). Member
    `sunrise-grizli-v7.5-fix.eazypy.zout.fits`, sha256
    `672d4bdf31926209633c1a7feeb436bab3e5b08e847a079521820f3409f548e0`.
  - Earendel: Welch et al. arXiv:2209.14866, "A highly magnified star at redshift 6.2".
- **RXC J0600−2007 (MACS J0600.1−2008)** (`arcs` only; accessed 2026-10-18).
  - Data: MAST 6796 `jw06796-o010_t001_nircam_clear-f200w_cat.ecsv` (PUBLIC; sha256 in
    `results/arcs/rxc0600-relics/summary.json`); i2d stamps by S3 byte range from `s3://stpubdata/jwst/public/jw06796/L3/t/o010/`.
  - DJA v7.2 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rxcj0600-grizli-v7.2-fix.photoz.tar.gz
    (33,760,489 B), member `rxcj0600-grizli-v7.2-fix.eazypy.zout.fits`, sha256
    `51d3eece0cfe0d86dc8e2be6c3a48c1fb21c316a58e86ecd651b960a851d93f7`.
  - RELICS models: https://archive.stsci.edu/hlsps/relics/rxc0600m20/models/ (DOI 10.17909/T9SP45): GLAFIC v3 and
    Lenstool v1 x/y-arcsec-deflect, kappa, z02-magnif (sha256 in `lensmodel.RXC0600_RELICS_*`).
  - z_lens 0.43: Fujimoto et al. 2021, arXiv:2101.01937 (ALCS [CII] at z = 6.0719 behind the cluster). Newer model:
    Fujimoto et al. 2024, arXiv:2404.03286 (maps not found public on 2026-10-18).
- **RXC J0018.5+1626 (MACS J0018.5+1626, Cl 0016+16)** (`arcs` only; accessed 2026-10-20).
  - Data: MAST 6882 `jw06882-o001_t001_nircam_clear-f200w_cat.ecsv` (PUBLIC; sha256 in
    `results/arcs/rxc0018-relics/summary.json`); i2d stamps by S3 byte range from `s3://stpubdata/jwst/public/jw06882/L3/t/o001/`.
  - DJA v7.5 photo-z: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rj0018p1626-grizli-v7.5-fix.photoz.tar.gz
    (31,118,096 B), member `rj0018p1626-grizli-v7.5-fix.eazypy.zout.fits`, sha256
    `2a7f9f8eb852dea00f0ff165acc493667d161887f5fee3d255efeca4ac8d604a`.
  - RELICS model: https://archive.stsci.edu/hlsps/relics/rxc0018p16/models/lenstool/v1/ (DOI 10.17909/T9SP45; the
    only model directory): x/y-arcsec-deflect, kappa, z02-magnif (sha256 in `lensmodel.RXC0018_RELICS_LENSTOOL`).
    On 2026-10-20 archive.stsci.edu sent no intermediate certificate (cloud session). Fix used: convert the leaf's AIA
    intermediate (http://repository.emsign.com/certs/EEEMIncommonOVG2C.crt, DER) to PEM, append it to a copy of
    the session CA bundle (`/root/.ccr/ca-bundle.crt`; it held the emSign root on 2026-10-20, not on 2026-10-21) and point `SSL_CERT_FILE` (and
    `curl --cacert`) at the copy; verification stays on.
    2026-10-21 (cloud session): the bundle held no emSign root; `curl --cacert` with the intermediate appended still
    verifies (curl accepts a chain ending at a trusted intermediate), Python urllib does not, so files were fetched
    with curl and their sha256 checked.
  - z_lens 0.546: Zitrin et al. 2011, arXiv:1002.0521 (as MACS J0018.5+1626); a JWST-era model was not searched for.
- **MAST stale product sizes** after the 2026-10-01 reprocessing (PR #29): `acquire` checks the Download
  service's Content-Length.
## HFF CATS lens models (accessed 2026-10-08; D-035)

- Hubble Frontier Fields lens models, CATS team (Lenstool), HLSP https://archive.stsci.edu/prepds/frontier/lensmodels/.
  Versions: MACS0416 v4.1, MACS1149 v4.1, Abell 370 v4, MACS0717 v4.1, Abell S1063 v4.1, Abell 2744 v4.1. The files
  (x/y-arcsec-deflect, kappa, z02-magnif, arcs.txt, params.txt where present) are pinned by sha256 in
  `lensmodel.HFF_CATS`. The release readmes quote rms 0.72″, 0.63″, n/a, 2.41″, 0.48″ and n/a.

- BUFFALO Abell 370 lens models, Niemiec et al. 2023, MNRAS 524, 2883 (arXiv:2307.03778), MAST HLSP
  https://archive.stsci.edu/hlsp/buffalo (DOI 10.17909/t9-w6tj-wp63); model A (Lenstool image-plane, quoted rms
  0.90″): maps and `input.par` listed in `download_scripts/hlsp_buffalo_hst_model_abell370_multi_niemiec-a_v1.0_lens-models-download.sh`
  (accessed 2026-10-22). Pinned 2026-10-22: the `main` x / y deflection and κ maps of `a/sl/`
  (`lensmodel.ABELL370_BUFFALO_A_MAPS`) and the gold strong-lensing list
  `catalogs/niemiec-lensing-dr1/hlsp_buffalo_hst_multi_abell370_multi_v1.0_sl-gold.dat` (104 images, sha256 in
  `scripts/a370_buffalo.py`). Lagattuta et al. 2019, MNRAS 485, 3738, CDS `J/MNRAS/485/3738` (accessed
  2026-10-22; `table1.dat` pinned by sha256 in `scripts/a370_buffalo.py`).

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
- VizieR `J/MNRAS/506/928` ("Stellar mass-size relation for low M* galaxies"; table `a1063clu`, HFF-DeepSpace
  SExtractor shapes), accessed 2026-10-16: HST axes of S1063 `arcs` rows 1316 and 2184
  (`results/arcs/abells1063-cats/vetting/hst_shapes.json`).

### MACS0717 field run (accessed 2026-10-08; D-041)

- ASTRODEEP Frontier Fields catalogues (Di Criscienzo et al. 2017, A&A 607, A30), VizieR `J/A+A/607/A30`, tables
  `m0717cla` (positions) and `m0717clz` (zbest, SIQR, relFlag, f_zbest), accessed 2026-10-17: 4,082 rows each; kept rows: relFlag = 1 or
  f_zbest = 1 (2,928)
  (`scripts/astrodeep_zout.py`, W1 shear photo-z; `results/w1_shear_macs0717/`).
- HFF-DeepSpace v3.9 (Shipley et al. 2018) download page http://cosmos.phy.tufts.edu/~danilo/HFF/Download.html:
  the catalogue tarballs returned HTTP 401 on 2026-10-17 (not used).

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
  - Pinned 2026-10-21 as `macs0717-sharon` / `macs0717-keeton` (`lensmodel.MACS0717_SHARON_V4COR`,
    `MACS0717_KEETON_V4`; the deflection and κ files above, re-downloaded and matching), plus
    `hlsp_frontier_model_macs0717_sharon_v4cor_z02-magnif.fits` 64,008,000 B `133be386832dda060afd180c716bb644e6ee000523add666bf5ad97edfffd436`
    and `hlsp_frontier_model_macs0717_keeton_v4_z02-magnif.fits` 64,039,680 B `26750c9f50083de610f7ab42ba8c02ee2d45310195a54cd9f484b3fcf97f8415`.
    Readmes (accessed 2026-10-21): `sharon/v4cor/hlsp_frontier_model_macs0717_sharon_v4cor_readme.pdf` (the shared
    four-cluster v4 document, sha256 `12a284d22368f29491670180d07f5776d35a28ed754743b866fdf7184af39636`, Table 7) and
    `keeton/v4/hlsp_frontier_model_macs0717_keeton_v4_readme.txt` (sha256
    `2a85d85ebb041c3f6fa1736bd88fe75bfe1afe2f3c104bdf5b87310c79c7cda6`; 29 images of 9 gold systems, rms 0.79″).

### MACS1149 field run (accessed 2026-10-08; D-037)

- JWST CANUCS program 1208, NIRCam level 3 `jw01208-o008_t004_nircam_clear-*`, 8 bands (F090W–F444W), jwst 3.0.0,
  S3 `s3://stpubdata/jwst/public/jw01208/L3/t/o008/`. Manifests: `data/manifests/macs1149*.ecsv`.
- HFF Sharon v4cor MACS1149 deflection maps (in `MODELS` as `macs1149-sharon` since 2026-10-20, with the kappa and
  z02-magnif maps; sha256 in `lensmodel.MACS1149_SHARON_V4COR`; readme
  `hlsp_frontier_model_macs1149_sharon_v4cor_readme.pdf`, sha256
  `12a284d22368f29491670180d07f5776d35a28ed754743b866fdf7184af39636`, Table 5 gold images; accessed 2026-10-20):
  `https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/sharon/v4cor/hlsp_frontier_model_macs1149_sharon_v4cor_{x,y}-arcsec-deflect.fits`,
  64,008,000 B each. sha256: x `44f1d21a78cead051c08295d13d7c9189638971cf1e9445bb2caa641a129324f`,
  y `54508deb502d9bcafca9733c3b28f3b4d794b64f11167589b56c1fc37d66a15c`.
- DJA v7 imaging index (dawn-cph GitHub Pages `/dja/imaging/v7/`): no MACS1149 mosaic.
- CANUCS DR1 (DOI 10.17909/18nv-np70) MACS1149 cluster-field photometry and EAzY photo-z:
  https://archive.stsci.edu/hlsps/canucs/dr1/macs1149/clu/hlsp_canucs_jwst-hst_multi_macs1149-clu_multi_v1_photometry-cat.fits.gz,
  29,779,449 B, sha256 `08ab67347f2c3dfe4f743cc1b2a9d4b77a1e40eb66ddeb611648bb296adc0739`. The same directory tree has
  CANUCS lens models (`model/`: deflection, κ, γ, best and samples). Used from 2026-10-22 (D-034 addendum f): the
  best-fit 100 mas `dpl1` / `dpl2` / `kappa` maps (Rihtarsic, preliminary v1; readme
  `hlsp_canucs_jwst-hst_multi_macs1149-lenstool-readme_multi_v1_model.txt`, accessed 2026-10-22; sha256 in
  `lensmodel.MACS1149_CANUCS_MAPS`) and the F200W 40 mas BCG-subtracted mosaic
  `.../macs1149/clu/hlsp_canucs_jwst_nircam_macs1149-clu-40mas_f200w_v1_bgsub-sci.fits.gz`, 104,672,989 B, sha256
  `34d8f2609fc49f0c88fc93d4ddcccc613fee27348204bdbeb86297d4e8339c26`.
- System 16 spec-z check (accessed 2026-10-20; D-037 addendum):
  - Treu et al. 2016, ApJ 817, 60, doi:10.3847/0004-637X/817/1/60; VizieR J/ApJ/817/60 (table3: multiple images with
    `z` and new `zspec` from HST-GO-13459/14041 grism and VLT-MUSE 294.A-5032), doi:10.26093/cds/vizier.18170060.
  - DJA NIRSpec emission-line table v4.4 (msaexp extractions):
    https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.4.csv.gz, 133,004,545 B,
    Last-Modified 2025-09-04, sha256 `728afcdd5ec4d164b9d816ea903e648b1d4937ae52c8a271056c9119cabf1ccf`. Per-spectrum
    files under `extractions/macs1149-v4/` returned 403 from the cloud session.
  - DJA NIRSpec emission-line table v4.6 (accessed 2026-10-21, streamed): same path with `_v4.6`, 171,878,708 B,
    Last-Modified 2026-09-15, sha256 `ee64bb460f417d3fcf96fc53818c06ca0e27d784833063ae8eb48feb18c6ddf3`, 99,478 rows.
    v4.5 (2026-02-23, 151,445,355 B) exists, not read.
  - CATS v4.1 `arcs.txt` (already pinned in `lensmodel.MODELS["macs1149-cats"]`).
### MACS0416 field run (accessed 2026-10-08; D-038)

- JWST CANUCS program 1208, NIRCam level 3 `jw01208-o004_t002`, 8 bands, jwst 3.0.0. Manifests:
  `data/manifests/macs0416*.ecsv`. Also available: PEARLS 1176 o211/o212/o213 `t009`.
- DJA v7 index https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html: no MACS0416 mosaic.
- CANUCS DR1 (DOI 10.17909/18nv-np70) MACS0416 cluster-field photometry and EAzY photo-z:
  https://archive.stsci.edu/hlsps/canucs/dr1/macs0416/clu/hlsp_canucs_jwst-hst_multi_macs0416-clu_multi_v1_photometry-cat.fits.gz,
  31,693,111 B, sha256 `339107a5c5041621d7bed4ecc8b4a51b5148a6913d4e513c3a5e783363f11bff`.
- MACS0417 (`arcs`, accessed 2026-10-16; `results/arcs/macs0417-relics/`):
  - JWST 1208 F200W pipeline catalogue `jw01208-o006_t003_nircam_clear-f200w_cat.ecsv` (MAST
    `mast:JWST/product/...`), 5,030,481 B, sha256 `e55dffb990a8e8124e37658933a47952a77b357b36dc5c8f523a5b487d423c24`;
  - CANUCS DR1 cluster photometry and EAzY photo-z
    https://archive.stsci.edu/hlsps/canucs/dr1/macs0417/clu/hlsp_canucs_jwst-hst_multi_macs0417-clu_multi_v1_photometry-cat.fits.gz,
    Last-Modified 2025-04-29, 27,820,129 B, sha256 `ca6b48aa6404444d2e0b8965a35376ff032917071ea64cbf6306d06ee57d90dd`;
  - CANUCS DR1 F200W BCG-subtracted 40 mas mosaic (accessed 2026-10-22, streamed, deleted after use;
    `results/lens_consistency/macs0417_canucsref/`):
    https://archive.stsci.edu/hlsps/canucs/dr1/macs0417/clu/hlsp_canucs_jwst_nircam_macs0417-clu-40mas_f200w_v1_bgsub-sci.fits.gz,
    104,841,074 B, sha256 `591228193a679d52ac523e8f5e0b7806327d7b2c382a57cbe8b0787d1eaea1cd`;
  - CANUCS DR1 Lenstool maps (Desprez, preliminary v1; readme
    `hlsp_canucs_jwst-hst_multi_macs0417-lenstool-readme_multi_v1_model.txt`): `dpl1`, `dpl2`, `kappa`
    `-dlsds1-best-100mas` under https://archive.stsci.edu/hlsps/canucs/dr1/macs0417/model/ (sha256 in
    `lensmodel.MACS0417_CANUCS_MAPS`; Last-Modified 2025-04-26);
  - RELICS Lenstool v1 maps (Mahler et al. 2019, MNRAS 485, 3954; DOI 10.17909/T9SP45)
    https://archive.stsci.edu/hlsps/relics/macs0417m11/models/lenstool/v1/ (sha256 in
    `lensmodel.MACS0417_RELICS_LENSTOOL`). GLAFIC v2/v3 and Lenstool v2 also exist there (not used).
- MACS0647 (`arcs`, accessed 2026-10-16; `results/arcs/macs0647-clash/`):
  - JWST 1433 (PI Coe) F200W pipeline catalogues (MAST `mast:JWST/product/...`):
    `jw01433-o010_t001_nircam_clear-f200w_cat.ecsv` 2,144,120 B, sha256
    `473e32d23c3fbc857e89d6b48196130d77670be2b9df2b4dfae5df82e71e7a7e`;
    `jw01433-o020_t001_nircam_clear-f200w_cat.ecsv` 2,761,197 B, sha256
    `29f74c91cf8984b3f44118073e183ef0bce9262addf7637ffed6dfbf0b642663`;
  - DJA v7.0 https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/macs0647-grizli-v7.0-fix.photoz.tar.gz
    (Last-Modified 2023-07-17, 41,494,640 B, sha256 `d56246b571ccfc78a54338a7ab16009fa8b73b2cdbcfefba4a2617c0058cdf0c`),
    member `macs0647-grizli-v7.0-fix.eazypy.zout.fits` (7,956 rows; sha256
    `ceb199c61f1ca67df6afc34df6ed7ad61267f59ccaf27e9a2fd1e64c2662c04e`, recorded 2026-10-17) extracted while streaming;
  - CLASH Zitrin v2 maps (Zitrin et al. 2009, MNRAS 396, 1985; 2013, ApJ 762, L30; CLASH, Postman et al. 2012,
    ApJS 199, 25): https://archive.stsci.edu/missions/hlsp/clash/macs0647/models/zitrin/{nfw,ltm-gauss}/v2/
    (sha256 in `lensmodel.MACS0647_CLASH_NFW`, `MACS0647_CLASH_LTM`; readme 2015-02, z_lens 0.591).
- MACS J2129.4-0741 (`results/lens_consistency/macs2129_clash_heldout/`, accessed 2026-10-25):
  - Cerny et al. 2026, "Strong LensIng and Cluster Evolution (SLICE) with JWST", ApJ 1001, 60, doi:10.3847/1538-4357/ae41b0,
    arXiv:2503.17498v2: e-print `main.tex` Table 8 (MACS J2129 block; e-print sha256 f7120d94…a9a7 in `slice_images.ecsv`
    meta). Its models are said to be in the LAM atlas, which answered "Restricted access" (security suspension) on 2026-10-25.
  - CLASH Zitrin v2 MACS2129 maps (PIEMD + eNFW, LTM-Gauss; readme 2015-02-18):
    https://archive.stsci.edu/missions/hlsp/clash/macs2129/models/zitrin/ (sha256 in `lensmodel.MACS2129_CLASH_*`;
    readme sha256 nfw 7ffa3f92…947d, ltm-gauss 3e1756ea…dae0da).
- SLICE × RELICS MACS J0358 (`results/lens_consistency/slice_relics_macs0358/`, accessed 2026-10-25):
  - SLICE Table 8 block of MACS J0358.8−2955 (same e-print; sha256 in `slice_images.ecsv` meta).
  - RELICS Lenstool v1 and GLAFIC v2 maps of Abell 3192:
    https://archive.stsci.edu/hlsps/relics/abell3192/models/ (HLSP DOI 10.17909/T9SP45; sha256 in
    `lensmodel.A3192_RELICS_{LENSTOOL,GLAFIC}`; no model readme or constraint list in those directories).
  - Hsu, Ebeling & Richard 2013, "The three-dimensional geometry and merger history of the massive galaxy cluster
    MACS J0358.8-2955", MNRAS 429, 833, doi:10.1093/mnras/sts379, arXiv:1209.2492 (from the SLICE e-print
    bibliography).
- SLICE × RELICS RX J0152 (`results/lens_consistency/slice_relics_rxj0152/`, accessed 2026-10-25):
  - SLICE Table 8 block of RX J0152.7−1357 (same e-print; sha256 in `slice_images.ecsv` meta).
  - RELICS Lenstool v1, GLAFIC v2 and Zitrin LTM-Gauss v1 maps:
    https://archive.stsci.edu/hlsps/relics/cl0152m13/models/ (HLSP DOI 10.17909/T9SP45; sha256 in
    `lensmodel.CL0152_RELICS_{LENSTOOL,GLAFIC,ZITRIN}`).
  - Acebron et al. 2019, "RELICS: High-resolution Constraints on the Inner Mass Distribution of the z = 0.83 Merging
    Cluster RXJ0152.7-1357 from Strong Lensing", ApJ 874, 132, doi:10.3847/1538-4357/ab0adf, arXiv:1810.08122
    (from the SLICE e-print bibliography).
- SLICE × RELICS (`results/lens_consistency/slice_relics_heldout/`, accessed 2026-10-25):
  - SLICE Table 8 blocks of RXC J0232.2-4420, RXC J2211.7-0349, MACS J0553.4-3342 (same e-print as above; sha256 in
    `slice_images.ecsv` meta).
  - RELICS Lenstool v1 maps (Coe et al. 2019, ApJ 884, 85; HLSP DOI 10.17909/T9SP45):
    https://archive.stsci.edu/hlsps/relics/{rxc0232m44,rxc2211m03,macs0553m33}/models/lenstool/ (sha256 in
    `lensmodel.{RXC0232,RXC2211,MACS0553}_RELICS_LENSTOOL`). The glafic directories of macs0553m33 and rxc2211m03 list no
    deflection maps; smacs2031m40 has no models.
- MACS0647 JWST-era image lists (`results/lens_consistency/macs0647_inputs/`, accessed 2026-10-23):
  - Meena et al. 2023, ApJL 944, L6, arXiv:2211.13334, DOI 10.3847/2041-8213/acb645: `lensed_systems.tex` of the
    arXiv source (sha256 in `summary.json` `images.tex_sha256`), 72 images in 23 systems.
  - Nishida et al. 2025, "The Jackknife method as a new approach to validate strong lens mass models",
    arXiv:2505.00553: `ms.tex` appendix tables `Table:MACS0647_1` / `_2` (sha256 in `summary.json`), 86 images from
    28 sources (JWST GLAFIC constraints).
  - DJA v7.0 `https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/macs0647-grizli-v7.0-fix_phot.fits`
    (25,427,520 B; headers only, HTTP range, 2026-10-23).
- RXCJ2211-0350 inputs (`results/lens_consistency/rxcj2211_inputs/`, accessed 2026-10-24):
  - Zenodo record 10.5281/zenodo.20312327 (https://zenodo.org/records/20312327, published 2026-05-20, CC0): VENUS
    GLAFIC v0.91 `RXCJ2211_lens_model_GLAFIC.zip` (md5 87a9bfab62c1e5357e676f93d44db4ac; member sha256 in
    `summary.json`); the record cites arXiv:2512.05180 (not read).
  - Other Zenodo lens-model records seen in the same search, not used: SPT-CL J0546-5345 LTM (10.5281/zenodo.16534337),
    Abell S1063 WSLAP+ JWST 2025 (10.5281/zenodo.17726321, 1.0 GB zip), CANUCS thesis data with MACS0416 and Bullet
    image lists (10.5281/zenodo.22883437).
- Limousin et al. 2016, "Strong-lensing analysis of MACS J0717.5+3745 from Hubble Frontier Fields observations: How
  well can the mass distribution be constrained?", A&A 588, A99, arXiv:1510.08077v2 (ar5iv HTML Tables 5–6 read via a
  summarising fetch, accessed 2026-10-25): ten spectroscopic systems 1, 3, 4, 6, 12, 13, 14, 15, 19, 68; used only to
  show that no spectroscopic image is held out from the MARS MACS0717 catalogue (CHANGELOG 2026-10-25).
- MARS HFF models (`results/lens_consistency/mars_hff_inputs/`, accessed 2026-10-24): Cha & Jee 2023, "Model-Independent
  Mass Reconstruction of the Hubble Frontier Field Clusters with MARS Based on Self-Consistent Strong Lensing Data",
  arXiv:2301.08765 (not read; README and code in the record read); Zenodo record 10.5281/zenodo.7804575
  (https://zenodo.org/records/7804575, published 2023-04-06, CC-BY-4.0): `MARS_HFF_clusters_result.zip` (2,713,478 B,
  md5 2d87f9a36aadb430f7a7b937fe70b8c2; sha256 in `summary.json`).
- Schuldt et al. 2024, MACS J1149.5+2223 VLT/MUSE spectroscopy and strong-lensing model, A&A 689, A42,
  arXiv:2405.12287v2 (HTML Table 2 read and transcribed to `results/lens_consistency/mars_m1149_heldout/schuldt24_new_images.ecsv`,
  accessed 2026-10-24); VizieR J/A+A/689/A42 (only the MUSE catalogue of the new pointing, tablea1). Exact title not
  recorded.
- SMACS 0723 Chow+2023 (`results/lens_consistency/smacs_chow23/`, accessed 2026-10-24): Chow et al. 2023, "A Lens
  Finder Map to check claimed High-z Galaxies behind SMACS J0723.3-7327", arXiv:2310.09790 (HTML read for the GLAFIC
  method, the five spectroscopic systems and the 0.492″ rms); Zenodo record 10.5281/zenodo.10158248
  (https://zenodo.org/records/10158248, published 2023-11-20): `Chow23+_LensModel_alpha_x.fits`, `..._alpha_y.fits`
  (42,514,560 B each; sha256 in `summary.json`), `readme.md`.
- Abell S1063 WSLAP+ JWST 2025 (`results/lens_consistency/s1063_wslap/`, accessed 2026-10-24): Zenodo record
  10.5281/zenodo.17726321 (https://zenodo.org/records/17726321, published 2025-11-26, CC-BY-SA-4.0), zip
  `LensModel_AS1063_WSLAPplus_JWST_Case3_2025.zip` (1,016,075,135 B, md5 6a6558557078c1dc2b0c8a2a13d17944); members
  `AlphaX_DdsDs1_arcsec.fits`, `AlphaY_DdsDs1_arcsec.fits` read by range (CRC-32 and sha256 in `summary.json`). No paper
  read; the record names BUFFALO identifications.
- SPT-CL J0546-5345 inputs (`results/lens_consistency/spt0546_inputs/`, accessed 2026-10-24):
  - Allingham et al. 2025, "SLICE: SPT-CL J0546-5345 -- A prominent strong-lensing cluster at z=1.07", ApJL 990, 1,
    arXiv:2507.08949v2 (HTML read for the rms ≈ 1.2″ and the z_s 3.5 anchor).
  - Zenodo record 10.5281/zenodo.16534337 (https://zenodo.org/records/16534337, published 2025-07-11, CC-BY-4.0):
    `alpha_x/y_sptcl0546_z3p5_20250611_LTM_Gc.fits`, `actcl0546_mulfile_v1.mul`, `actcl0546_spec-z.cat` (md5 checked,
    sha256 in `summary.json`).
- SN Eos follow-up epochs (`results/transients/epoch_inventory/`, MAST metadata accessed 2026-10-25): JWST programs
  9493 o001 (2025-10-08) and 12299 o001 (2026-09-28), NIRCam, target `SN_Eos_Image_101.1`; SN Eos spectroscopy:
  arXiv:2601.04156 (z 5.13; read from the search excerpt only, paper not opened).
- SN Eos flux-ratio literature check (`results/transients/sneos_ratio/LITERATURE.md`, accessed 2026-10-25):
  - Coulter et al. 2026, "A spectroscopically confirmed, strongly lensed, metal-poor Type II supernova at z = 5.13",
    arXiv:2601.04156v2 (2026-01-20), HTML https://arxiv.org/html/2601.04156; Tables S1 (photometry) and S2 (|μ|, τ)
    transcribed into `scripts/sneos_literature.py`. Names HST PID 17504 (2024-03) as the archival detection.
  - Allingham et al. 2026 (arXiv:2602.14074, HTML as served 2026-10-25), Table 3 (|μ|, τ).
- MACS J1931.8-2635 inputs (`results/lens_consistency/macs1931_inputs/`, accessed 2026-10-24):
  - Allingham et al. 2026, "VENUS: Strong-lensing model of MACS J1931.8-2635 -- revealing the farthest multiply imaged
    supernova", arXiv:2602.14074.
  - Zenodo record 10.5281/zenodo.18646525 (https://zenodo.org/records/18646525, published 2026-02-15): arcs list,
    MUSE z list, README, `chires.txt` (sha256 in `summary.json`); `Alpha_x/y.fits.zip` headers only (HTTP range).
  - Caminha et al. 2019, arXiv:1903.05103 (cited by the README for the spec-z of systems 1–7; not read).
  - Same record, used by `results/lens_consistency/macs1931_gatek/` (accessed 2026-10-24): `Alpha_x.fits.zip`,
    `Alpha_y.fits.zip`, `sourcepos.txt`, `chires.txt` (md5 checked against the record, sha256 in `summary.json`), and
    `Readme_Magnifications_summary.txt` (its three `k_dlds` values are constants in `scripts/macs1931_gatek.py`).
  - JWST GO 6882 (VENUS) level-3 tile `jw06882-o038_t041` (MAST CAOM metadata, accessed 2026-10-24; used by
    `results/lens_consistency/macs1931_counter/`): 10 NIRCam bands; F200W `_i2d.fits` / `_cat.ecsv` / `_segm.fits` under
    `s3://stpubdata/jwst/public/jw06882/L3/t/o038/` (sizes in `summary.json`; header only read).
- Bullet Cluster (1E 0657-56) inputs (`results/lens_consistency/bullet_inputs/`, accessed 2026-10-24):
  - Rihtaršič et al. 2026, "Mapping dark matter in the Bullet Cluster using JWST imaging and spectroscopy",
    A&A (2026), arXiv:2601.22245; image catalogue also as VizieR J/A+A/710/A207 (identifier from the GAVO registry record
    https://dc.g-vo.org/rr/q/lp/custom/CDS.VizieR/J/A+A/710/A207; the VizieR table itself was not read, it timed out).
  - CANUCS Bullet Cluster DR1 (GO 4598), https://niriss.github.io/data_release_bullet_cluster.html; files under
    https://www.canfar.net/storage/vault/file/canucs/BULLET_DR1/ (`lensing/bullet_gold_v4.dat`,
    `lensing/bullet_multim_full.dat`, `lensing/README.dat`: sha256 in `summary.json`; photometry catalogue
    `catalogs/hlsp_canucs_jwst-hst_multi_bullet-clu_multi_v1_photometry-cat.fits`, 51,721,920 B, header only).
- MACS0417 / MACS1423 multiple-image lists (`results/lens_consistency/macs0417_images/`, accessed 2026-10-21):
  - Jauzac et al. 2019, MNRAS 483, 3082, arXiv:1811.02505, DOI 10.1093/mnras/sty3312: Table 4 from
    https://arxiv.org/e-print/1811.02505 (`table4.tex`, sha256 `7e05be79da64fcfb62e71f081c1b18bebc249b6641e6f598930f01c6be53b4cc`).
  - Mahler et al. 2019, ApJ 873, 96, arXiv:1810.13439, DOI 10.3847/1538-4357/ab042b (RELICS; the gold/silver set
    Jauzac+19 reuse; `arcs_bronze.tex` in its e-print; not transcribed).
  - KALEIDOSCOPE, MNRAS 533, 4500 (2024), arXiv:2405.04577, DOI 10.1093/mnras/stae2069: table `tab:mulimg` from
    https://arxiv.org/e-print/2405.04577 (`strong_lensing.tex`, sha256 `86fd419c82d5d97c06c0542b1b621bf327728cd0ac5889918f5e28ab59f0c69d`).
  - Not in VizieR (checked J/ApJ/873/96, J/MNRAS/483/3082, J/MNRAS/533/4500). CANUCS v1 has no image catalogue
    (readme: with v2, Desprez et al. in prep.).
- MACS0417 / MACS1423 CANUCS extra predicted images (accessed 2026-10-21; `results/lens_consistency/macs0417_extra/`): the
  F200W catalogues above (sha256 re-checked) and their `_i2d.fits` by S3 byte ranges (o006 / o010); the CANUCS DR1
  model maps (`lensmodel.MACS0417_CANUCS_MAPS`, `MACS1423_CANUCS_MAPS`, sha256 as pinned) fetched through
  `https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/canucs/dr1/<field>/model/<file>` (archive.stsci.edu
  failed TLS from the cloud session).
- MACS0417 frame offset (accessed 2026-10-21; `results/lens_consistency/macs0417_frame/`): RELICS macs0417m11 HST
  catalogue v1, `https://archive.stsci.edu/hlsps/relics/macs0417m11/catalogs/hlsp_relics_hst_acs-wfc3ir_macs0417m11_multi_v1_cat.txt`,
  2,182,873 B, sha256 `26a210e599f37ae066cc3cfb5e8d53777f39ad322710cdc002f6372fc79f24b8` (header: produced 2017-06-19
  by Dan Coe, GO 14096); the o006_t003 F200W `_cat.ecsv` above via `https://stpubdata.s3.amazonaws.com/`, sha256
  `e55dffb990a8e8124e37658933a47952a77b357b36dc5c8f523a5b487d423c24`.
- MACS1423 (`arcs`, accessed 2026-10-16; `results/arcs/macs1423-canucs/`):
  - JWST 1208 F200W pipeline catalogue `jw01208-o010_t005_nircam_clear-f200w_cat.ecsv` (MAST
    `mast:JWST/product/...`), 4,987,129 B, sha256 `1a968574d0f206c72a74554c1b3cbbf71e9a5219a70a553a1458e80dd4d5bfa3`;
    its `_i2d.fits` read by S3 byte ranges (`s3://stpubdata/jwst/public/jw01208/L3/t/o010/`);
  - CANUCS DR1 cluster photometry and EAzY photo-z
    https://archive.stsci.edu/hlsps/canucs/dr1/macs1423/clu/hlsp_canucs_jwst-hst_multi_macs1423-clu_multi_v1_photometry-cat.fits.gz,
    25,278,098 B, sha256 `c257d8dff95319f4eea981cc838269c04aff0a79bac603c7ed64100deac622b0`;
  - CANUCS DR1 Lenstool maps (Desprez, preliminary v1, 5 spectroscopic systems; readme
    `hlsp_canucs_jwst-hst_multi_macs1423-lenstool-readme_multi_v1_model.txt`): `dpl1`, `dpl2`, `kappa`
    `-dlsds1-best-100mas` under https://archive.stsci.edu/hlsps/canucs/dr1/macs1423/model/ (sha256 in
    `lensmodel.MACS1423_CANUCS_MAPS`).
- JWST-era lens maps for other DJA photo-z fields (checked 2026-10-16, not used yet):
  - RXJ2129: six models updated with JWST NIRCam (Dhanasingham et al. 2026, arXiv:2604.11882; GLAFIC, Chen2020,
    HoliGRALE, WSLAP+, Zitrin-Analytic, Lenstool), Zenodo record 21460864, one 4.30 GB zip
    https://zenodo.org/api/records/21460864/files/RXJ2129+0005_Models.zip/content (range requests work; README
    says D_LS/D_S = 1, but the GLAFIC κ header has ZS = 3: check before use). Used 2026-10-16 for `arcs`: Lenstool
    `alphaX/alphaY/kappa_DLSDS1` and GLAFIC `x/y-arcsec-deflect`, `kappa` (sha256 in `lensmodel.RXJ2129_*_MAPS`),
    with JWST 2767 F200W catalogue `jw02767-o002_t001_nircam_clear-f200w_cat.ecsv`
    (https://stpubdata.s3.amazonaws.com/jwst/public/jw02767/L3/t/o002/, sha256 `2ed1a2c7666734483974e568e1c3b332eea223ea487b4cd279dd89b4e9261f53`)
    and DJA v7.0 https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rxj2129-grizli-v7.0-fix.photoz.tar.gz (member
    `…eazypy.zout.fits`, sha256 `651c0b24f09bafa27c36a2285f081d42884be3c6dcbbb1708014339afa726217`);
    image list: the paper's Appendix B table `Systems_RA_and_dec` read from the arXiv source
    (https://arxiv.org/e-print/2604.11882, accessed 2026-10-23), transcribed in `scripts/rxj2129_djaref.py`;
    DJA v7.0 photometry https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rxj2129-grizli-v7.0-fix_phot.fits
    (9,152,640 B, headers read 2026-10-23; v7.1–v7.5 return 403);
  - MRG0138 (Ertl et al. 2025, arXiv:2503.09694; Suyu et al. 2025, arXiv:2509.12319), RXCJ0600 (Furtak et al.
    2024, arXiv:2404.03286): JWST-era models without public maps; rj0018p1626: no JWST-era model found (RELICS
    HST-era maps only).
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
- D-060 addendum (line C, accessed 2026-10-21):
  - MAST CAOM TAP, `https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync`, `dbo.ObsPointing` (HST calib_level 3, F814W and
    F160W): coverage of the D-056 class, `results/w2_hst_blend/hst_coverage.ecsv`.
  - MAST HAPcut (`astroquery.mast.Hapcut`, `https://mast.stsci.edu/hapcut/api/v0.1/astrocut`): HAP-MVM F814W cutouts
    (programmes 9744, 14254, 15320, 15652, 17199 among the selected systems).
  - Gaia DR3 `gaiadr3.gaia_source` via astroquery.gaia: PSF stars.
  - Gilman et al. 2020a, MNRAS 491, 6077 (DESJ0405-3308 photometric z_d ∼ 0.3), as quoted in Nierenberg et al. 2023,
    arXiv:2309.10101 §3 (coordinator literature check, post hoc).

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
  undocumented; uncompressed) also answers; D-068 uses it only for vetting neighbours without a recorded tar offset.
- The other 21 field tars `bulk/gb{F}.tar` (D-068) are streamed with HTTP byte-range reads, never stored. Sizes and
  Last-Modified (HEAD, 2026-10-08) are in `moa.TAR_BYTES` / `moa.TAR_LAST_MODIFIED` (7.7–508.5 GB, ≈ 2.4 TB in
  all; Last-Modified 2023-10-10 … 13). No whole-file sha256 exists for a streamed tar: each is pinned by the sha256
  of every 64 MiB range (tracked per chunk in `results/w3_moa/prescreen/`) and a field digest over them
  (`moa_stream.range_digest`, rows in `data/manifests/moa_ii.ecsv`). The archive occasionally answers a range
  request with HTTP 200 (whole file); such replies are closed unread and retried.
- Nunota et al. 2024 Table 1 N_s (20 fields; `moa.NUNOTA_NS`) is the adopted N_s per field where published
  (D-068); gb6 and gb22 use the N_s-per-Cut-0-object model.
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
- Legacy Surveys DR10 Tractor + DR10 photometric redshifts (Data Lab TAP tables `ls_dr10.tractor`,
  `ls_dr10.photo_z`, joined on `ls_id`; https://www.legacysurvey.org/dr10/), same box queries, dereddened z < 22,
  accessed 2026-10-18 (D-070 addendum 2).

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

## S1 burst twins: CHIME/FRB Catalog 2 dynamic spectra (accessed 2026-10-21; D-071 addendum)

- CHIME/FRB Catalog 2 Stokes-I dynamic spectra, CANFAR doi:10.11570/25.0066,
  `data/dynamic_spectra/data/<TNS>_stokesi_dynamic_spectrum.h5` (4,536 files, 56.4 GB in the listing). Listing
  https://cadc-west-01.canfar.net/vault/nodes/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/dynamic_spectra/data?limit=10000 ;
  files https://cadc-west-01.canfar.net/vault/files/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/dynamic_spectra/data/<file>
  (`curl -L` to a minoc pre-signed URL). Route from issue #153 (2026-10-11); 30 files pinned in
  `data/manifests/s1_chime_dynspec.ecsv`. Format: arXiv:2601.09399 §3.3. Read with h5py 3.16.0 (PyPI).

## E1 + Einstein Probe (accessed 2026-10-18; D-074 addendum 7)
- GCN Circulars JSON archive, https://gcn.nasa.gov/circulars/archive.json.tar.gz (redirects to the GCN S3 bucket;
  rewritten daily; 45,863 circulars on 2026-10-18, `build_summary.json`), streamed; sha256 recorded in `data/manifests/e1_ep.ecsv`; the parsed table
  `results/e1_ep/ep_events.ecsv` keeps each row's circular id.
- GCN Circular 40906 (EP250702a likely associated with GRB 250702B,C,D,E), https://gcn.nasa.gov/circulars/40906.

## E1 + Swift BAT (accessed 2026-10-18; D-074 addendum 6)
- Swift GRB compilation by the burst advocates, HEASARC TAP table `swiftgrbba`
  (https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync; columns name, time, trigger_obs, ra, dec, error_radius),
  2,043 rows, 2004-12-17 to 2026-10-07 (session check, not stored), sha256 pinned in `data/manifests/e1_swift.ecsv`. HEASARC `swiftgrb`
  was checked and rejected: 872 rows, ending 2012-12-29 (session check, not stored).

## E1 causal event network: GBM × ICECAT-1 × GWTC × CHIME/FRB Cat 2 (accessed 2026-10-09; D-074)

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
- GWTC parameter-estimation sky maps (accessed 2026-10-10; D-074 addendum 4), Zenodo records, sha256 and sizes in
  `data/manifests/e1_gw_skymaps.ecsv`: GWTC-2.1 rec. 6513631 `IGWN-GWTC2p1-v2-PESkyMaps.tar.gz` (54 maps);
  GWTC-3 rec. 8177023 `IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz` (36); GWTC-4.1 rec. 20275769
  `IGWN-GWTC4p1-18965dda8_5-Archived_Skymaps.tar.gz` (88); GWTC-5.0 rec. 20348005
  `IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz` (104). Streamed into memory, reduced to nside 32, never
  stored. The 109 GWTC-4.1 / 5.0 list entries without a map also have no PE parameters in the GWOSC CSV.
  GW170817 has no map in these releases (a Gaussian at SSS17a stands in, positive control only).
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

## SN Refsdal light curves (accessed 2026-10-21; D-072 addendum)
- Kelly et al. 2023, ApJ 948, 93, arXiv:2305.06377: HST F105W / F125W / F160W difference-image photometry of
  S1–S4 and SX, CDS J/ApJ/948/93 (`table3`, 805 rows; `images`), streamed via TAPVizieR
  (https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync); SX − S1 delay 376.0 +5.6 −5.5 d from the abstract.
  `scripts/s3_refsdal.py`.

## D1 distance self-consistency (accessed 2026-10-09; D-073)

- H0LiCOW public distance posteriors, https://github.com/shsuyu/H0LiCOW-public (commit 57cf973, 2025-05-14):
  `h0licow_distance_chains/*`, `MontePython_cosmo_sampling/data/timedelay_6lenses/B1608_Dd_Ddt_params.dat` and the
  lens redshifts in `MontePython_cosmo_sampling/likelihoods/timedelay_6lenses/__init__.py`. Papers: Wong et al.
  2020, MNRAS 498, 1420 (arXiv:1907.04869); Suyu et al. 2010 (B1608 D_dt); Jee et al. 2019, Science 365, 1134
  (B1608 D_d); Chen et al. 2019, MNRAS 490, 1743 (HE0435, RXJ1131, PG1115); Birrer et al. 2019, MNRAS 484, 4726
  (J1206); Rusu et al. 2020 (WFI2033, arXiv:1905.09338).
- TDCOSMO 2025 public release, https://github.com/TDCOSMO/TDCOSMO2025_public (commit d7f38db, 2026-01-21):
  `TDCOSMO_sample/TDCOSMO_data/SDSS1206+4332/final_D_d.npy`, `final_D_dt.npy`, `TDCOSMO_sample/tdcosmo_sample.yaml`;
  `final_composite_D_dt.npy` and `final_power_law_D_dt.npy` (read 2026-10-10, D-073 addendum 3).
  hierArc kinematic likelihood pickles `TDCOSMO_sample/<lens>_const_processed.pkl` (8 lenses),
  `ExternalLenses/SLACS/slacs_{kcwi,sdss}_const_processed.pkl`, `ExternalLenses/SL2S/sl2s_const_processed.pkl`
  (read 2026-10-10 with a numpy-only unpickler, D-073 addendum 4).
- hierArc 1.2.0 (PyPI `hierarc`, https://github.com/sibirrer/hierArc; read 2026-10-10, not a dependency):
  `Likelihood/LensLikelihood/kin_likelihood.py`, `Likelihood/kin_scaling.py`, `Likelihood/transformed_cosmography.py`
  give the σ_v model and λ convention reproduced in `distance_consistency.kin_ln_ratio_loglike`.
  Paper: TDCOSMO Collaboration 2025, A&A 704, A63 (arXiv:2506.03023). Other per-lens files there (names only,
  listed 2026-10-09; read 2026-10-09 by `d1_distance.py tdcosmo`, sha256 in `data/manifests/d1_distance.ecsv`): power-law D_dt chains `HE0435-1223/he_powerlaw_Ddt.dat`,
  `PG1115+080/pg_powerlaw_Ddt.dat`, `RXJ1131-1231/rxj_powerlaw_Ddt.dat`, `WFI2033-4723/wfi2033_pl_dt_nokext.dat`,
  `B1608+656/B1608_Dtmod_n5e5.dat`, `DES0408-5354/power_law_dist_post_no_kext.txt`,
  `WGD2038-4008/desj2038_pl_nokext_nokin_dt_weight.csv` (weighted), `SDSS1206+4332/angular_diameter_pre_LOS_power_law.txt`
  (pickle); κ_ext
  files beside them (`*kext*`, `kappa_powerlaw_*.dat`, `kappahist_*`); kinematic likelihoods `*_const_processed.pkl`.
- FRBs/FRB repository, https://github.com/FRBs/FRB (commit 996fcda, 2026-05-06, BSD-3):
  `frb/data/Galaxies/public_hosts.csv` and `frb/data/FRBs/FRB*.json` (DM, DMISM). `DMISM` is NE2001 (Cordes &
  Lazio 2002, arXiv:astro-ph/0207156) from `frb/mw.py` `ismDM` (python `ne2001` package, `ElectronDensity().DM(l, b,
  100.)`), set by `frb/builds/build_frbs.py`; used as stored, not recomputed. Macquart et al. 2020, Nature 581,
  391 (arXiv:2005.13161) for ⟨DM_cosmic⟩, p(Δ) and the host log-normal; James et al. 2022, MNRAS 516, 4862
  (arXiv:2208.00819) for F ≈ 0.32.
- YMW16 Galactic electron-density model, Yao, Manchester & Wang 2017, ApJ 835, 29 (arXiv:1610.09448), through pygedm
  3.3.0 (https://pypi.org/project/pygedm/3.3.0/, sdist; Price, Flynn & Deller 2021, PASA 38, e038, arXiv:2106.15816):
  only its compiled `ymw16` extension and parameter files are used (D-073 addendum). Accessed 2026-10-09.

## Neutrino Frontier round-1 review (accessed 2026-10-10; docs/neutrino_frontier/round1_review.md)
- Stuttard 2021, "Neutrino signals of lightcone fluctuations resulting from fluctuating space-time", arXiv:2103.15313
- "Probing Lorentz Violation in Neutrino Propagation from a Core-Collapse Supernova", arXiv:1110.4848
- Perlman et al. 2015, "New Constraints on Quantum Gravity from X-ray and Gamma-Ray Observations", ApJ 805, 10, arXiv:1411.7262
- Beacom et al. 2004, "Pseudo-Dirac Neutrinos, a Challenge for Neutrino Telescopes", arXiv:hep-ph/0307151
- Rink & Sen 2022, "Constraints on pseudo-Dirac neutrinos using high-energy neutrinos from NGC 1068", arXiv:2211.16520
- Dixit, Miranda & Razzaque 2024, "Searching for Pseudo-Dirac neutrinos from Astrophysical sources in IceCube data", arXiv:2406.06476
- Martinez-Soler, Perez-Gonzalez & Sen 2022, "SN1987A still shining: A Quest for Pseudo-Dirac Neutrinos", PRD 105, 095019, arXiv:2105.12736
- Sen 2022, "Constraining pseudo-Dirac neutrinos from a galactic core-collapse supernova", arXiv:2205.13291
- Longo 1988, PRL 60, 173, doi:10.1103/PhysRevLett.60.173 (SN1987A nu-gamma Shapiro delay)
- Krauss & Tremaine 1988, PRL 60, 176, doi:10.1103/PhysRevLett.60.176 (SN1987A WEP test)
- Boran, Desai & Kahya 2019, "Constraints on differential Shapiro delay between neutrinos and photons from IceCube-170922A", EPJC 79, 185, arXiv:1807.05201
- "Multimessenger Tests of Einstein's Weak Equivalence Principle and Lorentz Invariance with a High-energy Neutrino from a Flaring Blazar", arXiv:1807.06504
- Minazzoli, Johnson-McDaniel & Sakellariadou 2019, "Shortcomings of Shapiro delay-based tests of the equivalence principle on cosmological scales", arXiv:1907.12453; Moriond summary arXiv:2203.11215
- LVK 2017, "Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A", arXiv:1710.05834
- IceCube 2018, "Neutrino emission from the direction of the blazar TXS 0506+056 prior to the IceCube-170922A alert", Science 361, 147, arXiv:1807.08794
- IceCube et al. 2018, "Multi-messenger observations of a flaring blazar coincident with high-energy neutrino IceCube-170922A", arXiv:1807.08816
- IceCube 2017, "Extending the search for muon neutrinos coincident with gamma-ray bursts in IceCube data", ApJ 843, 112, arXiv:1702.06868
- IceCube 2021, "Every Flare, Everywhere: An All-Sky Untriggered Search for Astrophysical Neutrino Transients Using IceCube Data", ICRC2021, arXiv:2107.12134
- IceCube 2021, "IceCube Data for Neutrino Point-Source Searches Years 2008-2018", arXiv:2101.09836
- IceCube 2023, "IceCat-1: the IceCube Event Catalog of Alert Tracks", arXiv:2304.01174
- Stein et al. 2021, "A tidal disruption event coincident with a high-energy neutrino", arXiv:2005.05340
- Päs, Pakvasa & Weiler 2005, "Sterile-active neutrino oscillations and shortcuts in the extra dimension", PRD 72, 095017, arXiv:hep-ph/0504096
- 2024, "The neutrino flavor oscillations in the static and spherically symmetric black-hole-like wormholes", arXiv:2412.02144
- LVK 2023, "Search for gravitational-lensing signatures in the full third observing run of the LIGO-Virgo network", arXiv:2304.08393

## Neutrino Frontier data audit (accessed 2026-10-10; docs/neutrino_frontier/data_audit.md)
- IceCube IceTracks-DR2, Harvard Dataverse doi:10.7910/DVN/MMIIZA (v3.1, 2026-10-05; CC0), 2008-04-06 to
  2022-05-23, 1,643,355 track events; paper arXiv:2605.19040. Events and uptime files (28, ~210 MB) streamed
  2026-10-10 via `https://dataverse.harvard.edu/api/access/datafile/<id>`; ids and sha256 of the served `.tab`
  in `data/manifests/nf_icetracks_dr2.ecsv` (Dataverse's md5 is of the original CSV, not the served file).
- IceTracks-DR2 IC86 IRFs and IC86_IV files for the TXS benchmark (D-080): Dataverse version 1.0 IRFs (file ids
  13597927, 13597978) and version 3.1 IRFs (14153506, 14153513), original CSVs; ids and md5 in
  `data/manifests/nf_skyllh_dr2.ecsv` (accessed 2026-10-10). Benchmark values: arXiv:2605.19040 Table 6 and section 5.
- IceCube IceTracks-DR1 (10-yr PS), doi:10.7910/DVN/VKL316 (v2.0) / data DOI 10.21234/CPKQ-K003, arXiv:2101.09836;
  same events as HEASARC TAP table `icecubepsc` (1,134,450 rows).
- IceCube HESE 12-yr DirectFit, doi:10.7910/DVN/PZNO2T (v2.0; `data.tab`, 164 events), PoS(ICRC2023)1030.
- IceCube HESE 7.5-yr, doi:10.21234/4EQJ-BB17 (zip 78,998,573 B), PRD 104, 022002 (arXiv:2011.03545).
- IceCube flavor composition 11.4 yr, doi:10.7910/DVN/CBNMEB (v2.0), arXiv:2510.24957.
- IceCube TXS 0506+056 2008-2017 events, doi:10.21234/B4QG92 (zip 31,520 B), Science 361, 147
  (doi:10.1126/science.aat2890).
- IceCube GW O3 joint search replication data, doi:10.7910/DVN/34B5AP (v1.0), arXiv:2601.07595.
- IceCat-2 preliminary (7 events), doi:10.7910/DVN/RX28YT (v1.0); proceedings arXiv:2507.06176.
- GCN AMON IceCube gold/bronze and cascade tables, https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html,
  https://gcn.gsfc.nasa.gov/amon_icecube_cascade_events.html (live pages; snapshot sha256 to be pinned on use).
- GCN Circulars, https://gcn.nasa.gov/circulars/<id>.json and archive.json.tar.gz (31.4 MB, daily).
- ANTARES 2007-2017 point-source tracks, https://opendata.km3net.de dataset "ANTARES 2007-2017 Point Source
  Analysis" v1.1 (DOI 10.5072/FK2/HZQTC5 is a DataCite test prefix, not persistent; CC BY 4.0).
- KM3-230213A event data, https://opendata.km3net.de v1.0 (test-prefix DOI 10.5072/FK2/JW72C9); Nature 638, 376
  (doi:10.1038/s41586-024-08543-1).
- Swift GRB Burst Advocate compilation, HEASARC TAP table `swiftgrbba` (2,043 rows on 2026-10-10).
- SkyLLH, PyPI `skyllh` 26.1.0 (2026-09-09), https://github.com/icecube/skyllh.
- SN1987A: Hirata+ 1988 PRD 38, 448; Bionta+ 1987 PRL 58, 1494; Bratton+ 1988 PRD 37, 3361; Alexeyev+ 1988
  PLB 205, 209; compilation Loredo & Lamb 2002 PRD 65, 063002 (astro-ph/0107260).

## Trial unit X1 feasibility (accessed 2026-10-12; docs/state/units/X1.md)
- GRAVITY 2018-07-22 flare astrometry, CDS J/A+A/618/L10 (https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/618/L10;
  `avg.dat`, `wais.dat`, `pfuhl.dat`). There is no CDS catalogue for A&A 677, L10 (ReadMe 404 on 2026-10-12).
- ESO pipeline gravity-kit 1.9.7, https://ftp.eso.org/pub/dfs/pipelines/instruments/gravity/; raw data in programmes
  0101.B-0576(C,E), 105.20B2.004 and 109.22ZA.002 (ESO archive).
- Gyoto, https://github.com/gyoto/Gyoto (GPL-3, tag 2.0.2); einsteinpy (MIT), https://github.com/einsteinpy/einsteinpy.

## Trial unit X2 screen (accessed 2026-10-13; docs/state/units/S0-2026-10-13.md, X2.md)
- GrayHawk-v2, Calzà, arXiv:2606.06339; https://github.com/marcocalza89/GrayHawk-v2 (commit 919126d; no licence file;
  Mathematica). GrayHawk v1: arXiv:2502.04041.
- Miani et al. 2023, "Constraints on the amplitude of gravitational wave echoes from black hole ring-down using
  minimal assumptions", arXiv:2302.12158 (Table I: per-event h_rss at 50 % DP, FAP 5 %).
- Lütfüoğlu 2026, arXiv:2606.08351 (black-bounce Hawking emission; not two-sided throat transmission).
- Krasnikov, arXiv:1910.07839; Dai & Stojkovic, arXiv:1912.07793 (comment and reply on arXiv:1910.00429).

## Trial unit X2 S1–S3 (accessed 2026-10-13; docs/state/units/X2.md)
- Miani et al., arXiv:2302.12158v2, Table II (95 % upper limits on echo h_rss and merger-ringdown h_rss), §III.5,
  §IV.1, §IV.4.
- Churilova & Stuchlík, arXiv:1911.11823, Eq. 6 (SV test-field potentials, s = 0, 1).
- GWTC-1, arXiv:1811.12907 (GW150914 final mass 63.1 M☉, z = 0.09, χ_f = 0.69); Berti, Cardoso & Will,
  arXiv:gr-qc/0512160 (Kerr QNM fits).

## Trial unit X2 S4 (accessed 2026-10-13; docs/state/units/X2_s4_matter.md)
- Yang & Lan, arXiv:2603.20594v1 (2026-03-21): axial gravitational perturbations of SV under anisotropic-fluid
  and NED + scalar sources; master equations Eqs. 33, 38–39, potentials Eqs. 55–56; supplementary QNM tables
  (wormhole branch a = 2.13–10 M; a = 2.20: AF 0.308238 − 0.000266i, NED GR-led 0.246140 − 0.000111i).
- Bronnikov & Walia, Phys. Rev. D 105, 044039 (2022), doi:10.1103/PhysRevD.105.044039, arXiv:2112.13198:
  SV sourced by a phantom scalar plus magnetic NED; L(F) for Schwarzschild-SV, Eq. L-fin.
- Nagar & Rezzolla, CQG 22 R167 (2005), arXiv:gr-qc/0502064v5: odd-parity Moncrief function Q^(o), CPM function
  Ψ^(o), ∂ₜΨ^(o) = −Q^(o) + source (Eq. Q_vs_Psi), h₊,× = −Ψ^(o)/r (§6.1), dE/dt ∝ |Ψ̇^(o)|² (Eq. power).
- Polar / radial stability survey for S4c (titles checked on arxiv.org 2026-10-13): 2607.10718 (Bronnikov et al.,
  NED in scalar-tensor gravity, ℓ = 0), 1205.2224 (Bronnikov, Konoplya & Zhidenko, PRD 86 024028), 0806.0608
  (Gonzalez, Guzmán & Sarbach, CQG 26 015010), 1911.13103 (Cremona, Pizzocchero & Sarbach, PRD 101 104061),
  2507.02258 (Wu, Wang & Miao), 2510.11406 and 2607.18399 (charged wormholes, radial and nonradial), 2608.20081
  (Alencar et al., de Sitter black bounces, axial).

## Trial unit X2 S6 rotation screen (accessed 2026-10-14; docs/state/units/X2_s6_rotation.md)
- Azad, Blázquez-Salcedo, Khoo & Kunz: arXiv:2301.05243v2 (PLB 848, 138349, doi:10.1016/j.physletb.2023.138349);
  arXiv:2403.08387v1 (PRD 109, 124051), Tables I-II, §V-VI.
- Blázquez-Salcedo et al., arXiv:2510.11406 v1 (13 Oct 2025, "Charged wormholes can be long-lived") and v2
  (13 Mar 2026; PLB 875, 140335).
- Dzhunushaliev, Folomeev, Kleihaus, Kunz & Radu, arXiv:1309.2448 (PRD 88, 124028), §V (perturbation ODEs, J_cr).
  X2 S6a used the arXiv source (accessed 2026-10-14): Eqs. (lineel), (eqa), (eqq), (Rll-expr), (pertlineel), (hop),
  (eqv00)-(eqv22), and values digitised from the gnuplot EPS of Figs. 1 and 5 (`om_vs_Jn.eps`, `M_vs_Jn.eps`).
  Citing literature: INSPIRE `refersto:recid:1253547` (47 records, 2026-10-14); erratum check via Crossref
  (10.1103/PhysRevD.88.124028, no update records). Log in `docs/state/units/X2_s6a_prior_work.md`.
- Kleihaus & Kunz, arXiv:1409.1503 (PRD 90, 121503); Chew, Kleihaus & Kunz, arXiv:1608.05253 (PRD 94, 104031).
- Khoo et al., arXiv:2401.02898 (PRD 109, 084013); Blázquez-Salcedo et al., arXiv:2312.10754 (PRD 109, 064028);
  Azad et al., arXiv:2509.22118 (Universe 11, 325); arXiv:2603.12748; arXiv:2610.13206 (not refereed).

## Trial unit X3 and S0-2026-10-14 (accessed 2026-10-14; docs/state/units/S0-2026-10-14.md, X3.md)
- Gaia DR3 Andromeda Photometric Survey (GAPS): `gaiadr3.gaia_source.in_andromeda_survey` (1,257,319 sources, all
  `has_epoch_photometry`). Gaia archive TAP https://gea.esac.esa.int/tap-server/tap. Epoch photometry via DataLink
  https://gea.esac.esa.int/data-server/data (`DATA_STRUCTURE=DATAMODEL_STANDARD`; `INDIVIDUAL` returns HTTP 500
  since at least 2026-10-14).
- Gaia DR3 bulk epoch photometry (accessed 2026-10-14): https://cdn.gea.esac.esa.int/Gaia/gdr3/Photometry/epoch_photometry/
  (listing via the S3 XML at https://gaia.eu-1.cdn77-storage.com/?prefix=Gaia/gdr3/Photometry/epoch_photometry/&delimiter=/;
  `_MD5SUM.txt`). Files named by level-8 nested HEALPix range. The five covering GAPS and their MD5s are in
  `scripts/x3_bulk.py` `FILES`. No `rejected_by_photometry` column (rebuilt from `g_transit_n_obs`).
- Reference variables: `gaiadr3.vari_classifier_result` (ECL, YSO, LPV; score > 0.5) in the M31 and Orion fields;
  SIMBAD and Gaia DR3 (VizieR I/355/gaiadr3) via CDS XMatch for vetting.
- X4 S1 tools (2026-10-15): lal 7.7.1 / lalsimulation 6.2.1 from PyPI (GPL-2.0-or-later; used in a separate
  environment, not a project dependency); IMRPhenomXAS and IMRPhenomD; aLIGO design PSD
  (correction 2026-10-21: "7.7.1" is the lal version; the PyPI package providing it is lalsuite 7.26.16, used again
  for S1g);
  `SimNoisePSDaLIGODesignSensitivityP1200087` (LIGO-P1200087).
- S0-2026-10-15b Kapitza screen (accessed 2026-10-15; read in the body, versions as listed):
  - Bukov, D'Alessio & Polkovnikov, arXiv:1407.4803v6, §5, Eqs. 92–106.
  - Krondorfer et al., arXiv:2507.02736v2, Eqs. 12–39.
  - Berenguer et al., "Floquet SYK wormholes", arXiv:2404.08394v3.
  - Journal references were not checked.
  - Gonzalez, Guzman & Sarbach, arXiv:0806.0608v2, Table I (Ellis unstable time 0.846 r_throat).
  - Azad et al., arXiv:2301.05243v2, Table 1; arXiv:2403.08387v1; Kleihaus & Kunz, arXiv:1409.1503v1.
  - Dzhunushaliev et al., arXiv:1309.2448v1; Friedman, Schleich & Witt, gr-qc/9305017v2, Theorem 1.
  - Poisson & Visser, gr-qc/9506083 (thin-shell stability; used for the toy model).
  - Not re-read this round: Isaacson 1968 (Phys. Rev. 166, 1263/1272); MTW §35.15; Brill & Hartle 1964.
  - Within four web queries only (dynamic/Kapitza/Floquet stabilisation of wormholes; Isaacson or Brill–Hartle
    averaging applied to wormholes), no paper was found. This does not mean none exists.
- S0-2026-10-16 unstable-wormhole message window (accessed 2026-10-16; body read unless noted):
  - Shinkai & Hayward, arXiv:gr-qc/0205041v2, PRD 66, 044005 (2002), doi:10.1103/PhysRevD.66.044005.
  - Shinkai, `wh_dual` program page https://www.oit.ac.jp/labs/is/system/shinkai/his/programs/wh_dual.html
    (archive `whdual_020430.tar.Z`, page updated 2007-04-19; no licence stated; archive not downloaded).
  - González, Guzmán & Sarbach, arXiv:0806.1370v2, CQG 26, 015011 (2009); with arXiv:0806.0608v2 (listed above).
  - Doroshkevich, Hansen, Novikov & Shatskiy, arXiv:0812.0702v2, IJMPD 18, 1665 (2009), doi:10.1142/S0218271809015230.
  - Calhoun, Fay & Kain, arXiv:2210.04905v2, PRD 106, 104054 (2022).
  - Xu, Chew & Yeom, arXiv:2503.07610v1 (2025).
  - Abstract only: Shirokov, arXiv:2604.00071v1; Koyama & Hayward, gr-qc/0406113; Hayward & Koyama, gr-qc/0406080.
- X4 S1d reuse check (accessed 2026-10-15):
  - SciPy 1.18.1 `optimize.differential_evolution` (BSD-3-Clause):
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html
  - Bilby 2.8.2 (MIT), https://pypi.org/project/bilby/
  - PyCBC 2.11.0 (GPLv3), https://pypi.org/project/pycbc/; its filter API, https://pycbc.org/pycbc/latest/html/pycbc.filter.html
  - cma 4.5.0 (BSD-3-Clause), https://pypi.org/project/cma/ (not installed)
  - The PyPI package "glow" is unrelated to GLoW.
- S0-2026-10-15 / X4 (accessed 2026-10-15):
  - LVK lensing searches: O3 arXiv:2304.08393; GWTC-4.0 arXiv:2512.16347 (isolated point mass; outlier
    GW231123_135430).
  - GW231123 lensing follow-ups: arXiv:2512.17631, 2512.19118, 2607.21834, 2604.08179.
  - Wave optics: Yoo, Harada & Tsukamoto arXiv:1302.7170 (Ellis wormhole, EM); GLoW arXiv:2409.04606 (licence not
    checked); wolensing arXiv:2410.19804.
  - aLIGO PSD fit attributed to Ajith, PRD 84, 084037 (2011), arXiv:1107.1267 (formula not checked in the body).
  - VASCO: arXiv:2507.15896, 2604.04810.
  - Quantum-processor wormhole: Jafferis et al., Nature 612, 51 (2022); Kobrin, Schuster & Yao comment
    arXiv:2302.07897; PhilSci-Archive 22120 (arXiv:2301.03522).
- X3 DR4 reopen check (accessed 2026-10-16):
  - https://www.cosmos.esa.int/web/gaia/dr4 ("Last update: 14 September 2026"; `epoch_astrometry`
    129,625,750,708 rows).
  - Draft data model:
    https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/gaia-dr4-prerelease-draft-data-model_2026-06-26.zip
    (last-modified 2026-06-26; PDF dated 2026-06-16; SHA-256 d807eae98acbeec0a0af2ce6a5d7d352298df6f270a1f13207cc8d1becf42c66).
    §3.1, §6.2 and §8.2 were read.
  - https://www.cosmos.esa.int/web/gaia/dr4-prerelease (no page date).
  - https://www.cosmos.esa.int/web/gaia/iow_20260703 (published 2026-07-03).
- Gaia DR4 content (accessed 2026-10-14): https://www.cosmos.esa.int/web/gaia/dr4 (`epoch_astrometry` for all
  sources in the main astrometric processing; 66 months, 2014-07-25 to 2020-01-20); image of the week
  https://www.cosmos.esa.int/web/gaia/iow_20260703 (2.8 billion sources, 142 billion detections).
- Gaia DR4 release date, 2 December 2026: https://www.cosmos.esa.int/web/gaia/release. Epoch-photometry coverage is
  not yet stated on https://www.cosmos.esa.int/web/gaia/data-release-4.
- PTA (S0 A):
  - NANOGrav 15-yr new-physics search, arXiv:2306.16219 (Sec. 6.2, Table 3);
  - Dror et al. arXiv:1901.04490;
  - Lee et al. arXiv:2104.05717;
  - Cherukupalli et al. arXiv:2607.03533;
  - Foster et al. arXiv:2606.18329;
  - Nakajima, Izumi & Asada arXiv:1404.2720 (negative time delay);
  - NANOGrav 15-yr data, Zenodo 10.5281/zenodo.16051178 (v2.1.0, CC-BY-4.0, not downloaded).
- FRB lensing (S0 B):
  - arXiv:2204.06001, 2204.06014, 2109.09251, 2112.03721, 2206.14310, 2002.12533, 1605.00008;
  - Yoo, Harada & Tsukamoto arXiv:1302.7170;
  - Takahashi & Asada arXiv:1303.1301.
- MACS0416 transients (accessed 2026-10-15): Yan et al. 2023, "JWST's PEARLS: Transients in the MACS J0416.1-2403
  Field", ApJS 269, 43, DOI 10.3847/1538-4365/ad0298, arXiv:2307.07579 (v2 source, Tables 4–5; no VizieR table).
  PEARLS NIRCam epochs: JWST program 1176 o211/o212/o213 (MAST), CANUCS 1208 o004.
- Abell 2744 two-epoch grid (accessed 2026-10-15): JWST DD program 2756 (PI W. Chen) NIRCam `jw02756-o003_t001`
  (2022-10-20) and `jw02756-o002_t001` (2022-12-06), level-3 `_i2d`/`_cat.ecsv` and level-2 `_cal` from
  `s3://stpubdata/jwst/public/jw02756/` (jwst 2.0.1, CRDS jwst_1535.pmap); UNCOVER `jw02561-o001_t003` byte-range
  stamps. Proposal: https://www.stsci.edu/jwst-program-info/download/jwst/pdf/2756/ (z = 3.47 SN, no public position).
  Pang et al. 2026, arXiv:2608.20781 v1 (submitted 2026-08-21; LaTeX source accessed 2026-10-15): PLATE-23a, a SN IIP
  candidate in Abell 2744. The printed "03h33m09.9216s, −30d24m01.2229s" is degrees in sexagesimal form:
  RA 3.552756°, Dec −30.400340° (verified on the ALT and UNCOVER F356W images, 2026-10-15).
- ALT (All the Little Things, JWST GO 3516; Naidu et al. arXiv:2410.01874) NIRCam direct imaging
  `jw03516-o001_t001` (2023-12-06/10) and `jw03516-o002_t001` (2023-11-27/12-05), `clear-f356w` `_i2d`/`_cat.ecsv` from
  `s3://stpubdata/jwst/public/jw03516/L3/t/` (jwst 3.0.0, CRDS jwst_1584.pmap; accessed 2026-10-15).

## HFF Abell S1063 second models (accessed 2026-10-16; `scripts/arcs_second_models.py`)
- GLAFIC v4, Sharon v4, Keeton v4 `hlsp_frontier_model_abells1063_<team>_v4_{kappa,gamma1,gamma2}.fits` under
  https://archive.stsci.edu/pub/hlsp/frontier/abells1063/models/<team>/v4/ (each map ~105.8 MB; GLAFIC kappa
  105,756,480 B). Read by HTTP byte ranges at 24 positions only; no whole-file sha256 (not downloaded).
- Chan et al. 2020, ApJ 888, 35 — redshifts of lensed systems toward RXC J2248.7−4431 (Abell S1063), VizieR
  `J/ApJ/888/35` table1 (121 rows), queried 2026-10-16 for the `arcs` anti positions. Streamed again 2026-10-21 as the D-035
  image list (`scripts/s1063_imagelist.py`, `results/lens_consistency/abells1063_images/chan20_images.ecsv`).
- ESO Science Archive, MUSE reduced cubes towards Abell S1063, accessed 2026-10-16 (`scripts/arcs_vet_spectra.py`):
  ObsCore TAP https://archive.eso.org/tap_obs (`ivoa.ObsCore`, instrument MUSE, cubes) and SODA sync cutouts
  https://dataportal.eso.org/dataPortal/soda/sync. Products used: ADP.2017-03-28T12:46:01.331 (095.A-0653, 5.78 GB,
  not downloaded) and ADP.2017-03-23T15:58:03.937 (60.A-9345, 5.73 GB, not downloaded), as 12–17 MB cutouts that
  were deleted after use. The coverage of every listed public cube is in
  `results/arcs/abells1063-cats/vetting/muse_coverage.ecsv`.
- ESO Science Archive, MUSE combined cube towards MACS J1931.8−2635, ADP.2017-03-27T15:26:07.707 (095.A-0525(A),
  8,127 s, ≈ 3.0 GB, not downloaded), accessed 2026-10-24 (`scripts/m1931_muse.py`): ObsCore TAP
  https://archive.eso.org/tap_obs/sync and SODA sync `BAND` cutouts (four × 105,131,520 B, deleted after use; URLs
  and DATA-HDU sha256 in `results/lens_consistency/macs1931_muse/*/summary.json`).
- Gaia DR3 `gaiadr3.gaia_source` box over that cube's footprint, ESA Gaia TAP https://gea.esac.esa.int/tap-server/tap,
  queried 2026-10-24 (10 rows, `results/lens_consistency/macs1931_muse/gaia_dr3.csv`), for the MUSE registration.
- GLIMPSE 3293 o001 F090W and F115W pipeline catalogues (`jw03293-o001_t001_nircam_clear-{f090w,f115w}_cat.ecsv`;
  sha256 `3eecb26c…`, `3775f1f9…`, which match `data/manifests/abells1063.ecsv`), read 2026-10-16 for continuum
  magnitudes.
- JWST VENUS (program 6882) observation 51 at El Gordo, 2026-08-05: `jw06882-o051_t051_nircam_clear-{f150w,f444w}_i2d.fits`
  and `_cat.ecsv` (jwst 3.0.0, CRDS jwst_1584.pmap), public level 3, read by S3 byte range on 2026-10-16 with PEARLS
  `jw01176-o241_t012_nircam_clear-{f150w,f444w}` (`results/transients/elgordo_arcdiff/`).

## Difference-imaging tool comparison (accessed 2026-10-17; `results/transients/sfft_trial/`)
- SFFT: Hu, L. et al. 2022, ApJ 936, 157, "Image Subtraction in Fourier Space". Code: https://github.com/thomasvrussell/sfft
  (MIT), PyPI `sfft` 1.7.3. Used: `CustomizedPacket` (Numpy backend).
- Hu, L. & Wang, L. 2024, AJ 167, 5, "Differencing and Coadding JWST Images with Matched Point Spread Function",
  arXiv:2309.09143, doi:10.3847/1538-3881/ad36cb. The B-spline method is `sfft.BSplineSFFT`, CuPy (GPU) only in 1.7.3;
  not run.
- MAST CAOM query (astroquery, 2026-10-17): VENUS 6882 NIRCam F150W level-3 targets and other programs within ±0.05°
  (`venus_pairs.json`).
- acstools 3.8.2 `findsat_mrt.TrailFinder` (STScI; BSD-3-Clause; https://github.com/spacetelescope/acstools,
  accessed 2026-10-17): median-Radon-transform satellite-trail finder, Stark et al. 2022, ACS ISR 2022-8
  (https://ui.adsabs.harvard.edu/abs/2022acs..rept....8S/abstract). Used with scikit-image 0.26.0 and photutils 3.1.1
  for the W1 streak masks (`scripts/w1_streak_mask.py`, `results/w1_streak_veto/`).

## W1 shear screen in the CANUCS flanking fields (accessed 2026-10-18; `results/w1_shear_ncf/`)
- JWST CANUCS 1208 F150W pipeline catalogues `jw01208-<obs>_nircam_clear-f150w_cat.ecsv`, public S3
  `s3://stpubdata/jwst/public/jw01208/L3/t/<o>/` (the NCF observations have no F200W); sha256:
  - o023_t002 (macs0416-ncf) `6ac987242db0435afb0aac549b825d12a7f67868e93cf23e63aef0dada8ec38c`
  - o029_t004 (macs1149-ncf) `e331b754f054cf3c31eb84096ff7dd3131b2aeb6f595f711953244ddad2c672d`
  - o020_t001 (abell370-ncf) `6be0bacd8cd437bbf557076a099888770c301ea2b81b835ddde84f12b7b9328c`
  - o026_t003 (macs0417-ncf) `b3d71250aec8bbe965c1c4908958e5d6366af7c8e25abf587006de65cee8383e`
  - o032_t005 (macs1423-ncf) `fb2f118625ae5a4a8c1107e546c5c19dfb7f2c4d1125d4c24d003ff4677898d2`
  - gate R cluster pointings: o006_t003 (MACS0417) `1d49cb6c579b0f2ae4923681a2be2b82db5fed8ebdb430fa64cf851362ccb2ac`,
    o010_t005 (MACS1423) `f165280389d6fbc192dd08c00208d9f4ed19f9f8f33875596b33706e611a4e6d`
  - the `_i2d.fits` of the same names were read by S3 byte range (streak masks, cutouts), not stored.
- CANUCS DR1 NCF photometry and EAzY photo-z: as `orphan_pairs.DEEP_FIELDS` (all five sha256 matched on download).

## W1 shear screen in GOODS-N (accessed 2026-10-18; `results/w1_shear_goodsn/`)
- JWST JADES 1181 F200W pipeline catalogues `jw01181-<obs>_nircam_clear-f200w_cat.ecsv`, public S3
  `s3://stpubdata/jwst/public/jw01181/L3/t/<o>/`, chosen from the MAST `s_region` footprints; sha256:
  - o001_t001 `28a3c404bd9a8d21828a23e4a9a34c0c56668b3d68678b9706a56af337a6c8e0`
  - o003_t007 `0563c8be132f3e9c147b54327389645f1df952a3b2498bd496b4abe1dda8e520`
  - o004_t008 `35bd5e3e7d21ef9a083d83e16f57ebde2c2f9be4ca698026cee0fe6a323cea06`
  - o005_t008 `27047942ed8fea19ef2a5466e0d5f26990ba6e1a09e9e3b1c50a29bc1b6b6ff1`
  - o009_t009 `818eb57775777d24dfba06a1cddb399c0a559aa59d399d8e387b7c619180f9b4` (gate P failed)
  - o010_t009 `a9060f9ce85856b73c4a4977a98ab58beec74e4936d36b24aae69c31b7e78574` (gate P failed)
- Photo-z: DJA v7.3 `gdn-grizli-v7.3-fix.eazypy.zout.fits`, streamed out of `gdn-grizli-v7.3-fix.photoz.tar.gz`
  (pins in `orphan_pairs.DEEP_FIELDS["goodsn-dja"]`; D-051). The F200W `_i2d.fits` SCI extensions were range-read
  for the streak masks and the sheet (not stored).

## W1 shear screen in CEERS (accessed 2026-10-18; `results/w1_shear_ceers/`)
- JWST CEERS (program 1345) F200W pipeline catalogues `jw01345-<obs>_nircam_clear-f200w_cat.ecsv`, public S3
  `s3://stpubdata/jwst/public/jw01345/L3/t/<o>/`, chosen from the MAST `s_region` footprints
  (`scripts/w1_footprints.py`; o052_t022 left out); sha256:
  - o001_t021 `878766359ac425b9510b0e02dba82f032f23a1e68d2deb097d51df9cc5f35f27`
  - o002_t022 `152648c84958a8bdc5e72287b55b70e3af292a0d0404917bf3fe65e325a11999`
  - o003_t023 `9a989e0c12144c65a332abab6694839a4ad5e8e3b00c28438fad957a573fc9d4`
  - o004_t024 `099d0beee716e6ba566a593908e1361e671d1fcf93d8957d435257f3b64f5dfa`
  - o061_t043 `b2c1b51ae6f8c3019716a7db1187560fe698852aef449ec5a9694d7aabba32b3`
  - o064_t044 `c3b2fce5200b0e6b690237d314c751687a38bfd3db029780d1d89411d8b9c5ff`
  - o066_t055 `93a96a25e0d9c8177044a2334cde57668759796f6f618b6e17c7d1a64bd5ba24`
  - o068_t048 `a897343adde44740e9031dd1cbde0340135a9af920d09f9caff9ab58800e919c`
  - o069_t046 `c8789ba9f7658db262c6432818990af0bf0318d5ed463818529060d984f28422`
  - o071_t047 `c775b3cd3dcaaf3698bb4a1c6b97dc0696aaf63825836537d8833a82c4f1ff93`
- Photo-z: DJA v7.4 `ceers-full-grizli-v7.4-fix.eazypy.zout.fits` (sha256
  `c50443d242c4c76b2d3c164b236911f14d926bf5c466ae90bfe0def299a48a50`, 78.4 MB), streamed out of
  https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/ceers-full-grizli-v7.4-fix.photoz.tar.gz (396.8 MB, Last-Modified
  2024-09-27; tarball sha256 not recorded; same DJA terms and citation, arXiv:2302.10936). The F200W `_i2d.fits` SCI
  extensions were range-read for the streak masks, the vetting and the sheet (not stored).

## W1 blank-field R simulation (accessed 2026-10-18; `results/w1_r_metacal/`)
- **GalSim** 2.8.5 (PyPI `galsim`), https://github.com/GalSim-developers/GalSim; Rowe et al. 2015, Astronomy and
  Computing 10, 121, doi:10.1016/j.ascom.2015.02.002. Sérsic profiles, Airy ⊗ Gaussian PSF, shears.
- **photutils** 3.1.1 (PyPI), `detect_sources` and `SourceCatalog` (`semimajor_axis` is 3.x's name for the
  catalogue's `semimajor_sigma`, per its deprecation warning).
## W1 shear in PRIMER (accessed 2026-10-18; `results/w1_shear_primer/`)
- **PRIMER** (JWST GO 1837, PI Dunlop) NIRCam F200W level-3 products, MAST: `jw01837-c1001_t000`, `c1002_t000`
  (UDS combined mosaics), `o003_t003`, `o004_t004` (COSMOS) `_cat.ecsv`; sha256 pins in `scripts/w1_stage_inputs.py`.
  Footprints from the MAST CAOM `s_region` (`footprints_*.json`). The `_i2d.fits` files are 31–69 GB (HTTP
  Content-Length, not stored).
- **DJA v7** photo-z tarballs for PRIMER, streamed 2026-10-18 (only the `.eazypy.zout.fits` member kept):
  `primer-uds-north-grizli-v7.2-fix`, `primer-uds-south-grizli-v7.2-fix`, `primer-cosmos-east-grizli-v7.4-fix`,
  `primer-cosmos-west-grizli-v7.4-fix` (`https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/<stem>.photoz.tar.gz`;
  the names without `-fix` return 403). Tarball and member sha256, and the content sha256 of the two stacked
  per-survey zouts: `PRIMER_DJA` / `PRIMER_STACK_SHA256` in `scripts/w1_stage_inputs.py`.
- The F200W `_i2d.fits` SCI extensions are 1.9–3.3 GB (NAXIS 25047 × 32984 for c1001; headers read
  2026-10-18, not stored); `w1_streak_mask.py` streams only that extension.
- Prior-work check 2026-10-18 (`docs/value_comparison.md`): DeCoursey et al. 2024, JADES Transient Survey,
  arXiv:2406.05060; Stein et al. 2021, "Self-supervised similarity search for large scientific datasets",
  arXiv:2110.13151; https://github.com/georgestein/ssl-legacysurvey (MIT); Palencia et al. 2026 catalogue,
  https://doi.org/10.5281/zenodo.19439926. sncosmo 2.x (BSD-3, https://github.com/sncosmo/sncosmo; isolated venv,
  templates hsiao, nugent-sn2p/sn2l/sn1bc and NIRCam bandpasses fetched by sncosmo 2026-10-18) for
  `scripts/achromatic_test.py simulate`; not a project dependency.
- MACS1423 h000 HST archive (queried 2026-10-18; `results/transients/macs1423_h000_hst/hst_obs.ecsv`): HST GO 16667
  (ACS F435W/F606W 2022–2023), GO 17504 (WFC3 F110W/F606W/F814W 2023-12 – 2024-04), CLASH 12790 (2012–2013), 13386
  (WFC3/UVIS F606W 2014), 10493 (ACS F814W 2006); pipeline association `_drc`/`_drz` products read from s3://stpubdata/hst/public/.
- ZTF DR24 light curves and objects, IRSA HATS collection (README version date 2026-05-13), DOI 10.26131/IRSA598;
  s3://ipac-irsa-ztf/ztf/enhanced/dr24/ (anonymous; read 2026-10-18); tutorial
  https://caltech-ipac.github.io/irsa-tutorials/ztf-lightcurves/; release notes
  https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/dr24/ztf_release_notes_dr24.pdf. Used by `scripts/w3_ztf.py`.
- ASAS-SN Sky Patrol v2 (Hart et al. 2023, arXiv:2304.03791; title not re-checked this round), official client
  `skypatrol` 0.6.21 on PyPI (module `pyasassn`; https://github.com/asas-sn/skypatrol), server
  http://asassn-lb01.ifa.hawaii.edu:9006 (deployment 0.6.20). Probed 2026-10-20: catalogue cones work, light-curve
  service down for maintenance (`results/w3_ztf/asassn/probe.json`). Installed in the session venv only; not a
  project dependency.
- W3 ZTF alternative-data inventory (`results/w3_ztf/altdata/`, queried 2026-10-22). The bibliographic details below
  are from memory and were not re-checked this round:
  - **ZTF science-exposure metadata:** IRSA IBE search https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci
    (`POS`, `WHERE filtercode='zg'`).
  - **NEOWISE-R single-exposure source table** `neowiser_p1bs_psd`: IRSA Simple Cone Search
    https://irsa.ipac.caltech.edu/SCS?table=neowiser_p1bs_psd. Epochs to MJD 60512.6 for these stars. Mainzer et al.
    2014, ApJ 792, 30; Explanatory Supplement https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/.
  - **TESS:** MAST CAOM via `astroquery.mast.Observations.query_region` (astroquery 0.4.11): SPOC FFIs and the HLSPs
    QLP, TESS-SPOC and TARS. Ricker et al. 2015, JATIS 1, 014003.
  - **Not queried:**
    - Gaia DR3 epoch photometry spans 2014-07-25 to 2017-05-28 (Gaia Collaboration, Vallenari et al. 2023, A&A 674,
      A1).
    - The Pan-STARRS1 3π survey ran 2010–2014 (Chambers et al. 2016, arXiv:1612.05560).
    - ATLAS forced photometry needs a registered account: https://fallingstar-data.com/forcedphot/ (Shingles et al.
      2021, TNS AstroNote 2021-7).
- **ATLAS forced photometry (used 2026-10-22, `results/w3_ztf/atlas/`):**
  - **Server and documentation:** https://fallingstar-data.com/forcedphot/ — API guide `apiguide/`, FAQ `faq/`, output
    description `resultdesc/`, OpenAPI schema `api/schema/`. The banner on 2026-10-22 read "2026-09-09 07:45 UTC: Normal
    operations have resumed".
  - **Shared access** was documented by NOIRLab's ANTARES devkit notebook "Lightcurves and Photometry Search", Chapter
    3: https://nsf-noirlab.gitlab.io/csdc/antares/devkit/notebooks/LightcurvesPhotometrySearch/ (read 2026-10-22). The
    token is not recorded.
  - **Required citations (server home page):** Tonry et al. 2018 (survey), Smith et al. 2020 (transients), Heinze et
    al. 2018 (variables), Shingles et al. 2021 (the service), plus the ATLAS acknowledgement text. Bibliographic details
    were not re-checked this round.
  - **Quality cuts:** the FAQ's common-sense filters, attributed there to Rest et al. (not re-checked).
- PyCBC 2.11.0 (PyPI, GPLv3; https://github.com/gwastro/pycbc): `pycbc.filter.matchedfilter.compute_max_snr_over_sky_loc_stat_no_phase`
  (Harry et al., arXiv:1709.09181), `pycbc.waveform.get_fd_waveform`, `pycbc.filter.matched_filter`; `bin/pycbc_optimize_snr`
  with `pycbc/live/snr_optimizer.py` (aligned-spin, 4-D); `bin/pycbc_banksim` (finite bank). Checked 2026-10-21 in a separate
  environment (`results/x4/handoff/pycbc_check.json`); not a project dependency.
- Surveyed, not installed (2026-10-21): Bilby relative binning `find_maximum_likelihood_parameters`
  (https://git.ligo.org/lscsoft/bilby/-/blob/master/bilby/gw/likelihood/relative.py); simple-pe `find_peak_snr`
  (https://git.ligo.org/stephen-fairhurst/simple-pe); ripple (JAX; no IMRPhenomXPHM listed). pesummary "maxL" and gwfast
  were assessed from memory, not checked.
- cogwheel 1.6.0 (PyPI `cogwheel-pe`, GPL-3.0-or-later; https://github.com/jroulet/cogwheel; docs
  https://cogwheel.readthedocs.io/en/stable/; paper arXiv:2207.03508, title not re-checked): `Posterior.from_event`,
  `Posterior.refine_reference_waveform`, `IASPrior`, `RelativeBinningLikelihood`, `ReferenceWaveformFinder`, `EventData`;
  read from the installed source 2026-10-21 in a separate environment (`results/x4/handoff/cogwheel_checks.json`,
  `cogwheel_runs.json`); not a project dependency.
- Metacalibration: Huff & Mandelbaum 2017, arXiv:1702.02600 ("Metacalibration: Direct Self-Calibration of Biases in
  Shear Measurement"); Sheldon & Huff 2017, arXiv:1702.02601 ("Practical Weak Lensing Shear Measurement with
  Metacalibration"); titles checked via the arXiv API 2026-10-18. JADES 1181 o001 F200W i2d:
  s3://stpubdata/jwst/public/jw01181/L3/t/o001/jw01181-o001_t001_nircam_clear-f200w_i2d.fits (anonymous, streamed
  2026-10-18). Used by `scripts/w1_r_realmetacal.py`.

## W2 DESI pairs (accessed 2026-10-23; `results/w2_desi_pairs/`)
- DESI DR1 (public release 2025): redshift catalogue `desi_dr1.zpix` via NOIRLab Astro Data Lab TAP
  (https://datalab.noirlab.edu/tap/sync, anonymous); spectra via SPARCL (https://astrosparcl.datalab.noirlab.edu/api,
  anonymous public access, data release `DESI-DR1`; SPARCL ids per spectrum in `fetch_log.json`).
- Hennawi et al. 2006, AJ 131, 1, binary quasars (VizieR J/AJ/131/1/binqso; pinned in
  `data/manifests/w12_niq_inputs.json`, as D-064).
- Route-choice literature (checked by the independent reviewer 2026-10-23):
  - Schmidt et al. 2023, "STRIDES: automated uniform models for 30 quadruply imaged quasars", MNRAS 518, 1260
    (arXiv:2206.04696);
  - Jing et al., "A Quasar Pair Sample Compiled from DESI DR1" (arXiv:2505.03103);
  - Lemon et al. 2023, "Gravitationally lensed quasars in Gaia – IV" (arXiv:2206.07714);
  - Chen et al., "Search for Quasar Pairs with Gaia Astrometric Data. III" (arXiv:2607.00093).

## W2 DESI pairs evidence audit (accessed 2026-10-23; `results/w2_desi_pairs/evidence/`)
- **Lemon et al. 2023** (arXiv:2206.07714): only v1 (2022-06-15), LaTeX source read (§4.2–4.3, table 1); VizieR
  J/MNRAS/520/3305 ReadMe (class notes). The MNRAS text was not read (OUP HTTP 403).
- **Jing et al.** (arXiv:2505.03103 v2): `DR1QP.csv` from https://github.com/astroliang/DESI_DR1_QP (1842 rows, last
  commit 2026-01-22, sha256 f3d05bf1fc8377c99848f88cf708569ecdeed41cae40a6e957fd647b1f789f5e; read, not stored).
- **Chen et al.** (arXiv:2607.00093) §3.1.1 and table 3. Papers I, II and IV of the series: arXiv:2504.17777,
  2605.09450, 2603.03042.
- **SQLS:**
  - Inada et al. 2008 (arXiv:0708.0828);
  - Inada et al. 2010 (arXiv:1005.5570, §4.1–4.2);
  - Inada et al. 2012 (arXiv:1203.1087, §4.1).
- **Binary catalogues:** Hennawi et al. 2006 (astro-ph/0504535, §3.2 and §5.1); Kayo & Oguri 2012 (arXiv:1203.6410).
- **Lens-candidate catalogues:**
  - Dawes et al. 2023 (arXiv:2208.06356; VizieR J/ApJS/269/61/table2; grades within 3″ of each pair stored in
    `identity.json`);
  - Ji et al. (arXiv:2412.08397), whose J0740+2926 'lensing candidate' cites Lemon et al. 2018 (not read).
- **Archives** (metadata only, `archives.json`):
  - MAST CAOM via `astroquery.mast.Observations`; HST GO 17308 (PI Lemon) and GO 17455 (PI Shen), titles and PIs
    stored in `archives.json`;
  - ESO ObsCore TAP https://archive.eso.org/tap_obs/sync;
  - KOA TAP https://koa.ipac.caltech.edu/TAP/sync;
  - CADC ObsCore TAP (HTTP 503, not checked).
- **Identity:** ESA Gaia TAP (https://gea.esac.esa.int/tap-server/tap/sync, `gaiadr3.gaia_source`); NOIRLab Data
  Lab TAP (`desi_dr1.zpix`, `ls_dr10.tractor`).

## W2 J2244−0550 HST test (accessed 2026-10-23; `results/w2_j2244_hst/`)
- **HST products:** pipeline drizzled products (calwf3 / calacs) read as byte-range cutouts from the public AWS copy
  `https://stpubdata.s3.amazonaws.com/hst/public/<ipppss>/<obs_id>/<obs_id>_{drz,drc}.fits`. ETags are in
  `products.json`.
  - Target: `if8a5k030` (GO 17455, PI Shen, WFC3/IR F160W, 2024-07-05) and `jf56c1010` (GO 17308, PI Lemon, ACS/WFC
    F814W, 2024-08-31).
  - Controls are listed in `runs/*.json`.
- **Programme field lists:** MAST CAOM (`astroquery.mast.Observations`), proposals 17455, 17199, 17308 and 9744.
  Star libraries: Gaia DR3 via ESA Gaia TAP (`library_*.json` lists fields and Gaia source ids).
- **Software** (isolated uv environment): lenstronomy 1.14.2 (Birrer & Amara 2018; Birrer et al. 2021, JOSS 6,
  3283), photutils 3.1.1 `EPSFBuilder`, astropy 8.0.2, fsspec 2026.9.0.

## W2 J2244 STScI PSF capability check (accessed 2026-10-23; `results/w2_j2244_stpsf/`)
- **Official entry:** "HST WFC3 Point Spread Function Modeling", hst_notebooks (M. Revalski; page updated 20 May
  2026), https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/point_spread_function/hst_point_spread_function.html
  (sha256 of the page in `sources.json`).
- **PSF model:** WFC3-team STDPSF grid `STDPSF_WFC3IR_F160W.fits` (hosted by J. Anderson; WFC3 ISR 2016-12),
  https://www.stsci.edu/~jayander/HST1PASS/LIB/PSFs/STDPSFs/WFC3IR/ (sha256 in `sources.json`).
- **Data:** 12 GO 17455 WFC3/IR F160W FLT exposures from the public AWS copy
  `https://stpubdata.s3.amazonaws.com/hst/public/<ipppss>/<rootname>/<rootname>_flt.fits` (sha256 in `sources.json`);
  Gaia DR3 positions via ESA Gaia TAP.
- **Software** (isolated uv environment, offline from cache): photutils 3.1.1 (`GriddedPSFModel`, `PSFPhotometry`),
  astropy 8.0.2, numpy 2.5.4, scipy 1.18.1.

## W3 / Kepler dense-cadence pilot (accessed 2026-10-24; `results/w3_kepler/`)
- **Kepler DR25 long-cadence light curves** (SOC 9.3, `DATA_REL` 25; Q1–Q17), MAST:
  https://archive.stsci.edu/missions/kepler/lightcurves/ (per-file sha256 in `files_<group>.json`).
- **Kepler DR25 stellar properties:** NASA Exoplanet Archive TAP table `keplerstellar`, delivery
  `q1_q17_dr25_stellar`, https://exoplanetarchive.ipac.caltech.edu/TAP/sync (query and sha256 in
  `sample.json`).

## W3 / Kepler hard tests (accessed 2026-10-25; `results/w3_kepler/hard/`)
- **Kepler Eclipsing Binary Catalog, third revision:** Kirk et al. 2016, AJ 151, 68,
  DOI 10.3847/0004-6256/151/3/68; VizieR J/AJ/151/68 (table 9, systems with no repeating events,
  `inputs/kirk2016_table9.csv`).
- **Kepler DR25 TCEs:** NASA Exoplanet Archive TAP table `q1_q17_dr25_tce`, rows with
  tce_depth ≥ 400,000 ppm (`inputs/tce_depth_ge_0.4.csv`).
- **Light curves:** Kepler DR25 long cadence, as in "W3 / Kepler dense-cadence pilot" (per-file
  sha256 in `files_<group>.json`).

## W3 / Kepler expansion search (accessed 2026-10-25; `results/w3_kepler/expand/`)
- **Kepler Eclipsing Binary Catalog, third revision:** Kirk et al. 2016, AJ 151, 68; the full KIC list with period
  and morphology from VizieR J/AJ/151/68/catalog (`inputs/kirk2016_catalog.csv`, 2,876 rows).
- **Light curves:** Kepler DR25 long cadence from MAST, as in the pilot. Per-file sha256 in `files_expand.json.gz`.

## W3 / Kepler 2,000-star search (accessed 2026-10-26; `results/w3_kepler/s2000/`)
- Kepler DR25 long-cadence light curves (Q1–Q17), MAST archive
  `https://archive.stsci.edu/missions/kepler/lightcurves/`; every file's name, size and sha256 per star in
  `files_b<b>.json.gz`.
- MAST AWS mirror `https://stpubdata.s3.amazonaws.com/kepler/public/lightcurves/`: one listing and one file fetched to
  test speed (sha256 equal to the archive copy); not used for the search.
- Cousins, R. D. & Highland, V. L. 1992, "Incorporating systematic uncertainties into an upper limit", Nucl. Instrum.
  Meth. A 320, 331–335 (cited as such in Conrad et al. 2003, arXiv:hep-ex/0202013): averaging a zero-event limit over
  the sensitivity's uncertainty, the idea behind `w3_kepler_2000.joint_limit` (which uses a flat rate prior instead).

## Lenstool in the cloud (accessed 2026-10-26; `results/lens_consistency/lenstool_cloud/`)
- **Lenstool 8.15.5**, conda-forge package `lenstool` (build `py314h03d985e_1`): https://anaconda.org/conda-forge/lenstool;
  project page https://projets.lam.fr/projects/lenstool (Jullo et al. 2007, NJPh 9, 447, doi:10.1088/1367-2630/9/12/447).
- **micromamba 2.9.0**: https://micro.mamba.pm/api/micromamba/linux-64/latest.
- **Model:** the pinned CANUCS DR1 MACS0416 Lenstool files (sha256 in `lensmodel.MACS0416_CANUCS`, D-044).
