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

## D-002 Archive access and acquisition (unit 1) (2026-10-07)

**Decision.**
- Query and download through `astroquery.mast.Observations` (0.4.11), wrapped by `query.py` and
  `acquire.py`: `query_criteria` (`obs_collection="JWST"`, `dataRights="PUBLIC"` server-side and
  re-checked client-side), one batched `get_product_list` per call, `filter_products(calib_level=[3],
  productSubGroupDescription=[...], dataRights=["PUBLIC"])`, and `download_file` per product. If
  `$MAST_API_TOKEN` is set, `Observations.login(token=...)` runs once; the token is never logged or stored.
- A level-3 product list also contains every level-2 member. For NIRCam F200W that is 1,158 rows,
  including 72 per-detector `_i2d.fits` of about 119 MB at `calib_level` 2. Members carry their own
  `obsID`, with `parent_obsid` pointing at the level-3 observation. `list_products` keeps a row only if
  its `obsID` is one of the requested observations, its file name starts with that `obs_id + "_"`
  (the `_` keeps F150W from matching F150W2), and its `calib_level` is 3.
- Downloads go to a `*.part` file next to the target. The size is checked against MAST, the sha256
  is computed, and `os.replace` moves the file into place, so a file at its final path is always
  complete. On a re-run, local files are hashed and skipped when the sha256 matches the manifest. A
  file that is present with the MAST size but has no manifest row is "adopted": it is hashed and
  recorded, not downloaded again. By default no single file over 200 MB is downloaded (CLAUDE.md).
- Manifest format: ECSV with `schema.MANIFEST_COLUMNS`, sorted by `dataURI`, LF line endings. It is
  written atomically and only when its content changes, so a re-run leaves tracked files
  byte-identical. `local_path` is POSIX and relative to `paths.data_root()`:
  `cache/mast/<obs_id>/<file>`.
- `pipeline_version` is read from each file. For `_cat.ecsv` it comes from the ECSV meta `version`
  (jwst, photutils, astropy); for FITS, from the primary header only (`CAL_VER`, `CRDS_CTX`). A
  warning is logged if it disagrees with MAST's `prvversion`. Versions differ per program: 2736
  uses jwst 2.0.1 and 1345 uses jwst 3.0.0.
- Cloud: `list_products(cloud_uris=True)` adds `cloud_uri` (`s3://stpubdata/...`) from MAST's
  `path_lookup` service via `astroquery.mast.utils.mast_relative_path` (up to 50 URIs per request, no
  boto3). The fetch script writes these for every listed level-3 product, including the NIRCam i2d
  files it skips, to `data/manifests/<config>_products.ecsv` for byte-range reads.

**Alternatives rejected.**
- `Observations.download_products`: it writes in place (no temp file), its cache check only compares
  size via HEAD, and it keeps no checksums. We use its per-file primitive `download_file` instead.
- `enable_cloud_dataset()` + `get_cloud_uris`: these need boto3, which is not a dependency. boto3
  1.43.108 needs `botocore>=1.43.108`, while aiobotocore 3.9.2 (pulled in by s3fs in the `cloud`
  extra) needs `botocore<1.43.107`, so the two would have to be pinned together. It also sends one S3
  HEAD per product. `path_lookup` is the same service without either cost.
- Downloading from S3 instead of MAST HTTPS: the bucket holds only public data. MAST HTTPS also
  accepts tokens and supports Range, and it is astroquery's default.
- Raw `requests` against the MAST API: astroquery already handles the service calls, auth and
  product filtering.
- pooch: it is built around hashes known in advance (MAST publishes none), its registry has no size,
  pipeline-version or provenance fields, and it would be a new dependency.
- DVC or DataLad: too heavy for about 20 files that MAST already serves at stable URIs.
- Other manifest formats: CSV has no types or meta, JSON gives noisy diffs, and Parquet cannot be
  diffed.

**Evidence.**
- Probes on 2026-10-07 with astroquery 0.4.11. Program 2736 has 12 PUBLIC level-3 image
  observations. Their combined product list has 5,176 rows; 72 are `calib_level` 3 and 24 are CAT or
  I2D.
- `path_lookup` keys were checked on S3 with s3fs (sizes match MAST) for the 2736 F1500W catalog,
  the F200W i2d and the 1345 F200W catalog.
- MAST's download endpoint sends `Accept-Ranges: bytes`.
- End-to-end run: `--catalogs-only` fetched 11 catalogs (20.7 MB). A second run downloaded nothing
  and left both manifests unchanged, and `sha256sum -c` agrees with the manifest.
- The HTTP 404 response for EXCLUSIVE_ACCESS products comes from the coordinator's MAST check and
  was not re-tested here.

**Revisit if.**
- astroquery drops `mast_relative_path`: switch to `get_cloud_uris` and boto3 once the botocore
  pins line up.
- MAST publishes checksums: verify against them, not just size.
- Hashing every file on every run gets slow: trust size plus mtime.
- We freeze a dataset to distribute: export the manifest as a pooch registry.
- Bulk downloads justify in-region S3.

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

**Decision.** Use official templates and GitHub-native or maintained tools for community and repository
hygiene, and write only the project-specific content (2026-10-07).
- `CITATION.cff` uses Citation File Format 1.2.0, which GitHub shows as "Cite this repository", and is
  validated with `cffconvert`. The author is recorded as the GitHub alias `DD-Ching` only. No real
  name, ORCID or DOI until the owner supplies them.
- `CODE_OF_CONDUCT.md` is Contributor Covenant 2.1, copied verbatim from contributor-covenant.org.
  Only the contact placeholder is filled in, with GitHub private reporting and GitHub "Report abuse"
  instead of an email address.
- `SECURITY.md` uses GitHub private vulnerability reporting as the only channel. It states that the
  project handles no user data and has a secrets policy (`MAST_API_TOKEN` comes from the environment only).
- Issue forms (GitHub YAML form schema): task, bug and anomaly candidate. Blank issues are disabled.
  The candidate form requires an outcome (ruled out / not ruled out / not tested) for each of steps
  1-3 in docs/methodology.md, plus provenance labels and integrity acknowledgements. The PR template
  mirrors the worker PR body plus a checklist.
- pre-commit uses `astral-sh/ruff-pre-commit`, restricted to CI's paths `src|tests|scripts`, and
  `pre-commit/pre-commit-hooks`: files ≤1 MB with `--enforce-all`, private keys, YAML/TOML syntax,
  whitespace and merge markers. Both are pinned to release tags. CI's unpinned ruff is authoritative.
- Dependabot runs weekly with one grouped PR per ecosystem and a 7-day cooldown, for `github-actions`,
  `pre-commit` and `pip`. For pip it uses `versioning-strategy: increase-if-necessary`. pyproject has
  `>=` floors only and no lock file, so this entry normally opens no PRs. It is kept so the ecosystem
  is configured; the real switch is to `uv` once a lock file exists.
- `network-tests.yml` runs `pytest --run-network` weekly and on manual dispatch (ubuntu, py3.12,
  read-only token). A failed scheduled run opens or updates the issue "Weekly network tests failing"
  from a separate job that has only `issues: write`, so agent cycles see it in `gh issue list`.
  PR CI stays offline.

**Alternatives rejected.**
- A hand-written or paraphrased code of conduct: it would drift from a maintained standard without
  adding value.
- Contributor Covenant 3.0 (published at contributor-covenant.org): not adopted because the bootstrap
  plan fixed 2.1. Moving to 3.0 is an owner choice.
- An email contact for reports: the repository is public and no address was provided. GitHub private
  reporting needs none.
- Markdown issue templates: forms enforce required fields and fixed options, which agents and humans
  fill in consistently. A pre-filled table in a required textarea would always count as filled in.
- Renovate: more configurable, but it needs a GitHub App install. Dependabot is built in and supports
  all three ecosystems.
- Dependabot's default pip strategy (`auto`). It applies `increase`, which raises lower bounds on
  every new release, to projects it classifies as apps. That has no reproducibility gain. A lock file
  is the right tool for reproducibility (see Revisit if).
- pre-commit.ci: it is a third-party app and pushes autofix commits to PR branches. Running
  `pre-commit run --all-files` as a step in ci.yml is instead recommended to the coordinator, because
  ci.yml is outside this unit. That step would enforce the 1 MB and private-key rules on the server and
  give CI and local runs one ruff version.
- Network tests on every PR: CI would depend on archive uptime and become flaky.

**Evidence.** `cffconvert --validate`: valid against CFF 1.2.0. `pre-commit run --all-files`: all
hooks pass and no files are modified. actionlint 1.7.12 is clean on both workflows. check-jsonschema
0.38.2 with its vendored SchemaStore schemas validates the issue forms, issue config, Dependabot config,
workflows and CFF. The `pre-commit` Dependabot ecosystem, `cooldown` and versioning strategies come
from the Dependabot options reference. Hook tags are the latest releases per the GitHub API on
2026-10-07. The Contributor Covenant text matches the official markdown line by line except for the
contact line. The official file's leading and trailing blank lines were trimmed, so the sha256 in
SOURCES.md is of the downloaded official file, not of CODE_OF_CONDUCT.md. Links are in SOURCES.md
"Open-source project tooling (unit 7)".

**Revisit if.** A `uv.lock` is committed: replace the pip entry with the `uv` ecosystem. The first
release is tagged: add `date-released` and a Zenodo DOI to CITATION.cff and keep its `version` equal
to pyproject's. More maintainers join, or the owner provides a contact: give conduct reports a
dedicated channel separate from the security form. Oversized files, keys or lint drift reach PRs: add
the pre-commit step to CI. The weekly network job fails repeatedly for upstream reasons: add retries
or narrow its scope. The repository has no activity for 60 days, which makes GitHub disable the
scheduled workflow: re-enable it from the Actions tab.

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9)

_Open._
