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

_Open._

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9)

_Open._
