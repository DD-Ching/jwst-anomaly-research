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

## D-002 Archive access and acquisition (unit 1)

_Open._

## D-003 Catalog ingestion and cross-band matching (unit 2)

_Open._

## D-004 Baseline features and anomaly ranking (unit 3)

_Open._

## D-005 Image cutouts and visualization (unit 4)

_Open._

## D-006 External catalog cross-checking (unit 5)

_Open._

## D-007 Candidate store and run provenance (unit 6)

**Decision.**
- Candidate store: one SQLite file through stdlib `sqlite3` (`<outputs>/candidates.sqlite`, gitignored)
  with tables `runs`, `candidates`, `status_history`, `vetting_notes`. Variable-shaped fields (per-method
  scores, per-band image flags, cross-match summary, cutout paths, run context) are JSON text columns,
  queryable with SQLite's built-in JSON functions. Schema version = `PRAGMA user_version`; a file from
  newer code is refused; migrations go in `CandidateStore._init_schema`.
- Status lifecycle `new → triaged → {artifact, known_object, explained, unexplained}`, `unexplained →
  followup`, re-open via `triaged`. Outcome statuses and `followup` require ≥1 vetting note (test,
  outcome `pass|fail|inconclusive`, evidence, provenance label, author, UTC), so every conclusion
  names its tests (docs/methodology.md). Candidates are keyed by `(run_id, sample_id, source_uid)`;
  statuses do not carry across runs, and a uid present in several runs or samples must be
  addressed with `--run`/`--sample` (uids come from per-run catalogs; the store never guesses).
- Run provenance with the standard library only: `git` via `subprocess` (commit, branch, dirty
  including untracked files, with `data/manifests/` reported separately; None unless the package
  sits at a checkout root), sha256 of the config bytes plus a verbatim copy in the run directory,
  `importlib.metadata` versions of key packages, Python/platform, UTC time. Run id
  `YYYYMMDDTHHMMSSZ-<8 hex>` sorts chronologically. Interrupted runs are recorded as `failed`.
- Runner: plain in-process Python that calls stage functions as module attributes and checks each
  output with `schema.validate`. A failing required stage aborts the run (recorded `failed`, report
  written); failing optional stages (cutouts, crossmatch) are recorded and the run continues. Only
  `CAT` products are fetched; i2d images are read via S3 by the cutout stage. Each sample is ranked
  separately. Intermediate tables are ECSV (Parquet optional) under `<outputs>/runs/<run_id>/`, plus
  `report.md`. New optional config block `outputs: {candidates_top_k, table_format}`.
- CLI: stdlib `argparse` (`jwst-anomaly run | runs list | candidates list|show|set-status|add-vetting`).

**Alternatives rejected.**
- DuckDB: a columnar engine for analytical scans and a new dependency. The store is small,
  transactional and updated row by row (status, notes). DuckDB can attach this SQLite file if
  analytics are needed later.
- SQLAlchemy + Alembic: too heavy for four tables; `user_version` plus in-code migrations suffices.
- Datasette / sqlite-utils as dependencies: good viewers of the same file, usable ad hoc, not needed by the code.
- GitPython: a handful of `git` subprocess calls do not justify a dependency. Click/Typer: argparse covers the CLI.
- Workflow engines (Snakemake, Prefect, Luigi): the run is one linear chain of in-process calls.
- Experiment trackers (MLflow, DVC): surveyed by unit 9 (D-010); the run context is a plain JSON dict
  that can be logged to a tracker later.

**Evidence.** SQLite positions itself as a local application file format that "does not compete with
client/server databases" (sqlite.org/whentouse.html); `user_version` pragma (sqlite.org/pragma.html);
JSON functions built in by default since SQLite 3.38.0 (sqlite.org/json1.html); CPython 3.12.10 ships
SQLite 3.49.1 (measured 2026-10-07). DuckDB attaches SQLite files (duckdb.org SQLite extension). Offline
tests (`tests/test_{candidates,provenance,pipeline,cli}.py`) run the full chain with contract-valid fake
stages; links in SOURCES.md "Candidate store and run provenance (unit 6)".

**Revisit if.** Several processes write the store concurrently (parallel cloud runs): use per-run files
merged later or a server DB. Candidates reach millions or analytics dominate: DuckDB/Parquet. Source
uids become stable across runs: carry vetting status across runs. Stages become long-running or fan
out over many programs: Snakemake.

## D-008 Open-source project tooling (unit 7)

_Open._

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9)

_Open._
