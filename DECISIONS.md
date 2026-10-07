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

## D-006 External catalog cross-checking (unit 5) (2026-10-07)

**Decision.**
- `crossmatch.py` composes astroquery 0.4.11 and pyvo 1.9.1. pyvo is already an astroquery
  requirement and is now declared directly because we import it. Our own code only normalizes
  columns and associates targets with returned objects via `astropy.coordinates.search_around_sky`.
  Each service has a backend chain:

  | service | primary (batched) | fallback |
  |---|---|---|
  | simbad | CDS XMatch `simbad`, one upload for all targets | SIMBAD TAP, `Simbad.query_region` (OR'd cones ≤300, upload above) |
  | gaia | CDS XMatch `vizier:I/355/gaiadr3`, one upload | ESA archive TAP (pyvo) on `gaiadr3.gaia_source`, 100 cones/request |
  | ned | NED TAP `NEDTAP.objdir` via pyvo, 100 OR'd cones/request | `Ned.query_region`, one request per target (≤100 targets) |

- No N-dependent crossover. The batched paths cost about one request per upload or per 100 cones
  (Evidence), so one backend serves 1 to 10⁴ targets and results don't depend on sample size. A
  per-object loop exists only as the NED fallback.
- Robustness.
  - Every HTTP request gets a 60 s client timeout and up to 2 retries with backoff (2 s, 4 s).
    Retries are per request, so one failed NED chunk doesn't repeat the others.
  - Deterministic failures (row-limit truncation, the fallback's target cap, bad input) are not
    retried.
  - Then the fallback backend. A service that still fails goes to
    `meta["services_failed"]`/`meta["service_errors"]` and the others continue;
    `CrossmatchError` is raised only if every service fails.
- Cache. Normalized results are written as ECSV to `data/cache/crossmatch/`, keyed by service,
  backend, radius and positions (30-day expiry; `refresh=True` forces a query). Each backend's entry
  is checked in chain order, so a cached fallback result is used only when the primary query
  fails. astroquery's own HTTP cache is bypassed, so `meta["query_utc"]` is the real query time.
- Output. `schema.XMATCH_COLUMNS` plus `n_simbad`/`n_ned`/`n_gaia` (−1 = not queried or failed),
  `is_lens_related` and `lens_types`. `query_matches()` gives the long format. Provenance `observed`.
- ASSUMPTIONS. `is_star` = the nearest SIMBAD/Gaia match is a star, so a galaxy with a Gaia star
  0.9″ away is not flagged. A star is a SIMBAD main otype in the non-candidate `*` branch minus
  PN/SN\*/Pl/out, or a Gaia DR3 source with parallax/σ ≥ 5 or proper-motion significance ≥ 5
  (χ² with the pmra–pmdec correlation). NED's `*` means "Star or Point Source" and is not used.
- Lens-related: any SIMBAD otype (main or other) in the `grv > gLS` branch, candidates included
  (`gLS LS? gLe Le? LeI LI? LeG LeQ`), plus `Lev`; or NED `G_Lens`/`Q_Lens`. Both otype sets come
  from SIMBAD `otypedef.path` via `derive_simbad_otype_sets` and are frozen in code; a network test
  checks for drift.
- Best match = smallest separation, ties broken by `services` order. Nearest is not the same as
  counterpart for extended sources.

**Alternatives rejected.**
- `Ned.query_region` as NED primary: no multi-position query; ≈0.9 s per target (Evidence).
- XMatch only above a size threshold, per-service cones below it: no speed gain at small N, and an
  N-dependent backend makes results depend on sample size.
- SIMBAD TAP as SIMBAD primary: equally fast. XMatch keeps SIMBAD and Gaia on one client and
  returns `other_types` directly. TAP stays as the fallback and gave identical results in the e2e.
- ESA Gaia archive as Gaia primary: 2× slower here, and astroquery warns that the archive "may be
  unstable" ahead of DR4.
- `astroquery.gaia` (TapPlus) for the Gaia fallback: its HTTP connections take no client timeout,
  and it adds `TOP 2000` to synchronous jobs. pyvo reaches the same ESA TAP service with both
  handled.
- One field-wide NED cone plus local association: 957 objects for a 3.6′ cone on SMACS 0723,
  versus only the matches with OR'd cones.
- An HTTP cache library (e.g. `requests-cache`): a new dependency for what ~40 lines of ECSV caching
  does. Our cache also keeps the normalized table and its query time as a reviewable record.

**Evidence.** Timings on 2026-10-07: radius 1″, random F200W-catalog positions in SMACS 0723.

| N | SIMBAD TAP | XMatch simbad | Gaia archive | XMatch Gaia | NED TAP | NED per-target |
|---|---|---|---|---|---|---|
| 5 | 0.58 s | 1.05 s | 1.68 s | 0.74 s | 1.34 s | 4.48 s |
| 50 | 0.73 s | 0.68 s | 2.04 s | 0.82 s | 1.80 s | — |
| 500 | 1.21 s | 0.75 s | — | 0.81 s | 10.87 s (5 chunks) | — |

SIMBAD TAP and XMatch returned the same row counts at N = 50 and 500.

E2E run: 55 real targets (5 in `tests/data/xmatch_targets.ecsv` plus 50 F200W-catalog sources).
- Primary chain: 4.5–4.7 s wall, no failures. Fallback chain: 41–49 s (NED per-target 31–39 s).
  A cached rerun took 0.04 s.
- One live test run hit a transient NED server error ("can't start the database search program").
  NED was recorded as failed and the other services completed.
- Primary and fallback agree on `is_known_object`, `is_star`, `is_lens_related`, `n_simbad` and
  `n_gaia` for 55/55 targets, and on every SIMBAD and Gaia match id.
- NED gave 12 matches either way, but only 8 preferred names agree. The TAP `objdir` table describes
  itself as release "N36.1" and can differ from NED's interactive search in names, positions and
  types (e.g. `*` vs `IrS` for one Gaia star).

**Revisit if.**
- NED TAP adds uploads or drifts further from NED's search.
- Samples reach ~10⁴ targets: NED TAP chunks run one after another at ~2 s each; run them
  concurrently or use one field-wide cone.
- Vetting needs epoch-propagated Gaia positions (J2016 vs JWST 2022 offsets up to 0.74″ were seen
  in this field).
- Gaia DR4 replaces `I/355`.
- Extended sources need a size-dependent radius.

## D-007 Candidate store and run provenance (unit 6)

_Open._

## D-008 Open-source project tooling (unit 7)

_Open._

## D-009 Agent operations: skills, permissions, scheduling (unit 8)

_Open._

## D-010 Tools for future milestones (unit 9)

_Open._
