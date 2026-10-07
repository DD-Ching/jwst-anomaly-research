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

**Decision.**
- Read `_cat.ecsv` with astropy's ECSV reader, which restores `sky_centroid`/`sky_bbox_*` as `SkyCoord`,
  the units and the pipeline meta. `load_pipeline_catalog` keeps every pipeline column and adds ICRS
  `ra`/`dec`. It records in `meta`: `band` (parsed from the MAST file name, since the ECSV has no
  filter keyword), `obs_id`, `input_sha256`, `jwst_version`, `photutils_version`, and a derived
  `pixel_scale_arcsec` (affine fit of sky offsets vs. centroids). Provenance `observed`.
- `merge_bands` uses astropy `search_around_sky`. Matching is one-to-one per band and greedy by
  increasing separation (ties broken by index); default radius 0.1″ (config). The reference band is
  matched first, then the other bands by wavelength. A row's position is its anchoring detection
  (per-row `anchor_band`); it is not averaged across bands. `ref_band` is the merge's reference band on
  every row. Provenance is `derived`, or `simulated` if any input is simulated (injection-recovery
  stays identifiable). Each input's meta, provenance included, is kept in `meta["inputs"]`.
- **Union by default** (`include_unmatched=True`). A detection without a counterpart starts a new row
  anchored on its own band, so dropouts and LW-only sources survive.
  `include_unmatched=False` returns reference-band rows only. Filtering downstream is cheap and
  reversible; dropping rows at ingestion is not. Flags: `<band>_detected`, `<band>_sep_arcsec`,
  `n_bands`, and `ref_nn_sep_arcsec` (distance to the nearest reference detection not in the row),
  which marks likely split or blended detections.
- `source_uid = <field_id>_<anchor band>_<anchor label>`, e.g. `jw02736-o001_t001_nircam_f200w_42`.
  `field_id` defaults to the reference `obs_id` minus its optics part and is restricted to
  `[A-Za-z0-9._-]` so it is safe in file names. The uid is stable for a fixed set of input files
  (sha256 in `meta["inputs"]`). It changes when MAST reprocesses the catalogs, because labels change.
- Non-detections are filled with NaN, -1 or False instead of being masked, so `np.asarray` and pandas
  never see hidden values; `<band>_detected` tells the cases apart. Values already masked in the input
  stay masked (floats become NaN).
- Inputs are validated: bands must be alphanumeric and agree with `meta["band"]`; `ra`/`dec` are
  converted to degrees from their units. Rows with non-finite positions are dropped with a warning
  and counted as `n_no_position`.
- **No bulk astrometric offset correction.** `meta["match_stats"]` reports, for matches to
  reference-band rows, the median separation and median offset. Warnings fire, given ≥10 possible
  pairs, when the offset exceeds 25% of the radius or fewer than 20% of possible matches are found
  (the case where the offset exceeds the radius). On program 2736 every offset is ≤13 mas, and removing
  the largest (F090W) adds 2 matches at 0.1″ (1615→1617). A warning is also raised when bands come from
  different jwst/photutils versions.

**Alternatives rejected.**
- Reference-only merge as the default: on 2736 it silently drops 2,109 of 5,254 rows. Among them are
  129 sources detected in F444W with no SW detection and no F200W source within 1″, which is exactly the
  dropout and very-red population.
- Mutual nearest neighbour: misses some pairs in crowded areas that greedy-by-separation still pairs
  uniquely. Both are deterministic.
- STILTS `tmatchn` (Java tool): equivalent for this task and adds a JVM; astropy is already a dependency.
- NWAY (Bayesian N-catalogue matching with positional errors and priors): built for catalogs with
  very different positional errors. Overkill at 8–20 mas median separations within one program.
- A per-row `ref_band` that varies with the anchor: one column name with two meanings. The anchor gets
  its own `anchor_band` column.
- CDS XMatch: a remote service for external catalogs (D-006), not for merging local bands.
- Multi-band averaged positions: they mix PSF-dependent centroids. The SW reference band is the
  sharpest single anchor.
- Masked columns for non-detections: `np.asarray(masked_column)` exposes arbitrary fill values downstream.
- Forced photometry (photutils on `_i2d`, or HLSP catalogs): the proper fix for independent detection,
  but out of scope for v0 (see D-001 "Revisit if").

**Evidence.** Program 2736 NIRCam level-3 catalogs (jwst 2.0.1, photutils 2.3.0), reference F200W,
radius 0.1″, astropy 8.0.1, run 2026-10-07:

Median separation and offsets are computed over the matches to F200W only:

| Band | Detections | Matched | Matched to F200W | Unmatched | Median sep | Median ΔRA·cosδ / ΔDec |
|---|---|---|---|---|---|---|
| F090W | 2652 | 1615 | 1615 | 1037 | 20.3 mas | −12.9 / −3.0 mas |
| F150W | 2962 | 2491 | 2381 | 471 | 8.1 mas | −1.0 / −0.1 mas |
| F277W | 1860 | 1636 | 1610 | 224 | 10.1 mas | −2.3 / +1.4 mas |
| F356W | 1983 | 1779 | 1614 | 204 | 11.3 mas | −0.4 / +0.3 mas |
| F444W | 1877 | 1704 | 1501 | 173 | 11.5 mas | +2.1 / +1.1 mas |

- The union has 5,254 rows (3,145 anchored on F200W). Rows with `n_bands` = 1…6:
  2362 / 679 / 540 / 256 / 387 / 1030.
- Single-band rows are mostly faint or fragments. F090W-only rows have median aper50 S/N 3.8, 61% are
  below 5, and their median `nn_dist` is half the catalog's.
- Pixel scales: NIRCam SW 0.03123″/pix, LW 0.06291″/pix.
- CEERS `jw01345-o001_t021_nircam_clear-f200w` (jwst 3.0.0, photutils 3.0.0) has the same columns but
  renames two meta keys: `apermask_method`→`aperture_mask_method`, `localbkg_width`→`local_bkg_width`.
- Tools and data: astropy matching (https://docs.astropy.org/en/stable/coordinates/matchsep.html) and
  JWST `source_catalog` (https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/main.html).
  Data sources are in SOURCES.md "Catalogs and cross-band matching".

**Limitations.**
- Every band is detected and deblended independently. A missing band means "no matching detection",
  not a flux limit.
- Colors are approximate because apertures are encircled-energy based per band and deblending differs.
- Pixel-unit columns (`nn_dist`, `semimajor_sigma`, `isophotal_area`) are in each band's own pixels;
  convert them with `pixel_scale_arcsec` before comparing bands.
- `aper_total_*` is valid for unresolved sources only, according to the pipeline docs.

**Revisit if.**
- Colors or dropout claims need to be quantitative: use forced photometry (D-001).
- A field triggers the offset warning: add a bulk offset correction.
- `n_contested` grows to a noticeable fraction of matches in denser or deeper fields.
- Multiple catalog versions per band need to coexist (uids would then need a version component).

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
