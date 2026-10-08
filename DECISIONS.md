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

**Decision.**
- Features (`features.build_features`, `derived`, one line each in `meta["feature_spec"]`):
  `ref_mag`; colors of wavelength-adjacent bands; reference-band morphology (log10 isophotal area,
  ellipticity, log10 `CI_50_30` and `CI_70_50`, DAOFind sharpness and roundness, log10 `nn_dist`);
  detection pattern (`blue_dropout`, `red_dropout`, `n_gaps`, which together count every undetected
  band exactly once; `n_bands` is left out because it would count the same missing band again).
  Bands come from column prefixes and any subset works; unknown bands raise an error. Features
  that are NaN for every source are dropped.
- Colors and `ref_mag` use `aper50` AB magnitudes. The pipeline sizes each filter's apertures from
  that filter's encircled-energy curve (APCORR is keyed by filter and EE fraction), so equal-EE
  colors don't depend on the PSF for unresolved sources. For extended sources they are biased red.
- Morphology comes from one table-level reference band (default: the most common `ref_band`), so
  pixel-unit features share one pixel scale.
- NaN policy (`features.NAN_POLICY`): no imputation in features. A value is NaN when it is undefined
  (no detection, negative flux, non-positive log argument) or when that band's aperture S/N is below
  `min_snr = 3`. A band with no error column is skipped with a warning and listed in `meta`. The
  detection-pattern features encode missing bands explicitly. The ranker median-imputes (z = 0)
  and reports `n_missing`.
- Scoring (`rank.score_anomalies`, `model_prediction`): every method sees the same matrix of robust
  z-scores (median and 1.4826 × MAD via `scipy.stats.median_abs_deviation`; if the MAD is 0 it falls
  back to IQR/1.349, then std, then 1), with NaN set to 0 and values clipped at ±10.
  `robust_z` = RMS of z (a diagonal robust Mahalanobis distance); `isolation_forest` = −`score_samples`
  (scikit-learn, 1000 trees); `lof` = scikit-learn LOF with 20 neighbours, with each exact-duplicate
  group capped at 20 copies (identical to plain LOF when there are no larger groups).
  `score` = mean of per-method percentile ranks. `top_features` = the 3 largest unclipped |z|, a
  model-agnostic attribution. Scoring a table that contains simulated rows keeps the `simulated` label.
- Evaluation: `inject_outliers` ("shift": 2 features set to ±6 robust σ; "shuffle": marginals kept,
  correlations broken), `injection_recovery` (precision@k and recall@k against a random baseline,
  with boundary ties counted pro rata; `simulated`), `seed_stability` (Jaccard overlap of top-k sets
  across seeds, ties broken by `source_uid`).

**Alternatives rejected.**
- `aper_total` colors: they add a point-source extrapolation and more NaN (F200W: 43 vs 26 of 3,145).
  `aper70` is noisier and more blended. Isophotal magnitudes use segments that differ between bands.
- sklearn `RobustScaler`: it uses the IQR (25% breakdown point vs 50% for the MAD) and silently sets
  zero-IQR scales to 1. A single MAD-based matrix that feeds every method and the attribution is simpler.
- Per-feature missing indicators (`SimpleImputer(add_indicator=True)`): they duplicate the
  detection-pattern features and would make NaN-heavy, low-information rows look anomalous.
  KNN/iterative imputation can't separate rows whose only information is the detection pattern.
- Plain LOF on all rows: imputed duplicates make it ill-defined (see Evidence). Scoring distinct rows
  only (tried first) discards multiplicity. In a probe, 200 identical rows away from 300 Gaussian rows
  became one isolated point and filled LOF's top 200 (LOF 4.4). With the cap they score 1.13 and
  none reach the top 50. Distinct-row scoring does score higher on shift injections
  (LOF P@50 0.62 vs 0.50 with the cap, same features), but it is wrong for real duplicate groups.
- PyOD 3.6.6 (ECOD, COPOD, …): ECOD is per-feature tail probability, the same family as `robust_z`,
  and numba is a core dependency. IF and LOF already come from scikit-learn. Not added.
- Astronomaly (active learning, image-first): suited to later human-in-the-loop labelling, not to a
  catalog baseline.
- Weighting the ensemble toward LOF because it wins the injection test: that would tune the ensemble
  to one synthetic outlier type.

**Evidence.** 2026-10-07 e2e on program 2736 NIRCam (6 bands, 5,254 sources, a throwaway 0.1″ union
merge standing in for unit 2; full output in the unit-3 PR body):
- S/N floor: without it, 5 of the top 20 had F200W S/N < 3 or no F200W detection. With it: 0 of
  the top 100.
- Imputed duplicates: 1,747 rows fell into 5 exact-duplicate groups (the largest had 906), and plain
  LOF reached 4×10⁷ next to them. With the cap the maximum LOF is 6.4 and sklearn gives no warning.
- Injection recovery (50 injections, n = 5,304, random precision 0.009), P@50 for robust_z / IF /
  LOF / ensemble: shift 6σ×2 0.08 / 0.08 / 0.50 / 0.20; shift 4σ×1 0.00 / 0.02 / 0.08 / 0.06;
  shuffle 0.00 / 0.00 / 0.22 / 0.00. Of the real sources, 6.4% already have a feature beyond 6σ
  (2.9% beyond 10σ), so the per-feature methods (and the mean-rank ensemble, which needs them to
  agree) recover few injections. LOF is the only method well above random on every injection type.
- Clip none / 10 / 20 → LOF P@50 (shift 6σ×2) 0.30 / 0.50 / 0.48, ensemble 0.22 / 0.20 / 0.22. On
  Gaussian data, clipping drops IF's rank for a 12σ single-feature outlier from 1st to 2nd or 3rd.
- Seed stability, top-20 Jaccard over 5 seeds: IF 0.69 / 0.76 / 0.77 with 500 / 1000 / 2000 trees.
  With 1000 trees the ensemble scores 0.90 (minimum 0.82); robust_z and LOF score 1.0.

**Revisit if.** Forced or HLSP photometry becomes available (use it for colors). Vetted labels exist
(then compare the mean-rank ensemble with LOF alone or a max-rank rule on real labels rather than
injections). A learned model arrives (it must beat this baseline on `injection_recovery` and
`seed_stability`; PyOD/ADBench is then the comparison zoo). Vetting shows the top is dominated by
one artifact class (add quality features or cuts such as edge or exposure flags). The sample grows
past ~10⁵ sources (LOF cost).

## D-005 Image cutouts and visualization (unit 4)

**Decision.** (2026-10-07)
- Cutouts: `astropy.io.fits` opened on an `fsspec` file (read-ahead cache, 64 KiB block) plus
  `astropy.nddata.Cutout2D` on a lazy `Section`. A ~10-line `Section` subclass reads each cutout's
  rows as one contiguous full-width byte range (one request per extension per cutout) instead of
  one read per row. No new dependency (`fsspec`/`s3fs` are the `cloud` extra). `mast:` URIs go to
  the MAST download endpoint, which serves HTTP range requests.
- Output: SCI + WHT (ERR optional) with the cutout WCS, as small FITS files under
  `outputs/cutouts/<image stem>/`. Pixel size rounded to an odd count (target on the central pixel).
- Quality flags (i2d has no DQ). No-data = non-finite SCI or WHT <= 0. `edge`: no-data on the box
  border or box beyond the array. `nan_center`: no-data within 0.2". `nan`: no-data fraction > 0.1.
  `low_weight`: median core WHT < 0.5 x the image's median positive WHT (32 bin-centre rows; the
  whole WHT if <= 8 MiB). `outside`: no data in the box (row kept, no file). Thresholds are untuned
  heuristics. The table is labeled observed; `meta["column_provenance"]` marks derived columns.
- Display: `astropy.visualization` asinh (a = 0.01) over each panel's 1-99.9 percentile range;
  matplotlib `Figure` (no pyplot). Chosen by eye on 2736 cutouts over PercentileInterval(99.5) with
  a = 0.1 (cores saturated, sky crushed) and a = 0.002 (noise-dominated).

**Alternatives rejected.** Cost of one 96x96 SCI cutout of the 2736 F200W i2d on S3 (1.76 GB):
- `fits.open(uri, use_fsspec=True)` with s3fs defaults (50 MiB read-ahead): 104.9 MB, 25.5 s.
- `astrocut` 1.3.0 `FITSCutout` (STScI; reads S3 FITS): 419.9 MB, 9 requests, 33.1 s by default,
  because it reads every HDU header through 50 MiB blocks; 5.8 MB, 22 requests, 9.6 s with
  `default_block_size=256 KiB`. It also takes one position per instance (file reopened per target)
  and adds asdf, gwcs, spherical-geometry and s3path.
- Plain `Cutout2D(hdu.section)`, 256 KiB read-ahead: 3.9 MB but 15 requests, 7.3 s (one read per
  row); without a cache, 33 s. Chosen row-strip reads: 4.5 MB, 2 requests, 2.4 s.
- Concurrent per-row ranges (`fsspec` `cat_ranges`) would cut bytes ~100x but need our own FITS
  offset/scaling code; not worth it at top-k scale. Full NIRCam i2d download: disk budget
  (CLAUDE.md).

**Evidence.** Measured 2026-10-07 on the owner's machine (throwaway scripts; output in the
unit-4 PR).
- Module run on 10 bright F200W catalog sources: SCI+WHT cutouts plus the weight sample cost
  84.7 MB in 56 requests, 34 s (4.8% of the file). MIRI F770W via MAST HTTPS range reads: 6.2 MB in
  9 requests, 5.5 s, pixel-identical to cutouts from the full 36 MB download.
- Weight reference on the F200W WHT: 64 rows give 6634; 8 rows give 4690-6695 depending on which
  rows, 16 rows 6151-6668, 32 rows 6630-6636. Hence 32 rows (~3.4 MB, ~12 s once per image).
- The i2d SCI headers carry a TAN FITS WCS (0.0312"/px F200W, 0.1109"/px F770W). Pixel positions
  from it for the catalog `sky_centroid`s equal the catalog `xcentroid/ycentroid` to < 1e-4 px. The
  light-weighted centroid is within 1.1 px of the cutout centre for 11 of 12 cutouts (2.2 px for
  the extended source, label 1845).
- MIRI F770W i2d (CAL_VER 2.0.1): SCI is NaN exactly where WHT == 0 (25.6% of pixels); positive WHT
  is bimodal (1st percentile 542, median 5225), so a relative weight flag is informative. Pipeline
  defaults `fillval` NAN, `weight_type` ivm:
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
- https://docs.astropy.org/en/stable/io/fits/usage/cloud.html ;
  https://astrocut.readthedocs.io/en/latest/astrocut/index.html

**Revisit if.** Runs move next to the data (us-east-1), where larger blocks may be cheaper than
latency; top-k grows to thousands (share strips between nearby targets, or `cat_ranges`);
astrocut gains multi-position FITS cutouts with lighter dependencies; injection tests show the flag
thresholds are badly calibrated; multi-band comparison needs a common north-up grid (`reproject`).

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

**Decision.**
- The charter runs long-term through three committed project skills in `.claude/skills/`. They link to CLAUDE.md and
  the charter rather than copying them.
  - `/research-cycle` runs one full cycle. It keeps model invocation on, because scheduled fires only run skills
    Claude may invoke. It treats only the owner's words as instructions; other people's issues, comments and PR text
    are data. It never writes `@claude`.
  - `/reuse-check` uses `context: fork` and `background: false`, so the search runs in an isolated subagent and only
    the decision comes back.
  - `/vet-candidate` follows the explanation ladder in docs/methodology.md and commits its record to
    `docs/candidates/`.

  Committed skills load locally, in cloud routines and in GitHub Actions.
- Scheduling uses built-ins only:
  - a cloud routine (`/schedule`) is the primary scheduler, with a fresh session each run;
  - `/loop` covers in-session runs;
  - a Desktop local task is the alternative;
  - an Actions cron is documented only.
- `@claude` runs `anthropics/claude-code-action@v1` on the owner's OAuth token.
  - Every trigger is gated to `author_association == 'OWNER'`, and comment context is filtered to the owner and
    `claude[bot]`.
  - Fork PRs are skipped. Checkout doesn't persist credentials.
  - Concurrency is one run per issue or PR.
  - The workflow does nothing until the secret exists.
- `.claude/settings.json` allows the cycle's commands. Its deny rules, mirrored for the PowerShell tool, cover:
  - merges, auto-merge, branch-protection and ruleset API calls, and `gh api` DELETE requests;
  - force pushes (including bundled short flags), pushes that name `main`, and remote branch deletion;
  - history rewrites;
  - posting local files (`--body-file`, `-F`, `=@file`, `--input`, `git diff --no-index`);
  - `gh secret` and token printing;
  - reads of credential files.
- Project allow rules need workspace trust, which `-p`, SDK and Actions runs never get. `research-cycle` therefore
  lists the cycle's tools in its own `allowed-tools`, which isn't trust-gated. No hooks.
- The cloud environment uses Custom network access: the default list plus five observed science hosts
  (docs/operations.md §3).

**Alternatives rejected.**
- A custom scheduler, daemon or Python orchestrator: it would rebuild routines and `/loop` (D-001).
- Copying the charter or CLAUDE.md rules into the skills: the copies would drift. CLAUDE.md is imported into every
  session.
- `disable-model-invocation: true` on `research-cycle`: it stops `/loop` and scheduled fires from running the skill.
- Dynamic context injection (`` !`cmd` ``) for orientation: a single failing `gh` call, such as the cloud GraphQL
  restriction, aborts the whole skill.
- A PreToolUse hook that parses `git push`: it would match more robustly, but it means a script on two shells for a
  risk that GitHub rulesets block server-side (`main` protection now; an all-branches ruleset is recommended).
- Actions `schedule` as the main scheduler: it runs only from the default branch, GitHub disables it after 60 days
  without activity, and it uses Actions minutes.
- `allowed_non_write_users`, bot triggers, or no OWNER gate on a public repository: prompt injection and spending
  the owner's subscription.
- No `concurrency` in `claude.yml`: one review with several `@claude` mentions would start parallel runs that race on
  one branch. The accepted cost is that GitHub collapses queued runs to the newest one; the cancelled runs are visible,
  and `/research-cycle` step 2 picks up unanswered owner `@claude` comments.
- `curl` in `/reuse-check`'s `allowed-tools`: the prefix grant would accept any extra flags, including uploads.
  WebFetch verifies the links instead.
- A deny rule for `git push * :<ref>` deletions: no clean pattern exists (`:*` is the legacy prefix syntax, and
  `:**` warns at every start). The recommended ruleset's "Restrict deletions" covers it server-side.
- "Full" network access for routines: five hosts plus the defaults are enough.

**Evidence.**
- Claude Code and claude-code-action docs, checked 2026-10-07 with Claude Code v2.1.292 (SOURCES.md "Agent tooling").
  The action embeds its own token in the remote URL and fetches fork PRs via `refs/pull` (`src/github/operations/`
  at tag `v1`), hence the fork skip.
- The allowlist hosts were observed by logging DNS lookups during live MAST, S3, SIMBAD, VizieR, XMatch, NED and Gaia
  calls (astroquery 0.4.11, 2026-10-07).
- `main` protection, read from the GitHub API on 2026-10-07: required checks and PRs, 0 approvals, no force pushes,
  enforced for admins.
- Headless tests in a throwaway repository with a local bare remote. With pushes and `gh` otherwise allowed, the deny
  rules blocked:
  - `push --force`, `-f`, `-uf`, `+ref`, `HEAD:main`, `<branch> main`;
  - `push --delete`, `-d`, `-qd`;
  - `gh pr merge`, `filter-branch`, `cat .env`;
  - `gh issue comment --body-file .env`, `gh api ... -F body=@.env`, `git diff --no-index ... .env`;
  - `gh api -X DELETE .../protection`, the `enablePullRequestAutoMerge` mutation, `gh auth status --show-token`.

  A plain push to `claude/*` went through. Without trust, the project allow rules were ignored, while the skill's
  `allowed-tools` still allowed the `claude/*` push.

**Revisit if.**
- Routines leave research preview with changed limits.
- Routines or the GitHub App can act as a non-admin identity. Then merges get a real server-side guard: required
  approval or a ruleset the agents can't bypass.
- A deny-rule bypass shows up in practice. Then add a PreToolUse hook or the sandbox.
- The WIP cap or cadence leaves the owner's review queue idle or overflowing.
- The cloud GitHub proxy blocks `gh` commands the cycle needs.
- claude-code-action changes major version.

**Revision (2026-10-07, owner decision).** The owner's agents now merge their own PRs. Plain
`gh pr merge` moved from deny to allow; `--admin` and `--auto` stay denied, so every merge waits for the
required checks. Guarded files and `needs-human` PRs still go to the owner. Other people's PRs need the
owner's approving review. See CLAUDE.md "Merge policy and version control".

## D-010 Tools for future milestones (unit 9, 2026-10-07)

The full survey, with versions, licenses and per-option verdicts, is in
[docs/landscape.md](docs/landscape.md).

**Decision.**
- **M2 photometry.** Start from the DJA grizli v7.4 catalogs, which cover both SMACS 0723 and the
  CEERS t021 data. They use matched apertures (not PSF-homogenised) and come with eazy photo-z.
  Validate CEERS colours against the PSF-matched CEERS DR1.0 catalog. Use photutils forced
  photometry (`SourceCatalog(detection_catalog=...)`) only where those colours prove inadequate.
  Don't make our own mosaics. The `ceers-full` DJA files are 250–400 MB, so state the reason
  before downloading them (CLAUDE.md).
- **M2 representations.** Start with Zoobot encoders loaded through `timm` (Apache-2.0 weights).
  Use DINOv2 as the generic-vision control. Use the Multimodal Universe JWST subset as unlabelled
  JWST data, but exclude CEERS and any other field being scored from it. Any embedding must beat
  the sklearn baseline on injection-recovery.
- **M2 detectors.** Keep the M1 sklearn baseline. Add PyOD when more detectors are needed, and
  coniferest for active learning.
- **M3 lensing.**
  - Test candidates against published SMACS 0723 models (Mahler+2022, CC0; RELICS HLSP; Caminha+2022)
    and the Noirot+2023 redshifts. Don't refit the cluster.
  - Lens finding: AnomalyMatch first, then a fine-tuned Zoobot.
  - Galaxy-scale modelling: lenstronomy, with PyAutoLens as a cross-check.
- **Calibration.** Use MAST products. Re-run `jwst` only when a product's `CAL_VER`/`CRDS_CTX` lags
  the operational build, when a top candidate needs a correction that is off by default (1/f,
  wisps, persistence, saturated cores), or when a comparison needs identical catalog settings.
  Record `CAL_VER`/`CRDS_CTX` in the manifests.
- **Tracking.** git + manifests now. Add MLflow (local backend only) at M2 if the run count
  outgrows a results table.

**Alternatives rejected.**
- **AstroCLIP, AstroPT, AION-1:** they are tied to Legacy Survey/HSC bands and pixel scales.
  AstroPT also has AGPL code and CC-BY-SA weights.
- **DINOv3:** gated download and a custom license.
- **anomalib:** built for industrial, normal-only training with pixel-level defect maps.
- **SEP:** duplicates photutils.
- **SourceXtractor++ at M2:** conda-only and heavy. Reconsider at M3 for blended arcs.
- **DVC, DataLad:** our inputs already live in a permanent public archive.
- **Refitting cluster lens models:** published models exist.
- **ssl-legacysurvey, CMU DeepLens:** ground-based training, or no weights.
- **ceers-nircam scripts:** no license, and stale.
- **lychee, markdown-link-check, LinkChecker** for `scripts/check_links.py`:
  - lychee is Rust, Apache-2.0, v0.24.2.
  - markdown-link-check is Node, ISC, v3.15.0.
  - LinkChecker is Python, GPL-2.0, v10.6.0.

  The checker had to be stdlib-only and run before any environment exists, without adding a
  Node/Rust toolchain or a GPL dependency. If CI link checking is wanted later, lychee-action is
  the upgrade path.

**Evidence.** docs/landscape.md (checked 2026-10-07; versions from the GitHub/PyPI APIs; arXiv IDs
checked against their abstract pages), DJA file listing
https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html, JWST build table
https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information,
and the MAST `s_region` of `jw01345-o001_t021_nircam_clear-f200w` (queried with astroquery).

**Revisit if.**
- A JWST-trained encoder or a lens finder with public weights appears.
- DJA colours disagree with CEERS DR1.0 beyond the photometric errors.
- MAST reprocesses program 2736 with build 13.0 (then redo the cross-program comparisons).
- A study needs more than ~50 tracked runs (MLflow).
- Derived datasets of several GB must be shared (DVC/DataLad).
- docs/landscape.md is more than 6 months old when a milestone starts.
- The link checker needs more than small fixes (then switch to lychee).

## D-011 Source-level quality gate before ranking (2026-10-07)

**Decision.**
- `quality.assess_sources` flags every merged source before ranking, and `rank.score_anomalies` sees
  only `quality_ok` sources. Flagged sources stay in `sources`/`quality` tables with a reason and are
  never deleted.
- Weight and edge come from a coarse WHT map of the reference-band i2d (`cutouts.sample_weight_map`):
  - Cells are 1″ (`grid_arcsec`, so SW, LW and MIRI behave the same).
  - Each cell takes the median WHT along one row through it; those rows are fetched as concurrent S3
    byte ranges, about 10 s for a 1.8 GB NIRCam mosaic and 5 s for MIRI.
  - Positions map through the SCI WCS. Edge distance is measured from the source's own fractional
    position, not snapped to cells.
- A map with no positive weight, a band mismatch, or under 10% of sources covered is a gate failure.
  So is a gate that leaves fewer than `min_ranked` sources. In all these cases everything is ranked
  and the run is marked "ungated"; the gate can never abort or empty a run.
- Flags (thresholds are ASSUMPTIONS, set in config):
  - `no_coverage`: off the image or on a zero-weight cell.
  - `low_weight`: below 0.5 × the median positive WHT. This matches unit 4's cutout threshold.
  - `edge`: within 1.0″ of a zero-weight region or the image border.
  - `sharper_than_psf`: ref-band `CI_50_30` in (0, 1.45), tested only for aper50 S/N ≥ 3 so noise is not
    called an artifact. A point source has ~0.5/0.3 ≈ 1.67, and real sources are at least that
    extended (program 2736 F200W: 1st/5th percentiles 1.50/1.58), so anything well below is a detector
    artifact.
- This is an optional stage. If it fails, everything is ranked and the report says "ungated".

**Alternatives rejected.**
- Per-source exact WHT/DQ lookup: i2d files have no DQ, and exact per-pixel reads for thousands of
  sources would read most of the mosaic. Cutouts already measure this exactly, but only for the top k.
- MAST `s_region` footprints: polygons miss internal gaps and low-depth areas, which were the main
  artifact source.
- DAOFind `sharpness` cuts: real compact blends reach 2–3.5, overlapping hot pixels.
- Downloading the full i2d: 1.8 GB per band.

**Evidence.**
- Run `20261007T033832Z-cfcf6032` (final implementation) vs `20261007T020124Z-4bfabaaa`: top-20 sources
  carrying a cutout image-quality flag fell from 20/60 to 2/60 (SMACS NIRCam 6→1, MIRI 10→0, CEERS 4→1).
- The coarse gate agrees with exact cutout flags on 19, 19 and 18 of the previous top 20. It also
  caught a CEERS hot pixel that the cutout flags missed.
- Excluded sources: 28% of SMACS NIRCam, 55% of MIRI (mostly the low-coverage mosaic outskirts) and
  28% of CEERS. A first version, with single-pixel cells and cell-snapped edges, flagged 553 CEERS
  sources as `edge` (now 10). It also called faint noise `sharper_than_psf` (SMACS 28, now 10).
  The code review caught both.

**Revisit if.**
- Science needs sources in shallow regions: then rank them as a separate stratum, or add local depth
  as a feature instead of excluding.
- A band other than the reference band drives colours: then gate per band and set that band's
  features to missing.
- MAST adds DQ to i2d products.

## D-012 Star/galaxy separation: rank stars as their own stratum (2026-10-07)

**Decision.**
- `classify.classify_sources` labels each merged source `star` or `other` from one bulk Gaia DR3 +
  SIMBAD cross-match of the whole catalog (`crossmatch.query_matches`, CDS XMatch, about 8 s per
  sample, 0.5″ radius). As in `crossmatch`'s `is_star`, the nearest SIMBAD/Gaia match decides.
  The source is a star when:
  - that match is a star (SIMBAD otype or Gaia DR3 astrometry), or
  - that match is a Gaia DR3 source within 0.3″ (ASSUMPTION: in extragalactic fields Gaia detections
    are mostly stars).

  A nearer SIMBAD galaxy wins, which keeps catalogued cluster members out of the star class.
- Stars are not discarded. The runner ranks them among themselves as a stratum `<sample>-stars`, with
  its own top k (default 10), cutouts, cross-match, candidates and report section. Unusual stars
  (e.g. brown dwarfs) stay findable, and they no longer crowd out galaxies.
- This is an optional stage, and it can never abort a run or drop sources. In each of these cases
  there is one ranking, as before, with all sources:
  - a SIMBAD or Gaia service failed, or the parameters are invalid;
  - the stratum id would collide with a configured sample;
  - fewer than `min_stars` (5) stars were found, or fewer than 2 other sources remain.

**Alternatives rejected.**
- The pipeline's `is_extended` flag: every bright, saturated star in the program 2736 top 20 has
  `is_extended = 1` (spikes and flat cores inflate CI), and 68% of all F200W detections have
  `is_extended = 0`.
- Excluding stars outright: that would lose unusual stars.
- A size–magnitude stellar locus as the first step: it needs care with saturated stars. It is the next
  step, for faint stars that Gaia misses.

**Evidence.** Runs `20261007T035338Z-8ea21f89` and `20261007T040938Z-01527ace` (nearest-match rule):
- Stars found: SMACS NIRCam 52 (46 astrometry/SIMBAD, 6 position-only), MIRI 12, CEERS 7.
- The galaxy-stratum top 20s contain no cross-matched stars (SMACS had 7 before) and 1/60 image-quality
  flags.
- Visual check (unvetted): about 4 faint PSF-like sources without a Gaia counterpart remain in the SMACS
  galaxy top 20. Most of the rest are interacting, clumpy or elongated galaxies, including three
  arc-like sources (`f200w_2925`, `f200w_2559`, `f200w_1096`).

**Known limitations.**
- No proper-motion propagation: SIMBAD positions are J2000, Gaia J2016, JWST about 2022. Fast movers
  (nearby M and brown dwarfs, >~100 mas/yr) can miss the 0.5″ match and stay in the galaxy stratum.
- Position-only Gaia stars (6 in SMACS) can be compact cluster galaxies.
- The top-k cross-match re-queries SIMBAD/Gaia that `classify` already fetched (about 4 s per stratum).

**Revisit if.**
- Faint stars or brown dwarfs keep reaching the galaxy top k: add a size–magnitude stellar locus.
- High-proper-motion stars show up as galaxy candidates: propagate Gaia positions to the JWST epoch.
- Gaia DR4 is released.
- A field has Gaia-detected compact galaxies in numbers (cluster cores): tighten `gaia_position` with
  Gaia's own classifiers.

## D-013 Colours from matched-aperture photometry; DAOFind statistics only for point-like sources (2026-10-07)

**Decision.**
- Colours come from the DJA v7.4 grizli catalog: the same 0.5″-diameter circular aperture in every band,
  defined on one combined detection image. The pipeline `_cat.ecsv` encircled-energy apertures are no longer
  used for colours.
  - `photometry.py` fetches the catalog by URL plus pinned sha256 into the cache.
  - It turns µJy into AB (23.9 − 2.5 log f). Measurements with SEP flags APER_TRUNC, ALLMASKED or
    NONPOSITIVE are dropped. APER_HASMASKED is kept, because grizli masks neighbours on purpose.
  - Magnitudes are joined to the merged sources one-to-one within 0.2″ (closest first), as
    `<band>_dja05_abmag`. Features then use `aperture: dja05` for colours. `ref_mag` and morphology
    keep the pipeline aperture and its S/N gate.
  - The catalog is verified against its sha256 before it enters the cache (downloads over 200 MB need
    `max_bytes`). The pinned URL and sha256 in the config are its reproducibility record.
  - The pipeline catalogs remain the source list. This is a per-sample, optional stage: if it fails, the
    pipeline colours are used and the report says so.
- `features.build_features(daofind_max_ci=1.8)` sets DAOFind `sharpness`/`roundness` to NaN unless the
  reference-band CI_50_30 ≤ 1.8 (point-like). Both statistics assume a point-source profile.

**Alternatives rejected.**
- Correcting the EE-aperture colours analytically: per-band EE radii differ in arcsec, and the correction
  depends on each source's profile.
- Isophotal colours: the isophotes are defined per band.
- Running our own forced photometry with photutils: it duplicates DJA, which already provides it. Revisit if
  PSF-matching is needed.
- Replacing the source list with DJA's: it would change the uids, the gate and the strata in one step.

**Evidence.** Run `20261007T122054Z-c52935ec` (SMACS NIRCam), `scripts/feature_size_bias.py --compare dja05`, on the same
1,423 rows (F200W−F277W, F200W aper50 S/N > 10). Spearman between log isophotal area and colour:
- aper50: +0.181 (p = 6.7e-12);
- isophotal: −0.312 (p = 1.5e-33);
- dja05: −0.012 (p = 0.66).

Re-measured after the exact pipeline S/N conversion (D-014 review) in run `20261007T124944Z-bbad4ab3`, on 1,418 rows:
+0.180 (p = 7.6e-12), −0.308 (p = 1.6e-32) and −0.010 (p = 0.71).

Only the matched-aperture colour is free of size dependence. 2,729 of 5,254 sources match one-to-one.
Nearest-neighbour matching would have given 191 duplicate assignments.

**Revisit if.**
- PSF differences matter for compact sources (DJA apertures are not PSF-homogenised): then validate against
  PSF-matched catalogs such as CEERS DR1.0.
- DJA releases v8.
- Other fields need it: CEERS `ceers-full` is 250–400 MB (state the reason first). MIRI colours in a fixed
  0.5″ aperture need their own evaluation.

## D-014 Rank only confirmed detections: best-band S/N floor and multi-band confirmation (2026-10-07)

**Decision.** The D-011 quality gate gains two configurable detection-confirmation tests (ASSUMPTION
thresholds):
- `low_snr`: the best aper50 S/N over the bands a source is detected in is below `min_detection_snr` (5).
  S/N is inverted exactly from the pipeline's `abmag_err = 2.5 log10(1 + 1/SNR)` (jwst `source_catalog`).
  The linear 1.0857/err approximation would admit a true S/N of about 4.5. DJA errors use the same
  convention.
- `single_band`: the source is detected in one band only, and no independent detection confirms it. A
  matched-photometry match counts as confirmation, because DJA detects on a stacked multi-band image.

Bands without an error column cannot fail a detection. With no error columns at all, the S/N test is
skipped, which `meta['thresholds']['snr_test']` records. If the two tests leave fewer than
`min_ranked` sources, ranking falls back to the D-011 image tests rather than going ungated. Flagged
sources stay in `quality.ecsv`. The methodology already treats single-band detections as artifact
candidates (snowballs, cosmic rays).

**Alternatives rejected.**
- A diffraction-spike mask around bright stars as the first step: it needs per-observation spike geometry.
  The 5–6 spike and stripe detections in run `20261007T122054Z-c52935ec` were all single-band with no DJA
  match, so confirmation removes them without geometry.
- An S/N floor in the reference band only (run `20261007T123127Z-3b78b3fc`): it excluded red dropouts.
  The published F150W-dropout candidate `f277w_829` (S/N 104–376 in F277W–F444W, undetected at F200W
  and bluer) left the top 20. 1,285 SMACS sources had S/N ≥ 5 in some band but not in F200W.

**Evidence.** Run `20261007T124944Z-bbad4ab3`, SMACS galaxy top 20, visual check (unvetted):
- spike and stripe detections fell from 5–6 to 1, and near-noise sources from about 6 to 0;
- `f277w_829` is rank 4, `f200w_1032` (`[MJR2023] 028.2`) rank 9 and `f200w_438` rank 15;
- 7 of the 20 are catalogued.

Gate pass counts, with the stars that are then ranked separately (D-012) in parentheses:
- SMACS NIRCam: 2,613 of 5,254 (33 stars), so 2,580 galaxies ranked;
- MIRI: 103 of 530 (8), so 95 ranked;
- CEERS: 2,757 of 7,582 (7), so 2,750 ranked.

**Known limitation.** Gates differ per sample. Only SMACS NIRCam has matched photometry, so MIRI and
CEERS drop every single-band detection, while SMACS keeps the DJA-confirmed ones. Science-vs-control
comparisons are therefore biased by sample construction until every sample has confirmation (TASKS).

**Revisit if.**
- Single-band sources matter scientifically, e.g. extreme emission-line objects: rank them as their own
  stratum.
- A MIRI-specific confirmation becomes available, e.g. a NIRCam counterpart.
- Spike or stripe detections that are multi-band appear: then add the spike-geometry mask.

## D-015 Stellar locus: point-like sources with ordinary stellar colours join the stars (2026-10-08)

**Decision.** After the catalogue classification (D-012), `classify.stellar_locus` adds sources that are
both point-like and stellar-coloured in the matched-photometry catalog (DJA). Thresholds are ASSUMPTIONS,
set in config:
- The detection-image half-light radius r50 is within 20% of `r50_psf`. `r50_psf` is the median over
  catalogued (Gaia/SIMBAD) stars with 20 < `mag_auto` < 22.5, which are unsaturated.
- `mag_auto` < 24. Fainter than that, size cannot separate stars from compact galaxies.
- F150W−F444W and F200W−F356W lie within the catalogued stars' 5–95% range, widened by 0.3 mag.
- At least 10 calibration stars are needed. Otherwise the locus is not applied and the report says so.

Point-like sources with unusual colours (brown dwarfs, compact high-z galaxies) are deliberately left in
the galaxy ranking, where they are legitimate anomaly candidates.

**Alternatives rejected.**
- Pipeline `CI_50_30`/`CI_70_50` or `semimajor_sigma`: bright stars' spikes and saturated cores
  inflate them, so stars overlap galaxies (catalogued stars: median CI_50_30 1.75, other sources 1.94,
  heavily overlapping).
- Size alone: 24 new members at mag < 24, with F150W−F444W from −1.43 to +0.71 against −0.85 to −0.52
  for catalogued stars, so compact galaxies get in.
- DJA eazy stellar-template χ²: a further 59 MB download. Revisit if colours plus size prove
  insufficient.

**Evidence.** Run `20261007T130207Z-c8ffc470` (SMACS NIRCam):
- `r50_psf` = 2.86 px from 22 catalogued stars; 16 stars were added (49 in the star stratum).
- Three of the six PSF-like sources in the previous galaxy top 20 (`2553`, `991`, `1456`) moved to the star
  stratum.
- The other three (`940`, a star pair; `2242` and `1571`, blends of a star and a galaxy) have r50
  inflated by neighbours and remain.
- `f277w_829`, `f200w_1032` and `f200w_438` are now ranks 3, 8 and 11; 8 of the top 20 are catalogued.

**Revisit if.**
- Blended stars dominate the galaxy top k: add image-based PSF-spike detection on cutouts.
- Samples without matched photometry need a locus: then use pipeline size and colour.
- Saturation sets in at another magnitude in other fields: adjust `calib_mag_range`.

## D-016 Stellar locus: one-sided size test (2026-10-08)

**Decision.** This supersedes D-015's size test. A source is point-like when `r50_floor · r50_psf ≤ r50 ≤
(1 + r50_tolerance) · r50_psf`, i.e. 0.5–1.2 × r50_psf. Previously the test was |r50/r50_psf − 1| ≤ 0.2. The
calibration (r50_psf from catalogued stars with 20 < `mag_auto` < 22.5), the colour box and `mag_max` 24 are
unchanged. Both thresholds are ASSUMPTIONS, set in config; `r50_floor: 0.8` restores the D-015 rule.

Why one-sided:
- Nothing real is smaller than the PSF.
- The calibration stars' r50 rises with brightness: 5–8 px at `mag_auto` < 20.5 (saturated cores and spikes),
  2.5–3.0 px at 20–21 and 2.1–2.2 px at 21–21.6. The median, 2.86 px, is therefore an upper envelope.
- The floor only rejects noise-like detections. D-011's `sharper_than_psf` already handles artifacts.

**Alternatives rejected.**
- Calibrating only on unsaturated catalogued stars (21–22.5 mag): there are 4–5 Gaia/SIMBAD stars there,
  near Gaia's limit, which is below `min_ref_stars`.
- A magnitude-dependent r50_psf(mag): more parameters, and not needed, since the one-sided bound with the colour
  box already separates the sequence (Evidence).
- Self-calibrating a two-sided band around the r50 mode of stellar-coloured sources: circular, and it would
  drop bright uncatalogued stars whose r50 is inflated (2.5–3.4 px).

**Evidence.** Run `20261007T144011Z-11342d69` (SMACS NIRCam), compared with `20261007T132538Z-96e911bb`:
- Among uncatalogued sources at `mag_auto` < 24 inside the colour box, r50 peaks at 1.9–2.1 px (65 sources).
  Galaxies start above about 2.5 px.
- Of the 75 uncatalogued sources at `mag_auto` < 24 with r50 < 2.3 px, 74 have stellar colours.
- The locus now adds 90 stars instead of 16 (76 past the quality gate). The star stratum grows from 49 to 109,
  and 2,504 galaxies are ranked.
- The star-stratum top 10 is all PSF-like with spikes on visual check, including new members `2054`, `1802`
  and `1005`. `f200w_1874` (TASKS) is now a star.
- In the galaxy top 20, `741`, `2915` and `1828` (not stars) were replaced by `1894`, `737` and `2043`, because the
  reference sample changed.
- Point-like sources left in the galaxy top 20 fall into three groups:
  - fainter than `mag_max`: `2317`, `365`, `432`, `2578`, `2043`, `737`;
  - non-stellar colours: `2804`, with F150W−F444W +3.7 and very red F356W−F444W;
  - the known blends: `1571`, `2242`, `940`.

**Revisit if.**
- Compact galaxies with stellar colours show up in the star stratum on visual check.
- Faint (mag > 24) point sources dominate the galaxy top k: then consider a fainter `mag_max` with a stricter
  colour test.
- A field with a different detection image or PSF is added.

## D-017 Two-epoch vetting from level-3 catalogs with a local frame tie (2026-10-08)

**Decision.** `scripts/epoch_compare.py` tests a candidate for proper motion and variability with two public
level-3 `_cat.ecsv` catalogs of the same filter from different programs:
- Mutual nearest matches within 0.3″ and 60″ of the target, with aper50 S/N ≥ 30 in both, fix the frame tie by
  their median offset. Galaxies do not move.
- The target's two rows must be mutual nearest neighbours. Its residual is compared with the tie's standard
  error plus its own centroid error (`semimajor_sigma` × pixel scale / S/N).
- Magnitudes get the same treatment, so zero-point differences between pipeline versions cancel.
- Thresholds are ASSUMPTIONS, set by arguments. Outputs are `derived`.

**Alternatives rejected.**
- Re-aligning the images (tweakreg/JHAT-style) and re-measuring centroids: it needs full `_i2d` or `_cal`
  downloads, which are GB for NIRCam (CLAUDE.md), for one target.
- An absolute tie to Gaia: there are few unsaturated Gaia stars per field, and they move.
- A global affine fit over the whole field: distortion residuals between programs grow with distance; a local
  median is simpler and enough for one target.

**Evidence.** SMACS `f200w_2804`, F444W, 2736 (2022-06-07) against VENUS 6882 o057 (2026-06-05):
- 64 references; frame offset (7.3, 20.6) mas; tie error (1.9, 2.3) mas;
- target residual 1.9 mas (0.95σ, including centroid errors of 0.4 mas per epoch), Δm 0.004 against a reference scatter of 0.037.

Unit tests recover an injected 30 mas motion and 0.5 mag change.

**Revisit if.**
- Many candidates need it: then batch per program pair and cache the frame tie.
- Targets sit near a chip edge, where local distortion dominates.
- Proper motions below about 1 mas/yr matter.
- Extended targets are tested: then measure motion on PSF-like cores, or across filters. For `2915`, a galaxy,
  the residual was 7 mas in F444W and 34 mas in F150W; the bands disagree, so it is not motion.

## D-018 Diffraction-spike flag on cutouts: orientation-free hexagonal harmonic (2026-10-08)

**Decision.** This is new code (tier 5): numpy and scipy, with no new dependency. `cutouts.spike_statistic` resamples
a NIRCam cutout on rings around the brightest smoothed peak within `search_arcsec` (0.3″) of the target. No-data
pixels there count as brightest, because saturated cores are NaN.
- Each ring loses its azimuthal median. The angular power spectra, summed over rings from 0.2″ to 0.8″, give
  `spike_s6 = P6 / sqrt(P4 · P8)`.
- JWST's six main spikes put power at m = 6 only. Elongated galaxies spread power over all even m, which the
  neighbours cancel. |F_m| does not depend on the spike angle, so no position angle is needed.
- Cutouts with `spike_s6` ≥ 3 get the `spikes` quality flag. The report lists them; ranks are unchanged.
- Radii, search radius and threshold are ASSUMPTIONS in `stages.cutouts.spike`. MIRI is skipped, because its
  cruciform PSF needs its own calibration.

**Alternatives rejected** (`/reuse-check`, 2026-10-08):
- STPSF templates: they need the spike angle and an 88 MB data set, and pin numpy < 2.4.
- photutils PSF fitting: it has no spike utility, and it fails on NaN cores.
- grizli `mask_IR_psf_spikes`: WFC3/IR only, and it works on level-2 exposures.
- DJA: it ships no star or spike masks.
- Catalog-driven geometric masks (LSST, JWST1PASS): they need per-exposure geometry and `_cal` files; D-014
  already deferred these.
- spike-psf: it needs GB of `_cal` files.
- MaxiMask: no JWST data in its training set, and it needs TensorFlow.
- i2d DQ: there is none.
- Centring on the catalog centroid (first implementation): blends scored only 0.4–2.6, because their star is
  offset from the target.

**Evidence.** Run `20261007T162133Z-2511977c` (SMACS NIRCam F200W cutouts), visually checked on both contact sheets:
- **Star-stratum top 10:** all flagged, `spike_s6` 3.9–28.3, including the saturated `1345` (11.6, `nan_center`).
- **Galaxy top 20:** exactly the three PSF-like blends are flagged, `940` 15.6, `2242` 10.9 and `1571` 5.7. The
  other 17 score ≤ 2.3, including faint point sources (`2043` 2.3, `432` 1.9) and arcs.
- **Radius scan** on the previous run's 30 cutouts: r_in = 0.2″ leaves the widest gap between the lowest star
  (3.85) and the highest unflagged galaxy (2.28).

**Revisit if.**
- Flagged sources should leave the galaxy ranking instead of being labelled, e.g. by backfilling the top k.
- Compact bright galaxies get flagged: then check the spike angle `arg(F6)/6` against the header PA.
- Close star pairs or crowded cluster cores are missed.
- MIRI is to be covered.
- MAST adds DQ to i2d products.

## D-019 Spike screening: drop spiky sources with stellar colours from the galaxy top k (2026-10-08)

**Decision.** With `stages.cutouts.spike.screen: true`, a galaxy stratum (not a star stratum) gets cutouts for
its top 2k sources.
- **Removed:** sources flagged `spikes` (D-018) whose matched-photometry colours lie in the D-015 stellar colour
  box (`classify.stellar_colour_mask`, size ignored). They are stars or star-dominated blends.
  - The next clean sources backfill the top k. Original rank numbers are kept.
  - Removed sources go to `screened.ecsv` and a report note.
- **Kept, with a note:** spiky sources with non-stellar colours, because a galaxy with a bright unresolved nucleus
  (an AGN, for instance) also shows spikes and is a legitimate anomaly.
- **Unchanged:** samples without stellar-locus colours are not screened, and the report says so.

**Alternatives rejected.**
- Removing every spike-flagged source: on run `20261007T170058Z-68af3efd` it also removed `2915`, an extended
  galaxy with a bright point-like nucleus and F150W−F444W +2.0 (not stellar), so it may be an AGN host.
- Flag only (D-018 alone): the blends kept 2 of the top 20 slots.
- Moving flagged sources into the star stratum: they were never ranked among stars, and the star stratum's top 10
  would then mix in blends.
- Iterative backfill with repeated cutout calls: a single 2k pool is simpler, and more than k flagged sources
  in 2k is unlikely (the report notes it if it happens).

**Evidence.** Run `20261007T170539Z-b369bfa2` (SMACS NIRCam), contact sheet checked:
- Removed: `1571` (#4) and `2242` (#9), star+galaxy blends with stellar colours.
- Kept and noted: `940` (#10), a saturated star with corrupted colours, and `2915` (#21), a bright nucleus.
- The top 20 now runs to rank 22. The cutout pool costs 20 extra S3 cutouts.

**Revisit if.**
- Saturated stars like `940` should also leave the ranking: use a pipeline-vs-matched magnitude mismatch
  (Δ = 5.1 mag for `940`, beyond the 99th percentile of 3.8) as a saturation test.
- AGN hosts turn out to have stellar-like colours in other fields.
- A sample without matched photometry needs screening.

## D-020 Host test for spiky sources: screen point sources without host light (2026-10-08)

**Decision.** `cutouts.host_ratio` is the azimuthal-median surface brightness in an annulus of 0.3–0.6″ around
the spike peak, divided by the peak, after subtracting the background measured beyond 0.9″. Per-ring medians
ignore spikes and small neighbours.
- With `stages.cutouts.spike.host_ratio_max` (0.004), D-019 screening also removes spike-flagged sources whose
  ratio is below that threshold: a point source with no host light is a star, whatever its colours.
- Spiky sources with host light and non-stellar colours stay ranked (a galaxy nucleus, possibly an AGN).
- The threshold and annulus are ASSUMPTIONS, set in config.

**Alternatives rejected.**
- Saturation tests (aper50 − isophotal magnitude; a no-data core within 0.3″): both failed for `940`, which is
  two stars in one segment, not a saturated star (CHANGELOG cycle 13).
- A centroid–peak offset rule: stars themselves show offsets up to 0.18″ (`1345` saturated, `2054` next to an arc),
  so it would not separate them.
- Fitting a PSF model: it needs STPSF or an ePSF and fails on NaN cores (D-018).

**Evidence.** Run `20261007T190843Z-db3c83fe` (SMACS NIRCam F200W):
- Star-stratum top 10: host ratio 0.0006–0.0021, except the saturated `1345` at 0.0045. Its NaN core lowers the
  peak, so it errs towards "host" (kept).
- Galaxy pool:
  - screened: `940` 0.0017, `2242` 0.0014, `1571` 0.0023;
  - kept: `2915` 0.0076, a catalogued quiescent galaxy with a compact core.
- Galaxy top 20 composition: spikes fall from 2 to 1 (only `2915`) and known objects rise to 40%. The contact
  sheet was checked visually.

**Revisit if.**
- A galaxy-stratum point source with host light under 0.004 is a real compact galaxy or AGN (vetting shows it).
- Saturated stars with NaN cores reach the galaxy stratum, since their ratio is biased high.
- Other fields or bands need their own calibration on the star stratum.

## D-021 Low-weight screening: re-apply the D-011 weight threshold on the cutout (2026-10-08)

**Decision.** With `stages.cutouts.screen_low_weight: true`, a top-k pool source is dropped when its cutout
carries the `low_weight` flag: median core WHT below `low_weight_frac` (0.5) of the image's typical WHT.
- This is the same threshold as D-011's `min_rel_weight`, measured at the source instead of on the 1″
  coarse map.
- It applies to every stratum except star strata, including single mixed rankings, and shares the D-019
  2 × top_k pool and backfill.
- `screened.ecsv` now has a `reason` column: D-019 stellar colours, D-020 no host light, or D-021 low
  weight.

**Alternatives rejected.**
- A finer D-011 weight map for every source: it reads more WHT rows for all sources (I/O), while the top-k
  cutouts already measure the core weight.
- Requiring DJA confirmation for faint short-wavelength-only detections: DJA's detection image may miss
  genuinely blue faint sources, so this risks dropping them. Revisit if low-weight screening proves
  insufficient.
- Raising `min_rel_weight` globally: the coarse map's smoothing, not the threshold, let `3034` through.

**Evidence.** Run `20261007T194733Z-ac5f1b59` (SMACS NIRCam):
- `3034` (#3) is a diagonal two-stripe feature, detected only in F150W and F200W (S/N about 3.6 and 7) and absent
  from DJA. Its coarse-map rel_weight is 0.558, which passed, but its cutout core weight is 0.468, which is
  screened.
- It was the only `low_weight` source in the 40-cutout pool.
- Galaxy top 20 now: 5% cutout-flagged (only `2915`, kept on purpose), 40% known objects, 4 screened.

**Revisit if.**
- Real sources in shallow regions (dither gaps, mosaic edges) are lost; check them with `screened.ecsv`.
- Short-wavelength-only artifacts appear in full-weight regions: then add the DJA-confirmation rule.

## D-022 CEERS control on DJA matched photometry; >200 MB download justified (2026-10-08)

**Decision.** `ceers_t021` gets the same DJA v7.4 matched-aperture photometry as SMACS (D-013): the
`ceers-full-grizli-v7.4-fix_phot.fits` catalog, aperture 0.5″, label `dja05`, radius 0.2″. Its colours,
detection confirmation (D-014), screening (D-019 to D-021) and top-k metric are then built the same way as
in the science field.

**Stated reason for the 250.5 MB download** (CLAUDE.md, >200 MB):
- The control field must share the science field's photometric system, or science-vs-control comparisons
  are invalid (TASKS since cycle 5).
- DJA's `ceers-full` catalog is the only public matched-aperture catalog in the same grizli system.
- A FITS binary table is stored row by row and is not sorted by position, so byte-range reads cannot extract the
  t021 rows.
- The file sits in the gitignored data cache with a pinned sha256. The config sets `max_bytes` explicitly.

**Alternatives rejected.**
- Keeping pipeline aper50 colours for CEERS: they are size-biased (D-013), so the samples would differ in
  construction.
- Our own forced photometry on the CEERS mosaics: GB of `_i2d` downloads and a re-implementation (D-013).
- CEERS DR1.0 alone: it is not on the SMACS system. It stays the validation target (D-013 "Revisit if").

**Evidence.** Run `20261007T202419Z-65728e80` (CEERS t021):
- 3,069 of 7,582 merged sources match DJA one-to-one.
- The D-014 gate passes 3,156, up from 2,757 without confirmation.
- 4 sources were screened from the top 20 (3 low weight, D-021; 1 star without host light, D-020).
- Top 20: 85% known objects (CEERS is heavily catalogued), 0% cutout-flagged.
- The stellar locus does not apply: only 2 catalogued stars at 20–22.5 mag, against the 10 needed.

**Known remaining differences between the samples.**
- Pipeline-catalog features come from jwst 3.0.0 in CEERS and 2.0.1 in SMACS (D-010).
- There is no stellar locus in CEERS.
- "Known object" fractions reflect literature coverage, not anomaly rates.

**Revisit if.**
- DJA colours disagree with CEERS DR1.0 beyond the errors.
- A newer DJA CEERS release appears.
- Disk space becomes tight: delete the cached file and re-fetch it by its hash.

## D-023 Owner priorities: lensing-violation search first, more clusters, faster cadence (2026-10-08)

**Decision** (owner, in chat): put the search for lensing that published models cannot explain at the top. It is
the most direct route to exotic-spacetime signatures, and it stays a hypothesis-generating search (charter).
- **M3 comes first.** Start on SMACS with the pinned Mahler+2022 model, and use only verified literature for
  exotic-lens signatures.
- **More lensing clusters, in parallel.** Abell 2744, El Gordo and Sunrise/WHL0137 run as worktree workers
  (label `batch-clusters`). MACS J0416 has no DJA v7 catalog, so Sunrise replaces it.
- **Time-domain:** use the VENUS second epoch for caustic-crossing transients and variable lensed images.
- **Process:**
  - no fixed pause between cycles; continue until usage limits near;
  - one `/code-review` per PR on its final diff, after all fixes are collected;
  - when the owner decides a needs-human candidate PR in chat, the agent records that on the PR, removes the
    `needs-human` label as instructed, and merges under the normal CLAUDE.md policy, which needs no label
    (PR #20).

**Alternatives rejected.** Continuing the M1/M2 contamination clean-up first: the SMACS galaxy top 20 is now
clean (D-019 to D-021), so further clean-up has diminishing returns.

**Revisit if.**
- The usage limits are reached often: then restore pauses.
- The lens-model comparison is dominated by model systematics: then use several published models per field.

## D-024 Lens-model stage: Lenstool dPIE port, Mahler+2022 ICLv2 validated; SMACS arc orientations (2026-10-08)

**Decision.**
- **`src/jwst_anomaly/lensmodel.py`** evaluates published Lenstool models. It parses `best.par` (French or English
  keywords) and sums its dPIE (profile 81) potentials in one lens plane. At sky positions, for any z_s, it returns
  convergence, shear, signed magnification, deflection, source-plane position and the expected arc orientation
  (`model_prediction`).
  - Anything else raises `UnsupportedModelError` and is never approximated: other profiles, `potfile` catalogs,
    several planes, non-flat cosmologies, relative coordinates.
  - The dPIE deflection and Hessian are ported to numpy from PyAutoGalaxy (`autogalaxy==2026.10.7.1`, MIT; licence
    text in the module).
- **Model: Mahler+2022 ICLv2**, the repository README's final model. It is pinned at commit `f36a41c`, with every
  file verified by sha256 (`SMACS0723_MAHLER22_ICLV2`).
- **`scripts/lens_consistency.py`** has two subcommands:
  - `validate`: compares the port with the published κ map and back-traces `arcs.dat`;
  - `arcs`: compares the orientation of elongated catalog sources with the predicted stretch at z_s = 1, 2, 4
    and at the DJA photo-z.

**Alternatives rejected.**
- PyAutoGalaxy as a dependency: it brings the JAX stack, heavy documentation dependencies and an older scipy pin,
  all for about 100 lines of closed-form maths.
- Lenstool itself: C code that needs compiling, not pip-installable. The port reproduces its outputs (see Evidence).
- Interpolating published maps: only a convergence map is pinned, at D_LS/D_S = 1. Deflection, shear and
  magnification at arbitrary z_s need the model itself.

**Evidence** (2026-10-08, `lens_consistency.py validate` / `arcs`).
- **κ map:** 3000² px of 0.0133″, i.e. the central 40″; every 5th pixel, 360,000 points.
  - Median |Δκ| = 1.6e-5; 95th percentile 1.2e-4; 99th 6.5e-4.
  - Maximum 1.22, in one galaxy core (κ = 11.2) at the map edge (−19.3″, +11.9″).
  - Median ratio 1.0000001. The map is therefore at D_LS/D_S = 1, and MCMC sample `0000` is `best.par`.
- **D_LS/D_S** matches the 16 values written in `best.par` to ≤ 2e-4. System 8's free redshift fits to 11.76; that
  is an unconstrained fit value (D_LS/D_S saturates at high z), not a redshift measurement.
- **Back-trace** of the 60 `arcs.dat` images (redshifts from `arcs.dat` or `z_m_limit`):
  - image-plane rms 0.319″, median 0.249″;
  - χ² at `sigposArcsec` 0.443″ is 31.18, against Lenstool's `Chi2pos` of 30.91 (1%);
  - no image lies beyond 3σ, and parities alternate within each system;
  - the worst systems are 4 (0.61″ rms) and 1 (0.56″).
- **Arc orientations**, from the F200W pipeline catalog (jw02736 o001) and DJA v7.4 eazy photo-z. ASSUMPTIONS:
  - selection: ellipticity ≥ 0.5, σ_major ≥ 2 px, S/N ≥ 10, r ≤ 50″, reduced shear ≥ 0.2 at every tested z;
  - classes: `aligned` ≤ 30° and `anti` ≥ 60° at every tested z. A background source is tested at its photo-z
    range (z16, z_phot, z84); other sources at z_s = 1, 2, 4;
  - background: z16 > z_lens + 0.1. Random orientations are aligned with probability 30/90.

  Results:
  - **Convention check:** 12 elongated (e ≥ 0.3) catalogued images match within 0.5″. Their median offset is 19°,
    and 7 are within 30°. All 4 with e ≥ 0.6 lie within 10°. The two beyond 60° (17.3, 19.1) are small
    (σ ≤ 2.8 px), with e ≈ 0.45 and μ 2.6–6.7.
  - **Strong-shear region, background sources:** 21 of 25 aligned (random expectation 1/3; binomial
    p = 2.6e-7), 1 anti. Members or foreground: 3/6 aligned, 2 anti. All 58 sources: 39 aligned, 8 anti.
  - **The 6 anti candidates** (background or no photo-z) were inspected in F200W cutouts:
    - `1159` and `1896` are blends whose centroid lies between two galaxies. The edge-on spiral in `1896` is
      itself aligned.
    - `1253` is a diffuse galaxy at z_phot 0.77: low lensing efficiency, so an intrinsic shape explains it.
    - `1791`, `1412` and `1804` are low-surface-brightness detections (S/N 23–37) on the BCG/ICL gradient or in
      noisy regions.
  - **Result:** at these thresholds, SMACS 0723 has no anti-tangential arc that ICLv2 fails to explain. A null
    result for this test.

**Revisit if.**
- A model uses other profiles, `potfile` catalogs or several planes: check El Gordo (Caminha+2023) and Abell 2744
  (Bergamini+2023) before porting.
- Shapes: reject blends (several peaks in one segment), require S/N ≥ 50 and compare bands; this cuts the noisy
  `anti` tail.
- Model uncertainty: only the best model is used. The MCMC samples would turn offsets into significances.

## D-025 Field robustness: extended-star veto and red exemption from the host test (2026-10-08)

**Decision.** Two rules found on Abell 2744 (PR #32). Both thresholds are ASSUMPTIONS, set in config.
- **`classify.extended_veto`.** This runs before the stellar locus calibrates.
  - A catalogue star (Gaia/SIMBAD) whose DJA r50 exceeds 4 × 2.5 px goes back to the galaxy ranking; for sources
    fainter than `mag_auto` 21 the limit is 2 × 2.5 px.
  - Gaia lists the cores of bright cluster galaxies, and those matches put whole galaxies into the star stratum.
    They also widened the locus colour box (Abell 2744: F150W−F444W up to +0.06).
  - Bright stars' spike wings reach r50 of about 8 px in SMACS, so the bright limit stays at 10 px.
  - Colour cannot separate stars from cluster galaxies here: cluster ellipticals have F150W−F444W ≈ −0.75, like
    stars.
- **`cutouts.spike.host_exempt_colour` [F150W, F444W, 1.0].** A spiky source without host light (D-020) that is
  redder than this is not a star, so it stays ranked. Abell 2744 `5904` (F150W−F444W +2.64, a catalogued
  emission-line galaxy) was screened from rank 1 before this rule.

**Alternatives rejected.**
- A colour veto: it fails for the reason above.
- Requiring Gaia astrometric significance: two of the extended "stars" have the astrometry basis.
- Dropping the D-020 host test: it correctly screens star pairs such as `940`.

**Evidence.**
- Run `20261007T211229Z-883a5066` (Abell 2744):
  - 8 extended catalogue stars went back to the galaxy ranking.
  - The stellar locus recalibrated to r50_psf 2.07 px from 10 stars (2.13 from 13 before) and added 14.
  - The star stratum fell from 48 to 39; its top 10 is now mostly point sources with spikes.
  - `5904` is kept ranked. 3 point sources without host light were screened.
- Run `20261007T210959Z-1ad9bca8` (SMACS): 1 star vetoed; the stellar locus (90 added) and the galaxy top 20 are unchanged.

**Revisit if.**
- Sources without DJA r50 are extended. Abell 2744 star rank 3, `7694`, is a smooth bright galaxy: fall back
  to the pipeline catalog's size.
- Real stars on galaxies (superpositions) are vetoed.

## D-026 Cluster fields: Abell 2744, El Gordo, WHL0137 (Sunrise) (2026-10-08)

**Decision.** Three lensing clusters run as separate configs with the reference stages and DJA matched
photometry. They were built by parallel worktree workers (PRs #30, #31, #32; label `batch-clusters`).

| Field | Association | Pipeline | DJA catalog | Notes |
|---|---|---|---|---|
| Abell 2744 | `jw02561-o001_t003`, 7 bands | jwst 3.0.0 | v7.2, 233.5 MB | >200 MB, reason stated in the config; core in footprint |
| El Gordo | `jw01176-o241_t012`, 8 bands | jwst 3.0.0, reprocessed 2026-10-01 | v7.0, 17.2 MB | stale MAST sizes fixed by #29; module A is a flanking field |
| Sunrise (WHL0137−08) | `jw02282-o010_t001`, 8 bands | jwst 2.0.1 | v7.5, 28.8 MB | later epochs 2282 o120 (same pipeline) and 6882 o052 |

MACS J0416 has no DJA v7 catalog, so Sunrise replaced it (D-023).

**Alternatives rejected** (per the worker PRs):
- Associations whose footprint misses the core or that lack wide bands (2561 o003_t006, 2756, 4111, 3516,
  6882 o051).
- Merging epochs of different depth or pipeline version (D-010).

**Evidence.** Galaxy top-20 composition on the first runs:
- Abell 2744: 95% known.
- El Gordo: 30% known, 15% lens-related.
- Sunrise: 10% known, 5% lens-related.

Field docs: `docs/fields/*.md`.

**Known differences.**
- DJA versions differ: v7.0 / v7.2 / v7.4 / v7.5.
- Pipeline versions differ (2.0.1 against 3.0.0).
- The stellar locus does not apply in El Gordo or Sunrise (too few catalogued stars).
- "Known" fractions reflect each field's literature coverage.

**Revisit if.**
- Newer DJA catalogs appear.
- Lens-model comparisons (D-023) need other associations or modules.
## D-027 Two-epoch transient search: catalog candidates, verified by forced photometry (2026-10-08)

**Decision.** Two scripts, both producing `derived` outputs with ASSUMPTION thresholds:
- **`scripts/transient_search.py`** compares two same-filter level-3 catalogs. It applies a local median position
  and zero-point tie, then lists `variable` (|Δm| ≥ 0.3 mag at ≥ 5σ), `appeared` and `disappeared` candidates.
  Footprint and depth are judged from catalog proxies.
- **`scripts/transient_forced.py`** checks candidates with **forced aperture photometry**:
  - 0.15″ aperture with a background annulus, at fixed sky positions;
  - both epochs' `_i2d` cutouts read by S3 byte ranges;
  - ERR-based errors plus a 0.05 mag systematic floor.
  Each catalog carries its own segmentation and deblending, and those changed between pipeline versions;
  forced photometry does not depend on them.

**Alternatives rejected.**
- Catalog-only comparison: in SMACS F444W it gave 56 "variables" and 311 "appeared" (F150W: 128 and 410),
  almost all from deblending differences between jwst 2.0.1 and 3.0.0.
- Image differencing: it needs full mosaics (GB) and PSF matching.
- Light-curve exotic tests: two epochs cannot sample them (docs/exotic_lensing.md).

**Evidence.** SMACS 2736 (2022-06-07) against VENUS 6882 o057 (2026-06-05):
- Catalog stage: 441 candidates in F444W and 628 in F150W; 144 coincide in both bands with the same kind.
- Forced photometry: the median Δm is 0.00 in both bands and the 5–95% range is ±0.2–0.3 mag. 4 pass in both bands.
- Visual check: all 4 are bright stars that look unchanged. Their spikes rotate between the epochs (different
  V3 PA), so a fixed aperture on a star's wing measures different spike light.
- After the review fixes (PR #34):
  - Matching follows a global frame tie: shifts of (0.001, 0.022)″ in F444W and (0.008, 0.019)″ in F150W.
  - Catalog candidates: 445 and 630.
  - Re-measuring with WCS-centred apertures: 6 pass in both bands. All look unchanged on visual check (5 bright
    stars and 1 compact source next to a brighter neighbour).
- **Result: no credible transient in the SMACS/VENUS overlap** at |Δm| ≥ 0.3 mag and ≥ 5σ in both F150W and
  F444W. This is a null result for caustic-crossing events at this depth.

**Amendment (2026-10-08, Sunrise 2282 o010/o120, same pipeline jwst 2.0.1).**
- **`scripts/transient_combine.py`** keeps positions with the same `kind` in ≥ 2 bands within 0.3″. It then
  excludes positions near Gaia DR3 sources (one VizieR cone), within `r = 1.5″ × 10^(0.2 (20 − G))`, clipped to
  1.5–12″ (ASSUMPTION). Gaia sources include compact galaxies, so the mask is conservative.
- **`transient_forced.py --recentre-arcsec`** moves the aperture to an iterated centroid in each band and
  measures both epochs at that sky position. `appeared` candidates are centroided in epoch 2, all others in epoch 1.
  Use it for positions taken from another frame or catalog.
- **Evidence:**
  - Catalog stage, same pipeline: 71/65/129/107 candidates in F090W/F115W/F277W/F356W, against hundreds in
    SMACS/VENUS. 58 coincide in ≥ 2 bands, and 6 of them lie near Gaia sources.
  - Forced photometry, recentred (0.15″ box): 0 of 52 candidates pass in all four bands (Earendel is measured as a separate reference). 2 pass in one or two bands, and
    both are ordinary on visual check: one has a linear streak through it in epoch 2, and the other lies on galaxy
    outskirts at the epoch-2 LW footprint edge.
  - Noise calibration on 150 ordinary sources (F356W 25.5–28, median 27.0): the robust std of the ERR-based flux
    significance is 1.18/1.23/1.49/1.36 (F090W/F115W/F277W/F356W). ERR therefore underestimates the noise by
    1.2–1.5×, and 3 of 600 control measurements exceed 5σ.
  - Earendel (known case) at the Scofield+2025 position, without recentring: F356W −0.76 mag at 5.7σ. Its
    centroid lies 0.15–0.16″ away. Recentred, |Δm| ≤ 0.18 mag at |σ| ≤ 1.5 in all four bands.
- **Failed approach:** a fixed aperture at a literature position. An offset of 0.1″ between the aperture and the
  source turns PSF-wing differences into a fake brightening, because the PSF rotates about 180° between the epochs.

**Revisit if.**
- More than two epochs make light curves possible.
- Thresholds are set on calibrated significances: divide by the control std (1.2–1.5) instead of trusting ERR.
- The search goes below catalog depth (image differencing or forced photometry on a grid).

## D-028 Cloud runs open and merge their own PRs with the session's GitHub MCP tools (2026-10-08)

**Decision.**
- Cloud routine sessions block GitHub GraphQL, so `gh pr` and `gh issue` fail. Cloud runs instead use:
  - the session's GitHub MCP tools: `create_pull_request`, `pull_request_read`, `merge_pull_request`;
  - REST (`gh api`) for everything else, labels included, because MCP `issue_write` replaces the whole label set.
- Merges follow CLAUDE.md's merge policy. The run first checks on GitHub:
  - the author and head repository;
  - the branch and labels;
  - every page of the PR's files;
  - the CI jobs (the skipped `claude` runs are not CI).
  It then merges with `merge_pull_request` (squash, `expectedHeadSha` = the reviewed head; docs/operations.md
  §3). When merging is unavailable, the PR gets `merge-ready`.

**Alternatives rejected.**
- A Bash allow rule for `gh api -X PUT .../pulls/*/merge -f merge_method=squash*` (#37, closed). Glob patterns span
  arguments, so the rule would also auto-approve other `pulls/...` writes (review dismissal, PATCH, DELETE).
- Leaving every cloud PR for the owner: the WIP cap of 3 stalls the hourly routine. The owner asked the agents not to
  wait for approval (2026-10-08).

**Evidence.** Routine run `cse_01GjmW75zAM9bsN2ebLZQ9Zh`:
- `gh pr list` returned 403 (GraphQL).
- `mcp__github__create_pull_request` opened #38 without a permission prompt.
- #38 had no check runs on its first three commits, because it conflicted with `main`.
- The github-mcp-server merge tool's head pin is `expectedHeadSha` (checked by the #39 review).

**Revisit if.**
- `.claude/settings.json` gains explicit `mcp__github__*` rules, e.g. denying the file-writing tools. These tools are
  outside the Bash rules today (docs/operations.md §8).
- The cloud proxy allows GraphQL.
- A merge through the MCP tool is refused.

## D-029 Counter-image prediction with an image-plane solver and forced photometry; SMACS null result (2026-10-08)

**Decision.**
- **`lensmodel.DeflectionGrid` and `lensmodel.find_images`** solve the lens equation in the image plane:
  - the D_LS/D_S = 1 deflection is computed once on a grid (default 0.1″ over ±60″) and cached by model sha256;
  - grid triangles are mapped to the source plane, and Newton steps on the analytic model refine each image
    (tolerance 1e-5″).
- **`lens_consistency.py images`** predicts every image of each catalogued system from its mean back-traced source.
  - Each image is classified against `arcs.dat` and a pipeline catalog.
  - With `--forced-image`, forced aperture photometry decides whether each non-catalogued image is there:
    `recovered`, `confused`, `absent`, `undetectable`, `ambiguous` or `no_reference`.
  - Forced-photometry settings (all ASSUMPTIONs): r = 0.2″, annulus 0.6–1.0″, ERR × 1.5, search 1″.
  - The reference is the least-magnified catalogued image at S/N > 5 with |μ| ≤ 50. `recovered` needs ≥ 5σ and
    1/3–3× the predicted flux. `absent` is predicted ≥ 10σ, with best < 3σ or < 1/3 of the predicted flux.
  - Apertures with < 80 % valid pixels are `off_image`. Offsets are West/North arcsec from the WCS.
- **Mahler+2022 κ tarball URL** moves from `github.com/.../raw/` to `raw.githubusercontent.com`. The sha256 is
  unchanged, and cloud runs get a 403 from the former.

**Alternatives rejected.**
- Catalog-only flux references: near cluster galaxies the pipeline segments of the catalogued SMACS images lie
  0.7–2.6″ away, with aper50 magnitudes of 29–32. They wrongly made system 9 `missing`.
- The source-plane back-trace χ² as the only check. It is not valid for models optimised in the image plane
  (El Gordo: 121.6 against Lenstool's 80.22; issue #41).

**Evidence** (SMACS 0723, ICLv2, F200W `jw02736-o001_t001`, run 2026-10-08).
- **75 images predicted for 21 systems.** All 60 catalogued images are reproduced within 0.04–0.91″, and no
  catalogued image is unpredicted. 4 central images are demagnified (|μ| < 0.5), as expected.
- **Forced photometry** of the 11 testable uncatalogued images:
  - **3 recovered:**
    - system 9's third image (μ 4.7) at 0.9″ from the prediction, predicted 11.1σ, flux ratio 1.7. Its
      F150W/F200W/F444W ratios match images 9.1 and 9.2;
    - system 8 (μ 3.4) at 1.0″, ratio 1.25;
    - system 17 (μ −10.2) at 0.7″, 7.8σ against a predicted 24.4σ (ratio 0.36, just above the 1/3 bound). It
      lies on the BCG's light gradient, so it is not a clean recovery.
  - **1 confused:** system 17 (μ 5.4), against a cluster galaxy.
  - **1 undetectable.**
  - **6 no_reference:** systems 11, 16 and 26. Their catalogued images have |μ| > 50 or S/N < 5.
  - **0 absent.**
- **Cutouts** (F150W/F200W/F444W) were inspected for every flagged position.
- **Result: no predicted image of ≥ 10σ is absent.** This is a null result for "missing image" lens-model
  violations in SMACS 0723 at this depth.
- **Ordinary findings:**
  - system 8's model-fitted z = 11.76 is contradicted by F090W detections of 8.2 (6.7σ) and by F150W detections of
    both images. It is a free model parameter, not a measurement;
  - system 26's μ ≈ 27 predicted image (no reference flux) shows at most 4.6σ within 1″ (scratch run, same
    aperture).

**Revisit if.**
- A BCG/ICL model is subtracted. That allows a clean test near the BCG (system 17).
- A reference flux is available for systems 11, 16 and 26 (DJA matched photometry, or a deeper image).
- The test runs on the other clusters: El Gordo and Abell 2744 (issue #41) once their parser bugs are fixed.

## D-030 Lenstool parser conventions and the El Gordo / Abell 2744 models; image-plane validation (2026-10-08)

**Decision.**
- **Image ids** (`lensmodel.image_system`):
  - a trailing lower-case letter names the image: `23a` → system `23`, `1.1a` → `1.1`, `A200.1a` → `A200.1`
    (Caminha+2023, Bergamini+2023b);
  - otherwise the last `.`-part does: `1.2` → `1` (Mahler+2022).
- **`z_m_limit`** takes one or more ids before the flag. Each id is a system (`4.0`) or one of its images (`7a`).
  The last three numbers fix the flag's position, so numeric ids are never read as the flag.
- **Radii given in both arcsec and kpc** may differ by max(2 %, 1e-6″). `best.par` prints 6 decimals, so tiny
  radii carry rounding error (Bergamini potential 37609: 0.000021″ against 2.142e-5″).
- **`load_lenstool_images`** keeps the `a` column as `err_arcsec`.
- **New `MODELS` entries** in `lens_consistency.py`:
  - `elgordo-caminha23`: CDS J/A+A/678/A3, best.par, image list and z = 2, 8 magnification maps;
  - `abell2744-bergamini23`: the authors' page, best.par and image list.
  - All files are pinned by sha256 in `lensmodel.py`.
- **`validate`** now always runs the exact image-plane check (`imageplane_check`, `find_images`). That is
  Lenstool's χ² for models optimised in the image plane.
  - It matches each image to its nearest prediction, Lenstool's convention.
  - A one-to-one pairing lists `shared_partner` images: two catalogued images that share one predicted image, so
    one of them is not reproduced separately (> 3σ). These are lens-model residuals to inspect.
  - It also compares magnification maps when a model publishes them. Pixels are selected on the published map only.
  - The source-plane back-trace is kept only for `sigposArcsec` models (SMACS).
- **Position errors per model:**
  - SMACS: `sigposArcsec` 0.44″;
  - Abell 2744: the per-image errors of its file;
  - El Gordo: 0.621″ for every image. Its CDS file lists "re-scaled" errors, 19 of 56 at 1.242″, which give χ² 52.0
    and do not reproduce best.par's Chi2pos.

**Alternatives rejected.**
- Rewriting the published files (the scratch workaround of issue #41): the pinned sha256 would then no longer
  identify the published product.
- Per-image errors for El Gordo: they do not reproduce the header χ² (52.0 against 80.22).

**Evidence** (`validate`, 2026-10-08; also `tests/test_lens_consistency.py::test_validate_reproduces_lenstool_image_plane_chi2`,
network). All results are `model_prediction` against the published products:

| Field | χ² (image plane) | Lenstool Chi2pos | rms | Other checks |
|---|---|---|---|---|
| SMACS ICLv2 | 30.87 | 30.91 | 0.318″ | κ map median 1.6e-5 |
| El Gordo | 82.53 | 80.22 | 0.754″ (paper 0.75″) | median \|Δμ\|/μ 3.8e-5 (z = 2), 6.7e-5 (z = 8) |
| Abell 2744 | 146.64 | 146.60 | 0.427″ | none |

- **`shared_partner` flags:**
  - SMACS and El Gordo have none.
  - Abell 2744 has two:
    - 34.1a: nearest 0.63″, one-to-one 2.42″. A pair near a critical curve, where the model's merging image has
      μ ≈ 45;
    - 700.1b: nearest 0.79″, one-to-one 14.7″. Its z = 1.217 is a model fit.
  - Both still need vetting (model and redshift error first).
- Grids (ASSUMPTIONs): 0.25″ over ±130″ (El Gordo) and ±190″ (Abell 2744); 0.1″ over ±60″ (SMACS).
- The solve takes 1–2 min per field, and the grid is cached.

**Revisit if.**
- A model uses `potfile`s, other profiles or several lens planes (`UnsupportedModelError` today).
- A published `sigposArcsec` or error convention changes.
