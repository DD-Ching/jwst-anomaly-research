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

## D-010 Tools for future milestones (unit 9)

_Open._
