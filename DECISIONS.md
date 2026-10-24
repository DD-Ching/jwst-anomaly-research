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

**Addendum 2026-10-13 (`arcs` on Abell 2744 and El Gordo).** The port handles both potfile models (D-030), so the
test ran unchanged; two outputs were added to `arcs`: a shape-noise expectation for `anti` background sources
(`anti_probability`: Gaussian intrinsic ellipticity, σ_ε = 0.26 per component, ASSUMPTION; no PSF or blending) and
`--chain-samples`, the fraction of `bayes.dat` rows in which each `anti` source stays anti.
- Abell 2744 (Bergamini+2023b): background 33/35 aligned, 0 anti (0.45 expected). El Gordo (Caminha+2023): 10/12, 1
  anti (0.17 expected; 742 anti in 55 % of chain rows, band-consistent shape). Catalogued images: 35/44 and 14/17
  within 30°. No anti background arc unexplained; docs/fields/abell2744.md and elgordo.md.
- Shapes are still single-band catalogue moments; the blend rejection and S/N ≥ 50 cut above remain open.

**Addendum 2 2026-10-13 (`arcs` on MACS1149 and MACS0416).** Unchanged code; CANUCS DR1 photo-z (zout columns as in
docs/fields/macs1149.md). MACS1149 (CATS v4.1 maps): background 20/24 aligned, 2 anti (0.61 expected); both looked at
and ordinary (4236 a diffuse segment beside a bright galaxy, 3.9σ centre peak after background removal; 1703 at a
critical curve, g 8). MACS0416 (CANUCS Lenstool): catalogued images 63/85 within 30°; background 26/27 aligned, 0 anti
(0.43 expected). Catalogue sources without a CANUCS match within 0.3″ are faint (median m ≈ 28) and near-random in
orientation; they carry no information. Null in both fields. Lesson: a catalogue S/N does not certify a compact
source near bright halos; look at a background-subtracted cutout before calling a shape real.

**Addendum 3 2026-10-13 (`arcs` on Abell 370, CANUCS Lenstool).** Background 47/56 aligned, 4 anti vs 0.94 expected
(P = 0.016). Two sit on critical curves (g_min 1.2 and 4.2), two are intrinsic shapes in weak shear (one is the
resolved spiral image 6.1). Five fields combined: 7 anti vs 2.6 with every source, 4 vs 2.13 (P = 0.17) for g_min < 1.
Rule: the shape-noise null (`anti_expected`) is valid only away from critical curves; quote it for g_min < 1 and
treat g_min ≥ 1 anti sources as direction flips to check against the curve position, not as shape-noise excess.

**Addendum 4 2026-10-13 (second band, `scripts/arcs_bands.py`).** F444W pipeline shapes for the stored F200W `arcs`
sources (nearest F444W source within 0.3″; usable if F444W e ≥ 0.3 and S/N ≥ 10, ASSUMPTIONS). Five fields,
background in strong shear: 90 of 154 usable (60 have no F444W source within 0.3″, mostly m ≈ 27–28; 4 fall below
e 0.3; no F444W source is shared by two F200W ones). Control: 81 of
83 F200W-aligned sources stay aligned, none turns anti; median band-to-band axis difference 1.7–3.9° per field.
F444W anti: 3 vs 1.30 expected (P = 0.14; shape-noise null at the F200W cut e ≥ 0.5): Abell 370 340 and 4540, El Gordo 742, the same axis in both bands
(Δ ≤ 4.5°), single sources in both bands' cutouts. The other four F200W anti (Abell 370 1565, 1788; MACS1149 4236,
1703) have no F444W source within 1″. Null. Rule: band agreement rules out colour blends but cannot separate an
intrinsic shape from a lensing stretch (both achromatic); it is a blend veto, not a lensing test.

**Addendum 5 2026-10-13 (arc curvature, `scripts/arcs_curvature.py`).** Curvature side of long background arcs
(strong shear, semimajor ≥ 10 px, e ≥ 0.5; 43 sources, five fields) from a weighted parabola in the major-axis frame
of the F200W i2d pixels (6″ S3 cutouts; 61 px median-filter halo removal; fit window 3 × semimajor; 3σ, 0.5″,
|c|/σ_c ≥ 3: ASSUMPTIONS), against the model mass direction −α. Raw counts 11 agree, 7 reverse, 24 flat, 1 no arc;
aligned control 9 / 7, so the raw sign is **not a lens-model test**. Reason: for a circular lens c_pred = 1/(2r) is
0.011–0.024 arcsec⁻¹ (sagitta L²/8r ≈ 0.02″ for L = 2″ at r = 25″, under a pixel), while the measured |c| median is
0.056 (4× c_pred): intrinsic bends, knots and blends dominate (looked at: contact sheets of all five fields). Only 5
sources are testable (c_pred ≥ 3 σ_c, formal error; ASSUMPTION): 3 agree (Abell 370 4614, Abell 2744 7247, El Gordo
2447), Abell 370 4613 flat (the window clips this 5″ arc), El Gordo 2527 reverse (|c| 1.9 c_pred, in the bimodal
merger where the circular-centre c_pred does not hold). Not a flag. Rule: a curvature test needs a forward-modelled
prediction (ray-trace an elliptical source at the arc's source position through the model and fit the same
parabola) and is limited to arcs ≳ 3″; the circular-centre sign test is a failed approach.

**Addendum 6 2026-10-13 (forward-modelled curvature, `scripts/arcs_curvature_model.py`).** Each model predicts the
curvature side by ray-tracing a round Gaussian source (σ 0.1″, ASSUMPTION) at the source position of the arc's
centre (photo-z) onto an 8″ grid and fitting the same parabola; checked on an SIS (side toward the lens, |c| = 1/2r).
Of the 18 curved arcs, 5 have a predicted bend ≥ 3σ of the measured error (ASSUMPTION): all 5 bend the way their
model says (Abell 370 3940, 4135; Abell 2744 7247; El Gordo 2447, 2527). **El Gordo 2527 is explained:** the
Caminha+23 model itself bends it away from −α (predicted 12.8σ, cos −0.97 vs observed −0.96), the saddle case
of addendum 5. The other 13 curved arcs have predicted bends < 3σ, so their measured bends are intrinsic and
untestable. Null. Limitation: a round source gives short model images (0.8–5.6″); intrinsic source shape is not
modelled.

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

**Amendment (2026-10-08, calibrated significances; Sunrise third epoch).**
- `transient_forced.py --controls <catalog>` measures a reproducible random set of ordinary sources (default
  200, `aper_total_abmag` 25.5–28, ≥ 1″ from candidates) with the candidates. Per band, the robust std
  (1.4826 MAD) of their significances divides the candidates' significances before thresholding, never by
  less than 1; raw values stay in `*_sigma_raw`.
- **Evidence:** Sunrise o010 against VENUS o052 gives scales 1.30 (F150W) and 1.18 (F444W), within the
  1.2–1.5 found by hand before; 0 of 57 candidates pass. Three controls change; two of them (`n0022`, `n0150`)
  are open candidates that the catalog-stage search did not list (docs/fields/sunrise.md).
- A check run with 40 requested controls (33 measurable in F444W) gave a scale of 0.54, below the 1st percentile
  (0.68) of 33-control subsamples of the 175-control run. The MAD of these heavy-tailed significances is
  unstable at small n, so a band needs ≥ 100 measurable controls (`MIN_CONTROLS`, ASSUMPTION) to be calibrated.
  Controls are selected on epoch-1 F150W magnitudes and applied to every band (a limitation for LW bands).
- **Failed approach:** candidates only from the catalog stage. Requiring the same kind in two bands drops
  transients seen in one band of the shared pair (blue sources in F150W/F444W).

**Addendum (2026-10-15, forced-photometry grid on every catalogued position; Sunrise o010/o120).**
- **`scripts/transient_grid.py`** measures every position of the union of both epochs' pipeline catalogs (all bands,
  duplicates within 0.2″ merged) in every band pair, with `transient_forced.py`'s aperture, annulus, recentring and
  zero point. The noise scale per band is the robust std over all measured positions (never below 1). A position is a
  hit when ≥ 2 bands pass (|Δm| ≥ 0.3, or Δm undefined with positive flux in one epoch; calibrated |σ| ≥ 5) with
  the same sign (ASSUMPTIONs); hits are grouped within 1.5″. Bands run in parallel processes. This
  removes the catalog stage and its two-band candidate rule, so a change in one module (F090W + F115W) is found.
- It reads whole local `_i2d.fits` mosaics (memory-mapped). A field-wide grid touches every row, so byte-range
  cutouts would read more than the file. Download reason (cloud disk rule): 8 mosaics, 6.4 GB in total (SW 1.29 GB
  each), streamed into the gitignored data root and deleted after the run.
- **Evidence:** 6,401 positions, 2,462 measured in at least one band pair (2,331 in all four; the rest lie outside
  the epoch-2 footprint or within ~1″ of an edge). Noise scales 1.48/1.36/1.47/1.43 (F090W/F115W/F277W/F356W), as
  in the control calibration above. 66 hits in 36 groups: 34 are bright-star diffraction spikes, streaks or PSF
  wings that move between epochs; `n0022` and `n0150` are the known persistence afterimages (D-039), recovered
  without any catalog-stage input. **0 candidates.** Contact sheets looked at.
- **Injection** (epoch 2 only, Gaussian FWHM 0.1″ at 500 random positions measured in all four bands, with the
  real grid's zero point and scale; flux max(frac × epoch-1 aperture flux, 10 × median aperture error)): at the
  nominal amplitude (positions where the fraction, not the floor, sets the flux) −0.75 mag is recovered at 100 %
  (198/198) and −0.44 mag at 53 % (131 positions); over all 500 positions 78 % and 39 %. Recovery counts the
  automatic thresholds only: the by-eye vetting and the catalog selection of positions are not emulated, so this is
  an efficiency for catalogued positions before vetting.
- **Failed approaches:** a first pass counted a difference of two negative fluxes (the background annulus
  outshining the aperture inside an extended galaxy) as an appeared/disappeared source; one epoch now needs
  positive flux. Caught by the injection test: recomputing the zero point on the injected positions; it
  absorbs a brightening common to all of them (5 % recovery). Injections must use the real grid's zero point and
  noise scale.

**Addendum (2026-10-15, MACS0416 PEARLS o211/o212; grid validated on published transients).**
- **Decision:** `transient_grid.py` now centroids each position once, in the image of the catalog it came from,
  and measures every band there (`--recentre-in origin`, the default). A hit also counts bands that pass with
  opposite signs (column `mixed`; `--same-sign` keeps the old rule). Both are ASSUMPTIONs.
- **Why (validation on known cases, all `derived`):** the 8 transients Yan et al. 2023 report between PEARLS
  epochs 1 and 2 (ApJS 269, 43; 6 caustic-crossing stars, SN01, SN02) were all missed by the original grid.
  - Their positions sit ≈ 0.20″ from the pipeline frame (offset (+0.196″ RA cos δ, −0.055″ Dec), measured from
    the F356W difference peaks of D21-W1, D21-W3, D21-S1 and SN01).
  - SN01 has a catalog position 0.03″ from it. Per-band recentring in epoch 1, where SN01 is absent in the LW
    bands, moved those apertures onto a neighbour. Its colour also changed (F115W −7.5σ fainter, F356W +6.4σ
    brighter), which the same-sign rule rejects. With both changes it is a hit (and `on_sky` in D-039's test, 31
    clean exposures). The new rules add 1 hit row on this field (SN01) and change no other; positions that
    centroid onto one peak are merged again after refinement.
  - The six caustic stars (F356W 27.1–28.6 on arcs) give |σ| ≤ 3.3 in 0.15″ apertures even at the corrected
    positions: below the grid's depth. SN02 is detectable at its position (F150W +6.2σ) but has no catalog
    entry (nearest position 0.34″, on the host nucleus): a position-set blind spot, as for SN Atalanta (#237).
- **Field result:** 8,322 catalog positions (8,110 after merging 212 that centroid onto one peak), 3,675 measured in
  F115W/F150W/F277W/F356W, 174 hits in 59 groups: 50
  bright-star spikes, streaks or wings (PA 293° → 33°), 6 footprint edges, 2 persistence afterimages (SW only:
  every detection suspect in `persistence_check.py`), SN01. **0 new candidates.** Contact sheets looked at.
  Injections (as above): −0.75 mag 95 % (198/208 nominal), −0.44 mag 46 %. Data: 8 mosaics, 6.0 GB, streamed and
  deleted.
- **o211/o213 (epoch 1 vs 3):** 104 hits in 50 groups: 47 spikes, wings or edges, 2 persistence, SN01; **0 new
  candidates**; injection −0.75 mag 95 %, −0.44 mag 40 %. Neither epoch's catalogs have an entry within 3.2″ of the
  Warhol arc, so its transients have no grid position (`results/transients/macs0416_grid_o213/`).
- **Limit of the grid:** catalogued positions only, 0.15″ apertures; transients on arcs or hosts without their own
  catalog entry, and point sources fainter than ≈ 28 mag on bright arcs, are outside it. Image differencing on the
  arcs would close this.

**Revisit if.**
- More than two epochs make light curves possible.
- The search goes below catalog depth (image differencing or forced photometry on a grid).


**Addendum 2026-10-15b (arc difference imaging).** `scripts/arc_diff.py` differences byte-range boxes of two
level-3 mosaics (resample, global + 3″ tiled registration, Gaussian blur and flux scale fitted on the brightest
10 % of pixels, matched filter). On MACS0416 o211/o212 it recovers 4 of 5 published Ep1→Ep2 arc transients
(Yan et al. 2023) but leaves positive ridges along thin bright arcs, so it is not yet a search. Rejected: a blur
chosen by the global rms (it over-blurs, because blurring also lowers the noise); a single global shift. Reused:
astropy, scipy.ndimage, fsspec/s3fs. Not adopted yet: a fitted pixel-kernel PSF match (HOTPANTS/SFFT-style);
revisit if the Gaussian+scale approach cannot remove the arc ridge.

**Addendum 2026-10-15c (pixel kernels).** The revisit condition held, so `arc_diff.py` now fits a 9 × 9 least-squares
pixel kernel on the box (Bramich 2008 style; numpy lstsq, no new dependency) plus a second, leave-one-stamp-out
kernel for compact-source stamps (`--kernel-px 4 --kernel-stamps`). On MACS0416 it keeps 4 / 5 (Ep1→Ep2) and 2 / 4
(Ep1→Ep3) published arc transients and cuts unflagged hits to the published ones plus 3 (o213). Rejected: one
kernel for the whole box (compact cores +11–28 % in SW); a stamp-only kernel (loses arc transients); one kernel on
stamps ∪ arc pixels. HOTPANTS/SFFT not adopted: the in-script fit passes the validation; revisit if a field with a
spatially varying PSF (box > 1′) leaves ridges.

**Addendum 2026-10-15d (Abell 2744, jw02756 o003/o002).** The first Abell 2744 pair is the DD 2756 visits (2022-10-20,
2022-12-06; six bands, jwst 2.0.1 for both, PA 30° vs 60°). UNCOVER 2561 o006 (2024-07) was rejected: its footprint
(Dec −30.43 to −30.55) misses the core, though MAST's centre coordinate suggests otherwise. Check mosaic corners, not
`s_ra`/`s_dec`. UNCOVER o001 SW mosaics are 5.2–7.7 GB, so `scripts/i2d_region.py` cuts SCI/ERR to a reference
footprint by byte-range rows (for a later UNCOVER pair). The grid is unchanged. The 30° roll makes rotated star
spikes the main hit class. `persistence_check.py` flagged the five epoch-1 SW-only hits as persistence. It returned on_sky / inconclusive for two epoch-2 SW-only hits, which are one NRCB3 pixel seen in two dithers (by hand). The tool tests each sky position alone, so it misses a detector-fixed feature that falls on several positions. Revisit
if a published Abell 2744 transient with coordinates becomes available for validation.

**Addendum 2026-10-15e (Abell 2744, UNCOVER o001 vs jw02756 o002).** UNCOVER 2561 o001 (2022-11-02/04, PA 40.8°) cut to
the o002 (2022-12-06, PA 59.5°) footprint with `i2d_region.py` (SW 373 MB, LW 90 MB per band, about 40 s each; the cut
is clipped at UNCOVER's lower-left corner, so the overlap is partial). The grid is unchanged (4,182 positions, noise
scales 1.22–1.44). Injections: −0.75 mag 94.6 %, −0.44 mag 82.6 %. Result: 46 hit groups, all ordinary, 0 candidates
(`results/transients/a2744_grid_uncover/hit_verdicts.ecsv`). One bright star near (3.5996, −30.4291) produces 21 of the
24 spike groups. Lessons:
- A saturated star's persistence can show in SW **and** LW at one sky position: both channels see the same sky, so
  colour does not separate it from a real source. g4 is detected in all six bands, but only in dither 4, on pixels
  where the star's core sat in dither 3. `persistence_check.py` identifies it (10 of 10 detections suspect); use it on
  every single-epoch compact hit, not only on SW-only hits.
- A single-band change on a star (g11: F115W −0.82 mag, the other bands within 0.02 mag) came from a suppressed core
  in the UNCOVER F115W mosaic. Compare bands from the same visit before reading a change as physical.
- The o002 detector-fixed artifacts (o003/o002 g16/g17) reappear here as g32/g33 against an independent first epoch,
  so both runs classify them the same way.

**Addendum 2026-10-15f (MACS0416 arc tiling, o211→o212; #245).** `arc_diff.py` (pixel kernels, unchanged) over 28
tiles of 20″ holding the 106 CATS v4.1 multiple images outside the Warhol box. Unflagged hits are triaged by the new
`scripts/arc_diff_morph.py` (edge < 0.5″, extended > 3.0″, rosette, dipole; ASSUMPTIONs). It is validated on the 11
Yan+2023 positions (11 / 11 pass in o211/o212) but not yet calibrated by injections through the chain. Every
morphology pass was looked at on contact sheets. The last open ones were vetted with `persistence_check.py` and
`arc_diff_vet.py` against CANUCS o004:
- t+2+2 h002/h003 are o212 NRCB3 persistence: one detection each, in the first dither after F115W exposures that
  saturated those pixels. h003 sits on the same detector pixel as the grid's g40.
- t+2+1 h063 is an F200W-only registration dipole on a 100 nJy knot (−6.7 ± 1.6 nJy; F150W −2.5 ± 2.0). The LW
  change is arc-coherent, and its S depends on the tile box (F150W −6.0 in t+2+1, −3.1 centred).

- t−1−3 h016 (a pass after the extent fix below) is o212 NRCB2 persistence: dither 4 only, after a saturated star in dither 3.
- t+0−3 h064/h067 and t+1−2 h041 lie on one straight o212 trail, which is visible in all five bands.

The extent test first rescaled the smoothed S as white noise (×2√π σ). S is a matched-filter output with correlated
noise, so the "3σ" level sat at about 1σ. Fixed: the smoothed map is now normalised by its robust spread. Review
fixes (same PR): only detected bands (|S_hit| ≥ 3) enter the edge, extent and rosette tests, and the extent test now
uses every detected band (it vetoes the h064/h067 trail by itself). The dipole ratio uses the opposite region's local
peak (within 0.8″). Border targets are judged, not dropped. The 25 stored morph JSONs were re-judged
(`--rejudge`) without new maps, giving 4 new passes. All 4 are ordinary: t+2+2 h007 is a detector-fixed NRCB3 pixel
(persistence) and t−2+0 h036 a detector-edge artifact in o211; t+1+0 h081 and t−2−3 h074 change sign between bands.

Result: **0 candidates over the o211→o212 arcs.** This is not a limit until the morphology veto is calibrated by
injection. Rejected: a 1.0″ extent cut (it removes the published D21-S1/S2), and a dipole test on bands with
|S_hit| < 3 (it flagged 3 / 11 Yan+2023 positions). Revisit if injections show the veto removes more than 10 % of
point transients.

**Addendum 2026-10-15g (MACS0416 morphology veto, injection calibration; #249).** `arc_diff.py --inject-grid`
puts Gaussian stars (27.5 / 28.0 AB) at common sky positions on bright epoch-1 pixels into epoch 2 before detection,
so `arc_diff_morph.py` is calibrated through the whole chain. o211→o212, 11 tiles (incl. the Warhol box), 403 stars: 402 hits, 326
unflagged; the original rosette rule passes 0.82 of them (maps kept: 207 / 251). The vetoing lobes exist without the
star (no-injection control of t−1−3): neighbouring residuals within 0.8″, mostly LW. **Decision:** a rosette counts
only lobes with peak ≥ 0.3 |S_hit| (`rosette_lobe_ratio`, ASSUMPTION): pass 0.92 (231 / 251), none of 244 real-hit
verdicts on the same maps changes, Yan+2023 11 / 11 pass. Rejected: ratio 0.5 (4 real rosettes change, +2 % only);
no ratio (rosettes alone remove 40 / 326 = 0.12 of unflagged injected hits). Old JSONs without `lobe_peaks` keep the old rule
on `--rejudge`. Not covered: ring residuals (t−1−4, 40 real passes, all on rings, not candidates); 7 tiles of the
o211→o212 run have no morph record. Revisit if: injections with real (WebbPSF) PSFs, or a new artifact class.

**Addendum 2026-10-15h (MACS0416 ring and sign-flip verdicts; #250).** t−1−4 re-run with maps: 46 unflagged passes,
almost all on ~0.3–1.2″ PSF-mismatch rings around epoch-1 galaxies. **Decision:** two new `arc_diff_morph.py`
verdicts (thresholds ASSUMPTIONs in `Params`): **ring** = azimuthal fill ≥ 0.5 (16 bins, bin max ≥ max(3, 0.3
|S_hit|), hit's own 0.3″ excluded) of the annulus through the hit around the brightest epoch-1 pixel within 1.2″,
in ≥ 2 detected bands; **sign_flip** = one band S_hit ≥ +5 and another ≤ −5 (SW and LW expose simultaneously, so a
new or vanished point source changes every band one way). t−1−4: 46 → 20 passes (ring 23, sign_flip 3). Injections
in t−1−4 (whole chain, 27.5–29.5 AB, 56 unflagged hits): ring removes 1, sign_flip 2. Yan+2023: 11 / 11 pass
(o211/o212), o211/o213 unchanged (10 pass, D21-S2 rosette). Rejected: sign level 3 (vetoes Yan+2023 D21-S1, F200W
−3.2 against LW +7); ring in ≥ 1 band (3 / 21 bright injected passes lost). Cost: a colour-evolving transient
present in both epochs, or one sitting on a galaxy's ring residual, is vetoed. Old JSONs have no `ring_fill` (= 0).
Revisit if: a transient with a published SED change between epochs is vetoed, or tiles with good PSF matching.

**Addendum 2026-10-15i (Sunrise and SMACS grids under the MACS0416 defaults).** The three grids run with per-band
recentring and the same-sign rule were re-run with `--recentre-in origin` and mixed-sign hits. **Result:** every
earlier hit group but one spike re-found; 5 new rows, all bright-star dipoles, spikes or streaks; SN Atalanta
recovered; 0 candidates (CHANGELOG 2026-10-15). **Decision:** the new defaults stand for all grid runs; the
earlier Sunrise/SMACS nulls hold under both rule sets. Contact sheets across epochs with different roll must
resample one epoch onto the other's WCS. Revisit if: a transient is found only under one rule set.

**Addendum 2026-10-15j (Abell 2744, UNCOVER o001 vs ALT 3516 o002, F356W only).** Single-band grid (`--min-bands 1`),
UNCOVER cut to the ALT footprint with `i2d_region.py`, 6,440 positions, noise scale 1.36; injections −0.75 mag 95.6 %,
−0.44 mag 50 %. Known case: PLATE-23a (Pang et al. 2026, a SN IIP candidate) at RA 3.552756°, Dec −30.400340°;
the paper prints degrees in sexagesimal form ("03h33m09.92s" = 3°33′09.92″), checked against its search centre and
the images. **The catalog grid misses it:** neither pipeline catalog has an entry for it (it is blended into its
20.95 mag host 2.2″ away), so no aperture is placed there. `arc_diff.py` on a 24″ box recovers it at +50σ, 0.06″
from the published position, the strongest hit in the box. **Result:** 45 hit groups; 43 ordinary (spikes and star
cores under an 18° roll, o002-only artifact blocks absent in ALT o001 a week later, UNCOVER F356W-only blocks absent
in the same visit's F277W/F410M/F444W, UNCOVER streaks, dipoles, an edge). Two real changes, both outside the
multiple-image region (≥ 45″ from any Bergamini+23 image), so neither is a lensing anomaly: g2 (3.613638,
−30.416380) is a stationary transient seen in all seven UNCOVER bands (12–43 nJy) and gone in both ALT visits;
g1 (3.562169, −30.400185) is a compact source that brightened by 0.55 mag in F356W in both ALT visits. 0 lensing
candidates. **Decision:** a catalog-driven grid is blind to transients blended into a catalogued host; every
field grid needs a difference-image pass (`arc_diff.py` tiles, `--min-bands 1` for single-band pairs) before it is
quoted as a null for SNe on hosts. Revisit if: a field-wide difference-image tool replaces the grid.

**Addendum 2026-10-16 (Sunrise arc, o010 → VENUS o052).** `arc_diff.py --kernel-px 4 --kernel-stamps` on a 24″ box
around Earendel (holds the arc's catalogued images), F150W + F444W (jwst 2.0.1 vs 3.0.0, 2.9 yr). Injections on bright pixels: 0.86 / 0.83
at 27.5 AB, 0.14 / 0.28 at 28.0. 9 hits; `arc_diff_morph.py` passes 4, all off the arc by eye (bright-core residuals,
an edge-on disk, a compact-core dipole at the box border; all ≥ 3.98″ from every Scofield+2025 multiple image, a post hoc rule): **0 candidates on the arc's
images** (injection recovery 0.86 / 0.83 at 27.5 AB is per band, before the morphology step; not a whole-box limit). Earendel is not a hit (S −5.2 / −1.5). The o010 → o120 pair is
still to run (`results/transients/sunrise_arcdiff/RESULT.md`). New `scripts/arc_diff_sheet.py` makes the hit sheets.

**Addendum 2026-10-16 (SMACS 0723 system 5 arc).** Same settings on a 24″ box between SN Atalanta's images. ERO o001 →
VENUS o057 (F150W + F444W): **Atalanta A and B recovered** as morphology passes (S +49 to +58), the known-case
validation; injections 0.82 / 0.76 at 27.5 AB per band. 73 hits; under the frozen rule (`classify` in
`scripts/smacs_arcdiff_images.py`: morph pass and ≤ 1″ from a Mahler+22 image) 0 to vet beyond Atalanta, **but the
1″ radius fails the known case** (Atalanta lies 1.59–1.69″ from its catalogued images). Post hoc at 2.0″: the one
extra pass (h063) shows no ≥ 2σ change in multi-epoch forced photometry; h061 (on 12.1, box-edge veto) is flat
2022 → 2026-06 but +35.6 ± 3.4 nJy in 2024 F444W, under a stripe by eye; follow-up: `undetermined` under a frozen stripe-subtraction rule (z 2.39), and post hoc matched in size by near-star controls and positions along the same 2024 spike of a bright star 4.7″ away: not a candidate. GO-4043 → DD 12774 (F115W + F444W)
**fails its control** (Atalanta not a hit: F115W S ≈ 2.4; injections 0.26 / 0.34 at 27.5), so it gives no null; its one pass is
AT 2024jfz (0.08″), a one-epoch point source in 2024. `results/transients/smacs_arcdiff/RESULT.md`. **Decision:** a
field's arc pass needs a known-transient (or injected) control that fires under `--min-bands 2`; a pair whose
shallower band misses it is reported as untrusted, not as a null; an "on an image" radius is checked on the
known case before use (images are catalogued at knots, not at every arc pixel).

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
- A reference flux is available for systems 11, 16 and 26 (DJA matched photometry, or a deeper image). Done for
  F200W/F444W peak fluxes in three models (addendum 2026-10-14): untestable at this depth.
- The test runs on the other clusters: El Gordo and Abell 2744 (issue #41) once their parser bugs are fixed.

**Addendum 2026-10-14 (uncatalogued bright arcs, `scripts/arcs_counter.py`).** The counter-image test now also runs
on background, strong-shear arcs of the stored `arcs` tables that lie > 2″ (ASSUMPTION) from every catalogued image:
the arc is traced at z16 / z50 / z84 of its photo-z; a counter-image is robust when all three agree within 3″; the
predicted flux is the arc's isophotal flux × |μ_c/μ_arc| with an encircled fraction of 0.5 (ASSUMPTION); a flag is a
robust D-029 `absent` with best S/N < 5 (ASSUMPTION, from the Abell 2744 control: catalogued image 7247 is `absent`
at 86σ by the 1/3 flux bound). Arcs are excluded near the model's own image list even where D-035/D-044 gate it
off for position tests, and near the HFF CATS list of the same cluster (identification needs only observed positions).
Five fields null (CHANGELOG 2026-10-14); the one earlier flag, Abell 370 3881, is CATS image 7.5 of system 7 with a
wrong photo-z. Rejected: a 1″ exclusion (knots of system 26 with a wrong photo-z made a false `absent`); the model's
own list alone (CANUCS A370 has no image within 7″ of 3881). Both false flags were photo-z errors on known systems,
so a flag needs a spectroscopic redshift or a wider independent-model image search before vetting. Revisit if: an
Abell 370 image list passes D-035, or spectroscopic redshifts for the tested arcs.

**Addendum 2026-10-14 (system 17 near the BCG, `--halo-filter-arcsec`, `scripts/smacs_s17_vet.py`).** `images` can
subtract a median filter (1.0″, ASSUMPTION; `halo_phot`'s `medfilt` route) from every forced stamp, with injected
compact images on a 1″ ring measuring what the filter keeps (`halo_recovery`, simulated). On SMACS F200W it turns
system 17's μ −10.2 image from a marginal recovery (7.8σ on the BCG's light) into `absent` (predicted 35.7σ, best
3.6σ, flux ratio 0.09, recovery 0.98; 0.6″ filter: 25.5σ / 2.7σ). Looked at: nothing within 1″. Every catalogued
sibling predicts ≥ 5× the measured flux (2.4 / 3.0 / 5.4 against 0.48). **Not a lens-model violation:** system 17
has no spectroscopic redshift (arcs.dat z = 0; model z 1.72 Pascale+22, 1.84–1.92 Caminha+22, 2.12 +0.11/−0.09
Mahler+22, photo-z 2.34), and the central pair is a near-caustic prediction: in ICLv2 itself it merges at z ≈ 1.92
(μ −36 / +31) and is gone at z ≤ 1.84; the RELICS Zitrin LTM model (its best z ≈ 1.5) predicts no central pair; RELICS
Glafic (z ≈ 1.5) predicts one 0.9–1.5″ from ICLv2's, also undetected; no published model lists a central image.
The absence constrains the inner mass profile (BCG) of the models that predict the pair. Rejected: the annulus-only
background of the original run (the BCG gradient made the 7.8σ). Revisit if: a spectroscopic redshift for system 17,
or a model with a fitted BCG profile that still predicts the pair at ≥ 10σ.

**Addendum 2026-10-14 (systems 11, 16 and 26, `scripts/smacs_ref_fluxes.py`).** Without the |μ| ≤ 50 cap, each
catalogued sibling predicts the uncatalogued images (peak flux / |μ_sib| × |μ|, F200W and F444W) in ICLv2 and the two
RELICS maps (each at the ICLv2-fitted z; ASSUMPTION). No prediction reaches 10σ from a sibling that passes the other
D-034 gates, except one that is recovered: system 11's pair is undetected (≤ 3.7σ); system 16's counter-images are
≤ 0.6σ (ICLv2, Glafic), while Zitrin's third image 2.5″ from the pair (gated 18σ) is `recovered` by the search (41σ,
×2.4, on a bright arc's edge, so the identification is uncertain) and its far image is 3.9σ; and system 26's
ICLv2 μ 26.6 image is 10.3σ / 10.7σ only from 26.2, whose peak is a neighbour's wing (2.3–3.2× the flux at its
position); from 26.1 it is 3.3σ / 2.2σ, with 4.5σ / 2.9σ found. Glafic and Zitrin give one image where 26.1/26.2 lie (the merging pair is not resolved by the 0.25″
solver grid, so an image count of one there is a grid limit) and no distant counter-image, which the grid would find.
16.1 and 16.2 share one 0.4″ peak (one knot), so their two references are not independent.
**Kept:** the |μ| cap and the gates; the gate-free numbers show the counter-images of these merging pairs lie 6–60×
below the pairs, so the systems are untestable at this depth (apart from Zitrin's recovered image), not
violations. Stamps are not halo-filtered: these systems lie 20–36″ from the BCG. Revisit if: deeper imaging, or a
spectroscopic redshift with a model that puts a counter-image of these systems at ≥ 10σ from a gated reference.

**Addendum 2026-10-14 (`--halo-filter-arcsec` on Abell 2744, MACS1149 and MACS0416).** `images` re-run with the forced
F277W stamps unfiltered and median-filtered at 1.0″ and 0.6″ (models `abell2744-bergamini23`, `macs1149-cats`,
`macs0416-cats`). **No new `absent` in any field:** the only `absent` rows are the ones already explained (Abell 2744
700.1 ×2, MACS1149 16, MACS0416 27 ×2), at the same S/N ratios. The filter makes four Abell 2744 references usable
(3.3, 26.3, A200.1, B200.2: all `recovered`) and moves others to `inconsistent_reference` (MACS1149 1192, MACS0416
122). **It does not decide the rows buried in cluster members:** member cores are compact against a 0.6–1.0″ filter,
so 16 rows stay `confused` (neighbour 4–900× the prediction). Halo-filtered photometry at their exact positions
(looked at, raw and filtered cutouts) is below the prediction in 14 of 16, but every one already has a recorded
ordinary explanation (model copies, CATS-only galaxy-scale caustics, near-critical μ, the 22.1 position residual,
MACS0416 51; MACS1149 13 is the model's copy of the unpredicted 13.2, 1.77″ away), and next to a bright core the
annulus is biased (at-position S/N down to −8σ). Kept: median filtering for smooth BCG/ICL gradients (SMACS 17);
member-buried images need galaxy-model subtraction (isophotes, as MACS0416 51), not a median filter. Also fixed: the
filter floor `3 × FORCED_R_AP` evaluated to 0.6000000000000001 and refused 0.6″. Revisit if: a member-subtraction
route (isophote or Sérsic models of the cluster members) exists for `forced_check`.

**Addendum 2026-10-14 (member-galaxy subtraction, `scripts/member_sub.py`).** Re-measured every `confused`,
`recovered` and `absent` forced row of the three fields (MACS1149 5, MACS0416 18, Abell 2744 8; F277W, same models and
catalogues as the previous addendum) after fitting one elliptical Sérsic (astropy.modeling, not PSF-convolved) per
catalogued neighbour brighter than the prediction, with the search disc masked except the neighbour cores. Each row
gets an empirical noise (ring of residual apertures around the nearest neighbour), injections (recovery ~1 by
construction while the disc is masked) and up to three control discs rotated about the neighbour, where no image is
predicted. Synthetic tests: an image 1.2″ from a Sérsic member is recovered (ratio 1.03), an empty position is
`absent`. **Result: no new `absent` that lacks an explanation.**
- Testable: one row. MACS1149 system 2 (μ −16.5) is `absent` (best residual 0.30 × prediction within 2.25″, at
  24σ_emp; the three control discs give 0.32 ×): the disc is empty after subtraction (looked at, with the prediction
  injected for comparison). Already explained (D-036 vetting): CATS merges 2.2/2.3 and puts this image on empty sky;
  Sharon puts it on the arc 2.6″ away. The margin to the 1/3 threshold is small, so this is not stronger evidence.
- `untestable` 14 rows: member-model residuals in the control discs reach 0.6–660 × the prediction, or the calibrated
  prediction is < 10σ_emp once the core residuals set the noise (MACS1149 21, 1192; MACS0416 15 ×2, 27, 34, 45, 51,
  122 ×2, 132, 133; Abell 2744 8.1, 18.1, 22.1). `blend` 4 rows (a core within 0.5″). `no_member`
  12 rows: no catalogued neighbour brighter than the prediction within search + 0.5″ (the `confused` ones among them,
  MACS0416 1, 27, 45, are confused by uncatalogued light or by sources counted only in the 0.2″ best-aperture search).
- One recorded `recovered` row is doubtful: MACS1149 21's forced best aperture (2.1″ off) lies on a member 2.0″ away.
Single unconvolved Sérsics leave core residuals as bright as the predictions, so this route decides almost nothing
in these fields. Rejected for now: GALFIT/imfit/pysersic with PSF convolution and multi-component members (a new
dependency and per-member tuning for ≤ 14 rows that all have recorded ordinary explanations). Revisit if: a PSF-
convolved member model is needed elsewhere (then rerun `member_sub.py` on these rows), or a member-buried row lacks
an ordinary explanation.

## D-030 Exact image-plane χ² in `validate`; El Gordo and Abell 2744 models in the repository (2026-10-08)

**Decision.**
- `lensmodel.imageplane_residuals` gives each catalogued image's distance to the nearest image that `find_images`
  (D-029) predicts from the family's mean back-traced source. `lens_consistency.py validate` reports this
  image-plane χ² next to the source-plane back-trace, for every model.
- `image_family` implements Lenstool's family rule: trailing letters (`23a` → `23`, `1.1a` → `1.1`), else the
  last `.N` (`4.1` → `4`). `z_m_limit` accepts several image ids on one line. The `_kpc` check allows 1e-6″ on top
  of 2 % (6-decimal rounding of tiny cores).
- Position errors per model (`MODELS[...]["sigpos"]`): SMACS `input.par` `sigposArcsec`; Abell 2744 the image
  list's error column (`forme -10`); El Gordo a uniform 0.621″ (ASSUMPTION: the CDS `best_fit.par` has an empty
  image section, and 0.621″ is the file's smallest error).
- New `MODELS`: `elgordo-caminha23` (CDS J/A+A/678/A3) and `abell2744-bergamini23` (authors' page); SOURCES.md
  has the sha256 values.

**Alternatives rejected.**
- Source-plane back-trace χ² alone: it fails for image-plane-optimised models (El Gordo 121.6, Abell 2744
  far off, against Lenstool's 80.22 and 146.60).
- The per-image error column for El Gordo: χ² 52.0, far from Lenstool's 80.22.

**Evidence** (`model_prediction`; grid 0.25″, cached by model sha256):

| Model | Images / families | Image-plane χ² | Lenstool Chi2pos | rms | max | > 3σ |
|---|---|---|---|---|---|---|
| SMACS Mahler+2022 ICLv2 | 60 / 21 | 30.87 | 30.91 | 0.318″ | 0.91″ | 0 |
| El Gordo Caminha+2023 | 56 / 23 | 82.53 | 80.22 | 0.754″ | 1.54″ (23c) | 0 |
| Abell 2744 Bergamini+2023b | 149 / 50 | 146.64 | 146.60 | 0.427″ | 1.62″ (22.1a) | 0 |

No image-position anomaly in any of the three fields: every catalogued image is reproduced within 3σ.

**Multiplicity residual (Abell 2744, `model_prediction`):** in 3 families, two catalogued images match the same
predicted image (`shared_match`): 3.2a/b (0.81″ apart, μ ≈ 50), 34.1a/b (1.27″, both odd parity, μ −64 and −30)
and 700.1a/b (1.72″, both odd parity, μ −4 and −2). The best model puts no critical curve between each pair, so
it predicts one image where two are catalogued. Each is still within 3σ. Ordinary explanations, not yet tested:
the critical curve's position uncertainty (test with `bayes.dat`), or two clumps of one galaxy catalogued as a
pair. Not a candidate.

**Revisit if.**
- Caminha+2023's σ is found stated otherwise (paper or Lenstool input file).
- A model needs multi-plane lensing or a non-dPIE profile (refused today).

## D-031 Exotic-lens screens: flux ratios (demagnification) and radial arcs around a dark centre (2026-10-08)

**Decision.** `scripts/exotic_screens.py` turns the two testable patterns of docs/exotic_lensing.md into screens.
Each hit is a candidate for `/vet-candidate`, never evidence.
- **`fluxratio`** (demagnification, Ellis-wormhole-like lenses):
  - forced aperture photometry (r = 0.2″, recentred) of every catalogued image in two bands;
  - `lum_ratio` = (flux / |μ|) over the median of the same quantity for the system's other images; its S/N
    includes the reference's error;
  - `underluminous` < 1/3 and `overluminous` > 3, in both bands at ≥ 5σ.
  - Gates:
    - images compact in both bands: f(0.2″)/f(0.4″) ≥ 0.6, from one peak and one background annulus;
    - |μ| ≤ 50;
    - achromatic: the bands agree within 0.5 mag, otherwise `chromatic`.
- **`radial`** (negative convergence):
  - inputs: elongated background sources (e ≥ 0.5, S/N ≥ 10; a photo-z, when present, must put them behind the
    lens) that are `anti` to the predicted stretch and that
    the model does not make radial (radial magnification 1/|1 − κ + γ| < 3 at every redshift of the class: the
    z = 1, 2, 4 grid or the photo-z range);
  - a grid search (0.5″ grid, 1″ line tolerance, 15″ line length) finds connected regions where at least 3
    major axes meet;
  - such a region with no catalog source within 1″ is a dark-centre candidate;
  - significance: the same search, with each arc's position angle redrawn inside its own `anti` window (200
    draws). A uniform 0–180° null would be biased, because the selected arcs point at the mass centre.
- All thresholds are ASSUMPTIONs (the script's module constants and CLI defaults).

**Alternatives rejected.**
- Fixed-aperture flux ratios on all images (**failed approach**, SMACS). A resolved arc's aperture flux follows its
  surface brightness, which lensing conserves, not |μ|. High-μ arcs therefore look underluminous and low-μ
  counter-images overluminous: systems 5 and 10 gave ratios of about 0.4 against about 4. The compactness gate
  removes this.
- Counting every grid point over threshold as a convergence centre: plateaus inflate the count (164 against a
  random 102). One peak per connected region is used instead.

**Evidence** (SMACS 0723, ICLv2, F150W/F444W `jw02736-o001_t001`, 2026-10-08; `derived`).
- **`fluxratio`:** 60 images:
  - 49 resolved in at least one band, 5 untestable, 6 consistent, **0 flags**;
  - with compactness tested in F150W only (first version), system 7 was flagged: 7.1 underluminous (μ 42.7),
    7.3 overluminous. Cutouts show 7.1 and 7.2 are single knots of a long, thin clumpy arc, and 7.3 is the
    compact whole counter-image. That is an ordinary knot-vs-whole aperture mismatch, and 7.1 is resolved in
    F444W;
  - **Result: no demagnification candidate.**
  - **Rule:** require compactness in both bands.
- **`radial`** (F200W catalog, DJA photo-z):
  - 224 elongated sources;
  - 26 dropped as cluster members or foreground (their photo-z is not behind the lens);
  - 34 `anti`, none predicted radial by the model;
  - 12 convergence centres (7 without a catalog source within 1″), against a mean of 9.0 (p95 15) for angles
    redrawn in each arc's anti window;
  - the strongest centre has 4 lines, and the null gives ≥ 4 in 94.5 % of draws.
  - **Result: consistent with chance; no dark-centre candidate.**

**Revisit if.**
- Total (deblended, model-subtracted) fluxes become available for resolved arcs. That would test the 38 resolved
  images.
- A field has many more `anti` arcs, or a centre reaches p_random < 0.01.

## D-032 Flux-ratio and colour test of catalogued multiple images on DJA photometry; no anomaly in SMACS or El Gordo (2026-10-08)

**Decision.**
- Complements D-031's `exotic_screens.py fluxratio` (forced photometry, compact images only): Kron totals
  test resolved images too, and the colour test needs no μ.
- `lens_consistency.py fluxratios` matches every `arcs.dat` image to DJA `fix_phot` (default 0.3″, after an
  optional `--offset-arcsec` frame shift). It computes the implied source magnitude `mag_auto + 2.5 log10 |μ|`
  (μ at the catalogued position and system z) and the 0.5″ aperture F150W−F444W colour, which needs no μ.
- Each image's residual is taken against the median of the *other* usable images of its system (leave-one-out).
  In a pair, both images carry the full pair difference.
- ASSUMPTIONs: `flux_outlier` |residual| > 0.75 mag; `colour_outlier` > 0.3 mag; S/N ≥ 10; |μ| ≤ 20 for the flux
  test.
- An image is untested when its DJA segment exceeds 20,000 px (it swallows host or ICL light), when one DJA source
  matches several images, or when the counterpart's 95 % photo-z interval, widened by 0.15 (1 + z), excludes the
  system redshift (`--photoz`, DJA eazy zout).
- Fluxes and colours come from `photometry.load_dja_catalog`, which applies SEP flags, masks and the unit check.
- Parity (the sign of μ) is reported, not tested.

**Alternatives rejected.**
- DJA `<band>_flux_aper_k × <band>_tot_corr` as a total flux. In v7, `<band>_tot_corr` is 1 and `tot_corr` is a
  point-source correction capped at 1.21. For extended arcs the aperture flux tracks surface brightness, which
  lensing conserves, so it cannot test μ. Only the detection-image Kron `mag_auto` is a total.
- A median that includes the image itself: it halves a pair's difference and hides a single bad image.
- No segment-size cut. SMACS 1.1, an arc on a cluster galaxy's halo with a 32,864 px segment, then looked 1.4 mag
  too bright. Visual check: the excess is host light.
- No photo-z gate. El Gordo 9a's counterpart has z_phot 0.89 (95 %: 0.73–0.99), a cluster-redshift object, for a
  z = 4.32 system, and it made 9a/9c a 1.2 mag colour pair. The 9c counterpart has z_phot 3.62.
- The bare 95 % eazy interval. It is too narrow and excluded plausible images: SMACS 3.3 (1.83–1.86 against
  z = 1.99) and El Gordo 5a (3.70–4.09 against 3.54).
- A wider match radius instead of a frame shift for El Gordo. DJA v7.0 sits at dRA +0.221″, dDec −0.018″ (median
  of 41 matches) from the image list. Matching at 0.5″ without the shift adds the 7b/7c pair, which is not a
  lensing effect (vetting below).

**Evidence** (`derived`; DJA v7.4 SMACS and v7.0 El Gordo, sha256 in SOURCES.md):

| Field (match) | Matched | Flux tested (systems) | Flux rms / max | Colour tested | Colour rms / max | Flagged |
|---|---|---|---|---|---|---|
| SMACS (0.3″) | 22 / 60 | 8 (4) | 0.60 / 1.14 mag | 9 | 0.11 / 0.22 mag | 6.1/6.3 (flux) |
| SMACS (0.5″) | 26 / 60 | 11 (5) | 0.52 / 1.14 mag | 12 | 0.14 / 0.22 mag | 6.1/6.3 (flux) |
| El Gordo (0.3″, shifted) | 37 / 56 | 23 (11) | 0.49 / 1.16 mag | 11 | 0.07 / 0.15 mag | 18b/18c (flux) |

- Photo-z excluded SMACS 8.1/8.2 (model z = 11.76 against z_phot 6.6, as D-029 found) and 11.2, and El Gordo 9a
  and 21b. **No colour outlier** in either field.
- Vetting of the flagged pairs (cutouts and forced 0.2″ photometry on the `_i2d`, S3 byte ranges):
  - **SMACS 6.1/6.3** (μ 15.3 / 3.1, pair difference 1.14 mag in `mag_auto`).
    - Both images are compact, with equal colours (Δ 0.01 mag).
    - Moving 6.1 by the model rms (0.32″) changes μ by only −0.14 to +0.17 mag.
    - Forced photometry reduces the difference to 0.61 / 0.66 / 0.78 mag in F150W / F200W / F444W. 6.2, which has
      no DJA match, agrees with 6.1 within 0.3–0.5 mag.
    - So 6.3 is about 0.6–0.7 mag brighter than the merging pair 6.1+6.2 predicts. That is below the threshold, in
      a pair whose μ ≈ 13–15 depends on the critical curve's position (no `bayes.dat` uncertainty yet).
    - **Not a candidate.**
  - **El Gordo 18b/18c and 7b/7c** (`mag_auto` 25.9–27.4, MUSE Lyα systems).
    - 18b has forced S/N 0.3–4.
    - The 7b/7c ratio changes from −0.16 to −1.48 mag between F200W and F277W, but lensing is achromatic.
    - These are faint-counterpart measurement failures. **Not candidates.**
- **Result: no flux-ratio or colour anomaly** in SMACS or El Gordo among the testable images.

**Limits.**
- Most images are untested. DJA has no segment for arcs inside cluster-galaxy or ICL light (SMACS: 38 of 60 images
  have no DJA source within 0.3″), and SMACS systems 7, 11, 16 and 26 have |μ| > 20.
- μ comes from the best-fit model only.
- `mag_auto` of faint sources (> 25.5 mag) next to bright neighbours is unreliable.
- The 0.5″ aperture colours are not PSF-matched here (the F444W PSF is wider), which may add about 0.1 mag of
  scatter between differently stretched images.
- Only each system's worst image is flagged, so a second discrepant image in a system of four would stay
  `consistent` (none of the tested systems has more than three usable images).
- `_loo_residual` duplicates the sibling reference of D-031's `luminosity_ratios`, with a different reference
  error (median/√n against 1.25 × mean/√n). Merge them into one helper when either changes.
- The S/N cuts invert each catalogue's own error definition: Pogson for SEP `magerr_auto`, and
  `snr_from_mag_err` for `load_dja_catalog`'s 2.5 log10(1 + 1/SNR).

**Revisit if.**
- A photometry with totals for arcs near cluster galaxies becomes available (BCG/ICL-subtracted; e.g. the DJA
  tarball's `_phot_apcorr.fits`, not yet inspected).
- `bayes.dat` μ uncertainties are added. Then use a χ² instead of a fixed threshold, and recheck SMACS 6.
**Addendum 2026-10-13 (Abell 2744).**
- Run on DJA v7.2 (233.5 MB `fix_phot` + 63.6 MB `eazypy.zout`, sha256 as in `configs/abell2744.yaml` and
  SOURCES.md). Reason for the > 200 MB download: the 2561 pipeline catalogues match only 20 / 149 images (core
  deblending), and the test needs the DJA Kron totals; no smaller DJA product carries them. Cloud session:
  downloaded to session disk, deleted after use.
- Result: 34 images flux-tested, 35 colour-tested, 7 flags, all ordinary (docs/fields/abell2744.md). **No anomaly.**
- New limitation: the 20,000 px segment cut misses BCG/ICL halo light for images in small segments on a halo
  (600.1b, 602.1c). A local-background or ICL-subtracted total would remove this class (see "Revisit if").
- Failed approach: forced photometry on a single UNCOVER tile (2561 t003) does not settle faint Lyα counterparts;
  its sign disagreed with DJA for 42.1a/c.

**Addendum 2026-10-13 (`bayes.dat` μ errors).**
- `fluxratios --mu-samples N` propagates the published chain (D-045) into each residual; `resid_z` divides by
  photometric ⊕ model errors ⊕ a scatter floor (median |z| = 0.6745) fitted on the other systems, and
  `flux_tension` flags a system's highest-|z| image at |z| > 3 (ASSUMPTION). The fixed 0.75 mag `flux_outlier`
  stays for continuity. `--inject-mag` calibrates it (flux and S/N scaled): 52–60 % recovery at 1.5 mag fainter,
  62–93 % at 2 mag.
- Evidence (derived, 100 rows): model error median 0.04 mag (Abell 2744) and 0.07 mag (El Gordo), max 0.14 / 0.21;
  σ_sys 0.36 / 0.32 mag (a pair counted once). The chain explains little of the scatter; photometry dominates. The
  tension flags equal the old flux flags (Abell 2744: 4, |z| 6.2–8.5; El Gordo 18b/c, |z| 4.6), all explained. No
  anomaly.
- Alternative rejected: a floor fitted on all systems. A single anomalous system raises it and hides itself (code
  review; a 1.5 mag deficit in a three-image field was not flagged).
- Alternative rejected: a χ² with photometric and model errors only. SEP `magerr_auto` (≈ 0.02 mag) and the chain
  errors leave a 0.3–0.4 mag floor, so it would flag most systems; the floor is fitted per field instead.
- Revisit if: a pinned SMACS chain appears (6.3), or ICL-aware totals lower σ_sys (then the 3σ tension test tightens).

**Addendum 2026-10-13 (halo-subtracted photometry).**
- `scripts/halo_phot.py` measures compact images on `_i2d` stamps after removing the smooth halo: a 1″ median filter
  (main route) and a 0.5–0.7″ annulus plane (second route), r = 0.2″, concentration ≥ 0.6 on the unfiltered
  light, an outlier only when both routes agree (ASSUMPTIONs). Blank apertures on a 1″ ring measure leftover halo light and add their scatter
  to the error; injected PSF-width Gaussians measure what a route removes from a compact image. Built
  from numpy/scipy and `lc.image_stamper`; photutils/SEP are not dependencies, and a two-line filter does not
  justify one.
- Evidence (derived, Abell 2744, UNCOVER t003 F150W/F444W): 602.1c F444W −0.27 / −0.32 mag (`mag_auto` −1.52),
  `consistent`; end-to-end injections (`--inject-mag`, the image's own filtered light scaled) flag 602.1a/b at
  ±1 mag and 602.1c at −1 mag; off-halo 602.1a/b within −0.04 to +0.18; leftover halo in blank apertures (≤ 0.11 mag) is subtracted. 600.1 (a
  resolved spiral) is untestable. Synthetic halo: +3 % / +10 % bias at halo ≈ source, +0.3 mag at 10×.
- Alternative rejected: concentration on the median-filtered image. The filter strips a resolved source's
  extended light, so the 600.1 spiral nuclei passed as compact.
- Alternative rejected: DJA `_phot_apcorr.fits` (inspected: El Gordo header; 405 MB for Abell 2744). Its
  `<band>_tot_<k>` are aperture fluxes times a detection-image Kron/aperture correction, so the halo light that
  inflates `mag_auto` enters the correction too.
- Revisit if: a halo/ICL-modelled catalogue (e.g. bCG-subtracted UNCOVER/MegaScience) is published.
- Three-field run (2026-10-13; CHANGELOG): Abell 2744 18 / 23, El Gordo 4 / 2, SMACS 3 / 6 images tested
  (F150W / F444W); scatter 0.10–0.27 mag on the compact subset; one flag (El Gordo 1c), a thin arc that passed the
  concentration gate. Failed approach: f(0.2″)/f(0.4″) alone does not reject thin knotty arcs; an elongation gate
  is required before the halo fluxes replace `mag_auto` in `fluxratios`.
- Elongation gate (2026-10-13; CHANGELOG): axis ratio q of the Gaussian-windowed (σ 0.15″) second moments of the
  filtered light about an iterated centroid, q ≥ 0.80 in F150W and F444W (ASSUMPTION: between the El Gordo 1a–c
  arcs, q 0.57–0.67, and Gaia stars, 5th percentile 0.91 / 0.94; `scripts/halo_phot_stars.py`,
  `results/halo_phot/stars/`). Tested after the gate (F150W / F444W): Abell 2744 2 / 11, El Gordo 0 / 0, SMACS
  0 / 2; no flag. Alternative rejected: moments about the peak pixel (a point source half a pixel off centre read
  q 0.68 at 0.04″ pixels; caught in review). Revisit if: a band beyond F150W/F444W is needed (calibrate it on stars
  first; `run` refuses uncalibrated bands).
- Halo fluxes into the tension test (2026-10-13; CHANGELOG): `fluxratios --halo-phot <json> --halo-band
  --halo-route` takes `halo_phot`'s tested images per band and route (its gates, instrumental magnitudes at the
  image-list positions) in place of DJA `mag_auto`, and keeps the chain μ errors, the leave-system-out floor and the
  injections. Abell 2744 F444W floor 0.13–0.17 mag (was 0.36). One flagged system (A2744-QSO1 = system 53, 0.55–0.68 mag,
  achromatic) is within published μ systematics (20–40 %, Furtak+2025) and AGN variability across 855 d / 7043 d
  delays: ordinary. Alternative rejected: a separate tension script, which would duplicate the floor and chain code.
  Revisit if: per-image μ systematics are modelled (then a floor from four systems stops understating them).
- Observational parity (2026-10-14; CHANGELOG; `scripts/parity_test.py`): image B's stamp fitted with image A's
  resampled through a free affine map in both handedness families (det > 0 and det < 0); the model predicts
  sign(μ_A μ_B). Per-pair power from synthetics drawn from a smoothed A (σ 2 px; drawn from A itself they share A's
  noise pattern and fake power on a symmetric source). Validated on synthetic pairs (handedness recovered both
  ways; symmetric source untestable). **Not validated on known systems:** on Abell 2744 (Bergamini+23, UNCOVER
  F200W, 1″ median halo filter) 6 of 41 decided pairs are `parity_residual` (31 `consistent`), including
  A2744-QSO1 (53a/b), a point source that carries no handedness, and system 1 (1.1a/b). Without the filter, 6 of 15
  decided pairs were flagged; the cutouts show ICL / BCG gradients and non-affine arc structure inside the 0.6″ fit
  disc. These flags are tool failures, not lens-model residuals; no parity conclusion is drawn.
  Alternatives rejected: peak-S/N template choice without halo subtraction (picks images on the BCG halo).
  Revisit with: member + ICL models (`member_sub.py`-type) before fitting, a fit mask from the source's own
  segmentation, and an acceptance test of ≥ 95 % `consistent` or `untestable` on catalogued, non-near-critical pairs
  before any field run counts.
- Parity acceptance test (2026-10-14; CHANGELOG): **failed; the tool stays unvalidated.** Added a per-stamp
  segmentation (S/N > 2 after 1 px smoothing; neighbours zeroed in A and masked in B), a resolution gate (major-axis
  second-moment σ ≥ 2 × the F200W PSF σ 0.028″, else `unresolved`), cross-family starts from 72 reflection axes,
  denser starts (24 angles, 8 refinements) and a process pool (2.8 min for the field on 4 cores). Abell 2744: QSO1
  and 42.1 become `unresolved`, 600.5a is cleared; of 56 fitted pairs, 27 `consistent`, 20 `untestable`,
  6 `undetermined`, 3 `parity_residual` (1.1b, 7b, 600.3b); acceptance (consistent + untestable + unresolved) 59 / 68 =
  87 % < 95 % (ASSUMPTION threshold). Two blockers: (1) **non-convergence**: between start sets two pairs swap
  `consistent` ↔ `parity_residual` (600.4a Δχ² +32 030 → −52 264 → +3 329; 2.5d −111 → +183), so a class can depend on
  the starts; (2) all 3 persistent flags are predicted-mirror pairs with compact, bright cores, which fits the
  **unlensed PSF** (the free map mirrors A's PSF; B's is not mirrored). A toy asymmetric PSF lowered mirror-pair power
  (Δχ² ≈ 780 → 440–580) but flipped none, so the PSF route is a hypothesis, not shown. Looked at the contact sheet:
  1.1b and 7b are elliptical, nearly mirror-symmetric sources with a core dipole residual in the predicted family.
  Revisit with: a PSF-commuting comparison (fit P ⊗ T(A) against T(P) ⊗ B with an empirical F200W PSF from mosaic
  stars, which is exact for an affine T) and a class only when Δχ² keeps its sign across two start densities.
- PSF-commuting parity statistic (2026-10-14; CHANGELOG): **implemented and rejected as the fix; the tool is paused.**
  P ⊗ T(A) fitted against K ⊗ B with K(d) ∝ P(M d) (exact for an affine map), P a Gaia DR3 star median stack
  (pedestal-subtracted) on the F200W mosaic, PSF-consistent synthetics, class only when 16- and 24-angle starts agree. Validated on a toy
  asymmetric PSF (χ² 1147 for 1185 px vs 34 483 plain). Abell 2744: acceptance 52 / 68 = 76 % after the review fixes (75 % before; plain 87 %);
  1.1b and 7b stay flagged (600.3b `unstable`), so the unlensed PSF is not their cause; 8 new flags, attributed to noise
  correlation from the convolutions (|Δχ²| ×4–8 against a fixed DCHI2_MIN) and non-affine arcs (2.4d). A χ²/n ≤ 2
  gate removes 2 of 12. K is truncated at 0.6″ (exact only for singular values ≥ 0.83). Alternatives rejected: Gaussian or toy PSFs (the empirical stack is cheap: one VizieR cone,
  12 byte-range stamps). Revisit if: a non-affine forward model through the published deflection field, or a
  correlated-noise-calibrated threshold, becomes available.
- MACS0416 halo fluxes (2026-10-14; CHANGELOG; `docs/fields/macs0416.md`): `halo_phot.py` + `fluxratios
  --halo-phot` with `macs0416-canucs` (unchanged code). 5 systems testable (F444W 11 images), max |r| 0.35 mag,
  σ_sys 0.25 / 0.26 mag (medfilt / annulus), no `flux_outlier`, `flux_tension` or colour outlier; 1 mag deficits
  recovered in 36–50 % of targets (simulated). No anomaly in a fourth cluster. Revisit if: a CANUCS MCMC chain is
  pinned (per-image μ errors) or a star calibration on this mosaic moves the q gate.

- Abell 370 halo fluxes (2026-10-18; `results/halo_phot/abell370/RESULT.md`): `abell370-canucs`' list is gated off
  as a whole (D-044), so `halo_phot` and `fluxratios` take `--gate-imageplane` (validate's `imageplane.ecsv`) and use
  only systems with every image solved and image-plane rms ≤ 0.5″ (`lens_consistency.systems_passing_rms`,
  `GATE_MAX_RMS_ARCSEC`, ASSUMPTION ≈ `macs0416-canucs`' global rms): 27 / 45 systems, 64 / 115 images. Only
  system 9 (2 images) passes the photometry gates; residuals ≤ 0.23 mag, no flag; injections recover a 1 mag
  change in ≤ 1 / 2 targets. No anomaly in a fifth cluster, with almost no power there. Alternative rejected:
  opening the whole list (global rms 2.32″, max 17.9″; μ near the worst systems is not the fitted model's).
  Revisit if: an image-plane Abell 370 model with fitted redshifts is published (it would open D-035 too).

## D-033 Lens models from published deflection maps; WHL0137 (Sunrise) RELICS Lenstool (2026-10-08)

**Decision.**
- **`lensmodel.MapLensModel`** evaluates a lens model from two published deflection maps (arcsec, D_LS/D_S = 1).
  - The maps must be on a north-up, east-left TAN grid. Rotated grids raise `UnsupportedModelError`.
  - Model-frame positions go to pixels through the maps' WCS (TAN), not a flat offset.
  - Deflection is interpolated bilinearly. The Hessian comes from centred finite differences in float64, so κ, γ
    and μ are resolution-limited at critical curves.
  - The frame origin is the map's reference pixel unless the `MODELS` entry gives a `centre`. The radial screen's
    `--max-radius` is measured from that origin.
  - Screens fetch only the two deflection maps; `validate` also fetches the κ and μ check maps.
  - It has the `LensModel` interface (`fields_xy`, `deflection_xy`, `kappa_xy`, `evaluate`), so `find_images`,
    `DeflectionGrid` and the exotic screens run unchanged.
  - Fields with only published maps (RELICS, HFF, UNCOVER) therefore need no Lenstool file.
- **`lens_consistency.py`:**
  - map models are `MODELS` entries with `kind: maps`;
  - `load_model` returns a Lenstool model or a map model;
  - `validate` compares a map model with the published κ map and magnification maps.
- **First map model:** `whl0137-relics-lenstool`, the RELICS Lenstool v1 maps of WHL0137-08.
  - The 4 maps are pinned by sha256 in `lensmodel.WHL0137_RELICS_LENSTOOL`, 100 MB each (5000² float32 at
    0.04″), each under the 200 MB limit.
  - z_lens 0.566, H0 70, Ωm 0.3.

**Alternatives rejected.**
- Interpolating the published κ/γ maps directly: they have no source-redshift scaling and no deflection, so they
  cannot solve the lens equation.
- Re-fitting a parametric model: that is new modelling, not the published model.

**Evidence** (`validate`, 2026-10-08; `model_prediction` against the published products).
- κ: median relative difference 3.5e-5 (p95 1.3e-4) over 17,991 pixels with 0.05 < κ < 2.
- μ at z = 6.2: median relative difference 2.0e-5 (p95 7.9e-5) over 17,133 pixels with |μ| < 10.
- These confirm the sign convention (+x along +i = West), the D_LS/D_S = 1 normalisation, z_lens and the
  cosmology.
- Unit test: maps sampled from an analytic dPIE reproduce its deflection (2e-3″), its Hessian (5e-3) and its image
  positions (0.01″).
- **Radial screen on Sunrise** (F200W `jw02282-o010_t001`, DJA v7.5 photo-z):
  - 265 elongated sources, 74 not behind the lens, 30 `anti`, 29 not model-radial;
  - 3 centres, all without a catalog source within 1″, against a null mean of 2.2 (p95 5);
  - max 3 lines, p = 0.885.
  - **Null result.**
- `fluxratio` needs a multiple-image list, which the RELICS HLSP lacks: not run.

**Revisit if.**
- A Sunrise image list becomes available (Welch+2022, Scofield+2025). It would allow `validate` χ² and
  `fluxratio`.
- Map resolution limits the Hessian near critical curves (0.04″ pixels).

**Addendum (2026-10-13, Scofield+2025 maps).** Not pinned as a `MapLensModel`: the four M5-G1 maps are 1.1 GB each
(11792², float64, 0.02″), and their deflection unit is undocumented (an empirical 0.0218–0.0223″ per unit reproduces
the κ map and the constraint images; docs/fields/sunrise.md). `scripts/sunrise_scofield.py` reads per-image patches
by byte range instead. Revisit when the authors state the unit (then pin a cutout as a map model).
**Addendum 2 (2026-10-13).** The κ-matched scale reads as 0.02″ pixels scaled to z_ref = 3.02 (16–84 % 2.86–3.16)
instead of the stated 6.2; `scripts/sunrise_zfit.py` uses s = 0.02″ × r(z)/r(z_ref) per source redshift (ASSUMPTION;
docs/fields/sunrise.md). Rejected: fitting the scale per system (it absorbs the redshift).
**Addendum 3 (2026-10-13).** All four RELICS v1 models (Glafic, Lenstool, WSLAP+, Zitrin LTM) are read per image by
byte range (`scripts/sunrise_relics.py`, not pinned: no full download, so no checksums). Their |μ| reproduces the
published z = 6.2 maps (median ratio 0.998–1.000), system 1 returns z = 5.8–6.3, and system 2 maps to one source at
z ≈ 2.6–3.1 (source-plane rms 0.11–0.23″ against M5-G1's ≥ 0.62″). System 2's M5-G1 misfit is that model's, not an
anomaly. Rejected: pinning the 100 MB maps per model (patches suffice; the Lenstool maps stay pinned for screens).
**Addendum 4 (2026-10-14, counter-images).** `scripts/sunrise_counter.py` reads each RELICS model's deflection maps
over a ±52″ box by byte range (in memory, nothing stored), solves every image of systems 1 and 2 with `find_images`
and applies D-029 forced photometry with an injection check and the D-041 cross-model rule. **0 flags:**
no `absent` counter-image is predicted by all four models. Glafic's two `absent` arc-end images (1.1, 1.2; |μ| 39, 79;
injections recovered) are predicted by 1–2 of 4 models and sit 3–5″ from the c images Glafic misses; the arc's
cross-cluster counter-image is predicted at S/N ≤ 1.7. Control through the same chain: 16 of 20 catalogued images
`recovered` leave-one-out, none `absent`. The D-041 copy rule is not applied as a class: with the re-predicted
images' rms it cannot fire beyond the 2″ known-image radius, and with the misses included it widens itself; the
distance to the nearest missed catalogued image is recorded instead. Rejected: a flag
from one model's near-critical image (the HST-era models re-predict only part of each system-1 knot's images: 1 of 4 to 5 of 5).
Revisit when a JWST-era model fits every knot (M5-G1 with a stated deflection unit).

## D-034 Counter-image and radial-screen rules from El Gordo (2026-10-08)

**Decision.** Four rules, each from an ordinary false flag in El Gordo (docs/fields/elgordo.md), so later fields skip
them automatically. All thresholds are ASSUMPTIONs.
1. **Frame offset:** a `MODELS` entry may give `frame_offset_arcsec` from its image list's frame to the JWST frame.
   `apply_frame_offset` shifts the model's reference point and the image list together at load, so `images`,
   `arcs` and both screens work in the JWST frame.
   El Gordo's is (+0.224″, −0.016″), the median offset of the Caminha image list from DJA v7.0 (issue #41). An
   unshifted 0.2″ aperture missed the images.
2. **Reference images** (`forced_check`) must have:
   - S/N > 5 at the catalogued position;
   - a 0.4″-recentred peak at most 2× brighter (else the peak is a neighbour);
   - compactness f(0.2″)/f(0.4″) ≥ 0.6 (a resolved arc or a galaxy wing does not scale with |μ|; D-031).
   When qualifying siblings differ in f/|μ| by more than 3×, the system is `inconsistent_reference`. That is a
   question for the flux-ratio screens, not a missing image.
3. **Search radius** = max(1″, 1.5 × the system's largest catalogued-image residual), and 2.25″ when one of the
   system's catalogued images is in the `unpredicted` list. The widening is capped at 3″; a larger user radius is
   kept. It used to be a fixed 1″.
4. **Radial screen:** elongated sources within the spike radius of a point source brighter than F200W 20 mag are
   dropped when their axis lies within 7° of the direction to it and that direction is one of the field's spike
   axes. The axes are a hexagonal set θ, θ + 60°, θ + 120° plus θ + 90°, with θ the mode of the aligned pairs'
   angles mod 60°. Spikes point radially at their star; genuine radial arcs off the spike axes survive.
   - Spike radius: 3″ × 10^(0.2 (20 − m)), clipped to 3–20″ (the D-027 mask form).

**Alternatives rejected.**
- A catalogue-counterpart requirement for references: it would tie the check to one catalogue's deblending.
- An S/N cut for spike stars: S/N saturates for bright stars; 553 of 613 sources were dropped.

**Evidence** (2026-10-08).
- **El Gordo forced photometry:** `absent` went from 3 to 0. All three were ordinary on vetting (23: model position
  error; 6: reference contamination plus far-image position error, with photo-z-consistent counterparts at 2.6–3.1″
  and flux ratios 0.9–1.0; 7: frame offset).
- **El Gordo radial screen:** the strongest centre went from 6 lines (p = 0.055, a mag 15.8 star) to max 4 lines,
  p = 0.64. 33 spike segments were dropped within the screened 110″.
- **SMACS re-run:** 0 absent, as in D-029. System 9 becomes `no_reference` (its catalogued images are resolved);
  radial max 4 lines, p = 0.965 (5 spike segments).
- **Tests:** `tests/test_lens_consistency.py` (neighbour peak, inconsistent references, residual radius) and
  `tests/test_exotic_screens.py` (spike segments; off-axis radial arcs kept). There are also tests for resolved
  references and the frame offset.

**Revisit if.**
- `bayes.dat` posteriors give model position errors for far images (system 6 needs about 3″).
- A field's image list carries its own frame solution.

**Addendum 2026-10-22 (MACS0417, `results/lens_consistency/macs0417_arcref/`).** Pipeline isophotal segment totals
(F = `isophotal_flux` of the segment whose 3σ ellipse contains the image, L = F / Σ|μ|) were tested as an arc-aware
reference and **rejected for this field**: frozen gate C2 failed, 10 / 27 testable controls recovered or confused
(0.37; arc-referenced 6 / 12) against ≥ 0.8. Post hoc (looked at): single segments of 3σ semi-major 9.2–65.1″ (three of them
16.2–17.6 mag) hold several lensed knots, or an image together with other light. Do not reuse pipeline
segment totals as references in cluster cores; an arc-aware reference needs a deblended arc catalogue
(e.g. CANUCS v2) or an arc-aperture sum on member-subtracted images.

**Addendum 2026-10-22b (MACS0417, `results/lens_consistency/macs0417_canucsref/`).** The CANUCS DR1 photometry
catalogue (detection on BCG-subtracted images, deblended; F200W Kron totals) **passes** the same frozen gate C2:
37 / 40 testable controls recovered or confused (arc-referenced 10 / 10). The pass is unchanged under four readings of the
catalogue's PA convention (post hoc, `pa_sensitivity.json`); the CANUCS-to-pipeline frame is assumed, not tested. Use it, not pipeline segments, as the
reference and search catalogue for the extra-image test on CANUCS fields (gate C2 per field first). Its two flags
(family 10) have post hoc ordinary readings: one lies on the BCG below the local aperture noise (field-median depth
does not apply there), the other has a photo-z-matched source of the predicted flux 1.49″ away. A local-noise
`undetectable` rule would be a new frozen rule for the next field, not a re-reading of this one.

**Addendum 2026-10-22c (MACS0417 family 10, `results/lens_consistency/macs0417_f10_zscan/`).** A frozen z scan in the
CANUCS v1 maps puts the third family-10 image within 1.0″ of CANUCS 1115938 (with 10.1 / 10.2 reproduced) only at
z 1.975–2.100 (closest 0.464″ at 2.025), above 1115938's z160–z840 interval 1.72–1.88 (offset there 1.96–3.95″). Frozen
outcome A′: the `candidate` flag of addendum b is explained by source redshift only at a z its photometry disfavours;
undetermined, not a candidate. A frozen colour test of 1115938 against 10.1 / 10.2, or a family-10 z_spec, decides it.

**Addendum 2026-10-22d (MACS0417 family 10, `results/lens_consistency/macs0417_f10_colour/`).** The frozen catalogue
colour test (CANUCS DR1 0.3″ colour fluxes, 8 JWST bands, 5 % floor, flux ratio free) passes its known-case gate
(8 / 9 other-family pairs consistent) and finds 1115938 consistent with 10.1 (p 0.90) and 10.2 (p 0.78): "H1 not
excluded". Post hoc: 7.5 % of similar-flux field sources pass against both; the fitted flux ratios (4.69, 2.78) match
the CANUCS-map |μ| ratios at z 2.025 (4.35, 2.78). Family 10 is closed as ordinary (third image ≈ 1115938 at z ≈ 2.0;
photo-z 0.1–0.2 low) unless a z_spec contradicts it. The colour test is reusable for other CANUCS counter-image checks.

**Addendum 2026-10-22e (MACS1423, `results/lens_consistency/macs1423_canucsref/`).** The same CANUCS DR1 reference
with a frozen local-noise rule (F_pred ≥ 5 σ_loc of 0.3″ apertures on the BCG-subtracted F200W, absent needs ≥ 10 σ_loc)
fails gate C2 on MACS1423: all 9 / 11, arc-referenced 2 / 4. The local-noise rule changed no class there. Post hoc:
3.2's own arc source has its centroid 1.02″ away (outside the 1.0″ radius) and 1.31× the predicted flux; 4.3's
reference is the central image beside the BCG. Decision: the CANUCS DR1 reference is not validated beyond MACS0417;
use it on another field only with a new, frozen counterpart rule gated on a third field.

**Addendum 2026-10-22f (MACS1149, `results/lens_consistency/macs1149_canucsref/`).** Third-field gate. The CANUCS v1
MACS1149 maps (`macs1149-canucs-maps`) pass gate K on the 23 Sharon gold images (rms 0.990″, 22 / 23). With the
CANUCS DR1 reference, the local-noise rule and a new containing-source rule (a source whose 3σ ellipse contains the
position also enters the comparison), gate C2 passes: all 15 / 15, arc-referenced 5 / 5. There are 0 flags among 6
targets. The containing-source rule changed no class on MACS1149. Post hoc it recovers MACS1423 3.2, but MACS1423
still fails (arc 3 / 4, BCG-adjacent 4.3). Decision: by the frozen plan, use the CANUCS DR1 reference + containing-source
rule on Abell 370 and MACS0416. The rule is validated only post hoc (MACS1423), and C2 passes on 2 of 3 fields.
Revisit if: a fourth field fails C2, or the rule turns a control `confused` on an unrelated bright source.

**Addendum 2026-10-22g (Abell 370, `results/lens_consistency/abell370_canucsref/`).** Gate K of the frozen plan on the
pinned `abell370-canucs` Lenstool model (115 images, z fixed): rms 2.34″ > 1.0″ (95.7 % within 3″), so it fails and no
prediction or flux was read. It reproduces D-044's 2.32″ (source-plane fit). Post hoc, five images (7.1, 5.3, 3a.2,
14a.3, 5.2) carry the rms (0.60″ without them); not used, because they were picked after the residuals were seen.
Decision: the CANUCS-reference extra-image route has no usable Abell 370 model; the line stops. Revisit if: an
image-plane-fitted Abell 370 model with public parameters or maps is pinned, or a family-restricted rule is frozen and
validated on MACS0416 / MACS1149 first.

**Addendum 2026-10-22h (Abell 370, BUFFALO model A, `results/lens_consistency/abell370_buffalo/`).** The pinned
BUFFALO Niemiec+2023 model A `main` deflection maps (`abell370-buffalo-a-maps`) pass gate L (κ from the maps vs the
published κ, median relative difference 3.4e-5). Gate K on Lagattuta+2019's 102 gold z_spec images (K1, frozen)
fails: rms 7.70″, 92.2 % within 3″. The added gate K2 on the model's own released gold list (104 images; a deviation,
frozen after K1) passes: rms 0.910″ (published 0.90″), 99.0 % within 3″ (28.1 at 3.36″). Decision: this model is
validated on its own constraints only and may serve the extra-image test on the BUFFALO gold systems; the K1 list is
not used with it. Revisit if: the extra-image stages need images outside the BUFFALO gold list, or K1's failure is
traced to the model rather than to list differences.

**Addendum 2026-10-23i (gate C2 shared-pair counting, frozen before the next field).** On Abell 370 (#468) three
BUFFALO systems used one CANUCS (reference, counterpart) segment pair, so one measurement counted as three gate C2
controls (#469, post hoc). Decision: every field this route runs on after Abell 370 uses
`m0417_arcref.gate_c2_shared`: arc-referenced controls of different systems that share one (reference, counterpart)
pair count once; the merged class is the first present in the order absent, recovered / confused, ambiguous, other
(an untestable member never hides a testable one). Non-arc controls and pairs within one system are not merged.
On the stored tables it gives the #469 numbers unchanged (`posthoc_c2/posthoc.json` byte-identical). It was found on
Abell 370, so it is not validated there, and Abell 370's flags stay unread (owner's call). Alternatives rejected:
dropping shared controls (loses a real absence); weighting by 1 / n (no frozen threshold for weights). Revisit if: a
field's controls share pairs whose members disagree between absent and recovered, which would mean the per-system
predicted flux, not the pair, drives the class.

**Addendum 2026-10-23j (Abell 2744, UNCOVER v2.0 maps, gate K).** The next field after the five CANUCS-reference
fields (`results/lens_consistency/abell2744_uncoverref/`, frozen plan first). Gate K (`m0417_images`, mean back-traced
source, 0.25″ grid) on the 172 spectroscopic images of the UNCOVER v2.0 table fails: rms 4.90″, 94.2 % within 3″.
Decision: the UNCOVER v2.0 maps are not used for the extra-image stages under this estimator. Post hoc: three families
(34.1, 80.1, 74) get one predicted image for three catalogued ones, consistent with a mean source outside the caustic;
D-041's WLS source reproduced all 181 own-list images (rms 0.84″). The UNCOVER DR3 SUPER catalogue (83.7 MB),
segmentation map and bCG-subtracted F200W mosaic (1.84 GB gzip) are reachable from cloud sessions. Revisit if: a frozen
gate K with the WLS estimator (a deviation, chosen after this failure) or with Bergamini+23 passes.

**Addendum 2026-10-23k (Abell 2744, UNCOVER v2.0 maps, gate K-WLS).** Frozen as a deviation after addendum j's failure
(`abell2744_uncoverref/PLAN_1b.md`): the same gate on the same 172 spectroscopic images and grid, with only the family
source swapped to D-041's WLS source. Passes: rms 0.835″, 100 % within 3″; no family has fewer predicted than
catalogued images, 31 have one more. Decision: the UNCOVER v2.0 maps go on to the extra-image stages **with the WLS
source** (own-constraint validation only; the pass is labelled post-failure wherever quoted). Revisit if: stage 2's
gate C2 fails, or the WLS and mean sources disagree on an extra-image flag.

**Addendum 2026-10-23l (Abell 2744, UNCOVER DR3 SUPER reference: column convention).** Frozen before the run
(`abell2744_uncoverref/PLAN_2a.md`; rule `a2744_uncoverref.theta_check`). The README calls `theta_J2000` degrees
clockwise from "the first image axis (along y)"; its values span ±π/2. On 500 sources' segments of the release segmap, the SEP reading (radians,
counter-clockwise from +x) passes (median residual 2.7°, 98.4 % within 20°); clockwise from +y fails (46.6°, 23.2 %), and post hoc clockwise from +x fails (43.4°, 24.0 %).
The catalogue grid is north-up (+y PA 0.006°, +x PA −90.0°) at 0.040″. Decision: stage 2b maps sky PA = deg(θ) + 90
(mod 180), i.e. `m0417_canucsref.to_cat`'s offset 90 and sign +1 on θ in degrees; fluxes `f_f200w` / `e_f200w` × 10
nJy (README, `PHOT_ZP` 28.9); local noise from the GO 2561 pipeline F200W mosaic by byte-range cutouts (the
bCG-subtracted mosaic has no byte ranges; a larger σ near members only moves classes away from `absent`). Revisit
if: a later UNCOVER release changes the column, or the 2b noise route fails its own check.

**Addendum 2026-10-23m (Abell 2744 stage 2b: DR3 SUPER reference passes gate C2).** Frozen plan
`abell2744_uncoverref/PLAN_2b.md`. One GO 2561 F200W tile (`o001_t003`) covers every image, so the stage runs on it
alone. The local noise reads the mosaic's SCI extension once into memory (≈ 430 MB range stream, cheaper than ~170
cutouts with 64 kB read-ahead) at the mosaic's own 0.031″ pixels (new `pix` argument; radii stay in arcsec). Gate C2
with shared pairs passes (arc-referenced 11 / 13, all 30 / 34); no flag; 4.2 and 78.3 are `ambiguous` at 9.7 σ /
8.7 σ. Decision: the DR3 SUPER + pipeline-noise reference counts as validated on Abell 2744 for follow-ups under the
same rules; the near-threshold targets get a frozen second-band or bCG-subtracted check, never a re-read of the
10 σ threshold. Alternatives rejected: the bCG-subtracted mosaic (1.84 GB gzip, no byte ranges). Revisit if: the
second-band check of 4.2 / 78.3 disagrees with F200W, or C2 fails on another band.

**Addendum 2026-10-23n (Abell 2744 stage 2c and the 78.3c vetting).** Stage 2c (`PLAN_2c.md`): gate C2 passes in
F277W and F356W; frozen flags 4.2 (= vetted 4.2c) and 78.3c. Vetting (`PLAN_vet78.md`): CATS v4.1 has no critical
curve between the 78 a/b images (test 1 "tool not validated"); in the UNCOVER maps the c images of 78.1 / 78.2
predicted from a+b miss their catalogued c by 1.455″ / 1.293″ and lie 0.07″ / 0.17″ from the predicted 78.3c (test
1b "explained"). Decision: 78.3c is closed as a model c-position error; a predicted image whose own system's
detected images are mispredicted by more than the search radius is not a flag. Alternatives rejected: μ
perturbation of the northern group first (needs a re-fit; position test is cheaper and decisive). Revisit if: a
model with a fold critical curve through the 78 arc predicts knot 3's c image ≥ 1″ from the 78.1c / 78.2c images.

**Addendum 2026-10-23o (El Gordo, Caminha+2023, gate K; reference inputs).** Frozen plan
`results/lens_consistency/elgordo_djaref/PLAN_1.md`: gate K passes on the 56 spectroscopic images (rms 0.754″, 56 /
56 within 3″; own constraints only); 15 families have more predicted than catalogued images (stage-2 inputs, not
read). Rule 2 fails as written: the DJA v7.0 photometry has segment shapes but no total-flux column; post hoc it has
aperture fluxes with a per-band `tot_corr`. Decision: the route continues only with a stage-2 plan that freezes the
reference flux from those columns and gates it with C2. Revisit if: a JWST-era El Gordo catalogue with segment
totals (e.g. a PEARLS release) is public.

**Addendum 2026-10-23p (El Gordo stage 2: DJA flux definition and unit check).** Frozen plan
`results/lens_consistency/elgordo_djaref/PLAN_2.md`. The grizli 1.14.3 source shows the total flux is
`{band}_flux_aper_k × (flux_auto / flux_aper_k) × {band}_tot_corr`; `tot_corr` alone is only a point-source
Kron-to-total factor, so the stage-1 handoff's `flux_aper_k × tot_corr` is not adopted. Rule 2b (grid) passes:
rotation 0.0002°, east left, 0.0400″ per pixel. Rule 2a (unit) fails as written: median `mag_auto + 2.5 log10(flux_auto)`
is 23.69 mag, neither 23.9 nor 28.9 within 0.05, so the run stopped before any flux row, pixel or prediction was
read. Post hoc: with the detection `tot_corr` the zero point is 23.9 mag in every row (5–95 %: 23.9–23.9), i.e.
µJy and `mag_auto` is a total magnitude; the header's `ZP` = 28.9 does not describe these columns. Decision: stage 2b
freezes rule 2a as `mag_auto + 2.5 log10(flux_auto × tot_corr)` (disclosed as already seen) and otherwise runs
PLAN_2.md unchanged. Revisit if: as addendum o.

**Addendum 2026-10-23q (El Gordo stage 2b: gate C2 passes, no flag).** Frozen plan
`results/lens_consistency/elgordo_djaref/PLAN_2b.md` (rule 2a with `tot_corr`, disclosed as already seen). Rule 2a
gives 23.9 mag (µJy); gate C2 (shared pairs) passes, 14 / 15 testable controls, arc-referenced 7 / 8; no flag.
Of 16 predicted extra images only 1 is testable (target 23, `missing`), and post hoc it is the model's 23c,
1.538″ from the catalogued 23c (`unpredicted` in the frame table): a model position error, as 78.3c in addendum n.
Contact sheet looked at (`sheet_2b.png`, `sheet_2b_t23.png`). Decision: El Gordo's reference-route test is closed
with no flag and low testability (1 / 16); the route moves to the next field. Revisit if: a deeper El Gordo
reference catalogue or a model with a better-fitting 23c is public.

**Addendum 2026-10-23r (SMACS 0723, Mahler+22 ICLv2, gate K; reference inputs).** Frozen plan
`results/lens_consistency/smacs_djaref/PLAN_1.md`. Rule 0 (ASSUMPTION): systems without `z_m_limit` in `input.par`
(1, 2, 3, 5, 7; 16 images) stand in for the spectroscopic set. Gate K passes (rms 0.321″, 16 / 16 within 3″; own
constraints only). The DJA v7.4 SMACS `fix_phot` header carries every column El Gordo stage 2b read (headers only).
9 families have more predicted than catalogued images (stage-2 inputs, not read). Decision: SMACS goes to stage 2
with El Gordo stage 2b's rules under a new frozen plan. Revisit if: a cross-checked spec-z list (Mahler+22 /
Noirot+23) differs from the rule-0 set.

**Addendum 2026-10-23s (SMACS 0723 stage 2: gate C2 fails).** Frozen plan
`results/lens_consistency/smacs_djaref/PLAN_2.md` (El Gordo stage 2b's rules on GO 2736 o001_t001 F200W and DJA
v7.4; deviations before any read: all 60 images with fitted z for free systems, z-edge family 8 not flaggable,
fixed / fitted z labels). Rules 2a / 2b pass (µJy; 0.028°). Gate C2 fails: 16 / 23 (0.70), arc-referenced 7 / 12
(0.58), against 0.8; it also fails without the containing-source or local-noise rule. No flag is read; the 11
targets would be untestable anyway. Post hoc contact sheet (`sheet_2.png`): sources sit inside the search circle
at 5 of 6 `absent` controls, all 6 `absent_sources_too_faint`, i.e. flux-ratio misses of the route, cause not
isolated. Decision: the reference route is closed for SMACS 0723 (model ICLv2, DJA v7.4); WHL0137 next only
with a public image list (D-033 addendum 4 already covers its systems 1–2).
Revisit if: a SMACS model with spectroscopic z for systems 12–25 or a total-flux reference built for arcs is public.

**Addendum 2026-10-23t (RXJ2129, Dhanasingham+2026 Lenstool, gate K; reference inputs).** WHL0137 is skipped: its
JWST-era image list (Scofield+2025) covers only systems 1–2, already tested in D-033 addendum 4, and its JWST-era model
has an undocumented deflection unit. Frozen plan `results/lens_consistency/rxj2129_djaref/PLAN_1.md`: the paper's
Appendix B table (31 images, 10 families, all spectroscopic, JWST/Gaia frame) against the pinned
`rxj2129-lenstool-maps`. Gate K passes (rms 0.449″, 31 / 31 within 3″; own constraints only). The DJA v7.0
`fix_phot` header carries every SMACS rule-2 column. 6 families have more predicted than catalogued images
(stage-2 inputs, not read). Decision: RXJ2129 goes to stage 2 with El Gordo stage 2b's rules under a new frozen
plan. Revisit if: a published erratum or data release changes the table (the text says 11 systems, the table has 10
families).

**Addendum 2026-10-23u (RXJ2129 stage 2: gate C2 fails).** Frozen plan
`results/lens_consistency/rxj2129_djaref/PLAN_2.md` (El Gordo stage 2b's rules on GO 2767 o002_t001 F200W and DJA
v7.0; deviation D1 before any read: overlapping targets of the one-galaxy group 3 / S3-1 / S3-2 are joined for flags,
`merge_knot_targets`, unit-tested). Rules 2a / 2b pass (µJy; rotation −2 × 10⁻⁵°). Gate C2 fails: 4 / 7 (0.57),
arc-referenced 2 / 5 (0.40), against 0.8 and minimum counts 8 / 3; it also fails without the containing-source or
local-noise rule. No flag is read. Post hoc: all 3 `absent` controls (1, 2, 7) take the same reference row, DJA
1506 (16.38 mag, 3a = 8.7″), so their predicted fluxes are 16.7–18.1 mag; the contact sheet (`sheet_2.png`) shows
faint sources in the circles of the system 1 and 2 controls. Without that row C2 would be 4 / 4 (arc 2 / 2; computed in the session, not stored), still
below the minimum counts, so the route cannot validate itself on this field. Decision: the reference route is
closed for RXJ2129 (Lenstool maps, DJA v7.0). Revisit if: a reference catalogue with arcs deblended from the bright
cluster light (e.g. a bCG-subtracted RXJ2129 catalogue) is public.

**Addendum 2026-10-23v (MACS0647 input check).** `results/lens_consistency/macs0647_inputs/RESULT.md`. Image lists:
Meena+2023 (72 images, 23 systems, photometric z only) and Nishida+2025 (86 images, 28 sources, the JWST GLAFIC
constraints); where they share images they agree to about 0.1″ (66 / 86 Nishida images within 0.3″). Reference: DJA v7.0 passes SMACS rule 2. Model maps:
no JWST-era model (Meena+2023 dPIEeNFW / Glafic, Nishida+2025 GLAFIC) has public maps in its arXiv source or on
Zenodo; only the HST-era CLASH Zitrin v2 maps are pinned. Decision: the fitted-model extra-image route does not run
on MACS0647. The open alternative is a held-out test of the CLASH maps on the JWST-first systems (Meena+2023 13–23),
only under a frozen plan with a redshift rule, a CLASH-to-JWST frame check and gate thresholds stated for a model
not fitted to those images. Alternatives rejected: transcribing a model from figures (no maps, no parameters);
running gate K with the CLASH maps under the fitted-model thresholds (most images were not its constraints).
Revisit if: Nishida+2025 or Meena+2023 maps (or GLAFIC parameter files) become public.

**Addendum 2026-10-23w (MACS0647 CLASH held-out test).** `results/lens_consistency/macs0647_clash_heldout/RESULT.md`,
frozen `PLAN.md`. Families at fitted z: CLASH LTM-Gauss passes gate K0 on the pre-JWST families (rms 1.66″, 37 / 38
within 3″) and frame check F, and reproduces 31 / 31 certain images of the JWST-first systems 13–23 (rms 1.41″); the
PIEMD + eNFW map fails K0 (2.34″ > 2.0″). No flag. Post hoc: at the family photo-z both models do worse on the
pre-JWST systems than on the held-out ones (LTM-Gauss 20 / 35 vs 25 / 31), so the photo-z are no redshift
constraint. Decision: MACS0647 gives no lensing-residual input from the CLASH maps; the fitted-z test only detects
images that miss by > 3″ at every z. Alternatives rejected: photo-z as the frozen redshift rule (inconsistent within
families, e.g. 6-3: 0.52 and 9.19). Revisit if: spectroscopic redshifts of systems 13–23, or a JWST-era model's maps.

**Addendum 2026-10-24x (Bullet Cluster input check).** `results/lens_consistency/bullet_inputs/RESULT.md`. The CANUCS
Bullet DR1 on CANFAR gives the Rihtaršič+2026 Lenstool gold image list (135 images, 49 systems, 27 galaxies with system
z), the fiducial model's best-fit deflection maps (50 mas, D_LS / D_S = 1) and Bayesian sample maps (200 mas), and a
CANUCS-format photometry catalogue that passes the frozen column mapping. Decision: the Bullet Cluster is the next
field for the fitted-model extra-image route (gate K, then gate C2), in the MACS0416 / MACS1149 CANUCS form.
Alternatives rejected: the paper's Dropbox folder link (one 24.5 GB zip, ranges ignored). Revisit if: gate K fails,
or the 50 mas maps cannot be fetched within the cloud disk rule (then the 200 mas sample maps).

**Addendum 2026-10-24y (Bullet Cluster gate K).** `results/lens_consistency/bullet_gatek/RESULT.md`. Gate K of the
CANUCS Bullet DR1 Lenstool best-fit maps (50 mas) on the 135 gold images fails as written: rms 1.196″ (≤ 1.0″), 97.0 %
within 3″. Two caustic-straddling sources (K1b, K1e) get one predicted image instead of three under the mean source;
post hoc, the other 129 images give rms 0.435″ against Lenstool's 0.419″, and a 0.1″ grid changes nothing. Gate K-WLS
(`PLAN_1b.md`, frozen as a post-failure deviation, the addendum k precedent unchanged) passes: rms 0.391″, 100 %
within 3″. Decision: the Bullet maps go on to gate C2 and the extra-image test **with the WLS source**; no flag read
yet. Alternatives rejected: exempting K1b / K1e by hand. Revisit if: gate C2 fails, or the WLS and mean sources
disagree on an extra-image flag.

**Addendum 2026-10-24z (Bullet Cluster extra-image test).** `results/lens_consistency/bullet_extra/RESULT.md`. The
CANUCS F200W `bgsub-sci` mosaic header carries `PC` without `CDELT` and a 0.08″ `CD`; the frozen deviation
`PLAN_1b.md` uses the header's `CRPIX` / `CRVAL` / `PC` with `CDELT` 0.04″, gated by the catalogue (W1, rms 1e-10″,
0-based `X`, `Y`) and by point-source centroids (W2, median offsets −0.035 / +0.205 px, 51 sources). Frame check F
passes (92 / 135, medians < 0.001″); gate C2 (shared) passes, 58 / 67 and arc-referenced 37 / 42; **no flag**
(0 / 22 targets absent; 13 untestable at local 5 σ). Decision: the Bullet Cluster extra-image test is closed as a null
with these limits; the route (CANUCS reference + local noise + containing-source rule) is validated here as on
MACS1423 / Abell 370. Alternatives rejected: a free TAN fit to catalogue `X`, `Y` (not needed once W1 agrees to
1e-10″). Revisit if: CANUCS republishes the mosaic with a different WCS, or a deeper / star-cleaned F200W image makes
N2a, N5a, N5b testable.

**Addendum 2026-10-24aa (MACS J1931 input check).** `results/lens_consistency/macs1931_inputs/RESULT.md`. The VENUS
AstroLensPy v2.1 model (Allingham+2026, Zenodo 10.5281/zenodo.18646525) gives a JWST-era image list (51 used images,
19 systems, 9 spec-z systems with 26 images) and 0.04″ deflection maps; its headers say `BUNIT = 'Pixels'` while the
README says arcsec, and the α normalisation redshift is not in the header. DJA v7 lists no MACS J1931 field.
Decision: MACS J1931 is a gate-K-only field: next a frozen plan that resolves the unit (the authors' description, or
both readings with gate K on the spec-z images deciding) before any residual is read; no gate C2 / extra-image test
until a deblended reference catalogue is public. Alternatives rejected (survey, not stored): PLCK G165.7+67.0,
MACS J0138, SPT-CL J0615, Sunburst (no public JWST-era maps found). Revisit if: a VENUS or DJA catalogue for the field
appears, or MAST hosts a VENUS HLSP with a catalogue.

**Addendum 2026-10-24ab (MACS J1931 gate K).** `results/lens_consistency/macs1931_gatek/RESULT.md`. The VENUS
`Alpha_x/y` maps are in pixels (0.04″; header `BUNIT`, not the README's arcsec) and normalised at z_s 5.079 (the
authors' `k_dlds` reproduced to 4e-6): unit check U accepts only the pixel reading (median 0.110″ to the authors'
sources) and gate K passes on the 26 spec-z images (rms 0.919″, 100 % within 3″). Decision: use these maps, with
that scaling, for the 9 spec-z systems only; the photo-z systems 8 and 10 are not reproduced as catalogued and the
mixed-z systems are not tested coherently by our single-z solver (reported 4.83″ over 44 images); they stay out until
their z is settled. Alternatives rejected: the
README unit (fails U by three orders of magnitude). Revisit if: a later VENUS release changes the maps or the z list.

**Addendum 2026-10-24ac (MACS J1931 imaging footprint).** `results/lens_consistency/macs1931_counter/RESULT.md`. VENUS GO 6882 tile
`jw06882-o038_t041` covers all 51 used images in 10 NIRCam bands; the 1.29 GB F200W i2d reads by S3 byte range
(0.0312″ / px). Decision: imaging is not the blocker; the reference catalogue is. Proposed route (not frozen, `PLAN.md`):
an own photutils deblended catalogue on a background-filtered F200W cutout, accepted for MACS J1931 only if it passes
the frozen MACS0417 gate C2 (where CANUCS passes and pipeline segments fail). Alternatives rejected: pipeline segments
(addendum 2026-10-22); waiting for a VENUS HLSP catalogue (unknown; archive.stsci.edu unreachable from cloud).
Revisit if: a public deblended MACS J1931 catalogue appears.

**Addendum 2026-10-24ad (MACS J1931 stage V: own deblended reference).** `results/lens_consistency/macs0417_ownref/RESULT.md`.
Reuse: photutils (`Background2D`, `detect_sources`, `deblend_sources`, `SourceCatalog`) builds the reference; it is
the library the JWST pipeline's catalogue step uses, pip-installable (new `deblend` extra; also in `dev` for tests).
Alternatives rejected: SEP (duplicates photutils, D-010); SourceXtractor++ (conda-only and heavy, D-010; this
need is met without it); pipeline segments (addendum 2026-10-22). On MACS0417 the frozen V1 run had no defined depth
(0 rows at Kron S/N 4–6), so its gate is not read; the V2 amendment (5σ median Kron error depth, containing rule,
local noise) passes gate C2 with shared pairs at the threshold, 28 / 35, arc-referenced 14 / 14 (CANUCS 37 / 40,
pipeline 10 / 27). Decision: the builder may be used for the MACS J1931 stage 3, with zero margin stated next to
any flag. Revisit if: a public deblended MACS J1931 catalogue appears, or a second validation field fails C2 with
this builder.

**Addendum 2026-10-24ae (MACS J1931 stage 3 stops at frame check F).** `results/lens_consistency/macs1931_stage3/RESULT.md`.
Gate K-WLS repeats (0.630″, 26 / 26) but frame check F fails as frozen: 18 / 51 images have an own-catalogue row within
0.5″ (rule ≥ 20; medians 0.045″ / 0.004″), and again under deviation 3b with every image inside the section (18 / 51). Post hoc only 15 / 51 positions lie inside any deblended segment; the
contact sheet (not stored) shows faint or unsegmented images (systems 5, 6, 7, 100, 101; arcs 8.2, 9.1, 10.1, 12.2). Decision: no
MACS J1931 counter-image reading from this builder; any retry is a new plan frozen before its run (e.g. a deeper or
multi-band detection image, or arc-aperture sums), validated on MACS0417 first, with the rule 3 depth defined (the pipeline S/N 4–6 proxy is NaN here). Alternatives rejected: lowering F's
minimum after seeing 18 (tuning on the outcome). Revisit if: a public deblended VENUS catalogue appears.

**Addendum 2026-10-24af (MACS J1931 stage 3c: wider background box; F fails again).** `results/lens_consistency/macs1931_stage3c/RESULT.md`, plan frozen in e1f21a8. The one builder change (background box 64 → 256 px, ≈ 8″) passes the MACS0417 gate (C2 shared 24 / 27, arc 10 / 10), but MACS J1931 frame check F fails: 17 / 51 within 0.5″ (rule ≥ 20), medians 0.033″ / −0.021″; the new chance control passes (mean shifted count 3.6). Contact sheet (not stored): about 20 of 51 positions show no visible F200W emission (by eye), and the wider box leaves cluster light as large segments. Decision: stop the single-band F200W own-reference route for MACS J1931; no counter-image reading. Alternatives rejected: further builder tuning against this F (each variant gave 17–18); lowering F's minimum. Revisit if: a frozen F restricted to images with a pre-declared continuum criterion (e.g. a pipeline row in any band) is validated on MACS0417, or a public deblended VENUS catalogue appears.

**Addendum 2026-10-24ag (MACS J1931 stage 3d: continuum-restricted frame check).** `results/lens_consistency/macs1931_stage3d/RESULT.md`, plan frozen in 1314a97, functions and tests in bce45c1. F_c counts only images with ≥ 5 σ F200W light in a 0.3″ aperture (local random-aperture noise). It passes the MACS0417 validation (41 / 46) and MACS J1931 (17 / 18 of 18 continuum images; 33 / 51 faint). Gate C2 then passes (11 / 11, arc 7 / 7), but 0 of 11 predicted extra images is testable: 10 have no flux reference (systems 5, 6, 7, 100, 101, whose 13 observed images have S/N −1.2 to 1.4) and 1 is predicted at 1.11 σ. Decision: the F200W own-reference route for MACS J1931 is finished with no flag and no null. F_c is the frame check for own references on fields with continuum-free images. Alternatives rejected: lowering stage 3's F minimum (not validated); a deeper continuum stack (the five systems show no F200W light at any observed image). Revisit if: a public MUSE cube or line map gives line-flux references for systems 5, 6, 7, 100 and 101, under a frozen plan.

## D-035 Six HFF clusters as CATS map models; image lists gated on their image-plane rms (2026-10-08)

**Decision.**
- **`lensmodel.HFF_CATS`** pins the CATS Lenstool v4/v4.1 HLSP products by sha256 for MACS0416, MACS1149, Abell 370,
  MACS0717, Abell S1063 and Abell 2744: x/y deflection (arcsec, D_LS/D_S = 1), κ, z = 2 magnification and
  `arcs.txt`, plus `params.txt` where it exists (MACS0416, MACS1149, MACS0717; it returns 404 for the others).
- **`lens_consistency.py`** has `MODELS` `<cluster>-cats` (`kind: maps`, through `MapLensModel`, D-033).
  - `validate` compares κ and μ(z = 2) with the published maps and measures the image-plane rms of `arcs.txt`, with
    system redshifts from `params.txt` (`read_z_m_limit`, which reads only `z_m_limit`).
  - An image list is used by `images` / `fluxratio` only where that rms is ≤ 1.5× the release's quoted rms
    (ASSUMPTION). `validate` reports this gate (`image_list_gate`) and whether it agrees with the `MODELS` setting.
  - For map models, χ² and the "> 3σ" list are null: there is no published position error.
- **Parser:** a parenthesised redshift in an image list, e.g. "(2.16)", is read as that value. `parse_lenstool_par` and
  `read_z_m_limit` share one `z_m_limit` line parser and key every id by `image_family`.

**Alternatives rejected.**
- Using every `arcs.txt` as published. Without the model's fitted redshifts, placeholder redshifts give 9–12″
  residuals (Abell 370, Abell S1063, Abell 2744), which would fake "absent" counter-images.

**Evidence** (`validate`, 2026-10-08; `model_prediction` against the published products).

| Cluster (z_lens) | κ median rel. diff | μ(z = 2) median rel. diff | image-plane rms (ours / quoted) | image list |
|---|---|---|---|---|
| MACS1149 (0.543) | 3.3e-3 | 5.6e-4 | 0.67″ / 0.63″ (145 images) | used |
| MACS0717 (0.545) | 3.2e-3 | 1.2e-3 | 3.21″ / 2.41″ (132) | used |
| MACS0416 (0.396) | 2.5e-3 | 3.7e-4 | 1.57″ / 0.72″ (116; system 26 off by 11″) | map only |
| Abell S1063 (0.348) | 1.8e-2 | 5.9e-3 | 11.8″ / 0.48″ (72) | map only |
| Abell 370 (0.375) | 3.0e-2 | 2.1e-3 | 10.7″ / not quoted (101) | map only |
| Abell 2744 (0.308) | 5.1e-4 | 8.8e-5 | 8.9″ / not quoted (87 of 109 solved) | map only |

- The lens redshifts reproduce every z = 2 map (median ratios 0.994–1.000), which confirms them.
- The κ differences of 0.3–3 % come from finite differences on 0.2–0.8″ pixels.

**Revisit if.**
- `params.txt` (or another source of fitted redshifts) becomes available for Abell 370, Abell S1063 and Abell 2744.
- MACS0416 system 26 is understood.

**Addendum 2026-10-21 (Abell S1063 with Chan+2020 redshifts, `results/lens_consistency/abells1063_d035/`).** The Chan
list with MUSE z_spec gives rms 18.5″ on the 12 z_spec systems (gate 0.72″), and no system reaches 0.5″. Post hoc:
the CATS `arcs.txt` z already equal those z_spec within 0.005 for the 33 one-to-one matched images (11 systems), so
for these images the residuals are not caused by missing redshifts (the other 39 CATS images were not compared). Five systems get a single predicted image from the mean traced
source, the other seven 0.80–1.83″. Cause not separated (mean-source estimator, map resolution, normalization).
S1063 stays map only. Revisit if a route with an image-plane source fit reproduces MACS1149 and passes on S1063's
z_spec systems.

**Addendum 2026-10-21b (image-plane source fit, `results/lens_consistency/imageplane_fit/`).** Fitting each system's
source in the image plane (one-to-one pairing) passes its frozen checks (synthetic 10 / 10; MACS1149 145 / 145 at
0.657″ vs quoted 0.63″). On S1063's z_spec systems it pairs all 36 images (the 18.5″ was the mean-source route) but
the rms is 0.978″ (> 0.72″) and no system reaches 0.5″ (0.633–1.460″): S1063 stays map only; the map, not the route,
is now the limit (cause not separated). Revisit if a public S1063 model with its own parameter file or image list
appears.
## D-036 Abell 2744 lens-model and exotic screens: null (2026-10-08)

**Decision.** The Abell 2744 screens (Bergamini+2023b) are a null result. Every flag has an ordinary explanation
(docs/fields/abell2744.md, 2026-10-08). The worker's single survivor of the cheap tests, 4.2c, was re-examined
under the D-034 rules: its reference images are resolved knots, and the prediction lies on the BCG halo. It is
untestable, not missing.

**Alternatives rejected.** Reporting the forced-photometry `absent` / `confused` classes directly. Near-critical
magnifications, the fitted-redshift system 700.1, resolved-arc references and the 1″ search made them misleading.
D-034's rules now remove these automatically.

**Evidence.**
- `validate`: χ² 146.64 against 146.60.
- `images`, after D-034: absent only for 700.1, which the model reproduces at none of the sampled z = 1, 1.5, 2, 3, 5 (finer scan: addendum 2026-10-20 at the end of this entry).
- `fluxratio`: 30 consistent, 0 under- or overluminous, 3 chromatic (a 602.1 blend).
- `radial`: 35 peaks against a random mean of 50.4; p ≥ 0.575.
- Cutouts were inspected for every flag.
- Tally: 176 predicted images, 149 flux-ratio images and 35 radial peaks screened; 16 flags; **0 surviving**.
- Family 4's c images (4.1c underluminous 4–8× after BCG subtraction; 4.2c undetected) were vetted
  (docs/candidates/abell2744-family4-c.md). They are explained by μ(4.1c) systematics next to member 34423:
  ±30 % changes give 3.9–28.7, and CATS v4.1 gives 7.3.
- **Rule:** an under- or overluminous image whose μ moves by more than 2× under ±30 % changes of the nearest member
  potential, or under an independent model, is untestable.

**Revisit if.**
- `bayes.dat` position spreads do not cover 34.1 / 700.1.
- An independent model (UNCOVER v2.0, or the CATS v4.1 maps of D-035) predicts a bright image where none is seen.
- **Addendum 2026-10-20 (700.1 fine z scan, `results/lens_consistency/a2744_s700_zscan/`):** frozen outcome B, "not
  reproduced at any z". Bergamini+23 gives 700.1a and 700.1b separate images at none of z 1.000–4.000 (step 0.025);
  at z 1.075–1.275 they share one image (μ ≈ −3) and 700.1c is reproduced; above z 1.3 the back-traced sources spread.
  The same rule reproduces 39 / 40 fixed-z systems (failure: the known 34.1 fold pair). So the residual is a missing
  critical curve between a and b in the infall region, not the fitted z. In the stamps a/b are the two ends of one
  thin 1.7″ feature (two images, or one image catalogued as two: untested). The system stays model-limited and is
  not a candidate. Revisit if a model with constraints near (3.579, −30.358) appears, or if UNCOVER v2.0 is scanned the
  same way.
- **Addendum 2026-10-20b (UNCOVER v2.0, `results/lens_consistency/a2744_s700_uncover/`):** frozen outcome "tool not
  validated" (the one-to-one rule reproduces 44 / 60 UNCOVER own systems at 1.5″, gate 0.8), so no z scan. Post hoc:
  700.1a/b/c are UNCOVER 78.21–78.23 (0.02″), z 2.3151 spectroscopic per the UNCOVER table ("Price+24"). At that z
  both models merge a and b into one predicted image, UNCOVER even with system 78 as a constraint. The nested knot
  pairs (midpoints within 0.05″) fit a fold arc or one symmetric galaxy. 700.1 is undetermined and not a candidate;
  next is `/vet-candidate` (perturber near the midpoint; one image or two). Revisit if that vetting or a model with a
  local perturber there changes the picture.
- **Addendum 2026-10-20c (vetting, `results/lens_consistency/a2744_s700_vet/`):** frozen outcome "undetermined
  (colour test not validated; member perturber present)". Pair colours fail on known UNCOVER pairs (10 / 20, gate
  0.8). Spec-z members at 2.11″ (DJA 45749 = Bergamini+23 member 20036) and 4.43″ (DJA 45736 = group halo O4/O9
  centre) are both in the model (post hoc). 700.1 is parked as model-limited near the northern group core, not a
  candidate. Revisit if a local re-fit with O4/O9/20036 free, or spectra of both arc ends, becomes available.

## D-037 MACS1149 screens: null; repeated-pair, fitted-redshift and model-dependent-μ rules (2026-10-08)

**Decision.** The MACS1149 screens (`macs1149-cats`, CANUCS program 1208 `jw01208-o008_t004`) are a null result
(docs/fields/macs1149.md). There were 3 forced-photometry flags. Two are ordinary model errors: system 16 (fitted
redshift / pair topology) and system 2 (topology). One is untestable: knot 1192 (model-dependent μ). New rules
(ASSUMPTION: thresholds):
- **Repeated pair.** An extra predicted image is untestable (`model_topology`) when its system's catalogued images
  share one predicted image (`validate`'s `shared_match_images`) and the image lies within 8″ of that pair. MACS1149
  system 16's flagged image is 7.2″ from 16.1 and 6.9″ from 16.2.
- **Fitted redshift.** The same applies when the system's redshift is model-fitted (`z_m_limit`, or no spectroscopic
  z in `arcs.txt`) and a redshift change within ±50 % moves the image onto a catalogued one. For system 16, CATS needs
  z 4.419 → 2.5 (−43 %) and Sharon v4cor needs → 3.0 (−32 %).
- **Model-dependent μ.** D-036's μ-stability rule is extended from flux-ratio flags to forced-photometry flags: a
  `confused` or `absent` image whose μ differs by more than 2× between independent models is untestable.
- **BCG annuli.** Where a forced-photometry background annulus crosses a BCG core, the flux goes negative. Use
  high-pass photometry there (0.6″ median filter).

**Alternatives rejected.** Treating system 16's third image as missing: the predicted image exists only at the
CATS-fitted z = 4.419.

**Evidence.**
- `validate`: κ median relative difference 3.3e-3; μ(z=2) 5.6e-4; image-plane rms 0.673″ over 145 images (the D-035
  gate passes). The frame offset is under 0.02″, so `MODELS["macs1149-cats"]` has no `frame_offset_arcsec`.
- `images`: 159 predicted, 16 forced-tested (F277W). `fluxratio`: 19 consistent, 106 resolved, 0 under- or
  overluminous. `radial`: 98 anti arcs; 32 peaks against a null mean of 30.8 (p95 41); p = 0.91. With CANUCS DR1 photo-z (134
  non-background sources dropped, including some blended lensed images): 62 anti arcs, 12 peaks against 12.8 (p95
  20), max 3 lines; null.
- System 16 (μ 15.8; empty sky at −0.4σ after high-pass): at z = 2.5 (CATS) or 3.0 (Sharon v4cor), both models move
  the image onto 16.2, so 16.1 and 16.2 become a merging pair. CANUCS DR1 photo-z for 16.2: z 2.25 (95 % 0.23–2.33; upper bound
  2.67 with the 0.1 × (1 + z) margin). This disfavours the fitted 4.419. In CATS the pair is close to merging at z
  2.0–2.5. (derived) The far image (μ ≈ 4, `recovered` at 0.95″ but 0σ at
  the fixed position) moves 3–5.5″ between z = 2.5 and 3.0, so it is untestable too.
- System 2: CATS merges 2.2 and 2.3; Sharon places the critical curve between them.
- Knot 1192 (SN Refsdal host) is 1.1″ from the BCG. CATS gives μ 7.5; Sharon gives 2.2 at its image 0.56″ away.
- Tally: 159 predicted images, 145 flux-ratio images and 98 anti arcs screened; 3 flags; **0 surviving**.

**Revisit if.**
- A spectroscopic z for system 16 is published (search the MUSE catalogues of Grillo+2016 and Treu+2016), or a source
  with 16.1's colour is found on the far-image track (z = 2.5–3.5 positions). Either would fix z.
- `images` is re-run with CANUCS DR1 photo-z (counterpart redshifts). DJA v7 has no MACS1149 mosaic.

**Addendum 2026-10-14 (the 11 `no_reference` rows).**
- **Decision.** The CATS forced rows without a reference flux are tested with sibling references (no |μ| cap,
  three bands) and Sharon v4cor solved per system (`scripts/m1149_refs.py`). They add no flag: 0 `absent`.
- **Evidence.** 13 is CATS's misplaced 13.2 (Sharon's image of it lies 1.5″ away; CATS does not reproduce 13.2);
  1002 and 1222 are CATS-only extras; 1007 is predicted by both, with μ 18.3 vs 0.1 (untestable); 14 is
  undetectable (2.4σ); six both-model rows lie in the resolved host images 1.2/1.3 (ungated `recovered`, `confused`
  or `undetectable`). Six rows stay untestable; no injection, so no limit. docs/fields/macs1149.md.
- **Revisit if** a JWST-era MACS1149 model with knot-level constraints appears.

**Addendum 2026-10-14 (system 16's redshift, `scripts/m1149_sys16.py`).**
- **Decision.** System 16's redshift is taken as z ≈ 2.5–2.9: the range where both CATS v4.1 (2.10–2.90) and Sharon
  v4cor (2.50–3.30) reproduce 16.1/16.2 within 1.0″ (ASSUMPTION). The predicted-image photometry did not constrain z
  (no extra image predicted ≥ 10σ). CANUCS 5106117 (z_phot 2.94, 2.49–3.16), 0.4–0.5″ from both models' far image at
  z 2.45–2.65 with pair-like colours, is a possible third image, not a confirmed one.
- **Evidence.** z scan 2.0–6.0; CATS's fitted 4.419 gives pair rms 4.3″ (3.9″ in Sharon); CANUCS photo-z of 16.2
  overlaps at 2.5–2.67. The other colour-matched track sources have photo-z 1.23, 1.68 and 4.59 (chance).
  docs/fields/macs1149.md.
- **Revisit if** a spectroscopic redshift for 16.1/16.2 or CANUCS 5106117 appears.
- **Addendum 2026-10-20 (spec-z check, `results/lens_consistency/macs1149_s16_specz/`):** frozen outcome *not
  decided*. No spectroscopic z for 16.1 / 16.2 in Treu et al. 2016 (their 27.1 / 27.2, compiling GLASS and MUSE
  redshifts) or in the DJA NIRSpec v4.4 table (16.2 has an ungraded CANUCS PRISM spectrum, `z_best` −1; 16.1 and
  5106117 have none). Post hoc: DJA's EAZY photo-z of 16.2 is 2.79, inside both models' merging
  ranges (CANUCS DR1's 2.25 is inside the CATS range only; exploratory). Revisit condition unchanged; a DJA grade of `1208_5114276` would meet it.

**Addendum 2026-10-20 (Sharon v4cor pinned, `scripts/m1149_sharon.py`).**
- **Decision.** Sharon v4cor is the `MODELS` entry `macs1149-sharon` (sha256-pinned maps, z_lens 0.543). Its frozen
  known-case gate G2 failed as written, because system 9 lies off the 200″ maps; it is used only inside the map
  footprint (post hoc condition).
- **Evidence.** G1: κ and z = 2 |μ| reproduce the published maps (median relative difference 0.00025 / 0.016). G2:
  readme Table 5 gold images, rms 0.226″, 7 / 7 in-footprint systems. CATS on the same images: system 13 needs its
  fitted z 1.301, not the spec-z 1.24. `results/lens_consistency/macs1149_sharon/RESULT.md`.
- **Revisit if** a wider Sharon map release appears, or a D-037 row lies outside the footprint.

**Addendum 2026-10-21 (D-041 rule with the pinned Sharon, `scripts/m1149_crossmodel.py`).**
- **Decision.** All 16 non-`observed` CATS rows are read with the D-041 cross-model rule against `macs1149-sharon`;
  D-037's tally (0 surviving) stands. The one frozen flag (system 16's `absent` image at the CATS-fitted z 4.419,
  Sharon extra image 0.90″ away) is D-037's fitted-redshift case: at 4.419 Sharon reproduces only 16.2, and at the
  adopted z 2.7 Sharon's nearest extra image is 20.2″ away (model-dependent).
- **Evidence.** Known-case gate 38 / 39 catalogued images (rule ≥ 0.8, ASSUMPTION); all rows inside the footprint;
  the 11 `no_reference` readings equal the scratch-Sharon ones. `results/lens_consistency/macs1149_crossmodel/RESULT.md`.
- **Revisit if** a spectroscopic z for system 16 is published (as above), or a JWST-era MACS1149 model appears.

## D-038 MACS0416: system 26 is a solver-grid artefact at a fold caustic; radial screen null (2026-10-08)

**Decision.** MACS0416 system 26's 11″ image-plane residual (D-035) comes from `find_images`. With the default 0.25″
solver grid, it misses the merging 26.1/26.2 pair, because the mean source lies 0.001–0.005″ from the fold caustic
(|μ| 142, 153, 69). It is not a redshift problem and not a model failure. `macs0416-cats` stays map-only (image list
gated off) until the solver refines its grid near high |μ|. Pairs whose catalogued images share one predicted image
(`shared_match_images`) are checked at a 0.1″ grid before they are called model topology. The radial screen is null
(docs/fields/macs0416.md).

**Alternatives rejected.**
- Dropping system 26 by hand (rms 0.809″): this hides a solver limitation.
- Opening the image list on the 0.1″ result alone: system 122 still shares one predicted image on both grids, and the
  gate should not depend on a hand-chosen step.

**Evidence.**
- Data: CANUCS 1208 `jw01208-o004_t002`, 8 bands, jwst 3.0.0.
- `validate`: κ 2.5e-3; μ(z=2) 3.7e-4. Image-plane rms:
  - 1.572″ on the 0.25″ grid (116 images; 2.18× the quoted 0.72″);
  - 0.811″ on a 0.1″ grid (1.13×), where system 26's residuals are 0.15″, 1.01″ and 1.12″.
- Redshift: `params.txt` has no `z_m_limit`. z = 3.238 is close to the source-plane rms minimum (0.071″ at z 3.30).
  All three images have F200W counterparts within 0.22″.
- Frame offset (arcs.txt → JWST, 54 images): dRA +0.208″, dDec −0.025″.
- MACS1149 control: its six shared pairs are identical at 0.25″ and 0.1″, so D-037's topology verdicts stand.
- `radial` (no photo-z): 139 anti arcs, 98 centres against a random mean of 99.6 (p95 120). The strongest centre has
  7 lines, p = 0.29. Its lines are noise segments: aper50 S/N 1.6–4.0, nothing in cutouts, forced S/N within ±1.5σ.
  With aper50 S/N ≥ 5, there are 7 centres against 6.8 and at most 3 lines. Re-run with CANUCS DR1 photo-z (77
  non-background sources dropped, some of them blended lensed images): 120 anti arcs, 77 centres against 85.0 (p95 103), max 7 lines, p 0.225; null.
  In the JWST frame (offset applied; D-040 fix): 101 against 100.2 without photo-z, 81 against 84.9 with; null.
- Tally: 137 anti arcs screened; 0 flags; **0 surviving**. `images` and `fluxratio` were not run (gated).

**Revisit if.**
- `find_images` gets adaptive refinement near |μ| > 50. Then open the MACS0416 gate and run `images` / `fluxratio`
  with `frame_offset_arcsec` (0.208, −0.025).
- `images` / `fluxratio` run with CANUCS DR1 photo-z once the image list opens.

## D-039 Persistence test on level-2 exposures; Sunrise `n0022` and `n0150` are afterimages (2026-10-08)

**Decision.** Before a single-epoch source counts as a transient, `scripts/persistence_check.py` measures it in every
level-2 `_cal` exposure that covers it, and measures the same *detector pixel* in the earlier exposures on that
detector (≤ 3 h). A detection is `suspect` when an earlier exposure put ≥ 20× its flux, or a saturated pixel, there;
earlier exposures that put the position itself on that pixel (< 2 px) are skipped. A position with detections
but none clean (clean = not suspect, with ≥ 1 earlier exposure checked) is `persistence`; two clean detections make
it `on_sky`. Thresholds are ASSUMPTIONs.

**Alternatives rejected.**
- Re-running calwebb_detector1's `persistence` step from `_uncal`: full raw downloads per exposure, and the archived
  `_cal` products already carry the afterimages with no DQ flag at their pixels (observed), so the test has to work
  on what the mosaic was built from.
- Inspecting the level-3 `_i2d` only: the mosaic hides which exposures contribute, and an afterimage that lands on
  the same sky in two dithers looks like a real source there (`n0022`).

**Evidence** (docs/fields/sunrise.md; all `derived`).
- `n0022`: seen only in o010 dithers 3 and 4 of each SW filter. A bright galaxy lit the same pixels in dithers 2 and
  1 (F150W 256 and 268 vs 5.0 and 2.1, i.e. 2.0 % and 0.8 %); the dither geometry puts both afterimages on one sky
  position.
- `n0150`: in each epoch it appears only in the exposure after a saturated star (19–31 saturated pixels) sat on that
  pixel: o010 d4 (SW and F277W), o120 d2, o052 d4, at 0.04–0.07 % of the star's flux. Its wandering position is the
  per-epoch dither vector.
- Controls: `n0153` (39 of 48 exposures detected, 0 suspect, 36 clean) and Earendel (4 detected, 4 clean; S/N ≈ 5 per
  exposure, a weak control) are `on_sky`.
- Synthetic test: a 1 % afterimage of a saturated star is `persistence`, a real faint source `on_sky`.

**Revisit if.** A candidate is detected in exposures whose earlier ones were dark at its pixel, but fades within one
visit (fast transient vs. a lookback that is too short), or afterimages appear > 3 h after the illumination.

## D-040 `find_images` refines grid cells on folds; frame offsets move map models; MACS0416 image list open (2026-10-08)

**Decision.**
- `find_images(..., refine_arcsec=0.02)` keeps the triangle scan on the deflection grid. It then subdivides, into
  sub-cells of at most 0.02″ (below the 0.05″ merge radius), every cell that lies within one cell of a critical curve
  and whose mapped bounding box, widened by its own size, contains the source. The model is evaluated directly at the
  sub-cell nodes. A critical curve is detected where the mapped triangles change orientation. `refine_arcsec=0`
  restores the plain scan.
- Frame offsets go through `LensModel.shift_frame`; `MapLensModel.shift_frame` also moves its WCS. Map models look
  their maps up by sky position, so moving only the reference point left the maps where they were. `radial` and
  `arcs` now apply the offset for every model, not only Lenstool ones.
- `macs0416-cats` gets `frame_offset_arcsec` (0.208, −0.025) and its image list is opened.

**Alternatives rejected.**
- A finer global grid (0.1″): about 6× the memory and time on every field, against the 2M-cell rule.
- Dropping system 26 by hand.

**Evidence.**
- MACS0416 `validate`: image-plane rms 1.572″ → 0.760″ (1.06× the quoted 0.72″; the gate passes). System 26 is
  solved as 3 images. Shared matches drop from 6 to 3 (only system 122 is left).
- The other eight models give the same rms with and without refinement (SMACS 0.318″, El Gordo 0.754″, Abell 2744
  0.427″, MACS1149 0.673″, MACS0717 3.21″, Abell S1063, Abell 2744 CATS), except Abell 370 (10.69″ → 10.66″).
- Grid alignment: a synthetic fold pair gives the same image count on 0.5″, 0.3″ and 0.25″ grids of different extents
  (an 8 × 8 split with 0.0625″ sub-cells did not).
- `validate` wall time: Abell 2744 29 s → 41 s, El Gordo 16 s → 18 s, MACS0416 2.6 s → 3.9 s, SMACS 14 s → 12 s.
- Tests: a synthetic fold pair about 0.1″ apart inside one 0.5″ cell is missed by the plain scan and found, with opposite
  parities, by the refined scan. A frame offset leaves a map model's model-frame deflection unchanged and moves it on
  the sky.
- No earlier result was affected: no map model had an offset before this change (MACS1149 and MACS0717 are under
  0.1″). MACS0416's radial screen was re-run in the JWST frame with this code after pinning its offset: still null
  (D-038).

**Revisit if.**
- A fold image is still missed with `refine=8` (cusps with three merging images may need recursion).
- Map models with rotated WCS grids appear (`MapLensModel` rejects them today).


## D-041 MACS0717 screens: null; model copies of unpredicted images; CATS-only extra images (2026-10-08)

**Decision.** The MACS0717 screens (`macs0717-cats`, VENUS program 6882 `jw06882-o029_t063`, 10 bands) are a null
result (docs/fields/macs0717.md). There were 51 flags, 0 surviving. Two rules (ASSUMPTION: thresholds):
- **Model copy.** A predicted image within 1.75× the model's image-plane rms (5.6″ here) of a catalogued, detected
  but unpredicted image of the same system is the model's copy of that image, not a missing counter-image. In
  MACS0717, 29 of the 51 flags are copies, 1.6–5.6″ from catalogued images.
- **Model-dependent extra image.** An extra image predicted by one model but by neither of two independent models
  solved from their deflection maps (here Sharon v4cor and Keeton v4) is model-dependent and untestable. This extends
  the D-036/D-037 μ rule to image existence.

**Alternatives rejected.** Reporting the raw `absent` / `confused` classes. Matching within 1.5″ for a model whose
image-plane rms is 3.2″: 89 of 132 catalogued images have no prediction within 1.5″. With `--match-arcsec 3.2` the forced-photometry flags
drop to 4 (plus the 2 system-65 flux-ratio flags, which do not depend on the radius), and all were vetted.

**Evidence.**
- `validate`: κ 3.2e-3; μ(z=2) 1.2e-3; image-plane rms 3.21″ against the quoted 2.41″ (the gate passes). With the
  D-040 fold refinement, rms and the 16 shared matches are unchanged, so the copies are not solver-grid misses (the
  `images` run itself used the 0.25″ grid before D-040). Frame offset
  (F200W, 63 matches): +0.023″, −0.061″, so no `frame_offset_arcsec`.
- `images` (forced F277W): 199 predicted. Classes: recovered 25, confused 38, absent 6, undetectable 13,
  no_reference 54, inconsistent_reference 5, off_image 3.
- `fluxratio`: system 65 under- and overluminous; the rest are 26 consistent, 82 resolved and 22 untestable.
  65.2's JWST source is 0.6″ from its catalogued position. Corrected, the flux ratio is 0.7–1.1. μ(65.1)/μ(65.2) is
  0.37–17 across six models (untestable).
- `radial` (no photo-z): 11 peaks against a random mean of 16.7; max 5 lines; p ≥ 0.70.
- Breakdown of the 51 flags: 29 model copies, 5 untestable μ, 6 CATS-only extra images, 6 below sensitivity, 1 with a
  counterpart inside the position uncertainty, and 4 rows for the system 65 pair. Cutouts were inspected for every flag.

**Revisit if.**
- A JWST-era MACS0717 lens model appears.
- A MACS0717 photo-z catalogue (DJA or a team release) appears: re-run `radial` with the background cut.
- `forced_check` ties its match and search radii to the measured image-plane rms when that rms exceeds 1″.

**Addendum 2026-10-14 (UNCOVER v2.0 as Abell 2744's second model).**
- **Decision.** Tracing fold pairs from the mean back-traced source fails the image-list gate spuriously.
  `a2744_uncover.py` uses the magnification-weighted least-squares source (`wls_source`, the first-order image-plane
  χ² source). UNCOVER v2.0 then passes (rms 0.84″ vs 1.5 × 0.60″, every traced image re-imaged) and is Abell 2744's second model for this rule.
- **Alternatives rejected.** The mean source (rms 4.8″: 34.1 and 80.1 traced outside the caustic); the best single
  image's source (1.01″); a finer solver grid (no change). Changing `lensmodel.imageplane_residuals` itself would
  shift every validated field's gate numbers; revisit if another map model fails the same way.
- **Evidence.** 4.2c is predicted by both models (μ 8.7 / 8.8, 0.93″ apart) and not seen at ≥ 0.5× (injections), so
  it is not model-dependent; its vetted explanation (member 34423 systematics) stands. CHANGELOG 2026-10-14.
- **Revisit if** a UNCOVER release refits the 4.2c region, or a third model removes 4.2c.

**Addendum 2026-10-14 (cross-model rule on every forced row).**
- **Decision.** Both models' extra images are force-photometered under the same rules. Two extra images of one
  family from the two models count as one when they are ≤ 1.5″ apart (`both`; 1.8× UNCOVER's rms) and as a loose
  match ≤ 3.0″ (1.3× its largest own-list residual). Both radii are ASSUMPTIONs. An unmet extra image is
  model-dependent unless the other model predicts an extra image within 3.0″. When the other model's image there is
  a catalogued image (≤ 3.0″) or is recovered, the unmet image is that image misplaced (ordinary).
- **Alternatives rejected.** 1.5″ alone: it calls three recovered Bergamini images (3.3, 28, 34.1) one-model
  predictions, although UNCOVER predicts them 1.5–2.7″ away.
- **Evidence.** Abell 2744: 4 `absent` rows. Two are Bergamini-only (700.1). Two are UNCOVER fold-pair merges (3.1,
  34.1) that Bergamini places ≤ 2.7″ away on observed images. No row is unmet in both models. CHANGELOG 2026-10-14.
- **Revisit if** a third Abell 2744 model is added, or a both-model `no_reference` row gets a reference flux.

**Addendum 2026-10-21 (Sharon v4cor and Keeton v4 pinned, `scripts/m0717_models.py`).**
- **Decision.** The scratch second models of this record are now `MODELS` entries `macs0717-keeton` (pinned for use:
  both frozen gates pass) and `macs0717-sharon` (G2 failed as written because image 6.2 lies off the 200″ maps; used
  only inside the map footprint, a post hoc condition as for `macs1149-sharon`). Deflection and κ sha256 are the
  2026-10-08 SOURCES.md values; the z = 2 magnification maps are new pins.
- **Evidence.** G1: κ / |μ(z = 2)| median relative differences 0.00002 / 0.00002 (Keeton) and 0.00025 / 0.024
  (Sharon). G2 on the 26 Table 7 spectroscopic images: Keeton rms 0.853″, 8 / 8 systems; Sharon 0.639″, 7 / 8 (0.662″
  on the 7 inside the maps). CATS on the same images: 2.77″. `results/lens_consistency/macs0717_models/RESULT.md`.
- **Revisit if** a forced row lies outside Sharon's footprint, or a JWST-era MACS0717 model is released.

## D-042 MACS0416 counter-image and flux-ratio screens: null; CANUCS Lenstool model as second model; two-plane check of foreground deflectors (2026-10-08)

**Decision.** With the image list open (D-040), the MACS0416 `images` and `fluxratio` screens are a null result
(docs/fields/macs0416.md). The JWST-era CANUCS Lenstool model is the independent second model for this field:
Rihtarsic et al. 2025, 222 potentials and 111 spectroscopic systems, loaded with `LensModel.from_par` from the
released best-fit parameter file. A flagged extra image that the CANUCS model does not predict is CATS-only, and
therefore model-dependent. This relaxes D-041's two-independent-models requirement to one JWST-era model (ASSUMPTION:
the JWST-era model, with spectroscopic constraints only, is the stronger check).

**Alternatives rejected.**
- HFF-era second models (Sharon, Keeton). The CANUCS model is fitted to JWST positions and uses spectroscopic
  constraints only, so it is the stronger independent check.
- Treating the system-27 extra images (predicted S/N 247–341) as missing: the CANUCS model reproduces system 27 with
  exactly its 3 catalogued images (within 0.2–0.66″) and predicts nothing at the 4 CATS extras (its nearest images are 4.5–12.8″ away).

**Evidence.**
- `images` (F277W forced, CANUCS photo-z, frame offset applied): 143 predicted. Classes: recovered 3, confused 13,
  absent 2, undetectable 3, no_reference 6, inconsistent_reference 3.
- Vetting (CANUCS model, cutouts inspected):
  - system 27: 2 absent and 2 confused, all CATS-only. Empty sky at the μ 26 and μ 19 positions in all bands;
    cluster members at z 0.3–0.4 lie 0.6–1″ away, so these are galaxy-scale caustics of the CATS model.
  - one system-55 row is CATS-only (5.3″ from any CANUCS image); the other two are reference-quality cases.
  - 5 confused rows (systems 34, 45, 47, 132, 133) are predicted by both models (within 0.25–0.9″) but lie next to
    neighbours 30–900× brighter, so they are untestable.
  - 4 confused rows (systems 1, 15, 122 ×2) are copies, under a second copy rule (ASSUMPTION). The flagged prediction
    is the only unmatched prediction of its system near exactly one unmatched catalogued image, it has the same
    parity, and it lies within 6.5 × the image-plane rms (1.6–4.7″ here, against 0.76″ rms). D-041's 1.75 × rms rule
    does not cover them.
  - System 45 (μ −9) lies 1.89″ from the unmatched 45.2 but with the opposite parity to CATS's μ +2.3 there. CANUCS
    gives μ +14.9 at 45.2, so 45.2 is near-critical and its parity is model-dependent. Untestable (D-036).
  - **System 51's fourth image is explained.** Both models predict it (μ 3.8 CATS, 5.5 CANUCS; 0.66″ apart). After
    isophote subtraction of the neighbour, nothing is seen; the expected signal, derived from 51.1–51.3 photometry,
    is 9–25σ. In CANUCS the image comes from potential 8757, the z_spec 0.268 galaxy 0.85–0.88″ away, modelled as a
    cluster member at z 0.396. CATS is a map model and cannot be decomposed; its fourth image lies 0.85″ from the same
    galaxy, consistent with the same origin (not tested). In a two-plane model (that galaxy at z 0.268, the rest of CANUCS at z 0.396), the
    image persists at σ 102 and 81 km/s (fitted, and rescaled to its luminosity distance). At σ ≤ 60 km/s the system
    has exactly 3 images and 51.3 is matched within 1.5″. That is 26 % below the rescaled σ, inside the assumed ±30 %
    scatter. Record: docs/candidates/macs0416-system51-fourth-image.md.
- `fluxratio` (F150W/F444W): 32 consistent, 58 resolved, 23 untestable.
  - 45.2 is overluminous: a 0.2″ aperture 0.5″ from a bright compact galaxy (the cutout), and its CANUCS match has
    spectroscopic z 2.544 but photo-z 0.38, a blend.
  - 45.1 is flagged as underluminous only relative to 45.2.
  - 38.1 is chromatic at S/N 2–3.5 (marginal).
- Tally: 143 predicted images and 116 flux-ratio images screened; 21 flags (18 forced rows, 3 flux-ratio images);
  **0 surviving**.
- **Rule:** before treating a missing extra image as a flag, check the potentials that produce it against their
  spectroscopic redshifts. A foreground or background galaxy modelled as a member is tested in a two-plane model with
  its σ rescaled and scanned over the scaling-relation scatter (ASSUMPTION: ±30 %).

**Revisit if.**
- The CANUCS model is pinned as a `MODELS` entry with its image list, and `images` is re-run on it directly.
- A multi-plane re-fit, or a measured velocity dispersion of CANUCS 3101008 (z 0.268), gives σ ≳ 65 km/s. System
  51's explanation then fails and the fourth image becomes a flag again.

## D-043 Abell 370 and Abell S1063 radial screens: null; Gaia-seeded spike veto (2026-10-08)

**Decision.**
- Only the `radial` screen applies to `abell370-cats` and `abells1063-cats`: their image lists are gated off (D-035;
  rms 10.7″ and 11.8″).
- Data: CANUCS 1208 `o002_t001` (Abell 370; 8 bands) and GLIMPSE 3293 `o001_t001` (Abell S1063; 9 bands, F200W 70 ks).
- Abell 370's frame offset (−0.121″, −0.015″) is pinned. S1063's 0.075″ is under the 0.1″ threshold and not pinned.
- Photo-z come from CANUCS DR1 (Abell 370) and DJA v7.5 eazy (S1063).
- `exotic_screens.py radial --spike-stars` (stars from `scripts/gaia_stars.py`) adds Gaia DR3 stars (G < 17) to
  `spike_segments`, with spikes up to 60″; catalogued stars keep the 20″ cap (D-034), and a Gaia star within 1″ of a
  catalogued one is seeded once (ASSUMPTIONs; G is used in the AB spike-length law without a colour term). The
  pipeline catalogue misses saturated and off-mosaic stars.
- Both fields are null (docs/fields/abell370.md, docs/fields/abells1063.md).

**Alternatives rejected.**
- Abell 370 data: MAGNIF 2883/3538 and JUMPS 5890 (medium bands only); 5324 (shallow).
- S1063 data: 1840 (about 1 ks per band).
- Taking the 15-line Abell 370 centre (p 0.0) at face value. It is a straight chain of diffraction-spike segments
  along the column axis (PA 62°), at a low-weight seam, pointing at a Gaia G = 13.7 star off the F200W mosaic.

**Evidence.**
- Abell 370 `radial`: 15 lines, p 0.0, with or without photo-z. Two flags, both instrumental (cutouts inspected):
  - the 15-line chain above;
  - an 8-line chain (p 0.01 with aper50 S/N ≥ 5) from a saturated Gaia G = 12.7 star that has no catalogue entry.
  - `spike_segments` missed both stars: it seeds only from catalogued point sources brighter than AB 20 and caps
    spikes at 20″, while these segments reach 12–37″.
  - With the committed Gaia-seeded veto: 82 segments dropped, 140 anti arcs, max 5 lines, p 0.945 (111 and p 0.495
    with photo-z). The scratch test used a fixed axis (63.4°) and removed 166 catalogue rows. The code estimates the
    axes and vetoes selected segments; the resulting anti counts and p-values are the same.
- S1063 `radial`: max 4 lines, p 0.435; with photo-z, max 3, p 0.95. The Gaia veto (11 stars) leaves it unchanged.
- Tally, counted the same way for both (anti arcs after all spike vetoes, no photo-z): Abell 370 140, S1063 51; 2
  flags; **0 surviving**.

**Revisit if.**
- A low-weight veto (relative WHT < 0.5) or an aper50 S/N floor is added to `radial`.
- Fitted redshifts for Abell 370 or S1063 appear, which would open their image lists.
- A photo-z catalogue for S1063 with fewer blends on arcs appears: DJA puts 13 of 46 matched images at z < 0.6.

## D-044 CANUCS Lenstool models pinned as `macs0416-canucs` and `abell370-canucs` (2026-10-08)

**Decision.**
- The CANUCS DR1 Lenstool best fits are now `MODELS` entries: `macs0416-canucs` (Rihtarsic et al. 2025, image-plane
  fit) and `abell370-canucs` (Gledhill et al. 2025, source-plane fit).
- Each entry pins the best-fit parameter file, the Lenstool multiple-image file and the input parameter file by
  sha256 (`lensmodel.MACS0416_CANUCS`, `ABELL370_CANUCS`). sigpos is read from the input file: 0.49″ and 0.3″.
- Frame offsets of the image lists against the JWST F200W catalogues, measured with the D-034 method:
  - MACS0416: (−0.008″, +0.052″), not pinned (under 0.1″);
  - Abell 370: (−0.148″, +0.002″) ± 0.017″, pinned.
- `abell370-canucs` has its image list gated off (`image_list_ok: False`).
- `validate` drops the Lenstool χ² reference from the image-plane comparison when the best fit is source-plane.
- `macs0416-canucs` is the independent JWST-era second model for the D-036/D-041/D-042 rules. It can be used
  directly by `validate` and `images`. The two-plane checks of D-042 still need scratch code (TASKS).
- MACS1149 has no released CANUCS Lenstool parameter file (only maps), so it is not pinned.

**Alternatives rejected.** Using the CANUCS deflection maps as map models: the parameter files allow galaxy-scale
tests (removing or rescaling one potential, D-042), which maps do not.

**Evidence.**
- `validate --model macs0416-canucs`: 303 of 303 images solved, rms 0.512″ (max 1.53″), χ²pos 330.8 against
  Lenstool's 344.30. The match is within 4 %, as for SMACS / El Gordo / Abell 2744 (D-030).
- `validate --model abell370-canucs`: 115 of 115 images solved, rms 2.32″ (max 17.9″). This model was fitted in the
  source plane (χ² 192.6 there), so its image-plane residuals are not comparable, and its image list is gated off.
- Network tests: the MACS0416 reproduction holds within 10 %, and both models load with finite redshifts for every
  image.

**Revisit if.**
- CANUCS releases a MACS1149 parameter file.
- The Abell 370 image list is needed: it requires an rms gate like D-035's.

## D-045 Lenstool MCMC posteriors (`bayes.dat`); the Abell 2744 multiplicity residual is model resolution, not an anomaly (2026-10-08)

**Decision.** `lensmodel.read_lenstool_bayes` / `posterior_par` rebuild a published model at any MCMC sample, and
`lens_consistency.py --model <m> posterior --systems ...` runs the image-plane solve for selected families over
best.par plus N random chain rows. Pinned chains: Abell 2744 Bergamini+2023b (70 MB, 133,690 rows) and El Gordo
Caminha+2023 (1.7 MB, 10,000 rows); kept apart from the model file sets so `validate` does not fetch them.
- `O<i> : <key> (<unit>)` columns replace keywords of the i-th best.par potential; units are checked.
- `Pot0 sigma` / `Pot0 rcut` rescale every potfile member (a best.par potential with `mag`) by the ratio to the
  best-fit value. That is exact for Lenstool's power-law scaling relations. The reference is best.par's own chain
  row (Abell 2744) or, for a thinned chain without it, the member at the input file's `mag0` (El Gordo: BCG 1758,
  `potfile_mag0` 17.9852 from CDS `to_sample.par`).
- `Redshift of <id>` columns set every family on that `z_m_limit` line.

**Alternatives rejected.**
- The CANUCS 100 MCMC deflection maps (MACS0416 only, 100 × 3 maps); a parametric chain covers every model we
  rebuild from best.par, at any number of samples.
- A full `images` run per sample: about 70 s per Abell 2744 sample at a 0.25″ grid, dominated by the deflection
  grid, so `posterior` restricts the solve to the named families.

**Evidence.**
- Abell 2744: best.par is chain row 117,869 (Chi2 146.60; parameter differences 0, family redshifts < 1e-3), and its
  rebuilt deflection field equals best.par's. Two random rows give χ²pos 173.5 and 173.0 against the chain's 178.1
  and 174.1 (all 149 images, 0.25″ grid). The rebuild is validated.
- El Gordo: the chain medians match best_fit.par (σ* 289.97 vs 289.48 km/s; O1 σ 1041 vs 1041), but the chain's
  `Chi2` column (54–77) does not track our χ²pos (93–106 at the 0.621″ sigpos of D-030), not even in rank. The
  column's definition is unknown, so the El Gordo chain is **not validated**; don't use it for conclusions yet.
  - **Resolved (2026-10-08, PR #69, issue #68):** the sampling run (`to_sample.par`, `forme -10`) uses an
    image-plane χ² with σ² = a·b from the image list (σ = 0.621″ for 37 images, 1.2421″ for 19), not D-030's uniform
    0.621″. Lenstool's `chi2_img` and `bayesapp.c` (git-cral.univ-lyon1.fr/lenstool, v8.15.6) give
    ln(Lhood) = −(Chi2 + Σ 2 ln(2π a b))/2; the file's Σ ln(2π a b) = 75.904 matches the chain's offset exactly. Our
    image-plane χ² with σ² = a·b reproduces `Chi2` for three random rows (60.05/60.00, 67.50/67.47, 72.45/72.36), and
    is 52.0 for best_fit.par (the chain minimum is 54.2). The El Gordo chain is **validated**. best_fit.par's
    `Chi2pos` 80.22 still corresponds to the uniform 0.621″ of D-030.
- Abell 2744 `posterior --systems 3.2,34.1,700.1,4.2 --samples 12 --seed 1` (`model_prediction`):
  - 3.2a/b, 34.1a/b, 700.1a/b stay a `shared_match` in 13 of 13 models. The MCMC spread does not split them.
  - Independent model, CATS v4.1 maps (scratch run with Bergamini's image list): 34.1a/b **split** (μ +24.9 / −21.7,
    residuals 2.8″ / 2.0″), so the 34.1 merger is model-dependent. 3.2a/b merge there too (μ 230): both models put
    the source on the caustic. The catalogue's 3.3a–3.1a–3.2a | 3.2b–3.1b–3.3b ordering is mirror-symmetric, so the
    critical curve must run between 3.2a and 3.2b. A point-source model with the source within its position error of
    the caustic predicts one merged image 0.28–0.53″ from each, which is within the model's 0.43″ rms. 700.1 has a
    free redshift (posterior 1.2–3.9); CATS does not reproduce it at Bergamini's z.
  - 4.2c (D-036) is present in 12 of 12 samples, 0.22–0.36″ (p16–p84) from the best-fit position, with μ 8.7–10.0.
    The statistical spread is far below the galaxy-scale systematics that explain it (μ 3.9–28.7 under ±30 %).
- Verdict: the D-030 multiplicity residual is model resolution at folds, not an anomaly. **0 surviving.**

**Revisit if.**
- ~~The El Gordo `Chi2` column is understood~~ (done: `forme -10` σ² = a·b; chain validated, PR #69).
- An independent model with 3.2's source well inside the caustic still merges 3.2a/b.
- Speed: `imageplane_residuals` exposes its predicted images (`posterior` currently solves each family twice).

## D-048 Orphan image pairs (blind dark-deflector screen): null in MACS0416, MACS1149 and Abell 370 (2026-10-08)

**Decision.**
- `scripts/orphan_pairs.py` ([docs/orphan_pairs.md](docs/orphan_pairs.md)) searches CANUCS DR1 catalogues for close
  pairs (0.3–3″) with matching SEDs that no published multiple-image system explains and that have no visible galaxy
  between them ("orphans"). It is a lens-model-independent dark-lens screen. The hypothesis is an unseen compact
  deflector; an SED match is not evidence of lensing.
- Result: **null**. SED-matched close-pair counts are explained by redshift clustering, and orphan counts match
  chance. Nothing goes to `/vet-candidate`.
- Thresholds (ASSUMPTIONs, in the script):
  - summed F277W+F356W+F444W S/N ≥ 10 and S/N ≥ 10 in at least 8 bands;
  - SED match: ≥ 8 shared bands, χ² probability ≥ 0.01 with a 3 % error floor and free normalisation, overlapping
    16–84 % photo-z intervals;
  - visible lens: a catalogued source at least 0.3″ from both members that lies within 0.3″ of the midpoint,
    within 0.3″ of the joining segment, or inside the circle with the pair as diameter.
- Inputs are pinned by URL and sha256 in `FIELDS` and downloaded on first use; the CANUCS Abell 370 image list is
  moved to the JWST frame with the `abell370-canucs` offset (D-044); its image-plane gate concerns model
  constraints, and here only positions are used (to mark pairs near a published image). A source with an invalid F277W, F356W or F444W
  measurement (NaN, or error ≤ 0) has no summed S/N and is dropped. Re-running with these rules (after review)
  reproduced every count below.

**Alternatives rejected.**
- A midpoint-only visible-lens rule: galaxies sitting between pair members were missed (first contact sheet).
- No per-band S/N cut: 25 % of random far-apart pairs pass the SED test.
- Random-position nulls alone: they ignore galaxy clustering at one redshift.

**Evidence** (`derived`; per-field `summary.json`, three contact sheets inspected by eye):

| | MACS0416 | MACS1149 | Abell 370 |
|---|---|---|---|
| Sources kept / catalogue rows | 1942 / 14149 | 1711 / 12851 | 1391 / 13567 |
| Pairs at 0.3–3″ / SED-matched | 1664 / 62 | 1418 / 98 | 1028 / 35 |
| SED-matched expected: null (a) shifted copies / (c) same-photo-z real pairs 10–30″ | 42.3 / 69.7 | 44.8 / 83.3 | 27.8 / 48.0 |
| Orphans observed / expected under null (c) | 11 / 13.5 | 18 / 12.9 | 9 / 8.1 |

- The SED-match excess over null (a) (significant in MACS1149) disappears under null (c): same-redshift groups.
- All 35 inspected top orphans sit where the CATS (or CANUCS) model gives |μ| ≈ 1–2.5 with no parity flip. They are
  knots of one galaxy, group members, or faint chance matches the null predicts.
- A dark deflector making the observed separations (θ_E 0.36–1.47″ at the cluster redshift) would need σ ≈ 130–300
  km/s; a normal galaxy of that mass would be m ≈ 17–20, 9–12 mag above the detection limit (`hypothesis`-level
  scaling, σ* = 180 km/s at m*(F160W) = 19, ASSUMPTION).
- Limits: the same-galaxy rule misses knots of large irregulars; pairs closer than ~0.6″ share aperture light; cluster
  cores are excluded by catalogue flags; 86–90 % of each catalogue is cut; no injection-recovery yet, so this is a
  count-level null, not an upper limit on dark deflectors.

**Revisit if.**
- Injection-recovery turns the count null into a limit (TASKS).
- CANUCS segmentation maps replace the same-galaxy rule, or deeper / core photometry is used.
- CANUCS v2 releases the MACS1149 image list, or spectroscopy targets an orphan pair.

**Addendum 2026-10-16 (Abell 2744).**
- **Run.** Same rules, frozen in `results/orphan_pairs/abell2744/PLAN.md` before the run, on DJA v7.2-fix (66,172
  rows, 9,385 kept). Published images come from CATS v4.1 and Bergamini+23; alignment is 0.02″.
- **Result: null.** 51 orphans observed against 43.8 / 42.4 expected under nulls (e) / (f); P = 0.16 / 0.11.
- **No positive control.** No SED-matched pair lies near a published image.
- **Injections** (`inject_pairs.py --fields abell2744`, efficiencies only): per-deflector efficiencies match the
  deep fields within noise (point 0.30 % / 0.95 % at 0.3″ / 0.7″). Completeness is low here as in every field. The
  cluster magnification of injected images is not modelled, and no pooled limit is written for a cluster field.
- **Nulls.** The excess over null (a) persists under (c) and (d) (P 0.016 / 0.021) and disappears under (e) and (f).
- **Escalation.** 11 of the top 15 orphans meet the frozen escalation rule. They are unvetted (step postponed) and
  not candidates until vetted. A sheet pass finds 5 companions or group members, 2 possible knots, 1
  artefact-affected pair and 3 featureless faint pairs. Details are in `RESULT.md`.
- **Per-pair vetting (2026-10-16, `vetting/RESULT.md`).**
  - Classes: 4 ordinary (#5, #6, #14 by E1 strip excess; #1 by a visual artifact flag), 7 insufficient data (#2, #3,
    #7, #9, #10, #11, #13), 0 anomalies.
  - E2 surface brightness could not be validated on published images, and no public spectra exist at any member.

## D-046 Multi-plane lens models: `LensModel.split_planes` and `MultiPlaneLensModel` (2026-10-08)

**Decision.**
- `lensmodel.MultiPlaneLensModel` holds several `LensModel` planes (one `z_lens` each) in one frame and cosmology and
  solves the standard multi-plane lens equation and its Jacobian recursion (Schneider, Ehlers & Falco 1992, ch. 9),
  with `D_ij / D_j = 1 − D_M(z_i) / D_M(z_j)` (flat ΛCDM).
- `LensModel.split_planes({name: z}, v_disp={name: σ})` moves named potentials to their own redshift, optionally
  with a new σ. A potential moved behind another plane is delensed: put where the ray through its fitted
  (observed) centre crosses its plane, so it is still seen where it was fitted; shapes and the other planes'
  positions are kept (ASSUMPTION). This turns the D-042 rule ("check every potential that produces an extra image against its
  spectroscopic redshift") into library code.
- `find_images`, `backtrace_images` and `imageplane_residuals` now go through `lens_map` / `source_points` /
  `source_grid`, so they take either model. The Jacobian is no longer assumed symmetric. Single-plane results are
  unchanged (tests).
- The ray positions θ_i do not depend on the source redshift; only the weights D_is/D_s do. So a `DeflectionGrid`
  of a multi-plane model stores α_i(θ_i) per plane (shape `(n_planes, n, n)`, cached by a content hash) and serves
  every source redshift, as for one plane. Computing it costs about one single-plane grid (≈70 s for the
  222-potential CANUCS MACS0416 model on a ±59″, 0.1″ grid).
- `MultiPlaneLensModel` has no `evaluate` and no scalar `z_lens` (`z_planes` instead), so single-plane-only code in
  `lens_consistency.py` fails loudly rather than silently. Sub-planes made by `split_planes` get their own hash;
  `find_images` refuses a grid computed for another model.

**Alternatives rejected.**
- Re-tracing every plane at every grid node per source redshift (first version): unnecessary, see above.
- lenstronomy `MultiPlane`: its dPIE-like profiles are not parametrised as Lenstool's, and our ported dPIE already
  reproduces Lenstool's χ² (D-030, D-044). The recursion itself is a few lines on top of it.
- Keeping two-plane checks as scratch code (D-042): not reproducible.

**Evidence** (`model_prediction`).
- Offline tests: `split_planes({})` reproduces the single-plane `find_images` and `backtrace_images`; the multi-plane
  Jacobian matches finite differences of `source_points` (and is not symmetric); a plane behind the source does not
  lens.
- MACS0416 system 51 (CATS positions, CANUCS model, potential 8757 moved to z 0.268; scratch
  `macs0416/mp_check2.py`): σ as fitted (102 km/s) keeps the fourth image (μ 6.1, 0.12″ from the single-plane
  position); σ 81 km/s gives 4 images; σ 70 and 60 km/s give 3, with 51.3 matched at 1.58″ and 1.49″. This
  reproduces the D-042 scratch result, except that the scratch found a fifth (faint) image at the fitted σ.

**Revisit if.**
- A multi-plane field run needs speed: interpolate per-plane deflection grids for the seeds.
- A model needs more than dPIE potentials on the extra planes.

**Addendum 2026-10-14 (`--plane NAME=Z[:SIGMA]` on `lens_consistency.py validate` and `images`).** Repeatable;
applies `split_planes` to the pinned Lenstool model (map models refused). `validate` then skips the published
single-plane κ and μ map comparisons and records the planes under `multi_plane`; outputs go to
`<model>/planes_<hash>/`, so single-plane tables are never overwritten. Use (`model_prediction`, other parameters
**not** refitted, so Δχ² is an upper bound on the cost): MACS0416 CANUCS, 303 images, σ_pos 0.49″. Single plane
χ² 330.81 (Lenstool 344.30), rms 0.512″. Potential 8757 (the z_spec 0.268 foreground galaxy behind the system-51
fourth image, D-042/D-046) at z 0.268: σ 102 km/s (fitted) χ² 330.55; 70 km/s 332.71; 60 km/s 336.44 (rms 0.516″).
Most affected: CANUCS D52a–c (1.3″ from 8757), first images +0.28–0.45″ at 60 km/s (1.11–1.18″ → 1.46–1.56″). The two-plane
σ ≤ 60 km/s solution that removes the fourth image (D-046) costs Δχ² ≤ +5.6 on the catalogued images: no tension.

**Addendum 2026-10-14b (off-plane potentials, `scripts/offplane_potentials.py`).** Every potential of a Lenstool
model is matched to the nearest CANUCS DR1 `Z_SPEC` source within 0.5″; off-plane at |Δz|/(1+z_cl) > 0.02; the
z_spec is trusted only when `Z_ML` agrees within 0.15 (1+z) (ASSUMPTIONs: compiled z_spec, e.g. Shipley+18 or MUSE
images, can belong to a blended lensed image). MACS0416 CANUCS: 12 / 222 off-plane, 3 rejected by photo-z
(misassociated). Seven near-cluster interlopers cost Δχ² +0.17. O5 (σ 136, z 0.114) costs Δχ² +39 to +73 for σ
136–175 km/s (best 150), +112 at the lensing-equivalent 123 and +924 at the luminosity-rescaled 67 km/s, unrefitted (`model_prediction`): a single-plane fit absorbs a
foreground galaxy's multi-plane deflection, so an unrefitted Δχ² is only an upper bound and is no tension at ~12 %.
Rejected: refitting the CANUCS model here (no Lenstool optimiser in the repo; not needed while no flag depends on O5).

**Addendum 2026-10-14c (MUSE redshifts for Abell 2744 and El Gordo).** For models with a published MUSE redshift
catalogue (`offplane_potentials.MUSE`: Richard+21 J/A+A/646/A83 for Abell 2744, Caminha+23 J/A+A/678/A3 for El
Gordo) the z_spec is trusted by its quality flag (zconf ≥ 2; q_zsp 2/3/9), matched in the model's HST frame.
Abell 2744 Bergamini+23: 72 / 180 potentials with z_spec, 0 off-plane. El Gordo Caminha+23: 188 / 265 with z_spec,
20 off-plane, all at z 0.6305–0.6355: exactly the paper's foreground group of 20 ellipticals, which the authors
modelled on the cluster plane with the member scaling (their Table: no-foreground χ² 103 vs reference 80, refitted).
Group moved to z 0.6325, unrefitted (`model_prediction`): χ² 82.5 → 100.3 (fitted σ), 98.2 (σ × 1.10, best of
0.82–1.40), 110.5 (× 0.913, lensing-equivalent at z_s 3), 127.8 (× 0.822, luminosity-rescaled); group removed
(σ × 0.01) 312.9. Per-galaxy z_spec instead of the group mean: 100.2. `images --plane` (σ × 1.10): the extra
predicted images of systems 3, 4, 6, 12, 15, 16, 19 persist within 0.7″; new are a system-17 extra (μ −1.9, F200W
29.2, below the 27.1 depth proxy) and a system-23 one beside 23a; no robust new extra image. **No tension**; the
foreground group is a known, modelled structure, so it is a positive control for the listing.

**Addendum 2026-10-14d (Abell 370 and SMACS 0723).** `MUSE` gains Richard+21 `a370.dat` for `abell370-canucs`
(median 0.03″ to the CANUCS potentials, so no offset; CANUCS DR1 z_spec stays available with `--field abell370`) and
Caminha+22 `muse.dat` for `smacs0723-iclv2`, with `frame_offset_arcsec` (−0.218, +0.262)″: Mahler+22's frame sits
0.34″ from Caminha+22's positions (median over 24 members, MAD 0.009″); unshifted, the matches all lie at 0.31–0.48″.
Results (`model_prediction`, unrefitted; details `docs/fields/abell370.md`, CHANGELOG 2026-10-14):
- Abell 370 CANUCS: 257 / 284 with z_spec (MUSE), 36 off-plane. CANUCS DR1 z_spec gives 206 / 32 and agrees on
  every common potential except 180 (MUSE 0.4649 vs GLASS 0.38; photo-z 0.33) and 61 (CANUCS 0.2845, MUSE member).
  The model is source-plane optimised, so the back-traced image-plane χ² (σ_pos 0.3″, 115 images; 1003.0 single
  plane) is a relative metric only. 23 foreground potentials at z 0.30–0.35 (σ up to 171 km/s): 985.3 (−18);
  11 background ones at z 0.41–0.47: 1005.2 (+2); 273 at z 2.80 (a lensed source on a member, σ 58): 1000.4;
  180 at z 0.4649: +46 at its fitted σ 154, +12 at σ 100 against the same σ on the plane (the σ change itself moves
  χ² by −164). `images` is gated off for this model (D-044), so no extra-image test.
- SMACS 0723 ICLv2: 24 / 149 with z_spec, 2 off-plane. 610 (σ 77) is the BCG potential; its z 0.9366 entry carries
  Caminha+22's amplification 3.6, a lensed source on the BCG: misassociated, the ASSUMPTION's known failure.
  489 (σ 102, F200W 18.4) is a z 0.746 background galaxy modelled as a member. Both moved: χ² 31.18 → 31.33;
  `images --plane`: the same 75 images, no class change, largest shift 0.29″ (system 9, |μ| 20). **Null.**

## D-047 Exotic-lens predictions: wormhole and negative-mass signatures are searchable, warp signatures are not (2026-10-08)

**Decision.** Phase 1 ("predictions first") gives the screens closed-form targets. docs/exotic_lensing.md has the
"Wormhole signatures" and "Warp signatures" sections; `src/jwst_anomaly/exotic_sim.py` holds the simulator.
- **Searchable with JWST data** (each one a hypothesis to test, never a result):
  - W1, negative-mass dark lens: a radially stretched pair on one side of an empty centre, and no images of sources
    within 2θ_E. It calibrates the `radial` screen (D-031).
  - W2, Ellis pair with no deflector: it calibrates the dark-lens search and `fluxratio`.
  - W3, inverted microlensing: a compact lensed source vanishes for 2 t_E √(4 − u₀²) between caustic spikes. It
    calibrates the two-epoch transient screen (D-027), which needs a dimming class and ≥ 3 epochs to see
    spike, dip, spike.
  - W5, count deficit: fewer sources inside about θ_E. It calibrates the counts screen, run only around
    `radial` centres.
- **Not searchable:**
  - W4, the Ellis gutter: it needs 0.8 % photometry per epoch;
  - W6, µas centroid shifts and femtolensing;
  - every warp-drive signature: none has a published imaging or photometric prediction for a distant observer
    (docs/exotic_lensing.md "Warp signatures"). That branch is stopped.
- **Model.** Kitamura, Nakajima & Asada (2013) power-law family, α = ε̄/bⁿ with ε of either sign; magnification from
  Izumi et al. (2013). It covers point mass, negative mass, Ellis and phenomenological n > 2 in one lens equation.
  - Solutions are closed form for n = 1, with bisection on the monotonic branches otherwise.
  - A uniform-disk finite source caps the caustic spikes. It is an exact 1-D radial integral, split at the caustic
    and the disk edges, and it agrees with inverse ray shooting to < 1 %.
  - `inject_images` and `inject_light_curve` return `simulated` tables for injection-recovery.
- **Recommended injections** (ASSUMPTIONs, to be tuned on recovery):
  - radial: n = 1, ε < 0, θ_E = 0.3″, 1″ and 3″ (|M| ≈ 2 × 10¹⁰ to 2 × 10¹² M☉), painted from catalogue
    sources with β ∈ [2, 4] θ_E.
  - transient: an umbra (lensed flux → 0, keeping the blend) of 2–4 t_E, plus spikes of × 7.0 (ρ = 0.01),
    × 2.35 (ρ = 0.1) and × 1.53 (ρ = 0.3). These are `simulated` peaks from `exotic_sim`, not paper values.
  - Ellis or n ≥ 3 dips of 4 %, 14 % and 59 %, as a sensitivity floor.

**Alternatives rejected.**
- lenstronomy point-mass or power-law profiles: not a dependency here. They are parametrised by a positive
  Einstein radius, and ε < 0 is not supported as far as we checked. The closed forms take about 100 lines of
  numpy and are tested against the papers.
- Inventing a warp-bubble lensing or flash model: the charter forbids it, and the Alcubierre exterior is flat.

**Evidence** (`tests/test_exotic_sim.py`, offline; `simulated`).
- n = 1, ε > 0 reproduces (u² + 2)/(u√(u² + 4)).
- n = 1, ε < 0 reproduces Safonova et al.'s umbra at u < 2, the caustic at u = 2 and (u² − 2)/(u√(u² − 4)).
- Ellis: inner images at −0.618 and −0.532 carrying 3.4 % and 1.3 % of the flux (Abe 2010), and a 4.2 % gutter
  ("about 4 %").
- n = 10: demagnification onset at β = 0.1875 (paper 0.187 numerically; 2/(n+1) = 0.182 is their leading-order
  estimate, now `demagnification_onset_approx`), and 59 % depletion at β ≈ 0.70 (paper ~60 % at ~0.7).
  The exact onsets are 1.111 (n = 2) and 0.643 (n = 3).
- n = 3: 14.3 % depletion at β ≈ 1.12. The paper's text says "~10 %", but its own Fig. 2c (read from the figure
  pixels) bottoms out at A ≈ 0.865 ± 0.01, i.e. 13–14 %. The text rounds; no model difference.
- Izumi κ and γ; Abe's Tables 1–2 (R_E, θ_E, t_E for a = 10³ and 10⁵ km); the docs' physical scales;
  injection round trips through the lens equation.
- Converged spike heights. The first version, with a centred 2-D disk quadrature, overestimated them near the
  caustic: 9.24 → 7.01 (ρ = 0.01), 2.56 → 2.35 (ρ = 0.1), 1.74 → 1.53 (ρ = 0.3).
- Every citation was fetched from the arXiv API or Crossref (SOURCES.md "Exotic-lensing predictions (D-047)").

**Revisit if.**
- A paper gives a distant-observer electromagnetic prediction for a warp bubble.
- Injection-recovery shows a screen is blind to W1 or W3 at the recommended amplitudes.
- Cluster macro-magnification needs a lens model with shear plus a microlens instead of an isolated lens (W3 in
  caustic-crossing arcs).

## D-049 W1 injection-recovery through the `radial` screen: blind below about 10¹² M☉, weak limits above (2026-10-08)

**Decision.**
- `scripts/inject_radial.py` turns the eight null `radial` screens into 95 % upper limits on the surface density of
  W1 negative-mass lenses (n = 1, ε < 0; D-047), per lens mass. The method, tables and caveats are in
  docs/exotic_limits.md.
- **Parameter.** The mass |M| is the parameter; the lens sits at the cluster redshift. Each lensed row gets
  θ_E(|M|, z_l, z_s), with z_s its photo-z, or 2 without one. θ_E at z_s = 2 is a label only.
- **Injection.** Each lens sits at a random point of the screened footprint. The lensed sources are the field's
  own background rows:
  - β < 2 θ_E: removed (umbra);
  - 2 ≤ β ≤ 4 θ_E: replaced by their images, with PSF-deconvolved moments mapped by the signed Jacobian,
    magnitude − 2.5 log₁₀|μ|, area × |μ| and S/N × √|μ|;
  - an image pair whose isophotes overlap (separation < sum of √(area |μ| / π)) is painted as one blended row;
  - rows below S/N 5 are dropped as undetected (bookkeeping only: the screen requires S/N ≥ 10);
  - painted images are flagged extended, so a magnified compact source never becomes a spike-veto "star".
- **Screen.** It is unchanged: `exotic_screens.radial_candidates`, `radial_grid`, `anti_window_draw` and
  `radial_defaults` are shared with `cmd_radial`, and the refactor gives byte-identical SMACS output.
- **Null.** Each batch of 10 trials has its own independent 200-draw null (trials within a batch share it, so the
  binomial σ is approximate). Each trial updates only the grid blocks
  its arcs touch, which gives results identical to a full recompute (tested).
- **Recovery.** A peak with p_random < 0.05 within 2″ of the injected centre. 200 lenses per field and mass.
- **Limit.** 2.996 / Σ ε_f A_f, with A_f the screened footprint: 1″ grid points with a catalogue source within
  4″. Its border excess is measured and given as an upper bound, with border-corrected limits alongside.
- All of these are ASSUMPTIONs.

**Alternatives rejected.**
- One θ_E for every source (the first version): θ_E depends on z_s, so the mass is the physical parameter.
- Painting the two images of a source as two rows even when they overlap: near the caustic, that invents two
  converging lines the catalogue would show as one blend.
- One fixed set of 200 null draws for all trials: it correlates every trial's p_random. Each batch now has its
  own.
- Re-running `cmd_radial` in full for every injection, or copying the whole null grid 200× per trial.
- Hard-coded copies of the screen defaults and a `nanmin(S/N)` detection floor, which removes nothing.
- Synthetic sources at random β: lensing the rows actually present keeps the real density, clustering and
  photo-z.
- Catalogue shapes used as intrinsic, with no PSF term: PSF-sized images would be over-elongated.
- Counting a lens as recovered at the ≥ 3-line peak: the screen's significance is p_random, and 3 lines is below
  every field's null (5–8 needed).
- Headline limits from all eight fields: without photo-z (MACS0717, Abell S1063) members and foreground galaxies
  get lensed and the efficiency is biased high; they enter only the optimistic set.

**Evidence** (`derived` from `simulated` injections; 8 fields, 51.2 arcmin²).
- The base screens reproduce the field docs: SMACS 31 arcs, 4 lines, p 0.965; MACS0416 120, 7, 0.225; Abell 370
  100, 5, 0.495; and so on.
  - Abell 2744 now gives 134 arcs, 5 lines, p 0.505: its doc predates the D-034 spike veto (42 segments).
- Recovered of 1,600 lenses per mass (200 in each field):

  | \|M\| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
  |---|---|---|---|---|---|
  | θ_E(z_s = 2) | 0.21–0.36″ | 0.65–1.14″ | 2.05–3.59″ | 4.11–7.18″ | 6.49–11.35″ |
  | recovered | 0 | 0 | 8 | 84 | 156 |

  - Headline 95 % limits, from the six fields with photo-z (38.0 arcmin²): none at 2 × 10¹⁰ and 2 × 10¹¹ M☉,
    < 6.1 × 10⁴ deg⁻² at 2 × 10¹², < 7.0 × 10³ at 8 × 10¹² and < 4.0 × 10³ at 2 × 10¹³. MACS0717 and Abell S1063
    have no photo-z, so their members and foreground galaxies get painted as W1 images and their efficiency is biased
    high; all eight fields (optimistic): 3.7 × 10⁴, 3.8 × 10³ and 2.1 × 10³ (border-corrected 4.1 × 10⁴, 4.2 × 10³,
    2.4 × 10³; headline border-corrected 6.6 × 10⁴, 7.6 × 10³, 4.3 × 10³).
  - The best headline limit is about 30× weaker than Takahashi & Asada's volume limit spread over 0 < z < 1 (about
    120 deg⁻²).
- History: #70 merged a first review round (photo-z-only headline, θ_E-parametrised); this record supersedes its
  numbers with the mass-parametrised, blend-aware, independent-null run.
- The loss is in the arc selection. An image reaches e ≥ 0.5 only for β ≲ 2.3 θ_E, and `anti` relative to the
  cluster keeps about a third of the images. A 2 × 10¹² M☉ lens therefore puts 0.5–2.7 arcs into the screen,
  while the null needs 5–8.
- Tests: `tests/test_inject_radial.py`, offline. They cover:
  - the windowed and incremental null against `line_counts` on the full grid;
  - θ_E(z_s, |M|) against D-047;
  - blends, the detection floor, painting, recovery of a dense synthetic lens, and the footprint border.

**Revisit if.**
- A W1-specific screen is built: collinear radial image pairs flanking an empty centre, with orientation measured
  relative to the candidate centre instead of the cluster, and a local rather than field-maximum null. This
  injection harness is its benchmark.
- The screen's thresholds change (e ≥ 0.5, 60° `anti` window, 15″ lines).
- Pixel-level injections (painted into cutouts and re-extracted) are needed to measure blending and
  incompleteness.

## D-050 W1 negative-tangential-shear screen: in-house catalogue aperture-mass map on scipy cKDTree; TreeCorr and lenspack rejected (2026-10-08)

**Decision.** Build the W1-specific screen (D-049 "Revisit if") as a catalogue aperture-mass map in
`scripts/exotic_screens.py` next to `radial_candidates`, with existing dependencies only (numpy, scipy
`cKDTree`). Reuse-check result; not implemented yet.
- **Statistic.** Schneider (1996) catalogue estimator on a grid of candidate centres (1″, the D-049 footprint):
  M_ap(θ₀) = Σ Q(|θ_i − θ₀|/R) e_t,i / Σ Q, with a Schirmer et al. (2007) shear-shaped filter and a top-hat
  option. Report −M_ap, because a W1 lens (ε < 0) gives *negative* tangential shear (radial alignment); report
  M_× as the B-mode/systematics check. Take Q's exact formula from the Schirmer et al. full text (only the
  abstract was checked).
- **Input.** Background rows as in D-049 (photo-z fields only), PSF-deconvolved second moments as in
  `inject_radial.lensed_shapes`, and the cluster model's reduced shear removed: e_int = (e − g)/(1 − g* e).
  Mask where |g| ≳ 0.5 (ASSUMPTION) instead of subtracting.
- **Null.** Random position-angle rotations with positions kept, ≥ 200 draws per field from cached neighbour
  lists; a local p-value per centre and a field-maximum p-value.
- **Adoption.** Benchmark with `scripts/inject_radial.py`; adopt only if it beats `radial`'s D-049 efficiencies.
- **Expectation (`derived`, rough).** γ_t ≈ (θ_E/θ)², so at θ_E = 1″ γ_t ≈ 0.11 at 3″; with about
  0.05–0.09 lensable rows per arcsec² (D-049) and σ_e ≈ 0.3, S/N ≈ 1–2 in a 5″ aperture. Gains are expected
  mainly at θ_E ≥ 2–3″; the injections decide.

**Alternatives rejected.**
- TreeCorr `NGCorrelation` / `calculateNMap` (5.1.4, BSD-3): ⟨N M_ap⟩ is stacked over all lens positions and
  returned against R only, with no per-centre map or local null. One correlation per grid centre would be slower
  than one cKDTree pass. No Windows wheels on PyPI or conda-forge (the owner's host would need an MSVC build).
- lenspack (1.0.0, 2020, MIT): `aperture_mass` filters a pixelised (binned, KS93) map, which loses the arcsecond
  scales W1 needs and adds edge/mask artefacts; its `gamma_tx` and `random_rotation` are ~10 lines each; no
  release since 2020.
- Weak-lensing peak finders on pixelised maps: same limitation.
- Keeping `radial`'s arc selection (e ≥ 0.5, `anti`): D-049 shows it discards most W1 images.

**Evidence.** TreeCorr NG docs (https://rmjarvis.github.io/TreeCorr/_build/html/ng.html) and its PyPI/conda-forge
file lists; lenspack source (https://github.com/CosmoStat/lenspack); Schneider 1996; Schirmer et al. 2007
(SOURCES "W1 shear screen").

**Revisit if.**
- Wide fields (≥ 1e5 sources) where tree-code speed matters and TreeCorr ships Windows wheels, or a stacked
  ⟨N M_ap⟩ around a list of `radial`/orphan-pair centres is wanted (exactly TreeCorr NG).
- A maintained catalogue-level aperture-mass map package with a per-centre null appears.

## D-051 Orphan-pair screen in deep fields, with injection-recovery: null; W2/W1/point-mass surface-density limits (2026-10-08)

**Decision.**
- **Catalogue adapter.** `scripts/orphan_pairs.py` reads every catalogue into one column layout (`as_standard`),
  so CANUCS DR1 and DJA grizli (catalogue plus eazy-py zout) both run through D-048's code. The three D-048
  cluster runs reproduce every pair and count bit for bit; only the column `mu_canucs` is renamed `mu_cat`.
- **Deep fields** (`DEEP_FIELDS`):
  - the five CANUCS NIRCam flanking fields (NCF: MACS0416, MACS1149, Abell 370, MACS0417, MACS1423; 23–38 MB
    each, same format as the cluster catalogues);
  - DJA v7.3 GOODS-North.

  They use these rules (ASSUMPTIONs):
  - no published images and no cluster model; the ordinary-lensing test is the catalogue |μ| (≤ 1.4, or 1);
  - redshift floor z_low > 0.5 (`Z_LENS_REF` 0.4 + 0.1);
  - S/N bands are those of F277W/F356W/F444W that exist (two NCFs lack F356W);
  - DJA: 0.36″ apertures, SEP flags, MIRI and `u` duplicates dropped, and a same_galaxy radius of
    3.3 × `flux_radius`.
- **Two new nulls.** (d) is the z-overlap match rate at 3–6″. (e) is null (c) conditioned on fainter-member S/N,
  larger-member size and LW/SW colour. Both are new keys; the D-048 keys are unchanged.
- **Injection-recovery.** `scripts/inject_pairs.py` paints simulated pairs into the real catalogue
  (docs/exotic_limits.md "W2"). Lensed rows are removed and replaced by |μ|-scaled copies of their SED (scatter max(|μ|, 1) ×
  the row's errors). Two images closer than the catalogue's 2nd-percentile nearest-neighbour separation (0.24–0.30″) merge
  into one row. The unchanged orphan rules then decide recovery.
  - Lens types: point mass (sanity), W2 Ellis (n = 2) and W1 negative mass (n = 1, ε < 0, β ∈ [2, 4] θ_E).
  - θ_E ∈ {0.15, 0.3, 0.7, 1.5}″.
  - Per-source (400) and per-deflector (2000) trials per field, type and θ_E.
- **Limits.** Two 95 % limits on the surface density, each over Σ ε_f A_f (108.1 arcmin² searched):
  - **no-candidate:** 2.996 / Σ ε_f A_f;
  - **background-aware:** s₉₅ = 72.3 (CLs) / Σ ε_f A_f, for 355 orphans observed against 315.4 expected.
- **Stated reason for the > 200 MB downloads** (CLAUDE.md):
  - **Catalogue, 223.8 MB.** The GOODS-N DJA catalogue is the only file with the matched-aperture fluxes, and
    FITS tables are row-major, so a byte-range read cannot skip columns. It is pinned with `max_bytes` 230 MB.
  - **Photo-z tarball, 371.1 MB.** It is streamed: `fetch_tar_member` hashes the whole archive in flight and keeps
    only the 67.6 MB zout, so nothing over 200 MB is written to disk.
  - Every other DJA deep field costs the same or more (SOURCES.md).
  - The CANUCS NCF catalogues give five fields with PSF-matched photometry for 160 MB in total, which is why they
    carry most of the area.

**Alternatives rejected.**
- **The standalone `gdn-grizli-v7.3-fix.eazypy.zout.fits` (60.6 MB).** It belongs to an older catalogue: 63,069
  rows, ids offset by a median 137″. It cannot be joined.
- **DJA Kron apertures for the same_galaxy rule.** `kron_radius` lies between 2.4 and 3.8, which gives
  8.7 × r₅₀ against 3.3 × in CANUCS. 0 of 2000 point lenses at θ_E = 0.3″ were recovered in GOODS-N, and the
  orphan count fell from 109 to 69.
- **A global null (c) only.** Flanking-field orphans exceed it at P = 0.04–0.08 per field. Conditioning on S/N,
  size and colour (null e) absorbs part of the excess, not all of it.
- **Re-running the whole field search per injected lens.** The rules are local (≤ 1.8″ from the pair), so the
  neighbourhood within β_max θ_E + 6″ gives the same classification at about 8 ms per trial.
- **Drawing synthetic sources.** Lensing the real rows keeps the real density, SEDs, depth and clustering (as in
  D-049).
- **CEERS, GOODS-S, PRIMER.** They are not run: 250–270 MB catalogues plus 350–410 MB tarballs each, for fields
  that add area but no new method.

**Evidence** (`derived`; summaries in `outputs/orphan_pairs/<field>/` and `outputs/inject_pairs/`; all six
contact sheets inspected).
- Orphans observed against null (e):

  | Field | Orphans / null (e) | P(≥ observed) |
  |---|---|---|
  | M0416-NCF | 51 / 41.4 | 0.083 |
  | M1149-NCF | 50 / 40.6 | 0.085 |
  | A370-NCF | 43 / 43.1 | 0.53 |
  | M0417-NCF | 49 / 38.0 | 0.049 |
  | M1423-NCF | 53 / 47.9 | 0.25 |
  | GOODS-N | 109 / 104.4 | 0.34 |
  | Total | 355 / 315.4 | 0.015 |

- **Vetting.** Of the 90 top orphans inspected, 50 are knots, companions, group members, satellites or
  artefact-affected. The other 40 are featureless faint pairs at θ_E-equivalent 0.4–1.4″ (`hypothesis`:
  4 × 10¹⁰–7 × 10¹¹ M☉ at z_l = 0.4). No pair goes to `/vet-candidate`.
- **The ~15 % excess is not read as lensing.** At the measured efficiency (~0.5 % per deflector), it would need
  ~140 dark galaxy-mass deflectors per arcmin², one per ~7 catalogued galaxies.
- **Per-deflector efficiency** (mean over 108 arcmin²):

  | Type | 0.15″ | 0.3″ | 0.7″ | 1.5″ |
  |---|---|---|---|---|
  | point | 0.02 % | 0.54 % | 0.84 % | 0 (pair separation > 3″) |
  | W2 | 0.01 % | 0.64 % | 0.83 % | 0 (pair separation > 3″) |
  | W1 | 0 | 0.10 % | 3.6 % | 18.6 % |

  - Per source, recovery is ≤ 3.5 %.
  - Only 1–12 % of lensed selected sources keep both images above the S/N cuts.
  - same_galaxy removes the 0.3–0.6″ pairs, and visible_lens removes ~73 % of the 0.7″ pairs.
  - At best, recovery reaches 11–15 % at m < 24.
- **Limits (deg⁻², no-candidate / background-aware).**
  - W1: < 5.4 × 10² / 1.3 × 10⁴ at θ_E 1.5″ (|M| 4.5 × 10¹¹ M☉) and < 2.8 × 10³ / 6.7 × 10⁴ at 0.7″
    (9.8 × 10¹⁰ M☉).
  - W2: < 1.2 × 10⁴ / 2.9 × 10⁵ at 0.7″ (a ≈ 10 pc). That is ~100× weaker than Takahashi & Asada's volume limit
    spread over z < 1 (~120 deg⁻²).
  - Point mass: < 1.2–1.9 × 10⁴ / 2.9–4.5 × 10⁵ at 0.3–0.7″.
  - Masses and throat radii are a `model_prediction` at z_l = 0.4 and z_s = 2 (Planck18).
- **Tests.** `tests/test_orphan_pairs.py` covers the adapter: CANUCS and DJA layouts give the same pairs; flags,
  misaligned zout, missing F356W, the deep-field lens check, tar streaming, footprints and polygons.
  `tests/test_inject_pairs.py` covers mass and throat inversion of `einstein_radius`, Poisson limits, recovery of
  a bright point-lens pair as an orphan, the merge rule, the W1 umbra and image side, and efficiency bounds.

**Revisit if.**
- A companion-aware null exists: spectroscopic close pairs, or a null matched in environment. The flanking-field
  orphan excess (P = 0.01) must be explained before any orphan is read as anything but chance.
- Pixel-level injections (painting into mosaics and re-extracting) are run. They would measure blending with
  neighbours and photo-z scatter of faint counter-images.
- The same_galaxy or visible_lens rule changes, for example to segmentation-map adjacency. Those rules set most of
  the efficiency loss.
- More deep fields are added (CEERS, GOODS-S, PRIMER), or a real-pair spectroscopic sample tests the SED-match
  power.

## D-052 W3 multi-epoch dimming / vanished-source screen: null in NEXUS, MACS0416 and Abell 2744; injection-calibrated rate limits (2026-10-08)

**Decision.** `scripts/dimming_screen.py` screens multi-epoch level-3 NIRCam catalogues for W3 (D-047). Results are in
docs/exotic_limits.md "W3 inverted microlensing / dimming (multi-epoch)"; epochs are in `configs/dimming_screen.yaml`.
- **Reuse.**
  - `epoch_compare.mutual_matches`, the D-027 global frame shift and `transient_search._neighbour_counts`.
  - `transient_combine.exclusion_radius` and `fetch_gaia`.
  - `transient_forced.measure`, `select_controls` and `robust_std`.
  - `cutouts.make_cutouts` (D-018/D-021) and `exotic_sim.inject_light_curve`.
- **Flags** (ASSUMPTIONs): `vanish`, achromatic in every testable band; `dim_achromatic`; `rise_dip_rise`.
- **Order of tests:** catalogue vetoes first, then forced confirmation, cutout tests, visual check and SIMBAD/NED.
  The catalogue vetoes are: Gaia star with the self-match excluded, a ≥ 100× brighter neighbour within 1.5″, edge,
  blend, single epoch, and sharper than the PSF.
- **Injection model.**
  - f = F f_obs + √max(0, 1 − F²) σ_n z, with error σ_n max(F, 1).
  - Light-curve vetoes are re-applied per injected copy.
  - Only baseline-unflagged sources count in the exposure.
  - The forced stage is emulated by the forced/catalogue fractional scatter measured on the same controls.
- **Calibration and limits.**
  - Headline limits come only from fields whose forced stage is calibrated (≥ 50 controls in some epoch pair and
    control zero points in every epoch). The others are reported separately.
  - Limit definitions: rate per source per year = 3 / Σ N ε (T + 4t_E); τ = rate × π t_E; per deg² per year and
    per deg² per epoch analogously (docs/exotic_limits.md).

**Alternatives rejected.**
- Image differencing, or forced photometry of whole mosaics: it needs full `_i2d` files (NEXUS o014 F200W is
  113 GB).
- Per-band vanish flags: they count single-band deblending misses.
- Catalogue non-detections as zero flux at any depth: MACS0416 gave 1,122 vanishes from shallower PEARLS epochs.
- `is_extended == False` as a point-source cut: 70 % of S/N ≥ 10 rows qualify.
- `transient_combine.near_bright` for the star mask: it returns the source's own faint Gaia match as the nearest
  star, so a saturated star 0.5–1″ away never masks. That let two MACS0416 spike artefacts through.
- `f_obs + (A − 1) f_ref` and `F f_obs + (1 − F) n`: the first leaves the bright-star epoch scatter, and the
  second's scatter (√(F² + (1 − F)²) σ) disagrees with the quoted error.
- Dividing the forced ERR scale by the catalogue scale for the inflation: the two are measured against different
  baselines.
- Positive-flux-only control zero points: they are biased against faded epochs.
- Controls drawn from the largest catalogue: NEXUS o014 barely overlaps the other epochs.
- The D-047 spike amplitudes ×7.5 / ×3.4, and the pre-merge quadrature (×9.2 / ×4.9): the merged simulator gives
  ×7.0 / ×2.35 / ×1.53.

**Evidence** (derived; `simulated` for injections).
- Catalogue flags: 3,793 / 854 / 530 (NEXUS / MACS0416 / Abell 2744).
- After catalogue tests (current code): 1,051 / 358 / 32.
- Last complete forced run: 100 / 375 / 32 measured, 0 / 6 / 0 confirmed. Every confirmed MACS0416 flag is a
  saturated star or lies within 1.5″ of one; the new `bright_neighbour` veto removes all six automatically
  (cutouts inspected). **0 surviving events.**
- Re-runs after the review fixes failed on S3 (s3fs "bucket does not exist" through the proxy, 2026-10-08 07:00 and
  08:29 UTC, for multi-target jobs; single reads worked). The last complete run was re-calibrated offline.
- MACS0416 is calibrated: noise 2.15 / 1.99, forced/catalogue scatter 1.00 / 1.98. NEXUS and Abell 2744 are not.
- Efficiency, MACS0416: W3 0.12–0.22 at t_E = 0.3–3 yr; dimming 0.09 / 0.59 / 0.35 for 20 / 50 / 100 %. The 100 %
  case is low because the per-copy `single_epoch` veto removes full vanishes of sources seen in only two F200W
  epochs.
- Headline 95 % limits (MACS0416, 34 sources): 0.16 per source per year at t_E = 0.3 yr, 0.056 at 1 yr; τ < 0.15
  at 0.3 yr. All fields, indicative: 0.025 and 0.015.

**Revisit if.**
- S3 cutout jobs work again: re-run `forced` for all fields with ≥ 300 overlapping controls to calibrate NEXUS and
  Abell 2744.
- The D-039 persistence test replaces the `single_epoch` veto for vanish flags. That restores sensitivity to full
  vanishes of two-epoch sources.
- JADES or more NEXUS epochs enlarge the compact sample.
- A W3 model with cluster macro-magnification is needed for caustic-crossing stars in arcs.

**Addendum 2026-10-19 (Abell 2744 calibrated).**
- S3 byte-range jobs work again, so `forced` was re-run on Abell 2744 with 300 controls, code unchanged
  (`results/dimming/abell2744_calib/RESULT.md`, plan frozen first).
- The gate passes on the o002–o006 pair only (247 / 256 controls; o001 pairs 1–39). Noise scale 4.05 / 7.21 and
  forced/catalogue inflation 1.00 / 3.60 (F200W / F444W), both higher than MACS0416's; not explained.
- 32 flags measured, 0 confirmed, 0 controls flagged: still null. Calibrated W3 efficiency 0.000–0.013, about
  4–7× below the 2026-10-08 uncalibrated values.
- Abell 2744 joins the headline: per-source rate limits tighten by 5–11 % (0.053 per source per yr at t_E = 1 yr,
  ρ = 0.1). NEXUS is still uncalibrated and is the field that could move the headline.

**Addendum 2026-10-19b (NEXUS calibrated).**
- `forced --workers N` (concurrent image reads) reproduces the Abell 2744 calibration exactly and cuts its wall time
  666 → 236 s; NEXUS `forced` (1,051 flags, 300 controls) took 1,756 s (`results/dimming/nexus_calib/RESULT.md`).
- The gate passes: noise 1.91 / 2.57, inflation 1.00 / 2.06 (F200W / F444W). 32 forced-confirmed, 5 survive the
  artefact tests, 10 / 300 controls flagged. By eye (post hoc, `/vet-candidate` pending): three single-band
  single-epoch F200W points, one extended-nucleus F200W-only variation, one SN-like brightening; none a W3 dimming.
- Calibrated W3 efficiency 0.026–0.218 (ρ = 0.1). If the survivors hold, the three-field headline tightens ×2.6–6.2:
  0.020 per source per yr at t_E = 1 yr, 0.029 at 0.3 yr (`headline.json`).

**Addendum 2026-10-19c (NEXUS survivors vetted).**
- `scripts/dimming_exposure_check.py` measures each flagged point on every level-2 `_cal` exposure that covers it
  (`results/dimming/nexus_vet/RESULT.md`).
- The three single-band points are single-exposure, JUMP-flagged events. Two of them are covered by one exposure
  only, and those two fail a post hoc chromatic veto: simultaneous F444W below F200W / 20 at 3σ. A lensing spike is
  achromatic, so the veto costs no efficiency.
- d73390's F200W variation is not reproduced on the `_cal` frames (epoch min / max 0.865). d79956 is a brightening
  (likely an ordinary transient; TNS not reachable).
- 0 W3 events: the ×2.6–6.2 three-field headline stands.
- Adopted: `chromatic_veto` is a cheap first test for single-band forced points in future `forced` runs. It is
  not yet in `dimming_screen.py`.

**Addendum 2026-10-19d (MACS0416 regenerated; joint `combine`).**
- MACS0416 `screen` / `forced --workers 16` / `inject` regenerated with unchanged code; screen counts identical,
  forced stage calibrated, 0 survivors; exposures within −3.7 % to +0.7 % of the D-052 values (frozen ±10 % gate
  passed; `results/dimming/macs0416_regen/`).
- The joint three-field `combine` is now the headline: < 0.44 … 0.016 per source per yr and < 3.0e3 … 118 per deg²
  per yr at t_E = 0.01 … 3 yr (ρ = 0.1; docs/exotic_limits.md).
- Catalogues read from S3 `cloud_uri` when the MAST API is unreachable, verified by the manifest's sha256.

**Addendum 2026-10-19e (chromatic veto in `forced`).**
- `dimming_screen.py forced` now writes `chromatic_veto` per source (`chromatic_vetoed`: some epoch meets
  `dimming_exposure_check.chromatic_veto` on its F200W / F444W forced fluxes) and removes those sources from
  `survives_artefact_tests`; `forced_summary.json` counts them. Fields without both bands are not tested.
- Efficiencies and limits are unchanged without a re-run: `inject` does not run `forced` (it emulates it by error
  inflation), and an injected W3 event scales both bands by the same factor, so its F200W / F444W ratio is the
  host's (unit test). Existing `forced.ecsv` outputs predate the column; the NEXUS reading of addendum c stands.

**Addendum 2026-10-19f (persistence-test swap scoped).**
- Replacing the `single_epoch` veto by the D-039 test is worth building
  (`results/dimming/persistence_scope/RESULT.md`, rule frozen first): veto off, the joint ρ = 0.1 exposure grows
  ×2.0–3.4 by t_E (simulated upper bound), and only 37 compact single-epoch vanish flags need the test.
- The swap waits for the persistence test's `on_sky` efficiency against mosaic S/N (measured on monitored
  sources); the injected copies take that efficiency, not 1.

**Addendum 2026-10-19g (persistence test run, partial).**
- On the 17 NEXUS compact single-epoch vanish flags the D-039 test gives 12 `on_sky`, 4 `inconclusive`, 1
  `undetected`, 0 `persistence` (`results/dimming/persistence_run/RESULT.md`). `on_sky` flags go to `forced` and
  `/vet-candidate`; the others stay vetoed.
- Partial controls: 61 / 61 monitored pairs `on_sky` (S/N 10–50). The swap waits for the full control set.

**Addendum 2026-10-19h (weight-map coverage term adopted).**
- `screen` and `inject` uncover proxy-covered epochs whose F200W i2d `WHT` (32-pixel grid, 3 × 3 node minimum) is
  below 0.5 × the median (`Params.wht_rel_min`, `wht_step_px`; ASSUMPTIONs), the rule gate V2 validated
  (`results/dimming/edge_coverage/RESULT.md`); `--no-weight-coverage` restores the catalogue-count proxy alone.
- Production reproduces the validated what-if exactly; the corrected three-field headline
  (`results/dimming/wht_coverage/limits_combined.ecsv`) replaces `macs0416_regen/limits_combined.ecsv`: ×1.07–1.33
  weaker per source at t_E ≥ 0.03 yr, ×0.82–0.92 at 0.01 yr (`results/dimming/wht_coverage/RESULT.md`).
- Alternative rejected: the per-stamp `stamp_edge` rule (post hoc in gate V2; needs `forced` cutouts of every
  monitored pair). Revisit if a finer grid or the stamp rule is validated on a new field.

**Addendum 2026-10-19i (persistence-test efficiency on the full control set, NEXUS partial).**
- Monitored pairs 438 (MACS0416 121, Abell 2744 127, NEXUS 190 from e0, e1, e3, e4): `persistence` 0 / 438. Rule 1
  passes in MACS0416 and Abell 2744 (247 / 248 `on_sky`) and fails in NEXUS (129 / 174 at S/N 20–50)
  (`results/dimming/persistence_run/RESULT.md` "Control completion").
- Post hoc: the NEXUS loss is set by exposures per source (≤ 4: 49 / 96 `on_sky`; > 4: 94 / 94), not mosaic S/N, so
  the swap will take a per-(source, epoch) efficiency by exposure count, to be frozen and checked on NEXUS e5–e7
  first. Alternative rejected: one pooled mosaic-S/N curve (passes rule 1 only because MACS0416 / Abell 2744
  dominate). Revisit if the efficiency model fails on e5–e7.

**Addendum 2026-10-20j (exposure-count efficiency checked on held-out NEXUS e5–e7).**
- Frozen model (`results/dimming/expcount_eff/PLAN.md`): swap efficiency 0.412 for `n_exposures` ≤ 4 and 0.961 for
  > 4 (training Wilson 95 % lower bounds). Held out e5–e7: ≤ 4 33 / 89 `on_sky` (Wilson 0.278–0.475), so the rule
  holds as written; > 4 untested (no pairs); `persistence` 0 / 89 (`results/dimming/expcount_eff/RESULT.md`).
- Post hoc: the held-out rate is lower (0.37); quote the swap with the pooled ≤ 4 lower bound 0.374 as a sensitivity
  case. Alternative rejected: a mosaic-S/N curve (S/N 10–33 throughout, cannot separate). Revisit if the 11
  held-out `undetected` pairs turn out to be a position effect, or NEXUS e2 controls fail the > 4 class.

**Addendum 2026-10-20k (persistence swap adopted in the injection chain).**
- `inject --persistence-verdicts` swaps the per-copy `single_epoch` veto for the D-039 test: a copy detected in one
  epoch that is a `vanish` flag counts ε by covering-exposure count (addendum j; pairs without a committed count
  count 0). The corrected three-field headline is `results/dimming/swap_inject/runs/limits_combined_swap.ecsv`
  (×1.54–2.51 joint exposure at ρ = 0.1 over `wht_coverage/limits_combined.ecsv`; gates G0 / G1 pass;
  `results/dimming/swap_inject/RESULT.md`). Every real flag the swap clears was tested and explained.
- Alternative rejected: a per-copy Bernoulli draw (adds noise; the expected weight gives the same mean).
  Revisit if NEXUS e2 controls or the 11 `undetected` stamps change ε, or a field's swap-cleared flags survive
  `forced`.

**Addendum 2026-10-20l (NEXUS e2 controls; k = 1 aperture not tested).**
- ε(> 4) = 0.961 holds on held-out NEXUS e2: 143 / 143 `on_sky` (Wilson 0.974–1.000), `persistence` 0 / 151
  (`results/dimming/nan_k1/RESULT.md`). The k = 1 NaN fill (G1 passed in `nan_aperture`) stays unadopted: e2 has
  8 ≤ 4 pairs (frozen minimum 10), all `on_sky` with or without it, so ε(≤ 4) stays 0.412. Revisit if a new wide
  NEXUS or JADES epoch gives ≥ 10 unseen ≤ 4 pairs (frozen value 0.584 in `results/dimming/nan_k1/PLAN.md`).
**Addendum 2026-10-20m (o015 Tile-6 checked).** The only unused public 5105 observation, `jw05105-o015_t005`,
  holds none of the 151 monitored sources (`results/dimming/o015_k1/RESULT.md`): no unseen ≤ 4 sample; k = 0 and
  0.412 stay. The revisit condition above is unchanged (a new Center wide epoch or JADES).

## D-053 W1 shear screen built and adopted for W1 limits at ≥ 2 × 10¹² M☉; four cluster fields null (2026-10-08)

**Decision.**
- D-050's catalogue aperture-mass map is `exotic_screens.py shear` (Schirmer Q_TANH, R_ap = 10″, 1″ grid, sources at
  1–10″, rotation null) with `scripts/inject_shear.py` for injection-recovery. It replaces `radial` for W1 limits at
  ≥ 2 × 10¹² M☉ (5–8× stronger); `radial` stays as an independent screen. Method and tables: docs/exotic_limits.md
  "W1 negative-mass lenses (shear screen)".
- **Responsivity.** The cluster shear removed is R g with R fitted per field (0.41–0.48): catalogue isophotal
  moments respond to shear by R, not 1. Injected images keep R of the lens-induced change of their measured moments (image minus source, so the cluster shear cancels), and painted images face the spike veto; only resolved sources are painted (unresolved ones have no measured shape), and not those with κ ≥ 1, |g| ≥ 1 or |R g| ≥ 1 or spike segments. Fewer than 20
  calibrating rows, or an R outside (0, 1.5] or below 3σ, is an error, not R = 1.
- **Spike veto.** Diffraction-spike segments (`spike_segments`, D-043, Gaia stars where the field uses them) are
  dropped: they point radially at their star, the W1 sign.
- **E/B rule.** A peak counts only if p_random < 0.05 against the rotation null *and* S exceeds the field's largest
  |S_×|. Adopted after seeing real B-mode extremes beyond the rotation null (conservative; post hoc).
- All thresholds are ASSUMPTIONs.

**Alternatives rejected.**
- Removing the full model g: leaves −(1 − R) g ≈ −0.55 g, a radial pattern of the W1 sign around every mass
  concentration (measured slope of ε along g: 0.41–0.48 in four fields).
- No spike veto (first run): Abell 370 E max 4.51 and B max 4.68 (both p < 0.005), Abell 2744 E 3.92 (p 0.01). With
  the veto: 3.48 / 3.46 and 3.39. Spikes, not lensing.
- The rotation null alone: B-mode extremes still reach p 0.005–0.04 in three fields after the veto.
- Point-mass 1/x² (4″, 10″) and top-hat 6″ filters: lower injection efficiency than Schirmer 10″ in MACS0416 and
  Abell 2744 (50 injections each; numbers in the doc).
- Hetterscheidt et al. (2005) as the source of the filter's cut-off: the exponential box E(x) is Schirmer et al.'s
  own (eq. 16), read from the full text.

**Evidence** (`derived`; four photo-z fields, 30.2 arcmin²).
- Unit tests: phases and PSF deconvolution, (ε − g)/(1 − g*ε) inverts the lens mapping, sign (radial ring S = +√2n,
  tangential −√2n), a strong synthetic ring beats the rotation null, R recovered from diluted shear, a massive
  W1 injection recovered on a synthetic field.
- Real data: ε along the model g rises with |g| in MACS0416 and Abell 2744 (e.g. +0.059 ± 0.010 at ⟨|g|⟩ = 0.14,
  +0.255 ± 0.024 at 0.42), so the screen sees the cluster's real shear; after removing R g the residual is
  consistent with 0.
- Real fields: S_max 3.39 / 3.79 / 3.40 / 3.48 (Abell 2744, MACS0416, MACS1149, Abell 370), p_rot 0.58 / 0.050 /
  0.19 / 0.20; only Abell 370 exceeds its max |S_×| (3.46), and not the rotation null: null.
- Recovered (of 800, four fields): 0, 1, 33, 197, 405 at 2 × 10¹⁰, 2 × 10¹¹, 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉
  (radial: 0, 0, 8, 84, 156 of 1,600). 95 % limits: 7.7 × 10³, 1.3 × 10³, 6.5 × 10² deg⁻² at the top three masses
  (radial headline 6.1 × 10⁴, 7.0 × 10³, 4.0 × 10³).
- Wall time 58–168 s per field (1,000 injections).

**Revisit if.**
- A proper weak-lensing shape catalogue (PSF-anisotropy-corrected, calibrated) exists for these fields: the E/B
  floor would drop and the limits improve.
- DJA photo-z for SMACS 0723 and El Gordo are reachable again (tarballs 404 on 2026-10-08).

**Addendum 2026-10-17 (Abell S1063 added).** The screen runs on Abell S1063 (CATS v4.1, DJA v7.5 photo-z, Gaia spike
veto): null (S_max 3.61, p_rot 0.07; max |S_×| 3.33). Shear-only inputs live in `inject_shear.SHEAR_FIELD_OVERRIDES`
so the radial FIELDS entry (no photo-z, D-049 headline) stays as it was. Five-field limits 5.8 × 10³ / 1.1 × 10³ /
5.7 × 10² deg⁻²; the four earlier fields reproduce exactly. Evidence: `results/w1_shear_s1063/RESULT.md`.

**Addendum 2026-10-17 (SMACS 0723 and El Gordo added).** The "Revisit if" on DJA photo-z holds: both tarballs answer
HTTP 200 and match SOURCES.md. The screen runs on both (Gaia spike veto at the model centre; D-053 defaults): null
(S_max 2.91 / 2.95, p_rot 0.77 / 0.62, both below the rotation-null median). R 0.438 / 0.486. Seven-field limits
5.3 × 10³ / 9.6 × 10² / 5.0 × 10² deg⁻². `--combine-only` reads screen options missing from older summaries as
the values those runs used (an Abell 2744 re-run reproduces its committed summary exactly). Evidence:
`results/w1_shear_smacs_elgordo/RESULT.md`.

**Addendum 2026-10-17 (aperture-coverage cut not adopted).** Dropping grid centres whose aperture annulus is
< 50 % inside the catalogue footprint (`--min-coverage 0.5`, default 0) was tested on the five fields under a
frozen rule (adopt only if the combined limit is lower at each of 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉): it fails
(5,765 → 5,672 / 1,082 → 1,088 / 567 → 570 deg⁻²); only MACS0416's b_max moves, and Abell 2744's 4.21 floor is at
full coverage. Revisit if a future field's floor sits at an edge centre. Evidence: `results/w1_edge_cut/RESULT.md`.

## D-054 Survey-agnostic exotic signatures; MulensModel for ordinary microlensing fits; OGLE-IV Mróz samples first for W3 (2026-10-08)

**Context.** Owner direction (2026-10-08): find observational evidence of traversable wormholes / negative-mass
objects or warp-drive spacetimes in any public dataset, not only JWST. Exotic physics stays a hypothesis; every hit
goes through `/vet-candidate`; nulls become limits; nothing is announced outside the repo without the owner.

**Decision.**
- **Signature layer** (`jwst_anomaly.signatures`): one `Signature` per D-047 signature (W1, W2, W3, W5) with its
  `exotic_sim` prediction and injector, the screens that implement it, the ordinary mimics vetting must rule out
  first, and where its limits live. A survey enters through a thin adapter satisfying `LightCurveSurvey`
  (`events`, `light_curve`, `efficiency`) or `CatalogueSurvey` (`catalogue`, `area_deg2`);
  `standard_light_curve` is the shared light-curve layout (`observed`). In-house dict registry, no plugin
  framework.
- **Ordinary microlensing fits:** MulensModel ≥ 3.12 (MIT; pulls VBMicrolensing, LGPL-3.0) as the optional extra
  `mulens`: point lens, finite source, binary lens, annual parallax, linear source/blend fluxes. The exotic models
  are not refitted in another library: `exotic_sim` magnification (n = 1, ε < 0 with the finite-source spike cap;
  n = 2 Ellis) is evaluated on the source trajectory from `MulensModel.Model.get_trajectory(times)`, so the
  geometry (and parallax) is identical, and source/blend fluxes come from the same linear solve. Comparison:
  Δχ² and BIC on identical data; one `derived` row per event and model.
- **W3 data, in order:** (1) Mróz et al. 2019 OGLE-IV bulge sample (5,836 fitted events, 48 MB calibrated I-band
  photometry, per-field detection efficiencies); (2) Mróz et al. 2020 OGLE-IV plane sample (630 events, with
  efficiencies); (3) Gaia DR3 `vari_microlensing` (363 events) with epoch photometry via `astroquery.gaia`;
  (4) KMTNet public seasons. OGLE-IV EWS seasons only after the owner decides on its terms (below).
- **W1/W2 in wide imaging:** reuse published lens catalogues before any finder: lenscat `catalog.csv` (32,838
  entries, MIT), Euclid Q1 Strong Lensing Discovery Engine CSV (CC-BY-4.0; not the 3 GB `lens.zip`), SuGOHI;
  pinned by sha256 in manifests.
- **Warp:** no imaging or photometric prediction for distant observers (D-047). One standing research task:
  monitor the literature, including warp-bubble-collapse gravitational waves (Clough, Dietrich & Khan 2024), and
  any detector band that could test it. No signature is invented.

**Alternatives rejected.**
- pyLIMA 1.9.8 (cleanest custom-model hook, but GPL-3.0-or-later and heavy dependencies); VBMicrolensing alone (no
  negative-mass or Ellis lens; MulensModel wraps it); VBBinaryLensing (superseded, GPL, win-only wheel); muLAn
  (unmaintained since 2020); the lenscat package (needs `ligo.skymap`; read its CSV); Master Lens Database (no
  response 2026-10-08; lenscat covers it); pluggy/stevedore (unneeded for < 10 in-repo modules).
- MOA alerts and ZTF/IRSA light curves for now: no verified bulk MOA endpoint; ZTF is ineffective toward the bulge.

**Evidence** (reuse-check, 2026-10-08; SOURCES "Signature layer and time-domain archives (D-054)").
- Prior art: no published archival search of OGLE, MOA, KMTNet, Gaia or MACHO light curves for negative-mass or
  Ellis events was found (arXiv API, `abs:wormhole AND abs:microlensing`, `abs:"negative mass" AND
  abs:microlensing`); existing limits come from SDSS quasar lensing (arXiv:1303.1301), BATSE (gr-qc/9805075,
  astro-ph/9904399) and femtolensing (arXiv:1302.7170). A W3 survey limit would be new: `needs-human` before any
  outside announcement.
- Selection caveat: the Mróz samples were selected with a PSPL-like finder and the EWS selects brightenings, so an
  inverted event with a long umbra may be excluded; use the published efficiencies only for events the finder
  could have selected, and state this limit.
- OGLE EWS terms ask users to contact the OGLE team before publishing (co-authorship may be required). The Mróz
  data products are published with their papers (cite them). Using EWS seasons in this public repo waits for the
  owner (`needs-human`).

**Revisit if.** MulensModel or VBMicrolensing add user-defined or negative-mass magnification; pyLIMA relicenses; a
paper publishes a survey light-curve limit on negative-mass or Ellis lenses (compare, do not repeat); OGLE or
KMTNet terms forbid a population reanalysis; > 10 signature modules or outside contributors need entry points.

## D-055 Orphan-pair null (e) made symmetric; companion-aware null (f); the D-051 deep-field excess was a null artefact (2026-10-08)

**Decision.** `orphan_pairs.pair_cells` bins colour per member and uses the unordered pair of bins; a non-finite
colour (no valid flux, or a non-positive mean flux) has its own bin. A cell needs ≥ 5 z-overlapping reference pairs
(`MIN_CELL_REF`, ASSUMPTION), else its S/N × size cell, else the global reference rate. A new null (f) applies the same cells to the 3–6″
z-overlapping pairs (null (d)'s annulus). Background-aware W1/W2/point-mass limits use null (e) as fixed.

**Alternatives rejected.**
- Pair mean colour: a red + blue pair would share a cell with two neutral members, which match far more easily.
- (f) as the limit background: it predicts more (347.5), so its limits would be less conservative.

**Evidence** (`derived`; `outputs/orphan_pairs/<field>/summary.json`, `outputs/inject_pairs/limits.json`; table
in docs/orphan_pairs.md).
- Orphans unchanged from D-048/D-051 in all nine fields. Deep fields, 355 orphans: null (e) fixed 334.0 (P = 0.13;
  D-051: 315.4, P = 0.015), null (f) 347.5 (P = 0.35). The D-048 clusters stay null (P ≥ 0.13).
- Injections re-run: no-candidate limits reproduce D-051; background-aware s₉₅ 72.3 → 54.9 (docs/exotic_limits.md).
- The symmetry, the non-finite bin and the cell floor were changed together; their shares were not measured.

**Revisit if.** A run gives P < 0.05 under both (e) and (f), or a segmentation-map same_galaxy rule changes the
orphan counts.

## D-056 W1/W2 in published lens catalogues: per-class deflector tests; weak limits; three open quasar pairs (2026-10-08)

**Decision.**
- `jwst_anomaly.lenscats` reads lenscat 1.1.3, the Euclid Q1 Discovery Engine tables and the SuGOHI list (pinned by
  sha256). It merges them by position (3″; `entries`, `refs` keep provenance) and exposes them as a
  `signatures.CatalogueSurvey` (`PublishedLensSurvey`). It holds the pure tests:
  - `brick_coverage`: footprint and depth from the brick summary;
  - `pair_images` and `quasar_pair_test`;
  - `bright_galaxy_near`;
  - SIS σ, the Faber–Jackson fit and the required lens magnitude;
  - position checks: signed J2000 name vs RA/Dec, and rounding.
- `scripts/w12_lenscats.py screen | vet`.
  - DR10 coverage and depth come from `ls_dr10.bricks_s`; Tractor boxes come from the Data Lab TAP, split at RA 0/360
    and full RA near the poles.
  - Only lensed-quasar and radio-interferometric systems are tested; galaxy-finder and sub-mm systems are
    insensitive.
  - "faint galaxy", "blended", "too close" and position or mask problems are undecided and removed from N.
  - Limits are s₉₅(k_class) / N_class, with test completeness assumed (no injection factor). Known biases of
    that assumption: a faint unrelated source in the circle removes a dark-lens system as "faint galaxy" (true
    efficiency near 0.98); a colourless pair takes a compact PSF-typed lens as an image (more "none").
- Thresholds are ASSUMPTIONs in `Params`.

**Result** (`derived`; run `summary.json`).
- Of 20,986 galaxy-scale systems, 17,555 are in the DR10 footprint: 17,102 galaxy-selected, 110 sub-mm, 325 quasar,
  18 radio.
- 29 quasar and radio systems are decided: 13 with a deflector, 16 "none". Of the 16 "none":
  - 13 have a literature lens galaxy (SIMBAD 6, published z_l 6, He et al. 2025 1);
  - 3 SuGOHI IX CHITAH pairs are open in the typical variant but explained by a lens below the LS depth in the
    conservative variant.
- Typical: f_dark < 0.48 (quasar, k = 3, N = 16), < 0.23 (radio, k = 0, N = 13).
- Conservative: < 0.50 (quasar, k = 0, N = 6); no radio limit (N = 0).
- Nothing goes to `/vet-candidate`. Earlier limits (1.5 × 10⁻⁴, then 0.13) are withdrawn (PR #81 reviews).

**Alternatives rejected.**
- Coverage from "a Tractor source within 5″": it drops the dark configuration itself, and it counted stray DECam
  detections north of the DR10 footprint.
- Counting galaxy-finder or sub-mm systems in N: they have zero sensitivity.
- Counting a galaxy fainter than required as an explanation: it is undecided.
- A recovery factor from deleting deflectors and re-running the same code: it is 1 by construction.
- Data Lab TAP uploads and `q3c_*` in ADQL (rejected by the service); per-object viewer calls.
- The Lemon lensed-quasar database (HTTP 500 on 2026-10-08) and HSC imaging (account required).

**Evidence.**
- FJ calibration on 605 galaxy-selected lenses (LS z): a = 20.41, k = 0.67, rms 0.88 mag.
- Median required m_z is 19.73 (typical) and 25.05 (conservative), against a median brick depth of 23.43.
- Catalogue defects: cluster-survey rows typed "galaxy" (14 references), AGEL declinations, SPT positions, 468
  rounded positions covered, rejected candidates kept (MJV16999), name-based merges.
- Tests: `tests/test_lenscats.py` and `tests/test_w12_lenscats.py`. They cover:
  - a dark pair reaching "none";
  - quad images never counted as the deflector, while a red PSF-typed lens is;
  - a fold-quad lens outside the brightest pair's circle (search circle holds every image; final review);
  - a deflector beyond image_radius inside the pair circle;
  - faint-galaxy, close, blended and insensitive cases;
  - empty inputs keeping the schema;
  - brick coverage independent of sources;
  - RA-wrap and pole boxes, signed declinations, rounding and Poisson limits.

**Revisit if.**
- HSC PDR or HST photometry is available for the 307 blended or close lensed quasars.
- The CHITAH lens models or spectra of the three open pairs are checked.
- A public list of rejected lensed-quasar candidates (binary or "nearly identical" quasars) appears.
- A lens list publishes image positions (W1 geometry).
- Tractor can be re-run on injected images (a measured completeness).

**Amendment (2026-10-08): the LS pair must be the catalogued pair (D-064 check); one entry per lens.**
- `lenscats.pair_match` (moved from `scripts/w12_niq.py`, tolerance `PAIR_SEP_TOL` = 0.5″, ASSUMPTION) is shared by
  both scripts. In `w12_lenscats.py`, `deflector_test` reports `used_pair`, and `pair_check` compares the LS pair
  with 2θ_E (SIS `model_prediction`; the catalogues give no image positions). A pair-based status with a mismatch
  is undecided. Without a catalogued θ_E the system stays decided and is counted as unchecked in `summary.json`.
- `dedup_same_lens` (vet) keeps one decided entry per lens. Entries match on the same designation
  (`lenscats.designation_key`, HHMM±DD after removing the prefix, J/B, spaces and suffixes) within 30″
  (`Params.same_lens_radius`, ASSUMPTION). A copy with a deflector is kept (a deflector seen at one catalogued
  position explains the lens), else the first copy. Merged ids: `summary.json` `vetting.same_lens_merged`. Three
  radio lenses are listed twice in lenscat about 11″ apart: MG0414+0534, B2114+022, B2319+052. For each, the
  "deflector" copy is kept and the "none" copy, on empty sky, is dropped.
- Re-run (`derived`): decided 29 → 25 (deflector 13 → 12, none 16 → 13). 115252+004733 (θ_E 1.67″, LS pair 4.18″)
  becomes undecided. Typical: quasar 3 / 15 / < 0.52, radio 0 / 10 / < 0.30, all 3 / 25 / < 0.31.
  Conservative: 0 / 5 / < 0.60. The three open CHITAH pairs are unchanged.
- Coverage is thin: 0 of the 15 pair-decided systems have a catalogued θ_E. A one-off check against the SQLS
  separations pinned for D-064 matched the three that have one (J1322+1052, J1349+1227, J1515+1511).
- Rejected: requiring a catalogued separation for every decided system (no pair-decided system has one, so the
  quasar class would have N = 0); importing the D-064 VizieR tables into the D-056 chain (3 of 15 matches, all
  consistent: a second pinned input set for no change); a hard-coded duplicate list; widening the 3″ catalogue merge
  (it would merge distinct close systems before any test).

## D-057 W3 in OGLE-IV Mróz samples: one fitter for ordinary and exotic models; disk sample null (2026-10-08)

**Decision.**
- `jwst_anomaly.ogle.OgleMrozSample` (a D-054 `LightCurveSurvey`) reads the published Mróz et al. 2019 (bulge,
  5,790 events) and 2020 (disk, 460 events of Table B1) products, pinned by sha256 (`data/manifests/ogle_mroz.ecsv`).
- `scripts/w3_microlensing.py fit` fits PSPL, FSPL, PSPL+parallax and the exotic `N1neg`, `E2pos`, `E2neg` models on
  the same straight trajectory with shared linear source/blend fluxes; ΔBIC = BIC(exotic) − min BIC(ordinary).
  ASSUMPTION: flag at ΔBIC < −10. `vet` then tests refits, robust errors, baseline variability, per-season offsets
  and drifts, binary source, binary lens, arXiv mentions and VSX / Gaia DR3 variable matches; a flag survives only
  if every test keeps the exotic preference.
- **Disk sample (all 460 events): null.** 6 flags (all `E2pos`, ΔBIC −12.6 … −33.8 on first fit); 0 survive. Every
  flag loses its preference once each season gets a free baseline offset and drift (`derived`).

**Alternatives rejected.** Refitting exotic models in a second library (D-054: geometry must match exactly);
treating a better exotic fit as a candidate without season-systematics tests (all 6 disk flags were systematics).

**Evidence.** CHANGELOG 2026-10-08 "W3 OGLE-IV disk sample"; contact sheet inspected (sparse peak coverage in
GD1279.14.87 and GD1081.21.615; post-peak dips below baseline in BLG979.24.9765 and BLG775.24.26593 that a
per-season drift absorbs). No rate limit yet: injections are built on bulge light curves only.

**Amendment (2026-10-08, same day): bulge sample, injections and the selection answer.**
- **Bulge sample (all 5,790 events): null.** 127 flags (E2pos 96, E2neg 22, N1neg 9; ΔBIC −10.2 … −2788); 0 survive.
  Free blend per season removes 73 → 14, season drifts → 9, binary source → 7, and two further tests, which need the
  survivors' fits and therefore run in `revet`, remove the rest: **feature coverage** (≥ 3 epochs where the exotic and
  the best ordinary model differ by > 3σ at the epochs, with the Δχ² coming from them) leaves 1; an **epoch jackknife**
  (drop up to 3 most influential epochs, keeping ≥ 3 inside the feature; review of PR #88: dropping 3 of 5 killed a
  synthetic t_E = 3 d W3 event) keeps it; **two unrelated PSPL bumps** explain it,
  BLG519.21.110304, whose N1neg spikes sit on the 2011 event and a 1-day bump in 2015 (ΔBIC +24.4 for the exotic).
- **The published samples cannot contain a W3 event.** 0 of 600 injected n = 1, ε < 0 events (t_E 3–300 d,
  ρ ∈ {0.01, 0.1}, u₀ ~ U[0, 2)) pass the emulated Mróz selection, in every t_E, ρ and u₀ bin, while PSPL controls on
  the same light curves pass 15–43 % and the fitter flags 42–97 % of the injections. They fail the one-bump, PSPL
  fit-quality, χ₃₊, three-consecutive-points and blend cuts, four at a time on average. The emulation is stricter than
  the published selection (63.9 % of the real events pass it), which cannot bridge that gap.
- **Therefore no W3 rate limit from these samples** (the measured efficiency is 0). The 95 % bound on the recovery
  fraction (3/60) says even the strongest limit they could give would be 0.8–3.9 × 10⁻⁸ per star per year, against an
  ordinary rate of 5–25 × 10⁻⁶. A W3 limit needs a search on the OGLE light curves before the PSPL selection, which
  the published products do not contain (and OGLE EWS terms are still an owner decision, D-054).
- Evidence: docs/exotic_limits.md "W3 in the OGLE-IV microlensing samples"; survivors inspected on a contact sheet
  before this entry. Wall time ≈ 5 CPU hours (bulge fit 10,413 s on 4 cores; vet 6,592 s; inject 1,974 s).
- **Params caveat:** this bulge run predates D-058, so its parallax fits are unbounded. A bound makes the ordinary
  family less flexible, so re-fitting under D-058 can only add flags, never remove one; the vetting chain and the
  injection/selection result are unaffected. The chunked tables of D-059 are the re-fit under the current Params.

**Revisit if.** A light-curve-level OGLE/KMTNet/MOA data set becomes usable (then a real W3 limit is possible, and the
selection emulation here is the baseline to beat); or the season-offset/drift test is shown to absorb injected W3
signals at small u₀ (it is applied to injections through `flag_vetted`, which stays at 43–86 %, so it does not now);
or the D-059 chunk re-fit under the bounded parallax of D-058 produces flags that pass the `revet` tests (this run,
with unbounded parallax, produced none); or the `revet` tests are run in the injection loop and remove injected W3
events (they are checked only on synthetic unit-test events now).

## D-058 W3 fitter: parallax bounded at |π_E| ≤ 5; disk sample still null (2026-10-08)

**Decision.** ASSUMPTION: `w3_microlensing.Params.pie_max = 5`; the PAR objective rejects larger |π_E| in `fit`,
`vet` and injections. Observed microlensing parallaxes are ≲ 1–2 even for nearby disk lenses, so 5 is generous.
Disk re-fit: 7 flags, 0 survivors (`derived`).

**Alternatives rejected.** Unbounded π_E (unphysical fits, π_E up to ~10³); a Gaussian prior on π_E (needs a
population model; a hard bound is enough for a flag screen).

**Evidence.** CHANGELOG 2026-10-08 "bounded parallax": simulated injections show no absorption either way
(25 / 25 flagged); the bound adds one marginal real flag (GD1217.10.8703), removed by season offsets.

**Revisit if.** Bulge injections with real cadences show the bound changes recovery, or a published event with
|π_E| > 5 appears in the samples.

## D-059 W3 bulge chunk fit tables are tracked in git (`results/w3_ogle/`) and joined by `merge-chunks` (2026-10-08)

**Decision.** `w3_microlensing.py fit --chunk K/N` (complete chunks, no `--limit`) also writes its `derived` table
as deterministic gzipped ECSV to `results/w3_ogle/fits_<sample>_chunkKofN.ecsv.gz`; `merge-chunks --n N` joins
chunks 1..N into the table `vet` reads and sets `chunk = ""` (the whole sample) only when every chunk exists,
was fitted with the current `Params` and holds exactly its own events (none skipped). `limit` keeps refusing anything else.
Raw light curves and all other outputs stay out of git.

**Alternatives rejected.** Fitting all 5,790 bulge events in one session (~4 h; cloud sessions end after ~40 min);
keeping chunk tables only under `$JWST_ANOMALY_DATA` (lost with each ephemeral session — chunk 1/12 of
2026-10-08 was lost this way and must be refitted); GitHub release assets (extra credentials and tooling for
~100 kB files); a compact column subset (the contact sheet and audits need the model parameters).

**Evidence.** A chunk table is ~480 rows × 57 columns; gzipped ~0.2 MB (CHANGELOG 2026-10-08 "chunk 2/12"),
so all 12 chunks stay ~2–3 MB, under the 1 MB per-file rule.

**Revisit if.** The tracked tables exceed ~10 MB in total, or a local session can fit the whole sample at once.

## D-060 W2 deflector test at HST resolution: the Hubble Source Catalog is not decisive (2026-10-08)

**Decision.** Do not use HSC v3 catalogue photometry to decide whether a lensed quasar lacks a deflector.
`scripts/w12_hsc_probe.py` stays as the reproducible probe and as the validation harness (known-lens efficiency)
for a pixel-level replacement.

**Alternatives rejected.** Counting HSC "none" systems toward f_dark (efficiency 0.46 on known lenses would need a
correction larger than the signal); widening the deflector radius or lowering the CI cut (the lens galaxy is
missing from the catalogue, not mis-typed, in the four inspected misses: H1413+117, HE1104−1805, SBS0909+532, HE2149−2745).

**Evidence** (`derived`, run 2026-10-08). 444 quasar/radio systems; 71 with HSC sources in ≥ 2 HSC images (91 with any); known-lens systems
13 deflector / 15 none / 38 undecided (sources in ≥ 2 HSC images; 15 / 20 / 48 without that cut); no-lens-z systems 1 none (HS0810+2554) / 4 undecided. Tests:
`tests/test_w12_hsc_probe.py`.

**Revisit if.** PSF-subtracted HST image models (or another deeper/sharper survey) reach an efficiency ≥ 0.9 on
the known-lens set; HSC v4 or a lens-aware HST catalogue appears.

**Addendum (2026-10-21, line C: PSF-subtracted HST models of the D-056 blended / too-close quasars).**
- **Decision.** The Moffat N-image residual chain (`scripts/w2_hst_blend_chain.py`) does not validate.
  - [pre-registered] It recovers 9/10 same-programme known lenses (0.9, one-sided 95 % lower bound 0.61).
  - [pre-registered] Its one negative control (a no-lens stamp with empirical-PSF images, SDSSJ0806+2006) gives S/N 13.0
    ≥ 5σ; the same stamp with Moffat images gives −0.3.
  - [owner decision 2026-10-21, after freezing] A negative control at threshold fails the method: target residuals are
    stored but not interpreted, and no candidate list is made.
  - A "not detected" from this chain is not a missing lens: [post hoc] DESJ0405-3308, "not detected", has a published
    deflector (photometric z_d ∼ 0.3, Gilman et al. 2020a, MNRAS 491, 6077, via Nierenberg et al. 2023,
    arXiv:2309.10101). "No catalogued lens z" ≠ "no lens".
- **Alternatives rejected.** Reading the "light" targets as deflectors; a post hoc flux-scaled false-positive threshold
  (stored only; a known control, WISEJ0259-1635, falls below it); a second PSF method or new thresholds inside this
  unit (owner decision); summing CAOM level-3 exposure times as depth (multiple products of the same exposures,
  over-counted 1.9–7×; fixed after /code-review, v1 kept as `*_v1`).
- **Evidence.** `results/w2_hst_blend/` (`PLAN.md` frozen first, `gate.json`, `psf_mismatch.json`, `verdict.json`,
  `selection_diff.ecsv`, `RESULT.md`).
- **Revisit if.** A per-image PSF matched to the HAP-MVM drizzle (STScI empirical PSF grids via photutils) or a joint
  lens-light + point-source model (lenstronomy) keeps no-lens negative controls below 5σ while recovering ≥ 0.9 of
  controls. More negative controls need a library PSF or star-rich fields (the Gaia ePSF found ≥ 3 stars in 1 of 5).

**Addendum (2026-10-23, J2244−0550 HST test, `results/w2_j2244_hst/RESULT.md`).**
- **Decision:** the revisit route named above, as implemented here, was tried on J2244−0550 (GO 17455 F160W, GO 17308
  F814W) and does not validate.
  - The implementation: lenstronomy 1.14.2 with two point sources, with or without a Sérsic, and PSF iteration from a
    programme-star ePSF, fitted separately for each model.
  - v3 (after two implementation fixes from the review): F814W detects 2 of 6 held-out close known lenses (≥ 80 %
    needed); F160W has 2 eligible positives (3 needed). The v2 run failed the same way.
  - Neither band has a false detection among 6 real-star negative pairs. The thresholds that achieve this (Δχ² 27203,
    10931), set from development star pairs, miss known lenses: real-star pairs reach Δχ² 13601.5 (F160W) and
    6468.4 (F814W), and 4 of 8 fitted F814W lenses at 0.93–1.60″ fall below the latter.
  - The target was not measured.
- **Rejected:**
  - raising or lowering T, or relaxing the 0.6″ centre rule, after validation;
  - a second PSF method inside the unit.
- **Untested reasons** (RESULT §7): PSF halo beyond the 41 px / 21 px kernel; two-field PSFs in the star pairs;
  non-nested M0 / M1 fits; brightness matched in Gaia G only.
- **Revisit if:** a different PSF route or kernel size (e.g. a larger kernel; STScI empirical PSF grids or forward
  modelling of the FLT / FLC exposures) keeps real-star pairs at the target's brightness below a threshold that still
  recovers ≥ 80 % of held-out ~1″ lenses. Or deeper / better-sampled imaging (JWST, AO) of J2244−0550 becomes public. Budget the
  traffic: use small fsspec blocks.

**Addendum 2 (2026-10-23, STScI PSF capability check, `results/w2_j2244_stpsf/RESULT.md`).**
- **Decision:** the self-built image-modelling line for J2244−0550 is paused (owner rule for a failed check). That
  covers the #474 lenstronomy chain and the official-PSF hook-up. No larger PSF, star library or budget is proposed.
- **Evidence:** the STScI default STDPSF model (hst_notebooks WFC3 PSF notebook §2.2.1, photutils
  `GriddedPSFModel`) was fitted on single F160W FLT exposures.
  - It hooks up: FLT positions agree across dithers to 1.3 mas (sub-pixel consistency).
  - It does not model a real G 18 star to the noise (χ²/ν 31–39; wings at 3–12 px 3.0–3.1 % above the model).
  - A no-galaxy ~1″ star pair leaves a residual of z 19 at the lens position. Its wing, core-shape and two-visit
    contributions are not separated.
- **Rejected inside the unit:** focus or spatial perturbation (hst1pass needs at least 10 bright isolated stars per
  exposure); stacked or drizzled PSFs (other methods; drizzled empirical models are not officially supported).
- **Revisit if:** the owner reopens the line, or an official PSF route models a real G ≈ 18 F160W star to the
  noise within 1.5″ (core and wings), on exposure-level or drizzled data.

## D-061 W3 in the Gaia DR3 microlensing candidates: the published selection also rejects W3 (2026-10-08)

**Decision.** Do not use `gaiadr3.vari_microlensing` (Wyrzykowski et al. 2023) to limit W3: its Sample A selection
passes almost no injected W3 event. `jwst_anomaly.gaia_mulens.GaiaDR3Microlensing` (a D-054 `LightCurveSurvey`:
TAP events, Table D.1 sample labels from the pinned arXiv source, DataLink G epoch photometry) and
`scripts/w3_gaia.py` (`fetch`, `fit`, `inject`, `summary`, `manifest`; the D-057 fitter) stay as the harness for
the next Gaia-cadence sample.

**Alternatives rejected.** Fitting the 363 candidates and quoting the flag count as a limit (the selection removes
the signal first, as in D-057); emulating the Extractor cuts with guessed definitions (they fail 126 of 163 real
Sample A events, so `selected` leaves them out; `selected_ext` keeps them for the record).

**Evidence** (`derived`, run 2026-10-08; tables in `results/w3_gaia/`, inputs in `data/manifests/gaia_dr3_mulens.ecsv`).
- Emulation audit: published cuts (score, u0, t_E, G0, skewness, amplitude, t_first, parallax χ²/dof and π_E,
  skew–Abbe with log10; the stricter first-year branch read in JD − 2450000, since Gaia time never reaches the published 6824.5–7189.5 window) pass **143 / 163** real Sample A events (149 with ln); 19 / 200 B-only events. G errors
  rescaled by Eq. 9–10; Level 0 refit reproduces the published t_E (checked on 6 events). Not emulated: colour and RP
  cuts, u0/t_E error cuts, Extractor cuts, visual inspection (ASSUMPTIONs in `SelParams`).
- Injections on 166 PSPL-subtracted candidate light curves (30 per cell; t_E 10/30/100/300 d, ρ 0.01/0.1, u0 < 2,
  t0 uniform in the DR3 window): **W3 2 / 240 selected (0.8 %)**, PSPL controls 17 / 120 (14 %); ε_W3/ε_PSPL ≈ 0.06.
  Both selected W3 events have u0 ≈ 1.9 (outside the caustic) and are not flagged; **0 / 240 are selected and
  flagged**. The fitter alone flags 139 / 240 (58 %). W3 events fail the skew–Abbe cut (213), the score (223), skewness
  < 0 (153) and the parallax χ²/dof (141).
- Fits of all 363 (0 errors): one flag, 4053892503992268288 (Sample B, ΔBIC −40.3, `E2neg`/`N1neg`, t_E ≈ 346 d,
  u0 ≈ 1.47). Light curve inspected: the event peaks at the end of the DR3 window (no post-peak data) and the
  baseline scatters by ~0.4 mag; the exotic fit uses the baseline wiggles. Not a candidate.

**Revisit if.** Gaia DR4 publishes epoch photometry for all sources (then a W3 screen before any microlensing
selection becomes possible); Gaia publishes the Extractor definitions; or a microlensing catalogue appears whose
selection does not require a single brightening.

## D-062 W3 in the MOA-II 9-year release: thin in-house adapter, notch pre-screen, gb22 pilot null and first W3 rate limit (2026-10-08)

**Decision.** Read the MOA-II 2006–2014 release (NASA Exoplanet Archive) with a thin adapter,
`jwst_anomaly.moa.MoaField`, a `signatures.LightCurveSurvey` over one field's bulk tar. Its light curves are read in
place through a tar-header index; they are difference fluxes in counts, returned by the new
`signatures.standard_flux_light_curve`. Its events are all Cut-0 objects of the field. The pinned files are the
metadata tarball (97 MB) and `gb22.tar` (3.5 GB; stated reason: light curves exist only as per-field tars, and gb22 is
the smallest). The W3 search is `scripts/w3_moa.py`:
- a notch pre-screen. A deficit must sit below both its flanks and the median flux, normalised by the light curve's
  own red noise. It must be single, and its epoch must not be shared by improbably many objects of the field or chip;
- the `w3_microlensing` fitter with the blend flux free (`linear_fluxes(f_min=inf)`), t_E ≤ 1,000 d and ρ ≤ 0.3;
- cheapest-first vetting that stops at the first failure. It adds MOA-specific ordinary tests (neighbour objects,
  seeing/airmass/sky regressors, a trapezoidal eclipse model, ≥ 3 nights in the feature, a night jackknife);
- injection-recovery through Cut-0 emulation, pre-screen, fit and vetting;
- a 95 % limit per monitored star per year. N_s comes from the N_s-per-Cut-0-object ratio of Nunota et al. 2024.
- `fit --chunk K/N` and `merge-chunks` (as D-059) for fields too large for one session; gb22 no longer needs them.
The 3.5 GB download is for cloud sessions; the owner's machine need not fetch it.

Pilot result for gb22: 18,599 light curves, 30 flags, **0 survivors**. Γ₉₅ ≈ 1.2–4.4 × 10⁻⁶ per star per year for
t_E = 10–300 d and 0.6–1.1 × 10⁻⁵ at 3 d (docs/exotic_limits.md "W3 in MOA-II (pilot: gb22)").

**Alternatives rejected.**
- astroquery `NasaExoplanetArchive`: TAP has no MOA table (2026-10-08).
- merida 0.3.2: it scrapes a temporary Firefly workspace URL with a spoofed User-Agent and has heavy runtime dependencies.
- qusi / ramjet `MoaDataInterface`: it reads internal feather files, not the public release.
- Per-object HTTP files (`data/Contributed/MOA/gb{F}/R/{C}/…ipac`): fine for a few events, but undocumented, and
  18,599 requests where a whole field is screened.
- Fitting every shape pass: 1,058 light curves at ~26–35 s each (~2.5 h on 4 cores; chunk 1/8 of them,
  `results/w3_moa/fits_gb22_chunk1of8.ecsv.gz`, flagged 131 / 133, so the flag carries no information there);
  97 % of them sit at epochs shared with other objects.
- The first pre-screen statistics:
  - a plain box mean against the median: 10,367 passes at S < −8, dominated by trends and season offsets;
  - a notch against the flanks with white-noise errors: 13,462 passes;
  - the notch normalised by its own spread: 4,855 passes, and it passed 44–76 % of Cut-0-detected PSPL controls (a baseline box
    between bright flanks). Requiring the box below the median flux too brought the PSPL controls to 0/250.
- Keeping parallax in the screen fit: it can only remove flags, so it is fitted in vetting. Unbounded t_E/ρ: one fit
  took 200 s with an exact finite-source integral on every epoch.

**Evidence.**
- `metadata.ipac` has 2,409,061 rows, the Cut-0 count of Nunota et al. 2024 (arXiv:2410.23553, Sect. 2.2).
- Cut-0 cuts: Koshimoto et al. 2023 (arXiv:2303.08279), Table 2. The archive page quotes the older Sumi et al. 2011
  values.
- N_s per Cut-0 object over 20 fields: median 188, range 131–240. gb22 is not in Nunota et al.'s Table 1.
- First vetting pass: 9 survivors, inspected; they are one- or two-night drops and spike-less flat dips. The eclipse
  model, the chip-level shared-epoch test and the ≥ 3-night rule were then added, and the injections were re-run
  after the change: 101/600 W3 injections recovered, 0/100 PSPL controls. The largest efficiency loss is the
  variable-baseline test (35 % of carriers have χ²/dof > 2 alone); a stricter carrier definition gave 171/600.
- D-054's "no verified bulk MOA endpoint" no longer holds.

**Revisit if.**
- NExScI publishes a TAP table or API for MOA, or a maintained reader appears.
- MOA publishes machine-readable W3-relevant efficiencies or image-level injections. The Cut-0 emulation is
  light-curve level.
- The other 21 fields are screened. The largest tars are 474 GB, so per-object HTTP or a cloud session is needed.
- A survivor appears: stop and report to the owner (/vet-candidate).
- Any W3 limit is quoted outside the repository: `needs-human` (D-054).

## D-063 W5 count-deficit screen: DR10 Tractor counts aggregated per nest4096 on Data Lab, astropy-healpix, cross-region null; 340.5 deg² null and first W5 limit (2026-10-08)

**Decision.**
- **Count maps, not catalogues:** `jwst_anomaly.countmap.LegacySurveysCountMap`, a `signatures.CountMapSurvey`
  (new protocol: `count_map()` per HEALPix pixel, `nside`, `area_deg2()`; it also answers `catalogue()`). Legacy
  Surveys DR10 `ls_dr10.tractor` is aggregated server-side on the Data Lab TAP with `GROUP BY nest4096`, three
  queries per 2° × 2° chunk (galaxies r < 23.5 extended unmasked; all sources with depth, nobs, E(B−V); sources
  with any maskbit → unmasked fraction w). Pixels split between chunks are summed (`combine_duplicates`).
- **Geometry:** astropy-healpix (now a core dependency; BSD-3-Clause, wheels everywhere).
- **Prediction:** `countmap.deficit_profile` = `exotic_sim.count_ratio` (n = 1, ε < 0) with measured counts of the
  same selection, |μ| capped at 30 at the critical curve.
- **Screen** (`scripts/w5_counts.py`): scipy-FFT local least-squares matched filter on a 0.25′ raster, θ_E = 2–32′;
  **null from a second, disjoint region** (each region's ordinarily-vetted peaks calibrate the other), with an
  exponential tail fit for the expected number of false peaks N_false; detection at N_false < 0.01.
- **Vetting**, cheapest first, the same code for flags and injections: mask, depth, depth_edge, dust, Gaia DR3
  bright star, HyperLEDA large galaxy, Wen & Han 2024 cluster, cosmic variance.
- **First run:** desA + desB (RA 20–40° and 50–70°, Dec −30° to −20°), 340.5 deg², 12.5 M galaxies: 40 flags,
  **0 survivors**; 95 % sky density of θ_E = 8–32′ lenses **n₉₅ ≈ 0.012–0.018 deg⁻²**; blind below θ_E ≈ 6′
  (docs/exotic_limits.md "W5 count deficits").

**Alternatives rejected.**
- Per-object catalogue downloads (~12 M rows) and DR10 random catalogues (~19.9 GB per file, unordered on the sky,
  so no byte-range region cut): counts, a mask proxy and depth are all the statistic needs.
- healpy 1.20.1 (GPL-2.0, no Windows wheels); healsparse 1.15.0 / hpgeom (GPL-3.0-or-later; only needed for DES/LSST
  mask files, not used); HSC-SSP PDR3 randoms (account needed, smaller area); DES Y6 Gold masks (healsparse tooling).
- Void finders (VIDE, REVOLVER, Pylians, 2-D tunnel finders): they find the under-densities that are the
  background W5 must beat, not the target; voids enter through the null.
- A Gaussian significance: the Poisson σ of A is ×1.1 (2′) to ×5.8 (32′) too small because of clustering.
- A null from the raw peaks of the other region: desB's artefacts (NGC 1398, deep tiles) raised desA's thresholds;
  peaks the ordinary tests remove are now left out of both samples.
- Veto radii growing with θ_E (a mimic anywhere inside θ_E): they vetoed 90 % of random 32′ positions. A mimic now
  vetoes only if it reaches the core (0.5 θ_E) and can empty ≥ 10 % of it.
- Data Lab ADQL: sub-selects, CASE, SIGN and GROUP BY on expressions are rejected (2026-10-08).

**Evidence.**
- Prior art: no published count-deficit search for negative-mass lenses was found (arXiv API, "negative mass" /
  wormhole with number counts, depletion, deficit, void; the 42 INSPIRE citers of Safonova, Torres & Romero 2001,
  astro-ph/0104075, which predicts the "central void"). Survey limits so far come from SDSS quasar lensing
  (Takahashi & Asada 2013, arXiv:1303.1301). A W5 limit would be new: `needs-human` before any outside use.
- Offline tests on synthetic Poisson maps: an injected deficit gives A = 1.10 ± 0.25 at the centre; the null map
  gives median Z ≈ 0; the chunk-edge artefact (pixels split between chunks gave a deficit along every chunk
  boundary in the pilot) is fixed and tested; pixels cut by the outer region border are dropped (code review).
- Injections (5,232 in total, 3 realisations per region and θ_E) through the full screen and vetting:
  ε = 0.50–0.73 at θ_E = 8–32′, 0.005 at 6′, 0 at ≤ 4′. The `depth_edge` test was added after inspecting the one
  first-pass survivor (a deficit of all sources along a deep-tile edge); the injections ran after the change.

**Revisit if.**
- Data Lab exposes DR10 randoms or a pixelised mask; or the w proxy shows systematics in a larger area.
- θ_E < 6′ becomes a priority (deeper counts, e.g. HSC or Euclid, or catalogue-level positions), or θ_E ≳ 1°
  (a larger contiguous area; the regions are 10° high).
- A survivor appears: stop and report to the owner (/vet-candidate).
- A published count-deficit search appears (compare; do not repeat). Any outside quotation: `needs-human`.

## D-064 W1/W2 in rejected lensed-quasar pairs: LS DR10 cannot decide them; measured control efficiency 0/5 (2026-10-08)

**Decision.** Test the pairs that lens searches rejected for lack of a lens galaxy with the unchanged D-056 chain
(`scripts/w12_niq.py`, importing `w12_lenscats`; D-056 Faber–Jackson calibration fixed). Add a control sample of real
lenses from the same tables to *measure* the test's efficiency, which D-056 had to assume.
- Inputs: VizieR tables pinned by the sha256 of their data lines (the ASU-TSV header carries the request time):
  - Lemon et al. 2023 table1: UQP / QSO pair rejected, lens / quad control;
  - SQLS DR3/DR5/DR7 candidate tables (Inada et al. 2008, 2010, 2012): "no lens(ing) object", "QSO pair" and
    "binary" rejected, "SDSS lens"/"known lens" control; SQLS QSO+star / different SED / not QSO rows describe
    another companion and are dropped without vetoing. A Lemon non-pair class (QSO + star, projected, …) vetoes a
    rejection of the same system, and Lemon classes with "?" are undecided;
  - Hennawi et al. 2006 binaries, for vetting.
- Mismatched or missing pins are refused.
- Sample rules (ASSUMPTIONs):
  - catalogued separation ≤ 3″, because positions are one image;
  - transitive 3″ merging, with a group containing a lens counted as a control;
  - the LS pair must match the catalogued separation within 0.5″.
- Vetting adds image colour (|Δ(g − z)| ≤ 0.5) and quoted-redshift agreement (|Δz| / (1 + z) ≤ 0.01).

Result: 0 of 5 decided control lenses (1.9–2.6″) show their lens galaxy. Rejected: 24 "none". Of these, 10 are
colour-mismatched, 2 are catalogued binaries and 1 has two redshifts; 11 remain untestable. No limit and no candidate
(docs/exotic_limits.md "W1/W2 in rejected lensed-quasar pairs").

**Alternatives rejected.**
- Deriving a dark-lens fraction from the rejected "none" count: control lenses give "none" too, so k carries no
  information.
- Separations up to 6″ (D-056's θ_E ≤ 3″): SQLS positions are one image, so 3–6″ pairs never fit the 3″ image search.
- Hashing the raw ASU-TSV: its header embeds the request time, so a fresh download never matches.
- The Lemon lensed-quasar database (HTTP 500 again on 2026-10-08).
- Gaia GraL invalidated candidates (Stern et al. 2021): mostly star pairs, with one quasar pair.
- Williams et al. 2018: no separation or redshift.
- Dawes et al. 2023: unconfirmed candidates, not rejections.
- NIQ tables that exist only in arXiv LaTeX (Lemon 2018/2019/2020, Anguita et al. 2018, Agnello et al. 2018):
  deferred. The Lemon 2023 UQPs and the SQLS rejections are the machine-readable superset.

**Evidence.**
- `results/w12_niq/summary.json` and `systems.ecsv`.
- The contact sheets show lens light blended into the control images, and blue+orange rejected pairs.
- Review bugs fixed before merge:
  - the Hennawi coordinates are sexagesimal and silently matched nothing; `binary_match` now raises if nothing
    parses;
  - wide pairs were tested on unrelated LS pairs;
  - the greedy dedup was not transitive;
  - Lemon's `z2` is a second quasar redshift only when `n_z2` is blank (or "zqso="). "z_lens=" and "zgal=" are
    other objects, and flagged values are not used;
  - a rejection that another catalogue classifies as a non-pair is now vetoed (J0947+0247);
  - SQLS pair-format tables put the quasar z on the primary row and θ, the comment and the companion's z on the
    next row. Both are now read, and flagged redshifts ("(") are not used.
- The Tractor rows the test ran on are pinned too (sha256 of the sorted rows), and so are the brick summary and
  every VizieR table. VizieR error, empty or truncated responses are refused.

**Addendum (2026-10-09, archival HST).** `scripts/w12_niq_hst.py` tests the pairs with HST F814W imaging: 2 of the 11,
and 2 controls. It fits two Moffat PSFs, removes the halo residual from the profile facing away from the other image,
and measures the residual flux between the images.
- Validation: both controls are detected (S/N 42 and 16), and injections run through the whole chain are recovered
  at ≥ 5σ above the baseline to F814W = 23. Errors are the larger of a drizzle-corrected pixel error and the
  empirical scatter of same-area sky apertures; the latter dominates (≈ 2.3×).
- Result:
  - J0130+0725 has no lens light to F814W ≈ 23 (injection-calibrated): ≈ 3.5 mag below a typical ordinary lens and
    ≥ 0.7 mag below a 2σ under-luminous one (F814W − z assumed 0–0.6). A binary quasar remains the untested
    ordinary explanation;
  - J0728+2607 is inconclusive with the Moffat model (PSF-core residuals). With an empirical PSF from 6 Gaia stars in
    the same cutout (`w12_niq_epsf.py`; control S/N 49) it has no light between the images (S/N 1.6; injections
    recovered to F814W = 23), ≥ 1.8 mag below any ordinary lens.
- Rejected alternatives:
  - the "combined_skycells" HAP cutout, whose WCS does not describe its pixels (separations of 10⁴″);
  - an aperture statistic without the halo correction, which gave S/N 18–24 from PSF mismatch alone;
  - white-noise aperture errors on drizzled pixels: S/N is ~3.5× too high once the drizzle correlation and
    large-scale sky structure are included;
  - injections that skip peak finding and the pair check (not the whole chain);
  - a 0.15″ core mask, which left core residuals.

**Revisit if** HST, Euclid or HSC PDR3 image models (PSF-subtracted) are available for the 11 colour-matched pairs;
spectra of both images can be compared (binary vs lens); or the LaTeX-only NIQ tables add pairs of 2–3″.

**Addendum (2026-10-23, DESI DR1 spectra, `results/w2_desi_pairs/RESULT.md`).**
- **Data:** DESI DR1 has QSO spectra of both components for 13 rejected pairs, among them J0130+0725, J0728+2607 and
  J1428+0500. It also has them for 9 catalogued lenses and 7 Hennawi binaries (θ ≤ 10″), all retrieved anonymously
  through Data Lab TAP and SPARCL.
- **Test:** a frozen spectral test, the mean per-pixel χ²_ν of continuum-normalised broad-line profiles with a free
  scale, T = 3.
  - It passes on the lenses: 9 / 9 with D ≤ T (8 of them at ≤ 2″, where fibre crosstalk may help). D ≤ T means
    no difference above threshold, not one source (wording corrected 2026-10-23); 6 / 7 binaries met it too.
  - It fails on the binaries: 1 / 7 "two spectra", with ≥ 70 % required.
  - So it is not validated. 11 of the 13 rejected pairs' spectra were not retrieved or read. The other 2 share their
    DESI TARGETIDs with two Hennawi controls and were read only in that role.
- **Rejected:** changing T, the lines or the statistic after the gate; a Δv rule seen post hoc on the same controls.
- **Revisit if:** a pre-registered statistic that flags binaries with small catalogued velocity differences (binned
  profiles, velocity offsets, narrow absorbers) is validated on fresh, non-overlapping controls.

**Addendum (2026-10-23, evidence audit, `results/w2_desi_pairs/evidence/RESULT.md`).**
- **Decision:** stop the DESI spectral-classifier line for the rejected pairs. Do not tune a statistic on the same
  DESI DR1 pairs as a default next step.
  - Jing et al. (arXiv:2505.03103) have already classified all 13 pairs visually from images and these spectra, and Chen et al.
    (arXiv:2607.00093) analysed J0130+0725.
  - A new statistic on them adds no observation.
- **Sample caveats** (recorded; the D-064 sample is not re-run):
  - VizieR "UQP" covers 13 of the 17 systems Lemon et al. 2023 discuss in §4.2 as likely lenses, J2244−0550 among
    them. D-064 counts them as rejected.
  - The "binary" labels of J0740+2926 and J1035+0752 (Hennawi et al. 2006, low-S/N spectra) and of J1511+3357 and
    J1552+3009 (SQLS IV, similar Subaru spectra) rest on not seeing a lens galaxy. They are not evidence against a
    dark deflector. SQLS IV §4.1: its process "would reject 'dark lenses'".
  - J2351+1343's DESI pair is not the catalogued pair (SLED: QUASAR+STAR), so #472's 13 holds 12 catalogued pairs.
  - Dropping the two disputed Hennawi controls (J0740, J1035) leaves the #472 gate at 1 / 5 binaries flagged, still a
    failure.
- **Rejected:**
  - re-running #472's statistic or a variant on the 11 unread pairs;
  - classifying from catalogue strings without their basis.
- **Revisit if:**
  - a pre-registered PSF test, with a known lens and a no-lens pair at ~1″ from the same HST programmes as controls,
    is affordable for J2244−0550 (HST GO 17308 F814W and GO 17455 F160W, public);
  - or the proposers publish those data;
  - or independent resolved spectra or a time delay appear for an open pair.

## D-065 W5 in Euclid Q1: no count screen; next is a radial-shear screen on Euclid Q1 shapes (2026-10-09)

**Decision.** Do not port the D-063 count screen to Euclid Q1. Its S/N at θ_E ≤ 4′ is limited by galaxy clustering,
not by the galaxy density, so 1.8× more galaxies buy ×1.0–1.4 in S/N on a fifth of the area
(docs/exotic_limits.md "Euclid Q1"). For θ_E < 6′ the next W5/W1 test is the sign of the tangential shear of
Euclid Q1 MER shapes around trial centres (radial for a negative-mass lens; forecast floor θ_E ≈ 26″).

**Alternatives rejected.**
- Euclid Q1 MER counts per pixel (IRSA TAP): ~3.5 h of row queries for 63 deg²; forecast in Evidence.
- HSC PDR3 counts: account needed (D-063), and the same clustering limit applies.
- Server-side HTM aggregation (`GROUP BY floor(htm20/65536)` works on IRSA for small boxes) would avoid the row
  fetch, but does not change the S/N argument.

**Evidence.** `results/w5_counts/euclid_q1_feasibility.json` (densities observed 2026-10-09; gains and the shear
floor are model_prediction with stated ASSUMPTIONs). IRSA TAP: a `ra BETWEEN` / `dec BETWEEN` box of 0.25 deg² did
not return in 5 min; the same selection with `CONTAINS(POINT, CIRCLE)` (r = 0.25°) returned 12,323 rows in 40 s;
`COUNT(*)` in a 0.1° disc takes a few seconds.

**Revisit if.**
- A count screen is needed at θ_E ≈ 4–6′ specifically (a larger Euclid release removes the area penalty), or a
  counts-in-cells measurement gives a Euclid/DR10 clustering-variance ratio well below 0.5 (ASSUMPTION range 0.5–1;
  at 0.3 counts would reach ≈ 4′).
- The shear screen finds Euclid Q1 MER moments unusable (PSF anisotropy) and no PSF-corrected shape catalogue is
  public.

## D-066 W5/W1 radial-shear screen on Euclid Q1 MER shapes: `ApertureMass` reused; pilot null; limits at θ_E = 1–2′ (2026-10-09)

**Decision.** Search for radial (negative-mass) shear with the D-050 catalogue aperture-mass statistic
(`exotic_screens.ApertureMass`, point-mass filter over 1.5–3 θ_E) on Euclid Q1 MER SExtractor moments, fetched per
0.3° disc from IRSA TAP. The catalogue `position_angle` is used as PA east of north. Pilot: EDF-F and EDF-S.

**Alternatives rejected.**
- Trusting the TAP column description (`position_angle` "CCW/x", THETA_IMAGE): image moments on 25 MER VIS cutouts
  contradict it (89° off); using it would flip tangential and radial, the sign under test.
- A PSF-corrected shear catalogue (none public for Q1); building one (KSB/metacal) before a pilot shows need.
- Grid-centred injections on deconvolved shapes: best case; injections are off-grid, applied to observed moments
  before the cuts, and detected at any centre within one step.

**Evidence.** docs/exotic_limits.md "Euclid Q1 radial-shear screen"; `results/w5_shear/`. Four SZ clusters show
tangential shear (S = −1.5 to −5.0); pilot null (0 flags); whole-chain injection efficiency 0.88–1.0 at 1′, 1.0 at
2′, ≤ 0.18 at 30″; n₉₅ ≈ 8.1 deg⁻² (1′), 12 deg⁻² (2′). *Superseded by D-067: these injections were boosted
~1.6× over R (efficiencies optimistic; `results/w5_shear/screen.json` is from that chain).*

**Revisit if.**
- A PSF-corrected Euclid shear catalogue becomes public (DR1), or cluster-calibrated R differs from 0.5 by > 30 %.
- PSF-anisotropy gradients (star ellipticity maps) show radial patterns on trial-centre scales.

## D-067 W5/W1 radial-shear survey of all Euclid Q1 Deep Fields: tile fetch by `tileid`, R calibrated on SZ clusters, 60 deg² null (2026-10-09)

**Decision.** Fetch the Q1 MER rows per tile with the indexed `tileid` (tile list from IRSA ObsCore), run the D-066
screen per Deep Field with per-tile PSF sizes and a coverage cut, and keep R = 0.5 after calibrating it on the four
SZ clusters against NFW haloes (in-house ~40 lines: Wright & Brainerd 2000, Duffy et al. 2008 c200).

**Alternatives rejected.**
- More `CONTAINS` discs: overlapping discs double-count rows; edge handling per disc wastes the aperture border.
- 16 concurrent IRSA queries: 504 gateway time-outs and no rate gain over 8 (0.3–0.5 queries/s).
- colossus / pyccl for the NFW shear: a dependency for textbook formulas checked by a unit test.
- Per-source photo-z (`euclid_q1_phz_photo_z`) for Σ_crit: one source plane bracketed by z_s = 0.8–1.2 already
  keeps R within ±20 %; worth it only if a flag needs a mass.

- Injecting by shearing the observed moments before the PSF deconvolution (the D-066 chain): it boosts the
  injected shear by tr_obs/tr_int (median 1.6×) over the R calibrated on deconvolved shapes. Injections now shear
  the deconvolved shapes by R g (`inject_shapes`); D-066's pilot efficiencies are superseded.
- Rule exception (scripts/CLAUDE.md "inject through the whole chain"): injections skip the resolved cut and the
  deconvolution, because R is defined on their output; a shear changes a galaxy's size by a fraction ≲ |g||ε|, so
  only galaxies that close to the 1.2 × PSF cut could change their resolved flag (kept as is). Revisit if a limit
  hinges on them.

**Evidence.** docs/exotic_limits.md "Euclid Q1 radial-shear survey"; `results/w5_shear/`. R = 0.56 ± 0.16
(z_s = 1, scaled by √(χ²/dof)); 0 flags over 60 deg² at θ_E = 1′, 2′, 4′; known clusters negative (8/9);
efficiency 0.43–0.58 at 1′ (no limit), 1.0 at 2′–4′; n₉₅ ≈ 0.049 / 0.051 deg⁻² at 2′ / 4′.

**Revisit if.**
- Euclid DR1 (≈ 2000 deg², PSF-corrected shapes) is public: the same code scales by tile.
- A flag appears, or per-tile star ellipticity gradients show radial patterns on 1–12′ scales.

## D-068 W3 MOA-II: calibrated variable-baseline test, LF-drawn injections and a streaming per-field pipeline; eleven fields, 0 candidates, limits withdrawn pending the corrected chain (2026-10-08, amended 2026-10-09)

**Decision.**
- The variable-baseline vetting test's threshold is the 95th percentile (ASSUMPTION, `BASELINE_Q`) of the field's
  quiet-light-curve χ²/dof about a constant (`calibrate_baseline`), not D-057/D-062's fixed χ²/dof > 2. gb22: 5.31.
- Injection magnitudes are drawn from the luminosity function (`--sampling lf`), 200 per t_E × ρ cell, 40 PSPL
  controls per t_E; uniform draws re-weighted afterwards left n_eff ≈ 20–25 per cell.
- Fields are streamed (`jwst_anomaly.moa_stream`): concurrent HTTP range reads of the uncompressed tar, member
  headers resynchronised at 512-byte blocks, each pre-screen process reading its own 64 MiB ranges. Per-chunk
  pre-screen tables are tracked in `results/w3_moa/prescreen/` (resumable, `merge-prescreen` checks the Cut-0
  member count); per-field limit and vetting summaries in `results/w3_moa/`; `combine` joins fields. Streamed tars
  are pinned by per-range sha256 in the manifest.
- gb22 re-run: same 30 flags as D-062, 0 survivors; Γ₉₅ ≈ 1.1–3.3 × 10⁻⁶ per star per year at t_E = 10–300 d,
  4–7 × 10⁻⁶ at 3 d (docs/exotic_limits.md "Calibrated re-run of gb22"). Supersedes D-062's limit numbers.
  **Amendment 2026-10-09: withdrawn** with every limit of this PR, see below.
- Six fields (gb11, gb16, gb19, gb20, gb21, gb22; 296,618 Cut-0 light curves) streamed and vetted: 234 flags,
  **0 survivors** (docs/exotic_limits.md "Six fields"). New vetting tests `exotic_in_domain`, `smooth_dip`,
  `feature_bracketed`, `step_ramp`.
- **Injections are vetted against the injected source magnitude** as the reference flux of `exotic_in_domain`'s
  source-flux bound (f_s ≤ 3 × reference). Before 2026-10-09 they got the lenient 14.2 mag default while real flags
  got their DoPHOT / Gaia magnitude, so efficiencies, and all per-field and combined limits, were biased high; the
  tables were removed. `CHAIN_VERSION` tags injection and limit tables; `limit` and `combine` refuse another
  version (bump it on any pre-screen, fit or vetting change).
- `fit`, `vet`, `inject` and `run-field` refuse to run without MulensModel (without it, the parallax refit and the
  binary-lens test were skipped silently). A failed Gaia XMatch is retried and, if it still fails, marks flags
  without a DoPHOT magnitude incomplete instead of passing them on the lenient default.

**Alternatives rejected.**
- Fixed χ²/dof > 2: 35 % of quiet carriers exceed it from the red noise of difference photometry alone; it was the
  largest injection loss (73 of 188 vetted injections, D-062).
- Downloading each tar: 3.5–508 GB per field (≈ 2.4 TB total); against the cloud-disk owner decision.
- A parent process that reads ranges and ships bytes to workers: OOM-killed on gb21.

**Evidence.** The HTTP stream reproduces the D-062 pre-screen of gb22 exactly (81 s vs 203 s from the local tar);
the fixed-width IPAC parser is bit-identical to the token parser on 15,461 gb21 light curves; gb21 test throughput
258 light curves/s (72 MB/s, 4 cores). 103 / 2,000 W3 injections recovered, 0 / 200 PSPL controls.

**Amendment (2026-10-09): vetting tests added post hoc, each after an automated survivor was seen.** Each change
was followed by a re-run of the vetting and the injections of every field already done, so every quoted limit
includes its cost. The candidate records give the verdicts and evidence.
- `exotic_in_domain` and `smooth_dip`, added after seeing gb20-R-4-0-49379: u₀ < 2 and f_s ≤ 3 × the
  DoPHOT or Gaia DR3 RP reference flux; a Gaussian dip. Cost on gb22: 0 injections.
- `feature_bracketed` and `step_ramp`, added after seeing gb19-R-4-4-31159: ≥ 20 epochs before the ingress and
  after the egress; a level change plus ramp. Cost on gb22: 103 → 93, mostly long t_E near the data edges.
- The night jackknife ranks nights against the eclipse model too (c63a6d5), changed after seeing
  gb7-R-8-6-94052. Its exotic preference over a trapezoid rested on one night. Paired re-run of gb22 (same seed,
  identical injections; W3 recovered, ρ summed, before → after): 3 d 8 → 6 (−25 %, small counts), 10 d 22 → 20,
  30 d 28 → 25, 100 d 22 → 21, 300 d 13 → 13; total 93 → 85 (−9 %). PSPL controls falsely recovered: 0 → 0.
  Kept, not tuned: a threshold set on the flag that motivated the change would be post hoc twice over.
  All seven paired fields (W3 recovered, before → after; PSPL controls 0 → 0 in each):

  | Field | 3 d | 10 d | 30 d | 100 d | 300 d | total |
  |---|---|---|---|---|---|---|
  | gb22 | 8 → 6 | 22 → 20 | 28 → 25 | 22 → 21 | 13 → 13 | 93 → 85 |
  | gb21 | 13 → 6 | 26 → 24 | 33 → 32 | 28 → 28 | 10 → 10 | 110 → 100 |
  | gb20 | 10 → 6 | 21 → 18 | 22 → 22 | 27 → 25 | 5 → 5 | 85 → 76 |
  | gb19 | 12 → 7 | 24 → 21 | 24 → 24 | 18 → 18 | 8 → 8 | 86 → 78 |
  | gb16 | 11 → 9 | 20 → 19 | 24 → 23 | 20 → 20 | 4 → 4 | 79 → 75 |
  | gb11 | 19 → 15 | 27 → 23 | 31 → 31 | 24 → 23 | 14 → 14 | 115 → 106 |
  | gb12 | 14 → 9 | 28 → 27 | 29 → 28 | 28 → 28 | 9 → 9 | 108 → 101 |
  | **all** | **87 → 58 (−33 %)** | **168 → 152 (−10 %)** | 191 → 185 (−3 %) | 167 → 163 (−2 %) | 63 → 63 | **676 → 621 (−8 %)** |

  The cost falls on short events, whose caustic spikes are sampled on one or two nights.
- **Follow-up (open, 2026-10-09).** gb18-R-9-4-24509 passed the whole chain and was explained only by a vetting note
  (docs/candidates/gb18-R-9-4-24509.md). The chain gap: season-offset comparisons are made only against single-lens
  models, and a slow (flat-bottomed) dip with per-season offsets is missing. The chain is left unchanged until the
  field queue finishes. Adding the model then would be post hoc and needs a full re-injection. Until then gb18 has an
  efficiency-only table and is not in the combined limit: its injections did not face that model, so a zero-event
  limit from them would overstate ε.
- **Follow-up (open, 2026-10-09), from gb17-R-6-1-3829** (a ~30-d periodic variable; docs/candidates/gb17-R-6-1-3829.md):
  (i) repeated deficits *inside* the feature window are not counted (the test looks only outside it);
  (ii) there is no periodicity / variability-model test;
  (iii) there is no cap on the exotic fit's own χ²/dof (6.7 here).
  gb17 is efficiency-only and stays out of the combined limit, as gb18. When the queue finishes, all four gaps
  (these three and the gb18 slow dip with season offsets) are added as one recorded post-hoc chain change, and
  every field is re-injected.

**Evidence (2026-10-09).** `/code-review` of the PR diff found the injection reference-flux gap; a unit test now
checks that injections carry `ref_mag = I_s`. A fresh cloud venv without the `mulens` extra re-vetted gb11 without
the binary-lens test (no error), which is how the silent skip was found.
The recovered counts in the amendment table above come from the old chain (injection reference flux at the 14.2 mag default) and are kept only as a record of the relative cost of each vetting change. Every limit and efficiency table was withdrawn (`CHAIN_VERSION` 2026-10-09.1); per-field limits wait for a re-injection under the corrected chain.

**Addendum (2026-10-10): the four open gaps closed as one recorded post-hoc chain change** (`CHAIN_VERSION`
2026-10-10.1; motivated by gb17-R-6-1-3829 and gb18-R-9-4-24509, both seen before the change). New tests in
`vet_one`, all thresholds ASSUMPTIONs in `Params` of `scripts/w3_moa.py`:
- `exotic_chi2_cap` (gap iii): the exotic fit's χ²/dof, caustic-spike epochs of a repulsive model left out (model
  above f_s + f_b by max(3 median errors, 0.1 f_s)), must not exceed the field's calibrated quiet-baseline threshold. Spikes are left out because a synthetic
  W3 event fitted at another u0 local optimum (u0 0 vs 0.5) had χ²/dof 1.69 against a white-noise 1.06, all of it on
  the spike nights.
- `residual_deficit` (gap i): `deficit_scan` z_min of the exotic residuals, spike epochs left out, must stay above
  −`REPEAT_S`; scanned over the whole curve so the red-noise normalisation never rests on a short window.
- `periodic_variable` (gap ii): Lomb–Scargle of the nightly mean exotic residuals, P = 2.5–100 d (above the
  one-sample-a-night Nyquist period; shorter periods are not covered); fails when a 2-harmonic term at that period
  gains BIC > 25 (errors rescaled by the harmonic model's χ²/dof, F-test-like) in **both** halves of the nights.
  One umbra crossing lies in one half; a periodic variable repeats in both. < 30 nights: `_untestable`.
- `slow_dip_seasons` (gap iv): the exotic model with one level per season against a generalized-Gaussian dip
  exp(−|Δt/σ|^p), p = 2–10 (Gaussian to flat-bottomed), with the same season levels.
Validation (synthetic, 6 W3 events in white noise, quiet threshold 1.06): all six pass every new test; unit tests
show a 30-d sinusoid fails in both halves, a single dip and white noise do not, and the slow-dip model recovers a
flat (p = 6) dip with three season levels at χ² ≈ N. Every field's real flags and injections must be re-run under
the new chain before any limit is quoted (`limit` and `combine` refuse 2026-10-09.1 tables).

**Revisit if.** A field's quiet χ²/dof distribution is bimodal or its 95th percentile exceeds ~10 (the quantile then
admits variables); the archive stops honouring range requests; a survivor appears (stop, /vet-candidate, owner).
## D-069 Owner direction 2026-10-09: System A (invention) / System B (destruction) hypothesis rounds; round 1 survivors (2026-10-09)

**Decision.** New hypotheses come from System A groups that work without literature and from minimal axioms, with
space, time, matter and causality not assumed fundamental ("Ontological Reset"). One System B reviewer per text then
attacks it with verified references. Only predictions that pass B get a data inventory and search compute. At least
30 % of exploration goes to quantities current searches never measure. Round 1 kept three conditional survivors,
in order:
- S1, burst twins (GRB/FRB pairs far apart on the sky at any delay);
- S2, flat-kernel SN residuals;
- S3, lensed-transient hybrid images (limit only).

W1–W5 screens and the MOA-II run continue to completion (PR #95), but new compute goes to S1 first.

**Alternatives rejected.**
- More literature-anchored lens signatures only: the owner asked for at least 30 % on unmeasured quantities.
- A blind causal-network timing search between unrelated events as the first test: B-on-A1 §7 shows it is testable
  in principle for pairs ≳ 0.1° apart, but its cost is chance coincidences and observer-side common causes, which S1
  controls with scrambled catalogues at lower cost.
- A2 P2/P3/P4 and all of A3 except P2b: see the failure list in docs/hypotheses/round-1/summary.md.

**Evidence.** docs/hypotheses/round-1/ (three System A texts, three System B reviews, summary). Key gap references
(opened 2026-10-09):
- arXiv:2006.07095, Fermi GBM lensed-GRB search: position-consistent pairs only.
- arXiv:2204.06014, CHIME/FRB lens interferometry: delays ≲ 100 ms.
- arXiv:2406.05047, DES SN lensing magnification.
- arXiv:0903.5303, disordered locality: the common renaming of the A1/A2/A3 links.

**Revisit if.**
- S1–S3 are run and null: start System A round 2 from the recorded failure modes (no Lorentz derivation,
  disordered-locality renaming, cosmic-opacity redundancy).
- The owner changes the 30 % budget or the process.

## D-070 S2 flat-kernel SN residuals: Pantheon+ × LS DR9 photo-z counts, α = 2 column, dilution-corrected limit (2026-10-09)

**Decision.** S2 (D-069) is tested with `scripts/s2_flat_kernel.py`. Pantheon+ residuals against flat ΛCDM
(Ω_m = 0.334) are regressed by WLS on [1, z, X_lens, X_flat]. The columns are per-SN counts of LS DR9 galaxies
(countmap selection, dereddened z < 21 mag, 2′ disc, z_g < z_s − 0.05, Δz = 0.05 shells) relative to the mean
shell densities per photometric region; discs with < 50 % of the region's median galaxy count are dropped.
X_flat uses A3's (1+δ)^α with α = 2 (Poisson-unbiased N(N−1)/E, only shells expecting ≥ 1 galaxy).
Errors come from χ²/dof scaling and a z-matched scramble. A 1 % leverage-trimmed refit is required, and
limits are divided by the shot-noise dilution λ from Poisson-only columns and by the recovery ratio of a whole-chain
injection on half-thinned counts. Result: null; γ < 0.025 mag per unit T/⟨T⟩ (one-sided 95 %).

**Alternatives rejected.**
- α = 2 over all shells: shells with E ≪ 1 turn 3–5 galaxies into X ≈ 20–54 and fake a tight limit.
- α = 1 flat column: 0.995-correlated with the lensing column at this depth, so it is kernel-degenerate.
- Scramble-carrier injections alone as the efficiency: they miss regression dilution.
- A halo-model κ (Shah+2024): more work for a test that cannot detect lensing in Pantheon+ anyway; it belongs to the
  DES-SN5YR follow-up.
- DES-SN5YR first: needs DES Y3 Gold access and setup; Pantheon+ and DR9 photo-z were reachable in one cycle.

**Evidence.** `results/s2_flat_kernel/fit_*.json`; CHANGELOG 2026-10-09.

**Revisit if.**
- DES-SN5YR + DES Y3 Gold are set up: repeat there with the lensing detection as the positive control.
- A3 fixes how its matter column maps to galaxy counts (bias), which sets the conversion of this limit.

**Addendum (2026-10-09, DES-SN5YR).** Repeated with `--sample des` on the DES-SN5YR Dovekie Hubble diagram (DES SNe,
MU/MUERR, P(Ia) ≥ 0.5, positions from the DES HEAD table) and the same LS DR9 z < 21 galaxies, so both samples share
one calibrated chain. 1,414 SNe fitted: γ_F = +0.0009 ± 0.0025, **γ < 0.019 mag per unit T/⟨T⟩** (chain recovery 0.57);
lensing control γ_L = −0.018 ± 0.011 (1.7σ, right sign). Rejected for this cycle: DES Y3 Gold (new access path and a
second galaxy calibration; it is the next step because the control needs deeper galaxies). Revisit if: a deeper
galaxy column is set up, or the full covariances change the errors by > 20 %.

**Addendum 2 (2026-10-18, LS DR10 z < 22).** The deeper column the addendum asked for was tried under a frozen plan
(`results/s2_flat_kernel/dr10/PLAN.md`): photo-z gate Q passed, positive-control gate C failed (γ_L = −0.0145 ±
0.0129, 1.13σ, p 0.25; `results/s2_flat_kernel/dr10/RESULT.md`). Decision: S2 is paused; γ < 0.019 stays the quoted
limit (the DR10 bound 0.0076 is not adopted without a detected control). Rejected: LS DR10 at z < 22 as the control
column (the error grows faster than λ_L). Revisit if: DES Y3 Gold or a spec-z / deeper galaxy column is accessible.

## D-071 S1 burst twins in Fermi GBM: bcat HDU 2 light curves, multi-band max cross-correlation, chi2 twin test; null and limit (2026-10-09)

**Decision.**
- Data: the HEASARC `fermigbrst` catalogue over TAP, and the per-burst bcat files streamed from the HEASARC FTP.
  Each burst is reduced to two matched bands from bcat HDU 2: `PHTFLUXB` (50–300 keV) and `PHTFLUX − PHTFLUXB`
  (10–50 plus 300–1000 keV). The window is T90 ± (0.25 T90 + 2 s) at 64 ms · 2^k, with at most 256 bins. The
  results are tracked in `results/s1_twins/lc_<YYYY>.ecsv.gz` (~150 kB per year).
- Eligibility: at least 2 pulses. A pulse is a summed-band peak whose prominence is ≥ 5σ of the local error,
  sqrt(σ_peak² + σ_base²), so the test includes source Poisson noise.
- Statistic: the maximum over lag of the multi-band normalised cross-correlation ρ, with s = 1 primary. Pairs are
  compared at a common resolution, |Δk| ≤ 1. A coarse s grid (0.5, 0.71, 1.41, 2) is secondary and runs on the top
  pairs only.
- Position cut: sep > 3·sqrt(σ1² + σ2² + 2·3.7°²).
- Chain: ρ > ρ* (95 % of the per-catalogue maximum over 100 pulse-shuffled surrogate catalogues), then a
  noise-consistency χ² test of a scaled, shifted copy with a 10 % flux systematic (p ≥ 10⁻³), then a re-trigger
  veto (delay ≥ 1 d). Only this primary chain gives the limit. A "deep" chain (ρ > 0.90) is post hoc and
  illustrative only: its threshold was chosen after seeing that it yields k = 0.
- Injections start from a smoothed template (Gaussian, 1.5 bins) of a real eligible burst. One member is re-noised
  as A. The other is a copy in another burst's slot, with that burst's background noise and independent
  realisations, at flux ratios 1–0.1, with a per-band gain mismatch of 15 % (see Evidence). The copy gets its own
  window and resolution k from a T90 proxy (A's catalogue T90 × the ratio of 5–95 % cumulative-fluence durations).
  A copy with |Δk| > 1 is lost. Both members pass through the whole chain.
- Limit units: f95 = μ95 / (N_eligible · ε) is the 95 % upper limit on twin PAIRS per eligible burst. The fraction
  of eligible bursts that have a twin is 2·f95.
- `scripts/s1_ingest.py`, `scripts/s1_twins.py`, `src/jwst_anomaly/burst_twins.py`. All thresholds are ASSUMPTIONs
  in `burst_twins.Params`.

**Alternatives rejected.**
- bcat HDU 1 `PHTCNTS` (per-detector deconvolved 8- or 128-channel spectra): in ~40 % of faint bins, single channels
  carry the fill value −9.9e36. Dropping those bins left only the positive fluctuations, which looked like spiky
  multi-pulse bursts (84 % of bursts appeared "multi-pulse"). Rejected after visual inspection.
- TTE/CTIME with our own background fits: these are ~10–100× more data and need per-burst background modelling,
  while bcat already holds rmfit's background-subtracted, deconvolved fluxes. They are kept as a follow-up for
  candidates only.
- GDT (official) for reading: no bcat reader is listed, and bcat is a plain FITS table, so `astropy.io.fits` is
  enough. GDT stays the tool for TTE/RSP follow-up.
- `astroquery.heasarc`: plain HEASARC TAP is one request (4 s) with no extra dependency.
- The 2006.07095 method (position, spectrum and duration cuts first, then cross-correlation): its position-first
  cut is exactly what S1 inverts. We keep its cross-correlation and its warning about single-pulse look-alikes, and
  replace the spectral cut with matched bands inside ρ plus the χ² test.
- Pulse-count threshold judged on the S/N with background error only: it counts Poisson fluctuations on bright
  pulses as pulses.
- Injecting a copy of the burst itself (with its own noise realisation) into the slot: the pair then shares noise,
  so faint copies reached ρ ≈ 1. That version overstated the ratio-1 efficiency (83 % instead of 37 %).

**Evidence.**
- The catalogue has 4,390 bursts and 4,389 bcat files reduced. 1,414 bursts are eligible. 998,991 eligible pairs
  give 510,294 pairs after the position and resolution cuts.
- Results are in `results/s1_twins/summary.json`, `null.json`, `pairs_top.ecsv` and `injections.ecsv.gz`; the PR
  body has the numbers.
- The surrogate null under-predicts the real high-ρ tail about 12×: 1,860 real pairs have ρ > 0.9, against 146
  per surrogate catalogue. It is also near-degenerate, because a 2-pulse burst has one alternative order. So ρ* is
  a pre-screen, not a calibrated false-alarm rate, and the χ² test decides.
- At ρ > 0.9, 0 of the 1,860 real pairs pass χ². At ρ > 0.8, 24 of 40,041 pass, all with a faint partner
  (S/N 5–10).
- Positive controls: 8 same-source re-trigger pairs (GRB 091024, 130925, 150201, 220627, 250702 ×3, and
  bn210925800/bn210926869 at 1.3°) are all removed by the position cut.
- Result: 6 pairs flagged at ρ > ρ* = 0.964. All 6 fail χ² (p < 10⁻¹⁰⁰): they are bright single-envelope
  look-alikes. 0 survivors (primary chain; the post-hoc deep chain also has 0).
- Band mismatch from data. In bright bins (summed S/N > 15) of the 522 eligible bursts with at least 8 such bins,
  the noise-subtracted scatter of ln(f_50–300 / f_rest) has percentiles 10/25/50/75/90 of
  0.075/0.13/0.21/0.32/0.42. The median 0.21 / √2 gives 0.15 per band, which is the injected value. It includes
  intrinsic spectral evolution within each burst, so it sits on the high side for a pure response mismatch.
- Error calibration. In off-T90 bins, std(flux/err) has a median of 1.00 in band 1 and 0.87 in band 2, so the
  quoted errors are about right and the injected noise equals them.
- Injection recovery, primary chain, 15 % jitter, 5,000 injections: 29.5 / 18.0 / 11.4 / 6.4 / 2.7 % at flux
  ratio 1 / 0.5 / 0.3 / 0.2 / 0.1.
  - By copy peak S/N: 72 % above 50, 48 % at 20–50, 7 % at 10–20, and 0 below 10.
  - The χ² test removes 13 % of flagged injections. The window/k emulation loses 0.9 % of copies, so the stored
    window and k (which depend only on catalogue T90) are a small effect.
- 95 % upper limit, primary chain only, at ratio 1: **7.2 × 10⁻³ twin pairs per eligible burst**, i.e. 1.4 % of
  eligible bursts have a twin. Averaged over ratios 1–0.1: 1.6 × 10⁻² pairs per burst (3.1 %).
- Sensitivity to the per-band mismatch (1,000 injections each; ratio-1 efficiency, f95 pairs per burst at
  ratio 1):
  - 5 %: 41.5 %, 5.1 × 10⁻³;
  - 10 %: 38.0 %, 5.6 × 10⁻³;
  - 15 %: 29.5 %, 7.2 × 10⁻³;
  - 20 %: 22.5 %, 9.4 × 10⁻³.
- Deep chain (post hoc, illustrative, not a limit): 49 % at ratio 1, 4.3 × 10⁻³ pairs per burst.
- Closest call: bn100528075/bn250207053 passes χ² only at s = 0.71 (p = 0.0016, ρ_s = 0.943 < ρ*), with a faint
  partner (S/N 12). Vetted as population similarity at low S/N (PR body).

**Revisit if.**
- A better null that keeps the population envelope becomes available (e.g. a pulse-parameter generative model
  fitted to the catalogue).
- TTE-based light curves are needed for short GRBs (64 ms is coarse for T90 < 1 s).
- CHIME/FRB Catalog 2 is added (the S1 second half).
- The χ² systematic (10 %) is calibrated on real same-burst data, e.g. GBM vs. Swift-BAT.

**Addendum (2026-10-21): S1 on CHIME/FRB Catalog 2, 30-event validation; the criteria do not carry over as frozen.**
- **Decision.** The D-071 chain (eligibility ≥ 2 pulses at 5σ, pulse-shuffled ρ\*, scaled-shifted-copy χ² with a 10 %
  systematic) is not used on CHIME dynamic spectra, and no CHIME search or limit is run. The conclusion is limited to
  "this frozen method did not transfer to CHIME"; it says nothing about whether CHIME holds same-signal replays.
  - Adapter `scripts/s1_chime_lib.py`: four 100-MHz sub-bands averaged in S/N units, band Jy scale = median fcf·σ.
  - The issue #153 χ² fixes (union window, detected per-band gain ≥ 0) plus a sub-bin lag (best p among δ with
    detected gains; corrected after /code-review, no outcome changed).
  - `burst_twins.py` is unchanged. Prior work not redone: issue #153 (route, format, 307-burst pilot).
- **Alternatives rejected.** Unweighted or 1/(fcf σ)² channel means (a few channels dominate the noise); integer lags
  (fail exact copies of ~1-bin pulses); fitburst models as noise-free templates (misfit χ²/dof up to 135);
  time-reversal null (#153).
- **Evidence** (`results/s1_chime/`; plan frozen first, `PLAN.md`):
  - 17 / 30 eligible; ρ\* 0.870; 0 / 133 chance and 0 / 3 same-repeater pairs pass.
  - Frozen recovery 6.9 % / 5.9 % / 3.1 % at ratio 1 / 0.5 / 0.3 (G_cal), against ≥ 50 % required (T1 fails).
  - Post hoc: 70.5 % for exact copies in Gaussian noise, 23.6 % in real host residuals, 10.0 % with +1 ms
    scattering at 600 MHz. The beam/calibration gain imprint (RMS 0.18 per band) costs ≤ 1 pp.
  - Same-source ρ median 0.411 (26 pairs) vs different-source 0.513 (380), over all 29 usable files (post hoc).
- **Revisit if.** Real same-signal CHIME pairs (multi-station or multi-beam data of one burst) calibrate the χ²
  systematic, or the owner approves a fitted scattering re-convolution in the twin model.

## D-072 S3 hybrid images: COSMOGRAIL XIX doubles, SIS leg-split lag window, template regression; null without a useful limit (2026-10-09)

**Decision.** S3 (D-069) is tested with `scripts/s3_hybrid.py` on the COSMOGRAIL XIX R-band curves (CDS
J/A+A/640/A105) of the 7 doubles whose published delay is not flagged uncertain and exceeds 3σ. Lags are a
model_prediction under an SIS (ASSUMPTION): the geometric delays of both images are equal, so with a fraction f of the
potential delay on the source leg, s_L − s_T = −Δt (D_ls/D_s + f). f is unobservable, so each direction scans the
window f ∈ [0, 1]; the copy into the trailing image arrives before the leading image when D_ls/D_s + f > 1. Image j is
fitted as a smooth cubic B-spline (3000-d knots) + m F_i(t − τ) + a F_i(t − τ − δ), F_i interpolated across gaps
≤ 40 d; r = a/m. The window maximum of r (|lag| ≥ 10 d) is compared with same-width windows at off-model lags,
and injections go through the same window maximum. A system is screened only if the main-term fit recovers its
published delay. Result: no copy (lowest p = 0.02 of 12 windows, on the main-image wing); sensitivity
r95 ≈ 0.38–2.4 per system and direction, i.e. only copies comparable to the main image are excluded.

**Alternatives rejected.**
- Microlensing B-splines with 120–730-d knots: the published delay is not recovered (6 of 7 systems off by
  50–250 d); only the 3000-d spline recovers 6 of 7 within 3σ + 4 d (J1620: wrong-sign minimum at −206 d).
- Quoting B-on-A1's ~10⁻² sensitivity: quasar variability is red (time-scales ≳ 100 d), so a copy at lags
  ≲ 100 d is nearly collinear with the main image and with a delay error. The same fit recovers r = 0.1 within 0.03
  on a synthetic curve with 20-d variability (`tests/test_s3_hybrid.py`), so the loss is the data, not the code.
- Treating the r95 values as limits: window-maximum injection efficiencies scatter from 0.67 to 4.7 over 8 trials,
  the r(δ) scans have a positive offset of 0.1–1 at all lags (template–microlensing degeneracy), and the null
  windows overlap and are not matched in |lag| (approximate p). They are sensitivities.

**Evidence.** `results/s3_hybrid/summary.json`, `scans.json`; CHANGELOG 2026-10-09.

**Revisit if.**
- A lensed transient with fast intrinsic variability (SN Refsdal, SN H0pe, a lensed FRB) gives a sharp template:
  then a copy at r ~ 10⁻² is reachable (B-on-A1).
- Lens models with image positions and a source position are wanted for quads (the SIS window does not apply).

**Addendum 2026-10-21 (SN Refsdal).** The first "Revisit if" was tested on Kelly et al. 2023's HST light curves
(`scripts/s3_refsdal.py`, `results/s3_hybrid/refsdal/`): S1 as template, copy lags 30–400 d, frozen two-filter
5σ rule. Gates pass; 0 flags; sensitivity r95 0.21–1.69, because the SN 1987A-like light curve is slow and the
F125W main-term fits have χ²/dof up to 3.85. A transient with weeks-long evolution is not sharp enough; the
condition is narrowed to day-scale lensed transients (lensed FRB, fast-rising SN) or a joint microlensing model.
SN H0pe not attempted (three images, fewer epochs; expected no better).

## D-073 D1 distance self-consistency: per-lens H0-free ratio R = D_dt/((1+z_d) D_d) and per-FRB DM-z predictive tails (2026-10-09)

**Decision.**
- Data (derived posteriors, pinned by sha256 in `data/manifests/d1_distance.ecsv`): H0LiCOW-public
  `h0licow_distance_chains/` (joint D_d–D_dt for RXJ1131, PG1115, J1206; D_dt for HE0435, WFI2033) and the
  B1608 analytic shifted-log-normal fits (`B1608_Dd_Ddt_params.dat`; D_d and D_dt independent, so no correlation).
  The TDCOSMO 2025 SDSS1206 `final_D_d.npy`/`final_D_dt.npy` are the same 67,570 samples as H0LiCOW
  `J1206_final.csv` (checked), so they add no lens. Redshifts from the H0LiCOW likelihood = TDCOSMO yaml.
- Statistics per lens (ln space, random-pair Δ = obs − pred; Gaussian pull plus an empirical tail check):
  A prior-predictive vs flat ΛCDM, Ωm ~ U(0.1, 0.5); B flat wCDM, w ~ U(−2, −0.5); C leave-one-out vs the other
  joint lenses with Ωm and a common ln-scale free (absorbs a sample-wide λ_MST); D leave-one-out ln D_dt vs the
  H0 + Ωm fit to the other five. Sidak over 18 (lens, statistic) pulls: 5σ global = 5.53σ local (ASSUMPTION).
- Null: all 23 non-identity permutations of the 4 joint lenses' (z_d, z_s); 100 random permutations of the 6.
  Injection: one lens's D_dt × f, f = 0.7–2.0, through the same statistics.
- FRBs: DM_obs vs the predictive distribution of DM_ISM (the JSON `DMISM`, NE2001 via FRBs/FRB `frb/mw.py`
  `ismDM`; Gaussian ±20 %, ASSUMPTION) + DM_halo U(10, 80) + ⟨DM_cosmic⟩(z)·Δ (Macquart+2020 p(Δ), σ = 0.32 z^−1/2,
  z floored at 0.02, Δ truncated at 20) + log-normal DM_host (e^μ = 68.2, σ = 0.88)/(1+z). The CDF is a deterministic
  grid convolution (ISM + halo analytic, cosmic on a uniform DM grid, host through its analytic CDF/SF; all sums of
  positive terms), so tails stay accurate far below 1e-10. A Monte Carlo floored at 1/N was rejected: it capped z at
  ~4.3σ, below the flag. One-sided tails both ways; Sidak over 2×94 with a global two-sided 5σ gives 5.81σ local
  one-sided. Per burst, the DM at which each tail reaches 5.81σ is recorded, and DM_obs = 0.9× / 1.1× those limits
  is injected through the same chain. All thresholds are in `distance_consistency.Params`/`FRBParams`.
- Code: `src/jwst_anomaly/distance_consistency.py`, `scripts/d1_distance.py`; outputs `results/d1_distance/`.

**Alternatives rejected.**
- hierArc / lenstronomy likelihoods: built for the hierarchical population fit in which λ_MST and anisotropy are
  population parameters, which is exactly what hides a per-lens outlier. The per-lens question needs only the
  published samples and D_s/D_ds, so astropy.cosmology (E(z), tabulated per (Ωm, w)) suffices. Revisit for
  TDCOSMO 2025 lenses whose kinematics exist only as `*_const_processed.pkl` likelihoods.
- `frb` package (FRBs/FRB) for ⟨DM_cosmic⟩: heavy dependency for one integral; the same formula is implemented
  with astropy (Planck18 parameters, f_d = 0.844) and tested against ~900–1000 pc cm⁻³ at z = 1.
- A Gaussian DM residual: the cosmic term is skewed and the host term log-normal (data-inventory §1.2).

**Evidence.** Lenses: max |pull| 1.3σ (RXJ1131, statistic D); every |pull| < 1.3σ; no flag. Shuffled-z null:
median max|pull| 3.8σ (R LOO) and 15σ (D_dt LOO) vs observed 0.84σ and 1.29σ. Injections (D_dt × f, f = 0.15–8, 21
values per lens, baseline pull included) give the factor at which the injected |z| crosses 5.53σ: R prior (A) ×0.21–0.43
down / ×2.0–6.8 up; R LOO (C) ×0.22–0.39 / ×2.5–7.9; D_dt LOO (D) ×0.54–0.78 / ×1.29–1.85 (per lens in
`results/d1_distance/lens_summary.json`). So only a D_dt/D_d inconsistency of a factor ≳ 2–8 (R) would be flagged;
D_d errors are 15–30 %. FRBs: no flag (threshold 5.81σ one-sided). Lowest: FRB 20220319D, DM 111 < its NE2001 DM_ISM
127, 3.26σ low (a low-latitude sightline; NE2001 is the likely overestimate). Highest: FRB 20190520B, 2.37σ (known
large DM_host). Injections just beyond each burst's limits flag 94/94 low and 94/94 high, with no wrong-side flags. The
low limit is above DM_ISM for most bursts (median 2.3× DM_ISM) because Macquart's p(Δ) has a sharp lower cutoff
(Δ ≳ 0.4); that cutoff is an ASSUMPTION and sets the low-side sensitivity.

**Revisit if.** Per-lens D_d posteriors for the other TDCOSMO 2025 lenses (or the 2025 λ_MST-free kinematics) are
released; JWST/KCWI spatially resolved kinematics shrink D_d errors below ~10 % (the R test then reaches ×1.5);
a YMW16 or empirical Galactic DM model replaces NE2001 at low latitude.

**Addendum (2026-10-09): YMW16 beside NE2001.** `scripts/d1_distance.py frb --ism ymw16` replaces the stored NE2001
`DMISM` by YMW16 (Yao, Manchester & Wang 2017) to 30 kpc at each burst's RA/Dec, everything else unchanged (±20 %,
halo, cosmic, host, 5.81σ). YMW16 comes from pygedm 3.3.0's compiled `ymw16` extension, called directly: pygedm's
NE2001 extension needs libf2c (not in the cloud image) and its package import fails on SciPy ≥ 1.14 (`simps`), so the
package is built without `ne21c` and only `ymw16` is imported (recipe in `results/d1_distance/README.md`). Known-case
check: FRB 20121102A gives 287 pc cm⁻³ (the published YMW16 value; NE2001 stored 158). **Null under YMW16 and under
either model:** no flag; max low pull FRB 20220319D 4.09σ (YMW16 211 vs NE2001 127 vs DM_obs 111), max high
FRB 20190520B 2.37σ; injections 94/94 per side. YMW16/NE2001 ratio median 0.84 (16–84 %: 0.71–1.31, range 0.36–1.91),
so the models differ by more than the ±20 % ISM error (ASSUMPTION) for 54 of 94 bursts; the ±20 % is therefore too narrow
as a model error. Widening it only broadens the predictive distribution and lowers the pulls, so the null stands, but
the per-burst low-side DM limits under either single model are optimistic; a burst
counts as flagged if either model flags it (`frb_ism_compare.ecsv`, `flag_either`). The 20220319D pull is model error
on a low-latitude sightline (both models exceed its total DM), not a sightline anomaly. Rejected: building NE2001 in
pygedm (needs a system f2c library; the stored NE2001 values are already the reference).

**Addendum 2 (2026-10-09): statistic D on the TDCOSMO 2025 power-law chains, 8 lenses.**
`scripts/d1_distance.py tdcosmo` runs the D_dt leave-one-out (H0 + Ωm free, flat ΛCDM) on the power-law,
kinematics-free D_dt^model chains behind the TDCOSMO 2025 likelihoods: the 6 H0LiCOW lenses plus DES0408 and WGD2038
(file pairing and readers as in TDCOSMO's `tdcosmo_sample.ipynb`; redshifts from `tdcosmo_sample.yaml`).
- κ_ext harmonised: every lens gets D_dt = D_dt^model/(1 − κ_ext) (hierArc convention, λ_int = 1) with an independent
  draw from its own TDCOSMO κ_ext distribution (`kext`); the `nokext` variant sets κ_ext = 0 for all as a diagnostic.
  The mass sheet is treated the same for all (power law, no internal MST; a sample-wide factor is absorbed by H0).
  WGD2038's chain is weighted (weights applied); DES0408 has no κ_ext in its chain, so its own PDF is applied like
  the others'.
- SDSS1206: the pre-LOS power-law pickle (D_d D_s/D_ds, D_d, κ_pert) is read with an unpickler that admits numpy
  arrays only (`distance_consistency.load_array_pickle`); D_dt = (1 + z_d) × first / (1 + κ_pert), as in the notebook.
- Sidak over the 8 `kext` pulls: 5.39σ local (ASSUMPTION; `nokext` is a diagnostic and never flags). Null: 100
  random redshift permutations. Injection: D_dt × f, f = 0.15–8 (baseline = the f = 1 run under the same settings).
- Checked: the RXJ1131, PG1115 and HE0435 κ_ext files have as many rows as their D_dt chains but are uncorrelated
  with them (|r| < 0.003; ln D_dt width identical paired or shuffled), so independent draws (TDCOSMO's own
  treatment) lose nothing; B1608's chain κ column is all zero. The WGD2038 histogram range [−0.2, 1], 2000 bins, is
  the one in TDCOSMO's `plotkappa_handpickedpaper2038.py`; its narrow κ_ext (16–84 %: −0.007 to 0.010) is as published.
- **Result, `kext`: null.** Max |pull| 1.12σ (RXJ1131 −1.12, PG1115 −1.05, B1608 +0.97; the rest < 0.7σ); LOO H0 of
  the others 72.2–74.4. Shuffled-z null median max|pull| 15.8σ. Detectable factors ×0.41–0.78 down / ×1.36–1.85 up
  for the seven tight lenses, ×0.26 / ×3.2 for WGD2038 (σ_ln D_dt ≈ 0.24).
- **`nokext`: RXJ1131 −5.67σ (beyond the threshold), DES0408 +3.9σ.** Explained by line-of-sight convergence: removing the
  measured κ_ext (RXJ1131 median 0.070, B1608 0.103, WFI2033 0.059, DES0408 −0.040) moves the others' H0 to 77.5–79.8
  while RXJ1131, the most precise chain (σ_ln = 0.023), stays at 73.4. With κ_ext the tension is gone (−1.1σ). This
  is a real-data positive control: the statistic resolves a ~7 % sightline convergence on the precise lenses, so an
  exotic per-sightline distance effect would have to be smaller than the κ_ext uncertainties (±3–5 %) or mimic κ_ext.
- Limitations: power-law models only (no composite); no kinematics, so no per-lens D_d (new D_d only via the hierArc
  `*_const_processed.pkl` likelihoods, not read: pickles from a download need a data-only reader first).

**Addendum 3 (2026-10-10): composite-model check on SDSS1206.** TDCOSMO 2025 publishes a composite-model D_dt
chain for SDSS1206 only (`final_composite_D_dt.npy`, κ_ext and κ_pert included; the other seven lenses have power-law
chains only), so the model-choice check is one lens. `scripts/d1_distance.py tdcosmo --j1206 final_composite` swaps it
into statistic D (`kext` only, everything else as addendum 2).
- Validation: `--j1206 final_power_law` (TDCOSMO's own final power-law chain) reproduces this script's pre-LOS power
  law + κ_ext treatment: D_dt 16/50/84 % 5231/5881/6540 vs 5236/5884/6540 Mpc; J1206 pull 0.606σ vs 0.608σ. In
  these runs the J1206 `kext_p*` columns are 0: κ_ext is inside the chain.
- **Result: null.** Composite D_dt median 5708 Mpc (−2.9 % against the power law, narrower); J1206 pull 0.54σ (power
  law 0.61σ); max |pull| 1.14σ (RXJ1131) against 5.39σ; LOO H0 of the others 72.2–74.1. Detectable factors for J1206:
  ×0.62 / ×1.49 (power law ×0.50 / ×1.74).
- Rejected: the pre-LOS composite pickle (`angular_diameter_pre_LOS_composite.txt`): its κ_pert array has 20,000
  rows against 400,000 distance rows, so the per-sample pairing the power-law reader uses is undefined there.

**Addendum 4 (2026-10-10): kinematic D_s/D_ds per lens, 76 lenses.** `scripts/d1_distance.py kinematic` reads the
hierArc kinematic likelihoods TDCOSMO 2025 ships as pickles (8 TDCOSMO lenses; SLACS KCWI 13, SLACS SDSS 41 not in
KCWI, SL2S 14: `ExternalLenses/`) with `distance_consistency.load_data_pickle` (numpy-only globals) and turns each
into a likelihood for D_s/D_ds, the H0-free ratio R of this record, without any D_dt.
- Model (hierArc 1.2.0 `KinLikelihood` + `_displace_lambda_mst`, read, not installed): σ_v = c √(J s(a_ani[, γ_pl])
  (D_s/D_ds) λ), covariance C_meas + C_√J s (D_s/D_ds) c², λ = λ_int (1 − κ_ext). Here λ = 1, a_ani flat over the
  published grid (31 nodes), γ_pl (RXJ1131 and the 13 SLACS KCWI lenses) from each lens's Gaussian prior truncated
  to the published grid [1.5, 2.5] (41 nodes; hierArc bounds it the same way; `gamma_prior_in_grid`, lowest
  SDSSJ0029−0055 0.23, SDSSJ1402+6321 0.59); flat prior in ln D_s/D_ds on [0.04, 40] (ASSUMPTIONs, all in
  `Params.kin_*`). The term is Gaussian in σ_v with an error proportional to the prediction, so its upper
  tail is a power law: the grid end is a prior bound and per-lens summaries are median and half 16–84 % width.
  Moving the bound to 400 (`--ratio-max 400`) changes the max pull by 2 × 10⁻⁴σ and δ, τ by < 10⁻³.
- Statistic: y = ln(D_s/D_ds)_kin − ln(D_s/D_ds)_ΛCDM; leave-one-out pull of each lens against the others' common
  offset δ (absorbs a population λ) and Gaussian intrinsic scatter τ (ML; absorbs per-lens λ_int, anisotropy and
  profile scatter), tail computed on the lens's own grid density. Ωm 0.3 (0.1 and 0.5 as checks). Sidak over 76:
  5.78σ local. Null: 100 shuffled redshift pairs. Injection: one lens's D_s/D_ds × f.
- **Result: null.** Max |pull| 2.49σ (SDSSJ2302−0840, SLACS SDSS, −0.63 in ln), 2.63σ / 2.54σ at Ωm 0.1 / 0.5; next
  SDSSJ1538+5817 (KCWI, +0.29) +2.48σ. Shuffled-z null median max|pull| 2.71σ (p ≥ observed 0.79): the redshift lever is weak against τ = 0.15,
  so this tests per-lens outliers, not the redshift dependence. Pooled δ = −0.06 ± 0.03. Detectable factors (median
  over lenses): × 0.21 down / × 5.6 up; TDCOSMO × 0.19–0.39 / × 2.9–9.0, SLACS KCWI × 0.30–0.45 / × 1.6–2.9,
  SLACS SDSS × 0.02–0.42 / × 3.7–52, SL2S × 0.02–0.22 / × 3.3–32 (two SL2S lenses not reached down, one up, within
  × 0.018–55).
- Sample offsets (δ ± sd, τ): KCWI +0.06 ± 0.04 (0.12), SLACS SDSS −0.15 ± 0.04 (0.03), TDCOSMO −0.11 ± 0.05
  (0.0), SL2S −0.36 ± 0.15 (0.35). KCWI vs SDSS differ by 0.2 on the same SLACS lenses (aperture and data
  systematics; hypothesis, not tested). With λ = 1 (diagnostic `z_lambda1`, never flags) the precise KCWI lenses sit
  up to 10σ high: λ_int ≠ 1 per lens, the mass-sheet degeneracy this test cannot break.
- Cross-check, not a validation: the TDCOSMO 2025 kinematic ratios differ from the H0LiCOW joint R (composite
  models, other anisotropy priors) by −0.43 to +0.30 in ln (PG1115 1.04 vs 1.61; J1206 2.48 vs 1.83), within
  their 1–1.6σ widths. Synthetic lenses recover the injected ratio to 0.01 in ln (tests).
- Limitations: λ and anisotropy are marginalised only through δ, τ and flat grids, so a real per-sightline distance
  anomaly smaller than ×2–5 is indistinguishable from a mass-profile difference; the Gaussian summaries of the
  others are approximate for the heavy-tailed SL2S lenses.

## D-074 E1 causal event network: lag × separation pair counts across GBM, ICECAT-1, GWTC and CHIME Cat 2 with sidereal-scrambled nulls (2026-10-09)

**Decision.**
- Hypothesis (owner idea 4): events from different directions are dependent at lags no ordinary path explains.
  Nodes:
  - 4,390 GBM bursts (`fermigbrst`);
  - 340 ICECAT-1 v4 tracks (8 `CR_VETO` dropped);
  - 391 GWTC events;
  - 3,641 CHIME/FRB Cat 2 sources (sub-bursts dropped, one node per repeater).

  All are pinned in `data/manifests/e1_events.ecsv`. The CANFAR download worked on the first try this time.
- Statistic: pair counts per channel (9 channels) × lag bin (0–10 s, 10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d, |Δt|) ×
  class.
  - "Same" means sep ≤ 3·σ_comb. "Wide" means sep > 3·σ_comb and > 0.1°.
  - σ_GBM = stat ⊕ 3.7°. Below an `error_radius` of 0.5° (a position from another instrument, usually 0) it is
    max(error_radius, 0.05°).
  - σ_ICECAT = mean 90 % error / 2.146. σ_CHIME = max(ra_err, dec_err).
  - GW events have no position in the CSV, so GW channels count all pairs (lag only).
- Nulls, 10⁴ each. All keep Dec and hour angle and permute times within catalogue × calendar year.
  - `perm`: the time multiset is unchanged, so this null tests only whether separation depends on lag.
  - `jit` (primary for lag excess): adds an independent shift per event. GBM gets k × 95.6 min orbit ± 5 min within
    ±3 d, which keeps orbit phase (SAA, occultation); the others get uniform ±3 d.
  - `jitday` (vetting): the ground instruments and GW shift by whole days.
- Trials (review finding 1 on PR #112): the observation is pooled with the scrambles. Every row's per-cell p is
  (#rows ≥ value)/(N + 1), and the global p is the rank of the observed minimum among all rows' minima. This is
  uniform under H0 (tested). The empirical floor is 1/(N + 1) = 10⁻⁴, which makes the family-wise floor about
  (cells at the floor)/(N + 1). A family-wise 3σ or 5σ is therefore not reachable empirically (`reachable()`,
  tested). Such claims use `analytic_p` (Poisson for null mean < 30, else Gaussian with the ensemble sd) ×
  Bonferroni, a model_prediction, quoted beside the pooled empirical p.
- Injection (review finding 2): a B event is moved to t_A ± lag (log-uniform in the bin) of a wide-separated A anchor
  within ±30 d.
  - Detection: the injected cell's analytic p × 75 cells ≤ 1.35 × 10⁻³ (family-wise 3σ), i.e. per-cell p ≤ 1.8 ×
    10⁻⁵. The earlier empirical rule (p ≤ 10⁻⁴) was only about 2.7σ family-wise.
  - Limits: the Poisson classical UL on the count when the null mean is < 30, else obs − q05(null). The rate is UL /
    min(eff, 1) / N_eligible.
  - N_eligible counts the A events with a B event within ±30 d, i.e. those inside B's live time (review finding 4).
  - A cell where no injected n reached 50 % detection has `valid` = False and no rate.
- Code: `src/jwst_anomaly/event_network.py`, `scripts/e1_events.py`, `scripts/e1_chime_exposure.py`, and
  `tests/test_event_network.py`. All thresholds are ASSUMPTIONs in `event_network.Params`.

**Alternatives rejected.**
- Reusing published coincidence-search code. Curtin+2023 (CHIME × GBM/BAT) and Masaoka+2026 (CHIME Cat 2 ×
  ICECAT-1) publish no code (checked 2026-10-09). Both test *same-direction* coincidences, which is the opposite of
  E1's wide channel. IceCube/GBM stacking analyses are likewise position-matched. Brainerd+1995 (BATSE
  time-dependent two-point correlation) is the closest method, but it is a same-direction repeater test. The new
  quantity is the wide-separation lag CCF, which is a few lines of numpy. No new dependency (`h5py` is used only
  for the one-off exposure reduction, from a scratch environment).
- A bootstrap with replacement for the uptime null: duplicated draws make self-pairs at lags < 2·jitter. Rejected.
- `perm` alone: it preserves every lag, so GW170817 × GRB 170817A gets p = 1 by construction.
- Ranking the observation with (1 + r)/(N + 1) but the scrambles with r/N. This dropped the trials factor whenever
  the observation lay beyond the ensemble: it reported 10⁻⁴ where the pooled p is 0.0036. Fixed after review.
- The in-day null (keep each event's UTC day, redraw the time of day) as an uptime test: it absorbs most of an
  injected real signal (`chime_vet.json`: 300 injected 1 h–1 d wide pairs give z = 5.6 under `jit` but 0.25 under
  in-day), so it cannot discriminate uptime from dependence.
- The CHIME Cat 2 exposure file as an uptime model. It holds two HEALPix nside-4096 maps of *time-integrated*
  exposure (upper and lower transit, 2018-09-04 to 2023-09-15) and has no time axis. It cannot give a per-day uptime
  series, and its sky part (Dec dependence) is already preserved by the Dec + hour-angle nulls. It was reduced to
  `chime_exposure_dec_profile.ecsv` and deleted. No CHIME Cat 2 file with a time-resolved uptime or injection
  series exists in the CANFAR release (checked: `table/`, `exposure/`, `localizations/`, `dynamic_spectra/`,
  `additional_figures/`).

**Evidence.** See `results/e1_events/`.
- Requested five channels (GBM–GBM, GBM–ICECAT, GBM–GW, ICECAT–ICECAT, ICECAT–GW), `jit`: min cell p = 0.014
  (ICECAT–ICECAT 1 h–1 d wide, 36 vs 25.4 ± 4.4). Trials-corrected p = 0.28. Null.
- All 9 channels, pooled empirical trials-corrected p:
  - all 75 cells: 0.0036 (`jit`) and 0.0047 (`jitday`);
  - wide and GW cells (45): 0.0058 (`jit`) and 0.0033 (`jitday`).

  Beside them, the analytic Bonferroni p over the 45 wide/GW cells is 4.5 × 10⁻³ (`jit`) and 4.5 × 10⁻⁴ (`jitday`,
  3.3σ). All 75 cells give an analytic 4 × 10⁻¹⁰, but that is CHIME–CHIME *same-direction* (see below).
- **CHIME–CHIME wide (and all-separation) clustering at 100 s–1 d.**
  - Status: unexplained under `jit` and `jitday`; the in-day null does not discriminate. **It is reproduced by a
    calibrated rate-modulated null** (test d below), so it is attributed to smooth (~week-scale) modulation of the
    CHIME detection rate. The cause of that modulation (uptime or sensitivity) is a hypothesis: no CHIME uptime
    series exists to confirm it.
  - Counts: 8,099 vs 7,793 ± 82 at 1 h–1 d (z = 3.7 `jit`, 4.3 `jitday`); 396 vs 329 ± 18 at 100 s–1 h. It is
    direction-independent: the wide fraction under `perm` is normal (z_perm = −0.7 and −2.5), and the
    all-separation counts show the same excess.
  - In the `excluded_flag` = 0 subset (post hoc, 2,000 scrambles) the per-cell z is 5.3 (`jit`) and 5.8 (`jitday`)
    at 1 h–1 d. The Gaussian Bonferroni over 45 cells is 5.2σ, or 4.96σ including the three post-hoc subsets.
  - Ordinary-explanation tests on the 1 h–1 d wide cell (`chime_flag_tests.json`, `scripts/e1_chime_flag.py`,
    1,000 scrambles unless stated):
    - (a) Busy days.
      - The 14 busiest days (top 1 %, 107 bursts) touch 578 cell pairs against 206 in scrambles. Dropping them
        lowers z from 3.9 to 2.9. Dropping the top 5 % lowers it to 2.1, and the top 10 % to 2.8.
      - A 30-day-block jackknife (63 blocks) keeps z between 3.4 and 4.8. No single month carries the excess.
    - (b) Epochs and seasons. The excess is present in most years (z = 3.0, 2.7, 1.5, 0.2, 2.2 and 3.0 for 2018–2023)
      and in both the Catalog 1 period (`catalog1_flag`, z = 3.2) and later (z = 2.9). It is strongest in MAM
      (z = 3.9; DJF 1.2, JJA 2.2, SON 2.7). It is not a commissioning-only effect.
    - (c) Property independence.
      - Cell pairs do not share DM, fluence or S/N more than scrambled pairs do (KS p = 0.55, 0.98 and 0.47).
      - Dec difference: p = 0.04, but the median differs by only 0.3°, and with 5 properties tested it is not
        significant.
      - Sidereal-phase difference: p = 0.09.
    - (d) Rate-modulated null. Times are redrawn from the catalogue's own 7-day running mean of daily counts
      (uniform within the day; Dec and hour angle kept).
      - Calibration: 300 injected 1 h–1 d wide pairs raise the observed count by 190, of which the null absorbs
        24. It keeps about 87 % of an injected dependent signal.
      - Data: z = −0.86 (8,099 vs 8,209 ± 128). The jit excess disappears.
      - A 3-day running mean over-predicts pairs (z = −5.5) because daily counts are under-dispersed (Fano 0.86).
  - Under the brief's STOP rule (> 5σ trials-corrected under a calibrated rate or exposure null) the full-catalogue
    cell has z = −0.9 under the calibrated null (d). No STOP.
- CHIME–CHIME same-direction at 10 s–1 h (12 vs 0.6 pairs at 100 s–1 h). 11 of the 12 pairs are one unflagged
  same-position episode: FRB20230825D–I, six bursts within 5 min at (347.35°, +48.75°) with DM 221–223 pc cm⁻³ and
  no `repeater_name`. This is a catalogue effect (one source, six nodes).
- Positive controls.
  - GW170817 × GRB 170817A: the only GBM–GW pair at 0–10 s (Δt = 2.07 s trigger − merger; the published onset is
    1.74 s). The null mean is 0.06, so the blind count gives p = 0.058: a single coincidence is not detectable
    without sky maps. 50 % detection at family-wise 3σ needs ≥ 3 injected pairs. The GBM row's position lies 1.2″
    from SSS17a.
  - GBM re-triggers (GRB 091024, 130925, 150201, 220627, 250702 ×3) appear in the same-direction channel.
  - CHIME with every repeater burst as a node: same-direction pairs 30 / 120 / 219 / 480 / 2,513 vs perm
    1.6 / 6.8 / 13.6 / 33 / 195.
- Orbit, sidereal-day and solar-day lag windows (k = 1–15 orbits, 1–7 days) show no excess (|z| < 2.5).
- 95 % limits, dependent wide (GW: any-separation) partners per eligible anchor, `limits.ecsv`:
  - lags ≤ 100 s: 6.8 × 10⁻⁴ (GBM–GBM), 1.0 × 10⁻³ (GBM–ICECAT, 2,910 eligible), 5.1 × 10⁻³ (GBM–GW, 935
    eligible), 8.8 × 10⁻³ (ICECAT–ICECAT), 4.5 × 10⁻² (ICECAT–GW, 66 eligible), and 0.8–2.8 × 10⁻³ (CHIME
    cross-channels);
  - longer lags: 4 × 10⁻³ to 0.75.

  ICECAT–ICECAT 100 s–1 h is not valid: no injected n reached 50 % detection.

**Revisit if.**
- GW sky maps are added (then GW channels get same/wide classes and GW170817 enters the same-direction channel).
- A time-resolved CHIME uptime or sensitivity series is published or obtained. It would confirm (or refute) that the
  rate modulation behind the CHIME flag is instrumental; an exposure-weighted null built from it must keep an
  injected signal. Sub-day (hour-resolved) uptime would also decide the 100 s–1 h residual (addendum).
- A CHIME Cat 2 revision assigns FRB20230825D–I to a repeater.
- IceTracks-DR2, Swift or Einstein Probe catalogues are added.
- Signed-lag (precursor) channels are wanted.

**Addendum (2026-10-09): the CHIME–CHIME 100 s–1 h wide cell under the rate-modulated null.**
- Hypothesis: the 100 s–1 h wide excess (396 vs 329.6 ± 17.8 under `jit`, z = 3.7) is lag dependence that a smooth
  detection-rate modulation does not explain. Ordinary explanations: week-scale rate modulation (as for 1 h–1 d),
  the daily duty cycle, and sub-day outages or sensitivity changes. A family-wise > 3σ residual under the calibrated
  null would have sent the cell to per-pair vetting.
- `scripts/e1_chime_flag.py --cell 100s-1h`, 2,000 scrambles; `results/e1_events/chime_rate_null_100s_1h.json`.
  The same D-074 (d) null, plus a `keep_tod` variant that keeps each event's UTC time of day.
- **Results (z; analytic p × 45 cells):**
  - 7-day running mean: 396 vs 348.9 ± 18.7, z = 2.52, p_fw = 0.26.
  - 7-day, time of day kept: 396 vs 349.6 ± 19.2, z = 2.42, p_fw = 0.35 (the daily duty cycle adds nothing).
  - 3-day running mean: z = 0.48. Not used as the primary null: it over-predicts the 1 h–1 d cell (z = −5.5, D-074 (d)).
- Calibration: 300 injected wide pairs in this bin raise the count by 300–305; the null mean moves by ≤ 8, so the
  7-day null keeps ≥ 97 % of an injected dependent signal (1 h–1 d: 87 %).
- **Decision: null at the family level.** The week-scale modulation explains about 30 % of the jit excess (19 of 66
  pairs); the rest is a local 2.5σ that a day-resolution null cannot test, because it needs hour-scale rate
  structure. Unlike 1 h–1 d (z = −0.9), the cell stays a weak, unexplained-at-low-significance residual; no STOP.
- Alternatives rejected: an hour-resolved running mean of the catalogue's own counts. With about 2 bursts per day
  it would absorb the very pairs it tests (as the in-day null did, `chime_vet.json`).
- Revisit if: a time-resolved CHIME uptime series is obtained (owner decision); then use an hour-resolved
  exposure-weighted null for both 100 s–1 h and 1 h–1 d.

**Addendum 2 (2026-10-09): antipodal lag channels.**
- Hypothesis (owner idea 4 follow-up): dependent events appear at the antipode of a first event (a "through the
  Earth / through the sky" link). Ordinary explanations: the same observer-side common causes as D-074; a GBM
  localization flip is not a concern because GBM errors set the antipodal radius.
- `scripts/e1_antipodal.py`: pairs with 180° − sep ≤ max(10°, 3σ_comb) (ASSUMPTION `ANTI_MIN_DEG`), five lag bins
  × six localized channels (GW has no positions) = 30 cells; `jit` null, 2,000 scrambles; pooled trials
  (`en.global_p`) and an analytic Bonferroni tail. Injection: B moved to t_A ± lag and to A's antipode (scattered by
  B's own error), 20 trials at n = 3, 10, 30, 100 per cell; detection = family-wise 3σ.
- **Result: null.** Pooled global p = 0.56 (min cell p 0.054); min analytic p × 30 = 1. Largest z: CHIME–GBM
  1 h–1 d, 258 vs 230.5 ± 16.4 (z = 1.7). `results/e1_events/antipodal.json` gives the per-cell 95 % upper limits.
- Sensitivity (`n50_injected`): 3 pairs at lags ≤ 100 s in most channels (10 for GBM–GBM 10–100 s and the CHIME–GBM
  10–100 s cell); 3–30 at 100 s–1 h; 10–100 at 1 h–7 d; GBM–GBM 1–7 d not reached at 100.
- Limitation: the CHIME–CHIME injections put B at A's antipode even where CHIME cannot see it (Dec < −11°). They
  calibrate the statistic, not the instrument.
- Revisit if: GW sky maps are added (GW–X antipodal channels), or signed-lag channels are built.

**Addendum 3 (2026-10-09): signed-lag ("which event comes first") channels.**
- Hypothesis: one catalogue's events lead another's (A → B) at lags no ordinary path explains, visible as a sign
  asymmetry that the pair counts (|Δt|) cannot see. Ordinary explanations: physical counterparts with a known order
  (GW170817 → GRB 170817A), follow-up chains, and uptime (sign-symmetric, so it cannot make D ≠ 0 on its own).
- `scripts/e1_signed_lag.py`: D = N(t_B > t_A) − N(t_B < t_A) for the six cross channels × five lag bins × class (same
  and wide; all for GW channels) = 45 cells; `jit` null, 2,000 scrambles; pooled two-sided trials (`en.global_p` on
  |D − μ|). Analytic tail: Gaussian, or Skellam(λ, λ) with λ = sd²/2 where the null sd < 1 (ASSUMPTION
  `GAUSS_SD_MIN`). The Gaussian on a one-pair cell gave z = −4.4 for GW170817 alone and was rejected. Injections
  move B to t_A + lag (one sign) in the tested class.
- **Result: null.** Pooled global p = 0.68 (min cell p 0.032); min analytic p × 45 = 1.
  - Largest |z|: GBM–GW 0–10 s, D = −1. This is GW170817 before GRB 170817A, the ordinary positive control (Skellam
    p = 0.05, not detectable alone).
  - Next: CHIME–GW 1 h–1 d, D = −45 vs +3.7 ± 22.6 (p = 0.03).
- Sensitivity (`n50_injected`): 3–10 one-sided pairs at ≤ 100 s, 10–100 at 100 s–1 d (2 cells not reached), and mostly not reached at
  1–7 d (≥ 100).
- Revisit if: GW sky maps are added (GW cells get same/wide classes).

**Addendum 4 (2026-10-10): GW channels with sky maps.**
- Hypothesis: GW events have same-direction, wide or antipodal partners (GBM, ICECAT-1, CHIME, GW) at lags no
  ordinary path explains. Before this the GW channels were lag-only. Ordinary explanations: as D-074, plus
  GW170817 × GRB 170817A (a known counterpart).
- Data: 282 PE sky maps (GWTC-2.1 / 3 / 4.1 / 5.0 Zenodo tarballs, 704 MB streamed in about 3.5 min with 12 range
  connections per tarball, one process per `--tar`; reduced to nside-32 NESTED; `scripts/e1_gw_skymaps.py`).
  The other 109 GWTC entries (GWTC-4.1: 52, GWTC-5.0: 57) have no PE parameters in the GWOSC CSV either and are
  left out. GW170817: Gaussian at SSS17a, σ 1.05° (ASSUMPTION).
- `scripts/e1_gw_directional.py`: same = partner within 3σ (+1.3° pixel slop, ASSUMPTION) of the 90 % region;
  wide = more than max(3σ, 0.1°) outside the 99 % region; antipodal = antipode within the same tolerance of the
  90 % region and not same; GW–GW by region overlap. Maps rotate in RA with the scrambled GW time (hour angle
  kept). 4 channels × 5 lags × 3 classes = 60 cells, `jit` null, 1,000 scrambles (266 s on 4 cores).
- The control pair GW170817 × GRB 170817A is classified `same` at 0–10 s and is left out of the family (both
  observation and null; with it the 0–10 s GW–GBM same cell would be 1 vs ≈ 0). The catalogue snapshot is pinned
  separately in `data/manifests/e1_gw_events.ecsv`, because the GBM TAP table grew since the D-074 pin.
- **Result: null.** Pooled global p = 0.83 (min cell p 0.053); min analytic p × 60 = 1. Largest z: GW–ICECAT
  1–7 d antipodal, 14 vs 8.3 ± 2.9 (z = 2.0); GW–GW 1 h–1 d antipodal, 17 vs 11.1 ± 3.4 (z = 1.8).
- Sensitivity (`n50_injected`, 10 trials, family-wise 3σ): 3 pairs at ≤ 10 s (GW–GBM, ICECAT, CHIME), 3–10 at
  10–100 s, 10–30 at 100 s–1 h; at 1 h–7 d mostly > 30 (not reached). GW–GW is not injected. 95 % upper limits per
  cell are in `results/e1_events/gw_directional.json` (3 extra pairs in every empty cell).
- Revisit if: new GWTC PE releases (more mapped O4 events), or signed-lag GW cells with directions are built.

**Addendum 5 (2026-10-10): signed-lag GW cells with sky-map classes.**
- Hypothesis: as addendum 3, with GW cells split by the addendum 4 sky-map classes. Ordinary explanations: known
  counterparts with a fixed order (GW170817 → GRB 170817A, left out of the family); catalogue edges (GW events
  after CHIME Cat 2 ends have only earlier partners). The per-year `jit` scramble keeps the edges, so its mean D is
  not 0 (GW–CHIME 1–7 d wide: −66.5 ± 52.7) and every cell is tested against the null mean, never against 0.
- `scripts/e1_gw_signed.py`: D = N(t_B > t_GW) − N(t_B < t_GW) for GW × {GBM, ICECAT, CHIME} × 5 lags × {same, wide,
  antipodal} = 45 cells; `jit` null, 1,000 scrambles; two-sided pooled trials and the addendum 3 Gaussian / Skellam
  tail. GW–GW is not a cell (D is antisymmetric within one catalogue). Injections put B after the GW event in the
  tested class (`e1_gw_directional.inject(..., sign=1)`), 10 trials at n = 3, 10, 30.
- **Result: null.** Pooled global p = 0.56 (min cell p 0.030); min analytic p × 45 = 0.78. Largest |z|: GW–CHIME
  1 h–1 d antipodal, D = 10 vs −0.3 ± 4.3 (z = 2.4).
- Sensitivity (`n50_injected`): 3 pairs at ≤ 10 s, 3–10 at 10–100 s, 10–30 at 100 s–1 h, and mostly not reached
  (> 30) at 1 h–7 d.
- Revisit if: new GWTC PE releases, or a time-resolved CHIME uptime series (owner decision).

**Addendum 6 (2026-10-18): Swift BAT GRBs as a fifth node set (`results/e1_swift/`).**
- Input: HEASARC `swiftgrbba` (burst-advocate compilation; the older `swiftgrb` table stops in 2012), on-board
  Swift triggers only, 1σ = BAT 90 % radius / 2.146. All five inputs re-pinned in `data/manifests/e1_swift.ecsv`.
- Family: BAT–BAT, BAT–ICECAT, BAT–CHIME (same, wide), BAT–GW (all) × 5 lags = 35 cells; BAT–GBM is the positive
  control (466 same-direction pairs at 0–100 s vs 0.10 expected: passes). BAT jit quantum = 5,700 s (ASSUMPTION).
- **Result: null.** Min cell BAT–CHIME 10–100 s wide (analytic p 0.047, × 35 = 1); pooled global p 0.54. Limits
  per eligible anchor in `results/e1_swift/limits.ecsv` (e.g. ≤ 1.8 × 10⁻³ dependent wide BAT–BAT pairs per burst
  at ≤ 100 s).
- Control-channel findings (post hoc): the 0.05° σ floor puts 8 same-GRB BAT–GBM pairs in "wide"; a BAT–GBM pair
  deficit at 1 h–4 d is carried by the 466 joint bursts (removed: −0.5σ), mechanism not determined. A null that
  keeps Sun-relative positions does not remove it.
- Revisit if: an XRT-position variant is needed for same-direction cells, or a satellite × satellite channel
  that contains joint detections enters a tested family (then remove joint events first).

**Addendum 7 (2026-10-18): Einstein Probe transients as a sixth node set (`results/e1_ep/`).**
- Input: no public EP event table exists; one row per EP name is parsed from its first WXT report in the GCN
  Circulars JSON archive (`scripts/e1_ep.py`; onset or first timestamp, J2000 position, 90 % radius / 2.146;
  trigger called "not a real source" dropped; streamed, not stored). 231 events. Rejected: tables in papers (PDF, curated subsets).
- Only GBM and BAT overlap EP in time. Family: EP–EP (5 lags) and EP–BAT, EP–GBM (1 h–1 d, 1–7 d), same and wide:
  18 cells; EP–GRB same-direction pairs at 0–1 h are the positive control (27 vs 0.63: passes). EP jit quantum
  5,760 s (ASSUMPTION).
- **Result: null.** Min cell EP–GBM 1 h–1 d same (p 2.2 × 10⁻³, × 18 = 0.039); pooled global p 0.017. 5 of its 27
  pairs are one source twice (EP250702a = GRB 250702B–E). The plan did not drop joint detections first
  (deviation from addendum 6); the post hoc re-count without them gives × 18 = 0.015. `jitday` not run: with no ground catalogue it equals `jit`.
- Revisit if: EP publishes an event catalogue (then replace the parsed table), or an EP–GBM test is planned (then
  define joint events by name and position, not by |lag| < 1 h).

**Addendum 8 (2026-10-18): IceTracks-DR2 northern tracks as a seventh node set (`results/e1_icetracks/`).**
- Input: the D-079 pinned DR2 files; Dec > −5°, log10 E ≥ 4.0, one track per (run, event), 1σ = AngErr: 7,583
  tracks (2008-04 to 2022-05). Null for DR2: the D-079 uptime-aware ±3 d jitter; other catalogues as addendum 6.
- Family: DR2–GBM, DR2–BAT, DR2–CHIME (same, wide), DR2–GW (all) × 5 lags = 35 cells. DR2–ICECAT same-direction
  pairs at 0–10 s are the positive control (270 one-to-one vs 0.0013 expected: passes). DR2–DR2 stays with D-079.
- **Result: null.** Min cell DR2–CHIME 1–7 d same (p 0.065, × 35 = 1); pooled global p 0.79. Limits per
  eligible anchor in `results/e1_icetracks/limits.ecsv` (e.g. ≤ 4.1 × 10⁻⁴ dependent wide DR2–GBM pairs per
  track at 0–10 s).
- Control-channel finding (post hoc): DR2–ICECAT wide pairs are short of the null at 1 h–1 d (z −4.5) and 1–7 d
  (z −6.5); six distinct null variants (downtime, DR2 end date, either side alone) leave it; cause not determined. No
  tested cell is below z −1.7.
- With this, every node set TASKS.md listed for E1 is done; E1 waits only on the CHIME uptime series (owner).
- Revisit if: f_astro per cut is computed (turns per-track limits into per-astrophysical-neutrino limits), or
  a southern / lower-energy DR2 variant is planned (atmospheric muons dominate; new control needed).

**Addendum 9 (2026-10-19): the DR2–ICECAT control deficit is the null displacing shared events (`results/e1_icetracks/selfpair/`).**
- **Decision.** In a channel whose two catalogues share events, the jit null moves each shared pair apart by up
  to a few days, so it over-predicts long-lag "wide" pairs. Shared events (here a DR2 track within 10 s of an
  alert; ASSUMPTION) are removed before such a channel's long-lag cells are read; the 0–10 s control cell keeps
  them. Addendum 6's "remove joint events first" now has this mechanism. Tested channels without shared events
  are unaffected; the shared null code is not changed.
- **Evidence.** Frozen rule "explained" (10³ scrambles per variant): dropping the 276 shared tracks gives wide z
  −1.24 (1 h–1 d) and +1.01 (1–7 d), from −4.53 and −6.48; dropping 276 random tracks leaves −4.63 and −6.98. The
  null mean falls by 122.38 / 450.22 pairs, the data by 60 / 240; the null-only parts (62.38 / 210.22) are within
  22 % / 11 % of the planned estimate (≈ 80 / 190).
- **Revisit if** a channel with shared events enters a tested family (remove them first), or the BAT–GBM deficit
  needs its own check.

**Addendum 10 (2026-10-19): the BAT–GBM control deficit at long lags behaves as the displacement predicts; frozen rule undetermined (`results/e1_swift/selfpair/`).**
- **Decision.** None on the tested family (no tested channel includes BAT–GBM). For a future joint-event rule
  (e.g. EP–GBM, addendum 7), pre-register a position match over a window longer than 100 s; the 1 h window
  used here is post hoc and is a hypothesis to freeze, not a settled parameter.
- **Evidence.** Frozen rule "undetermined": dropping the 476 BAT bursts with a GBM burst within 100 s gives
  all-separation z −0.33 (1 h–1 d) and −0.66 (1–7 d), from −4.78 / −5.04 (random drop: −5.83 / −5.31), but
  100 s–1 h z +4.2, where 25 same-direction pairs (23 within 0.1°, e.g. GRB 221009A at −3,197 s) look like joint
  GRBs with triggers 100 s–1 h apart. The post hoc 1 h variant is circular in that cell.
- **Revisit if** a channel with shared events enters a tested family.

## D-075 Hypothesis round 2: three worlds with derived Lorentz invariance; no survivor this project can test (2026-10-09)

**Decision.**
- Round 2 of the D-069 process ran with the round-1 failure modes as hard rules for System A.
- None of its predictions gets search compute here:
  - R2-A (perfect-reflector GW afterglow) is excluded by the ergoregion instability and the GWTC-4.0 remnant null,
    with the GW250114 reflectivity bound (preprint) as support.
  - R2-B (polarisation mirror at saddle images) is internally inconsistent and disfavoured by B0218+357.
  - R2-C (Σm_ν = 58.8 meV) is not novel. It is tracked externally through DESI DR3, KATRIN final and CMB birefringence.
- Round 3 targets a derived small effect in a large-N quantity that is measured but never analysed for it.

**Alternatives rejected.**
- Reprocessing the raw VLA data for B0218+357 now: R2-B already fails on internal consistency. Kept as an optional
  low-priority task.
- Re-running public DESI/CMB chains for R2-C: it would only reproduce the published collaboration results.
- A stacked GW echo search for R2-A: equivalent windows already tested (GWTC-4.0 cWB, O4a events, 2603.19021) and an
  energy-normalised bound set for GW250114 (2610.12429). R2-A's exact stack was not run; it is moot given the
  ergoregion exclusion.

**Evidence.** docs/hypotheses/round-2/ (texts, reviews, summary). Key references opened 2026-10-09:
- arXiv:2610.12429 (GW250114 near-horizon reflectivity);
- 2603.19021 (GWTC-4.0 remnant tests);
- 1706.06155 (echo recipe);
- 1802.10088 (B0218+357 VLA reanalysis);
- 2503.14744 (DESI DR2 neutrino constraints);
- 2605.21456 (DES-Dovekie: negative neutrino mass or negative dark energy).

**Revisit if.**
- Parity-resolved circular polarimetry of a lensed radio quasar is published, or a lensed FRB is confirmed (R2-B).
- DESI DR3 or a CMB analysis pushes the Σm_ν tension past 5σ, or relaxes it (R2-C).

## D-076 Hypothesis round 3: derived small effects in large-N data; no testable survivor; pause new rounds (2026-10-09)

**Decision.**
- Round 3 (brief: a derived small effect in a large-N quantity never analysed for it) gives nothing that earns
  search compute:
  - R3-A (universal same-exposure pair correlation) fails as written. Its decisive test is in a lab.
  - R3-B (condensate scalar charge) is mostly excluded by NS–WD dipole bounds (upper half at 6.7σ; self-consistent
    minimum at about 2σ; the no-feedback corner survives). Recorded: α_NS < 2.2 × 10⁻³ (95 %, reviewer's
    combination).
  - R3-C (SN Ia rate vs "clock depth") is not excluded but cannot be tested with about 5 × 10³ public SNe Ia under
    environment systematics.
- New System A rounds are paused. Compute goes to open searches limited by data (MOA-II, CHIME baseband, Euclid DR1,
  Gaia DR4) until the owner sets a new brief.

**Alternatives rejected.**
- A ZTF same-exposure correlation search for R3-A: about 3 TB to download, about 1.8σ at nominal D, and systematics bias it
  upward.
- A ZTF BTS Poisson GLM for R3-C: 0.5–2σ at N ≈ 5 × 10³; the kill needs N ≈ 8 × 10⁴.
- A stacked s(M) Ṗb template fit for R3-B: bound-setting only, and mass-resolved fits already exist.
- Continuing System A rounds at the same cadence: rounds 2–3 (6 worlds) produced no testable survivor.

**Evidence.** docs/hypotheses/round-3/ (texts, reviews, summary). Key references opened 2026-10-09:
- arXiv:1512.01216 (Holometer);
- 2605.01436 (J1738+0333 update);
- the cluster SN Ia rate papers listed in B-on-R3C.md.

**Revisit if.**
- The owner sets a new System A brief.
- Rubin/LSST public alerts give ≥ 10⁵ SNe Ia with hosts (R3-C).
- A lab cross-correlation of independent noise sources at 10⁻⁵ is published (R3-A).

## D-077 Owner brief 2026-10-10: Neutrino Frontier program; round 1 (four frameworks), data audit and E-NF1 ghost-pair null (2026-10-10)

**Decision.**
- Owner brief (DD-Ching, 2026-10-10) recorded here: a new, independent Neutrino Frontier program. Neutrinos are
  probes of emergent spacetime, locality and distance, mass and flavor, hidden degrees of freedom, causal
  structure, and wormhole or warp physics. The brief restarts hypothesis rounds **for this program only**; it
  supersedes D-076's pause here and nowhere else. It takes no resources from MOA-II, S1, E1/D1 follow-ups or
  data validation without a documented allocation decision. Merge policy, guarded files, WIP cap, cloud-disk
  and notification rules are unchanged.
- Layout: `docs/neutrino_frontier/` (README with the NF-H registry and experiments, round files, data audit);
  code in `src/jwst_anomaly/neutrino_frontier.py` and `scripts/nf_*.py`, reusing `event_network.py` and the
  E1 loaders and nulls.
- Round 1 (N-A literature-blind invention, N-B adversarial review with checked references):
  - NF-H01 granular locality (stochastic distance): KNOWN-REDUNDANT, UNTESTABLE NOW (needs ≥ 2 multi-event
    neutrino flares at different distances).
  - NF-H02 twin-sheet flavor = pseudo-Dirac partners: KNOWN-REDUNDANT; N-A's δm² window is excluded by public
    IceCube data (arXiv:2406.06476).
  - NF-H03 second causal cone = ν–γ Shapiro-delay EP test: KNOWN-REDUNDANT; robust |ε| ≲ 10⁻³–5 × 10⁻³
    (SN1987A); the LSS-scale bound is ill-defined (arXiv:1907.12453); D-074 GBM × ICECAT covers the rescan.
  - NF-H04 sparse nonlocal links ("ghost" neutrinos): CONDITIONALLY VIABLE as phenomenology (ladder B). The
    TXS-flare ghost test is rejected (all-sky flare trials swamp it).
- E-NF1 (NF-H04): ICECAT-1 v4 wide pairs (D-074 classes) in 7 lag bins (to 180 d), counts and Σ s_i s_j
  (signalness), 14 cells; `jit` ≤ 7 d, cyclic ±1 yr jitter beyond; 20,000 scrambles (the script refuses fewer
  than the per-cell threshold needs); empirical p; injection of ghosts with R_g × signalness; **95 % CLs** limits
  on R_g. ASSUMPTIONs: weighting, the long-bin null
  (stationary alert rate on ~1 yr), the R_g grid.

**Alternatives rejected.**
- F3 ε-line rescan of GBM × ICECAT (D-074 null; background ~2.4 pairs; no bound without a signal normalization).
- F1/F2 compute: no repo data can add to published work.
- A TXS 2014–15 ghost search in the 10-year track sample: IceCube's all-sky flare search finds chance hot spots at
  pre-trial p 9 × 10⁻⁶–3.5 × 10⁻⁷ (arXiv:2107.12134), so even R_g = 1 is not decidable.
- Classical (non-CLs) upper limits: a low fluctuation (1–7 d: 138 vs 149.5) "excluded" even R_g = 0 at 94 %.
- Testing injected catalogues against the observed catalogue plus ghosts: any ghost exceeds the observation, so
  limits came out too tight; the fix draws the background from the null and adds the ghost excess.
- Within-year permutation (D-074 `perm`) for the long bins: it keeps every pairwise lag, so it has no power.

**Evidence.** `results/nf/ghost_pairs.json`; reproduces D-074 (wide 36 at 1 h–1 d, 138 at 1–7 d). Pooled global
p = 0.058; min Bonferroni p = 0.094 (1 h–1 d, weighted, z = 2.7). The 36 pairs at 1 h–1 d were inspected:
spread over all 13 years; 6 are two alerts in one run and two are alert triplets (ordinary day-scale clustering
the ±3 d jitter keeps; without same-run pairs 30 vs 25.4, ~1σ, post hoc). 95 % CLs limits on R_g (per
astrophysical alert; the smallest grid R_g above which every grid point passes): ≤ 0.02 at Δt ≤ 10 s, ≤ 0.05
at 10 s–1 h, ≤ 0.2 at 1 h–1 d, ≤ 0.1 at 1–7 d, ≤ 0.5 at 7–30 d; at 30–180 d only R_g = 1 is excluded (R_g,50 not
reached). Round files: docs/neutrino_frontier/.

**Revisit if.**
- IceTracks-DR2 (doi:10.7910/DVN/MMIIZA) is streamed: R_g limits on the 1.6 M-event track sample, with its
  good-run lists as an uptime-aware null.
- A second public multi-event neutrino flare at a different distance appears (NF-H01).
- The owner sets a round-2 brief.

## D-078 Pipeline speed-ups (owner brief 2026-10-10, Part 0): numba kernel for GW-GW sky-map classes, durable derived-data store, 60-minute claim rule (2026-10-10)

**Decision.**
- S-2: `classify_gw_gw` (E1 GW sky maps) uses an exact ring-shift rotation (an RA rotation keeps a HEALPix pixel
  centre on its iso-latitude ring and moves it by floor(dr/w + 0.5) slots; the antipode sits on the mirror ring)
  and per-ring 128-bit masks. A numba kernel tests each pair ring by ring (window of B's doubled mask AND A's
  mask, early exit). numba is in the `dev` extra (so CI tests it) and a `fast` extra; without it a vectorized
  numpy bitmask path runs. The per-pair loop stays as the test reference. Unmapped GW events get no distance
  maps (never counted).
- S-1: `src/jwst_anomaly/derived_store.py`: reduced products (< 50 MB, never raw archives) are GitHub Release
  assets tagged `derived-data-YYYYMMDD`, pinned in `data/manifests/derived_data.ecsv` (sha256, size, sources,
  code commit). Scripts try the durable copy first and verify it, then the source (`e1_gw_skymaps.py` accepts it
  only when its content digest equals the pinned `maps_sha256`). `scripts/derived_publish.py` publishes.
- S-5: a claimed unit heartbeats at least every 10 min; takeover only after 60 min with no commit AND no
  heartbeat (docs/cloud-routine-prompt.md, research-cycle skill, docs/operations.md).
- S-3/S-4: profile before any run > 10 min; large fetches start in the background in the first minutes, with
  parallel range requests into memory, per-file caches and retried size probes (already in
  `e1_gw_skymaps.py`).

**Alternatives rejected.**
- Re-download each run: 704 MB from Zenodo at 0.15–2 MB/s per session (yesterday it did not finish in 40 min).
- git-lfs: bandwidth quota on a public repo and data in git history (CLAUDE.md).
- Owner-machine cache: not reachable from cloud sessions.
- numba as a core dependency (D-004 rejected PyOD partly for that); kept optional with a fallback.
- Precomputed per-pair shift tables (the pair set never saturates: 21k rows after 200 scrambles; 13×) and
  precomputed rotation tables (scattered gathers; 21×): slower than the kernel.
- A grid in RA rotation: not exact; the ring shift is exact (0 mismatches in 3.7 M rotated centres).

**Evidence.** CHANGELOG 2026-10-10 speed table. 30 fixed-seed scrambles, 19,313 real GW–GW pairs: loop, numpy
and numba classes identical; loop 808 ms, numpy 28 ms (×29), numba 1.79 ms (×452) per scramble.
- Release creation from a cloud session: HTTP 403 "Creating, editing, or deleting releases is not permitted for
  this session type". The store is built and tested; the first asset waits for a session that may publish.
- The live HEASARC GBM table no longer matches either pinned snapshot (D-074 5b42…, #122 0cb3…; live 8cc2…, 4,391
  rows): a pinned GBM snapshot cannot be re-fetched. It belongs in the durable store.

**Revisit if.** Cloud sessions gain release permission (publish the GW maps and the pinned GBM snapshot), or
another store is approved by the owner.

**Addendum (2026-10-10, local session on the owner's machine): first release; the pinned GBM snapshot is not
recoverable.**
- Published `derived-data-20261010` (pre-release, not latest): `gw_skymaps_nside32.npz`, 4,337,802 B, sha256
  ab033058…8e1ad5, code commit 67fffba (main; the tag itself sits on c92249c because `derived_publish.py` had no
  `--target` yet, fixed in the same PR), pinned in `data/manifests/derived_data.ecsv`. The npz is a fresh
  re-reduction (`e1_gw_skymaps.py --refresh`, 704 MB streamed into memory, 675 s); its content digest equals the
  pinned `maps_sha256` c5c2e4ff… and all four tarball sha256 equal their pins.
- Fresh-session check (`JWST_ANOMALY_DATA` = an empty directory): "durable copy verified in 1 s", 2.7 s wall, no
  tarball fetched.
- Neither pinned `fermigbrst.vot` (D-074 5b42…, 424,481 B; #122 0cb3…, 424,579 B) exists on the owner's machine:
  no file matched either name or either size in the data root, the old clone and its worktrees, the home folder
  (23.4 M files) or drives E: and G:. They were presumably fetched in ephemeral cloud sessions (inference, not
  checked), and the live HEASARC table has moved on, so neither can be re-fetched. **Decision:** E1 is re-pinned to a new snapshot in a separate PR that
  re-runs the affected E1 results and reports every change (TASKS "Part 0", cloud routine). No GBM snapshot was
  invented or downloaded here; `e1_events.fetch` is unchanged until a snapshot is published.

**Addendum 2 (2026-10-10, cloud routine): E1 re-pinned to one GBM snapshot.**
- `data/manifests/e1_events.ecsv` pins `fermigbrst.vot` sha256 8cc2823f… (424,672 B, 4,392 rows; the "live 8cc2…"
  above). ICECAT-1, GWTC and CHIME Cat 2 re-fetched with unchanged digests. `e1_gw_events.ecsv` (addendum 4 of
  D-074) is removed: the GW scripts read the same manifest.
- Relative to the D-074 pin: the new table minus bn261008763 and bn261007236 reproduces every tracked observed
  count (75 cells). All E1 scripts re-run with tracked seeds and scramble counts: D-074 core family p 0.278 → 0.273;
  addenda 2–5 global p 0.557 → 0.531, 0.676 → 0.568, 0.829 → 0.833, 0.555 → 0.551. No conclusion changes (CHANGELOG
  2026-10-10).
- The snapshot still exists only in an ephemeral session: a local session must publish it to the derived-data store
  (it can be re-fetched only while HEASARC serves this digest; otherwise re-pin there with `--refresh`).
- The null-ensemble cache is keyed by the input digests. Per-catalogue random streams are a follow-up.

**Addendum 3 (2026-10-10): MOA pre-screen chunks keyed on pre-screen inputs only.** `_chunk_done` and
`merge_prescreen` compared the whole `w3_moa.Params`, so vetting-only parameters (D-068 addendum: `period_*`,
`slow_dip_power`, `spike_*`) forced full re-streams with byte-identical rows (gb12 chunks 1, 4 and 7 diffed against
`main`: only metadata differs). Chunks now carry `prescreen_params` (`widths`, `n_min`, `prescreen_z/s/repeat`,
`COINC_Z`, `QUIET_TRACK_MOD`, `PRESCREEN_CODE`). Older chunks are judged on the same keys of their `params`
when `PRESCREEN_CODE` = 1. Bump `PRESCREEN_CODE` whenever `deficit_scan`, `baseline_chi2` or the tracked-row rule
changes. Rejected: keep the full-Params key (it costs 6–8 min and 24–30 GB of range reads per field re-run for no
change), and hash the scan source code (it breaks on comment edits and still needs a version for semantics).

## D-079 E-NF1b: IceTracks-DR2 short-lag ghost pairs: null; DR2 cannot reach R_g ~ 10⁻³ or improve lags > 1 h (2026-10-10)

**Decision.**
- Run NF-H04 on IceTracks-DR2 only at lags ≤ 1 h. Pre-registered (docs/neutrino_frontier/README.md, E-NF1b):
  northern tracks (Dec > −5°), log10 E ≥ 4.0 and ≥ 4.5, bins 0–10 s / 10–100 s / 100 s–1 h, D-074 wide pairs
  without pairs of one (run, event), uptime-aware ±3 d jitter null, 5,000 scrambles, 95 % CLs limits on R (ghosts
  per event above the cut). ASSUMPTIONs: cuts, bins, R grid, ghost passes the cut (ε = 1).
- Stream the 28 events/uptime files (sha256 of the served `.tab` pinned in `data/manifests/nf_icetracks_dr2.ecsv`),
  keep only a reduced `.npz` under the data root.
- Close the TASKS goal "reach R_g ~ 10⁻³ and 30–180 d with DR2": not reachable (Evidence).

**Alternatives rejected.**
- All lags on DR2: the optimistic floor max(3√B, 3)/N (every event astrophysical) is ≥ 0.055 at 1 h–1 d, ≥ 0.14 at
  1–7 d, ≥ 0.27 at 7–30 d and ≥ 0.64 at 30–180 d for every northern cut (log10 E ≥ 3.5–5.5): no gain over ICECAT-1.
- An ICECAT-alert-parent × DR2-ghost cross test: R_g per alert is bounded by ~3 / (astrophysical alerts ≈ 130)
  whatever the ghost sample, so it cannot beat the alert–alert limits by much.
- Southern tracks: dominated by atmospheric muons with a high energy cut (log10 E median 4.9).
- urllib downloads (403 from Dataverse); `requests` and curl work.

**Evidence.** `results/nf/ghost_pairs_dr2.json`, `results/nf/ghost_pairs_dr2_ic86.json` (identical on re-run, fixed
seeds). 1,643,355 events, no duplicate (run, event, subevent) (the paper's count); good-run union 4,963.4 d. TXS 0506+056 2014–15 box: 6 on
vs 1.11 expected (p = 1.0 × 10⁻³). Min Bonferroni p = 0.48; limits R ≤ 0.002 (log10 E ≥ 4, ≤ 100 s), ≤ 0.005
(≥ 4.5, ≤ 100 s), ≤ 0.01–0.02 (100 s–1 h). The 9 pairs at 0–10 s are all IC40/IC59; the rate above log10 E = 4 (north) is
×7–18 higher there than in IC79/IC86 (energy-proxy scale differs by configuration). IC86-only (post hoc): 0 pairs
at ≤ 100 s.

**Revisit if** NF-H04 gets a quantitative R_g or lag prediction; f_astro per cut is computed from the DR2 effective
areas (turns R into R_g); a public track sample with a much larger astrophysical count appears (IceCube-Gen2,
KM3NeT ARCA releases).

## D-080 Neutrino Frontier benchmark 1: SkyLLH reproduces the IceTracks-DR2 TXS 0506+056 2014–15 box fit (2026-10-10)

**Decision.**
- Adopt SkyLLH (PyPI `skyllh` ≥ 26.1, optional extra `nf`; IceCube's public-data likelihood framework, the tool
  arXiv:2605.19040 used) for point-source and time-window fits on IceTracks-DR2. No in-house likelihood.
- Benchmark before use: `scripts/nf_txs_benchmark.py` fits the TXS 0506+056 box (IC86_IV, catalogue position
  77.35°, 5.7°) with n_s and γ free, under the Dataverse 1.0 IRFs (what the paper used) and the ≥ 2.0 IRFs (binning
  fix). Two windows: Table 6 as printed ([T0 − ΔT/2, T0 + ΔT/2] = MJD 56927.5–57112.5) and the same window with
  each edge snapped to the nearest on-source event (< 1° from TXS, ≤ 1 d away; ASSUMPTIONs), because Table 6 rounds
  T0 and ΔT to whole days and a box flare fit puts its edges on events. Agreement tolerance |Δn_s| ≤ 1,
  |Δγ| ≤ 0.1 (ASSUMPTION).
- Raw files (IC86_IV events and uptime, IC86 IRFs: smearing matrix 817 MB v2 / 598 MB v1) are streamed into the data
  root, sha256-checked against `data/manifests/nf_skyllh_dr2.ecsv` and deleted after use (`--cleanup`). Stated reason
  for > 200 MB: the benchmark needs the full IC86 smearing matrix (SkyLLH builds the signal PDFs from it).

**Alternatives rejected.**
- The printed window as is: n_s = 11.1, γ = 2.22 (both IRF versions), 1.6 events short. The event at MJD 57112.653
  (0.38° from TXS, σ = 0.20°) lies 0.15 d after the printed stop; with the snapped edges (56927.860–57112.653:
  T0 = 57020.26, ΔT = 184.8 d, both round to Table 6) the fit matches.
- Multi-season fits (IC86_III–V): SkyLLH's time-dependent analysis supports one dataset only; the box lies inside
  IC86_IV (56757–57160), so nothing is lost.
- The whole-release zip that SkyLLH's dataset definition downloads: per-file streaming of the 4 needed files instead.

**Evidence.** `results/nf/txs_skyllh_benchmark.json`. Snapped window: n_s = 12.72, γ = 2.26, TS = 23.0 (v1 IRFs);
12.75 / 2.26 / 22.9 (v2). Paper: 12.7 / 2.3 (SkyLLH), 12.56 / 2.26 (internal). Printed window: 11.08 / 2.21
(v1), 11.11 / 2.22 (v2). The IRF binning fix moves n_s by 0.03. Fixed-window background (v2 IRFs), 50,000 SkyLLH scrambles per window: p < 2 × 10⁻⁵ (not comparable with the paper's
free-window pre-trial p).
The snap rule was chosen after the printed window missed (post hoc): the match is a consistency check that
the paper's box ends on these two events, not a blind test. n_s moves by ~1 event per edge event, so the tolerance
alone cannot separate a correct setup from a lucky window.

**Revisit if** a new SkyLLH release changes the PDFs (re-run the benchmark), or a test needs multi-season
time-dependent fits (SkyLLH single-dataset limit).

## D-082 W3 short-event (t_E = 3 d) recovery study: frozen sample, loss diagnosis on dev only (2026-10-11)

**Decision.** Owner directive V2.1, priority 1. The production chain (CHAIN_VERSION 2026-10-10.1, unchanged) runs
2,000 W3 injections at t_E = 3 d (ρ = 0.01 and 0.1, LF magnitudes, seed 3003; no limit uses this seed) on the
300 gb12 quiet carriers (`scripts/w3_short_te.py`, `results/w3_moa/short_te/`).
- **Split.** Each injection is assigned dev or validation by the SHA-256 of its own inputs (event_id, t_E, ρ,
  u0, t0, I_s) and the seed. This is independent of row order and of every outcome.
- **Analysis.** The loss is diagnosed on dev only. For each flagged injection, count spike nights (injected flux
  above `SPIKE_SIGMA` = 3 median errors) and umbra nights (injected flux below −3 median errors).
- **Validation.** It has been read only for its baseline recovery counts. The stage-2 classifier uses it once,
  as its final evaluation.

**Alternatives rejected.**
- A split keyed on the row index (the first version). `run_inject` collects rows with `imap_unordered`, so rows
  come in completion order, and slow (flagged) injections cluster late. `/code-review` found that this split
  was neither reproducible nor independent of the outcome: flagged injections divided 72/46, p ≈ 0.02.
- crc32 % 2: it is an affine function of the input bits, not a random draw.
- Re-using the production gb12 injections: only 23 flagged at t_E = 3 d, and they feed the published limits.
- Relaxing `jackknife_nights` directly: the test exists because one bright night carried gb7-R-8-6-94052's
  preference (w3-survey failed-approach rule).

**Evidence.** Baseline recovery of flagged injections:
- dev: 21 / 57 = 0.37 (Clopper–Pearson 95 % 0.24–0.51);
- validation: 27 / 61 = 0.44 (0.32–0.58).

First failing test on dev: `jackknife_nights` 14, `eclipse_dip` 12, others 10.

Spike nights on dev:
- **`eclipse_dip` losses (12):** 11 have ≤ 1 spike night, and 1 has spikes on both sides of the umbra. With the
  peaks unsampled, these signals cannot be told from a flat dip in the data. This is a sampling limit, not a rule
  defect.
- **`jackknife_nights` losses (14):** 1–4 spike nights; 9 have spikes on both sides.
- **Recovered (21):** 0–9 spike nights; 17 have spikes on both sides.

The jackknife losses and the recoveries overlap in spike nights, so this proxy cannot separate a sampling
limit from a rule loss for `jackknife_nights`. No bound on the rule-attributable loss is claimed. A first draft
gave "≤ 8 %"; it is withdrawn: it rested on the row-order split and on a proxy that differs from the
jackknife's own feature definition, |exotic − best ordinary model| > 3 median errors.

**Revisit if.** Stage 2 re-vets the flagged dev injections and records the jackknife's own feature nights,
dropped nights and ΔBIC. It then tries a night-robust likelihood term (a per-night systematic error instead of
dropping nights), tuned on dev, checked for false positives against the real gb12 flags and the PSPL
controls, and evaluated once on validation.

## D-083 Four-week high-risk exploration trial: split WIP counters, pending index, unit state files (owner-approved 2026-10-12)

**Decision.** The owner approved a trial from 2026-10-12 to 2026-11-09 (PR #174; owner instruction in session):
- a weekly split of run time: 45 % routine search, 30 % high-risk units (at most two in depth), 15 % verification
  and 10 % upkeep;
- two WIP counters: awaiting merge (cap 3) and awaiting owner review (backlog rule: with ≥ 10 open, at most one new
  needs-human PR per ISO week);
- `docs/state/pending.md`, which marks unmerged conclusions as provisional;
- per-unit state files with pre-written stage gates;
- a four-column outcome ledger;
- an idle rule.

Rules T1–T6 are in `docs/state/trial.md`, version trial-2026-10-12.1. First unit: X1 (far-side SV hot spots,
`docs/state/units/X1.md`). The routine prompt copy moves to version 2026-10-12.1.

**Alternatives rejected:**
- Rebuilding the agent system or adding a routine. Rejected by the owner.
- Counting all needs-human PRs toward the cap of 3. In practice this let only self-mergeable field sweeps advance:
  14 needs-human PRs were open, while draft claim PRs and issue-recorded research stayed allowed.
- A hard "< 10 % limit gain stops a field" rule. Replaced by batch judgement (T4).

**Evidence.** Audit `docs/audits/2026-10-12-exploration-trial.md`:
- The routine's activity was observed: claim PRs opened at minute :08–:11.
- Its configuration was not readable from a session.
- TASKS.md listed S1 work as "Next" that #167 and #169 had already done.

**Revisit if.** The trial ends (2026-11-09) and the review of finished unknowns, method capability, candidates and
verifications recommends continue, adjust or stop. Revisit earlier if the owner-review backlog keeps growing or a
unit overruns its budget.

## D-084 X3: W3 screen on unselected Gaia epoch photometry (GAPS dev → held-out → DR4); D-057 fitter plus an astrometric-quality veto (2026-10-14)

**Decision.**
- Test W3 (negative-mass inverted microlensing, D-047) on Gaia epoch photometry that was never selected by
  variability: GAPS now, DR4 later.
- Screen v1:
  - the D-057 fitter flag, ΔBIC(N1neg − best of PSPL/PAR) < −10, with D-061's G-error rescaling;
  - plus a veto of sources whose photometry depends on scan direction: no 5-parameter astrometric solution,
    RUWE > 1.4, or `ipd_frac_multi_peak` > 10.
- The DR3 variability flag is deliberately not used, because a real event would set it.
- Gaia DataLink calls use `DATA_STRUCTURE=DATAMODEL_STANDARD`.

**Alternatives rejected:**
- the PTA Shapiro sign test and FRB exotic lensing (S0-2026-10-14 A, B);
- `vari_microlensing` (D-061: the selection removes W3);
- vetoing `phot_variable_flag = VARIABLE` (it would veto real events);
- tuning ΔBIC instead of vetoing on independent astrometric quality: the false flags were blends, not marginal
  fits;
- DataLink for the 0.96 M held-out sources (about 13 days): use the bulk epoch-photometry files.

**Evidence** (`derived` from `simulated` injections into `observed` light curves; dev only):
- seed 2026: v0 efficiency 0.76 with real flags 0.027;
- seed 4242 (fresh, v1): efficiency 0.70 (t_E = 30–300 d), real flags 0.001, PSPL→W3 0, vetoed 7.5 %;
- tables in `results/x3/`.

**Revisit if:**
- DR4 epoch photometry does not cover unselected sources;
- G1's periodicity or achromaticity vetoes cost more than half the efficiency;
- the held-out false-flag rate exceeds the dev rate by more than 3×;
- Gaia changes DataLink again.

## D-085 X3 G1: frozen full W3 screen, Gaia bulk-file route, held-out GAPS run once (2026-10-14)

**Decision.**
- Screen, frozen before calibration (`scripts/x3_g1.py`): `MIN_POINTS` → v1 → pre-screen → v0 fit flag →
  repeat-dip veto (no period search: a sparse-scan periodogram peak can be an alias) → common-magnification colour
  test with per-band source and blend fluxes ("insufficient" never vetoes).
- Held-out light curves come from the DR3 bulk epoch-photometry files, streamed (`scripts/x3_bulk.py`). The split
  stays `MOD(source_id, 10) = 0`, which is 20 % of the sources.
- Flags are vetted in batch by ordinary rules (`scripts/x3_vet.py`). The rules are post hoc, so every rule is also
  applied to the dev W3 injections (event rules to the light curve, host rules to the host), and the efficiency
  includes their losses. v1-vetoed hosts are left out of the efficiency, because the exposure excludes them.

**Alternatives rejected:**
- re-splitting to a true 10 % (leaks already-read sources);
- DataLink for 423,598 sources;
- a period-search veto;
- requiring constant BP − RP (wrong for blended events);
- tuning the colour test after seeing that it barely separates W3 from dust (the screen was frozen; the fix goes
  to DR4);
- vetting the 196 flags one by one with the JWST-oriented skill.

**Evidence:**
- Bulk = DataLink on 309 / 309 dev sources.
- Dev calibration: W3 efficiency 0.32 (t_E 30–300 d), PSPL 0 / 1,995, real 0 / 1,995.
- RUWE inflation would veto an estimated 92–99.6 % of blended W3 events (model estimate under assumed blend and error
  priors; not a measured or strict bound).
- Held-out: 395,764 exposed sources, 196 flags excluded by rules (not each proven ordinary), no candidate;
  vetted efficiency 0.276.
- Γ₉₅ < 9.7 × 10⁻⁶ per star per year, photocentre-stable (`results/x3/g1_*`).

**Revisit if:**
- DR4 epoch astrometry is documented. The blend-ratio photocentre shift is NOT a signature: it is exactly
  degenerate with an ordinary blended variable. Only the lens image shift is new, measurable for G ≲ 17,
  θ_E ≳ 3 mas and t_E ≈ 100 d, under the assumed per-transit σ (X3.md, G2 prerequisite; addendum 2026-10-14 by owner correction);
- a joint G/BP/RP fit is built (it replaces the colour test);
- Gaia changes the bulk-file layout.

## D-086 X4: wave-optics GW lensing by repulsive lenses (S0-2026-10-15); X3 paused (2026-10-15)

**Decision.**
- Open unit X4 (`docs/state/units/X4.md`). Its S1 tests whether a free choice of a binary's intrinsic parameters
  absorbs a negative-mass lens's low-pass signature, before any strain analysis.
- X3 pauses until the DR4 epoch-astrometry content is confirmed.

**Alternatives rejected:**
- VASCO vanishing sources (no link to the hypothesis);
- the quantum-processor wormhole reanalysis (a simulation, not evidence);
- going straight to strain-level PE on GW231123 before the degeneracy question is answered.

**Evidence:**
- S0 A: closed-form F_s(w, y), checked two ways.
- Separable from a positive point mass at observed SNR ρ ≈ 3.1–16.5 for M_Lz ≥ 300 M☉, y ≤ 1, with intrinsic
  parameters fixed; those cases are demagnified (≳ 100× less detectable volume)
  (`results/s0/gw_neglens.json`).
- LVK searches use only positive lenses.

**Revisit if:**
- S1 shows the signature is absorbed by intrinsic parameters (stop X4);
- LVK or others publish a repulsive-lens GW analysis (reuse it).

## D-087 X4 S1: the negative-lens GW signature survives free ordinary fits in simulation; S2a pilot proposed (2026-10-15)

**Decision.**
- Record S1's gate as met.
- Propose the S2a pilot: a maximum-likelihood three-model comparison on GW231123 strain, with ≥ 50 off-source
  injections. Do not start it without the owner, and freeze its gate before the strain is read.

**Alternatives rejected:**
- going straight to full Bayesian PE on GW231123 (cost, before the cheaper pilot calibrates the statistic);
- using the S0 inspiral-only, fixed-intrinsic result alone (shown optimistic in principle; S1 now quantifies it).

**Evidence** (simulated only; `results/x4/s1_*.json`):
- 18 / 36 configurations keep 1 − FF² ≥ 9/ρ² at ρ ≤ 30 against free aligned-spin binaries (q ≤ 15) and a positive
  point lens, using the best fit over all stages (the frozen 6-parameter refinement was under-converged at 80
  evaluations). IMRPhenomD agrees within 0.0065. Bound-conditional: fits run to q = 15.
- Extreme mass ratios absorb part of the signature.

**Revisit if:**
- precessing, higher-mode or embedded/binary-lens fits absorb the residual;
- real noise or the detector network changes the picture (S2a).

## D-088 X4 S1b: with precessing, higher-mode fits the representative passing cases fall below the screening threshold; no real-data pilot (2026-10-15)

**Decision.**
- Do not propose the S2a pilot. D-087's "Revisit if" condition (precessing or higher-mode fits absorb the
  residual) is met for the two passing representative cases that could be checked.
- Keep X4 open only for R1-type (strongly demagnified, low-mass) cases. They first need a positive-lens optimiser
  that passes its own injection test, then XPHM. That step needs the owner's go-ahead.

**Alternatives rejected:**
- proposing S2a on GW231123 anyway (neither R2 nor R3 passes the screening threshold after the fuller ordinary fit; the
  other cases were not re-tested);
- declaring the GW route degenerate (R1-type cases are undetermined, not absorbed);
- lowering the screening threshold after seeing R2's ρ_req 30.8 against 30 (not allowed).

**Evidence** (simulated only; `results/x4/s1b_*.json`):
- Numerical stage: R2 and R3 keep a residual against aligned-spin models (fit errors ≤ 0.1 × mismatch). Strictly,
  they are undetermined, because the best P block descent was not converged. R1's P fit error is 0.13.
- IMRPhenomXPHM: best ordinary FF 0.99523 (R2) and 0.99647 (R3), so ρ_req 30.8 and 35.8, above 30. R2 is fitted
  through an edge-on inclination, R3 through in-plane spins, both at q 1.2–1.7.
- The convergence continuations and extra R2 starts were added post hoc and are labelled as such.
- R3 was already absorbed in the frozen run (ρ_req 31.4). R2 passed the frozen run (28.9) and is absorbed only
  after the post-hoc starts.

**Revisit if:**
- a validated P optimiser plus XPHM leaves an R1-type residual at ρ_req ≤ 30;
- priors on inclination and precession are added (a likelihood screen ignores them).

## D-089 X4 S1c: R1 is kept only against aligned-spin templates; undetermined with XPHM; still no pilot (2026-10-15)

**Decision.**
- Do not propose a real-data pilot. R1 meets outcome 2 (residual kept, validated) only for aligned-spin,
  dominant-mode templates. With IMRPhenomXPHM it is undetermined.
- Keep the shared fitter (`scripts/x4_s1c_r1.py`) for any later X4 fit. Validation and case fits must use the
  same code.

**Alternatives rejected:**
- treating the aligned-spin result as enough for a pilot (S1b showed precession and higher modes change verdicts);
- treating the best XPHM search value (ρ_req 9.2) as a bound (validation failed, no convergence);
- extending past the 2 h cap.

**Evidence** (simulated only; `results/x4/s1c_*.json`):
- Aligned-spin validation passes, including S1b's failed injection (ε 3.5e-6 ≤ 0.0145). P reaches FF 0.8545 from
  all 8 generic starts (ρ_req 5.8), at q 15 and |χ| 0.99. The result is provisional: the frozen saved-solution
  starts were cut for budget, although those solutions themselves score 0.8545.
- XPHM: unlensed validation ε 0.026 > 0.0125. The lensed fits do not agree across starts (0.9454 against 0.9441).

**Revisit if:** a precessing fitter passes the unlensed and lensed injection checks with agreeing starts, and R1
then stays at ρ_req ≤ 30.

## D-090 X4 S1d: SciPy differential evolution as the max-match reference; it fails validation, R1 stays undetermined (2026-10-15)

**Decision.**
- Use SciPy `differential_evolution` as the reference global search for XPHM max-match fits, with polarisation
  and amplitude maximised analytically.
- Do not fit R1 or propose a pilot while the tool fails its frozen injection checks.
- Next step, only on the owner's go: a joint (lens + intrinsic) DE with the lens snapped to the cached grid, and
  longer unlensed runs.

**Alternatives rejected:**
- Bilby + dynesty: not used this round, for task fit and cost. It does Bayesian inference, not max-match, at an
  estimated (not measured) 10⁵–10⁶ likelihood calls. It can use XPHM with phase marginalisation turned off.
  Kept for a later real-data stage. (Wording corrected in S1e.)
- PyCBC: no max-match optimiser for this objective (its complex overlap would give the real part, but the search
  would still be ours); GPLv3. Only the polarisation-maximisation idea is reused, in closed form.
- GLoW: it computes lens factors, which are already exact here.
- CMA-ES (`cma`): fallback if the joint DE fails.
- Treating S1c's best XPHM value (ρ_req 9.2) as a bound.

**Evidence** (simulated only; `results/x4/s1d_*.json`):
- Unlensed ε 5e-6 – 1.5e-3, against 0.026 in S1c; every run stopped at the 200-generation cap.
- Lensed ε 0.045–0.091 in 3 of 4 runs, from a weak-lens stall after the node-repeat stop; the polish also hit its
  round cap while still improving.
- The objective at the truth is ≤ 3.3e-16.

**Revisit if:** the joint DE passes all four injections on both seeds (then fit R1), or fails them (then try
CMA-ES once before any further tool work).

## D-091 X4 S1e: joint lens + intrinsic DE as the P search; one lensed check fails, R1 not fitted (2026-10-15)

**Decision.**
- Keep the joint 12-parameter DE (lens snapped to the cached grid) with exact continuous refinement as the P
  search, and the 600-generation DE as the U search.
- Do not fit R1 or propose a pilot while one validation run fails the frozen rule.
- Any continuation is the owner's call: more P budget, expert advice on a precession parametrisation, or an R1 U
  fit testing outcome 1 only.

**Alternatives rejected:**
- Running R1 anyway: the frozen plan put R1 after a full pass.
- Dropping or re-seeding the failing run, or loosening the 5e-4 rule.
- Adding CMA-ES this round (the owner excluded it).
- Treating the grid as exact: its floor is ε 9e-4 – 1.3e-3.

**Evidence** (simulated only; `results/x4/s1e_validate.json`):
- U: ε ≤ 2.1e-6 on all four runs.
- P: three runs at ε 2.4e-5 – 1.5e-4; V9 seed 11 at 1.5e-3, with its DE and its refinement both still improving
  at their caps.

**Revisit if:** the P search passes all four runs with a larger budget (then fit R1), or the owner authorises the
outcome-1-only U fit.

## D-092 X4 S1f: gate revision for constructive ordinary counter-examples; none found for R1 (2026-10-15)

**Decision (owner, 2026-10-15).**
- An exactly verified ordinary solution with ρ_req > 30 makes a case fail the screen without proof that the whole
  search tool is reliable.
- The absence of such a solution supports no claim; the case stays undetermined.
- For R1, the validated unlensed search found none (best FF 0.9135, ρ_req 7.4). R1 stays undetermined, with no
  pilot.

**Alternatives rejected:** treating the failure to find a counter-example as evidence that R1 is identifiable.

**Evidence** (simulated; `results/x4/s1f_*.json`): seeds 21–23 give 0.8923, 0.8923 and 0.9135. The independent
rebuild agrees to 1.2e-7. The best ordinary fit of any kind is 0.945 (S1c, lensed).

**Revisit if:** a validated positive-lens search exists, or the owner authorises a new R1 question.

## D-093 S0-2026-10-15b: Kapitza dynamic stabilisation of wormholes paused as a low-priority open item (2026-10-15)

**Decision.**
- Do not open a follow-up unit. Record that the only ordinary-matter "pivot shaker" found, an external high-frequency
  EM/GW radiation bath, acts through its averaged (Isaacson / Brill–Hartle) stress.
- In an exact thin-shell toy it can stabilise the radial mode, but only with radiation energy of the order of the
  exotic energy. For throats a₀ ≥ 4M
  that raises the exotic negative energy by 35–130 %, and it needs about 10⁵² W or a hypothetical
  confinement.
- The decisive next calculation, about 3–4 h: Ellis–Bronnikov plus a radiation fluid, the l = 0 eigenvalue. It
  waits for the owner.

**Alternatives rejected:**
- J → J(t), or an arbitrary sinusoidal term in the mode equation: no physical source.
- The internal phantom-field oscillation: it is still exotic, with no external control.
- Transferring the coupled-SYK Floquet results (2404.08394) to 4D.

**Evidence:**
- [lit] the Kapitza criterion; the Ellis growth rate 1.18 c/r; topological censorship (FSW).
- [der] the shell-theorem no-go for spherical external drives; the reduction to averaged stress; the first-order
  drive needing about 6–7 times the potential depth.
- [num] `results/s0/kapitza_shell.json`.

**Revisit if:** a field-supported calculation shows stabilisation at a small radiation fraction, or a mechanism
confining radiation at a throat is found.

**Addendum (2026-10-15, re-check):**
- The decisive next calculation named above is **cancelled**: the Ellis + static isotropic radiation-fluid
  background does not exist. ρ = C/N⁴ forces ρ ≡ 0 for asymptotic flatness, and there is no surface layer to
  truncate it.
- Two supporting claims are withdrawn: the general spherical no-go (it needs a vacuum gap), and "first-order drive
  is non-perturbative".
- A stability calculation will be proposed only with a fully specified radiation closure and a fixed matter
  theory.

## D-094 `arc_diff` per-stamp kernels: shared design by default (2026-10-17)
- **Decision:** `kernel_match --kernel-stamps` computes its leave-one-out stamp kernels with `kernel_fit_stamps`:
  - one weighted design and Gram matrix per band;
  - a per-fit downdate of the own stamp and the clipped rows, with a direct sum when the downdate dominates;
  - a Cholesky solve refined with the exact residual until converged, falling back to a direct sum.
  - It is the default (`STAMP_FIT_DEFAULT = "shared"`); `refit` stays available.
- **Alternatives rejected:**
  - dropping leave-one-out, fewer iterations or lower precision (they change the search; owner);
  - sharing one solution across fits;
  - plain normal equations without refinement (lose digits on bright stamps: review 2026-10-17);
  - a different algorithm or GPU (out of scope).
- **Evidence (`results/transients/kernel_loo/RESULT.md`):**
  - every B1 and B2 stamp kernel agrees to ≤ 1.2e-13 relative, with the same clipping (χ²);
  - on B2, S agrees to ≤ 3e-8, with identical hits, flags, morphology, review list and 60 injection outcomes;
  - the whole B2 run is ×14 (screen) and ×8–9 (injection run) faster, with the same peak memory.
- **Revisit if:**
  - a design is ill-conditioned (condition number ≫ 1e6; B1 and B2 are 320–1,800);
  - stamp fits exceed 200,000 rows (these already fall back to `refit`);
  - the clipping rule changes.

## D-095 Difference-imaging tool for batch searches: keep `arc_diff` (#291); SFFT not adopted (2026-10-17)
- **Decision:** the batch search keeps `arc_diff` with the D-094 stamp kernels.
- **Alternatives rejected:**
  - SFFT's JWST B-spline method (Hu & Wang 2024): GPU-only in sfft 1.7.3; a port is out of scope.
  - Ordinary SFFT without masking (corrected adapter): somewhat worse on B1 + B2 (fades 83/120 vs 92/120; unchanged
    false hits 5/60 vs 1/60), though about ×1.4 faster (indicative) with 1.49–2.51 × the peak memory.
  - SFFT's automatic masking: needs the SExtractor binary.
- **Evidence:** `results/transients/sfft_trial/RESULT.md`. Small development sample; not "generally worse".
- **Revisit if:**
  - a GPU is available, then test the B-spline method on the same B1/B2 injections;
  - SExtractor is installed for a masking-mode test;
  - a box larger than 1′ shows spatially varying PSF ridges (DECISIONS addendum 2026-10-15c).

## D-096 W1 PSF anisotropy: STPSF used as a diagnostic model; no PSF correction adopted (2026-10-17)
- **Decision:** STPSF 2.2.0 (data 2.2.0, streamed) models NIRCam F200W PSFs for the W1 floor tests; no constant
  PSF-χ correction is adopted for the shear screen. The model's isophotal |χ| depends on the isophote level and
  over-predicts the three Abell 2744 stars by 0.05–0.07; its axis is close to theirs (2 of 3 within 1.2°).
- **Alternatives rejected:**
  - a constant model χ per isophote level: only the two levels with |χ| 0.18–0.21 pass the 3.79 rule. That
    level fits the stars, not the fainter sources, and adds E-mode power;
  - a star-based map (D-entry of #298 / `results/w1_psf`: 3 stars).
- **Evidence:** `results/w1_stpsf/RESULT.md`. D-018's rejection of STPSF for spike templates is unchanged.
- **Revisit if:** a per-source forward model (galaxy ⊗ drizzled STPSF PSF, catalogue isophote) reproduces the
  three stars' χ, or the screen moves to PSF-weighted moments.

## D-097 W1 streak veto: STScI `acstools.findsat_mrt` trail masks, run through the injections (2026-10-17)
- **Decision:** the W1 shear screen drops shape rows (and painted injected images) within 1″ of persistent trails
  that `acstools.findsat_mrt` finds on the streamed F200W mosaic (`scripts/w1_streak_mask.py`,
  `inject_shear.py --streak-masks`). Adopted by the frozen rule of `results/w1_streak_veto/PLAN.md` (K1, K2 hold).
- **Alternatives rejected:**
  - the post hoc column test of `w1_peak_vet.py`: one peak only, columns only, not run through the injections;
  - catalogue collinearity: chains are common by chance (`results/w1_floor/RESULT.md`, 21 % of used sources);
  - a hand-written Hough transform: findsat_mrt is maintained STScI code (BSD-3) built for this (charter order 1);
  - zeroing > 5σ pixels before the finder: removes bright streaks (PLAN addendum A).
- **Evidence:** `results/w1_streak_veto/RESULT.md`. MACS1423 S 3.88 → 3.41, p_rot 0.01 → 0.22; no field passes.
  Vetoed rows are blanked inside `shear_screen` after the spike veto; dropping them earlier changed the spike veto.
- **Revisit if:** a faint-streak case (S1063's parallel streaks) needs a lower threshold; recalibrate through
  injections before using it.
- **Addendum (2026-10-18):** tried for the two flagged CANUCS NCF (`results/w1_streak_rowdiag/RESULT.md`):
  thresholds 3 and 2.5 give masks identical to 5 (the persistence check keeps the same trails); the features are a
  diffraction spike and detector row striping, not trails. A lower threshold is not the route; D-097 is unchanged.
- **Addendum 2 (2026-10-18):** a field-wide row/column-segment median veto (`w1_streak_mask.py --seg-len`,
  `results/w1_row_veto/RESULT.md`) covers the macs0417-ncf peak but keeps only 56–75 % of the NCF shape rows
  (cost rule ≥ 75 % in every field: 1 of 5 at L = 80, 0 of 5 at L = 40); not adopted, D-097 unchanged. A next stripe
  veto needs source footprints masked before the segment medians.

## D-098 W1 shear screen in blank fields: no lens model, R fixed from the same program's cluster pointings (2026-10-18)
- **Decision:** fields without a lens model run the unchanged D-053 shear screen with `inject_radial.NullLensModel`
  (κ = g = 0) and a fixed responsivity R. R cannot be fitted without shear; it is fitted on cluster pointings of the
  same program and band and the smaller value adopted (gate R, `results/w1_shear_ncf/PLAN.md`; CANUCS NCF in F150W:
  R = 0.425). The W1 volume headline adds only blank fields whose screen does not pass the E/B rule.
- **Alternatives rejected:**
  - R = 1: overstates the injected signal (fitted R 0.41–0.60 in ten of the eleven headline fields, 0.23 in MACS0647;
    `results/w1_volume/summaries/z0.3_*.json`);
  - the cluster fields' F200W R: the NCF have no F200W, and the band changes R by up to 0.07 (gate R);
  - DJA GOODS-N shapes: SEP moments, a different shape measurement whose R has no calibration yet.
- **Evidence:** `results/w1_shear_ncf/RESULT.md`. Fourteen fields: n₉₅ < 2.35 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉ (2.35×
  Takahashi & Asada, from 3.6×). Two NCF flags sit on linear artefacts: macs0416-ncf explained by the frozen V1/V2
  rule, macs0417-ncf only by a post hoc row-stripe test whose rotated-band control is weak (9 % of orientations).
- **Revisit if:** a blank field has a shape catalogue of another kind (calibrate its R first), or a row/diagonal
  streak veto is frozen and run through injections (then macs0416-ncf and macs0417-ncf can enter the headline).

## D-099 W1 in GOODS-N: R borrowed across programs, photo-z coverage gate; headline per D-098 (2026-10-18)
- **Decision:** JADES 1181 F200W pointings run the D-098 blank-field screen with R = 0.414, the smallest F200W
  cluster-fit R apart from MACS0647 (Abell 370; program 1181 has no cluster pointing to fit R on). A pointing enters
  only if ≥ 50 % of its used rows have a photo-z match (gate P, `inject_shear.passes_gate_p`). The headline follows
  D-098 (flagged pointings out, even when the frozen V1 ∪ V2 rule explains them).
- **Alternatives rejected:** R = 0.231 (MACS0647; unexplained and nearly blind, kept as the sensitivity case); the
  median cluster R (≈ 0.48; less conservative); the eastern pointings o009/o010 (photo-z share 0.22 / 0.02).
- **Evidence:** `results/w1_shear_goodsn/RESULT.md`. Sixteen fields: n₉₅ < 1.82 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉
  (1.82× Takahashi & Asada, from 2.35×); the two flags (o003, o004) sit on mosaic-edge artefacts and are explained
  by V1 ∪ V2 (eighteen fields would give 1.54×). With R = 0.231 the GOODS-N recoveries drop to 0.21–0.26 of the
  adopted ones (four-pointing totals at Z = 0.3).
- **Revisit if:** a blank-field R calibration exists (e.g. image simulations or a shear-calibrated catalogue of the
  same mosaics), or a photo-z catalogue covers o009/o010.
- **Addendum (2026-10-18, CEERS):** the same rule applied to CEERS (program 1345, F200W, 2.8 ks; DJA v7.4 zout;
  pointings with pairwise overlap ≤ 0.6 arcmin², o052 out). Ten pass gate P; o004 / o069 flags explained by V1 ∪ V2
  and kept out. Twenty-four fields: n₉₅ < 1.09 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉ (1.09× T&A); R 0.231 keeps 0.19–0.28
  of the CEERS recoveries (`results/w1_shear_ceers/RESULT.md`). The borrowed R now dominates the systematics.
- **Addendum 2 (2026-10-18, PRIMER):** mosaics larger than the 200″ screen radius are cut into non-overlapping 280″
  tiles (`w1_footprints.tile_grid`, each tile one mosaic's rows, share ≥ 0.3); same rule otherwise. Twenty tiles pass
  gate P; four flags kept out (x3y0 unresolved). Forty fields: n₉₅ < 6.38 × 10⁻⁵ h³ Mpc⁻³ at 2 × 10¹³ M☉ (0.638× T&A)
  under R = 0.414; R 0.231 keeps 0.15–0.28 of the PRIMER recoveries (`results/w1_shear_primer/RESULT.md`). Per the
  frozen plan, W1 turns to lower masses rather than more area.

## D-100 Same-source pairing: no model-free proposal step adopted; lens tools are pair testers, not pair finders (2026-10-18)
- **Decision:** the repository has no reliable model-free step that proposes which sources are images of one source
  (task 1). An identity-blind SED pairing at 2–60″ (`scripts/same_source_pairs.py`, reusing `orphan_pairs.match_table`)
  failed its frozen pass rule on Abell 2744; it is kept as a feasibility tool, not a search. lenstronomy, PyAutoLens and
  Lenstruction are not added as dependencies: each tests or reconstructs a pair or system already chosen (task 2),
  which published cluster models already cover here (`lensmodel.find_images`).
- **Alternatives rejected:**
  - affine morphology matching of either handedness (`parity_test.compare`) as the proposer: AUC 0.53 for true pairs
    vs matched decoys (6 + 6, development half), and the parity tool is paused (D-032);
  - integrating all three lens tools at once: none solves task 1; Lenstruction is unmaintained since 2022-11 and pinned
    to lenstronomy 1.10.3 (imports with 1.14.2 except one module needing photutils).
- **Evidence:** `results/same_source/a2744_feasibility/RESULT.md`: holdout recall@10 0.86 (relative ranking works) but
  only 10 of 33 true pairs pass the absolute SED rule (20 fail χ², 13 photo-z overlap); post hoc, 26,145 SED-matched
  pairs among 2,674,270 at 2–60″, of which 11 have both members within 0.3″ of a published image (10 same-system pairs pass; the rest are unverified, not shown to be false). Limits on later searches: `docs/same_source.md`.
- **Revisit if:** a ranking-based proposal with a calibrated threshold is frozen and tested on a new held-out cluster
  (e.g. MACS0416 CANUCS), or a specific pair needs a source-plane reconstruction (then lenstronomy in an isolated env).

## D-101 W1 blank-field R: GalSim image simulations of the pipeline-like shapes; borrowed R kept (2026-10-18)
- **Decision:** the blank-field W1 screens keep R = 0.414 (D-099). Image simulations of the catalogue shape
  measurement (`scripts/w1_r_sim.py`, GalSim + photutils, isolated Sérsic stamps) give a weighted measurement
  response ≥ 0.94 for the GOODS-N row distribution under every variant, so 0.414 understates injected W1 shear (conservative).
- **Alternatives rejected:** adopting R_sim now (the simulation has no blending, correlated noise or real PSF, and
  the frozen rule only allows "keep" or "re-run lower"); metacalibration on real i2d cutouts first (needs an STPSF
  deconvolution and S3 cutouts; the next stage); ngmix metacal (fits its own shapes, not the pipeline's isophotal
  moments).
- **Evidence:** `results/w1_r_metacal/RESULT.md`, `r_sim.json`, `r_sim_rawthr.json` (simulated).
- **Revisit if:** a real-image calibration gives R < 0.414, or the cluster-fit gap (≈ ×2) is traced to the
  measurement rather than the lens models.

## D-102 Achromaticity as the discriminating test of a two-epoch fade; G165 h002 cannot be decided (2026-10-18)
- **Decision:** for a fade seen in the same bands at two bright epochs, the test between an intrinsic transient and
  magnification by a compact deflector (point mass, Ellis wormhole n = 2, negative-mass n = 1; D-047 family) is the
  band-independence of the calibrated epoch ratio (`scripts/achromatic_test.py`: per-band matched-control zero point
  and scatter, χ² of a constant, p < 0.01). Its power is simulated first with SN templates (sncosmo, isolated env).
  It does not separate the deflector types from each other. It is a band-consistency test under stated assumptions
  (two epochs, three bands, a uniform-colour source, geometric optics), not a general supernova-versus-lensing
  classifier (wording narrowed 2026-10-18).
- **Alternatives rejected:** comparing one SN light curve with one wormhole light curve (two ideal curves, no
  ordinary lens, no errors); a light-curve-shape test (W3 spike–umbra–spike): needs ≥ 3 well-sampled epochs, none of
  the open objects has them.
- **Evidence:** `results/transients/confirm_r2/RESULT.md`: B1-PLCKG165 h002 not decided (too few usable controls; h002
  is near the edge of the 4446 coverage); power 0.15–0.62 against SN templates at its real epochs and errors (two
  noise models); post hoc p 0.03 / 0.30.
- **Revisit if:** a new public epoch or a spectrum of a fading point source exists while it is bright, with ≥ 3 bands
  calibrated to a few per cent; or a local-error photometry of the 4446 edge region recovers usable controls.

## D-103 W1 lower masses: the frozen Schirmer 10″ shear filter stays; re-weighted filters do not open ≤ 2 × 10¹² M☉ (2026-10-18)
- **Decision:** the W1 shear screen keeps its frozen filter (Schirmer, 10″, r_min 1″; D-053) for lower masses.
- **Alternatives rejected:** the point-mass (1/θ²) filter at 10″ (N(1e12 + 2e12) 0.33× / 0.43× of V0 on the two
  pilot pointings); a 5″ Schirmer aperture (0.93× / 0.57×); the frozen adoption threshold was ≥ 1.5× in both.
- **Evidence:** `results/w1_low_mass/RESULT.md` (simulated injections, ceers-o001 and goodsn-o004, Z = 0.3).
  V0 recovers 3–4 / 100 at 1 × 10¹² M☉, so a forty-field 1 × 10¹² grid would sit at roughly 30–45× Takahashi & Asada.
- **Revisit if:** a filter with r_min near θ_E, deeper shape catalogues, or a source-density increase (e.g. a
  second band's shapes) is proposed.

## D-104 W3 on ZTF DR24 light curves: detections-only data; feasibility trial failed its frozen validation twice (2026-10-18)
- **Decision:** ZTF DR24 light curves (IRSA HATS on S3, read column-wise) support a W3 (negative-mass n = 1) umbra
  search only through non-detections inferred from quadrant coverage. The frozen pipeline (`scripts/w3_ztf.py`, 1,000
  isolated stars near RA 272.65°, Dec +64.95°) failed its validation in v1 and, after bug fixes with unchanged
  thresholds, in v2 (W3 recovery 0.24 at ρ = 0.1, t_E = 100 d; 9 real development survivors > 5). No holdout search,
  no limit.
- **Alternatives rejected:** IRSA's per-object light-curve API (a 50-s request returned nothing here) in favour of
  one HATS partition read; a full-vanish flux calibration (the release has no upper limits); changing thresholds
  after the validation; a third version within this round.
- **Evidence:** `results/w3_ztf/RESULT.md`: v2 controls pass vetting at 0–1.5 % (eclipses, fades, PSPL; common mode
  by construction); the 9 survivors fall in two windows where ~10 % of development stars are mostly undetected.
- **Revisit if:** a sample-wide shared-epoch test and a correct comparison dim fraction are frozen first and
  validated on a new development set, an edge test that allows ρ t_E-long edges is defined without the injected
  results, or ZTF forced photometry (fluxes at every epoch) replaces the detections-only release.
- **Addendum 2026-10-18 (v3, `results/w3_ztf/v3/RESULT.md`):** exposures from IRSA's `ztf_current_meta_sci`, unique
  per-exposure counts, the target excluded and an in-play denominator (regression-tested). Survivors 9 → 8, recovery
  at ρ 0.1, t_E 100 d 0.24 → 0.30; still failing. Images of 23 exposures (1.34 MB) of 3 of the 8 survivors: the 11
  catalogue non-detections are 6 no-image/invalid/edge and 5 present in a neighbouring quadrant without an objectid;
  0 real drops. For these, a catalogue non-detection was not a disappearance; the rest is not established. Next, only as a new
  independent validation: footprint containment from the metadata corners, cross-quadrant objectid merging, and IBE
  forced photometry of survivors' non-detections.
- **Addendum 2026-10-18 (v4, `results/w3_ztf/v4/RESULT.md`):** coverage from each exposure's quadrant footprint
  (≥ 10″ inside IRSA's WCS corners, gnomonic) added to the time-span rule; ambiguous groups dropped; expid dedupe.
  On 1,000 new stars: 0 flags (2,565 under the v3 mapping; images of 5 removed flags: 0 dimming in 38 non-detections). Still fails the
  frozen rule at ρ 0.1, t_E 100 d (0.405; sharp-edge test in a post hoc diagnostic). Revisit if: the owner accepts a slow-edge-aware vetting
  test calibrated on injections before any new data are read.
- **Addendum 2026-10-18 (v5, `results/w3_ztf/v5/RESULT.md`):** an exploratory holdout search with the #339 code
  and the capability limited to the 10 cells at ≥ 0.5. 486 groups, 0 flags. The holdout was not strictly blind
  (aggregate comparison use in v1–v3). No abundance or limit. Revisit if: the owner authorizes another batch or a
  slow-edge vetting change.
- **Addendum 2026-10-18 (v6, `results/w3_ztf/v6/RESULT.md`):** a cross-sky validation in Npix 7134 (fields
  disjoint from v1–v5) fails: all 10 calibrated cells < 0.5 (48–90/200). In the one diagnosed t_E 30 d cell
  (ρ 0.01, blend 1), 82 of 112 unflagged injections have < 3 umbra nights (post hoc). The v4 capability did not
  reproduce in Npix 7134. The search group was not read. Revisit if: the
  owner wants capability conditioned on cadence (e.g. per-star umbra coverage), decided before new data are read.
- **Addendum 2026-10-18 (v7, `results/w3_ztf/v7/RESULT.md`):** the route is **paused**. The v4→v6 drop is about
  654–704 of 756 injections from sampling: DR24 light curves lack most 2022-10 to 2023-10 exposures in Npix 7134, and
  2024 is sparse. About 52–102 come from pre-screen and vetting given eligibility. No bug changes the outcome. Reopen if: a later
  release covers those epochs (`month_script.py.txt` check on calibration data first), or the owner authorizes
  image-level forced photometry for the post-event baseline. Not reopened by another batch or area of DR24.
- **Addendum 2026-10-18 (v8, `results/w3_ztf/v8/RESULT.md`):** for field 672 in Npix 7134 the gap exposures are absent from DR24
  in both HATS and the IRSA light-curve API (13 objectids checked; consistent with release-notes advisory 10; the
  reason for each exposure is not established). The #338/#339 cutout photometry fails its control
  (0.54–0.75 × DR24 flux, no field stars), so it cannot fill baselines as is. The exposure mid-time bug (exptime) is
  fixed, with no effect on stored results. **Paused.** Reopen if a newer release fills the field 672 gap (checked
  with `entry_script.py.txt`), image photometry passes a control against DR24 (owner go-ahead needed), or the owner
  revisits the infobits rule under a new plan.
- **Addendum 2026-10-19 (v9, `results/w3_ztf/v9/RESULT.md`):** for 3 stars in field 672 / CCD 8 / q2 (Npix 7134),
  IRSA's single-exposure psfcat catalogues pass the frozen overlap test against DR24 (13 pairs; per-star Δ −0.015 to
  +0.020 mag). They give ≥ 20 usable gap observations per star. This shows consistency with DR24, which is built
  from the same PSF photometry; it is not an independent validation, and it is not tested beyond this quadrant.
  Next-step preference for this gap: psfcat before ZFPS or cutout aperture photometry. Revisit if: the owner
  authorizes a psfcat data path and a calibration re-run.
- **Addendum 2026-10-19 (v10, `results/w3_ztf/v10/RESULT.md`):** psfcat was used as the data path for all
  validated DR24-gap exposures of Npix 7134 (32 of 33 combinations pass the overlap test). On the same 2,000 v6
  events, recovery rises from 711 to 872, but only 2 of 10 cells reach 0.5, so the frozen validation still fails and
  the search was not run. 527 of the 1,128 lost events (47 %) have fewer than 3 official infobits-0 r-band exposures
  in the umbra. Decision: keep psfcat as the gap-fill path if W3/ZTF is reopened; do not retune thresholds. Revisit
  if: a longer-baseline release or another survey adds epochs inside these gaps, or the owner changes the model range
  (e.g. longer t_E only).
- **Addendum 2026-10-19 (v11, `results/w3_ztf/v11/RESULT.md`):** in Npix 3839, chosen by official observing
  conditions first (worst-cell fraction 0.7025), a new 1,000-group sample passes the
  frozen validation (10 of 10 cells ≥ 0.5), and one search of 488 groups gives 0 flags. Decision: for any further W3 /
  ZTF batch, run the metadata check (`obs_script.py.txt`) before reading light curves, and use only areas whose
  worst-cell fraction is ≥ 0.5. Npix 7134 stays excluded. Revisit if: the owner wants a limit (not derived here) or
  another model range.
- **Addendum 2026-10-19 (v12, `results/w3_ztf/v12/RESULT.md`):** a new-field area chosen the same way (Npix
  2382; sampling indicator 0.6965) fails the frozen validation in one cell. At ρ 0.1, blend 0.7, t_E 30 d it recovers
  84/200; 80 of the flagged losses fail the sharp-edge test first, with 171/200 flagged. In Npix 2382 the
  observing-condition check did not predict recognition: the failing cell has a necessary fraction of 0.828 but
  recovers 0.42. Decision: no further W3 / ZTF batches under the current rules. Revisit
  if: the owner authorizes a re-specified edge test for finite-source events, under a new plan and a new validation.
- **Addendum 2026-10-19 (v13, `results/w3_ztf/v13/RESULT.md`):** the edge test is re-specified as edge
  *duration* ≤ 0.25 × run span. The threshold comes from the W3 finite-source model (`edge_frac.json`), and the rule is
  selected by `--edge time`; the count rule stays the default for earlier versions. It passes a new validation in
  Npix 2382, with fade-control contamination 8.5 % (limit 10 %). One search of 2,000 groups gave 0 vetted survivors.
  All 17 auto-classified "dimming" non-detections coincide with a dark detector column within 2″ of the source,
  against 0 of 20 detections (post hoc). Revisit if: the fade-control rate exceeds 10 % on another set, or a column-mask check is added to the
  catalogue stage (it would remove the spurious runs earlier).
- **Addendum 2026-10-19 (v14, `results/w3_ztf/v14/RESULT.md`):** the frozen v13 rule reproduces its capability on
  exposure-independent data (Npix 8170: cells ≥ 0.57; fade 5.0 %, Wilson 3.8–6.5 %, n = 1,000). One search of 5,000
  groups gave 1 vetted survivor. Images show it is a near-limit g-band catalogue drop, not a dimming. The 3 newly checked #370 class-4 NDs sit on official mask DEAD + LOW RESPONSIVITY pixels; the other #370 NDs were not masked.
  Revisit if: the owner wants the pre-screen to treat an ND as dim only when the source is well above the exposure
  limit (a method change; not made here).
- **Addendum 2026-10-19 (v15, `results/w3_ztf/v15/RESULT.md`):**
  - **3837:** the frozen v13 method passes validation and is searched: 4,778 groups, 0 vetted survivors.
  - **183:** fails validation. Its short-event cells are sampling-limited, and fade passes vetting at 15.4 %.
  - **Fade rate depends on the area:** 5.0–8.5 % in the three passing areas, 15.4 % in 183.
  - **Decision:** keep the per-area validation gate before any search; the metadata sampling margin (183: 0.547) is
    a weak predictor of passing.
  - **Revisit if:** a larger batch is wanted, then choose areas with a worst-cell margin ≳ 0.65 first.
- **Addendum 2026-10-19 (v16, `results/w3_ztf/v16/RESULT.md`):**
  - A frozen 8-area queue ran with the v13 method: 3 areas passed and 5,995 groups were searched. The 2 survivors are
    consistent with a low-response detector column and with poor-seeing images (no confirmed dimming).
  - **Post hoc** (v15–v16, 10 areas): inventory worst-cell sampling fraction (a sampling estimate, not the injected
    recovery rate) ≥ 0.593 passed and ≤ 0.578 failed. Not a calibrated rule.
  - **Decision:** as a heuristic, future queues may rank by that sampling fraction (≥ 0.59) before cost; the
    per-area validation gate stays. The remaining identified supply is
    4 partitions (17,702 searchable groups, 5.76 GB, shared fields).
  - **Revisit if:** a new partition list or a later DR adds unused areas.
- **Addendum 2026-10-20 (v17, `results/w3_ztf/v17/RESULT.md`):**
  - Post hoc evidence check of the two v16 survivors: both are **undetermined**.
    - 1624: a data problem is supported, not proven (a LOW RESPONSIVITY column within 1.70 px of the source in all
      9 NDs, ≥ 2.11 px away in the detections; split post hoc).
    - 1625: consistent with a wide-seeing blend, not tested to a decision.
  - **Decision:** neither is a candidate. No rule changes, and the column relation is not made a veto.
  - **Revisit if:** a mask-aware PSF fit or difference-image forced photometry becomes available. That would turn
    the unknown NDs into measurements.
- **Addendum 2026-10-20 (v18, `results/w3_ztf/v18/RESULT.md`):**
  - **Reuse:** a mask-aware PSF fit, using the official ZTF per-exposure DAOPHOT PSF (reader adapter) with
    photutils `PSFPhotometry`.
  - **Validated NDs:** 1 of 9 (1624) and 1 of 3 (1625). Both are normal and exclude a 1-mag drop in that exposure.
    1625's validation used its own PSF (not independent) and has no PSF-mismatch term.
  - **Not validated:** the rest stay unknown.
    - Masked-core injections (injected PSF from another exposure, fitted with the own PSF) are off by −17 % to +7 %.
    - Blend fits with a seeing-mismatched independent PSF give recovered/true from −4.6 to 3.0.
  - **Decision:** both sources stay undetermined and this follow-up is closed. Rejected: tuning the PSF choice
    after the validation failed (owner rule).
  - **Revisit if:** a seeing-matched independent PSF, or difference-image forced photometry, becomes available.
- **Addendum 2026-10-20 (v19, `results/w3_ztf/v19/RESULT.md`):**
  - **Finding:** v18's validation mostly tested cross-exposure PSF mismatch. It is relabelled a stress test.
  - **Decision:** validate with the measurement's own conditions. That means:
    - the own-exposure PSF;
    - real controls at a matched sub-pixel phase;
    - real-stamp injections with their noise counted;
    - a reference from validated detections only;
    - a cross-exposure repeatability term.
  - **Result:** NDs passing the corrected validation are 1 of 9 (1624) and 2 of 3 (1625), all normal. The 1625 /
    214647008 pass is conditional (dim-level ρ 0.92). The v18 1624 / 206251129 claim is suspended: it fails a V1b′
    criterion whose injection σ is narrower than the measurement's (a known defect). Both sources stay
    undetermined.
  - **Rejected:** iterating further methods on these two sources (owner rule).
  - **Revisit if:** difference-image forced photometry or the official PSF-star lists become available, or a rerun
    fixes the V1b′ σ and the per-level ρ.
- **Addendum 2026-10-20 (v20, `results/w3_ztf/v20/RESULT.md`):**
  - **Run:** the 9 remaining v16-inventoried partitions, with the frozen v13 method.
  - **Result:** 2 of 7 started pass (8126, 2900), 7,724 groups searched, 1 undetermined survivor (a near-limit
    catalogue miss is the untested leading explanation).
  - **Decision:** 1624 and 1625 stay archived as undetermined; reopen only with new data or an already-validated
    method. Failed areas stay unsearched. The #385 worst-cell pattern is not used as a gate.
  - **Revisit if:** the budget allows 8124 and 8127, or a new inventory is authorized.
- **Addendum 2026-10-20 (v21, `results/w3_ztf/v21/RESULT.md`):**
  - **8124 and 8127:** both pass; 9,978 groups searched.
  - **The v20 survivor:** its NDs fail the frozen depth and completeness criteria (0.32–0.60 mag below the
    baseline; 71–89 % catalogue completeness for comparable stars).
  - **The new survivor:** it sits on an official low-responsivity column (consistent with it, not proven).
  - **Decision:** both survivors stay undetermined, with no new photometry. The column relation is still not made a
    global veto.
  - **Revisit if:** a new inventory is authorized.
- **Addendum 2026-10-20 (v22, `results/w3_ztf/v22/RESULT.md`):**
  - **Inventory:** a new one of 31 partitions; 25 qualify, none sharing fields.
  - **Queue:** of 8 areas only 1 passed (1,882 groups searched, 0 survivors). 7 failed the frozen fade control at
    10.2–13.8 %.
  - **The 8124 neighbour:** in the shared NDs it also sits on a low-responsivity column (each within 0.5 px; column
    identity not stored); 8124 stays undetermined.
  - **Decision:** no change to the rule this round. The fade control was the only failing condition in 3 of the 7
    failed areas; the other 4 also failed short-t_E cells. A similar excess occurred in v20.
  - **Revisit if:** the owner wants a method decision on the fade control, such as calibrating it per field or
    diagnosing the excess. That is not done here.
- **Addendum 2026-10-20 (v23, `results/w3_ztf/v23/RESULT.md`):**
  - **Diagnosis:** 129 / 138 (4991) and 76 / 81 (2344) of the passing fades rely on an edge the data do not
    bound. `edge_time()` records it as 0 d.
  - **The fix tested:** the pre-registered "bounds" rule (`edge_bounds()`, `--edge bounds`) removes every fade pass
    but recovers only 43–114 / 200 W3 per cell.
  - **Decision:** v13 (`--edge time`) stays the search rule; bounds is kept as an option, not adopted. No second
    rule this round.
  - **Revisit if:** the owner accepts a lower recovery target for a measured-edge rule, or edge information is
    added (e.g. psfcat fills for unusable exposures).
- **Addendum 2026-10-20 (v24, `results/w3_ztf/v24/RESULT.md`):**
  - **Finding:** the undecided edge sides are mostly without a further covering r-band exposure (64 %) or quality-flagged
    (27 %). Candidate psfcat
    DR24 gaps are 8 % (W3) and 12 % (fades).
  - **Gain:** the optimistic gain moves no cell across 100 / 200 and leaves the 4991 fade > 100 / 1,000.
  - **Decision:** existing-product restoration (the v9–v10 psfcat route) is closed for the v13 and bounds rules.
  - **Revisit if:** new edge information exists (another survey or band covering the same stars), or the owner
    decides on using infobits ≠ 0 exposures or on a lower recovery target.
- **Addendum 2026-10-20 (ASAS-SN, `results/w3_ztf/asassn/RESULT.md`):**
  - **Trial:** ASAS-SN Sky Patrol, via the official `skypatrol` client, is the edge-information source tried.
  - **Outcome:** its light-curve service was down (maintenance), so the trial did not run.
  - **Kept for the re-run:** 40 stars selected by catalogue (bright, clean 16″ apertures) and their ZTF-only
    baseline.
  - **Revisit if:** `get_block_servers` is non-empty. Re-run the committed plan unchanged.
- **Addendum 2026-10-22 (alternative data, `results/w3_ztf/altdata/RESULT.md`):**
  - **Compared:** ZTF g (IRSA metadata + DR24), NEOWISE-R single exposures and TESS, on the same 40 stars and gap
    windows. W3 sides covered: 36, 3 and 0 of 88. Not queried: Gaia DR3 / PS1 (time), ATLAS (account), ASAS-SN
    (brief).
  - **Tested:** ZTF g as added edge points (`edge_bounds_alt`) on the replayed v23 events. 3 of 63 insufficient W3
    events restored (0.048; required ≥ 0.50); fade passes 0 / 79; success false; coverage ceiling 0.30.
  - **Decision:** ZTF g edge information is not adopted, and this route stops. Rejected: lowering the 0.50 target;
    own TESS FFI photometry (21″ pixels; no window covered); NEOWISE (3 sides).
  - **Revisit if:** a nightly-cadence source without an account becomes reachable, or the owner wants a rule in which
    other-band dim points extend the run. That would be a method change with a new validation; the post hoc
    diagnostic puts 59 / 66 window g epochs inside the event.
- **Addendum 2026-10-22 (g + r jointly, `results/w3_ztf/joint/RESULT.md`):**
  - **Tested (worth check only):** one candidate in which both bands bound the event (`joint_extent`,
    `edge_bounds_joint`). The g dims inside the simulated events are dim in the injected model (64 / 64; simulated-injection measurements, not
    real events) and decide 14 / 88 sides, against 5 for #459.
    Of the 9 new: 2 from the bracket only, 3 from the span only, 1 from either, 3 need both.
  - **Ceiling:** 7 / 63 insufficient W3 events restored (noise-free replay). A generous post hoc oracle gives 16 / 63.
    Both are below 0.50.
  - **Decision:** stop. No rule is frozen; bounds and v13 are unchanged.
  - **Rejected:** a second window or threshold, and using non-adjacent g epochs (a different method).
  - **Revisit if:** a source with epochs inside the shared r / g gaps exists, or the owner wants a model-fitting edge
    test (new plan and validation).
- **Addendum 2026-10-22 (ATLAS forced photometry, `results/w3_ztf/atlas/RESULT.md`):**
  - **Access correction:** ATLAS forced photometry is reachable without our own account through the Astro Data Lab
    token, which NOIRLab's ANTARES devkit notebook documents for use (shared queue). It worked on 2026-10-22, with 4 / 4
    tasks. The earlier "needs an account" (#459) is narrowed to: shared access exists; bulk use needs a personal
    account.
  - **Measured on 4 stars and 19 ZTF r / g common-gap sides (17 distinct gaps):**
    - reduced-mode photometry passes the FAQ cuts at 91–94 %;
    - 11 gaps have a passing precise visit in the decisive window;
    - 6 sides (5 distinct seasonal gaps) have none.
  - **Not measured:** the effect on decisions (no paired injection; no ZTF replay this round).
  - **Decision:** no rule change and no 40-star run without the owner. Rejected: difference-mode flux as a brightness
    ratio (template changes).
  - **Revisit if:** the owner authorizes ATLAS for the 40 stars and a paired injection with a ZTF replay and a frozen
    ATLAS noise rule.
- **Addendum 2026-10-22 (ATLAS paired injections on the 40 stars, `results/w3_ztf/atlas40/RESULT.md`):**
  - **Run:** the shared queue ran one task at a time (2–19 min each plus 2–22 min wait). By the frozen deadline 17 / 40
    stars had data, so the trial is **not completed**.
  - **Result:** with ATLAS points added as edge points, 6 / 63 insufficient W3 events are restored, 4 recovered W3 are
    lost, and fade passes are 0 / 79. Even if every insufficient event on the 23 missing stars were restored, the
    best case is 27 / 63 < 0.50. Verdict: **fails regardless of the missing stars.** (#470 scope: 27 / 63 is the
    best-case upper bound under this fixed method, data and A decisions, not a general ceiling; 0 / 79 holds with
    ATLAS on 17 stars, and 0–4 / 79 is the range for all 40 under this method.)
  - **Mechanism (post hoc diagnostic):** the restorations come from added baseline epochs. 7 of the 8 W3 → "slow"
    moves come from single observed ATLAS dips (32–79 % below baseline) that pass the FAQ cuts. Quoted σ and the 9999
    rule change nothing.
  - **Decision:** ATLAS as added edge points does not rescue the bounds rule on this sample under this criterion. No
    independent validation; no seed, threshold or second-method retry.
  - **Rejected:**
    - finishing the remaining 23 stars, which cannot change the verdict;
    - a post hoc two-epoch intermediate rule, which would be a new method.
  - **Revisit if:** a new pre-registered method that is robust to single-epoch survey dips is proposed with an ideal
    bound that can reach the success condition.
- **Addendum 2026-10-22 (ATLAS dip check, post hoc, `results/w3_ztf/atlas_dips/RESULT.md`):**
  - **Found:** 7 of the 8 simulated W3 mis-decisions of #462 trace to 8 low ATLAS rows from 6 exposures (4
    unit-nights). The 8th traces to a normal row where the model itself was intermediate.
    - One exposure holds 3 of them; there all 10 stored area-4991 stars are 32–34 % low (exposure-level photometric
      clue).
    - One night shows several stars low.
    - One single-exposure dip is undetermined.
    - One group is unchecked: cap of 3, and a selection-count deviation that would have chosen it over G4.
  - **Decision:** no real-variability evidence and no source to follow up in the 3 checked groups. The line is archived
    as is, with no new rule, gap filling or photometry request.
  - **Note:** a common-mode test across stars in the same exposure would flag G1. Adding it would be a new
    pre-registered method. Check the cleaning guide the ATLAS FAQ cites (Rest et al.) first.
  - **Revisit if:** the owner asks for the unchecked G2 group (done: post-review G2 check below), opens a new method
    round, or ATLAS exposure images become reachable without an account.
  - **Post-review G2 check (2026-10-22, PLAN §3):** G2 is the group the selection deviation left out.
    - One exposure shows common mode.
    - Two rows are undetermined: their comparison stars scatter > 5 σ in both directions on that night.
    - No individual variability is supported, and no source is worth follow-up from these rows. All 4 groups are
      now checked.
    - The line is closed with no further check of this kind. Revisit only on the owner's request for a new method
      round, or if exposure images become reachable without an account.

## D-105 W1 blank-field R from real-image metacalibration: borrowed R = 0.414 kept (2026-10-18)
- **Decision:** keep R = 0.414 for the blank-field W1 screens. GalSim metacalibration of six GOODS-N F200W i2d tiles
  with the pipeline-like shapes gives R = 0.953 ± 0.239, "conservative" under the frozen rule (R ≥ 0.414 − 2σ).
- **Alternatives rejected:** adopting the measured R (its 2σ lower bound is 0.475); ngmix metacal (fits its own
  shapes, not the screen's isophotal moments; D-101); STPSF as PSF_d (an analytic Airy ⊗ Gaussian was frozen first
  for speed; an empirical drizzled PSF is the better next step); the post hoc tile exclusion as the headline.
- **Evidence:** `results/w1_r_realmetacal/RESULT.md`, `r_real.json`, `validate.json` (simulated gate passed, weak),
  `validate_pooled_posthoc.json`, `r_real_posthoc_excl_tile.json`, `sheet.png`.
- **Revisit if:** an empirical-PSF, masked metacal on more tiles gives R < 0.414 − 2σ (re-run the limits lower), or
  its 2σ lower bound exceeds ~0.8 (a frozen plan to re-run the limits with it).
- **Addendum (2026-10-18, v2):** noise-filled mask, PSF σ floor 0.97 px, unique matching, a frozen tile gate and all
  22 coverage-passing tiles (4 excluded, all ringing after deconvolution): R = 1.056 ± 0.156 (tile bootstrap; 804
  sources); R − 2σ = 0.744 < 0.8, still "conservative"; R = 0.414 kept. The PSF_d × 1.15 variant (0.875) shifts R by
  as much as the error, so the empirical PSF is the next step before more area (`results/w1_r_realmetacal/v2/RESULT.md`;
  post hoc output-PSF floor 1.37 px: 1.114 ± 0.168, same call).
- **Addendum 2 (2026-10-18, v3):** an empirical PSF_d from Gaia DR3 stars of the same mosaic fails its frozen gate
  (2 usable stars of 6; ≥ 5 needed), so v3 gives no R. Post hoc, the two-star ePSF fails the simulated gate (V1) and
  gives R = 0.707 ± 0.169 against 1.006 ± 0.134 for the analytic PSF_d in the same run. R = 0.414 kept; the route is
  paused, because the PSF model moves R by about 0.3 (`results/w1_r_realmetacal/v3/RESULT.md`). Revisit with a
  validated drizzled STPSF or a field with ≥ 5 usable stars.

## D-106 W3 inside a macro model: in-house numpy inverse ray shooter; saddle images show no W3 vanishing in the sampled windows (2026-10-20)
- **Decision:** `exotic_sim.macro_magnification_map` (numpy inverse ray shooting, fine grid near the lens, coarse
  grid scaled to the macro stretch, ray-area weights) models a point lens of either sign inside a macro (κ, γ).
  Result (simulated): minimum images keep a zero-flux umbra, narrowed to an on-axis half-width 2√λ1 along the
  stretch (not elliptical off axis); saddle images have none in the sampled windows.
- **Alternatives rejected:** lenstronomy (point mass parametrised by θ_E², no negative sign without a custom
  profile; a new dependency for ~60 lines); the isolated-lens formulas for arcs (wrong for saddle images, too
  wide along the stretch for minima).
- **Evidence:** `results/w3_macro_shear/RESULT.md`, `result.json`, `maps.png`; gates G1, G2, G3′ pass (G3 as written
  was mis-specified, PLAN addendum 1); `tests/test_w3_macro_shear.py`.
- **Revisit if:** a W3 search targets caustic-transient stars (then add a microlens population κ_* and λ2 ≠ 1), or
  a published treatment of negative-mass microlensing with external shear disagrees.

## D-107 W3 on Kepler DR25 long-cadence SAP: pilot validation passed, 60-star search null; shoulder test needs hardening before any scale-up (2026-10-24)
- **Decision:**
  - For W3 (finite-source negative-mass point lens), Kepler long-cadence SAP light curves are a
    working route for t_E 0.1–10 d at the pilot level: injections are recovered at 78–91 / 100 per
    cell, and recovered edges are bounded with a median of 0.041 d.
  - The vetting on this route adds a **shoulder test** (a brightening just outside the umbra, on
    both sides). It is the only test that separates W3 from a single box eclipse.
  - Edges are bounded per quarter: a cross-quarter edge is not measured.
- **Alternatives rejected:**
  - The ZTF chain unchanged (single eclipses then pass 85–91 / 100).
  - PDCSAP as the primary flux (the owner specified SAP; PDCSAP is a diagnostic from the same
    observation).
  - t_E ≥ 30 d with the per-quarter join (the umbra spans quarters).
- **Evidence:** `results/w3_kepler/RESULT.md`.
  - The 0 / 100 control rates follow from the control models' design (RESULT §2).
  - 4 of 100 real SAP light curves already meet the shoulder criterion within a quarter
    (`posthoc2.json`).
  - ZTF DR24 under the bounds rule: 43–114 / 200 per cell, gate failed (D-104).
- **Revisit if:** before any larger Kepler search, harden the shoulder test: several cadences per
  side, a local in-quarter baseline instead of the quarter median, and a control with real
  in-quarter brightenings and long dips. Also revisit if long t_E is needed (a different quarter
  join first), or if a second band or pixel-level injection becomes available.
