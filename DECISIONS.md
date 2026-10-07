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
