# Architecture

## Pipeline (v0, catalog-first)

```
MAST (CAOM metadata)            query.py      query_observations / list_products
  → reproducible acquisition    acquire.py    fetch_products → data/manifests/*
  → science-ready products      level-3 _cat.ecsv (photometry/shape) + _i2d.fits (images)
  → source representation       catalog.py    load_pipeline_catalog / merge_bands
  → matched-aperture colours    photometry.py join_matched_photometry (DJA, per sample, D-013)
  → derived features            features.py   build_features
  → quality gate (D-011/D-014)  quality.py    assess_sources: image tests (coarse WHT map via
                                              cutouts.sample_weight_map) + detection confirmation
                                              (best-band S/N, multi-band or matched-photometry match)
  → star/galaxy strata (D-012)  classify.py   classify_sources (bulk Gaia/SIMBAD) + stellar_locus (D-015/016,
                                              DJA r50 + colours); stars ranked as <sample>-stars
  → baseline anomaly ranking    rank.py       score_anomalies
  → image evidence (top-k)      cutouts.py    make_cutouts (S3 byte-range reads; `spikes` flag D-018, screening D-019/D-021, host test D-020); viz.contact_sheet
  → time domain (vetting)       scripts/transient_search.py + transient_combine.py + transient_forced.py (D-027)
  → external cross-check        crossmatch.py crossmatch (SIMBAD / NED / Gaia / CDS XMatch)
  → candidate store             candidates.py CandidateStore (SQLite) + provenance.py
  → orchestration               pipeline.py   run(config) ; cli.py `jwst-anomaly`
  → interpretation              humans + /vet-candidate skill (never automatic)
```

Why catalog-first: a NIRCam level-3 `_cat.ecsv` is ~3 MB while its `_i2d.fits` is ~1.8 GB
(program 2736, measured 2026-10-07). Ranking on pipeline catalogs and pulling image cutouts only
for the top-k keeps the first slice cheap. Known limitation: pipeline catalogs detect each band
independently (no forced photometry), so cross-band colors are approximate. See DECISIONS.md.

## Contracts

- Every stage takes/returns `astropy.table.Table`; required columns are in `src/jwst_anomaly/schema.py`.
- `meta["provenance"]` ∈ {observed, derived, simulated, model_prediction, assumption, hypothesis};
  `meta["source"]` names the inputs. `schema.validate()` checks both.
- Per-band columns: `schema.band_column(band, quantity)` → `f200w_aper50_abmag`.
- Locations: `paths.data_root()` (`$JWST_ANOMALY_DATA`), `paths.manifests_dir()` (tracked),
  `paths.outputs_dir()` (`$JWST_ANOMALY_OUTPUTS`).
- Public stage signatures are fixed; extend with keyword arguments that have defaults.
- One tracked manifest per config: `data/manifests/<config stem>.ecsv` (downloads) and
  `<config stem>_products.ecsv` (all level-3 products with S3 URIs), shared by
  `scripts/fetch_reference_sample.py` and `pipeline.run`.
- Field footprints come from MAST `s_region`, not `s_ra/s_dec` (which can be another
  instrument's prime target, e.g. CEERS t021).

## Stage ownership (bootstrap batch, 2026-10-07)

| Stage | Module(s) | Unit |
|---|---|---|
| Archive query + acquisition | `query.py`, `acquire.py`, `scripts/fetch_reference_sample.py` | 1 |
| Catalog ingestion | `catalog.py` | 2 |
| Features + baseline ranking | `features.py`, `rank.py` | 3 |
| Cutouts + visualization | `cutouts.py`, `viz.py` | 4 |
| External cross-check | `crossmatch.py` | 5 |
| Candidate store + runner + CLI | `candidates.py`, `provenance.py`, `pipeline.py`, `cli.py` | 6 |
| OSS hygiene | community files, templates, pre-commit | 7 |
| Agent operations | `.claude/`, `.github/workflows/claude.yml`, `docs/operations.md` | 8 |
| Reuse landscape | `docs/landscape.md` | 9 |
