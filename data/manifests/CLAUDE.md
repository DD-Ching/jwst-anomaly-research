# data — data root, downloads and provenance

Root CLAUDE.md applies.

- Data root is `$JWST_ANOMALY_DATA` (default `<repo>/data`). Everything under it except `data/manifests/` is
  gitignored. Never commit data; no tracked file over 1 MB.
- Manifests (URI + sha256 + size + retrieval time + pipeline version) are the reproducibility record; every archive
  product a run uses gets a row, and SOURCES.md says how to recover the source (URL, version or DOI, access date).
- Owner's machine: single downloads > 200 MB need a stated reason; disk is limited. Don't download a full NIRCam
  `_i2d.fits` (~1.8 GB) when an S3 byte-range cutout suffices.
- Cloud sessions (owner decision 2026-10-08): stream data, never store a whole archive tar, log the reason for any
  large transfer, delete raw data after use.
- Prefer server-side aggregation (TAP, byte ranges) over bulk downloads; batch network I/O (one query per program,
  CDS XMatch for many sources).
- A MAST token comes only from env `MAST_API_TOKEN`. Never commit secrets.
