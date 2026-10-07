# Data manifests

Tracked, machine-readable records of every archive product the pipeline used: data URI, file name,
size, sha256, retrieval time (UTC) and pipeline version. The files themselves live under
`$JWST_ANOMALY_DATA` and are never committed; re-fetch them with `scripts/` and verify against these manifests.
