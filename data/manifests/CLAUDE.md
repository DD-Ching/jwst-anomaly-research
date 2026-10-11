# data/manifests — provenance records

Root CLAUDE.md holds the data-root, download-size and cloud-disk rules; `scripts/CLAUDE.md` holds the streaming
topology. This file covers the manifests themselves.

- Every archive product a run uses gets a manifest row: data URI, file name, size, sha256, retrieval time (UTC) and
  pipeline version. The manifests, not the data, are the reproducibility record; SOURCES.md says how to recover each
  source (URL, version or DOI, access date).
- Adapters pin their inputs by sha256 against these manifests and refuse a mismatch.
- A streamed archive (cloud disk rule) is recorded by URI, size and the sha256 of each member or chunk read, plus
  the reason for the transfer.
- Prefer server-side aggregation (TAP, byte ranges) over bulk downloads; batch network I/O (one query per program,
  CDS XMatch for many sources).
