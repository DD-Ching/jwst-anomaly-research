# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.

## 2026-10-07 — M0 bootstrap started
- Created repository skeleton: interfaces (`schema.py`, stage stubs with fixed signatures), state files,
  CI, BSD-3-Clause license, verbatim agent charter (`docs/agent-charter.md`) imported by `CLAUDE.md`.
- Verified on MAST: program 2736 level-3 NIRCam/MIRI imaging is public; catalogs ~3 MB vs NIRCam i2d
  ~1.8 GB, hence the catalog-first slice (D-001). The S3 mirror key pattern is confirmed.
- Owner decisions: public repo; every change via PR and the owner merges; chat in Chinese.
- **Handoff:** nine bootstrap units are being implemented in parallel (`batch/*` branches). Next comes integration
  on real data (TASKS.md "Next").
