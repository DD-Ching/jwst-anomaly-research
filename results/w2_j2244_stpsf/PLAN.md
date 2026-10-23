# W2 J2244−0550 follow-up: minimal STScI official PSF route, capability check (WFC3/IR F160W only) — plan

Status: claim and scope (2026-10-23 15:55 UTC). The test plan (cases, settings, comparison quantities, success and
stop conditions) is added here and committed before any comparison is run.

- **Question:** is there a directly reusable official STScI PSF route that improves the point-source modelling that
  failed in #474 (`results/w2_j2244_hst/RESULT.md`)?
- **Type:** development feasibility test. Not an independent validation; J2244−0550 is not measured.
- **Hard limits:** wall clock 2 h (to 17:55 UTC); 2 cores; total download 1 GB including calibration files and new
  dependencies; WFC3/IR F160W only; no ACS; no multi-field star library; at most 3 small cases.
- **Download control before any image request:** a persistent, shared ledger (reserve before each request, record the
  bytes actually received, retries and repeats counted, survives restarts); verified with one small request first.
  No download if the remaining budget cannot be metered or reserved.
- **One official entry, read for at most 30 minutes and pinned to a version:**
  https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/point_spread_function/hst_point_spread_function.html
- **Out of scope:** #474 numbers and verdicts (unchanged); DESI classifier, X4, W3 / ATLAS and other paused units; no
  outside contact, accounts or schedules.
