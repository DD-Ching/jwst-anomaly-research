# MACS1149 system 16: DJA NIRSpec table v4.6 re-check (plan, frozen before the table is read)

Why: `../RESULT.md` (2026-10-20) read DJA v4.4 (2025-09-04), where the 16.2 PRISM spectrum
`macs1149-v4_prism-clear_1208_5114276` had no grade. Newer tables exist on the same bucket: v4.5 (2026-02-23,
151,445,355 B) and v4.6 (2026-09-15, 171,878,708 B) (observed, HTTP HEAD 2026-10-21). Only v4.6 is read (the newest).

- **Hypothesis under test (unchanged):** the flagged extra image of system 16 is a model artefact of a wrong fitted
  redshift (D-037; CATS fits z 4.419).
- **Data:** `dja_msaexp_emission_lines_v4.6.csv.gz`, streamed and filtered on the fly (not stored; cloud disk
  rule), rows within 0.5″ (ASSUMPTION, as before) of CATS v4.1 16.1 (177.40971, 22.387662) / 16.2 (177.40989,
  22.387828), plus any row whose file name contains `5114276` or `5106117` (D-037's possible far image).
- **What counts as a spectroscopic redshift (ASSUMPTION):** DJA visual `grade` 3 (robust). Grade 2 is reported as
  tentative and decides nothing; grade ≤ 1, missing, or automated fits only → none.
- **Outcomes (frozen; the first two as in `../PLAN.md`):**
  - grade-3 z ≤ 3.5 for 16.1 or 16.2 → "fitted z wrong": the extra-image flag is explained (D-037 addendum).
  - grade-3 z with |z − 4.419| ≤ 0.1 → "fitted z confirmed": the flag reopens; re-run forced photometry for the
    predicted image at that z.
  - grade-3 z elsewhere → "undetermined: new z"; re-run the model z track at that z (`../../macs1149-cats/sys16_track.json`).
  - no grade-3 z → "not decided"; D-037 unchanged; record the v4.6 state of the rows.
