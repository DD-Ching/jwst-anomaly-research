# VizieR CfA mirror fallback (plan, 2026-11-21)

- **Hypothesis:** `vizier.cfa.harvard.edu` returns the same rows as `vizier.cds.unistra.fr` for the catalogues the
  live tests read (J/A+A/607/A30 M0717, I/355/gaiadr3 cone), so an ordered CDS → CfA fallback restores those
  callers while CDS answers 404 (#819).
- **Ordinary explanations for a mismatch:** mirror lag (older catalogue version), different row limits or column
  defaults, different server-side cone rounding.
- **Gate (ASSUMPTION):** identical row counts and identical IDs / positions to the stored CDS-derived test
  expectations; otherwise the mirror is not used and the PR reports that.
- **Not covered:** cdsarc FTP files (no reachable mirror found), TAPVizieR, CDS XMatch.

## Result (2026-11-21)
Gate passed: on `vizier.cfa.harvard.edu`, M0717 4082 / 4082 rows, 2928 zout rows (CDS-recorded: 4082 / 2928);
the A370 Gaia cone test passes. CDS VizieR still 404 at 09:08Z.
