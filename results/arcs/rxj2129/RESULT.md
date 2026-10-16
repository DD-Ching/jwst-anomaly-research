# `arcs` on RXJ2129: screen result (2026-10-16)

Plan: `PLAN.md` (committed 20:13Z, `bbe0306`, before the run at 20:21Z). Labels: shapes `observed`; predicted axes
`model_prediction`; classes, counts and P values `derived`. All thresholds are the D-024 defaults (ASSUMPTIONs).

## Inputs (`summary.json` `inputs`; map sha256 in `lensmodel.RXJ2129_*_MAPS`)
- F200W pipeline catalogue `jw02767-o002_t001_nircam_clear-f200w_cat.ecsv` (program 2767, 2,061 s; 1,030,198 B,
  sha256 `2ed1a2c7…261f53`), read from the public S3 bucket (the MAST download API returned an internal error).
- DJA grizli v7.0 `rxj2129-grizli-v7.0-fix.eazypy.zout.fits` (sha256 `651c0b24…726217`).
- Maps read from the Zenodo 21460864 zip by HTTP range requests (no full download): Lenstool (primary, D_LS/D_S = 1)
  and GLAFIC (second, normalised at z_s 3). `validate`: κ from the deflections vs the published κ, median relative
  difference 2.9 × 10⁻⁵ (Lenstool) and 1.2 × 10⁻⁵ (GLAFIC) (`*/validate.json`): both pass.

## Result (`summary.json`, `*/arcs.json`)
| Model | Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) | P(≥ anti) |
|---|---|---|---|---|---|
| Lenstool (primary) | 12 | 8 (67 %) | 2 | 0.22 | 0.021 |
| GLAFIC (second) | 12 | 8 (67 %) | 1 | 0.176 | 0.162 |

**Frozen outcome: no result** (both controls fail, 8 / 12 < 70 %). The P values above do not count: the rule gives
them no weight when the control fails.

## Post hoc (not in the plan; do not change the outcome)
- Lenstool anti rows (`rxj2129-lenstool-maps/arcs_anti.ecsv`): 447 (m 23.76, z_phot 2.87, g_min 0.35) and 464
  (m 23.84, z_phot 1.91, g_min 0.27), 8.1″ apart (not stored). GLAFIC calls 447 anti and 464 mixed.
- The sample is small: 52 selected rows, 12 in strong shear (`*/arcs.ecsv`), against 55 strong-shear rows on the
  6,399 s MACS1423 catalogue. The 2,061 s exposure is a plausible reason; not tested.
- F200W cutouts of 447, 464, two aligned and two mixed rows looked at (not stored): 447 is a compact elongated source
  on the halo of a bright galaxy; 464 is compact, about 1.5″ from a brighter compact neighbour; mixed 495 is a large
  clumpy spiral. The sheet's axis overlay used an unchecked PA convention, so no orientation is read from it.

## Reading
No result on the S1063 replication question from RXJ2129 under the frozen rule; the tally stays: MACS0417 does not
replicate, MACS0647 no result, MACS1423 explained by the clump-pair census, RXJ2129 no result. Not a candidate.
