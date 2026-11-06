# SMACS 0723 jwst 3.0.0: 3 of 139 star-stratum sources have a NED galaxy-like match; no star rule changes

Plan: `PLAN.md` (pushed before the run). Run `20261106T192505Z-3470b4cf`, `configs/smacs_jwst3_simbad.yaml`, main at
`ea3fba8`. Script: `scripts/smacs_star_ned.py`. Tables: `matches.csv` (one row per star–NED pair, stars without a
match included), `summary.json`, `followup.json`.

## Run (observed / derived)
- The star stratum reproduces the SIMBAD-mirror run (`../simbad_mirror/RESULT.md`): 139 ranked stars, 36
  SIMBAD/Gaia astrometry, 6 Gaia position only, 97 stellar locus.
- **Deviation:** ESA Gaia TAP answered 503 most of the time (2 of 12 probes 200, not stored), so two pipeline
  runs failed classify by design. The Gaia `field` result was fetched once in a retry loop and the third run read it
  from the cross-match cache (log line "gaia: using cached field result"). CDS XMatch now answers but rejects
  `simbad` and `vizier:I/355/gaiadr3` as unknown tables (`ValueError`), so the chains fall through to the mirror and
  ESA TAP. Not investigated.

## Result (NED TAP, radius 0.5″, ASSUMPTION; galaxy-like types listed in the script)
| Star | Basis | Star rank | NED object | Type | Sep | NED z |
|---|---|---|---|---|---|---|
| f090w_1561 | stellar locus | 2 | `SMACS J0723:[NDA2023] 00713` | G | 0.030″ | 0.3828 |
| f200w_3670 | stellar locus | 24 | `SMACS J0723-73:[MJR2023] 002.3` | G_Lens | 0.369″ | 1.3715 |
| f200w_3588 | SIMBAD/Gaia astrometry | 68 | `WISEA J072311.11-732757.7` | G | 0.064″ | 0.24 |

Per basis (`summary.json`): astrometry 1 / 36, Gaia position 0 / 6, stellar locus 2 / 97 (1 within 0.3″). 18
NED `*` matches. The two 2026-10-07 examples of TASKS item 4 reappear with new jwst 3.0.0 ids (same NED objects,
`f090w_1719` → `f090w_1561`; `f200w_2054` → `f200w_3670`, now 0.369″ instead of the 0.23″ recorded in TASKS).

## Per source (ordinary explanations)
- **f200w_3588 (astrometric star, NED galaxy):** Gaia DR3 parallax 0.606 ± 0.039 mas (15.7 σ) and pmRA −6.84 ±
  0.05 mas/yr at 0.052″ (`followup.json`); F200W 16.3 mag. A significant parallax and proper motion make it a
  Galactic star; the NED "G" type of a WISE source is the weaker label. **Plan trigger** ("any astrometric star with
  a NED galaxy match → re-examine the star rule"): fired, and the re-examination keeps the rule (NED does not feed
  `is_star`, `crossmatch.py`). Image not viewed (rank 68).
- **f090w_1561 (locus star, NED galaxy at 0.03″):** NED z 0.3828 (redshift type not checked) is close to the cluster
  redshift (0.388, as used in `results/transients/smacs_at2024jfz/PLAN.md`); no Gaia or SIMBAD object within 0.5″.
  Contact sheet (stars top 10, #2) seen: round and compact, without the hexagonal PSF structure of the other nine
  panels. Reading: **likely a compact cluster member mislabelled by the stellar locus** (hypothesis; not re-ranked).
- **f200w_3670 (locus star, G_Lens at 0.369″):** outside 0.3″, so a neighbour rather than the same object is the
  narrowest reading; the recorded spikes and adjacent arc (TASKS item 4) not re-checked. Image not viewed.

## Reading
Neither plan threshold for a rule change holds: the astrometric case is a Gaia star, and 1 (< 10) locus star has a
NED galaxy within 0.3″. The locus mislabels at least one compact galaxy (f090w_1561); this is a known D-015 risk,
not a candidate. No re-rank was done (the plan's paired re-rank needed ≥ 10).

## Limitations
- NED type and redshift are taken as catalogued; NED redshift types and the `NDA2023` reference were not checked.
- Only f090w_1561 was viewed; f200w_3670 and f200w_3588 were not.
