# SMACS 0723: Gaia DR3 proper motion does not change any Gaia-rule star label

Plan: `PLAN.md` (written before the query). Script: `scripts/gaia_pm_epoch.py`. Numbers: `summary.json`,
`fast_movers.ecsv`. Run 2026-11-06 on the jwst 3.0.0 F200W pipeline catalogue
`jw02736-o001_t001_nircam_clear-f200w_cat.ecsv` (S3 `stpubdata`, sha256 in `summary.json`; not stored in git).

- **Input (observed / assumption).** One Gaia DR3 TAP cone (centre and radius in `summary.json`) holding 103 Gaia
  sources, 96 with proper motions. JWST epoch 2022.43 (ASSUMPTION: the ERO observation date, 2022-06),
  Δt 6.43 yr. Linear propagation, parallax ignored (derived).
- **Frame check (derived).** The 22 slow Gaia sources (shift < 0.05″) within 0.5″ of a catalogue source sit a median
  0.028″ from it. The three fast stars inside the footprint fall to 0.009–0.019″ after propagation and were
  0.15–0.34″ before, so the pipeline catalogue is on the observation-epoch frame and the propagation is right.
- **Fast movers (derived).** 4 sources shift > 0.15″, 1 shifts > 0.3″, 0 shift > 0.5″. The one > 0.3″
  (Gaia DR3 5262899199049945728, G 17.85, 54.5 mas/yr, 0.350″) lies 0.343″ from catalogue source 4025 at J2016:
  outside the 0.3″ `gaia_position` radius, inside the 0.5″ query radius, and an astrometric star (PM S/N 480 ≥ 5),
  so the Gaia part of the D-012 rule labels it a star with or without propagation. A fourth fast source lies 64″
  from any catalogue source (outside the footprint; the cone is a circle around a square mosaic).
- **Decision (frozen rule).** No star with shift > 0.3″ loses its label, so no propagation is built into
  `crossmatch` now. A label could be lost only for μ > 0.5″ / 6.43 yr = 78 mas/yr; this field has none.
- **Limits.** Only the Gaia part of the rule is checked: a SIMBAD galaxy nearer than 0.343″ to source 4025 would
  win under D-012 and was not queried (post hoc caveat). Only F200W positions; one field. Revisit for fields at
  low Galactic latitude, with many high-proper-motion stars, or observed > 10 yr after J2016 (shift > 0.3″ already
  at 30 mas/yr).
