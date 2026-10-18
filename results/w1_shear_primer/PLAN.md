# W1 shear screen in PRIMER (UDS and COSMOS, program 1837, no lens model): plan (frozen 2026-10-18, before any screen)

**Question.** The twenty-four-field W1 shear volume limit (`results/w1_shear_ceers/RESULT.md`) is 1.09× weaker than
Takahashi & Asada (2013) at 2 × 10¹³ M☉ over 199.9 arcmin² of screened area. PRIMER has four large F200W level-3
mosaics in MAST plus two combined ones. They are the next blank-field W1 area (TASKS "JWST focus" W1, after the
blank-field R calibration of #329, which found the borrowed R = 0.414 conservative).

**Hypothesis being limited.** As `results/w1_volume/PLAN.md` (`hypothesis`: W1 lenses uniform in comoving volume).
A screen plus a limit; a peak that passes the E/B rule is a flag to vet, not a candidate.

**What would change the plan.** If the combined V_eff brings n₉₅ at 2 × 10¹³ M☉ below 10⁻⁴ h³ Mpc⁻³ (the Takahashi &
Asada value), W1 has matched the published bound and the next step is lower masses, not more area. If the PRIMER
tiles add < 20 % to V_eff at every mass, blank-field additions stop (T4: cost per gain).

## Inputs (frozen)
- **Footprints** (`footprints_uds.json`, `footprints_cosmos.json`; MAST `s_region`, rasterised at 1.2″, queried
  2026-10-18 with `w1_footprints.py --program 1837 --ra-range`). MAST holds, per field, two large mosaics and
  shallow single-visit products (837 s, all inside the mosaics; not used):
  - UDS: `c1001_t000` (129.2 arcmin²) and `c1002_t000` (127.9) are the combined north and south mosaics; they
    contain `o002_t002` and `o001_t001` (overlaps 115.7 and 119.0 arcmin²) and overlap each other by 3.6 arcmin².
    **Kept: c1001, c1002.**
  - COSMOS: `o003_t003` (92.4) and `o004_t004` (87.4), overlapping by 38.7 arcmin². **Kept: o003, o004.**
- **Tiling (new; the mosaics exceed the screen's 200″ radius).** `w1_footprints.tile_grid`: square tiles of
  half-width 140″ (21.8 arcmin² each; half-diagonal 198″ < the 200″ screen radius) on one RA/Dec grid per field over
  the union of the kept mosaics. Each tile goes to the kept mosaic covering most of it and is kept when that share
  is ≥ 0.3. Only that mosaic's catalogue rows inside the tile are used (`inject_radial.in_box`, applied in
  `load_field` through the spec key `box_half`). Tiles do not overlap, so no sky is screened twice; sky of a tile
  that only the other mosaic covers is lost (conservative). Result: **12 UDS tiles** (share-weighted 204.8 arcmin²)
  and **8 COSMOS tiles** (129.8 arcmin²), listed under `tiles` in the two JSON files. Each tile is a field
  `primer-<uds|cosmos>-<tile>` with `NullLensModel` at the tile centre.
- **Catalogues.** The F200W `_cat.ecsv` of each kept mosaic (MAST; 21.5–31.3 MB), sha256 pinned in
  `w1_stage_inputs.py` at the first fetch.
- **Photo-z.** DJA eazy zouts streamed out of the tarballs, never stored whole (stated reason: the only photo-z
  catalogues of these mosaics in the format the loader reads; cloud disk rule; D-051 precedent):
  `primer-uds-north-grizli-v7.2` (393.0 MB), `primer-uds-south-grizli-v7.2` (344.9 MB),
  `primer-cosmos-east-grizli-v7.4` (407.9 MB, SOURCES.md), `primer-cosmos-west-grizli-v7.4` (536.4 MB). Sizes are
  HTTP Content-Length on 2026-10-18 (not stored). A tile whose rows fall in two DJA mosaics gets both zouts
  stacked before the 0.3″ match (`attach_photoz`); tarball and member sha256 are recorded at the first fetch.
- **Vetoes.** Gaia spike veto (D-043) at each tile centre (4′ radius, which covers the tile); streak masks from
  `w1_streak_mask.py` (frozen rule, `results/w1_streak_veto/PLAN.md`), read only for the tile's rows of the
  mosaic's `_i2d.fits` SCI extension: the files are 31–69 GB, so the whole-SCI stream the tool now does is not
  used for PRIMER (a tile-range read is needed first).

## Gate P (computable; before the screen)
`inject_shear.passes_gate_p` per tile: ≥ 50 % of the rows that pass the screen's S/N and resolution cuts have a
photo-z match. Fail → the tile is dropped and reported.

## Responsivity (ASSUMPTION)
R = 0.414 as in D-099 (Abell 370 F200W fit), shown conservative for blank fields by `results/w1_r_metacal`
(R_sim ≥ 0.94). Sensitivity (reported, not used in the headline): the injections at Z = 0.3 repeated with R = 0.231.

## Screen (computable)
As `results/w1_shear_ceers/PLAN.md`: Schirmer 10″, 1″ grid, rotation null, E/B rule, streak veto, screen radius
200″ (the tile's rows bound the area). **Screen rule:** unchanged E/B rule (p_random < 0.05 and S > the tile's max
|S_×|). A flagged tile gets a contact sheet and the frozen vetting rule (`w1_shear_more` addendum A, V1 ∪ V2), and
stays out of the headline either way (D-098). Peaks within 20″ of a tile edge are flagged as such (edge rows are cut
by the tile, not by the data).

## Injections (frozen)
`inject_shear.py --lens-z Z` on each tile that passes gate P, Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}, masses 2 × 10¹²,
8 × 10¹², 2 × 10¹³ M☉, 100 lenses per cell, seed 50. Lenses are placed inside the tile's screened footprint, so the
tile-edge loss is in the efficiency. **Combined limit:** `w1_volume.py` over the twenty-four headline summaries plus
the unflagged PRIMER tiles, reported beside the twenty-four-field values.

## ASSUMPTIONs
Every threshold above (tile half-width, 0.3 share, R); single-plane lenses as in `results/w1_volume`; flat
tangent-plane tile offsets; R from a cluster fit in another program applies to PRIMER (same band and pipeline).
