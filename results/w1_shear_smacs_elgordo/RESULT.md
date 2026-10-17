# W1 shear screen (D-053) on SMACS 0723 and El Gordo: result (2026-10-17)

Plan: `PLAN.md` (committed before the run). Labels: shapes `observed`; ε, S maps, R `derived`; cluster shear
`model_prediction` (SMACS: Mahler+22 ICLv2; El Gordo: Caminha+23); injections `simulated`; limits `derived` from
`simulated`; every threshold is an ASSUMPTION (D-053 defaults, unchanged).

## Inputs (sha256-checked against SOURCES.md and `data/manifests/`)
- F200W pipeline catalogues `jw02736-o001_t001` (SMACS, 3,145 rows) and `jw01176-o241_t012` (El Gordo, 3,108 rows).
- DJA eazy zout: SMACS v7.4 (`82db8910…`), El Gordo v7.0 (`efd051e4…`), each extracted from its tarball (HTTP 200 on
  2026-10-17; tarball sha256 matched SOURCES.md; tarballs deleted).
- Gaia DR3 G < 17 within 4′ of the model centre (`scripts/gaia_stars.py`): 25 stars (SMACS, sha256 `242d41bb…`), 7
  (El Gordo, `cf507cd4…`); `screen_shear_*.json`. Not stored.
- Commands: `exotic_screens.py --model <model> shear --catalog <cat> --photoz <zout> --spike-stars <gaia>
  --max-radius <60|110>`; injections `inject_shear.py --fields <field>` (seed 50, 200 per mass) with the new
  `SHEAR_FIELD_OVERRIDES` entries (Gaia veto and centre; photo-z was already in FIELDS).

## Screen (`screen_shear_smacs0723.json`, `screen_shear_elgordo.json`)
| Field | spike segments | sources used | R | S_max (p_rot) | max S_× / max −S_× (p_rot) | rule |
|---|---|---|---|---|---|---|
| SMACS 0723 | 72 | 1,902 | 0.438 ± 0.038 | 2.91 (0.765) | 2.97 (0.66) / 3.19 (0.25) | not passed |
| El Gordo | 59 | 1,114 | 0.486 ± 0.058 | 2.95 (0.62) | 3.25 (0.15) / 3.01 (0.49) | not passed |

- Known-case check: R lies in the 0.41–0.48 range of the other five fields (El Gordo 0.486 is 0.006 above it,
  within its ± 0.058), so the screen sees each cluster's shear.
- Neither S_max reaches its field's rotation-null median (`null_max_p50_p95`: 3.06 and 3.01); no grid centre is
  above it (`n_peaks_above_null_median` 0). **D-053 rule not passed in either field: null.** Nothing for
  `/vet-candidate`.
- `sheet_smax.png` (F200W, 24″, S3 byte-range read; looked at; not in the plan, post hoc): the SMACS maximum's
  aperture holds a bright star with spikes about 6″ from its centre and spike streaks from a star at the box edge;
  El Gordo's holds one long linear streak in the outer aperture. By-eye reading, not measured.

## Injections (`inject_summary_*.json`; recovered / 200)
| \|M\| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| SMACS 0723 | 0 | 0 | 11 | 62 | 106 |
| El Gordo | 0 | 0 | 0 | 14 | 48 |

Screened area 2.99 arcmin² (SMACS, 60″ grid half-width) and 4.83 arcmin² (El Gordo).

## Seven-field limits (`combined_7fields.json`; 43.7 arcmin²)
- The five earlier fields enter from their committed summaries (`results/w1_shear_s1063/`). They predate two screen
  options (`psf_chi`, `min_coverage`); `inject_shear.py --combine-only` now reads a missing option as the value
  those runs used (`LEGACY_SCREEN_OPTIONS`: isotropic PSF, no coverage cut). `DEFAULT_FIELDS` now holds all seven, so
  `inject_shear.py --combine-only --out <dir>` rebuilds `combined_7fields.json` (checked, equal).
  Check (deviation from the plan, added): Abell 2744 re-run with the current code reproduces its committed
  recoveries, S_max and area exactly (`rerun_abell2744_summary.json`); `min_coverage` 0 had already reproduced all
  five fields (`results/w1_edge_cut`).
- Σ(W1) < 5.3 × 10³ / 9.6 × 10² / 5.0 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (five fields: 5.8 × 10³ /
  1.1 × 10³ / 5.7 × 10²); border-corrected 5.8 × 10³ / 1.05 × 10³ / 5.4 × 10². Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- Same caveats as D-053 (catalogue isophotal moments, no PSF-anisotropy model, one B map per field).
- El Gordo (z_l = 0.87) has the lowest efficiency of the seven fields; its cause is not measured.
- SMACS uses the jwst 2.0.1 catalogue (older pipeline than the other fields; D-010).
