# `arcs`: public spectral coverage of the open anti rows (metadata only)

Plan: `PLAN.md` (committed before any query). Script: `scripts/arcs_spec_coverage.py`; rule functions
(footprint containment, instrument groups, class, gate Q) unit-tested before the run
(`tests/test_arcs_spec_coverage.py`). Every number below is in `summary.json` (queried 2026-11-08 19:12Z).
Archives: ESO ObsCore TAP (`INTERSECTS` on the server), MAST CAOM (3″ cone, footprint tested here). No data downloaded.

## Gate Q (validation): **passed**
- Positive control, Abell 2744 core: ESO MUSE 13 public cubes; MAST NIRISS WFSS 10, NIRCam grism 13, WFC3/IR 9,
  ACS/WFC 2 (slitless), NIRSpec MSA 21 (all `observed` metadata).
- Negative control, empty sky (RA 90.5, Dec −35.0): no row, class `none`.

## Frozen classes (derived)
| Row | Class | Public covering observations (optical / NIR) |
|---|---|---|
| RXC J0600 581 | **`none`** | none (ESO ALMA cubes only, 12, group `other`; no MAST spectroscopic row within 3″) |
| MACS0308 2157 | **`spectrum_available`** | 1 VLT/MUSE cube: `ADP.2020-01-04T23:56:39.902`, programme 0104.A-0801(A), 2629 s, 475–935 nm, public since 2019-11-17 |
| WHL0137 2194 (explained; no label) | — | none (ESO ALMA cubes only, 4) |

## Reading (narrowest statement)
- 2157's membership can be decided from public data: the MUSE range covers [O II] 3727 Å, Ca H+K and Hα at the
  cluster redshift 0.356, and [O II] / Ca H+K at the stored photo-z 0.756 (wavelength arithmetic, not stored). Whether
  the spectrum has the S/N for a redshift is not known until it is extracted.
- 581 has no public optical / NIR spectroscopic coverage in ESO or MAST; its membership stays undecidable from
  public spectra. Its line rests on model tests only (addenda 8, 13).
- Not a lensing anomaly claim and not a candidate.

## Deviations
- ALMA cubes (ESO ObsCore `cube`) are put in group `other`, not counted as coverage: the plan named optical / NIR IFUs
  and slitless modes; this grouping was written into the tested function before the run, not after.
- The MAST query asked for `dataproduct_type` spectrum and cube (the plan said spectrum); no target row was returned
  either way.

## Limitations
- Footprint containment is an ASSUMPTION for coverage (edge spaxels, masked regions, gaps are not checked).
- Only ESO and MAST were queried (no Keck KOA, Gemini, or survey-specific releases); a `none` is for these two archives.
- MSA footprints would not show whether a slit was placed on the target (none occurred here).
