# Abell 2744 700.1 (UNCOVER 78) vs public NIRSpec MSA coverage — result

Plan: `PLAN.md` (frozen and committed before the run). Numbers: `summary.json` (run), `side_checks.json` (post hoc).
Script: `scripts/a2744_s700_nirspec.py`; `decide` unit-tested in `tests/test_a2744_s700_nirspec.py` before the run.

**Outcome (frozen rule): `no new input`.** 700.1 stays parked (D-036 addendum c).

- Gate K passes: programme 2561 is returned on every knot, and DJA v4.6 (99,478 rows read) has a grade-3 row at
  z 2.3108 within 0.5″ of a 78 knot (|Δz| within 0.01 (1 + z) of 2.3151). (observed)
- End a (78.11 / 78.21 / 78.31) and end b (78.12 / 78.22 / 78.32): no DJA row within 0.5″, any grade. (observed)
- Programmes with an MSA footprint on a knot: 2561 (released 2023-08-08) and the unexecuted 6793 and 9645
  (`calib_level` −1, no release date). None was released after 2026-10-20. (observed, archive metadata)

**What the data say beyond the rule.**
- The system 78 z_spec rests on one spectrum, UNCOVER 60018 (PRISM, grade 3), which matches 78.13 at 0.012″ and
  78.23 at 0.112″: the **c** knots. No public spectrum measures the a/b ends, so the z used for 700.1a/b in D-036 is
  assigned through the system identification, not measured at a or b. (derived) This does not change the earlier
  readings: Bergamini+23 reproduces 700.1 at no z in 1.0–4.0 (D-036 addendum 2026-10-20, `a2744_s700_zscan/summary.json`), so the z assignment is not what
  keeps it open.
- Footprint ≠ shutter: being inside a 2561 / 6793 / 9645 footprint does not mean a shutter was or will be placed on
  a or b. (assumption stated)

**Post hoc side checks** (`side_checks.json`, deviation: added after the plan, not in its rules):
- Post hoc: programme 8204 (NIRSpec MSA, 803 rows, released 2026-11-05) has pointing centres 1.52′ from the a/b knot mean;
  700.1 lies outside every 8204 `s_region` polygon (nearest vertex 8.8″), as the per-knot footprints in `summary.json`
  also show. So the new release does not cover 700.1. Within 6″ of the a/b mean only 2561, 6793 and 9645 have MSA
  footprints. (Correction before review: a first computation at a mistyped position gave 0.17′ / 8.4″; superseded.)
- Post hoc, DJA: no emission-line table newer than v4.6 found by name (v4.7, v4.8, v5.0, v5.1, v4.6.1: 403; the bucket listing is
  also 403, so this is absence by name only). v4.6 is unchanged (Last-Modified 2026-09-15). The recheck conditions of
  MACS1149 16.2 and El Gordo 6c / 13c / 20c ("a DJA table newer than v4.6") are not met.
- Post hoc, El Gordo: MAST has 0 JWST NIRSpec rows (MSA / IFU / fixed slit, any status) within 0.06° of the cluster. The
  same query returns NIRCam 1176 (8 rows) and 6882 (4 rows). The "El Gordo NIRSpec programme" recheck condition is
  not met.

**Failed approach.** CAOM `s_ra` / `s_dec` of NIRSpec MSA rows are the pointing centre (every 9645 row
returned had the same position; not stored), not the source, so they cannot test shutter placement. Use footprints (`s_region`) for
coverage and DJA / `_msa.fits` source tables for shutters.

**Next.** Recheck when 6793 or 9645 is released (MAST `t_obs_release`), or when a DJA table newer than v4.6 appears.
A local re-fit of the northern group (O4 / O9 / 20036) remains the other route (owner machine).
