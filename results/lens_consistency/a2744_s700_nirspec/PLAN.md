# Abell 2744 700.1 (UNCOVER 78) vs public NIRSpec MSA coverage — plan (frozen before the run)

**Question.** 700.1a/b (UNCOVER 78.x1 / 78.x2) is parked "until a local re-fit or spectra of both arc ends"
(D-036 addendum c, CHANGELOG 2026-10-20). Has a NIRSpec MSA spectrum of both ends become public since the
2026-10-20 vetting? NIRSpec MSA programme 8204 in Abell 2744 became public on 2026-11-05.

**Hypothesis (H-new).** A programme released after 2026-10-20 has an MSA footprint on a 78 knot.
**Ordinary readings.** If both ends have grade-3 redshifts, equal z is expected for both "two images" and "one
galaxy" (does not decide); different z would rule out two images of one source.
**What changes the plan.** `both ends measured` → compare z next run; `new input, not in DJA` → wait for a DJA table
newer than v4.6, or read the programme's own products; `no new input` → 700.1 stays parked.

**Inputs.**
- UNCOVER v2.0 image table (`a2744_uncover.UNCOVER_V2["images"]`): the 8 knots 78.11–78.32.
- MAST CAOM via `astroquery.mast.Observations.query_criteria`: `obs_collection` JWST, `instrument_name`
  NIRSPEC/MSA, radius 1″ around each knot (s_region intersection = footprint on the knot, not a shutter).
- DJA NIRSpec emission-line table v4.6 (pinned in `scripts/macs0647_specz.py`), streamed: rows within 0.5″ of a knot.

**Frozen rules (ASSUMPTIONs).**
- `RELEASED_AFTER` = 2026-10-20 (MJD of `t_obs_release`); planned rows (`calib_level` −1) are listed, never inputs.
- Gate K (known case, combined: both must hold): programme 2561 (Price+24, the source of z_spec 2.3151 for system 78)
  is returned on at least one knot, **and** DJA has at least one grade-3 row within 0.5″ of a 78 knot with
  |z − 2.3151| ≤ 0.01 (1 + z).
- An "end" is the set {78.11, 78.21, 78.31} (a) or {78.12, 78.22, 78.32} (b); an end is measured if any of its knots
  has a grade-3 DJA row within 0.5″.
- `decide(gate_k, end_a, end_b, new_public)`: not gate_k → `not validated`; end_a and end_b → `both ends measured`;
  else new_public non-empty → `new input, not in DJA`; else `no new input`.
