# W1 faint-streak veto for the CANUCS flanking fields: result (2026-10-18)

Plan: `PLAN.md` (frozen before any mask). Code: `scripts/w1_streak_mask.py --mrt-threshold` (default 5 unchanged).
Labels: pixels `observed`; masks `derived`; every threshold an ASSUMPTION.

## K1n: failed for both variants, so the route stops
| Variant | Field | status-2 trails | all candidates | mask pixels | covers S_max (1″) |
|---|---|---|---|---|---|
| D-097 (threshold 5) | macs0416-ncf | 2 | 67 | 6,393 | no |
| A (threshold 3) | macs0416-ncf | 2 | 86 | 6,393 | no |
| B (threshold 2.5) | macs0416-ncf | 2 | 87 | 6,393 | no |
| D-097 (threshold 5) | macs0417-ncf | 2 | 53 | 20,242 | no |
| A (threshold 3) | macs0417-ncf | 2 | 90 | 20,242 | no |
| B (threshold 2.5) | macs0417-ncf | 2 | 105 | 20,242 | no |

(`masks_t3/`, `masks_t2.5/`: JSON and mask per field; sheets in `masks_t2.5/`.) A lower threshold adds candidates,
but TrailFinder's persistence check keeps the same two status-2 trails, so the masks are identical pixel for pixel.
The threshold is not what hides these features. Per the plan: no injection run, the headline stays at fourteen
fields (`results/w1_shear_ncf/RESULT.md`: n₉₅ < 2.35 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉, 2.35× Takahashi & Asada).
**Deviation:** the plan said all five NCF; only the two flagged NCF were masked, because K1n is decided by them alone
and no variant passed.

## What the sheets show (looked at, `masks_t2.5/*.png`; post hoc reading, not a test)
- macs0416-ncf: the S_max centre (block-4 pixel ≈ 1518, 619) sits beside the bright star at the centre of the
  right-hand detector block; the "diagonal streak" through the aperture lines up with that star's long diffraction
  spikes. That is the V2 (spike) case of `w1_peak_vet.py`, which already explains this flag under the frozen V1/V2
  rule; a trail finder is the wrong tool for it.
- macs0417-ncf: the centre (≈ 1914, 713) is at the left edge of the upper-right detector, which shows many parallel
  horizontal bands (1/f-type row striping) over the whole detector, plus a chip gap. These are not single long trails,
  so a median Radon transform with a persistence check is not built to find them.

## Failed approach and next step
Lowering the findsat_mrt threshold (D-097 "Revisit if") does not catch spike streaks or detector row striping.
A veto for these needs a different detector: a field-wide row-segment median test (V1r of `w1_peak_vet.py` made
field-wide on fixed-length row segments, frozen and run through the injections), and a check of why the
Gaia spike veto (D-043) leaves the macs0416-ncf star's spike sources in (not checked here). Either needs its own frozen plan.
