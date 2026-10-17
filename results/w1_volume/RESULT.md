# W1 shear limits as a volume density over lens redshift: result (2026-10-17)

Plan: `PLAN.md` (frozen before any run; addendum A after the first gate attempt, before the grid). Code:
`scripts/inject_shear.py --lens-z` (new option, default path unchanged), `scripts/w1_volume.py`. Labels: injected
lenses `simulated`; efficiencies and limits `derived`; the comparison value is a published limit; every threshold an
ASSUMPTION.

## Gates
- **Gate Z, attempt 1: failed as written** (`gate_attempt1.json`): a three-mass run cannot reproduce the stored
  five-mass run, because seeds are spawned per `--mass` entry (plan error, addendum A).
- **Gate Z, five masses: passed** (`gate_z.json`): Abell 2744 at `--lens-z 0.3072` gives 0 / 0 / 14 / 72 / 146 of
  200, equal to `results/w1_streak_veto/inject_summary_abell2744.json`.
- **Inputs** (`inputs_check.json`): on all eleven fields the regenerated streak masks (pixel counts), the vetoed
  rows, the real-field S_max and the screened area equal the stored `results/w1_streak_veto` values. The real
  screens are therefore unchanged (null in every field; their sheets were looked at in that run). No new sheet.

## Grid (eleven fields, 66.3 arcmin²; 100 lenses per field, mass and Z; `summaries/z<Z>_<field>.json`)
Recovered lenses summed over the eleven fields, out of 1,100 per cell (`limits.json`, `n_recovered_by_z`):

| \|M\| (M☉) | Z = 0.1 | 0.3 | 0.6 | 1.0 | 1.5 |
|---|---|---|---|---|---|
| 2 × 10¹² | 269 | 72 | 13 | 2 | 0 |
| 8 × 10¹² | 598 | 314 | 147 | 40 | 9 |
| 2 × 10¹³ | 580 | 529 | 304 | 141 | 25 |

θ_E(z_s = 2) at 2 × 10¹³ M☉ falls from 19.5″ (Z 0.1) to 3.8″ (Z 1.5) (`theta_e_zs2_arcsec_by_z`).

## Volume density limit (95 %, zero detections; Planck18; `limits.json`)

| \|M\| (M☉) | V_eff (Mpc³) | n₉₅ (h³ Mpc⁻³) | / Takahashi & Asada 10⁻⁴ |
|---|---|---|---|
| 2 × 10¹² | 1.48 × 10³ | 6.5 × 10⁻³ | 65 |
| 8 × 10¹² | 1.13 × 10⁴ | 8.5 × 10⁻⁴ | 8.5 |
| 2 × 10¹³ | 2.65 × 10⁴ | 3.65 × 10⁻⁴ | 3.6 |

**Reading.** The eleven JWST cluster fields do not reach Takahashi & Asada's (2013) SDSS-quasar bound at any mass
tested: 3.6× weaker at 2 × 10¹³ M☉, 8.5× at 8 × 10¹², 65× at 2 × 10¹². This replaces the hand-spread comparison in
docs/exotic_limits.md (~3.2× at the cluster redshift only). By the plan's own criterion (weaker by ≥ 10× closes the
line), the line is not closed at 2 × 10¹³ and 8 × 10¹² M☉, but no mass reaches the bound either: at similar efficiency, n₉₅ scales as 1 / area, so reaching it at
2 × 10¹³ M☉ needs about 3.6× the area (ASSUMPTION; not run).

## Limitations
- Single-plane ASSUMPTION (PLAN.md): the cluster's lensing of the W1 images (Z < z_cl) or of the W1 position
  (Z > z_cl) is neglected. Not measured.
- ε is zero outside 0.1 ≤ Z ≤ 1.5 and linear between five grid points (conservative outside the grid; the
  interpolation error between grid points is not measured).
- 100 lenses per cell: binomial errors per field and cell are in the summaries (`efficiency_err`); not propagated.
- Deviations (addendum A): RXJ2129 left out (MAST catalogue HTTP 404, 2026-10-17); 100 instead of 200 lenses.
- The comparison treats Takahashi & Asada's |M| > 10¹² M☉ bound as applying at each mass (ASSUMPTION).
