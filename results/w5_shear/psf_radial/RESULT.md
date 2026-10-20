# Result: star-ellipticity (PSF) radial pattern vs the Euclid Q1 radial-shear null (D-067)

Plan: `PLAN.md`, committed at 17:11Z (its header's "17:22Z" is a typo) before any star rows were read. Numbers: `summary.json`. Script:
`scripts/w5_euclid_psf_radial.py`. Labels: star shapes observed; ē_r, ρ derived; the signal reference is a
model_prediction; β, the star-count floor and the |e| clip are ASSUMPTIONs.

## Data
- 181,816 stars (VIS 19–22, `point_like_prob > 0.9`) on 343 of 344 D-067 tiles (EDF-F 20,391 / 72; EDF-S 63,839 /
  148; EDF-N 97,586 / 123), from ESA's Euclid archive (EAS TAP), not IRSA (**deviation**: IRSA answered about one
  tile query per minute this run). Identity check on the 12 tiles also fetched from IRSA: 2,031 / 2,031 rows matched,
  max |Δe| 5 × 10⁻¹² (`irsa_identity`), after dropping EAS rows with a NULL `point_like_prob` (or shape column).
- Grid (**deviation**): D-067's step θ_E and annulus 1.5–3 θ_E, but centres laid from the stars' extent (not the
  resolved galaxies'), so they are offset from D-067's screened centres by up to one step. Coverage cut as D-067
  (annulus galaxy count ≥ 0.8 × field median), counted on a server-side 0.5′ map of D-067 galaxies before the
  resolved cut (D-067 counted resolved galaxies).
- Per-centre error: field rms of one component / √n (as `w5.star_psf`). The first committed run used the
  per-centre sample std; /code-review found it unstable at n ≈ 5–9 and it was replaced before this result.

## Frozen outcome
ρ_max(β = 1) 1.19, ρ_max(β = 2.3) 2.73, both at EDF-F θ_E = 1′, where χ²/dof is 0.98 →
**"no PSF pattern above noise where rho is largest (rules 2-3 not read)"** (`decide`, with the plan's χ² > dof
condition). There the maximum |ē_r| 0.122 (8.6σ against the field rms) comes from 41,747 centres with a median of 9
stars, i.e. from single heavy-tailed stars, not from a map-wide pattern. The frozen rule therefore does not
close D-067's PSF revisit condition at 1′.

## By θ_E (same stored numbers; the split per θ_E is post hoc)
| θ_E | max ρ(β = 1) | max ρ(β = 2.3) | χ²/dof (F / S / N) |
|---|---|---|---|
| 1′ | 1.19 (EDF-F) | 2.73 | 0.98 / 0.99 / 1.01 |
| 2′ | 0.42 (EDF-F) | 0.96 | 1.03 / 1.08 / 1.18 |
| 4′ | 0.17 (EDF-F) | 0.40 | 1.20 / 1.44 / 1.66 |

## Post hoc diagnostics (not in the plan)
- At 4′ a PSF radial / tangential pattern rises above noise: |z| > 3 at 25 / 75 / 106 covered centres, against 7 /
  16 / 13 expected for a Gaussian (overlapping annuli are correlated, so the counts are not independent). After
  subtracting each tile's mean star e, it remains at 25 / 77 / 107 (χ²/dof 1.16 / 1.41 / 1.63). So it comes from
  gradients within tiles, not from a constant PSF in partly covered annuli. Its largest amplitude is 0.018 raw
  (ρ(β = 1) 0.17) and 0.017 tile-subtracted.
- Typical (median) 95 % upper bound |ē_r| + 2σ: 0.015–0.027 at 1′, 0.008–0.014 at 2′, 0.004–0.007 at 4′.
- Clipping stars at |e| < 0.2 lowers the 1′ maximum to 0.087 / 0.049 / 0.032 (ρ(β = 1) ≤ 0.85).

## Reading
- PSF anisotropy gradients are measured. At 2′–4′ the largest radial term is ≤ 0.42 (2′) and ≤ 0.17 (4′) × the
  annulus-mean injected radial shear at β = 1. At 1′ the star sample cannot bound it below the signal.
- This test does not change D-067's 2′ / 4′ limits. Their injection efficiency was 1.0, and the worst PSF term is
  below the signal. At β = 2.3 (galaxies at the resolved cut) the 2′ term could reach 0.96 × signal at one centre;
  that bound is not propagated through the injections.
- D-067 "Revisit if (PSF gradients)" is not closed; it narrows to θ_E = 1′ and to galaxies near the resolved cut. A
  per-galaxy anisotropy correction (star-interpolated e_psf × β) before a DR1 run would settle it.

## Failed approach
- First implementation without D-067's coverage cut: partly covered annuli turn a constant PSF ellipticity into a
  radial term (4′ maxima 0.04–0.14, not stored). Any star or galaxy ring statistic needs the
  coverage cut.
- A max-over-centres rule on a heavy-tailed star statistic cannot pass where annuli hold ~10 stars. Future plans
  should use a significance-aware, outlier-robust bound.
- A per-centre sample std of ~5 values as the error (unstable z; replaced by the field rms).
