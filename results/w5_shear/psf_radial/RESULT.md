# Result: star-ellipticity (PSF) radial pattern vs the Euclid Q1 radial-shear null (D-067)

Plan: `PLAN.md`, committed at 17:11Z (its header's "17:22Z" is a typo) before any star rows were read. Numbers: `summary.json`. Script:
`scripts/w5_euclid_psf_radial.py`. Labels: star shapes observed; ē_r, ρ derived; the signal reference is a
model_prediction; β, the star-count floor and the |e| clip are ASSUMPTIONs.

## Data
- 181,816 stars (VIS 19–22, `point_like_prob > 0.9`) on 343 of 344 D-067 tiles (EDF-F 20,391 / 72; EDF-S 63,839 /
  148; EDF-N 97,586 / 123), from ESA's Euclid archive (EAS TAP), not IRSA (**deviation**: IRSA answered about one
  tile query per minute this run). Identity check on the 12 tiles also fetched from IRSA: 2,031 / 2,031 rows matched,
  max |Δe| 5 × 10⁻¹² (`irsa_identity`), after dropping EAS rows with NULL `point_like_prob`.
- Grid: D-067's (step θ_E, annulus 1.5–3 θ_E) with D-067's coverage cut (annulus galaxy count ≥ 0.8 × field median)
  from a server-side 0.5′ count map of D-067 galaxies before the resolved cut (**deviation**: D-067 counted resolved
  galaxies only).

## Frozen outcome
ρ_max(β = 1) 1.19, ρ_max(β = 2.3) 2.73 → **"PSF could fake or hide a flag"** (rule 3; `decide`). The maximum is EDF-F
at θ_E = 1′: |ē_r| 0.122 at 1.3σ, with a median of 9 stars per annulus. The plan's noise caveat applies: at
1′ the max over 44,016 centres is set by star-shape noise, not a resolved PSF pattern, so the frozen rule has no power
there. D-067 quotes no 1′ limit (efficiency below 0.5) but did screen 1′ for flags.

## By θ_E (same stored numbers; the split per θ_E is post hoc)
| θ_E | max ρ(β = 1) | max ρ(β = 2.3) | max |ē_r| significance | per-θ_E rule reading |
|---|---|---|---|---|
| 1′ | 1.19 (EDF-F) | 2.73 | 1.2–1.4σ | rule 3, noise-limited |
| 2′ | 0.42 (EDF-F) | 0.96 | 1.1–2.3σ | rule 2 |
| 4′ | 0.17 (EDF-F) | 0.40 | 2.7–5.4σ | rule 2 |

## Post hoc diagnostics (not in the plan)
- A real PSF radial / tangential pattern exists at 4′: |ē_r / σ| > 3 at 18 / 103 / 116 covered centres against 7 /
  16 / 13 expected for a Gaussian (EDF-F / S / N; overlapping annuli are correlated, so these counts are not
  independent); χ²/dof 1.31 / 1.59 / 1.76. Its largest amplitude is 0.008–0.018, i.e. ρ(β = 1) ≤ 0.17.
- Typical (median) 95 % upper bound |ē_r| + 2σ: 0.013–0.020 at 1′, 0.007–0.012 at 2′, 0.004–0.007 at 4′
  (ρ(β = 1) ≤ 0.19 everywhere).
- Clipping stars at |e| < 0.2 lowers the 1′ maximum to 0.087 / 0.049 / 0.032 (ρ(β = 1) ≤ 0.85).

## Reading
- PSF anisotropy gradients are measured. Where they can be resolved (2′–4′), the largest radial term is ≤ 0.42 ×
  (2′) and ≤ 0.17 × (4′) the annulus-mean injected radial shear at β = 1. At 1′ the star sample cannot bound it
  below the signal.
- This does not change D-067's 2′ / 4′ limits: their injection efficiency was 1.0, and the worst-centre PSF term
  is below the signal. At β = 2.3 (galaxies at the resolved cut) the 2′ term could reach 0.96 × signal at one
  1.1–2.3σ centre. This bound has not been propagated through the injections.
- D-067 "Revisit if (PSF gradients)": not closed by the frozen rule; narrowed to θ_E = 1′ and to galaxies near the
  resolved cut. A PSF-anisotropy correction (star-interpolated e_psf × β per galaxy) before a DR1 run would settle it.

## Failed approach
- First implementation without D-067's coverage cut: partly covered annuli turn a constant PSF ellipticity into a
  radial term (4′ maxima 0.04–0.14, |z| up to 10.6, not stored). Any star or galaxy ring statistic needs the
  coverage cut.
- A max-over-centres rule on a noisy star statistic cannot pass at scales where annuli hold ~10 stars. Future plans
  should use a significance-aware bound (e.g. the 95 % upper bound at the worst significant centre).
