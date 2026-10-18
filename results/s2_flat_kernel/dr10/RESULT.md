# S2 with LS DR10 z < 22 galaxies: gate Q passed, gate C (positive control) failed; S2 paused

Plan: `PLAN.md` (frozen before any DR10 galaxy was read). Code: `scripts/s2_flat_kernel.py --sample des --survey dr10
--mag-z-max 22`. Outputs: `../fit_des_dr10_mz22_alpha2.json` (planned), `../fit_des_dr10_mz22_alpha1.json` (planned
α = 1 companion), `../fit_des_dr10_mz21_alpha2.json` (**post hoc**: same DR10 fetch cut at z < 21, to separate the
release change from the depth change). Exploratory: the same DES-SN5YR SNe as the DR9 run (`../fit_des_alpha2.json`).

## Data (observed)
- DR10 Tractor + DR10 photo-z, Data Lab TAP, same selection as DR9 (brick_primary, maskbits = 0, no PSF/DUP), 2′ discs
  around the 1,518 DES-SN5YR SNe: 116,190 unique galaxies at z < 22 (47,828 at z < 21; DR9 z < 21 had 47,228).
  Fetch 7.9 MB in the session scratch data root (not stored in git), 26 TAP queries of ≈ 25–36 s.

## Gate Q (photo-z quality, function `photoz_gate`): passed
| Bin (dereddened z mag) | n with z_spec | σ_NMAD | outlier fraction |
|---|---|---|---|
| < 21 | 6,306 | 0.020 | 0.008 |
| 21–22 | 2,021 | 0.035 | 0.031 |
Thresholds 0.10 / 0.20 (ASSUMPTION). Caveat: spectroscopic objects are a biased, brighter-for-their-bin subset, so the
faint-bin photo-z quality of the full sample is likely worse than shown (not measured).

## Gate C (positive control): failed
| Column | n SNe fit | γ_L (mag per unit X) | |γ_L| / err | scramble p_L | λ_L |
|---|---|---|---|---|---|
| DR9 z < 21 (D-070 addendum) | 1,414 | −0.0183 ± 0.0110 | 1.67 | 0.095 | 0.603 |
| **DR10 z < 22 (planned)** | 1,418 | **−0.0145 ± 0.0129** | **1.13** | **0.252** | 0.649 |
| DR10 z < 21 (post hoc) | 1,393 | −0.0098 ± 0.0115 | 0.85 | 0.387 | 0.614 |
C needs |γ_L| / err ≥ 3 and p_L < 0.01; both fail (gate evaluated from the JSON values above, not by a committed
function: deviation from the owner rule, recorded here). The deeper column does **not** raise the control's S/N:
the error grows (0.0110 → 0.0129) while λ_L moves only 0.60 → 0.65. At fixed depth (z < 21) the release change moves γ_L
by 0.0085 (0.7σ; the fitted SN sets differ by 21), so the control value is not stable against the galaxy catalogue at this level.

## Flat-kernel term (reported, not calibrated by a detected control)
- α = 2, DR10 z < 22: γ_F = −0.0011 ± 0.0031 (p 0.72), trimmed −0.0073 ± 0.0059, colour term on X_F +0.0008 ± 0.0023;
  chain recovery 0.83, λ_F 0.88; computed one-sided 95 % bound 0.0076 mag per unit T/⟨T⟩ (DR9: 0.0193). Because gate C
  failed, the plan does **not** adopt it as the calibrated S2 limit; the D-070 addendum value stays the quoted one.
- α = 1: corr(X_L, X_F) = 0.991, kernel-degenerate as in D-070 (no limit).

## Verdict and next
Per the plan: **S2 paused.** LS-depth photo-z columns (DR9 z < 21, DR10 z < 22) do not detect the ordinary lensing
control on DES-SN5YR. Reopen only with DES Y3 Gold (Shah+2024 detect lensing with it) or a spec-z / deeper column, and
the full DES-SN5YR covariance.
