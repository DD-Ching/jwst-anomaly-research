# S1 deformed-twin recovery test: pre-registration

Written 2026-10-11 ~21:20 UTC, before any deformed injection or new-method score was computed. The code reading
and the band-ratio calculation below were done first.

**Question.** Does the S1 burst-twin screen (D-071; `scripts/s1_twins.py`) miss same-source signals that arrive
**deformed**? "Deformed" means time-stretched with the matching energy shift, rescaled in brightness, and seen
through a different response. If it does, which step loses them, and does a minimal change recover more at the
same false-alarm rate?

**Facts established from the code (to be confirmed by tests):**
- Delay: every position-inconsistent pair over 2008–2026 is scored, so multi-year delays are covered.
- Pairs whose stored resolutions differ by more than ×2 (|Δk| > 1) are not scored.
- Frequency axis: there are two bands. The statistic is the multi-band normalised cross-correlation, so the two
  bands' brightness ratio must agree. An energy shift is not modelled. Injections only jitter the band gains by
  15 %.
- Flagging and χ² vetting use only s = 1. The stretch grid is applied, for display only, to the top 300 pairs and
  the flagged ones. A stretched twin can therefore never be flagged through the stretch grid.
- Injections are unstretched copies, so the published efficiency and limit (D-071) apply to s = 1 twins only.
- The pulse-shuffle surrogate null is near-degenerate: 618 / 1,414 eligible bursts have two pulses, whose
  shuffle is a fixed swap, and the 100 surrogate maxima take two values.

**Deformation model** (physics fixes the coupling; ASSUMPTIONs flagged):
- One relative factor D = s ties time and energy, as for a relative Doppler or (1 + z) factor between two paths:
  t → s t and E → E / s. The copy's band fluxes come from shifting a Band spectrum (α = −1, β = −2.3;
  E_peak ∈ {100, 200, 400} keV, ASSUMPTION: typical GBM values) by 1/s.
- Range s ∈ [0.5, 2] (ASSUMPTION: relative Doppler |β| ≤ 0.6 between exits, or a (1 + z) ratio ≤ 2; it equals
  the existing S1 grid range).
- Treated separately:
  - **brightness:** a free ratio from {1, 0.5, 0.3, 0.2, 0.1}, as in D-071;
  - **response:** per-band gain jitter of 0.15, as in D-071;
  - **delay:** the real time between slot bursts.
- Ordinary lensing has s = 1 with no spectral change; it is the s = 1 control. Same-position pairs are excluded
  by the S1 design.

**Methods compared** (identical eligibility; scores at the coarser resolution of each pair):
- **M0:** the existing statistic. ρ(s = 1), total normalisation, |Δk| ≤ 1.
- **M1:** the minimal change. max over s ∈ {0.5, 0.63, 0.79, 1, 1.26, 1.59, 2} of ρ(y stretched by s), |Δk| ≤ 2.
- **M2:** M1 plus a free per-band scale in the statistic, i.e. the band-ratio constraint is dropped.

**False-alarm control.**
- The thresholds come from the real-pair background of a **calibration half** of the eligible bursts (random
  split, seed 1).
- Each method's threshold is set so that a fraction f = 10⁻⁴ of all calibration position-inconsistent eligible
  pairs exceeds it. The maximum over s is part of the statistic, so the multiplicity of deformations is counted.
- χ² vetting (p ≥ 10⁻³, at the best s; per-band scale for M2) and the ≥ 1-day delay are applied identically
  afterwards.
- The real background may contain true twins. That only makes the thresholds conservative.

**Validation.** Injections use only bursts of the other (validation) half. The grid is s0 ∈ {0.5, 0.6, 0.71,
0.85, 1, 1.2, 1.41, 1.7, 2} × E_peak × ratio, with fixed seeds. The real validation-half pairs are scored once
with the fixed methods, as the small hold-out trial.

**Decision rules.**
- **Blind spot confirmed:** M0 recovery at |ln s0| ≥ ln 1.41 is below half of M0's recovery at s0 = 1, for ratio 1.
- **Improvement:**
  - M1 or M2 recovers more deformed injections than M0 at the same f;
  - its recovery at s0 = 1 is not lower than M0's by more than 2σ (binomial);
  - its χ²-passing false pairs in the calibration half stay ≤ M0's + 2.
- Real hold-out pairs passing a method are **candidates for vetting only**, listed separately from injections.
  Recovering injections shows search capability, not evidence for wormholes.
