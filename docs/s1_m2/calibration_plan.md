# S1 M2: consistency-test calibration and the freeze rule

Written 2026-10-11 ~22:05 UTC. It follows the first calibration pass (`results/s1_m2/calibration.json`) and comes
before the refined pass and before any full-catalogue score was computed.

**What the first pass showed** (failed approach, kept on record):
- On simulated true pairs with independent noise, the per-band χ² at the best grid stretch rejects 7.2 % of the
  pairs above threshold (p < 10⁻³). All of them are stretched copies (s0 = 1: 0 %).
- The p distribution is not uniform: KS p = 3 × 10⁻¹⁰⁶. Its median is 0.998 (the 10 % systematic is
  conservative), but its tail is heavy (1st percentile 4 × 10⁻³⁹).
- The rule written into the code, c = min(10⁻³, 1st percentile of the true-pair p), therefore gave c = 4 × 10⁻³⁹.
  That disables the test: 370 / 2,000 top unrelated real pairs pass.
- Cause (hypothesis): grid discretisation. At step 1.05 the best grid point is up to 2.5 % off the true stretch,
  which leaves large residuals for bright copies. Linear interpolation without covariance adds to it.

**Minimal fix.** Refine s locally (ln s within ±0.03 in steps of 0.002, around the best grid value) and re-select
the lag before computing the per-band χ². This changes only the consistency step; the M2 score and its threshold
stay as they are.

**Freeze rule, declared before the refined pass:**
- The consistency statistic is the refined per-band χ² p-value, used as a **calibrated score**, not as an exact
  probability.
- The cut c is the 5th percentile of the refined p over simulated true pairs above the M2 threshold, i.e. 95 % of
  true pairs are kept. The simulation uses seed 777 and all eligible bursts as templates. If that percentile
  exceeds 10⁻³, c = 10⁻³.
- **False-candidate estimate.**
  - Measure the pass rate f_u of unrelated real pairs at c, on the 2,000 most similar pairs of the calibration
    half (M2 range 0.922–0.984). The expected number of false candidates is N_flagged × f_u.
  - Pairs share bursts, so this is an expectation, not a family-wise probability.
  - If it is ≥ 1, the full search is labelled **exploratory**: every candidate needs individual vetting, and no
    global significance is claimed.
- **Frozen afterwards** (written to `results/s1_m2/frozen.json` before the full search):
  - grid (step 1.05, 29 values on [0.5, 2]);
  - M2 threshold (calibration half, per-pair FAR 10⁻⁴);
  - the refinement and c;
  - delay ≥ 1 d, ≥ 2 pulses, position inconsistency at 3σ with the 3.7° systematic, |Δk| ≤ 2;
  - de-duplication by connected components of shared bursts;
  - the efficiency seed.
- **Thresholds are not loosened after the search.**
- **Data already used for tuning:**
  - calibration half: thresholds, unrelated pass rate;
  - hold-out half: looked at with M0/M1/M2 and with the fine grid (#167).
  - The full search therefore re-uses both halves and is not an independent validation.
- **Candidate review rules:**
  - de-duplicate;
  - data quality (bad-bin fraction, valid fraction);
  - re-trigger and same-day checks;
  - known populations (single-envelope similarity, low S/N partner);
  - only then a CSPEC/TTE spectral check of whether E → E/s matches the fitted s.
