# S1 TTE round 2: do TTE time and energy information improve same-source discrimination?

Written 2026-10-12 01:15 UTC, before any of the checks below was run. Download cap: 1 GB for this round.

**Step 2: are the five new multi-pulse events robust?** (bn240905186, bn101129726, bn210528586, bn081113230,
bn250330911)
- **Variants, 12 per event.** One run uses each of:
  - 4 bin-start offsets (0, ¼, ½, ¾ of the bin) at every resolution;
  - 2 background treatments: [−20, −2] s with [T90 end + 5, + 30] s, linear; and [−40, −5] s with
    [T90 end + 10, + 60] s, constant;
  - 2 detector sets: the bcat mask, and its two brightest NaI by net counts.
  - That is 4 + 2 + 2 = 12 runs, each varying one factor from the nominal.
- **Robust** = ≥ 2 pulses in the nominal run and in ≥ 10 of the 12 variants.
- **Exposure check.** Both background windows and the burst window must lie inside the TTE GTIs, with no gaps
  longer than 1 bin. Otherwise the event is flagged "exposure-limited" and excluded.
- **Negative control.** One smooth pulse (FRED with rise/decay in 1 : 3, width equal to the event's T90) plus
  independent Poisson noise. The peak and background count rates match each pilot burst. Ten realisations per
  pilot burst (300 in all) go through the full multi-resolution counter with all 4 offsets, counting max over
  offsets.
  - The false rate of "≥ 2 pulses" must be ≤ 5 %. If it is higher, raise the prominence threshold until it holds,
    and re-count the pilot and the positive control with that threshold.
- **Positive control.** The 42 eligible short bursts (bcat ≥ 2 pulses), through the same counter. The rate is
  reported with its binomial interval, and it is not used as a general efficiency.
- Reduced light curves (npz) and file versions or hashes are saved under `results/s1_blind/tte/`.

**Step 3: joint time–energy constraint.** Official GBM catalogue spectral fits (fermigbrst; response- and
background-folded by the GBM team), no new download.
- **Model:** one factor s maps t → s t and E → E/s. The catalogue fluence-spectrum E_peak of a true copy is then
  E_peak,A / s.
- **Measurement:**
  - s comes from the time-domain M2 fit (the refined stretch);
  - E_peak comes from the catalogue, with the best-fitting-model value (Band or Comptonized) and its asymmetric
    errors;
  - a systematic of 0.10 in ln E_peak is added per burst (ASSUMPTION: cross-detector or model scatter).
- **Statistic:** z = [ln E_peak,B − ln E_peak,A + ln s] / σ, where σ² = σ_A² + σ_B² + 2 × 0.10² + σ_ln s².
  Pass if |z| ≤ 2.
- **Null (unrelated pairs that look alike in time):**
  - the 64 full-catalogue pairs above the M2 threshold (`results/s1_m2/search.json`);
  - the 2,000 most similar calibration-half pairs, recomputed with their s;
  - only pairs where both members have a catalogue E_peak.
- **True pairs (simulated):** for each null pair's first member, a copy with s0 drawn from the frozen grid. Its
  E_peak = E_peak,A / s0 is "measured" with the catalogue error of a burst of the copy's brightness: σ scales as
  ratio^(−1/2) (ASSUMPTION), capped at 1.0. The copy's s estimate carries σ_ln s = 0.01.
- **Hard-to-soft evolution** does not change the fluence-weighted E_peak scaling under this model. The null test
  uses real bursts, so ordinary evolution is included there.
- **Discrimination** = (true-pair pass rate) / (null pass rate), with binomial intervals.
  - The test helps if the null pass rate is ≤ 0.5 and the true-pair pass rate is ≥ 0.85.
  - With 64–2,000 null pairs, the test cannot support false-alarm rates below about 1 / n_null. That is stated,
    not extrapolated.

**Step 4: candidate check.** Only for events that pass Step 2.
- **Partners:** catalogue bursts at an inconsistent position (3σ, 3.7° systematic) with a T90 ratio in [0.5, 2]
  (the stretch range) and bcat ≥ 2 pulses. Partners need TTE light curves for a common representation, so they are
  limited to the eligible short bursts already downloaded for the positive control.
- **Count reported:** compatible pairs with all bursts (by T90 ratio and position, any data product), so the
  coverage of the eventual TTE search is stated.
- **Scoring:** M2 on TTE two-band light curves at the pair's common resolution, refined s, and per-band χ²
  (c = 10⁻³). Then the E_peak z test.
- **Status:** the sample is too small for a calibrated false-alarm rate, so the check is exploratory. A pair that
  passes goes to vetting (`/vet-candidate`); it is not evidence of anything.
