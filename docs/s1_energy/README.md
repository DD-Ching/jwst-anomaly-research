# S1 round 3: catalogue E_peak in the full search at equal false-candidate burden

**Files:**
- Plan: `plan.md`, written before any A/B number. Section 4a was amended after the literature check and still
  before any result.
- Code: `scripts/s1_energy.py`; tests: `tests/test_s1_energy.py`.
- Results in `results/s1_energy/`:
  - `error_checks.json`;
  - `compare.json`;
  - `gbm_spectral.ecsv.gz`: the HEASARC fermigbrst spectral columns, 4,392 rows; the TAP query, access time and
    VOTable sha256 are in its header.
- No TTE was downloaded. The run used the existing bcat light curves (`results/s1_twins`) and one catalogue query.

**Labels:**
- [obs]: catalogue values and real-pair scores.
- [sim]: injected copies, a catalogue-level simulation only.
- [assump]: error models.

## 1. Answer
**Adding the catalogue E_peak to the existing M2 shape method did not reach the pre-declared practical gain that is
robust across the error scenarios. Under the pre-declared rule, the S1 energy upgrade is paused.**
- **Small positive gains:** in the optimistic and central scenarios, total recovery rises by +1.8 percentage
  points (95 % CI 1.2–2.5) and +1.4 percentage points (0.8–2.1). These gains are real within the simulation.
  They fall short of the pre-declared minimum of +2 percentage points.
- **No gain** in the sensitivity scenario based on the peak/fluence window spread: −0.1 percentage points
  (−0.9 to +0.6).
- **Negative** in the literature-based scenario of two independent analyses: −1.0 percentage points (−2.0 to
  −0.1).
- **Scope of what was tested:** catalogue-level fluence-spectrum E_peak added to the existing shape thresholds,
  in a catalogue-level simulation. Not tested:
  - full time-resolved spectroscopy with the instrument response;
  - other observables;
  - other messengers.
  The result does not rule out wormhole or other exotic-lensing mechanisms. It only bounds what this particular
  catalogue-level upgrade adds to S1.
- No candidate. No evidence of new physics; simulated recovery, tests and CI are not physical evidence.

## 2. Error checks (catalogue only) [obs]
- **Availability:**
  - Round 2's best-model E_peak covered 64 % of bursts.
  - The same-model Comptonized fit with the quality rule (σ ≤ 0.5 in ln, 10 keV–10 MeV) covers 80 %.
  - It depends on brightness: 65 % in the faintest fluence quartile, 87 % in the brightest.
  - It depends on type: 67 % of short and 82 % of long bursts.
  - It depends on the best-fitting model: 63 % for power-law-best bursts.
  - Among S1-eligible (multi-pulse) bursts it is 87–88 %.
  - Missing values are never counted as failures.
- **Model choice:** for the same data, ln E(Band) − ln E(Comp) has median −0.12 and robust sd 0.18. Mixing models
  biases z, so both members of a pair use the Comp fit.
- **Integration window and evolution:** for the same burst and model, ln E(peak 1 s) − ln E(fluence interval) has
  median +0.08 (harder at the peak) and robust sd 0.33. The median statistical error is 0.19, so the excess is about
  0.27 and the pull sd 1.86.
  - **How this number is used:** only to set a sensitivity scenario (0.29). It is **not** a proven lower or upper
    bound on the measurement error of a same-source copy over the same physical interval. The two catalogue
    windows share photons and differ in definition, so they are not two independent measurements of one interval.
- **Literature:**
  - Collazzi, Schaefer & Moree 2011 (arXiv:1101.1285): the 1σ uncertainty from analyst choices is "28 %, i.e.
    0.12 in log10(E_peak)". In ln this is 0.25 (from ln 1.28) to 0.28 (from 0.12 × ln 10); the paper's two
    forms differ slightly.
    - The value is inferred from pairs of published E_peak values for the same burst made by different
      analyses.
    - Our 0.35 scenario assumes two such **independent** analysis errors (√2 × 0.25). Within one homogeneous
      catalogue, part of the analyst-choice error is shared, and the paper itself says identical definitions
      reduce it. The applicable difference systematic could therefore be smaller. This is not verified.
    - The 0.35 run is kept as an independent-analyses literature scenario, not as the established value.
  - Gruber et al. 2014 (arXiv:1401.5069): background selection up to 30 % for weak bursts.
  - Goldstein et al. 2012 (arXiv:1201.2981): Band and Comp E_peak differ beyond 1σ for 25 % of bursts.
  - No published GBM number exists for "same spectrum, different brightness". **The relevant systematic is not
    established by any of these. It is kept as a sensitivity range, 0.10–0.35, not chosen.**
- **Stretch uncertainty** [sim]: on 188 calibration-half injections with M2 > 0.9, ln(s_prog · s_phys) has robust
  sd 0.011 and no case with |Δ| > 0.1. The test uses 0.02. Stretch error is not the limiting factor.
- **Time and energy checks are correlated** [obs]: among unrelated real pairs, the E_peak pass rate rises with
  shape similarity.
  - It is 38–39 % for M2 0.80–0.93, 43 % for 0.96–0.97 and 56 % (10 / 18) above the M2 threshold.
  - Similar-looking bursts tend to have similar spectra. Multiplying separate pass rates would overstate the gain.
  - Round 2's 20 % pass rate used σ_sys = 0.10. With the 0.20 adopted here it is about 40 %.

## 3. Fair comparison [obs + sim]
**Method convention:**
- The program stretches b by s_prog to match a. The physical model is t → s_phys·t and E → E/s_phys, so
  s_prog = 1/s_phys and E_b = E_a·s_prog.

**The two methods:**
- **A** is the frozen M2 method: M2 > 0.9716, χ² p ≥ 1e-3, delay ≥ 1 d.
- **B** is staged, using the same cuts plus energy:
  - for pairs where both members have a quality E_peak: M2 > thr_B and |z| ≤ 2;
  - otherwise it uses A's rule.

**Calibration and validation** (split by burst, seed 1):
- **Calibration half** (185,392 real pairs, 141,554 with E_peak for both):
  - A flags 18 pairs.
  - B at the same flagged count gives thr_B = 0.9664, range 0.9656–0.9684 for 18 ± √18.
- **Validation half** (188,055 real pairs):
  - flagged E-testable pairs: A 6, B 4;
  - candidates after χ² and delay: A 0, B 0.
  - The flagged counts were **matched at calibration**. The validation sample does not show B adding burden.
    False-alarm-rate equivalence is **not established**: the counts are small, so the Poisson intervals are
    wide.
- **Validation injections:**
  - 1,800 copies from 643 validation-half templates, with independent noise;
  - brightness ratio 1–0.1, s_phys 0.5–2;
  - the template's real two-band evolution is kept;
  - the stretch is estimated by the pipeline;
  - copy E_peak availability and errors come from catalogue bursts of similar fluence and E_peak, and 74 % of
    copies are E-testable;
  - a window systematic is resampled from the empirical peak − fluence residuals;
  - the test side always uses σ_sys = 0.20, so the truth model differs from the test model.

| truth systematic (difference) | true copies pass E_peak | R_A | R_B | ΔR (95 % CI, cluster bootstrap) | gained / lost |
|---|---|---|---|---|---|
| optimistic 0.10 | 98 % | 0.181 | 0.199 | +0.018 (0.012–0.025) | 33 / 0 |
| central 0.20 | 95 % | 0.181 | 0.196 | +0.014 (0.008–0.021) | 32 / 6 |
| window-spread sensitivity 0.29 | 89 % | 0.181 | 0.180 | −0.001 (−0.009–0.006) | 29 / 31 |
| literature, independent analyses 0.35 | 85 % | 0.181 | 0.171 | −0.010 (−0.020 to −0.001) | 29 / 47 |

- The gain comes only from E-testable injections (central: 0.198 → 0.217). The 469 untestable ones are
  unchanged by construction.
- Moving thr_B by ±√18 changes ΔR by less than 0.008.
- Differences are in percentage points of total recovery.
- **Why the gain is small:** at equal burden the energy cut removes only about 60 % of shape look-alikes, so the
  M2 threshold drops only from 0.9716 to 0.9664. Most of A's losses are elsewhere (pulse eligibility and the χ²
  test, about 50 % of selected copies), and energy information does not reach them.

## 4. Decision
- **Paused.** No frozen B and no exploratory B screen of the 747,709 pairs. That step was conditional on the
  pre-declared gain and is not run.
- No threshold was changed after seeing these numbers.
- **What would justify reopening (minimum inputs):**
  - A same-source E_peak systematic of about 0.10 or less in the ln difference, demonstrated on real data. For
    example: matched-interval spectral fits of the same burst at artificially reduced S/N, or repeated bursts
    from one source measured by GBM and an independent instrument. In this simulation the 0.10 scenario gave
    +1.8 percentage points. That is a single scenario result, not an upper limit on what energy information can
    add.
  - Or an energy observable that sharply separates look-alikes while keeping true copies, i.e. passing well
    under 60 % of shape look-alikes at the real systematic. Candidates are time-resolved E_peak(t) tracks fitted
    over stretch-matched intervals with the GBM response (GBM Data Tools). That is a new spectral pipeline and is
    out of scope.
  - A larger lever is the shape side: eligibility and the χ² test lose about half of selected true copies, which
    is outside this question.
- **Scope:** catalogue-level E_peak plus the existing shape thresholds, in a catalogue-level simulation. It does
  not test:
  - time-resolved spectra with the instrument response;
  - other observables or messengers;
  - any wormhole or exotic mechanism as such.
- **Owner decision 2026-10-12:** the pause is accepted. #169 is closed out and S1 is paused; other searches
  continue as scheduled.
