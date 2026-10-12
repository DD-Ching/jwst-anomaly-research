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
**At equal false-candidate burden, catalogue E_peak does not recover meaningfully more same-source pairs, so the
S1 energy upgrade is paused.**
- The gain is +1.4 % of injected pairs in the central scenario (95 % CI 0.8–2.1 %). That is below the
  pre-declared minimum of +2 %.
- It vanishes at the measured window systematic (−0.1 %).
- It turns negative at the published systematic (−1.0 %, CI −2.0 to −0.1 %).
- Even the optimistic scenario (+1.8 %) misses the minimum.
- No candidate. There is no evidence of new physics.

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
  - A fainter copy's catalogue interval lies between these two windows. This brackets the copy-vs-original
    systematic at about 0–0.29 in the difference.
- **Literature:**
  - Collazzi, Schaefer & Moree 2011 (arXiv:1101.1285): about 28 % per measurement from analysis choices, i.e.
    about 0.35 in the difference.
  - Gruber et al. 2014 (arXiv:1401.5069): background selection up to 30 % for weak bursts.
  - Goldstein et al. 2012 (arXiv:1201.2981): Band and Comp E_peak differ beyond 1σ for 25 % of bursts.
  - No published GBM number exists for "same spectrum, different brightness". **The systematic is therefore not
    constrained below about 0.29–0.35. It is kept as a sensitivity range, not chosen.**
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
  - candidates after χ² and delay: A 0, B 0. The burden is matched.
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
| measured bracket 0.29 | 89 % | 0.181 | 0.180 | −0.001 (−0.009–0.006) | 29 / 31 |
| literature 0.35 | 85 % | 0.181 | 0.171 | −0.010 (−0.020 to −0.001) | 29 / 47 |

- The gain comes only from E-testable injections (central: 0.198 → 0.217). The 469 untestable ones are
  unchanged by construction.
- Moving thr_B by ±√18 changes ΔR by less than 0.008.
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
    from one source measured by GBM and an independent instrument. Even then, the catalogue-level gain is at
    most about +2 %.
  - Or an energy observable that sharply separates look-alikes while keeping true copies, i.e. passing well
    under 60 % of shape look-alikes at the real systematic. Candidates are time-resolved E_peak(t) tracks fitted
    over stretch-matched intervals with the GBM response (GBM Data Tools). That is a new spectral pipeline and is
    out of scope.
  - A larger lever is the shape side: eligibility and the χ² test lose about half of selected true copies, which
    is outside this question.
- **Scope:** catalogue-level simulation only. It does not test instrument-level spectral fitting.
