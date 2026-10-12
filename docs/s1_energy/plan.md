# S1 round 3 plan: does energy information raise full recovery at equal false-candidate burden?

Written 2026-10-12, after the catalogue error checks in section 1 and **before** any A/B recovery,
threshold or validation-half number was computed. Changes after this point are labelled post hoc.

## 0. Question
At the same false-candidate burden, does a time + energy method (B) recover more simulated
same-source pairs than the existing time-shape method (A, frozen M2; `results/s1_m2/frozen.json`)?

## 1. Inputs and error facts already checked (catalogue only, no twin scores)
Source: HEASARC `fermigbrst`, 4,392 rows (`results/s1_energy/gbm_spectral.ecsv.gz`).
- Round 2 used the best-fitting-model flnc E_peak. That was available for 64 % of bursts: 38 % in
  the faintest fluence quartile, 75–78 % in the brightest two, 49 % of short and 67 % of long bursts.
  Missingness depends on brightness and type.
- Model choice matters. For the same burst and data, ln E_peak(Band) − ln E_peak(Comp) has median
  −0.12 and robust sd 0.18. Comparing a Band value with a Comp value would add a bias.
  **Rule:** both members use the same model, the flnc Comptonized fit. It is fitted for 3,773 bursts.
- Integration window and spectral evolution. For the same burst and model, ln E_peak(peak-flux
  1 s) − ln E_peak(fluence interval) has median +0.07 to +0.08 (harder at the peak) and robust sd
  0.33–0.36. The median statistical error is 0.19–0.21, so the excess scatter is about 0.28 and
  the pull sd about 1.85.
  A fainter copy's catalogue interval lies between these two windows. This bounds the window
  systematic at roughly 0 to 0.29 plus a bias up to +0.08. It is not measured more precisely.
- Asymmetric errors: the median upper/lower log-error ratio is 0.87, so the mean of the two log
  errors is used.
- Quality rule: σ(ln E) ≤ 0.5 and 10 keV ≤ E_peak ≤ 10 MeV. Otherwise the value counts as
  **missing**, never as failing.

## 2. Methods
Convention: the program stretches burst b (y) by `s_prog` to match burst a (x).
- The physical model says b is a's source with times t → s_phys·t and energies E → E/s_phys.
- Hence s_prog = 1/s_phys and E_b = E_a / s_phys = E_a · s_prog.
- z = (ln E_a − ln E_b + ln s_prog) / sqrt(σ_a² + σ_b² + 2·σ_sys,test² + σ_lns²).

**A:** M2 > thr_A (0.9716), refined χ² p ≥ 1e-3 and delay ≥ 1 d. This is unchanged.

**B (staged, fixed now):**
- For a pair where both members have a quality Comp E_peak (an "E-testable" pair): M2 > thr_B and
  |z| ≤ 2, plus the same χ² and delay cuts.
- For any other pair, A's rule applies. A missing E_peak never rejects a pair.
- thr_B ≤ thr_A is calibrated so that, among E-testable calibration-half real pairs, the count
  passing (M2 > thr_B and |z| ≤ 2) equals the count passing M2 > thr_A. Equal burden on the
  untestable pairs holds by construction.
- Burden is matched at the flagged level, before χ² (counts are larger there). The candidate-level
  burden (with χ² and delay) is reported as a check.
- σ_sys,test = **0.20**, the middle of the 0–0.29 bracket, rounded. It is fixed now.
- σ_lns,test = the robust sd of ln(s_prog·s_phys) on calibration-half injections with M2 > 0.9,
  rounded up to 0.01.

## 3. Simulation (catalogue level; it does not validate instrument-level gain)
- Templates and noise slots are from the **validation half** (split seed 1). Calibration-half
  injections are used only for σ_lns.
- Existing light-curve injection (`s1_deform.inject_one`) with independent noise, a brightness
  ratio of 1 to 0.1, s_phys from 0.5 to 2 and 15 % gain jitter. The band factors use the
  template's own catalogue Comp E_peak (200 keV if it has none).
- The template's two-band light curves keep its real spectral evolution.
- The copy's E_peak is measured or missing like a real catalogue burst:
  - the copy's fluence is F_a·ratio·s_phys;
  - a donor is drawn at random from the 20 catalogue bursts nearest in (log fluence,
    log E_peak,true);
  - the copy takes the donor's availability (quality rule) and its σ(ln E);
  - its measured value is ln E_a − ln s_phys − ε_a + ε_c + w, with ε ~ N(0, σ);
  - w is the window systematic, resampled from the empirical (peak − fluence) residuals, scaled
    to sd ∈ {0.10, 0.20, 0.29} with bias {0.025, 0.05, 0.075}. These are the optimistic, central
    and pessimistic truth scenarios; the pessimistic one is the measured window bracket edge.
- The test-side σ_sys stays 0.20 in every scenario. The truth model and the test model are
  therefore not the same.
- The stretch estimate is the injection's own M2 grid value, so its errors, including
  wrong-peak tails, are real for this pipeline.

## 4. Deliverables and decision
- Main number: total recovery R = recovered / injected for A and B on the same injections, ΔR =
  R_B − R_A, with a 95 % interval from a cluster bootstrap over template bursts. Also reported:
  - split into E-testable and untestable injections;
  - split by brightness ratio;
  - the validation-half burden for A and B (E-testable val–val real pairs, flagged and candidate
    level), with Poisson intervals.
- **Minimum worthwhile improvement:** ΔR ≥ +0.02 absolute (about +12 % of A's 0.167) **and**
  a 95 % lower bound > 0, in **both** the central and the pessimistic scenarios. The
  validation-half flagged burden of B must also not exceed A's by more than its Poisson
  interval.
- **If it holds:** freeze B and run it once on the existing 747,709-pair screen, the same
  eligible subset. List pairs that B passes and A does not. This is exploratory; the data were
  already seen in #167/#169.
- **If not, or only in the optimistic scenario:** pause the S1 energy upgrade and record the
  minimum input needed to reopen it. No threshold retuning.

## 4a. Amendment before any result (literature check, same day)
Published systematics are larger than the catalogue bracket:
- Collazzi, Schaefer & Moree 2011 (arXiv:1101.1285): analyst choices give about 28 % (1σ) on E_peak
  for the same burst, i.e. about 0.25 in ln E_peak per measurement and about 0.35 for a difference.
- Gruber et al. 2014 (arXiv:1401.5069): background selection gives up to 30 % for weak GBM bursts.
- Goldstein et al. 2012 (arXiv:1201.2981): for 25 % of bursts, Band and Comp E_peak differ by more
  than their combined 1σ errors.

A fourth truth scenario is therefore added, **literature**: sd 0.35, bias 0.075.

The decision rule now requires the gain in the central, pessimistic **and** literature scenarios.
A gain that holds only at a difference systematic ≤ 0.29 counts as depending on an unverified
narrow assumption. σ_sys,test stays 0.20.

## 5. Budget
About 2 h. Catalogue/TAP queries and the existing reduced light curves only; no TTE downloads.

## 6. Post-result notes (2026-10-12, owner review; thresholds and runs unchanged)
- The 0.29 peak/fluence figure in sections 1 and 3 defines a sensitivity scenario only. It is not a proven bound
  on same-interval, same-source measurement error.
- Collazzi et al. quote "28 %, i.e. 0.12 in log10". In ln that is 0.25–0.28. Their value comes from independent
  analyses, so √2 × 0.25 = 0.35 assumes independent errors. A homogeneous catalogue may share part of them. The
  0.35 run is kept as that scenario, not re-run.
- The burden was matched at calibration by flagged count. Validation did not show B adding burden; false-alarm
  equivalence is not established.
