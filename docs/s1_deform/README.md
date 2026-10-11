# S1: does the burst-twin screen miss deformed same-source signals?

Owner brief (2026-10-11). Two signals may share a source but reach us by different paths. They can then differ by
a delay of years, a time stretch with a matching frequency shift, magnification or dimming, and different
instrument responses. Does S1 (D-071) find such pairs, and if not, where does it lose them?

**Pre-registration:** `docs/s1_deform/preregistration.md`, written before any deformed injection or new score was
computed.
**Code:** `scripts/s1_deform.py`. **Results:** `results/s1_deform/summary.json` (pre-registered run) and
`results/s1_deform/summary_posthoc_grid105_holdout.json` (post-hoc grid). **Tests:** `tests/test_s1_deform.py`.

**Labels:** [obs] real GBM data; [sim] injected copies; [derived] computed here; ASSUMPTION marks a stated
assumption. A recovered injection measures search capability only; it is not evidence of a wormhole.

## 1. What S1 and E1 actually covered (checked in code; confirmed by tests and the injections below)

| Question | Answer | Where |
|---|---|---|
| Multi-year delays? | **Yes in S1.** Every position-inconsistent eligible pair over 2008–2026 is compared. **No in E1:** lag bins stop at 7 d, and its null permutes times within a calendar year. | `scripts/s1_twins.py` `score_all`; `src/jwst_anomaly/event_network.py` `lag_edges` |
| Frequency shift? | **No.** S1 has two bands (50–300 keV; 10–50 + 300–1000 keV). The statistic is the multi-band normalised cross-correlation with one overall scale, so the band ratio must stay the same. Injections only jitter band gains by 15 %. A shift E → E/s changes the band ratio by up to 40 % for s = 0.5–2 (E_peak = 100 keV; `band_factors`). | `burst_twins.xcorr_max`, `inject_twin` |
| Stretch tested only after an undeformed pre-screen? | **Worse than that.** Flagging and χ² use s = 1 only. The stretch grid is computed for display on the top 300 pairs and the flagged ones, and never flags anything. | `s1_twins.main` (flag = ρ(s = 1) > ρ*) |
| Which signals do the injections and the null cover? | **Undeformed copies only** (s = 1, flux ratio 0.1–1, 15 % gain jitter). The D-071 limit (< 7.2 × 10⁻³ twin pairs per eligible burst at ratio 1) applies to s = 1 twins. | `injection_one` |
| Other cuts | Pairs whose stored resolutions differ by more than ×2 are not scored. That is 51 % of 998,991 eligible pairs (`summary.json`: 510,294 scored). | `allowed` in `s1_twins.main` |
| Null validity | **The pulse-shuffle null is close to degenerate.** 618 / 1,414 eligible bursts have exactly two pulses, and shuffling them is a fixed swap. The 100 surrogate maxima take only two values: 0.9569 (×82) and 0.9636 (×18). ρ* = 0.9636 is therefore not a well-estimated family-wise threshold. | `results/s1_twins/null.json` |

The external comparison paper, Ahlgren & Larsson (arXiv:2006.07095, ApJ), searches 11 years of GBM data for
lensed GRBs. It allows delay and flux scaling only, with no stretch and no energy shift, and pre-selects on T90,
time-averaged spectra and position. It gives no efficiency figure and states that null results from cross-correlation
searches give unreliable constraints. S1 shares its "undeformed" assumption, so this check is not a duplicate of it.

## 2. Deformation model (the physics ties the parameters)

- **Time and energy.** One relative factor s acts on both, t → s t and E → E/s, as a relative Doppler or (1 + z)
  factor between two paths does.
- **Spectrum.** The copy's band fluxes come from shifting a Band spectrum (α = −1, β = −2.3; E_peak = 100, 200,
  400 keV; ASSUMPTION, typical GBM values) by 1/s.
- **Range** s0 ∈ {0.5, 0.6, 0.71, 0.85, 1, 1.2, 1.41, 1.7, 2}. ASSUMPTION: a relative exit speed |β| ≤ 0.6, or a
  (1 + z) ratio ≤ 2.
- **Handled separately:**
  - brightness: ratio {1, 0.5, 0.3, 0.2, 0.1};
  - response: per-band gain jitter 0.15 (as D-071);
  - delay: the real time between the slot bursts;
  - noise: the slot burst's background plus Poisson-like source variance.
- **Ordinary lensing** is the s0 = 1 case with no spectral change; it is the control. S1 excludes same-position
  pairs by design.

## 3. Methods and false-alarm control

- **The three methods:**
  - **M0:** the existing statistic, ρ(s = 1) with |Δk| ≤ 1.
  - **M1:** max over an s grid of ρ(stretched), with |Δk| ≤ 2.
  - **M2:** M1 with a free per-band scale. The band-ratio constraint is dropped, and χ² is computed per band.
- **Equal false-alarm rate.**
  - Each method's threshold is set on the real-pair background of a **calibration half** of the 1,414 eligible
    bursts (239,690 pairs). Exactly 10⁻⁴ of those pairs exceed it, i.e. 23 pairs for every method.
  - Taking the maximum over s is part of the statistic, so the multiplicity of deformations is counted.
  - After the χ² and delay cuts, **0** of those 23 pass for any method.
- **Validation.** Injections and the hold-out real pairs use only the other half of the bursts.

## 4. Results [sim on real light curves; derived]

**Pre-registered run** (grid 0.5, 0.63, 0.79, 1, 1.26, 1.59, 2; 1,620 injections; recovered = threshold + χ² +
≥ 1 d):

| Injected s0 | n | M0 | M1 | M2 |
|---|---|---|---|---|
| 1 (undeformed control) | 180 | **29** | 26 | 26 |
| all deformed (s0 ≠ 1) | 1,440 | **0** | 74 | 117 |
| on or near the grid (0.5; 2; 1.2) | 540 | 0 | 68 | 98 |
| off the grid (0.6; 0.71; 0.85; 1.41; 1.7) | 900 | 0 | 6 | 19 |

- **Where M0 loses them** (ratio 1, eligible copies): at the ρ(s = 1) threshold.
  - The median ρ(s = 1) is 0.70–0.90 for deformed copies, against 0.96 at s0 = 1.
  - The |Δk| cut loses only 1 copy, and χ² loses a few.
  - Selection (≥ 2 pulses for both members and an emulated window) passes 56–78 % of copies at ratio 1, but only
    26–48 % over all brightness ratios. The JSON field `eligible` is the **pass** fraction. Selection is the same
    for all three methods, so the differences between methods arise after it (§8 gives selection rate, recovery
    after selection and total recovery separately).
- **Pre-registered decision rules:**
  - "Blind spot confirmed": M0 recovers 0 / 36 at |ln s0| ≥ ln 1.41 versus 12 / 36 at s0 = 1 (ratio 1). **Met.**
  - "Improvement": M1 and M2 recover more deformed copies at the same rate (74 and 117 vs 0); s0 = 1 is 26 vs 29,
    within 2σ; false pairs after χ² are 0 vs 0. **Met.**
  - The coarse grid misses off-grid stretches. Recovery needs |ln(s·s0)| ≲ 0.05.

**Post hoc** (chosen after seeing the pre-registered result; fresh injection seed): a log-uniform grid with step
1.05 (29 values).

| Injected s0 | n | M0 | M1 | M2 |
|---|---|---|---|---|
| 1 (control) | 180 | 39 | 36 | 37 |
| all deformed | 1,440 | **0** | 185 | **259** |
| per s0 (all ratios), M2 | 180 each | — | — | 26, 18, 31, 29, 32, 37, 38, 48 (0.5, 0.6, 0.71, 0.85, 1.2, 1.41, 1.7, 2) |

- At the same false-alarm rate, the thresholds are M0 0.948, M1 0.964 and M2 0.970.
- M2 with the fine grid recovers deformed copies at **18 %**, against **21 %** for undeformed ones. M0 recovers 0 %.
- Brightness ratio 1, M2: 11–19 of 36 per s0, against 16 / 36 at s0 = 1.
- Because this grid was chosen after the pre-registered result, the gain is a method result on independent
  injections. It is not a pre-registered claim.

**Hold-out trial on real data** [obs] (validation half, 240,104 pairs; methods fixed beforehand):
- Pairs above threshold: 16 (M0), 11 (M1), 10 (M2) with the pre-registered grid; 16, 11, 12 with the fine grid.
- **None passes χ² (p ≥ 10⁻³) and the ≥ 1-day delay. No new candidates.**
- The calibration half was used only for thresholds; its real pairs were not searched for candidates.

**What the frequency axis can and cannot show.**
- With two coarse bands, M2 only **tolerates** a band-ratio change. It cannot **test** that the energy shift
  matches the stretch (E → E/s).
- Testing that link needs time-resolved spectra: GBM CSPEC (128 channels) or TTE for the two members of a
  candidate pair. Fitting E_peak(t) of both, a true deformed twin predicts E_peak,2(s t) = E_peak,1(t)/s.
- This is cheap **per candidate** (two TTE files). It does not need the whole archive, and it is not needed until
  a candidate exists.

## 5. Conclusions

- **Search capability (method result, not physics):**
  - **Code-level finding** (holds generally): S1's candidate selection never used the stretched similarity. Flagging
    and χ² are s = 1 only, and the stretch grid is display-only.
  - **In this injection set:** the existing chain recovered **0 / 1,440** deformed copies (s0 ≠ 1, under this
    stretch, brightness, spectral, noise and selection model). The D-071 limit applies only to undeformed twins.
  - With a stretch grid in the flagging statistic, a free per-band scale and a 5 % grid step, about 18 % of the
    deformed injections were recovered (21 % undeformed) at the same per-pair false-alarm rate. **That 18 % holds
    only for these assumptions:**
    - s0 ∈ [0.5, 2];
    - a Band spectrum with α = −1, β = −2.3 and E_peak ∈ {100, 200, 400} keV;
    - brightness ratios 1–0.1;
    - 15 % gain jitter;
    - the slot-burst noise;
    - S1's selection.
    It is not a general efficiency for every same-source signal.
- **Real data:** no candidate in the hold-out half.
- **Not done:**
  - the calibration half's real pairs were not searched;
  - the full catalogue was not re-screened (the owner asked not to);
  - no rate limit for deformed twins was derived.
- **Separate issue:** the S1 pulse-shuffle null is close to degenerate (§1). D-071's ρ* should be re-derived
  with a null that randomises two-pulse bursts, for example time-reversal or pulse-amplitude permutation.

## 6. Owner's hypothesis: microscopic links, short physics check

The claims: links may be common; nearby mouths show little effect; external fields may change the links; and we
may see the other side's scene rather than recognising "light in, light out".

**Closest literature** (read in full by a subagent; key points re-checked):
- **Savelova & Kirillov, arXiv:1405.3171** (gr-qc; no journal found).
  - Spacetime foam is modelled as a gas of **virtual (Euclidean) Planck-size wormholes**. These are tunnelling
    configurations, not solutions that carry signals.
  - They compute the Euclidean scalar Green function (G = N(k)/k², Eqs. 10, 17) and its shift in an external
    current (Eq. 22).
  - Within their model, a dilute gas of Euclidean virtual wormholes coupled to a massless scalar through the bias
    kernel σ² of Eq. 22 (left unspecified), with an external current J, they find δn < 0 (Eq. 23): the external
    field **reduces** the virtual-wormhole density. That result holds for that model and that coupling only. It is
    not a general statement about fields and wormholes, but it does not support "a field creates a link".
  - They give no falsifiable formula; their remarks on observables (a shift in Λ, a speed of light slightly above
    c in a field) are qualitative.
- **James, von Tunzelmann, Franklin & Thorne, arXiv:1502.03809** (Am. J. Phys. 83, 486).
  - A **conditional optical calculation for a specified geometry**: null geodesics through a given classical
    traversable wormhole metric (Ellis / "Dneg" with chosen ρ, a and W). Within that geometry, images inside the mouth
    show the **other side's sky**, and an Einstein ring appears.
  - It assumes the metric exists and says nothing about how such a geometry is sourced or how common it is.
- **Already in this repository:** macroscopic-mouth signatures (W2 demagnification, W4 gutters;
  `docs/exotic_lensing.md`) and the Ellis/W3 limits.

**Which part can become a computable, falsifiable prediction?** Only the "we see the other side" claim, applied to
a **population** of mouths.
- Suppose mouths with cross-section σ(λ) and number density n reroute photons to an exit drawn from an exit
  kernel K(x_out | x_in).
- Rerouting is **not** pure absorption. Photons leave a line of sight (mouths on it send them elsewhere) and also
  **enter** it (photons from the other side's sources exit into it). The direction and frequency maps of the exit
  kernel decide the net effect.
- A treatment therefore has to balance removal against injection along each sightline, including the other
  side's radiation field and the angular and frequency mapping. Applying e^(−nσL) alone is wrong. In the limit
  of an isotropic, statistically identical other side, the losses and gains can cancel for diffuse light while
  imaging is still blurred.
- **Possibly falsifiable against existing data** (not computed): image blurring or decorrelation of distant
  compact sources, distance-duality and cosmic-opacity tests, and diffuse-background limits. Each needs the full
  transfer equation, not an absorption factor.

**What is missing for a definition.** Three quantities, without which nothing is testable: σ(λ), the exit kernel
K, and n.
- "Links everywhere but invisible when the mouths are close" corresponds to K being concentrated near x_in. Then
  rerouting looks like a short-range diffusion of photon direction: image blurring that grows with distance.
  That too is computable, but only once K's width is defined.
- "External fields change the link" needs a stated dependence n(F) or σ(F). The only computed case (1405.3171)
  has n decreasing.
- "It is all waves" and "the displacements cancel" are not mechanisms. They specify none of σ, K or n.
- Virtual (Euclidean) wormholes do not carry signals. A population of traversable mouths needs a source of
  exotic stress-energy (rounds 4–7).

**One next step** (not started, owner: on hold): write σ(λ) and K for the Ellis mouth from 1502.03809's ray map,
then the radiative-transfer balance (removal plus injection) along a sightline. Only after that can it be compared
with blurring or opacity data.

## 7. Next actions (owner decides)

1. Re-derive D-071's ρ* with a non-degenerate null, then re-screen the **full** catalogue with M2 and the 1.05
   grid. About 15 min of CPU, no download.
   - Stop if no pair passes the chain at a properly estimated family-wise rate.
   - Any pair that does pass goes to `/vet-candidate` and a TTE spectral check (§4).
2. A deformed-twin rate limit from the full screen, with these efficiencies. Cheap once step 1 is done.
3. The owner's mouth-population idea (§6) is on hold. It needs a transfer-equation treatment, not an absorption
   factor.

## 8. Calibration, frozen method and the full-catalogue M2 re-screen (owner review, 2026-10-11)

**Files:**
- `docs/s1_m2/calibration_plan.md`: the freeze rule, written before the refined pass and the search.
- `scripts/s1_m2_search.py`, `tests/test_s1_m2_search.py`.
- Results: `results/s1_m2/` (`calibration_pass1_unrefined.json`, `calibration_refined.json`, `frozen.json`,
  `search.json`, `efficiency.json`).

**8.1 Consistency test, calibrated on simulated true pairs** (independent noise; seed 777; all eligible bursts as
templates):
- **Failure found.** The per-band χ² at the best grid stretch rejected **7.2 %** of true pairs above threshold
  (p < 10⁻³), all of them stretched copies (s0 = 1: 0 %).
  - Cause: grid discretisation (up to 2.5 % off the true stretch at step 1.05) on bright copies. Linear
    interpolation without covariance adds to it.
  - The cut rule then written in the code (1st percentile of the true p) degenerated to c = 4 × 10⁻³⁹, which
    disables the test: 370 / 2,000 top unrelated pairs passed.
- **Minimal fix.** Refine s locally (ln s ± 0.03 in steps of 0.002) and re-select the lag before χ².
  - True-pair rejection at 10⁻³ falls to **1.4 %**: 0.8 % at brightness ratio 1, rising to 9 % at 0.1.
  - The p distribution is still not uniform (median about 1, because the 10 % systematic is conservative). It is
    used as a calibrated score, not as an exact probability.
  - By the declared rule, c = 10⁻³ (the 5th percentile of the true-pair p is 0.37).

**8.2 Does the old surrogate null represent unrelated bursts?** No (`results/s1_twins/null.json`).
- Above ρ = 0.8, real pairs (almost all unrelated) outnumber the mean pulse-shuffled catalogue 7.5–14×. At
  ρ ≥ 0.95 there are 56 real pairs against 4.8 surrogate pairs.
- So the shuffle is less similar than real unrelated bursts. The old ρ* = 0.9636 therefore **underestimated** the
  background false-alarm rate; its two-valued maxima were a symptom of the same problem.
- This study's background is the real-pair distribution of the calibration half. It keeps real shapes, noise and
  selection, and any true twins in it would only make it conservative.

**8.3 False-alarm rates, kept separate:**
- **Per-pair pre-screen:** 10⁻⁴ by construction (M2 threshold 0.9716 on 185,392 calibration-half pairs).
- **False candidates over the full catalogue** (after stretch, lag, per-band scale and consistency):
  - unrelated real pairs pass the consistency and delay cuts at **2.4 %** (48 of the 2,000 most similar
    calibration-half pairs; 2 of the top 100);
  - with about 75 pairs above threshold, **about 1.8 false candidates are expected**.
  - Pairs share bursts (the top 2,000 involve only 291 distinct bursts), so this is an expectation, not a
    family-wise probability. Extreme tails were not estimated.
  - **The search is exploratory. No global significance is claimed.**

**8.4 Frozen method** (`results/s1_m2/frozen.json`, frozen 2026-10-11T22:06:22Z and committed in 4b4d235 before
any full-catalogue score was computed):
- data hash `0403d265d7c3d47f` (results/s1_twins at main eed9b36);
- grid step 1.05 on [0.5, 2] (29 values) and the M2 threshold 0.9716;
- the refinement, c = 10⁻³, delay ≥ 1 d, ≥ 2 pulses, position inconsistency at 3σ with the 3.7° systematic,
  |Δk| ≤ 2;
- de-duplication by shared bursts; calibration seed 777; efficiency seed 4242.

**Data already used for tuning:**
- **Calibration half:** the threshold and the unrelated pass rate.
- **Validation half:** seen with M0/M1/M2 and again with the fine grid. It is **not** an untouched validation
  sample; a fresh noise seed is not fresh data.
- **Simulated true pairs:** the consistency cut.

**8.5 Full-catalogue re-screen** [obs; `search.json`; 414 s on 4 cores]:

| | count |
|---|---|
| catalogue bursts / reduced | 4,390 / 4,389 |
| eligible bursts (≥ 2 pulses) | 1,414 (2,975 excluded: < 2 pulses, mostly single-pulse or faint) |
| eligible pairs | 998,991 |
| excluded: positions consistent (same-source or ordinary-lensing region) | 38,632 |
| excluded: resolution mismatch \|Δk\| > 2 (durations differing by more than about ×4–8) | 212,650 |
| **scored** | **747,709** (calibration–calibration 185,392; validation–validation 188,055; cross-half 374,262) |
| above the M2 threshold | 64 (74.8 expected from the per-pair rate) |
| **pass consistency (c = 10⁻³) and delay ≥ 1 d** | **0** |

- **No candidates.** Of the 64 pairs above threshold, 63 have χ² p ≤ 4 × 10⁻¹³.
  - The closest call has p = 2.2 × 10⁻⁶, still more than two orders of magnitude below c (pair listed in
    `search.json`: bn160724444 / bn260107831, s = 0.505 at the grid edge, two-pulse bursts, 9.5 yr apart).
  - 11 of the 64 sit at the grid edges (s = 0.505 or 1.98).
  - High similarity with strongly inconsistent residuals is what generic shape similarity between unrelated
    bursts looks like; a stretched copy would show consistent residuals.
- No threshold was changed after the search.

**8.6 Efficiency of the frozen chain** (`efficiency.json`; seed 4242; templates from both halves):

| injected | n | selected (pass) | recovered given selected | total |
|---|---|---|---|---|
| deformed, s0 ≠ 1 | 1,440 | 34.0 % | 49.9 % | **16.9 %** |
| undeformed, s0 = 1 | 180 | 35.6 % | 42.2 % | 15.0 % |
| brightness ratio 1 / 0.5 / 0.3 / 0.2 / 0.1 (all s0) | 324 each | 68 / 47 / 31 / 19 / 7 % | 54 / 48 / 42 / 44 / 48 % | 37 / 23 / 13 / 8 / 3 % |

**8.7 Conditional upper limit.**
- **The population it assumes.** Twins distributed like the injection set:
  - s0 uniform over the nine values in [0.5, 2];
  - a Band spectrum (α = −1, β = −2.3, E_peak ∈ {100, 200, 400} keV);
  - brightness ratio uniform over {1, 0.5, 0.3, 0.2, 0.1};
  - partner at a position-inconsistent sky location, with GBM-like noise and 15 % gain jitter.
- **Result.** Zero candidates give f95 < 3.0 / (1,414 × 0.169) = **1.3 × 10⁻² deformed-twin pairs per eligible
  burst** (5.8 × 10⁻³ at brightness ratio 1, where the total efficiency is 0.367).
- **Scope.** This limit is **not** valid for:
  - other stretch, spectral or brightness distributions;
  - single-pulse bursts;
  - pairs excluded by the position or resolution cuts;
  - any wormhole population in general.

**8.8 What remains invisible to this search:**
- **Single-pulse and faint bursts:** 2,975 of 4,389 reduced bursts. The search requires ≥ 2 resolved pulses on
  both members.
- **Same-position pairs:** the S1 design excludes them. Ordinary lensing searches cover them (Ahlgren & Larsson
  2020).
- **Large duration mismatch:** stretches beyond about ×2, or |Δk| > 2; 212,650 eligible pairs were not scored.
- **Faint copies.** Total efficiency is 3 % at brightness ratio 0.1, and selection loses most of them.
- **Strong spectral evolution, or frequency shifts beyond what two bands can absorb.** The two-band data can
  tolerate a shift but cannot test E → E/s. That needs CSPEC/TTE per candidate.
- **Time-dependent stretch** (non-uniform s), dispersion, or partial copies (one pulse of a burst).
- **Copies hidden in data gaps**, Earth occultation or SAA passages, i.e. no GBM trigger at all.
- **Other messengers.** FRBs (CHIME Catalog 2 is still open), neutrinos and GW are not in this screen.
