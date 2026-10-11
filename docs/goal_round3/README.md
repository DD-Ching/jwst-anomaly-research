# Goal round 3: which finite-speed hidden-influence criterion is worth an experiment?

Owner brief 2026-10-11 (continues #152 and #154). Question: which criterion most plausibly becomes an experiment
that can discriminate new physics?
- Bancal et al. 2012, S ≤ 7 (Nature Phys. 8, 867; arXiv:1110.3795v2).
- Li, Hu, Deng & Scarani 2026, LC4: S4 ≤ 6 and LC5: S5 ≤ 10 (arXiv:2608.05271v1).

Both papers were read in full, including App. E and the Supplemental Material.

Labels: **[paper]**, **[derived]** (this round), **[literature]** (experimental capability, read by a subagent this
round; see SOURCES.md), **[ASSUMPTION]**, **UNKNOWN**. Nothing here is evidence for new physics. Standard quantum
mechanics (QM) is the no-signalling baseline.

Code: `scripts/goal_hic_witnesses.py`; tests: `tests/test_goal_hic_witnesses.py`; output:
`results/goal_round3/hic_witnesses.json` (about 30 s). The round-2 report revisions requested in the same brief
are in `docs/goal_round2/README.md` §3.1, §5 and §8.

## 1. Verdict

- **Most promising criterion: LC4.** It is the most noise-robust under every realistic noise model computed, needs
  ≥ 50× fewer trials than Bancal (3.5 × 10³ against 2.0 × 10⁵ at 5σ for ideal states), uses 6 low-weight terms, and fits Bancal's 4-site geometry unchanged.
  - LC5 is better only under pure white noise.
  - The Bancal criterion is not a practical target: it is violated only above 97.2 % visibility and with less than
    0.5 % readout error.
- **What an LC4 experiment can discriminate.** Bancal Experiment v1 tests QM against **no-signalling** v-causal
  models. These models cannot carry messages. A null result removes them, so any surviving finite-speed explanation
  must signal.
- **The communication-relevant test is out of reach.** It is the direct observation of signalling in B ∼ C
  ("R") trials, and with an unknown frame it needs ~10¹³ trials (§6).
- **Decision (§8):**
  - Keep the LC4 Experiment-v1 specification as the only design in this line.
  - Wait for one capability: four spacelike-separated parties with no postselection, visibility ≳ 0.95 and detection
    ≳ 0.97.
  - Pause further expansion of this line.

## 2. Do the criteria answer the same question? [paper]

| | Bancal 2012 | Li 2026 LC4 | Li 2026 LC5 |
|---|---|---|---|
| Model set | NS ∩ local(BC given A, D) (App. C, Lemma 1) | HIC_AD = NS ∩ CL(BC\|AD) (Eq. 1) | HIC_ADE (Eq. 1) |
| Spacetime premise | A < D < (B ∼ C), with App. E i)–iii) | the same, with E = AD (Fig. 1) | A ≺ D ≺ E ≺ (B ∼ C) with five separate sites and a 2-D layout (S8–S10); or Bancal's 4-site layout with D and E as one system (main text after Eq. 6) |
| How QM marginals reach the B ∼ C case | App. D: random B and C timing; QM holds in every v-connected chain | S18–S19 delay lemma (the same argument) | the same |
| Measured data | ABD and ACD marginals only (no BC term) | ABD and ACD (supports AB, ABD, CD, ACD) | ABDE and ACDE |
| State | 4-qubit state of App. B | linear cluster LC4 | linear cluster LC5 |
| Extra stated assumptions | measurement independence; continuous random timing or λ-independence (App. E) | measurement independence, no postselection, stable shared marginals across blocks (SM I.D) | the same |

- **They test the same model class under the same premises.** Their predictions differ only through the state, the
  measurements and the inequality.
- Li et al. note that S4 ≤ 6 is also the ordinary local bound (main text, below Eq. 5). A violation measured without
  the v1 timing structure is therefore just a Bell violation. It says nothing about v-causal models unless the
  App. E conditions hold.
- **LC5 needs a 2-D layout only with five separate sites.** In 1-D, a middle early party never has a collection
  point (S21). Li et al. also allow LC5 in Bancal's 4-site layout, with D and E treated as one system that has two
  four-outcome measurements and the same bound of 10.
- **Frame-free protocol.** Bancal's eight rotated configurations cover every frame direction. The construction uses
  only the roles A, D early and B, C late, so it carries over to LC4 [derived]. With D and E merged it plausibly carries
  over to LC5 as well; that has not been derived.

## 3. Three different outcomes and what each supports [derived from the papers]

| Observation (Experiment v1: ABD from B < C trials, ACD from C < B trials, random B and C times) | Supports | Does not support |
|---|---|---|
| **S above the bound, at the QM value** (QM's prediction) | Excludes **no-signalling** v-causal models for every finite v, in the frames the configurations cover. | Anything about signalling v-causal models: they reproduce QM marginals there too (App. D), so they predict the same S. It is not a channel. |
| **S at or below the bound in v1, and also in other measurements of the same ABD and ACD marginals** | A deviation from QM in the four-party marginals. Under λ-independence of the timing (App. D, Eq. 10), the ABD marginal does not depend on when C is measured. So a no-signalling v-causal model has to break QM even in v-connected configurations, i.e. break Bancal's QM-in-chains premise. This would be new physics only if every loophole of §6 is closed. It would also appear in ordinary lab measurements of the same marginals; it is not a v1-only signature. | Communication: these models do not signal. |
| **S at or below the bound only with genuine continuous random B and C timing, but QM values with pre-set timing** | A model whose hidden variables use prior knowledge of the measurement times (the loophole in App. E). Otherwise retrocausality or an instrument artifact. | A v-causal model obeying λ-independence: it cannot produce this pattern. |
| **Setting dependence of BCD on x, or of ABC on w, in B ∼ C trials only, absent in v-connected control trials** | Direct, controllable superluminal signalling: the communication-relevant result. | Anything when it also appears in control trials: then it is an instrument effect (§6). |

**Consequence [derived].** Under λ-independence, v1's added value over a lab measurement of the same marginals lies
in replacing assumptions with measurement:
- spacelike separation removes ordinary communication;
- continuous random timing removes timing knowledge.

Under Bancal's Experiment-v3 hypothesis, one high-fidelity device would already decide the no-signalling question
(LC4 marginals at visibility > 0.879). No such published LC4 dataset was found (UNKNOWN). Either way, the outcome
would not test the signalling models that matter for communication.

## 4. Noise robustness [derived; positive controls reproduce the papers]

Each entry is the largest noise at which the witness is still violated, with everything else ideal. For visibility
and detection efficiency it is the smallest value that still violates.

| | Bancal | **LC4** | LC5 |
|---|---|---|---|
| Quantum value / bound (gap) | 7.2014 / 7 (2.9 %) | 6.8284 / 6 (13.8 %) | 11.657 / 10 (16.6 %) |
| Terms (max body) | 23 (3) | 6 (3) | 14 (4) |
| White-noise visibility (Li et al.: 0.9720 / 0.8787 / 0.8579) | 0.9720 | 0.8787 | **0.8579** |
| Local dephasing p (Z flip per qubit) | 0.0078 | **0.0416** | 0.0396 |
| Local depolarising q | 0.0098 | **0.0506** | 0.0440 |
| Amplitude damping γ | 0.0128 | **0.0713** | 0.0601 |
| Symmetric readout flip | 0.0049 | **0.0253** | 0.0220 |
| Detection efficiency, no-click → fixed outcome (no postselection) | 0.986 | **0.949** | 0.956 |
| Readout flip at visibility 0.98 / 0.95 | 0.0014 / none | **0.021 / 0.015** | 0.019 / 0.015 |
| Detection efficiency at visibility 0.98 / 0.95 | 0.996 / none | **0.957 / 0.969** | 0.961 / 0.970 |

- **LC4 versus LC5.** LC4 beats LC5 under every local noise model computed. LC5's lower white-noise threshold does
  not survive its four-body terms, as Li et al. anticipate in SM VIII.B.
- **Bancal is out.** Its criterion fails at the visibilities that the best published sources reach (§7).
- **Assigning no-clicks keeps the bound valid.** Mapping a no-click to a fixed outcome is a local relabelling, so the
  bound stays valid. Discarding no-clicks (postselection) does not keep it valid.

## 5. The exact linear programme: best threshold, forced signal, the maximal Bancal state [derived]

LP over full distributions, with:
- the ABD(E) and ACD(E) marginals equal to the observed ones;
- BC local given every early event;
- either full no-signalling (feasibility) or the smallest forced signal.

For the Bancal case it reproduces round 2 to 10⁻⁸.

**Best possible threshold for these data.** The lowest white-noise visibility at which the observed marginals still
force signalling:
- Bancal state: 0.9720;
- LC4: 0.8787;
- LC5: 0.8579.

These equal each witness's own threshold. No other no-BC witness does better *for the same state and settings*.

**The maximal Bancal state is not needed.** Its threshold, 0.9526 (Li et al. p.2), still loses to LC4 (0.879)
under white noise. The 23 terms of S cannot plausibly reverse the order under local noise. Computing that state
would not change the choice, so it was not computed.

**Forced signal** (smallest largest mean total variation over the early senders, in R trials; any v-causal model
that reproduces these marginals):

| Marginals | Bancal | LC4 | LC5 |
|---|---|---|---|
| ideal | 0.0063 | **0.0259** | 0.0194 |
| visibility 0.99 / 0.97 / 0.95 / 0.90 | 0.0040 / 0 / 0 / 0 | 0.0238 / 0.0195 / 0.0152 / 0.0046 | 0.0181 / 0.0153 / 0.0126 / 0.0058 |
| dephasing 0.02 | 0 | 0.0133 | 0.0096 |
| readout flip 0.01 | 0 | 0.0154 | 0.0103 |

- For the ideal 4-party states the forced signal equals (S_Q − B)/32. LC4 forces 4.1× more signal than Bancal, so a
  direct signalling test needs ~17× fewer trials.
- The LP is a relaxation of the model class, so these are lower bounds.

## 6. What can make a fake signal, and how the design removes it [derived]

**Postselection and detection loss**
- Coincidence postselection voids the HIC bound: the proofs assume no postselection (Li SM I.D).
- It can also turn setting-dependent loss into apparent signalling.
- The design rule: every trial is kept, with no-clicks assigned as in §4. That needs per-detector efficiency
  ≥ 0.949 (LC4, ideal state), or ≥ 0.969 at visibility 0.95.

**Choosing which trials to analyse**
- v1 records ABD only from B < C trials and ACD only from C < B trials. The classification uses the recorded times
  with Bancal's 2δ margins, never the outcomes.
- B and C times must be drawn continuously at random, or the λ-independence hypothesis must be stated (App. E).
  Otherwise a model could behave classically exactly when B ∼ C is possible.

**Ordinary signal crosstalk**
- An electronic, optical or clock link from A's or D's setting to B, C or D's outcome produces real setting
  dependence at ≤ c.
- Removed by lab spacelike separation of every setting choice from every other party's outcome event. Bancal's
  geometry already requires this for condition iii.
- Lab-scale single-device runs (ions, superconducting chips) cannot meet it. There, a violation tests only Bancal's
  Experiment-v3 hypothesis.

**Instrument drift and shared marginals**
- The AD marginal enters both families. Li et al. require stable shared marginals across acquisition blocks.
- Interleave B- and C-family trials randomly, and test the equality of the shared AD marginal as a pre-registered
  control.

**Finite statistics**
- The per-trial v1 estimator (unbiased; exact variance from the distribution) needs the trial counts below.
- Gaussian 5σ assumes i.i.d. trials. The Azuma–Hoeffding count at p = 10⁻⁶ holds with memory between trials, given
  the per-trial range.
- Both are sufficient for these tests, not minimal.

| Visibility | Bancal 5σ / Azuma | **LC4** 5σ / Azuma | LC5 5σ / Azuma |
|---|---|---|---|
| 1.00 | 2.0 × 10⁵ / 2.5 × 10⁶ | 3.5 × 10³ / 4.1 × 10⁴ | 5.9 × 10³ / 4.7 × 10⁴ |
| 0.98 | 2.5 × 10⁶ / 3.0 × 10⁷ | 5.2 × 10³ / 5.9 × 10⁴ | 8.1 × 10³ / 6.3 × 10⁴ |
| 0.95 | no violation | 1.1 × 10⁴ / 1.2 × 10⁵ | 1.4 × 10⁴ / 1.1 × 10⁵ |
| 0.92 | no violation | 3.3 × 10⁴ / 3.6 × 10⁵ | 3.2 × 10⁴ / 2.4 × 10⁵ |

Frame-free coverage multiplies these by the 8 configurations of App. E.

**Direct signalling test**
- The signal is only expected in R trials. Controls are the timelike-ordered trials, where QM, and hence no
  signal, is predicted. The instrument floor measured in the controls must be below f·δ.
- With β_max at the CMB speed [ASSUMPTION] and r = 5 × 10⁶, the worst-case R fraction is f = 8.1 × 10⁻⁵
  (`docs/goal_round2` §5). At δ = 0.015 (LC4, visibility 0.95), a 5σ detection needs ~25/(fδ)² ≈ 1.7 × 10¹³ trials
  [estimate].
- Without a bound on the frame speed there is no positive guarantee.

## 7. Existing capability [literature; numbers as read by a subagent; anything not found is UNKNOWN]

| Platform (source) | Fidelity or visibility | Rate | Spacelike? | Postselected? |
|---|---|---|---|---|
| SPDC 4-photon GHZ (Zhang et al. arXiv:1702.04130, Fig. 3) | F = 0.981 | 0.42 /s | no | yes |
| SPDC 4-photon GHZ Bell test (Mao et al. PRL 129, 150401 (2022)) | F = 0.974, visibility 0.972 | 33,252 events over ~14,741 × 20 s cycles | no; settings switched every ~20 s | yes |
| 4-photon cluster (Kiesel et al. PRL 95, 210502 (2005); Walther et al. Nature 434, 169 (2005)) | F = 0.74 / 0.63 | 150 /h; ~0.2 /s | no | yes |
| 3-photon GHZ, locality closed (Erven et al. Nat. Photon. 8, 292 (2014)) | Mermin M = 2.77 | 0.5 /s at the stations | yes (~0.7–0.8 km) | yes (link efficiencies 14–33 %) |
| 3-node NV GHZ (Pompili et al. Science 372, 259 (2021)) | F = 0.54 | ~1 per 90 s | no (~30 m) | heralded |
| Loophole-free 2-party: Hensen 2015; Rosenfeld 2017; Storz 2023 | S = 2.42; S = 2.22 (detection ≥ 0.98); S = 2.075 (readout 97.6–99.05 %) | 245 events in 220 h; 1–2 /min; 12.5 kHz | yes (1.3 km; 398 m; 30 m) | no |
| 4-ion GHZ in one trap (Monz et al. PRL 106, 130506 (2011), Table I) | F = 0.957 | UNKNOWN | no | no |

UNKNOWN: superconducting 4/5-qubit cluster fidelity, 5-ion GHZ fidelity, any photonic 4-fold rate ≥ 1 kHz since 2023
(the search was not exhaustive). No experiment of any finite-speed multipartite test was found, which agrees with
Li et al. p.2 (three searches; not a full crawl of the citing papers).

**Gaps for LC4 Experiment v1**
1. **No platform combines the requirements**: four parties spacelike-separated from each other, no postselection with
   per-detector efficiency ≥ 0.969, and visibility ≥ 0.95.
   - Photonic sources have the visibility but 14–76 % system efficiency.
   - Heralded matter nodes have the efficiency, but their multi-node fidelity (0.54) is far below threshold.
2. **Rate.** At visibility 0.95, LC4 needs 1.2 × 10⁵ trials (Azuma) per configuration.
   - At the 3-node GHZ rate (1 per 90 s) that would be ~4 months per configuration.
   - At the postselected photonic 0.42 /s it would be ~3 days, but that run would not be a valid test.
3. **Per-trial random settings and continuous random B and C timing** at four stations. These are demonstrated
   separately (Pockels cells, as in Erven and Shalm), but not together at four stations.
4. **Timing is not the bottleneck for v1.** v1 needs no synchronisation (App. E). The windows that matter are the
   spacelike margins: 4.27 µs at 1.3 km (Hensen) and 109 ns at 33 m (Storz).

## 8. Decision and what would change it

**Decision.**
1. **Continue: the LC4 Experiment-v1 specification** (§2–§6) is the only design in this line worth keeping.
2. **Wait for capability.** No new experiment design work is started until a platform reports one of the trigger
   conditions below.
3. **Pause expansion of the low-latency line through v-causal models.** No more theory units on it:
   - the Bancal criterion is dominated;
   - the maximal state cannot change the choice;
   - the direct signalling test needs ~10¹³ trials, or a known frame.

**Triggers that change this decision**
1. **A published platform with four spacelike-separated parties**, no postselection, detection efficiency ≥ 0.97
   and 4-party visibility ≥ 0.95 (equivalently F ≳ 0.95). Then write the full LC4 v1 protocol, with rates and the
   8 configurations.
2. **Any report of S4 ≤ 6 (or S ≤ 7) in a v1 configuration while timelike-ordered controls give the QM value.**
   That has top priority: vet loopholes, then design the direct signalling test.
3. **Independent evidence fixing the preferred frame** (speed and direction), or a bound β_max far below the CMB
   speed. The worst-case R fraction ≈ 1/(2 r β_max) then rises, and the trial count falls as β_max². With a fully known
   frame and B–C timing better than d/v, R trials can be scheduled (f → 1) and a 5σ direct test needs ~25/δ² ≈ 10⁵ trials. Recompute it.
4. **A model that fixes the branch and the signal amplitude per direction.** Then a one-way channel test with the
   receivers clustered (round 2 §3.2) becomes designable.
5. **A photonic source with system efficiency ≥ 0.97 and a 4-fold rate ≥ 1 kHz.** Re-cost the experiment.
6. **A theory paper giving a no-BC witness with a lower threshold under local noise than LC4**, for example from a
   full facet list (Li et al. Discussion). Recompute §4.

## 9. Bookkeeping

- **Independent review** (separate agent, own code, against both papers). Confirmed: the states, the observables and
  the terms; every threshold in §4; the LP encoding and the values of §5 (with its own LP); the estimator and the trial
  counts; the R fraction (exact integration and Monte Carlo). Fixed from it:
  - §3 row 2 was logically wrong. Under λ-independence, a deviation cannot be specific to v1; the row is now split.
  - LC5 also fits Bancal's 4-site layout.
  - With no bound on the frame speed, the worst-case R fraction is exactly 0 for windows shorter than d/c.
- **Time.** 2026-10-11 05:40 → about 06:10 UTC, including the round-2 revisions.
- **Compute.** About 30 s for the full script (LPs of up to 1024 variables). The tests take about 10 s.
- **Data.** None. Papers only (arXiv PDFs in the scratchpad, not committed).
- **Not done.**
  - Noise correlated across qubits.
  - The frame-free 8-configuration protocol for LC5, with D and E merged.
  - The detection-loss LP: the witness thresholds are computed, but the LP threshold under loss is not.
  - A full literature crawl of the experiments citing Bancal et al.
