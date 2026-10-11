# Goal round 2: when could a remote receiver read the sender's bit in a Bancal-type v-causal model?

Owner brief 2026-10-11 (continues #152, `docs/goal_round1/README.md`). Labels: **[paper]** Bancal, Pironio, Acín,
Liang, Scarani & Gisin, Nature Phys. 8, 867 (2012), arXiv:1110.3795v2, read in full (main text, Fig. 3, App. B–E);
**[derived]** this round; **[estimate]**; **[ASSUMPTION]**. Nothing here is evidence for new physics. Standard
quantum mechanics (QM) is the no-signalling baseline throughout.

Code: `scripts/goal_bancal_channel.py`; tests: `tests/test_goal_bancal_channel.py`; outputs:
`results/goal_round2/bancal_channel.json` (L = 20 ly) and `bancal_channel_100ly.json`. Run time about 2 s.

## 1. The model and what the theorem actually provides [paper]

**Model and premises**
- A preferred frame exists. Influences propagate inside v-cones (v = r c > c).
- The model reproduces the QM predictions whenever each later measurement lies inside the v-cones of the earlier
  ones. This holds for A < D < B < C and A < D < C < B.

**Configuration R = A < D < (B ∼ C), the conditions of App. E**
- i) t_A, t_D < t_B, t_C.
- ii) B and C lie outside each other's v-cones.
- iii) Two points exist:
  - A′, inside the future light cones of A, B and C but outside D's;
  - D′, inside those of D, B and C but outside A's.

**What follows in R** (App. C, D; Lemma 1)
- The BC correlations conditioned on (a, x, d, w) are local.
- The ABD and ACD marginals equal the QM ones.
- So S = 7.2014 > 7 forces a violation of no-signalling: the BCD marginal depends on A's setting x, or the ABC
  marginal depends on D's setting w, or both (Fig. 3 caption).

**What the theorem does not provide**
- Which of the two branches occurs.
- How large the dependence is.
- Any rate.
- Any statement about layouts that break i)–iii).

## 2. Positive control

| Quantity | This round | Paper |
|---|---|---|
| ⟨Ψ\|Ψ⟩ | 1.000 | — |
| S (Lemma 1, Eq. 4) from the App. B state and observables | 7.20139 | 7.2014 |
| QM: BCD marginal dependence on x; ABC marginal dependence on w | exactly 0 | — (QM is no-signalling) |
| LP with no-signalling imposed on top of the R constraints | infeasible | Lemma 1 |

## 3. Event timeline, counted from the sender's free choice [derived]

Sender A at Earth, x = 0. D at the star, x = L. The privileged frame is taken as the rest frame of both
[ASSUMPTION]. Times are in L/c.

| # | Event | Fig. 3, r → ∞ | Needs |
|---|---|---|---|
| E0 | A chooses bit x and measures | t = 0 | A's share of a 4-party state |
| E1 | D measures with its own setting w | t = 1/r (inside A's v-cone) | D's share |
| E2 | B and C measure | t = 2/(1+r), with B ∼ C | B's and C's shares; \|t_B − t_C\| < \|x_B − x_C\|/v |
| E3 | B and C broadcast outcome and setting at light speed | from E2 | ordinary radio |
| E4 | the B, C and D records are together at the star (D′) | max(t_B + L − x_B, t_C + L − x_C) | E3 |
| E5 | x decoded at the star | E4 + statistics (§6) | many trials |
| Back | if instead ABC depends on w: A, B and C records together at Earth (A′) | max(t_B + x_B, t_C + x_C) | E3 |

**Fig. 3 itself**
- Positions [paper]: d_B = ¼(1 + 1/r) + 1/(1+r), d_C = ¾(1 + 1/r) − 1/(1+r), d_D = 1.
- Times: t_A = 0, t_B = t_C = 2/(1+r), t_D = 1/r.
- At L = 20 ly the stations sit at 5 and 15 ly for r → ∞, and E4 comes at **15.00 yr** in both directions. That
  confirms the owner's estimate. Light takes 20 yr, so the effective speed is only 1.33 c, however large v is.

| Layout (L = 20 ly) | r = 2 | r = 10 | r = 10⁴ | r = 6.3 × 10⁸ |
|---|---|---|---|---|
| Fig. 3: A → star / star → A (yr) | 19.17 / 19.17 | 16.32 / 16.32 | 15.00 / 15.00 | 15.00 / 15.00 |
| B and C at the midpoint, 0.04 ly apart | 25.0 (not superluminal) / 15.0 | 13.0 / 11.0 | 10.02 / 10.02 | 10.02 / 10.02 |
| best layout keeping both directions superluminal (grid) | 17.6 | 12.0 | 10.1 | 10.1 |
| B and C next to D (0.02–0.04 ly) | 10.1 / 20.0 | 2.0 / 20.0 | 0.04 / 19.98 | 0.04 / 19.98 |

At 100 ly the same layouts give 75.0, 50.1 and 50.5 yr, and 0.2 / 99.9 yr for B and C next to D, at r = 6.3 × 10⁸.

### 3.1 No layout guarantees better than half the light time [derived]

The theorem allows either branch, so a guaranteed channel needs both collection points (condition iii). Put B and
C anywhere in space, and let L be the A–D distance.
- The star collects the B and C records no earlier than max_k(t_k + |x_star − x_k|).
- Earth collects them no earlier than max_k(t_k + |x_Earth − x_k|).
- By the triangle inequality, |x_Earth − x_k| + |x_star − x_k| ≥ L for each station k.
- So the worse of the two latencies is ≥ L/2 + O(L/v).

The grid search reaches 0.505 L (B and C at the midpoint). For 20 ly that is ≥ 10 yr; for 100 ly, ≥ 50 yr. So,
**under the theorem's guarantee alone, no arrangement gives low latency.** That includes putting all the receiving
stations on one planet: as the last table row shows, that keeps one direction superluminal only marginally
(19.98 yr against 20).

### 3.2 Concentrating the receiving stations needs an extra assumption [derived]

With B, C and D next to the star, the A → D′ latency falls to the v-cone limit: B and C must lie inside A's future
v-cone, so t_B ≥ L/v (1 s at 20 ly for v = 6.3 × 10⁸ c). It also includes the local record collection.

That helps only if the realised branch is "BCD depends on x" in this layout, and the theorem does not supply that
(§4). **The hidden-influence speed v bounds the communication speed from above. It is not the communication speed.**

## 4. How large is the signal, and in which direction? [derived, exact LP]

**Linear programme** over P(abcd|xyzw) ≥ 0 (256 variables).
- Constraints:
  - normalisation;
  - ABD = QM for every z, ACD = QM for every y (App. D);
  - for each (a, x, d, w), the conditional BC correlations are local (the eight CHSH facets in joint form; App. C).
- Signal measure: total-variation distance between the receiver's record distributions under the two values of
  the sender's setting, averaged over the 8 setting combinations of the other parties.

| Quantity | Value |
|---|---|
| smallest A → D′ signal (BCD dependence on x) | 0 |
| smallest D → A′ signal (ABC dependence on w) | 0 |
| smallest max(A → D′, D → A′) | **0.0063** |
| smallest sum | 0.0126 |
| smallest A → D′ signal if D → A′ is silent | 0.0160 |
| smallest D → A′ signal if A → D′ is silent | 0.0126 |

**What this means [derived]**
- Any v-causal model that reproduces the QM marginals of this state signals with a mean total variation of at
  least 0.0063 in at least one direction.
- Either direction may carry nothing.
- A two-way channel is therefore **not guaranteed**, even in a single R configuration.
- Two mirrored configurations, one with A on Earth and one with A at the star, give a two-way channel only if the
  model's branch choice follows the parties' roles rather than a spatial direction. That is an extra assumption
  [ASSUMPTION], and the theorem does not supply it.
- The numbers hold for this state and these observables only. Other states give other bounds; the maximal state
  (S = 7.3481) was not computed.

## 5. Sync, deployment, record arrival and statistical time, kept separate [derived / estimate]

**Deployment.** Every trial consumes one 4-party entangled state.
- The B, C and D shares must reach their stations at ≤ c before use: about 10 yr to the midpoint and 20 yr to the
  star for 20 ly.
- They must stay coherent until they are used.
- Nothing can top the stock up faster than light.

**Single-record arrival.** This is E4 in §3: 15 yr for Fig. 3 and ≥ 10 yr for any guaranteed layout. With an extra
directional assumption (§3.2) it falls to about L/v.

**Trials per bit.** With a per-trial total variation δ, a threshold test on the best event decides a bit at error
ε after n = 2 ln(1/ε)/δ² trials (Hoeffding).
- At δ = 0.0063: n ≈ 3.5 × 10⁵ for ε = 10⁻³ and 7.0 × 10⁵ for ε = 10⁻⁶.
- An 8 kbit/s voice link at ε = 10⁻³ would consume about **2.8 × 10⁹ pre-positioned 4-party states per second**
  [estimate].
- Trials can run in parallel, so the statistical *time* is n divided by the trial rate. It is limited by detector
  rates and the size of the stock.

**B ∼ C has to occur.** Signalling appears only in R. In B < C or C < B configurations the model reproduces QM, which
does not signal.
- So B and C must be measured within |x_B − x_C|/v of each other in the privileged frame. For 0.04 ly and
  v = 6.3 × 10⁸ c, that window is 2.0 ms.
- Arranging this needs the frame to be known: the simultaneity ambiguity at β ≈ 1.23 × 10⁻³ (CMB) is
  β·d/c ≈ 1,560 s for that separation.
- With random timing (App. E) over a window T_w, only a fraction ~(2|x_B − x_C|/v)/T_w of trials are in R, and the
  signal is diluted by that factor.
- For B and C 1 km apart at v = 6.3 × 10⁸ c the window is 5.3 fs.

## 6. Missing physical assumptions

None of these is in the theorem. All are needed for a rate or an error rate in a chosen direction:
1. Which branch the model realises in a given layout and frame, and whether it depends on geometry. This is
   needed for any one-way guarantee and for concentrating the receiving stations.
2. The signal amplitude per direction. The LP gives only the joint bound 0.0063 (at least one direction).
3. The value of v and the privileged frame. Both are needed to schedule B ∼ C, and v is unknown above the bipartite
   bounds (≳ 5 × 10⁶ c in near-CMB frames; round 1).
4. The state dependence: the bound changes with the state and the observables.
5. Decoherence of stored shares over decades, and the stock needed per bit.
6. Whether a signalling model is consistent with the rest of physics. Signalling with a privileged frame avoids
   paradoxes, but it contradicts QM's no-signalling, which is tested to high precision in ordinary settings.

## 7. Answers

**Which model and layout deserve an experimental design?**
- The Bancal v-causal class in configuration R, with B and C adjacent (the central layout is the best guaranteed
  one).
- The interstellar layouts are not the target: under the theorem they cannot beat L/2 (10 yr at 20 ly).
- A laboratory-scale R experiment is what can separate the model class from QM.

**What does it predict?** QM predicts exactly zero setting dependence in the BCD and ABC marginals. Any v-causal
model that reproduces the QM marginals predicts a mean total variation of at least 0.0063 in at least one of them,
**in R trials only**.

**Minimal next verification steps**
1. **Bancal's Experiment v1** (App. E): S > 7 with random B and C timing.
   - It needs no synchronisation.
   - It excludes **no-signalling** v-causal models at any finite v. It cannot test signalling ones.
   - It has never been done (Li et al. 2026, p.2). It needs S > 7 after noise, from a gap of 0.2014.
2. **A direct test for signalling in R trials.**
   - It needs B ∼ C, i.e. a B–C timing better than d_BC/v.
   - To reach models not already excluded by the bipartite bounds (v ≳ 5 × 10⁶ c), d_BC/(cδ) must exceed 5 × 10⁶.
     For d_BC = 100 km that means δ < 67 ps; reaching 6.3 × 10⁸ c needs δ < 0.53 ps [derived].
   - A 1 km, 10 ps set-up reaches only v ≲ 3 × 10⁵ c, which is already excluded in those frames.
   - If the frame is unknown, the simultaneity ambiguity β·d_BC/c ≈ 0.41 µs at 100 km (CMB β) far exceeds that
     window. With random timing over ≳ 0.4 µs, only ~3 × 10⁻⁴ of trials fall in R, and the signal is diluted by that
     fraction.

   So the decisive direct test needs picosecond B–C timing over ~100 km plus either knowledge of the frame or a very
   large number of trials.

**Scope.** This round is limited to the Bancal class with this 4-qubit state, the theorem's premises, and 1-D or 3-D
layouts with point stations. It does not exclude other new physics, other signalling mechanisms, or models outside
the v-causal premises. Within that scope, **no layout found gives guaranteed low-latency or two-way
communication**. Low latency needs the additional assumptions in §6, items 1–3.

## 8. Bookkeeping

- **Time.** 2026-10-11 05:09 → about 05:30 UTC.
- **Compute.** About 2 s per run (7 LPs and the grids). The tests take under 1 s.
- **Data.** The arXiv PDF of 1110.3795v2 (in the scratchpad, not committed). No other downloads.
- **Not done.**
  - The maximal-violation state (7.3481) and other states.
  - The LP with the per-setting maximum instead of the mean.
  - Noise models.
  - Noise-robust witnesses (Li et al. 2026).
