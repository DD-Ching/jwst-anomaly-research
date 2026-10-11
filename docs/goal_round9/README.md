# Goal round 9: is there a next executable physics unit for the communication goal?

Owner brief (2026-10-11). The Einstein–æther line stays paused. Choose one worthwhile physics problem from the whole
project, and run a minimal check only if a clear prediction, adequate data and a light test all exist. Otherwise
report that there is no basis for a choice, and say which input is missing.

This round adds no new search, no G1, no change to existing searches or schedules, and no merges.

**Labels:**
- [lit]: stated in a cited paper.
- [rec]: recorded in this repository; the file is named.
- [obs]: an observational result already in the repository.
- [derived]: reasoned in this round.
- [conj]: a conjecture or hypothesis.

## 1. What was read

- Latest `main` (eed9b36): TASKS.md, DECISIONS.md, CHANGELOG.md, `docs/hypotheses/`, `docs/neutrino_frontier/`,
  `docs/exotic_limits.md`.
- Open PRs: #152, #154–#158, #161, #163 (goal rounds 1–8), #145, #150, #160 (MOA-II gb3, claimed by another
  session) and #61.
- Goal rounds 1–8 from their branches.

A read-only survey subagent compiled the inventory. The decision-relevant facts were re-checked here.

**Inventory in brief.** Every data program is null or limited by access:
- W3 MOA-II, 14 fields [obs; `results/w3_moa/limits_combined.ecsv`];
- S1, S2 and S3 (D-070–D-072);
- D1 (D-073);
- E1 (D-074, with a CHIME excess reproduced by a rate-modulated null);
- the neutrino ghost test NF-H04 (D-077, D-079).

None of them records a map from a positive result to a controllable superluminal signal [rec; goal round 1 §§1.5,
5.3].

The goal rounds have each paused their communication candidate:
- v-causal hidden influences (round 3);
- the beyond-Horndeski wormhole (round 4);
- Scharnhorst and Lentz (round 5);
- Krasnikov (rounds 6–7);
- Einstein–æther (round 8).

## 2. Three unresolved, concrete problems

### P1. Do finite-speed hidden influences exist (multipartite v-causal test)?

| Question | Answer |
|---|---|
| Missing evidence | A four-party experiment with spacelike-separated parties, no postselection and continuous random timing that reaches S4 ≥ 6.2 within ≲ 10⁶ Azuma trials per configuration [rec; round 3 §6a, §8]. No such dataset or platform was found (round 3; the search was not exhaustive). |
| Prediction different from standard physics | QM: S4 = 6.83 in configuration R. No-signalling v-causal models with finite v: S4 ≤ 6 [lit: Li, Hu, Deng & Scarani 2026; Bancal et al. 2012]. Signalling v-causal models: setting dependence in B ∼ C trials only (round 3 §3). |
| Decision-changing result | An S4 ≤ 6 counts only against a QM prediction computed for the actual state, noise, detection efficiency and statistics of that run (round 3 §§4, 6a); a low S4 alone is not evidence that QM fails. A low S4 that stays below that prediction once all of these are accounted for, or B ∼ C-only setting dependence absent in control trials (round 3 trigger list), would reopen the line toward a channel design. A clean QM-valued violation excludes the no-signalling class for all finite v, but says nothing about the signalling class that matters for communication [rec; round 3 §3]. |
| Minimal verification | Data from such an experiment. **Not available.** A single-device LC4 dataset would test only Bancal's Experiment-v3 hypothesis, and "would not test the signalling models that matter for communication" [rec; round 3 §3]; no such dataset was found either. |
| Gaps to communication | A model that fixes the signalling branch and amplitude per direction. The preferred frame. The ≥ L/2 latency guarantee (10 yr at 20 ly) and ~10¹³ trials per bit in the unknown-frame case [rec; round 2]. |

### P2. Does negative or exotic mass exist in nature (inverted microlensing, W3)?

| Question | Answer |
|---|---|
| Missing evidence | Any confirmed event. Fourteen MOA-II fields are null, with Γ₉₅ ≈ 2.1 × 10⁻⁸–1.2 × 10⁻⁷ per star per year [obs; TASKS.md, D-068]. There is also no abundance model to turn the limit into a decision [rec; goal round 1 §5.3]. |
| Prediction different from standard physics | The templates differ by lens type [rec; `docs/exotic_lensing.md`]. **W3** is a negative-mass (repulsive) microlens, or a repulsive n ≥ 2 variant: the source flux drops to zero between two caustic spikes. Standard lensing never dims a source below baseline. An **Ellis wormhole** (attractive n = 2) gives a different light curve: shallow, time-symmetric ~4 % "gutters" (W4), which the record calls not practical to detect, plus static demagnification (W2). |
| Decision-changing result | A W3 event that survives every vetting test would be a **candidate**, not proof. It would still need vetting by `/vet-candidate`, independent photometry and an alternative-explanation review before it could count as evidence of exotic stress-energy. Even then it would not by itself supply a source for any specific shortcut geometry. Further nulls change nothing for the communication goal without an abundance model. |
| Minimal verification | MOA-II fields gb3–gb6, gb9, gb10, gb14, plus gb13 after the owner's decision on #145. Data **are** available, but this is an existing, scheduled search (gb3 is claimed in #160). Selecting it here would duplicate it or change its schedule. **Not selected.** |
| Gaps to communication | Even a confirmed exotic lens would not supply a controllable source, a stable shortcut, or the negative energy that a specific channel model needs. The ~10³⁸ M☉ figure belongs only to the round-6 Krasnikov model (20 ly, 1 µm radius, under its QI assumptions); it is not a general requirement. |

### P3. Does the Hořava / khronometric "instantaneous" term give a controllable remote response?

| Question | Answer |
|---|---|
| Changed premise | Local Lorentz invariance of the gravity sector: a dynamical preferred foliation T (the round-5 §5 "not covered" class). |
| How verified physics is recovered | The khronometric (BPSH IR) action is Einstein–æther with u_a = N∂_aT [lit: Jacobson, arXiv:1001.4823v3, the published version with its erratum; v3 fixes a typo in Eq. 10]. With small couplings it restores GR, PPN and GW170817 [lit: 1802.04303; 1311.7144]. |
| First experimental constraint | PPN α1 and α2, binary pulsars and GW170817, as in round 8. |
| Missing evidence | Settled at linear order by this round's check (§2b). The non-linear elliptic part is not settled (§2b). |
| Prediction different from standard physics | At linear order with ordinary matter: none beyond round 8 (§2b). |
| Decision-changing result | A demonstration that **local, conserved matter**, with the prior state and the outer boundary conditions held fixed, drives the slice-wise elliptic (instantaneous) part to a non-zero tidal response at a distant receiver. |
| Minimal verification | §2b: done at linear order; the non-linear step is costed there. |
| Gaps to communication | Round 8's amplitude trade-off at linear order, plus the non-linear question below. |

## 2b. P3 check (this round): khronometric instantaneous term, ordinary matter, tidal readout

**Three different things** [lit]:
1. **The propagating khronon (spin-0) mode** has speed c_χ² = (β + λ)/α (Blas–Sibiryakov arXiv:1110.2195v2 Eq. 9),
   the same mode as EA's spin-0.
2. **An extra instantaneous (elliptic) part.**
   - Perturbatively, it is the second term of the retarded propagator of the decoupled khronon (BS Eq. 14). It
     builds up from zero and falls as 1/r².
   - Non-perturbatively, it is the slice-wise elliptic equation for the lapse, with a free function f_im(T) fixed
     by boundary or regularity conditions on every slice (Bhattacharyya, Coates, Colombo & Sotiriou
     arXiv:1512.04899v2 §IV, Eqs. 34–39).
   - It is **not** the c_S → ∞ limit of EA: it exists at finite c_χ.
3. **The EFT range.** BS assume α, β, λ ≪ 1 and neglect the khronon's back-reaction on the metric ("decoupling":
   khronon on a fixed background; BS §2). The UV completion with higher spatial derivatives is not used.

**Which source produces the instantaneous term** [lit]. BS couple the khronon to an *external source vector*,
S_source = ∫√−g S^μu_μ (Eq. 10). Ordinary matter couples only to the metric, so this term is not ordinary matter,
and it is not assumed here.
- Its instantaneous piece is transverse and traceless, so it "will vanish if the source S^μ is a gradient of a
  localized scalar configuration" (BS §2).
- In EA it is cancelled by the transverse æther modes (BS App. A).

**Controls reproduced** [derived; `scripts/goal_khronometric.py`, `bs_controls`]:
- The partial-fraction split 1/(p²(αp0² − Bp²)) = (1/αp0²)[1/p² + B/(αp0² − Bp²)] gives the instantaneous piece
  p_ip_j/(αp0²p²). In position space, ∂_i∂_j(1/r) has zero trace and zero divergence away from the source.
- A gradient source leaves no 1/p² term, so there is no long-range instantaneous field.
- In EA, the longitudinal instantaneous piece (α → α + α0) plus the transverse one sum to zero.

**New: ordinary matter, linear order** [derived; `khronometric_vs_ea`].
- **Setup.** The khronometric theory is linearised about flat space with the æther at rest. Matter is coupled to
  the metric only, with a general conserved source. The transverse covariant æther is set to zero in the action
  (u_⊥ = v_⊥ + h_0⊥ = 0, so the h_0⊥ equation becomes E[h_0⊥] − E[v_⊥]), in gauges h_0i = 0 and unitary khronon.
- **Result.** The tidal tensor E_ij = R_{0i0j} at a receiver equals EA's exactly in the scalar and tensor sectors,
  for every c_i.
  - The vector sector differs only by EA's c_V pole, which is proportional to c13²; khronometric has no pole there.
  - With c13 = 0 (GW170817) the two theories give the **identical** linear response.
- **Consequences:**
  - At linear order, the BS-type instantaneous term does not appear in the tidal response to ordinary matter.
  - Round 8's results carry over exactly to khronometric theory at linear order: the causal c_S cone and the
    coupling κ ≈ c14/4 ∝ 1/c_S².
  - The two cases are **the same** under: linear order, matter coupled to the metric only, c13 = 0, tidal readout.
  - They **differ** only in the vector sector when c13 ≠ 0 (EA has the c_V pole, khronometric does not).
- **Interpretation** [derived, not computed]. This is consistent with BS: at linear order, matter feeds the
  khronon only through metric scalar potentials, i.e. gradient-type couplings, for which their instantaneous
  piece vanishes.

**Not settled: the non-linear elliptic part** [lit + open].
- **Spherical symmetry.** In spherical symmetry f_im(T) = 0 when the centre (or a universal horizon) is regular,
  and the khronometric equation then reduces to EA's (Bhattacharyya et al. Eqs. 34–35 and text). With the
  boundary conditions fixed, a local spherical source cannot switch it on.
- **Beyond spherical symmetry.** Bhattacharyya et al. argue only qualitatively. The lapse obeys a slice-wise
  elliptic equation (their Eq. 39) with two integration functions of T, one fixed by reparametrisation and one
  that they relate to the instantaneous mode and that is set by boundary or asymptotic conditions. Eq. 39 also holds in
  EA, which is causal. So ellipticity alone does not show that a local source drives a distant response.
- **Steps still missing:**
  - (a) a second-order (or fully non-linear) calculation of whether a localised, conserved matter source with
    fixed asymptotics and a fixed prior state produces a non-zero instantaneous component of E_ij at distance;
  - (b) its magnitude against the round-8 scale;
  - (c) whether it survives in the EFT range.
- **Cost:** step (a) is a second-order perturbative computation of the khronometric field equations about flat
  space (symbolic, several times the linear one) or a 3+1 numerical evolution in the preferred foliation (heavy).
  Neither was started.
- **Also not covered:** the UV (higher-spatial-derivative) Hořava dispersion, which has unbounded group velocity at
  high momentum.

**Answer to the P3 question.** At linear order, no: with ordinary matter coupled only to the metric, the khronometric
"instantaneous" term does not enter the receiver's tidal tensor, and the response equals EA's. Non-linearly the
question is open; it is not shown to exist and not shown to be absent. Using the word "instantaneous" in the
literature is not evidence of a usable channel.

## 3. Comparison and choice

| | Different prediction | Data or model available | Light test | Changes a communication decision |
|---|---|---|---|---|
| P1 v-causal | yes (S4 below the realistic QM prediction; B ∼ C signalling) | **no** (no four-party spacelike data) | — | yes, if positive |
| P2 exotic lenses | yes (W3: flux to zero) | yes, but already scheduled in another program | the existing pipeline | a survivor is a candidate only; nulls need an abundance model |
| P3 khronometric | not at linear order (§2b); non-linear part open | linear: done; non-linear: needs a second-order calculation | linear done; non-linear not light | only if the non-linear part gives a matter-driven response |

**Result.** Not enough basis to select a further executable unit beyond the P3 linear check done here.
- The next P3 step is a second-order or non-linear computation, which is beyond a light check. It is costed in §2b
  for the owner to decide.
- P1 lacks data.
- P2 is an existing, scheduled search.

**New inputs that would make a unit executable:**
1. **P1:** a public dataset or platform with four spacelike-separated parties, no postselection, and the round-3
   §6a capability. Or an S4 that stays below the realistic QM prediction. Or a published model that fixes the
   v-causal signalling branch and amplitude.
2. **P2:** an unexplained W3 survivor (then vetting), or a quantitative abundance or rate model for exotic lenses.
3. **P3:** an owner decision to fund the second-order khronometric calculation (§2b), or a published non-linear
   result with a matter-only source.

## 4. Plain-language summary

- **New evidence?** No new observational evidence. One new derived result: at linear order, the Hořava/khronometric
  "instantaneous" term does not reach a receiver's tidal measurement when the sender only moves ordinary matter.
  The theory then behaves exactly like round 8's Einstein–æther case.
  - [lit]: the controls (Blas–Sibiryakov instantaneous term and its EA cancellation; Bhattacharyya et al.'s
    elliptic mode).
  - [derived]: the khronometric = EA linear response.
  - [conj]: the non-linear question.
- **Judgement changed:** P3 no longer rests on "a recomputation would only repeat round 8" or "needs a strong-
  coupling completion first". The linear part is now checked; only the non-linear part remains, with a stated cost.
  The earlier remark that the MOA-II searches were the only near-term route to new evidence is withdrawn: no full
  comparison of all routes was made.
- **Still missing for communication** [derived]:
  - a mechanism that (a) a sender can drive with ordinary matter or a lab device, (b) produces a readable change at
    a receiver 20–100 ly away before light, (c) keeps a usable amplitude when made fast, and (d) works in both
    directions without closed loops;
  - every candidate examined so far fails at least one of (a)–(c), or lacks the evidence to decide it.

## 5. Bookkeeping

- Inputs: repository files and PRs; three papers read for P3 (arXiv PDFs in the scratchpad, not committed):
  Jacobson 1001.4823v3 [66fd0b22…], Blas & Sibiryakov 1110.2195v2 [d1a83782…] (§2, App. A), Bhattacharyya et al.
  1512.04899v2 [b8343d3b…] (§§II, IV).
- Computation: `scripts/goal_khronometric.py` (~20 s), results `results/goal_round9/khronometric.json`, tests
  `tests/test_goal_khronometric.py`.
- Corrections in this revision (owner review): the Jacobson citation (1003.4987 → 1001.4823); W3 survivors are
  candidates only; the Ellis and negative-mass templates are separated; the ~10³⁸ M☉ figure is scoped to round 6's
  model; a low S4 needs a realistic QM prediction; the "only route" remark is withdrawn.
- State files: CHANGELOG entry and a proposed DECISIONS record (D-TBD), with no change to TASKS priorities.
