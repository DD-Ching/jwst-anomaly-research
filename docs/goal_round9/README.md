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
| Decision-changing result | Any S4 ≤ 6 where QM predicts a violation, or B ∼ C-only setting dependence (round 3 trigger list), would reopen the line toward a channel design. A clean QM-valued violation excludes the no-signalling class for all finite v, but says nothing about the signalling class that matters for communication [rec; round 3 §3]. |
| Minimal verification | Data from such an experiment. **Not available.** A single-device LC4 dataset would test only Bancal's Experiment-v3 hypothesis, and "would not test the signalling models that matter for communication" [rec; round 3 §3]; no such dataset was found either. |
| Gaps to communication | A model that fixes the signalling branch and amplitude per direction. The preferred frame. The ≥ L/2 latency guarantee (10 yr at 20 ly) and ~10¹³ trials per bit in the unknown-frame case [rec; round 2]. |

### P2. Does negative or exotic mass exist in nature (inverted microlensing, W3)?

| Question | Answer |
|---|---|
| Missing evidence | Any confirmed event. Fourteen MOA-II fields are null, with Γ₉₅ ≈ 2.1 × 10⁻⁸–1.2 × 10⁻⁷ per star per year [obs; TASKS.md, D-068]. There is also no abundance model to turn the limit into a decision [rec; goal round 1 §5.3]. |
| Prediction different from standard physics | The source flux drops to zero between two caustic spikes (Ellis or negative-mass lens) [rec; `docs/exotic_lensing.md`]. Standard lensing never dims a source below baseline. |
| Decision-changing result | A W3 event that survives every vetting test would establish that exotic stress-energy exists in nature. That is the missing input behind "no self-consistent source found" for shortcut geometries (rounds 6–7) [derived]. Further nulls change nothing for the communication goal without an abundance model. |
| Minimal verification | MOA-II fields gb3–gb6, gb9, gb10, gb14, plus gb13 after the owner's decision on #145. Data **are** available, but this is an existing, scheduled search (gb3 is claimed in #160). Selecting it here would duplicate it or change its schedule. **Not selected.** |
| Gaps to communication | Even a confirmed exotic lens would not supply a controllable source, a way to assemble ~10³⁸ M☉ of negative energy (round 6), or a stable shortcut. |

### P3. Can a preferred-frame theory with an unbounded scalar speed (khronometric / Hořava IR) avoid round 8's trade-off?

| Question | Answer |
|---|---|
| Changed premise | Local Lorentz invariance of the gravity sector: a preferred foliation (round 5 §5 "not covered" class). |
| How verified physics is recovered | The same way as in Einstein–æther. Small couplings restore GR, PPN and GW170817 [lit: 1802.04303; 1311.7144 treats Hořava jointly]. |
| First experimental constraint | PPN α1 and α2, binary pulsars, and GW170817, as in round 8. |
| Missing model condition | Khronometric theory is Einstein–æther with a hypersurface-orthogonal æther, which removes the vector mode [lit, not re-read this round: Jacobson, arXiv:1003.4987]. Its spin-0 sector then has the same structure as round 8, where the matter coupling falls as 1/c_S² [derived inference, not computed]. Escaping the trade-off would need the infinite-speed (instantaneous) limit, which is strongly coupled; no controlled theory of it was found [rec; round 8 §1, "EFT range: open"]. |
| Prediction different from standard physics | The same breathing-mode precursor as round 8 [rec; round 8 §4], unless a strong-coupling completion changes it. Nothing new without that completion. |
| Decision-changing result | A controlled strong-coupling or UV-complete khronometric theory in which matter sources the instantaneous mode with an amplitude that does not vanish. |
| Minimal verification | A theory input that does not exist. A round-8-style recomputation would only reproduce the same trade-off. |
| Gaps to communication | All of round 8's gaps, plus the strong-coupling regime. |

## 3. Comparison and choice

| | Different prediction | Data or model available | Light test | Changes a communication decision |
|---|---|---|---|---|
| P1 v-causal | yes (S4 ≤ 6; B ∼ C signalling) | **no** (no four-party spacelike data) | — | yes, if positive |
| P2 exotic lenses | yes (flux to zero) | yes, but already scheduled in another program | the existing pipeline | only if positive; nulls need an abundance model |
| P3 khronometric | not beyond round 8 without a strong-coupling completion | **no** (no controlled theory) | — | only with that completion |

**Result: there is not currently enough basis to select a next executable unit.** No candidate has all four of: a
prediction that differs, data in hand, a light test, and a decision-changing outcome. P2's data are in hand, but
P2 is an existing, scheduled search, and the owner asked not to duplicate it or change its schedule.

**New inputs that would make a unit executable:**
1. **For P1:** a public dataset or platform with four spacelike-separated parties, no postselection, and the round-3
   §6a capability. Or any reported S4 ≤ 6. Or a published model that fixes the v-causal signalling branch and
   amplitude.
2. **For P2:** an unexplained W3 survivor from the running MOA-II fields (after #145 is decided). Or a quantitative
   abundance or rate model for exotic lenses, so that nulls become decisions.
3. **For P3 (or the other round-5 uncovered classes: Gödel/Gott geometries, non-local quantum gravity, tachyons):**
   a concrete published claim with equations, a matter coupling, and a prediction of controllable early reading.

## 4. Plain-language summary

- **New evidence?** No [derived]. This round compared the existing records; it did not run a new measurement or
  calculation.
- **Judgement changed:** "what to do next". No candidate on hand justifies a new physics unit. The useful inputs are
  external (an experiment, a confirmed lens event, or a new theory with a controllable coupling). The ongoing W3
  searches and the owner's decision on #145 are the only routes that could produce new evidence soon.
- **Still missing for communication** [derived]:
  - a mechanism that (a) a sender can drive with ordinary matter or a lab device, (b) produces a readable change at
    a receiver 20–100 ly away before light, (c) keeps a usable amplitude when made fast, and (d) works in both
    directions without closed loops;
  - every candidate examined so far fails at least one of (a)–(c), or lacks the evidence to decide it.

## 5. Bookkeeping

- Inputs: repository files and PRs only. No downloads, no computation beyond reading.
- State files: CHANGELOG entry and a proposed DECISIONS record (D-TBD), with no change to TASKS priorities.
