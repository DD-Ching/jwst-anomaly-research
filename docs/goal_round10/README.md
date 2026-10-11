# Goal round 10: external search, de-duplication and one minimal verification

Owner brief (2026-10-11). Search for external, reusable results (GitHub, papers, corrections, public data) that
could close a specific gap in the communication chain: a sender acting on ordinary matter, a propagating
influence, a receiver reading it, and the physical channel itself, over 20–100 ly with low latency in both
directions; later, transport of matter. Compare up to three opportunities. Then run one cheap verification
directly. No second-order, non-linear or strong-coupling work and no large scans were run. Paused lines
(round 9) were not reopened. The observational searches and their schedule are unchanged.

**Labels:** [lit] stated in a cited source; [obs] measured data; [derived] computed or reasoned here; [hyp] our
hypothesis; ASSUMPTION marks a stated assumption.

## 1. What was searched

Three read-only survey subagents ran in parallel: tools and code, theory since 2023, and experiments and public
data. They read arXiv abstracts, PDFs, GitHub pages and PyPI metadata. Facts used below were re-checked where they
matter; anything only second-hand is marked.

**De-duplication against the record.** Rounds 1–9, the System A/B hypotheses, MOA-II W3, S1–S3, D1/E1 and the
neutrino program already cover the following: Lentz, Bobrick–Martire, Fell–Heisenberg, Santiago–Schuster–Visser
(generic warp drives violate the NEC), Krasnikov, Scharnhorst, Einstein–æther and khronometric preferred-frame
gravity, v-causal hidden influences, beyond-Horndeski wormholes, and the NIST 2015 data **as a timing bound on V**
(round 1 §5). Only what is new is listed below. Forks, ports and AI restatements of one paper count as one
source.

| Item | Kind | What it adds here | Verdict |
|---|---|---|---|
| WarpFactory (github.com/NerdsWithAttitudes/WarpFactory; MIT; MATLAB; no tags, last push 2024-05-24); paper arXiv:2404.03095 (CQG) | tool | Metric → stress-energy, plus energy conditions on **sampled** observers (100 directions × 10 speeds). A Python port exists (JohnsterID/WarpFactory, MIT; unreleased; it claims "within 2 %" of the paper's figures, which I did not run). | reusable tool; not a physics candidate |
| **warpax** 1.5.1 (github.com/anindex/warpax; PyPI; MIT; JAX; 2026-09-24); paper arXiv:2602.18023 (v6; README says provisionally accepted by CQG) | tool + derivation | Energy conditions over **all** causal directions (LMI / S-lemma, interval certificates). Eulerian-frame checks miss 54–74 % of WEC/DEC violations for Van den Broeck and Rodal (Rb = 1, σ = 8, v_s = 0.5) [lit]. Lemma 4: the violation scales exactly as v_s² for a unit lapse, flat slice and linear shift [lit]. | best tool to reuse **if** a metric needs checking; no new route |
| Fuchs et al., "Constant velocity physical warp drive", arXiv:2405.02709 (CQG 41, 095013) | numerical derivation | Positive-energy shell: R1 = 10 m, R2 = 20 m, M = 4.49 × 10²⁷ kg, β_warp = 0.02. Energy conditions reported satisfied on sampled observers [lit]. warpax finds the NEC margin positive on the sampled slice [lit]. Barzegar, Buchert & Vigneron (arXiv:2602.16495) argue the construction is not an Einstein-equation solution (TOV equations not actually solved) [lit, second-hand summary]. | **subluminal**; disputed; no latency gain |
| Maldacena–Milekhin humanly traversable wormholes, arXiv:2008.06618 | theory | "Time through the wormhole is always longer than through the outside" [lit]. | no latency gain by construction |
| Einstein–Dirac–Maxwell wormholes (2010.07317) and corrections (2108.13361; Kain 2305.11217, 2308.00049; Weinbaum 2607.28738) | theory | Not solutions, not traversable, or only partial wormholes [lit]. | closed |
| Garattini–Zatrimaylov 2502.13153; Rodal 2512.18008; Le 2606.22531 ("steering without exotic matter") | theory | Energy conditions hold only up to divergence terms; NEC violation reduced but not removed; a radiating photon rocket [lit]. | no positive-energy superluminal solution |
| Jafferis et al. "wormhole teleportation" (Nature 612, 51; 2022); critique arXiv:2302.07897 | quantum simulation | 9 qubits, sparsified SYK [lit]. | **not a spacetime channel**; excluded as a physics candidate |
| Rybotycki et al., arXiv:2409.11348 (no-signalling violation on IBM qubits) | measurement | The authors cannot exclude technical imperfections or subluminal communication [lit]. | example of an artefact class to rule out |
| Broz et al., arXiv:2206.12976 (PRL 130, 200201) | measurement | Causal nonlinear QM (Kaplan–Rajendran) parameter < 5.4 × 10⁻¹² [lit]. | bounds a model built not to signal |
| LVK GWTC-4 tests of GR (2603.19019/20); LHAASO GRB 221009A Lorentz-invariance violation (2402.06009) | measurement | Tighter null results [lit]. | no route |
| NIST 2015 loophole-free Bell data (s3 `nist-belltestdata`; public, attribution) | **data** | Whole-run counts N(o_A, o_B \| x, y) per run and pulse grouping. They allow a direct test of whether one site's free setting choice moves the other site's click rate. Round 1 used these data only for timing. | **chosen for §4** |

The tools mature enough to reuse are warpax, OGRePy 2.0.0 (symbolic; MIT) and WarpFactory. None is needed this
round.

## 2. Three opportunities

**O1. No-signalling in loophole-free Bell data (NIST 2015; Delft as a later cross-check).**
- *Done by others:* measurements. NIST computed two-proportion no-signalling z values in diagnostic spreadsheets
  (supplement Table S-III). Gill (arXiv:2209.00702) used no-signalling deviations to sharpen Bell statistics.
  Adenier–Khrennikov (1606.00784) and Bednorz et al. (1906.05503) argued that apparent signalling appears in
  some data.
- *What it adds:* the only public data in which a **free, spacelike-separated choice** is the input and a
  **remote detector** is the output. That is exactly the send/receive primitive of an entanglement channel.
  Round 1 used the data only for a timing bound on V.
- *Unproven assumption:* that no-signalling holds at the level these data can test, for this apparatus.
- *Chain link:* sender and receiver of a quantum-correlation channel.
- *Next check:*
  - continue if a family-wise significant shift with the same sign appears in ≥ 2 runs;
  - stop if not, recording the bound.

**O2. Positive-energy "physical warp" shells (Fuchs et al.) with observer-robust checking (warpax).**
- *Done by others:* a numerical metric construction and a re-analysis with a tool.
- *What it adds:* the only claimed positive-energy moving shell with a payload region. Rounds 5–6 covered only
  superluminal solitons and tubes.
- *Unproven assumption:* that the matter shell solves Einstein's equations with a physical equation of state
  (disputed by 2602.16495). Acceleration is unsolved.
- *Chain link:* the **matter-transport** goal only. It is **subluminal**, so it does not address latency.
- *Next check:* re-solve the shell's TOV equations independently and test whether the published density and
  pressure profiles satisfy them, at a cost of about 1 h of sympy/scipy.
  - continue if they do and β_warp can rise without violating the DEC under warpax's all-observer check;
  - stop if the profiles fail TOV.
  - This does not move the communication goal.

**O3. Signalling through nonlinear quantum mechanics.**
- *Done by others:* measurement (Broz et al. 2023) and theory (Bielińska–Eckstein–Horodecki 2412.20854: some
  global nonlinear dynamics permit statistical superluminal transfer).
- *What it adds:* the textbook mechanism by which a quantum channel could signal (Gisin, Polchinski).
- *Unproven assumption:* that any nonlinearity exists. The measured bound concerns a model that is causal by
  design.
- *Chain link:* sender and receiver.
- *Next check:* map a signalling-capable nonlinear model's parameter onto a Bell-type marginal shift. The O1 bound
  then constrains it directly.
  - continue only if a published signalling model predicts a shift above the O1 bound;
  - stop otherwise.

**O1 was chosen.** It is the only one with public data, a direct link to the send/receive step, and a check that
fits this round. O2 does not touch latency. O3 needs the O1 bound first.

## 3. Minimal verification: no-signalling in the NIST 2015 loophole-free Bell data

**Pre-registration:** `docs/goal_round10/preregistration.md`. It fixes the hypotheses, the ordinary explanations,
the primary family, the Bonferroni threshold and the continue/stop rules. It was written before the inference
was run.

**Data** [obs]:
- **Source:** NIST belltestdata, `Peter's Data Calculations/Diagnostic Statistics/*.xlsx`. Six runs (19_45, 02_54,
  03_43, ClassicalXOR 1–3) of Shalm et al., PRL 115, 250402 (2015). Accessed 2026-10-11; sha256 prefixes in the
  results.
- **Counts:** whole-run counts for the five pulse groupings, about 4.4 × 10⁷ trials per setting pair per run.
- **Check:** the trial totals implied by the 16 counts match the sheet's totals in all 30 tabs.

**Method** [derived]:
- `scripts/goal_nist_nosignal.py`: binomial likelihood-ratio (G) tests and Fisher exact tests, written
  independently of the spreadsheet formulas. NIST's pooled z is read back only for comparison.
- **Primary family:** grouping 2345678 (the largest reported as spacelike) × 6 runs × 4 comparisons, so 24 tests,
  α = 0.05/24.
- Comparisons: Alice's click rate given her setting a or a′ as Bob's setting varies, and the converse for Bob.

**Results** [derived; `results/goal_round10/nist_nosignal.json`]:
- **Reproduction.** The independent pooled z equals NIST's in every one of 120 cells (max |Δz| = 3 × 10⁻⁴). The G and
  Fisher p-values agree.
- **Primary family:**
  - One test is below the Bonferroni threshold: run XOR_2, Alice | a, z = 3.13, p = 0.0017 (relative shift 3.7 %).
  - A second test is below 0.05 (03_43, B | b, z = 2.30).
  - Fisher's combined p = 0.39; KS test of the 24 p-values against uniform, p = 0.095.
  - **Neither pre-registered rule applies (owner review, 2026-10-11).**
    - "Continue" needs the shift in ≥ 2 runs. No other run shows one in that comparison; the largest is XOR_1,
      z = 1.58.
    - "Stop (null)" needs no family-wise significant shift. XOR_2 is one.
    - The rules did not cover a single, unreplicated hit. The result is therefore **undetermined under the
      pre-registered rules**. It is not a pre-defined null, and the rules are not changed after the fact.
    - The combined p-value is not used to erase the single-run deviation.
- **Indicative pooled limits** (fixed-effect inverse-variance over the six runs; Wald 95 % interval on the
  setting-induced shift |δ| of the remote click probability). These are **not rigorous capacity bounds**:
  - they assume one common δ for all runs, which XOR_2 does not obviously share (heterogeneity p = 0.12);
  - the interval is a large-sample Wald interval;
  - the "mutual information" is computed for a uniform binary input switching a Bernoulli output. It is not
    maximised over input distributions, and it ignores memory and timing.
  - XOR_2 alone has δ̂ = 1.2 × 10⁻⁵ (3.7 % of p; Wald 95 % upper ≈ 1.9 × 10⁻⁵), which is above the pooled limit.
    The pooled numbers therefore do not bound a run-specific effect.

| Comparison | click prob. p | δ̂ (z) | \|δ\|₉₅ | \|δ\|₉₅ / p | indicative MI, uniform input (bits/trial) | indicative trials per bit |
|---|---|---|---|---|---|---|
| Alice \| a (Bob's choice) | 3.06 × 10⁻⁴ | 3.2 × 10⁻⁶ (2.09) | 6.2 × 10⁻⁶ | 2.0 % | 2.3 × 10⁻⁸ | 4.4 × 10⁷ |
| Alice \| a′ | 8.94 × 10⁻⁴ | (0.78) | 7.2 × 10⁻⁶ | 0.8 % | 1.0 × 10⁻⁸ | 9.6 × 10⁷ |
| Bob \| b (Alice's choice) | 3.06 × 10⁻⁴ | (1.08) | 4.7 × 10⁻⁶ | 1.5 % | 1.3 × 10⁻⁸ | 7.8 × 10⁷ |
| Bob \| b′ | 9.73 × 10⁻⁴ | (−1.25) | 8.8 × 10⁻⁶ | 0.9 % | 1.4 × 10⁻⁸ | 7.0 × 10⁷ |

  - At 99,100 trials/s these correspond to about 2 × 10⁻³ bits/s, under the assumptions above (indicative only).
  - The four combined z values were not pre-registered as a family. With 4 comparisons, Alice | a's z = 2.09 is not
    significant.
- **Post hoc** (not pre-registered; reported for completeness):
  - Within XOR_2, the Alice | a shift appears in both the first 85 % of the run (z = 2.41) and the last 15 %
    (z = 2.29). This comes from NIST's time-segmented spreadsheet for that run (00_25 ClassicalRNGXOR_2), whose
    totals agree with the diagnostic counts to within a few counts.
  - Without XOR_2, the combined Alice | a is z = 0.90.
  - XOR_2 differs from the rest by 2.5σ; Cochran's heterogeneity test across all six runs gives p = 0.12.
  - So the hit is either a fluctuation (family-wise p ≈ 0.04) or a run-specific systematic. It does not behave like
    a law that holds across runs.

**Conclusion** [derived]:
- **Undetermined under the pre-registered rules** (see above). Five of the six runs show no setting-dependent
  shift. One run (XOR_2) shows an unreplicated 3.1σ deviation that is not explained.
- It is **not** evidence of superluminal signalling: it does not replicate, and a run-specific systematic has not
  been excluded. It is **not** erased either: the combined p-value does not remove a single-run deviation.
- The pooled effect is ≤ 0.8–2.0 % of the click probability. That is an indicative fixed-effect limit, not a
  capacity bound.
- This agrees with NIST's own diagnostics, which it reproduces exactly, and with quantum mechanics.
- **It is not new physics evidence.**
- **What it adds:**
  - an independent, pre-registered reproduction;
  - an indicative (not rigorous) limit on the send/receive step of entanglement-based signalling on real
    spacelike data;
  - one localised, unexplained 3.1σ run-level shift (XOR_2). It is flagged for an optional raw-time-tag crosstalk
    check, not claimed.

**Scope.** This is 184.9 m, one apparatus and 2015 data. It constrains signalling models only to the extent that
they predict a marginal shift in this geometry. Mapping it to 20–100 ly needs a model.

## 4. Our hypothesis (labelled [hyp]; not tested this round)

**H-R10 (frame-gated signalling).**
- *Statement:* if a signalling v-causal influence exists, it acts only when the two setting events are close to
  simultaneous in its preferred frame. Then a whole-run average dilutes any shift, and the shift is modulated by
  sidereal time.
- *Minimal operational definition.* Let Δt_pref(φ) be the preferred-frame time difference between the setting
  events at sidereal phase φ, computed from the NIST geometry (round 1) and an assumed frame (CMB, ASSUMPTION).
  Then δ(φ) = P_A(+ | a, b) − P_A(+ | a, b′) is predicted to be non-zero only where |Δt_pref(φ)| < L/V.
- *Distinguishing prediction:* δ depends on φ. QM and every ordinary drift predict no dependence on φ in the
  difference between Bob's settings.
- *Cheapest falsification:* compute δ in 1-minute bins from the NIST processed hdf5 files and fit a sidereal
  template.
  - Cost: about 1–2.5 GB per run, streamed and deleted, about 2–3 h for six runs.
  - Limitation: the runs cover only a few sidereal hours on 18–19 Sept 2015, and the Delft data (another site and
    epoch) would add coverage.
  - **Stop** if no bin-level shift exceeds the per-bin bound scaled from §3.
- *Status:* not started. Proposed for the owner's decision.

## 5. Collaboration (draft only; not sent)

- **To:** NIST Bell test team (belltest@nist.gov, the contact given in the data description).
- **We would need:**
  - the analysis YAML (click window delays and radii);
  - the UTC start time of each run to the second;
  - the fixed Alice–Bob detection offset;
  - any record of Pockels-cell/RNG crosstalk tests for run 00_25 (ClassicalRNGXOR_2).
- **We offer:** a reproducible, pre-registered no-signalling analysis (`scripts/goal_nist_nosignal.py`) that
  reproduces their diagnostics to |Δz| ≤ 3 × 10⁻⁴, the capacity bounds above, and the flagged XOR_2 shift, for
  their assessment.
- **Draft text:** "We reproduced the no-signalling diagnostics of your 2015 data (six runs, five pulse groupings)
  with an independent likelihood-ratio analysis, pre-registered before inference. The results agree with your
  spreadsheets. One run (00_25, ClassicalRNGXOR_2) shows a 3.1σ dependence of Alice's click rate on Bob's setting
  at the 7-pulse grouping, present in both time segments. It is not significant across runs. Do you know of a
  hardware cause (e.g. RNG/Pockels-cell crosstalk) in that run, and could you share the UTC run start times and
  the analysis YAML so we can test it at the time-tag level?"

## 6. One-page answer

- **New physics evidence?** No.
- **External results worth continuing?**
  - Yes as tools: warpax (all-observer energy-condition checks) is the one to reuse if a metric must be checked.
  - Yes as data: the NIST 2015 loophole-free Bell data, and Delft once its downloads return.
  - No external physics result found gives a latency advantage over light at 20–100 ly:
    - "physical" warp shells are subluminal and disputed;
    - known-matter wormholes fail or are slower than the outside path;
    - "wormhole teleportation" is a quantum simulation;
    - new data tightened Lorentz-invariance bounds.
- **What was verified and the result.**
  - The NIST 2015 no-signalling diagnostics were reproduced exactly and extended with a pre-registered test.
  - Outcome: **undetermined** under the pre-registered rules. There is one unreplicated 3.1σ run-level deviation
    (XOR_2); the other five runs are null.
  - Indicative pooled limit: |δ| ≤ 0.8–2.0 % of the remote click probability (fixed-effect, uniform-input MI
    ~10⁻⁸ bits/trial). It is not a rigorous capacity bound.
- **Which gap moved?** Only a little. The send/receive step now has an indicative measured limit on real
  spacelike data, plus one open run-level deviation. It closes nothing by itself.
- **Next concrete action, cost and stop condition** (owner decides; nothing started):
  - (a) Optional: a time-tag crosstalk check of run XOR_2. One ~1.4 GB hdf5 file, about 1 h. Stop if Alice's
    shift tracks Bob's Pockels-cell drive timing or vanishes under a stricter click window.
  - (b) H-R10 sidereal test: **kept as a hypothesis only**. Before any multi-GB download, these must be shown
    from metadata alone:
    - sidereal coverage of the six runs;
    - the per-bin detectable shift given the counts per minute;
    - controls for ordinary drift and device periodicities (fridge, laser mode-locking, PPS jumps).
    Stop if coverage or power is insufficient.
  - (c) The Delft 2015 cross-check once 4TU downloads return: about 1 h. Stop if null.
  - Second-order, non-linear and strong-coupling work stays unstarted.
