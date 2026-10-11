# Goal round 1: low-latency two-way contact and fast return over 20–100 light-years

Owner brief, 2026-10-11. The goal is a real, verifiable and controllable mechanism that lets two ends 20–100
light-years apart talk with low latency in both directions, and later lets people travel between them quickly,
without the family on Earth or the traveller waiting decades. Warp drives, wormholes, new propagation rules and new
physics are candidates. None is assumed feasible.

Labels: **[paper]** read in the full text, with the section, equation or figure given; **[derived]** worked out in
this round; **[estimate]** order of magnitude; **[ASSUMPTION]**; **[hypothesis]**. Nothing here is evidence for new
physics.

Round started 2026-10-11 03:43 UTC from `main` d2d9484. Code: `scripts/goal_vcausal_bounds.py`,
`scripts/goal_wormhole_timeline.py`. Tests: `tests/test_goal_*.py`. Outputs: `results/goal_round1/`.

## 1. Corrections to the previous assessment (chat report, 2026-10-11 ~01:30 UTC)

Each item is kept, modified or withdrawn.

### 1.1 "Raising the hidden-influence speed bound closes the communication path": withdrawn

A bipartite test that keeps violating a Bell inequality while the two measurements are simultaneous, to precision
ρ, in a candidate frame excludes only v-causal models with influence speed V below the bound, for the frames the run
covers (Salart et al. 2008 Eqs. 2–3; Yin et al. 2013 Eq. 2; Cocciaro et al. 2018 Eq. 5) [paper]. Such a test:
- does not exclude faster models;
- does not test frames outside its geometry (Cocciaro: χ < π/10 or χ > 9π/10; large β weakens every bound);
- says nothing about whether a model allows controllable signalling.

| Observation | What it supports | What it excludes | What it does not show |
|---|---|---|---|
| Bipartite: Bell violation persists at frame simultaneity | — | v-causal models with V < bound in the covered frames | anything about V > bound; controllability |
| Bipartite: violation vanishes reproducibly at a fixed sidereal phase, not at a solar phase, on two baselines | a v-causal model with V near the bound in one frame (to be vetted against drifts) | QM as the complete description of that run | a usable channel |
| Multipartite (Bancal configuration R, B and C timed at random, App. E v1): S > 7, as QM predicts (7.2014) | QM | **no-signalling** v-causal models at any finite V in that frame (Bancal's claim covers all finite v), given λ does not know the measurement timings (App. E) | signalling v-causal models: they reproduce S > 7 and must be caught through their marginals |
| Multipartite: S ≤ 7 in R | a contradiction with QM, which must first be vetted as an artefact; or a model whose λ knows the timings | QM as the complete description | a channel |
| Multipartite: a three-party marginal depends on a remote setting (direct check) | controllable signalling at speed ≲ V | QM's no-signalling | a rate or practical interstellar geometry |

### 1.2 Bancal et al. 2012 (Nature Phys. 8, 867; arXiv:1110.3795): modified and completed

**Assumptions [paper]**
- A preferred frame, or a foliation, exists; influences propagate inside "v-cones" (Fig. 1, p.2–3).
- Settings are freely chosen, and screening-off holds (Fig. 2).
- The model reproduces the quantum predictions whenever each successive measurement lies inside the v-cone of the
  earlier ones (p.3).
- No-signalling is defined as the setting-independence of every tripartite marginal (Eqs. 1, 3).
- For tests: λ must not depend on future measurement timings (App. E).

**Arrangement [paper].** Four parties with A < D < (B ∼ C). B and C lie outside each other's v-cones. Fig. 3 gives
explicit positions and times in units of r = v/c, and the configuration is called R (Fig. 4, App. C).
- Lemma 1 (App. B, Eq. 4): all inputs and outputs are binary, and S ≤ 7 when the BC correlations conditioned on AD
  are local and the distribution is no-signalling.
- A specific 4-qubit state gives S = 7.2014; the maximum is 7.3481. S contains only ABD and ACD terms.

**Logic [paper].**
- The model reproduces the quantum ABD marginal when B < C and the ACD marginal when C < B. B cannot know when C
  is measured, so these marginals must stay the same when B ∼ C (App. D).
- That gives S > 7 in R, which violates no-signalling.
- The theorem is therefore a disjunction: the model either deviates from QM in an ordered configuration, or it
  signals. Equivalently, finite-v models are inconsistent with QM plus no-signalling.

**What a receiver gets [paper].** The predicted signal is a dependence of the BCD marginal on A's setting, or of the
ABC marginal on D's setting.
- Reading it needs B's and C's outcomes brought to the receiver (D′ or A′) by ordinary light-speed communication, so
  it is not a purely local marginal. The theorem does not fix the direction.
- All tests in App. E (v1–v3) pool data from every station after all runs are completed.
- Barnea et al. 2013 (PRA 88, 022123; arXiv:1304.1812) give a tripartite version with a tiny violation, −2.00015
  against the bound −2 (Sec. III B), and call it "too fragile to noise to test" (Sec. IV).

**Status.** No experiment has tested the multipartite prediction. Li, Hu, Deng & Scarani, arXiv:2608.05271 (theory,
Aug 2026), p.2: "no experiment has reported even a proof-of-principle demonstration yet". They propose cluster-state
witnesses (4 qubits: bound 6, quantum 4 + 2√2) whose certificates the audit could not check.

**Interstellar use [derived; open].** The Fig. 3 geometry scales with the station separations. Whether B, C and the
receiver can sit on one planet while A sits 20 ly away is not established by the paper. Until that is shown, the
"receiver reads locally" step is an open item.

**Previous report.** The statement "a finite v allows signalling, so the model itself is a communication channel" is
withdrawn and replaced by the disjunction above.

### 1.3 Carried wormhole, event by event: kept (no loop for one link), extended

See section 2.

### 1.4 Energy-condition theorems: modified

Each hypothesis is labelled G (geometric or causal), F (uses the field equations) or M (matter energy condition).

| Result | Statement (where) | Hypotheses | Not covered |
|---|---|---|---|
| Olum 1998, PRL 81, 3567 | "superluminal" (Condition 1, p.4: a local, necessary condition) + null generic ⇒ WEC violated on the path (p.4–6, Eqs. 3–8) | G: Condition 1, generic; F: WEC ⇒ R_ab K^a K^b ≥ 0; M: WEC (only its null projection is used); 4D | a wormhole longer than the ambient path (Condition 1 fails; Graham–Olum §IV.A); CTC regions (p.6); signals not bounded by g's light cones. The Casimir effect satisfies Condition 1 locally (Eq. 11), so the theorem does not by itself forbid superluminal travel |
| Gao & Wald 2000, CQG 17, 4999 | NEC (= R_ab k^a k^b ≥ 0, Eq. 2) + null generic + completeness ⇒ no shortcut through a compact region (Th. 1); AdS boundary signals are fastest (Th. 2) | purely geometric (R_kk), any dimension (§1) | non-generic (exact AdS, de Sitter) or incomplete spacetimes. A statement about **matter** needs the field equations, which the paper does not use |
| Friedman, Schleich & Witt 1993, PRL 71, 1486 | asymptotically flat, globally hyperbolic + ANEC on every inextendible null geodesic ⇒ topological censorship (Th. 1, p.3) | G + M (ANEC, no achronality restriction); F implicit; 4D | non-asymptotically-flat or non-globally-hyperbolic spacetimes. The semiclassical status is open: ANEC fails in generic curved spacetimes (p.9) |
| Graham & Olum 2007, PRD 76, 064001 | Condition 1 (a conjecture, §I): no self-consistent semiclassical solution violates ANEC on a complete achronal null geodesic. With it: no complete achronal null geodesics (Lemma 1), topological censorship for **simply connected** spacetimes (Th. 1), no compactly generated time machines (Th. 2–3) | F (self-consistent Einstein equations) + M; sub-Planckian curvature | **§IV.D: achronal ANEC does not rule out Olum-type superluminal travel** (it would need ANEC on half-geodesics, which the Casimir effect violates). A wormhole joining a region to itself with a throat longer than the exterior distance is not excluded (§IV.A) |
| Kontou & Olum 2015, PRD 92, 124009 | achronal ANEC proven for a free, massless, minimally coupled scalar in a Hadamard state, to first order in curvature, in a background obeying null convergence (Th. 1) | as stated | back-reaction, non-minimal coupling ("remains an open question", §VIII), higher spin, interacting fields |
| Hartman, Kundu & Tajdini 2017; Faulkner et al. 2016 | ANEC in flat-space QFT (from causality / relative entropy) | interacting d > 2 CFT-like theories (HKT); a formal continuum limit (Faulkner et al. §5.1) | curved spacetime and gravity (HKT §1); general achronal geodesics (Faulkner et al. §5.2) |
| Maldacena, Milekhin & Popov 2018, arXiv:1807.04726 | traversable 4D wormholes, constructed in a specific setup, are "long": going through takes longer than the ambient path. Short ones "are not allowed by the Einstein equations combined with the achronal average null energy condition" (§1; Fig. 7, §5.5) | Einstein equations + achronal ANEC | — |
| Gao, Jafferis & Wall 2016, arXiv:1608.05687 | null geodesics through the wormhole are chronal, so the wormhole gives no faster-than-light travel over long distances (§5) | — | — |

**Corrections to the previous report**
- **Modified.** "Every shortcut violates achronal ANEC, so Graham–Olum rule it out." This holds only for a
  configuration that needs a complete achronal null geodesic. Olum-type local superluminal travel is explicitly not
  excluded (Graham–Olum §IV.D), and simple connectedness is needed for their topological-censorship theorem.
- **Modified.** "Shortcuts need NEC violation." This is a statement about R_kk (geometry) in Gao–Wald. It becomes a
  statement about matter only through the Einstein equations. In modified gravity the geometric condition can fail
  without exotic matter, and each such theory then faces its own constraints (e.g. GW170817; not reviewed here).
- **Kept.** In GR with semiclassical matter, the published constructions of traversable wormholes give no shortcut,
  and making or holding a short throat open is unsolved.

### 1.5 Over-strong statements in the previous report

| Statement | Verdict | Reason |
|---|---|---|
| "The carried wormhole is the only configuration that meets the goal" | **withdrawn** | It was the only one found among the GR configurations the round looked at; other topologies, modified gravity and Krasnikov-type tubes were not exhausted. Reworded: "the one configuration examined that meets both conditions on paper" |
| "Bell (v-causal) is the only direction decidable in the lab" | **withdrawn** | The search was not systematic. Other laboratory tests (e.g. of non-minimally coupled ANEC violation, or analogue tests) were not assessed |
| "No remaining search changes any decision" | **modified** | True only in the absence of a model that predicts a rate. Example: the W3 MOA-II combined limit, Γ₉₅ ≈ 2.8 × 10⁻⁸–1.5 × 10⁻⁷ per star per year (docs/exotic_limits.md), would be decisive against any abundance model predicting more. None exists yet. The remaining fields would improve the limit by a factor I estimated at ≲ 2–3 without checking N_s, so that number is unverified |
| "Cocciaro 2018 reaches ≈ 5 × 10⁶ c" | **kept, now verified** | Fig. 7, β ≈ 10⁻³, χ = π/2. It is reproduced from the paper's own parameters in §3 |
| "Exotic mass ~ c²a/G ≈ 1.35 × 10²⁷ kg per metre of throat" | kept as [estimate] | dimensional estimate, not a literature value |
| "A v-causal model would itself be a channel" | **withdrawn** | §1.2 |

## 2. Carried wormhole: event accounting [derived; premise unsolved]

`scripts/goal_wormhole_timeline.py`. Model [ASSUMPTION]:
- Flat ambient spacetime and the Earth frame S. Mouth A stays on Earth; mouth B rides the ship.
- The throat is short, rigid and static inside, so B's proper time τ is identified with A's Earth time t = τ
  (the Morris–Thorne–Yurtsever idealisation). Traversal time is neglected.
- The ship coasts at speed β, or accelerates and decelerates at 1 g.

| Trip | Earth-frame arrival | Traveller τ | Shift s = t − τ | Loop margin D − s | Mother's elapsed time at first call (through the throat) | Mother's elapsed time at reunion via the throat (1 yr stay) | Reunion by ship instead (0.99c) |
|---|---|---|---|---|---|---|---|
| 20 ly, 0.99c | 20.2 yr | 2.85 yr | 17.35 yr | +2.65 yr | 2.85 yr | 3.85 yr | 41.4 yr |
| 20 ly, 0.999c | 20.0 | 0.90 | 19.12 | +0.88 | 0.90 | 1.90 | 41.2 |
| 20 ly, 1 g | 21.85 | 6.04 | 15.81 | +4.19 | 6.04 | 7.04 | 43.1 |
| 100 ly, 0.999c | 100.1 | 4.48 | 95.62 | +4.38 | 4.48 | 5.48 | 202.1 |
| 100 ly, 1 g | 101.9 | 9.02 | 92.90 | +7.10 | 9.02 | 10.02 | 203.9 |

**Events [derived]**
1. **Departure.** t = 0, τ = 0.
2. **Arrival.**
   - Earth frame: t = D/β.
   - Through the throat, the mother sees the arrival at her own time τ.
   - By ordinary light she sees it at D/β + D.
3. **First two-way call.** Earth event t ↔ destination event t + s, with near-zero round-trip latency. The mother's
   clock and the ship's clock read the same value (τ at arrival). Destination clocks synchronised with Earth in S
   read τ + s.
4. **Return through the throat** after a stay Δ: the traveller exits on Earth at t = τ + Δ, which is the mother's
   elapsed time.
5. **Loops.**
   - Earth → light → mouth B → throat → Earth comes back D − s later. Since s < D for every β < 1,
     (1 − √(1−β²))/β < 1, so there is never a loop. The tests check this on a 400-point β grid and for 1 g.
   - Adding any second short connection between the same two places with a different shift s₂ closes a loop of size
     |s − s₂| immediately. A slowly moved wormhole (s₂ ≈ 0) is one example.
   - The traveller flying home by ship creates no loop, but the mother then waits more than 40 years.

**What stays unsolved.**
- Making a traversable, short throat: §1.4, the open premise.
- Keeping it rigid and short during a relativistic trip.
- Transporting a mouth at γ ~ 7–22.
- The chronology-protection status of a second link or of a returning mouth (section 2.1). A single one-way link is
  not near s → 0.

### 2.1 Literature on the time-shift and loop conditions

**Read in full or in the relevant sections**
- Morris, Thorne & Yurtsever 1988 (PRL 61, 1446), Caltech scan.
- Morris & Thorne 1988 (Am. J. Phys. 56, 395).
- Kim & Thorne 1991 (PRD 43, 3929).
- Visser: hep-th/9202090, hep-th/9303023, gr-qc/9311026 and gr-qc/9702043.

**Not opened:** Frolov & Novikov 1990 and Hawking 1992. Hawking's position is cited as Kim & Thorne report it.

**Identification and time shift**
- **Identification by equal mouth proper time** [paper]:
  - Morris, Thorne & Yurtsever p. 1447 and Fig. 2 (the metric's t is each mouth's proper time). The idealisation
    needs a rigid throat with g_max·S ≪ 1 (peak acceleration times the distance through the throat).
  - Kim & Thorne Eq. 28.
- **Time shift after transport** [paper]: T = ∫(γ − 1) dt in the frame of the resting mouth (Visser
  hep-th/9202090 Eq. 2.4), i.e. T = t_E − τ, as in `Trip.shift_years`.

**Loop conditions**
- **One wormhole** [paper]: a closed timelike curve exists once ℓ < T, and the closed geodesic through the throat
  has length √(ℓ² − T²) (Visser hep-th/9202090 Sec. II.C, Eq. 2.6; gr-qc/9702043 Eq. 2). Here ℓ = D after arrival.
- **No loop for a one-way trip** [derived; also derived by the audit; in no source]:
  - Along a trajectory that never backtracks, √(1 − v²) ≥ 1 − v at every instant, so t − τ(t) ≤ x(t) during the
    trip and after it. Hence T < D strictly: a one-way trip never makes a CTC, and only approaches the chronology
    horizon as v → 1.
  - At constant v, D − T ≈ D/γ and the closed geodesic has length s ≈ D√(2/γ). For 20 ly at 0.999c that is about
    6 ly.
  - No source treats a one-way transport. Morris, Thorne & Yurtsever get CTCs only after the mouth turns round and
    comes back (pp. 1447–1448).
- **Two or more wormholes** [paper]:
  - Morris & Thorne Sec. III.I (p. 407): two wormholes in relative motion give backward time travel.
  - Visser gr-qc/9311026 Sec. VI (the "Roman configuration") gives the condition |T₁ + T₂| > ℓ₁→₂ + ℓ₂→₁, for the
    normal-space separations between the mouths of the two wormholes.
  - When both wormholes join the same two places, ℓ₁→₂ = ℓ₂→₁ ≈ 0. Any difference in shift then closes a loop, as
    `loop_margin_years(..., other_shift_years)` reports [derived from the paper's condition].
  - Rings of N wormholes can approach a CTC when no subset does (gr-qc/9702043 Eq. 2).

**Quantum back-reaction** [paper]
- The renormalised stress-energy diverges only as the closed geodesic length s → 0 (Kim & Thorne Sec. II; Visser
  hep-th/9303023 Eq. 16: ⟨T⟩ ~ Σ Δ^{1/2} t/(π² s⁴)).
- The effect needs a short closed spacelike geodesic, not a CTC.
- A single-wormhole time machine is disrupted at s ≈ Δ^{1/8}√(π l_P R) (gr-qc/9311026 Eq. 5).

**Inference for a one-way trip** [inference]: s ≈ 6 ly at 0.999c, so ⟨T⟩ ~ ħ/s⁴ is negligible and no chronology
horizon forms. The protection arguments apply only to configurations driven towards s → 0, e.g. a second link or a
return of the mouth.

## 3. Selected unit: can existing public data inform the v-causal channel question, and what would?

**Candidates after the corrections**

| Candidate | Falsifiable question | Observable | Data / feasibility this round | Decision value |
|---|---|---|---|---|
| A. v-causal influences (bipartite and multipartite) | Are finite-V influences present in some frame? Do they signal? | Bell violation at frame simultaneity; three-party marginals | Public NIST 2015 data (checked); published bounds | settles whether any public data bears on it, and fixes a lab target |
| B. Carried wormhole | Does one carried mouth close a causal loop? | none (theory) | done (§2) | bookkeeping only; the premise stays open |
| C. Modified-gravity escape from Gao–Wald | Is there a theory with R_kk < 0 and no exotic matter that survives GW170817 and lab tests? | c_T, gravitational-wave decay, solar-system bounds | literature-heavy; no data step this round | high if found, but not a minimal unit |
| D. More MOA-II / S1 / E1 / neutrino scans | rate limits | existing pipelines | owned by other sessions | depends on a rate model that does not exist |

A was chosen. It is the only candidate with a decidable data question in one round (does this dataset carry any
information?), a positive control (three published bounds) and a direct link to the communication goal.

**Results [derived from paper parameters]**

1. **Positive control.** The formula, with each paper's own ρ, window and β = 10⁻³, gives:

   | Experiment | Reproduced | Published |
   |---|---|---|
   | Salart 2008 (ρ = 5.4 × 10⁻⁶, T = 360 s) | 5.40 × 10⁴ c | 5.4 × 10⁴ c |
   | Yin 2013 (ρ = 6.84 × 10⁻⁶, T = 1800 s) | 1.38 × 10⁴ c | 1.38 × 10⁴ c |
   | Cocciaro 2018 (ρ = 1.83 × 10⁻⁷, δt = 0.494 s) | 4.97 × 10⁶ c | ≈ 5 × 10⁶ c |

   `tests/test_goal_vcausal_bounds.py` asserts agreement within 2 %. Notes from the independent review:
   - The orientation factors cos 5.8° (Salart) and A = cos 18° (Cocciaro Eq. 8) are dropped. With them the values
     are 5.42 × 10⁴ and 5.00 × 10⁶.
   - Salart's 5.4 × 10⁴ c is at χ = 90°. His frame-independent claim is "> 10⁴ c", with the minima near
     χ ≈ 0°/180°.
   - By Cocciaro's own argument (Eq. 6), the Salart and Yin windows should be doubled. That gives 3.2 × 10⁴ c and
     7.2 × 10³ c (Cocciaro Fig. 7 caption). The published values were reproduced as published.

2. **Public NIST 2015 data has no discriminating power** for V. Shalm et al. 2015 (PRL 115, 250402) give:
   - an Alice–Bob distance of 184.9 m (Fig. 4a);
   - a timing uncertainty of 3.7 ns (Table I margins, 63.5–25.7 ± 3.7 ns);
   - a violation of only p = 2.5 × 10⁻³ (1 pulse) to 5.9 × 10⁻⁹ (5 pulses) over the whole dataset, which is about
     2 days of runs.

   The repository is public. It holds two-party time tags at 78.125 ps with a GPS PPS channel (100 ns jumps
   possible) and run start times to the minute (`File_Folder_Descriptions.pdf`; NIST licence: free use with
   attribution).

   The Alice–Bob detection offset t_AB is a fixed fibre-and-geometry delay that the paper does not give. Both photons
   come from one pump pulse, so aggregating slots does not change it. Assuming a zero offset, which can only flatter
   the data, the envelopes are:

   | Case | Bound |
   |---|---|
   | timing uncertainty 3.7 ns | V ≤ 1.7 × 10² c |
   | time-tag floor 78 ps | V ≤ 7.9 × 10³ c |
   | tag floor + rotation over one ~30 min run (Table I: 1.77 × 10⁸ trials at 99,100 trials/s gives the p = 5.9 × 10⁻⁹ result; ASSUMPTION: no shorter window is usable; CMB frame) | V ≤ 4.8 × 10³ c |

   - Every envelope is ≥ 600× below the published 5 × 10⁶ c.
   - The data is bipartite, so it cannot test the Bancal prediction at all.
   - A first draft of this argument used the slot spread as ρ. The independent review showed that is the wrong
     quantity. The conclusion is unchanged.
   - No scan was run, by design. Nothing was downloaded except the two description PDFs (0.45 MB).

3. **Requirements for the goal** (frame = CMB, β = 1.23 × 10⁻³ [ASSUMPTION]; the budget is split equally between the
   lab offset and the rotation residual):

   | Target | V/c | ρ | window T | timing σ_t at 1.2 / 18 / 100 km |
   |---|---|---|---|---|
   | 20 ly, 1 s one-way | 6.3 × 10⁸ | ≤ 7.9 × 10⁻¹⁰ | ≤ 18 ms | 3.2 / 48 / 264 fs |
   | 20 ly, 0.3 s | 2.1 × 10⁹ | ≤ 2.4 × 10⁻¹⁰ | ≤ 5.3 ms | 1.0 / 14 / 79 fs |
   | 100 ly, 1 s | 3.2 × 10⁹ | ≤ 1.6 × 10⁻¹⁰ | ≤ 3.5 ms | 0.6 / 9.5 / 53 fs |
   | 100 ly, 0.3 s | 1.05 × 10¹⁰ | ≤ 4.8 × 10⁻¹¹ | ≤ 1.1 ms | 0.2 / 2.9 / 16 fs |

   - The gap to the best published bound is ×130 to ×2,100 in V.
   - "Window" T is the half-width convention (Yin Eq. 5). If T is the acquisition time per measurement, Cocciaro's
     Eq. 6 halves it: 8.8, 2.6, 1.8 and 0.53 ms.
   - Statistics [estimate]: var(S) ≈ 16(1 − E²)/N with E = S/4, so a 5σ separation of S = 2.6 from 2 needs about
     650 coincidences per window. Windows at one sidereal phase recur every sidereal day, so they can be stacked.
   - The multipartite test is harder. The Bancal gap (7.2014 − 7)/7 = 2.9 % is 0.069 of the CHSH relative gap, so it
     needs ~200× more events at equal noise [estimate]. It also needs near-ideal 4-qubit fidelity (Barnea Sec. IV;
     Li et al. 2026).

**Limitations**
- The formula is the published bipartite one. A dedicated derivation of the multipartite timing tolerance for
  configuration R was not done.
- The CMB frame is one candidate frame; others need the multi-orientation scan (Bancal App. E: eight planar
  configurations).
- The NIST bounds are optimistic envelopes, with a zero offset assumed. Any real offset only lowers them.
- The literature search for re-analyses of loophole-free data was rate-limited and found none. That is not proof
  that none exists.

## 4. Round bookkeeping

- **Time.** 2026-10-11 03:43 → 04:05 UTC, from `main` d2d9484.
- **Compute.**
  - New calculations: well under 1 CPU-minute. Both scripts run in under a second.
  - The full offline test suite takes about 1 minute.
  - Four literature and review subagents ran for about 18 minutes of wall time in total, all of it reading
    papers.
  - The 2-hour budget was not needed.
- **Data.**
  - NIST belltestdata description PDFs (0.45 MB; S3 bucket `nist-belltestdata`).
  - arXiv and journal PDFs, read by the subagents into the session scratchpad. Not committed.
  - No dataset was scanned.
- **Done.**
  - §1 corrections, with paper citations.
  - §2 event accounting, with tests.
  - §3 calibration, data-adequacy check and requirement table, with tests.
  - An independent numerical and logical review. It found that two rows of the §1.1 table were inverted and that
    the NIST argument used the wrong quantity; both are fixed.
- **Not done.**
  - A dedicated timing-tolerance derivation for Bancal configuration R.
  - Noise and visibility modelling of the 4-party and cluster-state tests.
  - A scan of other loophole-free datasets (Delft, Vienna, Munich). Their two-party geometry excludes the
    multipartite test in the same way, and their bipartite reach was not computed.
  - Candidate C (modified-gravity escape from Gao–Wald).
- **Unchanged.** The schedules of other work, the selection thresholds, permissions, and the MOA-II, S1, E1 and
  neutrino work.

## 5. What this round supports next

1. **The decisive v-causal experiment is multipartite and has never been done** (Li et al. 2026, p.2). The cheapest
   new information is theoretical:
   - derive the timing tolerance of configuration R;
   - find the noise threshold of the 4-qubit (or cluster-state) witness at a stated visibility;
   - check whether R can be laid out with B, C and the receiver co-located, which an interstellar channel needs.

   Only then does an experimental design make sense. Two-party public Bell data cannot help.
2. **Carried wormhole.** The accounting is closed for one link. The open premise is a short throat, which
   Maldacena–Milekhin–Popov say "is not allowed by the Einstein equations combined with the achronal average null
   energy condition". The question worth effort is candidate C: is there a gravity theory with R_kk < 0 on a
   throat-crossing geodesic, no exotic matter, and consistency with GW170817 and the solar-system bounds? That needs
   a literature unit before any data.
3. **Existing searches** (MOA-II, S1, E1, neutrinos) remain bound-setting. Their decision value depends on the first
   abundance or rate model that predicts a number in their reach. No such model exists yet.
