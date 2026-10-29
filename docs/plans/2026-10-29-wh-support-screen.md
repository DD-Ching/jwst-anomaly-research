# Wormhole support mechanisms and finite signal loads: feasibility screen (2026-10-29)

Owner-requested bounded round (≤ 90 min). **No solver, evolution or scan was run**, and D-111's emitter was not
built. #592's results are not re-run or re-classified. This note asks whether a published wormhole support mechanism
with a better physical basis than #592's assumed ghost field is worth testing for its capacity to carry a
finite-energy signal. It is a literature screen plus a few order-of-magnitude estimates from published formulas
(`scripts/wh_support_scales.py` → `results/wh_support_screen/scales.json`, with sources, assumptions and units).
Nothing here is a communication result. No independent review was run in this round.

**Labels.** [lit] published, read in the paper body (section or equation given); [abstract] read in the abstract
only; [der] our arithmetic from published formulas (in `scales.json`); [assumption]; [inference] our reading, not
stated by the authors. Two research subagents extracted the papers from their arXiv sources; every quotation below
was re-checked against the extracted text, except the items marked [abstract].

## 0. What the repository already answers (cited, not redone)

The goal rounds are on the unmerged stacked branch `claude/goal-round10` (PRs #152–#166, `docs/goal_roundN/README.md`).

| Question | Answer already recorded | Where |
|---|---|---|
| Does a published traversable wormhole beat the exterior light time between two points of one universe? | No published GR + semiclassical construction checked does. Maldacena–Milekhin–Popov (MMP) 4D wormholes are "long"; Gao–Jafferis–Wall (GJW) paths are chronal; Fu–Grado-White–Marolf (FGM) have transit d + logs, and the generalized second law forbids wormholes from being the fastest causal curves; Maldacena–Milekhin 2020: "always longer" | round 1 §1.4; round 5 §0; round 10 |
| Which theorems frame this? | Gao–Wald (a geometric NEC statement), Friedman–Schleich–Witt topological censorship (ANEC), Graham–Olum's self-consistent achronal ANEC (a conjecture; with it, no wormholes between different asymptotically flat regions) | round 1 §1.4; round 6 |
| Does GJW need anything outside the spacetime? | Yes: a non-local coupling between the two boundaries | round 5 §0 |
| Can quantum negative energy be thick and long-lasting? | Quantum energy inequalities (Ford–Roman; Fewster–Smith) restrict free fields; non-minimal coupling has no state-independent bound (Fewster–Osterbrink); the one interacting example obeys ANEC (Olum–Graham); vacuum polarisation of realistic field content supports only Planck-scale walls. Analysed for Krasnikov tubes, not for wormholes | rounds 5–7 |
| Classical support without exotic matter (Einstein–Dirac–Maxwell)? | Closed: not solutions, not traversable, or only partial wormholes | round 10 |
| Do rotation or charge stabilise phantom-supported wormholes? | Rotation is not an established stabiliser; charged wormholes are long-lived near extremality but still unstable, with an extra polar ℓ = 2 instability; the 5D result of X2 S6a is unchanged and not applied to 4D | X2 S6 screen; X2 S6a |
| What do ghost-supported evolutions with signals show? | Collapse on the 0.846 r/c e-folding scale; test rays can return; finite-energy two-way propagation in one configuration (#592); a reply chosen after receipt needs an idealised device (D-111) | S0-2026-10-16 §1; #592; D-111 |
| Is the quantum-processor "wormhole" a channel? | No: a quantum simulation, not a spacetime channel | round 10; S0-2026-10-15 C |
| How much signal can a quantum-supported wormhole carry? | **Not analysed.** #592 §1 only names the AdS information-transfer bounds as "a different setting" | this note |

So the open part is the last row: for supports built from established quantum field theory, are the support, a
finite-energy signal and its back-reaction handled in one model, and what limits the signal?

## 1. Three published support methods

### M1. Holographic double-trace coupling (the GJW family)

Papers: Maldacena, Stanford & Yang, arXiv:1704.05333 (MSY); Freivogel, Galante, Nikolakopoulou & Rotundo,
arXiv:1907.13140v2 (FGNR); Maldacena & Qi, arXiv:1804.00491 (MQ); Lu, Yang & Zheng, arXiv:2603.26051 [abstract].

1. **Source.** An explicit interaction between the two boundary systems of an eternal AdS black hole,
   e^{igV} with V = (1/K) Σ_j O_L^j O_R^j [lit: MSY (2.2); MQ (2.11)], makes the negative null energy. Established
   within its framework: semiclassical gravity with quantum fields. **Assumed:** asymptotically AdS spacetime with
   a holographic dual (not the asymptotics of our universe), and the coupling itself, which MSY's protocol can
   replace by measuring and sending K classical bits [lit: MSY (2.3)].
2. **One model?** Yes, at the probe level. The coupling, the signal and the signal's gravitational back-reaction
   enter one eikonal formula: the signal's momentum shifts the coupled pair, "This is the total effect of the
   backreaction in this approximation" [lit: MSY (2.14)–(2.15)]. MQ add bulk energy to the Schwarzian-plus-matter
   system [lit: MQ (4.41)]. Beyond the probe regime FGNR write: "we do not have a persuasive bulk argument
   explaining why the information that can be sent is bounded by these probe regime calculations" [lit: §4].
3. **Compatibility conditions.**
   - Number of quanta: "N_send ∼ p₊^total/p₊^each ≤ g" with g ≲ K, the number of coupled operators [lit: MSY
     (2.21), (2.25)], so at most about one signal quantum per coupled operator; total momentum −p₊^total G_N e^t ≲ 1
     [lit: (2.24)].
   - Window: the wormhole is open for a proper time that is "indeed Planckian" [lit: FGNR (2.27)]; with K fields,
     N ≲ K r_h/ℓ ≲ S_BH [lit: FGNR (3.23)–(3.26)]; "For black holes with horizons of order the AdS radius,
     information cannot be reliably sent through the wormhole" [lit: FGNR abstract].
   - Energy: "If we add a large amount of energy in the bulk, then we could remove the possibility of having bound
     states" [lit: MQ §4.1, at (4.41)]; then only black-hole-like configurations remain [lit: MQ §4]. That the
     threshold is about the ground-state binding energy is our [inference].
   - Capacity: the quantum channel capacity of the GJW protocol is governed by the time derivative of an
     out-of-time-ordered correlator and bounded above by the Einstein-gravity limit [abstract: 2603.26051].
4. **Ends.** Two boundaries joined by the coupling: "This does not violate any causality constraint since we are
   adding direct interactions between the two sides" [lit: MQ §6]. The wormhole depends on a pre-existing
   non-local coupling or classical channel and cannot beat it (rounds 1 and 5).
5. **New against the repository.** Quantitative load limits (≲ g ≤ K quanta, ≲ S_BH, a Planckian window) and a
   2026 capacity paper. They do not change "needs a coupling" or "no shortcut".

### M2. Casimir energy of massless charged fermions on magnetic field lines (MMP 2018; Maldacena–Milekhin 2020)

Papers: Maldacena, Milekhin & Popov, arXiv:1807.04726 (MMP); Maldacena & Milekhin, arXiv:2008.06618 (MM20).

1. **Source.** About q two-dimensional lowest-Landau-level modes of massless charged fermions run along closed
   magnetic field lines that pass through the throat and return through the outside; their negative Casimir
   energy holds the throat open [lit: MMP §5.2, (5.25)–(5.26)]. Established ingredients: Einstein–Maxwell gravity,
   charged fermions, Landau levels, the Casimir effect. **Assumed:** a pair of near-extremal magnetically charged
   black holes with integer charge q ≫ 1, and massless charged fermions. With Standard Model fields alone the
   mouths must be closer than the electroweak scale ("say 1/TeV"), with r_e much smaller [lit: MMP §8], i.e.
   d < 2 × 10⁻¹⁹ m [der]. MM20 reach human scale only with a Randall–Sundrum dark sector (two-dimensional central
   charge c₂ ≈ 7 × 10⁷²), which they call "science fiction" [lit: MM20 §1, (3.27)]. Neither paper gives a
   formation route: "we did not describe a simple procedure" [lit: MMP §9.3]; "We have not given any plausible
   mechanism for their formation" [lit: MM20 §4].
2. **One model?** Support and geometry: yes, solved together semiclassically [lit: MMP §5.3.2]. **The signal is a
   probe only:** "the wormhole can be safely explored by sending low energy waves" [lit: MMP §1]; MM20 use an energy
   budget, a tidal estimate and a boost argument [lit: §3]. No paper found in this search evolves a finite-energy
   signal with its back-reaction through such a wormhole. New since round 1 (2026 preprints): on a fixed MMP
   background, low-frequency scalars are mostly reflected or trapped on light-crossing times, while charged
   massless fermions pass "with essentially unit probability", and probes must keep r_e ω ≪ 1, "otherwise, we would
   destroy it" [lit: Freivogel, Fumagalli & Tomašević, arXiv:2606.12528, abstract and §2]; first-order linear
   perturbations find no purely exponential growth in the sectors treated [abstract: Sadhukhan, arXiv:2609.33511,
   2610.09847].
3. **Compatibility conditions.**
   - Signal energy: sending "an energy that is larger than the binding energy of the throat (5.31)" makes a
     near-extremal black hole [lit: MMP §5.4]. With (2.3) and (5.31), E_max r_e/(ħc) ≈ g²/(256π) ≈ (1.1–1.6) × 10⁻⁴
     [der; g = 0.30–0.36 assumed]: one signal may carry about 10⁻⁴ of the energy of a single quantum whose
     wavelength is the mouth radius.
   - Information: the binding energy is "larger than the gap by a factor of q" [lit: MMP §5.4; gap ∝ 1/ℓ, (5.21)],
     so at most about q gap-sized quanta per transit [inference]; MM20: "∼ c₂ qubits per ℓ" [lit: §2].
   - Lifetime and range: the orbiting mouths radiate, T_lifetime ∝ d³/(q² l_p²) [lit: MMP (6.54)]; the solution
     needs q^{5/3} l_p ≪ d ≪ q^{5/2} l_p [lit: (6.52), (5.49)].
   - MM20's example [lit: (3.26)–(3.27), (B.32); der]: r_e ≈ 1.5 × 10⁷ m, ℓ ≈ 3.1 × 10³ ly, binding energy
     ≈ 5.2 × 10⁹ kg c², each mouth ≈ 2.0 × 10³⁴ kg (≈ 1.0 × 10⁴ M☉). The outside transit time is πℓ ≈ 9.8 × 10³ yr
     while the mouths must be more than 0.5 ly apart, so at that separation the wormhole takes ≈ 2.0 × 10⁴ times
     the exterior light time.
4. **Ends.** One universe; the mouths orbit each other; the coupling between them "is generated automatically by
   the exchange of massless fields in the bulk" [lit: MMP §1], so no external non-local channel is needed. Transit
   πℓ > d (round 1).
5. **New against the repository.** Round 1 recorded only "long". New here: the per-signal energy cap
   (≈ 10⁻⁴ ħc/r_e), the information cap (≈ q or c₂ quanta per ℓ), the Standard Model size limit (sub-attometre), the
   probe dependence of transmission, and first linear-stability preprints. None changes "no shortcut".

**Also checked (same class, not counted as a fourth method): FGM, arXiv:1908.03273.** Hartle–Hawking fluctuations
of N compact cosmic strings make two charged black holes in one asymptotic region traversable at first order
[lit: §2]. The traversability "is thus extremely fragile, as an exponentially small positive-energy perturbation
will negate this effect", which "includes possible perturbations associated with any signal one might attempt to
send" [lit: §4, fn 6]. The authors leave that back-reaction undone; it "will place fundamental limits on the amount
of information that can be transmitted" [lit: fn 6]. Transit d + logs (round 5); the holes merge on a time ∼ d^{3/2}
[lit: §1]; footnote 8 notes that the logarithmic Shapiro delay makes rigorous 4D bounds on transit times hard.

### M3. Casimir-sourced and vacuum-polarisation-supported static throats in 4D

Papers: Garattini, arXiv:1907.03623 (EPJC 79, 951); Butcher, arXiv:1405.1283 (PRD 90, 024019); Hochberg, Popov &
Sushkov, gr-qc/9701064 (PRL 78, 2050; HPS).

1. **Source.** Garattini: the parallel-plate Casimir energy density with the plate gap promoted to the radial
   coordinate [lit: (4)] and an assumed equation of state p_r = ωρ; no plates or fields are in the model
   [assumption]. Butcher: the renormalised stress tensor of a massless, conformally coupled scalar on a fixed
   infinite R × S² throat [lit: (59)], with a length a₀ that "can only be determined by experiment". HPS: a
   self-consistent one-loop scalar on a static wormhole joining "two asymptotically spatially flat regions" [lit].
2. **One model?** Static (Garattini, HPS) or quasi-static on a fixed background (Butcher); no finite-energy signal.
   Butcher traces test null rays only: "It remains to be seen whether the light pulse can actually escape" [lit].
   Citers 2020–2026 found only test-field quasinormal modes and a first-order, not self-consistent, back-reaction
   study [abstract: 2302.09141, 2610.06893, 2603.11724]; a finite-energy signal with back-reaction: not found in
   this search.
3. **Compatibility conditions.** Garattini: traversable only if ω = r₀²/r₁² [lit: (29)]; ω = 3 gives
   r₀ ≈ 1.02 l_P ≈ 1.6 × 10⁻³⁵ m [lit: (41)–(42), fn 6; der], "traversable in principle but not in practice", and
   "with a Planckian wormhole nothing can traverse it" [lit: fn 9]. HPS: throat r(0) ≈ 0.02 l_P in the numerical
   solution and ≈ 67 l_P in an analytic example [lit], i.e. 3 × 10⁻³⁷ to 1.1 × 10⁻³³ m [der]. Butcher: it is
   "impossible to solve the Einstein equations for the static wormhole … using only its Casimir energy-momentum
   … and ordinary matter"; a long throat collapses slowly [lit: (99)].
4. **Ends.** Two asymptotic regions, no exterior path to compare; topological censorship and the Graham–Olum
   conjecture apply (rounds 1 and 6).
5. **New against the repository.** Nothing for signal loads. It agrees with round 6's estimate that vacuum
   polarisation of realistic field content supports only Planck-scale structures.

## 2. Summary of the three

| Method | Basis of the support | Support + finite signal + back-reaction in one model | Load limit | Ends |
|---|---|---|---|---|
| M1 GJW family | explicit boundary coupling, AdS + holography | yes, probe / eikonal and Schwarzian level | ≲ g ≤ K quanta; Planckian window; ≲ S_BH | needs the coupling (a channel); no shortcut |
| M2 MMP family | Casimir energy of charged massless fermions; magnetically charged black-hole pair assumed | support yes; signal probe only | ≈ 10⁻⁴ ħc/r_e per signal [der]; ≈ q or c₂ quanta per ℓ | one universe, no extra channel; slower than the exterior path (≈ 2 × 10⁴ × in MM20's example [der]) |
| M3 Casimir static | effective equation of state or one-loop scalar | static only | Planck-scale throats | two regions; no exterior path |

M2 has the best physical basis of the three: the support follows from established field theory once its objects
exist. But those objects (magnetically charged black-hole pairs with q ≫ 1, massless charged fields beyond the
Standard Model at any size above 10⁻¹⁹ m, or a dark sector) are assumed, no formation route is given, and the channel
is slower than light outside.

## 3. Recommendation: zero

**No test is recommended.** For the communication goal, a support mechanism matters only if it gives a channel
that is faster than the exterior light path between the same events, or one that links regions that have no other
path, without needing a channel that already exists. M1 needs a coupling that is itself a channel. M2's published
solutions are slower than the exterior path (πℓ > d; round 1, MMP §5.5); FGM's transit is d plus logarithmic
terms, approaching the exterior light time without being shown to beat it (round 5; FGM §4, fn 8). M3 is
Planck-scale and static.

The one open physics question these papers name is real: the back-reaction of a finite signal on an FGM- or
MMP-type throat, and the information limit it sets [lit: FGM fn 6]; no such computation was found in this search.
Its two possible answers are (a) the throat carries a message below its binding energy and stays open, or (b) any
useful signal closes it. Neither changes our next step: neither gives a faster or independent channel, and neither
supplies the assumed objects. It therefore fails the owner's criterion of two answers that change our next step,
and it is not recommended.

**Missing physical conditions** (any one would reopen this line):
1. A published, self-consistent construction whose transit beats the exterior light time between the same events.
   None of the constructions checked does; the generalized-second-law argument forbids it in D ≥ 5, and FGM's
   footnote 8 notes that rigorous 4D bounds are hard.
2. A support with an established physical basis that links otherwise disconnected regions without a direct
   coupling. Topological censorship and the Graham–Olum conjecture stand against it; #592's two-universe wormhole
   rests on the assumed ghost field.
3. A formation route for MMP-type objects, or evidence for magnetically charged black holes or for new massless
   charged fields.

A limited search that did not find something is not evidence that it does not exist.

## 4. Plain answers

1. **What we know.** No published wormhole built from established quantum field theory that this or earlier
   rounds checked beats light outside.
   - The best grounded one (MMP) needs a pair of magnetically charged black holes. With Standard Model fields it is
     smaller than 10⁻¹⁹ m.
   - It can carry a signal only below about 10⁻⁴ of the energy of a quantum as long as its mouth. In its
     published solutions the transit is slower than the exterior path.
   - The human-size version needs a hypothetical dark sector and takes about 10⁴ years from the outside for mouths
     half a light year apart.
   - The GJW family needs a direct coupling between the two sides, which is already a channel, and carries at most
     about one quantum per coupled field.
   - Casimir-sourced static throats are Planck-sized.
   - #592's ghost-supported model stays the only finite-energy, back-reacted two-way calculation we have; no paper
     found does the same for a quantum-supported wormhole.
2. **Premise without a basis.** A wormhole support that gives a channel faster than, or independent of, ordinary
   light paths. It is not just unbuilt: the published arguments (the generalized second law, topological censorship
   and the achronal ANEC conjecture) stand against it. Next to it, the objects every grounded construction assumes
   (magnetically charged black-hole pairs, cosmic strings, a dark sector) have no formation route.
3. **One new question, or zero?** Zero. The only open question in these papers, the back-reaction of a finite
   signal on an MMP- or FGM-type throat, has two possible answers that both leave our next step unchanged.
4. **What success would have bought.** Even a complete load calculation would only say how many bits a channel
   slower than light, built from assumed objects, can carry. It would not move any link of the communication chain
   (sender, faster or independent path, receiver) closer.
