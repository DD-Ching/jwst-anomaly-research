# Goal round 4: does the beyond-Horndeski wormhole give a stable, usable channel?

Owner brief 2026-10-11. The question: can a concrete modified-gravity wormhole model provide a stable channel
that carries messages and gives both ends a time advantage?

Labels used below:
- **[paper]**: arXiv:2212.05969v1, Mironov, Rubakov & Volkova, PRD 107, 104061 (2023), and arXiv:2404.06297v2,
  Mironov & Volkova (v2 of 2025-01-17, IJMPA ICNFP2022 issue). Both were read in full; the formulas were read from
  page images.
- **[derived]**: this round.
- **[abstract only]**, **[ASSUMPTION]**, **UNKNOWN**.

Nothing here is a new-physics result. Reproducing the papers is not a discovery, and a failure of this one model
says nothing about other mechanisms.

Code: `scripts/goal_bh_wormhole.py`. Tests: `tests/test_goal_bh_wormhole.py`. Output:
`results/goal_round4/bh_wormhole.json` (about 18 s).

## 0. What the project already had, and what this round adds

| Earlier work | What it covered | Reused here |
|---|---|---|
| Goal round 1 (#152) §1.4, §2 | Energy-condition theorems; carried-wormhole bookkeeping under an *assumed* short rigid throat; "candidate C: a modified-gravity escape from Gao–Wald" left open | The bookkeeping and the theorem table. The open premise is the one tested here. |
| D-047, W2/W3 (Ellis lenses), docs/exotic_limits.md | Observational searches for Ellis-type (1/b² deflection) and negative-mass lenses | Only to check whether they apply to this model (§5). |
| Hypothesis rounds A/B (D-069, D-075, D-076), D1 (D-073), E1 (D-074) | New-physics hypotheses, distance self-consistency, event networks | None treats a modified-gravity wormhole, so nothing is reused. |

**The unsolved assumption this round tests.** Round 1 assumed a short traversable throat can be held open without
NEC-violating matter. The question here is whether one specific modified-gravity solution does that *stably*.

## 1. The fixed model [paper]

**Theory.** Quadratic beyond Horndeski, action (1) of 2212.05969, with the functions F(π, X), G4(π, X) and
F4(π, X).

**Solution** (Sec. 4, Eqs. 70, 72, 73, 80, 81 and 83):
- metric: A = 1, B = 1 + sech(r/τ), J = τ ln(1 + 2 cosh(r/τ)), for r ∈ (−∞, ∞);
- scalar: π = tanh(r/τ) − 1;
- parameters: τ = 100, w = 2, u = 0.

**Lagrangian.** It is reconstructed from the solution ("reverse engineering"):
- H = G = 1, Ξ = 0 and Γ = 0;
- H − 2F4B²π′⁴ = 1 − 2 sech(r/τ);
- c²_s2 = cosh/(1 + cosh).

The functions are given only as plots and asymptotics (Eqs. 84–86).

**Units [ASSUMPTION].** G4 → ½, so these are reduced Planck units. The throat radius is J(0) = τ ln 3 = 110
l_Pl ≈ 9 × 10⁻³³ m.

**Consequence [derived].** Because Ξ = Γ = 0 and H = G = 1, every quantity in the 2404.06297 stability set is
fixed by the metric. That includes F, which follows from their Eq. (71): F = J²[B′J′/(2J) + B(J′²/J² + J″/J)].
So the complete set can be checked without the unpublished Lagrangian functions.

**Literature status** (INSPIRE citations and an arXiv search, 2026-10-11):
- **No erratum or second version of 2212.05969 exists.**
- **2404.06297 v2, by the same authors, states that this solution "does not satisfy the complete set of stability
  conditions and, hence, cannot claim to be fully stable"** (Conclusion). It does not say which condition fails.
- Of the five papers citing 2212.05969, none constructs a fully stable version. These include no-go and
  speed-of-gravity papers by the same group: 2411.05416, 2408.01480 [abstract only].
- The 2609.00664 obstruction concerns black holes with static scalar hair, not wormholes [abstract only].

## 2. Step by step: what holds, what fails, what is unknown

| Step | Status | Evidence |
|---|---|---|
| **Solution exists** | **Holds, as published** (not re-derived here) | Reconstruction in 2212.05969 §4.1–4.2. The background equations were not re-checked: the Lagrangian functions are not tabulated. |
| No ghosts, no radial or angular gradient instabilities (high ℓ) | **Holds** [paper, reproduced] | F_min = 1.0000001 > 1. K11, det K, G11, det G, M11^(ℓ²), det M^(ℓ²) > 0 on the grid, as in their Fig. 2. 2P1 − F > 0 and N > 0. |
| Odd-parity tachyons, condition (6d) | **Holds** [derived] | V + 6H/J² > 0 at every grid point for ℓ = 2–50. |
| Even-parity tachyons and low-ℓ angular modes, (39b) and (39c) | **Sufficient condition fails** [derived, new]. Whether a real tachyon exists is **UNKNOWN**. | (39b) M22 > 0 fails for ℓ = 2 at \|r\| < 0.455τ. (39c) fails near the throat and the Θ = 0 crossing, and also for r ≳ 6τ when ℓ ≲ 20. The set is sufficient, not necessary (2404.06297 §1, §4). |
| (Sub)luminality | **Fails** [derived, quantitative; agrees with 2212.05969 §4.2] | One even-parity mode has angular speed² up to **1.83** at \|r\| ≈ 1.4τ. It stays above 1 at every grid point and tends to 1 far away. The other mode lies in 0.16–0.82. Radial speeds are ≤ 1. |
| Nonlinear stability | **UNKNOWN** | Not addressed by any source found. |
| Validity of the effective theory | **UNKNOWN**, with two warning signs | The throat is 110 Planck lengths. The reconstructed functions grow as (−π)^(−7/2) and (−π)^(−11/2) as π → 0 at large r (Eq. 86), so the couplings of perturbations there are not shown to be weak. |
| Formation or control | **UNKNOWN** | No mechanism is given. The solution is static and eternal. |

**Scale.** Every condition above is unchanged when τ is rescaled: the τ = 10⁴ check gives the same speeds to 10⁻⁹.
So making the throat macroscopic does not cure the failures. That would in any case be a different reconstructed
Lagrangian.

## 3. The new calculation [derived; positive controls first]

**Unsolved question chosen.** Does the published solution satisfy the *complete* stability set of 2404.06297 that
followed it, and if not, which part fails? This decides whether the model can even be called linearly stable,
which is a prerequisite for any channel.

**Method.**
- The model's choices reduce every coefficient to the metric (§1).
- sympy builds Eqs. (5), (22), (24), (27), (33)–(37), (41)–(43) and (129)–(133) symbolically.
- numpy evaluates them on 4,000 points in −10 ≤ r/τ ≤ 10, avoiding the point r = 0, where J′ = 0.
- An mpmath (60-digit) evaluation at eight points agrees with the float64 values.

**Positive controls.**
- F > 1, matching their statement below Eq. 72.
- All six Fig. 2 quantities are positive.
- The ℓ → ∞ limit of the full-ℓ Eq. (33) matches Appendix B Eqs. (132) and (133) to 7 × 10⁻⁷ and 2 × 10⁻⁵
  (relative) at ℓ = 10⁴. M11 matches Eq. (41) except near its zeros.
- The tensor-like angular eigenvalue tends to 1 in the flat far zone, which confirms the normalisation of the speeds.

**Inconsistencies found in 2404.06297.**
- det K from Appendix B (129)–(130) is exactly 4 × Eq. (22).
- det M^(ℓ²) from (41), (132) and (133) is not proportional to Eq. (42).

The speeds use the Appendix B components, because those agree with Eq. (33) and give a luminal tensor mode far
away. This choice changes no conclusion: Eqs. (22) and (42) also give a positive det K and det M^(ℓ²).

**Results.**
1. **Superluminality is confirmed and quantified.** The maximum angular speed² is 1.83, at |r| ≈ 1.4τ.
2. **The low-ℓ sufficient conditions fail.** This is new: no earlier source evaluated them on this solution.
   - (39b) fails for ℓ = 2 around the throat.
   - (39c) fails over 47 % of the inner zone (|r| < 3τ) for ℓ = 2. The failing fraction falls with ℓ, to 2 % at
     ℓ = 50, where only the neighbourhood of the Θ = 0 crossing fails.
   - (39c) also fails in a far zone, r ≳ 6τ, for ℓ ≲ 20. That zone is nearly flat (B − 1 ≲ 5 × 10⁻³). The failure
     there may come from the criterion treating v and v′ as independent rather than from a real mode. **UNKNOWN.**
3. **What it settles.** The model is *not shown* to be linearly stable, because the only published complete criterion
   fails, and it *is* superluminal in one sector. **It does not show that a tachyon exists.** Only the low-ℓ
   eigenvalue problem can decide that (§7).

## 4. Clocks, signal time and the channel [derived]

**The two ends.** The ends are the two asymptotically flat regions r → +∞ and r → −∞. **The solution connects two
separate asymptotic regions**: there is no exterior path between them, so "exterior light time" and "time
advantage" are undefined in this solution.
- A shortcut inside one universe would need the two far regions glued into one exterior. No such solution of this
  theory was found in the literature searched; it is UNKNOWN.

**Clock comparison.**
- A = 1 everywhere, so static clocks at both ends tick at the Killing-time rate. Their ratio is exactly 1.
- They can be synchronised by Einstein exchange through the throat.
- **There is no relative time shift between the ends.** Round 1's carried-mouth time shift needs one mouth to move.
  That is a non-static configuration outside this solution, and its stability is UNKNOWN.

**Signal time.**
- Radial light from r = −R to r = +R takes Δt = 2R − τI, with I = ∫(1 − 1/√(1 + sech x))dx = 1.1066.
- At τ = 100 the throat saves 111 l_Pl ≈ 3 × 10⁻⁴¹ s against the coordinate length 2R. The light time equals the
  proper radial distance, since A = 1.
- The one-way time is about (proper distance)/c, and a round trip twice that.

**Deployment, continued use and round trips.**
- Deployment, i.e. creating a throat or moving its ends into place, is not described: UNKNOWN.
- Continued use needs the throat to survive finite-amplitude signals. That requires nonlinear stability: UNKNOWN.

**Message transmission at τ = 100 [estimate].** Waves much longer than the 9 × 10⁻³³ m throat are mostly reflected.
A carrier that fits would have a quantum energy of order 10¹⁷ GeV, outside the regime where the effective theory is
known to hold. For a macroscopic member of the family (larger τ) this objection disappears classically, but §2 and
§3 apply unchanged.

**Answer to the brief's channel question.** This model gives no channel with a time advantage for both ends:
- it links two separate asymptotic regions with equal clock rates and no exterior path to beat;
- its stability is not established;
- it is superluminal in one sector.

## 5. Observational limits: what applies to this model, background and parameters

- **Lensing.** The ADM mass is zero and B − 1 = sech(r/τ) decays exponentially. Deflection is exponentially small for
  impact parameters ≫ τ.
  - The project's Ellis-lens limits (W2/W3, D-047) assume a 1/b² deflection law, so they **do not constrain this
    model**.
  - At τ = 100 l_Pl, no astronomical observation applies.
- **GW170817 (c_T = c).** That limit constrains tensor speeds on the *cosmological* scalar background of a theory.
  - Here the Lagrangian is reconstructed only along this solution (π ∈ (−2, 0)), so its cosmological background is
    undefined.
  - The radial tensor speed² on *this* background, G/F = 1/F, is below 1 near the throat (down to 0.68). That is
    not excluded by GW170817, which says nothing about a local static background.
  - Whether a cosmological completion exists, and lies in the c_T = c class, is **UNKNOWN**. 2408.01480 [abstract
    only] relates the two for time-dependent backgrounds.
- **Lab and solar-system tests.** The modification is confined to r ≲ few τ around a Planck-size object. None
  applies.

**No prediction found** that an observation or experiment could distinguish at the published parameters.

## 6. Answers

**Where it holds.** The solution exists, as published. It is free of ghosts and of high-ℓ gradient
instabilities, and of odd-parity tachyons.

**Where it fails or is unknown.**
- It fails the complete sufficient stability set: low-ℓ (39b) and (39c), and subluminality, with c_a² up to 1.83.
- Real tachyons, nonlinear stability, validity of the effective theory, formation and control are all UNKNOWN.
- It provides no shortcut and no time shift: separate asymptotic regions, equal clock rates.

**What the new calculation settled.**
- The authors' unexplained statement that the solution fails their complete set is confirmed.
- The failing conditions are located (low-ℓ even parity near the throat, superluminal angular modes) and quantified.
- Two internal inconsistencies of 2404.06297 are flagged.

**Distinguishable prediction.** None at the published parameters.

**Next minimal verification step.** Solve the ℓ = 2 even-parity eigenvalue problem, ω² from K⁻¹(G, Q, M), in a
formulation that is regular at the Θ = 0 crossing. That decides whether a real tachyon exists:
- a negative ω² would make the model linearly unstable (a FAIL);
- none would leave superluminality and the missing nonlinear and EFT analysis as the remaining obstacles.

**Decision proposed.** Do not extend this line toward communication. Even a stable version of this solution has no
exterior path and no time shift. A concrete candidate worth a next round would need an *intra-universe* traversable
solution of a stated theory with a computed stability spectrum. If none can be named from the literature, pause the
line. This proposal is for the owner to decide.

## 7. Not done

- The low-ℓ eigenmode calculation (§6).
- Re-deriving the background equations, which needs the reconstructed functions.
- The 2018 predecessor (1812.07022). 2404.06297 also says it fails the complete set.
- Wave transmission through the throat (only an estimate is given above).

## 8. Bookkeeping

- **Time.** 2026-10-11 06:15 → about 06:50 UTC.
- **Compute.** About 18 s for the script. A 60-digit cross-check takes about 1 min.
- **Data.** Only the two papers (arXiv PDFs, sha256 `2a8ea802…` and `22eb223e…`, kept in the scratchpad and not
  committed), the INSPIRE citation lists and the arXiv search results.
