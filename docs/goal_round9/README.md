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
| Missing evidence | Linear order: settled for c13 = 0 (§2b); for c13 ≠ 0 the linear response has an instantaneous part (§2c), whose usability is not settled. The non-linear elliptic part is not settled (§2b). |
| Prediction different from standard physics | At linear order with ordinary matter: none beyond round 8 for c13 = 0 (§2b). For c13 ≠ 0: a transient tidal term with no retardation in the æther frame, ∝ c13²/c14 (§2c; model-internal). |
| Decision-changing result | A demonstration that **local, conserved matter**, with the prior state and the outer boundary conditions held fixed, drives the slice-wise elliptic (instantaneous) part to a non-zero tidal response at a distant receiver. |
| Minimal verification | §2b and §2c: done at linear order (c13 = 0 and c13 ≠ 0); the non-linear step is costed in §2b and not started. |
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
  for every c_i. **Scope (owner review 2): the conclusions below hold for c13 = 0 only; c13 ≠ 0 is §2c.**
  - The vector sector differs only by EA's c_V pole, which is proportional to c13²; khronometric has no pole there.
  - With c13 = 0 (GW170817) the two theories give the **identical** linear response.
- **Consequences (c13 = 0 only):**
  - At linear order, the BS-type instantaneous term does not appear in the tidal response to ordinary matter.
  - Round 8's results carry over exactly to khronometric theory at linear order: the causal c_S cone and the
    coupling κ ≈ c14/4 ∝ 1/c_S².
  - The two cases are **the same** under: linear order, matter coupled to the metric only, c13 = 0, tidal readout.
  - They **differ** only in the vector sector when c13 ≠ 0 (EA has the c_V pole, khronometric does not). GW170817
    bounds |c13| ≲ 10⁻¹⁵ but does not set it to zero, so this difference must be checked in position space (§2c).
- **Interpretation** [conj], for c13 = 0 only. This is consistent with BS: at linear order, matter then feeds the
  khronon only through gradient-type couplings, for which their instantaneous piece vanishes. §2c shows that this
  fails for c13 ≠ 0.

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

**Answer to the P3 question (scope corrected, §2c).** At linear order with c13 = 0, no: with ordinary matter coupled
only to the metric, the khronometric "instantaneous" term does not enter the receiver's tidal tensor, and the response
equals EA's. With c13 ≠ 0 it does enter (§2c). Non-linearly the
question is open; it is not shown to exist and not shown to be absent. Using the word "instantaneous" in the
literature is not evidence of a usable channel.

## 2c. P3 with c13 ≠ 0: the linear khronometric response has an instantaneous part (owner review 2)

**Scope correction.** GW170817 bounds |c13| ≲ 10⁻¹⁵ [lit: 1802.04303, via round 8] but does not set c13 = 0. §2b's
"the instantaneous term does not enter at linear order" is established for c13 = 0 only. This section checks
c13 ≠ 0 in position space. It stays at linear order: no second-order or non-linear step was started.

**3D kernel** [derived; `scripts/goal_khronometric.py` `c13_check`, `results/goal_round9/khronometric.json`]:
- **Ansatz.** E_ij = a k̂_ik̂_jρ + b k̂_ik̂_jS^T + c P_ijρ + d P_ijS^T + v(k̂_iV_j + k̂_jV_i) + t S^TT_ij, with
  P = 1 − k̂k̂, V = Pj, S^T = P_mnS_mn and S^TT = PSP − PS^T/2. It reproduces every helicity-frame component of EA
  and of khronometric exactly (residual 0).
- **Only v differs.** ΔE_ij = E_khronometric − E_EA = A ω (k_iP_jm + k_jP_im) j_m / (ω² − c_V²k²), with
  A = 4πG c13² / [c14(1 − c13)²]. This is the owner's form; the residual is 0 in every direction.
  - Khronometric alone: v = −8πGω / ((1 − c13)k), t = 8πGω² / (k² − (1 − c13)ω²). The scalar coefficients equal EA's.
- **Reduction check.** Khronometric depends on c1 only through c13 and c14 (∂/∂c1 = 0 at fixed c13 and c14), as a
  hypersurface-orthogonal æther requires [lit: Jacobson 1001.4823v3].

**Support test** [analytic split; the full-tensor test is a program check at exact rational parameter points, not a general proof]:
- **The source is local.** The input is the conserved T^{μν} itself, not Pj. Any conserved source is generated by
  its stresses: j_i = k_mS_im/ω, ρ = k_mk_nS_mn/ω².
- **Algebraic split.** With k·j = ωρ and ω²/[k²(ω² − v²k²)] = 1/k² + v²/(ω² − v²k²):
  - ΔE_ij = [Aω(k_ij_j + k_jj_i) − 2Ac_V²k_ik_jρ] / (ω² − c_V²k²) + I_ij,
  - I_ij = −2A k_ik_jρ / k² (residual 0).
  - This split alone proves nothing, so the full tensor was tested.
- **Full-tensor test.** A bare k² = k_x² + k_y² + k_z² in a denominator of the 3D kernel is a spatially non-local
  term with no retardation. Factors ω^n (time integrals of a local source) and ω² − c²k² (retarded propagators) are
  causal. Exact rational parameter points:

| Kernel | Bare k² |
|---|---|
| GR limit (c_i → 0); EA at c13 = 0, 3/20 and 1/1000; khronometric at c13 = 0 | none |
| khronometric at c13 = 3/20 and 1/1000 | **yes** |
| khronometric − I_ij at c13 = 3/20 and 1/1000 | none |

- **No further cancellation.** Nothing else in the full tidal tensor cancels I_ij, and I_ij is the whole
  instantaneous part.
- **Controls.** GR and EA show none; EA is causal, as round 8 assumed.

**Position space** [derived]. Retarded prescription, a quiet prior state and fields that decay at infinity, the
same for both theories:
- **Instantaneous term.** I_ij(t, x) = −2A ∂_i∂_j∇⁻²ρ(t, x) = −κ_inst ∂_i∂_jΦ_N[ρ(t)], with
  κ_inst = 2c13²(1 − c14/2) / [c14(1 − c13)²]. It is the Newtonian tidal tensor of the source's **current** mass
  distribution, taken on the æther slice t.
- **Validation source.** A smooth, local, conserved, non-spherical source with the stresses of its drive:
  - T^{ij} = q̈_ij f, T^{0i} = −∂_j(q̇_ij f), ρ = ∂_i∂_j(q_ij f);
  - f is a compact C³ radial bump; q is a generic symmetric matrix switched on from 0.
  - Outside the source, I_ij = κ′ G I_mn(t) ∂_i∂_j∂_m∂_n(1/r), with κ′ = c13² / [c14(1 − c13)²] and I_mn = ∫ρx_mx_n.
  - An FFT of the kernel agrees at three exterior points to ≤ 4.4 % (`position_space_check`; box periodicity and
    grid). This is a numerical cross-check of the formula, not an experimental verification.
- **Window.** Every other term is retarded and vanishes at the receiver until the fastest front arrives, at
  (r − R)/max(c_S, c_T). A source switched on at t = 0 therefore gives the receiver **exactly I_ij** for
  0 < t < (r − R)/max(c_S, c_T).
- **Transient, not a static change.** The static response equals EA's (ΔE(ω = 0) = 0; `pole_bookkeeping` checks
  two components at one k). The c_T and c_S fronts that follow restore the static field. How the compensation
  splits between them was not computed.
- **Multipoles.** Energy and momentum conservation of the whole sender, drive included, fix ∫ρ and make ∫ρx move
  uniformly. The controllable part therefore starts at the quadrupole and falls as r⁻⁵.

**Status.** This is a **linear response inside the model** (khronometric, linearised about flat space with the æther
at rest, matter coupled to the metric only). It is not a demonstrated channel. It is symmetric between sender and
receiver. It is simultaneous on the æther foliation, so a single global foliation gives no closed loops [derived].

**Dependences** [derived]:
- **Amplitude:** κ′ G ΔI × (angular factor) / r⁵. The angular factor is 36 on axis for ΔI in I_xx − I_yy.
- **Arrival:** no retardation in the æther frame, so it does not depend on r.
- **Duration:** at most (r − R)/max(c_S, c_T).
- **Coupling:** κ′ ∝ c13²/c14, which grows as c14 → 0 at fixed c13. The c14 denominator is kept.

**Size at a fixed receive time τ** (`link_budget_c13`; round 8's generous receiver, ASD 10⁻²⁰ s⁻²/√Hz and SNR 5;
sender: 10¹² kg moved 10 m at 1 km non-spherically, ΔI = 2 × 10¹⁶ kg m², ASSUMPTION):
- SNR = |I_xx|√τ/ASD, valid while τ fits inside the window, i.e. max(c_S, c_T) ≤ r/(cτ).
- Required values:

| r | τ | κ′ needed | c14 needed at \|c13\| = 10⁻¹⁵ | window needs max(c_S, c_T) ≤ |
|---|---|---|---|---|
| 20 ly | 1 s / 1 d / 1 yr | 2.5 × 10⁵⁹ / 8.6 × 10⁵⁶ / 4.5 × 10⁵⁵ | 4 × 10⁻⁹⁰ / 1 × 10⁻⁸⁷ / 2 × 10⁻⁸⁶ | 6.3 × 10⁸ / 7,305 / 20 |
| 100 ly | 1 s / 1 d / 1 yr | 7.9 × 10⁶² / 2.7 × 10⁶⁰ / 1.4 × 10⁵⁹ | 1 × 10⁻⁹³ / 4 × 10⁻⁹¹ / 7 × 10⁻⁹⁰ | 3.2 × 10⁹ / 36,525 / 100 |

**What decides usability** (first pass; §2d supersedes the window-versus-coupling bullet with khronometric PPN bounds):
- **Window versus coupling.** c_S² = (c13 + c2)(2 − c14) / [c14(1 − c13)(2 + c13 + 3c2)] grows as c14 → 0.
  - With c2 ≥ 0 and c13 > 0, c_S² ≳ c13/c14, so κ′ ≲ 1.2 c13 c_S² ≲ 5 × 10⁻¹³ for the 20-ly, 1-yr window
    (c_S ≤ 20). That is ~68 orders below need [derived].
  - Only a tuned c13 + c2 ≈ 0 (c13 < 0, c2 ≈ |c13|, tuned to ~c14 (r/cτ)², i.e. ~10⁻⁸⁴ in absolute terms) escapes
    this.
  - Whether such a point is allowed in **khronometric** theory was not checked. The region (i)/(ii) bounds quoted
    in round 8 are EA fits (1802.04303), and khronometric has no vector mode.
- **EFT range.** c14 → 0 is the limit where the æther equation loses its time derivatives (round 8). No
  strong-coupling scale was read. Perturbative validity of the linear khronon response at κ′ ≫ 1 was not checked.
- **Existing data.** For any accelerated mass, a κ′ ≫ 1 transient would exceed the Newtonian tidal field of the
  same mass change by that factor. Whether existing data already exclude this was not checked.
- **Decision-changing result.** An allowed khronometric parameter point with κ′ ≳ 10⁵⁵ and
  max(c_S, c_T) ≤ r/(cτ), inside a controlled EFT. It is neither shown nor excluded here.

**Answer.** With c13 ≠ 0, the linear khronometric tidal response to ordinary, local, conserved matter has an
instantaneous, spatially non-local part, I_ij = −2A k_ik_jρ/k². Nothing else in the full tidal tensor cancels it.
It is a model-internal linear response, not evidence of a usable channel. Its size and window are tied to c13²/c14
and c_S as above. For the fastest-front scheme with the designated devices, §2d closes it within the controlled
domain.

## 2d. PPN bounds × receive window: the fastest-front scheme with the designated devices (owner review 3)

**Scope.** This section covers only one scheme: receiving **before the fastest front arrives**, with round 8's
designated sender and receiver (§2c table). It does not exclude other schemes, such as reading after the fronts
arrive or repeated modulation. Those need the full time response (every pole); the I_ij amplitude cannot be
integrated past the window. No second-order, non-linear or large parameter scan was run.

**Notation.** Khronometric only; EA parameter regions are not used. α = c14, β = c13, λ = c2, and
κ′ = β²/[α(1 − β)²].

**Identity** [analytic: derived by hand, then checked symbolically with general symbols (`ppn_identity`, residual 0)].
- **Inputs** from Yagi, Blas, Barausse & Yunes, arXiv:1311.7144v4 [sha256 445976e4…]:
  - Eq. 48: α1 = 4(α − 2β)/(β − 1).
  - Eq. 49: α2 = (α − 2β)[−β² + β(α − 3) + α + λ(−1 − 3β + 2α)] / [(β − 1)(λ + β)(α − 2)].
  - Eq. 118: c_s² = (α − 2)(β + λ) / [α(β − 1)(2 + β + 3λ)]. This equals the c_S² of §2c.
- **Conventions.** The paper uses signature (+,−,−,−) and defines α, β, λ in the khronometric action as used here.
  The formulas are ratios of couplings, so the signature does not change them.
- **Steps:**
  - α2 − α1/8 = (α − 2β)²(2 + β + 3λ) / [2(α − 2)(β − 1)(β + λ)]. The bracket 2N − (λ + β)(α − 2) factorises to
    (α − 2β)(2 + β + 3λ).
  - κ′ − (β − α/4)/(1 − β)² = (α − 2β)² / [4α(1 − β)²].
  - Therefore **κ′ = (c_s²/2)(α2 − α1/8) + (β − α/4)/(1 − β)²** for every finite point, either sign of β and any λ.
    This is the owner's form.
- **Singular and degenerate boundaries:**
  - α = 0: κ′ and c_s² diverge.
  - λ + β = 0: c_s² = 0, and α2 diverges unless α = 2β.
  - α = 2: c_s² = 0, and α2 diverges.
  - β = 1.
  - 2 + β + 3λ = 0: c_s² diverges.
  - c_s² < 0: gradient instability; no bound is derived there.
  - α = 2β: α1 = α2 = 0 and κ′ = β/[2(1 − β)²] ≲ 10⁻¹⁵.

**Bound** [analytic]:
- **Assumptions:** α > 0 and c_s² ≥ 0; |β| ≲ 10⁻¹⁵ of either sign; |α1| ≲ 10⁻⁴ and |α2| ≲ 10⁻⁷ (as quoted in
  1311.7144v4 §II).
- **Result:** κ′ ≤ (c_s²/2)(|α2| + |α1|/8) + |β|/(1 − β)² = **6.3 × 10⁻⁶ c_s² + 1.0 × 10⁻¹⁵**.
  - The −α/4 term only lowers κ′.
  - Tuning λ ≈ −β acts only through c_s², and the window bounds c_s.
- **Window:** a receive time τ inside the window requires max(c_S, c_T) ≤ X = D/(cτ). So c_s ≤ X, and
  κ′ ≤ 6.3 × 10⁻⁶ X² + 10⁻¹⁵.

**Numbers** [arithmetic on the formulas above; `ppn_window_check`]:

| D | τ | X = D/(cτ) | κ′ needed | κ′ allowed (PPN) | shortfall (orders) | α2 implied by the need | α2 v²Φ at 1 AU (v = 30 km/s) |
|---|---|---|---|---|---|---|---|
| 20 ly | 1 s | 6.3 × 10⁸ | 2.5 × 10⁵⁹ | 2.5 × 10¹² | 47 | 1.3 × 10⁴² | 1.2 × 10²⁶ |
| 20 ly | 1 d | 7,305 | 8.6 × 10⁵⁶ | 3.4 × 10² | 54 | 3.2 × 10⁴⁹ | 3.1 × 10³³ |
| 20 ly | 1 yr | 20 | 4.5 × 10⁵⁵ | 2.5 × 10⁻³ | 58 | 2.2 × 10⁵³ | 2.2 × 10³⁷ |
| 100 ly | 1 s | 3.2 × 10⁹ | 7.9 × 10⁶² | 6.3 × 10¹³ | 49 | 1.6 × 10⁴⁴ | 1.5 × 10²⁸ |
| 100 ly | 1 d | 36,525 | 2.7 × 10⁶⁰ | 8.4 × 10³ | 57 | 4.0 × 10⁵¹ | 3.9 × 10³⁵ |
| 100 ly | 1 yr | 100 | 1.4 × 10⁵⁹ | 6.3 × 10⁻² | 60 | 2.8 × 10⁵⁵ | 2.7 × 10³⁹ |

- **Formal edge.** The allowed κ′ would reach the need only for τ ≲ 5 × 10⁻³² s (20 ly) or 2 × 10⁻³³ s (100 ly).
  At such τ the flat-ASD receiver model has no meaning.

**Q1. Do the PPN derivation and bounds apply at tiny α and tiny λ + β?** [derived]
- **What the formulas assume.** They are exact in the couplings, with no small-coupling expansion. They come from
  a weak-field, slow-motion (1PN) expansion with the æther/khronon perturbed linearly about the solar-system
  background. They therefore apply at small α and small λ + β only if all of these hold:
  - (a) no singular boundary above is reached;
  - (b) c_s ≫ v ≈ 10⁻⁴–10⁻³, so the khronon is quasi-static on solar-system scales (PN ordering). The |α2| bound
    from the Sun's spin axis is derived in this regime.
  - (c) the perturbations stay small:
    - the preferred-frame metric term, α2 v²Φ ≪ 1;
    - the æther's O(v) coefficient, ≈ (α − 2β)/(λ + β) ≈ 2α2/(α − 2β) (from Yagi Eq. 53), times vΦ ≪ 1;
  - (d) the EFT strong-coupling scale is not exceeded (not computed).
- **At the points the device needs, (c) fails.** The identity forces α2 ≥ 2(κ′_needed − 10⁻¹⁵)/X² − |α1|/8 ≥ 10⁴²
  in every row. Then α2 v²Φ at 1 AU is ≥ 10²⁶ with v at Earth's orbital speed (a floor over a year), and ≥ 10²⁸
  at the CMB dipole speed (ASSUMPTION for the æther frame).
- **Consequence.** There the solar-system PPN expansion loses control. That alone does not show that every
  flat-background linear solution fails. It does mean that, at the extreme parameters required, no controlled
  prediction has been established for the actual sender and receiver environments, which are themselves inside
  gravitating, moving systems. These points are labelled **"no controlled prediction yet"** and cannot be treated
  as a feasible exception.
- **Answer to Q1.** The bounds apply wherever the expansion is under control. At the points the device would
  require they do not apply, and no controlled prediction for the actual sender and receiver environments has
  been established there either.

**Q2. Do the existing inequalities close the tuned scheme?** [derived]
- **Inside the controlled domain: yes.**
  - Every point that satisfies the PPN bounds has κ′ ≤ 6.3 × 10⁻⁶ X² + 10⁻¹⁵. That is 47–60 orders of magnitude
    below the need.
  - This holds for either sign of β, for the λ ≈ −β tuning (which acts only through c_s ≤ X) and for c_s < c_T,
    where the window is set by c_T and c_s ≤ X still holds.
  - Perturbativity alone gives the same answer without the observed α2 value: α2 v²Φ < 1 at 1 AU caps κ′ at
    X²/(2v²Φ) + 10⁻¹⁵ ≈ 2 × 10¹⁸–5 × 10³⁴, still 26–39 orders below need.
- **Outside it:** the points with α2 v²Φ ≫ 1 are "no controlled prediction yet" (as above). Settling them would
  need a controlled treatment of the real environments, possibly non-linear or strong-coupling. It was not
  started.
- **Stop.** The bound is enough for this scheme, so no strong-coupling scale was computed.

**Status (owner decision, 2026-10-11):**
- **Paused:** the designated-device, fastest-front reception scheme. This round's closing conclusion is accepted.
- **Undetermined, not started:**
  - reception after the fronts arrive;
  - repeated modulation;
  - the non-linear and strong-coupling questions (including the §2b elliptic part).
  None is started automatically.
- **Kept, with its standing:** the model-internal instantaneous linear response at c13 ≠ 0 (§2c). It is a derived
  property of the linearised model. It is **not** an experimental finding and **not** a usable communication
  channel.

**Kinds of evidence in §§2c–2d:**
- **Analytic:** the helicity ansatz and ΔE form (symbolic in all c_i), the split, the PPN identity and the bounds.
- **Program checks at exact rational parameter points** (not general proofs): the bare-k² tests (§2c table) and
  the pole bookkeeping (one k, two components).
- **Numerical cross-check:** the FFT agreement (≤ 4.4 %) checks the position-space formula. It is not an
  experimental verification.

## 3. Comparison and choice

| | Different prediction | Data or model available | Light test | Changes a communication decision |
|---|---|---|---|---|
| P1 v-causal | yes (S4 below the realistic QM prediction; B ∼ C signalling) | **no** (no four-party spacelike data) | — | yes, if positive |
| P2 exotic lenses | yes (W3: flux to zero) | yes, but already scheduled in another program | the existing pipeline | a survivor is a candidate only; nulls need an abundance model |
| P3 khronometric | c13 = 0: not at linear order (§2b); c13 ≠ 0: yes, a linear instantaneous transient ∝ c13²/c14 (§2c); non-linear part open | linear: done; usability needs the allowed khronometric (c13, c14, c2) region and an EFT range; non-linear: needs a second-order calculation | linear done; non-linear not light | the fastest-front scheme with the designated devices is closed within the controlled domain (§2d: 47–60 orders short); other schemes need the full time response; or the non-linear part gives a matter-driven response |

**Result.** Not enough basis to select a further executable unit beyond the P3 linear checks done here (§2b, §2c).
- The c13 ≠ 0 linear term exists in the model (§2c). For the fastest-front scheme with the designated devices, the
  khronometric PPN bounds close it inside the controlled domain; the points it would need are outside control (§2d).
  Schemes after the fronts arrive are not covered.
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
   result with a matter-only source. For the linear c13 ≠ 0 term: a scheme other than fastest-front reception, analysed
   with the full time response (§2d scope).

## 4. Plain-language summary

- **New evidence?** No new observational evidence. Two derived results:
  - At linear order with c13 = 0, the Hořava/khronometric "instantaneous" term does not reach a receiver's tidal
    measurement when the sender only moves ordinary matter. The theory then behaves exactly like round 8's
    Einstein–æther case.
  - With c13 ≠ 0 (allowed: GW170817 gives |c13| ≲ 10⁻¹⁵, not 0), it does. Moving ordinary matter gives a tidal
    transient that arrives with no delay in the æther frame. Its size is c13²/c14 times the Newtonian tidal
    field of the quadrupole change, and it lasts until the khronon and tensor fronts arrive. This is a
    model-internal linear response. For a lab-scale sender it needs κ′ ≳ 10⁵⁵ (c14 ≲ 10⁻⁸⁶ at |c13| = 10⁻¹⁵).
  - The khronometric PPN identity κ′ = (c_s²/2)(α2 − α1/8) + (β − α/4)/(1 − β)² with the PPN bounds and the
    receive window caps κ′ at 47–60 orders below that, including the λ ≈ −β tuning.
  - At the points that would reach it, the solar-system PPN expansion fails. No controlled prediction for the
    actual sender and receiver environments has been established there, so they cannot count as a feasible
    exception (§2d).
  - This covers only the fastest-front scheme with the designated devices, which is now paused. Reception after
    the fronts, repeated modulation and the non-linear questions are undetermined and not started.
  - None of this is an experimental finding or a usable channel.
  - [lit]: the controls (Blas–Sibiryakov instantaneous term and its EA cancellation; Bhattacharyya et al.'s
    elliptic mode).
  - [derived]: the khronometric = EA linear response (c13 = 0) and the c13 ≠ 0 instantaneous part.
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
- Computation: `scripts/goal_khronometric.py` (~80 s), results `results/goal_round9/khronometric.json`, tests
  `tests/test_goal_khronometric.py`.
- Corrections in this revision (owner review): the Jacobson citation (1003.4987 → 1001.4823); W3 survivors are
  candidates only; the Ellis and negative-mass templates are separated; the ~10³⁸ M☉ figure is scoped to round 6's
  model; a low S4 needs a realistic QM prediction; the "only route" remark is withdrawn. Owner review 2: the
  linear conclusion is scoped to c13 = 0, and c13 ≠ 0 is checked in §2c (3D kernel, full-tensor support test,
  position space, signal size). Owner review 3: the khronometric PPN identity and bounds against the receive
  window (§2d), read from 1311.7144v4 Eqs. 48, 49 and 118.
- State files: CHANGELOG entry and a proposed DECISIONS record (D-TBD), with no change to TASKS priorities.
