# B on R2-A: destruction report for "Co-Witness Consistency" (CWC)

System B ("Physics Destruction"), round 2. Target: `R2A-consistency.md`. Labels: **[lit]** = from a
reference I opened (arXiv IDs at the end); **[B-est]** = my own estimate (scratch script, not committed); **[std]** = textbook
result.

**Verdict.** P1 is the renamed |R|² = 1 Planck-surface exotic compact object (ECO); its closure amplitude is a
published theorem, its premise dies to the ergoregion instability, and its own kill criterion is already met. P2 and
P3 restate GR.

## 1. Internal consistency

1. **A4 does not forbid black-hole event horizons, so "all carriers must come back out" is a non sequitur.** An event
   horizon is the boundary of J⁻(𝓘⁺), not a pair with disjoint futures. A4 asks only for *some* common witness, and an
   outside observer can supply one by falling in. In Kruskal coordinates J⁺(r) ∩ J⁺(s) is empty only when
   max(U)·max(V) ≥ 1, i.e. beyond the singularity [std]. A4 forbids **singular endings and cosmological event
   horizons**. A future-complete trapped region satisfies A4 and still swallows everything. §4.3 needs an unstated
   [HYP]: the witness must be reachable from 𝓘⁺.
2. **The closure rule picks one target arbitrarily.** If every carrier that crosses the would-be horizon must return,
   that includes all the collapsed stellar matter and the two progenitor holes. The remnant would then re-emit ~Mc²,
   not E_hor. CWC limits the rule to the perturbative ringdown flux without saying why. That restriction is [HYP], not
   [DER].
3. **Dense records versus counting.** A1 makes every interval infinite, so chains have no "tick count", and §2.2's
   radar coordinates and A5's counts are undefined. Order fixes geometry only up to a conformal factor (Malament 1977;
   Hawking–King–McCarthy 1976); causal sets recover volume by counting, which A1 forbids. CWC has no proper time, mass
   scale or k: its "derived" kinematics is conformal.
4. **The exclusion of κ = 0 is wrong.** The Galilean order (s ≻ r iff t_s > t_r) is acyclic: simultaneous events are
   unrelated, so no cycle r ≺ s ≺ r arises. Only the k-calculus premise of a finite cone excludes it. The Lorentz
   result is Zeeman (1964) plus radar with the conclusion in the inputs: circular, though harmless.
5. **The rational-points model is a model, not a derivation.** Q⁴ with Minkowski causal order does satisfy A0, A1
   and A4. But it is built from real Minkowski space, so the metric is an input. Its automorphisms (rational Poincaré
   plus dilations) form a countable group, A3 transitivity holds only for rational velocities, and "unaccelerated"
   needs a metric.
6. **Noether and Lovelock need structure CWC lacks.** Noether needs an action and a Lie symmetry (SO⁺(3,1;ℚ) is
   not one). Lovelock needs a metric. §3.3–3.4 import GR wholesale.
7. **§3.4 contradicts §4.2.** In unimodular gravity Λ is an integration constant. For a given solution it is constant
   forever, so any Λ > 0 gives a future de Sitter event horizon and violates A4. CWC therefore needs Λ_int = 0. With
   A6 forbidding every dimensionful dark-energy parameter, CWC has **no source for the observed acceleration**. The
   claim that "the acceleration must end" assumes an acceleration that CWC cannot produce.
8. §4.2 numbers check out [B-est]: χ_eh = 5.08 Gpc, χ_p = 14.17 Gpc, 2χ_eh at z = 12.8.

## 2. Conservation, stability and causality

**EHT shadows: no kill.** The shadow is set by the light ring, which a Planck-surface ECO has. For zero absorption
the Sgr A* luminosity bound only needs μ ≡ 1 − 2M/R* ≲ 10⁻¹⁶ [lit 2205.13555]. A proper distance ℓ_P gives
μ ≈ ℓ_P²/(16M²) ≈ 5×10⁻⁹¹ [B-est], since trapped radiation needs ~1/μ bounces to escape. CWC passes; so does M87
[1906.11238, 1906.11243].

**Kerr ringdown: no kill.** The prompt ringdown probes only the light ring [lit 1602.07309]; GW250114's (2,2,0)/(2,2,1)
modes are within ±30 % of Kerr [lit 2509.08054]. CWC's delay for it (68.4 M☉ detector frame, χ_f = 0.68) is 0.292 s
[B-est], exactly the Planck-surface delay Prasad et al. probe [lit 2610.12429].

**Ergoregion instability: fatal.** A perfectly reflecting Kerr-like object is unstable on dynamical timescales
[lit 1703.03696]. Gravitational perturbations need **≥ 0.3 % absorption at χ ≲ 0.7, ≥ 6 % at χ ≲ 0.9, ≈ 60 % for any
spin** [lit 1807.08840]. CWC's closure fixes absorption at exactly 0. The "phase-scrambled" escape fails: ω_I ≈ |𝒜|²/t₀
is a bounce-and-amplify argument on *energy* [lit 1805.08229]; each round trip the band ω < mΩ_H gains (1+Z) and A5
returns all of it.

Numbers [B-est], from eq. (2) of 1805.08229 with t₀ = Δt_CWC/2 ≈ 430–600 t_M:

| Remnant | Unstable-mode frequency f_R | e-folding τ | Rotational energy (fraction of M) |
|---|---|---|---|
| GW250114 | 184 Hz | ≈ 4×10³ s | 6.9 % (≈ 4.3 M☉c²) |
| 20 M☉, χ = 0.7 | 650 Hz | ≈ 1×10³ s | 7.4 % |
| 150 M☉, χ = 0.7 | 87 Hz | ≈ 9×10³ s | 7.4 % |

The seed is an O(1) merger perturbation, so each remnant would shed more than its merger radiated (≈ 3.1 M☉ for
GW250114, from m₁ + m₂ − M_f [2509.08054]) as an ~hour-long quasi-monochromatic signal near 180 Hz [B-est].
At population level the missing background already gave **X < 50 % of all compact objects in O1, even at low natal
spin**, and X < 1 % at design [lit 1805.08229]. O3 gives Ω_GW ≤ 5.8×10⁻⁹ [lit 2101.12130] (I did not re-derive X:
different band). CWC needs X = 1. **P1's premise is excluded before any echo analysis.**

X-ray-binary spins are a weaker test [lit 1805.08229], and GW remnant spins (measured within ms) do not conflict.

## 3. Already excluded? (P1)

**Claims.** Abedi et al., 2.5σ in O1 [1612.00266]; GW170817 [1803.10454]; GW190521 [2201.00047]. Abedi 2023 finds
nothing significant in 47 events and A < 0.4 (90 %) against a canonical A ~ 1 [2301.00025].

**Null results.** Ashton et al. [1612.05625]; Westerweck et al., "entirely consistent with noise" [1712.09966];
Nielsen et al. [1811.04904]; Lo et al. [1811.07431]; Uchikata et al. [1906.00838, 2309.01894]; Tsang et al.,
smallest p = 0.03 [1906.11168]; LVK testing-GR papers [1903.04467, 2010.14529, 2112.06861]. Also:
- GWTC-4.0 remnants paper: four echo tests and 42 O4a events, "no evidence for post-merger echoes" [2603.19021]. Its
  cWB test sums coherent energy in **[t_echo − 0.05 s, 4t_echo + 0.05 s], 16–1024 Hz**, which is P1's window
  "[Δt, a few Δt]", and it is morphology-independent, so phase is irrelevant.
- Wu et al., phase-marginalised long-lived QNMs in GW150914, GW231226 and GW250114 [2512.24730].

**Did anyone constrain energy relative to the horizon-absorbed energy?**
- LVK did not. Its amplitude A is relative to the truncated IMR waveform.
- Miani et al. (cWB, < 1 s after merger) give h_rss(echo)/h_rss(merger–ringdown) ≤ 0.11–0.34 for 15 GWTC-3 events.
  That is 0.21 for GW150914, so **E_late ≲ 1–12 % of E_MR** [lit 2302.12158].
- **Prasad, Sarkar & Sarkar (2026-10-08, preprint)** do exactly what CWC calls "never measured". They calibrate the
  inward plunge radiation of GW250114 with its Direct Wave and bound the reflected energy fraction to
  **|R_s|² < 3.5×10⁻³ (90 %)** for surfaces from 0.1 M_f down to ℓ_P. |R_s| = 0.2 would already have given network
  SNR ≈ 11 [lit 2610.12429]. Caveats: the template is coherent with one global phase, it covers only the first return,
  and it relies on a point-particle plunge model.

CWC's "quantity 1" is therefore not new. Mark et al. had already proved that for |R| = 1 the echo energy **equals**
the energy that would have crossed the BH horizon, with "> 97 %" radiated for x₀ > 20M [lit 1706.06155]. That is
CWC's closure rule, published in 2017.

**How much energy does GR send into the horizon?** In nonlinear GR this is not gauge-invariant. The event horizon is
teleological and eventually swallows the whole binary, so CWC's claim that it is "deterministic, not tunable" holds
only after a choice of perturbative start time. The estimates I could verify:

- Inspiral tidal absorption is tiny: ≤ 0.1 GW cycle to b = 2M [lit gr-qc/0107080].
- Around merger, ingoing to outgoing flux is "around 1:1" (Wang & Afshordi [1803.02845], reading the NR of
  Gupta et al. [1801.07048]; Gupta et al. themselves show normalised flux shapes).
- For R = 1, total ringdown-plus-echo energy is up to ≈ 38× the ringdown, and total SNR ≈ 18 ρ_ringdown
  [lit 1806.04253].

Taking f ≡ E_hor/E_post-peak ≈ 0.3–1 is therefore conservative. Nothing I found supports CWC's "percent-level" worry.

**Predicted SNR [B-est].** GW250114 has post-peak SNR 26 at t ≥ 6 t_M [lit 2509.08054]. An incoherent excess-power
statistic over T ≈ 0.88 s and B ≈ 100 Hz (N ≈ 88 pixels) gives S ≈ fρ²/√N:

| f | GW250114 | GW150914 (ρ_RD ≈ 8) |
|---|---|---|
| 0.05 | 3.6σ | — |
| 0.3 | 22σ | 2σ |
| 1 | 72σ | 7σ |

If phase survives, ρ ≈ 26√f; scaling Prasad's ρ ≈ 11 at |R_s| = 0.2 linearly gives ≈ 55 for the first return alone. Stacking the GWTC-3/4 events
with ringdown SNR ≳ 8 adds roughly another 1–2× in ρ² [B-est]. **The predicted SNR is ≥ 5 from GW250114 alone for any
f ≳ 0.07.** The observed limits are E_late/E_MR ≲ 0.04 (Miani, GW150914), the null cWB windows, and |R|² < 0.0035
(Prasad). By CWC's own kill criterion (observed < ⅓ of predicted at 95 %, predicted SNR ≥ 5), **P1 is killed, not
untested.**

The only escape, f ≲ 0.05, needs a GR horizon flux 10× below every estimate above. That nobody has run CWC's exact
stacked statistic is a bookkeeping gap, not a reprieve.

## 4. Renaming check

| CWC element | Existing name | Verdict |
|---|---|---|
| A4 "no horizons" + Planck turning surface | Horizonless ECO, "ClePhO" / canonical perfectly reflecting model [1805.08229, 1904.05363] | Renamed |
| Δt formula | Cardoso et al. [1602.07309, 1608.08637]; Abedi et al. [1612.00266] | Admitted |
| Energy closure Σ E_n = E_hor | \|R\| = 1 echo-energy theorem [1706.06155]; measured as \|R_s\|² [2610.12429] | Renamed |
| "Phase not fixed" | Complex or frequency-dependent reflectivity [1907.03091]; free φ in BHP template [2603.19021]; morphology-independent searches [1906.11168, 2302.12158] | Renamed |
| Λ as integration constant | Unimodular / trace-free Einstein [1008.1196] | Renamed, and contradicts A4 (§1.7) |
| "No eternal de Sitter" | Swampland dS conjecture [1806.08362] | Resemblance only |
| Dense causal order | Zeeman (1964); Malament (1977); HKM (1976); escapes Bombelli–Henson–Sorkin [gr-qc/0605006] only by giving up volume | Known |

## 5. Novelty verdict

- **P1:** not novel. The closure rule is in Mark et al. (2017), and the energy-normalised test was posted on 2026-10-08
  (2610.12429). The incoherent window test is LVK's cWB echo test.
- **P2:** not novel; the SR/GR null. c_g − c ∈ [−3×10⁻¹⁵, +7×10⁻¹⁶] c [1710.05834]; E_QG,1 > 10 E_Pl
  [2402.06009]; GRB 090510 [0908.1832]; lab isotropy 10⁻¹⁸ [1412.6954]. A pass cannot separate CWC from GR.
- **P3:** not novel; the standard distance-duality / opacity test [1004.2053]. Keil et al. find a parametrised
  **6σ "violation"** with SH0ES + BBN and ≤ 1σ model-independently [2504.01750]: the H₀ calibration already mimics a
  kill.
- **Δη_∞:** an extrapolation, as CWC itself admits, and inconsistent with its own Λ treatment (§1.7).

## 6. Triage

| Prediction | Verdict | Reason |
|---|---|---|
| **P1** energy-closure afterglow | **FAIL (excluded)** | Premise excluded by the ergoregion instability for \|R\|² = 1 (needs ≥ 0.3–6 % absorption). O1 stochastic background caps X at 50 % while CWC needs X = 1. Predicted SNR ≫ 5 (GW250114 alone, f ≳ 0.07). Observed: \|R\|² < 0.0035 (GW250114), E_late/E_MR ≲ 0.01–0.12 (GWTC-3 cWB), LVK GWTC-4 cWB null in P1's window. Renamed R = 1 ECO. |
| **P2** zero Lorentz / c_g violation | **FAIL (non-discriminating)** | Same prediction as GR. Already measured. Survival means nothing. |
| **P3** η = 1 exactly | **FAIL (non-discriminating, mimicked)** | Same as GR/ΛCDM. H₀-calibration systematics already produce published "violations". |
| Δη_∞ (cosmology) | **FAIL (not a prediction)** | No present-day amplitude. Contradicts §3.4. |

**No prediction earns PASS or conditional PASS.** For the record, P1 data: GWOSC strain plus public M_f, χ_f
posteriors. Mimics: glitches (e.g. the ~30–40 Hz post-merger excesses in GW190701 and GW200224 [2302.12158]),
non-Gaussian noise, imperfect primary subtraction. Sensitivity already published (E_late/E_MR ~ 1–10 % per loud event,
|R|² ~ 10⁻³ for GW250114) excludes CWC.

**Rescue** needs ≥ 0.3–6 % surface absorption, which breaks A5 closure and turns CWC into the known partially
reflecting ECO family. CWC also still needs a dark-energy mechanism compatible with A4 and A6.

## References (all opened on 2026-10-09; arXiv IDs as given; titles checked on the abstract page)

- 1612.00266 Abedi, Dykaar, Afshordi, "Echoes from the Abyss: Tentative evidence for Planck-scale structure at black hole horizons" (PRD 96, 082004)
- 1612.05625 Ashton et al., "Comments on: 'Echoes from the abyss…'"
- 1712.09966 Westerweck et al., "Low significance of evidence for black hole echoes in gravitational wave data" (PRD 97, 124037)
- 1811.04904 Nielsen et al., "Parameter estimation for black hole echo signals and their statistical significance" (PRD 99, 104012)
- 1811.07431 Lo, Li, Weinstein, "Template-based Gravitational-Wave Echoes Search Using Bayesian Model Selection" (PRD 99, 084052)
- 1906.00838 Uchikata et al., "Searching for black hole echoes from the LIGO-Virgo Catalog GWTC-1" (PRD 100, 062006)
- 2309.01894 Uchikata et al., "Searching for gravitational wave echoes from black hole binary events in the third observing run…"
- 1906.11168 Tsang et al., "A morphology-independent search for gravitational wave echoes in data from the first and second observing runs…" (PRD 101, 064012)
- 1803.10454 Abedi & Afshordi, "Echoes from the Abyss: A highly spinning black hole remnant for the binary neutron star merger GW170817" (JCAP 1911, 010)
- 2201.00047 Abedi, Micchi, Afshordi, "GW190521: Search for Echoes due to Stimulated Hawking Radiation from Black Holes"
- 2301.00025 Abedi, "Search for echoes on the edge of quantum black holes"
- 2302.12158 Miani et al., "Constraints on the amplitude of gravitational wave echoes from black hole ring-down using minimal assumptions"
- 2512.24730 Wu et al., "Model-agnostic search of gravitational wave echoes in LVK data" (PRD 113, 124023)
- 2610.12429 Prasad, Sarkar, Sarkar, "Listening to the Horizon: Probing Near-Horizon Reflectivity with GW250114" (preprint, 2026-10-08)
- 1903.04467 LVC, "Tests of General Relativity with the Binary Black Hole Signals from the LIGO-Virgo Catalog GWTC-1"
- 2010.14529 LVC, "Tests of General Relativity with Binary Black Holes from the second LIGO-Virgo Gravitational-Wave Transient Catalog"
- 2112.06861 LVK, "Tests of General Relativity with GWTC-3"
- 2603.19021 LVK, "GWTC-4.0: Tests of General Relativity. III. Tests of the Remnants"
- 2509.08054 LVK, "GW250114: testing Hawking's area law and the Kerr nature of black holes"
- 1602.03841 LVC, "Tests of general relativity with GW150914"
- 1602.07309 Cardoso, Franzin, Pani, "Is the gravitational-wave ringdown a probe of the event horizon?" (PRL 116, 171101)
- 1608.08637 Cardoso et al., "Echoes of ECOs: gravitational-wave signatures of exotic compact objects and of quantum corrections at the horizon scale" (PRD 94, 084031)
- 1904.05363 Cardoso & Pani, "Testing the nature of dark compact objects: a status report" (Living Rev. Rel. 22, 4)
- 1706.06155 Mark, Zimmerman, Du, Chen, "A recipe for echoes from exotic compact objects" (§V energy theorem)
- 1803.02845 Wang & Afshordi, "Black Hole Echology: The Observer's Manual" (PRD 97, 124044)
- 1801.07048 Gupta, Krishnan, Nielsen, Schnetter, "Dynamics of marginally trapped surfaces in a binary black hole merger…" (PRD 97, 084028)
- 1806.04253 Testa & Pani, "Analytical template for gravitational-wave echoes…" (PRD 98, 044018)
- 1907.03091 Maggio et al., "Analytical model for gravitational-wave echoes from spinning remnants" (PRD 100, 064056)
- gr-qc/0107080 Alvi, "Energy and angular momentum flow into a black hole in a binary" (PRD 64, 104020)
- 1703.03696 Maggio, Pani, Ferrari, "Exotic Compact Objects and How to Quench their Ergoregion Instability" (PRD 96, 104047)
- 1807.08840 Maggio et al., "Ergoregion instability of exotic compact objects: electromagnetic and gravitational perturbations and the role of absorption" (PRD 99, 064007)
- 1805.08229 Barausse et al., "The stochastic gravitational-wave background in the absence of horizons" (CQG 35, 20)
- 2101.12130 LVK, "Upper Limits on the Isotropic Gravitational-Wave Background from Advanced LIGO's and Advanced Virgo's Third Observing Run"
- 2205.13555 Carballo-Rubio et al., "Constraints on horizonless objects after the EHT observation of Sagittarius A*"
- 1906.11238 and 1906.11243 EHT Collaboration, "First M87 EHT Results. I. The Shadow…" and "VI. The Shadow and Mass of the Central Black Hole"
- 1710.05834 LVC + Fermi-GBM + INTEGRAL, "Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A"
- 2402.06009 LHAASO, "Stringent Tests of Lorentz Invariance Violation from LHAASO Observations of GRB 221009A"
- 0908.1832 Fermi LAT/GBM, "Testing Einstein's special relativity with Fermi's short hard gamma-ray burst GRB090510" (Nature 462, 331)
- 1412.6954 Nagel et al., "Direct Terrestrial Test of Lorentz Symmetry in Electrodynamics to 10⁻¹⁸" (Nat. Commun. 6, 8174)
- 1004.2053 Avgoustidis et al., "Constraints on cosmic opacity and beyond the standard model physics from cosmological distance measurements"
- 2504.01750 Keil et al., "Probing the Distance Duality Relation with Machine Learning and Recent Data"
- 1008.1196 Ellis et al., "On the Trace-Free Einstein Equations as a Viable Alternative to General Relativity" (CQG 28, 225007)
- 1806.08362 Obied, Ooguri, Spodyneiko, Vafa, "De Sitter Space and the Swampland"
- gr-qc/0605006 Bombelli, Henson, Sorkin, "Discreteness without symmetry breaking: a theorem"
- Not on arXiv (verified via Crossref DOI): Zeeman, "Causality Implies the Lorentz Group", J. Math. Phys. 5, 490 (1964), doi:10.1063/1.1704140; Malament, "The class of continuous timelike curves determines the topology of spacetime", J. Math. Phys. 18, 1399 (1977), doi:10.1063/1.523436; Hawking, King & McCarthy, "A new topology for curved space–time which incorporates the causal, differential, and conformal structures", J. Math. Phys. 17, 174 (1976), doi:10.1063/1.522874.
