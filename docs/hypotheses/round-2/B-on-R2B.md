# System B on R2-B ("Rate-Cone World"): round 2

Target: `R2B-spectrum.md`. Labels: **[lit]** = verified in a paper listed below;
**[calc]** = my arithmetic on published numbers; **[arg]** = my theoretical argument.

## Verdict in one paragraph

The scalar sector is GR renamed. The one new claim (saddle images mirror χ and flip V) does not follow
from the axioms: A6 plus superposition **derives the GR answer**. The author's kill criterion is not
formally met (no published ≥3 saddles with 5° intrinsic EVPAs and eigen-axes), but every relevant
dataset points the same way. The strongest is B0218+357: the polarisation-angle *variations* of the saddle
image track the minimum image with the same sign at the 11-d lens delay, and the mirror requires the
opposite sign. B1422+231 and Q0957+561 disfavour the mirror at roughly the 2–3σ level, limited by
systematics. A cheap archival reanalysis would make the kill formal.

---

## 1. Internal consistency

**1a. Jordan algebra [arg].** Rank-2 Euclidean Jordan algebras are the spin factors, and their cones of
squares are Lorentz cones. Extreme rays ≅ S^{n−2}. That part is correct. Three errors:
- The linear automorphism group of the Lorentz cone is ℝ₊ × **O**⁺(1,n−1), which includes spatial
  reflections. SO⁺ is only its identity component. This matters because the Section-4 rule is exactly an
  improper (reflection) element. A3's complex structure on CP¹ is not preserved by reflections, since a
  reflection acts antiholomorphically and sends Chern number c → −c. So the theory has to *choose* how
  parity acts, and it never does.
- "S⁶ is only almost complex" is the open Hopf problem, not a theorem. Atiyah's claimed proof
  (arXiv:1610.09366) is not accepted. The selection n = 4 still survives on different grounds:
  H²(S^k; ℤ) ≠ 0 only for k = 2, so S² is the only ray sphere that admits non-trivial U(1) bundles.
  The stated argument is wrong but redundant.
- "Herm(2,ℂ) cone = light cone" is textbook; deriving 3+1 dimensions and Lorentz from the qubit cone is
  prior work (§4).

**1b. "Gravity changes only the phase function" (A6) implies GR polarisation [arg, decisive].** If the
mass label only alters the scalar phase φ, and amplitudes superpose (A4), then the observed field is
E_obs(ω) = F(ω)·E_src with the scalar amplification factor F = ∫ e^{iωφ(θ)} d²θ. Stationary phase
gives F = Σ_j |μ_j|^{1/2} e^{iωT_j − iπn_j/2}. Every image carries the **same** Jones vector up to a
helicity-independent scalar, including the Morse phase. That is exactly GR's χ_sad = χ_min and
V_sad/I_sad = V_min/I_min. The reflection appears only through A5's separate rule, "frame = polar factor
of ∂β/∂θ". That rule assigns gravity a second, non-phase action, which contradicts A6. The author's
weak point 2 concedes that the rule is imposed. I go further: it is **inconsistent** with A6.

**1c. The transport rule is not canonical even on its own terms [arg].** For det A < 0 the polar factor
is canonical in O(2), but the choice of O(2) over SO(2) is put in by hand. In 2-D every rotation is
exactly equidistant from diag(1, −1) (‖D − R_φ‖² = 4 for all φ), so "nearest rotation" is undefined,
and "nearest orthogonal map" is simply the decision to allow reflections. The lens map A also acts on
**angular displacements** between neighbouring rays, which is why image *shapes* are mirrored (a fact
known since Kemball et al. 1999: "source structural position angles may be changed"). The polarisation
of a single ray is a different object, and no axiom ties the two together. Near a fold caustic a
smooth wave field cannot interpolate between the two components of O(2), so the rule also fails in wave
optics.

## 2. Conservation, stability, causality

- **Angular momentum [arg].** Item 4 ("lens absorbs 2ħ per photon like a half-wave plate") needs a
  birefringent coupling, meaning a polarisation-dependent optical path difference of λ/2. A static
  weak-field metric acts as an isotropic, impedance-matched medium (ε = μ ≈ 1 − 2Φ/c²), so it has no
  birefringence at any order in geometric optics. Where a static lens's field is axisymmetric about the
  saddle ray's lens axis, J_z conservation would force an l = ∓2 orbital vortex onto the saddle image to
  compensate ΔS_z = ±2ħ. A scalar Fermat phase (A6) has no such winding. The torque must therefore go
  either to the lens, which A6 forbids ("no local vertex"), or to an OAM imprint, which A6's scalar phase
  cannot produce. Either way A6 is contradicted.
- **GR [lit].** Polarisation is parallel-transported along null geodesics. Dyer & Shaver (1992) showed
  with symmetry arguments that astrophysical lenses induce no rotation, except for extreme rotators.
  Faraoni (1993, astro-ph/9211012) confirmed this post-Newtonianly. The only helicity-dependent GR
  effect is the gravitational spin Hall effect (Oancea et al. 2020, arXiv:2003.04553). It scales as
  ~ λ/b times the bending angle. For λ = 3.6 cm, b ~ 1 kpc and α ~ 10⁻⁵ that gives ~10⁻²⁶ rad [calc],
  some 25 orders of magnitude below an O(1) flip.
- **Causality.** Nothing new beyond SR. c_gw = c follows from any single-cone theory, so it is not a
  prediction.

## 3. Already excluded? (key step)

**Parity assignments.** The minimum image leads in time.
- B0218+357: B lags A by 11.3 d [lit, Biggs & Browne 2018], so **B is the saddle**.
- Q0957+561: B lags A, so **B is the saddle**.
- B1422+231: **A, C are minima and B, D saddles** (cusp; standard labelling, as in Mao & Schneider 1998
  and arXiv:1409.1326). Some papers swap A and B, so check the labels.
- PKS 1830−211: NE leads, so **SW is the saddle** [lit, 1810.11012].

For a near-isothermal lens the saddle's positive eigen-axis ê₊ is radial (lens → image).

**(a) B0218+357: variability sign test (the strongest).** Biggs & Browne (2018, arXiv:1802.10088)
reanalysed 1996/97 VLA monitoring at 15 and 8.4 GHz. They found polarisation-angle "variations common
to both A and B", and plotted both images with the *same* ΔPPA scale. Pearson r(τ) for the PPA curves
has its global maximum near 11 d. The PPA-only CCF delays are 11.4 d (15 GHz) and 9.6 d (8.4 GHz).
Under R2-B, χ_B(t) = 2θ₊ − χ_A(t − τ), so dχ_B/dt = −dχ_A/dt and r should reach a **minimum**
(anti-correlation) at the lens delay. A constant RM, a constant lens screen and the unknown θ₊ all drop
out, so this test has no free parameter.
- Caveats: the authors call the PPA curves the weakest delay constraint ("shallow r curves"). No r value
  is published for the PPA alone, only the sign and location of the peak.
- Mimics: lens-plane RM variations would not correlate at the lens delay. Source-internal RM changes are
  intrinsic changes, which the mirror would also flip.
- Kemball, Patnaik & Porcas (1999) add that the A and B cores have "parallel PAs … at 43 GHz (where
  Faraday rotation is negligible)". They give no numbers. The mirror is consistent with that only if
  χ ≡ θ₊ (mod 90°) by coincidence.
- **Status: contradicts R2-B for one saddle. Qualitative and published, not yet quantified.**

**(b) B1422+231.** Patnaik et al. (1999, astro-ph/9909329), VLA from 1.4 to 43 GHz, give RM-corrected
intrinsic χ_A = 90±10°, χ_B = 57±10° and χ_C = 59±10°. The RMs are −4230, −3440 and −3340 rad m⁻².
The CASTLES lens position puts the radial axis of B at PA 48.5° or 131.5° (the RA sign convention is
unstated). By cusp symmetry, B's Jacobian eigen-axes are radial and tangential, and since A–B–C is
stretched tangentially, ê₊ is radial [calc]. Results:

| comparison | GR residual | mirror residual (θ₊ = 48.5° / 131.5°) |
|---|---|---|
| B vs C | 2° | 19° / 33° |
| B vs A | 33° | 50° / 64° |
| B vs mean(A,C) = 74.5° | 17.5° | 34.5° / 48.5° |

The two minima disagree by 31±14°, which both theories forbid, so ~30° systematics exist. The mirror
is disfavoured at about 2σ. The kill criterion is not met: errors are ±10°, and |χ_C − θ₊| is only
10.5° or 17.5°.

**(c) Q0957+561.** Patnaik et al. (1999) report RM_A = −61, RM_B = −91 rad m⁻² "with equal intrinsic
PAs". No values are given, so this cannot be scored without the archival data. Optical polarimetry
(Popović et al. 2021, arXiv:2101.07154) at a single epoch, with the 417-d delay not removed, gives
χ_A = 130±3° (g) and 117±3° (r), and χ_B = 153±2° and 151±2°. The radial axis of B from G1 is at
PA 9.8° or 170.2° [calc].
GR residuals are 23° (g) and 34° (r); mirror residuals are 57–83° [calc]. For the mirror to fit, ê₊
would have to sit 34–41° from radial. The authors attribute the residuals to extended-source or
microlensing effects. Disfavours the mirror; not decisive.

**(d) Others.**
- PKS 1830−211: ALMA image polarisations vary wildly by epoch (~5% vs ~0%; 1810.11012). Not usable
  without delay-matched epochs. MG J0414 is unpolarised (0.2%); B1600+434 has RMs only.
- **No parity-resolved Stokes V** for any lensed quasar was found, so the circular half of P1 is
  untested.

**(e) FRBs.** The CHIME/FRB lens-interferometry search (Kader et al. 2022, arXiv:2204.06014; Leung et
al. 2022, arXiv:2204.06001) correlates the X and Y feeds separately. It *requires* ε_X = ε_Y with the
same sign (Table II, the "polarization condition").
A reflection R(θ₊) acting on unpolarised voltages gives ε_X ∝ +cos 2θ₊ and ε_Y ∝ −cos 2θ₊ (tr R = 0)
[arg], so an R2-B echo is **rejected by construction**. P2's claim that existing searches would miss it
is correct. No lensed FRB has been confirmed anyway.

**Kill criterion status:** **not formally met.** Its precision demands (5° errors, ≥3 saddles, model
θ₊) have never been published together. The B0218 sign test is a stronger, parameter-free kill: reanalyse public VLA
programmes AH593/AB809 and correlate χ_B(t) with ±χ_A(t − τ).

## 4. Renaming check

- Lorentz group from the qubit cone, and 3-D space from the bit: Höhn & Müller (1412.8462); Müller &
  Masanes (1206.0630); see also 1004.1483 and 0911.0695.
- Helicity as monopole charge on CP¹: Wigner little group, Penrose spinors.
- Positions as duals of momenta: relative locality (1101.0931).
- The scalar lensing sector (Morse sign, Etherington, Δd/d = 0) is standard (astro-ph/0407232;
  2008.12814).
- The parity-dependent polarisation mirror itself is **not in the literature**. It is new only because GR
  forbids it, and it is inconsistent with A6 (§1b).

## 5. Novelty per prediction

- **P1-linear:** novel, but contradicted in sign by B0218 and disfavoured by B1422 and Q0957.
- **P1-circular:** novel and untested.
- **P2:** a novel test template that existing pipelines would veto.
- **P3:** novel, but partly degenerate with the GR Morse phase; no confirmed lensed GW pair (2304.08393).
- **Point-lens V/I ∝ 1/μ:** novel, no feasible target.
- **Scalar sector:** GR.

## 6. Triage

| Prediction | Verdict | Data | Mimics | Sensitivity |
|---|---|---|---|---|
| P1-linear: χ_sad = 2θ₊ − χ_min | **FAIL** (disfavoured; formal kill pending a 1-day reanalysis) | B0218 VLA monitoring (AH593/AB809); B1422 VLA; Q0957 | Faraday rotation (∝ λ², removable); differential depolarisation (image A in B0218); source structure at VLA resolution (~30° minimum–minimum scatter in B1422); time delay vs variability | Sign test is parameter-free. Static test needs ≤5° intrinsic χ from VLBI at ≥43 GHz. |
| P1-circular: V_sad = −V_min | **Conditional PASS** | Full-Stokes VLBA/ALMA of B0218, PKS 1830 (SW saddle), B1422 at delay-matched epochs | Faraday conversion (steep in λ); instrumental V leakage ~0.1%, equal to the signal; source V variability over the 11–26 d delays | Needs V/I ≳ 0.3% at ≥5σ per image. Marginal today. Kill only if (a) is not already decisive. |
| P2: FRB echo with flipped V and mirrored χ | **Conditional PASS** (cheap add-on, low yield) | CHIME/FRB baseband (172+ bursts), CHIME Outriggers | Scintillation and multipath (no handedness flip); instrumental reflections <300 ns; Faraday rotation within the burst | Add an X↔−Y template to Kader et al. Expected lens rate is tiny under existing PBH limits, so a null result is uninformative. |
| P3: GW ι₂ = π − ι₁ | **FAIL** (now) / conditional for O5 | LVK O3/O4 lensing candidates | Degenerate with Morse phase and polarisation angle for 22-only signals | Needs higher modes or precession and a confirmed lensed pair. None exists. |
| Point-lens V/I ∝ 1/μ | **FAIL** (no feasible target) | none | Limb polarisation; finite source | Stellar V ≪ 10⁻³. No pulsar or maser microlensing with polarimetry. |
| Cross-image anti-diagonal coherence | Same as P2 | — | — | — |

**Bottom line.** Hand to the pipeline only one task, a reanalysis rather than a search: the B0218+357
χ_A/χ_B sign-correlation test on public VLA data. If r(+) > r(−) at τ = 11.3 d at >3σ, R2-B is dead.

---

## Verified references

Each was opened on 2026-10-09. arXiv IDs were checked against the abstract-page metadata.

1. Dyer, C. C. & Shaver, E. G. 1992, ApJ 390, L5. Not on arXiv. Bibliographic data verified from
   Faraoni's reference list; the ADS page was not reachable (HTTP 405). The title is as given in the
   brief and was not independently seen.
2. Faraoni, V. 1993, A&A 272, 385, "On the rotation of polarization by a gravitational lens",
   arXiv:astro-ph/9211012.
3. Kemball, A. J., Patnaik, A. R. & Porcas, R. W. 1999, "VLBI Polarisation Images of the Gravitational
   Lens B0218+357", arXiv:astro-ph/9909330 (Boston Univ. lensing conference proceedings).
4. Patnaik, A. R., Menten, K. M., Porcas, R. W. & Kemball, A. J. 1999, "Polarisation Observations of
   Gravitational Lenses", arXiv:astro-ph/9909329.
5. Biggs, A. D. & Browne, I. W. A. 2018, "A revised lens time delay for JVAS B0218+357 from a
   reanalysis of VLA monitoring data", MNRAS, arXiv:1802.10088.
6. Popović, L. Č. et al. 2021, "Spectroscopy and polarimetry of the gravitationally lensed quasar
   Q0957+561", arXiv:2101.07154.
7. Martí-Vidal, I. & Muller, S. 2019, "Submillimeter polarization and variability of quasar
   PKS 1830-211", A&A 621, A18, arXiv:1810.11012.
8. Ros, E. et al. 2000, "VLBI imaging of the gravitational lenses B1422+231 and MGJ0414+0534",
   arXiv:astro-ph/0010650. Images A–B–C stretched along their joining direction.
9. CASTLES lens database, B1422+231, Q0957+561 and Q2237 pages
   (lweb.cfa.harvard.edu/castles/Individual/), accessed 2026-10-09. The RA sign convention is unstated.
10. Kader, Z., Leung, C. et al. 2022, PRD 106, 043016, arXiv:2204.06014 (CHIME/FRB lens
    interferometry; polarisation condition ε_X = ε_Y).
11. Leung, C., Kader, Z. et al. 2022, PRD 106, 043017, arXiv:2204.06001.
12. Oancea, M. A. et al. 2020, "Gravitational spin Hall effect of light", arXiv:2003.04553.
13. Ezquiaga, J. M., Holz, D. E. et al. 2021, "Phase effects from strong gravitational lensing of
    gravitational waves", arXiv:2008.12814.
14. LIGO-Virgo-KAGRA 2023, "Search for gravitational-lensing signatures in the full third observing
    run of the LIGO-Virgo network", arXiv:2304.08393.
15. Müller, M. P. & Masanes, Ll. 2013, "Three-dimensionality of space and the quantum bit",
    arXiv:1206.0630.
16. Höhn, P. A. & Müller, M. P. 2016, "An operational approach to spacetime symmetries: Lorentz
    transformations from quantum communication", arXiv:1412.8462.
17. Masanes, Ll. & Müller, M. P. 2011, "A derivation of quantum theory from physical requirements",
    arXiv:1004.1483.
18. Dakić, B. & Brukner, Č. 2009, "Quantum Theory and Beyond: Is Entanglement Special?",
    arXiv:0911.0695.
19. Amelino-Camelia, G., Freidel, L. et al. 2011, "The principle of relative locality",
    arXiv:1101.0931.
20. Atiyah, M. 2016, "The Non-Existent Complex 6-Sphere", arXiv:1610.09366 (claimed proof, not
    accepted).
21. Kochanek, C. S. 2004, "The Saas Fee Lectures on Strong Gravitational Lensing",
    arXiv:astro-ph/0407232. Seen in search results only; cited for standard image-parity facts.
22. arXiv:1409.1326 ("Constraints on warm dark matter from weak lensing in anomalous quadruple
    lenses"). B1422 parity A, C minima and B, D saddles, seen via a search snippet only.

Not verified, and used only as background: Mao & Schneider 1998 (B1422 parity); Jordan, von Neumann
& Wigner 1934, Ann. Math. 35, 29; Greenfield et al. 1985 (Q0957 RM; not found online); Kronberg et al.
1991 (only seen cited as assuming lensing preserves polarisation).
