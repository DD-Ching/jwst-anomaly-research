# System B on R3-C "Resonant Partner Density"

*Labels: [pub] published (arXiv ID verified 2026-10-09); [calc] my calculation; [est] estimate; [data] measured by me today.*

## Verdict in one paragraph

Not excluded, but the headline test cannot work. (1) A2 makes λ a *kinematic shell-crossing* variable, not a potential depth: −Φ/c² is the wrong regressor and infallers predict **zero** effect. (2) A5 leaves nothing for A3 to gate. (3) Redoing the white-dwarf slope with *central* clock depths gives K_Ia ≈ −2.4×10³ to −4×10³ [calc]., the band's weak edge, where the kill cannot fire. (4) The published environmental SN Ia rate anomalies range from ×0.6 to ×3–8 [pub], 10–40× the predicted −6 % to −20 %. Unadjusted, the record points to the *positive* K that the author's own criterion counts as a kill. Novelty: a renamed khronometric preferred-frame rate coupling.

## 1. Internal consistency

**1a. Is λ well defined, and is it Lorentz invariant?** τ_max(x) = sup over the past of the Lorentzian distance. This is the Andersson–Galloway–Howard *cosmological time function* [pub, gr-qc/9709084]. It is locally Lipschitz but *not* C¹ everywhere, and every point is reached by a maximising "rest curve" (a timelike geodesic). Where τ_max is differentiable, g(∇τ_max, ∇τ_max) = −1, so dτ_max/dτ = −u·∇τ_max = γ_rel and

  **λ = 1 − 1/γ_rel ≈ |v − v_rest|²/2c²** [calc]

exactly, with no Φ term. λ is a scalar, but it is the Lorentz factor relative to a hypersurface-orthogonal unit-timelike field. — a khronon [pub, 0909.3525]. The −Φ/c² form comes back only for matter that is *not* on a rest curve, for example a static white-dwarf centre reached by a rest curve that fell in from outside (v_rest² /2 ≈ Φ_out − Φ). Three consequences:

- **A first-infall galaxy is itself close to a rest curve, so λ_ext ≈ 0.** λ_ext becomes non-zero only after orbit crossing, when geodesics stop maximising. So the "infallers at 50–70 % of central depth" lever arm (§6, item 2) predicts *no* signal under A2. λ_ext effectively tracks time since accretion, the same clock that drives quenching and the DTD age effect: the degeneracy is built into the definition.
- **The v²/2c² term is not "negligible at 10⁻⁶".** For a virialised satellite, |v − v_rest|² ~ σ² with σ ≈ 300–1000 km/s, giving 0.5–5.6×10⁻⁶. That is the same size as the group Δλ in the R3-C table (2–8×10⁻⁶). The cross term 2v·v_rest/c² for radial-in versus radial-out cluster members is ~4×10⁻⁵ at 1000 km/s [calc], larger than the cluster-core Δλ of 3×10⁻⁵.
- **∇τ_max at Earth is the incoming rest-curve direction, not the CMB frame**, because the Solar System and Milky Way are bound and phase-mixed. So the "CMB-dipole phase" of the annual modulation is not derived. The amplitude scale is v_⊕·v_rest/c² with an unknown direction.

**1b. The single-quantum exemption (A5).** Every collective runaway is made entirely of exempt parts: individual reactions (A5) plus continuous evolution (A5). Carbon-ignition hot spots, recoil-nucleated bubbles and Geiger discharges are all single quanta amplified by deterministic dynamics. The gate has no physical locus.

Worse, "links only change the timing of the release" (§2) contradicts a −20 % rate change. A pure delay Δt in a population with DTD ∝ t^−1.1 gives R_obs(t) = DTD(t − Δt), which is an *enhancement* in old populations and of order Δt/t [calc]. A suppression needs runaways that permanently fizzle, so the white dwarf never explodes. That is an energy-release change, not a timing change.

**1c. Annual modulation already excluded?** Not cleanly. Precise decay nulls use non-runaway detectors, which A5 exempts: HPGe at Gran Sasso, 3.4–3.5×10⁻⁶ (1σ) [pub, 1802.09373]; TDCR liquid scintillation [pub, 1407.2493]; Norman et al. [pub, 0810.3265]. The claimed modulations (BNL gas counter, PTB ionisation chamber [pub, 0808.3283]; Parkhomov's *Geiger–Müller* ⁹⁰Sr series, disproved by LSC [pub, 1407.2493]) carry documented humidity artefacts [pub, 2108.00116]. So GM data neither confirm nor exclude ≲10⁻³. I found no annual-modulation analysis of bubble-chamber nucleation rates. With sign (λ_LSS,here) and phase (1a) both free, kill (c) cannot fire.

**1d. w = A_s^{1/2}.** The number is right: Planck A_s ≈ 2.1×10⁻⁹ [pub, 1807.06209] gives 4.58×10⁻⁵. The identification fails inside the framework it cites, though:

- In a Poisson sprinkling, longest chains track geodesics [pub, gr-qc/0512073].
- Longest-chain fluctuations scale sub-linearly: ∝ L^{1/3} in the 2-D (Ulam) case [pub, math/9810105].
- With L ~ t₀/t_P ≈ 8×10⁶⁰, the fractional fluctuation is ~L^{−2/3} ≈ 2.5×10⁻⁴¹ [calc]. The 2-D exponent is used as an indicator; d = 4 is not proven.

Off by ~36 orders of magnitude. Bombelli–Henson–Sorkin [pub, gr-qc/0605006] show that sprinkling picks out no frame, so a frame-dependent λ cannot come from the discreteness itself.

**1e. K_Ia ≈ −8×10³.** The author mixes surface and central clock depths. I used Nauenberg radii (μ_e = 2) and polytropic centre-to-surface potential ratios (2.35 for n = 1.5, 4.42 for n = 3, interpolated) [calc, est]:

| M (M_☉) | R (km) | λ_surf | λ_centre |
|---|---|---|---|
| 0.6 | 8750 | 1.0×10⁻⁴ | 2.7×10⁻⁴ |
| 1.0 | 5560 | 2.7×10⁻⁴ | 8.8×10⁻⁴ |
| 1.2 | 3970 | 4.5×10⁻⁴ | 1.7×10⁻³ |
| 1.37 | 2200 | 9.2×10⁻⁴ | 3.8×10⁻³ |

- **Central values.** The trigger sits at the centre (§4), where Δλ(0.6→1.2) = 1.4×10⁻³, not 4×10⁻⁴.
- **Jacobian.** ρ_b(λ) ∝ n(M)·M·dM/dλ adds a factor 0.28 between 1.2 and 0.6 M_☉ [calc]. A 30× fall in n(M) therefore becomes ~107× in ρ_b.
- **Result.** K_centre ≈ −ln(30…300)/1.4×10⁻³ = **−2.4×10³ to −4.1×10³** [calc].
- **Observed mass function.** The input fall-off is itself uncertain. The 100 pc sample shows a narrow peak at 0.59 M_☉, a 0.7–0.9 M_☉ shoulder and a near absence of cool DA white dwarfs above 1 M_☉ [pub, 2006.00323]. 5–9 % of massive white dwarfs are cooling-delayed by ²²Ne distillation [pub, 2407.04827]. A luminosity-selected snapshot is not the record-mass census that A6 requires.
- **K is not one number.** It varies with λ*, from ~9×10⁻⁴ for sub-Chandrasekhar to ~4×10⁻³ for near-Chandrasekhar ignitions. So R3-C *predicts* that the channel mix, and with it the x₁ distribution, shifts with Δλ_ext. This contradicts side prediction (i) "rates shift, properties don't".
- **Partner count:** ~10⁷⁵–10⁷⁶ in the window [est], not 10⁵⁰–10⁶⁰ (harmless, but unnormalised).

## 2. Already excluded?

| Study | Result | Relevance |
|---|---|---|
| Mannucci+08 [0710.1094] | Cluster early-type SNuM 0.066 vs field 0.019 (×3.5, 98 % CL); cluster CC rate ≈ field | Opposite sign, ~2σ |
| Dilday+10 SDSS-II [1003.1521] | Cluster/field early-type ratio 1.94 (+1.31/−0.91) for C4 and 3.02 (+1.31/−1.03) for maxBCG; hint of enhancement in cores | A ratio of 0.8 is ~2σ low vs maxBCG |
| Sharon+07 [astro-ph/0610228]; Sand+12 MENeaCS [1110.1632]; Graham+08 SNLS [0801.4968] | Cluster or red-sequence rate per mass ≈ field ellipticals (±30–50 %) | Cannot resolve −20 % |
| Friedmann & Maoz 18 [1803.04421]; Freundlich & Maoz 21 [2012.00793] | Cluster DTD amplitude 2–3× field (3.8σ), even with extended star-formation histories | +100–200 % "environment effect" |
| Toy+23 DES [2302.05184] | Per-galaxy cluster/field ratio 0.594 ± 0.068 at 10 ≤ log M* ≤ 11.25, attributed to age. Passive SNuM 0.0386 cluster vs 0.0625 field (≤1.7σ) | Only result on R3-C's side, and age-explained |
| Cooper+09 [0901.4338] | Star-forming hosts: SNe Ia favour *lower* density at fixed M* and SFR; attributed to metallicity | Mimic with R3-C's sign |
| Lordet+26 ZTF DR2 [2604.12714] | 12–14σ excess of SNe Ia vs uniform rate; rates ×3–8 in specific clusters (Perseus, Coma, Hercules); not a linear tracer of 2M++/Manticore density | Closest thing to a potential regression; strongly positive |
| McGee & Balogh 10 [0912.3455] | 19/59 group SNe Ia hostless | Hostless events have no host covariates |

**Satellite vs central and halo-mass regressions.** I found none on arXiv or the web. **Regression on potential depth:** none. The nearest is Lordet+26, which uses the density field.

**Conclusion.** The predicted −6 % to −20 % is *not formally excluded*, because every rate comparison carries 30–50 % errors and an unmodelled astrophysical baseline. However:

- Measured "environment effects" span ×0.6 to ×8 between studies. The systematic floor is a few % per 5×10⁻⁶ of Δλ, which already means |K_sys| ≳ 10⁴ > |K_pred| [calc].
- The un-adjusted record (Lordet, Friedmann/Freundlich, Mannucci, Dilday) gives an apparent positive K of ~+5×10⁴ to +10⁵ [calc: ln 3–8 over Δλ ≈ 1–3×10⁻⁵]. The author's own kill ("positive K at ≥3σ") fires unless confounders are invoked; with them, the test has no power (§5).

**Core collapse.** Graham+12: 7 cluster SNe II, 0.026 (+0.085/−0.018) SNuM [1205.0015]; Mannucci+08: cluster ≈ field. Known to ~×3 only.

## 3. Renaming check

- **Khronometric / Einstein-aether** [0909.3525; gr-qc/0410001]: τ_max *is* a khronon. "Derived, not postulated" changes nothing observable; matter couples to u·∇τ.
- **Mach:** generic causal-past dependence. **Chameleons** [1306.4326]: screened in deep wells, so the sign is reversed; superficial. **Fischbach–Jenkins** [0808.3283]: dodged by A5 fiat.

## 4. Novelty

The specific claim (SN rate log-slope in clock depth set by the white-dwarf mass function) is new as far as I found; the mechanism is a preferred-frame rate coupling. **Novel but ill-posed.**

## 5. Triage

| # | Prediction | Verdict |
|---|---|---|
| P1 | SN Ia rate ∝ exp(K·Δλ_ext) at fixed host, K ∈ [−3×10³, −3×10⁴] | **Conditional PASS** (not excluded; public data; underpowered) |
| P2 | Cluster cores −20 % | **FAIL** as a test: ~50–100 core SNe Ia available (Larison+23: 102 in X-ray clusters within 2r₅₀₀ [2306.01088]); the published baseline scatters ×0.6–×8 |
| P3 | Large-scale-potential term ∓4–8 % | **Conditional PASS**: the only lever arm partly orthogonal to age, but λ ≠ −Φ (1a) |
| P4 | Infaller lever arm breaks the age degeneracy | **FAIL**: A2 gives λ_ext ≈ 0 for first infall |
| P5 | Rates change, light curves don't | **FAIL**: K depends on λ*, so the channel mix shifts; x₁ is already strongly environment-dependent (inner-cluster ≥75 % fast decliners, Larison+23), so this cannot discriminate |
| P6 | K_CC/K_Ia ∈ [0.2, 5] | **FAIL**: λ*_CC is arbitrary; cluster CC rates are known only to ~×3 |
| P7 | K grows more negative with cosmic time | **FAIL**: unmeasurable |
| P8 | Lab annual modulation ≲10⁻³, CMB phase | **FAIL** as a falsifier: sign is free, phase not derived; not excluded (nulls use non-runaway detectors) |
| P9 | No event–event coincidences | Trivial PASS (no content) |

**Data for P1/P3:**

- **ZTF BTS public table** (sites.astro.caltech.edu/ztf/bts/explorer.php). Today's statistical-quality CSV is 0.87 MB with **7072 SNe: 5326 Ia (5020 at z < 0.1) and 1651 CC**, peaks from 2018-05 to 2026-08 [data]. Completeness is 97/93/75 % at <18/18.5/19 mag [pub, 2009.01242].
- ZTF SN Ia DR2: 3628 events, a ZTF subset [2409.04346]. ASAS-SN: 1776 SNe Ia from 2014–2024, mostly overlapping [2602.00223].
- TNS: no selection function, unusable for rates.
- **Realistic clean N_Ia ≈ 5–6×10³, not 2×10⁴.**
- **Hosts:** Legacy Survey SED fits. **Groups:** Lu+16 2MRS (z ≤ 0.08, whole sky, halo mass ±0.35 dex) [1607.03982]; Yang+07 [0707.4640] and Tempel+17 [1704.04477] (SDSS footprint only, roughly a third to half of BTS hosts [est]).
- **Potentials:** 2M++ [1105.6107] / Carrick+15 [1504.04627]; CF4 [2209.11238]. Total ≲300 MB.

**Mimics:** DTD/age via time since infall (which *is* λ_ext under A2); quenching; metallicity (Cooper+09); cluster DTD normalisation (×2–3); detection inefficiency on bright cluster hosts; hostless intracluster SNe (~10–30 %); host-mass and IMF M* biases; peculiar-velocity distance errors correlated with Φ_LSS; the unexplained ZTF cluster excess.

**Sensitivity** [calc]: σ_K = 1/(0.74 σ_Δλ √(1−R²) √N), with σ_Δλ = 5×10⁻⁶ and R² the collinearity of Δλ with the local-density covariate.

| N | R² | σ_K | Signif. at K = −8×10³ | Signif. at K = −3×10³ |
|---|---|---|---|---|
| 5×10³ | 0 | 3.8×10³ | 2.1σ | 0.8σ |
| 5×10³ | 0.6 | 6.0×10³ | 1.3σ | 0.5σ |
| 2×10⁴ | 0.6 | 3.0×10³ | 2.6σ | 1.0σ |

- For the main kill to fire on a null (95 % upper bound > −3×10³), N ≈ **8×10⁴** is needed.
- The 5σ detection at my K ≈ −3×10³ needs N ≈ 2×10⁵ to 5×10⁵, which means the LSST era with photometric classification.
- The systematic floor (|K_sys| ~ 10⁴) exceeds the signal at any N.

**Feasibility.** A Poisson GLM with N ≈ 5×10³ events over ~10⁵ census galaxies takes minutes on a laptop. The hard part is the per-host control time, which needs the BTS selection function. **Feasible to run; unable to confirm or kill R3-C.**

## Verified references (arXiv IDs checked against the arXiv API, 2026-10-09)

- Andersson, Galloway & Howard 1998, The Cosmological Time Function — gr-qc/9709084
- Blas, Pujolas & Sibiryakov 2009, A healthy extension of Hořava gravity — 0909.3525
- Eling, Jacobson & Mattingly 2004, Einstein-Aether Theory — gr-qc/0410001
- Khoury 2013, Chameleon Field Theories — 1306.4326
- Bombelli, Henson & Sorkin 2006, Discreteness without symmetry breaking — gr-qc/0605006
- Ilie, Thompson & Reid 2005, paths in causal sets vs geodesics — gr-qc/0512073
- Baik, Deift & Johansson 1998, longest increasing subsequence — math/9810105
- Planck 2018 VI — 1807.06209
- Jenkins et al. 2008 — 0808.3283; Norman et al. 2008 — 0810.3265; Kossert & Nähle 2014 — 1407.2493; Bellotti et al. 2018 — 1802.09373; Pommé & Pelczar 2021 — 2108.00116
- Kilic et al. 2020, 100 pc WD sample — 2006.00323; Jewett et al. 2024, massive WDs — 2407.04827
- Mannucci, Maoz & Sharon 2007/08 — 0710.1094; Sharon, Gal-Yam & Maoz 2006/07 — astro-ph/0610228; Dilday et al. 2010 — 1003.1521; Graham et al. 2008 — 0801.4968; Sand et al. 2011/12 — 1110.1632; Graham et al. 2012 — 1205.0015; Friedmann & Maoz 2018 — 1803.04421; Freundlich & Maoz 2020/21 — 2012.00793; Toy et al. 2023 — 2302.05184; Cooper, Newman & Yan 2009 — 0901.4338; McGee & Balogh 2009/10 — 0912.3455; Larison et al. 2023 — 2306.01088; Lordet et al. 2026 — 2604.12714
- Perley et al. 2020, BTS II — 2009.01242; Rigault et al. 2024, ZTF SN Ia DR2 — 2409.04346; Desai et al. 2026, ASAS-SN Ia rates III — 2602.00223
- Lu et al. 2016, 2MRS groups — 1607.03982; Yang et al. 2007 — 0707.4640; Tempel et al. 2017 — 1704.04477; Lavaux & Hudson 2011, 2M++ — 1105.6107; Carrick et al. 2015 — 1504.04627; Tully et al. 2022, CF4 — 2209.11238
- Data: ZTF BTS explorer CSV (statistical-quality SN sample), https://sites.astro.caltech.edu/ztf/bts/explorer.php, accessed 2026-10-09.
