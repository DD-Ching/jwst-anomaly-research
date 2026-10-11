# B on A3: "Bridged Markov substrate" — adversarial review

Reviewer: System B (Physics Destruction). Every reference was checked against the arXiv API or its abstract page on 2026-10-09. Labels: **[derived]** = my own calculation, **[obs]** = published measurement, **[est]** = order-of-magnitude estimate.

## 1. Internal consistency

**1.1 Complex amplitudes in a positive measure.** A1 and A2 define `P` as a strictly positive Markov random field. Hammersley–Clifford requires positive clique potentials, so `ψ_b = 1 + ε e^{iφ_b}…` with `φ_b ≠ 0` is not a probability weight. The whole `sinφ_b`/`cosφ_b` split, and with it `Q`, therefore has no definition in the axioms. A phase could only enter after a Wick rotation to a Lorentzian transfer operator, which needs reflection positivity.

**1.2 Massless emergent physics vs the Dobrushin regime.** A light cone is derived from Dobrushin uniqueness, where correlations decay as `λ^d`. In that regime every emergent field is massive, with correlation length `~1/|ln λ|` steps. A massless `□_g`, as needed for photons and gravitational waves, requires criticality (power-law correlations), which is the opposite regime. In addition, Gromov–Hausdorff convergence applies to metric spaces, so it cannot produce a Lorentzian `g`.

**1.3 The optical-theorem analogy undermines "universal Q at every wavelength" [derived].** For a dilute set of scatterers with forward amplitude `f(0)`:
- the extinction per unit length is `Nσ_ext = (4πN/k) Im f`;
- the refractive index is `n−1 = (2πN/k²) Re f`.

Achromatic dimming needs `Im f ∝ k`. An achromatic time delay (constant `n−1`) needs `Re f ∝ k²`. Then `cot φ_b = Re f / Im f ∝ k`, so `Q = −½cot φ_b` changes with frequency. From radio to optical it shifts by about 10⁴. From LIGO-band gravitational waves (λ ~ 10³ km) to optical it shifts by about 10¹³. If instead `φ_b` is held fixed, either the dimming or the delay becomes chromatic. That contradicts P1's "same β in every band" and the claim that `Q` is messenger-universal.

**1.4 "d_A is clean" contradicts a differential δ_T [derived].** P1 needs the bridge delay to differ between images, `Δt_res ∝ Δτ_ij`. Any extra delay field `t_b(θ)` on the sky enters Fermat's principle and deflects rays by `∇_θ t_b`. A differential delay at O(ε) therefore moves image positions at O(ε). In a lens fit this is absorbed as fake surface density `∝ ∇² t_b`, so most of the predicted delay residual disappears into the macro model. Saying there is "no centroid change at O(ε)" is inconsistent with what P1 predicts.

**1.5 The κ-linearity contradicts α ≳ 2 [derived].** With `n_b ∝ (1+δ)^α`, the column through a lens is `Δτ ∝ ∫ρ^α dl`, not `Σ`. For an isothermal profile (`ρ ∝ r⁻²`), `∫ρ² dl ∝ R⁻³` while `κ ∝ R⁻¹`, which gives `Δτ ∝ κ³`. The relation `Δm = βΔκ` therefore holds only for `α = 1`, and that value is the one the author's own CMB argument rules out.

**1.6 The CMB "shielding" argument is wrong.** Raising α does not put CMB photons outside collapsed structure: they cross the same z < 2 halos that SN light crosses. So `τ(1100) ≥ τ(z≈2) ≥ τ(1)` for any α. The fiducial SN dimming of 1–4% in flux therefore implies that CMB anisotropy power is damped by about 2–8%. That is degenerate with `A_s e^{−2τ}`. Planck measures `τ_reio = 0.054 ± 0.007` from low-ℓ polarisation (1807.06209), and this extra damping, unlike reionisation, acts at all ℓ. The inferred `A_s` and σ₈ would drop by roughly 1–4% **[est]**. That pressure on the upper fiducial range comes from CMB-lensing and A_L consistency checks; it is pressure, not an exclusion.

**1.7 Mass-sheet sign argument, checked.** Under `κ → λκ + (1−λ)`, delays scale as `Δt → λΔt` and magnifications as `μ → μ/λ²` (1306.0901, §2). A3 states this correctly: λ < 1 gives shorter delays and brighter images. However:
- (a) `μ_i/μ_j` and `Δt_ij/Δt_ik` are invariant under a mass-sheet transform, so in quads (relative fluxes only) the transform has no flux signature to be "opposite" to. The sign test needs absolute magnifications, which only lensed SNe Ia provide.
- (b) The sign of `δ_T` is set by the free parameter `φ_b`. For `sinφ_b < 0` the bridge signal has the same sign as a mass-sheet transform, so the test kills only half of parameter space.
- (c) TDCOSMO marginalises an internal mass sheet `λ_int` per lens, and a uniform `δ_T` per lens is fully degenerate with it.

**1.8 Smaller errors.**
- P1 writes `Δt_res/Δt ≈ (2δ_T/δ_L)`, but `1/Q = δ_T/δ_L`. The factor of 2 is unexplained.
- In P4, `Γ` is a width in `Ω ∝ ln ρ` units, so `Γ⁻¹` is not a time. "Lorentzian in density" is really Lorentzian in ln ρ.

## 2. Conservation, stability, causality

- **Local non-conservation is incompatible with any metric gravity.** Where flux "leaves at one endpoint and reappears at the other", the Bianchi identity (`∇_μG^{μν} = 0 ⇒ ∇_μT^{μν} = 0`) fails. Bridge endpoints act as time-variable mass sources, so the "static" potentials in P3 are not static.
- **The no-signalling rule does not work.** "Norm but no controllable phase" still allows signalling by amplitude: if a patch can be pushed into or out of the resonance set R, then switching norm transfer on and off sends bits across `ℓ_b`. Frozen bridges carry energy across spacelike separations. The preferred foliation avoids paradoxes only by breaking Lorentz invariance, and the author then has to evade Lorentz-violation bounds, which are not addressed.

## 3. Already excluded? (numbers)

| Prediction | Existing constraint | Status |
|---|---|---|
| Mean η_E > 0 (δ_L ≈ 0.005–0.02 at z=1) | Avgoustidis+2010: ε = −0.04 +0.08/−0.07 (2σ) for d_L = d_A(1+z)^{2+ε}, which gives Δln d_L < 0.03 at z=1 **[derived]**. Martinelli+2020 tighten current data by ×2.5. | Upper fiducial end pressured; lower end allowed |
| GW universality | GW170817 (1801.08160): d_L^GW/d_L^EM exponent γ = 1.01 +0.04/−0.05 at 40 Mpc, where τ ≈ 0. GWTC-3: H₀ = 68 +8/−6 (≈10%). Speed: \|c_g−c\|/c < 3×10⁻¹⁵ (1710.05834). | Not testable at 1–4%. Universality makes Ξ₀ = 1 trivially, so EM–GW tests are blind to it. |
| Time-delay residuals | TDCOSMO-2025: 8 lenses, H₀ = 71.6 +3.9/−3.3, with the mass sheet marginalised. TDCOSMO I: the 7 lenses are mutually consistent. | Degenerate with λ_int; not excluded |
| P1 quads | Keeley+2024: 9 MIRI + 5 narrow-line quads; the survey targets 31 (2309.10101), with 3% flux-ratio precision; Gilman+2020: 8 quads | About 14 analysed quads are public; ≥30 is not yet available |
| Lensed SNe Ia | iPTF16geu, SN Zwicky (θ_E = 0.167″, delays < 1 day), SN H0pe (H₀ = 75.7 +8.1/−5.5, absolute magnification used) | n ≈ 2–3 usable systems |
| SN lensing | Smith+2014: 608 SNe, 1.4σ; Shah+2024 (DES-5YR): lensing magnification detected at 6.0σ | An opposite-sign flat-kernel term would partly cancel this; not yet quantified |
| P3 static deficit | Applegate+2016 (relaxed): M_lens/M_HSE = 0.96 ± 9% ± 9% at r₂₅₀₀. WtG: M_Planck/M_WL = 0.688 ± 0.072. CLASH-VLT: η(r₂₀₀) = 1.01 +0.31/−0.28. ESO 325-G004: γ_PPN = 0.97 ± 0.09 at kpc scales. Chameleon tests: Terukina+2014, Wilcox+2015 (58 clusters). | Not excluded. The r² scaling predicts only 2–4% at r₂₅₀₀ ≈ 0.45 r₅₀₀ **[derived]**, which Applegate cannot resolve. The kill criterion (±5% at r₅₀₀) is beyond current precision. |

## 4. Renaming check

- **Bridges** are "disordered locality": nonlocal links in an emergent-geometry graph (Markopoulou & Smolin 2007). Their cosmology was worked out by Prescod-Weinstein & Smolin (2009). "Shortcut edges" is small-world networks (Watts & Strogatz 1998).
- **δ_L** (photon loss into an isotropic background) is the cosmic-opacity / distance-duality framework (Avgoustidis+2010), with grey dust (Aguirre 1998) and photon–axion mixing (Csáki, Kaloper & Terning 2001) as concrete mechanisms.
- **The messenger-universal version** is GW "leakage" (Deffayet & Menou 2007; Pardo+2018) applied to every field.
- **P3** is scale-dependent gravitational slip Ψ/Φ ≠ 1 with a Yukawa-like range ℓ_b. The parametrisation and the test design (static vs lensing mass) already exist: Pizzuti+2016, Wilcox+2015, Terukina+2014. The only new element is the sign: a deficit, where chameleon models predict an excess.
- **The correlation of delay with κ** is degenerate with the mass sheet (A3 admits this).

## 5. Novelty per prediction

| Prediction | Novelty |
|---|---|
| P1 Δm = βΔκ, parity- and band-blind | Partially known: flux-ratio anomaly vs κ/parity studies exist; linking it to a delay residual is new |
| Universal Q = δ_L/δ_T | Novel as a statistic, but ill-defined (§1.1, 1.3) |
| P2a messenger-universal η_E | Partially known (opacity + GW leakage) |
| P2b flat-kernel vs lensing-kernel SN regression | Novel as a test (no standard pipeline builds a ρ^α-weighted column) |
| P3 aperture² static deficit | Known (gravitational slip / Yukawa) |
| P4 (1+z_l)-scaled achromatic dips | Novel phrasing, but ill-defined |

## 6. Survivor triage

| Prediction | Verdict | Reason |
|---|---|---|
| P1 Δm = βΔκ per image | **FAIL** | Inconsistent with α ≳ 2 (§1.5) and with "clean d_A" (§1.4). Within quads, κ correlates with parity (saddles sit at slightly higher κ), so the known demagnification of saddles by substructure already produces β > 0. Sensitivity **[est]**: the residual scatter after the macro model is ≈ 0.1–0.2 mag per pair and σ(Δκ) ≈ 0.1. With 30 quads × 3 pairs, σ_β ≈ 0.15 / (√90 × 0.1) ≈ 0.16 mag per unit κ, so reaching 0.02 needs about 2,000 quads. Only a weak limit (≈ 0.3 mag at 2σ) is possible now. |
| Universal Q | **FAIL** | Undefined in the axioms (§1.1). Not universal across wavelength (§1.3). Unmeasurable: flux ratios and delay ratios are mass-sheet invariant, and absolute magnifications exist for only about 3 SNe with σ_int ≈ 0.1 mag plus microlensing. |
| P2a GW = EM η_E | **FAIL** (not testable now) | Bright and dark sirens give ≳ 10% at z > 0.05, while 1–4% is needed. Needs ET/LISA (cf. 1712.08108). |
| P2b flat-kernel SN regression | **PASS (conditional)** | Testable now. See below. |
| P3 aperture² deficit | **FAIL** | A renamed scale-dependent slip/Yukawa (§4). The sign is the only new element, and it is fully mimicked by non-thermal pressure plus WL calibration bias. Existing tests are the right way to examine it, and they find nothing. |
| P4 transient dips | **FAIL** | Duration is dimensionally undefined (§1.8). Microlensing Einstein times also scale with (1+z_l), and large-source caustic crossings are nearly achromatic. The flash–dip fluence balance cannot be tested. The no-signalling failure (§2) makes the mechanism acausal. |

### P2b: what a search would have to do

**Mimics to exclude:**
1. Weak-lensing magnification. Model it with a halo model, as in Shah+2024, rather than with a single κ proxy.
2. Grey or weakly chromatic dust in foreground halos. Fit colour excess and SALT c per SN against T_flat.
3. Selection effects. SNe on dimmed sightlines drop out near the detection limit, which biases γ toward 0. Use simulated bias corrections.
4. Photometry systematics from blending and host contamination on lines of sight with bright foreground galaxies.
5. Host-environment correlations at low z, where the "foreground" overlaps the host group. Restrict to z > 0.2.
6. Photo-z scatter in the foreground catalogue. This dilutes T_flat and κ unequally.

**Data:**
- DES-SN5YR public distances plus DES Y3 Gold (exactly Shah+2024's setup);
- Pantheon+ public tables (2112.03863) with Legacy Survey photo-z galaxies for the foreground;
- DESI DR2 BAO (2503.14738) for the mean-η_E anchor.

**Sensitivity [est]:**
- With N ≈ 1,500 SNe at σ ≈ 0.17 mag, `σ_γ ≈ 0.17/√1500 ≈ 0.0044` mag per unit normalised column. Inflating by 1/√(1−ρ²) ≈ 1.7 for collinearity gives ≈ 0.007.
- With α ≳ 2, the sightline-to-sightline scatter of τ is comparable to its mean, so γ ≈ 0.01–0.04 mag per unit column. That would be a 1.5–6σ effect in existing data.
- The author's kill threshold (γ < 0.005) is not reachable now. A **new limit at ≈ 0.015 mag (2σ)** is, and it would cut into the fiducial range. LSST will reach the threshold.

## References (all verified, 2026-10-09)

1. 0903.5303: Prescod-Weinstein & Smolin, *Disordered Locality as an Explanation for the Dark Energy*, 2009
2. gr-qc/0702044: Markopoulou & Smolin, *Disordered locality in loop quantum gravity states*, 2007
3. hep-ph/0111311: Csáki, Kaloper & Terning, *Dimming Supernovae without Cosmic Acceleration*, 2001
4. astro-ph/9811316: Aguirre, *Dust Versus Cosmic Acceleration*, 1998
5. 1004.2053: Avgoustidis et al., *Constraints on cosmic opacity and beyond the standard model physics from cosmological distance measurements*, 2010
6. 2007.16153: Martinelli et al., *Euclid: Forecast constraints on the cosmic distance duality relation…*, 2020
7. 0709.0003: Deffayet & Menou, *Probing Gravity with Spacetime Sirens*, 2007
8. 1801.08160: Pardo et al., *Limits on the number of spacetime dimensions from GW170817*, 2018
9. 1712.08108: Belgacem et al., *The gravitational-wave luminosity distance in modified gravity theories*, 2017
10. 1710.05834: LVC + Fermi + INTEGRAL, *Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A*, 2017
11. 2111.03604: LVK, *Constraints on the cosmic expansion history from GWTC-3*, 2021
12. 1306.0901: Schneider & Sluse, *Mass-sheet degeneracy, power-law models and external convergence…*, 2013
13. 2506.03023: TDCOSMO Collaboration, *TDCOSMO 2025: Cosmological constraints from strong lensing time delays*, 2025
14. 1912.08027: Millon et al., *TDCOSMO. I. An exploration of systematic uncertainties in the inference of H0…*, 2019
15. 2309.10101: Nierenberg et al., *JWST lensed quasar dark matter survey I*, 2023
16. 2405.01620: Keeley et al., *JWST Lensed quasar dark matter survey II*, 2024
17. 1908.06983: Gilman et al., *Warm dark matter chills out… with 8 quadruple-image strong gravitational lenses*, 2019
18. 1611.00014: Goobar et al., *iPTF16geu: A multiply imaged, gravitationally lensed type Ia supernova*, 2016
19. 2211.00656: Goobar et al., *Uncovering a population of gravitational lens galaxies with magnified standard candle SN Zwicky*, 2022
20. 2403.18902: Pascale et al., *SN H0pe: The First Measurement of H0 from a Multiply-Imaged Type Ia Supernova…*, 2024
21. 1307.2566: Smith et al., *The Effect of Weak Lensing on Distance Estimates from Supernovae*, 2013
22. 2406.05047: Shah et al. (DES), *Detection of weak lensing magnification of supernovae and constraints on dark matter haloes*, 2024
23. 2112.03863: Scolnic et al., *The Pantheon+ Analysis: The Full Dataset and Light-Curve Release*, 2021
24. 2202.04077: Brout et al., *The Pantheon+ Analysis: Cosmological Constraints*, 2022
25. 2503.14738: DESI, *DESI DR2 Results II: BAO and Cosmological Constraints*, 2025
26. 1807.06209: Planck, *Planck 2018 results. VI. Cosmological parameters*, 2018
27. 1509.02162: Applegate et al., *Relaxed galaxy clusters IV: calibrating hydrostatic masses with weak lensing*, 2015
28. 1402.2670: von der Linden et al., *Robust Weak-lensing Mass Calibration of Planck Galaxy Clusters*, 2014
29. 1602.03385: Pizzuti et al., *CLASH-VLT: Testing the Nature of Gravity with Galaxy Cluster Mass Profiles*, 2016
30. 1312.5083: Terukina et al., *Testing chameleon gravity with the Coma cluster*, 2013
31. 1504.03937: Wilcox et al., *The XMM Cluster Survey: Testing chameleon gravity using the profiles of clusters*, 2015
32. 1806.08300: Collett et al., *A precise extragalactic test of General Relativity*, 2018
33. 2002.05736: Millon et al., *COSMOGRAIL XIX: Time delays in 18 strongly lensed quasars…*, 2020
34. Watts & Strogatz, *Collective dynamics of 'small-world' networks*, Nature 393, 440 (1998); not on arXiv, cited from memory.
