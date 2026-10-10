# Neutrino Frontier — Round 1: Adversarial review (team N-B)

Review of `docs/neutrino_frontier/round1_invention.md` (team N-A, branch `claude/neutrino-frontier-r1`).
Literature checked 2026-10-10 (arXiv API, web search). Every reference below was opened or its arXiv ID
was resolved through the arXiv API; none comes from memory alone. Labels: **[O]** observed / published
measurement, **[D]** derived here, **[A]** assumption, **[P]** model_prediction. Nothing here is a detection.

## Summary

| ID | Verdict | Ladder | Key constraint | Prior art |
|----|---------|--------|----------------|-----------|
| F1 stochastic distance | KNOWN-REDUNDANT + UNTESTABLE NOW | B (hidden preferred frame) / C in region | TXS 2014–15 flare width ≈ 110 d ⇒ κ(E/100 TeV)^{2n} ≲ 2×10⁻¹⁰; no second resolved multi-event ν flare exists | lightcone fluctuations (2103.15313), stochastic LV with SN ν (1110.4848), path-length foam / image blurring (1411.7262) |
| F2 hidden partners | KNOWN-REDUNDANT (it is the pseudo-Dirac neutrino) | C, D already done by others | N-A's "testable window" 10⁻¹⁹–10⁻¹⁷ eV² is **already excluded** (2.1×10⁻²¹–2.0×10⁻¹⁶ eV², ≥ 90 % CL, public IC86 tracks) | hep-ph/0307151, 2211.16520, 2406.06476, 2105.12736, 2205.13291 |
| F3 second causal cone | KNOWN-REDUNDANT (ν–γ differential Shapiro delay / WEP, PPN Δγ) | B (absolute-Φ dependence ill-posed) / C | SN1987A: Shapiro delay equal to 0.2–0.5 % ⇒ \|ε\| ≲ few×10⁻³; N-A's LSS-based 10⁻⁴ is not valid (1907.12453) | PRL 60, 173 and 176 (1988); 1807.05201; 1807.06504 |
| F4 nonlocal ghosts | CONDITIONALLY VIABLE as phenomenology only; TXS-ghost test REJECTED (undecidable) | B (fragile; A2 underived) | short-delay alert ghosts already bounded by D-074 ICECAT–ICECAT wide nulls; R_g ≲ 0.04 for Δt_g ≤ 100 s [D] | no ghost-image paper found; nearest: ν "shortcuts" through extra dimension (hep-ph/0504096), GW lensing repeat searches (2304.08393) |

None reaches READY FOR QUANTITATIVE INVESTIGATION in the sense of a test with a plausible discovery.
One cheap bound-setting experiment (F4, re-using D-074) is specified at the end.

## Corrections to N-A's common inputs

| Input | N-A | Checked value | Source / note |
|-------|-----|---------------|---------------|
| TXS 2014–15 flare | ≈ 13 events over ≈ 110 d | 13 ± 5 events, Gaussian width 110 d, box 158 d; 5.8 bkg events in 1° over the box; ≈ 3.5σ **at the fixed a-priori position** | 1807.08794 (secondary summaries agree) |
| TXS distance | 1.4 Gpc | z = 0.3365: comoving ≈ 1.35 Gpc (use this for path length), luminosity ≈ 1.8 Gpc | standard ΛCDM [D] |
| SN1987A ν–optical offset | ν ~hours early | ≈ 3 h early; Galactic Shapiro delay toward LMC 1–6 months depending on halo model | PRL 60, 176 |
| SN1987A Shapiro test | — | ν and γ delays equal to 0.2 % (Longo) / ≥ 0.5 % (Krauss–Tremaine) | PRL 60, 173; 60, 176 |
| ICECAT-1 | ~350 tracks / 12 yr | 348 in v4, 340 after CR_VETO (D-074); 275 in the v1 paper (2011–2020) | 2304.01174, D-074 |
| GBM × ICECAT overlap | ≈ 3,500 | 2,910 GBM anchors with an ICECAT event within ±30 d | D-074 limits.ecsv |
| GWTC × ICECAT | ≈ 150 | 66 eligible ICECAT anchors in ICECAT–GW channel | D-074 |
| GBM × ICECAT spatial match | combined 5° radius | D-074 uses 3σ_comb with σ_GBM ≥ 3.7° (sys) ⇒ typical radius 11–50°; observed σ_comb of matched pairs 4–38° | D-074 same_dir_pairs |
| 10-yr track sample | ~1.1×10⁶ | 1.13×10⁶ (2008–2018 release) — **not in this repo** | 2101.09836 |
| Pseudo-Dirac SN1987A window | "compatible 10⁻²⁰–10⁻¹³" | SN1987A excludes [2.55, 3.01]×10⁻²⁰ eV² (Δχ² > 9), mild hint at 6.3×10⁻²⁰ (Δχ² ≈ 3) | 2105.12736 |

## F1 — Granular locality (stochastic distance)

**Prior art.** Stochastic propagation of neutrinos through fluctuating spacetime (lightcone fluctuations;
arrival-time spread and flavor decoherence) is worked out in Stuttard, arXiv:2103.15313. Stochastic
(non-systematic) Lorentz violation for SN neutrinos: arXiv:1110.4848. Random-walk path-length foam and the
resulting image degradation, with the random-walk model (α = 1/2) excluded for photons by Chandra/Fermi/IACT:
Perlman et al., arXiv:1411.7262. The D² delay from small-angle multiple scattering is the textbook
formula (also used for cosmic-ray deflection and FRB scattering). F1's only new element is the assumption
that photons are exempt (A3), which converts the photon exclusion of 1411.7262 into a free parameter.

**Consistency.**
- A1 says kicks are isotropic about the direction of motion, with a fixed variance per unit length. A
  stochastic medium with a rate per length has a rest frame; a boosted observer sees a different rate, so
  "Lorentz invariant on average" is not true. F1 is a Lorentz-violating neutrino-only medium. That is
  allowed but must then face neutrino LV bounds, not just timing.
- A3 (photons exempt) contradicts "spacetime is the substrate": a geometric fluctuation is universal (EP).
  A neutrino-only effect is a new interaction (a dark-sector scattering medium), not geometry.
- Energy conservation with direction kicks requires the medium to absorb momentum; fine for a medium,
  again breaks translation invariance of the vacuum.
- Path jitter δL ~ cΔt ≈ 3×10¹² km at the TXS edge gives oscillation phase jitter ~10⁵ rad at 100 TeV [D]:
  full decoherence, but astrophysical flavor is already averaged, so no observable flavor effect. N-A's
  "lab δφ ≈ 0" is right.

**Constraints (re-checked).** TXS: the 2014–15 multi-event flare has width 110 d ≈ 9.5×10⁶ s [O]; this
bounds the *spread* without any photon timing. N-A's τ_TXS ≲ 10⁷ s and κ(E/100 TeV)^{2n} ≲ 2×10⁻¹⁰ stand
(with D_C = 1.35 Gpc the bound tightens by ≈ 8 %) [D]. SN1987A bound for n < 0 as N-A wrote.

**Distinguishable prediction vs sensitivity.**
- Angle: ≲ 4″ vs 0.5–1° PSF. Undecidable (N-A correct).
- Time: needs the *width* of ≥ 2 multi-event neutrino flares at different distances. Single alerts measure
  no width; delays of single alerts behind photons (e.g. TDE AT2019dsg, ν ≈ 150 d after optical peak,
  arXiv:2005.05340) are dominated by astrophysical modelling. Only one multi-event flare (TXS) exists, and
  NGC 1068 is steady. The GBM × ICECAT loss-of-signal argument cannot be separated from the null result of
  IceCube GRB stacking (arXiv:1702.06868), which already finds no GRB neutrinos at any window tested.
- Background estimate in N-A (GWTC × ICECAT 2.6×10⁻³) uses a 5° radius; GW maps are 10²–10³ deg² and the
  repo has no GW positions (D-074), so the realistic background is 10–100× higher.

**Verdict.** KNOWN-REDUNDANT (mechanism) and UNTESTABLE NOW (distinct D² prediction). Ladder B/C.
Revisit if: a second resolved multi-event neutrino flare with known distance (IceCube-Gen2, KM3NeT).

## F2 — Twin-sheet flavor (hidden partner states)

**Prior art.** This is exactly the pseudo-Dirac neutrino: each mass state is a maximally mixed pair of
active and sterile Weyl states with tiny δm². Proposed for neutrino telescopes by Beacom et al.,
arXiv:hep-ph/0307151; constrained with NGC 1068 by Rink & Sen, arXiv:2211.16520; with public IC86 PSTracks
(single sources and stacking: NGC 1068, NGC 4151, PKS 1424+240, TXS 0506+056) by Dixit, Miranda & Razzaque,
arXiv:2406.06476; with SN1987A by Martinez-Soler, Perez-Gonzalez & Sen, arXiv:2105.12736; future Galactic SN
reach in arXiv:2205.13291. The "second sheet" wording adds no prediction beyond the pseudo-Dirac one.

**Consistency.** Fine (SM + right-handed neutrinos with tiny Majorana terms; unitary; Lorentz invariant).
N-A's survival formula P = 1 − sin²(1.27 δm² L/E) is the maximal-mixing pseudo-Dirac formula. Arithmetic
re-checked: TXS phase 5.5×10¹⁷ δm² at 100 TeV and E₁ ≈ 70 TeV for δm² = 2×10⁻¹⁸ eV² are right [D];
for z = 0.34 the redshift-weighted effective baseline is ≈ 15–20 % shorter than D_C [D, rough].

**Constraints (corrections).**
- **N-A's target window 10⁻¹⁹–10⁻¹⁷ eV² is already excluded**: 2406.06476 excludes
  δm² ∈ [2.1×10⁻²¹, 2.0×10⁻¹⁶] eV² (0.5 TeV–1 PeV; ≥ 90 % CL; [1.1×10⁻²¹, 3.0×10⁻¹⁶] with 0.1 TeV), using
  public IceCube data. N-A's "undecidable with public data" is therefore wrong as stated; the published
  exclusion relies on source-flux model assumptions (single power laws), which is its main caveat.
- SN1987A: excluded band [2.55, 3.01]×10⁻²⁰ eV²; weak preference at 6.3×10⁻²⁰ eV² (Δχ² ≈ 3) [O, 2105.12736].
- Remaining open: δm² ≲ 10⁻²¹ eV² (diffuse/cosmological baselines) and 2×10⁻¹⁶ ≲ δm² ≲ ~10⁻¹¹ eV²
  (upper end from solar data, as summarized in 2406.06476).
- ICECAT-1 cannot add: one or two alerts per source and no flavor. Diffuse index hardening Δγ ≈ 0.2 is
  degenerate with astrophysical spectral freedom (N-A correct).

**Verdict.** KNOWN-REDUNDANT. Ladder C (and D, done by others). No compute here.
Revisit if: Galactic SN (2205.13291 reach ~10⁻²⁰ eV²) or KM3NeT/Gen2 per-source spectra.

## F3 — Second causal cone (neutrino-specific Shapiro delay)

**Prior art.** An energy-independent, flavor-universal difference between neutrino and photon Shapiro
delays is the standard ν–γ weak-equivalence-principle test, parametrized by the PPN difference
Δγ = γ_ν − γ_γ; with Δt = (1+γ)∫U dl/c³, N-A's ε ≈ Δγ/2 [D]. SN1987A: Longo, PRL 60, 173 (1988);
Krauss & Tremaine, PRL 60, 176 (1988). TXS / IC-170922A: Boran, Desai & Kahya, arXiv:1807.05201
(|Δγ| < 5.5×10⁻² for a 175 d window and ≈ 6300 d delay); arXiv:1807.06504 (~10⁻⁶ using the
Laniakea potential). Direction dependence ε T_S(n̂) is built into every such test (each computes T_S for
its own line of sight); it is not a new prediction.

**Consistency.**
- g̃ depends on the absolute Newtonian potential Φ and on the matter frame u_μ. Φ is defined only up to a
  constant locally and diverges on cosmological scales; Minazzoli et al. (arXiv:1907.12453, 2203.11215)
  show cosmological Shapiro-delay bounds have no conservative lower bound without extra assumptions. So
  N-A's extragalactic rms T_S ≈ 10¹¹ s and the resulting |ε| ≲ 10⁻⁴ are **not valid**; only the Milky Way
  term is usable.
- u_μ selects a frame wherever Φ ≠ 0: local Lorentz violation of size ε Φ/c² ~ ε×10⁻⁶ at Earth. SN1987A
  speed bound |v−c|/c ≲ 10⁴ s / 5.3×10¹² s ≈ 2×10⁻⁹ ⇒ |ε| ≲ 2×10⁻³ [D], the same scale as the Shapiro bound.
- A3 (photons and GW on g) is consistent with GW170817 (|c_GW − c|/c ≲ 10⁻¹⁵; arXiv:1710.05834).

**Constraints (corrected).** SN1987A (MeV ν̄e): delays equal to 0.2–0.5 % of a 1–6 month delay ⇒
|ε| ≲ 10⁻³–5×10⁻³ (model range) [D from O]. TXS (TeV): |ε| ≲ 3×10⁻² with 1807.05201's potential; the
10⁻⁶ of 1807.06504 inherits the 1907.12453 objection. With A2 (energy-independent) SN1987A dominates.

**Distinguishable prediction vs sensitivity (corrected numbers).**
- Lags ε T_S(MW) ≲ 2×10⁻³ × 10⁷ s ≈ 2×10⁴ s: inside D-074's GBM–ICECAT "same" bins (0–1 d).
- **Background was underestimated ~20×.** N-A: 0.12 chance pairs over ±10⁴ s with a 5° radius. D-074's
  scrambled null with realistic 3σ_comb radii gives 0.84 pairs at 100 s–1 h and 20.2 pairs at 1 h–1 d
  (2.4×10⁻⁴ s⁻¹) ⇒ ≈ 2.4 chance pairs for |Δt| ≤ 10⁴ s [D]. "Two pairs on one ε line ⇒ p ≈ 10⁻²" fails:
  the Galactic T_S model is uncertain by a factor 2–6 (1–6 months for one line of sight), so an "ε line"
  constrains lag ratios only to that factor.
- Data already in hand (D-074): GBM–ICECAT same-direction, 0–10 s: 0 (null 0.0); 10–100 s: 0 (0.02);
  100 s–1 h: 2 (0.84); 1 h–1 d: 19 (20.2). The two 100 s–1 h pairs (bn201120402–IC201120A, +350 s,
  33° sep, σ_comb 16°; bn221223997–IC221224A, +3543 s, 22°, σ_comb 9°) both have ν after γ, but their
  lag ratio ≈ 10 exceeds the T_S model range for one ε; signed-lag test p = 0.09 per cell, global 0.68.
- Signal normalization is unknown and IceCube GRB stacking (1702.06868) is null, so a null scan yields no
  bound on ε, only a joint (ε, N_assoc) statement.

**Verdict.** KNOWN-REDUNDANT; the repo-level scan is effectively done (D-074, null). Ladder B/C.
Revisit if: a Galactic transient (SN) or a confirmed GRB/GW neutrino; then the Galactic-only ε T_S fit is
well posed.

## F4 — Sparse nonlocal links (neutrino ghost images)

**Prior art.** No paper found proposing neutrino-only ghost images at unrelated sky positions. Nearest:
active–sterile "shortcuts" through an extra dimension, where only gauge-singlet states take the short path
(Päs, Pakvasa & Weiler, arXiv:hep-ph/0504096; this is the only derived mechanism for A2-like selectivity,
and it acts on the sterile admixture, not on active neutrinos); neutrino oscillations near a
wormhole-like lens (arXiv:2412.02144); GW repeated-image searches (LVK O3, arXiv:2304.08393) use sky-overlap
statistics and would *not* flag a graviton ghost at an unrelated position — a genuinely untested corner,
but GWTC in this repo has no positions or PE samples (D-074), so not testable here.

**Consistency.**
- A2 (only neutrinos transmit) has no derivation. A wormhole is geometry; by the EP every particle
  follows it. Active neutrinos carry weak isospin and hypercharge, so a "neutral-only" channel would act
  on sterile components (hep-ph/0504096), making R_g ∝ sterile admixture, which oscillation data bound.
- A1: mouths at rest in the CMB frame still sit in different potentials and move with peculiar velocities
  ~300–600 km/s; accumulated clock offsets generically produce closed timelike curves for short links.
  N-A's own ladder B (fragile) stands.
- Flux bookkeeping (R_g = P_h/(1 − P_h)) is unitary; the isotropic ghost component is compatible with
  diffuse isotropy.

**Distinguishable predictions vs sensitivity (corrected).**
- *TXS 2014–15 ghost* — REJECTED as a test. The real TXS flare (13 ± 5 events) reached only ≈ 3.5σ at a
  fixed, a-priori position, with full energy weighting [O, 1807.08794]. A ghost at an unknown position
  pays all-sky trials: IceCube's untriggered 10-yr all-sky flare search found background hot spots at
  pre-trial p = 9.2×10⁻⁶ (north) and 3.5×10⁻⁷ (south), post-trial 0.69 and 0.06 [O, 2107.12134]. Even
  R_g = 1 (a full TXS copy) would be weaker than the background hot spots. N-A's "energy weighting cuts
  background 10× ⇒ R_g ≳ 0.45 decidable" double-counts: energy weighting is already in the 3.5σ. Since
  P_h ≤ 0.5 ⇒ R_g ≤ 1, the test is undecidable for every allowed R_g. It also needs the 10-yr track
  sample, which this repo does not hold.
- *Alert ghosts at short delay* — **already bounded by D-074** [D]. A ghost of an astrophysical alert-
  quality neutrino is another alert-quality track at a random direction with lag Δt_g: exactly the
  ICECAT–ICECAT *wide* channel. D-074 95 % limits on dependent wide partners per anchor (340 anchors):
  8.8×10⁻³ (Δt ≤ 100 s), invalid (100 s–1 h), 0.15 (1 h–1 d), 0.039 (1–7 d). Dividing by mean
  astrophysical signalness (~0.4, Gold/Bronze mix) and by the alert acceptance of a random-direction
  ghost relative to the parent (~0.5–0.7; alert effective area is declination-dependent) gives
  **R_g ≲ 0.03–0.04 for Δt_g ≤ 100 s**, ≲ 0.15–0.2 for 1–7 d, ≲ 0.5–0.7 for 1 h–1 d [D, rough; the
  correction factors are assumptions until computed per event].
- *Transient ghosts of GRBs*: undecidable (direct GRB–ν signal is null, 1702.06868); N-A correct.

**Verdict.** CONDITIONALLY VIABLE only as a phenomenological parameter space (P_h, r_M); the idea itself
stays at ladder B (A2 underived, A1 fragile). A null constrains (R_g, Δt_g), not the framework.

## Recommended first experiment (< 30 min, data in repo)

**F4 alert-ghost limit: R_g(Δt_g) from the ICECAT-1 wide channel.** It re-uses D-074 code
(`src/jwst_anomaly/event_network.py`, `scripts/e1_events.py`) and the pinned ICECAT-1 v4 table. It is the
only F-prediction this repo can bound that is not already covered by published work. Expect a bound, not
a discovery; label it a parameter limit on an ad-hoc model.

- **Statistic.** Number of ICECAT–ICECAT pairs with separation > 3σ_comb (and > 0.1°) and lag (ghost
  after parent, signed) in log bins 0–10 s, 10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d, 7–30 d, 30–180 d
  (extend D-074's bins to cover r_M/c up to the TXS-flare scale). Weight each anchor by its ICECAT-1
  signalness (column name to verify in the Dataverse file).
- **Null.** D-074 `jit` scrambles (Dec and hour angle kept; times permuted within year plus ±3 d jitter),
  10⁴ scrambles; for the 7–180 d bins use a year-block time permutation (jitter ±3 d is too short).
  Pooled-rank trials correction as in D-074.
- **Injection.** For each anchor, with probability R_g × signalness, add a synthetic track at Δt_g
  (log-uniform within the bin), RA uniform, Dec drawn from the ICECAT-1 Dec distribution (alert
  acceptance), error radius resampled from the catalogue. Find the R_g giving 50 % / 90 % detection at
  family-wise 3σ; the 95 % UL on R_g per bin is the limit.
- **Known case to reproduce first.** D-074's ICECAT–ICECAT counts (wide 1 h–1 d: 36 vs 25.4 ± 4.4;
  1–7 d: 138 vs 149.4) and `limits.ecsv` (8.8×10⁻³ at ≤ 100 s) must be reproduced exactly before the new
  bins are added. Sanity: the TXS 2014–15 events are not alerts, so no TXS ghost is expected in ICECAT.
- **Rejection rule.** A family-wise > 3σ excess in any signed wide bin that survives `jitday` and the
  year-block null and is not one alert multiplet or a reconstruction-correlated pair (same run / same
  day split events) would send the pairs to `/vet-candidate`. Otherwise the result is the R_g(Δt_g) UL
  table; the parameter region R_g above it is rejected for links with Δt_g in that bin.
- **Cost.** One run of the existing script with new bins and an R_g injection loop: 340 events, 10⁴
  scrambles, minutes on a laptop.

Not recommended: F3 ε-line rescan (D-074 already null, background ~2.4 pairs, T_S model factor 2–6,
no bound without signal normalization); F1 and F2 (no data in repo can add to published work).

## SOURCES.md lines (accessed 2026-10-10)

```
## Neutrino Frontier round-1 review (accessed 2026-10-10; docs/neutrino_frontier/round1_review.md)
- Stuttard 2021, "Neutrino signals of lightcone fluctuations resulting from fluctuating space-time", arXiv:2103.15313
- "Probing Lorentz Violation in Neutrino Propagation from a Core-Collapse Supernova", arXiv:1110.4848
- Perlman et al. 2015, "New Constraints on Quantum Gravity from X-ray and Gamma-Ray Observations", ApJ 805, 10, arXiv:1411.7262
- Beacom et al. 2004, "Pseudo-Dirac Neutrinos, a Challenge for Neutrino Telescopes", arXiv:hep-ph/0307151
- Rink & Sen 2022, "Constraints on pseudo-Dirac neutrinos using high-energy neutrinos from NGC 1068", arXiv:2211.16520
- Dixit, Miranda & Razzaque 2024, "Searching for Pseudo-Dirac neutrinos from Astrophysical sources in IceCube data", arXiv:2406.06476
- Martinez-Soler, Perez-Gonzalez & Sen 2022, "SN1987A still shining: A Quest for Pseudo-Dirac Neutrinos", PRD 105, 095019, arXiv:2105.12736
- Sen 2022, "Constraining pseudo-Dirac neutrinos from a galactic core-collapse supernova", arXiv:2205.13291
- Longo 1988, PRL 60, 173, doi:10.1103/PhysRevLett.60.173 (SN1987A nu-gamma Shapiro delay)
- Krauss & Tremaine 1988, PRL 60, 176, doi:10.1103/PhysRevLett.60.176 (SN1987A WEP test)
- Boran, Desai & Kahya 2019, "Constraints on differential Shapiro delay between neutrinos and photons from IceCube-170922A", EPJC 79, 185, arXiv:1807.05201
- "Multimessenger Tests of Einstein's Weak Equivalence Principle and Lorentz Invariance with a High-energy Neutrino from a Flaring Blazar", arXiv:1807.06504
- Minazzoli, Johnson-McDaniel & Sakellariadou 2019, "Shortcomings of Shapiro delay-based tests of the equivalence principle on cosmological scales", arXiv:1907.12453; Moriond summary arXiv:2203.11215
- LVK 2017, "Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A", arXiv:1710.05834
- IceCube 2018, "Neutrino emission from the direction of the blazar TXS 0506+056 prior to the IceCube-170922A alert", Science 361, 147, arXiv:1807.08794
- IceCube et al. 2018, "Multi-messenger observations of a flaring blazar coincident with high-energy neutrino IceCube-170922A", arXiv:1807.08816
- IceCube 2017, "Extending the search for muon neutrinos coincident with gamma-ray bursts in IceCube data", ApJ 843, 112, arXiv:1702.06868
- IceCube 2021, "Every Flare, Everywhere: An All-Sky Untriggered Search for Astrophysical Neutrino Transients Using IceCube Data", ICRC2021, arXiv:2107.12134
- IceCube 2021, "IceCube Data for Neutrino Point-Source Searches Years 2008-2018", arXiv:2101.09836
- IceCube 2023, "IceCat-1: the IceCube Event Catalog of Alert Tracks", arXiv:2304.01174
- Stein et al. 2021, "A tidal disruption event coincident with a high-energy neutrino", arXiv:2005.05340
- Päs, Pakvasa & Weiler 2005, "Sterile-active neutrino oscillations and shortcuts in the extra dimension", PRD 72, 095017, arXiv:hep-ph/0504096
- 2024, "The neutrino flavor oscillations in the static and spherically symmetric black-hole-like wormholes", arXiv:2412.02144
- LVK 2023, "Search for gravitational-lensing signatures in the full third observing run of the LIGO-Virgo network", arXiv:2304.08393
```
