# Exotic-lens signatures: what the literature predicts and what our data can test

Survey of 2026-10-08 (D-023). Every reference was opened on arXiv; titles are exact. This page only says what
would distinguish an exotic lens from ordinary lensing. Any detection would be a hypothesis to test, never a
conclusion (charter).

## Signatures

| Lens | Predicted signature | Testable here? |
|---|---|---|
| Ellis wormhole (the "1/r^n" lens family, n = 2) | **Demagnification**: sources beyond about 2/(n+1) Einstein radii appear fainter than unlensed. Images are still stretched tangentially. Light curves dip by ~4% just outside each Einstein-ring crossing. The image centroid shift is a few µas and traces an oval path. | Static demagnification against a cluster model: in principle, from one epoch. Light curves and µas shifts: no. |
| Negative mass / negative convergence | **Radially stretched images**; an umbra (eclipsed region); light curves with a central dip flanked by spikes; an early, not delayed, arrival time | Radial stretching plus demagnification around one dark centre that no positive-mass perturber can fit: yes, in principle |
| Alcubierre-type warp bubble | No published imaging or light-curve prediction for a distant observer, only views from or near the bubble and gravitational waves from bubble collapse | No; a model would have to come first (re-checked for D-047, see "Warp signatures") |

Key references:
- Abe, arXiv:1009.6084, "Gravitational Microlensing by the Ellis Wormhole".
- Kitamura, Nakajima & Asada, arXiv:1211.0379, "Demagnifying gravitational lenses toward hunting a clue of exotic
  matter and energy".
- Izumi et al., arXiv:1305.5037, "Gravitational lensing shear by an exotic lens object with negative convergence
  or negative mass".
- Cramer et al., arXiv:astro-ph/9409051, "Natural Wormholes as Gravitational Lenses".
- Safonova, Torres & Romero, arXiv:gr-qc/0105070, "Microlensing by natural wormholes: theory and simulations".
- Asada, arXiv:1711.01730, "Gravitational lensing by exotic objects" (review).
- Warp drives:
  - Clark, Hiscock & Larson, arXiv:gr-qc/9907019, "Null geodesics in the Alcubierre warp drive spacetime: the view
    from the bridge";
  - Clough, Dietrich & Khan, arXiv:2406.02466 (gravitational waves from warp-drive collapse).

## Where the signatures overlap with ordinary effects

The literature itself names these degeneracies, and each must be excluded first.
- **Radial images:** they also form near a cluster's radial critical curve, and voids give negative convergence
  (Izumi et al.).
- **Demagnified images:** also from saddle-point images and from microlensing by intracluster stars, e.g. Kelly
  et al., arXiv:1706.10279, "Extreme magnification of a star at redshift 1.5 by a galaxy-cluster lens".
- **Flux-ratio anomalies:** also from dark-matter substructure (Mao & Schneider, arXiv:astro-ph/9707187) or dust.
- **Two-epoch brightness changes:** not diagnostic of exotic lenses. They are degenerate with microlensing by
  stars, intrinsic variability with time delays, and supernovae.

## Existing observational limits

- Takahashi & Asada, arXiv:1303.1301 (ApJL 768, L16): no exotic lensing among about 50,000 SDSS quasars.
  - Negative-mass objects: n < 10⁻⁸ h³ Mpc⁻³ above 10¹⁵ M☉, and n < 10⁻⁴ h³ Mpc⁻³ above 10¹² M☉.
  - Ellis wormholes: n < 10⁻⁴ h³ Mpc⁻³ for throat radii 10–10⁴ pc.
- Anchordoqui et al., arXiv:astro-ph/9904399, "In search for natural wormholes": 631 BATSE gamma-ray-burst
  profiles.

**Consequence for this project.** The realistic outcome is a null result with better limits, which is still
publishable science. The search therefore tests lensing consistency in cluster fields against published models:
- First, ordinary checks: image positions, parities, flux ratios, and arc orientation against the predicted
  shear.
- Only then the exotic-specific patterns from the table, after ordinary explanations fail (TASKS "Now (M3)").

## Wormhole signatures (D-047)

Predictions for JWST imaging and time-domain data, from the verified papers in SOURCES.md "Exotic-lensing
predictions (D-047)". Every number below is a **model_prediction** of the cited lens model, or **simulated** with
`src/jwst_anomaly/exotic_sim.py` (tests in `tests/test_exotic_sim.py` reproduce the papers' statements). Physical
scales are **derived** for z_l = 0.4, z_s = 2 (Planck18) and v⊥ = 1000 km/s (**assumption**). The existence of
any such lens is a **hypothesis**.

**Model.** Weak-field deflection α = ε̄/bⁿ (Kitamura, Nakajima & Asada 2013). n = 1 is a point mass, n = 2 the Ellis
wormhole (ε̄ = πa²/4, a = throat radius; Abe 2010), and ε < 0 a negative mass. In Einstein-radius units the lens
equation is β = x − sign·sgn(x)/|x|ⁿ, and μ = 1/(λ_t λ_r) with λ_t = 1 − sign/|x|ⁿ⁺¹ and λ_r = 1 + sign·n/|x|ⁿ⁺¹
(Izumi et al. 2013). θ_E = (|ε̄| D_LS / (D_S D_Lⁿ))^(1/(n+1)) (Kitamura et al. 2014).
- Negative point mass: θ_E = 2.2″ (|M|/10¹² M☉)^½; resolved by NIRCam (θ_E ≳ 0.05″) for |M| ≳ 10⁹ M☉.
- Ellis wormhole: θ_E = 0.03″, 0.15″, 0.7″, 3.3″ for a = 0.1, 1, 10, 100 pc.
- Microlensing time scale: t_E ≈ 12 yr (|M|/M☉)^½ for a negative point mass. JWST epochs (months to years)
  therefore probe |M| ≲ 1 M☉. An Ellis lens needs a ≲ 10⁶ km for t_E ≲ 10 yr.

| # | Signature (label) | What JWST data can see | Amplitude (formula, parameters) | Time scale | Abundance (verified limits) | Ordinary mimics | Screen |
|---|---|---|---|---|---|---|---|
| W1 | Negative mass, static (model_prediction) | Two images on the source's side, both radially stretched, around a dark centre. No image of a source with β < 2θ_E. | x = (β ± √(β² − 4))/2; μ = x⁴/(x⁴ − 1); total μ = (u² − 2)/(u√(u² − 4)) (Safonova et al. 2002). At β = 2.05: μ = 1.69 and −0.69, radial/tangential stretch 4.6. At β = 3: μ = 1.02 and −0.02, stretch 1.3. | static | n < 10⁻⁸ h³ Mpc⁻³ for \|M\| > 10¹⁵ M☉; n < 10⁻⁴ h³ Mpc⁻³ for \|M\| > 10¹² M☉ (Takahashi & Asada 2013) | radial arcs at a cluster's radial critical curve; voids (Izumi et al. 2013); edge-on or irregular galaxies; diffraction spikes (D-034, D-043) | radial (D-031); counts (W5) |
| W2 | Ellis wormhole, static (model_prediction) | Two tangential images on opposite sides (like a point mass) with no luminous deflector. | Total μ < 1 for β > 1.11; minimum 0.958 at β = 1.67. Inner image at x = −0.618 (β = 2) carries 3.4 % of the flux (Abe 2010). | static | n < 10⁻⁴ h³ Mpc⁻³ for a = 10–10⁴ pc (Takahashi & Asada 2013) | dark subhalo or faint galaxy lens; chance pairs | dark-lens search; flux ratio (D-031). The 4 % deficit is not measurable on one object. |
| W3 | Negative mass, microlensing (model_prediction, simulated) | A compact source (a caustic-crossing star in a cluster arc, or an AGN disk) **vanishes**, between two caustic spikes. | Umbra β < 2: the lensed flux goes to 0 for 2 t_E √(4 − u₀²). Spike peak (uniform disk of radius ρ θ_E): ≈ 7.5 for ρ = 0.01, 3.4 for ρ = 0.1, 1.7–2.1 for ρ = 0.3. For n = 2 (n = 3) repulsive lenses the caustic is at β = 1.89 (1.76). | t_E ≈ 12 yr (\|M\|/M☉)^½ | negative-mass density ≲ O(10⁻³⁶) g cm⁻³ from GRB data (Torres, Romero & Anchordoqui 1998) | demagnified saddle-point images and microlensing by intracluster stars (Kelly et al.); SN fading; AGN variability; persistence (D-039); subtraction residuals | transient / dimming (D-027) |
| W4 | Ellis or n > 1, microlensing (model_prediction) | Shallow, time-symmetric "gutters" either side of the peak. | 4.2 % deep for n = 2 (Abe: about 4 %), still 4 % for ρ = 0.3. 14 % for n = 3 (the paper says ~10 %) and 59 % for n = 10 (Kitamura et al. 2013). A finite source makes them shallower (Tsukamoto & Gong 2018). | few t_E | none | any 4 % variability | Not practical: 5σ needs about 0.8 % photometry per epoch. Only n ≳ 3 is testable. |
| W5 | Image-plane counts (simulated) | A deficit of sources above a flux limit inside about θ_E of a dark centre, surrounded by radial arcs (ε < 0). | N_obs/N = N(>S/\|μ\|)/(\|μ\| N(>S)) (`count_ratio`). At x = 0.3, \|μ\| = 0.008. | static | as W1 | masks around bright stars and galaxies; deblending; cosmic variance | counts, around radial-screen centres only (without a lens position, stacking is impossible) |
| W6 | Centroid shifts, wave optics (model_prediction) | Nothing: µas shifts (Toki et al. 2011; Kitamura et al. 2014) and femtolensing of GRBs (Yoo et al. 2013) | — | — | n ≲ 10⁻⁹ AU⁻³ for a ~ 1 cm (Yoo et al. 2013) | — | none |

Notes.
- For n > 1 and ε < 0 the inner image at x < 1 has |λ_r| > λ_t, so its shape is tangential. Izumi et al.'s
  "everywhere radially elongated" compares the signed eigenvalues. The outer image is always radial.
- The isolated-lens formulas ignore cluster shear and macro-magnification (**assumption**). Inside a cluster arc, a
  negative-mass microlens sits on a macro-caustic network. Treat the W3 amplitudes as indicative, and the source-plane
  velocities and t_E as lower limits on how fast things change.
- Injection helpers: `inject_images` (image positions, μ, stretch and position angle for catalogue painting) and
  `inject_light_curve` (fluxes per epoch, with a blend fraction). Both return `simulated` tables.

## Warp signatures (D-047)

**Result: no published warp-drive model gives an imaging or photometric prediction usable with JWST data.** This
branch stops here until a paper supplies one. We do not invent one.
- **Alcubierre, constant velocity.** Outside the bubble the spacetime is exactly flat: no ADM mass and no
  gravitational waves (Clough, Dietrich & Khan 2024, citing Schuster, Santiago & Visser 2023). Only rays that cross
  the wall are deflected. A bubble of radius R at distance D subtends R/D, about 3 × 10⁻¹⁵ rad for R = 100 m at
  1 kpc (**derived**). That is far below any resolution or lensing cross-section.
- **Views from or near the bubble** (Clark, Hiscock & Larson 1999; Müller & Weiskopf 2011): aberration, horizons and
  redshifts seen by an observer **inside** the bubble. Not applicable to a distant telescope.
- **Collapse** (Clough et al. 2024): gravitational waves only, at f ~ 300 kHz and h ~ 10⁻²¹ at 1 Mpc for a 1 km bubble
  (**simulation**). The paper makes no electromagnetic prediction.
- **Positive-energy or "physical" warp drives** (Lentz 2021; Fell & Heisenberg 2021; Bobrick & Martire 2021): a warp
  drive is a shell of matter. If it has an ADM mass, it lenses as an ordinary compact mass of that mass
  (**assumption**). No distinct signature has been published. Fell & Heisenberg's example energies are about 10⁻⁴ M☉,
  which gives a θ_E of µas or less.
- **Particles swept up by a decelerating bubble** (McMonigal, Lewis & O'Byrne 2012) and **atmospheric glow** of a
  bubble faster than 0.1c in Earth's air (Fell & Loeb 2026): both effects are local to the destination or to Earth.
  Neither paper gives a flux, spectrum or rate for a distant source.
- The warp row of the "Signatures" table above stands. No JWST screen is calibrated for warp signatures.

## Screens and results (D-031)

`scripts/exotic_screens.py fluxratio` (demagnification) and `radial` (radial arcs around a dark centre). Every
hit goes to `/vet-candidate`.

| Field | fluxratio (compact images) | radial (convergence centres) | Verdict |
|---|---|---|---|
| WHL0137 / Sunrise (RELICS Lenstool maps, D-033) | no image list | 3 against a null mean of 2.2 (p95 5); max 3 lines, p 0.885 | null |
| SMACS 0723 (ICLv2) | 6 compact images consistent, 0 flags (49 resolved) | 12 against a null mean of 9.0 (p95 15); max 4 lines, p 0.945 | null |
| MACS1149 (HFF CATS, D-037) | 19 consistent, 0 flags (106 resolved; no photo-z) | with CANUCS photo-z: 12 against a null mean of 12.8 (p95 20); max 3 lines, p 1.0 (without: 32 against 30.8, p 0.91) | null |
| MACS0416 (HFF CATS, D-038, D-042) | 32 consistent, 3 flags explained: 45.1/45.2 blend, 38.1 marginal (58 resolved) | JWST frame, with CANUCS photo-z: 81 against a null mean of 84.9 (p95 105); max 7 lines, p 0.225 (without: 101 against 100.2, p 0.245) | null |
| MACS0717 (HFF CATS, D-041) | 26 consistent; system 65 flags explained (catalogue offset; μ model-dependent) | 11 against a null mean of 16.7; max 5 lines, p ≥ 0.70 (no photo-z) | null |
| Abell 370 (HFF CATS, D-043) | not run (image list gated) | with a Gaia spike veto: max 5 lines, p 0.945 (0.495 with CANUCS photo-z); two raw flags were spike chains | null |
| Abell S1063 (HFF CATS, D-043) | not run (image list gated) | max 4 lines, p 0.435; with DJA photo-z max 3, p 0.95 | null |
