# Plan: star-ellipticity (PSF) radial pattern test for the Euclid Q1 radial-shear null (D-067)

Frozen 2026-10-20 17:22Z, before any star rows of this run were read.

**Question.** D-067 "Revisit if": do per-tile star-ellipticity (PSF anisotropy) gradients show radial patterns on
1–12′ scales? Galaxy shapes in `w5_euclid_shear.py` are deconvolved for PSF size only, so PSF anisotropy leaks into
them with factor β = T_psf / T_int (≈ tr_obs/tr_int − 1; median ≈ 0.6 from D-067's 1.6; ≤ 2.3 at the 1.2 × PSF
resolved cut). A radial PSF pattern could fake a flag (radial) or hide one (tangential).

**Hypothesis under test (ordinary):** PSF anisotropy contributes a coherent radial/tangential term comparable to
the injected signal.

**Data.** Observed: Euclid Q1 MER star rows (same selection as D-067: `point_like_prob > 0.9`, VIS 19–22) of all
344 Deep-Field tiles, IRSA TAP by `tileid` (cached; not committed).

**Measurement (derived).** Per Deep Field, on the D-067 grid for θ_E = 1′, 2′, 4′ (annulus 1.5–3 θ_E, step θ_E):
the star-weighted annulus mean radial component ē_r = −⟨Re(e · e^{−2iφ})⟩ over stars, and its error from the star
rms / √n (centres with n_star ≥ 5). Signal reference (model_prediction): the injected radial shear averaged over
the annulus, R · 2 ln 2 / 6.75 = 0.103 for R = 0.5 (area-weighted ⟨(θ_E/θ)²⟩, scale-free).

**Statistic.** ρ = β · max over centres of |ē_r| / 0.103 (field and θ_E), with β = 1 fiducial and 2.3 bound
(ASSUMPTIONs). Also χ² of ē_r / σ over centres against zero (independent-centre approximation, labelled).

**Frozen decision rule.**
- ρ(β = 2.3) < 0.2 at every field and θ_E: "PSF gradients cannot fake or hide a flag"; D-067 stands, its PSF
  revisit condition is closed.
- 0.2 ≤ ρ(β = 2.3) and ρ(β = 1) < 1: "PSF term not negligible at the edge of the resolved cut"; D-067 stands,
  correct anisotropy before Euclid DR1.
- ρ(β = 1) ≥ 1 anywhere: "PSF could fake or hide a flag"; re-screen with an anisotropy correction (new PR).

Noise caveat frozen with the rule: ē_r at n_star ≈ 5 is noise-dominated; max |ē_r| is then an upper bound on the
PSF pattern, so rule 1 is conservative, rules 2–3 need the χ² to exceed its dof before they are read as a pattern.
