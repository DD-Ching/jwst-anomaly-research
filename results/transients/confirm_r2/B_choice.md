# B. One discriminating test: is the B1-PLCKG165 h002 fade achromatic? (fixed before new data, 08:47Z session clock)

## Why this object and this test
- Among the open real changes, B1-PLCKG165-box2-h002 (171.805899, +42.481802) is the only one with **two epochs in
  the same three bands while it was bright**: JWST 1176 o251 (MJD 60033.2, 2023-03-30) and 4446 o002 (MJD 60057.0,
  23.8 d later) in F090W, F150W, F200W, plus late non-detections (4744 2025, VENUS 2026) and a 2016 HST limit
  (`../batch_venus/boxes/review/lightcurves.ecsv`, `../batch_venus/followup/h002_hst2016.json`). MACS1423 h000 and
  Bullet h019 have one bright epoch each in JWST; their extra HST epochs are unchecked and would need new reads.
- **Exploratory, not blind:** these h002 fluxes were measured and looked at in #296 (e.g. F090W 271 → 173 nJy, F150W
  240 → 214). New here are the matched-control calibration of the 1176 → 4446 flux scale per band and the
  identifiability simulation; the decision rule below is fixed before either is computed.

## The two explanations compared
- **H1, intrinsic transient** (e.g. a supernova): the source's own spectrum changes, so the flux ratio between the two
  epochs differs between bands (colour evolution).
- **H2, gravitational magnification of a persistent compact source** by a compact deflector of any of these kinds:
  an ordinary point mass (star, compact object), an **Ellis wormhole (Kitamura–Nakajima–Asada n = 2, ε > 0)** or a
  **negative-mass point lens (n = 1, ε < 0; the W3 "umbra" fade of D-047)**. In geometric optics the magnification of
  a point or uniform-colour source does not depend on wavelength, so the ratio is the same in every band.
- What the test cannot do: it does not separate the three H2 deflectors from each other (that needs the light-curve
  shape: W3 spike–umbra–spike vs a single peak, ≥ 3 well-sampled epochs). A source with a colour gradient (limb
  darkening, an accretion disc) crossing a caustic can be chromatic at the few-per-cent level; the test only rejects
  chromaticity larger than its calibrated noise. H2 also needs the source to persist; the late non-detections
  (F150W 2025: −2.0 ± 5.3 nJy) imply a magnification ≳ 15 at the first epoch if H2 holds (recorded, not decisive).

## Measurement (fixed)
- Controls: `hit_lightcurve_controls.py pick` on the 1176 o251 F150W mosaic, `--half-arcsec 16 --flux-range 120 480
  --n 80 --seed 2`; fallback `--half-arcsec 30` if < 20; drop controls within 0.3″ of h002, de-duplicate at 0.5″
  (as `PLAN.md` A). Photometry: `hit_lightcurve.py` 0.15″, `--query-only h002 --group-arcsec 90 --only-program 01176
  --only-program 04446`.
- Per band b ∈ {F090W, F150W, F200W}: h002 ratio R_b = F(4446 o002)/F(1176 o251); control median ratio m_b and robust
  scatter s_b = 1.4826 × MAD of the control ratios (controls with S/N ≥ 5 in both epochs). Calibrated ratio
  ρ_b = R_b / m_b; σ_b² = (stat error of R_b)² + s_b² (conservative: s_b already contains control noise).
- **Rule:** χ² = Σ_b (ρ_b − ρ̄)² / σ_b² with ρ̄ the weighted mean, 2 dof. p < 0.01 ⇒ "the fade is chromatic beyond the
  calibrated noise: H2 (achromatic magnification of a uniform-colour source, including the Ellis and negative-mass
  cases) is disfavoured for this event"; p ≥ 0.01 ⇒ "achromatic within errors: not distinguished".
  Fewer than 20 usable controls in any band ⇒ not decided.

## Identifiability simulation (run after the controls, before computing ρ_b; `simulated`)
- Same epochs, bands and errors (h002's measured flux errors at both epochs plus s_b). 2,000 draws each:
  - H2: equal ratio in all bands, r ~ U(0.3, 1.5), t1 fluxes as measured → false-rejection rate at p < 0.01;
  - H1: sncosmo 2.x templates hsiao (Ia), nugent-sn2p, nugent-sn2l, nugent-sn1bc; z ~ U(0.3, 1.5); first epoch at
    rest-frame phase ~ U(−10, +40) d from peak; no dust; scaled to the measured t1 F150W flux → power = fraction
    rejected at p < 0.01.
- If the H1 power is < 0.5, a "not distinguished" outcome is uninformative and is reported that way.
- Stop for B: 11:30Z.

## Deviation (08:44Z session clock; after the control calibration, before any target ratio or simulation)
- The 16″ pick gave 26 positions (23 after the entry/duplicate filter) but only 10 / 15 / 15 controls usable in
  F090W / F150W / F200W (S/N ≥ 5 at both epochs; `g165_h002/calibration.json`), below the 20 the rule needs. The
  frozen result for that set is "not decided" and is kept. As a labelled deviation, the 30″ pick (seed 2, same flux
  range; `g165_h002/pos30.txt`) is run through the same filter, photometry, calibration and rule; its outcome is
  reported as the deviation's result. No threshold changes.
