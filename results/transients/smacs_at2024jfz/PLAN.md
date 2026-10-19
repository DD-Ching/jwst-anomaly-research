# Plan: lens-model counter-images and delays for AT 2024jfz (SMACS 0723; frozen 2026-10-19 before any prediction)

**Why.** Open item (TASKS "SMACS 0723 arc `arc_diff`"; #237 / `results/transients/smacs_arcdiff/RESULT.md` h000):
AT 2024jfz is a point source seen only in the 2024-05-11 GO-4043 mosaics (F090W 2579 nJy, F115W 2688, F444W 258;
`../smacs_arcdiff/o4043_o12774/vet/lightcurves.ecsv`), 1.81″ from catalogued image 6.1 of Mahler+22 ICLv2. If its host
is a background source, the cluster model predicts counter-images that repeat the event at model time delays.

**Hypotheses (labels).**
- H_bg (hypothesis): the transient sits in a lensed background source; then its counter-images are a lens-model
  prediction that existing or future epochs can test.
- Ordinary alternatives: a cluster-member or foreground transient (z ≤ 0.388, one image only); a background source
  at another redshift than system 6 (z 1.71 is Lenstool's fitted value for system 6, `z_m_limit`, a model value).
  The redshift of the host is not known here (no spectrum checked), so the prediction is made on a z_s grid.
- **What changes the plan:** if any counter-image is *testable* in a public JWST epoch (rule below), the next step is
  forced photometry at its predicted position in that epoch (a missing counter-image would argue against H_bg at that
  z_s, or against the model there); if none is testable, the result is the list of future windows.

**Method (`scripts/smacs_at2024jfz.py`, unit-tested).**
- Model: `smacs0723-iclv2` best model (Lenstool `best.par`, analytic dPIE; repo `lensmodel.LensModel`). No MCMC chain
  is pinned for this model, so delays carry no model spread (stated limitation).
- Source point β = θ − α(θ) at the TNS position of AT 2024jfz (110.83416, −73.451752; recorded in
  `../smacs_arcdiff/o4043_o12774/h000_tns.json`); images with `lensmodel.find_images` on a 0.1″ grid of ±60″.
- z_s grid (ASSUMPTION): 0.6, 1.0, 1.425 (system 5 / Atalanta's arc), 1.712766 (system 6 model z), 2.5, 4.0.
- Fermat potential φ(θ) = ½|θ − β|² − ψ(θ); ψ differences from the line integral of α along a straight path
  (Simpson, 2000 steps; α = ∇ψ). Delay Δt_k = (1 + z_l) D_l D_s / (c D_ls) [φ(θ_k) − φ(θ_0)], cosmology of the model
  file. θ_0 = the image nearest the TNS position. Positive Δt_k: image k trails the 2024 event.
- **Validation before use:** (a) unit test on a singular isothermal sphere (ψ = θ_E |θ|, analytic delays);
  (b) path independence on the real model: straight path vs an L-shaped path, |difference| ≤ 1 % of |Δψ| (else stop).

**Frozen testability rule (function `testable`, unit-tested before the run).** For each z_s and counter-image k
(k ≠ 0) and each public JWST NIRCam epoch e covering the field (MJD 59737.0, 60441.8, 61196.9, 61227.7, from the
lightcurve file above): testable if
`−2 ≤ t_e − (t_0 + Δt_k) ≤ 30` observed days (ASSUMPTION: the event, seen on one day only, stays near its 2024 flux
for ≤ 30 d) **and** the predicted flux 2579 nJy × |μ_k / μ_0| ≥ 10 × 17 nJy (the largest forced-photometry error of
any covering epoch's blue band in that file, rounded; ASSUMPTION).
Per z_s: outcome "testable" if any (k, e) pair passes, else "not testable"; for each counter-image with
|μ_k / μ_0| ≥ 0.066 the next window [t_0 + Δt_k − 30 d, t_0 + Δt_k + 2 d] is listed. Exploratory; no candidate can
come from this step.

**Data:** model files from the pinned Mahler+22 commit (sha256 in `lensmodel.SMACS0723_MAHLER22_ICLV2`); no images.
