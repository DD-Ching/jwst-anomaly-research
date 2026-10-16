# Vetting 1316 and 2184, step 2: deblend 1316, arc-fragment test for 2184 (plan, 2026-10-16)

Written before any measurement in this step. Parent: `RESULT.md` "Next". Script: `scripts/arcs_vet_deblend.py`
(rule functions unit-tested). Data: GLIMPSE `jw03293-o001_t001` level-3 i2d, SW bands F090W–F200W (LW bands are too
coarse for a 0.16″ pair and are used only for 2184's colours), 8″ S3 byte-range cutouts, 1″ median filter subtracted
as in `arcs_bgsub_shapes.py`. Labels: pixels observed, fits and colours derived, injections simulations.

## 1316: one elongated source or two knots?
- **Ordinary explanation tested:** the anti axis is the separation vector of two blended knots (two sources, or two
  clumps of one source), not a single lensed shape.
- **Models (F200W, pixels within 0.4″ of the catalogue position, ERR as σ):** M1 one elliptical Gaussian
  (x, y, flux, σ, q, θ); M2 two round Gaussians (x, y, flux, σ each). ΔBIC = BIC(M1) − BIC(M2).
  ASSUMPTIONs: Gaussian profiles that include the PSF; the fit radius; i2d pixel correlations ignored (χ² too small,
  which is why the threshold below is large and the rule is calibrated by injection).
- **Gate (method validity, frozen):** 20 injections each of (a) one elliptical Gaussian with 1316's catalogue flux,
  semimajor and ellipticity and (b) two round knots at 0.16″ with equal halves of that flux, at empty positions
  1.5–3.5″ from 1316 on the same image; the rule must call (a) not `two-knot` and (b) `two-knot` in ≥ 80 % each;
  otherwise the test is invalid and supports nothing.
- **Rule `deblend_class` (frozen):** `two-knot` if ΔBIC ≥ 10, both knots S/N ≥ 3 and separation ≥ 0.10″;
  `single` if ΔBIC ≤ −10; else `undecided`.
- **If `two-knot`:** (i) separation PA vs the z_s = 2 prediction (anti ≥ 60°, as `arcs.json`); (ii) knot SEDs: fixed
  positions and σ from F200W, fluxes fitted in F090W, F115W, F150W, F200W; `different SED` if the flux ratio f1/f2 is
  inconsistent across bands at χ² p < 0.003, else `consistent SED`.
- **Outcome classes for 1316:** `two-knot` + separation anti + `different SED` → **ordinary: blend of two sources**
  (the `arcs` test assumes one source), closed; `two-knot` + `consistent SED` → **clump pair**: the axis is a clump
  separation, stays open with that note; `single` or `undecided` → blend not supported, stays open.

## 2184: a fragment of the giant arc (z 1.2282)?
- **Ordinary explanation tested:** 2184 is a piece of the nearby giant tangential arc, so its shape is not an
  independent background galaxy's.
- **Measurement:** 0.2″-radius aperture fluxes in all nine bands (LW bands on their own pixels, no PSF matching;
  ASSUMPTION: colours of 0.2″ apertures are PSF-limited only at ≲ 0.1 mag in LW) at 2184 and at two control points on
  the arc ridge (the brightest 1″-filtered F200W pixels of the arc within 3″ of 2184, ≥ 0.5″ apart; found
  automatically and shown on the figure).
- **Control (frozen):** the two arc points must be `SED-consistent` with each other; otherwise the test is invalid.
- **Rule `sed_match` (frozen):** normalise each SED to its error-weighted mean; χ² of the difference over bands with
  S/N ≥ 2 in both; `SED-consistent` if p ≥ 0.003, else `SED-different`. 2184 `SED-different` from both arc points →
  the fragment explanation is excluded; else it stays possible.
- Recorded alongside, not a rule: the DJA photo-z interval (z16–z84 1.80–1.85) against the arc's 1.2282.

## What changes the plan
1316 closed as a blend and 2184 not a fragment → the S1063 sample has two survivors fewer or one, and the next step
is 2184's redshift (MUSE) or a second photo-z; both ordinary explanations excluded or undecided → the targets stay
open and the next step is a second photo-z code and the ESO MUSE cube.
