# Vetting 1316 and 2184 (Abell S1063 `arcs` anti): measurement-noise injections (plan, 2026-10-16)

Written before the run. Parent: `../RESULT.md` "Next".

- **Hypothesis to falsify (ordinary):** 1316 and 2184 have true axes aligned with the lensing prediction, and the
  anti axis comes from photon and ICL noise in the shape measurement at m ≈ 27.5. The earlier null (0.524 anti
  expected) models intrinsic shape noise only, not measurement noise.
- **Other ordinary explanations, not tested here:** photo-z failure (unlensed member or foreground), fragment of
  nearby extended light (2184 next to the giant arc), blend.
- **Test (`scripts/arcs_vet_noise.py`):** for each of the 24 strong-shear background rows, inject an elliptical
  Gaussian with the row's catalogue flux, σ_major = `semimajor_px`, ellipticity = catalogue `ellipticity`, major
  axis along `pa_pred_z2` (aligned), at N = 60 empty positions 1.5–6″ from the row on the same GLIMPSE F200W i2d
  (S3 byte range); re-measure with the `arcs_bgsub_shapes.py` method (1″ median filter, σ_w 3 px, peak within 3 px).
  ASSUMPTIONs: Gaussian profile (no PSF convolution, catalogue sizes are post-PSF), N, the annulus, "empty" =
  pre-injection |filtered pixel| peak within 4 px below 3 robust σ.
- **Frozen rules (functions `gate_ok`, `noise_class`, `excess_p` with tests):**
  - gate (method validity): the same injections at m = 24 come back `aligned` (≤ 30°) in ≥ 95 % of valid
    injections, per row, for every row; otherwise the test is invalid and supports nothing;
  - per target (1316, 2184 separately): P_anti = anti fraction (≥ 60°) of valid injections at the row's own
    magnitude; `noise-explicable` if P_anti ≥ 0.05, else `not noise-explicable`;
  - sample: Poisson-binomial P(≥ n_obs) with the rows' P_anti, n_obs = rows `anti` under `bgsub_shapes.ecsv`;
    `excess explained by measurement noise` if ≥ 0.05, else `excess remains`.
- **What changes the plan:** both targets noise-explicable → the S1063 anti excess is closed as a measurement-noise
  effect and the next step is the same calibration on the other `arcs` fields; a target not noise-explicable → the
  remaining tests (MUSE z, deblend) run on it next.
