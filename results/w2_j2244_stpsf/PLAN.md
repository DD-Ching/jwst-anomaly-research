# W2 J2244−0550 follow-up: minimal STScI official PSF route, capability check (WFC3/IR F160W only) — frozen plan

Claim and scope committed 2026-10-23 15:55 UTC (37f60ce). This test plan is committed before any fit is run.

- **Question:** is there a directly reusable official STScI PSF route that improves the point-source modelling that
  failed in #474 (`results/w2_j2244_hst/RESULT.md`)?
- **Type:** development feasibility test on controls already looked at in #474. Not an independent validation.
  J2244−0550 is not measured, and no exposure of its visit is read.
- **Hard limits:** wall clock 2 h (to 17:55 UTC); 2 cores; total download 1 GB including calibration files and new
  dependencies; WFC3/IR F160W only; no ACS; no multi-field star library; at most 3 small cases.

## 0. Download control (in force before any image request)
- `scripts/net_ledger.py`: a persistent, shared JSONL ledger (`net_ledger.jsonl`, cap 10⁹ B).
  - Every request reserves its expected size first; the bytes actually received (body, response header and an 8 kB
    TLS allowance) are recorded.
  - Retries and repeats are separate reservations. A reservation that is never closed keeps counting, so a restart
    does not reset the total.
  - Transfers are cut off before the cap; downloads made by other tools enter as measured `external` records.
- **Checked before any image request:**
  - unit tests on a local server (`tests/test_net_ledger.py`);
  - one 2880 B live request, counted 11651 B;
  - a live stop test on a 40 kB ledger (`net_ledger_stoptest.jsonl`): refused before the request, and cut mid-stream.
- **Dependencies:** installed with `uv --offline` from the local cache; 0 B downloaded.

## 1. Official entry and the route chosen (read 2026-10-23, ≤ 30 min)
- **Entry:** https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/point_spread_function/hst_point_spread_function.html
  - Pinned as the rendered page fetched 2026-10-23 (sha256 in `sources.json`), "Updated: 20 May 2026" (version
    history: photutils 3.0.0, hst1pass 2026.03.24_v2b).
  - The GitHub API commit lookup was refused (HTTP 403), so no commit hash is recorded.
- **Routes the notebook offers, and why one is chosen:**
  - **2.2.1 Default empirical models (chosen):** the WFC3-team STDPSF grid, read with photutils
    `GriddedPSFModel.read(format="stdpsf")` and fitted with `PSFPhotometry` on single exposures (FLT/FLC).
  - **2.2.2 / 2.2.3 Focus or spatial perturbation (hst1pass PERT):** needs at least 10 bright isolated stars per
    exposure. The GO 17455 F160W library holds 91 stars over 144 fields (#474), so the requirement is not met.
  - **2.3 Stacks (own or MAST cutouts):** build a PSF by stacking stars. That is a different method from the one
    tested here, and the notebook reports that MAST stacks fit only "marginally".
  - **3 Drizzle models:** stack stars in the drizzled image. This is what #474 did with photutils `EPSFBuilder`.
  - **Drizzling the empirical models into FLCs:** "under development and not yet officially supported".
- **What the chosen route supports** (from the notebook and from the file itself, `STDPSF_WFC3IR_F160W.fits`):
  - **Frame:** the exposure frame only. The notebook says the models "do not accurately represent the PSF in drizzled
    image mosaics". So the fits here run on the 4 FLT exposures of each visit, not on the `_drz`.
  - **Detector position:** a 3 × 3 grid at x, y = 0, 507, 1014 of the 1014-pixel science frame. Positions are
    interpolated.
  - **Sub-pixel position:** 4× oversampled, 101 × 101 sub-pixels (25 detector pixels across). The sub-pixel phase
    comes from the fitted position.
  - **Subarray:** IRSUB512 FLTs have LTV1 = LTV2 = −251, so full-frame x = x_sub − LTV1.
    - If the physical frame were the 1024-pixel raw frame instead, the position would be off by 5 px. That changes
      the model by at most 4.3 × 10⁻⁴ of the peak at the detector centre (derived), which is negligible.
  - **Normalisation and units:**
    - the sum over the 25 × 25 px model is 1.027–1.037 across the 9 grid PSFs; the sum within r ≤ 6 px is 1.008,
      and within r ≤ 3 px 0.943 (centre PSF, derived);
    - so the fitted "flux" is about the flux within r ≈ 6 px (0.8″), not the total;
    - FLT units are e⁻/s.
  - **Halo:** the models are "tapered to zero flux at large radii … at the expense of PSF wing accuracy". The centre
    PSF is zero beyond r ≈ 12.5 px (1.6″). The halo that #474 listed as a possible reason is therefore not modelled
    by this route, and the outer annuli below test how much that matters.
  - **Extra data or code:** none beyond the 371 kB model file and the FLTs (`sources.json`). No PSF solver is
    written: the model is used as is, with no PSF iteration and no fitting to the controls.

## 2. Data (observed; `sources.json` holds URLs, sizes and sha256)
- 12 FLT exposures of GO 17455 (3 visits × 4, WFC3/IR F160W, IRSUB512, 105.518 s, ELECTRONS/S):
  - if8a06030: if8a06cuq, cvq, cwq, cxq;
  - if8a4j030: if8a4jf4q, f5q, f6q, f7q;
  - if8a38030: if8a38guq, gvq, gwq, gxq.
- `STDPSF_WFC3IR_F160W.fits` from https://www.stsci.edu/~jayander/HST1PASS/LIB/PSFs/STDPSFs/WFC3IR/.
- Gaia DR3 positions (ESA TAP) of the three stars below; they are used only as starting positions.

## 3. Cases (all development: these stars were controls in #474)
1. **Real isolated star:** Gaia DR3 2755194341830999040 (G 18.03; visit if8a06030; star a of #474 development pair 0).
2. **No-galaxy double:** the two stars of #474 validation pair 3.
   - A: Gaia 2834388281729993472 (G 18.29, if8a4j030), the frame host;
   - B: Gaia 299642795747798272 (G 17.70, if8a38030).
   - Exposure k of B (k = 0–3 in rootname order) is background-subtracted, cut out ±15 px around its peak, and added
     to exposure k of A at an integer offset of (dx, dy) = (−3, +7) px from A's peak pixel. That is a 7.6 px (≈ 1.0″)
     separation at PA ≈ 111° in pixel axes, as in #474 pair 3.
   - The error is sqrt(ERR_A² + ERR_B²) where B is added, and the masks are combined. Fitting B with the PSF of A's
     detector region is an approximation; the size of the PSF change is reported.
3. **Pre-fixed extended-light injection (simulated) into case 2:** Sérsic n = 4, q = 0.8, r_eff 0.2″ (1.56 px at
   0.1283″/px), PA 30°, F160W AB 20, 21 and 22 mag (ZP 25.9311, #474 `runs/inj.json`).
   - Centre P = the fainter star + 0.4 × (brighter − fainter), from the case-2 initial peak positions.
   - Rendered on the 4× sub-grid and convolved with the STDPSF at P (the ePSF already contains the pixel response).
   - Poisson noise of the galaxy counts at 105.518 s, seed 1000 + 10 k + mag index.
   - Evidence strength: a simulation in real noise next to real stars. It tests absorption by the point-source fit,
     not a real lens galaxy. No known lens is used.

## 4. Fit and comparison quantities (identical in every case)
- **Per exposure:** `PSFPhotometry` with the STDPSF `GriddedPSFModel`, `fit_shape` 11 × 11 px, free x, y and flux
  per point source; the two stars of a pair are fitted as one group.
  - Background: the 3σ-clipped median of pixels 10–15 px from the star (pair: 12–18 px from the midpoint) is
    subtracted first.
  - Error: the ERR array. Masked: DQ bits 1, 2, 4, 8, 16, 32, 256, 512.
  - No Sérsic or other extended component is fitted. No PSF iteration.
- **Case 1:**
  - χ²/ν over unmasked pixels within r ≤ 5 px of the fitted centre;
  - radial residual z = Σ resid / sqrt(Σ ERR²) in annuli [0,1), [1,2), [2,3), [3,5), [5,8), [8,12) px, plus
    Σ resid / flux;
  - across the 4 exposures: the rms of the fitted positions on the sky (each FLT's WCS; mas and px at 0.128″);
  - the fitted-flux scatter against the median flux error;
  - aperture flux (r ≤ 6 px) / fitted flux.
- **Case 2:**
  - the same per-star χ²/ν and annuli;
  - the extended-light proxy at P: S = Σ resid and N = sqrt(Σ ERR²) in an r = 2 px aperture, combined over the 4
    exposures as z_P = Σ S / sqrt(Σ N²).
- **Case 3:**
  - z_P with the injection;
  - recovered fraction R = (S_inj − S_noinj) / F_ap, where F_ap is the noiseless injected flux in the same aperture.
    1 − R is what the point-source fit absorbs.
- FLT pixels are not resampled, so no drizzle correlation enters. Whether ERR is a calibrated per-pixel error is
  checked by case 1's z values.
- **Not compared:**
  - #474's Δχ² values: they come from a different pixel grid (drz) and a different statistic;
  - a model with versus without galaxy light, because no galaxy model is fitted.
  - The only same-pixel comparison is the same pair with and without the injection.

## 5. Decision rules (fixed here; computed by `scripts/w2_j2244_stpsf.py`)
- **H — hook-up:** all 12 FLTs read, and the STDPSF evaluated at full-frame positions.
  - Case-1 sky-position rms ≤ 0.1 px: hook-up passes.
  - > 0.3 px: stop (mapping problem). In between: passes, flagged.
- **S — star modelled to the noise:** in every case-1 exposure, χ²/ν(r ≤ 5) ≤ 2 and |z| ≤ 3 in every annulus.
- **P — no false extended-light signal:** |z_P| < 3 for case 2 (one pair; not a false-positive rate).
- **I — injection recovered and not mostly absorbed:** z_P ≥ 5 and R ≥ 0.5 at 20 mag. z_P and R are reported at 21
  and 22 mag.
- **Worth a formal validation:** H, S, P and I all pass, and z_P ≥ 5 at 21 mag.
  - #474 F160W recovered no injection at 19–23 mag, but with a different statistic, so this is a pointer, not a
    like-for-like improvement.
- **If S fails:** the default official model does not reach the noise on real stars at the target's brightness in
  these exposures. The failing annuli locate the problem. No perturbation, stack or other route is tried in this
  round.
- **Stop:** the download cap; a mapping failure (H); an unreadable model; or a need to write a PSF solver.

## 6. Out of scope
- #474 numbers and verdicts (unchanged); no new thresholds for #474.
- J2244−0550.
- DESI classifier, X4, W3 / ATLAS and other paused units.
- Outside contact, accounts or schedules.

## Errata after the run (from the independent review; no setting or rule changed)
- **"No exposure of its visit is read" (header):** no pixel of J2244−0550's visit was read. But the ledger's live
  tests (§0) used its drz product `if8a5k030`: a 2880 B primary-header read, and an un-ranged request cut after the
  HTTP header (RESULT §7).
- **§1 "zero beyond r ≈ 12.5 px":** the model is a 25 × 25 px square (±12 px on the axes) and is tapered at large
  radii (`posthoc.json`).
- **§3 case 3:** the injected magnitude is normalised within the ±12 px box before PSF convolution.
- **§5 H** checks sub-pixel consistency across dithers. A constant subarray offset cancels in it.
- **§5 I:** the no-galaxy residual alone (z_P 19) satisfies the z condition, so I rests on R.
- **Not recorded:** photutils' fit flags.
- **§0 (ledger), corrected in RESULT §1:**
  - the TLS diagnosis of www.stsci.edu ran outside the ledger (an estimated `external` record);
  - error and redirect responses were counted at the TLS allowance only until the post-review fix;
  - the live "cut mid-stream" happened before the first body read.
