# W2 J2244 follow-up: the STScI default WFC3/IR F160W PSF hooks up to our FLT exposures, but a no-galaxy ~1″ pair of G ≈ 18 stars still leaves a residual at the lens position (z 19); not worth a formal validation, and the image-modelling line is paused

- **Plan:** `PLAN.md`, frozen at 908ecca before any fit.
  - One analysis run (5307add) and an identical rerun for the recomputation check (§6). No setting changed.
  - Errata after the review are at the end of the PLAN.
- **Type:** development feasibility test on stars already used as controls in #474. Not an independent validation.
- **J2244−0550 was not measured**, and no pixel of its visit was read. The ledger's live tests did touch its drz
  product header (deviation, §7).
- **Labels:**
  - FLT pixels, the STDPSF model and Gaia positions: observed;
  - fits, residuals, z and R: derived (`result.json`; post hoc additions in `posthoc.json`);
  - the Sérsic injections: simulated;
  - thresholds: ASSUMPTIONs (`scripts/w2_j2244_stpsf_rules.py` `Params`).

## 1. Download control (PLAN §0; `net_ledger.jsonl`)
- **Through the ledger:** every data request went through `scripts/net_ledger.py`. Each one reserves first, counts
  the bytes received (header, body and an 8 kB TLS allowance per request) and is cut off before the cap. The JSONL
  ledger is shared and survives restarts.
- **One exception:** the TLS diagnosis of www.stsci.edu (curl / openssl) ran outside the ledger. It was entered
  afterwards as an estimated 20 kB `external` record.
- **Checked before any image request:**
  - unit tests;
  - one 2880 B live request (counted 11651 B);
  - a live stop test on a 40 kB ledger (`net_ledger_stoptest.jsonl`). It was refused before the request (that
    refusal left no record at the time), and an underestimated request was cut, before its first body read.
- **Total: 105.3 MB of the 1 GB cap.**
  - First fetch 52.8 MB: 12 FLTs, the PSF model, 3 drz headers, the notebook page and Gaia.
  - 52.1 MB for the recomputation re-fetch (§6).
  - 0.38 MB for the post hoc model re-fetch (`posthoc.json`).
  - Dependencies: 0 B (`uv --offline` from the local cache).
- **Ledger gaps found in the review, fixed after the run** (offline tests added, no data re-read):
  - error and redirect responses were counted at the 8 kB allowance only, which here means two GitHub 403s, a few kB.
    They are now counted in full, and redirects are not followed.
  - refusals are now recorded;
  - a cut after body data has arrived is now tested;
  - the Windows lock now waits at most 300 s.
- **Not covered:** `external` records are not capped, and bytes the operating system buffers after a cut-off are
  not counted.
- **www.stsci.edu certificate:** the host's chain ends at "emSign Root CA − G1", which the container bundle lacks.
  Verification used the proxy bundle plus certifi's Mozilla roots (`NET_LEDGER_CAFILE`); it was never disabled.

## 2. The official entry (PLAN §1)
- **Pinned page:** the rendered page of 2026-10-23, "Updated: 20 May 2026" (sha256 in `sources.json`). The GitHub
  commit lookup was refused (HTTP 403).
- **Route chosen:** §2.2.1 "Default Models". This is the WFC3-team STDPSF grid, read by photutils
  `GriddedPSFModel.read(format="stdpsf")` and fitted with `PSFPhotometry`.
- **How it fits our data** (model properties in `posthoc.json`):
  - It is an exposure-frame model, so it was fitted on the 4 FLT exposures of each visit, not on the drizzled image.
  - F160W is covered by a 3 × 3 detector grid, 4× oversampled; the sub-pixel phase follows the fitted position.
  - IRSUB512 maps to the full frame through LTV1 = LTV2 = −251. A 5 px error in the frame would change the model by
    4.3 × 10⁻⁴ of its peak.
  - The centre PSF sums to 0.943 within r ≤ 3 px, 1.008 within 6 px and 1.026 within 10 px. The 9 grid PSFs sum to
    1.027–1.037. So the fitted flux is about the flux within r ≈ 6 px (0.8″), in e⁻/s.
  - The model is a 25 × 25 px square (±12 px, 1.5″ on the axes) and is tapered at large radii, per the notebook.
    Light beyond it is not modelled.
  - No PSF solver was written and no extra data were needed. The model was used as is (no iteration, nothing fitted
    to the controls).
- **Routes not used:**
  - focus or spatial perturbation (hst1pass), which needs at least 10 bright isolated stars per exposure;
  - stacks, a different method;
  - drizzled empirical models, which are not yet officially supported.

## 3. Cases (PLAN §3; `result.json`, `cases.png`)

| case | what | measured | rule → outcome |
|---|---|---|---|
| H | case-1 star on the sky across 4 dithered FLTs | rms 1.29 mas = 0.010 px (offsets 0.4–1.6 mas) | ≤ 0.1 px → **pass** |
| 1 | real isolated star, Gaia G 18.03 (fitted 8074–8254 e⁻/s) | χ²/ν (r ≤ 5 px) 30.9–38.6; flux scatter 0.95 % against 0.15 % formal | **S fails** (χ²/ν ≤ 2 and \|z\| ≤ 3 needed) |
| 2 | no-galaxy double: two real stars ~1″ apart (fitted 7.3–8.7 px) | z_P = 19.0 at the lens position (S 12.5–29.5 e⁻/s per exposure, N 1.8–2.3) | **P fails** (\|z_P\| < 3 needed) |
| 3 | case 2 + Sérsic n = 4, r_eff 0.2″ at P | 20 / 21 / 22 mag: z_P 107.1 / 55.7 / 34.5, R 0.909 / 0.892 / 0.918 | **I passes** (z_P ≥ 5 and R ≥ 0.5 at 20 mag) |

- **What rule I rests on:** the no-galaxy residual alone gives z_P = 19, so its z condition holds at any magnitude.
  - Rule I therefore rests on R.
  - Without the baseline, z is 90.1 / 37.6 / 15.8 at 20 / 21 / 22 mag (post hoc, `posthoc.json`).
- **H checks sub-pixel consistency only.** A constant subarray offset would cancel in it, because the fitted position
  returns to the sky through the same offset.
- **Case 1, residual by annulus** (z, then residual / fitted flux, ranges over the 4 exposures):

| annulus (px) | 0–1 | 1–2 | 2–3 | 3–5 | 5–8 | 8–12 |
|---|---|---|---|---|---|---|
| z | −0.5 to 7.3 | −4.7 to 3.0 | −4.8 to 0.3 | 10.0 to 18.2 | 21.9 to 26.1 | 32.3 to 34.8 |
| resid / flux | −0.1 to 1.0 % | −0.3 to 0.2 % | −0.2 to 0.0 % | 0.3 to 0.6 % | 0.9 to 1.1 % | 1.6 to 1.7 % |

- **Reading the table:**
  - In the core (r < 3 px) the sums are near zero, but the pixels are not: χ²/ν 31–39 within 5 px.
  - From 3 px (0.4″) outward the data are brighter than the model in every exposure, by 3.0–3.1 % of the fitted flux
    within 12 px (sum of the three outer annuli, post hoc).
  - The background ring (10–15 px) lies inside the star's halo, so the background is over-subtracted. That pushes
    the residuals down, so 3.0–3.1 % is a low estimate, and the failures of S and P are not caused by it.
  - Aperture flux within 6 px / fitted flux is 0.92–1.01.
- **Contact sheet** (looked at, `cases.png`; colour scales are not stored in JSON):
  - the star leaves a structured core residual (colour scale ±12.35 e⁻/s) and a faint positive halo along the
    diffraction-spike directions;
  - the pair leaves the same core residuals and a positive residual between the stars;
  - the injected Sérsic sits between them, as placed.
- **What makes up the pair residual at P is not separated:**
  - its sign matches the star's wing excess;
  - but P lies 3.1–3.8 px from star A, inside A's structured residual zone, and S varies 2.4× across the four
    exposures (12.5–29.5 e⁻/s), while the star's outer excess is steady;
  - the pair combines two visits (#474's untested "two PSFs in one pair" reason). The reported 0.94–0.95 %-of-peak
    PSF change covers only the detector-position difference, not focus.
  - Wing excess, core-shape or sub-pixel mismatch, and the two-visit difference may each contribute.
- **Case 3 against case 2:**
  - The point-source fit absorbs only about 10 % of the injected light in the aperture (1 − R = 0.08–0.11). The PSF
    is fixed here, and no PSF iteration takes up galaxy light.
  - The no-galaxy residual at P (S 12.5–29.5 e⁻/s per exposure) is about equal to the aperture flux of the 22 mag
    injection (17.8–18.8 e⁻/s). The 21 mag aperture flux is 44.8–47.3 e⁻/s.

## 4. Answers
1. **Does the official route fit our data?**
   - Mechanically, yes: F160W model, IRSUB512 subarray, sub-pixel position, single FLT exposures. Positions agree
     across dithers to 1.3 mas.
   - The subarray offset itself is not tested by that check; a 5 px error would change the model by ≤ 4.3 × 10⁻⁴ of
     its peak.
   - It does not cover two of #474's untested possible reasons: it does not apply to drizzled images, and light
     beyond its ±12 px square is not modelled.
   - Its documented remedy for star-to-star shape differences (hst1pass focus perturbation) needs at least 10 bright
     isolated stars per exposure, which these fields lack.
2. **What improved and what still fails?**
   - **Measured here:**
     - positions are consistent to 0.010 px;
     - injected galaxy light is mostly not absorbed by the point-source fit (R 0.89–0.92 at 20–22 mag).
   - **Uncorrelated noise holds by construction, not by measurement:** FLT pixels are not resampled. The planned check
     of the error calibration through case 1's z values cannot decide this, because model mismatch dominates them.
   - **Still failing:**
     - a real star matched to the target in Gaia G (the F160W match is not checked) is not modelled to the noise
       (χ²/ν 31–39);
     - its wings at 3–12 px are 3.0–3.1 % brighter than the model;
     - a no-galaxy ~1″ pair leaves a residual at the lens position of z 19, about the aperture flux of a 22 mag
       compact injection.
   - **Not like-for-like with #474:** the statistic here is an aperture residual, not Δχ², and the pixel grid differs.
     The two are not compared.
   - **Not enough for a formal validation** under the frozen rule: S and P fail.
3. **Is the problem located?**
   - **Only partly described; the cause is not known.**
   - The single star shows structured core residuals and a wing excess. The pair residual at P matches the wing excess
     in sign but not in stability, and it sits inside the companion's core zone.
   - Wing excess, core-shape or sub-pixel mismatch, and the two-visit PSF difference are not separated, and no test
     in this round could separate them.
4. **What it takes to reproduce:**
   - **Inputs:** `sources.json` (URLs, sizes, sha256).
   - **Scripts:** `scripts/net_ledger.py`, `scripts/w2_j2244_stpsf_fetch.py`, `scripts/w2_j2244_stpsf.py`,
     `scripts/w2_j2244_stpsf_rules.py`.
   - **Settings:** the PLAN, and `settings` in `result.json`.
   - **Environment:** photutils 3.1.1, astropy 8.0.2, numpy 2.5.4, scipy 1.18.1, matplotlib, in an isolated
     environment. The CA note for www.stsci.edu is in §1.
   - **Cost:** about 52 MB of download and about 11 s on 2 cores per run (`cost.json`).

## 5. Consequence (owner rule for a failed check)
- **The self-built image-modelling line for J2244−0550 is paused.** That covers the #474 lenstronomy chain and this
  official-PSF hook-up.
- No larger PSF, star library or budget is proposed.
- #474's numbers and verdicts are unchanged.

## 6. Recomputation check before clean-up (added after the run, for the owner's clean-up rule)
- `scripts/w2_j2244_stpsf_fetch.py` re-fetched all 14 inputs through the ledger. All sha256 matched the recorded ones,
  including the Gaia response.
- A rerun in a fresh directory gave a `result.json` and `cases.png` byte-identical to the committed ones.
- The scratch data were deleted afterwards: FLTs, model file and isolated environments.

## 7. Deviations
- **The target's product was touched by the ledger tests.** Both used J2244−0550's own drz product (`if8a5k030`):
  - the 2880 B live request read its FITS primary header;
  - the stop test's un-ranged request was cut after the HTTP header.
  - No pixel data were read or used.
- **Post hoc:**
  - the recomputation rerun (§6);
  - the model re-fetch and `posthoc.json` (model properties and the baseline-subtracted z);
  - the sum of the three outer annuli;
  - the comparison of the pair residual with the injected aperture fluxes.
  - None of these enters a decision.
- **Not in the plan text:** the injected magnitude is normalised within the ±12 px box before PSF convolution (R is
  self-consistent).
- **Not recorded:** photutils' fit flags.
- **Cosmetic:** the titles in `cases.png` row 3 overlap.
- The GitHub commit of the notebook is not recorded (HTTP 403). The page's sha256 and date stand in for it.
