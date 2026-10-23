# W2 J2244 follow-up: the STScI default WFC3/IR F160W PSF hooks up to our FLT exposures, but at the target's brightness it leaves wing residuals that make a no-galaxy pair look like extended light (z 19); not worth a formal validation, and the image-modelling line is paused

- **Plan:** `PLAN.md`, frozen at 908ecca before any fit. One run (5307add); no setting was changed after it.
- **Type:** development feasibility test on stars already used as controls in #474. Not an independent validation.
  J2244−0550 was not measured, and no exposure of its visit was read.
- **Labels:**
  - FLT pixels, the STDPSF model and Gaia positions: observed;
  - fits, residuals, z and R: derived (`result.json`);
  - the Sérsic injections: simulated;
  - thresholds: ASSUMPTIONs (`scripts/w2_j2244_stpsf_rules.py` `Params`).

## 1. Download control (PLAN §0; `net_ledger.jsonl`)
- Every request went through `scripts/net_ledger.py`: reserve first, count the bytes received (header, body and an
  8 kB TLS allowance per request), cut off before the cap. The JSONL ledger is shared and survives restarts.
- **Checked before any image request:** unit tests; one 2880 B live request (counted 11651 B); a live stop test on
  a 40 kB ledger, refused before the request and cut mid-stream (`net_ledger_stoptest.jsonl`).
- **Total: 104.9 MB of the 1 GB cap.**
  - First fetch 52.8 MB: 12 FLTs, the PSF model, 3 drz headers, the notebook page and Gaia.
  - The rest is the re-fetch for the recomputation check (§6).
  - 20 kB is estimated: a TLS diagnosis with curl / openssl.
  - Dependencies: 0 B (installed with `uv --offline` from the local cache).
- **Not counted:** bytes the operating system buffers after a cut-off; none occurred outside the stop test.
- **www.stsci.edu certificate:** the host's chain ends at "emSign Root CA − G1", which the container bundle lacks.
  Verification used the proxy bundle plus certifi's Mozilla roots (`NET_LEDGER_CAFILE`); it was never disabled.

## 2. The official entry (PLAN §1)
- **Pinned page:** the rendered page of 2026-10-23, "Updated: 20 May 2026" (sha256 in `sources.json`). The GitHub
  commit lookup was refused (HTTP 403).
- **Route chosen:** §2.2.1 "Default Models". This is the WFC3-team STDPSF grid, read by photutils
  `GriddedPSFModel.read(format="stdpsf")` and fitted with `PSFPhotometry`.
- **How it fits our data:**
  - It is an exposure-frame model, so it was fitted on the 4 FLT exposures of each visit, not on the drizzled image.
  - F160W is covered by a 3 × 3 detector grid, 4× oversampled; the sub-pixel phase follows the fitted position.
  - IRSUB512 maps to the full frame through LTV1 = LTV2 = −251. The 5 px ambiguity of the frame changes the model by
    ≤ 4.3 × 10⁻⁴ of its peak.
  - The model sums to 1.008 within r ≤ 6 px and 1.03 over its 25 × 25 px footprint, so the fitted flux is about the
    flux within 0.8″ (e⁻/s).
  - It is tapered to zero beyond r ≈ 12.5 px (1.6″), so the halo is not modelled.
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
| 1 | real isolated star, G 18.03 (fitted 8074–8254 e⁻/s) | χ²/ν (r ≤ 5 px) 30.9–38.6; flux scatter 0.95 % against 0.15 % formal | **S fails** (χ²/ν ≤ 2 and \|z\| ≤ 3 needed) |
| 2 | no-galaxy double: two real stars ~1″ apart (fitted 7.3–8.7 px) | z_P = 19.0 at the lens position (S 12.5–29.5 e⁻/s per exposure, N 1.8–2.3) | **P fails** (\|z_P\| < 3 needed) |
| 3 | case 2 + Sérsic n = 4, r_eff 0.2″ at P | 20 mag z_P 107.1, R 0.909; 21 mag 55.7, 0.892; 22 mag 34.5, 0.918 | **I passes** (z_P ≥ 5 and R ≥ 0.5 at 20 mag) |

- **Case 1, residual by annulus** (z, then residual / fitted flux, ranges over the 4 exposures):

| annulus (px) | 0–1 | 1–2 | 2–3 | 3–5 | 5–8 | 8–12 |
|---|---|---|---|---|---|---|
| z | −0.5 to 7.3 | −4.7 to 3.0 | −4.8 to 0.3 | 10.0 to 18.2 | 21.9 to 26.1 | 32.3 to 34.8 |
| resid / flux | −0.1 to 1.0 % | −0.3 to 0.2 % | −0.2 to 0.0 % | 0.3 to 0.6 % | 0.9 to 1.1 % | 1.6 to 1.7 % |

- **Reading the table:**
  - In the core (r < 3 px) the sums are near zero, but the pixels are not: χ²/ν 31–39 within 5 px.
  - From 3 px (0.4″) outward the data are brighter than the model in every exposure, by 3.0–3.1 % of the fitted flux
    within 12 px (sum of the three outer annuli).
  - Aperture flux within 6 px / fitted flux is 0.92–1.01.
- **Contact sheet** (looked at, `cases.png`):
  - the star leaves a structured core residual (colour scale ±12.35 e⁻/s) and a faint positive halo along the
    diffraction-spike directions;
  - the pair leaves the same core residuals and a positive residual between the stars;
  - the injected Sérsic sits between them, as placed.
- **Case 3 against case 2:**
  - The point-source fit absorbs only about 10 % of the injected light in the aperture (1 − R = 0.08–0.11). The PSF
    is fixed here, and no PSF iteration takes up galaxy light.
  - But the no-galaxy residual at P (S 12.5–29.5 e⁻/s per exposure) is as large as the aperture flux of the 22 mag
    injection (17.8–18.8 e⁻/s). Without a no-galaxy version of the same pair, which a real target never has, a
    compact galaxy fainter than about 21–22 mag cannot be told apart from this residual.
- **Approximation in case 2:** star B is fitted with the PSF of its new detector position, which differs from the PSF
  at its own position by 0.94–0.95 % of the peak. The case-1 star shows the same wing excess without any relocation.

## 4. Answers
1. **Does the official route fit our data?**
   - Mechanically, yes: F160W model, IRSUB512 subarray, detector and sub-pixel position, single FLT exposures.
     Positions agree across dithers to 1.3 mas.
   - It does not cover what #474 needs most. It does not apply to drizzled images, and it is tapered beyond 1.6″,
     so the halo is not modelled.
   - Its documented remedy for star-to-star shape differences (hst1pass focus perturbation) needs at least 10 bright
     isolated stars per exposure, which these fields lack.
2. **What improved and what still fails?**
   - **Measured here:**
     - positions are consistent to 0.010 px;
     - galaxy light is mostly not absorbed by the point-source fit (R ≈ 0.9 at 20–22 mag);
     - the noise is per-pixel and uncorrelated, because FLT pixels are not resampled.
   - **Still failing:**
     - a real star at the target's brightness is not modelled to the noise (χ²/ν 31–39);
     - its wings at 3–12 px carry 3.0–3.1 % more light than the model;
     - for a ~1″ double that excess appears between the images as a false extended-light signal of z 19, about the
       aperture flux of a 22 mag compact galaxy.
   - **Not like-for-like with #474:** the statistic here is an aperture residual, not Δχ², and the pixel grid differs.
     The two are not compared.
   - **Not enough for a formal validation** under the frozen rule (H, S, P and I all needed): S and P fail.
3. **Is the problem located?**
   - **Partly.** The false signal comes from light the default model lacks at r ≥ 3 px. That light is present for a
     single isolated star, so it is not a pair-construction effect.
   - Which mechanism dominates is not separated: the documented taper, which removes the outer wings, or a
     focus / shape mismatch that the notebook fixes with perturbation (not available here). Neither was tested in
     this round.
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

## 6. Recomputation check before clean-up
- `scripts/w2_j2244_stpsf_fetch.py` re-fetched all 14 inputs through the ledger. All sha256 matched the recorded ones,
  including the Gaia response.
- A rerun in a fresh directory gave a `result.json` and `cases.png` byte-identical to the committed ones.
- The scratch data were deleted afterwards (52 MB of FLTs, the model file and the isolated environment).

## 7. Deviations
- None from the frozen plan's cases, settings or rules.
- **Post hoc, descriptive only:** the comparison of the false residual with the 22 mag aperture flux (§3), and the sum
  of the three outer annuli. Neither enters a decision.
- The GitHub commit of the notebook is not recorded (HTTP 403). The page's sha256 and date stand in for it.
