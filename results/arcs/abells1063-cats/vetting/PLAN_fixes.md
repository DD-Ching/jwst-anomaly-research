# Plan: fix the four review issues and recompute affected results (owner-approved 2026-10-16; written before any re-run)

Scope: correct existing results only. No new exposure, no 116.28VK stacking, no whole cubes, no new line of work.
Budget: 2 h. If time runs out, the affected numbers are marked invalid or pending re-computation, not current.

**Data.**
- Re-fetch the two MUSE SODA cutouts already used: ADP.2017-03-28T12:46:01.331 at 1557 and
  ADP.2017-03-23T15:58:03.937 at 1316, each with a 2.4″ radius.
- Steps 2–4 also need the JWST stamps they used before: the same GLIMPSE 3293 i2d products, bands, positions and
  8″ size, read by S3 byte range.
  - Those stamps were never stored, because earlier sessions were ephemeral.
  - Re-reading the identical stamps is the only way to recompute the deblend and the census.
  - It adds no product, band or region.

**Fixes (code, decided now; none is tuned on the 6000 Å peak).**
1. *Spectra noise rescale* (`arcs_vet_spectra.py`): the MAD rescale of the boxcar S/N is computed only over valid,
   unmasked channels. These are the same channels that are tested; channels flagged by fix 2 and by the
   sky-line flag are excluded.
2. *Per-channel validity* (`arcs_vet_spectra.py`):
   - A channel is valid only when every aperture spaxel has finite DATA and STAT > 0.
   - Invalid channels get no flux, and every boxcar that touches one is excluded, like a sky-line channel.
   - The sky-line flag (variance > 1.5 × its 101-channel running median) is evaluated on valid channels only.
   - The running-median continuum uses valid channels only (NaN-aware median).
3. *Deblend flux scale* (`arcs_vet_deblend.fit_models`):
   - Old rule: flux start and upper bound came from the filtered aperture sum, floored at 1e-6, which pins faint
     rows whose sum is ≤ 0.
   - New rule: the start is that sum if it is > 0, else the sum of the positive pixels. The upper bound is
     10 × max(sum, sum of positive pixels).
   - Everything else is unchanged: models, BIC, the class rule `deblend_class`, the gate rule.
4. *Census sample consistency* (`arcs_vet_census.py`): a count that removes a row also removes that row's
   probability from the Poisson-binomial null. `k2_core_flux` (2184 removed) becomes P(≥ 2) over the 23 other rows.
   `k4_screen` and `k3_bgsub` are counts over all 24 rows and stay as they are.

**Recompute (only what the fixes touch), keeping before/after.**
- Fixes 1–2: `spectra.json`.
- Fix 3: everything that calls `fit_models` (observed fits and gate injections):
  - `deblend.json` (1316; the 2184 SED part is re-run with it but does not use the fit);
  - `deblend_1557.json`;
  - `clumps_1557.json`, if it reads the step-3 fit;
  - `census.json` / `.ecsv` (classes, gates, f̂, P).
  The noise injections (moments, no fit) are not affected.
- Fix 4: `census.json` sample block.
- Old files are kept as `*_before_fixes.*`, and `fixes_diff.json` lists every changed number.

**Verification (unit tests, written with the fixes):**
- invalid channels are excluded and cannot produce a detection;
- a sky-line-contaminated channel does not inflate the rescale of clean channels;
- a faint source whose filtered sum is ≤ 0 gets a non-pinned fit, with flux not at the bound and ΔBIC not
  trivially 0;
- removing a row removes its probability from the tail.

**6000 Å check for 1316 (decided now, run after fixes 1–2).**
- **Validity:** report whether the peak channel range is valid and unflagged after the fixes.
- **Narrow-band image:**
  - Map = sum of the 6 channels at the peak minus the mean of two 30-channel side bands (±20–50 channels).
  - Is the map peak within 0.5 × FWHM (0.72″) of the target? This is a pass or fail.
- **Neighbouring background:**
  - Six apertures of the same radius, centred 1.5″ from the target at 60° steps, inside the 2.4″ cutout.
  - Their boxcar S/N at the same channels, with each aperture's own rescale.
  - The residual is local to the target only if none of the six has S/N ≥ 3 there.
- **False-alarm rate:**
  - The observed S/N is a local threshold. It is calibrated empirically: count how many boxcar positions in the
    six offset apertures, over the valid unmasked channels of the C III] window, reach the target's peak S/N.
  - This count gives the expected number of chance peaks in one aperture per window search (look-elsewhere
    included).
  - The window is searched once per target and line; there are 3 windows (1557 [O II]; 1316 C III] and C IV).
- **Reading (fixed now):**
  - A peak passes only if all of these hold: it is ≥ 5 after the fixes, valid, spatially coincident, not in the
    offset apertures, and its expected chance count is < 0.05.
  - A peak that passes is only a **line to be checked**. One line cannot identify the species, so no redshift is
    claimed. It is not evidence for unusual lensing.
  - Otherwise it is recorded as not significant under the calibration.
