# El Gordo arc difference imaging, PEARLS 1176 o241 (2022-07-29) → VENUS 6882 o051 (2026-08-05): result

Plan: `PLAN.md`, committed before any run (6ba237c). Addenda 1–2 were written before the controls they
define were measured. This is an owner-approved exploratory pilot: one box and one epoch pair. Earlier fields fixed
the tool settings, so this is **not a blind test**. Every number is `derived` unless labelled; thresholds are
ASSUMPTIONs.

## Data
- **Images (`observed`):** F150W and F444W of both epochs, S3 byte ranges.
- **Pipeline:** both epochs have `CAL_VER` 3.0.0 and CRDS `jwst_1584.pmap` (`o241_o051/headers.json`), so
  pipeline-version
  differences are excluded.
- **Orientation:** PA_V3 is 276.25° and 277.09° (`headers.json`), a difference of 0.84°.
- **Box:** 24″ × 24″ at (15.7300, −49.2530). It holds 11 Caminha+2023 images of 6 systems.
- **Data-quality gate passed:** epoch-1 images are 100 % finite in the box; difference maps 99.2 % / 99.0 %
  (`headers.json` `maps_finite_fraction`).
- **Registration (`reg6882/`):** catalogue offsets between the epochs at the box are ≤ 0.02″ (F150W
  (0.016, −0.016)″, M0; F444W (0.004, −0.012)″, M1).

## Sensitivity (what this search can see)
- **Point injections per band (`o241_o051/summary.json`; 100 positions on epoch-1 pixels above the 90th
  percentile):**
  - F150W 0.99 / 0.99 / 0.92 at 27.5 / 28.0 / 28.5 AB (0.74 / 0.76 / 0.70 away from bright cores);
  - F444W 0.94 / 0.71 / 0.36 (0.69 / 0.43 / 0.18 away from bright cores).
- **Full-chain control (`inject_grid_seed0/inject_review.json`; 80 stars, seed 0):**
  - Each star is a Gaussian point source added to epoch 2 of both bands on galaxy, arc or ICL pixels, 1.5″ apart.
  - Recovered = hit, unflagged and morph `pass`, over the whole box.
  - **26.0 AB: 23 / 27 = 0.85 (gate ≥ 0.80: passed)**; 27.5 AB: 21 / 27 = 0.78; 28.0 AB: 14 / 26 = 0.54.
  - These apply to this box and these positions only.
- **Not covered:**
  - **Fadings are not calibrated.** All injections are brightenings, while a W3-type signal (a compact image that
    dims or vanishes) is a fading. The chain detects either sign, but no fading injection was run.
  - Hits within 0.8″ of bright cores (32 of 50) are not reviewed (chain stage 2, unchanged). Changes on bright
    galaxy cores and stars are therefore not searched.
  - Only F150W and F444W, one 24″ box, and two dates 4 years apart. A transient that rose and faded between
    2022-07 and 2026-08 is invisible.
  - Sources fainter than about 28 AB, and extended changes, are largely missed.

## Screen
- 50 hits; 18 not `near_bright_core`.
- All-hits morphology: 28 rosette, 10 dipole, 6 edge, 1 sign_flip, 5 pass.
- **Review list (`review.json`):** 4 entries (h024, h031, h043, h049), all priority A (in the arc mask). The cap is
  15, so all were vetted. h024 and h049 are negative in both bands (`hits.ecsv` sign `-`); h031 and h043 are
  `mixed`, with the largest |S| negative in F150W.

## Vetting (sheets looked at: `sheet_review.png`, `vet/stamps_*.png`)
- Forced photometry: `vet/lightcurves.ecsv`, 0.15″ aperture, conservative random-aperture errors.
- **Spike geometry (`spike/spike.json`, `scripts/elgordo_spike_check.py`):**
  - All four lie 2.6–3.9″ from the box's brightest point source (15.72722, −49.25394).
  - Three are on its diffraction spikes (spikes = rays with an opposite partner): h024 in both bands, h043 in both,
    h049 in F150W (F444W 2.0° off). h031 is not on a spike.
  - After the final review, S is sampled across the matched spike line. On that line no entry meets the
    rotation-dipole rule. The first version sampled across the entry's own direction and reported a dipole for
    h043 in F150W; that version is superseded.
  - With PA_V3 0.84° apart, a spike moves
  sideways by about 0.05″ at that distance.

| Entry | Nearest image | Forced photometry 2022 → 2026 | Check | Class |
|---|---|---|---|---|
| h024 | 1.36″ from 11b (z 4.33) | F150W 502 → 473 nJy, F444W 296 → 283 (−1.7σ each) | on a spike in both bands; no change ≥ 2σ | no significant change (spike residual in S) |
| h031 | 0.06″ from 11b | F150W 133 → 141, F444W 91 → 91 (+0.5σ, 0.0σ) | no change | no significant change |
| h043 | 0.68″ from 9b (z 4.32) | F150W 110 → 85 (−1.4σ); F444W 655 → 624 (−4.1σ in `vet/lightcurves.ecsv`) | F150W: on the vertical spike, no change ≥ 2σ. F444W: on-spike controls 1 / 8 (`undetermined`), then bright-source controls: raw z −5.3 (ratio 0.954) against a control median z −9.8 and robust std 5.9 → calibrated z 0.77, tail fraction 0.45 (`vet/bright_ctl30/control_z.ecsv`) | **explained:** program-to-program photometric scatter of bright compact sources (F444W); no significant F150W change |
| h049 | 3.07″ from 11b | F150W 34 → 2 (−1.9σ); F444W 37 → 10 (−3.5σ in `vet/lightcurves.ecsv`) | on the horizontal spike (F150W); on-spike controls: 6 / 8 fade at ≥ 2σ in F444W (`spike/onspike_controls.json`, `onspike-classify`) | **explained by the spike shift** |

- **σ values differ between runs.** `hit_lightcurve.py` errors are the robust spread of random apertures in each
  run's boxes, so the same flux change has slightly different σ in each run:
  - h049 F444W −3.5σ in `vet/` and −3.0σ in `spike/ctl/`;
  - h043 F444W −4.1σ in `vet/` and −3.5σ in `spike/ctl/`; `compare` gives raw −5.3 with its own errors.
  The classes come from the rule run named in each row.
- **Deviations:**
  - The spike check first chose a different source in F150W; this was fixed before its result was used.
  - After the final review: spikes must have an opposite partner, which drops a 121.5° F444W ray, so h031 is no
    longer on a spike. S is now sampled across the matched spike line. The on-spike rule is now committed code
    (`onspike-grid` reproduces the measured positions exactly; `onspike-classify` reproduces both classes).
  - The bright-source control box was widened from 16″ to 30″ because 17 controls are fewer than the tool's
    minimum of 20; this was done before the control z values were seen.
  - The addendum 1 wording first left out h043; it was corrected before measuring.
- **Bright compact sources fade systematically between programs:** the 20 controls with F444W 300–1300 nJy give a
  median z of −9.8 in 0.15″ apertures. Any single-aperture fading of a bright compact source in this pair must be
  compared with that.

## Outcome
- **No real, sky-fixed flux change without an ordinary or instrumental explanation** was found in this box.
- **Nothing was sent to `/vet-candidate`.**
- **No SN or variable was found either.**
- This is a limited-sensitivity search of one box: brightenings are recovered at about 0.85 / 0.78 / 0.54 at
  26.0 / 27.5 / 28.0 AB, and fadings are not calibrated. **It gives no exclusion and no abundance limit.**

## Regenerate
- `arc_diff.py`: the command is in `o241_o051/summary.json`; the inject-grid command is in
  `inject_grid_seed0/summary.json`. The maps were written to the session scratchpad and are not stored.
- `arc_diff_morph.py --hits … --maps … --all-hits`, then `arc_chain.py review --field elgordo …` and
  `inject-review --field elgordo …`.
- `reg_offsets.py --pair 6882:… --ref 1176:… --bands f150w f444w --offsets-at 15.7300 -49.2530`.
- `hit_lightcurve.py` with the four entries, then the on-spike positions (`spike/ctl`) and the
  `hit_lightcurve_controls.py pick --half-arcsec 30 --flux-range 300 1300 --n 60 --seed 1` positions, with
  `compare` on the F444W pair.
- `elgordo_spike_check.py headers --uri <4 i2d> --maps …`; `spike --maps … --review … --dpa 0.842`;
  `onspike-grid --source 15.727222 -49.253939 --entry h049 … --entry h043 …`, then `hit_lightcurve.py` on those
  positions and `onspike-classify --lc spike/ctl/lightcurves.ecsv --entry h049 --entry h043`.
