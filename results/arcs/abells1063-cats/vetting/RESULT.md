# Vetting of the Abell S1063 `arcs` anti rows (2026-10-16)

## Summary (owner-corrected wording, 2026-10-16; the step records below are kept unchanged)
- **Sample statistic:** 4 background anti rows vs 0.52 expected from intrinsic shape noise (screen, P 0.0021;
  Poisson-binomial 0.00082). Under the step-4 clump-aware null, P(≥ 4) = 0.026. This number is an **exploratory
  result that depends on a simplified ordinary model**:
  - the clump-pair fraction is uncertain (f̂ 0.238 from 21 rows; the mixture null at its 95 % upper value f = 0.472 gives P = 0.093);
  - point clumps with isotropic separations are assumed, and the κ < 1 formula is also used for g > 1 rows;
  - 2004's κ ≈ 1 degeneracy of the predicted axis is not in the model; it is counted as anti with a generic row
    probability (0.191). Counting only the other three anti rows (1316, 1557, 2184), P(≥ 3) = 0.116
    (`census.json` `k3_bgsub`, which still includes 2004's row probability).
  The clump idea came from looking at 1316 and 1557, so the test is post hoc.
- **What is not excluded:**
  - The tests did not exclude ordinary morphology (clumps, blends, intrinsic shapes), measurement noise or lens-model
    error in general. They only showed that the specific models tested do not reproduce the anti axes.
  - GLAFIC, Sharon, Keeton and CATS agree, but their errors are not fully independent: they share the same image
    constraints and similar parametrisations.
- **Per row:**
  - **2004:** at κ ≈ 0.95–1.00 in all three team models, where the predicted axis is degenerate.
  - **2184:** compatible with the measurement-noise model tested (core-flux injections, P_anti 7 / 60; post hoc).
  - **1316, 1557:** **pending vetting** (clump pairs along the anti axis). Not anomaly candidates.
- **Spectroscopy (step 5, below; numbers after the review fixes of 2026-10-16):**
  - **No reliable redshift:** the spectra examined did not confirm a redshift for 1316 or 1557.
  - **Peaks above the threshold:** one exists. 1316 has a boxcar S/N 6.72 peak at 6000.3 Å after the fixes. The
    pre-registered checks classify it as **not significant under the calibration**:
    - the narrow-band peak lies 2.79″ from the target (sky PA 224°, at the edge of the valid map);
    - the offset aperture 1.5″ W shows the same feature more strongly (S/N 22.9), and the one at PA 150° shows
      S/N 8.9;
    - the planned chance metric (positions at or above 6.72 in the offset apertures) gives 1.83 per window search.
  - Deviation: excluding the feature's own channels, the offset apertures have no chance run at that level. So
    the 6000 Å feature looks like real emission in the field, not a noise fluctuation, but it is centred off 1316.
  - It is not a line of 1316, and no redshift is derived from it.
  - Ordinary explanations are not excluded. Pausing further spectroscopy is a resource decision (owner,
    2026-10-16), not a scientific result.
- **Owner decision 2026-10-16:**
  - 1316 and 1557 stay as pending-vetting flags with insufficient data.
  - No stacking of single-OB cutouts, no whole-cube download, no new spectroscopy tool unit.
  - Other authorized searches continue.
  - Reopen only for a reliable target redshift, a usable 1-D spectrum, or materially better data. The arrival of
    data does not trigger an automatic download or re-run; it triggers an assessment.

## Step records

## Vetting 1316 and 2184 (Abell S1063 `arcs` anti): result (2026-10-16)

Plan: `PLAN.md` (committed before the run of `scripts/arcs_vet_noise.py`). Labels: injected sources are
simulations; re-measured axes and rates `derived`; predicted axes `model_prediction`; HST shapes `observed`.

## Measurement-noise injections (frozen plan; `noise_injections.json`, `.ecsv`)
- Gate passed: the m = 24 copies come back aligned in ≥ 0.9667 of valid injections for all 24 rows.
- **1316: P_anti = 0 / 60; 2184: P_anti = 0 / 60** → both `not noise-explicable` under the frozen rule (< 0.05).
- Sample: 3 bgsub anti rows observed, 0.017 expected from measurement noise alone → `excess remains`.
- Limitation: the injected ellipticity is each row's catalogue ellipticity (0.58 and 0.55 for the targets); a
  Gaussian profile without PSF convolution; only the bgsub (moment) measurement is calibrated, not the pipeline
  segmentation axis that defined `anti`.

## Sensitivity (post hoc; deviation from the plan)
- Injected e = 0.2 (`noise_injections_e0.2.json`): the all-row gate **fails** (min 0.45; the failing rows are
  large galaxies, `semimajor_px` 14–24 px against the 3 px weight, `.ecsv`). For the two targets the bright copies
  are aligned in 60 / 60 each, and P_anti = 2 / 60 (1316) and 1 / 60 (2184), both below 0.05 (per-target gate,
  post hoc).
- Injected e = 0 (`noise_injections_e0.0.json`): no true axis, so the gate is meaningless (P_anti 26 / 60 and
  25 / 60, near a random axis's 1/3), so this test cannot exclude a near-round true source.
- Reading: if the true sources are as elongated as observed (or e ≥ 0.2), measurement noise at their flux does not
  turn an aligned axis anti. A near-round true source remains a possible explanation; the intrinsic-shape null
  (`../arcs.json`, P = 0.0021 for the sample) is the test for that.

## Independent HST axis (post hoc; `hst_shapes.json`)
- HFF-DeepSpace SExtractor shapes (VizieR `J/MNRAS/506/928`), separate data and noise: 1316 PA 49.4°, 65.9° from
  the prediction (`anti`); 2184 PA 65.3°, 86.1° from it (`anti`). Both are at S/N ≈ 1 in that catalogue (σ_mag 1.07
  and 1.26), so this agreement is weak. ASSUMPTION: the HFF mosaics are north-up.

## Looked at (F200W 6″ raw, 1″-filtered and 0.75″ zoom; figure not stored)
- 1316: two compact knots about 5 px (0.16″) apart along the measured axis (by eye): a blend of two faint sources
  or a two-clump galaxy is a leading ordinary explanation; a deblend has not been run.
- 2184: a faint irregular blob about 2″ from the giant tangential arc, whose direction matches the prediction.

## Status
Both stay **open, unvetted, not candidates**. Tested and passed: measurement noise (observed elongation), ICL,
model error (four models). F444W axis: untestable. Untested: a blend or clumps (1316 by eye), photo-z failure (an all-VizieR
query within 1.5″ found no spectroscopic redshift for either; the nearest is the giant-arc image at z 1.2282, 1.34″
from 2184; not stored), fragment of the giant arc (2184).

## Next
Deblend 1316 into its two knots and measure each knot's colour and axis; a spectroscopic redshift for 2184 (MUSE
cube of the field, if public through the ESO archive) or a photo-z from a second code.

# Step 2 (2026-10-16): deblend 1316, arc-fragment test for 2184, core-flux noise injections

Plan: `PLAN_deblend.md` (before the run) and its addendum (after the deblend run, before the core-flux run). Script:
`scripts/arcs_vet_deblend.py` → `deblend.json`, `deblend.png`; `scripts/arcs_vet_noise.py --flux core` →
`noise_injections_core.{json,ecsv}`. Labels as above; core fluxes `derived`.

## 1316 (frozen rule)
- Two round knots beat one elliptical Gaussian: ΔBIC 134.2, knot S/N 12.8 and 10.6, separation 0.208″ →
  `two-knot`. The frozen gate passed (20 / 20 of each injected kind classified correctly).
- The separation axis is 78.5° from the z_s = 2 prediction (the one-component axis 74.6°): the anti axis is the
  knot separation.
- Knot flux ratio across F090W–F200W: χ² p = 0.146 → `consistent SED`. **Frozen outcome: clump pair; stays open.**
- **Post hoc (deviation):** the fitted F200W flux (0.52 image units, two knots) is a third of the catalogue flux the
  frozen gate injected (1.52). At the fitted fluxes and sizes the gate gives 0.75 / 0.75 (`gate_fitted_post_hoc`),
  below 0.8, so at the source's own flux the `two-knot` call is not validated (5 of 20 single injections were
  called `two-knot`). The observed ΔBIC (134) is not compared with those injections' ΔBIC (not stored).

## 2184 (frozen rule)
- Control failed: the two arc-ridge points (1.60″ and 1.93″ from 2184; `deblend.png`) are `SED-different` at
  S/N ≈ 1000 with ERR-only errors (p = 0). **Frozen outcome: test invalid.**
- **Post hoc (deviation):** with a 10 % per-band floor (the plan's own 0.1 mag PSF ASSUMPTION) the control passes
  (p 0.999) and 2184 differs from both arc points (p 8.8 × 10⁻⁶ and 2.7 × 10⁻⁶): 2184's SED peaks at F150W and falls
  to F444W, the arc's is flat (`deblend.png`). Under that post hoc rule a fragment of the z 1.2282 giant arc is
  excluded.

## Core-flux noise injections (addendum)
- The 0.3″ core magnitude minus the catalogue `mag`: 1.24 mag (1316), 2.49 mag (2184), from −0.47 to 3.09 mag over
  the 24 rows (`noise_injections_core.ecsv`, `mag` vs `mag_injected`). Step 1 injected every row at its catalogue
  flux.
- Gate passed (min 0.9667, 24 rows, none dropped).
- **1316: P_anti 0 / 60 → `not noise-explicable`. 2184: P_anti 7 / 60 (0.117) → `noise-explicable`.**
- Sample: 3 bgsub anti rows (1316, 1557, 2184), 0.333 expected; P(≥ 3) = 0.00039 → `excess remains`.
- Exploratory runs with one offset for all rows (`noise_injections_magoff1.24`, `_magoff2.49`; not in the plan,
  superseded by the per-row core run): 1316 P_anti 0 and 0.133; 2184 0 and 0.117.

## Status after step 2
- **2184: explained by measurement noise at its measured core flux** (post hoc core-flux injections). A fragment of
  the giant arc is excluded under a post hoc rule. Not a candidate.
- **1316: open, unvetted, not a candidate.** Not explained by measurement noise at core flux (0 / 60). Its axis is
  the separation of two knots with consistent colours (frozen), but the two-knot call is not validated at its own
  flux (post hoc gate 0.75). A clump pair of one source is a leading ordinary reading that the data here do not
  confirm or exclude.
- **Sample:** the anti excess on S1063 remains at core flux (1316 and 1557 at P_anti 0). 1557 (earlier reading:
  likely an intrinsic shape or blend in weak shear) has no injection test of its own.

## Next
- 1557: the same step-2 tests (it is a bgsub anti row with P_anti 0 at core flux).
- 1316: a redshift (MUSE cube through the ESO archive) or a second photo-z code; a PSF-convolved two-component fit
  whose gate passes at the source's own flux.
- Use `--flux core` as the default injected flux for the other `arcs` fields (catalogue fluxes can be up to ≈ 3 mag
  brighter than the core).

# Step 3 (2026-10-16): deblend 1557

Plan: `PLAN_1557.md` (before the run; the gate that counts is the one at the source's own fitted flux and observed
separation). Script: `scripts/arcs_vet_deblend.py --target 1557` → `deblend_1557.json`, `deblend_1557.png`. Post hoc
checks: `scripts/arcs_vet_clumps.py` → `clumps_1557.json`. Labels as above.

## 1557 (frozen rule)
- Two round knots beat one elliptical Gaussian: ΔBIC 227.4, knot S/N 15.0 and 14.1, separation 0.147″ → `two-knot`.
- Gate at the fitted fluxes and the observed separation (0.15″): 0.85 single, 0.95 pair → **passed**. The
  catalogue-flux gate (reported, not used): 0.95 / 0.90. The fitted F200W flux (0.76–0.77 image units) is close to
  the catalogue flux (0.85), unlike 1316.
- The separation axis is 75.2° from the z_s = 2 prediction (one-component axis 77.6°): the anti axis is the knot
  separation.
- Knot flux ratio across F090W–F200W: χ² p = 0.0193 (4 bands) → `consistent SED` (threshold 0.003).
- **Frozen outcome: clump pair; stays open.**

## Post hoc (not in the plan; `clumps_1557.json`)
- Looked at (`deblend_1557.png`): besides the fitted body and compact knot, a PSF-sized point (σ 0.032″) sits 0.294″
  from the body, at PA 111.7° from it, i.e. along the same axis; it is inside 1557's catalogue segment (the nearest
  other F200W catalogue source is 2.39″ away; catalogue `a750a902…`, not stored). The 0.4″ two-component fit does
  not model it.
- Three-component forced fluxes: the point is `SED-consistent` with the body (p 0.081) and with the compact knot
  (p 0.0069); body vs knot p 0.035. With ERR-only errors (the 2184 control showed these can be too small) all three
  share one colour within the frozen threshold. Reading: one clumpy (chain) source or a same-colour group whose
  catalogue axis is the chain direction; the data here do not tell which.
- **Clump-orientation model** (model_prediction; ASSUMPTIONs: isotropic point-pair separations, point clumps, κ < 1):
  a separation vector lensed by reduced shear g lands ≥ 60° from the shear axis with probability 0.144 / 0.137
  (g 0.43 / 0.45, 1557), 0.084–0.113 (g 1.6–1.9, 1316), 0.073 (g 1.5). The catalogue shape-noise null
  (`../arcs.json`) draws from the elongation distribution of all sources and does not model clump pairs, whose
  intrinsic "elongation" is maximal. Under this model each of the two clump-pair rows is anti with probability 0.08–0.14;
  how many of the 24 strong-shear rows are clump pairs is not measured, so no sample probability is computed.

## Status after step 3
- **1557: open, unvetted, not a candidate.** Frozen: clump pair. Post hoc: a three-component same-colour chain along
  the anti axis.
- **1316: unchanged** (open; clump pair, two-knot call not validated at its flux).
- **Sample:** both remaining anti rows have axes set by clump separations; the step-2 `excess remains` (P 0.00039) is
  computed against a null that does not model clumps. The ordinary reading "intrinsic clump orientation" is not
  excluded; it is not tested either until the clump-pair fraction among the aligned rows is measured.

## Next
- Run the frozen `deblend_class` (with the own-flux gate) on all 24 strong-shear background rows, and re-compute the
  anti expectation with clump-pair rows given the clump-orientation probability (a plan first: the sample test is
  then post hoc, so label it exploratory).
- Redshifts for 1316 and 1557: no catalogued z_spec (VizieR, step 1); the ESO MUSE cube of the field (size not
  checked) needs a download decision under the 200 MB rule.

# Step 4 (2026-10-16): clump-pair census of the 24 strong-shear background rows (exploratory)

Plan: `PLAN_census.md` (committed before the run, 02d9133). Script: `scripts/arcs_vet_census.py` →
`census.json`, `census.ecsv`. The sample test is post hoc (the clump idea came from looking at 1316 and 1557), so
the result is **exploratory**. Labels: pixels observed; classes derived; injections simulations; P_cp, P_sn
model_prediction.

## Census (frozen rule, own-flux gate per row)
- Gate passed for 21 / 24 rows; failed for 1404, 2004 and 2184 (the faint, low-S/N rows; `census.ecsv`
  `gate_single_ok`, `gate_two_ok`).
- `two-knot` among the 21: **5** (1061, 1316, 1452, 1557, 2449) → f̂ = 0.238, 95 % upper 0.472. Three of the five
  are screen-aligned, two are the screen-anti 1316 and 1557.
- Looked at (F200W 1.2″ 1″-filtered stamps of all 24, figure not stored): 1061 and 1557 are two compact knots;
  1452 and 2449 are knots inside or beside a larger body. The large aligned galaxies (1994, 3551, 3763, 3816) are
  clumpy chains along their axis but are `single` at the 0.4″ fit radius: the census measures pairs at ≲ 0.4″ only.

## Sample test (frozen outcome)
| Null (`census.json` "sample") | expected anti | P(≥ 4), screen count | P(≥ 3), bgsub | P(≥ 2), 2184 removed from count and null (after fix 4; [before]) |
|---|---|---|---|---|
| Shape noise only (Poisson-binomial check) | 0.521 | 0.00082 | 0.011 | 0.0885 [0.0885] |
| **Clump-aware, per row (primary)** | 1.243 | **0.026** | 0.116 | 0.354 [0.359] |
| Mixture, f̂ = 0.238 | — | 0.023 | 0.102 | 0.318 [0.323] |
| Mixture, f = 0.472 (95 % upper) | — | 0.093 | 0.259 | 0.536 [0.544] |

- **Frozen outcome: `excess remains under the clump-aware null`** (primary P 0.026 < 0.05).
- Reading (narrow): letting clump pairs point anywhere raises the expected anti count from 0.52 to 1.24 and the
  primary P from 0.00082 to 0.026. P ≥ 0.05 holds for the secondary counts (k = 3, k = 2) and for the primary at the
  upper end of f; those are not the frozen outcome.
- The k = 4 count includes 2004, explained by κ ≈ 1 (step 1), which no null here models. Its row probability in the
  table is the gate-failed mixture (0.191), not a κ ≈ 1 model.
- The screen quoted P 0.0021 with a Poisson approximation; the Poisson-binomial over the same rows gives 0.00082
  (`p_ge_k_shape_noise_only_check`). The comparison above uses the Poisson-binomial throughout.

## Limits
- Point-pair model (ASSUMPTIONs as step 3): isotropic separations, point clumps, κ < 1 formula also used for
  g > 1 rows (2004, 2167, 2215 have g_min 3–10).
- Pairs wider than 0.4″ and chains of ≥ 3 clumps are not counted as pairs; 1557's third point (step 3) is not modelled.
- f̂ comes from 21 rows of one field.

## Status after step 4
- **1316, 1557: open, unvetted, not candidates.** The clump-aware null reduces but does not remove the anti excess
  at the screen count (P 0.026, exploratory).
- **Next (per plan):** a spectroscopic redshift for 1316 and 1557. No catalogued z_spec exists (VizieR, step 1);
  the ESO MUSE cube of the field needs a download decision under the 200 MB rule (owner). Without it: a second
  photo-z code, or the same census on another `arcs` field with an anti excess, to test whether clumps explain it
  there.

# Step 5 (2026-10-16): can public spectroscopy settle 1316 and 1557? (owner request; feasibility only)

Plan: `PLAN_spectra.md` (committed before measuring, d443354). Script: `scripts/arcs_vet_spectra.py` →
`spectra.json`. Labels: cube pixels observed; limits derived; expected fluxes model_prediction.

## 1. Catalogued redshifts and 1-D spectra
- No published redshift lies within 2.4″ of either target. This was checked with an all-VizieR cone of 6″ around
  1316 and 1557, and each table with a redshift column was inspected (not stored). The step-1 queries (Chan+2020,
  1.5″ all-VizieR for 1316 and 2184) found none either.
- The nearest rows are 2.5–2.7″ away and are other objects, e.g. cluster-member photo-z rows in
  `J/MNRAS/470/95` and `J/ApJ/848/37`. Their distances agree with the owner's earlier check of the CLASH-VLT
  November 2021 AS1063 table (5.67″ and 3.47″, not re-downloaded here).
- Richard+2021 (`J/A+A/646/A83`, MUSE atlas of 12 clusters) does not include AS1063.
- No public extracted 1-D spectrum for either target was found.
- **Conclusion:** no existing redshift or spectrum can be used.

## 2. MUSE cubes in the ESO archive (ObsCore query 2026-10-16; `access_estsize`, `s_resolution` = PSF FWHM)
| Product | Program | Exposure | FWHM | Size | Covers |
|---|---|---|---|---|---|
| ADP.2017-03-28T12:46:01.331 (combined) | 095.A-0653 | 4.41 h | 1.08″ | 5.78 GB | 1557 |
| ADP.2017-03-23T15:58:03.937 (combined) | 60.A-9345 | 2.19 h | 1.45″ | 5.73 GB | 1316 (listed for 1557 too, but no data there) |
| single-OB cubes, public, with data at the target | 116.28VK.001 | 6 OBs, 3.79 h (1316); 6 OBs, 4.14 h (1557) | 0.62–1.65″ | 5.1–5.2 GB each | `muse_coverage.ecsv` |
| single-OB cubes, proprietary until 2026-11 to 2027-03 | 116.28VK.001 | 7.6 h (1316), 10.7 h (1557), nominal | — | ~5.2 GB each | — |

- Whole cubes are 25–30 times the 200 MB rule; none was downloaded.
- ESO's SODA service returns spatial cutouts: a 2.4″-radius cutout over all wavelengths is 12–17 MB. This route was
  used below.
- ObsCore footprints are not reliable at the cube edges. Every public cube listed as covering a target was checked
  with a 0.5″ SODA cutout (`muse_coverage.ecsv`). For 1316, 11 of 16 have data; for 1557, 13 of 26 have data.
  Examples without data: ADP.2025-11-20T15:47:55.061 (both targets) and the 60.A combined cube at 1557.
- Public depth with data at the target (sums from `muse_coverage.ecsv`, not stored):
  - 1316: 60.A combined 2.19 h + one more 60.A OB 0.75 h + 116.28VK 3.79 h ≈ 6.7 h;
  - 1557: 095.A combined 4.41 h + 116.28VK 4.14 h ≈ 8.6 h.

## 3. Measured line limits (`spectra.json`; aperture radius 0.5 × FWHM; values after the 2026-10-16 fixes, before-fix values in brackets)
Wording corrected 2026-10-16 (owner): every limit and equivalent width below is a **sensitivity estimate**. It
depends on the assumed line width (one 7.5 Å boxcar), on the per-channel noise (STAT rescaled by the MAD spread),
on the continuum (the F090W isophotal magnitude) and on the aperture (0.5 × FWHM, with an assumed point-source
enclosed fraction).
| Target | Cube | Window | Detection | median 5σ limit over usable channels (aperture, 7.5 Å) | Expected total flux, rest EW 30 / 100 / 300 Å |
|---|---|---|---|---|---|
| 1557 | 095.A combined | [O II], 8463–9059 Å (59 % not usable) | none (max S/N 2.31) [2.33] | 1.34 × 10⁻¹⁸ [1.33] | 0.33 / 1.1 / 3.3 × 10⁻¹⁸ |
| 1316 | 60.A combined | C III], 5478–6413 Å (26 % not usable) | peak S/N 6.72 at 6000.3 Å [4.77]; not significant under the calibration (below) | 3.37 × 10⁻¹⁸ [4.75] | 0.82 / 2.7 / 8.2 × 10⁻¹⁸ |
| 1316 | 60.A combined | C IV, 4750–5205 Å (5 % not usable) | none (max S/N 2.21) | 3.89 × 10⁻¹⁸ [5.51] | 1.25 / 4.2 / 12.5 × 10⁻¹⁸ |

- **Expected fluxes:** continuum from F090W isophotal AB 28.69 (1557) and 28.86 (1316); pipeline catalogue, sha256
  `3eecb26c…`, which matches the manifest.
- **Aperture loss (ASSUMPTION, not stored):** a point source puts about half its flux into a 0.5 × FWHM aperture
  (Gaussian PSF), so the effective total-flux limits are about twice those in the table.
- **Scale of the sensitivity** (from the table values and the ×2 aperture assumption; not stored):
  - The median limits correspond to rest EW of order 240 Å for [O II] in 1557 and of order 250 Å for C III] in
    1316 (was ~350 Å before the fixes). Typical C III] EWs are a few Å to ~20 Å.
  - A median limit does not exclude a line anywhere in its redshift window. The per-channel limit varies with sky
    residuals, and the flagged channels are not covered at all.
- **Whole-range check (after the fixes):**
  - 1557: among usable channels (valid data, outside the sky-line mask), no position reaches the threshold used
    (6-channel boxcar S/N ≥ 5 after the usable-channel MAD rescale); the maximum is 3.04.
  - 1316: the only position at or above it is the 6000.3 Å peak, classified below as not significant.
  - Not covered:
  - the positions that are not usable: 36.6 % of 4750–9350 Å for 1557 and 41.6 % for 1316 (`spectra.json`
    `whole_range.frac_not_usable`; no invalid channels, 17–20 % sky-line channels), including 59 % of the [O II]
    window of 1557;
  - wavelengths outside 4750–9350 Å;
  - lines broader than the boxcar, or fainter than the per-channel threshold;
  - flux outside the 0.5 × FWHM aperture;
  - 1557 in the 60.A cube and both targets in the public 116.28VK OBs, which were not used (`muse_coverage.ecsv`).
  So the check gives no hint of a wrong photo-z, but it cannot exclude one.
- **Deviation from the plan (implementation, made before any result was used):**
  - The planned sky-annulus flag failed: it sat on cluster light around 1316, and it missed the 5577 Å sky line
    next to 1557.
  - It was replaced by continuum subtraction (101-channel running median) and a variance-based sky-line flag.
  - The first-run output was overwritten.

## 4. Can MUSE separate the clumps?
- No. The clump separations are 0.147″ (1557) and 0.208″ (1316). The best public MUSE seeing is 0.62–0.68″ FWHM
  and the deep combined cubes have 1.08–1.45″, so the clumps are 3–10 times smaller than the PSF.
- At the present resolution and S/N, MUSE cannot reliably separate the two clumps' spectra, so these data cannot
  test whether the clumps share a redshift. (A blended spectrum can in principle contain lines at more than one
  redshift; that would need far more S/N than these targets give.)

## 5. What each outcome would mean (stated in the plan)
- **Line confirming a position in front of the cluster (z < 0.348):** the orientation test, which treats the row
  as a background source of this cluster, would not apply to it. This would remove the row from the anti count;
  it is an ordinary explanation of that row's flag. A line at about the cluster redshift would mark a member,
  which the background cut also excludes.
- **Line inside the photo-z window:** only the photo-z-failure explanation is removed. The clump-orientation and
  shape explanations stay open. It is not support for unusual lensing.
- **No line (the current result):** undetermined. It supports neither reading. These data could not measure
  whether the clumps share a redshift in any case.

## 6. Options and cost
- Options as assessed before the owner chose (b) on 2026-10-16; (a), (c) and (d) are not being pursued.
- **(a) Stack the public single-OB cutouts with valid data.** These are the 6 + 1 extra OBs for 1316 and the 6
  116.28VK OBs for 1557, on top of the combined cubes used here. Each cutout is about 15 MB, so the
  transfer is about 0.2–0.4 GB in pieces under 200 MB. The work is about 1–2 h: per-OB sky residuals, astrometric
  alignment, variance-weighted sum. The noise gain would be about √(8.6 / 4.4) ≈ 1.4 for 1557 and √(6.7 / 2.2) ≈ 1.75 for
  1316 if the noise scales with exposure time (a sensitivity estimate; it ignores sky residuals and OB-to-OB
  systematics, and seeing could change it either way). That corresponds to rest EW of order 130–170 Å for [O II]
  in 1557 and 150–200 Å for C III] in 1316 (median-limit estimates, not stored). This is still far above typical C III] (≲ 20 Å), so
  it does not help 1316.
- **(b) Wait for 116.28VK.001.** Its remaining OBs become public between 2026-11 and 2027-03, and a combined
  product or a published redshift catalogue may follow. This costs nothing.
- **(c) Whole-cube downloads (5.7–5.8 GB each).** They add nothing over the SODA cutouts for this question.
- **(d) Not spectroscopic:** a second photo-z code tests only the photo-z explanation, not the anti axis.

## Step 5 status (2026-10-16, closed by owner decision)
- **Conclusion:** the spectra examined did not confirm a redshift; ordinary explanations are not excluded; the pause
  is a resource decision.
- **Kept for a later reassessment:** product IDs and coverage (`muse_coverage.ecsv`), the cutout method (ESO SODA
  sync, `scripts/arcs_vet_spectra.py`), and the current limits (`spectra.json`).
- **Reassess, without automatic download or re-run, when one of these appears:**
  - a reliable redshift for either target;
  - a usable 1-D spectrum;
  - materially better data, e.g. the full 116.28VK.001 data set or a published redshift catalogue from it.

## Review fixes and re-run (2026-10-16; plan `PLAN_fixes.md`, committed before the re-run)
Before-fix outputs: `*_before_fixes.*`. Every changed number: `fixes_diff.json` (relative changes < 1e-6 ignored).
Data: the two MUSE SODA cutouts used before, and the same JWST 8″ stamps re-read by byte range. No other data.

**Fixes:**
1. The spectra noise rescale uses usable channels only.
2. Channels with any invalid aperture spaxel are excluded.
3. `fit_models` takes its flux start from the positive-pixel sum when the filtered sum is ≤ 0 or below 10 % of
   it, and its bound from the larger of the two (the 10 % clause was added after the final review).
4. Census counts that drop a row drop it from the null as well.

All four have unit tests (`tests/test_arcs_vet_fixes.py`).

**What changed:**
- **Spectra:**
  - Noise rescale 2.44 → 1.73 (1316) and 1.51 → 1.52 (1557). No invalid channels were found in either aperture.
  - 1316 C III] median limit 4.75 → 3.37 × 10⁻¹⁸; C IV 5.51 → 3.89 × 10⁻¹⁸. 1557 [O II] 1.33 → 1.34 × 10⁻¹⁸.
  - The 1316 peak at 6000.3 Å went from S/N 4.77 to 6.72.
- **Deblend (fix 3):**
  - Observed-fit quantities changed by ≤ 3e-5 relative (`fixes_diff.json` `max_relative_change_observed_fits`); no
    observed class changed.
  - Injection gates:
    - 1316 post hoc own-flux gate: two-knot recall 0.75 → 0.95, single side unchanged at 0.75, so it still
      fails (< 0.8). The two-knot call for 1316 is still not validated at its own flux.
    - 1557 gates: two-knot recall 0.90 → 1.00 (catalogue flux) and 0.95 → 1.00 (fitted flux); still passed.
- **Census:**
  - Classes and gate outcomes are unchanged (21 passed, 5 two-knot, f̂ 0.238). The gate still fails for 1404, 2004
    and 2184 (2004 two-knot recall 0.20 → 0.35; 2184 0.50 → 0.55), so **the flux bound was not the cause of those
    gate failures**, contrary to the earlier guess.
  - P(≥ 4) = 0.026 and P(≥ 3) = 0.116 are unchanged.
  - With 2184 removed from both the count and the null, P(≥ 2) = 0.354 (was 0.359 with 2184 kept in the null);
    mixtures 0.318 / 0.536; shape-noise only 0.0885.
- **Clumps (`clumps_1557.json`):** unchanged.

**6000 Å check for 1316 (`spectra.json` `peak_6000_1316`, `peak_6000.png`; rules from `PLAN_fixes.md`):**
| Check | Result | Pass |
|---|---|---|
| Valid and not sky-masked | yes | yes |
| S/N ≥ 5 after the fixes | 6.72 | yes |
| Narrow-band (5996.5–6002.8 Å) peak within 0.72″ of the target (NaN-aware smoothing) | 2.79″ at sky PA 223.8°, on the edge of the valid map (`narrowband_peak_at_valid_map_edge`) | no |
| No offset aperture (1.5″, 6 positions) with S/N ≥ 3 at the same channels | S/N 22.9 at sky PA 270° (W) and 8.9 at PA 150°; the one at PA 210° is not usable there | no |
| Expected chance positions ≥ 6.72 per window search < 0.05 (offset apertures, usable positions; planned metric) | 1.83 | no |

- **Deviation (reported, not used for the class):**
  - The planned chance metric counts the tested feature itself in the offset apertures that contain it.
    Counting separate runs gives 0.83, and excluding ±12 channels around the feature gives 0, i.e. no chance
    peak this strong elsewhere in the window.
  - So the 6000 Å signal is not a noise fluctuation. It is real emission in the field.
- **Class (frozen rule): not significant as a line of 1316.** The emission is stronger 1.5″ W and peaks
  2.79″ away, so the target aperture sits on the wing of nearby emission. A contribution from 1316 itself at the
  same wavelength is not excluded.
- Possible origin (not tested, not stored): earlier VizieR rows 5.3–5.6″ from 1316 list a z ≈ 0.611 galaxy, whose
  [O II] would fall near 6004 Å.
- **Evidence level:** not a line of 1316. No redshift is derived from it, and it is not evidence for unusual
  lensing.
- **Effect on the decision:** none. The decision to wait for 116.28VK.001 stands. There is still no reliable
  redshift, and the only threshold-crossing peak is explained as neighbouring emission at the level tested.

**Other changes after the final review (2026-10-16):**
- `window_result`'s `detected` is renamed `above_threshold`: a position above the threshold is not a detected
  line until the peak checks pass.
- `peak_check` runs only on cutouts with data.
- Direction values are stored as sky PAs (`sky_pa_deg`, `narrowband_peak_sky_pa_deg`).
- Census counts come from the screen and bgsub sets, and every null uses `p_tail`.
- `fixes_diff.json` uses a relative tolerance only.
- **Deviation:** the 116.28VK single-OB cube ADP.2025-11-20T15:47:55.061 was not re-requested, per the owner's
  limit of two cutouts. Its earlier `no_data` records are kept in `spectra_before_fixes.json` and
  `muse_coverage.ecsv`.

**Not changed (cleanup, no effect on numbers):** duplicated helpers (`axis_class`, the empty-position injection
loop), census thresholds hard-coded at the `arcs.json` values, and tests of two scripts living in
`tests/test_arcs_faint.py`.
