# Vetting 1316 and 2184 (Abell S1063 `arcs` anti): result (2026-10-16)

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
