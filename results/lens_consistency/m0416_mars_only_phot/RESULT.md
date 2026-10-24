# MACS0416: forced photometry at the MARS-only extra image (D-034 addendum as)

Plan, rule functions and tests frozen in fc8b209d before the run. Numbers from `summary.json` (frozen run),
`controls.ecsv` / `targets.ecsv` and `posthoc.json` (post hoc). Labels: CANUCS image positions and catalogue fluxes
are observed; MARS positions, |μ| and F_pred are model_prediction / derived; the 3.0″ pairing radius, the 1.5″ search
radius and the classify / local-noise thresholds are ASSUMPTIONs (PLAN.md).

## Frozen result
- Gate C2 (controls at MARS positions with MARS |μ|) **passes**: 39 / 41 testable recovered or confused (0.951);
  arc-referenced 12 / 13 (0.923). Control classes: 37 recovered, 2 confused, 2 absent, 3 undetectable,
  27 no_reference, 3 inconsistent_reference (40 families).
- One target, as in addendum ar: D115a at (64.044742, −24.060901), |μ_MARS| 29.2, F_pred 412 nJy (F200W), local σ
  4.79 nJy (F_pred / σ 86.1). Class **recovered** (4 compared sources; best F / F_pred 0.55), reading
  `seen_where_only_mars_predicts`.

## Post hoc (not frozen; `posthoc.json`)
- The two sources that make it "recovered" are the CANUCS counterparts of catalogued images: 3121427 (D115a.1's
  counterpart, 0.90″, F / F_pred 0.55) and 3121428 (D115b.1's counterpart, 0.38″, 0.35). The other two are faint
  (F / F_pred 0.009, 0.004).
- The target lies 0.85″ from observed D115a.1 and 0.21″ from observed D115b.1 (D115a / D115b are two knots of one
  z 2.2999 source), on the same arc (`targets.png`, looked at). Observed D115a.2 is 3.24″ away and is the only
  D115a image not paired to a MARS image (pairing radius 3.0″).

## Reading (narrow)
- The flux at the MARS-only position is the already-catalogued D115a.1 / D115b.1 arc, not an uncatalogued image.
  Two ordinary readings fit the stored numbers and are not separated here: the target is MARS's counterpart of
  D115a.2 displaced by 3.24″ (just beyond the pairing radius), or MARS splits the D115a.1 arc into a fold pair
  that the catalogue lists as one image. Either way it is model geometry near a critical curve (1.515″ map pixels),
  not a new image.
- The D115 two-model route stops here as frozen (PLAN.md "what changes the plan"); the D115b record stays
  `inconclusive` (its open step is the Lenstool re-fit with 8515's σ free).

## Limits
- Single band (F200W); no colour or morphology test of the arc, because the post hoc check ties the flux to catalogued
  images.
- MARS-only images of families without a CANUCS extra prediction were not searched (as addendum ar).
