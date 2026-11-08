# `arcs`: model shear error at the single-component anti rows (plan, 2026-11-08)

Frozen before any realisation map or new model was read for the tested rows (only the RELICS directory listings
and file sizes were looked at). Script: `scripts/arcs_shear_error.py`; tests: `tests/test_arcs_shear_error.py`
(rule functions written and tested before the run). Stored `arcs.ecsv` rows of the deciding runs plus published
RELICS maps (streamed into the cache, sha256 recorded, deleted after use).

**Question.** Addendum 8 found the 3 single-component anti rows (WHL0137 2194, RXC J0600 581, MACS0308 2157) anti
under a second published model. Is that robust to (a) the models' own statistical error (MCMC realisation maps)
and (b) further independent models?

**Hypothesis H_err (ordinary).** The anti class of these rows comes from model error in the predicted stretch
direction or |g|. Prediction: in realisation maps, or in further models of the same field, the rows often lose the
anti class.

**Inputs.**
- RXC J0600 581: RELICS Lenstool v1 `range/` (100 realisations, `lenstool-mapNNN`, x/y arcsec deflections, 167 kB
  each). This is the second model of addendum 8 (deciding model GLAFIC has no realisations).
- WHL0137 2194: two further independent RELICS v1 models, Zitrin LTM-Gauss and WSLAP (best fits, arcsec
  deflections), beyond Lenstool (deciding) and GLAFIC (addendum 8).
- MACS0308 2157: **not tested here** (deviation from the TASKS wording, stated now): its only realisations are
  Zitrin LTM-Gauss `range/` maps of 36 MB each in pixel units (≈ 7.2 GB for 100 x/y pairs); the next route is a
  byte-range read of the rows around the source. Reported `untested`.

**Re-classification.** `lens_consistency.orientation_table` unchanged on the deciding run's stored rows (observed
axis, photo-z range), default thresholds (aligned ≤ 30°, anti ≥ 60°, z margin 0.1). Maps are taken as scaled to
D_LS/D_S = 1 (RELICS readme v3, ASSUMPTION for WSLAP / Zitrin as for the other RELICS maps).

**Rules (functions in the script, unit-tested before the run).**
1. Map gate per model (and for the range set: on map000, map049, map099): κ from our deflection-map model against the
   published κ map, `lens_consistency.map_check`, median relative difference ≤ 0.05 (ASSUMPTION). Failing → that
   model / set is `invalid`.
2. Control gate per model: the deciding run's `aligned` strong-shear background rows of the field, re-classed;
   retention r = share still `aligned`; passes with r ≥ 0.70 (ASSUMPTION, the `arcs` control threshold). For the
   range set, r is the median over realisations.
3. Per row: WHL0137 2194 → `robust` if every valid, gate-passing further model gives `anti`; `model_sensitive`
   if at most half of them do; else `partly`; no valid gate-passing model → `undetermined`.
   RXC J0600 581 → f = share of realisations in which it is `anti`: `robust` f ≥ 0.68; `model_sensitive`
   f ≤ 0.32; else `partly`; failed gates → `undetermined` (thresholds ASSUMPTIONs).
4. Outcome over the tested rows: every tested row `robust` → `model error does not remove the tested rows`; any
   `model_sensitive` → `model error explains <rows>`; any `undetermined` and no `model_sensitive` → `undetermined`;
   otherwise `mixed`.
Secondary (reported, not rules): 16 / 50 / 84 % of the realisations' predicted axis offset and g_min at 581; per
model offset and g_min at 2194.

**What would change the plan.** `model error explains`: the row leaves the open set and the pooled excess shrinks
further. `robust`: the remaining ordinary routes are a spectrum (membership) and MACS0308 2157 via byte ranges.
