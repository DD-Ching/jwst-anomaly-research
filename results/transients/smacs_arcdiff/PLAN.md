# SMACS 0723 arc difference imaging, ERO o001 → VENUS o057 — plan (written before the run)

**Why.** The Sunrise handoff (CHANGELOG 2026-10-16) names SMACS 0723 next. The SMACS two-epoch catalog grid found
one lensed-transient pair (SN Atalanta, system 5 arc; vetted on #237), but a catalog grid misses transients blended
into an arc (PLATE-23a, D-027 addendum 2026-10-15j). `scripts/arc_diff.py` has not been run on SMACS.

**Box.** 24″ × 24″ centred between Atalanta images A (110.825136, −73.452178) and B (110.823747, −73.452975), i.e.
(110.824442, −73.452577); it holds the Mahler+22 system 5 arc. Pair: jw02736-o001 (2022-06-07, ERO) → jw06882-o057
(2026-06-05, VENUS), F150W + F444W, as the catalog grid.

**Validation first (known case).** Atalanta A and B (0.22 / 0.28 μJy F150W, ≈ 25.4 AB, on the arc's light) must be
recovered as hits passing `arc_diff_morph.py`. If either is missed, the tool is not trusted on this box and the run
is reported as a tool failure, not a null. Plus `--inject 100` per-band point injections (recovery per magnitude).

**Hypothesis (to test).** Apart from Atalanta, no point-like flux change on the box's arcs beyond instrument,
processing and known populations. Ordinary explanations first: kernel residuals on bright cores, persistence,
spikes, pipeline-version differences, registration dipoles, caustic microlensing (ordinary).

**What changes the plan.** A non-Atalanta hit passing the morphology step within 1″ of a catalogued Mahler+22 image
(ASSUMPTION, post-hoc-free rule fixed here) gets `/vet-candidate`. A clean null gives a depth limit on the arc.

## Addendum 1 (2026-10-16T11:12Z, before the run): h061 stripe test
- **Hypothesis to test:** h061's 2024 GO-4043 F444W excess (+35.6 ± 3.4 nJy over 2022) is a straight stripe
  (diffraction spike) crossing the 0.15″ aperture but only part of the annulus. Ordinary alternatives if it
  survives: a real brightening of image 12.1 (it would then need a counterpart in its counter-images 12.2/12.3),
  or another processing artifact.
- **Method:** `scripts/stripe_phot.py`: stripe PA found in the 2024 F444W mosaic, the same sky-fixed profile
  subtraction applied to every covering mosaic, the same aperture/annulus as `hit_lightcurve.py`.
- **Validation first:** 40 Gaussian injections of 35 nJy on the stripe 0.8–2.0″ from h061 in the reference mosaic.
- **Frozen rule** = `classify_excess` (unit-tested): `invalid` unless the injection median recovered / alone is in
  [0.7, 1.3]; else z = (2024 corrected F444W − inverse-variance mean of the other corrected F444W epochs) /
  combined error; `explained` z < 2, `open` z ≥ 3, `undetermined` otherwise (ASSUMPTIONs). Applied to F444W only
  (the band of the excess); F090W/F115W reported, not classified.
- **What changes the plan:** `open` → `/vet-candidate h061` and a counter-image check; `explained` → h061 closed.
