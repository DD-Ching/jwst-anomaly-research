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
