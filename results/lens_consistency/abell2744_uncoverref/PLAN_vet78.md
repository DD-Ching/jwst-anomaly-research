# Vet 78.3c, test 1: does an independent model (CATS v4.1) put knot 78.3's c image where UNCOVER does? (plan)

Frozen and pushed before any CATS deflection value at system 78 is read (2026-10-23, session_01SsXPbQ).
Script: `scripts/a2744_vet78_cats.py`. Output: `vet78_cats.json` in this directory.

**Hypothesis under test (ordinary, H_pos).** The stage-2c flag 78.3c is a position error of the shared c image,
not a missing image. Reasons on record: 78.3c is a *predicted* image (UNCOVER v2.0 maps, class `candidate`); its
search circle (1.0″) lies 1.5″ from the catalogued 78.1c / 78.2c (0.13″ apart), which the same arc reference
predicts and which are `recovered`; the UNCOVER maps miss the catalogued 78.1c / 78.2c by 1.05″ / 0.99″ even when
fitted with them (`imageplane_wls_*.ecsv`); the a/b midpoints of the three knots agree within 0.05″, so the three
knots lie close together in the source plane and their c images should lie close together too (`hypothesis`).

**Alternative (H_miss).** Knot 78.3 has a c image near the UNCOVER position that is fainter than predicted.

**Measurement** (`model_prediction`; positions of 78.1a/b/c, 78.2a/b/c, 78.3a/b `observed` from the UNCOVER v2.0
list, z 2.3151 catalogue spec-z). Model `abell2744-cats` (D-035 `MapLensModel`, sha256-pinned deflection maps).
For each knot k in 78.1, 78.2, 78.3: back-trace its a and b images (only a and b, so all three knots get the same
treatment), source = `a2744_uncover.wls_source`, `lensmodel.find_images` on a 0.25″ grid (D-040 refinement). The
c-type image is the predicted image farthest from the a/b midpoint (ASSUMPTION; a fold pair plus one far image).
Report for each knot its c-type position, |μ|, and the distance to the catalogued 78.1c / 78.2c and to the
UNCOVER-predicted 78.3c (3.575198, −30.355756).

**Frozen decision (in order):**
1. *Tool check:* fails if CATS gives no c-type image for 78.1 or 78.2, or if either lies > 1.5″ (D-041 match radius,
   ASSUMPTION) from its catalogued c → "tool not validated", 78.3c stays unvetted for this test.
2. *Ordinary, position:* CATS 78.3c lies ≤ 1.0″ (the stage-2c search radius) from 78.1c or 78.2c → "explained by a
   model position error of the shared c image (H_pos supported by an independent model)".
3. *Models disagree:* CATS 78.3c lies > 1.0″ from both the catalogued c images and the UNCOVER 78.3c position →
   "untestable (D-036: an independent model disagrees)".
4. Otherwise (CATS agrees with UNCOVER within 1.0″, away from the catalogued c) → "flag survives test 1"; next is
   μ(78.3c) under ±30 % of the nearby member and group components.

**What changes the plan.** Outcome 2 closes 78.3c as ordinary; outcomes 3–4 keep it open for tests 1b / 2.
**Data.** CATS x/y deflection maps (2 × ≈ 16 MB) to the session scratchpad, deleted after the run. Cloud TLS:
`archive.stsci.edu` needs the certifi bundle (CHANGELOG 2026-10 entry).

## Test 1b (added after test 1 returned "tool not validated"; frozen before the UNCOVER maps are read)

Test 1 could not run as intended: CATS v4.1 puts no critical curve between the a and b images of any 78 knot
(`vet78_cats.json`). Test 1b asks H_pos inside the model that made the flag. Same script, `--model uncover2`
(UNCOVER v2.0 maps, `a2744_s700.load_uncover2`, ≈ 250 MB of maps to the scratchpad, deleted after the run; reason:
the flag is defined in these maps, and no smaller product holds them). Same per-knot procedure (a and b only), so
the c-type images of 78.1 and 78.2 are predictions here, not fits. Output `vet78_uncover2.json`.

**Frozen decision 1b (in order; `decide_1b`):**
1. no c-type image for any of 78.1, 78.2, 78.3 → "tool not validated";
2. the c-type images of 78.1 **and** 78.2 both lie > 1.0″ (search radius) from their own catalogued c → "explained:
   this model's c positions for system 78 miss by more than the search radius, so a non-detection inside the 78.3c
   circle is not informative";
3. else, the 78.3 c-type image lies ≤ 1.0″ from 78.1c or 78.2c → "explained: model position error of the shared c
   image";
4. else → "flag survives test 1b".
