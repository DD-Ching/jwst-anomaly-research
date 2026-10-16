# Full-chain known-case check for arc difference imaging: plan (written before any run of this unit, 2026-10-16)

Owner request 2026-10-16: before more arcs are searched, check that known real signals survive the whole chain. One
unit of ≤ 2 h. It uses existing tools (`arc_diff.py`, `arc_diff_morph.py`) plus one small analysis script,
`scripts/arc_chain.py`. No large tool development.

## The chain (each stage computed, never judged by eye)
1. **Detection.** An `arc_diff.py` hit: |S| ≥ 5 in ≥ `--min-bands` bands, with the settings of the original run.
2. **Morphology.** `arc_diff_morph.py --all-hits` verdict `pass`.
3. **Position.** The hit lies inside the **arc mask** (definition below).
4. **Multi-band.** Already required at stage 1 (`--min-bands 2`, or the original run's value).

A known signal "passes the full chain" only if it passes all four stages. Being detected is stage 1 only.

## Arc mask (definition fixed now; function `arc_mask` in `scripts/arc_chain.py`, unit-tested)
- **Image.** The epoch-1 image (`e1img`) of the reddest band of the run, from the `--maps` file.
- **Processing:**
  1. subtract a 3″ median filter (removes ICL and cluster-member halos);
  2. smooth with a Gaussian of σ = 0.06″;
  3. keep pixels above 3 × the robust sky σ of the filtered image (1.4826 × MAD).
- **Components.** Keep the 8-connected components that contain at least one catalogued multiple-image position
  (within 0.3″), then dilate by 0.3″.
- **Image lists:**
  - SMACS: Mahler+22 ICLv2 `arcs.dat`, as in `smacs_arcdiff_images.py`;
  - MACS0416: the CANUCS image list (`macs0416-canucs`);
  - Sunrise: Scofield+2025 (`sunrise_scofield.IMAGES`).
- **Status: post hoc.** The rule is written after #275 showed that SN Atalanta fails the 1″ radius rule. Its
  parameters are fixed here, before any mask is computed. They are not tuned to Atalanta.

## Mask validation (decided by these conditions, written now)
- **Area.** The mask must cover ≤ 25 % of the box; otherwise it does not localise arcs and is rejected for that box.
- **Positives.** Each must lie inside its mask:
  - SMACS: SN Atalanta A and B (ERO o001 → VENUS o057, F150W + F444W);
  - MACS0416: the Yan+2023 arc transients that the stored o211 → o212 run detected (D21-W1, D21-W2, D21-S1, D21-S2);
  - Sunrise: the steady arc star Earendel (position only).
- **Independence.** Atalanta motivated the rule, so it is not an independent test. The MACS0416 transients and
  Earendel were not used to design it; they are the independent checks.
- **Negative check (reported, no threshold).** Count the artifact hits, i.e. the non-`pass` verdicts, inside each mask.

## Full-chain injection efficiency (one box, SMACS pair 1)
- **Run.** `arc_diff.py --inject-grid` at 27.5 and 28.0 AB, same settings as #275, with `--maps`.
- **Count.** For injected stars that fall inside the arc mask, count those that are:
  1. a hit;
  2. hit and morph `pass`;
  3. hit, morph `pass` and in the mask.
- **Result.** The full-chain efficiency is the third count over all injected stars inside the mask.
- **Comparison.** The pre-morphology detection figures quoted earlier (Sunrise 0.97–1.00, SMACS 0.76–0.82) are
  detection only.

## Decision
- **Validated.** If every independent positive passes the full chain and the mask passes its area rule, the chain
  counts as validated on known cases, and the next new data may be searched with it.
- **Not validated.** Otherwise the failing case and its stage are recorded, and the search does not expand.
