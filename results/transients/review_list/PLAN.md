# Position as a label: review-list chain for arc difference imaging (plan, written before any new count)

**Owner decision (2026-10-16).** Position becomes a label and a review priority, no longer a veto. Every other
criterion stays as it is. The arc mask of #276 is not changed. This unit has one budget of ≤ 2 h, starting at
12:30Z. It adds a `review` subcommand to `scripts/arc_chain.py` and no new tool.

## The chain (each stage computed by a function, not by eye)
1. **Detection.** An `arc_diff.py` hit, with the settings of the run (|S| ≥ 5 in ≥ `--min-bands` bands; multi-band is
   part of this stage).
2. **Bright cores.** The hit must not be `near_bright_core`. This is unchanged: the default morphology run judges only
   unflagged hits.
3. **Morphology.** `arc_diff_morph.py` verdict `pass` with default parameters.
4. **Review list.** Every hit that passes stages 1–3 is listed. Nothing is dropped for its position.

## Position labels (kept on every entry; they only order the list)
- **Fields:**
  - `in_mask`, `mask_distance_arcsec` and `arc_level_sigma` (the #276 arc mask, unchanged);
  - `nearest_image_arcsec` (the field's catalogued image list, as in #276);
  - `old_rule` = within 1.0″ of a catalogued image. That was the frozen rule of SMACS (#275) and Sunrise o120 (#274).
    Sunrise o052 had no frozen position rule, so the same 1.0″ is applied there.
- **Priority:**
  - **A** = `in_mask`;
  - **B** = not in the mask and ≤ 2.0″ from a catalogued image;
  - **C** = everything else.
- **Priority sets the order of review only.** A C entry is never called an artifact because of its position.

## Development controls (not independent validation)
SN Atalanta A and B (SMACS ERO o001 → VENUS o057) and Yan+2023 D21-W1, D21-W2, D21-S1 and D21-S2 (MACS0416 o211 →
o212). All six have been seen in earlier units, so they count as development controls only. No threshold is changed
for any of them.

## Acceptance (fixed now)
- **R1 (recovery of known cases).** All 6 development controls are on the review list.
- **R2 (injection, full chain without position).** SMACS uses the 3 stored `--inject-grid 80` runs from #276
  (`results/transients/arc_chain/runs/smacs_inj*`, seeds 0–2, 27.5 and 28.0 AB). The count is hit, unflagged and
  morph `pass`, over **all** injected stars in the box.
  - **Threshold:** efficiency ≥ 0.50 at 27.5 AB. The 28.0 AB value is reported with no threshold.
  - **Reported strata, no threshold:**
    - in the mask / out of the mask;
    - epoch-1 filtered arc level < 3σ / ≥ 3σ at the star.
  - **Already seen.** The in-mask stratum was published in #276 (0.50 / 0.31). The whole-box and out-of-mask values
    have not been computed.
  - **What the injections cover:**
    - Gaussian point sources added to epoch 2;
    - only on pixels above the epoch-1 90th percentile, 1.5″ apart (galaxies, arcs, ICL);
    - none on blank sky and none on the bright cores' PSF.
  - **Scope.** The efficiencies apply to these positions in this box only. In particular, an in-mask efficiency is
    never extended to out-of-mask positions.
  - **MACS0416.** The stored measurement (`inject_morph/summary_o211_o212.json`: 0.79 / 0.74 at 27.5 / 28.0 AB, hit,
    unflagged, morph pass) already has no position stage. It is cited, not re-run.
- **B1 (review burden).** For each stored box-pair used for acceptance:
  - SMACS ERO → VENUS;
  - Sunrise o010 → o052;
  - Sunrise o010 → o120;
  - MACS0416 Warhol box o211 → o212.

  Two counts are made per box-pair:
  - **added entries** = on the review list but outside `old_rule`; the threshold is **≤ 10 per box-pair**;
  - **total list**, reported.

  SMACS GO-4043 → DD 12774 failed its control in #275. Its list is reported but not used for acceptance.
- **Already seen (B1).** The morphology pass counts of these boxes are already in #274, #275 and #276, so B1 is not
  blind for them. It is fixed here so that the next box is judged by the same number.

## Classifying the review entries (step 4; labels, not rules)
Each entry gets one label from the existing records:
- `known real` (Atalanta, Yan+2023);
- `vetted ordinary` (an existing vetting record, cited);
- `artifact by eye` (the contact sheet is cited; this is a by-eye judgement and is labelled as one);
- `unvetted`.

An ordinary supernova or variable star is called that and nothing more. No entry is called a wormhole or exotic
candidate.

## Decision
- **R1, R2 and B1 all hold.** Within the remaining budget, run the next unsearched comparison from TASKS.md: "the rest
  of the SMACS arcs (more boxes)", same ERO o001 → VENUS o057 pair and settings as #275. Its box is fixed by a rule:
  1. lay a 12″ grid of 24″ box centres within the area covered by both epochs;
  2. choose the box holding the most Mahler+22 catalogued images that are not already in box 1;
  3. on a tie, take the box nearest the cluster centre (110.8270, −73.4545).

  Every entry on its review list is reported. Unreviewed entries are listed as unreviewed.
- **Any of R1, R2 or B1 fails, or the budget runs out.** Save the failing case and the review list, and start no new
  box.
- **In every outcome**, an empty or unreviewed list is not written as "zero candidates" or as an exclusion.
