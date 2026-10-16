# Full-chain known-case check for arc difference imaging: result (2026-10-16)

Plan: `PLAN.md` (committed before any run). Code: `scripts/arc_chain.py` (mask function unit-tested in
`tests/test_arc_chain.py`). Every number below is `derived` and comes from the JSON files named next to it.

## Answer in one paragraph
Of the known real transients tested, the four MACS0416 Yan+2023 transients that the detection step finds pass
every stage (detection, morphology, arc mask; multi-band is part of detection). SN Atalanta in SMACS 0723 does
**not** pass. It is detected and passes morphology, but fails the position stage under both the frozen 1″ rule
(#275) and the new arc mask. The frozen decision rule of `PLAN.md` counts the chain as validated on its independent
positives. Those are 4 transients in one MACS0416 pair, 3 of them non-trivial tests of the mask. The validation
covers only transients on bright, catalogued arc light. Atalanta-like transients on faint host light outside the
dilated, seeded arc components are still missed (see "What the chain still misses").

## Detection vs full chain (owner item 1)
| Box | Detection only | Full chain (hit, morph `pass`, in mask) | Source |
|---|---|---|---|
| Sunrise o010 → o120 LW | 1.00 / 0.97 at 27.5 AB (F277W / F356W, per-band injections, before morphology) | not measured | `results/transients/sunrise_arcdiff/o010_o120/summary.json` (#274) |
| SMACS ERO → VENUS, whole box | 0.82 / 0.76 at 27.5 AB (F150W / F444W, per-band injections, before morphology) | — | `results/transients/smacs_arcdiff/o001_o057/summary.json` |
| SMACS, injected stars inside the arc mask | 15 / 26 = 0.58 (27.5 AB); 9 / 29 = 0.31 (28.0 AB) | **13 / 26 = 0.50** (27.5); **9 / 29 = 0.31** (28.0) | `smacs_inject.json` |
| MACS0416 o211 → o212, 11 tiles (no mask) | 152 / 152 (27.5), 149 / 150 (28.0) | 120 / 152 = 0.79 (27.5), 111 / 150 = 0.74 (28.0): hit, unflagged, morph `pass` (tiles with maps kept; `pass_r0.3`, the default rosette lobe ratio) | `results/transients/macs0416_arcdiff/inject_morph/summary_o211_o212.json` |

- The Sunrise 0.97–1.00 is a **pre-morphology** number from per-band injections. Sunrise was not re-run (owner).
- The SMACS full-chain numbers come from 3 `--inject-grid 80` runs (seeds 0, 1, 2; `runs/smacs_inj*`), pooled. They
  count only stars that land inside the mask, so they are small samples. On this arc, most of the losses happen at
  detection (11 of 26 missed at 27.5 AB) rather than at morphology (2 of 15 hits vetoed). No in-mask injected hit was
  flagged `near_bright_core` (`chain_unflagged` = `chain`).

## SMACS: where Atalanta stops (owner item 2)
`smacs_known.json` (rerun of #275 with `--maps`; same settings; 73 hits, 5 morph passes, as in #275):

| Signal | Hit | Morph | Nearest catalogued image | Distance to mask | Arc level at the position | Full chain |
|---|---|---|---|---|---|---|
| Atalanta A | h021 | pass | 1.66″ | 0.86″ outside | 1.78 σ | no |
| Atalanta B | h024 | pass | 1.59″ | 0.76″ outside | 0.30 σ | no |

- **The SMACS mask fails the positive check in `PLAN.md`** ("Positives": Atalanta A and B inside the mask). Atalanta
  is not an independent case, so the PLAN "Decision" rule does not count this failure.
- **Failing stage: position.** Atalanta is stopped there under the frozen 1″ rule and also under the post-hoc mask.
- **Why the mask misses it** (`mask_smacs.png`): both supernova images lie about 0.8″ from the bright F444W ridge of
  the system 5 arc. The filtered epoch-1 level at their positions is below the mask's 3σ threshold (1.78 σ and
  0.30 σ). In this image, the host light under Atalanta is not part of any arc component that contains a catalogued
  image.
- **No retuning.** The parameters were not changed after this result. Widening the seed radius, the dilation or the
  threshold until Atalanta is inside would be fitting to the one case that motivated the rule.

## Arc mask: definition and checks (owner item 3)
- **Definition** (`PLAN.md`; post hoc with respect to #275, fixed before any mask was computed):
  1. take the epoch-1 image of the reddest band;
  2. subtract a 3″ median filter;
  3. smooth with a Gaussian of σ = 0.06″;
  4. keep pixels above 3 robust σ;
  5. keep the 8-connected components within 0.3″ of a catalogued image;
  6. dilate by 0.3″.
- **Area.** The mask covers 0.120 of the SMACS box, 0.067 of the MACS0416 box and 0.029 of the Sunrise box. All are
  within the 0.25 limit.
- **Independent positives: MACS0416 Yan+2023** (`macs0416_known.json`, `mask_m0416.png`).
  - All 11 published positions are inside the mask.
  - The 4 that this pair detects (D21-W1, D21-W2, D21-S1, D21-S2) pass the full chain.
  - D21-S2 is 1.43″ from the nearest catalogued image, so the frozen 1″ rule would have rejected it; the mask keeps
    it. D23-S3 (1.80″) is also inside the mask (position only; it is not detected in this pair).
  - D21-W1 is 0.19″ from a catalogued image, inside the seed radius, so it is a weak test of the mask.
- **Earendel is not an independent check (plan error).** Earendel is catalogued image 1.1e in
  `sunrise_scofield.IMAGES`, so it seeds its own component. Being inside the mask (`sunrise_known.json`) is
  automatic. `PLAN.md` wrongly listed it as independent. The Sunrise mask follows the arc along all its catalogued
  images (`mask_sunrise.png`), but that was judged by eye.
- **Artifact hits inside the mask** (reported; no threshold):
  - MACS0416: 27 hits are in the mask. 23 are `rosette` and `near_bright_core`; the other 4 are the 4 Yan+2023
    passes. No artifact hit passes.
  - SMACS: one hit is in the mask, h061, verdict `edge`.
  - Sunrise: no hits are in the mask.
- **Visible weaknesses of the mask:**
  - when a catalogued image sits on a cluster-member halo, the whole halo enters the mask (SMACS upper left; MACS0416
    upper right);
  - faint host light below 3σ never enters the mask.

## Status of each statement
- **Validated on independent known cases:** detection + morphology + mask for transients on bright, catalogued arc
  light. The cases are 4 / 4 in one MACS0416 pair; D21-W1 is 0.19″ from a catalogued image, so only 3 of them test
  the mask non-trivially.
- **Deviation from the PLAN decision (stricter):** the PLAN says that a validated chain may be used on the next new
  data. Because a known real transient (Atalanta) still fails, this unit does not start new data and leaves that
  call to the owner.
- **Not independently validated:**
  - the mask for transients offset from the bright arc ridge. The only such case, Atalanta, fails.
  - full-chain efficiency in boxes other than SMACS (the MACS0416 figure has no mask stage).
- **Post hoc:** the mask rule itself (motivated by #275), the `chain_unflagged` count and the multi-seed pooling of
  the SMACS injections. These were added during this unit and were not changed after the numbers were seen.

## What the chain still misses
- **Transients on host light outside the dilated, seeded components** (Atalanta type). They are missed by
  construction. Faint host light within 0.3″ of a seeded component is kept.
- **Faint transients on arc light:** about half at 27.5 AB and two thirds at 28.0 AB in the SMACS mask (small samples,
  above). Most are lost at detection.
- **Transients near bright cores.** They are flagged `near_bright_core` and are not judged by the default
  morphology run.
- **Bands not run:** Sunrise SW (owner: not pursued).

A null in any box is a statement about detectable point-like flux changes on its masked arc light at the
efficiencies above. It does not exclude the arc as a whole, dark lenses or wormholes.

## Reproduce
- **Inputs.** The detection runs are the stored commands with `--maps`:
  - SMACS: `runs/smacs/summary.json` `command`;
  - MACS0416: `runs/m0416/summary.json`;
  - injections: `runs/smacs_inj*/summary.json` (`--inject-grid 80 --inject-grid-mags 27.5 28.0`, seeds 0–2);
  - Sunrise: the #274 run (`results/transients/sunrise_arcdiff/o010_o120/hits.ecsv`, its maps; morphology
    re-run with `--all-hits` in `runs/sunrise_o120/morph.json`).
- **Morphology.** `arc_diff_morph.py --hits … --maps … --all-hits --out runs/<run>/morph.json`.
- **Known check.**

  ```
  arc_chain.py known --field {smacs|macs0416|sunrise} --maps … --hits … --morph … --known … [--frame-offset 0.196 -0.055] --out <field>_known.json --figure mask_<field>.png
  ```

- **Injection check.** `arc_chain.py inject --field smacs --maps m0 m1 m2 --injected … --morph … --out smacs_inject.json`.
- **Not committed.** The maps files (`*.npz`, 8–30 MB each) were deleted after use.
