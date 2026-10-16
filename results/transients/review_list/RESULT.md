# Position as a label: review-list chain, result (2026-10-16)

Plan: `PLAN.md`, committed at 12:32Z before any new count. Code: `scripts/arc_chain.py review` and `inject-review`,
unit-tested in `tests/test_arc_chain.py`. All numbers are `derived` and come from the JSON files named next to them.
The detection and morphology inputs are re-runs with `--maps` of the stored commands; their hits and morphology are
identical to the stored ones.

## Outcome against the frozen acceptance rules
| Rule | Threshold | Measured | Status |
|---|---|---|---|
| R1 known cases on the list | 6 / 6 | 6 / 6 | pass |
| R2 SMACS whole-box review-chain efficiency, 27.5 AB | ≥ 0.50 | 25 / 74 = **0.34** | **fail** |
| B1 added entries per box-pair | ≤ 10 | 4, 1, 0, 1 | pass |

R2 fails, so per the plan no new box is searched. The failing case and the review lists are saved here.

## R1: development controls (not independent validation)
All six controls are on the review list:
- **SMACS ERO → VENUS** (`smacs_o001_o057.json`): Atalanta A (h021) and B (h024), both priority B (not in the mask;
  1.69″ and 1.59″ from a catalogued image).
- **MACS0416 o211 → o212 Warhol box** (`macs0416_o211_o212_warhol.json`): D21-W1 (h061), D21-W2 (h137), D21-S1
  (h129) and D21-S2 (h138), all priority A. D21-S2 is the one entry that the old 1″ rule would have dropped.

No threshold was changed for any of them.

## R2: injections through the same chain (`smacs_inject_review.json`)
The chain counted is hit, not `near_bright_core`, morph `pass`, with no position stage. The injections are the 3
stored `--inject-grid 80` runs of #276 (seeds 0–2). These are Gaussian point sources added to epoch 2, placed only on
pixels above the epoch-1 90th percentile, 1.5″ apart; none sit on blank sky.

| Stratum | 27.5 AB: n, hit, listed (eff.) | 28.0 AB: n, hit, listed (eff.) |
|---|---|---|
| whole box | 74, 44, 25 (0.34) | 72, 26, 16 (0.22) |
| in the mask | 26, 15, 13 (0.50) | 29, 9, 9 (0.31) |
| outside the mask | 48, 29, 12 (0.25) | 43, 17, 7 (0.16) |
| arc level ≥ 3σ | 70, 41, 24 (0.34) | 66, 24, 15 (0.23) |
| arc level < 3σ | 4, 3, 1 (0.25) | 6, 2, 1 (0.17) |

- **Where recovered stars are lost** (`hit_outcomes`). Outside the mask at 27.5 AB, the 29 hits break down as:
  - 12 `pass`;
  - 8 `rosette`;
  - 5 `dipole`;
  - 3 `extended`;
  - 1 `near_bright_core`.

  Inside the mask, the 15 hits are 13 `pass` and 2 `edge`.

  Detection loses 19 of 48 outside the mask and 11 of 26 inside it.
- **Reading.** Off the arcs, morphology vetoes most of the recovered stars, mainly as rosettes and dipoles. That is
  where the review-list chain is least complete. The cause of the vetoes (for example lobes from neighbouring
  sources, as found in MACS0416, CHANGELOG 2026-10-15) was **not tested here**.
- **Scope.** These numbers hold for this box and for these injected positions. The in-mask figure does not apply
  outside the mask. Fewer than 10 stars per magnitude fall below 3σ arc level, so that stratum is too small to read.
- **MACS0416** (cited, not re-run): 0.79 / 0.74 at 27.5 / 28.0 AB on 11 tiles (`inject_morph/summary_o211_o212.json`;
  hit, unflagged, morph `pass`, no position stage). The SMACS box falls far below it. This unit does not test why.

## B1: review burden on the stored test boxes
| Box-pair | Hits | Unflagged | Listed | Added (outside 1″) | Priority A / B / C | Source |
|---|---|---|---|---|---|---|
| SMACS ERO → VENUS | 73 | 27 | 4 | 4 | 0 / 3 / 1 | `smacs_o001_o057.json` |
| Sunrise o010 → o052 | 9 | 3 | 1 | 1 | 0 / 0 / 1 | `sunrise_o010_o052.json` |
| Sunrise o010 → o120 (LW) | 12 | 1 | 0 | 0 | — | `sunrise_o010_o120.json` |
| MACS0416 Warhol box o211 → o212 | 145 | 4 | 4 | 1 | 4 / 0 / 0 | `macs0416_o211_o212_warhol.json` |
| SMACS GO-4043 → DD 12774 (failed control; not used) | 38 | 23 | 0 | 0 | — | `smacs_o4043_o12774.json` |

In total, 6 entries are added across the four acceptance pairs. These pass counts were already known from #274–#276
(see PLAN.md), so this check is not blind.

## The entries: real variability or artifact? (`labels.json`; sheets `sheet_*.png`, looked at)
| Pair | Entry | Label | Basis |
|---|---|---|---|
| SMACS | h021, h024 | known real: SN Atalanta (an ordinary lensed supernova) | published, #275 |
| SMACS | h063 | vetted, no change measured | F444W-only ± pair; forced photometry shows no ≥ 2σ change (#275) |
| SMACS | h039 | artifact by eye | F150W ± pair across a compact core; not vetted by photometry |
| Sunrise o052 | h005 | artifact by eye | F444W ± pair on a galaxy core near the map edge; not vetted by photometry |
| MACS0416 | h061, h137, h129, h138 | known real: Yan+2023 transients | published |

- **Added entries.** Of the 6 added entries, 3 are known real (Atalanta A and B, D21-S2). The other 3 are a
  no-change vetting and two by-eye artifact calls. The two by-eye calls are judgements, not rule outputs.
- **What these labels mean.** The real entries are ordinary transients. None of the entries is a wormhole or exotic
  candidate.
- **Open items.** h039 and h005 have not been vetted by photometry, so their review is not closed.

## What this does not show
- **No new data.** R2 failed, so no new observation was searched.
- **No claim about unsearched data.** Nothing here is a "zero candidates" result or an exclusion of any arc, dark lens
  or wormhole.
- **Untested ground.** The lists cover only the stored test boxes. The chain's efficiency off the arcs (0.25 at
  27.5 AB in SMACS) means that a short list there does not mean that little changed.

## Reproduce
`arc_chain.py review --field … --maps … --hits … --morph <all-hits morph> [--known … --frame-offset …] --out …`;
`arc_chain.py inject-review --field smacs --maps m m m --injected <3 injected.csv> --morph <3 morph> --out …`
(epoch 1 is the same in all three seeds, so one maps file gives the mask). The all-hits morphology of the two pairs
not stored before is in `morph/`. The maps files (8–30 MB) were deleted after use.
