# Per-pair vetting of the 11 escalated Abell 2744 orphans: result (2026-10-16)

## Plan, timing and reproducibility
- **Timing.** The main plan (`PLAN.md`) was committed at 06:33 UTC, before any measurement. The addendum (E2 controls,
  `classify`, `e2_valid`, the neighbour rule) was committed at 06:53 UTC. That was **after** the first orphan
  measurements (06:38) and before the controls, so it is **post hoc for the orphans**.
- **Code.** `scripts/orphan_vet.py` computes everything:
  - measurement: `measure_pair`;
  - controls: `--controls abell2744-cats`;
  - classes: `--classify`, which also runs the VizieR queries and the E4/SIS comparison.

  `classify` and `e2_valid` are unit-tested.
- **Files.** Every number below is in `measurements.json`, `controls.json`, `e2_validation.json`,
  `classification.json`, `vizier.json` or `pair7_uncover_dr3.json` in this folder, unless marked "not stored".
- **Exploratory.** All numbers are derived. The pairs were seen in #269.

## Deviations and fixes after review (all post hoc)
1. **E2 validation on published images** was added at the owner's request (addendum).
2. **Monte Carlo error patches** are now restricted to source-free sky (|pixel| < 4 × sky std).
3. **LW stack.** The first run's stack was mislabelled. Bands with slightly different pixel grids (F356W) were summed
   without alignment. The rerun sums a band only when both members fall on the same pixel coordinates to ≤ 0.1 px.
   - Because the grids differ, E1, E2 and E4 are measured in **F277W alone**, except for #7 and #10, where F444W is
     added. The bands are listed per pair in `lw_stack_bands`.
   - **The fix changed #3:** its E1 strip excess fell from 3.2σ to 1.2σ.
4. **E4** is now a stored conservative upper bound: flux + 3σ in a 0.2″ aperture at the brightest strip point,
   converted to AB magnitude in the code. The aperture includes the members' PSF wings, so it is an upper bound on any
   deflector light, not a non-detection limit.
5. **E4 scaling.** E4 is compared with the D-048 F160W Faber–Jackson scaling. A real galaxy at z ≈ 0.3 is brighter
   in F277W than in F160W, so the margin is conservative.

## E2 validation (11 CATS pairs ≤ 4″ apart with known z; `e2_validation.json`)
- **Failures.** The tool fails on these controls in three ways:
  - **detection:** in every control, at least one member reaches 3σ in at most 1 band (0 bands in 10 of the 11);
  - **E1:** for 7 of them the random-strip placement fails (n_random 0–3), because 2.5–4″ strips do not fit in the
    10″ box;
  - **E3:** no control has a band with both fluxes usable.

  The controls lie on bright cluster and ICL light, and several catalogue positions sit off the visible arc in F277W
  (`stamps_lw_controls.png`).
- **Rule result.** Seven controls have both members resolved: median sb_z −0.28, robust std 0.64, 1 of 7 (14 %) with
  |sb_z| > 3 (the value is −3.37). The rule allows at most 10 %, so **E2 is not valid.** With n = 7, that single
  value decides.
- **Reading.** The controls cannot validate E2 in either direction. **E2 is "undetermined" for every pair**; it
  neither supports nor excludes any pair.

## E5 and E6 public catalogues (`vizier.json`, all queries stored, empty ones included)
- **E5 spectroscopy.** No Price+25 (J/ApJ/982/51) row lies within 1.5″ of any member. Wang+24 (J/ApJS/270/12)
  `zsp` is empty at every member. E5 is unavailable for all 11 pairs.
  - For context, Price+25 has 633 rows within 4′ of the field (not stored).
- **E6 same-z neighbours.** These are Wang+24 sources within 1.5″ of the midpoint, more than 0.2″ from both members,
  with z50 within 0.15 (1 + z) of the pair's photo-z:
  - #1: 1
  - #2: 3
  - #5: 2
  - #6: 3
  - #13: 1
  - all others: 0

  Under the plan, a group never decides on its own.

## #7, member a (`pair7_uncover_dr3.json`)
- **Position.** UNCOVER DR3 has a source (id 23249) 0.006″ from our aperture centre. DJA catalogues it too (id 26610).
  The position is correct.
- **DR3 fluxes** (units of 10 nJy):
  - member a: F356W 2.45 ± 0.11, i.e. about 24.5 nJy at S/N ≈ 22 (not stored, derived from the stored values);
  - member b: 49.12 ± 0.28;
  - flux ratio about 20 (not stored).
- **Our photometry.** Our 0.15″ aperture gives 17.6 ± 6.0 nJy (2.9σ) in F356W. This is consistent with the DR3 total
  flux. Our random-aperture error is the conservative one.
- **Classification.** Under the frozen detection rule (≥ 3σ in ≥ 3 bands of our photometry), #7 counts as
  "insufficient" on detection as well as on E2.

## Classification (`classification.json`, produced by `orphan_vet.py --classify`)
| # | class | condition | E1 strip excess (σ) | E6 neighbours | E4 upper bound vs SIS-scaled (AB) |
|---|---|---|---|---|---|
| 1 | ordinary | E6 artifact flag; **visual** (star spike, #269 sheet) | 0.4 | 1 | 27.5 vs 19.0 |
| 5 | ordinary | E1 strip excess ≥ 3σ | 3.9 | 2 | 27.8 vs 19.1 |
| 6 | ordinary | E1 strip excess ≥ 3σ | 3.3 | 3 | 27.6 vs 20.0 |
| 14 | ordinary | E1 strip excess ≥ 3σ | 4.8 | 0 | 27.2 vs 19.7 |
| 2 | insufficient | E2 not valid | 1.6 | 3 | 28.5 vs 20.3 |
| 3 | insufficient | E2 not valid | 1.2 | 0 | 28.3 vs 21.0 |
| 7 | insufficient | E2 not valid; member a < 3σ in our photometry | 1.6 | 0 | 25.6 vs 18.8 |
| 9 | insufficient | E2 not valid | 2.2 | 0 | 27.1 vs 18.4 |
| 10 | insufficient | E2 not valid | 0.8 | 0 | 27.5 vs 18.5 |
| 11 | insufficient | E2 not valid | 1.2 | 0 | 27.6 vs 18.4 |
| 13 | insufficient | E2 not valid | 0.5 | 1 | 27.4 vs 19.3 |

- **Totals:** 4 ordinary (3 measured, 1 visual), 7 insufficient, 0 anomaly worth independent checking.
- **E1 strip excess.** This is the plan's "bridge" statistic: a 0.1″ strip starting 0.15″ from each member, against
  40 random strip placements in the same 10″ box.
  - An excess of ≥ 3σ is what the frozen rule counts as ordinary.
  - It can also come from overlapping PSF or profile wings of close members, which the random null does not control.
  - Light from an extended lensed source along a fold would also produce it.
  - So it is an ordinary indicator under the rule, not proof of one physical system.
- **E3 colours.** No pair has p < 0.01; the range is 0.046–0.999. The ratios are not PSF-matched, and the empirical
  errors are conservative, so E3 has little power.
- **E4.** The upper bounds are 6.8–9.2 mag fainter than the SIS-scaled magnitudes. This states only how dark a
  deflector would have to be; it is not evidence of one.

## What would decide the 7 "insufficient" pairs
- **Spectra of both members** (e.g. NIRSpec MSA): redshifts equal to < 0.01 (1 + z) and matching lines would be
  needed for the lensing reading. A ≥ 3σ difference would make the pair ordinary (E5). No public spectrum exists at any
  member.
- **Without new data:** a surface-brightness and morphology comparison that passes on published multiple images, i.e.
  PSF-matched forward modelling on aligned multi-band data. The present aperture-moment E2 fails its controls. This
  would be new tool work.

## Three separate statements
1. **Population (#269).** 51 orphans against 43.8 / 42.4 expected under nulls (e) / (f); P = 0.16 / 0.11. Null under
   the frozen rule.
2. **Individual pairs (this record).**
   - 3 pairs have a measured ordinary indicator (E1 strip excess ≥ 3σ). 1 has a visual artifact flag.
   - 7 cannot be decided with the existing data and tools.
   - None qualifies as an anomaly worth independent checking. Nothing here is a confirmed anomaly or a candidate.
3. **Completeness (#269 injections).** Per-deflector efficiencies are about 0.3–1 % at θ_E 0.3–0.7″. The screen would
   miss most dark lenses, so these checks constrain only a small part of such a population.
