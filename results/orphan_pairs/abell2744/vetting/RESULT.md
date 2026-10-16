# Per-pair vetting of the 11 escalated Abell 2744 orphans: result (2026-10-16)

The plan is `PLAN.md`: the main rules plus a 06:55 UTC addendum, both committed before the measurements they govern.
- Every number below is in a committed file of this folder: `measurements.json`, `controls.json`,
  `e2_validation.json`, `classification.json`, `vizier.json` and `pair7_uncover_dr3.json`.
- Classes are computed by `orphan_vet.classify` and E2 validity by `orphan_vet.e2_valid`; both are unit-tested.
- All numbers are derived. The pairs were seen in #269, so this vetting is exploratory.

## Deviations from the frozen plan
1. **E2 control validation was added** (owner request; addendum committed before the controls were measured).
2. **Monte Carlo error patches** are restricted to source-free sky (|pixel| < 4 × sky std) after the first run.
   Patches containing galaxies had inflated the E2 errors. This is a post-run fix of the tool (`STATUS.md`).
3. **The E4 comparison** uses the D-048 F160W Faber–Jackson scaling against a 3σ limit in the LW stack. A real galaxy
   at z ≈ 0.3 is brighter in the LW bands than in F160W, so the stated margin is conservative.

## E2 validation (11 CATS pairs ≤ 4″ apart with known z; `e2_validation.json`)
- **Detection.** No control member reaches 3σ in 3 bands with this tool. All 11 have `min(detected bands) = 0`.
  They lie on bright cluster and ICL light, and on the LW stamps several catalogue positions are off the visible arc
  (`stamps_lw_controls.png`).
- **Rule result.** Eight controls have both members resolved, with sb_z = −0.45, −3.42, −0.21, 0.02, 0.07, −0.38, −1.47
  and −0.16: median −0.29, robust std 0.36, 1 of 8 (12.5 %) with |sb_z| > 3. The rule allows at most 10 %, so
  **E2 is not valid.**
- **Reading.** The controls cannot validate the method in either direction: their members are not detected with this
  tool. **E2 is therefore "undetermined" for every pair.** It neither supports nor excludes any pair; the sb_z values in
  `measurements.json` are not used.

## E5 and E6 public catalogues (`vizier.json`)
- **E5 spectroscopy.** Price+25 (J/ApJ/982/51) has 633 rows within 4′ of the field (count not stored), but none within
  1.5″ of any member. Wang+24 (J/ApJS/270/12) `zsp` is empty at every member. E5 is unavailable for all 11 pairs.
- **E6 neighbours.** These are Wang+24 sources within 1.5″ of the midpoint, more than 0.2″ from both members, with
  z50 within 0.15 (1 + z) of the pair's mean photo-z:
  - #1: 1
  - #2: 3
  - #5: 2
  - #6: 3
  - #13: 1
  - all others: 0

  Under the plan, a group never decides on its own.

## #7, member a (`pair7_uncover_dr3.json`)
- **Position.** UNCOVER DR3 (J/ApJS/270/7) has a source 0.006″ from our aperture centre (DR3 id 23249). DJA catalogues
  it too (id 26610). Our position is correct.
- **Detection.** DR3 lists F356W = 2.45 ± 0.11 (catalogue units; S/N 22) for member a, against 49.12 ± 0.28 for
  member b. That is a flux ratio of 20.
- **Our photometry.** Our 0.15″ aperture gives 17.6 ± 6.1 nJy in F356W (2.9σ), below 3σ in every band. The nJy
  values are converted from `measurements.json` with the pixel area and are not stored as nJy. Our
  random-aperture error is the conservative one.
- **Under the frozen detection rule**, #7 is "insufficient" on detection as well as on E2.

## Classification (`classification.json`)
| # | class | condition that fired | E1 bridge z | E6 same-z neighbours | E4 3σ limit (LW, AB) vs SIS-scaled mag |
|---|---|---|---|---|---|
| 1 | ordinary | E6 artifact (star spike; #269 sheet, a visual judgement) | 0.8 | 1 | 28.0 vs 19.0 |
| 3 | ordinary | E1 bridge | 3.2 | 0 | 28.7 vs 21.0 |
| 5 | ordinary | E1 bridge | 4.7 | 2 | 28.4 vs 19.1 |
| 6 | ordinary | E1 bridge | 4.0 | 3 | 28.6 vs 20.0 |
| 14 | ordinary | E1 bridge | 6.1 | 0 | 27.8 vs 19.7 |
| 2 | insufficient | E2 not valid | 1.6 | 3 | 28.9 vs 20.3 |
| 7 | insufficient | E2 not valid; member a below 3σ in our photometry | 1.3 | 0 | 25.8 vs 18.8 |
| 9 | insufficient | E2 not valid | 1.8 | 0 | 27.7 vs 18.4 |
| 10 | insufficient | E2 not valid | 0.3 | 0 | 28.0 vs 18.5 |
| 11 | insufficient | E2 not valid | 2.9 | 0 | 28.1 vs 18.4 |
| 13 | insufficient | E2 not valid | 0.0 | 1 | 27.9 vs 19.3 |

- **Results:** 5 ordinary, 6 insufficient, 0 anomaly worth independent checking.
- **E1 bridge.** It is measured in the stacked F277W + F356W + F444W image, against 40 random strip placements in the
  same 10″ box.
- **E3 colours.** No pair has p < 0.01 (range 0.047–0.999), so no pair is ordinary on colour alone.
- **E4.** The 3σ limits are 7–10 mag fainter than the SIS-scaled magnitudes. This only states how dark a deflector
  would have to be; it is **not** evidence of one (plan, E4).

## Minimal new data or tools for the 6 "insufficient" pairs
- **Spectra of both members**, with redshifts equal to < 0.01 (1 + z) and matching line ratios. This is the deciding
  test. With a ≥ 3σ difference it would also settle a pair as ordinary (E5). Example: NIRSpec MSA prism or gratings.
  No public spectrum exists at any member.
- **Without new data:** a surface-brightness/morphology comparison that passes on published multiple images. That
  means PSF-matched forward modelling of both members that recovers known equal-surface-brightness pairs; the
  current aperture-moment E2 does not. This is new tool work, not a rerun.

## Three separate statements
1. **Population (#269):** 51 orphans against 43.8 / 42.4 expected under nulls (e) / (f); P = 0.16 / 0.11. Null under
   the frozen rule.
2. **Individual pairs (this record):**
   - 5 of the 11 have a measured ordinary indicator.
   - 6 cannot be decided with the existing data and tools, because E2 cannot be validated and no spectra exist.
   - None qualifies as an anomaly worth independent checking.
   - Nothing here is a confirmed anomaly or a candidate.
3. **Completeness:** per-deflector efficiencies are about 0.3–1 % at θ_E 0.3–0.7″ (#269, as in the deep fields). The
   screen would miss most dark lenses even if present. The population null and these individual checks constrain
   only a small part of such a population.
