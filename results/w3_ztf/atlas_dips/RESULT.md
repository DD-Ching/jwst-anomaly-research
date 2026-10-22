# W3 ZTF × ATLAS dip check: #462's 8 W3 simulated mis-decisions trace to 7 ATLAS exposures (8 low rows in 6); all 4 groups now checked (G2 after review, §7): exposure- or night-level photometric clues in G1, G2, G4, one undetermined single exposure (G3); no reliable evidence of real variability among these rows

> **Post-review update (§7):** G2, left unchecked in #470 because of a selection-count deviation, was checked after the review with the same code and rules. Sections 2–5 are kept as #470's record; §7 gives the status of all four groups.


- **Plan:** `PLAN.md`, committed in three steps:
  - the scope and de-duplication rule (b2893a8);
  - the selection, check and classification rules, before any comparison flux was read (a1bc291);
  - the selection by counts only (7448bf5).
- **Scope:** post hoc diagnostic of #462. Its thresholds, outcomes and verdict are unchanged. No new ATLAS task, no
  injection rerun, no new rule, no download (0 MB).
- **Inputs:**
  - the stored ATLAS tables (`../atlas40/phot/`);
  - `../atlas40/calib.json` and `diag_*.json`;
  - the ZTF rows in the local cache that #462 already used (`data/cache/w3_ztf/<npix>/lc_rows.parquet`);
  - the official ATLAS result description and FAQ saved in #461.
- **Labels:**
  - the 8 mis-decisions are **simulated** (injections);
  - every flux below is an **observed** ATLAS or ZTF measurement without injection;
  - nothing here is a wormhole or exotic-physics candidate.

## 0. Scope correction of #462 (wording only; numbers and verdict kept)
- **27 / 63** is the best-case upper bound under #462's fixed method, data and A decisions, if every remaining event
  were restored. It is not a lower bound, and not a general ceiling for ATLAS or for other methods.
- **0 / 79** fade passes holds with ATLAS on 17 of the 40 stars. Under the fixed method the 40-star value lies in
  0–4 / 79; it was not measured.
- **"All three routes failed"** narrows to: the tested options (ZTF g edge points, the g + r joint worth check, ATLAS
  edge points) did not meet the requirement on the specified sample under the specified criterion.

## 1. De-duplication (`dedupe.json`)
- **Counts:**
  - 8 W3 mis-decisions (4 lost, 4 insufficient → slow) have 9 triggering ATLAS rows; event 35 has one on each side.
  - The 9 rows are 9 distinct (star, filter, exposure) keys on 7 stars, from 7 distinct exposures.
  - No row is reused by two mis-decisions. Each row is the "L" point of exactly one of the 49 diagnosed events.
  - Each row was available to every flagged injection on its star: 3–11 injections (6–11 for the low rows), yet each
    decided only one. Reuse beyond the 49 diagnosed events is not counted; that needs a replay, which the brief excludes.
- **8 of the 9 rows are low in the un-injected data** (rel 0.21–0.68):
  - With the injected model at or above baseline (model dm −0.37 to −0.003), the injection rule keeps them low.
  - Each then becomes a significant intermediate point (dm ≥ 0.3 + 2σ) more than 0.25 span from the ZTF run edge, so
    the side is "slow".
- **The 9th row** (632216400001580, `02a58567o0511o`) is normal (rel 1.046). There the injected W3 model itself was
  intermediate (model dm 0.61). It is not an observed anomaly.
- **Exposures:** the 8 low rows come from 6 exposures in 4 (unit, night) groups:

  | case | exposure(s) | stars (rows) | rel | mis-decisions |
  |---|---|---|---|---|
  | G1 | `02a59790o0824c` | 3 (3) | 0.658–0.677 | 3 (#35 pre, #186, #188) |
  | G2 | `01a60480o0694o`, `…0722o`, `…0749o` | 2 (3) | 0.417–0.677 | 3 (#35 post, #161, #132) |
  | G3 | `01a59632o0033o` | 1 (1) | 0.209 | 1 (#134) |
  | G4 | `01a60416o0482o` | 1 (1) | 0.460 | 1 (#16) |

- **Reading checks, all 9 rows pass:**
  - the low flux is in the stored un-injected row;
  - uJy matches the baseline's unit;
  - no 9999 code;
  - each `Obs` is published once;
  - the `Obs` filter letter and night match F and MJD;
  - HMJD − MJD is −292 to +334 s (light travel time + 15 s).

## 2. Selection (`selection.json`; counts only, before any comparison flux)
Cases were ranked by the comparison measurements available (other stars in the same exposure, same-star same-night
exposures, other ATLAS filter ±2 d, ZTF g / r / i ±2 d; PLAN named g / r, and i was also counted, which changes no
rank):

| case | total | checked |
|---|---|---|
| G1 | 66 | yes |
| G3 | 17 | yes |
| G4 | 8 | yes |
| G2 | 6 | **not checked** |

- **G2's count:** 6 other stars (17 rows) share its 3 exposures. There are no other-filter or ZTF data within ±2 d.
- **Deviation found at review:**
  - PLAN count (b) says "other exposures of the same star in the same night".
  - The script excluded all of the case's own exposures, not just the row's own (`check_script.py.txt`, stage `select`).
  - Under PLAN's wording, G2's 3 rows add 6 same-night rows (total 12), so G2 would have ranked above G4 (8), and G2
    rather than G4 would have been checked.
  - The selection was committed before any comparison flux was read. Re-selecting now would be post hoc, so G2 stays
    unchecked. Checking it is the owner's call.

## 3. Checks and classes (`check.json`; frozen rules C1–C6, first matching class)
| case | C3: other stored stars in the same exposure | C4: same star, other exposures that night | C5: other bands ±2 d | class |
|---|---|---|---|---|
| G1 | **all 9 other stars also low**, rel 0.658–0.683 (−20 to −43 σ); 10 / 10 stars of area 4991 with data | 4 c exposures 2.2–2.6 h earlier: rel 0.98–1.01 | ATLAS o 0.9–1.9 d later: rel 0.98–1.02; ZTF g / r / i 0.2–1.2 d earlier: no unflagged epoch dimmer than its band median (−0.10 to −0.01 mag) | **2 data-quality clue** (exposure-level) |
| G3 | 8 of 9 other stars **normal** (rel 0.988–1.020); 601213300007095 **high**, rel 1.84 (+19.7 σ) | 3 o exposures 5–30 min later: rel 0.997–1.004 | no ATLAS c; ZTF r 1.1 d earlier, all flagged (catflags 32768) | **4 undetermined** |
| G4 | 2 of 4 other stars > 5 σ low (rel 0.678, 0.468); 2 at 0.849, 0.902 | 1 o exposure 16 min later also low (0.432, −5.8 σ); that exposure is common-mode too (`common_mode` true) | ZTF i 1 h earlier: +0.014 mag; ZTF g 2 d earlier: normal | **2 data-quality clue** (night-level) |

- **No checked case reaches class 3:** no low point has an independent exposure within ±2 d that is also low and free
  of common mode. So no known-variable lookup was made.
- **Script gap:** the support test reads same-night exposures (C4) and the other filter (C5). Same-band exposures on
  other nights within ±2 d were not tested by the script.
  - The reviewer's separate count found none below −2.4 σ for G1 and G4, and none at all for G3 (not stored).
  - No class changes.
- **Quality fields (C2):** none is outside its star's 0.5–99.5 % range in the documented failure direction.
  - G3 is near it: mag5sig 17.62 at the 1.3rd percentile and Sky 17.80 at the 0.6th, close to the FAQ cuts (17).
  - The FAQ says such points come from obscured or moonlit observations and "a failed astro/photometric fit because of
    clouds" is one source of nonsense values.
- **Post hoc note on G3** (seen after the checks; not used for the class):
  - The high star in the same exposure (601213300007095, chi/N 299) sits at detector x = 3005.6.
  - The low star (601209200002983) sits at x = 3009.5, 1,131 px away along y.
  - A shared column artefact (e.g. a trail raising the local sky estimate) would fit both, but it is not shown.
  - Deciding it needs the exposure image or difference-image photometry for `01a59632o0033o`.
    - The image needs an ATLAS account.
    - Difference-mode forced photometry is reachable through the shared Data Lab token (#461), but this round's brief
      excludes new forced-photometry tasks.

## 4. What is settled and what is missing
- **Confirmed reading or processing errors:** none. All 9 rows are read correctly.
  - #462's rule cannot tell an exposure-level zero-point failure from a stellar dip. That is a known method limit, not a
    reading error.
- **Specific data-quality clue, not proven:**
  - G1, 3 rows and 3 mis-decisions: one exposure, every stored star 32–34 % low;
  - G4, 1 row and 1 mis-decision: one night, several stars low.
  - Missing: the per-exposure zero point or image, which is not in the forced-photometry output used here.
- **Flux anomaly undetermined:** G3, 1 row (development star 601209200002983), 1 mis-decision.
  - It is a single exposure: the star is normal 5 min later, and 8 of 9 other stars are normal in the same exposure.
  - Missing: an image or difference flux for that exposure, or an independent same-time measurement.
- **Not checked** (cap of 3; see the selection deviation in §2): G2, 3 rows on 2 stars, 3 mis-decisions. The first
  check would be C3, on the 6 other stored stars in those exposures.
- **Variability supported by independent observations:** none in the 3 checked groups. G2 is unchecked (checked after
  review: §7).

## 5. Answer
- **Counts:** the 8 simulated mis-decisions come from 9 rows. 8 of the rows are low, and they come from 6 exposures on 4
  nights; the 9th row is not an observed anomaly.
- **Exposure-level clue:** the single biggest source is one exposure (G1). There, all 10 stored stars of area 4991 are
  32–34 % low. That points to an exposure-level photometric problem rather than 3 stellar dimmings; it is a clue, not
  proven.
- **New evidence:** in the 3 checked groups there is no reliable new evidence of real variability, and no specific
  source worth follow-up. G2 is unchecked (see §7 for the post-review check).
- **G3:** stays undetermined. Deciding it needs the exposure image (an ATLAS account) or a new difference-mode task
  (excluded this round).

## 6. Cost
- About 0.65 h of the 1.5 h wall-clock budget (20:56 to the merge, ≈ 21:35 UTC), 0 MB downloaded.
- One local ZTF cache read (filtered to the group members of 4 stars).

## 7. Post-review check of G2 (PLAN §3; `selection_plan_wording.json`, `check_g2.json`)
- **Label:** a post-review supplementary check, not a blind test or validation.
  - #470's selection, checks and deviation record are unchanged.
  - Same code (`check_script.py.txt`, stage `post_review`), baselines, quality cut and windows.
  - The `select` and `check` outputs regenerate byte-identical.
  - 0 MB downloaded; about 0.15 h of analysis.
- **Identity:**
  - With PLAN §2's wording, G2 totals 12: 6 other stars + 6 same-night rows.
  - It ranks 3rd, above G4 (8), so it is the group that should have been checked.
- **What it covers:** G2's 3 simulated mis-decisions (#35 post, #161, #132) are **3 raw rows of 2 stars in 3 exposures
  on 1 night** (unit 01a, MJD 60480). They are not 3 astronomical events.
- **Per star and exposure:** `check_g2.json` → `per_star_per_exposure` covers all 10 stored area-4991 stars in each
  exposure:
  - uJy, duJy, baseline and its (filter, unit) group, k, rel, calibrated σ, deviation and quality fields;
  - 2 stars have no row in these exposures (unknown);
  - 601213300007095 fails the quality cut in `…0694o` (excluded there by the rule).

| row (mis-decision) | target rel | other stars in the same exposure (passing rows) | target's other exposures | class (frozen rule) |
|---|---|---|---|---|
| 601213200006434 `…0694o` (#161) | 0.417 (−5.7 σ) | 5 of 6 low (rel 0.28–0.77; 2 beyond −5 σ), 1 high (1.71) | same night: 1.327 (+12.4 σ), 0.677 (−16.7 σ); next nights in o: 0.997–1.12 | **2 common-mode clue** |
| 601213200006434 `…0749o` (#35 post) | 0.677 (−16.7 σ) | 1 of 7 beyond −5 σ (0.832); **5 of 7 beyond +5 σ** (1.12–1.18; chi/N 5–37); 1 normal | same night: 0.417 (common mode), 1.327 (+12.4 σ); next nights: 0.986–1.12 | **4 undetermined** |
| 601213400002567 `…0722o` (#132) | 0.492 (−13.0 σ) | 1 of 7 beyond −5 σ (0.695), 1 at −4.8 σ (0.81); **4 of 7 beyond +5 σ** (1.20–1.33); 1 normal | same night: 0.338 (common mode), 1.122 (+9.3 σ); next night 02a: 0.69 (common mode), 0.999–1.007 | **4 undetermined** |

- **Other bands:** no ATLAS c and no ZTF epoch within ±2 d of this night (unknown).
- **Common photometric shift?** Partly.
  - `…0694o` has a common low shift under the frozen rule.
  - `…0749o` and `…0722o` do not meet the rule's "≥ 2 other stars > 5 σ low". Their comparison stars are still not
    normal: most deviate by > 5 σ, mostly **high**. Within 39 min, 601213200006434 reads 0.42 → 1.33 → 0.68 and 601213400002567 reads
    0.34 → 0.49 → 1.12.
  - This is a night-level photometric instability clue (observed; not used for the class, and not proven). The
    quality fields of the target rows are inside their stars' 0.5–99.5 % ranges.
- **Individual variability supported by independent exposures?** No.
  - 601213200006434 has two low exposures that night. The frozen rule's class-3 test is formally met for the `…0694o`
    row, by `…0749o`, but class 2 matches first.
  - The owner's condition, comparison stars normal, fails in both exposures: common mode in one, two-sided > 5 σ
    scatter in the other. Its exposure between them is 33 % high.
  - 601213400002567's next exposure is 12 % high, and its other low points are in common-mode exposures.

## 8. Status of the four groups (after §7)
| group | rows / exposures / stars | status | missing to decide |
|---|---|---|---|
| G1 | 3 / 1 / 3 | exposure-level common mode: all 10 stored stars 32–34 % low (clue, not proven) | exposure zero point or image |
| G2 | 3 / 3 / 2 | one exposure common mode (clue); two undetermined under the rule, with comparison stars scattering > 5 σ both ways (night-level clue, not proven) | the night's exposure images or zero points; no other band within ±2 d |
| G3 | 1 / 1 / 1 | undetermined single exposure (development star); 8 of 9 other stars and its own exposures 5–30 min later normal | the exposure image (ATLAS account) or a difference-mode task (excluded) |
| G4 | 1 / 1 / 1 | night-level common mode (clue, not proven) | exposure zero point or image |

- **Overall** (unchanged by §7):
  - no confirmed reading error;
  - no variability supported by independent observations among these 8 low rows;
  - no source worth follow-up from **these rows**.
- **Scope:** this covers only the ATLAS rows behind #462's mis-decisions. It says nothing about the stars' other
  epochs or other data.
- **Closure:** the check line is closed as the owner's brief requires, with no further check of this kind proposed.
  #462's injection results and verdict are unchanged.
