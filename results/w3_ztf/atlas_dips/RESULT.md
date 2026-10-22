# W3 ZTF × ATLAS dip check: #462's 8 W3 simulated mis-decisions trace to 7 ATLAS exposures; the largest shows exposure-level common mode (10 / 10 stars 32–34 % low; a clue, not proven); in the 3 checked groups no reliable evidence of real variability

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
- **Variability supported by independent observations:** none in the 3 checked groups. G2 is unchecked.

## 5. Answer
- **Counts:** the 8 simulated mis-decisions come from 9 rows. 8 of the rows are low, and they come from 6 exposures on 4
  nights; the 9th row is not an observed anomaly.
- **Exposure-level clue:** the single biggest source is one exposure (G1). There, all 10 stored stars of area 4991 are
  32–34 % low. That points to an exposure-level photometric problem rather than 3 stellar dimmings; it is a clue, not
  proven.
- **New evidence:** in the 3 checked groups there is no reliable new evidence of real variability, and no specific
  source worth follow-up. G2 is unchecked.
- **G3:** stays undetermined. Deciding it needs the exposure image (an ATLAS account) or a new difference-mode task
  (excluded this round).

## 6. Cost
- About 0.65 h of the 1.5 h wall-clock budget (20:56 to the merge, ≈ 21:35 UTC), 0 MB downloaded.
- One local ZTF cache read (filtered to the group members of 4 stars).
