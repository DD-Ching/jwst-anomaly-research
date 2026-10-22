# W3 ZTF × ATLAS: dip check on the stored measurements (post hoc diagnostic of #462)

Owner brief 2026-10-22 after #462 (session start 20:56 UTC):
- **Budget:** 90 min wall clock, ≤ 100 MB download.
- **Inputs:** only the stored measurement tables (`../atlas40/phot/`, `../atlas40/calib.json`, `../atlas40/diag_*.json`).
- **Excluded:**
  - no new ATLAS forced-photometry task;
  - no large ZTF partition read;
  - no injection rerun;
  - no two-point rule;
  - no completion of the 23 missing stars.
- **Scope:** a post hoc diagnostic. #462's thresholds, outcomes and verdict are unchanged. Nothing here can be a
  wormhole candidate.

## 0. Scope correction of #462 (wording only; numbers and verdict unchanged)
- **27 / 63:** the best-case upper bound under #462's fixed method, the data acquired and the existing A decisions,
  if every remaining event were restored. It is not a lower bound, and not a general ceiling for ATLAS or for other
  methods.
- **0 / 79 fade passes:** holds with ATLAS on 17 of the 40 stars. It does not say what all 40 would give. Under
  #462's fixed method the 40-star value lies in 0–4 / 79.
- **"All three routes failed"** narrows to: the tested options (ZTF g edge points, the g + r joint worth check,
  ATLAS edge points) did not meet the requirement on the specified sample under the specified criterion.

## 1. De-duplication (step 1; frozen before any per-exposure detail is read)
- **Cases:** the 8 W3 simulated mis-decisions of #462 (`diag_*.json`):
  - 4 A-recovered → slow;
  - 4 A-insufficient → slow.
- **Mapping:** each triggering ATLAS point is mapped back to its stored row in `phot/<star>.csv.gz` by carrier and
  exact relative flux and σ. The key is (star, filter, exposure `Obs`).
- **Listed per unique exposure:**
  - star, filter, unit, `Obs`;
  - time (MJD as published, and HMJD as used);
  - uJy, duJy, baseline, k;
  - quality fields (err, chi/N, x, y, maj, min, apfit, mag5sig, Sky);
  - how many of the 8 mis-decisions use it;
  - how injection turned it into an intermediate point.
- **Checks on each mapped row:**
  - the low flux is in the un-injected table;
  - the units (uJy) match the baseline's;
  - the time conversion matches;
  - the row is not a missing-value code (9999) or a duplicate publication (same `Obs` twice).
- **Coverage limit:** reuse is counted only over the 49 events diagnosed in #462. Counting it over all events needs a
  replay, which the brief excludes.

## 2. Case selection and check rules (frozen after `dedupe.json`, before any comparison flux is read)
The step-1 table (`dedupe.json`, committed with this section) gives:
- 9 trigger rows on 7 stars, from 7 distinct exposures.
- One row (632216400001580, `02a58567o0511o`) is not low: rel 1.046. There the injected model itself was intermediate.
  It is listed, not checked.
- The 8 low rows fall in 4 **(unit, night) groups**, each one independent case:
  - G1 `02a` night 59790: exposure `02a59790o0824c`, on 3 stars;
  - G2 `01a` night 60480: 3 exposures on 2 stars;
  - G3 `01a` night 59632: 1 exposure, 1 star;
  - G4 `01a` night 60416: 1 exposure, 1 star.
- **Selection** (at most 3 cases):
  - Rank by the number of comparison measurements available in the stored data, counted only (no flux read):
    - (a) other stored stars with a row in the same exposure(s);
    - (b) other exposures (distinct `Obs`) of the same star in the same night;
    - (c) exposures of the same star in the other ATLAS filter within ±2 d;
    - (d) ZTF g / r epochs of the same star within ±2 d, from the local `data/cache/w3_ztf/<npix>/lc_rows.parquet` that
      #462 already used (no download, no partition read).
  - Ties go to the earlier night. The top 3 are checked; the 4th is listed as not checked.
- **Checks per selected case** (stored data and saved official pages only):
  - **C1 reading:** the row is the stored, un-injected value. Units, baseline group, time conversion, 9999 code and
    duplicate publication are already tested in step 1.
  - **C2 quality fields:** the meaning of err, chi/N, mag5sig, Sky, maj / min and apfit, from the official result
    description and FAQ saved in #461. Each field is ranked against the star's passing rows in the same (filter, unit):
    an outlier is outside the 0.5–99.5 percentile range.
  - **C3 common mode across stars:** every other stored star with a passing row in the same exposure gets its rel
    (own (filter, unit) baseline) and its deviation in calibrated σ.
  - **C4 same star, same night:** the other exposures, each band and unit separately.
  - **C5 nearby other bands:** ATLAS c / o within ±2 d; ZTF g / r within ±2 d (catalogue magnitudes, catflags
    reported; each band against its own median).
  - **C6 baseline stability:** the star's median rel in the same (filter, unit) over ±30 d, its row count, and the units
    seen in that window (camera or product changes).
- **Classification (first match):**
  1. **Confirmed reading or processing error:** C1 fails.
  2. **Specific data-quality clue, not proven:**
     - C3: at least 2 other stars low by > 5σ in the same exposure (exposure-level photometric problem); or
     - C2: a quality field outside the star's 0.5–99.5 % range in the documented failure direction (high chi/N; low
       mag5sig or Sky).
  3. **Variability supported by independent observations, worth follow-up:**
     - at least one other independent exposure of the same star (another `Obs`, another night or band) is also low
       by > 5σ, within ±2 d;
     - and C3 shows no common mode for those exposures.
  4. **Flux anomaly undetermined:** everything else.

  A single low point is never counted as real dimming. Normal neighbours alone never make a point an artefact. Only
  class 3 leads to a known-variable lookup (VSX / Gaia DR3 variability, small queries). Nothing changes #462.

## 3. Post-review check of G2 (owner brief 2026-10-22 after #470; written before G2's comparison flux is read)
- **Why:** #470's review found that the selection counts deviated from §2's wording.
  - Count (b) excluded all of a case's own exposures, not only the row's own.
  - Under §2's wording G2 totals 12, above G4's 8 (RESULT §2), so G2 should have been checked.
  - The `post_review` stage recomputes the counts with §2's wording and asserts this before checking anything.
- **Labels:**
  - This is a **post-review supplementary check**, not a blind test and not an independent validation.
  - The original selection (`selection.json`), the completed checks (`check.json`) and the deviation record stay as
    they are.
  - Output goes to `selection_plan_wording.json` and `check_g2.json`.
- **Inputs:**
  - the stored tables only, with the same baselines, quality cut and windows: night ±0.5 d, other bands ±2 d, baseline
    ±30 d;
  - the same C1–C6 code, moved into one function so that the `check` stage output stays byte-identical;
  - no download, no ATLAS task, no injection rerun, no new photometry.
  - Analysis time cap: 30 min.
- **What is checked for G2** (3 rows, 3 exposures `01a60480o0694o`, `…0722o`, `…0749o`, 2 stars):
  - C3: every other stored star with a passing row in each exposure. All stars are included, with none picked.
  - C4: the target's other exposures that night.
  - C5: the target's other ATLAS filter and ZTF within ±2 d.
  - C6: baseline stability.
  - The support test follows §2's text, which the original script had narrowed in two ways:
    - any other exposure of the same star counts, including the case's other exposures;
    - same-band exposures on other nights within ±2 d are also tested.
  - Classes and their order are §2's, unchanged.
- **Interpretation limits** (owner brief):
  - several stars low together supports a common photometric problem but does not prove a cause;
  - a single low point of the target stays undetermined, and normal neighbours do not make it an artefact;
  - only several independent exposures of the target low, with normal comparison stars, would be a variability lead
    worth follow-up. It would still need ordinary explanations checked, and it would not be a wormhole candidate.
- **No further round:** if nothing new is found, this check line is closed and no further check of this kind is
  proposed. #462's injection results and verdict stay as they are.
