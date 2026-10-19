# W3 ZTF v12 — bounded batch search of new sources in areas chosen by observing conditions

Plan fixed 2026-10-19 before any objects table, exposure metadata or light curve of this round is read. Owner
authorization: the round after #365 (4 h, ≤ 3 GB, W3 / ZTF reopened for this round only).
- Reused unchanged: the #365 code (`scripts/w3_ztf.py`, `../v11/obs_script.py.txt` logic), the W3 negative-mass
  point-lens model (n = 1, ε < 0), the 10 calibrated cells, every frozen threshold, the controls and the burden rule.
- Hypothesis scope: catalogue-stage dimming runs of this one model. This is not a general wormhole search, no
  abundance limit is derived, and a null is not an exclusion of wormholes.

## 1. Dedupe list (`dedupe.json`, written before any new read)
- Sources used before as calibration or search groups: v1–v3 (`../sample.json`), v4, v6 and v11 samples. Each is
  stored with its position and member objectids.
- **A new group is not new** if any of its objectids is in that list, or it lies within 3″ (2 × 1.5″ group radius) of
  a listed group.
- Comparison (control) stars of earlier rounds were read but never searched. They are defined by rule (refmag 14–18,
  non-transient, in the area) in Npix 3838, 3839 and 7134, and the list records them by that rule. A group in those
  partitions may have served as a comparison star; such groups are counted as newly searched and labelled "read
  before as a comparison star".

## 2. Candidate areas (`candidates.json`, from the HATS partition lists only)
- **Rule:**
  - Norder-5 light-curve partitions whose order-4 parent has an objects partition and is not 959 or 1783 (the areas
    used so far);
  - 300k–700k rows;
  - centre dec +25° to +65°, |b| ≥ 30°.
  - Three parents evenly spaced in Npix order; their 2, 2 and 1 largest children (one parent had only one eligible
    child each, so 4 partitions); plus Npix 3837 (parent 959, old ZTF fields, never read) as the fallback.
- **Result:** Npix 183, 2380, 2382 and 8170, plus fallback 3837. That is 5 candidates (≤ 6).

## 3. Observing-condition and product check (metadata only; objects table, TAP, IBE HEAD)
- **Same indicator as v11** (`../v11/PLAN.md` section 2):
  - 200 positions per area from that area's seed-7 draw of eligible new groups;
  - r-band official infobits-0 exposures inside the footprint (10″);
  - each group's own exposure span, seasonal gaps kept;
  - 2,000 draws per cell, seeds 2100 + cell;
  - the noise-free model must pass the frozen `flagged_runs` over the umbra.
- **Product availability:** per position, the summed DR24 `ngoodobs` of its r-band objectids divided by the number of
  counted official exposures. This catches release gaps such as field 672 (v8–v10). The psfcat HEAD check from v11 is
  kept as information.
- **Fields:** ZTF fields of each area's positions, compared with the fields of all earlier samples.
- **An area qualifies if all hold:**
  1. its worst-cell necessary fraction is ≥ 0.5;
  2. its median product ratio is ≥ 0.8;
  3. it has ≥ 1,000 eligible new groups.
- **Choice:** up to 2 areas.
  - First the qualifying areas that share **no** ZTF field with any earlier sample, by worst-cell fraction (ties: more
    eligible groups).
  - If fewer than 2 such areas qualify, the remaining place goes to qualifying areas that share fields, including the
    fallback 3837. These are labelled as sharing observing systematics.
  - If none qualifies, the round stops before any light curve.
- The necessary fraction is used only to choose; it is never quoted as recovery.

## 4. Sample per chosen area (seed 7 draw, as before)
- Up to 3,000 eligible new groups per area (`select()`, N_MAX 3,000; cuts unchanged).
- **Split by rank of crc32("v12:" + key):**
  - the lowest min(500, n // 2) form the calibration group;
  - the rest (≤ 2,500) form the search group.
  - Across 2 areas that is ≤ 5,000 search groups.
- **Control (comparison) group:** stars of the area in the comparison range, without any v12 sample member; it is
  separate from both groups.
- Groups bind all objectids of a source (fields, quadrants, filters) within 1.5″.
- **Coverage check (reported, not a gate):** refmag and r-band epoch-count percentiles (5 / 50 / 95) of the
  calibration and search groups.
- **Separation of the search group:** it is never loaded before its area's validation passes (`source_usage.jsonl`).

## 5. Validation per area (frozen v6 rule; `decision_script`)
- `inject --cells calibrated --n 200`, plus controls and the real calibration light curves (`search --half dev`).
- **Pass if all hold:**
  - every one of the 10 cells recovers ≥ 0.5, with all injections in the denominator;
  - every control class passes vetting ≤ 10 %;
  - ≤ 5 calibration survivors.
- **On failure:** the area stops, with no retry, no new seed, no threshold change and no partial pass.

## 6. Search and image check
- **Search:** a passing area is searched once: `search --half holdout`, its frozen search group, nothing added.
- **Image check (only if there are flags):**
  - **Order, fixed now:** runs that pass vetting first (by key), then the other flagged groups by the number of
    non-detections in their longest run, across both areas.
  - **Limits:** ≤ 5 flags and ≤ 64 cutouts (≤ 12 per flag), ≤ 50 MB.
  - **Steps:**
    1. Exposure, footprint, cross-quadrant and common-mode checks first (catalogue).
    2. Then cutouts before, during and after the run, with same-exposure comparisons (`../v6/image_script.py.txt`
       method).
    3. Only a confirmed dimming is then compared with eclipse, dust, ordinary microlensing and the W3 model.
  - Anything beyond the limits is recorded as pending, not as excluded.
- A missing catalogue record is unknown: never zero flux, total disappearance or a limit.

## 7. Budget and stop
- **Downloads ≤ 3 GB:**
  - objects tables ≈ 0.1 GB per parent, 4 parents;
  - light-curve partitions ≈ 0.3–1.0 GB each (from parquet metadata before the fetch);
  - TAP;
  - cutouts ≤ 50 MB.
  - If two areas would exceed the cap, only the first is used.
- **Time:** 4 h in total, with ≥ 30 min kept for checks and records. Each file is downloaded once; caches are
  deleted after the round. There is no further round.
