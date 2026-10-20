# W3 ZTF v21 — limited check of the v20 survivor and the last two inventoried partitions (frozen v13 method)

Plan fixed 2026-10-20, before any light curve of 8124 or 8127 is read and before the 8126 check is computed.
Owner authorization: after #393, 4 h, ≤ 6 GB. No new inventory, no new photometric method, no account, no external
contact. 1624 and 1625 stay archived as in v19. Other paused units are untouched.

## A. 8126 / 634210200007847: one limited check (≤ 30 min, ≤ 50 MB; `survivor_script.py.txt` → `survivor_check.json`)
- **Inputs:**
  - the committed v20 files `8126/sample.json.gz`, `search_holdout_v20_time.json`, `image_check.json` and
    `mask_check.json`: coordinates, g band, the run, the 4 NDs and the detections;
  - the DR24 light-curve API (by objectid);
  - the official psfcat of the 4 NDs and the 2 after-run detections.
- **Question:** were the ND images deep and normal enough to distinguish the rule's dimming (≥ 1 mag)?
- **Baseline:** the median g magnitude of clean (catflags 0) DR24 epochs of the survivor's own objectid, and their
  robust σ.
- **An ND exposure is "able to discriminate" only if all of these hold:**
  1. infobits 0;
  2. no mask bit 2, 8, 9 or 10 at the source (from `mask_check.json`);
  3. maglimit ≥ baseline + 0.75 mag, so an unchanged source would be well above the 5σ limit;
  4. completeness: of the isolated psfcat stars (flags 0, no other source within 15″) within ±0.25 mag of the
     baseline in an after-run detection catalogue, ≥ 90 % are found (within 1.5″) in that ND's psfcat.
- **Otherwise** the ND is "not able to discriminate": a catalogue miss is then expected even for an unchanged
  source.
- **What this check does not do:**
  - It measures **no** target flux, since there is no validated method for this source; the target flux in the NDs
    stays **unknown**.
  - A missing psfcat row is never zero flux. A limiting magnitude alone never shows that the target did not dim.
  - The result is not a global veto.

## B. Frozen queue (`queue_script.py.txt` → `queue.json`; from `../v20/remaining.json`)

| pos | Npix | unused (verified at sampling) | calibration / search | LC MB | shared fields (v1–v20) | new fields | seed offset |
|---|---|---|---|---|---|---|---|
| 1 | 8124 | 5,604 | 500 / 5,104 | 1,704.6 | 634, 635, 1631, 1679 | — | 30000 |
| 2 | 8127 | 5,374 | 500 / 4,874 | 1,734.0 | 634, 635, 681, 1678, 1679 | — | 31000 |

- **Order:** unused eligible groups per GB, as before.
- **Fields:** every field of both areas already has earlier sample data (v20 added 1631 and 1678), so **no exposure
  is independent of earlier data**. The sources are new.
- **Dedupe:** against all 53,529 earlier sample groups, v1–v20 (`dedupe.json.gz`, git-ignored, rebuilt by
  `dedupe_script.py.txt`), by position within 3″ **or** shared objectid. 0 groups were dropped for a shared objectid
  after the position cut.
- **Samples:** drawn and committed before any light curve (`8124/`, `8127/` → `sample.json.gz`, `sample_info.json`).
  - **Split:** the lowest 500 by crc32("v21:" + key) are calibration, the rest search.
  - **References:** the comparison range minus all sample members (21,159 and 16,386 stars).
- **Method:** the frozen v13 method as in v16/v20.
  - `inject --cells calibrated --n 200 --n-fade 1000 --edge time --seed-offset 30000 + 1000 × (position − 1)`;
  - `search --half dev`.
  - All injections are in the denominator.
  - **Pass:** every cell ≥ 0.5, every control ≤ 10 %, ≤ 5 calibration survivors.
  - **On a pass:** search the frozen search group at once.
  - **On a failure:** record it and leave the search group unread. No seed change, no retry, no other area.
- **Stop rules** (fixed v20 driver, 2 workers in queue order):
  - no new area once the projection (objects 97.3 MB + LC + 30 MB TAP per area) would pass 5,900 MB. The
    projection is about 3,596 MB, so both fit;
  - or after 09:10 UTC.
  - Each LC cache is deleted when no longer needed.

## C. Image check of final survivors
- **Caps:** ≤ 10 groups and ≤ 100 cutouts, in queue order and then by key; ≤ 50 MB together with part A.
- **Method:** as in v20 (`image_script.py.txt`, `mask_check_script.py.txt`). Coverage, masks, saturation,
  neighbours and local depth come before any reading of a drop.
- **Classes:**
  - a catalogue miss is not a disappearance;
  - data that are too weak to judge mean undetermined;
  - anything beyond the caps is pending, not excluded.
- **If a reliable real dimming appears:** start no new area, and use the remaining budget for its independent check
  (ordinary variables, occultation, known events, via existing catalogues).

## D. Records
- `timeline.jsonl` holds measured times.
- The last 30 min are kept for review and records.
- Zero confirmations are written "no confirmed target event"; an undetermined item is not an artefact; an unsearched
  item is not a zero; no uncalibrated limit.
