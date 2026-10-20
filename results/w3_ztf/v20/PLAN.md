# W3 ZTF v20 — frozen batch search with the v13 method on the remaining v16-inventoried partitions

Plan fixed 2026-10-20, before any new light curve. Owner authorization: after #391, 4 h, ≤ 6 GB (calibration, search,
image check and records). The v13 W3 rule, the 10 cells, every threshold and all quality and dedupe cuts are frozen.
No new classifier, no limit, no wormhole statement. Other paused units are untouched.

## 0. Archived follow-up (no work this round)
- **1624 / 708205400009336 and 1625 / 709212400017577:** "no reliable evidence of a real dimming; some exposures
  still cannot be judged".
  - Per-exposure normal, unknown and conditional results stay as recorded in `../v19/RESULT.md`, including the known
    V1b′ σ defect and the post hoc checks. Nothing is re-run; no further photometry.
  - **Reopen only with:** new data, or a method that has already passed validation. Undetermined status alone is no
    reason to continue.
- No ZTF account is registered and no one outside is contacted.
- Old flagged diagnostics and the two sources are not re-checked.

## 1. Queue (`queue_script.py.txt` → `queue.json`; from `../v16/inventory.json`, no new inventory)
- **Candidates:** the 9 sampling-qualified partitions that v16 inventoried but did not queue. Qualification is
  unchanged: worst-cell necessary fraction ≥ 0.5, product ratio ≥ 0.8, ≥ 1,000 unused eligible groups.
- **Order:** unused eligible groups per GB, descending. The #385 "worst cell ≈ 0.59" pattern is **not** used as a
  gate or a sort key.

| pos | Npix | worst cell | unused | search if 500 cal | LC MB | shared fields (v1–v16) | fields without earlier data |
|---|---|---|---|---|---|---|---|
| 1 | 1356 | 0.563 | 1,221 | 721 | 289.2 | 713, 1711, 1712 | 669, 670 |
| 2 | 2490 | 0.569 | 1,901 | 1,401 | 471.6 | — | 760, 1759, 1760, 1796, 1797 |
| 3 | 2491 | 0.5395 | 1,828 | 1,328 | 462.3 | — | 760, 1796, 1797 |
| 4 | 8126 | 0.659 | 5,108 | 4,608 | 1,299.7 | 634, 1679 | 1631, 1678 |
| 5 | 2488 | 0.5365 | 1,860 | 1,360 | 490.3 | — | 760, 1759, 1760, 1796, 1797 |
| 6 | 2489 | 0.569 | 1,933 | 1,433 | 515.6 | — | 760, 1760, 1797 |
| 7 | 2900 | 0.6885 | 3,616 | 3,116 | 1,021.4 | 848, 1856 | 825 |
| 8 | 8124 | 0.6355 | 5,604 | 5,104 | 1,704.6 | 634, 635, 1679 | 1631 |
| 9 | 8127 | 0.6565 | 5,374 | 4,874 | 1,734.0 | 634, 635, 681, 1679 | 1678 |

- **Two kinds of unit:** "Npix" is the DR24 HATS catalogue partition, which is sky area; ZTF fields are telescope
  pointings. Exposures in shared fields are not independent of earlier samples (the same images). Those in the
  listed new fields are.
- **Totals:** 28,445 unused groups and 23,945 searchable if 500 calibrate. The light curves total 7,988.7 MB, more
  than the budget allows (section 4).

## 2. Per area (frozen v13 method, as v16; scripts copied from `../v16` with only paths, version and salt changed)
- **Sample** (`area_sample_script.py.txt`):
  - every eligible group under the unchanged `select()` cuts;
  - deduplicated against **all** earlier sample groups, v1–v16 (36,062 groups in `dedupe.json.gz`), by position
    within 3″ **and** shared objectid, i.e. by physical source;
  - groups bind all objectids and filters of a source.
- **Split:** the lowest 500 by crc32("v20:" + key) are calibration; the rest are search.
  - Reference (control) stars are the comparison range minus every sample member.
  - Search groups take no part in reference counts or validation; they are not read before a pass.
- **Validation:**
  - `inject --cells calibrated --n 200 --n-fade 1000 --edge time --seed-offset 20000 + 1000 × (position − 1)`;
  - `search --half dev --edge time`.
  - All injections are in the denominator.
  - **Pass:** every cell ≥ 0.5, every control ≤ 10 %, ≤ 5 calibration survivors (`decision_script.py.txt`).
- **On a pass:** search at once (`--half holdout`), over every frozen search group of that area.
- **On a failure:** record the reason and move to the next queued area. No seed change, no tuning, no added area,
  no deep diagnosis.

## 3. Image check of final (vetted) survivors only
- **Caps:** ≤ 10 groups, ≤ 100 cutouts, ≤ 50 MB (counted in the 6 GB). **Order:** queue position, then key.
- **Method:** `image_script.py.txt` and `mask_check_script.py.txt` as in v16.
- **Checks:** coverage, mask, saturation, blending, image depth and same-exposure comparisons, before any reading of
  a flux drop.
- **Classes:**
  - a catalogue miss is not a disappearance;
  - bad pixels or blending that prevent a measurement mean **undetermined**, with no new method this round;
  - anything beyond the caps stays pending.
- **If an image-confirmed dimming without an adequate ordinary explanation appears:** no new area is started, and
  the remaining budget goes to that event.

## 4. Stop rules (`driver_script.py.txt`, 2 workers, strictly in queue order)
- **No new area is started when:**
  - the projected download (measured objects MB + each started area's LC MB + 30 MB TAP) would pass 5,900 MB;
  - or the time is after 08:00 UTC (start ≈ 04:55 UTC + 3 h 05 min).
- Running areas finish.
- **Expected** (projection, not a result): with about 0.35 GB of objects tables, positions 1–7 fit (13,967
  searchable groups if all pass); 8124 and 8127 do not.
- **Image check capacity:** 10 groups.
- **Cloud disk:** each area's LC cache is deleted when it is no longer needed. Areas with vetted survivors keep
  theirs until the image check.
- **Records:** `timeline.jsonl` holds UTC times, seconds and projections. At least 30 min are kept for review and
  records.
