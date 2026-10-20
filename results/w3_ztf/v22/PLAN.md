# W3 ZTF v22 — new-area inventory and frozen batch search (v13 method), plus one small 8124 neighbour check

Plan fixed 2026-10-20, before the inventory reads any data and before any new image. Owner authorization: after #395,
4 h, ≤ 6 GB. No new photometric method, no account, no external contact. Other paused units are untouched.

## 0. Undetermined items (scope stated explicitly)
- **Currently tracked, from the latest batches (v20–v21):**
  - 8126 / 634210200007847: NDs fail the frozen depth and completeness criteria (v21); stays undetermined.
  - 8124 / 634205400002525: on an official low-responsivity column (v21); stays undetermined.
- **Archived (v19), not reopened:** 1624 / 708205400009336 and 1625 / 709212400017577.
- These are the W3 ZTF items only. Other units of the project have their own undetermined items; this list is not a
  project total.

## A. 8124 neighbour (≤ 15 min, ≤ 20 MB; `neighbour_script.py.txt` → `neighbour_check.json`)
Fixed check items for the 8124 calibration survivor 634205400002482, 9.9″ from 634205400002525:
1. **A different source?** Both have their own sample group. Both must be matched separately (each within 1.5″, no
   ambiguity) in the psfcat of one bracket detection.
2. **Common NDs:** list each object's run NDs, in the same band and field-quadrant, from the stored run windows.
   Exposures in a run window with no DR24 detection (LC API) are NDs. Report the shared expids.
3. **The same detector problem?** In each shared ND, the official mask bits at each source's position and the
   quadrant x distance of each from the LOW RESPONSIVITY column. The neighbour gets ≤ 12 small mask cutouts; the
   target's come from v21.
- **Report only what these items show.**
  - A shared miss is not by itself an instrument problem.
  - Visibility is not normal brightness.
  - No flux is measured, no global veto is added, and the search is not held up.

## B. Inventory (metadata only, ≤ 30 min; `candidates_script.py.txt` → `candidates.json`, `inventory_script.py.txt` → `inventory.json`)
- **Candidates:** 31 partitions under 8 new parents.
  - Rule: the v16 per-partition rule; parents never read before (v1–v21); Npix < 10000; parents with ≥ 3 eligible
    children.
  - The 8 parents are taken in farthest-point order from all used parents, so new ZTF fields are more likely.
  - Up to 4 children per parent. The cap of 50 is not reached.
  - **Scope:** this is the inventory scope of this round only, not the remaining supply of ZTF or the sky.
- **Method:** the v16 inventory, copied with only the dedupe list (v1–v21, 64,507 groups) and parallel I/O changed.
  Per partition it records:
  - unused eligible physical sources (the unchanged `select()` cuts, minus any within 3″ of or sharing an objectid
    with an earlier group);
  - the necessary-sampling fraction for the 10 cells (official r exposures);
  - the product ratio;
  - fields shared with v1–v21 samples;
  - the LC MB from parquet metadata, and searchable = unused − 500.
- **Qualification is unchanged:**
  - worst cell ≥ 0.5, product ratio ≥ 0.8, ≥ 1,000 unused;
  - the #385 0.59 pattern is not used;
  - no area that failed calibration before is retried; none is in the candidate set.
- **Queue rule** (`queue_script.py.txt`, before any light curve):
  - qualifying areas **with no shared field** come first, ordered by searchable groups per GB;
  - then areas with shared fields, by the same key;
  - at most 8 areas.
  - Seeds 40000 + 1000 × (position − 1). Split salt "v22:". About 500 calibration groups each; references are the
    comparison range minus all sample members.
  - Samples are drawn and committed before any light curve.
  - The queue, samples, costs and stop rules are committed in one commit before the driver starts.

## C. Batch (frozen v13 method; the fixed v21 driver)
- **Validation:**
  - `inject --cells calibrated --n 200 --n-fade 1000 --edge time`;
  - `search --half dev`;
  - all injections are in the denominator.
  - **Pass:** every cell ≥ 0.5, every control ≤ 10 %, ≤ 5 calibration survivors.
- **On a pass:** search the frozen search group at once.
- **On a failure:** record it and leave the search group unread.
- **Not allowed:** a seed change, a retry or a replacement area.
- **On a crash without a verdict:** the same inputs and seed may be run once more, logged.
- **Stop rules:**
  - no new area once the projection (measured inventory bytes + per-area LC MB + 30 MB TAP) would pass 5,850 MB;
  - or 2 h 50 min after the start of this round (≈ 09:40 UTC).
- **Cap:** at most 20,000 new searched groups. Nothing is relaxed to reach it.

## D. Image check of final survivors
- **Caps:** ≤ 10 groups and ≤ 100 cutouts, in queue order and then by key.
- **Before any photometry,** each entry lists:
  - the physical position;
  - the event band;
  - that band's objectid, which is **not** the group key;
  - the ND expids;
  - the same-band baseline (DR24 clean median).
- **Method:** the v20/v21 image and mask check. A catalogue miss is not a dimming; unvalidated photometry neither
  confirms nor excludes anything.
- **If a reliable real dimming appears:** stop starting new areas, and check the other band, independent exposures
  and ordinary explanations.
- The last 30 min are kept for review and records.
