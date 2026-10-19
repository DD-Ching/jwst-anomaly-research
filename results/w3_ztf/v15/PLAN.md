# W3 ZTF v15 — bounded batch search with the frozen v13 method (≤ 2 partitions, ≤ 20,000 groups)

Plan fixed 2026-10-19, before any objects table or light curve of this round is read. Owner authorization: the round
after #380 (4 h, ≤ 6 GB). The W3 model, the 10 cells, `--edge time` (EDGE_FRAC 0.25, EDGE_NSIG 2) and every
threshold are frozen. No new classifier, no limit.

## 1. Partition choice (metadata, product availability and download cost only; `../v12/obs_conditions.json`)

| Npix | worst-cell necessary fraction | product ratio | eligible groups | unused eligible | LC partition (MB) | fields shared with earlier samples |
|---|---|---|---|---|---|---|
| 8170 | 0.6615 | 1.05 | 6,238 | ≈ 738 (5,500 used in v14) | 2,207.9 | — (v14 itself) |
| 3837 | 0.6975 | 1.02 | 5,278 | 5,278 | 1,233.2 | 826, 1857 (v1–v5) |
| 183 | 0.547 | 1.01 | 2,303 | 2,303 | 270.2 | none |

- **8170 is not used.** A new 500-group calibration set would leave ≤ 238 search groups for a 2.2 GB re-read
  (≈ 9 MB per searched group, against ≈ 0.26 MB for 3837 and ≈ 0.15 MB for 183). Its unused groups stay unused.
- **Chosen: Npix 3837 and Npix 183,** both qualified in the #368 metadata check.
  - 3837 shares fields 826 and 1857 with v1–v5: **new sources, but partly shared exposures and systematics.** It is
    not an independent-area validation.
  - 183 has the lowest worst-cell margin (0.547).
- The 20,000-group target cannot be reached with two qualifying partitions. The expected search total is about 6,600
  groups, and nothing is relaxed to raise it.

## 2. Sources (`sample_script.py.txt`)
- **Cuts:** every eligible group of each partition (`select()`, unchanged; N_MAX lifted).
- **Dedupe:**
  - against all earlier calibration, development and search groups (v1–v3, v4, v6/v10, v11, v12/v13 2382, v14
    8170): no shared objectid, and none within 3″;
  - earlier comparison (reference) stars were all restricted to their own partitions (3838, 3839, 7134, 2382, 8170),
    which are disjoint sky from 3837 and 183, so no earlier reference star can be a new source.
- **Groups** bind all objectids of a source.
- **Split per partition:** the lowest 500 by crc32("v15:" + key) form the calibration group; the rest form the search
  group. **Control stars:** the comparison range minus every sample member.
- **Coverage:** refmag percentiles are reported before the read; series-epoch percentiles after the search.

## 3. Validation per partition (independent; a failure stops that partition only)
- `inject --cells calibrated --n 200 --n-fade 1000 --edge time`, with seed offsets 6000 (3837) and 8000 (183);
  `search --half dev --edge time` for the real calibration light curves.
- **Pass:** every cell ≥ 0.5 (all injections counted), every control ≤ 10 %, ≤ 5 calibration survivors.
- **Fail:** stop that partition, with no new seed, no threshold change and no replacement partition.

## 4. Search and image check
- **Search:** each passing partition is searched once (`search --half holdout --edge time`).
- **Image check of final (vetted) survivors only:**
  - **Order, fixed now:** 3837 before 183, then by key.
  - **Limits:** ≤ 10 groups, ≤ 100 cutouts, ≤ 50 MB. The `../v14/image_script.py.txt` method, plus official mask
    cutouts at the source for the run NDs when the cutout cap allows.
- **Questions per survivor:**
  1. coverage and mask at the target pixel;
  2. whether a missing record discriminates (the source's usual brightness against that image's depth and the
     same-exposure comparisons);
  3. whether the source is present and the flux drops reliably before / during / after;
  4. saturation, low-response or dead pixels, pairing, or common-mode effects.
- **Classes:** near-limit or poor images are unknown, never dimming. Ordinary explanations versus W3 come only after a
  confirmed dimming. If an event deserves independent follow-up, no new search is started and the rest of the budget
  goes to it.
- None of these image facts becomes a global veto.

## 5. #380 pending items
- 656 rebuilt NDs and the undetermined 6th group stay pending. They come from flagged runs that the frozen rule
  rejected at vetting (0 survivors in #370), so they are diagnostic only.
- They do not block this batch and are not re-downloaded.

## 6. Budget and stop
- **Downloads ≈ 1.75 GB:**
  - 3837 LC 1,233.2 MB;
  - 183 LC 270.2 MB;
  - objects 959 and 45, ≤ 0.2 GB;
  - TAP ≈ 30 MB;
  - cutouts ≤ 50 MB.
- **Each partition is read once.** If time runs short, no new stage is started and the unfinished list is saved.
- At least 30 min are kept for review and records.
