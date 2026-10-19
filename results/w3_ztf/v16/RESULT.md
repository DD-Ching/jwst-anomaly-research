# W3 ZTF v16 — inventory of 25 partitions, frozen queue of 8: 3 pass, 5,995 groups searched, 2 survivors, both data problems

- **Commits:**
  - `PLAN.md` and the scripts were committed before the inventory results;
  - `queue.json` (5e80009) before any light curve;
  - `timeline.jsonl` holds every stage's UTC times.
- **Labels:** DR24, TAP, IBE images and masks `observed`; counts `derived`; injections `simulated`.
- **Scope:**
  - the frozen v13 W3 rule, the 10 calibrated cells, catalogue stage only;
  - no limit derived; nothing here excludes wormholes.
  - The inventory covers the 25 listed partitions only; it is not an all-sky census.

## 1. Inventory (`inventory.json`, `queue.json`; metadata only)
- **Candidates:** 25 partitions under 8 parents never used before (rule in `candidates.json`).
- **Inventory time:** 48.7 min, longer than the planned 30 min, mostly TAP for 246,499 exposure rows and the
  sampling computation.
- **Downloads:** objects 666.1 MB and TAP 47.6 MB.
- **Three counts, kept apart:**

| count | partitions | unused groups | searchable after 500 calibration each |
|---|---|---|---|
| catalogue-eligible (unchanged cuts, deduplicated against 20,153 earlier groups) | 25 | 54,080 | — |
| sampling-qualified (worst cell ≥ 0.5, product ratio ≥ 0.8, ≥ 1,000 groups) | 17 | 44,354 | 35,854 |
| queued (top 8 by unused groups per GB) | 8 | — | 11,909 |
| **calibration passed** | **3** | — | **5,995** |

## 2. Queue results (`summary.json`; each area's `decision.json`)

| position | Npix | worst-cell fraction | verdict | cells < 100/200 | fade passing / 1,000 | searched |
|---|---|---|---|---|---|---|
| 1 | 187 | 0.578 | failed | 4 (three t_E 10 d cells, one t_E 100 d) | 97 | — |
| 2 | 1626 | 0.551 | failed | 2 (ρ 0.1, t_E 10 d) | 73 | — |
| 3 | 1627 | 0.593 | **passed** | 0 | 70 | 2,064 groups, 9 flagged runs, 0 survivors |
| 4 | 2274 | 0.508 | failed | 3 (t_E 10 d) | 78 | — |
| 5 | 1624 | 0.6385 | **passed** | 0 | 89 | 1,887 groups, 16 runs, **1 survivor** |
| 6 | 1358 | 0.5105 | failed | 2 (ρ 0.1, t_E 10 d) | 39 | — |
| 7 | 1625 | 0.6615 | **passed** | 0 | 64 | 2,044 groups, 5 runs, **1 survivor** |
| 8 | 1357 | 0.517 | failed | 2 (t_E 10 d) | 61 | — |

- **All failures are in the cell condition.** Every control (eclipse, PSPL, common, fade ≤ 9.7 %) and every
  calibration burden is within its limit. The failing cells are mostly t_E 10 d.
- **Post hoc pattern, not a calibrated rule** (`remaining.json`): in v15–v16, areas with a worst-cell fraction ≥ 0.593
  passed (4) and those ≤ 0.578 failed (6).
- Failed areas' search groups were never loaded (their `source_usage.jsonl`). No retries; nothing was added to the
  queue.

## 3. Final survivors and image check (2 groups, 19 science cutouts and 11 mask cutouts, 1.75 MB)
- **1624 / 708205400009336** (r band, 9 epochs over 2 nights, MJD 59815.5–59817.5):
  - The source is bright: about 2,450 µJy, r ≈ 15.4 AB. It is **visible in all 8 checked NDs**; one more ND is pending
    (cap). The image limits are 19.1–19.8.
  - The official mask flags **bit 2 LOW RESPONSIVITY** within the 2″ aperture in all 8 NDs (centre value 4 in 5 of
    them).
  - A dark column lies 0.0–1.75″ from the source in the NDs; in the 4 detections it lies ≥ 2.2″ away.
  - **Reading:** consistent with a low-response detector column under the source (a data problem). Not a dimming.
- **1625 / 709212400017577** (g band, 3 epochs, MJD 59896.4–59900.5):
  - All 3 NDs are poor images (seeing 4.05–4.64″; class 2), and the source is visible in each.
  - The mask carries bit 1 (SExtractor detection) at the source.
  - 2″-aperture flux is 101–147 µJy against a 154 µJy reference. That is consistent with aperture loss in wide seeing
    (not separated).
  - **Reading:** a catalogue miss in poor images. Not a dimming.
- **Result:** 2 final survivors, 0 image-confirmed dimmings, 0 undetermined among the survivors. One ND exposure of
  708205400009336 is pending.
- Flagged runs removed by vetting (30 in these 3 areas) are **not** shown by image to be artefacts; they were not
  image-checked.
- Older pending items (#380 / #383) are unchanged.

## 4. Speed and volume (measured; `summary.json`, `timeline.jsonl`)
- **Queue execution** (8 areas, 2 workers): 23:25–23:39 UTC, 14 min.
- **Stage seconds per searched area** (sample, fetch, exposures, inject and dev, search, image): 1,142 s in total for
  the 3 searched areas, i.e. 190 s per 1,000 searched groups (concurrent stages are summed).
- **The whole round,** including the inventory and the 5 failed areas: 66.7 min to the last logged event, i.e. about
  11 min per 1,000 searched groups.
- **Downloads:**
  - 3,561.5 MB from objects, TAP and 8 partitions (2,800 MB of light curves);
  - 239 MB of light curves per 1,000 searched groups in passing areas;
  - plus 1.75 MB of cutouts.

## 5. What remains (`remaining.json`, from the inventory only)
- **9 sampling-qualified partitions were not queued:** 28,445 unused groups, 23,945 searchable.
  - **4 of them have a worst-cell fraction ≥ 0.59:** 2900, 8124, 8126 and 8127. Together 19,702 unused groups (17,702
    searchable), 5.76 GB of light curves. All share ZTF fields with earlier samples.
  - **The other 5** (1356, 2488–2491; 8,743 unused) sit at 0.536–0.569, where the queued areas failed.
- **Not computed:** pass rates for these partitions.

## Deviations
- Inventory time: 48.7 min against the planned 30 min.
- The image-check `timeline.jsonl` entries were written by hand (start and end only).
