# W3 ZTF v21 — the v20 survivor's NDs cannot discriminate; 8124 and 8127 both pass and are searched (9,978 groups); 1 final survivor, undetermined (a low-response column under the source)

- **Plan, queue and samples:** committed before any light curve and before the survivor check (ea80798).
- **Labels:** DR24, TAP, IBE images, masks and psfcat are `observed`; counts are `derived`; injections are
  `simulated`.
- **Scope:**
  - the frozen v13 W3 rule, the 10 cells, the v16 thresholds and cuts;
  - no limit; nothing here bears on wormholes;
  - the cumulative result is "no confirmed target event", not "no star dimmed".
- **Unchanged:** 1624 and 1625 stay archived as in v19.

## A. v20 survivor 8126 / 634210200007847 (`survivor_check.json`; 7 min, about 1.4 MB, not stored as one number)
- **Correction to v20:** the group key 634210200007847 is the **r-band** objectid. The run was in g, on member
  634110200003853 (field 634, rcid 37, g).
  - The first run of the check used the key's r light curve. That was fixed before recording.
  - **The g baseline is 18.92 ± 0.05** (median of 1,553 clean DR24 epochs). v20's "g ≈ 19.3" came from an
    uncorrected 2″ aperture and is superseded.
- **The 4 NDs** (frozen criteria, PLAN A):

| ND | maglimit | seeing | comparable stars found (of 323) | depth ≥ baseline + 0.75 | completeness ≥ 0.9 | able to discriminate |
|---|---|---|---|---|---|---|
| 45450718 | 19.31 | 3.27″ | 230 (71 %) | no | no | **no** |
| 45546604 | 19.52 | 3.58″ | 268 (83 %) | no | no | **no** |
| 45548609 | 19.51 | 3.60″ | 286 (89 %) | no | no | **no** |
| 45636370 | 19.24 | 2.67″ | 267 (83 %) | no | no | **no** |

- **Reading:** none of the 4 NDs could distinguish the rule's dimming.
  - The images reach only 0.3–0.6 mag below the baseline.
  - In the same exposures the catalogue misses 11–29 % of isolated stars of the same brightness.
  - So a catalogue miss of an unchanged source is expected there.
- **What this does not settle:** no target flux was measured (no validated method), so the NDs' flux stays
  **unknown**. The survivor remains **undetermined**, and this is not proof that it did not dim. Infobits are 0 and
  there are no bad mask bits.

## B. Queue (`summary.json`, each area's `decision.json`)

| pos | Npix | verdict | lowest cell / 200 | fade passing / 1,000 | calibration survivors | searched |
|---|---|---|---|---|---|---|
| 1 | 8124 | **passed** | 126 | 74 | 1 | 5,104 groups (10,199 series), 29 runs in 14 groups, **1 survivor** |
| 2 | 8127 | **passed** | 118 | 54 | 0 | 4,874 groups (9,745 series), 24 runs in 5 groups, 0 survivors |

- Eclipse, PSPL and common-mode controls were 0 / 200 in both areas.
- **Exposure independence:** every field of both areas already had earlier sample data, so no exposure is
  independent of earlier data. The groups are in no earlier sample group: they were deduplicated against 53,529
  groups by 3″ or a shared objectid, and 0 were dropped for a shared objectid.
- **Deviation:**
  - The first 8127 injection job was killed by the system (rc −9, no output written; most likely memory, with 4
    heavy jobs running at once). The driver then stopped that area without a verdict.
  - The same validation (seed 31000) was run once more, alone, at 06:15–06:19 (`timeline.jsonl`, hand-logged
    note). That is not a retry after a verdict.
  - `summary.json` stage-seconds count only the second run. The first used about 174 s: inject 82 s and dev
    search 92 s (not stored).

## C. Final survivor and image check (`8124/image_check.json`, `mask_check.json`, `image_check.png`; 7 science and 3 mask cutouts, 0.62 MB)
- **Source:** 8124 / 634205400002525, r band, 3 NDs (HMJD 59897.10–59899.08), field 634. A bright source: reference
  3,300 µJy, refmag r 14.65.
- **Automated classes:** 2 × class 4 ("credible dimming", 2″ ratio 0.03) and 1 × poor image. The support test
  reports `real_sky_change` true.
- **Official mask:** bit 2 **LOW RESPONSIVITY** at the source pixel in all 3 NDs (centre value 4). The dark column
  lies 0.27–0.50″ from the source in the NDs and 1.74–15.63″ away in the 4 detections.
- **Visual reading** (`image_check.png`, separate from the classes): the star is bright and clearly present in all
  3 NDs, with the black column through or beside its core.
- **Reading:**
  - The automated class 4 is not reliable here, because the 2″ aperture includes official low-responsivity pixels.
  - The catalogue miss is consistent with the column, as in v16's 1624.
  - **Undetermined**, with a data problem supported; **not** a confirmed dimming. No validated measurement on the
    masked core exists, so the flux stays unknown.
  - No variable-star or known-event check was opened.

## D. Totals and cost (measured or estimated, as marked)
- **Searched this round:** 9,978 new groups (19,944 series).
- **Cumulative with the v13 rule:** 35,475 groups (25,497 + 9,978; not stored as one number). No confirmed target
  event.
- **Final survivors this round:** 1, plus the v20 survivor carried over.
  - Reliable real dimmings: 0.
  - Undetermined: 2 (8124 / 634205400002525 new; 8126 / 634210200007847 carried, now with NDs shown unable to
    discriminate).
- **Downloads:** 3,558.8 MB, from the objects table (97.3 MB), light curves from parquet metadata, and TAP 23 MB.
  Add about 1.4 MB for part A and 0.62 MB of cutouts (not stored).
- **Time:**
  - the queue, 06:05:07–06:26, is 21.4 min including the image check and the 8127 rerun;
  - part A took about 7 min (not stored);
  - **per 1,000 searched groups:** 167 stage-seconds, 344.6 MB of light curves, 2.1 queue minutes.
- **Review load:** 1 image-checked group (cap 10) and 10 cutouts (cap 100).

## E. Inventory status (`remaining.json`)
- **All 17 sampling-qualified partitions of the v16 25-partition inventory have now started.** Nothing is left
  unstarted.
- The 6,243 search groups of the v20 failed areas stay unsearched under the frozen rule.
- Any further search needs a new inventory, which was not authorized this round.
