# W3 ZTF v4 — independent new-source validation: failed (one W3 cell); 0 flags on 1,000 new stars

Plan: `PLAN.md` (frozen in commit 7fcd456 before any v4 light curve was read). Every number below comes from a JSON
file in this folder. Labels: injections are `simulated`; everything else `derived` from observed ZTF DR24 data and
IRSA metadata.

## Verdict (frozen rule, applied by code: `decision.json`, script `decision_script.py.txt`)
- **Failed.** The W3 cells of the pass rule must each reach 0.5 (t_E 30 and 100 d, ρ 0.01 and 0.1, blend 1). They
  reached 0.845, 0.865, 0.705 and **0.405** (ρ 0.1, t_E 100 d).
- The other two conditions hold:
  - controls passing vetting: eclipse 0/200, fade 9/200, PSPL 0/200, common mode 0/200;
  - survivors: 0 of 1,000 groups (`search_all_v4.json`; 1,991 series searched, 0 flagged runs).
- This is a source-independent sample from the neighbouring partition (Npix 3838). It shares the same fields,
  cameras and observing windows as v1–v3, so shared systematics remain possible. The v1–v3 holdout is still unread.

## 1. Did the fixes solve the known mapping problem?
- **On v3 data** (`oldcheck.json`): the v4 footprint test (≥ 10″ inside the exposure's quadrant corners) puts all
  11 image-classified v3 catalogue non-detections outside the member's quadrant. These are the 4 with no image,
  the 5 with the star in the neighbouring quadrant, the infobits case (the star was in the neighbouring quadrant)
  and the non-finite-aperture case. This covers those 11 only; the other v3 flags were not re-checked (one partition
  per round).
- **On the v4 sample:**
  - Under the v3 mapping, the same 1,000 groups give 2,565 flagged runs in 144 groups, 10 of which pass vetting
    (`search_all_v3rule_posthoc.json`, post hoc).
  - Under the v4 mapping: 0 flagged runs.
- **Image check of what the fix removed** (post hoc deviation; v4 had no flag to check, so the frozen selection rule
  was applied to the v3-mapping flags; `image_check_posthoc.json`, `stamps_posthoc.png`):
  - Coverage: 5 flags (all pass vetting under the v3 mapping, all in MJD 58655–58665), 44 cutouts, 38 catalogue
    non-detections checked and 6 left pending.
  - Classes of the 38 non-detections:
    - 21 class 3: the star is present in a quadrant that has no objectid for it. Where the ratio is defined it is
      0.96–1.04 with deficits < 3σ.
    - 17 class 1: 14 not covered (position outside every image of the exposure, IBE) and 3 non-finite aperture
      (image edge).
    - 0 class 4 (credible dimming) and 0 undetermined.
  - For all 38, the v4 footprint test says "outside the member quadrant", so the fix removed only exposures that
    the images show as uncovered or off-quadrant.
- **Merge across quadrants and filters** (`merge_check.json`, v4 sample):
  - Members lie ≤ 0.56″ from the anchor; the nearest non-member object is ≥ 2.07″ away. No neighbour is inside the
    1.5″ tolerance.
  - 2 ambiguous groups (two objectids in one quadrant and filter) were dropped; 0 of 1,000 v3 groups are ambiguous.
  - The class-3 cases have no objectid in the other quadrant, so merging cannot recover them; the footprint test
    excludes them instead.

## 2. Injections on the new sample (`injections_v4.json`, 200 per cell, 1,000 carriers)
- Recovered / injected, by event class (`decision.json`):
  - **Full disappearance (ρ 0.01):**
    - blend 1: 90, 145, 169, 173 at t_E 3, 10, 30, 100 d;
    - blend 0.7: 90, 141, 165, 148.
  - **Short (t_E 3, 10 d):**
    - ρ 0.01: 90, 145 (blend 1); 90, 141 (blend 0.7);
    - ρ 0.1: 84, 137 (blend 1); 80, 128 (blend 0.7).
  - **Slow edges (ρ 0.1):**
    - blend 1: 84, 137, 141, 81 at t_E 3, 10, 30, 100 d;
    - blend 0.7: 80, 128, 120, 60.
- **Short-event limit is flagging:** at t_E 3 d only 83–92 of 200 are flagged at all (fewer than 3 dim epochs on
  3 nights).
- **Slow-edge limit is vetting** (`slowedge_posthoc.json`, post hoc, same draws). Of the flagged but rejected
  injections, the sharp-edge test is failed by every hit run in 84 of 84 (ρ 0.1, blend 1, t_E 100 d) and 94 of 102
  (blend 0.7, t_E 100 d).
- These are light-curve-table efficiencies. They do not include the image step or the footprint margin. The
  fraction of exposures lost to the 10″ margin was not stored.
- v3 → v4 differences mix the code change with a new sample: 514 dev carriers before, 1,000 new carriers now.

## 3. Image check of v4 flags
None: the v4 search flagged nothing. No credible dimming was found in any checked image (section 1).

## Cost (logged separately)
- Archive reads: objects partition Npix 959 (11 columns of a 574.2 MB file) and light-curve partition Npix 3838
  (7 columns of an 846.0 MB file). The upper bound is about 1.42 GB, from file sizes; exact transferred bytes were
  not stored.
- Queries: IRSA TAP 6.47 MB (`exposures_summary.json`) plus 0.18 MB (`oldcheck.json`); IBE search 4.6 kB.
- Cutouts: 3.10 MB, 44 cutouts (`image_check_posthoc.json`).
- Wall time (not stored): sample 68 s, fetch 54 s, search and injections a few minutes each.

## Decision (pre-declared, PLAN section 5)
- Partial pass, so the overall verdict is **failed**. Capability per class:
  - full disappearance at t_E ≥ 10 d: 0.70–0.87;
  - t_E 3 d: 0.40–0.45;
  - ρ 0.1, t_E 100 d slow edges: 0.30–0.41, limited by the frozen sharp-edge test.
- No credible dimming. W3 stays a restricted hypothesis; nothing is a wormhole candidate. Nothing expands this round.
