# W3 ZTF v5 — exploratory holdout search within the calibrated cells: 0 flags in 486 groups

- Plan: `PLAN.md`, frozen in commit 24cdcbb before any holdout light curve was screened.
- Method: the #339 code, unchanged; the only addition is output plumbing.
- #339's verdict (**failed**) stands. This is a narrower exploratory search; no abundance or limit is derived.
- Labels: `derived` from observed ZTF DR24 data and IRSA metadata.

## 1. What was searched, and with what capability
- **Searched** (`search_holdout_v5.json`, `search_stats.json`):
  - 486 holdout groups;
  - 971 light-curve series in the r and g filters;
  - 1,438,841 series epochs (star-exposures); median span 2,767 d.
- **What the search saw:**
  - 531 catalogue non-detections that passed the footprint test;
  - 0 detections ≥ 1 mag below the median;
  - 16 dim runs of ≥ 2 epochs, none flagged by the frozen run rule.
- **Calibrated capability** (W3 n = 1, ε < 0; recovered / injected on the v4 carriers, `../v4/injections_v4.json`):

  | ρ | blend | t_E 10 d | t_E 30 d | t_E 100 d |
  |---|---|---|---|---|
  | 0.01 | 1 | 147/200 | 170/200 | 173/200 |
  | 0.01 | 0.7 | 142/200 | 165/200 | 148/200 |
  | 0.1 | 1 | 124/200 | 142/200 | — |
  | 0.1 | 0.7 | 128/200 | 120/200 | — |

  These are catalogue-stage efficiencies measured on other stars of the same fields. They are not image-to-verdict
  efficiencies, and they were not re-measured on holdout carriers.
- **Blind spots:**
  - t_E 3 d;
  - ρ 0.1 at t_E 100 d;
  - any model or parameter outside the table;
  - full disappearance as such: a missing catalogue record is never counted as "vanished".

## 2. Flags and image check
- 0 flagged runs, so there was nothing to check. No cutouts were downloaded (0 MB). The image script is committed
  (`image_script.py.txt`) but was not run.
- **Functional check** (post hoc; also the `dev` block of `search_stats.json`, development half, already seen in v1–v3; `search_dev_v5.json`): the same code
  gives 0 flagged runs on 514 development groups (1,023 series). With the v3 mapping these groups gave 1,520 flagged
  runs and 8 survivors (`../v3/search_dev_v3.json`). The change is the footprint test of #339; no rule was altered.

## 3. Real dimming
None was proposed by the catalogue stage, so none was checked. This is not evidence of absence, and no limit is
derived.

## 4. Independence and caveats (`independence.json`)
- No objectid is shared with the development half or the v4 sample. The nearest development group is 9.96″ away, the
  nearest v4 group 77.2″.
- The holdout was not strictly blind: v1–v3 downloaded its light curves, and the holdout stars in the comparison
  magnitude range entered aggregate per-exposure comparison counts. None was screened or inspected individually.
- The search uses the same fields, cameras and observing windows as v1–v4, so shared systematics are possible.

## 5. Cost (logged separately)
- Light-curve partition Npix 3839: 7 columns, 944.9 MB compressed (parquet metadata, read before fetching;
  `cost.json`).
- Objects catalogue Npix 959: 5 columns, 79.9 MB compressed (`cost.json`). Total 1,024.8 MB.
- IRSA TAP: 7.90 MB (`exposures_summary.json`).
- Cutouts: 0 MB.
- Wall time: fetch 39.6 s (`fetch_info.json`); search and statistics about 1 min each (not stored).
