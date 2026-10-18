# W3 ZTF v4 — independent validation on new sources (plan, fixed before any v4 light curve is read)

Written 2026-10-18 before the v4 light-curve partition was read. Code version: the commit that adds this file
(`scripts/w3_ztf.py`, `VERSION = 4` when run with `--pix 3838`). Owner authorization: the round after #338.

## 1. Data-mapping fixes (checked on v3 data before freezing)
- **Coverage from each exposure's footprint.** A non-detection counts only if the star is ≥ `EDGE_MARGIN` = 10″
  (ASSUMPTION) inside that exposure's quadrant polygon (IRSA `ztf_current_meta_sci` corners `ra1..dec4`, computed by
  IRSA from each science image's WCS). Comparison stars are "in play" at an exposure only under the same test (or a
  detection). Detections are kept regardless of the margin.
- **One physical star across quadrants and filters.** Groups merge every objectid (any field, quadrant, filter)
  within `GROUP_RADIUS` = 1.5″ of the anchor (unchanged since v1). v4 adds: a group with two objectids in the same
  (field, quadrant, filter) is ambiguous (two sources or a duplicate) and is dropped from the sample; per group the
  farthest member and the nearest non-member object within 5″ are recorded; within one filter each exposure (expid)
  enters a light curve once, a detection beating a non-detection.
- **Kept from v3:** each star counts once per exposure; the target is not its own comparison; official exposure IDs;
  the denominator counts only comparison stars in play (now also inside the footprint).
- **Check on v3 data** (`oldcheck.json`, script `oldcheck_script.py.txt`; the v3 light curves were not re-read: one
  partition per round): all 11 catalogue non-detections that the v3 image check classified fall outside the
  member's quadrant footprint under the 10″ rule; 0 of 1,000 v3 groups are ambiguous. This explains those 11 only;
  it does **not** clear the other v3 flags, which stay as they were.

## 2. Sample (metadata only)
- Light-curve partition Norder 5 Npix **3838** (rule: the smallest unread sibling of 3839 in the same order-4
  pixel 959; 846.0 MB; objects partition Npix 959, 574.2 MB, the same file as v1–v3; never re-read).
- v1–v3 selection cuts unchanged (r band, refmag 14.5–17.5, |sharp| ≤ 0.3, ngoodobs ≥ 300, non-transient,
  isolated), seed 7, ≤ 1,000 groups; v4 drops ambiguous groups and groups within 3″ of a v1–v3 group position.
  No dev/holdout split: all v4 groups are the validation set. The v1–v3 holdout stays unread.
- Same fields, cameras and observing windows as v1–v3 are likely (neighbouring sky): shared systematics
  (e.g. the MJD 58657–58661 and 59106–59109 windows) may recur. This is an independent *source* sample, not an
  independent survey.

## 3. Validation (unchanged thresholds and injection set from #336/v3)
- Pre-screen, vetting, injection cells (ρ 0.01/0.1 × blend 1/0.7 × t_E 3/10/30/100 d, n = 200 each), control classes
  (eclipse, fade, PSPL, common mode; n = 200 each), seeds (W3 cell c: 1100 + c; controls 1200 + k), hit rule
  (overlaps the umbra, not a run flagged in the uninjected carrier), noise model and pass conditions are the v3
  ones. Carriers: all v4 groups.
- **Pass (all must hold):** W3 recovery ≥ 0.5 in the t_E = 30 and 100 d cells for both ρ (blend 1); each control
  class passes vetting in ≤ 10 % of injections; the real v4 light curves leave ≤ 5 survivors after vetting (5 kept
  although the sample is up to twice the v3 development half). Computed by code from the JSON outputs.
- Missing records are never zero flux and never a measured dimming.
- Reported separately: full disappearance (umbra, ρ 0.01), short events (t_E 3–10 d), slow-edge events (ρ 0.1, blend
  0.7). The recovery is light-curve-table efficiency only, not image-to-verdict efficiency.

## 4. Image check (fixed before any flux)
- At most 5 flags, in this order: groups with a run passing vetting (by key), then other flagged groups ordered by
  the number of non-detections in their longest flagged run (descending, then key); one run per group (the
  longest; ties: earliest).
- Per flag ≤ 12 exposures: up to 8 catalogue non-detections inside the run (evenly spaced in index if more), and
  the 2 nearest detected exposures before and the 2 after the run in the same filter. Total ≤ 64. Unchecked
  non-detections stay "pending", never "excluded".
- Forced photometry as v3 (IRSA IBE cutouts, 120″, 2″ aperture, 4–6″ annulus, random-aperture error, 5 field stars
  of the same cutout as same-exposure comparisons; ≤ 50 MB of cutouts; 60-s timeouts, ≤ 4 retries).
- Classes: (1) not covered (position outside every image of the exposure) or invalid (infobits ≠ 0, non-finite
  aperture); (2) poor image (seeing > 4″ or maglimit < 19.5); (3) source present: flux ratio to the median of the
  bracketing detections (field-star corrected) ≥ 0.5 and the deficit < 3σ (combined error), or target position in
  a quadrant without the group's objectid; (4) credible dimming / non-detection: field-corrected ratio < 0.5 and
  the deficit ≥ 3σ, with field stars within 0.8–1.25 of their own reference; (5) undetermined: anything else.
  A ratio < 1 alone is never dimming.

## 5. Decision (pre-declared)
- All pass conditions hold → report the supported signal ranges and the size/cost of a next batch; nothing expands
  this round.
- Some fail → overall verdict "failed"; per-class capability and limits reported.
- Flags still dominated by missing records or common anomalies → keep a minimal reproducible case and stop the
  route; no further tuning round.
- Any class-4 dimming → ordinary explanations first; real dimming, model compatibility and exotic evidence kept
  apart. Nothing is called a wormhole candidate.

## 6. Budget and stop
- Downloads ≤ 1.5 GB (object partition columns + v4 light-curve partition + TAP), cutouts ≤ 50 MB, logged apart.
  Stop computing by 14:10Z; records and review after.

## 7. Sample as drawn (metadata only; `sample.json`, `merge_check.json`)
- 743,134 objects in the partition area; 4,781 pass the r-band cuts; 4,635 isolated; 4,633 groups; 2 ambiguous
  groups dropped; 0 near the v1–v3 sample; 1,000 drawn.
- Merge check: farthest member ≤ 0.56″ in every group (99th percentile 0.23″); nearest non-member object ≥ 2.07″
  (5 groups have one within 3″). No neighbour falls inside the 1.5″ tolerance in this sample.
