# W3 ZTF v6 — cross-sky validation, then one fixed search (plan, frozen before any light curve is read)

Written 2026-10-18 before any light curve in the new sky area was read. Owner authorization: the round after #342.
Both stages are authorized: if validation passes, the search group is searched once; if it fails, the round stops.

## 1. Scope
- The full validation in #339 failed, and that verdict stands. This round keeps #342's narrow model range: the 10
  injection cells that reached 0.5 in v4 (W3 negative-mass point lens, n = 1, ε < 0).
- Short events (t_E 3 d) and slow edges (ρ 0.1, t_E 100 d) were **tested and found weak**, not untested. They stay
  blind spots and are not needed for this limited search.
- #342's holdout entered aggregate comparison counts, so it was not a strict blind test. Its zero flags do not show
  that the method failed, and they rule nothing out.
- All efficiencies are catalogue-stage only; none is an image-level efficiency for full disappearance.

## 2. Sky area (chosen from metadata only)
- **Rule, applied to IRSA's HATS `partition_info.csv`:**
  - light-curve partitions at Norder 5 whose order-4 parent has an objects partition;
  - centre at dec +15 to +35 and |b| > 20°;
  - the one with the row count closest to that of Npix 3839 (508,684), ties broken by Npix.
- **Result:** Npix **7134** (RA 184.2, Dec +32.8, b +80.6; 508,709 rows). Objects come from Norder 4 Npix 1783.
- **Independence from v1–v5** (`independence.json`):
  - new ZTF fields 671, 672, 1669, 1713, 1714; old fields 826, 848, 1831, 1856, 1857;
  - no shared field or quadrant, so no shared exposure ID.
- **Still shared:** the ZTF camera, DR24 processing, the 2018–2024 survey period, the IRSA metadata and this code.
  Systematics common to the instrument or processing are not made independent.

## 3. Sources (`sample.json`, fixed before any light curve)
- Cuts as in v1–v5: r band, refmag 14.5–17.5, |sharp| ≤ 0.3, ngoodobs ≥ 300, non-transient, isolated. All objectids
  of a star (any field, quadrant or filter) within 1.5″ form one group. Ambiguous groups and groups near the v1–v3
  sample are dropped. Up to 5,000 groups, seed 7.
- **Drawn:** 1,072 groups (all that pass; none are added). There are 0 ambiguous groups.
- **Three fixed sets:**
  - **search group:** crc32("v6:" + key) even, 539 groups;
  - **calibration group:** the other 533 groups;
  - **control group:** 5,591 objects in the comparison range (count from the cache; not stored) (refmag 14–18, non-transient) that belong to neither
    of the other two. It is the only source of per-exposure comparison counts.
- **Enforcement in code:** `context()` loads only control rows plus the members of the groups in use. Each stage
  appends a line to `source_usage.jsonl` (rows loaded, search-member rows). During calibration that count must be 0.

## 4. Method (frozen; the #342 code plus v6 plumbing)
- Unchanged from #339:
  - exposure IDs; one count per star and exposure;
  - the target is not its own comparison; in-play denominator;
  - gnomonic footprint test with a 10″ margin; objectid merging; expid dedupe;
  - flagged exposures unknown; common-mode test.
- Thresholds unchanged: DIM_MAG 1.0, MIN_RUN 3, MIN_NIGHTS 3, MIN_BRACKET 20, COMMON_FRAC 0.2, INTER 0.3–1.0,
  MAX_INTER 1.
- Injections: the 10 cells, n = 200 each, with their original seeds (1100 + cell index). Noise, missing-epoch and
  limit handling are as before. Controls (eclipse, fade, PSPL, common mode) use n = 200 and seeds 1200 + k. Carriers
  come from the calibration group only.
- Each cell reports recovered / injected and the number of distinct stars used. Injections on one star are not
  independent sky samples.

## 5. Pass conditions (all required; computed by `decision_script.py.txt`)
- Each of the 10 cells recovers ≥ 0.5.
- Each control class passes vetting in ≤ 10 % of injections.
- The real calibration light curves leave ≤ 5 survivors after vetting (the original burden rule).
- If any condition fails, the failure is recorded and the round stops: the search group is not read and nothing is
  re-tuned.

## 6. Search and image check (only if section 5 passes)
- One search of the 539 search groups (`search --half holdout`), with no added sources and no rule changes.
- The image check uses `image_script.py.txt`, the v5 script configured for v6:
  - **order:** runs passing vetting first (by key), then other flagged groups by the number of non-detections in
    their longest run;
  - **capacity:** ≤ 5 flags and ≤ 12 exposures each (≤ 64 in total);
  - **sampling:** 2 detected exposures before the run, up to 8 inside it (evenly spaced), 2 after;
  - anything beyond the cap is archived as pending.
- Classes and the per-source rule are as in v5 PLAN section 5:
  - **real sky change:** ≥ 2 class-4 exposures on ≥ 2 nights, plus a detection before and a recovery after;
  - **W3-compatible:** timing and shape match a calibrated cell at the actual epochs;
  - **ordinary explanations:** eclipse, dust, variable star and instrument are each recorded as open or closed.
  - A missing record is never zero flux. A missing key phase means undetermined.

## 7. Limits
- Downloads: light-curve partition 7134, 7 columns, 299.9 MB compressed; objects Npix 1783, 11 columns, 86.8 MB
  (both from parquet metadata). Total 386.7 MB, under the 1.5 GB cap. Cutouts ≤ 50 MB, logged separately.
- Network: TAP batched per field and cached; 60 s timeouts for cutouts, ≤ 4 retries.
- Computation stops by 21:40Z.

## Deviations
- After a TAP connection reset, `stage_exposures` gained a bounded retry (4 tries, 10–30 s back-off). This is
  network handling only; the stage was re-run before any light-curve statistic.
- Added after the failure: `failure_script.py.txt`, a post hoc diagnostic run on the calibration group only.
