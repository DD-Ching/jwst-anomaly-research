# E1 + Einstein Probe: result (2026-10-18)

Plan: `PLAN.md` (frozen first, committed before the count). Code: `scripts/e1_ep.py` (`build`, `count`; frozen
rules `control_passes`, `outcome`, `test_mask`, unit-tested in `tests/test_e1_ep.py`), post hoc diagnostic
`scripts/e1_ep_diag.py`. Inputs: GCN Circulars archive pinned in `data/manifests/e1_ep.ecsv`; GBM and BAT reused
from `data/manifests/e1_swift.ecsv`.

**Node table (derived, `ep_events.ecsv`, `build_summary.json`):** 231 EP/WXT transients, MJD 60384.8–61328.7, one
row per name with its circular id; 218 with an arcmin (WXT) radius, 13 with only an arcsec radius; time from a
stated transient onset for 76, else the first timestamp (155). 8 names with a WXT mention do not parse or were
called "not a real source" (listed in `build_summary.json`). The archive is streamed, never stored.

**Positive control (frozen rule passes):** EP–BAT + EP–GBM same-direction pairs at 0–1 h: 27 vs a jit-null mean
of 0.63 (sd 0.80), analytic p 2.2 × 10⁻³⁴ (`summary.json`).

**Tested family (18 cells, 10⁴ jit scrambles; `counts.ecsv`): outcome "null".**
- Smallest cell: EP–GBM 1 h–1 d same, 27 pairs vs 14.53 ± 3.76; analytic p 2.2 × 10⁻³; × 18 = 0.039 (rule:
  ≤ 1.35 × 10⁻³).
- Pooled global p 0.017 (jit) (`summary.json`).
- **Deviation:** the plan's `jitday` vetting null was dropped: no catalogue here is ground-based, so its day
  quantum applies to none and it equals `jit` (review finding).
- EP–EP: no same-direction pair within 7 d (`same_dir_pairs_EP-EP.ecsv` is empty); wide cells within 1σ of the null.

**Limits (`limits.ecsv`, D-074 method, EP–EP only, all 5 rows valid):** 95 % upper limits on dependent wide EP–EP
pairs per EP transient (231 anchors): 1.3 × 10⁻² (0–10 s, 10–100 s), 3.4 × 10⁻² (100 s–1 h), 5.8 × 10⁻² (1 h–1 d),
9.7 × 10⁻² (1–7 d).

**The smallest cell, looked at (post hoc; `same_dir_pairs_EP-GBM.ecsv`, 27 rows at 1 h–1 d):**
- 5 of the 27 are one source seen twice: EP250702a with the three GBM triggers bn250702548 / 581 / 682 at 0.011°
  (GCN 40906: EP250702a "likely associated with GRB 250702B,C,D,E"); EP250919a / bn250919020 (0.002°) and
  EP261007a / bn261007236 (0.054°), where the EP time is 1–2 h after the GBM trigger. Without them the cell is
  22 vs 14.53 (not stored).
- The other 22 have separations 4.7–94.6° at combined σ 6.1–31.8°: GBM's large errors make "same direction" a
  wide region.
- **Deviation:** D-074 addendum 6 asks that a satellite × satellite channel drop joint detections first; the frozen
  plan did not. The post hoc re-count without them (`diag_no_joint.json`: 24 EP, 9 BAT, 18 GBM events removed,
  10⁴ jit scrambles) gives EP–GBM 1 h–1 d same 26 vs 12.85 ± 3.59, p 8.2 × 10⁻⁴, × 18 = 0.015: still above the
  family-wise rule. The joint-event rule (|lag| < 1 h) keeps EP250702a, EP250919a and EP261007a, whose EP times are
  more than 1 h from GBM.
- **Failed approaches:** the first run's `control_sum` mixed axes on the null ensemble (numpy advanced-index
  order; control null mean 3,175); fixed before any cell was read. The first build took "the observation started
  at" as an onset and accepted WXT mentioned up to 400 characters after the position; both fixed after review
  (archive sha256 unchanged), and every number here is from the re-run.

**Limitations:** EP times are trigger, onset or observation start (up to ~20 min apart); the table holds only
transients reported in circulars (a selected subset of WXT detections); positions with arcsec radii are FXT or
refined ones. The circular archive changes daily, so the parsed table, not a re-parse, is the record.
