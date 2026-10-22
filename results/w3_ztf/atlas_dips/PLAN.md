# W3 ZTF × ATLAS: dip check on the stored measurements (post hoc diagnostic of #462)

Owner brief 2026-10-22 after #462 (session start 20:56 UTC):
- **Budget:** 90 min wall clock, ≤ 100 MB download.
- **Inputs:** only the stored measurement tables (`../atlas40/phot/`, `../atlas40/calib.json`, `../atlas40/diag_*.json`).
- **Excluded:**
  - no new ATLAS forced-photometry task;
  - no large ZTF partition read;
  - no injection rerun;
  - no two-point rule;
  - no completion of the 23 missing stars.
- **Scope:** a post hoc diagnostic. #462's thresholds, outcomes and verdict are unchanged. Nothing here can be a
  wormhole candidate.

## 0. Scope correction of #462 (wording only; numbers and verdict unchanged)
- **27 / 63:** the best-case upper bound under #462's fixed method, the data acquired and the existing A decisions,
  if every remaining event were restored. It is not a lower bound, and not a general ceiling for ATLAS or for other
  methods.
- **0 / 79 fade passes:** holds with ATLAS on 17 of the 40 stars. It does not say what all 40 would give. Under
  #462's fixed method the 40-star value lies in 0–4 / 79.
- **"All three routes failed"** narrows to: the tested options (ZTF g edge points, the g + r joint worth check,
  ATLAS edge points) did not meet the requirement on the specified sample under the specified criterion.

## 1. De-duplication (step 1; frozen before any per-exposure detail is read)
- **Cases:** the 8 W3 simulated mis-decisions of #462 (`diag_*.json`):
  - 4 A-recovered → slow;
  - 4 A-insufficient → slow.
- **Mapping:** each triggering ATLAS point is mapped back to its stored row in `phot/<star>.csv.gz` by carrier and
  exact relative flux and σ. The key is (star, filter, exposure `Obs`).
- **Listed per unique exposure:**
  - star, filter, unit, `Obs`;
  - time (MJD as published, and HMJD as used);
  - uJy, duJy, baseline, k;
  - quality fields (err, chi/N, x, y, maj, min, apfit, mag5sig, Sky);
  - how many of the 8 mis-decisions use it;
  - how injection turned it into an intermediate point.
- **Checks on each mapped row:**
  - the low flux is in the un-injected table;
  - the units (uJy) match the baseline's;
  - the time conversion matches;
  - the row is not a missing-value code (9999) or a duplicate publication (same `Obs` twice).
- **Coverage limit:** reuse is counted only over the 49 events diagnosed in #462. Counting it over all events needs a
  replay, which the brief excludes.
