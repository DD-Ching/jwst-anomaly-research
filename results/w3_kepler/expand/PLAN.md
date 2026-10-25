# W3 / Kepler fixed-rule limited expansion search: frozen plan

Owner authorization 2026-10-25. Start 21:49 UTC. Limits: 4 h (to 01:49 UTC), 6 GB through
`scripts/net_ledger.py` (`net_ledger.jsonl` here), 4 cores. This is an exploratory search for real
events worth checking. It is not a claim that the method can identify wormholes.

**Committed before any new flux is read:** this plan, `scripts/w3_kepler_expand.py`, its tests and
`sample_expand.json`. Before this commit, the ledger holds only the Kirk et al. 2016 EB catalogue
query (metadata).

## 0. Orientation and previous counts
- **Repository:** `main` at cada17c. No open PR or branch works on Kepler; `claude/w3-kepler-pilot`
  and `claude/w3-kepler-hard` are merged (#500, #529).
- **The "425 stars" figure:** it came from my chat answer after #529 and is not in the repository.
  It was wrong. From `../sample.json` and `../hard/sample_hard.json`, deduplicated:
  - development 20;
  - validation 45 (#500 20 + #529 25);
  - search 260 (#500 60 + #529 200);
  - ordinary-event controls 46 (EB stars);
  - interface check 1.
  - In total 372 distinct stars, of which 325 are FGK pool stars (development + validation +
    search). Stored as `previous_rounds` in `sample_expand.json`.

## 1. Sample (`sample_expand.json`; metadata only)
- **Pool:** the same as #500 / #529: NExScI `keplerstellar` delivery `q1_q17_dr25_stellar`, logg ≥ 4.0,
  4500–6500 K, Kp 12–14, all 17 quarters. The pinned sha256 holds: 18,570 stars.
- **Exclusions:** every one of the 372 stars read before, whatever its role. That leaves 18,240
  candidates.
- **Selection:** sort by sha256(`w3-kepler-expand-2026-10-25:<kepid>`) and take the **first 500**. The
  list and its order are fixed.
  - No star is swapped for light-curve shape, result or download failure.
  - A failed download is recorded in `files_expand.json.gz` and counted as "selected, not read".
- **Cost estimate:** 500 × ~7.7 MB ≈ 3.9 GB of light curves, leaving ~2 GB for vetting. That is within
  the 6 GB cap, so no reduction is needed.

## 2. Processing (unchanged from #529)
- **Product:** Kepler DR25 long-cadence `_llc.fits` files from MAST, Q1–Q17, `SAP_FLUX`.
  `PDCSAP_FLUX` is a diagnostic only.
- **Rule:** `w3_kepler.screen_star(..., local=True)`, unchanged. The parts:
  - the quality mask (`QMASK`);
  - per-quarter median join;
  - 1-mag dim threshold;
  - ≥ 3 consecutive dim cadences, with ≥ 20 non-dim cadences on each side;
  - common mode < 0.2 (leave-one-out over the 500 read stars; `frac_loo`, equal to
    `w3_kepler.frac_for` by test);
  - bounded fast edges (v23 bounds, per quarter);
  - one-off;
  - the #529 local shoulder rule.
  - A run "survives" if `pass_local` is true.
- **Changes:** no threshold, no second method, no injection, no new photometry or detrending.

## 3. Validity of the existing recovery numbers
- **Source:** recovery comes from #529's validation (25 stars, `../hard/validation.json`: 55–86 / 100
  per cell, weakest 55). It is **not re-measured on this batch**. The weakest cell does not show a
  population recovery above 50 %.
- **Applicability (per star, stored):** all 17 quarters present, and median σ_rel ≤ 2 × the #500
  maximum (3.96 × 10⁻⁴; `../posthoc2.json`).
- **Stop rule:** if fewer than 90 % of the read stars are applicable, the recovery numbers are not
  quoted for this batch. The layer counts are still reported. Stars outside the range are screened
  and counted separately.

## 4. Stored layers (`layers.json`, `screen_expand.json.gz`, `runs_expand.json.gz`)
1. Selected (500), read, screened, applicable; effective observing time = Σ good cadences × 29.4 min.
2. a. Stars with any good cadence ≥ 1 mag below the quarter median, and their total cadences.
   b. Flagged deep runs and their sources.
3. Runs left after each test in order (common mode → edges → one-off → local shoulder), and each test
   on its own.
4. Survivor runs and sources. For every flagged source, the first test its best run failed, and
   whether it is in the Kirk et al. 2016 EB catalogue (`inputs/kirk2016_catalog.csv`).

If layer 2a is empty, the result is written as "this batch did not trigger the pre-screen". It is
not credited to the shoulder test.

## 5. Vetting (only if there are survivors)
- **Order:** survivor sources in sample order (then run start time). At most **10** are checked;
  the rest are listed as "not vetted", not as excluded.
- **Per survivor, with existing data and tools only:**
  - `SAP_QUALITY` (all bits) within ±1 d of the run;
  - flux-weighted centroids `MOM_CENTR1/2` in the run and its shoulders against the quarter;
  - same-exposure comparison: the other read stars on the same channel and quarter at the same
    cadences;
  - quarter and aperture changes (run near a quarter boundary; aperture header);
  - repeated events elsewhere in the light curve (other runs, and dips ≥ 0.3 mag);
  - PDCSAP, never an independent confirmation;
  - if needed, the official long-cadence target pixel file for the run's quarter.
- **Specific checks:** whether each shoulder could be a flare (shape: fast rise, slower decay,
  duration), a baseline drift or an instrumental feature.
- **Classification:** "confirmed real dimming", "artefact with evidence" or "undetermined /
  not vetted". No new rule is invented during vetting.
- **Interpretation:** a real dimming is not W3, and not a wormhole.

## 6. Budget and stop
- **Time:** if the 4 h would be exceeded, the finished layers and the unchecked list are reported.
  Nothing is dropped from the denominators.
- **Data:** if the ledger would pass 6 GB, the fetch stops. The unread stars stay "selected, not
  read".
- **Fixes:** an implementation error is fixed and the affected stage re-run, with both outputs kept.
- **Labels:** light curves and flags are observed; counts are derived; there are no injections in
  this round.
