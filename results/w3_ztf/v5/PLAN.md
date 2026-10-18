# W3 ZTF v5 — exploratory search of the v1–v3 holdout half within the calibrated cells (plan)

Written 2026-10-18 before any holdout light curve is screened. Owner authorization: the round after #339. Code: the
commit that adds this file (`scripts/w3_ztf.py`, the #339 method; `--pix 3839 --vdir v5`).

## 1. Scope (a new, narrower question)
- #339's verdict, **failed**, stands and is not revised. This is an exploratory search for catalogue-proposed dimming
  events. It is not a pass of the full W3 search, and no abundance or limit is derived.
- Capability is stated only for the injection cells that reached the original 0.5 threshold. Numbers are recovered /
  injected from `results/w3_ztf/v4/injections_v4.json` (W3 n = 1, ε < 0; v4 carriers, not holdout carriers):

  | ρ | blend | t_E 10 d | t_E 30 d | t_E 100 d |
  |---|---|---|---|---|
  | 0.01 | 1 | 147/200 | 170/200 | 173/200 |
  | 0.01 | 0.7 | 142/200 | 165/200 | 148/200 |
  | 0.1 | 1 | 124/200 | 142/200 | — |
  | 0.1 | 0.7 | 128/200 | 120/200 | — |

- **Blind spots (not calibrated):**
  - t_E 3 d (79–90/200);
  - ρ 0.1 at t_E 100 d (81/200 and 60/200);
  - every other model, t_E, ρ or blend not in the table;
  - full disappearance as such, because "no catalogue record" is never counted as "vanished".
- Efficiencies cover the catalogue stage only, not image-to-verdict.

## 2. Sample and independence (`independence.json`, script `independence_script.py.txt`)
- The v1–v3 holdout half: 486 groups (2,788 objectids) in light-curve partition Npix 3839. All of them; no
  substitutes are added.
- Independence:
  - no shared objectid with the v1–v3 development half or with the v4 sample;
  - nearest development group 9.96″, nearest v4 group 77.2″;
  - 0 ambiguous groups;
  - no committed file names a holdout id.
- **Exposure history (stated, not hidden):**
  - v1–v3 downloaded the holdout members' light curves together with the development sample.
  - Holdout stars in the comparison magnitude range entered the aggregate per-exposure comparison counts.
  - No holdout light curve was screened, injected into or plotted. This is therefore not a strict blind test.
- **Shared systematics:** the same fields, cameras and observing windows as v1–v4 (for example the MJD 58655–58665
  window that dominated the v3-mapping flags). Common systematics can recur.

## 3. Frozen method
- `scripts/w3_ztf.py` exactly as merged in #339:
  - exposure IDs; one count per star and exposure;
  - target not its own comparison; in-play denominator;
  - gnomonic footprint test with a 10″ margin; objectid merging within 1.5″; expid dedupe;
  - flagged exposures unknown; common-mode test.
- Thresholds unchanged from #336 / v3:
  - DIM_MAG 1.0; MIN_RUN 3; MIN_NIGHTS 3; MIN_BRACKET 20;
  - COMMON_FRAC 0.2; INTER 0.3–1.0; MAX_INTER 1.
- The only code change is output plumbing: a `--vdir` folder and a `comparisons` stage that rebuilds the comparison
  list with the same criterion as `stage_sample`.
- No rule changes after the search. If a bug affecting the outcome is found, the results are kept and the round is
  marked affected.

## 4. Data and time limits
- Light curves: one HATS partition, Npix 3839, 7 columns, 944.9 MB compressed (parquet metadata).
- Comparison list: 5 metadata columns of the objects catalogue Npix 959, 79.9 MB compressed. This is a second file
  but contains no light curves.
- Total about 1.02 GB, under the 1.5 GB cap.
- TAP exposure metadata: batched per field, cached.
- Cutouts: ≤ 50 MB, logged separately.
- Timeouts: 300 s for TAP, 60 s for cutouts, ≤ 4 retries.
- Computation stops by 16:10Z; records and review follow.

## 5. Image check (selection fixed now)
- **Order:** groups with a run passing vetting (by key); then other flagged groups, by the number of non-detections in
  their longest flagged run (descending, then key). One run per group (the longest; ties: earliest). At most 5 groups.
- **Per group, at most 12 exposures (≤ 64 in total):**
  - before: the 2 nearest detected exposures before the run;
  - during: up to 8 catalogue non-detections or dim detections inside the run, evenly spaced in index if more;
  - after: the 2 nearest detected exposures after the run.
- Groups or exposures beyond the cap stay **pending** and are archived in the selection file; they are never counted
  as excluded.
- **Photometry:** forced photometry as in v4 (`image_script.py.txt`):
  - 120″ IBE cutouts; 2″ aperture with a 4–6″ annulus;
  - random-aperture error;
  - 5 field stars as same-exposure comparisons (field ratio must lie in 0.8–1.25).
- **Exposure classes:**
  - (1) not covered or invalid;
  - (2) poor image;
  - (3) source present, or in a quadrant without the group's objectid;
  - (4) credible dimming or non-detection: field-corrected ratio < 0.5 and deficit ≥ 3σ;
  - (5) undetermined.
  - A ratio < 1 alone is never dimming.
- **Per-source result:**
  - **real sky change:** ≥ 2 class-4 exposures on ≥ 2 nights, with a detection before and a recovery after (ratio
    ≥ 0.8, deficit < 3σ);
  - **model-compatible:** additionally the timing and shape agree with a calibrated cell at the actual epochs, with
    no curve filled in where there are no epochs;
  - **ordinary explanations open:** eclipse (repeats in the full light curve), dust (colour where the other filter
    has epochs), ordinary variable (catalogue lookups not done in this round: open), instrument.
  - Missing key phases are reported as **undetermined**.
- **Overall result classes:**
  - no flags: report the number searched and the capability;
  - flags all explained: report which class each falls in;
  - a real sky change: report it with its open ordinary explanations; nothing is called a wormhole candidate.
