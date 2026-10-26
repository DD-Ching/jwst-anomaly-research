# W3 / Kepler fixed-rule search of new stars: 1,500 of 2,000 read and searched; no star reaches the 1-mag pre-screen, 0 survivors; 3,710–4,855 effective star-years per W3 cell (4,959–6,475 with #546)

- **Pre-registration:** `PLAN.md`, `sample_2000.json`, `scripts/w3_kepler_2000.py` and its tests were committed at
  04:55:06 UTC (89d041f), before any new flux. The batch-1 fetch started at 04:55:22 (`cost.json` stage log).
- **What it is:** a fixed-rule exploratory search with the #529 / #542 screen and the #546 injection-recovery
  measurement, unchanged. A survivor would not be wormhole evidence; a zero is not an exclusion.
- **Labels:** light curves and flags are observed; injections simulated; counts, efficiencies and limits derived.
- **Deviation:** batch 4 (500 stars) was not fetched (`batch4_decision.json`).
  - The plan allowed a fetch start until 09:15 UTC. At 08:49, batch 3's last 30 minutes ran at 466,733 B/s, which
    projected batch 4's download to end at 11:07, after the 6-h limit (about 10:45).
  - Its 500 stars are unread, not replaced.
- **Notes:**
  - **Implementation error, fixed:** the first batch-1 screen call crashed before writing anything. The shared
    `w3_kepler_expand.layers` reads an `applicable` field that the new screen rows lacked. The field was added with
    #542's definition (descriptive only) and the screen re-run.
  - **V9 code:** the pixel-check code (procedure fixed in PLAN §6) was committed at 04:56:08 (2980edc), 46 s after
    the batch-1 fetch began and before any vetting. It was never used.
  - **Download connections:** batch 3 used 16 download connections instead of 8 (`scripts/CLAUDE.md` allows 8–16),
    because files were latency-bound: the median per-file time was 3.1 s in batch 2 with 8 connections (ledger; not
    stored). Throughput did not improve: the median was 10.5 s in batch 3 (not stored). The environment's download
    rate fell from batch to batch: 1,771.6 s, 4,621.1 s and 7,369.4 s per batch for the same about 3.7 GB
    (`cost.json`).
  - **Mirror test:** one file from the MAST AWS mirror was fetched to test its speed. Its sha256 equals the archive
    copy. It took 12.6 s, against about 8 s for archive files at the same time (ledger; not stored), and it was not
    used.
  - **Ledger snapshot:** the ledger snapshot is split in two gzip parts (`net_ledger.part1/2.jsonl.gz`) to keep each
    file under 1 MB.

## 1. Search results (`layers_b<b>.json`; every star in `screen_b<b>.json.gz`)

| batch | selected | read | unread | observing star-yr | stars with a ≥ 1-mag cadence | flagged runs | flagged sources | survivors |
|---|---|---|---|---|---|---|---|---|
| 1 | 500 | 500 | 0 | 1,731.3 | 0 | 0 | 0 | 0 |
| 2 | 500 | 500 | 0 | 1,732.3 | 0 | 0 | 0 | 0 |
| 3 | 500 | 500 | 0 | 1,731.6 | 0 | 0 | 0 | 0 |
| 4 | 500 | 0 (not fetched) | 500 | — | — | — | — | — |

- **Pre-screen:** no read star has any good cadence 1 mag or more below its quarter median.
  - The lowest relative flux of any star is 0.51 (batch 1), 0.6243 (batch 2) and 0.7306 (batch 3) (per-star
    `min_rel`).
  - So no run was flagged, nothing reached the common-mode, edge, one-off or shoulder tests, and there was nothing
    to vet. The contact sheet is empty.
- **Files:** every read star has all 17 quarters, except KIC 6350306 (16 quarters; its listing has one file fewer).
  It was screened like the others.
- **EBs:** the Kirk 2016 EBs in the sample were kept: 36 of the 2,000 selected (`sample_2000.json`), 31 of them in
  batches 1–3 (from the `batch` field; not stored). None of the 31 read ones reached 1 mag in a good cadence.
- **Downloads:** 0 failures. Per batch, 8,500, 8,499 and 8,500 Q1–Q17 files (3,719,520,000, 3,719,013,120 and
  3,719,520,000 B; `files_b<b>.json.gz`). Ledger total: 11,391,544,013 B of the 18 GB cap.

## 2. Vetting
None. With no flagged run there was no survivor:
- confirmed real dimmings: 0;
- supported artifacts: 0;
- undetermined: 0;
- not vetted: 0.

The 2 GB reserve was not used.

## 3. Search capability per batch (`summary_b<b>.json`; every event in `events_b<b>.json.gz`)

| batch | stars | events | recovery rate (range over cells) | E per cell (star-yr) | spurious |
|---|---|---|---|---|---|
| 1 | 500 | 20,000 | 0.633–0.803 | 1,264.7–1,604.3 | 0 |
| 2 | 500 | 20,000 | 0.61–0.811 | 1,218.7–1,620.3 | 0 |
| 3 | 500 | 20,000 | 0.61–0.84 | 1,218.7–1,678.3 | 0 |

- **Comparison with #546:** per-batch recovery ranges are 0.61–0.84, against 0.625–0.811 for #546 on its 500 stars.
  The weakest cells are again t_E = 10 d.
  - The highest value, batch 3's 0.84, is cell 2 (ρ 0.01, blend 1, t_E 1 d). Batches 1 and 2 have 0.788 and 0.802
    there, and #546 has 0.811.
  - Against the mean of the other three, that is 2.7σ (binomial, 1,000 events each). Over all 80 cell × batch values
    (3 batches + #546), χ²/dof ≈ 0.79 and the largest deviation from a cell mean is 2.4σ (not stored).
  - So the batches scatter about as binomial counts would. No cell was dropped and no seed was changed.

## 4. Round totals (`combined.json`)
- **This round:** 1,500 stars; observing 5,195.2 star-years (good cadences × 29.4 min); sampling 5,993.8
  star-years (windows, gaps included).
  - Effective search volume per cell: 3,710.2–4,855.0 star-years.
  - Joint limits: 6.173–8.08 × 10⁻⁴ per star-year.
- **This round plus #546:** combined, because #546 used the same rule, event distribution, injection method,
  recovery definition and batch-of-500 common-mode set. 2,000 stars; observing 6,927.5; sampling 7,991.7
  star-years.
  - Effective search volume per cell: 4,958.9–6,475.3 star-years (cell 19 ρ 0.1, blend 0.7, t_E 10 d lowest;
    cell 2 ρ 0.01, blend 1, t_E 1 d highest).
  - Joint limits: 4.628–6.044 × 10⁻⁴ per star-year.
- **Not added:** the earlier 60 (#500 rule) and 200 (#529 search, other comparison set, no #546-style
  injection).

## 5. Limits: method and the #546 check (`limit_check_546.json`)
- **Joint limit (PLAN §7):**
  - a flat prior on the rate;
  - the efficiency integrated over its Beta posterior;
  - the count and the efficiency uncertainty handled in one step, not two stacked 95 % steps.
  - For an exactly known efficiency it equals the classical 2.996 / E. It is a Bayesian credible limit, and it
    treats the draws as independent given the stars.
  - It follows the idea of Cousins & Highland 1992 (averaging over the sensitivity's uncertainty), not their exact
    formula.
- **Correction of #546:** #546 quoted Γ95 = 2.996 / E_low, i.e. the Poisson 95 % step at the low end of a 95 %
  efficiency interval.
  - That stacks two 95 % steps. It is not a rigorous overall 95 % statement; it is an approximate, conservative
    one.
  - It is 1.031–1.049 × the joint limit on the same stored data. #546's numbers stay valid as conservative,
    approximate upper bounds; the joint values are the ones to quote.
- **Conditions of every limit here:**
  - zero survivors;
  - events independent and Poisson at one constant rate over stars and time;
  - all events in that one cell, with t0 uniform in time and u0 uniform on [0, 2);
  - at most about 1 % of events sharing a star (λ < 0.01, #546).
  - They are product-level: W3-shaped signals as they would appear in Kepler SAP of these stars. The sky-to-SAP
    response, including the archive's quality flags, is untested.
  - The 20 per-cell limits are separate statements, not simultaneous.

## 6. What the zero supports
- **Supports:** in 2,000 FGK Kepler stars (1,500 here + 500 in #542 / #546), no good cadence in this round's 1,500
  dims by ≥ 1 mag. Under the conditions above, W3 signals of any one tested shape occur in SAP at fewer than about
  4.6–6.0 per 10,000 star-years (95 % credible, flat rate prior, per cell).
- **Does not support:**
  - a physical lens rate, wormhole or negative-mass abundance, or FTL statement;
  - a limit for events shallower than 1 mag in SAP, t_E values other than the five tested (0.1, 0.3, 1, 3, 10 d), u0 ≥ 2, repeating events, other stars or
    other surveys;
  - a joint limit across cells;
  - "W3 excluded".
- **The unread batch 4:** if it were searched with similar efficiency and the result stayed zero, it would add about
  500 × the round's per-star means, 2.47–3.24 star-years each, to each cell (`combined.json` `per_star_mean_yr`).
  That is about 1,240–1,620 effective star-years (not stored), and it would lower each joint limit combined with #546 by about a fifth (factor about 0.80; each round-only limit by about a quarter, 0.75)
  (not stored).

## 7. Answers
| question | answer |
|---|---|
| 1. stars planned / read / searched / not done | 2,000 / 1,500 / 1,500 / 500 (batch 4 not fetched: the projected download would end after the 6-h limit) |
| 2. pre-screen sources / survivors / confirmed real dimming / artifacts / undetermined | 0 / 0 / 0 / 0 / 0 (no star reached the 1-mag pre-screen) |
| 3. effective search volume added | 3,710–4,855 effective star-years per W3 cell (5,195 observing, 5,994 sampling star-years); with #546, 4,959–6,475 |
| 4. new observation unexplained by ordinary astrophysics? | No. None of the 1,500 light curves has a good cadence ≥ 1 mag below its quarter median, so nothing reached the shape tests or vetting |
| 5. what the zero supports / does not | Supports per-cell, product-level, conditional limits of about 4.6–6.0 × 10⁻⁴ per star-year (95 % credible, flat rate prior, efficiency uncertainty included; with #546). Does not support a physical rate, abundance, a joint limit across cells, or an exclusion |
