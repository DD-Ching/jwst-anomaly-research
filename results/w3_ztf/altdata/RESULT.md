# W3 ZTF × alternative public data — ZTF g restores 3 of 63 insufficient W3 events (needs ≥ 50 %); fade passes stay 0 / 79; NEOWISE and TESS do not reach the gaps; not continued

- **Pre-registration:**
  - PLAN sections 0–3 (gap windows, inventory, selection rule) were committed at 05:32:14 (780b528), before the
    inventory output was read.
  - **Deviation:** the inventory queries had already started at 05:30:34 (`cost.json`). The rule text was written
    while they ran, so it is pre-registered against the output, not against the start of the queries.
  - Sections 3b–6 (selection outcome, trial, criteria, stop) and `edge_bounds_alt` with its unit test were committed at
    05:45:24 (7b7554d). That was after the inventory (coverage only) and before the data re-read and any flux.
- **Labels:**
  - ZTF DR24 light curves, IRSA exposure metadata, NEOWISE-R epochs and flags, and TESS product time ranges are
    `observed`.
  - Injected events are `simulated`. Counts are `derived`.
  - The diagnostic in section 4 uses injection truth after classification and is **post hoc**.
- **Scope:**
  - A feasibility trial on an already-used calibration sample (the 40 stars of `../asassn/`, the v23 events). It is
    **not** an independent validation.
  - No search group was read (`../v23/dev*/source_usage.jsonl`: `search_member_rows` 0).
  - ASAS-SN, the ZTF r-band psfcat route, X4 and the other paused units were not touched.

## 1. Inventory (`inventory.json`; metadata and epoch times only)
- **Gap windows:** 106 bounds-insufficient sides on the 40 stars: 88 sides of 63 W3 events and 18 sides of 12 fades.
  The windows span MJD 58359.4–60821.7 (2018-08-29 to 2025-05-26).

| source (access) | band, resolution | W3 sides with a usable epoch in the window | fade sides | flux, errors, flags | status |
|---|---|---|---|---|---|
| ZTF g (IRSA IBE metadata; DR24 light curves in the repo's reader) | g, same camera as r | 36 / 88 (all W3 sides of 19 / 63 events) | 8 / 18 | PSF mag + error, infobits, catflags | **observed** (infobits-0 exposure with the star ≥ 10″ inside) |
| NEOWISE-R single exposures (IRSA SCS, 4″) | 3.4 / 4.6 µm, ~6″ | 3 / 88 (1 event) | 0 / 18 | profile-fit mag + error, frame and source flags | observed in 3 sides; epochs end MJD 60512.6 (2024-07-21) |
| TESS (MAST CAOM: FFIs, QLP / TESS-SPOC / TARS products) | 600–1000 nm, 21″ pixels | 0 / 88 | 0 / 18 | all 40 stars fall on FFIs; light-curve products exist for 6 stars, all in 2344 | **confirmed no coverage**: no sector overlaps any window; the pixels also fail the 16″ resolution condition |

- **Not queried** (sources in SOURCES.md):

| source | reason | coverage of the windows |
|---|---|---|
| Gaia DR3 epoch photometry | the data end 2017-05-28, before every window | confirmed none |
| Pan-STARRS1 3π | 2010–2014 | confirmed none |
| ASAS-SN Sky Patrol | excluded by the brief (paused, `../asassn/`) | unknown |
| ATLAS forced photometry | needs a registered account (brief: no new account) | unknown |
| Gaia DR4 epoch photometry | not released at 2026-10-22 | unknown |

- **NEOWISE quality cuts:** fixed in the script, and modelled on the Explanatory Supplement's single-exposure advice
  (not re-checked against its text this round). They decide only 3 sides, so they do not change the choice.

## 2. Selection (PLAN §3)
ZTF g has the most "observed" W3 sides (36 > 3 > 0). It has 19 events whose insufficient sides are all observed, so
the stop rule (≥ 1) did not trigger. All 40 stars are used: every one has a g series (33–2,002 g detections per star,
`result.json`), and none is excluded.

## 3. Paired comparison on the same events (`trial_<npix>.json`, `trial_events_<npix>.json.gz`, `result.json`)
- **Check:** A (bounds on ZTF r alone, replayed) equals the stored v24 outcome for all 3,000 events of each area,
  under both rules (asserted per event; the runs completed). Fetched rows were identical to v23 (13,497 and 9,404).
- **The post hoc diagnostic** (section 4) saw the same A-insufficient sides as the trial and matched each to an
  inventory side (`result.json`: `diag_consistency_checks`).
- **New observations:**
  - The 40 stars carry 40,584 g detections.
  - In the 88 W3 windows, 22 sides hold at least one g detection (40 detections in total) and 5 sides a normal one
    (dm < 0.3, S/N ≥ 3).
  - All 36 sides with an infobits-0 g exposure in the window hold a DR24 g epoch there (detection or non-detection;
    post hoc, `diag_posthoc_<npix>.json`).

| events on the 40 stars | A: ZTF r | B: ZTF r + g |
|---|---|---|
| W3 bounds-insufficient | 63 | 60 |
| W3 recovered | 48 | 51 |
| W3 recovered in A, lost in B | — | 0 |
| W3 sides decided (of the 88): fast / slow | — | 5 / 0 |
| fade passes (of 79 fade injections) | 0 | 0 |
| fade insufficient | 12 | 10 (2 became slow) |

- **Restored W3 events:** 3, all in 4991 (cells ρ 0.01 / blend 1 / t_E 100 d, ρ 0.1 / 1 / 10 d and ρ 0.1 / 0.7 / 10 d),
  on stars 601213400002567, 601213200006434 and 601213200002823. Two more sides became fast in events that keep another
  undecided side.
- **Decision (PLAN §5, computed in `result_script.py.txt`):**
  - restored / A-insufficient = 3 / 63 = 0.048, below the required 0.50;
  - B fade passes / fade injections = 0 / 79 = 0.0, within the limit of 0.10;
  - **success: false.**
- **Coverage ceiling (inventory):** only 19 of the 63 events have every insufficient side covered by an infobits-0 g
  exposure. So even if each such exposure gave a normal point, at most 19 / 63 = 0.30 could be restored, still below
  0.50.
- **Cost of B:** no W3 recovery was lost, and no fade passes. B made 2 fades decidable as slow (correct rejections),
  through significant g intermediates beyond 0.25 × span.

## 4. Why so few windows are bounded (POST HOC, `diag_posthoc_<npix>.json`; injection truth after classification)
- **What the g epochs show:** the 88 W3 windows hold 66 g epochs.
  - By observed state: 26 non-detections, 28 dim detections (dm ≥ 1), 5 significant intermediates and 7 normal.
  - By truth: at 59 of the 66 the injected model is ≥ 0.3 mag below baseline. These epochs fall inside the event,
    not after its edge.
- **By side class** (W3 sides with a g epoch / with a normal g point):
  - A (no further r exposure): 26 / 2 of 52;
  - B (quality-flagged r): 6 / 2 of 29;
  - C1 (DR24 r gap): 4 / 1 of 7.
- **Reading:** the r gap next to a run is where the true edge lies. The g epochs there mostly see the event still in
  progress. Under bounds, the run and its span come from r alone, and a dim g point outside the r run bounds nothing.
  Letting other-band dim points extend the run would be a method change. It is not made here.

## 5. Deviations from PLAN (labelled)
- **Inventory rerun:** the first run stopped at 05:42, because 4 of 40 NEOWISE responses were IRSA server errors
  (ORA-24459 pool timeouts; count hand-logged in `cost.json`).
  - The four error files were deleted by hand. Retry-on-error was added, and the inventory rerun reusing the valid
    cached responses.
  - After review, the script also refetches a cached error response itself.
  - This is not a method change.
- **Review refactor (rerun 06:02–06:04, `cost.json`):** `edge_bounds` now calls `edge_bounds_alt` without alt points, so there is one
  implementation. The trial was rerun once with it:
  - both areas passed the A = v24 assert;
  - the event records were identical (decompressed contents equal; only the run time differed; check not stored);
  - the committed trial files are those of the first run.
- **PLAN §2 erratum:** PLAN gives Pan-STARRS1 3π as 2009–2014. The survey ran 2010–2014 (SOURCES.md). The exclusion
  is unchanged.
- **Timing of §3b–§4:** these sections were written after the inventory output had been read (coverage only), as
  PLAN §3b states.
- **Post hoc diagnostic:** section 4 was not in the plan.
- **v23 files:** `fetch_info.json` was restored after the re-read (only `read_s` differed). `source_usage.jsonl` keeps
  the appended usage lines.

## 6. Limits of the rule as applied
- **σ is small:** σ for g (and r) points comes from the maglimit model (S/N 5 at the limit). For stars several
  magnitudes above the limit, as here, σ is of order 10⁻³ mag (typical value not stored).
  - So the S/N ≥ 3 condition and the 2σ margin do not bind, and "significant intermediate" is in effect
    0.3 ≤ dm < 1.0.
  - This is the same as for r epochs in bounds (v23–v24). It is kept, not changed post hoc.
- **Dim alt points:** dim alt detections (dm ≥ 1) outside the run are not used. A normal g point followed by a dim one
  still bounds the edge, as for r epochs in bounds.

## 7. Cost (`cost.json`)
- About 1,246 MB downloaded, under the 2 GB budget:
  - the two calibration re-reads, the same files as v24 (1,038.6 MB light curves, 168.4 MB objects, 14.3 MB TAP);
  - 25 MB of inventory responses (MAST queries not measured).
- About 31 min from the start (05:22) to the end of the diagnostic (05:53).

## 8. Answers
- **New information that changes a decision:** yes, but little. ZTF g decides 5 of 88 undecided W3 sides (3 events
  restored) and 2 fades (rejected as slow). NEOWISE reaches 3 sides; TESS reaches none.
- **Where it applies, and its cost:**
  - It applies to these 40 bright stars (ZTF r 14.5–14.8), 3 of 63 insufficient W3 events.
  - It cost no lost W3 recovery and no fade pass; 1.25 GB of re-reads and no new photometry.
  - The pre-registered success condition fails (0.048 < 0.50), and the inventory ceiling (0.30) is also below it.
- **Next step:** stop this route; ZTF g is not continued.
  - No source checked here covers the edge moments well enough.
  - Reopen only with a source of nightly cadence that is reachable without an account, or with an owner decision on a
    rule in which other-band dim points extend the run. Such a rule needs a new plan and validation.
  - The ASAS-SN plan stays paused and unchanged.
