# W3 ZTF × alternative public data — can another source supply the missing edge information for the 40 stars? (one feasibility trial)

Owner brief 2026-10-22 (session start 05:22 UTC): 2 h, ≤ 2 GB; skip ASAS-SN maintenance and retries; at most 3 directly
accessible public sources compared, at most 1 tested; no account, no paid resource, no proxy workaround, no outside
contact. X4 and the other paused units stay as they are. A feasibility trial on an already-used calibration sample:
**not** an independent validation; no search group is read. The question is only: can other public data supply the
missing edge moments of these stars' events, so that some events become decidable?

## 0. Prior state (checked 05:22–05:30)
- `main` at 51bf60a; no open PR or `claude/*` branch claims W3 / ZTF or the 40 stars (open: #458 MACS0416 claim,
  dependabot #442, `needs-human` backlog). X4 untouched.
- Reused, not redone: the 40 stars (`../asassn/selection.json`), their ZTF-only baseline (`../asassn/ztf_baseline.json`),
  the v23 events and seeds and the v24 per-event records (`../v24/step12_events_*.json.gz`), the bounds rule
  (`scripts/w3_ztf.py: edge_bounds`). v24 closed the ZTF r-band psfcat route; ASAS-SN is paused (light-curve service
  down 2026-10-20; blocked by this environment 2026-10-22).
- The repository has never examined TESS, NEOWISE or ZTF g-band exposures as edge information. In the W3 ZTF chain the
  g band enters only the achromatic test (epochs inside the run), never the edges.

## 1. Gap windows (no new data)
From the v24 records, every side of a bounds-insufficient event on the 40 stars (run passing the other tests, side class
"insufficient"): 106 sides, 88 of 63 W3 events and 18 of 12 fades. Window: pre [edge_t − 0.25 span, edge_t), post
(edge_t, edge_t + 0.25 span]; a bounding normal point must fall there.

## 2. Inventory (`inventory_script.py.txt` → `inventory.json`; metadata and epoch times only, no target flux)
Sources compared (all official services, no account):
1. **ZTF g band** (IRSA IBE exposure metadata; the repo's footprint rule `in_quad`, 10″): observed = an infobits-0 g
   exposure with the star inside in the window; possible = only flagged exposures; none.
2. **TESS** (MAST CAOM): possible = a timeseries product's time range overlaps; ffi_only = on an FFI of an overlapping
   sector without a photometry product; none.
3. **NEOWISE-R single exposures** (IRSA SCS, 4″): observed = an epoch in the window passing the fixed quality cuts in
   the script; possible = only failing epochs; none.

Excluded without a query (sources in `exclusions` of RESULT): Gaia DR3 epoch photometry (ends 2017-05-28, before every
window), Pan-STARRS1 3π (2009–2014), ASAS-SN (brief), ATLAS forced photometry (needs an account).

## 3. Selection rule (fixed before the inventory output was read)
- **Eligible:** per-epoch flux or magnitude with errors and quality flags from an existing product (no own photometry),
  angular resolution no coarser than the 16″ clean radius of the selection, no account.
- **Choose** the eligible source with the most W3 sides "observed". It must have at least one W3 event whose
  insufficient sides are all "observed"; otherwise stop and deliver the exclusion table.
- **Tie:** the repository's existing interface first.
- Light-curve shapes of the targets are not looked at for the choice.

## 3b. Inventory outcome (05:44 UTC; coverage only, `inventory.json`)
- W3 sides "observed": ZTF g 36 / 88 (19 of 63 events with every insufficient side observed), NEOWISE 3 / 88 (1 event),
  TESS 0 / 88 (no sector overlaps any window; its 21″ pixels also fail the resolution condition).
- **Selected: ZTF g.** "Observed" here means an infobits-0 exposure exists; whether DR24 holds a usable g detection
  there is measured in the trial.

## 4. Trial with the selected source (fixed before any flux is read; written for ZTF g, the source the rule selects)
- **Sample:** the same 40 stars. Every event of the v23 replay whose carrier is one of them: the 10 calibrated W3 cells
  × 200 and the 1,000 fades per area, with carriers drawn by the v23 seeds (offsets 45000 for 4991, 46000 for 2344). No
  re-selection. A star without g detections adds nothing (B = A for its events) and stays in every denominator.
- **Data:** the ZTF DR24 g-band series of the same groups, as built by `series(g, rows, cov, 1)`: good exposures, star in
  the footprint, flagged epochs unknown, non-detections inferred from coverage. These are the files of v23 / v24 (two
  calibration re-reads, about 1.22 GB). No psfcat, no other product, no new photometry.
- **Injection:** the same physical event (same dm(t), same blend) on the g times, as in every version since v2. The g
  noise is the draw the v23 replay already makes: same generator, after the r draw, one independent normal draw per g
  epoch with σ from that epoch's maglimit (SNR 5 at the limit; `apply_dmag`). The r series and its noise stay
  identical. No new random number is drawn.
- **Measurements:** only g detections enter the edges. A g non-detection is neither a measurement nor zero flux: it is
  left out of the edges and keeps its existing role in the achromatic test. DR24 has no upper limits; none is used. A
  missing epoch is never zero.
- **Cross-band handling (fixed now):** W3 lensing is achromatic and the model applies the same dm(t) with the cell's
  blend in both bands. A g point's dm is relative to the g median of its quadrant (the series' own baseline). Normal:
  dm < 0.3 and S/N = 1.0857 / σ ≥ 3 (the ASAS-SN plan's condition). Significant intermediate: 0.3 + 2σ ≤ dm < 1.0. A
  colour change of an ordinary fade is not modelled (the fade control applies the same dm in both bands).
- **Rules:**
  - **A** is bounds on ZTF r alone, replayed. It must equal the stored v24 outcome for every event of both areas under
    both rules (v13 and bounds), else stop.
  - **B** is bounds with the g detections as added edge points: `edge_bounds_alt` in `scripts/w3_ztf.py` (unit test
    `test_altdata_edge_bounds_alt`, written before the run). Run finding, span, EDGE_FRAC and the other three tests are
    taken from A unchanged.
- **Script:** `trial_script.py.txt` → `trial_<npix>.json`, `trial_events_<npix>.json.gz`.

## 5. Criteria (from `../asassn/PLAN.md` §3, unchanged)
- **Success (development):** restored = W3 events on the 40 stars that are bounds "insufficient" in A and "recovered" in
  B. Success if restored / A-insufficient ≥ 0.50 **and** B fade passes / fade injections on the 40 stars ≤ 0.10. Both
  areas combined, the 10 calibrated cells pooled. Per-area numbers are reported, not gated.
- **Reported:** g detections per star; g points inside the A-insufficient windows; sides with a normal g point there;
  sides turning fast or slow; A-recovered W3 events that B loses; fade passes in A and B; the A → B transition table.
- **Scope:** a feasibility trial on the already-used calibration sample (development data), not an independent
  validation. It applies to these 40 bright stars (ZTF r 14.5–14.8, clean 16″ apertures) and these simulated events
  only. No search group is read.

## 6. Stop and cost
- Stop computing at 07:00 UTC or at 2 GB downloaded in total. If A differs from the stored outcomes, stop and report.
- One run per area. A crash without a verdict may be rerun once with the same inputs. No threshold or rule change.
