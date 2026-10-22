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
