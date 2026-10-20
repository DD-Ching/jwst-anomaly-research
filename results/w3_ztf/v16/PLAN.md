# W3 ZTF v16 — inventory of up to 30 partitions, then a frozen queue of up to 8 with the v13 method

Plan fixed 2026-10-19, before the inventory results and before any light curve of this round. Owner authorization:
the round after #383 (4 h, ≤ 6 GB, ≤ 8 partitions, up to 20,000 new groups). Frozen and unchanged: source cuts,
quality flags, `--edge time` (EDGE_FRAC 0.25, EDGE_NSIG 2), the 10 cells, all thresholds. No new classifier, no
limit, no wormhole statement.

## 1. Inventory (metadata only: DR24 objects, IRSA TAP, parquet metadata; `inventory_script.py.txt`)
- **Candidates** (`candidates.json`, HATS partition lists only):
  - Norder 5; the parent has objects; the parent was not used before (959, 1783, 595, 2042, 45);
  - 300k–700k rows; centre dec +25° to +65°, |b| ≥ 30°;
  - 8 parents evenly spaced in Npix order, up to 4 children each. That gives **25 candidates**.
- **Per candidate:**
  - **unused eligible groups:** all groups passing the unchanged `select()` cuts, minus any within 3″ of, or sharing
    an objectid with, the 20,153 earlier sample groups (`dedupe.json.gz`). Earlier reference stars lie in other
    partitions;
  - **searchable if 500 go to calibration;** reference (control) stars are other objects and do not reduce this;
  - the v11/v12 necessary-sampling fraction for the 10 cells (200 positions, 2,000 draws per cell);
  - the product ratio (DR24 `ngoodobs` / counted exposures);
  - ZTF fields shared with earlier samples;
  - the light-curve partition size (7 columns, from parquet metadata).
- **Three counts kept apart:** catalogue-eligible unused; sampling-qualified; calibration passed (only after
  section 3). The first two are not reliable search counts. The inventory covers these 25 partitions only; it is not
  an all-sky census.

## 2. Frozen queue (`queue_script.py.txt` → `queue.json`, committed before any light curve)
- **Qualifies:** worst-cell fraction ≥ 0.5, product ratio ≥ 0.8, ≥ 1,000 unused eligible groups.
- **Order:** unused eligible groups per GB of partition, descending. Ties (to 0.1) put areas with no shared field
  first, then lower Npix. At most 8 areas.
- **Execution** (`driver_script.py.txt`): two workers take areas strictly in queue order. A failed area is recorded
  and the next queued area starts; nothing is added and the order does not change.
- **Stop rules for starting a new area:**
  - the projected download (objects + TAP + each started area's partition + 30 MB TAP) would pass 5,900 MB;
  - or the time is after the no-new-area limit (start + 3 h 05 min).
  - Running areas finish.

## 3. Per area (frozen v13 method)
- **Sample** (`area_sample_script.py.txt`):
  - every unused eligible group;
  - the lowest 500 by crc32("v16:" + key) are calibration, the rest search;
  - controls are the comparison range minus all sample members;
  - groups bind all objectids and filters of a source;
  - one `configure()` per process; each area has its own folder `v16/<npix>` and cache `data/cache/w3_ztf/<npix>`.
- **Validation:**
  - `inject --cells calibrated --n 200 --n-fade 1000 --edge time --seed-offset 10000 + 1000 × (position − 1)`;
  - `search --half dev --edge time`.
  - **Pass:** every cell ≥ 0.5 (all injections counted), every control ≤ 10 %, ≤ 5 calibration survivors
    (`decision_script.py.txt`).
- **Search:** a passing area is searched at once (`search --half holdout --edge time`). On failure the area is
  recorded and the search group stays unread. No new seed, no tuning.

## 4. Image check of final (vetted) survivors only
- **Limits:** ≤ 10 groups, ≤ 100 cutouts, ≤ 50 MB over the whole batch. **Order:** queue position, then key.
  Method: `../v14/image_script.py.txt` (survivors only), plus official masks when the cap allows.
- **Questions per survivor:**
  - coverage, mask and saturation at the source;
  - same-exposure comparisons and the local image depth against the source's usual brightness;
  - the flux before / during / after.
- **Classes:** near-limit or poor images are unknown, never dimming. Anything beyond the limits stays pending.
- **If an image-confirmed dimming without an ordinary explanation appears:** no new area is started, and the rest of
  the budget goes to it.
- #380 / #383 pending items stay as they are.

## 5. Records
- `timeline.jsonl` holds UTC start and end and seconds for every stage, plus download projections. Rates per 1,000
  searched groups are measured from it; the #383 1.9 min figure stays an estimate.
- At least 30 min are kept for review and records.
