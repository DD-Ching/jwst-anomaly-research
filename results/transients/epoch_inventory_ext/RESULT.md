# Two-epoch NIRCam inventory on 15 further lensing clusters: result

Plan `PLAN.md` (committed before any query). Tool `scripts/epoch_inventory_ext.py`, rules tested in
`tests/test_epoch_inventory_ext.py`, committed before the run. Outputs: `inventory.json` (centres, epochs, new pairs
with per-filter overlap and core areas), `upcoming.json`. Every row is MAST metadata (`observed`); areas are
`derived`. Nothing was downloaded or read as science.

## Result (frozen rules)
- 15 / 15 field names resolved. 30 NIRCam level-3 epochs within 3′, 29 public (`summary.json`); none is named
  under `results/transients` outside the two inventory directories (`used` false for all).
- **8 `core` pairs on 4 clusters** (pair footprint overlap > 1 arcmin² inside the 1′ disc in at least one shared
  filter), 2 `overlap`, 0 `none`. The other 11 fields have no new pair (fewer than two public epochs, or only
  pairs failing rule (b) or (c)).

| Field | Pair (early → late) | Δt (d) | Shared filters with core area (arcmin²) |
|---|---|---|---|
| MACS J0308.9+2645 | 5293 o005 → 6882 o015 | 200.2 | F115W, F150W, F200W, F300M, F410M: 3.142 each |
| MACS J0308.9+2645 | 5293 o005 → 9478 o002 | 230.8 | F115W 2.919 |
| MACS J0308.9+2645 | 6882 o015 → 9478 o002 | 30.6 | F115W 2.919, F444W 2.920 |
| Abell 1689 | 1840 o019 → 5782 o001 | 869.8 | F115W, F150W, F200W 2.811; F277W, F356W 2.802 |
| Abell 1689 | 5782 o001 → 6882 o048 | 342.9 | F150W 2.012, F356W 1.904 |
| Abell 1689 | 1840 o019 → 6882 o048 | 1212.8 | F150W 1.681, F356W 1.564, F444W 1.566 |
| MACS J2129.4-0741 | 4903 o004 → 6882 o040 | 392.1 | F115W 2.270, F200W 2.271, F277W 2.227, F444W 2.224 |
| MACS J0138.0-2155 | 2345 o005 → 6549 o003 | 52.0 | F150W 1.004 (F444W 0.983) |

- `overlap`: MACS J0138 6549 o001 → o003 (33.3 d, core 0.872–0.891) and 2345 o005 → 6549 o002 (36.6 d, 0.958–0.978).
- Not yet public: MACS J0138 8799 o003 (release MJD 61412.9, `upcoming.json`).
- Outcome per the plan: ≥ 1 `core` pair, so a D-027 grid on one of them is the next two-epoch input (plan first).

## Reading (post hoc, not in the plan)
- 6882 is VENUS, the same programme as the MACS J1931 grid input. Grids with a VENUS late epoch reuse that run's
  pipeline settings most directly.
- Most core coverage and filters: MACS J0308 5293 o005 → 6882 o015 (full 1′ disc in five shared filters, 200 d).
  Longest single-year baseline with LW and SW: MACS J2129 4903 o004 → 6882 o040 (392 d, four filters). Abell 1689
  5782 o001 → 6882 o048 is 343 d but shares only F150W / F356W.
- MACS J0308 9478 o002 (30.6 d after VENUS, F115W / F444W) and MACS J0138 6549 (three epochs, o001 → o003 33.3 d apart) look like
  transient follow-up visits (the MACS J0138 ones fit SN Encore / Requiem follow-up); not checked here. A grid there
  must first name the known transient as its positive control.
- None of the four clusters has a lens model validated in this repository except MACS J0138 (`encore_mu`,
  `requiem_dt` used published models). A grid is a transient search through a lensing core; a lens-model test
  needs a pinned model first.

## Limits
`s_region` is the mosaic outline: chip gaps and masks inside it count as covered, so areas are upper bounds. One
centre per field (the MAST name resolution, `centre_deg`), one core radius (ASSUMPTION 1′), 2″ grid cells. "Used"
counts only epochs named under `results/transients`; epochs used by lens-consistency runs elsewhere count as unused.
The 3′ cone bounds which epochs were seen.
