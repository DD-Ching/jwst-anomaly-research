# Two-epoch NIRCam inventory refresh on all 34 cluster fields (2026-11-08)

Re-run of the frozen inventories (`results/transients/epoch_inventory/PLAN.md`,
`results/transients/epoch_inventory_ext/PLAN.md`) with unchanged rules, to confirm currency before choosing the next
D-027 two-epoch input (docs/efficiency-rules.md §4.1: one minimal query when currency must be confirmed). Metadata
only (`observed`); nothing downloaded.

- **Question:** has any NIRCam level-3 epoch within 3′ of the 34 fields been added or become public since the
  snapshots of 2026-10-26 (19 fields, MJD 61338.85) and 2026-10-27 (15 fields, MJD 61340.05)?
- **Decision rule (written before the diff):** any added or newly public epoch → plan a D-027 grid on its best pair;
  none → no new two-epoch input; next check at the first `upcoming` release.

## Result
- Queried MJD 61352.22 (2026-11-08 05:11 UTC). `core19/diff.json`, `ext15/diff.json` (`scripts/epoch_inventory_diff.py`):
  **0 / 34 fields with an added, removed or newly public epoch** → no new two-epoch input.
- New-pair counts fall only because epochs have since been named under `results/transients` ("used": 96 → 109).
  Pairs still unused on `main` (Abell 2744 60, MACS0647 1433 o010 → o020, MACS J0138 2, Abell 1689 2, MACS J0308 2)
  were all in the October snapshots; 1433 o010 → o020 and 2561 o002 → o006 were gridded on open PRs #556 / #553.
- **Upcoming** (unchanged release dates): Abell 370 9371 o001 2026-12-29 (MJD 61403.8), 7345 o101 61404.1,
  o102 61406.3, o103 61414.8, o104 61433.2, 5324 o015 61451.0, o105 61628.4; MACS J0138 8799 o003 2027-01-07
  (MJD 61412.9).

## Deviation (tool)
- `mastresolver.stsci.edu` fails TLS verification from this cloud session (curl and requests, CA bundle set), while
  `mast.stsci.edu` works. Both scripts now resolve through `resolve_centre`: astroquery first, else MAST's
  `Mast.Name.Lookup` service (here: SIMBAD-CfA). The core-19 script now queries by coordinates instead of
  `objectname`. Four ext-15 centres differ from October by 8.6–28.4″ (`ext15/diff.json`: RX J0152 28.43, RXC J2211
  19.02, SPT0546 10.11, MACS0358 8.63) with identical epoch sets; the 3′ cone and 1′ core test are therefore only
  approximately the same as in October for those four (core areas not compared).
