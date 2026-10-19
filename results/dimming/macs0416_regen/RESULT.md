# W3 dimming screen: MACS0416 regenerated; joint three-field `combine` (D-052 addendum 2026-10-19d)

Plan frozen first (`PLAN.md`). Code and thresholds unchanged (`Params` defaults): `screen`, `forced --n-controls 300
--workers 16` (410 s wall), `inject`; then `combine --field macs0416,abell2744,nexus` over this `limits.ecsv` and the
committed `abell2744_calib/limits.ecsv` and `nexus_calib/limits.ecsv`. Outputs here: `screen_summary.json`,
`forced_summary.json`, `inject_summary.json`, `limits.ecsv`, `repro_gate.json`, `limits_combined.ecsv`,
`limits_combined_all_fields.ecsv` (derived; injections simulated).

## Data route (deviation)
`fetch` failed: the MAST API (`mast.stsci.edu/api/v0/invoke`) timed out or dropped the proxy connection three times
(13:15–13:33 UTC). The 9 catalogues were read instead from their S3 `cloud_uri` (listed in
`data/manifests/dimming_macs0416_products.ecsv`); all 9 match the committed `data/manifests/dimming_macs0416.ecsv`
in size and sha256 (not stored beyond this statement; the check compares against the committed manifest).

## Result
- **Screen (derived):** reproduces D-052 exactly: 7,708 sources, 44 monitored compact, 854 flags, 358 after the
  catalogue vetoes, 10 compact (`screen_summary.json`).
- **Forced stage: calibrated** (`forced_summary.json`): 358 flags + 300 controls; 1 forced-confirmed (rise-dip-rise),
  0 cutout-clean, 0 survive the artefact tests; 0 / 300 controls flagged. The confirmed flag is `d03349`, the
  saturated star D-052 already named (spikes visible in `forced_sheet.png`, looked at; not committed).
- **Reproduction gate (`repro_gate.json`, decided by `scripts/dimming_repro_gate.py`): passed.** Regenerated MACS0416
  exposures differ from the D-052 values by −3.7 % to +0.7 % (limit ±10 %, ASSUMPTION).
- The forced noise scale (3.11 / 4.57, F200W / F444W) and the F444W inflation (2.62) differ from the D-052 values
  (2.15 / 1.99 and 1.98), whose forced fluxes were re-calibrated offline after an S3 failure (docs/exotic_limits.md);
  the cause of the difference was not investigated. With the gate passed, the exposure moves by ≤ 3.7 %.
- **Joint headline, three calibrated fields** (`limits_combined.ecsv`; ρ = 0.1; 0 events, Poisson 95 % = 3.0;
  235 monitored compact sources):

  | t_E (yr) | exposure (source yr) | per source per yr | τ | per deg² per yr | per deg² per epoch |
  |---|---|---|---|---|---|
  | 0.01 | 6.85 | 0.44 | 0.014 | 3.0e3 | 94 |
  | 0.03 | 26.8 | 0.11 | 0.011 | 776 | 73 |
  | 0.1 | 69.1 | 0.043 | 0.014 | 295 | 93 |
  | 0.3 | 103.5 | 0.029 | 0.027 | 201 | 190 |
  | 1 | 145.5 | 0.021 | 0.065 | 157 | 493 |
  | 3 | 191.2 | 0.016 | 0.15 | 118 | 1.1e3 |

  The per-source column agrees with the summed headline of D-052 addendum b (0.43 / 0.11 / 0.043 / 0.029 / 0.020 /
  0.016) within the MACS0416 regeneration difference. The per-deg² columns are new: 4.0–12× tighter than the
  MACS0416-only ones (2.9e4 / 7.2e3 / 3.5e3 / 1.8e3 / 635 / 545 per deg² per yr), at these fields' compact-source
  densities.
- `limits_combined_all_fields.ecsv` equals `limits_combined.ecsv` (all three fields are now calibrated); the older
  "including the two uncalibrated ones" table in docs/exotic_limits.md is superseded.

## Limits of this result
- Same caveats as D-052: isolated lens, no macro-magnification; catalogue-level injection; the `single_epoch` veto
  blinds full vanishes of two-epoch sources.
- `chromatic_veto` (#371) is still not in `dimming_screen.py forced`; these limits do not use it.

## Next
- Fold `chromatic_veto` into `forced` and re-run the injections of all three fields; replace the `single_epoch`
  veto for vanish flags with the D-039 persistence test; add JADES and new NEXUS epochs.
