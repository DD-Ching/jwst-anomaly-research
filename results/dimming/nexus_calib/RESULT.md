# W3 dimming screen, NEXUS: forced stage calibrated (D-052 addendum 2026-10-19b)

Plan frozen first (`PLAN.md`). Thresholds unchanged (`Params` defaults); `fetch`, `screen`, `forced --n-controls 300
--workers 16` (all 1,051 flags; 1,756 s wall), `inject`. Outputs here: `forced_summary.json`, `inject_summary.json`,
`limits.ecsv`, `headline.json` (derived; injections simulated), `survivors_and_flagged_controls.ecsv`,
`survivors_sheet.png`.

## Speed option validated first
`forced --workers N` reads the (band, epoch) images concurrently; nothing after the reads changed. On Abell 2744
`--workers 6` reproduces `abell2744_calib/forced_summary.json` exactly (all fields except wall time; 666 → 236 s;
`workers_validation.json`). Gate passed, so the option was used.

## Result
- **Screen (derived):** reproduces D-052: 3,793 catalogue flags, 1,051 after the catalogue vetoes.
- **Forced stage: calibrated** (D-052 gate). Noise scale 1.91 (F200W) / 2.57 (F444W); forced / catalogue
  inflation 1.00 / 2.06 — close to MACS0416's (2.15 / 1.99; 1.00 / 1.98), unlike Abell 2744's.
- 1,051 flags measured: 32 forced-confirmed (11 vanish, 22 rise-dip-rise), 5 survive the cutout and catalogue
  artefact tests. 10 of 300 controls are flagged in forced photometry (all vanish-type).
- **The 5 survivors, by eye (post hoc; the plan sends them to `/vet-candidate`, which has not run — deviation):**
  - d53823 (e0), d60156 (e3), d79928 (e6): one sharp F200W point in a single epoch (4.3, 10.4, 24.4 µJy), with
    nothing in the simultaneous F444W exposure (≤ 0.05 µJy). A single-band, single-epoch point is the signature of a
    detector or cosmic-ray residual; not a W3 dimming.
  - d73390: nucleus of a two-nucleus galaxy (NED: WISEA J175433.33+651306.6); F200W 0.43–0.86 µJy across epochs
    while F444W stays 3.49–3.92 µJy. Varies in one band only; the 0.15″ aperture sits on an extended nucleus. Not a
    W3 dimming.
  - d79956: a point source on a faint galaxy that brightens in e5–e7 in both bands (F200W 0.02 → 0.16 µJy, F444W
    0.01 → 0.056 µJy). It looks like an ordinary transient such as a supernova (not checked against TNS); a
    brightening, not a W3 dimming.
- **0 W3 events among the survivors, by eye (post hoc).**
- **Injection (simulated):** W3 efficiency 0.026–0.218 at ρ = 0.1 (0.218 at t_E = 0.3 yr), against 0.07–0.24
  uncalibrated in D-052.
- **Headline, three fields** (`headline.json`; ρ = 0.1; exposures add; 95 % Poisson for zero events):

  | t_E (yr) | NEXUS (source yr) | MACS0416 + Abell 2744 | per source per yr (was) | tighter by |
  |---|---|---|---|---|
  | 0.01 | 5.60 | 1.29 | 0.43 (2.3) | ×5.3 |
  | 0.03 | 21.47 | 5.26 | 0.11 (0.57) | ×5.1 |
  | 0.1 | 58.04 | 11.08 | 0.043 (0.27) | ×6.2 |
  | 0.3 | 82.82 | 21.05 | 0.029 (0.14) | ×4.9 |
  | 1 | 89.99 | 56.57 | 0.020 (0.053) | ×2.6 |
  | 3 | 123.19 | 68.78 | 0.016 (0.044) | ×2.8 |

## Limits of this result
- The headline holds only if none of the 5 survivors is a W3 event; that rests on the by-eye reading above until
  `/vet-candidate` runs (cheap: the cutouts are on S3).
- The 10 / 300 flagged controls (3.3 %, all vanish-type; Abell 2744 0 / 300) were not examined; a vanish flag in
  controls is not "confirmed" because controls carry no catalogue flag.
- Same caveats as D-052 (isolated lens, no macro-magnification; catalogue-level injection; the `single_epoch` veto).
- `combine` was not run (MACS0416's `limits.ecsv` must be regenerated first, ≈ 32 min sequential, ≈ 10 min with
  `--workers`); the per-deg² headline is not recomputed.

## Next
- `/vet-candidate` on d79956 (TNS / SIMBAD at the position) and the three single-band points (exposure-level
  cosmic-ray / snowball check); look at the 10 flagged controls.
- Regenerate MACS0416 `forced` with `--workers`, then `combine` all three fields; update docs/exotic_limits.md.
