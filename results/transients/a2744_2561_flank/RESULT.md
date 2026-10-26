# Abell 2744 flanking pair 2561 o002 → o006: D-027 two-epoch grid (result, 2026-10-26; vetting pending)

Plan `PLAN.md` (pushed before any pixel was read; same rules and gates as #552). Tool `scripts/transient_grid.py`,
unchanged. Numbers: `summary.json`, `hits.ecsv` (92 hits, 42 groups). Imaging (observed): 2561 o002_t001 and
o006_t007 level-3 `_i2d` in F277W / F356W / F444W, 365.6 d apart (#547), outside the 1′ core; 4.5 GB streamed, deleted
after the run. Every number is derived; thresholds are the D-027 ASSUMPTIONs.

## Frozen result
- **G1 passes:** −0.75 mag injections recovered 0.918 (500; 0.989 where the fraction set the flux); −0.44 mag 0.402.
- **G2 passes:** noise scales F277W 1.45, F356W 1.37, F444W 1.29.
- 9,782 positions, 4,764 measured in ≥ 1 band, 4,465 in all three; **92 hits in 42 groups**.
- **Outcome C (pending vetting):** 41 groups have an instrumental explanation; 1 does not (below).

## Hit groups (contact sheets looked at: e1, e2 resampled onto e1, difference, three bands, 3″; not stored)
| Verdict (by eye) | Groups |
|---|---|
| diffraction spikes / star wings (roll change between epochs) | 0, 1, 3–14, 16–18, 20–26, 29, 30, 33–35, 38–41 (33) |
| bright point source with a dipole residual (registration or proper motion) | 15, 27, 28 (27 and 28 are one star) (3) |
| mosaic edge in epoch 2 (NaN region in the aperture) | 19, 31, 32, 37 (4) |
| galaxy-core dipole (registration) | 2 (1) |
| **no instrumental explanation seen** | **36** (1) |

- **Group 36** (`g09725`, RA 3.61170, Dec −30.45765): a faint compact source in a clean area (no spike or edge in
  the 3″ stamps) that brightens in all three bands: F277W 6.0σ, Δm −1.03; F356W 8.7σ, −1.32; F444W 7.6σ, −1.52
  (o002 → o006). It is an unvetted screen hit, not a candidate (trial rule T2). Next: `/vet-candidate` on it
  (persistence and detector-fixed checks on the o006 exposures, other bands and epochs, catalogue and AGN / SN
  host checks).

## Limits
Only catalogued positions; LW bands only; recovery by the automatic thresholds. `grid.ecsv` (2.5 MB) is not kept,
as on the earlier grids.
