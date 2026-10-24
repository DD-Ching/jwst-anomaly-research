# Bullet Cluster CANUCS extra-image test: status (incomplete, 2026-10-24)

Plan `PLAN.md` (frozen in 2a7e851; script and CANUCS sha256 pins committed in 1ae40dd, both before the run).
**The frozen run has not finished; no gate C2 outcome and no flag exist.**

## What ran
- Rules 1–2: the run went past rule 1, so frame check F passed (by the script's control flow; its numbers were held
  only in the unwritten top-level summary: **not stored**). Rule 2 repeated gate K-WLS: rms 0.391″, 135 / 135 within
  3″ (`summary_wls.json`, identical to `../bullet_gatek/summary_wls.json`).
- Rule 3 (frame tables, WLS source, GO 4598 `o001_t001`): `frame/summary.json`, `frame/predicted_bullet.ecsv`,
  `frame/control_bullet.ecsv`. Reported, not used (PLAN.md rule 3): pipeline forced-test gate C 58 / 60 testable;
  image classes 136 observed, 6 candidate, 7 demagnified, 7 faint, 8 no_flux_ref, 1 missing; 0 pipeline flags with
  the gate ignored.
- Rule 5 stopped at its input check: the CANUCS F200W `bgsub-sci` mosaic (sha256 pinned) fails the frozen 0.04″
  pixel-scale check. astropy reads its WCS as 3600″ per pixel: the header carries both `CD` keywords (|CD1_1| = 2.22e-5°,
  i.e. 0.08″, and `PIXSCALE` = 0.08) and `PC` keywords without `CDELT`. Session check, not stored: the photometry
  catalogue's `X`, `Y` against `RA`, `DEC` give 0.0400″ per pixel for 9 source pairs, the image (10240 × 8240) matches
  that grid's extent, and a `CD`-only WCS misplaces catalogue sources by hundreds of pixels. So neither header reading is
  the catalogue's 0.04″ grid.

## Failed approach (operations)
The first attempt ran as a background shell job with a 10-minute limit and was killed while building the frame tables; the
second (same plan, unchanged) took about 16 min (file times, not stored) for rules 1–3 (0.1″ deflection grid of half-width 219″, then `find_images`
per family and S3 stamps of a 1.9 GB mosaic). Run it detached (`nohup`) from the start.

## Next
A frozen deviation note (`PLAN_1b.md`) before any CANUCS flux or pixel is read: derive the mosaic WCS from the catalogue
(e.g. a TAN WCS fitted to catalogue `X`, `Y` → `RA`, `DEC`, with a stated residual gate) or use the CANUCS F200W
`sci` / `rms` files if their headers pass the frozen check; then rules 4–8 unchanged on the stored frame tables.
