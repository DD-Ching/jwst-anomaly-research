# VENUS batch 2: the three "real change" entries (vetting record)

Run: `results/transients/batch_venus2` (#304). Positions and photometry: `boxes/review/review.json`,
`boxes/review/lightcurves.ecsv` (0.15″ forced apertures, random-aperture errors; host light adds a constant offset per
band, so compare within a band). Checks: `boxes/review/post_review_checks.json`, `boxes/review/final_classes.json`,
`boxes/review/crossmatch.json`. Sheets and stamps looked at (`boxes/review/sheet_*`, `stamps_*`). Vetted 2026-10-17.
Every number is `derived` unless labelled. An ordinary transient is not wormhole evidence; none of these is an
exotic candidate.

## Common results
- **Instrumental tests** (all three): inside all four footprints, finite fraction ≥ 0.9999 (`results/transients/batch_venus2/boxes.json`);
  not within 1″ of the box edge; F150W and F444W S peaks ≤ 0.06″ apart (sky-fixed); a point source in ≥ 3 bands of
  the same visit (not a cosmic ray or single-band artifact); sheets show no spike through the position.
- **Catalogues:** no Gaia DR3 source within 3″ (not a star); no reported transient found within 1″ in the sources
  searched. **TNS could not be queried** (HTTP 403), so "no TNS record" is not established. Finding none does not
  mean the public record was exhausted.

## MACSJ1423+2404_2:h000 (215.949099, +24.086948): fading
- 1208 (2023-02-26): a point source on a galaxy's outskirts in all 8 bands (e.g. F150W 551.9 ± 4.4 nJy including host
  light). VENUS (2026-03-05): F150W 275.1 ± 8.4 nJy, z −29.1; F444W 268.9 → 178.1 nJy, z −12.5; no point source on the
  stamps, i.e. **not detected as a point source in later observations**.
- Host: GLASS MACS1423-00487 (SIMBAD, 0.72″, z 0.525; NED: SSTSL2 J142347.77+240514.0, 1.08″, z 0.53, redshift type
  not determined). Catalogue position near the source does not by itself establish the host.
- **Status:** real change; an ordinary transient (supernova-like) is a reasonable hypothesis; identity and type not
  confirmed (no spectrum, one detected epoch).

## MJ0417.5-1154_2:h002 (64.396407, −11.907040): brightening
- 1208 (2022-10-13): no point source. VENUS (2025-10-06): a point source in F150W, F210M, F300M and F444W; F150W
  299.0 ± 26.3 → 867.3 ± 32.3 nJy (z +13.7, host light included), F444W 289.0 → 532.5 nJy (z +3.2).
- Position: on image 1.1 of the multiply-imaged lensed system 1 (z = 0.8718; `[MSF2019] MACS J0417 1c.1` 0.26″,
  Mahler+2019 / Jauzac+2019 via SIMBAD/NED). Its counter-images 1.2 lie 6–7″ away.
- **Status:** real change; an ordinary transient in a lensed background galaxy is a reasonable hypothesis; identity
  and type not confirmed. If it is a transient in system 1, the counter-images should show it with lensing time
  delays: **not checked** (needs a lens model's image positions for that knot and delays; not done in this budget).

## BULLET_2:h019 (104.655304, −55.945271): brightening (faint)
- 4598 (2025-01-20): not detected in 8 bands (F150W −4.1 ± 2.5, F444W 0.7 ± 1.7 nJy). VENUS (2026-01-21): F150W
  38.0 ± 7.0, F210M 44.7 ± 7.7, F300M 23.5 ± 4.7, F444W 25.2 ± 3.6 nJy (S/N 5–7).
- No catalogued object within 3″; nearest lensed image `[PKR2016] E.3` (z ≈ 2.96, likely photometric) at 3.7″ (weak
  link). A horizontal streak 0.5″ away is present in both epochs at the same orientation (a fixed sky feature).
- The VENUS–4598 brown-dwarf paper (arXiv:2604.23668) reports proper motions of 24–49 mas/yr; such a mover would stay
  inside the 0.15″ aperture over one year, and nothing is detected in 4598, so a known moving dwarf does not explain it
  (positions in that paper not checked).
- **Status:** real change at modest S/N; an ordinary transient is plausible; identity and type unknown.

## Not explained / open
- None of the three is "unexplained": each has an ordinary transient hypothesis that the data do not contradict.
  Their identities and types are open (no spectra, single detected epochs, TNS not queried).
