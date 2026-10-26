# g09725 short-wavelength difference, 2561 o002 → o006: result (2026-10-26)

Plan `PLAN.md` (pushed in 7e8d668f before any SW pixel). Tool `scripts/g09725_sw.py`, tests
`tests/test_g09725_sw.py`. Numbers from `summary.json`; stamps in `contact_sheet.png`. Fluxes are observed
(difference images, µJy in a 0.20″ aperture, not aperture-corrected); thresholds are ASSUMPTIONs.

## Frozen result: outcome U
- **Gate C passes:** F356W (Δf 0.0255 µJy, 16.4σ, shift 0.0007″) and F444W (0.0255 µJy, 11.1σ, shift 0.029″)
  are registered and `brighter`. F277W measures 7.7σ but its residual offset 0.044″ exceeds the 0.03″ rule, so it
  is `unregistered`.
- **Outcome U:** two SW bands are `unregistered` (F115W 0.170″, F150W 0.138″), which rule 6 maps to U.
- Per SW band (frozen classes): F115W `unregistered` (−0.15σ), F150W `unregistered` (+1.9σ), F200W `none`
  (+1.35σ; 5σ limit 0.0102 µJy against F356W's Δf 0.0255 µJy).

## Post hoc (not in the plan)
- **Failed approach:** the rule-1 registration check is a cross-correlation of a 6″ stamp with g09725 masked.
  The stamps hold almost no bright sources (contact sheet), so the correlation peak follows noise; the SW
  "offsets" of 0.14–0.17″ (4–6 SW pixels) are not credible as mosaic misregistration while F356W gives 0.0007″.
  A usable check needs a larger stamp or the mosaics' own WCS alignment; it was not re-run (no second look at the
  rule after seeing the data).
- Read as measured despite the registration flag, no SW band reaches 2σ (F115W −0.15σ, F150W +1.9σ, F200W
  +1.35σ), while all three LW bands do (7.7–16.4σ). Contact sheet viewed: the source is visible in both epochs in
  every band; the LW differences show a compact positive residual in the aperture, the SW differences show none
  standing out of the noise. So the change looks red-only at this depth (exploratory; not a frozen outcome).

## Reading
- No instrumental cause appears in the SW stamps (no dipole, edge or NaN in the aperture; `n_bad_pix` 0 in every
  band). The verdict of `docs/candidates/a2744_2561_flank_g09725.md` stays `inconclusive`; a red-only change does
  not by itself separate a supernova from AGN variability (both `hypothesis`).

## Deviations
- The first run executed with the driver uncommitted (a lint failure stopped the commit in the same command); the
  committed re-run (a variable rename and per-panel display scaling of the direct images only) wrote an identical
  `summary.json` (compared key by key).
