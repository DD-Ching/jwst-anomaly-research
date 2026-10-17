# PLCKG165+67.0_2:h002 (VENUS batch, #296): vetting record

## Source
- Run: `results/transients/batch_venus` (box `PLCKG165+67.0_2`, screen hit `h002`).
- Program: `main` 7c975ac with `arc_diff.py` unchanged (#291, D-094, D-095); settings in
  `results/transients/batch_venus/PLAN.md` §3.
- Position: 171.805899, +42.481802 (11:27:13.42 +42:28:54.5), in PLCK G165.7+67.0.
- Pair: PEARLS 1176 o251 (2023-03-30) → VENUS 6882 o050 (2026-03-26), F150W + F444W.
- Score (`model_prediction`): matched-filter S of −35.9 (F150W) and −27.9 (F444W), a fading. It is rank 1 of 5 in
  the batch review queue.
- Vetted 2026-10-17. Queries: SIMBAD, NED and Gaia via astroquery, 2026-10-17.

## Tests
| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge or low coverage | A | does not explain | inside all four footprints; finite fraction 1.0; 1″-from-edge rule not triggered (`boxes.json`, `review.json`) | derived |
| Diffraction spike or bright-star neighbour | A | does not explain | brightest epoch-1 pixel within 5″ is 0.59 × (F150W) and 0.32 × (F444W) the box's 99.9th percentile (`post_review_checks.json`); no spike on the sheet | derived |
| Saturation or DQ | A | does not explain | 240 nJy in F150W, far below saturation; finite pixels | observed |
| Persistence | A | does not explain | detected in two programs (1176 on 2023-03-30, 4446 on 2023-04-22) and in 8 bands | observed |
| Cosmic ray or snowball (one band) | A | does not explain | present in F090W, F115W, F150W, F200W, F277W, F356W, F410M and F444W (S/N 14–41; `lightcurves.ecsv`) | observed |
| Registration or matching residual | A | does not explain | F150W and F444W S peaks 0.03″ apart (sky-fixed); a compact point-like negative lobe, not a dipole (`sheet_*h002.png`) | derived |
| Real flux change | — | yes | forced photometry with 0.15″ apertures and random-aperture errors: F150W 240 ± 7 → −6 ± 6 nJy (z −27), F444W 74 ± 4 → −2 ± 4 nJy (z −13); every band falls at z −14 to −28 between some pair of epochs | derived |
| Star (point-like, Gaia or SIMBAD) | C | does not explain | no Gaia source within 3″; no SIMBAD object within 3″; the source disappears completely | observed |
| Catalogued galaxy at the position | C | consistent with a host | NED `[FPP2024] NS_0481`, type G, z = 1.3236, 0.06″ away; reference 2024ApJ...961..171F (Frye et al. 2024). Redshift type and origin not verified (NED lists one "Uncertain origin" row). Whether the catalogued object is the host or the transient itself, measured in 2023, is not known. | observed |
| Light curve | C | explains, as a class | present on 2023-03-30 and still on 2023-04-22 (F090W 271 → 173 nJy, F150W 240 → 214, F200W 187 → 175); absent in 2025 (4744) and 2026 (VENUS). A transient lasting weeks to months that fades over two years is typical of a supernova. | derived; the class is a hypothesis |
| SN H0pe identification | C | inconclusive | SN H0pe's catalogued position (11:27:15.60 +42:28:33.7) is about 32″ away; image coordinates were not checked | observed |
| Transient-catalogue match (TNS) | C | not run | no TNS access in this session | — |

## Verdict
**known population: extragalactic transient, supernova-like (hypothesis).**
- **Basis:** a real, sky-fixed point source seen at two epochs in 2023 and gone in 2025–2026. It sits on a catalogued
  galaxy with z = 1.32 (Frye et al. 2024; redshift type not verified).
- **Not unexplained.** It is a fading that the search was calibrated to find; it is not wormhole evidence.
- **Open questions:**
  - TNS or literature identity: one of the transients reported in G165 by PEARLS or the SN H0pe follow-up?
  - Whether `NS_0481`'s redshift comes from the host or from the transient.

## Regenerate
- `scripts/batch_review.py` (sheets, `hit_lightcurve.py` forced photometry); `post_review_checks.json`
  (session script, stored output).
- Cross-matches: astroquery `Simbad.query_region` (3″, 30″), `Ned.query_region` (5″), `Ned.get_table(...,
  'references')`, `Gaia.cone_search_async` (3″), all on 2026-10-17.
