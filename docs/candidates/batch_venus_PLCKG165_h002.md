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
| Catalogued galaxy near the position | C | no catalogued host at the position | NED `[FPP2024] NS_0481` (type G) is **3.77″ away**, not 0.06″ (corrected 2026-10-17: NED's separation column is in arcmin, 0.063′). Its z = 1.3236 ± 0.0002 is spectroscopic (Frye et al. 2024, arXiv:2309.07326, Table 2, NIRSpec MSA ID NS 481 from DDT 4446; grating not determined). At 3.8″ (about 31 kpc at z = 1.32, not stored) it is a weak host association; no host is seen directly under the transient (F150W −6 ± 6 nJy in 2026 in the 0.15″ aperture). | observed |
| Light curve | C | explains, as a class | present on 2023-03-30 and still on 2023-04-22 (F090W 271 → 173 nJy, F150W 240 → 214, F200W 187 → 175); absent in 2025 (4744) and 2026 (VENUS). A transient lasting weeks to months that fades over two years is typical of a supernova. | derived; the class is a hypothesis |
| SN H0pe identification | C | does not explain | the three SN H0pe images (Frye et al. 2024 Table 2, NS 2/3/4) are 24.9″, 31.9″ and 37.9″ away (subagent check 2026-10-17, not stored); no image within 1″ | observed |
| Transient-catalogue match (TNS) | C | not run | wis-tns.org returned HTTP 403 (its own load balancer, not the session proxy) on 2026-10-17; Frye et al. 2024 discusses only SN H0pe; two web searches found no other G165 transient report | — |

## Verdict
**known population: extragalactic transient, supernova-like (hypothesis).**
- **Basis:** a real, sky-fixed point source seen at two epochs in 2023 and gone in 2025–2026. No catalogued galaxy lies
  under it; the nearest catalogued galaxy with a redshift (NS_0481, spectroscopic z = 1.3236) is 3.77″ away, so
  the host and redshift are not determined (an earlier version said "0.06″", a unit error).
- **Not unexplained.** It is a fading that the search was calibrated to find; it is not wormhole evidence.
- **Open questions:**
  - TNS identity (not reachable 2026-10-17); no literature report found (Frye et al. 2024 only lists SN H0pe).
  - The host: none catalogued at the position; NS_0481 (3.77″) is a weak association.

## Revision 2026-10-17 (identity follow-up; the record above is kept)
- **Status: real fade; an ordinary transient is a reasonable explanation; identity and physical type not confirmed.**
- Reported transient: none matched in the sources that could be queried; TNS could not be queried (HTTP 403 ×2).
- SN H0pe images: 24.9″, 31.9″, 37.9″ (Frye et al. 2024 Table 2; now stored in
  `results/transients/batch_venus/followup/h002_identity_followup.json`); further lens-model images not checked.
- NS_0481: separate galaxy at 3.77″; its z = 1.3236 is the galaxy's own NIRSpec redshift; not shown to be the host.
- HST F160W 2016: 4.8 ± 9.5 nJy (3σ 28.5 nJy), so the source most likely appeared after 2016
  (indicative: different instrument, band and aperture; `followup/h002_hst2016.json`).
- The 23.8-day decline is significant only in F090W (z −4.0; F150W −1.7, F200W −1.3).
- Missing for a type: transient spectrum, host and redshift, classification report, better light-curve sampling.
- Identity unknown ≠ first discovery; not an exotic candidate.

## Regenerate
- `scripts/batch_review.py` (sheets, `hit_lightcurve.py` forced photometry); `post_review_checks.json`
  (session script, stored output).
- Cross-matches: astroquery `Simbad.query_region` (3″, 30″), `Ned.query_region` (5″), `Ned.get_table(...,
  'references')`, `Gaia.cone_search_async` (3″), all on 2026-10-17.
