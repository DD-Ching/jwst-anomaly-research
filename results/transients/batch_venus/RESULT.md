# Bounded VENUS batch search with fixed `arc_diff`: result

Plan: `PLAN.md`, committed before box selection. `boxes.json` was committed before any difference image. Owner
instruction 2026-10-17; budget 06:01Z → 10:01Z. Used 06:01Z → about 07:00Z, including two review rounds (session
clock, not stored). Every number is `derived` unless labelled;
injected sources are simulations. An ordinary fading is not wormhole evidence; zero entries is not an exclusion.

## What was searched (new data)
- **Data:** 20 boxes of 24″ × 24″, 3.2 arcmin² in all, in 10 cluster fields. Each pairs an earlier public NIRCam
  epoch with VENUS (6882) in F150W + F444W (`boxes.json`).
- **Boxes per field:** two each, at the cluster core and its west neighbour (east for El Gordo):
  - Abell 370 (1208 o002, 2022-12 → 2026-01);
  - Abell S1063 (3293, 2024-09 → 2025-09);
  - Bullet (4598, 2025-01 → 2026-01);
  - El Gordo (1176 o241, 2022-07 → 2026-08; two new boxes beside B1–B3);
  - MACS1423 (1208 o010, 2023-02 → 2026-03);
  - MACS0416 (1176 o212, 2022-12 → 2026-01);
  - MACS0417 (1208 o006, 2022-10 → 2025-10);
  - MACS1149 (1176 o231, 2024-01 → 2026-05);
  - MS0451 (5058, 2025-02 → 2025-09);
  - PLCK G165 (1176 o251, 2023-03 → 2026-03).
- **Pipeline versions:** jwst 3.0.0 in every epoch, except Bullet and MS0451 epoch 1 (2.0.1).
- **Selection:** coverage, finite fraction ≥ 0.95 in all four images, and depth only.
  - Every pair's shallower F150W exposure is the same VENUS 1,245 s, so the depth ranking ties. The stable sort then
    took pairs in alphabetical order of target. This tie-break is labelled **post hoc**; it uses no difference
    information.
  - RXCJ0018, RXJ2129, SPT0615, SMACS and Sunrise fell beyond the 20-box cap and were not searched.
- **Already searched and excluded:** El Gordo B1–B3, the Sunrise box and the SMACS system-5 box (3 grid boxes
  dropped).

## Calibration (per box, fixed rule; `boxes/batch.json`, `boxes/*/inj_eval.json`)
- **12 boxes passed:** vanish 26.0 + 27.0 AB on the review list with sign ≥ 0.70, and ≤ 2 unchanged-control false
  hits. That is 1.92 arcmin².
  - A370 ×2, Bullet ×2, El Gordo 1, MACS1423 2, MACS0416 2, MACS0417 2, MACS1149 ×2, MS0451 2, PLCK G165 2.
- **8 boxes are undeterminable:**
  - too many control false hits: AS1063 ×2 (4/10 and 12/20), El Gordo 2 (3/20), MACS0416 1 (4/17);
  - too few vanish recoveries: MACS1423 1 (3/8), MACS0417 1 (3/7), MS0451 1 (4/10), PLCK G165 1 (9/14).
  - They are not reported as searched with zero candidates.
- **Sample size:** controls per box range 9–20 (fewer positions where the box has few bright pixels), so single-box
  pass decisions rest on small numbers.
- **Pooled over the 12 passing boxes** (review list with sign; 68 % Wilson intervals; `totals.json`):

  | Change | 26.0 AB | 27.0 AB | 27.5 AB |
  |---|---|---|---|
  | vanish (all flux gone) | 62/71 = 0.87 [0.83, 0.91] | 56/68 = 0.82 [0.77, 0.87] | 50/63 = 0.79 [0.74, 0.84] |
  | partial fade (50 % drop) | 54/68 = 0.79 [0.74, 0.84] | 47/63 = 0.75 [0.69, 0.80] | 37/71 = 0.52 [0.46, 0.58] |
  | appear (brightening) | 48/63 = 0.76 | 59/73 = 0.81 | 55/70 = 0.79 |
  | unchanged controls, any hit | 8/207 = 0.04 [0.03, 0.05] | | |

- **Where these efficiencies apply:**
  - point sources injected on galaxy and arc pixels;
  - not within 0.8″ of bright cores, which are flagged and not reviewed;
  - not extended changes;
  - not fainter than 27.5 AB.
  - The unchanged-control fraction is not a real-sky false-alarm rate.

## Screen and review (`boxes/review_queue.json`, `boxes/review/`)
- **Flags:** 113 hits in the 20 screens; 13 not near a bright core; **5 review-list entries.**
- **After the post-review checks:** 1 ordinary transient (h002), 3 with no significant change, 1 undetermined
  (h003). 0 pending, 0 unexplained.
  - 3 are in passing boxes, 2 in undeterminable ones.
  - All 5 were reviewed in detail (capacity 10). Pending: 0.
- **Each detailed review:** the contact sheet, looked at; forced photometry in every covering public mosaic (180
  measurements); the coded class rule.

| Rank | Entry | Box status | z F150W / F444W (search epochs) | Class | Note |
|---|---|---|---|---|---|
| 1 | PLCK G165 2, h002 (171.805899, +42.481802) | passed | −27.1 / −13.4 | **real change** | point source in all 8 PEARLS bands on 2023-03-30 (S/N 14–41) and on 2023-04-22 (programme 4446, F090W 271 → 173 nJy, F150W 240 → 214 nJy); absent in 2025 (4744) and 2026 (VENUS) |
| 2 | MACS0416 2, h003 | passed | −0.35 / −2.87 | **undetermined** (post-review) | F444W ridge along a bright elongated arc; band peaks 0.49″ apart; across programmes ≥ 3σ differences of opposite sign (F277W +5.6, F356W −9.1, F444W −7.3) |
| 3 | MACS1423 2, h006 | passed | −0.89 / −1.03 | no significant change | weak S on a compact source |
| 4 | MACS0416 1, h004 | undeterminable | −1.10 / −2.82 | no significant change | F444W ridge along a thin arc |
| 5 | AS1063 1, h005 | undeterminable | −0.78 / −1.56 | no significant change | weak S on a compact source |

- **h002 is vetted** (`docs/candidates/batch_venus_PLCKG165_h002.md`): known population, an extragalactic
  transient, supernova-like (hypothesis).
  - It is sky-fixed: the band peaks are 0.03″ apart.
  - It is not a spike, persistence, star or cosmic ray: no bright source within 5″, 8 bands, 2 programmes, no
  Gaia or SIMBAD match within 3″.
  - The nearest catalogued galaxy with a redshift, NED `[FPP2024] NS_0481` (spectroscopic z = 1.3236, Frye et al.
  2024), is 3.77″ away (corrected 2026-10-17: an earlier "on the galaxy, 0.06″" read NED's arcmin column as arcsec).
- **It is not an unexplained change:**
  - Its multi-epoch light curve is that of a transient: present over at least 23 days in 2023 and fading in the blue
    bands, then gone by 2025. That is the light curve of a supernova or similar ordinary transient.
  - It lies about 1.2″ from a compact galaxy; no host is visible directly under it.
  - It is about 32″ from the published SN H0pe position in this cluster (11:27:15.60 +42:28:33.7). Whether it is one
    of SN H0pe's images or another transient is **not established**. The per-image coordinates were not found, and
    no transient-catalogue cross-match was done.
  - The frozen rule sent it to `/vet-candidate`, and the vetting found an ordinary explanation (a known population,
    as a hypothesis).
  - It is a fading recovered by the search, as the calibration says the search can. It is not wormhole evidence.
### Revision 2026-10-17 (follow-up, `followup/`; the text above is kept as first written)
- **h002 status: real fade; an ordinary transient is a reasonable explanation; its identity and type are not
  confirmed** (`followup/h002_identity_followup.json`).
  - **Reported transient:** no match found in the sources that could be queried (Frye et al. 2024 full text, NED,
    SIMBAD, Gaia, web search). TNS could not be queried (HTTP 403 twice), so "no TNS record" is not established.
  - **SN H0pe:** the three published image positions (Frye et al. 2024 Table 2) are 24.9″, 31.9″ and 37.9″ away; the
    lens model's image-position RMS is 0.65″. This replaces "about 32″ ... per-image coordinates not found" above.
    Lens-model positions of any further image were not checked.
  - **NS_0481** (3.77″, spectroscopic z = 1.3236 from a DDT 4446 MSA filler slit) is a separate galaxy; the spectrum
    is not of the transient, and a host relation is not established.
  - **Appearance:** HST WFC3/IR F160W in 2016 (GO-14223) gives 4.8 ± 9.5 nJy (3σ < 33 nJy, 0.2″, no aperture
    correction; `followup/h002_hst2016.json`) against F150W 240 ± 7 nJy in 2023: it appeared after 2016.
  - **"Fading in the blue bands" narrowed:** over 23.8 d only F090W declines significantly (271 → 173 nJy, z −4.0);
    F150W (z −1.7) and F200W (z −1.3) are consistent with constant.
  - **Missing for a type:** a spectrum of the transient, a host and redshift, a classification report, and light-curve
    sampling (two detected epochs, no rise or peak). "Supernova-like" stays a hypothesis.
  - Identity unknown does not mean first discovery, and this is not an exotic candidate.
- **h003 status: the cross-programme differences are explained by between-mosaic photometric scatter of matched arc
  pixels (one pre-fixed control test); the specific mechanism is not isolated** (`followup/h003_recheck.json`).
  - The large differences are same band, same fixed sky position and same 0.15″ aperture; they occur inside
    programme 1176 (o211 → o213, F356W 359.8 → 308.3 nJy) as well as between programmes, and the mosaics differ in
    orientation.
  - 31 arc-pixel controls (`pick-arc`, parameters fixed in `followup/PLAN.md` before the run) on the decisive pair
    F356W 1176 o211 → 1208 o004: calibrated z −0.76, control tail fraction 0.39 (rule: ≥ 0.05 ⇒ explained). The other
    four pairs give tail fractions 0.26–0.52.
  - The raw z differs from #296 (−17.7 vs −9.1) because the random-aperture error is measured in a larger box when
    controls are grouped; the fluxes are identical.
  - This changes h003 from "undetermined" to "explained by measurement scatter" under the fixed rule; the owner's
    instruction expected "undetermined" unless a concrete check supported a processing effect (flagged in the PR).
- **Review queue vs open science questions:** the review queue is empty (5 of 5 reviewed, 0 pending). Detailed review
  finished does not mean every object is explained: h002's identity and physical type remain open. No evidence for
  special physics was found; that does not mean every ordinary explanation was proven.
- **Not reviewed (by design):** 108 hits that are near bright cores or fail morphology. They are saved in each
  box's `screen/hits.ecsv` and `morph_all.json`, not deleted.

## Deviations and corrections (labelled)
- **Post hoc tie-breaks:**
  - Across targets, equal depths were taken in input order (alphabetical).
  - Within a target, the pair chosen among equal depths followed MAST row order. The code now takes the lowest
    obs_id pair; `boxes.json` is not re-run.
  - Grid ties at equal distance go to the lowest RA, which gives the west neighbour.
- **The first review used a narrower coded class rule** than PLAN §5: only the two search epochs; no sky-fixed,
  spike or persistence test.
  - The post-review checks (`boxes/review/post_review_checks.json`) apply the full rule: all covering mosaics,
    sky-fixed offset and bright neighbour.
  - This moves h003 from "no significant change" to "undetermined".
- **PA_V3:** `boxes.json` recorded PA_V3 from the primary header, where level-3 files lack it. The SCI-extension
  values are in `headers_sci.json`.
- **Query pattern:** box selection queried MAST once per target and read boxes serially, not in batches (an
  efficiency limitation only).

## Time used
- **Box selection** (MAST metadata, footprints, finite fractions, headers): about 5 min.
- **Computation:**
  - 40 `arc_diff` runs plus samples and morphology, 3,849 task-seconds;
  - 22 min of wall time on 4 workers (06:10–06:32Z);
  - peak memory 1.53 GB per task.
- **Review:**
  - 2.4 min of computation (sheets and forced photometry);
  - looking at the 5 sheets and the h002 stamps;
  - the literature and catalogue checks and the h002 vetting;
  - about 20 min in all, including the second review round (not stored).
- **Records:** about 15 min (not stored).

## Pending
- Identify h002 in TNS (403 on 2026-10-17). Literature: Frye et al. 2024 reports only SN H0pe; NS_0481's z is
  spectroscopic (NIRSpec) and 3.77″ away (host not determined).
- h003 (MACS0416 box 2) stays undetermined: it needs a PSF-matched multi-programme check on the arc. *Superseded 2026-10-17: explained by matched arc-pixel scatter (see Revision).*
- Boxes not searched in this batch: the undeterminable 8, which need a decision on whether to re-calibrate with more
  injections, and pairs beyond the cap.
- No follow-on batch is started.

## Regenerate
- `python scripts/batch_pairs.py --pairs ../sfft_trial/venus_pairs.json --exclude … --out boxes.json`
  (the exclusions are listed in `boxes.json` `rule`).
- `python scripts/batch_run.py run --boxes boxes.json --out boxes --scratch <dir> --workers 4`, then
  `batch_run.py totals` (writes `totals.json`).
- `python scripts/batch_review.py --batch boxes --boxes boxes.json --scratch <dir> --capacity 10
  --sheet-verdicts boxes/review/sheet_verdicts.json`.
