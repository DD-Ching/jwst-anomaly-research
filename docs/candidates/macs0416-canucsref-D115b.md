# MACS0416 D115b: predicted extra image not seen (CANUCS DR1 reference test)

## 2026-10-22: vetting record

- **What was flagged.** `results/lens_consistency/macs0416_canucsref/` (PR #458, frozen plan `PLAN.md`): the pinned
  CANUCS Lenstool model (`macs0416-canucs`, D-044) predicts an extra image of D115b (second knot of source D115,
  z 2.2999) at 64.048551, −24.064825, |μ| 26.6, m_pred 24.45. Rule 8 reads it `absent` (max F / F_pred 0.049 over
  3 CANUCS DR1 sources). D115a has three catalogued images and D115b two; the model predicts six images of each knot.
- **source_uid / run.** `macs0416-canucs:D115b-extra` (no catalogue source); the #458 run at e24031c.
- **Score.** None: a lens-model consistency flag. Predictions are `model_prediction`; photometry `observed` /
  `derived`.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge, detector artefacts | A | does not explain | Inside the CANUCS F200W mosaic; first-look cutout (05:40Z run, not stored) shows BCG-subtraction residuals by a cluster galaxy, no edge or spike | observed |
| Bright neighbour light, deblending | A/B | inconclusive, not needed | Cluster members 9087 and 8515 lie 1.11″ and 1.34″ away; the companion D115a extra prediction 0.6″ away is `recovered` at F / F_pred 129, i.e. on member light | derived |
| Fold position (critical curve) | C | consistent | Fiducial: both D115a and D115b have an extra pair at |μ| 15–21 within 0.5″ of the flag (0.1″ grid, scratch `pert.py`); the image sits on a galaxy-scale critical curve | model_prediction |
| Which potential makes it | C | **8515** | σ × 0.7 or removal of 9087 (σ 96.6 km/s) moves the D115b image 0.5–1.1″ but keeps it. σ × 0.7 (84.4 km/s) or removal of 8515 (σ 120.6 km/s, r_cut 5.5″) leaves 3 images per knot, none within 7″ of the flag; σ × 1.3 adds a seventh. Potential 8908 (σ 32.7) changes nothing. Scan of 8515 (`pert_scan.py`): the D115a extra pair is gone at σ × ≤ 0.9; the D115b extra image stays within 0.7″ of the flag down to × 0.8 (|μ| 15–18) and is gone at × 0.7 | model_prediction |
| Catalogued images without the extra image | C | still matched | 8515 σ × 0.7: D115a 0.26 / 0.51 / 0.59″, D115b 0.15 / 0.33″ (fiducial 0.16–0.38″; model rms 0.512″) | model_prediction |
| Member redshift (two-plane) | C | members | NED (positional, 1″, 2026-10-22): 8515 = CLASH-VLT J041611.69−240354.4, z 0.3907 (SUN); 9087 = CLASH-VLT J041611.69−240352.5, z 0.4015 (SUN); cluster z 0.396. No foreground/background deflector | observed |
| χ² cost of the lower mass | C | **inconclusive** | `lens_consistency.py --model macs0416-canucs validate [--plane 8515=0.3907[:84.4]]` (0.25″ solver grid, 303 images): χ²_pos 330.81 (as fitted) → 330.73 (8515 on its z_spec plane, σ as fitted) → 5261 (σ 84.4; rms 2.04″, max 25.4″), the increase from D115a.1–3 all sharing one predicted image. That disagrees with the 0.1″ single-plane solve above (D115a matched within 0.59″); the cause (grid step, plane z 0.3907 vs 0.396, or a real loss) is not isolated | model_prediction |
| Scaling-relation scatter | — | — | ±30 % in σ for an individual member (ASSUMPTION, as D-042 and the system-51 record) | assumption |

**Verdict:** `inconclusive: needs the χ² cost of a lower 8515 σ on a fine solver grid`. Leading ordinary reading
(`hypothesis`): lens-model systematics, i.e. the galaxy-scale mass of cluster member 8515 (z_spec 0.3907), as in the
MACS0416 system-51 record (`macs0416-system51-fourth-image.md`). The extra D115a / D115b images exist only while 8515's
σ is above × 0.7–0.8 of its fitted value (within the assumed ±30 % scatter), and on a 0.1″ grid the catalogued images
stay matched without them. But the 0.25″ `validate` run at σ × 0.7 loses D115a's three-image match (χ² 330.8 → 5261),
so the cost of the lower mass is not established. Not a candidate for anything unusual; not closed either.

### Open questions
- **Next (cheapest):** `validate --plane 8515=0.396:84.4 --grid-step 0.1` and the same at 8515=0.3907, i.e. the cost
  at σ × 0.7 on the fine grid, single-plane and on the z_spec plane. If χ² stays within the system-51 precedent
  (Δχ² ≲ +6), the verdict becomes the known-population reading above; if D115a still breaks, scan σ × 0.75–0.85.
- A full re-fit with 8515's σ free would test the lower σ against all constraints; here the rest of the model is
  fixed.
- Rule for later fields (extends D-042's rule): before reading a high-|μ| extra prediction next to a cluster member
  as a flag, run the ±30 % σ test on the nearest potentials.

### Regeneration
Scratch scripts (not committed): `pert.py` (nearest potentials, ±30 % and removal for the three nearest, images of
D115a / D115b at z 2.2999 by `a2744_s700.solve` on a 0.1″ grid, half-width 58″) and `pert_scan.py` (8515 σ × 0.8–0.95).
Inputs: the pinned `macs0416-canucs` files (sha256 in `lensmodel.MACS0416_CANUCS`) and
`canucs_m0416_images.ecsv`. In the cloud, `SSL_CERT_FILE` must point to certifi + the proxy CA for archive.stsci.edu.
