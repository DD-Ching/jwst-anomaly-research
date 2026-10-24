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

## 2026-10-22: follow-up (07:20Z run): fine-grid validate and a fitted source

Frozen plan `results/lens_consistency/macs0416_canucsref/d115_refit/PLAN.md` (committed before the refit);
`scripts/m0416_d115_refit.py`, `d115_refit/result.json`; post hoc `d115_refit/posthoc_families.json`.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Grid step / plane z as cause of the D115a break | C | ruled out | `validate --grid-step 0.1`: χ²_pos 330.81 (fiducial), 330.73 (8515 at z 0.3907), 5262.01 (σ 84.4, z 0.396), 5261.44 (σ 84.4, z 0.3907); D115a.1–3 share one predicted image in both σ 84.4 runs (`validate_runs.json`) | model_prediction |
| Source definition as cause | C | **confirmed** | With the source fitted to the image plane (Nelder–Mead, one-to-one assignment), σ 84.4 at z 0.3907 gives 3 D115a images at 0.35 / 0.46 / 0.58″ (family χ² 2.77, against 4924.8 with the mean back-traced source) and 3 D115b images at 0.20 / 0.30″; no predicted image within 0.7″ of the flag (nearest 7.08″) | model_prediction |
| χ² cost of σ × 0.7 | C | **above the frozen bound** | χ²_fit 339.21 (σ 84.4) vs 326.62 (σ as fitted), Δχ² +12.6 > 6 (ASSUMPTION, system-51 precedent). Post hoc: +10.6 of it is in other families with `validate`'s source (D70 +5.6, D69 +4.6, D15 +4.0, D72b +2.7; D34 −8.3); nothing else in the model was refitted | model_prediction |

**Verdict (frozen rule):** still `inconclusive`. The flag's dependence on member 8515's mass is established: the
extra image disappears at σ × 0.7 with every catalogued D115 image matched. But that σ costs Δχ² +12.6 in a model
where only the two D115 sources were refitted, above the 6 accepted for system 51. Whether a re-fit (8515 σ free,
other parameters free, all sources fitted) absorbs the cost is not tested. Not a candidate for anything unusual.

### Open questions (updated)
- **Next:** a Lenstool re-fit with 8515's σ free (owner machine; Lenstool not in the cloud environment), or a
  σ × 0.75–0.85 scan here with every source fitted, to find the smallest σ change that removes the image and its cost.

## 2026-10-23: σ(8515) scan with every source fitted

Frozen plan `results/lens_consistency/macs0416_canucsref/sigma_scan/PLAN.md` (+ pre-run amendment);
`scripts/m0416_d115_sigma_scan.py`; `sigma_scan/result.json`, `sigma_scan/fits.jsonl` (555 / 555 tasks).

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Smallest σ factor that removes the flag image, all 111 sources fitted | C | f* = 0.75 (90.4 km/s) | flag image present at f 1.00 / 0.85 / 0.80, absent at 0.75 / 0.70; D115a / D115b and all 111 families matched one-to-one at every f; gates G1, G2, G3 pass | model_prediction |
| χ² cost at f* | C | **above the frozen bound** | Δχ² +12.28 > 6 (ASSUMPTION, system-51 precedent); +5.96 / +0.70 / +15.85 at f 0.85 / 0.80 / 0.70. Post hoc: the f 0.80 value is one family (D24, start-point gap); without it +8.88, monotonic | model_prediction |

**Verdict (frozen rule):** still `inconclusive` (H1). Removing the extra image by lowering 8515's σ alone costs
more than the accepted bound even with every source refitted; the cost is spread over D70, D69, D72b, D15 and D115a.
A model tension of a fixed-parameter Lenstool model, not evidence of anything unusual.

### Independent model, 2026-10-24 (`results/lens_consistency/m0416_mars_extra/RESULT.md`, D-034 addendum ar)
The MARS free-form maps (gate H passed on MACS0416) do not predict the flagged image: nearest free MARS D115b image
7.29″ away (frozen radius 3.0″, ASSUMPTION). Adds a model-specific reading; verdict unchanged (`inconclusive`), since
the 1.515″ MARS pixels are coarse at a galaxy-scale critical curve.

### Open questions (updated 2026-10-23)
- **Only next step:** a Lenstool re-fit with 8515's σ and the other parameters free (owner machine).

### Regeneration
Scratch scripts (not committed): `pert.py` (nearest potentials, ±30 % and removal for the three nearest, images of
D115a / D115b at z 2.2999 by `a2744_s700.solve` on a 0.1″ grid, half-width 58″) and `pert_scan.py` (8515 σ × 0.8–0.95).
Inputs: the pinned `macs0416-canucs` files (sha256 in `lensmodel.MACS0416_CANUCS`) and
`canucs_m0416_images.ecsv`. In the cloud, `SSL_CERT_FILE` must point to certifi + the proxy CA for archive.stsci.edu.
