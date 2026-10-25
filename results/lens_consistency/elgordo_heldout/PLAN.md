# El Gordo: Caminha+2023 (pre-JWST) model vs JWST-era Diego+2023 images (gate H) — frozen plan

Frozen 2026-10-25 before any prediction for this list is computed. Script `scripts/elgordo_heldout.py`.

## Hypothesis and ordinary explanations
- H: the Caminha+2023 dPIE model (`elgordo-caminha23`, D-030; fitted to 56 HST/MUSE images of 23 systems) predicts
  images it was not fitted to well enough to serve as a held-out-validated model on El Gordo (D-034 addendum ap).
- Ordinary reasons for a failure: model precision (Caminha's own rms 0.754″), knot-to-knot position differences
  inside one arc, wrong image identifications in the newer list, and wrong redshifts (geometric z).
- Plan change: a pass makes Caminha+2023 a held-out-validated model on El Gordo (a second model for an extra-image
  stage needs a JWST-era model, Diego+2023 maps are not public in this repository). A fail closes this route; no image
  is flagged by this stage either way.

## Inputs (observed unless stated)
- Held-out candidates: Diego et al. 2023, A&A 672, A3 (arXiv:2210.06514), Table `tab_arcs` (`TableArcs_NEW3.tex`
  in the arXiv source): ID, RA, Dec, z_spec, z_geo. Families 1–23 follow Caminha's naming and carry Caminha's MUSE
  z_spec on their first row (observed). Families ≥ 24 use z_spec when given, otherwise z_geo, a model_prediction of
  Diego's own model (ASSUMPTION; not independent of a lens model).
- Seen list: the model's own 56 catalogued images (`lens_consistency.image_list`).

## Rules (addendum ap rules 1–5 unchanged, via `m1149_hst_heldout.model_block`)
1. Frame offset: median offset of pairs within 1.0″ (ASSUMPTION), applied to the Diego list.
2. An image is seen when a Caminha image lies within 0.5″ (ASSUMPTION); all others are held out.
3. Per family with a held-out image: source from the seen images (type E, ≥ 2 seen) or from all images (N / M).
   Family = `lensmodel.image_family` of the Diego ID (knots `1.1`, `1.2` are separate families).
4. Gate H (ASSUMPTION): rms ≤ 1.0″ and ≥ 80 % within 3.0″ on tested held-out images inside the map;
   at least 5 tested images (rule 5), else `insufficient`.
5. **Primary gate:** held-out images in families with spectroscopic z (z_spec, MUSE). **Reported only:** the same
   numbers over all held-out images (geometric-z families included).

## Readings fixed in advance
- Images > 3″ in the primary set are listed; post hoc checks (if any) are labelled post hoc.
