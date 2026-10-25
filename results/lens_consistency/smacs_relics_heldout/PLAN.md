# RELICS HST-era SMACS 0723 models vs the JWST-era Mahler+22 images (held-out; plan, frozen)

Frozen 2026-10-25 (commit time of this file), before any map value is read. Follows D-034 addenda ap / at / au
(gate H on MACS0416, Abell 370, MACS1149). SMACS 0723 has no held-out test yet: Chow+2023 (`../smacs_chow23/`) and
ICLv2 (`../smacs_djaref/`) were fitted to the JWST images. Script `scripts/smacs_relics_heldout.py`; tests
`tests/test_smacs_relics_heldout.py` (written before the run).

## Inputs
- model_prediction: the three RELICS SMACS 0723 HLSP models (DOI 10.17909/T9SP45), all built before the JWST ERO
  release: `glafic/v2`, `zitrin-ltm-gauss/v2` (Golubchik+22; its readme: "JWST-blind") and `lenstool/v1` best map,
  arcsec deflections at D_LS/D_S = 1, z_lens 0.39, flat H0 70 Ωm 0.3, loaded as `smacs_s17_vet.relics` does
  (`MapLensModel.from_fits`). The lenstool sha256 is pinned at first download.
- observed: held-out list = Mahler+22 ICLv2 `arcs.dat` (pinned, `smacs0723-iclv2`), 21 systems, 17 of them first
  found in JWST ERO data (Mahler+22 abstract).
- ASSUMPTION (z): systems without a spectroscopic z take ICLv2 `best.par`'s optimised `z_m_limit` z (a Mahler+22
  model value, not a measurement); system 8 (11.76) is unconstrained (D-024). Reported per family.
- observed: seen list = Golubchik+22 Table 1 (`golubchik22_images.ecsv`, 17 rows, five HST-era families).
  ASSUMPTION for glafic v2 and lenstool v1: their own constraint lists are not public; Mahler+22 §5.6 says the
  RELICS-Lenstool model used seven systems, so a few images counted here as held-out may have been constraints of
  that model. A pass for those two is therefore optimistic; a fail stands.

## Question, ordinary explanations, what changes the plan
- Do HST-era (JWST-blind) models predict the JWST-found SMACS 0723 images? A pass makes that model a validated
  second model on SMACS 0723 for a later two-model extra-image stage against ICLv2. A fail is expected from ordinary
  causes first: model incompleteness away from the five HST families (the western diffuse component of Mahler+22),
  model-z assumption, map pixels, frame. No large residual is a flag in this stage.

## Frozen rules
Rules 1–5 of `../mars_hff_heldout/PLAN.md` with the same functions (`frame_offset`, `shift`, `split_seen`,
`heldout_residuals`, `gate_h`, `outcome` in `scripts/mars_hff_heldout.py`) and thresholds (offset pairs within
1.0″, ≥ 5 pairs; seen within 0.5″; gate H rms ≤ 1.0″ and ≥ 80 % within 3.0″; ≥ 5 tested images), per model, with
the Golubchik+22 list as the seen list (rule 1 moves the Mahler positions into the HST/RELICS frame). Grid step
0.25″.

## Readings fixed in advance
- zitrin-ltm-gauss v2 passes: a validated JWST-blind second model on SMACS 0723.
- glafic / lenstool pass: candidate second model, with the constraint-overlap caveat above.
- All fail with shared > 3″ images: likelier the families (z, identification) or a shared HST-era gap (the west).
- No flag, flux or extra-image reading; no image pixels read.

## Stop rules
- No rule change after the first residual; any departure is labelled deviation or post hoc.
- Maps (≈ 40–60 MB each) are fetched into the data root and deleted after use (cloud disk rule); only
  `summary.json` and per-model residual tables are stored.
