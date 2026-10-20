# W3 dimming: `forced` + footprint/depth test on the 20 MACS0416 / Abell 2744 `on_sky` persistence targets

Frozen 2026-10-20 03:10 UTC, before any forced flux was read (handoff 2 of #386, still open after #388;
`results/dimming/persistence_run/RESULT.md` "Next" 2). Same method as `results/dimming/nexus_o014/PLAN.md`;
screen code, thresholds and `forced` unchanged from `main` (c4ee84a).

- **Targets.** The 20 compact single-epoch vanish flags whose D-039 persistence verdict was `on_sky`
  (`persistence_run/{macs0416,abell2744}/verdicts.ecsv`, `kind == flag`): MACS0416 12, all catalogue-detected only
  in epoch 2 (`jw01208-o004`, CANUCS); Abell 2744 8, only in epoch 1 (`jw02561-o002`).
- **Screen.** Rebuilt with the production chain as on `main`, i.e. with the weight-map coverage term on (D-052
  addendum h; the persistence run used the screen before it). Catalogues from S3 by sha256
  (`persistence_run/fetch_cats.py`), Gaia mask from the ESA TAP (`persistence_scope/gaia_esa.py`). Targets are
  matched by uid (amended 03:14 UTC, before any forced flux: the verdict tables carry no positions, so the
  check is that each uid is catalogue-detected in F200W in its detection epoch only). A target that is no longer a flag under the term is
  explained by `weight_edge` (the committed term removes it) and is not measured.
- **Hypotheses.** (H0) Each is a real source present in every epoch that the other epochs' pipeline catalogues
  miss (footprint, depth or deblending), i.e. not a vanish. (H1) It is fainter by the screen's vanish criterion in
  the other epochs in forced photometry at a fixed position.
- **Ordinary explanations tested here:** (a) footprint / image quality; (b) catalogue depth / deblending. Not tested
  here: SN / AGN / variable star, which go to `/vet-candidate` for survivors.
- **Measurement.** `dimming_screen.py forced` per field with the targets as the only flags (every other flag marked
  `catalogue_veto` in a copy of `lightcurves.ecsv`; the targets' `catalogue_veto` cleared because D-039 `on_sky`
  replaces `single_epoch` for them) and the default 300 matched controls.
- **Decision rule** (the `nexus_o014/summarise.py` rule with the detection epoch per field): survivor iff
  `survives_artefact_tests`; else `catalogue_depth_or_deblend` if `forced_vanish` is False and some other F200W
  epoch has a finite forced flux; else `footprint_or_image_quality` if no other epoch is finite or a cutout flag is
  set; else `other_forced_test`. Survivors go to `/vet-candidate` before any D-039 swap.
- **Calibration requirement.** Per field, `calibrated` true and `controls_flagged` ≤ 1 % of controls (ASSUMPTION;
  the committed MACS0416 and Abell 2744 `forced` runs had 0 / 300). A field failing it does not count as a
  calibrated null; its footprint readings (coverage only) still stand.
- **What would change the plan:** any survivor (vetting first, the swap waits); a calibration failure (record and
  stop for that field). If all 20 are explained, the next step is the swap in `inject` (#388 handoff 1).
