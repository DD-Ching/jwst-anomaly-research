# SMACS MIRI colours from DJA fixed apertures: plan (frozen 2026-11-07 15:20Z, before measuring)

TASKS Next M2 "MIRI: evaluate fixed-aperture colours". D-013 moved NIRCam colours to DJA 0.5″ fixed apertures because
the pipeline encircled-energy (EE) apertures have different angular sizes per band. MIRI colours in the ranking still
come from the pipeline `aper50` EE apertures, whose radii grow ~2.3× from F770W to F1800W.

**Data.** DJA v7.4 SMACS `fix_phot.fits` (sha256 361a4eb9…, the D-013 catalogue) carries SEP fixed-aperture fluxes in
F770W, F1000W, F1500W, F1800W (2686–2719 finite rows per band). Aperture 2 (0.7″ diameter) is used, chosen before
measuring because 0.5″ is smaller than the F1800W PSF FWHM (≈ 0.59″, ASSUMPTION from the JDox PSF table).
Pipeline run: `configs/smacs_jwst3_miri_dja.yaml` = `smacs_jwst3_miri_confirm.yaml` plus a `matched_photometry`
block (aperture 2, label `dja07`, 0.2″ join), classify / cutouts / cross-match as that config.

**Hypothesis H.** The pipeline MIRI colours carry a source-size bias that fixed apertures remove.
**Ordinary explanations.** A real colour–size trend in the galaxy population (e.g. low-z large galaxies with PAH
emission in F770W); aperture-correction differences between bands for point sources (constant, harmless for ranking).

**Measurements** (`scripts/feature_size_bias.py --compare dja07`, common rows, F770W aper50 S/N > 10):
F770W−F1000W and F1000W−F1500W; Spearman ρ between log F770W isophotal area and colour for `aper50`, `isophotal`,
`dja07`.

**Frozen rules (ASSUMPTIONs).**
- V (join and units): among the common rows in the lowest quartile of F770W isophotal area, the MAD of
  (`dja07` colour − `aper50` colour) for F770W−F1000W ≤ 0.15 mag. If V fails, stop: no adoption.
- C (coverage): ≥ 0.80 of the ranked (gate-passing) MIRI sources get a finite `dja07` F770W−F1000W colour.
- B (bias): F770W−F1000W |ρ_dja07| ≤ |ρ_aper50| − 0.10, or |ρ_aper50| < 0.10 (then no bias to remove: H rejected).
- Adopt `dja07` as the MIRI colour aperture only if V, C and B all pass; otherwise not adopted, ranking unchanged.
- Report the top-20 Jaccard against the confirmation run (`smacs_jwst3_miri_confirm.yaml`) either way; look at the
  contact sheet of new entries before writing any reading.

**What would change the plan.** B passing with C failing → a coverage fallback is needed (pipeline colours where DJA
has no row), queued, not built here.
