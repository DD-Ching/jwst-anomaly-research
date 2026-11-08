# MACS J0308.9+2645: SN Ares image order in the pre-JWST RELICS maps (2026-11-08)

Plan frozen in `PLAN.md` (2972398) before the run; `decide` and its unit test committed before the run.
Script `scripts/macs0308_ares.py`; numbers from `summary.json` (frozen rule) and `posthoc.json` (post hoc).
Labels: SN position and z (observed, Larison et al. 2026, arXiv:2610.33161v1); image positions, magnifications and
delays (model_prediction, RELICS best models, no MCMC); thresholds are ASSUMPTIONs (PLAN.md).

**Result (frozen rule).** Both public RELICS models give `agrees`:

| Model | n images | Ares image μ | Later arrivals (yr) | Path check |
|---|---|---|---|---|
| GLAFIC v2 | 3 | +6.434 | 57.21, 65.0 | 1.27e-05 |
| Zitrin LTM-Gauss v1 | 3 | +4.233 | 76.92, 85.44 | 3.51e-05 |

In both, the SN position re-images onto itself (0.0″) and is the first image to arrive; no image arrives earlier, so
neither model predicts an earlier appearance of the SN in archival imaging. The published JWST-era models give
72–79.1 yr (paper text; `published_later_arrivals_yr`) and μ 4.6–5.6 for Ares (`published_ares_mu`); the HST-only RELICS models bracket
those delays (57.21–85.44 yr). The image order is not a lens-model residual.

**Post hoc** (no verdict depends on it). The two later images lie near the cluster centre, as the paper states.
Against the planned future-image targets of JWST programme 9478 (`Ares_Image2`, `Ares_Image3`, from
`results/transients/macs0308_grid/PLAN.md`; not a refereed source): Zitrin 0.326″ / 3.389″, GLAFIC 1.525″ /
5.023″. The image nearest `Ares_Image2` is demagnified in both RELICS models (μ 0.359, 0.395), while the paper
quotes ~1.5× for one later image (`published_later_mu_approx`); which later image that is cannot be read from the text, and a central image's
magnification depends on the unconstrained core profile, so this is a model difference, not an anomaly.

**Contact sheet** `images.png` (looked at): predicted images over each model's κ map; Ares NE of the centre, the
two later images near the BCG.

**Limits.** None of the JWST-era models is public, so no map-level comparison. No MCMC spread; H0 70, Ωm 0.3
(map cosmology). The deflection maps (RELICS HLSP, sha256-pinned in `lensmodel.MACS0308_RELICS_*`) were read from
the data root cache; nothing stored in git.
