# RXC J0018.5+1626: SN Helios image order in the pre-JWST RELICS Lenstool map (2026-11-08)

Plan frozen in `PLAN.md` (claim commit) before the run; `decide` and its unit test committed before the run.
Script `scripts/rxc0018_helios.py`; numbers from `summary.json` (frozen rule) and `posthoc.json` (post hoc).
Labels: SN position and z (observed, Fujimoto et al. 2026, arXiv:2609.30440v2); published next-image positions,
delays and magnifications (model_prediction, their Table 5, models not public); RELICS image positions, μ and delays
(model_prediction, best model, no MCMC); thresholds are ASSUMPTIONs (PLAN.md).

**Result (frozen rule): `agrees`.** The SN position re-images onto itself (0.0″, μ +5.554); path check 3.3e-05.

| Image | Delay vs Helios (yr) | μ | From the published next-image centre (″) |
|---|---|---|---|
| at the SN position | 0 | +5.554 | 18.589 |
| past | −29.615 | +4.081 | 42.196 |
| first later | +6.753 | −2.201 | 2.384 |
| later | +9.963 | −1.61 | 12.854 |
| central | +12.046 | +0.001 | 6.615 |

No image arrives within 1 yr before Helios, so the map predicts no second image in the discovery mosaic, matching the
single reported image. The first later image arrives in 6.753 yr, 2.384″ from the mean of the five published next-image
positions (window 1–10 yr, radius 5″). The published JWST-era range is 2.32–6.0 yr; 6.753 yr exceeds the longest
(LENSTOOL 6.0) but lies inside that model's quoted upper error (paper Table 5; error not stored). Not a lens-model
residual.

**Post hoc** (`posthoc.json`; no verdict depends on it).
- Post hoc: the first later image lies 1.2″ from the published LENSTOOL next image and 2.198–6.679″ from the other four.
- Post hoc: the past image (29.61 yr ago, μ +4.081) lies 4.238″ from the published LENSTOOL past image (42.7 yr ago; ALP 46 yr).
  Position agrees within a few arcsec; the delay is shorter than the published 42.7 / 46 yr.
- **Amplitude sensitivity (post hoc; deviation: added after the run).** RXC J0018's Lenstool v1 map shares its map
  maker with the RELICS maps whose amplitude is offset (D-034 addenda bi–bk, which predicted a factor of about two
  thirds for offset maps; not stored here), and its amplitude has not been tested against spectroscopic families.
  Post hoc, with the deflection scaled by 0.85 the outcome
  is `differs` (5 images; first later image +10.708 yr, 19.966″ from the centre); at 0.70 Helios is singly imaged
  (`undetermined`, post hoc). So the `agrees` holds only at the map's nominal normalisation (D_LS/D_S = 1 as released). That the
  nominal map reproduces the JWST-era prediction is consistent with no offset for this map, but it is not a test of it.

**Contact sheet** `images.png` (looked at): RELICS images over the κ map; Helios NNE of the core, the first later
image among the published next-image positions (green), the past image to the SE.

**Deviations.** After the run and the review, `decide` was changed to count a non-finite μ ratio as bright
(conservative; all ratios in this run are finite, so the outcome is unchanged), and the case "only earlier images"
(not covered by the plan; returns `differs`) is now commented; it did not occur in the frozen run.

**Limits.** One public model only (no GLAFIC or Zitrin RELICS map for this cluster); no JWST-era map is public, so no
map-level comparison; best model, no MCMC; H0 70, Ωm 0.3 (map cosmology); amplitude untested (above). The maps (RELICS
HLSP, sha256-pinned in `lensmodel.RXC0018_RELICS_LENSTOOL`) were cached for the run and deleted after it; nothing
stored in git.
