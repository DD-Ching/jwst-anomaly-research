# El Gordo: Caminha+2023 vs the JWST-era Diego+2023 images — gate H fails (2026-10-25)

Plan `PLAN.md` (frozen in 64f923ef; parser, primary-gate rule and tests committed before the run). Script
`scripts/elgordo_heldout.py`; tests `tests/test_elgordo_heldout.py`. Numbers from `summary.json` (plan) and
`posthoc.json` (post hoc). Positions and z_spec are observed; z_geo and all residuals are model_prediction;
thresholds are ASSUMPTIONs (addendum ap).

## Result (frozen rules)
- Input: Diego et al. 2023 Table `tab_arcs`, 218 image rows; 23 families carry a MUSE z_spec.
- Frame offset dx −0.235″, dy +0.003″ (69 pairs). All 56 Caminha images are matched (seen); 162 held out
  (10 E, 150 N, 2 M), none outside the map.
- **Primary (spec-z families): gate H fails on rms**: 46 images, rms **1.64″**, 93 % within 3″. Images > 3″: **6c,
  13c, 20c**.
- All held-out images (geometric z included, reported only): rms 8.98″, 88 % within 3″.
- No held-out-validated second model on El Gordo. No image is flagged by this stage (plan).

## Post hoc (not in the plan)
- 6c, 13c and 20c are among the 7 tested `\dag` rows: third-image candidates first proposed from the JWST images.
  Those 7 give rms 3.55″, 57 % within 3″; the other 39 primary images give 0.961″, 100 % within 3″. By anchor type,
  E (source from seen images only, 10) 3.11″ and N (source includes the image itself, 36) 0.88″.
- Caminha's own predicted images of systems 6, 13 and 20 (`elgordo_djaref/predicted_elgordo.ecsv`) have m_pred
  28.5–29.2 and were classed `undetectable` there.
- Contact sheet (F200W, 6″ stamps, 0.5″ circles; session script, not stored): faint compact sources sit at Diego's
  6c, 13c and 20c; Caminha's predicted positions for 6 (`no_flux_ref`), 13 and 20 show no visible source, as
  expected at m_pred ≈ 29; Caminha's 6 `candidate` position falls on a bright galaxy.

## Reading
Caminha+2023 reproduces the JWST-era knots and images of its own MUSE families except the newly proposed third
images. Two ordinary explanations remain and are not separated here: the candidates are other faint galaxies
(misidentifications), or the model's counter-image positions, far from its constraints, are off by 3–6″. Spectra
of 6c / 13c / 20c (or a JWST-era model with public maps) would decide. No exotic reading is considered.

## Failed approaches
- First run stopped in the parser: z columns carry LaTeX footnote marks (`2.64\ddag`, `3.1\P`); stripped, tested and
  committed before any prediction was computed.
