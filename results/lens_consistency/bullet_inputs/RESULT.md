# Bullet Cluster reference-route input check (2026-10-24)

Input check named in the MACS0647 handoff (CHANGELOG 2026-10-23): a public JWST-era image list, maps of a model
fitted to it, and a reference catalogue that deblends arcs from cluster light. No model prediction, map value or
catalogue flux is read here. Script `scripts/bullet_inputs.py`; numbers from `summary.json`; image list
`bullet_gold_images.ecsv` (observed); tests `tests/test_bullet_inputs.py`.

| Input | Found | Status |
|---|---|---|
| Image list | CANUCS Bullet DR1 `lensing/bullet_gold_v4.dat` (Rihtaršič+2026, arXiv:2601.22245; the Lenstool constraints) | 135 images, 49 systems, 27 galaxies, 27 distinct z (0.897–6.64), as the paper states; the graded full list adds 22 silver and 42 bronze images, and its gold rows equal the Lenstool list (135 / 135 ids, max separation 0.0″, max Δz 0) |
| JWST-era model maps | `lensing/best_fit_model/best_fit_maps/` (Lenstool fiducial, scaled to D_LS / D_S = 1 per the README) | deflection `dpl1` / `dpl2`, κ, γ, γ1, γ2 at 50 mas (gzip, 205–221 MB each, not read); ψ at 200 mas; 200 mas Bayesian sample maps in `sample_model/sample_maps/` (13.8 MB each; the first listing page has 35 entries, so the listing may be paginated) |
| Reference catalogue | `catalogs/hlsp_canucs_jwst-hst_multi_bullet-clu_multi_v1_photometry-cat.fits` (51,721,920 B; header only, HTTP range) | 20,138 rows, 323 columns; every column of the frozen CANUCS DR1 mapping (`../macs0417_canucsref/PLAN.md`) present, `FLUX_KRON_TOTAL_F200W` in nJy: **passes** |

**Reading.** All three inputs of the fitted-model route exist for the Bullet Cluster, in the same CANUCS
format as MACS0416 / MACS1149, where the route already ran (gate K, gate C2). Unlike MACS0647 (D-034 addendum v),
the model was fitted to the JWST image list, so gate K applies as written. Next, under a frozen plan: gate K with
the 50 mas best-fit deflection maps (two downloads of ≈ 0.2 GB each, over the 200 MB per-file rule: the reason goes
into the plan and DECISIONS.md), or the 200 mas sample maps if the plan judges that resolution sufficient; then
gate C2 and the extra-image test with this catalogue. Not checked here: whether the image list and the catalogue
share a frame (frame check F of the plan).

**Failed approach (data handling).** The paper-page Dropbox link serves the whole lens-model folder as one
24.5 GB zip (sizes from the transfer log, not stored) and ignored an HTTP range request: about 5.0 GB reached the session
scratchpad within 60 s before the time limit cut it. Only its local-header index (file names, sizes) and the
README were read; the partial file was deleted and nothing from it is used here. Use the CANFAR vault (per-file URLs, ranges honoured), never the Dropbox folder link.
