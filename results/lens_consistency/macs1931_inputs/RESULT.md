# MACS J1931.8-2635 fitted-model route input check (2026-10-24)

Input check named in the Bullet Cluster handoff (CHANGELOG 2026-10-24): a public JWST-era image list, maps of a model
fitted to it, and a reference catalogue that deblends arcs from cluster light. No map value, model prediction or
catalogue flux is read here. Field chosen by a literature survey (research subagent, not stored): of PLCK G165.7+67.0,
MACS J0138, SPT-CL J0615, Sunburst, Abell 1689 and others, only MACS J1931 had public JWST-era deflection maps. Script
`scripts/macs1931_inputs.py`; numbers from `summary.json`; image list `m1931_images.ecsv`; tests
`tests/test_macs1931_inputs.py`.

| Input | Found | Status |
|---|---|---|
| Image list | VENUS AstroLensPy v2.1 `m1931_AstroLensPy_v2.1_arcs_WCS.txt` (Allingham+2026, arXiv:2602.14074; Zenodo 10.5281/zenodo.18646525; md5 checked against the record) | 54 rows; 51 used images in 19 systems, as the README states (`c`-prefixed c11.3, c18.1, c18.2 not used). 26 images in 9 spec-z systems (1–7, 100 = SN Eos host, 101 = SN Eos), one z each. 10 other systems: 5 carry different z per image (11, 13, 15, 16, 17), 7 images have z 0 (free) |
| Spec-z cross-check | the record's `m1931_redshifts_MUSE.txt` (50 entries) | **fails as run**: only 100 and 101 match (6 images, |Δz| ≤ 8e-5); for systems 1–7 the list holds no entry within 0.5″ at the same z. Post hoc (not stored): its lensed-z entries are few (z 2.195, 3.541, 3.999, 5.133); the README attributes systems 1–7 to Caminha+19 (arXiv:1903.05103), not checked here |
| JWST-era model maps | `Alpha_x.fits.zip`, `Alpha_y.fits.zip` (67.7 / 67.0 MB zipped; headers read from the first 200 kB by HTTP range) | 3001 × 3001 px, 0.04″ / px, RA---TAN / DEC--TAN, FK5, model `vA_52`; ±1σ and Std maps also public. **Unit conflict:** header `BUNIT = 'Pixels'`, README "in arcsec" (a factor 25 apart). The header comment says "model v2.0", the README v2.1. The source-z normalisation of α is not in the header (README default z 5.079, not stored) |
| Model fit quality | `chires.txt` | rms over the 51 used images: 0.429″ (Avg), 0.444″ (Bary; README 0.44″, not stored), 0.662″ (ByEye) |
| Reference catalogue | DJA v7 imaging index (35 fields listed) | **absent**: no `macs1931` / `j1931` field. The record has only a VENUS photo-z table (12.8 kB, not stored) of lensed images. MAST (archive.stsci.edu) did not answer from this session (connection failure), so a VENUS HLSP was not checked |

**Review fixes (/code-review, before the PR):** COMMENT cards joined without the 72-character wrap break, full DJA
field names matched, unsupported or encrypted zip members and an empty `chires` match refused, MUSE counts stored;
the verdicts are unchanged.

**Reading.** Two of three inputs exist: a JWST-era image list and the maps of the model fitted to it, so gate K can
run once a frozen plan fixes (a) the deflection unit, (b) the α normalisation redshift and (c) the source-redshift rule
for the photo-z / free systems. (a) is decidable before any residual is read only from the authors' description;
otherwise the plan states both readings and lets gate K on the 26 spec-z images decide (0.04″ × px vs arcsec differ
by 25×, so one reading fails grossly). Gate C2 and the extra-image test cannot run: no deblended reference catalogue
is public for this field. D-034 addendum aa.

**Data:** only the four small text files (≤ 4 kB, not stored) and 2 × 200 kB zip heads were read, in the session; nothing stored.
