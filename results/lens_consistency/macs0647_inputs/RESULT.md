# MACS0647 reference-route input check (2026-10-23)

Input check named in the RXJ2129 handoff (CHANGELOG 2026-10-23): a public JWST-era image list, model maps, and a
reference catalogue that deblends arcs from cluster light. No model prediction or residual is read here. Script
`scripts/macs0647_inputs.py`; numbers from `summary.json`; tables `meena23_images.ecsv`, `nishida25_images.ecsv`
(observed, parsed from the arXiv TeX sources; the parsers refuse any row or redshift cell they cannot read; tests
`tests/test_macs0647_inputs.py`).

| Input | Found | Status |
|---|---|---|
| Image list 1 | Meena+2023 (arXiv:2211.13334) Table `tab:lensed` | 72 images in 23 systems (25 families: MACS0647-JD's three components each a family); photometric z only (13 rows printed `--`); 3 ids marked `?` (8.3, 16.3, 16.4); systems 13–23 marked "New system" (`new_in_jwst`) |
| Image list 2 | Nishida+2025 (arXiv:2505.00553) appendix tables, the JWST GLAFIC constraints | 86 images from 28 sources, as the paper states; z printed for 26 families (14 with ±, 12 without; the paper does not say what a value without ± means), none for families 4 and 22 |
| Frame of the two lists | nearest-neighbour match | 66 / 86 Nishida images within 0.3″ of a Meena image (64 distinct Meena images; 69 within 1″). Over the 66 pairs within 0.5″ only, median offset +0.12″ RA, +0.02″ Dec, so where the lists share images they agree to about 0.1″; the other 17 Nishida images have no Meena image within 1″ (not compared) |
| JWST-era model maps | none public found | Meena+2023 (dPIEeNFW, Glafic) and Nishida+2025 (GLAFIC; the paper quotes χ²/DoF 136/95, not stored) give no map or data link in their arXiv sources; a Zenodo search for MACS0647 finds no model record |
| HST-era model maps | CLASH Zitrin v2 PIEMD+eNFW and LTM-Gauss (`lensmodel.MACS0647_CLASH_NFW` / `_LTM`, pinned 2026-10-16) | fitted to the pre-JWST systems in the CLASH frame; the MAST directory was not reachable from this session (curl TLS error), the pinned files were not re-read |
| Reference catalogue | DJA v7.0 `macs0647-grizli-v7.0-fix_phot.fits` (headers only, HTTP range; v7.1–v7.5 returned 403 in the session, not stored) | SMACS rule 2 **passes**: 7,956 rows, 455 columns, nothing missing, `APER_1` 12.5 px |

**Review fixes (/code-review, before commit):** the Nishida z carry now follows the last printed family, not the
last kept row; unreadable rows or z cells raise; outputs are written only after the network check; the hash covers
every column; the offset is stated over its matched subset.

**Reading.** The route as run on Abell 2744, El Gordo, SMACS 0723 and RXJ2129 needs a model fitted to the JWST
image list; for MACS0647 no such model's maps are public, so that route cannot run here. What can run is a
different test: the HST-era CLASH maps against the JWST lists, where the systems first identified in JWST data
(Meena+2023 systems 13–23, `new_in_jwst`) were not among the model's constraints, so their positions are a held-out prediction.
That needs, before any read, a frozen plan with a redshift rule (no list gives spectroscopic redshifts as such),
a CLASH-to-JWST frame check, and gate K's thresholds stated for a model not fitted to the JWST-first systems.
D-034 addendum v.

**Data:** the two arXiv source archives (≈ 8.8 MB and 4.4 MB, not stored) were read in the session scratchpad and deleted
(the first came in two range reads after the transfer was cut at 8 MB); only the DJA headers (first blocks) were
read. Nothing else was downloaded.
