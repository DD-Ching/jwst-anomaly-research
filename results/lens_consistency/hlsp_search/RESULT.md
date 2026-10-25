# MAST HLSP search for JWST-era cluster lens models (2026-10-25)

**Question** (TASKS Supporting, "a MAST HLSP / paper-data search for JWST-era maps or new releases"): does MAST
host a cluster lens-model HLSP that this repository has not used yet? The Zenodo half was done on 2026-10-24
(`../mars_hff_inputs/search.json`).

**Method:** `hlsp_scan.py` queries `obs_collection = HLSP` within 0.05° of 18 lensing-cluster centres (observed:
MAST listings). Unfamiliar provenances were opened with `get_product_list` and their file names checked for
model / kappa / deflection / magnification products. Seven guessed `archive.stsci.edu/hlsps/<name>/` directories
(GLIMPSE, PEARLS, VENUS, UNCOVER, ALCS, Technicolor, MegaScience) return 404. `search.json` holds the script output (`provenance_names_by_field`; all 18 queries
succeeded) and the hand-made product checks.

**Result:** the only lens-model HLSPs at these fields are BUFFALO, CANUCS, RELICS, CLASH and Frontier, all already
used. JWST-ERO (SMACS 0723) holds four MIRI mosaics; GLASS-JWST (Abell 2744) mosaics and NIRSpec spectra;
TEQUILA / EVEREST are TESS / K2 light curves that fall in the cone. **No new lens-model HLSP.**

**Limits:** cone search by cluster centre only (a model registered under offset coordinates would be missed);
GLASS-JWST product names were checked for its first 50 observations only; no HLSP index page was readable (the
MAST HLSP table is JavaScript-loaded). Do not repeat this search before a new MAST HLSP announcement.
