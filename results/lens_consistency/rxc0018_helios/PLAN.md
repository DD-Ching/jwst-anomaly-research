# RXC J0018.5+1626: SN Helios image order and next-image prediction in the pre-JWST RELICS Lenstool map (plan, frozen before the run)

**Input (observed, literature).** SN Helios (Fujimoto et al. 2026, arXiv:2609.30440v2, HTML read 2026-11-08):
spectroscopic Type II SN at z 3.34 behind RXC J0018.5+1626, discovered in VENUS NIRCam imaging of 2026-07-03 at
RA 00:18:31.50, Dec +16:26:18.0 (J2000); no counterpart in the 2024 December (SLICE) imaging; DD imaging on
2026-08-17. Only one image is reported. Published JWST-era models (their Table 5; none public) put the **next** image
2.32–6.0 yr after discovery, at positions that agree within a few arcsec: ALP 00:18:31.54 +16:26:00.3 (+2.32 yr),
zitrin-analytic 00:18:31.84 +16:26:04.3 (+2.70), glafic 00:18:31.78 +16:25:59.3 (+3.09), WSLAP+ 00:18:31.96
+16:25:57.0 (+4.18), LENSTOOL 00:18:31.73 +16:25:58.2 (+6.0). ALP and LENSTOOL also give a past image 46 / 42.7 yr ago.

**Question.** Does the public pre-JWST RELICS Lenstool v1 map (already pinned as `lensmodel.RXC0018_RELICS_LENSTOOL`,
the `arcs` primary) predict the same image order: no image that should have been visible next to Helios at
discovery, and a next image within years at the published place? An extra image predicted to be bright at
discovery would be a "multiple image where none is seen"; a next image far from the published region or decades away
is a model disagreement.

**Method** (reuses `scripts/macs0308_ares.py` / `smacs_at2024jfz.py`): trace the SN position to the source plane at
z_s 3.34, re-image with `lensmodel.find_images`, Fermat-potential delays with the straight vs L-shaped path check.

**Frozen rule** (`decide`, unit-tested before the run):
- gate P: path check ≤ 1 %, else `stopped: path`; gate S: an image within 0.3″ of the SN, else `stopped: no SN image`;
- `conflicts`: some other image arrives 0–1.0 yr (ASSUMPTION: still on the SN II plateau at discovery, ≤ 112 rest-frame
  days past maximum at z 3.34) **before** Helios with |μ| ≥ 0.5 |μ_Helios| (ASSUMPTION: ≤ 0.75 mag fainter, so it
  would have been seen in the discovery mosaic);
- `agrees`: no conflict, and the first later image arrives 1–10 yr (ASSUMPTION; published 2.32–6.0) after Helios
  within 5″ (ASSUMPTION; "a few arcsec") of the mean of the five published next-image positions;
- `differs`: no conflict, a later image exists, but the first one misses the delay window or the 5″ radius;
- `undetermined`: no other image.

**Ordinary explanations of a `conflicts` or `differs`.** The RELICS map predates JWST (HST constraints only), has no
z > 3 constraint and is a single best model; a different image order follows from mass-profile and member-halo
differences. Either outcome is a model disagreement, not an anomaly, unless the JWST-era models share it. A
`conflicts` would be checked against the VENUS discovery mosaic at the predicted position.

**Not done here.** No JWST-era model is public, so no map-level comparison; no MCMC spread; H0 70, Ωm 0.3 (map
cosmology). Post hoc (labelled): past-image delays vs the published 42.7–46 yr.
