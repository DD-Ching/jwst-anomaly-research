# `arcs`: public spectral coverage of the open anti rows (plan, frozen before any query)

**Question.** The two `robust` anti rows (RXC J0600 581, MACS0308 2157; D-024 addenda 13–14) can only be decided by
membership, and membership only by a spectrum (photo-z test `undetermined`, addendum 10; NED and DJA v4.6 have no
spec-z). Does any public spectroscopic observation outside NED / DJA cover either position?

**Hypothesis (H).** A released IFU or slitless observation (VLT/MUSE or other ESO IFU; JWST NIRISS / NIRCam WFSS;
HST grism) covers a target, so a spectrum can be extracted next cycle.
**Ordinary expectation.** RELICS clusters at z 0.36–0.43 have little public spectroscopy beyond the RELICS
programmes; none may exist.

**Route (metadata only, no data download).**
- ESO: `ivoa.ObsCore` TAP at archive.eso.org, `dataproduct_type` cube or spectrum, `INTERSECTS(s_region, POINT)`.
- MAST: `astroquery.mast.Observations.query_criteria` 3″ cone, `dataproduct_type` spectrum, then the footprint
  (`s_region`) must contain the point (`regions`-free polygon test on the STC-S string).
- WHL0137 2194 (explained, addendum 13) is queried too and reported, without a label.

**Classes per target (ASSUMPTION: footprint containment = coverage).**
- `spectrum_available`: ≥ 1 public (release date passed) IFU or slitless observation whose footprint contains the
  target.
- `msa_footprint_only`: only NIRSpec MSA / fixed-slit footprints contain it; a slit on the target needs the MSA
  metafile, so not decided here.
- `none`: neither.

**Gate Q (validation, before reading the targets).** Positive controls must be returned: Abell 2744 core
(RA 3.5862, Dec −30.4002; MUSE GTO mosaic, Mahler+18; GLASS NIRISS WFSS, Treu+22) → ESO ≥ 1 MUSE cube and MAST
≥ 1 NIRISS WFSS. Negative control: an empty-sky point (RA 90.5, Dec −35.0) → `none`. If gate Q fails, no target is
labelled (`undetermined`).

**What changes the plan.** `spectrum_available` → next cycle extracts the spectrum and measures z (decides
membership). `msa_footprint_only` → read the metafiles. `none` for both → the membership step is closed until new
public spectroscopy; the anti line rests on model tests only.
