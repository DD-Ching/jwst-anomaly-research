# MARS HFF free-form maps: input check C / U and gate K on six fields (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the deflection maps was read. Search and claim
note: no unused JWST-era cluster map record turned up in a Zenodo search on 2026-10-24 (queries and hits listed in
`search.json`); the one unused cluster-map record with image lists is HST-era: Cha & Jee 2023, MARS free-form models
of the six Hubble Frontier Fields (arXiv:2301.08765; Zenodo 10.5281/zenodo.7804575, CC-BY-4.0, one 2.7 MB zip, md5
checked). Script `scripts/mars_hff_inputs.py`; tests `tests/test_mars_hff_inputs.py`.

## Inputs (seen when this section was written: file listing, README, code, headers, image catalogue)
- README: `deflection_angle_*_{x,y}_w_header.fits` are 100 × 100 maps "in the unit of arc second", scaled to
  D_ds / D_s = 1; the catalogue has positions, z (`^{**}` = z fitted by MARS), class gold / silver and a per-image
  `Delta_rms`. Released code: flat ΛCDM Ωm 0.27, ΩΛ 0.73; deflection at catalogue positions is bicubic (Hermite)
  interpolation of the map; z_lens per field in `load_parameter_torch.py` (A370 0.375, A2744 0.308, AS1063 0.348,
  MACS0416 0.396, MACS0717 0.545, MACS1149 0.543).
- Headers: RA---TAN / DEC--TAN, no BUNIT, pixel 1.010″ (A2744), 1.515″ (A370, AS1063), 1.263″ (MACS0416), 2.020″
  (MACS0717), 0.758″ (MACS1149).
- Catalogue `MARS_HFF_Multiple_image_catalog.txt`: 887 rows, 672 gold.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the released MARS maps, loaded by our map-model code (bilinear), reproduce the model's own
  spectroscopic images? A pass gives a second, free-form (non-parametric) model on six fields for later stages,
  above all a held-out test on JWST-era images absent from the 2023 catalogue. A fail stops this record. Nothing
  exotic is read: residuals here are residuals of a fit to these very images; no flag is read.
- Ordinary explanations of a fail: our loading (unit, orientation, WCS), the coarse maps (0.76–2.02″ pixels) with
  bilinear instead of the authors' bicubic interpolation, the map extent, or the model.

## Frozen rules (computable; each a function with a unit test written before the run)
Per field, independently. Cosmology flat H0 70, Ωm 0.27 (authors' code). Maps are the deflection at D_ds / D_s = 1.
"Spec images" = gold rows whose z carries no `**`; a system enters only with ≥ 2 such images.
1. **C (coverage)** `coverage_field`: fraction of spec images inside the map. Passes when ≥ 0.9 (ASSUMPTION);
   images outside are dropped and listed. Fails → that field stops.
2. **U (unit)** for reading A (map value in arcsec) and P (map value in pixels × the header pixel): per system at its
   catalogue z, rms of `backtrace_images` `dtheta_arcsec`; `spt0546_inputs.unit_check` with threshold 0.5″
   (ASSUMPTION): passes when exactly one reading has median ≤ 0.5″. Fails → that field stops.
3. **K** `m0417_images.gate_k` on the spec images at catalogue z, mean source (`imageplane_residuals`, grid step
   0.25″): rms ≤ 1.0″ **and** ≥ 80 % within 3.0″. Reported, not a gate: K-WLS (`a2744_uncover.imageplane`).
4. **Record outcome** `record_outcome`: the record is usable for a held-out stage on the fields whose C, U and K all
   pass; "usable on N / 6 fields".
Reported (fixed now): per-system rms and image counts; the authors' median `Delta_rms` of the same images.

## Stop rules
- No rule changes after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No extra-image, flux or flag reading in this stage. Data: the zip lives in the session scratchpad and is deleted
  after the run; only `summary.json`, per-field image-plane tables and `search.json` are stored.
